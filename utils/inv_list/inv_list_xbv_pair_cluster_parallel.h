/** @file inv_list_xbv_pair_cluster_parallel.h
    @brief Cross-bvector (xbv) research support for inv_list.
    @note Development-only; included by inv_list.cpp after shared options and I/O.
*/
/*
Copyright(c) 2002-2026 Anatoliy Kuznetsov(anatoliy_kuznetsov at yahoo.com)

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0
*/

// One BM thread-pool task per RS cluster. Included after inv_list_xbv_pair_shard.h.
#ifndef INV_LIST_XBV_PAIR_CLUSTER_PARALLEL_H
#define INV_LIST_XBV_PAIR_CLUSTER_PARALLEL_H

namespace
{

struct pair_cluster_job
{
    std::atomic<bm::id64_t> progress{0};
    std::vector<pair_match> matches;
    bm::id64_t category_pairs[3] = {};
    bm::id64_t actual_memory = 0;
    double seconds = 0;
    std::string error;
};

void analyze_pairs_cluster_parallel()
{
    if (bv_in_file.empty() || pair_csv.empty() || pair_index_in.empty() || rs_cluster_in.empty())
        throw std::runtime_error("Parallel cluster search requires -bvin, -csv, -pair-index-in and -pair-cluster-in");
    if (pair_max_ratio != 0 || pair_inspect)
        throw std::runtime_error("Parallel cluster search requires -pair-max-ratio 0 and does not support -pair-inspect");
    if (!pair_cluster_memory)
        throw std::runtime_error("-pair-cluster-memory must be positive");

    const std::string best_name = pair_csv + ".best.csv";
    const std::string report_name = pair_csv + ".report.txt";
    for (const auto& name : std::vector<std::string>{pair_csv, best_name, report_name})
    {
        std::ifstream f(name.c_str(), std::ios::binary);
        if (f) throw std::runtime_error("Output already exists: " + name);
    }
    bool source_complete = false;
    bm::id64_t indexed_bytes = 0;
    auto index = pair_read_shard_index(source_complete, indexed_bytes);
    const std::size_t indexed_vectors = index.size();
    auto entries = pair_apply_rs_clusters(index, indexed_bytes);
    std::ifstream source(bv_in_file.c_str(), std::ios::binary | std::ios::ate);
    if (!source || source.tellg() < std::streampos(0) ||
        bm::id64_t(source.tellg()) < indexed_bytes ||
        (source_complete && bm::id64_t(source.tellg()) != indexed_bytes))
        throw std::runtime_error("Pair index input size does not match collection");
    source.close();

    std::vector<std::pair<std::size_t, std::size_t> > clusters;
    std::vector<bm::id64_t> estimates, pair_counts;
    bm::id64_t planned_pairs = 0;
    for (std::size_t begin = 0; begin < entries.size(); )
    {
        std::size_t end = begin + 1;
        bm::id64_t memory = entries[begin].memory;
        while (end < entries.size() && entries[end].shard == entries[begin].shard)
        {
            if (entries[end].memory > (std::numeric_limits<bm::id64_t>::max)() - memory)
                throw std::runtime_error("Cluster memory estimate overflow");
            memory += entries[end++].memory;
        }
        if (memory > pair_cluster_memory)
            throw std::runtime_error("Cluster " + std::to_string(clusters.size() + 1) +
                                     " exceeds -pair-cluster-memory; increase the budget or replan clusters");
        clusters.emplace_back(begin, end);
        estimates.push_back(memory);
        pair_counts.push_back(pair_choose2(end - begin));
        planned_pairs += pair_counts.back();
        begin = end;
    }
    const bm::id64_t all_pairs = pair_choose2(entries.size());

    std::ofstream csv(pair_csv.c_str()), best_csv(best_name.c_str()), report(report_name.c_str());
    csv.exceptions(std::ios::failbit | std::ios::badbit);
    best_csv.exceptions(std::ios::failbit | std::ios::badbit);
    report.exceptions(std::ios::failbit | std::ios::badbit);
    report << "Input: " << bv_in_file << "\nPair shard index: " << pair_index_in
           << "\nRS cluster manifest: " << rs_cluster_in
           << "\nIndexed vectors: " << indexed_vectors << "; searched vectors: " << entries.size()
           << "\nStatus: incomplete\nSerializer level: " << c_level
           << "\nPair evaluation threads: " << pair_threads
           << "\nExecution: one BM thread-pool task per cluster"
           << "\nActive cluster memory budget bytes: " << pair_cluster_memory
           << "\nMinimum structural gain: " << pair_min_gain
           << "; minimum target gain percent: " << pair_min_pct << '\n' << std::flush;
    const char* size_columns = ",lead_bytes,target_bytes,smaller_bytes,larger_bytes,size_ratio";
    csv << "lead_id,target_id,lead_cost,target_cost,xor_cost,gain_cost,gain_pct,xor_count"
        << size_columns << '\n' << std::setprecision(12);
    best_csv << "lead_id,target_id,lead_cost,target_cost,xor_cost,gain_cost,gain_pct,xor_count"
             << size_columns << '\n' << std::setprecision(12);

    std::vector<pair_best> best(indexed_vectors);
    std::vector<bool> participants(indexed_vectors, false);
    bm::id64_t completed = 0, qualified = 0, targets = 0;
    bm::id64_t category_pairs[3] = {}, category_qualified[3] = {};
    bm::id64_t bins[5] = {}, ratio_candidates[4] = {};
    bm::id64_t reserved = 0, peak_reserved = 0, max_actual = 0;
    unsigned active = 0;
    const auto start = pair_clock::now();
    auto last = start;
    typedef bm::thread_pool<bm::task_descr*, bm::spin_lock<bm::pad0_struct> > pool_type;
    typedef bm::task_batch<pair_bv::allocator_type> batch_type;
    pool_type pool;
    pool.start(pair_threads);
    batch_type batch;
    batch.get_task_vector().reserve(clusters.size());
    std::vector<std::unique_ptr<pair_cluster_job> > jobs;
    jobs.reserve(clusters.size());
    for (std::size_t cluster_no = 0; cluster_no < clusters.size(); ++cluster_no)
    {
        jobs.emplace_back(new pair_cluster_job);
        pair_cluster_job* job = jobs.back().get();
        const std::size_t begin = clusters[cluster_no].first;
        const std::size_t end = clusters[cluster_no].second;
        bm::task_function_t task([job, &entries, begin, end](void*)
        {
            try
            {
                const auto start_job = pair_clock::now();
                std::ifstream input(bv_in_file.c_str(), std::ios::binary);
                if (!input) throw std::runtime_error("Cannot open pair input in cluster worker");
                bm::streams_deserializer<pair_bv> reader;
                bm::serializer<pair_bv> serializer;
                serializer.set_compression_level(c_level);
                bm::serializer<pair_bv>::buffer buffer;
                std::vector<std::unique_ptr<pair_item> > items;
                items.reserve(end - begin);
                for (std::size_t p = begin; p < end; ++p)
                {
                    auto item = pair_load_indexed_item(input, entries[p], reader, serializer, buffer);
                    pair_bv::statistics st;
                    item->bv.calc_stat(&st);
                    job->actual_memory += st.memory_used;
                    items.push_back(std::move(item));
                }
                bm::id64_t done = 0;
                for (std::size_t i = 0; i < items.size(); ++i)
                    for (std::size_t j = i + 1; j < items.size(); ++j)
                    {
                        const auto& ei = entries[begin + i];
                        const auto& ej = entries[begin + j];
                        ++job->category_pairs[unsigned(ei.super_sparse) + unsigned(ej.super_sparse)];
                        const bool first_lead = items[i]->cost < items[j]->cost ||
                            (items[i]->cost == items[j]->cost && ei.id < ej.id);
                        const std::size_t target = first_lead ? j : i;
                        const auto residual = bm::calc_complexity_xor(items[i]->bv, items[j]->bv);
                        const auto cost = pair_cost(residual), baseline = items[target]->cost;
                        if (cost < baseline)
                        {
                            const auto gain = baseline - cost;
                            const double pct = 100.0 * double(gain) / double(baseline);
                            if (gain >= pair_min_gain && pct >= pair_min_pct)
                                job->matches.push_back(pair_match{i, j, cost, residual.count});
                        }
                        if ((++done & 255u) == 0)
                            job->progress.store(done, std::memory_order_relaxed);
                    }
                job->progress.store(done, std::memory_order_relaxed);
                job->seconds = std::chrono::duration<double>(pair_clock::now() - start_job).count();
                return 0;
            }
            catch (const std::exception& e) { job->error = e.what(); return 1; }
            catch (...) { job->error = "Unknown cluster worker failure"; return 1; }
        });
        batch.add(task, 0);
    }
    std::vector<std::size_t> admission_order(clusters.size());
    for (std::size_t i = 0; i < admission_order.size(); ++i) admission_order[i] = i;
    std::stable_sort(admission_order.begin(), admission_order.end(),
                     [&](std::size_t a, std::size_t b) { return estimates[a] > estimates[b]; });
    std::vector<unsigned char> state(clusters.size(), 0); // pending, running, done, emitted
    std::size_t emitted = 0;
    std::string failure;
    struct pool_guard
    {
        pool_type& pool;
        bool joined = false;
        ~pool_guard() { if (!joined) join(); }
        void join()
        {
            pool.set_stop_mode(pool_type::stop_when_done);
            pool.join();
            joined = true;
        }
    } guard{pool};

    auto progress = [&](bool force)
    {
        const auto now = pair_clock::now();
        if (!force && now - last < std::chrono::seconds(10)) return;
        bm::id64_t observed = completed;
        for (std::size_t i = 0; i < jobs.size(); ++i)
            if (state[i] == 1) observed += jobs[i]->progress.load(std::memory_order_relaxed);
        const double seconds = std::chrono::duration<double>(now - start).count();
        const double rate = seconds > 0 ? double(observed) / seconds : 0;
        if (!is_silent)
            cout << "Clusters " << emitted << '/' << clusters.size() << "; pairs " << observed
                 << '/' << planned_pairs << " (" << (planned_pairs ? 100.0 * double(observed) / double(planned_pairs) : 100.0)
                 << "%), pairs/s=" << rate << ", elapsed_s=" << seconds
                 << ", ETA_s=" << (rate > 0 ? double(planned_pairs - observed) / rate : 0)
                 << ", active=" << active << ", reserved_MiB=" << double(reserved) / (1024.0 * 1024)
                 << ", qualifying=" << qualified << ", targets=" << targets << endl;
        last = now;
    };

    while (emitted < clusters.size())
    {
        bool changed = false;
        for (std::size_t i = 0; i < clusters.size(); ++i)
        {
            if (state[i] != 1) continue;
            bm::task_descr* task = batch.get_task(batch_type::size_type(i));
            if (!task->done.load(std::memory_order_acquire)) continue;
            state[i] = 2;
            --active;
            reserved -= estimates[i];
            completed += pair_counts[i];
            max_actual = (std::max)(max_actual, jobs[i]->actual_memory);
            for (unsigned k = 0; k < 3; ++k) category_pairs[k] += jobs[i]->category_pairs[k];
            if (task->err_code || !jobs[i]->error.empty())
                failure = "Cluster " + std::to_string(i + 1) + " failed: " + jobs[i]->error;
            changed = true;
        }
        while (emitted < clusters.size() && state[emitted] == 2)
        {
            const std::size_t begin = clusters[emitted].first;
            const std::size_t match_count = jobs[emitted]->matches.size();
            for (const auto& match : jobs[emitted]->matches)
            {
                const auto& ei = entries[begin + match.first];
                const auto& ej = entries[begin + match.second];
                const bool first_lead = ei.cost < ej.cost || (ei.cost == ej.cost && ei.id < ej.id);
                const auto& le = first_lead ? ei : ej;
                const auto& te = first_lead ? ej : ei;
                const auto gain = te.cost - match.residual_cost;
                const double pct = 100.0 * double(gain) / double(te.cost);
                ++qualified;
                ++category_qualified[unsigned(le.super_sparse) + unsigned(te.super_sparse)];
                participants[le.id] = participants[te.id] = true;
                csv << le.id << ',' << te.id << ',' << le.cost << ',' << te.cost << ','
                    << match.residual_cost << ',' << gain << ',' << pct << ',' << match.residual_count;
                pair_write_sizes(csv, le.bytes, te.bytes);
                csv << '\n';
                ++ratio_candidates[pair_ratio_bucket(pair_size_ratio(le.bytes, te.bytes))];
                auto& b = best[te.id];
                if (!b.found) ++targets;
                if (!b.found || gain > b.gain || (gain == b.gain && le.id < b.lead))
                { b.found = true; b.lead = le.id; b.residual = match.residual_cost;
                  b.gain = gain; b.count = match.residual_count; }
            }
            jobs[emitted]->matches.clear();
            report << "Cluster " << emitted + 1 << ": pairs=" << pair_counts[emitted]
                   << "; matches=" << match_count << "; actual decoded MiB="
                   << double(jobs[emitted]->actual_memory) / (1024.0 * 1024)
                   << "; worker seconds=" << jobs[emitted]->seconds << '\n';
            state[emitted] = 3;
            ++emitted;
            changed = true;
        }
        if (failure.empty())
        {
            for (const auto i : admission_order)
            {
                if (active == pair_threads) break;
                if (state[i] != 0 || estimates[i] > pair_cluster_memory - reserved) continue;
                pool.get_job_queue().push(batch.get_task(batch_type::size_type(i)));
                state[i] = 1;
                ++active;
                reserved += estimates[i];
                peak_reserved = (std::max)(peak_reserved, reserved);
                changed = true;
            }
        }
        progress(false);
        if (!failure.empty() && !active) break;
        if (!changed) std::this_thread::sleep_for(std::chrono::milliseconds(5));
    }
    guard.join();
    if (!failure.empty()) throw std::runtime_error(failure);
    progress(true);

    std::vector<const pair_index_entry*> by_id(indexed_vectors, nullptr);
    for (const auto& entry : entries) by_id[entry.id] = &entry;
    for (std::size_t target = 0; target < best.size(); ++target)
    {
        const auto& b = best[target];
        if (!b.found) continue;
        const auto& le = *by_id[b.lead];
        const auto& te = *by_id[target];
        const double pct = 100.0 * double(b.gain) / double(te.cost);
        ++bins[pct < 25 ? 0 : pct < 50 ? 1 : pct < 75 ? 2 : pct < 100 ? 3 : 4];
        best_csv << b.lead << ',' << target << ',' << le.cost << ',' << te.cost << ','
                 << b.residual << ',' << b.gain << ',' << pct << ',' << b.count;
        pair_write_sizes(best_csv, le.bytes, te.bytes);
        best_csv << '\n';
    }
    csv.close();
    best_csv.close();
    auto summary = [&](std::ostream& out)
    {
        out << "Sharded search complete for " << (source_complete ? "whole collection" : "loaded prefix")
            << ". Indexed vectors=" << indexed_vectors << "; searched vectors=" << entries.size()
            << "; shards=" << clusters.size()
            << "\nAll searched-vector pairs=" << all_pairs << "; internal-shard pairs=" << planned_pairs
            << "; deliberately omitted cross-shard pairs=" << all_pairs - planned_pairs
            << "; coverage pct=" << (all_pairs ? 100.0 * double(planned_pairs) / double(all_pairs) : 100.0)
            << "\nPairs considered=" << completed << "; completed=" << completed
            << "\nPair evaluation threads=" << pair_threads
            << "\nExecution: one BM thread-pool task per cluster"
            << "\nActive cluster memory budget bytes=" << pair_cluster_memory
            << "; peak reserved estimate MiB=" << double(peak_reserved) / (1024.0 * 1024)
            << "\nQualifying pairs=" << qualified << "; candidate targets=" << targets
            << " (" << (entries.empty() ? 0 : 100.0 * double(targets) / double(entries.size()))
            << "%); participating vectors=" << std::count(participants.begin(), participants.end(), true)
            << "\nPair categories normal-normal,mixed,sparse-sparse: evaluated="
            << category_pairs[0] << ',' << category_pairs[1] << ',' << category_pairs[2]
            << "; qualifying=" << category_qualified[0] << ',' << category_qualified[1] << ',' << category_qualified[2]
            << "\nBest-target gain distribution [0,25), [25,50), [50,75), [75,100), 100: "
            << bins[0] << ',' << bins[1] << ',' << bins[2] << ',' << bins[3] << ',' << bins[4]
            << "\nPeak loaded shard allocation MiB=" << double(max_actual) / (1024.0 * 1024)
            << "\nStructural gains are heuristic units, not serialized bytes. Cross-shard outcomes are unknown.\n"
            << "Size ratio buckets for qualifying pairs [1,2),[2,4),[4,8),[8,infinity): "
            << ratio_candidates[0] << ',' << ratio_candidates[1] << ',' << ratio_candidates[2] << ',' << ratio_candidates[3]
            << "\nElapsed seconds=" << std::chrono::duration<double>(pair_clock::now() - start).count() << '\n';
    };
    summary(cout);
    summary(report);
    report << "Status: complete\n";
}

} // namespace

#endif
