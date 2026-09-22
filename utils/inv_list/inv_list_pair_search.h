/*
Copyright(c) 2002-2026 Anatoliy Kuznetsov(anatoliy_kuznetsov at yahoo.com)

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

For more information please visit:  http://bitmagic.io
*/

// Pairwise collection analysis for inv_list.cpp.
// Included after the utility options; not a public library header.
#ifndef INV_LIST_PAIR_SEARCH_H
#define INV_LIST_PAIR_SEARCH_H

namespace
{
using pair_bv = bm::bvector<>;
using pair_clock = std::chrono::steady_clock;

// Same structural heuristic as sample28; these units are not serialized bytes.
bm::id64_t pair_cost(const bm::bvector_complexity_statistics& s)
{
    return s.count ? s.bit_blocks * 8192u +
        s.gap_words * sizeof(bm::gap_word_t) + s.zero_spans + s.full_spans : 0;
}
struct pair_item
{
    pair_bv bv;
    bm::id64_t cost = 0, bytes = 0;
};
struct pair_best
{
    bool found = false;
    std::size_t lead = 0;
    bm::id64_t residual = 0, gain = 0, count = 0;
};

struct pair_coord
{
    std::size_t first = 0, second = 0;
};

struct pair_match
{
    std::size_t first = 0, second = 0;
    bm::id64_t residual_cost = 0, residual_count = 0;
};

struct pair_job
{
    std::vector<pair_coord> pairs;
    std::vector<pair_match> matches;
    bm::id64_t considered_end = 0;
};

// Ratio is larger/smaller; serialized empty vectors still have a BLOB header.
double pair_size_ratio(bm::id64_t a, bm::id64_t b)
{
    return double((std::max)(a, b)) / double((std::min)(a, b));
}
unsigned pair_ratio_bucket(double ratio)
{
    return ratio < 2 ? 0u : ratio < 4 ? 1u : ratio < 8 ? 2u : 3u;
}
void pair_write_sizes(std::ostream& out, bm::id64_t lead, bm::id64_t target)
{
    out << ',' << lead << ',' << target << ',' << (std::min)(lead, target)
        << ',' << (std::max)(lead, target) << ',' << pair_size_ratio(lead, target);
}
struct pair_measured_best
{
    bool found = false;
    std::size_t lead = 0;
    bm::id64_t xor_bytes = 0, gain = 0;
};

void analyze_pairs()
{
    if (bv_in_file.empty() || pair_csv.empty())
        throw std::runtime_error("Pair search requires -bvin FILE and -csv FILE");
    // Never truncate an existing output (including an input alias).
    const std::string best_name = pair_csv + ".best.csv";
    const std::string report_name = pair_csv + ".report.txt";
    const std::string measured_name = pair_csv + ".measured.csv";
    std::vector<std::string> output_names = {pair_csv, best_name, report_name};
    if (pair_inspect) output_names.push_back(measured_name);
    for (const auto& name : output_names)
    {
        std::ifstream existing(name.c_str(), std::ios::binary);
        if (existing.good())
            throw std::runtime_error("Output already exists: " + name);
    }
    std::ifstream input(bv_in_file.c_str(), std::ios::binary);
    if (!input) throw std::runtime_error("Cannot open pair input");
    input.seekg(0, std::ios::end);
    const auto file_size = input.tellg();
    input.seekg(0);
    if (!input || file_size <= std::streampos(0))
        throw std::runtime_error("Empty or inaccessible pair input");
    std::ofstream csv(pair_csv.c_str()), best_csv(best_name.c_str()), report(report_name.c_str());
    csv.exceptions(std::ios::failbit | std::ios::badbit);
    best_csv.exceptions(std::ios::failbit | std::ios::badbit);
    report.exceptions(std::ios::failbit | std::ios::badbit);
    report << "Serializer level: " << c_level << "; maximum size ratio: "
           << pair_max_ratio << " (0 means disabled; heuristic prefilter)\n";
    report << "Pair evaluation threads: " << pair_threads
           << "; accepted pairs per task: 128\n";
    report << "Input: " << bv_in_file << "\nStatus: incomplete\n"
           << "Vector limit: " << pair_limit << " (0 means all)\n"
           << "Minimum structural gain: " << pair_min_gain
           << "; minimum target gain percent: " << pair_min_pct << '\n' << std::flush;
    const char* size_columns = ",lead_bytes,target_bytes,smaller_bytes,larger_bytes,size_ratio";
    csv << "lead_id,target_id,lead_cost,target_cost,xor_cost,gain_cost,gain_pct,xor_count"
        << size_columns << '\n';
    bm::serializer<pair_bv> serializer;
    serializer.set_compression_level(c_level);
    bm::serializer<pair_bv>::buffer buffer;
    csv << std::setprecision(12);
    std::vector<std::unique_ptr<pair_item> > items;
    bm::streams_decoder decoder(input);
    bm::streams_deserializer<pair_bv> reader;
    bm::id64_t memory = 0;
    const auto start = pair_clock::now();
    auto last = start;
    while (input.tellg() != file_size && (!pair_limit || items.size() < pair_limit))
    {
        const auto before = input.tellg();
        std::unique_ptr<pair_item> item(new pair_item);
        if (!reader.deserialize(item->bv, decoder))
            throw std::runtime_error("Error reading pair input");
        const auto after = input.tellg();
        if (before < std::streampos(0) || after <= before || after > file_size)
            throw std::runtime_error("Invalid pair input boundary");
        item->bv.optimize();
        item->cost = pair_cost(bm::calc_complexity(item->bv));
        serializer.serialize(item->bv, buffer);
        item->bytes = buffer.size();
        if (!item->bytes) throw std::runtime_error("Empty serialized BLOB");
        pair_bv::statistics stat;
        item->bv.calc_stat(&stat);
        memory += stat.memory_used;
        items.push_back(std::move(item));
        const auto now = pair_clock::now();
        if (!is_silent && now - last >= std::chrono::seconds(10))
        {
            cout << "Loaded " << items.size() << " vectors; input " << after << '/'
                 << file_size << "; bvector allocation MiB=" << double(memory) / (1024.0 * 1024) << endl;
            last = now;
        }
    }
    const bool limited = input.tellg() != file_size;
    const bm::id64_t n = items.size();
    if (n > 1 && n > (std::numeric_limits<bm::id64_t>::max)() / (n - 1))
        throw std::runtime_error("Pair count overflow");
    const bm::id64_t total = n * (n - 1) / 2;
    cout << "Loaded " << n << " vectors (" << (limited ? "prefix only" : "whole collection")
         << "); bvector allocation MiB=" << double(memory) / (1024.0 * 1024)
         << "; pairs=" << total << endl;
    std::vector<pair_best> best(items.size());
    std::vector<bool> participants(items.size(), false);
    bm::id64_t considered = 0, scheduled = 0, completed = 0;
    bm::id64_t skipped = 0, qualified = 0, targets = 0;
    bm::id64_t bins[5] = {};
    bm::id64_t ratio_candidates[4] = {}, ratio_inspected[4] = {}, ratio_gains[4] = {};
    double ratio_positive_bytes[4] = {};
    auto search_start = pair_clock::now();
    last = search_start;

    typedef bm::thread_pool<bm::task_descr*,
                            bm::spin_lock<bm::pad0_struct> > pair_pool_type;
    typedef bm::task_batch<pair_bv::allocator_type> pair_task_batch_type;
    pair_pool_type pool;
    pool.start(pair_threads);

    const std::size_t pairs_per_task = 128;
    const std::size_t tasks_per_window = std::size_t(pair_threads) * 4;
    std::size_t next_i = 0, next_j = 1;

    auto next_pair = [&](pair_coord& coord) -> bool
    {
        while (next_i < items.size())
        {
            if (next_j >= items.size())
            {
                ++next_i;
                next_j = next_i + 1;
                continue;
            }
            const std::size_t i = next_i, j = next_j++;
            ++considered;
            // Keep the boundary inclusive; skip before any XOR complexity work.
            if (pair_max_ratio != 0 &&
                pair_size_ratio(items[i]->bytes, items[j]->bytes) > pair_max_ratio)
            {
                ++skipped;
                continue;
            }
            coord.first = i; coord.second = j;
            return true;
        }
        return false;
    };

    auto write_progress = [&](bool force)
    {
        const auto now = pair_clock::now();
        if (!force && now - last < std::chrono::seconds(10))
            return;
        const double seconds = std::chrono::duration<double>(now - search_start).count();
        const double rate = seconds > 0 ? double(considered) / seconds : 0;
        if (!is_silent)
            cout << "Pairs " << considered << '/' << total << " ("
                 << (total ? 100.0 * double(considered) / double(total) : 100.0)
                 << "%), pairs/s=" << rate << ", elapsed_s=" << seconds
                 << ", ETA_s=" << (rate > 0 ? double(total - considered) / rate : 0)
                 << ", scheduled=" << scheduled << ", completed=" << completed
                 << ", pending=" << scheduled - completed
                 << ", ratio_skipped=" << skipped
                 << ", qualifying=" << qualified << ", targets=" << targets << endl;
        last = now;
    };

    bool input_pairs_done = false;
    while (!input_pairs_done)
    {
        std::vector<pair_job> jobs;
        jobs.reserve(tasks_per_window);
        while (jobs.size() < tasks_per_window && !input_pairs_done)
        {
            jobs.emplace_back();
            pair_job& job = jobs.back();
            job.pairs.reserve(pairs_per_task);
            pair_coord coord;
            while (job.pairs.size() < pairs_per_task)
            {
                if (!next_pair(coord))
                {
                    input_pairs_done = true;
                    break;
                }
                job.pairs.push_back(coord);
            }
            job.considered_end = considered;
            if (job.pairs.empty())
                jobs.pop_back();
        }
        if (jobs.empty())
            break;

        pair_task_batch_type batch;
        batch.get_task_vector().reserve(jobs.size());
        for (pair_job& job : jobs)
        {
            job.matches.reserve(job.pairs.size());
            pair_job* job_ptr = &job;
            bm::task_function_t task([job_ptr, &items](void*)
            {
                for (const pair_coord& coord : job_ptr->pairs)
                {
                    const std::size_t i = coord.first, j = coord.second;
                    const std::size_t lead = items[i]->cost <= items[j]->cost ? i : j;
                    const std::size_t target = lead == i ? j : i;
                    const auto residual = bm::calc_complexity_xor(items[i]->bv, items[j]->bv);
                    const auto cost = pair_cost(residual);
                    const auto baseline = items[target]->cost;
                    if (cost < baseline)
                    {
                        const auto gain = baseline - cost;
                        const double pct = 100.0 * double(gain) / double(baseline);
                        if (gain >= pair_min_gain && pct >= pair_min_pct)
                            job_ptr->matches.push_back(
                                pair_match{i, j, cost, residual.count});
                    }
                }
                return 0;
            });
            batch.add(task, 0);
            scheduled += job.pairs.size();
        }

        auto& queue = pool.get_job_queue();
        pair_task_batch_type::size_type submitted = 0;
        try
        {
            for (; submitted < batch.size(); ++submitted)
                queue.push(batch.get_task(submitted));
        }
        catch (...)
        {
            // Submitted tasks capture jobs in this window. Keep those jobs
            // alive until all submitted workers release them.
            for (pair_task_batch_type::size_type k = 0; k < submitted; ++k)
                while (!batch.get_task(k)->done.load(std::memory_order_acquire))
                    std::this_thread::sleep_for(std::chrono::milliseconds(1));
            throw;
        }

        std::exception_ptr window_error;
        for (std::size_t k = 0; k < jobs.size(); ++k)
        {
            bm::task_descr* task =
                batch.get_task(pair_task_batch_type::size_type(k));
            while (!task->done.load(std::memory_order_acquire))
            {
                write_progress(false);
                std::this_thread::sleep_for(std::chrono::milliseconds(1));
            }
            pair_job& job = jobs[k];
            completed += job.pairs.size();
            if (task->err_code && !window_error)
                window_error = std::make_exception_ptr(
                    std::runtime_error("Pair evaluation task failed"));
            if (!window_error)
            {
                try
                {
                    for (const pair_match& match : job.matches)
                    {
                        const std::size_t i = match.first, j = match.second;
                        const std::size_t lead =
                            items[i]->cost <= items[j]->cost ? i : j;
                        const std::size_t target = lead == i ? j : i;
                        const auto baseline = items[target]->cost;
                        const auto gain = baseline - match.residual_cost;
                        const double pct =
                            100.0 * double(gain) / double(baseline);
                        ++qualified;
                        participants[i] = participants[j] = true;
                        csv << lead << ',' << target << ','
                            << items[lead]->cost << ',' << baseline << ','
                            << match.residual_cost << ',' << gain << ','
                            << pct << ',' << match.residual_count;
                        pair_write_sizes(csv, items[lead]->bytes,
                                         items[target]->bytes);
                        csv << '\n';
                        ++ratio_candidates[pair_ratio_bucket(pair_size_ratio(
                            items[lead]->bytes, items[target]->bytes))];
                        auto& b = best[target];
                        if (!b.found) ++targets;
                        if (!b.found || gain > b.gain ||
                            (gain == b.gain && lead < b.lead))
                        {
                            b.found = true; b.lead = lead;
                            b.residual = match.residual_cost;
                            b.gain = gain; b.count = match.residual_count;
                        }
                    }
                }
                catch (...)
                {
                    window_error = std::current_exception();
                }
            }
            write_progress(false);
        }
        if (window_error)
            std::rethrow_exception(window_error);
        csv.flush();
    }
    pool.set_stop_mode(pair_pool_type::stop_when_done);
    pool.join();
    write_progress(true);
    best_csv << "lead_id,target_id,lead_cost,target_cost,xor_cost,gain_cost,gain_pct,xor_count";
    best_csv << size_columns;
    if (pair_inspect)
        best_csv << ",measured_lead_id,measured_target_id,independent_bytes,xor_bytes,lead_plus_xor_bytes,gain_bytes,pair_gain_pct,target_gain_pct";
    best_csv << '\n' << std::setprecision(12);
    bm::id64_t measured_gains = 0, inspected = 0;
    std::vector<pair_measured_best> measured_best(items.size());
    last = pair_clock::now();
    for (std::size_t target = 0; target < items.size(); ++target)
    {
        const auto& b = best[target];
        if (!b.found) continue;
        const double pct = 100.0 * double(b.gain) / double(items[target]->cost);
        ++bins[pct < 25 ? 0 : pct < 50 ? 1 : pct < 75 ? 2 : pct < 100 ? 3 : 4];
        best_csv << b.lead << ',' << target << ',' << items[b.lead]->cost << ','
                 << items[target]->cost << ',' << b.residual << ',' << b.gain << ','
                 << pct << ',' << b.count;
        pair_write_sizes(best_csv, items[b.lead]->bytes, items[target]->bytes);
        if (pair_inspect)
        {
            const auto a = items[b.lead]->bytes, t = items[target]->bytes;
            const auto bucket = pair_ratio_bucket(pair_size_ratio(a, t));
            ++ratio_inspected[bucket];
            pair_bv residual(items[target]->bv);
            residual ^= items[b.lead]->bv;
            residual.optimize();
            serializer.serialize(residual, buffer);
            const bm::id64_t independent = a + t;
            const bm::id64_t referenced = (std::min)(a, t) + buffer.size();
            const auto measured_lead = a < t ? b.lead : t < a ? target : (std::min)(b.lead, target);
            const auto measured_target = measured_lead == b.lead ? target : b.lead;
            const double gain = double(independent) - double(referenced);
            if (independent > referenced)
            {
                ++measured_gains; ++ratio_gains[bucket];
                ratio_positive_bytes[bucket] += gain;
                auto& m = measured_best[measured_target];
                const auto saved = independent - referenced;
                if (!m.found || saved > m.gain || (saved == m.gain && measured_lead < m.lead))
                {
                    m.found = true; m.lead = measured_lead;
                    m.xor_bytes = buffer.size(); m.gain = saved;
                }
            }
            best_csv << ',' << measured_lead << ',' << measured_target << ',' << independent
                     << ',' << buffer.size() << ',' << referenced << ',' << gain << ','
                     << (independent ? 100.0 * gain / double(independent) : 0) << ','
                     << 100.0 * gain / double((std::max)(a, t));
        }
        best_csv << '\n';
        ++inspected;
        const auto now = pair_clock::now();
        if (now - last >= std::chrono::seconds(10))
        {
            if (!is_silent) cout << "Best candidates written/inspected: " << inspected << '/' << targets << endl;
            best_csv.flush(); last = now;
        }
    }
    csv.close(); best_csv.close();
    std::vector<std::size_t> measured_rank;
    for (std::size_t id = 0; id < measured_best.size(); ++id)
        if (measured_best[id].found) measured_rank.push_back(id);
    std::sort(measured_rank.begin(), measured_rank.end(), [&](std::size_t a, std::size_t b)
        { return measured_best[a].gain != measured_best[b].gain ?
                 measured_best[a].gain > measured_best[b].gain : a < b; });
    if (pair_inspect)
    {
        std::ofstream measured_csv(measured_name.c_str());
        measured_csv.exceptions(std::ios::failbit | std::ios::badbit);
        measured_csv << "measured_lead_id,measured_target_id,lead_bytes,target_bytes,smaller_bytes,larger_bytes,size_ratio,xor_bytes,independent_bytes,lead_plus_xor_bytes,gain_bytes,pair_gain_pct,target_gain_pct\n"
                     << std::setprecision(12);
        for (auto id : measured_rank)
        {
            const auto& m = measured_best[id];
            const auto a = items[m.lead]->bytes, t = items[id]->bytes;
            measured_csv << m.lead << ',' << id;
            pair_write_sizes(measured_csv, a, t);
            measured_csv << ',' << m.xor_bytes << ',' << a + t << ',' << a + m.xor_bytes
                         << ',' << m.gain << ',' << 100.0 * double(m.gain) / double(a + t)
                         << ',' << 100.0 * double(m.gain) / double(t) << '\n';
        }
        measured_csv.close();
    }
    std::vector<std::size_t> ranked;
    for (std::size_t id = 0; id < best.size(); ++id)
        if (best[id].found) ranked.push_back(id);
    const auto top_count = (std::min)(std::size_t(10), ranked.size());
    std::partial_sort(ranked.begin(), ranked.begin() + std::vector<std::size_t>::difference_type(top_count), ranked.end(),
        [&](std::size_t a, std::size_t b)
        { return best[a].gain != best[b].gain ? best[a].gain > best[b].gain : a < b; });
    auto summary = [&](std::ostream& out)
    {
        out << "Search complete for " << (limited ? "loaded prefix" : "whole collection")
            << ". Vectors=" << n << "; pairs considered=" << considered
            << "; pairs scheduled=" << scheduled
            << "; pairs completed=" << completed << "; ratio skipped=" << skipped
            << "\nMaximum size ratio=" << pair_max_ratio << " (0 means disabled); skipped outcomes unknown"
            << "\nPair evaluation threads=" << pair_threads
            << "\nQualifying pairs=" << qualified << "; candidate targets=" << targets
            << " (" << (n ? 100.0 * double(targets) / double(n) : 0) << "%)"
            << "; participating vectors=" << std::count(participants.begin(), participants.end(), true)
            << "\nBest-target gain distribution [0,25), [25,50), [50,75), [75,100), 100: "
            << bins[0] << ',' << bins[1] << ',' << bins[2] << ',' << bins[3] << ',' << bins[4]
            << "\nStructural gains are heuristic units, not serialized bytes. Overlapping gains are not additive.\n";
        if (pair_inspect)
            out << "Inspected best pairs=" << inspected << "; measured improving pairs=" << measured_gains
                << "; distinct measured targets=" << measured_rank.size()
                << "; serializer level=" << c_level
                << "\nMeasured sizes include BLOB headers, exclude reference metadata; no reference graph constructed.\n";
        out << "Size ratio buckets: qualifying_pairs,inspected_best_pairs,improving_pairs,positive_gain_bytes_sum\n";
        const char* labels[] = {"[1,2)", "[2,4)", "[4,8)", "[8,infinity)"};
        for (unsigned k = 0; k < 4; ++k)
        {
            out << labels[k] << ':' << ratio_candidates[k];
            if (pair_inspect)
                out << ',' << ratio_inspected[k] << ',' << ratio_gains[k] << ',' << ratio_positive_bytes[k];
            else out << ",not inspected";
            out << '\n';
        }
        out << "Buckets cover heuristic candidates only; positive byte sums overlap and are not collection savings.\n";
        if (pair_inspect)
        {
            out << "Top measured candidates, deduplicated by target (lead,target,gain_bytes):\n";
            for (std::size_t k = 0; k < (std::min)(std::size_t(10), measured_rank.size()); ++k)
            {
                const auto id = measured_rank[k];
                out << measured_best[id].lead << ',' << id << ',' << measured_best[id].gain << '\n';
            }
        }
        out << "Bvector allocation MiB=" << double(memory) / (1024.0 * 1024) << '\n';
        out << "Top candidates by structural gain (lead,target,gain_cost):\n";
        for (std::size_t k = 0; k < top_count; ++k)
        {
            const auto id = ranked[k];
            out << best[id].lead << ',' << id << ',' << best[id].gain << '\n';
        }
        out << "Elapsed seconds=" << std::chrono::duration<double>(pair_clock::now() - start).count() << '\n';
    };
    summary(cout); summary(report);
    report << "Status: complete\n";
    report.close();
}
} // namespace

#endif
