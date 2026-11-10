/** @file inv_list_xbv_ref_eval.h
    @brief Cross-bvector (xbv) research support for inv_list.
    @note Development-only; included by inv_list.cpp after shared options and I/O.
*/
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

// Offline reference-serializer experiment for inv_list.cpp.
// Included after inv_list_xbv_pair_search.h; not a public library header.
#ifndef INV_LIST_XBV_REF_EVAL_H
#define INV_LIST_XBV_REF_EVAL_H

namespace
{

struct ref_candidate
{
    std::size_t id = 0;
    long long structural_gain = 0;
};

struct ref_measurement
{
    bm::id64_t bytes = 0;
    bm::id64_t equal_blocks = 0, xor_blocks = 0, chain_blocks = 0;
    double model_us = 0, serialize_us = 0, deserialize_us = 0;
    bool model_ok = false, roundtrip_ok = false;
};

struct ref_eval_result
{
    std::size_t target = 0, plain_lead = 0;
    std::vector<std::size_t> candidates;
    bm::id64_t standalone_bytes = 0, plain_xor_bytes = 0;
    ref_measurement ref1, refk;
    std::string error;
};

struct ref_eval_job
{
    std::size_t target = 0;
    ref_eval_result result;
};

long long ref_signed_delta(bm::id64_t smaller, bm::id64_t larger)
{
    const bm::id64_t limit = bm::id64_t((std::numeric_limits<long long>::max)());
    if (smaller > limit || larger > limit)
        throw std::runtime_error("Serialized size exceeds signed report range");
    return static_cast<long long>(smaller) - static_cast<long long>(larger);
}

double ref_gain_pct(bm::id64_t baseline, bm::id64_t result)
{
    return baseline ? 100.0 * double(ref_signed_delta(baseline, result)) /
                          double(baseline) : 0.0;
}

ref_measurement ref_serialize_target(
    const std::vector<std::unique_ptr<pair_item> >& items,
    std::size_t target, const std::vector<std::size_t>& candidates,
    std::size_t candidate_count)
{
    typedef bm::serializer<pair_bv> serializer_type;
    ref_measurement result;
    serializer_type::bv_ref_vector_type refs;
    refs.add(&items[target]->bv, pair_bv::size_type(target));
    for (std::size_t i = 0; i < candidate_count; ++i)
        refs.add(&items[candidates[i]]->bv,
                 pair_bv::size_type(candidates[i]));

    serializer_type serializer;
    serializer.set_compression_level(c_level);
    serializer.set_ref_vectors(&refs);
    serializer_type::xor_sim_model_type model;
    const auto model_start = pair_clock::now();
    result.model_ok = serializer.compute_sim_model(model, refs,
                                                    bm::xor_sim_params());
    const auto model_end = pair_clock::now();
    result.model_us = std::chrono::duration<double, std::micro>(
                          model_end - model_start).count();

    serializer_type::buffer buffer;
    const auto serialize_start = pair_clock::now();
    if (result.model_ok)
    {
        serializer.set_sim_model(&model);
        serializer.set_curr_ref_idx(0);
        serializer.serialize(items[target]->bv, buffer);
    }
    else
    {
        serializer.set_ref_vectors(0);
        serializer.serialize(items[target]->bv, buffer);
    }
    const auto serialize_end = pair_clock::now();
    result.serialize_us = std::chrono::duration<double, std::micro>(
                              serialize_end - serialize_start).count();
    result.bytes = buffer.size();
    const auto* stats = serializer.get_compression_stat();
    result.equal_blocks = stats[bm::set_block_ref_eq];
    result.xor_blocks = stats[bm::set_block_xor_ref32];
    result.chain_blocks = stats[bm::set_block_xor_chain];

    // The decoder needs the external row IDs actually encoded in the BLOB.
    // Supplying all candidates is a conservative dependency superset.
    serializer_type::bv_ref_vector_type decode_refs;
    for (std::size_t i = 0; i < candidate_count; ++i)
        decode_refs.add(&items[candidates[i]]->bv,
                        pair_bv::size_type(candidates[i]));
    pair_bv restored;
    const auto deserialize_start = pair_clock::now();
    const std::size_t consumed = bm::deserialize(
        restored, buffer.data(), 0, result.model_ok ? &decode_refs : 0);
    const auto deserialize_end = pair_clock::now();
    result.deserialize_us = std::chrono::duration<double, std::micro>(
                                deserialize_end - deserialize_start).count();
    result.roundtrip_ok = consumed == buffer.size() &&
                          restored == items[target]->bv;
    return result;
}

void ref_evaluate_one(ref_eval_job& job,
                      const std::vector<std::unique_ptr<pair_item> >& items,
                      const std::vector<std::vector<ref_candidate> >& candidate_sets)
{
    try
    {
        ref_eval_result& result = job.result;
        result.target = job.target;
        result.standalone_bytes = items[job.target]->bytes;
        const auto& selected = candidate_sets[job.target];
        bm::serializer<pair_bv> serializer;
        serializer.set_compression_level(c_level);
        bm::serializer<pair_bv>::buffer buffer;
        std::vector<std::pair<bm::id64_t, std::size_t> > measured;
        measured.reserve(selected.size());
        for (const auto& candidate : selected)
        {
            pair_bv residual(items[job.target]->bv);
            residual ^= items[candidate.id]->bv;
            residual.optimize();
            serializer.serialize(residual, buffer);
            measured.push_back(std::make_pair(bm::id64_t(buffer.size()),
                                               candidate.id));
        }
        std::sort(measured.begin(), measured.end(),
                  [](const std::pair<bm::id64_t, std::size_t>& a,
                     const std::pair<bm::id64_t, std::size_t>& b)
                  { return a.first != b.first ? a.first < b.first :
                                                a.second < b.second; });
        result.plain_xor_bytes = measured.front().first;
        result.plain_lead = measured.front().second;
        for (const auto& value : measured)
            result.candidates.push_back(value.second);
        result.ref1 = ref_serialize_target(items, job.target,
                                           result.candidates, 1);
        result.refk = ref_serialize_target(items, job.target,
                                           result.candidates,
                                           result.candidates.size());
        if (!result.ref1.roundtrip_ok || !result.refk.roundtrip_ok)
            result.error = "Reference serialization round-trip failed";
    }
    catch (const std::exception& ex)
    {
        job.result.error = ex.what();
    }
    catch (...)
    {
        job.result.error = "Unknown reference evaluation failure";
    }
}

void ref_write_candidate_ids(std::ostream& out,
                             const std::vector<std::size_t>& ids)
{
    for (std::size_t i = 0; i < ids.size(); ++i)
    {
        if (i) out << ';';
        out << ids[i];
    }
}

void ref_write_result(std::ostream& out, const ref_eval_result& r)
{
    const auto plain_gain = ref_signed_delta(r.standalone_bytes,
                                             r.plain_xor_bytes);
    const auto ref1_gain = ref_signed_delta(r.standalone_bytes, r.ref1.bytes);
    const auto refk_gain = ref_signed_delta(r.standalone_bytes, r.refk.bytes);
    out << r.target << ',';
    ref_write_candidate_ids(out, r.candidates);
    out << ',' << r.candidates.size() << ',' << r.standalone_bytes
        << ',' << r.plain_lead << ',' << r.plain_xor_bytes
        << ',' << plain_gain << ','
        << ref_gain_pct(r.standalone_bytes, r.plain_xor_bytes)
        << ',' << r.ref1.bytes << ',' << ref1_gain << ','
        << ref_gain_pct(r.standalone_bytes, r.ref1.bytes)
        << ',' << r.refk.bytes << ',' << refk_gain << ','
        << ref_gain_pct(r.standalone_bytes, r.refk.bytes)
        << ',' << ref_signed_delta(r.plain_xor_bytes, r.refk.bytes)
        << ',' << r.ref1.model_us << ',' << r.ref1.serialize_us
        << ',' << r.ref1.deserialize_us << ',' << unsigned(r.ref1.model_ok)
        << ',' << unsigned(r.ref1.roundtrip_ok)
        << ',' << r.refk.model_us << ',' << r.refk.serialize_us
        << ',' << r.refk.deserialize_us << ',' << unsigned(r.refk.model_ok)
        << ',' << r.refk.equal_blocks << ',' << r.refk.xor_blocks
        << ',' << r.refk.chain_blocks << ',' << unsigned(r.refk.roundtrip_ok)
        << '\n';
}

void analyze_references()
{
    if (bv_in_file.empty() || ref_pairs_in.empty() || ref_csv.empty())
        throw std::runtime_error(
            "Reference evaluation requires -bvin FILE, -pairs-in FILE and -ref-csv FILE");
    const std::string report_name = ref_csv + ".report.txt";
    for (const auto& name : std::vector<std::string>{ref_csv, report_name})
    {
        std::ifstream exists(name.c_str(), std::ios::binary);
        if (exists.good())
            throw std::runtime_error("Refusing to overwrite output file: " + name);
    }
    std::ifstream input(bv_in_file.c_str(), std::ios::binary);
    if (!input)
        throw std::runtime_error("Cannot open bvector collection");
    input.seekg(0, std::ios::end);
    const auto file_size = input.tellg();
    input.seekg(0);
    if (!input || file_size <= std::streampos(0))
        throw std::runtime_error("Empty or inaccessible bvector collection");

    std::ofstream csv(ref_csv.c_str()), report(report_name.c_str());
    csv.exceptions(std::ios::failbit | std::ios::badbit);
    report.exceptions(std::ios::failbit | std::ios::badbit);
    report << "Input: " << bv_in_file << "\nPair CSV: " << ref_pairs_in
           << "\nStatus: incomplete\nVector limit: " << pair_limit
           << " (0 means all)\nReference candidates: " << ref_top
           << "\nSerializer level: " << c_level
           << "\nEvaluation threads: " << pair_threads << '\n' << std::flush;

    std::vector<std::unique_ptr<pair_item> > items;
    bm::streams_decoder decoder(input);
    bm::streams_deserializer<pair_bv> reader;
    bm::serializer<pair_bv> serializer;
    serializer.set_compression_level(c_level);
    bm::serializer<pair_bv>::buffer buffer;
    bm::id64_t memory = 0;
    const auto start = pair_clock::now();
    auto last = start;
    while (input.tellg() != file_size &&
           (!pair_limit || items.size() < pair_limit))
    {
        const auto before = input.tellg();
        std::unique_ptr<pair_item> item(new pair_item);
        if (!reader.deserialize(item->bv, decoder))
            throw std::runtime_error("Error reading bvector collection");
        const auto after = input.tellg();
        if (before < std::streampos(0) || after <= before || after > file_size)
            throw std::runtime_error("Invalid bvector collection boundary");
        item->bv.optimize();
        item->cost = pair_cost(bm::calc_complexity(item->bv));
        serializer.serialize(item->bv, buffer);
        item->bytes = buffer.size();
        pair_bv::statistics stat;
        item->bv.calc_stat(&stat);
        memory += stat.memory_used;
        items.push_back(std::move(item));
        const auto now = pair_clock::now();
        if (!is_silent && now - last >= std::chrono::seconds(10))
        {
            cout << "Loaded " << items.size() << " vectors; input " << after
                 << '/' << file_size << "; bvector allocation MiB="
                 << double(memory) / (1024.0 * 1024) << endl;
            last = now;
        }
    }
    const bool limited = input.tellg() != file_size;
    cout << "Loaded " << items.size() << " vectors ("
         << (limited ? "prefix only" : "whole collection")
         << "); bvector allocation MiB="
         << double(memory) / (1024.0 * 1024) << endl;

    std::ifstream pairs(ref_pairs_in.c_str());
    if (!pairs)
        throw std::runtime_error("Cannot open pair CSV");
    std::string line;
    if (!std::getline(pairs, line))
        throw std::runtime_error("Empty pair CSV");
    if (!line.empty() && line.back() == '\r') line.pop_back();
    const auto header = inv_list_csv::split(line);
    const std::size_t lead_id_col = inv_list_csv::column(header, "lead_id");
    const std::size_t target_id_col = inv_list_csv::column(header, "target_id");
    const std::size_t lead_cost_col = inv_list_csv::column(header, "lead_cost");
    const std::size_t target_cost_col = inv_list_csv::column(header, "target_cost");
    const std::size_t xor_cost_col = inv_list_csv::column(header, "xor_cost");
    const std::size_t lead_bytes_col = inv_list_csv::column(header, "lead_bytes");
    const std::size_t target_bytes_col = inv_list_csv::column(header, "target_bytes");
    const std::size_t required_col = (std::max)(
        (std::max)((std::max)(lead_id_col, target_id_col),
                   (std::max)(lead_cost_col, target_cost_col)),
        (std::max)(xor_cost_col, (std::max)(lead_bytes_col,
                                            target_bytes_col)));
    std::vector<std::vector<ref_candidate> > candidate_sets(items.size());
    bm::id64_t rows = 0, eligible_rows = 0;
    for (std::size_t line_no = 2; std::getline(pairs, line); ++line_no)
    {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty()) continue;
        ++rows;
        const auto fields = inv_list_csv::split(line);
        if (fields.size() <= required_col)
            throw std::runtime_error("Short row in pair CSV at line " +
                                     std::to_string(line_no));
        const auto csv_lead = std::size_t(inv_list_csv::parse_u64(
            fields[lead_id_col], line_no, "lead_id"));
        const auto csv_target = std::size_t(inv_list_csv::parse_u64(
            fields[target_id_col], line_no, "target_id"));
        if (csv_lead >= items.size() || csv_target >= items.size())
            continue;
        const auto lead_cost = inv_list_csv::parse_u64(
            fields[lead_cost_col], line_no, "lead_cost");
        const auto target_cost = inv_list_csv::parse_u64(
            fields[target_cost_col], line_no, "target_cost");
        const auto xor_cost = inv_list_csv::parse_u64(
            fields[xor_cost_col], line_no, "xor_cost");
        const auto lead_bytes = inv_list_csv::parse_u64(
            fields[lead_bytes_col], line_no, "lead_bytes");
        const auto target_bytes = inv_list_csv::parse_u64(
            fields[target_bytes_col], line_no, "target_bytes");
        if (lead_cost != items[csv_lead]->cost ||
            target_cost != items[csv_target]->cost ||
            lead_bytes != items[csv_lead]->bytes ||
            target_bytes != items[csv_target]->bytes)
        {
            std::ostringstream msg;
            msg << "Pair CSV does not match loaded vectors/serializer level at line "
                << line_no;
            throw std::runtime_error(msg.str());
        }
        const auto lead_key = std::make_pair(lead_bytes, csv_lead);
        const auto target_key = std::make_pair(target_bytes, csv_target);
        const std::size_t candidate = lead_key < target_key ? csv_lead : csv_target;
        const std::size_t target = candidate == csv_lead ? csv_target : csv_lead;
        const bm::id64_t baseline = items[target]->cost;
        long long gain;
        if (baseline >= xor_cost)
            gain = static_cast<long long>(baseline - xor_cost);
        else
            gain = -static_cast<long long>(xor_cost - baseline);
        candidate_sets[target].push_back(ref_candidate{candidate, gain});
        ++eligible_rows;
    }

    std::vector<std::size_t> targets;
    for (std::size_t target = 0; target < candidate_sets.size(); ++target)
    {
        auto& candidates = candidate_sets[target];
        if (candidates.empty()) continue;
        std::sort(candidates.begin(), candidates.end(),
                  [](const ref_candidate& a, const ref_candidate& b)
                  { return a.id != b.id ? a.id < b.id :
                         a.structural_gain > b.structural_gain; });
        candidates.erase(std::unique(candidates.begin(), candidates.end(),
                         [](const ref_candidate& a, const ref_candidate& b)
                         { return a.id == b.id; }), candidates.end());
        std::sort(candidates.begin(), candidates.end(),
                  [](const ref_candidate& a, const ref_candidate& b)
                  { return a.structural_gain != b.structural_gain ?
                           a.structural_gain > b.structural_gain : a.id < b.id; });
        if (candidates.size() > ref_top)
            candidates.resize(ref_top);
        targets.push_back(target);
    }
    cout << "Read " << rows << " qualifying pair rows; eligible rows="
         << eligible_rows << "; reference targets=" << targets.size() << endl;

    csv << "target_id,candidate_ids,candidate_count,standalone_bytes,"
           "plain_xor_lead,plain_xor_bytes,plain_xor_gain_bytes,plain_xor_gain_pct,"
           "ref1_bytes,ref1_gain_bytes,ref1_gain_pct,"
           "refk_bytes,refk_gain_bytes,refk_gain_pct,refk_vs_plain_gain_bytes,"
           "ref1_model_us,ref1_serialize_us,ref1_deserialize_us,ref1_model_ok,ref1_roundtrip_ok,"
           "refk_model_us,refk_serialize_us,refk_deserialize_us,refk_model_ok,"
           "refk_equal_blocks,refk_xor_blocks,refk_chain_blocks,refk_roundtrip_ok\n"
        << std::setprecision(12);

    typedef bm::thread_pool<bm::task_descr*,
                            bm::spin_lock<bm::pad0_struct> > ref_pool_type;
    typedef bm::task_batch<pair_bv::allocator_type> ref_task_batch_type;
    ref_pool_type pool;
    pool.start(pair_threads);
    const std::size_t window_size = std::size_t(pair_threads) * 4;
    std::size_t completed = 0;
    bm::id64_t plain_improving = 0, ref1_improving = 0;
    bm::id64_t refk_improving = 0, refk_beats_plain = 0;
    long long plain_gain_sum = 0, ref1_gain_sum = 0, refk_gain_sum = 0;
    last = pair_clock::now();
    std::string worker_error;
    for (std::size_t begin = 0; begin < targets.size() && worker_error.empty();
         begin += window_size)
    {
        const std::size_t count = (std::min)(window_size,
                                              targets.size() - begin);
        std::vector<ref_eval_job> jobs;
        jobs.reserve(count);
        ref_task_batch_type batch;
        batch.get_task_vector().reserve(count);
        for (std::size_t i = 0; i < count; ++i)
        {
            jobs.emplace_back();
            jobs.back().target = targets[begin + i];
        }
        for (auto& job : jobs)
        {
            ref_eval_job* job_ptr = &job;
            bm::task_function_t task([job_ptr, &items, &candidate_sets](void*)
            {
                ref_evaluate_one(*job_ptr, items, candidate_sets);
                return job_ptr->result.error.empty() ? 0 : 1;
            });
            batch.add(task, 0);
        }
        auto& queue = pool.get_job_queue();
        ref_task_batch_type::size_type submitted = 0;
        try
        {
            for (; submitted < batch.size(); ++submitted)
                queue.push(batch.get_task(submitted));
        }
        catch (...)
        {
            for (ref_task_batch_type::size_type i = 0; i < submitted; ++i)
                while (!batch.get_task(i)->done.load(std::memory_order_acquire))
                    std::this_thread::sleep_for(std::chrono::milliseconds(1));
            throw;
        }
        for (std::size_t i = 0; i < jobs.size(); ++i)
        {
            bm::task_descr* task = batch.get_task(
                ref_task_batch_type::size_type(i));
            while (!task->done.load(std::memory_order_acquire))
            {
                const auto now = pair_clock::now();
                if (!is_silent && now - last >= std::chrono::seconds(10))
                {
                    const double seconds = std::chrono::duration<double>(
                                               now - start).count();
                    cout << "Reference targets " << completed << '/'
                         << targets.size() << "; elapsed_s=" << seconds << endl;
                    last = now;
                }
                std::this_thread::sleep_for(std::chrono::milliseconds(1));
            }
            const auto& r = jobs[i].result;
            if (task->err_code || !r.error.empty())
            {
                if (worker_error.empty()) worker_error = r.error.empty() ?
                    "Reference evaluation task failed" : r.error;
                continue;
            }
            ref_write_result(csv, r);
            ++completed;
            const long long pg = ref_signed_delta(r.standalone_bytes,
                                                   r.plain_xor_bytes);
            const long long r1g = ref_signed_delta(r.standalone_bytes,
                                                   r.ref1.bytes);
            const long long rkg = ref_signed_delta(r.standalone_bytes,
                                                   r.refk.bytes);
            if (pg > 0) { ++plain_improving; plain_gain_sum += pg; }
            if (r1g > 0) { ++ref1_improving; ref1_gain_sum += r1g; }
            if (rkg > 0) { ++refk_improving; refk_gain_sum += rkg; }
            if (r.refk.bytes < r.plain_xor_bytes) ++refk_beats_plain;
        }
        csv.flush();
    }
    pool.set_stop_mode(ref_pool_type::stop_when_done);
    pool.join();
    if (!worker_error.empty())
        throw std::runtime_error(worker_error);
    csv.close();

    auto summary = [&](std::ostream& out)
    {
        out << "Reference evaluation complete for "
            << (limited ? "loaded prefix" : "whole collection")
            << ". Vectors=" << items.size() << "; qualifying pair rows="
            << rows << "; eligible rows=" << eligible_rows
            << "; evaluated targets=" << completed
            << "\nCandidate limit=" << ref_top
            << "; evaluation threads=" << pair_threads
            << "; serializer level=" << c_level
            << "\nImproving targets: plain XOR=" << plain_improving
            << "; one-reference=" << ref1_improving
            << "; reference-table=" << refk_improving
            << "; reference-table beats plain XOR=" << refk_beats_plain
            << "\nPositive target-byte gain sums: plain XOR=" << plain_gain_sum
            << "; one-reference=" << ref1_gain_sum
            << "; reference-table=" << refk_gain_sum
            << "\nPlain XOR and reference BLOB sizes exclude dependency storage. "
               "Reference BLOBs include encoded row IDs and masks; candidate_ids "
               "is the supplied dependency superset.\n"
            << "Candidate edges descend by (standalone_bytes,id), so the supplied "
               "dependency graph is acyclic. Gain sums exclude collection framing.\n"
            << "Bvector allocation MiB=" << double(memory) / (1024.0 * 1024)
            << "\nElapsed seconds="
            << std::chrono::duration<double>(pair_clock::now() - start).count()
            << '\n';
    };
    summary(cout); summary(report);
    report << "Status: complete\n";
    report.close();
}

} // namespace

#endif
