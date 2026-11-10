/** @file inv_list_xbv_pair_shard.h
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

// Size-sorted, greedy sharded pair search for inv_list.cpp.
// Included after inv_list_xbv_pair_search.h and inv_list_xbv_csv.h.
#ifndef INV_LIST_XBV_PAIR_SHARD_H
#define INV_LIST_XBV_PAIR_SHARD_H

namespace
{

struct pair_index_entry
{
    std::size_t id = 0, shard = 0;
    bm::id64_t offset = 0, source_bytes = 0, bytes = 0;
    bm::id64_t cost = 0, count = 0, memory = 0;
    bool super_sparse = false;
};

struct pair_shard_range
{
    std::size_t begin = 0, end = 0, band = 0;
    bm::id64_t bytes = 0, memory = 0;
};

bm::id64_t pair_choose2(bm::id64_t n)
{
    if (n > 1 && n > (std::numeric_limits<bm::id64_t>::max)() / (n - 1))
        throw std::runtime_error("Pair count overflow");
    return n * (n - 1) / 2;
}

std::vector<pair_shard_range> pair_make_shards(
    const std::vector<pair_index_entry>& entries)
{
    std::vector<pair_shard_range> shards;
    std::size_t band_id = 0, band_begin = 0;
    while (band_begin < entries.size())
    {
        std::size_t band_end = band_begin + 1;
        if (pair_max_ratio == 0)
            band_end = entries.size();
        else
            while (band_end < entries.size() &&
                   pair_size_ratio(entries[band_begin].bytes,
                                   entries[band_end].bytes) <= pair_max_ratio)
                ++band_end;

        const std::size_t first_shard = shards.size();
        std::size_t begin = band_begin;
        while (begin < band_end)
        {
            pair_shard_range shard;
            shard.begin = begin; shard.band = band_id;
            while (begin < band_end)
            {
                const auto& e = entries[begin];
                const bm::id64_t count = begin - shard.begin;
                const bool count_full = count && count >= pair_shard_vectors;
                const bool bytes_full = count &&
                    (shard.bytes > pair_shard_bytes ||
                     e.bytes > pair_shard_bytes - shard.bytes);
                if (count_full || bytes_full) break;
                shard.bytes += e.bytes; shard.memory += e.memory; ++begin;
            }
            shard.end = begin;
            shards.push_back(shard);
        }
        if (shards.size() >= first_shard + 2 && pair_shard_tail_pct > 0)
        {
            pair_shard_range& tail = shards.back();
            pair_shard_range& prev = shards[shards.size() - 2];
            const std::size_t tail_count = tail.end - tail.begin;
            const std::size_t prev_count = prev.end - prev.begin;
            if (prev_count && 100.0 * double(tail_count) <
                              pair_shard_tail_pct * double(prev_count))
            {
                prev.end = tail.end;
                prev.bytes += tail.bytes; prev.memory += tail.memory;
                shards.pop_back();
            }
        }
        band_begin = band_end; ++band_id;
    }
    return shards;
}

void build_pair_shard_index()
{
    if (bv_in_file.empty() || pair_index_out.empty())
        throw std::runtime_error("Pair indexing requires -bvin FILE and -pair-index-out FILE");
    const std::string report_name = pair_index_out + ".report.txt";
    for (const auto& name : std::vector<std::string>{pair_index_out, report_name})
    { std::ifstream f(name.c_str(), std::ios::binary); if (f) throw std::runtime_error("Output already exists: " + name); }
    std::ifstream input(bv_in_file.c_str(), std::ios::binary);
    if (!input) throw std::runtime_error("Cannot open pair input");
    input.seekg(0, std::ios::end); const std::streampos file_size = input.tellg(); input.seekg(0);
    if (!input || file_size <= std::streampos(0)) throw std::runtime_error("Empty or inaccessible pair input");
    std::ofstream index(pair_index_out.c_str()), report(report_name.c_str());
    index.exceptions(std::ios::failbit | std::ios::badbit);
    report.exceptions(std::ios::failbit | std::ios::badbit);
    report << "Input: " << bv_in_file << "\nStatus: incomplete\nSerializer level: " << c_level
           << "\nMaximum size ratio: " << pair_max_ratio
           << "\nShard vector limit: " << pair_shard_vectors
           << "\nShard serialized-byte limit: " << pair_shard_bytes
           << "\nSmall-tail merge percent: " << pair_shard_tail_pct << '\n' << std::flush;

    std::vector<pair_index_entry> entries;
    bm::streams_decoder decoder(input); bm::streams_deserializer<pair_bv> reader;
    bm::serializer<pair_bv> serializer; serializer.set_compression_level(c_level);
    bm::serializer<pair_bv>::buffer buffer;
    const auto start = pair_clock::now(); auto last = start;
    bm::id64_t indexed_bytes = 0, total_serialized = 0, total_memory = 0, sparse_count = 0;
    while (input.tellg() != file_size && (!pair_limit || entries.size() < pair_limit))
    {
        const std::streampos before = input.tellg(); pair_bv bv;
        if (!reader.deserialize(bv, decoder)) throw std::runtime_error("Error reading pair input");
        const std::streampos after = input.tellg();
        if (before < std::streampos(0) || after <= before || after > file_size)
            throw std::runtime_error("Invalid pair input boundary");
        bv.optimize(); serializer.serialize(bv, buffer);
        pair_bv::statistics stat; bv.calc_stat(&stat);
        pair_index_entry e;
        e.id = entries.size(); e.offset = bm::id64_t(before);
        e.source_bytes = bm::id64_t(after - before); e.bytes = buffer.size();
        e.cost = pair_cost(bm::calc_complexity(bv)); e.count = bv.count();
        e.memory = stat.memory_used; e.super_sparse = is_super_sparse(bv);
        entries.push_back(e); indexed_bytes = bm::id64_t(after);
        total_serialized += e.bytes; total_memory += e.memory; sparse_count += e.super_sparse;
        const auto now = pair_clock::now();
        if (!is_silent && now - last >= std::chrono::seconds(10))
        { cout << "Indexed " << entries.size() << " vectors; input " << after << '/' << file_size << endl; last = now; }
    }
    const bool complete = indexed_bytes == bm::id64_t(file_size);
    std::sort(entries.begin(), entries.end(), [](const pair_index_entry& a, const pair_index_entry& b)
    { return a.bytes != b.bytes ? a.bytes < b.bytes : a.id < b.id; });
    std::vector<pair_shard_range> shards = pair_make_shards(entries);
    bm::id64_t planned_pairs = 0, max_shard_count = 0, max_shard_bytes = 0, max_shard_memory = 0;
    std::size_t bands = 0;
    for (std::size_t s = 0; s < shards.size(); ++s)
    {
        auto& shard = shards[s]; bands = (std::max)(bands, shard.band + 1);
        const bm::id64_t count = shard.end - shard.begin;
        planned_pairs += pair_choose2(count);
        max_shard_count = (std::max)(max_shard_count, count);
        max_shard_bytes = (std::max)(max_shard_bytes, shard.bytes);
        max_shard_memory = (std::max)(max_shard_memory, shard.memory);
        for (std::size_t i = shard.begin; i < shard.end; ++i) entries[i].shard = s;
    }
    index << "shard_id,vector_id,source_offset,source_bytes,serialized_bytes,cost,count,memory_bytes,super_sparse,serializer_level,indexed_input_bytes,source_complete\n";
    for (const auto& e : entries)
        index << e.shard << ',' << e.id << ',' << e.offset << ',' << e.source_bytes
              << ',' << e.bytes << ',' << e.cost << ',' << e.count << ',' << e.memory
              << ',' << unsigned(e.super_sparse) << ',' << c_level << ',' << indexed_bytes
              << ',' << unsigned(complete) << '\n';
    index.close();
    const bm::id64_t all_pairs = pair_choose2(entries.size());
    auto summary = [&](std::ostream& out)
    {
        out << "Pair shard index complete for " << (complete ? "whole collection" : "loaded prefix")
            << ". Vectors=" << entries.size() << "; bands=" << bands << "; shards=" << shards.size()
            << "\nIndexed input bytes=" << indexed_bytes << "; summed level-" << c_level
            << " bytes=" << total_serialized << "; estimated decoded allocation MiB="
            << double(total_memory) / (1024.0 * 1024)
            << "\nSuper-sparse vectors=" << sparse_count
            << "\nAll collection pairs=" << all_pairs << "; planned internal-shard pairs=" << planned_pairs
            << "; deliberately omitted cross-shard pairs=" << (all_pairs - planned_pairs)
            << "; search coverage pct=" << (all_pairs ? 100.0 * double(planned_pairs) / double(all_pairs) : 100.0)
            << "\nLargest shard: vectors=" << max_shard_count << "; serialized bytes=" << max_shard_bytes
            << "; estimated decoded allocation MiB=" << double(max_shard_memory) / (1024.0 * 1024)
            << "\nElapsed seconds=" << std::chrono::duration<double>(pair_clock::now() - start).count() << '\n';
        for (std::size_t s = 0; s < shards.size(); ++s)
        {
            const auto& sh = shards[s];
            out << "Shard " << s << ": band=" << sh.band << "; vectors=" << sh.end - sh.begin
                << "; size range=" << entries[sh.begin].bytes << '-' << entries[sh.end - 1].bytes
                << "; serialized bytes=" << sh.bytes << "; estimated memory MiB="
                << double(sh.memory) / (1024.0 * 1024) << "; pairs=" << pair_choose2(sh.end - sh.begin) << '\n';
        }
    };
    summary(cout); summary(report); report << "Status: complete\n";
}

std::vector<pair_index_entry> pair_read_shard_index(bool& source_complete,
                                                     bm::id64_t& indexed_bytes)
{
    std::ifstream in(pair_index_in.c_str());
    if (!in) throw std::runtime_error("Cannot open pair shard index");
    std::string line;
    if (!std::getline(in, line)) throw std::runtime_error("Empty pair shard index");
    if (!line.empty() && line.back() == '\r') line.pop_back();
    const auto h = inv_list_csv::split(line);
    const std::size_t c_shard = inv_list_csv::column(h, "shard_id");
    const std::size_t c_id = inv_list_csv::column(h, "vector_id");
    const std::size_t c_offset = inv_list_csv::column(h, "source_offset");
    const std::size_t c_source = inv_list_csv::column(h, "source_bytes");
    const std::size_t c_bytes = inv_list_csv::column(h, "serialized_bytes");
    const std::size_t c_cost = inv_list_csv::column(h, "cost");
    const std::size_t c_count = inv_list_csv::column(h, "count");
    const std::size_t c_memory = inv_list_csv::column(h, "memory_bytes");
    const std::size_t c_sparse = inv_list_csv::column(h, "super_sparse");
    const std::size_t c_serializer_level = inv_list_csv::column(h, "serializer_level");
    const std::size_t c_indexed = inv_list_csv::column(h, "indexed_input_bytes");
    const std::size_t c_complete = inv_list_csv::column(h, "source_complete");
    const std::size_t required = (std::max)({c_shard,c_id,c_offset,c_source,c_bytes,c_cost,c_count,c_memory,c_sparse,c_serializer_level,c_indexed,c_complete});
    std::vector<pair_index_entry> entries; std::vector<bool> ids;
    std::size_t previous_shard = 0; bool first = true;
    for (std::size_t line_no = 2; std::getline(in, line); ++line_no)
    {
        if (!line.empty() && line.back() == '\r')
        {
            line.pop_back();
        }
        if (line.empty())
        {
            continue;
        }
        const auto f = inv_list_csv::split(line);
        if (f.size() <= required) throw std::runtime_error("Short pair index row at line " + std::to_string(line_no));
        pair_index_entry e;
        e.shard = std::size_t(inv_list_csv::parse_u64(f[c_shard], line_no, "shard_id"));
        e.id = std::size_t(inv_list_csv::parse_u64(f[c_id], line_no, "vector_id"));
        e.offset = inv_list_csv::parse_u64(f[c_offset], line_no, "source_offset");
        e.source_bytes = inv_list_csv::parse_u64(f[c_source], line_no, "source_bytes");
        e.bytes = inv_list_csv::parse_u64(f[c_bytes], line_no, "serialized_bytes");
        e.cost = inv_list_csv::parse_u64(f[c_cost], line_no, "cost");
        e.count = inv_list_csv::parse_u64(f[c_count], line_no, "count");
        e.memory = inv_list_csv::parse_u64(f[c_memory], line_no, "memory_bytes");
        e.super_sparse = inv_list_csv::parse_u64(f[c_sparse], line_no, "super_sparse") != 0;
        const auto level = inv_list_csv::parse_u64(f[c_serializer_level], line_no, "serializer_level");
        const auto ib = inv_list_csv::parse_u64(f[c_indexed], line_no, "indexed_input_bytes");
        const auto sc = inv_list_csv::parse_u64(f[c_complete], line_no, "source_complete");
        if (sc != 0 && sc != 1) throw std::runtime_error("Invalid pair index metadata");
        if (level != c_level) throw std::runtime_error("Pair index serializer level does not match -level");
        if (first) { indexed_bytes = ib; source_complete = sc != 0; }
        else if (ib != indexed_bytes || (sc != 0) != source_complete)
            throw std::runtime_error("Inconsistent pair index metadata");
        if (first && e.shard != 0)
            throw std::runtime_error("Pair index must begin with shard zero");
        if (!first && (e.shard < previous_shard || e.shard > previous_shard + 1))
            throw std::runtime_error("Non-contiguous shard IDs in pair index");
        if (!first && e.bytes < entries.back().bytes) throw std::runtime_error("Pair index is not size sorted");
        previous_shard = e.shard; first = false;
        if (e.id >= ids.size()) ids.resize(e.id + 1);
        if (ids[e.id]) throw std::runtime_error("Duplicate vector ID in pair index");
        ids[e.id] = true; entries.push_back(e);
    }
    if (entries.empty()) throw std::runtime_error("Pair shard index has no vectors");
    if (ids.size() != entries.size() || std::find(ids.begin(), ids.end(), false) != ids.end())
        throw std::runtime_error("Pair index vector IDs are not dense");
    std::vector<const pair_index_entry*> source_order(entries.size(), nullptr);
    for (const auto& e : entries) source_order[e.id] = &e;
    bm::id64_t expected_offset = 0;
    for (const auto* e : source_order)
    {
        if (e->offset != expected_offset ||
            e->source_bytes > (std::numeric_limits<bm::id64_t>::max)() - expected_offset)
            throw std::runtime_error("Pair index source offsets are not contiguous");
        expected_offset += e->source_bytes;
    }
    if (expected_offset != indexed_bytes)
        throw std::runtime_error("Pair index source byte total is inconsistent");
    return entries;
}

std::unique_ptr<pair_item> pair_load_indexed_item(
    std::ifstream& input, const pair_index_entry& e,
    bm::streams_deserializer<pair_bv>& reader,
    bm::serializer<pair_bv>& serializer,
    bm::serializer<pair_bv>::buffer& buffer)
{
    input.clear(); input.seekg(std::streamoff(e.offset));
    if (!input) throw std::runtime_error("Cannot seek pair input");
    bm::streams_decoder decoder(input);
    std::unique_ptr<pair_item> item(new pair_item);
    if (!reader.deserialize(item->bv, decoder)) throw std::runtime_error("Cannot load indexed bvector");
    if (bm::id64_t(input.tellg()) != e.offset + e.source_bytes)
        throw std::runtime_error("Indexed bvector boundary mismatch");
    item->bv.optimize(); item->cost = pair_cost(bm::calc_complexity(item->bv));
    serializer.serialize(item->bv, buffer); item->bytes = buffer.size();
    if (item->cost != e.cost || item->bytes != e.bytes || item->bv.count() != e.count)
        throw std::runtime_error("Pair index does not match input collection");
    return item;
}

std::vector<pair_index_entry> pair_apply_rs_clusters(
    const std::vector<pair_index_entry>& index, bm::id64_t indexed_bytes)
{
    std::ifstream status(rs_cluster_in + ".report.txt");
    std::string line;
    bool complete = false;
    while (std::getline(status, line))
        if (line == "Status: complete") complete = true;
    if (!complete) throw std::runtime_error("RS cluster manifest report is not complete");
    std::ifstream in(rs_cluster_in);
    if (!in || !std::getline(in, line))
        throw std::runtime_error("Cannot open RS cluster manifest");
    const auto h = inv_list_csv::split(line);
    const auto c_cluster = inv_list_csv::column(h, "cluster_id");
    const auto c_id = inv_list_csv::column(h, "vector_id");
    const auto c_offset = inv_list_csv::column(h, "source_offset");
    const auto c_bytes = inv_list_csv::column(h, "source_bytes");
    const auto c_indexed = inv_list_csv::column(h, "indexed_input_bytes");
    const auto required = (std::max)({c_cluster,c_id,c_offset,c_bytes,c_indexed});
    std::vector<const pair_index_entry*> by_id(index.size(), nullptr);
    for (const auto& e : index) by_id[e.id] = &e;
    std::vector<bool> seen(index.size(), false);
    std::vector<pair_index_entry> selected;
    std::size_t previous_cluster = 0;
    for (std::size_t line_no = 2; std::getline(in, line); ++line_no)
    {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty()) continue;
        const auto f = inv_list_csv::split(line);
        if (f.size() <= required) throw std::runtime_error("Short RS cluster row");
        const auto group = std::size_t(inv_list_csv::parse_u64(f[c_cluster], line_no, "cluster_id"));
        const auto id = std::size_t(inv_list_csv::parse_u64(f[c_id], line_no, "vector_id"));
        const auto offset = inv_list_csv::parse_u64(f[c_offset], line_no, "source_offset");
        const auto bytes = inv_list_csv::parse_u64(f[c_bytes], line_no, "source_bytes");
        const auto total = inv_list_csv::parse_u64(f[c_indexed], line_no, "indexed_input_bytes");
        if (id >= by_id.size() || !by_id[id] || seen[id] || total != indexed_bytes ||
            by_id[id]->offset != offset || by_id[id]->source_bytes != bytes)
            throw std::runtime_error("RS cluster manifest does not match pair index");
        if ((selected.empty() && group != 0) ||
            (!selected.empty() && (group < previous_cluster || group > previous_cluster + 1)))
            throw std::runtime_error("RS cluster IDs must be contiguous and ordered");
        seen[id] = true;
        auto e = *by_id[id]; e.shard = group; selected.push_back(e);
        previous_cluster = group;
    }
    if (selected.size() < 2) throw std::runtime_error("RS cluster manifest has fewer than two vectors");
    return selected;
}

void analyze_pairs_sharded()
{
    if (bv_in_file.empty() || pair_csv.empty() || pair_index_in.empty())
        throw std::runtime_error("Sharded pair search requires -bvin, -csv, and -pair-index-in");
    if (pair_inspect)
        throw std::runtime_error("-pair-inspect is not implemented for sharded search; run -ref-eval or -xor-build for measured serialization");
    const std::string best_name = pair_csv + ".best.csv", report_name = pair_csv + ".report.txt";
    for (const auto& name : std::vector<std::string>{pair_csv,best_name,report_name})
    { std::ifstream f(name.c_str(), std::ios::binary); if (f) throw std::runtime_error("Output already exists: " + name); }
    bool source_complete = false; bm::id64_t indexed_bytes = 0;
    std::vector<pair_index_entry> entries = pair_read_shard_index(source_complete, indexed_bytes);
    const std::size_t indexed_vectors = entries.size();
    if (!rs_cluster_in.empty())
    {
        if (pair_max_ratio != 0)
            throw std::runtime_error("RS cluster search requires -pair-max-ratio 0");
        entries = pair_apply_rs_clusters(entries, indexed_bytes);
    }
    std::ifstream input(bv_in_file.c_str(), std::ios::binary);
    if (!input) throw std::runtime_error("Cannot open pair input");
    input.seekg(0, std::ios::end); const std::streampos file_size = input.tellg(); input.seekg(0);
    if (!input || bm::id64_t(file_size) < indexed_bytes ||
        (source_complete && bm::id64_t(file_size) != indexed_bytes))
        throw std::runtime_error("Pair index input size does not match collection");
    std::vector<std::pair<std::size_t,std::size_t> > shards;
    for (std::size_t begin = 0; begin < entries.size(); )
    {
        std::size_t end = begin + 1;
        while (end < entries.size() && entries[end].shard == entries[begin].shard) ++end;
        shards.push_back(std::make_pair(begin,end)); begin = end;
    }
    bm::id64_t planned_pairs = 0;
    for (const auto& sh : shards)
    {
        if (pair_max_ratio != 0 &&
            pair_size_ratio(entries[sh.first].bytes,
                            entries[sh.second - 1].bytes) > pair_max_ratio)
            throw std::runtime_error(
                "Pair shard exceeds -pair-max-ratio; rebuild the shard index with the requested ratio");
        planned_pairs += pair_choose2(sh.second - sh.first);
    }
    const bm::id64_t all_pairs = pair_choose2(entries.size());

    std::ofstream csv(pair_csv.c_str()), best_csv(best_name.c_str()), report(report_name.c_str());
    csv.exceptions(std::ios::failbit | std::ios::badbit); best_csv.exceptions(std::ios::failbit | std::ios::badbit);
    report.exceptions(std::ios::failbit | std::ios::badbit);
    report << "Input: " << bv_in_file << "\nPair shard index: " << pair_index_in
           << "\nRS cluster manifest: " << (rs_cluster_in.empty() ? "none" : rs_cluster_in)
           << "\nIndexed vectors: " << indexed_vectors << "; searched vectors: " << entries.size()
           << "\nStatus: incomplete\nSerializer level: " << c_level
           << "\nPair evaluation threads: " << pair_threads
           << "\nMinimum structural gain: " << pair_min_gain
           << "; minimum target gain percent: " << pair_min_pct << '\n' << std::flush;
    const char* size_columns = ",lead_bytes,target_bytes,smaller_bytes,larger_bytes,size_ratio";
    csv << "lead_id,target_id,lead_cost,target_cost,xor_cost,gain_cost,gain_pct,xor_count" << size_columns << '\n' << std::setprecision(12);
    best_csv << "lead_id,target_id,lead_cost,target_cost,xor_cost,gain_cost,gain_pct,xor_count" << size_columns << '\n' << std::setprecision(12);
    std::vector<pair_best> best(indexed_vectors); std::vector<bool> participants(indexed_vectors, false);
    bm::id64_t considered = 0, scheduled = 0, completed = 0, qualified = 0, targets = 0;
    bm::id64_t category_pairs[3] = {}, category_qualified[3] = {};
    bm::id64_t bins[5] = {}, ratio_candidates[4] = {};
    bm::id64_t peak_memory = 0;
    const auto start = pair_clock::now(); auto last = start;
    typedef bm::thread_pool<bm::task_descr*, bm::spin_lock<bm::pad0_struct> > pool_type;
    typedef bm::task_batch<pair_bv::allocator_type> batch_type;
    pool_type pool; pool.start(pair_threads);
    const std::size_t pairs_per_task = 128, tasks_per_window = std::size_t(pair_threads) * 4;
    bm::serializer<pair_bv> serializer; serializer.set_compression_level(c_level);
    bm::serializer<pair_bv>::buffer serial_buffer;
    bm::streams_deserializer<pair_bv> indexed_reader;

    auto progress = [&](bool force, std::size_t shard_no)
    {
        const auto now = pair_clock::now(); if (!force && now - last < std::chrono::seconds(10)) return;
        const double seconds = std::chrono::duration<double>(now - start).count();
        const double rate = seconds > 0 ? double(considered) / seconds : 0;
        if (!is_silent) cout << "Shard " << shard_no << '/' << shards.size() << "; pairs " << considered
            << '/' << planned_pairs << " (" << (planned_pairs ? 100.0 * double(considered) / double(planned_pairs) : 100.0)
            << "%), pairs/s=" << rate << ", elapsed_s=" << seconds
            << ", ETA_s=" << (rate > 0 ? double(planned_pairs - considered) / rate : 0)
            << ", scheduled=" << scheduled << ", completed=" << completed
            << ", pending=" << scheduled - completed
            << ", qualifying=" << qualified << ", targets=" << targets << endl;
        last = now;
    };

    for (std::size_t shard_no = 0; shard_no < shards.size(); ++shard_no)
    {
        const std::size_t begin = shards[shard_no].first, end = shards[shard_no].second;
        std::vector<std::unique_ptr<pair_item> > items; items.reserve(end - begin);
        bm::id64_t shard_memory = 0;
        for (std::size_t p = begin; p < end; ++p)
        {
            auto item = pair_load_indexed_item(input, entries[p], indexed_reader,
                                               serializer, serial_buffer);
            pair_bv::statistics st; item->bv.calc_stat(&st); shard_memory += st.memory_used;
            items.push_back(std::move(item));
        }
        peak_memory = (std::max)(peak_memory, shard_memory);
        std::size_t next_i = 0, next_j = 1;
        auto next_pair = [&](pair_coord& coord) -> bool
        {
            while (next_i < items.size())
            {
                if (next_j >= items.size()) { ++next_i; next_j = next_i + 1; continue; }
                const std::size_t i = next_i, j = next_j++; ++considered;
                const auto& ei = entries[begin+i]; const auto& ej = entries[begin+j];
                const unsigned cat = unsigned(ei.super_sparse) + unsigned(ej.super_sparse);
                ++category_pairs[cat];
                coord.first = i; coord.second = j; return true;
            }
            return false;
        };
        bool done = false;
        while (!done)
        {
            std::vector<pair_job> jobs; jobs.reserve(tasks_per_window);
            while (jobs.size() < tasks_per_window && !done)
            {
                jobs.emplace_back(); auto& job = jobs.back(); job.pairs.reserve(pairs_per_task);
                pair_coord coord;
                while (job.pairs.size() < pairs_per_task)
                { if (!next_pair(coord)) { done = true; break; } job.pairs.push_back(coord); }
                if (job.pairs.empty()) jobs.pop_back();
            }
            if (jobs.empty()) break;
            batch_type batch; batch.get_task_vector().reserve(jobs.size());
            for (auto& job : jobs)
            {
                job.matches.reserve(job.pairs.size()); pair_job* jp = &job;
                bm::task_function_t task([jp,&items,&entries,begin](void*)
                {
                    for (const auto& coord : jp->pairs)
                    {
                        const std::size_t i=coord.first,j=coord.second;
                        const bool first_lead = items[i]->cost < items[j]->cost ||
                            (items[i]->cost == items[j]->cost &&
                             entries[begin+i].id < entries[begin+j].id);
                        const std::size_t target=first_lead?j:i;
                        const auto residual=bm::calc_complexity_xor(items[i]->bv,items[j]->bv);
                        const auto cost=pair_cost(residual), baseline=items[target]->cost;
                        if (cost < baseline)
                        {
                            const auto gain=baseline-cost; const double pct=100.0*double(gain)/double(baseline);
                            if (gain>=pair_min_gain && pct>=pair_min_pct)
                                jp->matches.push_back(pair_match{i,j,cost,residual.count});
                        }
                    }
                    return 0;
                });
                batch.add(task,0); scheduled += job.pairs.size();
            }
            auto& queue=pool.get_job_queue();
            for (batch_type::size_type k=0;k<batch.size();++k) queue.push(batch.get_task(k));
            for (std::size_t k=0;k<jobs.size();++k)
            {
                bm::task_descr* task=batch.get_task(batch_type::size_type(k));
                while (!task->done.load(std::memory_order_acquire))
                { progress(false,shard_no+1); std::this_thread::sleep_for(std::chrono::milliseconds(1)); }
                if (task->err_code) throw std::runtime_error("Sharded pair evaluation task failed");
                completed += jobs[k].pairs.size();
                for (const auto& match:jobs[k].matches)
                {
                    const std::size_t li=match.first,lj=match.second;
                    const bool first_lead = items[li]->cost < items[lj]->cost ||
                        (items[li]->cost == items[lj]->cost &&
                         entries[begin+li].id < entries[begin+lj].id);
                    const std::size_t lead_local=first_lead?li:lj;
                    const std::size_t target_local=lead_local==li?lj:li;
                    const auto& le=entries[begin+lead_local]; const auto& te=entries[begin+target_local];
                    const auto baseline=items[target_local]->cost, gain=baseline-match.residual_cost;
                    const double pct=100.0*double(gain)/double(baseline); ++qualified;
                    const unsigned cat=unsigned(le.super_sparse)+unsigned(te.super_sparse); ++category_qualified[cat];
                    participants[le.id]=participants[te.id]=true;
                    csv<<le.id<<','<<te.id<<','<<items[lead_local]->cost<<','<<baseline<<','
                       <<match.residual_cost<<','<<gain<<','<<pct<<','<<match.residual_count;
                    pair_write_sizes(csv,le.bytes,te.bytes); csv<<'\n';
                    ++ratio_candidates[pair_ratio_bucket(pair_size_ratio(le.bytes,te.bytes))];
                    auto& b=best[te.id]; if(!b.found)++targets;
                    if(!b.found||gain>b.gain||(gain==b.gain&&le.id<b.lead))
                    { b.found=true;b.lead=le.id;b.residual=match.residual_cost;b.gain=gain;b.count=match.residual_count; }
                }
                progress(false,shard_no+1);
            }
            csv.flush();
        }
        if (!is_silent) cout << "Completed shard " << shard_no+1 << '/' << shards.size()
                             << "; vectors=" << items.size() << "; allocation MiB="
                             << double(shard_memory)/(1024.0*1024) << endl;
    }
    pool.set_stop_mode(pool_type::stop_when_done); pool.join(); progress(true,shards.size());
    std::vector<const pair_index_entry*> by_id(indexed_vectors, nullptr);
    for (const auto& entry : entries) by_id[entry.id] = &entry;
    for (std::size_t target=0;target<best.size();++target)
    {
        const auto& b=best[target]; if(!b.found)continue;
        const auto& le=*by_id[b.lead]; const auto& te=*by_id[target];
        const double pct=100.0*double(b.gain)/double(te.cost); ++bins[pct<25?0:pct<50?1:pct<75?2:pct<100?3:4];
        best_csv<<b.lead<<','<<target<<','<<le.cost<<','<<te.cost<<','<<b.residual<<','<<b.gain<<','<<pct<<','<<b.count;
        pair_write_sizes(best_csv,le.bytes,te.bytes); best_csv<<'\n';
    }
    csv.close(); best_csv.close();
    auto summary=[&](std::ostream& out)
    {
        out<<"Sharded search complete for "<<(source_complete?"whole collection":"loaded prefix")
           <<". Indexed vectors="<<indexed_vectors<<"; searched vectors="<<entries.size()<<"; shards="<<shards.size()
           <<"\nAll searched-vector pairs="<<all_pairs<<"; internal-shard pairs="<<planned_pairs
           <<"; deliberately omitted cross-shard pairs="<<(all_pairs-planned_pairs)
           <<"; coverage pct="<<(all_pairs?100.0*double(planned_pairs)/double(all_pairs):100.0)
           <<"\nPairs considered="<<considered<<"; completed="<<completed
           <<"\nPair evaluation threads="<<pair_threads
           <<"\nQualifying pairs="<<qualified<<"; candidate targets="<<targets
           <<" ("<<(entries.empty()?0:100.0*double(targets)/double(entries.size()))<<"%); participating vectors="
           <<std::count(participants.begin(),participants.end(),true)
           <<"\nPair categories normal-normal,mixed,sparse-sparse: evaluated="
           <<category_pairs[0]<<','<<category_pairs[1]<<','<<category_pairs[2]
           <<"; qualifying="<<category_qualified[0]<<','<<category_qualified[1]<<','<<category_qualified[2]
           <<"\nBest-target gain distribution [0,25), [25,50), [50,75), [75,100), 100: "
           <<bins[0]<<','<<bins[1]<<','<<bins[2]<<','<<bins[3]<<','<<bins[4]
           <<"\nPeak loaded shard allocation MiB="<<double(peak_memory)/(1024.0*1024)
           <<"\nStructural gains are heuristic units, not serialized bytes. Cross-shard outcomes are unknown.\n"
           <<"Size ratio buckets for qualifying pairs [1,2),[2,4),[4,8),[8,infinity): "
           <<ratio_candidates[0]<<','<<ratio_candidates[1]<<','<<ratio_candidates[2]<<','<<ratio_candidates[3]
           <<"\nElapsed seconds="<<std::chrono::duration<double>(pair_clock::now()-start).count()<<'\n';
    };
    summary(cout);summary(report);report<<"Status: complete\n";
}

} // namespace

#endif
