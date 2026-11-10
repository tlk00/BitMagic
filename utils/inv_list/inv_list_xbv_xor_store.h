/** @file inv_list_xbv_xor_store.h
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

// Experimental XOR collection builder/verifier for inv_list.cpp.
// Included after inv_list_xbv_ref_eval.h; not a public library header.
#ifndef INV_LIST_XBV_XOR_STORE_H
#define INV_LIST_XBV_XOR_STORE_H

namespace
{

const unsigned xor_store_version = 1;
const unsigned xor_store_slots = 4;
const bm::id64_t xor_payload_header_size = 40;
const bm::id64_t xor_refs_header_size = 48;
const bm::id64_t xor_index_header_size = 48;

enum xor_mode { xor_standalone = 0, xor_plain = 1, xor_reference = 2 };

struct xor_plan
{
    xor_mode mode = xor_standalone;
    std::vector<std::size_t> refs;
    bm::id64_t bytes = 0;
};

struct xor_loaded
{
    std::vector<std::unique_ptr<pair_item> > items;
    std::vector<std::vector<ref_candidate> > candidates;
    std::vector<std::size_t> targets;
    bm::id64_t memory = 0, standalone_bytes = 0, source_bytes = 0;
    bm::id64_t total_ints = 0;
    bm::id64_t csv_rows = 0, eligible_rows = 0;
    bool prefix = false;
    bool indexed = false;
};

void xor_write_raw(std::ostream& out, const void* ptr, std::size_t size)
{
    out.write(static_cast<const char*>(ptr), std::streamsize(size));
    if (!out) throw std::runtime_error("Error writing XOR store");
}

void xor_write_u32(std::ostream& out, std::uint32_t value)
{
    unsigned char b[4];
    for (unsigned i = 0; i < 4; ++i) b[i] = (unsigned char)(value >> (i * 8));
    xor_write_raw(out, b, sizeof(b));
}

void xor_write_u64(std::ostream& out, std::uint64_t value)
{
    unsigned char b[8];
    for (unsigned i = 0; i < 8; ++i) b[i] = (unsigned char)(value >> (i * 8));
    xor_write_raw(out, b, sizeof(b));
}

void xor_read_raw(std::istream& in, void* ptr, std::size_t size)
{
    in.read(static_cast<char*>(ptr), std::streamsize(size));
    if (!in) throw std::runtime_error("Truncated XOR store");
}

std::uint32_t xor_read_u32(std::istream& in)
{
    unsigned char b[4]; xor_read_raw(in, b, sizeof(b));
    std::uint32_t v = 0;
    for (unsigned i = 0; i < 4; ++i) v |= std::uint32_t(b[i]) << (i * 8);
    return v;
}

std::uint64_t xor_read_u64(std::istream& in)
{
    unsigned char b[8]; xor_read_raw(in, b, sizeof(b));
    std::uint64_t v = 0;
    for (unsigned i = 0; i < 8; ++i) v |= std::uint64_t(b[i]) << (i * 8);
    return v;
}

void xor_write_magic(std::ostream& out, const char magic[9])
{ xor_write_raw(out, magic, 8); }

void xor_expect_magic(std::istream& in, const char magic[9])
{
    char actual[8]; xor_read_raw(in, actual, sizeof(actual));
    if (std::memcmp(actual, magic, 8))
        throw std::runtime_error("Wrong XOR store file type");
}

bm::id64_t xor_file_size(std::ifstream& in)
{
    in.seekg(0, std::ios::end);
    const std::streampos pos = in.tellg();
    in.seekg(0);
    if (!in || pos < std::streampos(0))
        throw std::runtime_error("Cannot determine file size");
    return bm::id64_t(pos);
}

std::vector<unsigned char> xor_read_blob(std::istream& in, bm::id64_t size)
{
    if (size > bm::id64_t((std::numeric_limits<std::size_t>::max)()))
        throw std::runtime_error("XOR store BLOB is too large for this process");
    std::vector<unsigned char> blob((std::size_t(size)));
    if (size) xor_read_raw(in, blob.data(), std::size_t(size));
    return blob;
}

xor_loaded xor_load_analysis_input()
{
    xor_loaded data;
    std::ifstream input(bv_in_file.c_str(), std::ios::binary);
    if (!input) throw std::runtime_error("Cannot open bvector collection");
    const bm::id64_t file_size = xor_file_size(input);
    if (!file_size) throw std::runtime_error("Empty bvector collection");
    bm::streams_decoder decoder(input);
    bm::streams_deserializer<pair_bv> reader;
    bm::serializer<pair_bv> serializer;
    serializer.set_compression_level(c_level);
    bm::serializer<pair_bv>::buffer buffer;
    const auto start = pair_clock::now();
    auto last = start;
    while (bm::id64_t(input.tellg()) != file_size &&
           (!pair_limit || data.items.size() < pair_limit))
    {
        const std::streampos before = input.tellg();
        std::unique_ptr<pair_item> item(new pair_item);
        if (!reader.deserialize(item->bv, decoder))
            throw std::runtime_error("Error reading bvector collection");
        const std::streampos after = input.tellg();
        if (before < std::streampos(0) || after <= before ||
            bm::id64_t(after) > file_size)
            throw std::runtime_error("Invalid bvector collection boundary");
        item->bv.optimize();
        item->cost = pair_cost(bm::calc_complexity(item->bv));
        serializer.serialize(item->bv, buffer);
        item->bytes = buffer.size();
        data.standalone_bytes += item->bytes;
        data.total_ints += item->bv.count();
        pair_bv::statistics stat; item->bv.calc_stat(&stat);
        data.memory += stat.memory_used;
        data.items.push_back(std::move(item));
        data.source_bytes = bm::id64_t(after);
        const auto now = pair_clock::now();
        if (!is_silent && now - last >= std::chrono::seconds(10))
        {
            cout << "Loaded " << data.items.size() << " vectors; input "
                 << after << '/' << file_size << "; bvector allocation MiB="
                 << double(data.memory) / (1024.0 * 1024) << endl;
            last = now;
        }
    }
    data.prefix = data.source_bytes != file_size;
    cout << "Loaded " << data.items.size() << " vectors ("
         << (data.prefix ? "prefix only" : "whole collection")
         << "); bvector allocation MiB="
         << double(data.memory) / (1024.0 * 1024) << endl;

    std::ifstream pairs(ref_pairs_in.c_str());
    if (!pairs) throw std::runtime_error("Cannot open pair CSV");
    std::string line;
    if (!std::getline(pairs, line)) throw std::runtime_error("Empty pair CSV");
    if (!line.empty() && line.back() == '\r') line.pop_back();
    const auto header = inv_list_csv::split(line);
    const std::size_t lead_id_col = inv_list_csv::column(header, "lead_id");
    const std::size_t target_id_col = inv_list_csv::column(header, "target_id");
    const std::size_t lead_cost_col = inv_list_csv::column(header, "lead_cost");
    const std::size_t target_cost_col = inv_list_csv::column(header, "target_cost");
    const std::size_t xor_cost_col = inv_list_csv::column(header, "xor_cost");
    const std::size_t lead_bytes_col = inv_list_csv::column(header, "lead_bytes");
    const std::size_t target_bytes_col = inv_list_csv::column(header, "target_bytes");
    const std::size_t required = (std::max)(
        (std::max)((std::max)(lead_id_col, target_id_col),
                   (std::max)(lead_cost_col, target_cost_col)),
        (std::max)(xor_cost_col, (std::max)(lead_bytes_col, target_bytes_col)));
    data.candidates.resize(data.items.size());
    for (std::size_t line_no = 2; std::getline(pairs, line); ++line_no)
    {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty()) continue;
        ++data.csv_rows;
        const auto f = inv_list_csv::split(line);
        if (f.size() <= required)
            throw std::runtime_error("Short row in pair CSV at line " + std::to_string(line_no));
        const std::size_t a = std::size_t(inv_list_csv::parse_u64(f[lead_id_col], line_no, "lead_id"));
        const std::size_t b = std::size_t(inv_list_csv::parse_u64(f[target_id_col], line_no, "target_id"));
        if (a >= data.items.size() || b >= data.items.size()) continue;
        const bm::id64_t ac = inv_list_csv::parse_u64(f[lead_cost_col], line_no, "lead_cost");
        const bm::id64_t bc = inv_list_csv::parse_u64(f[target_cost_col], line_no, "target_cost");
        const bm::id64_t xc = inv_list_csv::parse_u64(f[xor_cost_col], line_no, "xor_cost");
        const bm::id64_t ab = inv_list_csv::parse_u64(f[lead_bytes_col], line_no, "lead_bytes");
        const bm::id64_t bb = inv_list_csv::parse_u64(f[target_bytes_col], line_no, "target_bytes");
        if (ac != data.items[a]->cost || bc != data.items[b]->cost ||
            ab != data.items[a]->bytes || bb != data.items[b]->bytes)
            throw std::runtime_error("Pair CSV does not match vectors/serializer at line " + std::to_string(line_no));
        const bool a_is_dep = std::make_pair(ab, a) < std::make_pair(bb, b);
        const std::size_t dep = a_is_dep ? a : b, target = a_is_dep ? b : a;
        const bm::id64_t baseline = data.items[target]->cost;
        const long long gain = baseline >= xc ?
            static_cast<long long>(baseline - xc) :
            -static_cast<long long>(xc - baseline);
        data.candidates[target].push_back(ref_candidate{dep, gain});
        ++data.eligible_rows;
    }
    for (std::size_t target = 0; target < data.candidates.size(); ++target)
    {
        auto& c = data.candidates[target];
        if (c.empty()) continue;
        std::sort(c.begin(), c.end(), [](const ref_candidate& x, const ref_candidate& y)
        { return x.id != y.id ? x.id < y.id : x.structural_gain > y.structural_gain; });
        c.erase(std::unique(c.begin(), c.end(), [](const ref_candidate& x, const ref_candidate& y)
        { return x.id == y.id; }), c.end());
        std::sort(c.begin(), c.end(), [](const ref_candidate& x, const ref_candidate& y)
        { return x.structural_gain != y.structural_gain ? x.structural_gain > y.structural_gain : x.id < y.id; });
        if (c.size() > ref_top) c.resize(ref_top);
        data.targets.push_back(target);
    }
    cout << "Read " << data.csv_rows << " qualifying pair rows; eligible rows="
         << data.eligible_rows << "; reference targets=" << data.targets.size() << endl;
    return data;
}

std::vector<xor_plan> xor_make_plans(const xor_loaded& data)
{
    std::vector<xor_plan> plans(data.items.size());
    for (std::size_t i = 0; i < plans.size(); ++i) plans[i].bytes = data.items[i]->bytes;
    typedef bm::thread_pool<bm::task_descr*, bm::spin_lock<bm::pad0_struct> > pool_type;
    typedef bm::task_batch<pair_bv::allocator_type> batch_type;
    pool_type pool; pool.start(pair_threads);
    const std::size_t window = std::size_t(pair_threads) * 4;
    const auto start = pair_clock::now(); auto last = start;
    std::size_t completed = 0; std::string error;
    for (std::size_t begin = 0; begin < data.targets.size() && error.empty(); begin += window)
    {
        const std::size_t count = (std::min)(window, data.targets.size() - begin);
        std::vector<ref_eval_job> jobs(count); batch_type batch;
        for (std::size_t i = 0; i < count; ++i) jobs[i].target = data.targets[begin + i];
        for (auto& job : jobs)
        {
            ref_eval_job* p = &job;
            bm::task_function_t task([p, &data](void*)
            { ref_evaluate_one(*p, data.items, data.candidates); return p->result.error.empty() ? 0 : 1; });
            batch.add(task, 0);
        }
        auto& queue = pool.get_job_queue();
        for (batch_type::size_type i = 0; i < batch.size(); ++i) queue.push(batch.get_task(i));
        for (std::size_t i = 0; i < count; ++i)
        {
            bm::task_descr* task = batch.get_task(batch_type::size_type(i));
            while (!task->done.load(std::memory_order_acquire))
            {
                const auto now = pair_clock::now();
                if (!is_silent && now - last >= std::chrono::seconds(10))
                {
                    cout << "XOR plans " << completed << '/' << data.targets.size()
                         << "; elapsed_s=" << std::chrono::duration<double>(now - start).count() << endl;
                    last = now;
                }
                std::this_thread::sleep_for(std::chrono::milliseconds(1));
            }
            const ref_eval_result& r = jobs[i].result;
            if (task->err_code || !r.error.empty()) { error = r.error.empty() ? "XOR planning task failed" : r.error; continue; }
            xor_plan& p = plans[r.target];
            if (r.plain_xor_bytes < p.bytes)
            { p.mode = xor_plain; p.bytes = r.plain_xor_bytes; p.refs.assign(1, r.plain_lead); }
            if (r.ref1.bytes < p.bytes)
            { p.mode = xor_reference; p.bytes = r.ref1.bytes; p.refs.assign(1, r.candidates[0]); }
            if (r.refk.bytes < p.bytes)
            { p.mode = xor_reference; p.bytes = r.refk.bytes; p.refs = r.candidates; }
            ++completed;
        }
    }
    pool.set_stop_mode(pool_type::stop_when_done); pool.join();
    if (!error.empty()) throw std::runtime_error(error);
    return plans;
}

xor_loaded xor_load_indexed_analysis_input(std::vector<pair_index_entry>& by_id)
{
    bool complete = false; bm::id64_t indexed_bytes = 0;
    auto index = pair_read_shard_index(complete, indexed_bytes);
    if (pair_limit && pair_limit < index.size())
        throw std::runtime_error("Indexed XOR build requires a full index or a separately indexed prefix");
    by_id.resize(index.size());
    xor_loaded data; data.indexed = true; data.prefix = !complete;
    data.source_bytes = indexed_bytes; data.items.resize(index.size());
    data.candidates.resize(index.size());
    for (const auto& e : index)
    {
        by_id[e.id] = e;
        data.standalone_bytes += e.bytes;
        data.total_ints += e.count;
    }
    std::ifstream pairs(ref_pairs_in.c_str());
    if (!pairs) throw std::runtime_error("Cannot open pair CSV");
    std::string line;
    if (!std::getline(pairs,line)) throw std::runtime_error("Empty pair CSV");
    const auto h = inv_list_csv::split(line);
    const auto c_lead = inv_list_csv::column(h,"lead_id");
    const auto c_target = inv_list_csv::column(h,"target_id");
    const auto c_lcost = inv_list_csv::column(h,"lead_cost");
    const auto c_tcost = inv_list_csv::column(h,"target_cost");
    const auto c_xcost = inv_list_csv::column(h,"xor_cost");
    const auto c_lbytes = inv_list_csv::column(h,"lead_bytes");
    const auto c_tbytes = inv_list_csv::column(h,"target_bytes");
    const auto required = (std::max)({c_lead,c_target,c_lcost,c_tcost,c_xcost,c_lbytes,c_tbytes});
    for (std::size_t line_no=2; std::getline(pairs,line); ++line_no)
    {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty()) continue;
        ++data.csv_rows;
        const auto f = inv_list_csv::split(line);
        if (f.size() <= required) throw std::runtime_error("Short pair CSV row at line " + std::to_string(line_no));
        const auto a = std::size_t(inv_list_csv::parse_u64(f[c_lead],line_no,"lead_id"));
        const auto b = std::size_t(inv_list_csv::parse_u64(f[c_target],line_no,"target_id"));
        if (a >= by_id.size() || b >= by_id.size()) continue;
        const auto ac = inv_list_csv::parse_u64(f[c_lcost],line_no,"lead_cost");
        const auto bc = inv_list_csv::parse_u64(f[c_tcost],line_no,"target_cost");
        const auto xc = inv_list_csv::parse_u64(f[c_xcost],line_no,"xor_cost");
        const auto ab = inv_list_csv::parse_u64(f[c_lbytes],line_no,"lead_bytes");
        const auto bb = inv_list_csv::parse_u64(f[c_tbytes],line_no,"target_bytes");
        if (ac != by_id[a].cost || bc != by_id[b].cost ||
            ab != by_id[a].bytes || bb != by_id[b].bytes)
            throw std::runtime_error("Pair CSV does not match indexed vectors at line " + std::to_string(line_no));
        const bool a_dep = std::make_pair(ab,a) < std::make_pair(bb,b);
        const auto dep = a_dep ? a : b, target = a_dep ? b : a;
        const auto baseline = by_id[target].cost;
        const long long gain = baseline >= xc ?
            static_cast<long long>(baseline-xc) : -static_cast<long long>(xc-baseline);
        data.candidates[target].push_back(ref_candidate{dep,gain});
        ++data.eligible_rows;
    }
    for (std::size_t target=0; target<data.candidates.size(); ++target)
    {
        auto& c = data.candidates[target];
        if (c.empty()) continue;
        std::sort(c.begin(),c.end(),[](const ref_candidate& x, const ref_candidate& y)
        { return x.id != y.id ? x.id < y.id : x.structural_gain > y.structural_gain; });
        c.erase(std::unique(c.begin(),c.end(),[](const ref_candidate& x, const ref_candidate& y)
        { return x.id == y.id; }),c.end());
        std::sort(c.begin(),c.end(),[](const ref_candidate& x, const ref_candidate& y)
        { return x.structural_gain != y.structural_gain ? x.structural_gain > y.structural_gain : x.id < y.id; });
        if (c.size() > ref_top) c.resize(ref_top);
        data.targets.push_back(target);
    }
    cout << "Read " << data.csv_rows << " qualifying pair rows; eligible rows="
         << data.eligible_rows << "; reference targets=" << data.targets.size() << endl;
    return data;
}

std::vector<std::size_t> xor_load_indexed_items(
    xor_loaded& data, const std::vector<pair_index_entry>& by_id,
    std::ifstream& input, bm::streams_deserializer<pair_bv>& reader,
    bm::serializer<pair_bv>& serializer, bm::serializer<pair_bv>::buffer& buffer,
    std::size_t target, const std::vector<std::size_t>& refs)
{
    std::vector<std::size_t> loaded;
    bm::id64_t memory = 0;
    auto load = [&](std::size_t id)
    {
        if (data.items[id]) return;
        if (by_id[id].memory > xor_memory - memory)
            throw std::runtime_error("Indexed XOR working set exceeds -xor-memory");
        data.items[id] = pair_load_indexed_item(input,by_id[id],reader,serializer,buffer);
        loaded.push_back(id); memory += by_id[id].memory;
    };
    load(target);
    for (const auto id : refs) load(id);
    data.memory = (std::max)(data.memory,memory);
    return loaded;
}

std::vector<xor_plan> xor_make_plans_indexed(
    xor_loaded& data, const std::vector<pair_index_entry>& by_id)
{
    std::vector<xor_plan> plans(by_id.size());
    for (std::size_t i=0;i<plans.size();++i) plans[i].bytes=by_id[i].bytes;
    std::ifstream input(bv_in_file.c_str(),std::ios::binary);
    if (!input) throw std::runtime_error("Cannot open indexed XOR input");
    bm::streams_deserializer<pair_bv> reader;
    bm::serializer<pair_bv> serializer; serializer.set_compression_level(c_level);
    bm::serializer<pair_bv>::buffer buffer;
    auto last=pair_clock::now();
    for (std::size_t k=0;k<data.targets.size();++k)
    {
        const auto target=data.targets[k];
        std::vector<std::size_t> refs;
        for (const auto& c : data.candidates[target]) refs.push_back(c.id);
        const auto loaded=xor_load_indexed_items(data,by_id,input,reader,serializer,buffer,target,refs);
        ref_eval_job job; job.target=target;
        ref_evaluate_one(job,data.items,data.candidates);
        const auto& r=job.result;
        if (!r.error.empty()) throw std::runtime_error("Indexed XOR plan target " + std::to_string(target) + ": " + r.error);
        xor_plan& p=plans[target];
        if (r.plain_xor_bytes < p.bytes)
        { p.mode=xor_plain; p.bytes=r.plain_xor_bytes; p.refs.assign(1,r.plain_lead); }
        if (r.ref1.bytes < p.bytes)
        { p.mode=xor_reference; p.bytes=r.ref1.bytes; p.refs.assign(1,r.candidates[0]); }
        if (r.refk.bytes < p.bytes)
        { p.mode=xor_reference; p.bytes=r.refk.bytes; p.refs=r.candidates; }
        for (const auto id : loaded) data.items[id].reset();
        const auto now=pair_clock::now();
        if (!is_silent && now-last >= std::chrono::seconds(10))
        { cout << "Indexed XOR plans " << k+1 << '/' << data.targets.size() << endl; last=now; }
    }
    return plans;
}

void xor_serialize_plan(const xor_loaded& data, std::size_t target,
                        const xor_plan& plan,
                        bm::serializer<pair_bv>::buffer& buffer)
{
    bm::serializer<pair_bv> serializer; serializer.set_compression_level(c_level);
    if (plan.mode == xor_standalone)
        serializer.serialize(data.items[target]->bv, buffer);
    else if (plan.mode == xor_plain)
    {
        pair_bv residual(data.items[target]->bv);
        residual ^= data.items[plan.refs[0]]->bv; residual.optimize();
        serializer.serialize(residual, buffer);
    }
    else
    {
        bm::serializer<pair_bv>::bv_ref_vector_type refs;
        refs.add(&data.items[target]->bv, pair_bv::size_type(target));
        for (std::size_t id : plan.refs)
            refs.add(&data.items[id]->bv, pair_bv::size_type(id));
        bm::serializer<pair_bv>::xor_sim_model_type model;
        serializer.set_ref_vectors(&refs);
        if (!serializer.compute_sim_model(model, refs, bm::xor_sim_params()))
            throw std::runtime_error("Reference model failed during XOR store write");
        serializer.set_sim_model(&model); serializer.set_curr_ref_idx(0);
        serializer.serialize(data.items[target]->bv, buffer);
    }
    if (buffer.size() != plan.bytes)
        throw std::runtime_error("XOR plan size changed during store write");
}

template<class SV>
std::vector<unsigned char> xor_serialize_sv(SV& sv)
{
    BM_DECLARE_TEMP_BLOCK(tb)
    sv.optimize(tb);
    bm::sparse_vector_serial_layout<SV> layout;
    bm::sparse_vector_serialize(sv, layout, tb);
    return std::vector<unsigned char>(layout.data(), layout.data() + layout.size());
}

std::vector<unsigned char> xor_serialize_rsc(sparse_vector_u32& sv)
{
    BM_DECLARE_TEMP_BLOCK(tb)
    rsc_sparse_vector_u32 rsc; rsc.load_from(sv); rsc.optimize(tb);
    bm::sparse_vector_serial_layout<rsc_sparse_vector_u32> layout;
    bm::sparse_vector_serialize(rsc, layout, tb);
    return std::vector<unsigned char>(layout.data(), layout.data() + layout.size());
}

double xor_pct(bm::id64_t part, bm::id64_t whole)
{ return whole ? 100.0 * double(part) / double(whole) : 0.0; }

void build_xor_store()
{
    if (bv_in_file.empty() || ref_pairs_in.empty() || xor_out_file.empty() ||
        xor_refs_file.empty() || xor_index_file.empty())
        throw std::runtime_error("XOR build requires -bvin, -pairs-in, -xorout, -xor-refs, and -xor-index");
    if (!ref_top || ref_top > xor_store_slots)
        throw std::runtime_error("XOR store requires 1 to 4 reference candidates");
    if (!xor_memory) throw std::runtime_error("-xor-memory must be positive");
    const std::string report_name = xor_out_file + ".report.txt";
    for (const auto& name : std::vector<std::string>{xor_out_file, xor_refs_file, xor_index_file, report_name})
    { std::ifstream f(name.c_str(), std::ios::binary); if (f) throw std::runtime_error("Refusing to overwrite output file: " + name); }
    const auto start = pair_clock::now();
    std::vector<pair_index_entry> by_id;
    xor_loaded data = pair_index_in.empty() ? xor_load_analysis_input() :
                      xor_load_indexed_analysis_input(by_id);
    if (data.items.size() > (std::numeric_limits<unsigned>::max)() / 2u)
        throw std::runtime_error("Too many vectors for packed 32-bit reference IDs");
    std::vector<xor_plan> plans = data.indexed ? xor_make_plans_indexed(data,by_id) :
                                              xor_make_plans(data);
    const std::uint64_t build_id = std::uint64_t(std::chrono::high_resolution_clock::now().time_since_epoch().count());

    std::ofstream payload(xor_out_file.c_str(), std::ios::binary);
    if (!payload) throw std::runtime_error("Cannot create XOR payload file");
    xor_write_magic(payload, "BMXORDB1"); xor_write_u32(payload, xor_store_version);
    xor_write_u32(payload, unsigned(xor_payload_header_size)); xor_write_u64(payload, build_id);
    xor_write_u64(payload, data.items.size()); xor_write_u32(payload, c_level);
    xor_write_u32(payload, data.prefix ? 0u : 1u);

    typedef bm::sparse_vector<bm::id64_t, pair_bv> offset_sv_type;
    offset_sv_type offsets;
    std::vector<std::unique_ptr<sparse_vector_u32> > slots;
    for (unsigned i = 0; i < xor_store_slots; ++i)
    { slots.emplace_back(new sparse_vector_u32(bm::use_null)); slots.back()->resize(pair_bv::size_type(data.items.size())); }
    bm::serializer<pair_bv>::buffer buffer;
    bm::id64_t payload_bytes = 0, standalone_payload = 0, plain_payload = 0, ref_payload = 0;
    bm::id64_t standalone_count = 0, plain_count = 0, ref_count = 0, dependency_edges = 0;
    std::array<bm::id64_t, xor_store_slots> slot_values = {{0, 0, 0, 0}};
    offset_sv_type::back_insert_iterator offset_bi = offsets.get_back_inserter();
    offset_bi.add(xor_payload_header_size);
    auto last = pair_clock::now();
    std::ifstream indexed_input;
    bm::streams_deserializer<pair_bv> indexed_reader;
    bm::serializer<pair_bv> indexed_serializer;
    indexed_serializer.set_compression_level(c_level);
    bm::serializer<pair_bv>::buffer indexed_buffer;
    if (data.indexed)
    {
        indexed_input.open(bv_in_file.c_str(),std::ios::binary);
        if (!indexed_input) throw std::runtime_error("Cannot reopen indexed XOR input");
    }
    for (std::size_t i = 0; i < plans.size(); ++i)
    {
        std::vector<std::size_t> loaded;
        if (data.indexed)
            loaded=xor_load_indexed_items(data,by_id,indexed_input,indexed_reader,
                                          indexed_serializer,indexed_buffer,i,plans[i].refs);
        if (data.indexed && plans[i].mode == xor_standalone)
        {
            buffer = indexed_buffer;
            if (buffer.size() != plans[i].bytes)
                throw std::runtime_error("Indexed standalone size changed during store write");
        }
        else
            xor_serialize_plan(data, i, plans[i], buffer);
        for (const auto id : loaded) data.items[id].reset();
        xor_write_raw(payload, buffer.data(), buffer.size()); payload_bytes += buffer.size();
        offset_bi.add(xor_payload_header_size + payload_bytes);
        if (plans[i].mode == xor_standalone)
        { ++standalone_count; standalone_payload += buffer.size(); }
        else if (plans[i].mode == xor_plain)
        {
            ++plain_count; plain_payload += buffer.size(); ++dependency_edges;
            slots[0]->set(pair_bv::size_type(i), unsigned(plans[i].refs[0] << 1));
            ++slot_values[0];
        }
        else
        {
            ++ref_count; ref_payload += buffer.size(); dependency_edges += plans[i].refs.size();
            slots[0]->set(pair_bv::size_type(i), unsigned((plans[i].refs[0] << 1) | 1u));
            ++slot_values[0];
            for (std::size_t k = 1; k < plans[i].refs.size(); ++k)
            {
                slots[k]->set(pair_bv::size_type(i), unsigned(plans[i].refs[k]));
                ++slot_values[k];
            }
        }
        const auto now = pair_clock::now();
        if (!is_silent && now - last >= std::chrono::seconds(10))
        { cout << "Wrote XOR payloads " << (i + 1) << '/' << plans.size() << endl; last = now; }
    }
    offset_bi.flush(); payload.close();
    if (!payload) throw std::runtime_error("Error closing XOR payload file");
    const bm::id64_t payload_file_bytes = xor_payload_header_size + payload_bytes;

    std::array<std::vector<unsigned char>, xor_store_slots> slot_blobs;
    for (unsigned i = 0; i < xor_store_slots; ++i) slot_blobs[i] = xor_serialize_rsc(*slots[i]);
    std::ofstream refs(xor_refs_file.c_str(), std::ios::binary);
    if (!refs) throw std::runtime_error("Cannot create XOR reference file");
    xor_write_magic(refs, "BMXREFS1"); xor_write_u32(refs, xor_store_version);
    xor_write_u32(refs, unsigned(xor_refs_header_size)); xor_write_u64(refs, build_id);
    xor_write_u64(refs, data.items.size()); xor_write_u32(refs, xor_store_slots);
    xor_write_u32(refs, c_level); xor_write_u64(refs, payload_file_bytes);
    for (const auto& blob : slot_blobs)
    { xor_write_u64(refs, blob.size()); xor_write_raw(refs, blob.data(), blob.size()); }
    refs.close();
    if (!refs) throw std::runtime_error("Error closing XOR reference file");
    bm::id64_t refs_file_bytes = xor_refs_header_size + xor_store_slots * 8;
    for (const auto& blob : slot_blobs) refs_file_bytes += blob.size();

    std::vector<unsigned char> offset_blob = xor_serialize_sv(offsets);
    std::ofstream index(xor_index_file.c_str(), std::ios::binary);
    if (!index) throw std::runtime_error("Cannot create XOR index file");
    xor_write_magic(index, "BMXOFFS1"); xor_write_u32(index, xor_store_version);
    xor_write_u32(index, unsigned(xor_index_header_size)); xor_write_u64(index, build_id);
    xor_write_u64(index, data.items.size()); xor_write_u64(index, payload_file_bytes);
    xor_write_u64(index, offset_blob.size()); xor_write_raw(index, offset_blob.data(), offset_blob.size());
    index.close();
    if (!index) throw std::runtime_error("Error closing XOR offset index");
    const bm::id64_t index_file_bytes = xor_index_header_size + offset_blob.size();
    const bm::id64_t store_bytes = payload_file_bytes + refs_file_bytes + index_file_bytes;
    const long long net_gain = ref_signed_delta(data.standalone_bytes, store_bytes);
    const long long payload_gain = ref_signed_delta(data.standalone_bytes, payload_file_bytes);
    const long long source_gain = ref_signed_delta(data.source_bytes, store_bytes);
    std::ofstream report(report_name.c_str());
    if (!report) throw std::runtime_error("Cannot create XOR report file");
    auto summary = [&](std::ostream& out)
    {
        out << std::setprecision(12)
            << "XOR store build complete for " << (data.prefix ? "loaded prefix" : "whole collection")
            << ". Vectors=" << data.items.size() << "; serializer level=" << c_level
            << "; candidate limit=" << ref_top << "; evaluation threads=" << (data.indexed ? 1u : pair_threads)
            << "\nPair CSV input=" << ref_pairs_in
            << "\nPair index input=" << (data.indexed ? pair_index_in : "none")
            << "\nInput bytes consumed=" << data.source_bytes
            << "; standalone level-" << c_level << " bytes=" << data.standalone_bytes
            << "\nTotal ints=" << data.total_ints
            << "; input CM6 bits/int=" << (data.total_ints ? 8.0 * double(data.source_bytes) / double(data.total_ints) : 0)
            << "; standalone level-" << c_level << " bits/int="
            << (data.total_ints ? 8.0 * double(data.standalone_bytes) / double(data.total_ints) : 0)
            << "\nPayload file bytes=" << payload_file_bytes << " (header=" << xor_payload_header_size
            << ", BLOBs=" << payload_bytes << ")"
            << "\nPayload modes: standalone=" << standalone_count << " vectors/" << standalone_payload
            << " bytes; plain XOR=" << plain_count << " vectors/" << plain_payload
            << " bytes; reference group=" << ref_count << " vectors/" << ref_payload << " bytes"
            << "\nReference metadata bytes=" << refs_file_bytes << " (header+lengths="
            << (xor_refs_header_size + xor_store_slots * 8) << ", RSC BLOBs="
            << (refs_file_bytes - xor_refs_header_size - xor_store_slots * 8) << ")"
            << "\nOffset index bytes=" << index_file_bytes << " (header=" << xor_index_header_size
            << ", sparse-vector BLOB=" << offset_blob.size()
            << ", ratio vs dense uint64="
            << (data.items.empty() ? 0 : double(offset_blob.size()) / double((data.items.size() + 1) * sizeof(bm::id64_t))) << ")";
        for (unsigned i = 0; i < xor_store_slots; ++i)
            out << "\nReference RSC slot " << i << ": values=" << slot_values[i]
                << "; BLOB bytes=" << slot_blobs[i].size()
                << "; ratio vs dense uint32="
                << (data.items.empty() ? 0 : double(slot_blobs[i].size()) / double(data.items.size() * sizeof(unsigned)));
        out
            << "\nDependency edges=" << dependency_edges
            << "\nTotal XOR store bytes=" << store_bytes
            << "; XOR store bits/int=" << (data.total_ints ? 8.0 * double(store_bytes) / double(data.total_ints) : 0)
            << "; compression ratio vs standalone=" << (data.standalone_bytes ? double(store_bytes) / double(data.standalone_bytes) : 0)
            << "; net saving bytes=" << net_gain << "; net saving pct="
            << (data.standalone_bytes ? 100.0 * double(net_gain) / double(data.standalone_bytes) : 0)
            << "\nSize ratio vs input CM6=" << (data.source_bytes ? double(store_bytes) / double(data.source_bytes) : 0)
            << "; compression factor input/store=" << (store_bytes ? double(data.source_bytes) / double(store_bytes) : 0)
            << "; saving vs input bytes=" << source_gain << "; saving vs input pct="
            << (data.source_bytes ? 100.0 * double(source_gain) / double(data.source_bytes) : 0)
            << "\nPayload-only saving bytes=" << payload_gain << "; payload-only saving pct="
            << (data.standalone_bytes ? 100.0 * double(payload_gain) / double(data.standalone_bytes) : 0)
            << "\nMetadata+index overhead bytes=" << (refs_file_bytes + index_file_bytes)
            << "; overhead pct of standalone=" << xor_pct(refs_file_bytes + index_file_bytes, data.standalone_bytes)
            << "; overhead pct of payload saving=" << (payload_gain > 0 ? 100.0 * double(refs_file_bytes + index_file_bytes) / double(payload_gain) : 0)
            << "\nPair CSV rows=" << data.csv_rows << "; eligible rows=" << data.eligible_rows
            << "\nBvector allocation MiB=" << double(data.memory) / (1024.0 * 1024)
            << (data.indexed ? " (peak indexed working set)" : " (all loaded vectors)")
            << "\nElapsed seconds=" << std::chrono::duration<double>(pair_clock::now() - start).count() << '\n';
    };
    summary(cout); summary(report); report << "Status: complete\n";
}

struct xor_store_info
{
    std::uint64_t build_id = 0, vectors = 0, payload_bytes = 0;
    unsigned level = 0, flags = 0;
};

void verify_xor_store()
{
    if (bv_in_file.empty() || xor_in_file.empty() || xor_refs_file.empty() || xor_index_file.empty())
        throw std::runtime_error("XOR verification requires -bvin, -xorin, -xor-refs, and -xor-index");
    const auto start = pair_clock::now();
    std::ifstream payload(xor_in_file.c_str(), std::ios::binary);
    if (!payload) throw std::runtime_error("Cannot open XOR payload file");
    const bm::id64_t payload_file_size = xor_file_size(payload);
    xor_expect_magic(payload, "BMXORDB1");
    if (xor_read_u32(payload) != xor_store_version || xor_read_u32(payload) != xor_payload_header_size)
        throw std::runtime_error("Unsupported XOR payload version/header");
    xor_store_info info; info.build_id = xor_read_u64(payload); info.vectors = xor_read_u64(payload);
    info.level = xor_read_u32(payload); info.flags = xor_read_u32(payload); info.payload_bytes = payload_file_size;
    if (!info.vectors || info.vectors > std::uint64_t((std::numeric_limits<std::size_t>::max)()))
        throw std::runtime_error("Invalid XOR vector count");
    const std::size_t n = std::size_t(info.vectors);
    if (n > (std::numeric_limits<unsigned>::max)() / 2u)
        throw std::runtime_error("XOR vector count exceeds metadata address range");

    std::ifstream refs_file(xor_refs_file.c_str(), std::ios::binary);
    if (!refs_file) throw std::runtime_error("Cannot open XOR reference file");
    const bm::id64_t refs_size = xor_file_size(refs_file);
    xor_expect_magic(refs_file, "BMXREFS1");
    if (xor_read_u32(refs_file) != xor_store_version || xor_read_u32(refs_file) != xor_refs_header_size ||
        xor_read_u64(refs_file) != info.build_id || xor_read_u64(refs_file) != info.vectors ||
        xor_read_u32(refs_file) != xor_store_slots || xor_read_u32(refs_file) != info.level ||
        xor_read_u64(refs_file) != payload_file_size)
        throw std::runtime_error("XOR reference file does not match payload");
    std::vector<std::unique_ptr<rsc_sparse_vector_u32> > slots;
    for (unsigned i = 0; i < xor_store_slots; ++i)
    {
        const bm::id64_t len = xor_read_u64(refs_file); auto blob = xor_read_blob(refs_file, len);
        slots.emplace_back(new rsc_sparse_vector_u32);
        BM_DECLARE_TEMP_BLOCK(tb)
        if (bm::sparse_vector_deserialize(*slots.back(), blob.data(), tb) != 0)
            throw std::runtime_error("Cannot deserialize XOR reference metadata");
        // RSC serialization may trim a nullable tail. The collection header
        // carries the logical vector count; missing tail entries are NULL.
        if (slots.back()->size() > n) throw std::runtime_error("Wrong XOR reference vector size");
    }
    if (bm::id64_t(refs_file.tellg()) != refs_size) throw std::runtime_error("Extra data in XOR reference file");

    typedef bm::sparse_vector<bm::id64_t, pair_bv> offset_sv_type;
    std::ifstream index_file(xor_index_file.c_str(), std::ios::binary);
    if (!index_file) throw std::runtime_error("Cannot open XOR offset index");
    const bm::id64_t index_size = xor_file_size(index_file);
    xor_expect_magic(index_file, "BMXOFFS1");
    if (xor_read_u32(index_file) != xor_store_version || xor_read_u32(index_file) != xor_index_header_size ||
        xor_read_u64(index_file) != info.build_id || xor_read_u64(index_file) != info.vectors ||
        xor_read_u64(index_file) != payload_file_size)
        throw std::runtime_error("XOR offset index does not match payload");
    const bm::id64_t offset_len = xor_read_u64(index_file); auto offset_blob = xor_read_blob(index_file, offset_len);
    if (bm::id64_t(index_file.tellg()) != index_size) throw std::runtime_error("Extra data in XOR offset index");
    offset_sv_type offsets;
    BM_DECLARE_TEMP_BLOCK(tb)
    if (bm::sparse_vector_deserialize(offsets, offset_blob.data(), tb) != 0 || offsets.size() != n + 1)
        throw std::runtime_error("Cannot deserialize XOR offsets");
    if (offsets[0] != xor_payload_header_size || offsets[pair_bv::size_type(n)] != payload_file_size)
        throw std::runtime_error("Invalid XOR payload boundaries");

    std::vector<std::vector<std::size_t> > deps(n), followers(n);
    std::vector<xor_mode> modes(n, xor_standalone);
    std::vector<unsigned> indegree(n, 0);
    bm::id64_t standalone_count = 0, plain_count = 0, ref_count = 0, edges = 0;
    for (std::size_t i = 0; i < n; ++i)
    {
        unsigned packed = 0;
        if (!slots[0]->try_get(pair_bv::size_type(i), packed)) { ++standalone_count; continue; }
        const bool is_ref = (packed & 1u) != 0; const std::size_t first = packed >> 1;
        modes[i] = is_ref ? xor_reference : xor_plain;
        deps[i].push_back(first);
        if (is_ref)
        {
            ++ref_count;
            for (unsigned k = 1; k < xor_store_slots; ++k)
            { unsigned id = 0; if (slots[k]->try_get(pair_bv::size_type(i), id)) deps[i].push_back(id); }
        }
        else
        {
            ++plain_count;
            for (unsigned k = 1; k < xor_store_slots; ++k)
            { unsigned unused = 0; if (slots[k]->try_get(pair_bv::size_type(i), unused)) throw std::runtime_error("Plain XOR has extra reference slots"); }
        }
        std::sort(deps[i].begin(), deps[i].end());
        if (std::unique(deps[i].begin(), deps[i].end()) != deps[i].end())
            throw std::runtime_error("Duplicate XOR dependency");
        for (std::size_t dep : deps[i])
        {
            if (dep >= n || dep == i) throw std::runtime_error("Invalid XOR dependency ID");
            followers[dep].push_back(i); ++indegree[i]; ++edges;
        }
    }
    std::queue<std::size_t> ready;
    for (std::size_t i = 0; i < n; ++i) if (!indegree[i]) ready.push(i);
    std::vector<std::size_t> order; order.reserve(n);
    while (!ready.empty())
    {
        const std::size_t id = ready.front(); ready.pop(); order.push_back(id);
        for (std::size_t f : followers[id]) if (!--indegree[f]) ready.push(f);
    }
    if (order.size() != n) throw std::runtime_error("Cycle in XOR dependency graph");

    std::ifstream oracle_file(bv_in_file.c_str(), std::ios::binary);
    if (!oracle_file) throw std::runtime_error("Cannot open oracle collection");
    const bm::id64_t oracle_size = xor_file_size(oracle_file);
    bm::streams_decoder oracle_input(oracle_file); bm::streams_deserializer<pair_bv> oracle_reader;
    std::vector<std::unique_ptr<pair_bv> > decoded(n);
    auto last = pair_clock::now(); bm::id64_t decoded_blob_bytes = 0;
    bm::id64_t total_ints = 0;
    bm::id64_t expected_oracle_bytes = 0;
    if (pair_index_in.empty())
    {
        for (std::size_t count = 0; count < order.size(); ++count)
        {
            const std::size_t id = order[count];
            const bm::id64_t begin = offsets[pair_bv::size_type(id)];
            const bm::id64_t end = offsets[pair_bv::size_type(id + 1)];
            if (begin < xor_payload_header_size || end <= begin || end > payload_file_size)
                throw std::runtime_error("Invalid XOR BLOB offset range");
            payload.seekg(std::streamoff(begin)); if (!payload) throw std::runtime_error("Cannot seek XOR payload");
            auto blob = xor_read_blob(payload, end - begin); decoded_blob_bytes += blob.size();
            decoded[id].reset(new pair_bv);
            std::size_t consumed = 0;
            if (modes[id] == xor_reference)
            {
                bm::serializer<pair_bv>::bv_ref_vector_type refv;
                for (std::size_t dep : deps[id]) refv.add(decoded[dep].get(), pair_bv::size_type(dep));
                consumed = bm::deserialize(*decoded[id], blob.data(), 0, &refv);
            }
            else
            {
                consumed = bm::deserialize(*decoded[id], blob.data());
                if (modes[id] == xor_plain) *decoded[id] ^= *decoded[deps[id][0]];
            }
            if (consumed != blob.size()) throw std::runtime_error("XOR BLOB boundary mismatch");
            const auto now = pair_clock::now();
            if (!is_silent && now - last >= std::chrono::seconds(10))
            { cout << "Decoded XOR vectors " << (count + 1) << '/' << n << endl; last = now; }
        }
        for (std::size_t i = 0; i < n; ++i)
        {
            pair_bv oracle;
            if (!oracle_reader.deserialize(oracle, oracle_input)) throw std::runtime_error("Cannot read oracle vector");
            total_ints += oracle.count();
            if (!oracle.equal(*decoded[i]))
                throw std::runtime_error("XOR verification mismatch at vector " + std::to_string(i));
            const auto now = pair_clock::now();
            if (!is_silent && now - last >= std::chrono::seconds(10))
            { cout << "Compared oracle vectors " << (i + 1) << '/' << n << endl; last = now; }
        }
    }
    else
    {
        if (!xor_memory) throw std::runtime_error("-xor-memory must be positive");
        bool source_complete = false; bm::id64_t indexed_bytes = 0;
        auto index = pair_read_shard_index(source_complete,indexed_bytes);
        if (index.size() != n || indexed_bytes > oracle_size ||
            source_complete != bool(info.flags & 1u))
            throw std::runtime_error("XOR offset index does not match oracle collection");
        expected_oracle_bytes = indexed_bytes;
        std::vector<pair_index_entry> by_id(n);
        for (const auto& e : index) by_id[e.id] = e;
        std::vector<unsigned char> active(n,0);
        std::vector<std::size_t> loaded;
        bm::id64_t active_memory = 0;
        std::function<void(std::size_t)> decode = [&](std::size_t id)
        {
            if (active[id] == 2) return;
            if (active[id] == 1) throw std::runtime_error("Cycle in XOR decode");
            active[id] = 1;
            for (const auto dep : deps[id]) decode(dep);
            if (by_id[id].memory > xor_memory - active_memory)
                throw std::runtime_error("XOR dependency closure exceeds -xor-memory at vector " + std::to_string(id));
            const auto begin = offsets[pair_bv::size_type(id)];
            const auto end = offsets[pair_bv::size_type(id+1)];
            if (begin < xor_payload_header_size || end <= begin || end > payload_file_size)
                throw std::runtime_error("Invalid XOR BLOB offset range");
            payload.clear(); payload.seekg(std::streamoff(begin));
            if (!payload) throw std::runtime_error("Cannot seek XOR payload");
            auto blob=xor_read_blob(payload,end-begin); decoded_blob_bytes+=blob.size();
            decoded[id].reset(new pair_bv);
            std::size_t consumed=0;
            if (modes[id] == xor_reference)
            {
                bm::serializer<pair_bv>::bv_ref_vector_type refv;
                for (const auto dep : deps[id]) refv.add(decoded[dep].get(),pair_bv::size_type(dep));
                consumed=bm::deserialize(*decoded[id],blob.data(),0,&refv);
            }
            else
            {
                consumed=bm::deserialize(*decoded[id],blob.data());
                if (modes[id] == xor_plain) *decoded[id] ^= *decoded[deps[id][0]];
            }
            if (consumed != blob.size()) throw std::runtime_error("XOR BLOB boundary mismatch");
            active[id]=2; loaded.push_back(id); active_memory+=by_id[id].memory;
        };
        for (std::size_t i=0;i<n;++i)
        {
            decode(i);
            pair_bv oracle;
            if (!oracle_reader.deserialize(oracle,oracle_input))
                throw std::runtime_error("Cannot read oracle vector");
            total_ints += oracle.count();
            if (!oracle.equal(*decoded[i]))
                throw std::runtime_error("XOR verification mismatch at vector " + std::to_string(i));
            for (const auto id : loaded) { decoded[id].reset(); active[id]=0; }
            loaded.clear(); active_memory=0;
            const auto now=pair_clock::now();
            if (!is_silent && now-last >= std::chrono::seconds(10))
            { cout << "Verified indexed XOR vectors " << i+1 << '/' << n << endl; last=now; }
        }
    }
    const bool oracle_consumed = bm::id64_t(oracle_file.tellg()) == oracle_size;
    const bm::id64_t oracle_prefix_bytes = bm::id64_t(oracle_file.tellg());
    if (expected_oracle_bytes && oracle_prefix_bytes != expected_oracle_bytes)
        throw std::runtime_error("XOR verification did not consume indexed oracle prefix");
    if ((info.flags & 1u) && !oracle_consumed) throw std::runtime_error("Whole-collection XOR store has an oracle tail");
    const bm::id64_t total_store = payload_file_size + refs_size + index_size;
    cout << std::setprecision(12)
         << "XOR verification complete. Vectors=" << n << "; exact matches=" << n
         << "; source=" << ((info.flags & 1u) ? "whole collection" : "prefix")
         << "\nModes: standalone=" << standalone_count << "; plain XOR=" << plain_count
         << "; reference group=" << ref_count << "; dependency edges=" << edges
         << "\nPayload bytes=" << payload_file_size << "; reference metadata bytes=" << refs_size
         << "; offset index bytes=" << index_size << "; total store bytes=" << total_store
         << "\nOracle CM6 bytes=" << oracle_prefix_bytes
         << "; total ints=" << total_ints
         << "; input CM6 bits/int=" << (total_ints ? 8.0 * double(oracle_prefix_bytes) / double(total_ints) : 0)
         << "; XOR store bits/int=" << (total_ints ? 8.0 * double(total_store) / double(total_ints) : 0)
         << "; size ratio XOR/oracle=" << (oracle_prefix_bytes ? double(total_store) / double(oracle_prefix_bytes) : 0)
         << "; compression factor oracle/XOR=" << (total_store ? double(oracle_prefix_bytes) / double(total_store) : 0)
         << "\nDecoded payload BLOB bytes=" << decoded_blob_bytes
         << (pair_index_in.empty() ? "" : " (including dependency reloads)")
         << "; elapsed seconds=" << std::chrono::duration<double>(pair_clock::now() - start).count() << endl;
}

} // namespace

#endif
