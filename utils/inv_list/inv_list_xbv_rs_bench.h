/** @file inv_list_xbv_rs_bench.h
    @brief Cross-bvector (xbv) research support for inv_list.
    @note Development-only; included by inv_list.cpp after shared options and I/O.
*/
/* Copyright(c) 2026 Anatoliy Kuznetsov
   Licensed under the Apache License, Version 2.0.
   http://www.apache.org/licenses/LICENSE-2.0 */
#ifndef INV_LIST_XBV_RS_BENCH_H
#define INV_LIST_XBV_RS_BENCH_H

#include <cassert>
#include <cstdint>
#include <limits>
#include <mutex>
#include <condition_variable>

namespace
{
double rs_seconds(pair_clock::time_point start)
{ return std::chrono::duration<double>(pair_clock::now() - start).count(); }

bm::id64_t rs_abs(bm::id64_t a, bm::id64_t b)
{ return a > b ? a - b : b - a; }

unsigned rs_bucket(bm::id64_t n)
{ unsigned b = 0; while (n >>= 1) ++b; return b; }

struct rs_profile
{
    std::vector<std::uint16_t> regions; // three disjoint region counts per block
    bm::id64_t count = 0;
};

rs_profile rs_extract(const pair_bv& bv)
{
    pair_bv::rs_index_type idx;
    bv.build_rs_index(&idx);
    rs_profile p; p.count = idx.count();
    pair_bv::size_type last = 0;
    if (!bv.find_reverse(last)) return p;
    const unsigned blocks = unsigned(last / bm::gap_max_bits) + 1;
    p.regions.reserve(std::size_t(blocks) * 3);
    bm::id64_t total = 0;
    for (unsigned b = 0; b < blocks; ++b)
    {
        const auto packed = idx.sub_count(b);
        const unsigned a = unsigned(packed & 65535u);
        const unsigned c = unsigned((packed >> 16) & 65535u);
        const unsigned count = idx.count(b);
        if (a + c > count) throw std::runtime_error("Invalid RS regional populations");
        const unsigned third = count - a - c;
        if (third > (std::numeric_limits<std::uint16_t>::max)())
            throw std::runtime_error("RS regional population exceeds 16 bits");
        p.regions.push_back(std::uint16_t(a));
        p.regions.push_back(std::uint16_t(c));
        p.regions.push_back(std::uint16_t(third)); total += count;
    }
    if (total != p.count) throw std::runtime_error("RS profile count mismatch");
    return p;
}

// This is a heuristic compression screen, not a bound on compression cost.
// The accumulated value IS a lower bound on Hamming distance.
bool rs_accept(const rs_profile& a, const rs_profile& b, double threshold)
{
    const double limit = threshold * double(a.count + b.count);
    if (double(rs_abs(a.count, b.count)) > limit) return false;
    bm::id64_t distance = 0;
    const std::size_t a_size = a.regions.size(), b_size = b.regions.size();
    const std::size_t common = (std::min)(a_size, b_size);
    std::size_t r = 0;
    for (; r < common; r += 3)
    {
        distance += rs_abs(a.regions[r], b.regions[r]) +
                    rs_abs(a.regions[r+1], b.regions[r+1]) +
                    rs_abs(a.regions[r+2], b.regions[r+2]);
        if (double(distance) > limit) return false;
    }
    if (a_size != b_size)
    {
        const auto& tail = a_size > b_size ? a.regions : b.regions;
        for (; r < tail.size(); r += 3)
        {
            distance += unsigned(tail[r]) + unsigned(tail[r+1]) + unsigned(tail[r+2]);
            if (double(distance) > limit) return false;
        }
    }
    return true;
}

struct rs_distances
{
    bm::id64_t block_l1 = 0, region_l1 = 0;
    double shape = 0;
};

void rs_distance(rs_distances* result, const rs_profile& a, const rs_profile& b)
{
    assert(result);
    const std::size_t a_size = a.regions.size(), b_size = b.regions.size();
    if (a_size % 3 || b_size % 3)
        throw std::runtime_error("RS profile has an incomplete block");
    auto& d = *result;
    d = {};
    const bool has_shape = a.count && b.count;
    const std::size_t common = (std::min)(a_size, b_size);
    std::size_t r = 0;
    for (; r < common; r += 3)
    {
        const unsigned a0 = a.regions[r], a1 = a.regions[r+1], a2 = a.regions[r+2];
        const unsigned b0 = b.regions[r], b1 = b.regions[r+1], b2 = b.regions[r+2];
        d.region_l1 += rs_abs(a0, b0) + rs_abs(a1, b1) + rs_abs(a2, b2);
        d.block_l1 += rs_abs(a0 + a1 + a2, b0 + b1 + b2);
        if (has_shape)
        {
            d.shape += std::abs(double(a0)/double(a.count) - double(b0)/double(b.count));
            d.shape += std::abs(double(a1)/double(a.count) - double(b1)/double(b.count));
            d.shape += std::abs(double(a2)/double(a.count) - double(b2)/double(b.count));
        }
    }
    if (a_size != b_size)
    {
        const bool a_tail = a_size > b_size;
        const auto& tail = a_tail ? a.regions : b.regions;
        const double count = double(a_tail ? a.count : b.count);
        for (; r < tail.size(); r += 3)
        {
            const unsigned t0 = tail[r], t1 = tail[r+1], t2 = tail[r+2];
            const unsigned block_count = t0 + t1 + t2;
            d.block_l1 += block_count;
            d.region_l1 += block_count;
            if (has_shape)
            {
                d.shape += double(t0)/count;
                d.shape += double(t1)/count;
                d.shape += double(t2)/count;
            }
        }
    }
    d.shape = has_shape ? d.shape / 2 : -1; // undefined for empty vectors
}

// Join all workers before captured data can go out of scope, including failures.
template<class Function>
void rs_parallel(std::size_t begin, std::size_t end, Function fn,
                 const std::string& label = std::string())
{
    std::atomic<std::size_t> next(begin), done(0);
    std::atomic<bool> failed(false);
    std::mutex error_mutex; std::exception_ptr error;
    std::condition_variable finished;
    unsigned workers_left = pair_threads;
    auto worker = [&]()
    {
        try
        {
            while (!failed.load())
            {
                const auto row = next.fetch_add(1);
                if (row >= end) break;
                fn(row); ++done;
            }
        }
        catch (...)
        {
            std::lock_guard<std::mutex> guard(error_mutex);
            if (!error) error = std::current_exception();
            failed = true;
        }
        {
            std::lock_guard<std::mutex> guard(error_mutex);
            --workers_left;
        }
        finished.notify_one();
    };
    std::vector<std::thread> workers;
    try
    {
        for (unsigned t = 0; t < pair_threads; ++t) workers.emplace_back(worker);
    }
    catch (...)
    { failed = true; for (auto& t : workers) t.join(); throw; }
    auto last = pair_clock::now();
    for (;;)
    {
        std::unique_lock<std::mutex> lock(error_mutex);
        if (finished.wait_for(lock, std::chrono::seconds(1), [&](){return workers_left == 0;})) break;
        lock.unlock();
        if (!label.empty() && !is_silent && rs_seconds(last) >= 10)
        { cout << label << ": rows " << done.load() << '/' << end-begin << endl; last = pair_clock::now(); }
    }
    for (auto& t : workers) t.join();
    if (error) std::rethrow_exception(error);
}

struct rs_pair_result
{
    bm::id64_t cost = 0, count = 0;
    rs_distances d;
};

bool rs_qualifies(bm::id64_t cost, bm::id64_t baseline)
{
    return cost < baseline && baseline-cost >= pair_min_gain &&
           100.0 * double(baseline-cost) / double(baseline) >= pair_min_pct;
}

struct rs_totals
{
    bm::id64_t evaluated = 0, qualified = 0, checksum_cost = 0, checksum_count = 0;
    void add(const rs_totals& b)
    { evaluated += b.evaluated; qualified += b.qualified;
      checksum_cost += b.checksum_cost; checksum_count += b.checksum_count; }
};

void run_rs_benchmark()
{
    if (rs_out.empty() || bv_in_file.empty() || pair_index_in.empty())
        throw std::runtime_error("RS benchmark requires -rs-out, -bvin, -pair-index-in");
    if (rs_set != "shard" && rs_set != "stratified")
        throw std::runtime_error("-rs-set must be shard or stratified");
    if (rs_size < 2 || rs_repeats < 1 || rs_repeats > 100)
        throw std::runtime_error("RS size must be >=2; repeats must be 1..100");
    std::vector<double> thresholds;
    for (const auto& s : inv_list_csv::split(rs_thresholds))
    {
        std::size_t pos = 0; const double t = std::stod(s, &pos);
        if (pos != s.size() || !std::isfinite(t) || t < 0 || t > 1)
            throw std::runtime_error("RS thresholds must be comma-separated numbers in [0,1]");
        thresholds.push_back(t);
    }
    std::sort(thresholds.begin(), thresholds.end());
    thresholds.erase(std::unique(thresholds.begin(), thresholds.end()), thresholds.end());
    const std::vector<std::string> suffixes = {".vectors.csv", ".profiles.csv", ".pairs.csv", ".timings.csv", ".histogram.csv", ".report.txt"};
    for (const auto& s : suffixes)
    { std::ifstream exists(rs_out+s); if (exists) throw std::runtime_error("Output already exists: " + rs_out+s); }
    std::ofstream report(rs_out+".report.txt"); report.exceptions(std::ios::failbit|std::ios::badbit);
    report << "Status: incomplete\nSet=" << rs_set << "\nSource=" << bv_in_file
           << "\nIndex=" << pair_index_in << "\nSeed=" << rs_seed
           << "\nRequested size=" << rs_size << "\nMinimum bytes=" << rs_min_bytes
           << "\nByte coverage=" << rs_byte_coverage << "\nThreads=" << pair_threads
           << "\nSerializer level=" << c_level << "\nMinimum gain=" << pair_min_gain
           << "\nMinimum percent=" << pair_min_pct << "\nThresholds=" << rs_thresholds
           << "\nRepeats=" << rs_repeats << '\n' << std::flush;
    const auto overall_start = pair_clock::now();
    bool complete = false; bm::id64_t indexed_bytes = 0;
    auto entries = pair_read_shard_index(complete, indexed_bytes);
    bm::id64_t collection_bytes = 0;
    for (const auto& e : entries) collection_bytes += e.bytes;
    for (unsigned pct : {1u, 5u, 10u, 25u})
    {
        const std::size_t top = (entries.size()*pct+99)/100;
        bm::id64_t bytes = 0;
        for (std::size_t k = 0; k < top; ++k) bytes += entries[entries.size()-1-k].bytes;
        report << "Largest " << pct << "% vectors: count=" << top << "; bytes=" << bytes
               << "; byte share pct=" << 100.0*double(bytes)/double(collection_bytes) << '\n';
    }
    // Keep a descending size prefix contributing the requested fraction of bytes.
    std::size_t cutoff = entries.size(); bm::id64_t coverage_bytes = 0;
    while (cutoff && double(coverage_bytes) < double(collection_bytes)*rs_byte_coverage/100)
        coverage_bytes += entries[--cutoff].bytes;
    std::vector<std::size_t> eligible;
    std::map<unsigned, std::vector<std::size_t> > buckets;
    std::map<unsigned, bm::id64_t> bucket_weight, all_bytes, all_count;
    bm::id64_t eligible_bytes = 0;
    for (std::size_t i = 0; i < entries.size(); ++i)
    {
        const auto& e = entries[i]; const unsigned b = rs_bucket(e.bytes);
        all_bytes[b] += e.bytes; ++all_count[b];
        if (i < cutoff || e.bytes < rs_min_bytes || e.super_sparse) continue;
        eligible.push_back(i); buckets[b].push_back(i); bucket_weight[b] += e.bytes;
        eligible_bytes += e.bytes;
    }
    if (eligible.size() < 2) throw std::runtime_error("Fewer than two eligible vectors; relax RS sampling options");
    std::mt19937_64 rng(rs_seed);
    std::vector<std::size_t> selected;
    if (rs_set == "shard")
    {
        std::map<std::size_t, bm::id64_t> weights;
        for (auto i : eligible) weights[entries[i].shard] += entries[i].bytes;
        std::size_t chosen = std::size_t(rs_shard);
        if (rs_shard == (std::numeric_limits<bm::id64_t>::max)())
        {
            auto best = weights.begin();
            for (auto it = weights.begin(); it != weights.end(); ++it)
                if (it->second > best->second) best = it;
            chosen = best->first;
        }
        for (auto i : eligible) if (entries[i].shard == chosen) selected.push_back(i);
        report << "Selected shard=" << chosen << '\n';
        if (selected.size() > rs_size)
        { std::shuffle(selected.begin(), selected.end(), rng); selected.resize(std::size_t(rs_size)); }
    }
    else
    {
        // One sample per eligible log2 size bucket, then draw buckets proportional
        // to original byte contribution; draw members without replacement.
        for (auto& b : buckets) std::shuffle(b.second.begin(), b.second.end(), rng);
        for (auto& b : buckets)
            if (selected.size() < rs_size)
            { selected.push_back(b.second.back()); b.second.pop_back(); }
        while (selected.size() < (std::min)(std::size_t(rs_size), eligible.size()))
        {
            std::vector<unsigned> keys; std::vector<double> weights;
            for (const auto& b : buckets) if (!b.second.empty())
            { keys.push_back(b.first); weights.push_back(double(bucket_weight[b.first])); }
            std::discrete_distribution<std::size_t> draw(weights.begin(), weights.end());
            auto& bucket = buckets[keys[draw(rng)]];
            selected.push_back(bucket.back()); bucket.pop_back();
        }
    }
    if (selected.size() < 2) throw std::runtime_error("Selected shard has fewer than two eligible vectors");
    std::sort(selected.begin(), selected.end(), [&](std::size_t a, std::size_t b){return entries[a].id < entries[b].id;});
    const std::size_t n = selected.size(); const auto pair_count = pair_choose2(n);
    std::map<unsigned, bm::id64_t> sample_count, sample_bytes;
    for (auto i : selected) { ++sample_count[rs_bucket(entries[i].bytes)]; sample_bytes[rs_bucket(entries[i].bytes)] += entries[i].bytes; }
    std::ofstream histogram(rs_out+".histogram.csv"); histogram.exceptions(std::ios::failbit|std::ios::badbit);
    histogram << "log2_size,collection_vectors,collection_bytes,eligible_bytes,sample_vectors,sample_bytes\n";
    for (const auto& b : all_bytes)
        histogram << b.first << ',' << all_count[b.first] << ',' << b.second << ',' << bucket_weight[b.first]
                  << ',' << sample_count[b.first] << ',' << sample_bytes[b.first] << '\n';
    histogram.close();
    report << "Indexed collection complete=" << complete << "\nIndexed vectors=" << entries.size()
           << "\nCollection standalone bytes=" << collection_bytes << "\nEligible vectors=" << eligible.size()
           << "\nEligible bytes=" << eligible_bytes << "\nExcluded target byte ceiling=" << collection_bytes-eligible_bytes
           << "\nSelected vectors=" << n << "\nPairs=" << pair_count
           << "\nNo pair-size ratio filter. No high-density exclusion. Labels are structural, not measured bytes.\n" << std::flush;
    cout << "RS set=" << rs_set << "; vectors=" << n << "; pairs=" << pair_count << endl;
    std::ifstream input(bv_in_file, std::ios::binary);
    if (!input) throw std::runtime_error("Cannot open RS source");
    input.seekg(0, std::ios::end); const auto file_size = input.tellg(); input.seekg(0);
    if (file_size < std::streampos(0) || bm::id64_t(file_size) < indexed_bytes ||
        (complete && bm::id64_t(file_size) != indexed_bytes)) throw std::runtime_error("RS input/index size mismatch");
    bm::streams_deserializer<pair_bv> reader; bm::serializer<pair_bv> serializer;
    serializer.set_compression_level(c_level); bm::serializer<pair_bv>::buffer buffer;
    std::vector<std::unique_ptr<pair_item> > items;
    auto stage = pair_clock::now(); auto last = stage;
    for (auto i : selected)
    {
        items.push_back(pair_load_indexed_item(input, entries[i], reader, serializer, buffer));
        if (!is_silent && rs_seconds(last) >= 10)
        { cout << "Loaded RS vectors " << items.size() << '/' << n << endl; last = pair_clock::now(); }
    }
    const double load_s = rs_seconds(stage);
    std::vector<rs_profile> profiles(n); stage = pair_clock::now();
    rs_parallel(0, n, [&](std::size_t i){profiles[i] = rs_extract(items[i]->bv);}, "Build RS profiles");
    const double profile_s = rs_seconds(stage);
    bm::id64_t profile_bytes = 0, bv_memory = 0;
    std::ofstream vectors(rs_out+".vectors.csv"); vectors.exceptions(std::ios::failbit|std::ios::badbit);
    std::ofstream profile_file(rs_out+".profiles.csv"); profile_file.exceptions(std::ios::failbit|std::ios::badbit);
    profile_file << "id,block,first_count,second_count,third_count\n";
    vectors << "id,shard,serialized_bytes,popcount,structural_cost,memory_bytes,stored_bit_blocks,stored_gap_blocks,runs,block_runs,zero_blocks,full_blocks,predicted_gap_blocks,gap_words,predicted_bit_blocks,zero_spans,full_spans,rs_regions,rs_bytes\n";
    for (std::size_t i = 0; i < n; ++i)
    {
        const auto& e = entries[selected[i]]; pair_bv::statistics st; items[i]->bv.calc_stat(&st);
        const auto c = bm::calc_complexity(items[i]->bv);
        const auto bytes = profiles[i].regions.capacity()*sizeof(std::uint16_t); profile_bytes += bytes; bv_memory += st.memory_used;
        vectors << e.id << ',' << e.shard << ',' << e.bytes << ',' << e.count << ',' << e.cost << ',' << st.memory_used
                << ',' << st.bit_blocks << ',' << st.gap_blocks << ',' << c.runs << ',' << c.block_runs
                << ',' << c.zero_blocks << ',' << c.full_blocks << ',' << c.gap_blocks << ',' << c.gap_words
                << ',' << c.bit_blocks << ',' << c.zero_spans << ',' << c.full_spans << ',' << profiles[i].regions.size() << ',' << bytes << '\n';
        for (std::size_t r = 0; r < profiles[i].regions.size(); r += 3)
            profile_file << e.id << ',' << r/3 << ',' << profiles[i].regions[r] << ','
                         << profiles[i].regions[r+1] << ',' << profiles[i].regions[r+2] << '\n';
    }
    vectors.close(); profile_file.close();
    // Independent timing passes have identical traversal and aggregation and NO pair CSV I/O.
    std::ofstream timings(rs_out+".timings.csv"); timings.exceptions(std::ios::failbit|std::ios::badbit);
    timings << "repeat,order,mode,threshold,pairs,evaluated,rejected,qualified,seconds,load_seconds,profile_seconds,profile_bytes,bvector_bytes,checksum_cost,checksum_count\n" << std::setprecision(17);
    std::vector<rs_totals> observed(thresholds.size()+1);
    for (bm::id64_t repeat = 0; repeat < rs_repeats; ++repeat)
    {
        std::vector<std::size_t> modes(thresholds.size()+1);
        for (std::size_t k = 0; k < modes.size(); ++k) modes[k] = k;
        if (repeat % 2) std::reverse(modes.begin(), modes.end());
        for (std::size_t order = 0; order < modes.size(); ++order)
        {
            const auto mode = modes[order]; const double threshold = mode ? thresholds[mode-1] : -1;
            std::vector<rs_totals> rows(n);
            cout << "Timing repeat " << repeat+1 << '/' << rs_repeats << ": "
                 << (mode ? "RS threshold=" : "baseline ") << threshold << endl;
            stage = pair_clock::now();
            rs_parallel(0, n, [&](std::size_t i)
            {
                auto& r = rows[i];
                for (std::size_t j = i+1; j < n; ++j)
                {
                    if (mode && !rs_accept(profiles[i], profiles[j], threshold)) continue;
                    const auto c = bm::calc_complexity_xor(items[i]->bv, items[j]->bv);
                    const auto cost = pair_cost(c); ++r.evaluated;
                    r.qualified += rs_qualifies(cost, (std::max)(items[i]->cost, items[j]->cost));
                    r.checksum_cost += cost; r.checksum_count += c.count;
                }
            }, "Timing pairs");
            const double seconds = rs_seconds(stage);
            rs_totals totals; for (const auto& row : rows) totals.add(row);
            if (repeat && (totals.evaluated != observed[mode].evaluated || totals.qualified != observed[mode].qualified ||
                totals.checksum_cost != observed[mode].checksum_cost || totals.checksum_count != observed[mode].checksum_count))
                throw std::runtime_error("RS repeated timing results differ");
            observed[mode] = totals;
            timings << repeat+1 << ',' << order+1 << ',' << (mode ? "rs" : "baseline") << ',' << threshold
                    << ',' << pair_count << ',' << totals.evaluated << ',' << pair_count-totals.evaluated << ',' << totals.qualified
                    << ',' << seconds << ',' << load_s << ',' << profile_s << ',' << profile_bytes << ',' << bv_memory
                    << ',' << totals.checksum_cost << ',' << totals.checksum_count << '\n' << std::flush;
            cout << "Timing seconds=" << seconds << "; evaluated=" << totals.evaluated << "; qualifying=" << totals.qualified << endl;
        }
    }
    timings.close();
    std::ofstream pairs(rs_out+".pairs.csv"); pairs.exceptions(std::ios::failbit|std::ios::badbit);
    pairs << "a_id,b_id,target_id,a_bytes,b_bytes,a_count,b_count,a_cost,b_cost,xor_cost,xor_count,qualified,block_l1,region_l1,rs_distance,shape_distance\n" << std::setprecision(17);
    std::vector<rs_totals> expected(thresholds.size()+1);
    stage = pair_clock::now(); last = stage;
    for (std::size_t begin = 0; begin < n; begin += 16)
    {
        const auto end = (std::min)(n, begin+16);
        std::vector<std::vector<rs_pair_result> > rows(end-begin);
        rs_parallel(begin, end, [&](std::size_t i)
        {
            auto& row = rows[i-begin]; row.reserve(n-i-1);
            for (std::size_t j = i+1; j < n; ++j)
            {
                const auto c = bm::calc_complexity_xor(items[i]->bv, items[j]->bv);
                rs_pair_result r; r.cost = pair_cost(c); r.count = c.count;
                rs_distance(&r.d, profiles[i], profiles[j]);
                if (r.d.block_l1 > r.d.region_l1 || r.d.region_l1 > r.count)
                    throw std::runtime_error("RS Hamming lower-bound violation");
                row.push_back(r);
            }
        });
        for (std::size_t i = begin; i < end; ++i)
            for (std::size_t j = i+1; j < n; ++j)
            {
                const auto& r = rows[i-begin][j-i-1];
                const auto& a = entries[selected[i]]; const auto& b = entries[selected[j]];
                const auto baseline = (std::max)(a.cost,b.cost);
                const bool qualified = rs_qualifies(r.cost, baseline);
                const auto denominator = a.count+b.count;
                const double distance = denominator ? double(r.d.region_l1)/double(denominator) : 0;
                const auto target = a.cost > b.cost ? a.id : b.id; // IDs ascend; equal-cost lead is smaller ID
                pairs << a.id << ',' << b.id << ',' << target << ',' << a.bytes << ',' << b.bytes
                      << ',' << a.count << ',' << b.count << ',' << a.cost << ',' << b.cost << ',' << r.cost
                      << ',' << r.count << ',' << unsigned(qualified) << ',' << r.d.block_l1 << ',' << r.d.region_l1
                      << ',' << distance << ',' << r.d.shape << '\n';
                for (std::size_t k = 0; k < expected.size(); ++k)
                    if (!k || double(r.d.region_l1) <= thresholds[k-1]*double(denominator))
                    {
                        auto& e = expected[k]; ++e.evaluated; e.qualified += qualified;
                        e.checksum_cost += r.cost; e.checksum_count += r.count;
                    }
            }
        if (!is_silent && rs_seconds(last) >= 10)
        { cout << "Detailed pair CSV: rows of vectors " << end << '/' << n << endl; last = pair_clock::now(); }
    }
    pairs.close();
    for (std::size_t k = 0; k < expected.size(); ++k)
        if (expected[k].evaluated != observed[k].evaluated || expected[k].qualified != observed[k].qualified ||
            expected[k].checksum_cost != observed[k].checksum_cost || expected[k].checksum_count != observed[k].checksum_count)
            throw std::runtime_error("RS screening/timing disagrees with full pair labels");
    report << "Load seconds=" << load_s << "\nRS construction and extraction seconds=" << profile_s
           << "\nRS cached population bytes=" << profile_bytes << "\nBvector bytes=" << bv_memory
           << "\nDetailed pair collection seconds=" << rs_seconds(stage)
           << "\nTotal experiment seconds=" << rs_seconds(overall_start)
           << "\nTiming counters and checksums agree with full pair labels.\nStatus: complete\n";
    report.close(); cout << "RS benchmark complete: " << rs_out << endl;
}
} // namespace
#endif
