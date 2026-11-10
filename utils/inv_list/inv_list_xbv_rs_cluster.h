/** @file inv_list_xbv_rs_cluster.h
    @brief Cross-bvector (xbv) research support for inv_list.
    @note Development-only; included by inv_list.cpp after shared options and I/O.
*/
/* Copyright(c) 2026 Anatoliy Kuznetsov
   Licensed under the Apache License, Version 2.0.
   http://www.apache.org/licenses/LICENSE-2.0 */
#ifndef INV_LIST_XBV_RS_CLUSTER_H
#define INV_LIST_XBV_RS_CLUSTER_H

#include <numeric>

namespace
{

struct rs_cluster_profile
{
    rs_profile fine;
    std::vector<unsigned> blocks8;
};

void rs_cluster_finish_feature(rs_cluster_profile& p, bool full,
                               std::size_t& dimensions, bm::id64_t& profile_bytes)
{
    if (full)
    {
        dimensions = (std::max)(dimensions, p.fine.regions.size());
        profile_bytes += p.fine.regions.capacity()*sizeof(std::uint16_t);
        return;
    }
    const auto& regions = p.fine.regions;
    p.blocks8.reserve((regions.size()+23)/24);
    for (std::size_t r = 0; r < regions.size(); r += 24)
    {
        unsigned total = 0;
        for (std::size_t t = r; t < (std::min)(regions.size(),r+24); ++t)
            total += regions[t];
        p.blocks8.push_back(total);
    }
    dimensions = (std::max)(dimensions, p.blocks8.size());
    profile_bytes += p.blocks8.capacity()*sizeof(unsigned);
    p.fine.regions.clear(); p.fine.regions.shrink_to_fit();
}

void rs_cluster_read_profiles(const std::vector<pair_index_entry>& index,
                              std::vector<pair_index_entry>& entries,
                              std::vector<rs_cluster_profile>& profiles)
{
    std::ifstream status(rs_cluster_profiles_in + ".report.txt");
    std::string line; bool complete = false;
    while (std::getline(status, line)) if (line == "Status: complete") complete = true;
    if (!complete) throw std::runtime_error("RS benchmark profile report is not complete");
    std::vector<const pair_index_entry*> by_id(index.size(), nullptr);
    for (const auto& e : index) by_id[e.id] = &e;
    std::ifstream vectors(rs_cluster_profiles_in + ".vectors.csv");
    if (!vectors || !std::getline(vectors, line))
        throw std::runtime_error("Cannot open RS benchmark vectors CSV");
    const auto vh = inv_list_csv::split(line);
    const auto v_id = inv_list_csv::column(vh, "id");
    const auto v_bytes = inv_list_csv::column(vh, "serialized_bytes");
    const auto v_count = inv_list_csv::column(vh, "popcount");
    const auto v_memory = inv_list_csv::column(vh, "memory_bytes");
    const auto v_regions = inv_list_csv::column(vh, "rs_regions");
    const auto v_required = (std::max)({v_id,v_bytes,v_count,v_memory,v_regions});
    std::vector<std::size_t> id_position(index.size(), index.size());
    std::vector<std::size_t> expected_regions;
    for (std::size_t line_no=2; std::getline(vectors,line); ++line_no)
    {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty()) continue;
        const auto f = inv_list_csv::split(line);
        if (f.size() <= v_required) throw std::runtime_error("Short RS benchmark vector row");
        const auto id = std::size_t(inv_list_csv::parse_u64(f[v_id],line_no,"id"));
        if (id >= by_id.size() || !by_id[id] || id_position[id] != index.size())
            throw std::runtime_error("Invalid RS benchmark vector ID");
        const auto& e = *by_id[id];
        if (e.bytes != inv_list_csv::parse_u64(f[v_bytes],line_no,"serialized_bytes") ||
            e.count != inv_list_csv::parse_u64(f[v_count],line_no,"popcount") ||
            e.memory != inv_list_csv::parse_u64(f[v_memory],line_no,"memory_bytes"))
            throw std::runtime_error("RS benchmark vectors do not match pair index");
        id_position[id] = entries.size(); entries.push_back(e);
        expected_regions.push_back(std::size_t(inv_list_csv::parse_u64(f[v_regions],line_no,"rs_regions")));
    }
    if (entries.size() < 2) throw std::runtime_error("RS benchmark has fewer than two vectors");
    profiles.resize(entries.size());
    std::ifstream in(rs_cluster_profiles_in + ".profiles.csv");
    if (!in || !std::getline(in,line)) throw std::runtime_error("Cannot open RS benchmark profiles CSV");
    const auto h = inv_list_csv::split(line);
    const auto c_id = inv_list_csv::column(h,"id"), c_block = inv_list_csv::column(h,"block");
    const auto c_first = inv_list_csv::column(h,"first_count");
    const auto c_second = inv_list_csv::column(h,"second_count");
    const auto c_third = inv_list_csv::column(h,"third_count");
    const auto required = (std::max)({c_id,c_block,c_first,c_second,c_third});
    for (std::size_t line_no=2; std::getline(in,line); ++line_no)
    {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty()) continue;
        const auto f = inv_list_csv::split(line);
        if (f.size() <= required) throw std::runtime_error("Short RS benchmark profile row");
        const auto id = std::size_t(inv_list_csv::parse_u64(f[c_id],line_no,"id"));
        if (id >= id_position.size() || id_position[id] == index.size())
            throw std::runtime_error("Unknown RS benchmark profile ID");
        auto& p = profiles[id_position[id]];
        const auto block = inv_list_csv::parse_u64(f[c_block],line_no,"block");
        if (block != p.fine.regions.size()/3)
            throw std::runtime_error("RS benchmark profile blocks are incomplete or out of order");
        for (const auto column : {c_first,c_second,c_third})
        {
            const auto count = inv_list_csv::parse_u64(f[column],line_no,"region_count");
            if (count > (std::numeric_limits<std::uint16_t>::max)())
                throw std::runtime_error("RS benchmark region count exceeds 16 bits");
            p.fine.regions.push_back(std::uint16_t(count));
            p.fine.count += count;
        }
    }
    for (std::size_t i=0;i<entries.size();++i)
        if (profiles[i].fine.regions.size() != expected_regions[i] ||
            profiles[i].fine.count != entries[i].count)
            throw std::runtime_error("RS benchmark profile is incomplete");
}

std::vector<unsigned> rs_cluster_features(const rs_cluster_profile& p, bool full,
                                          std::size_t dimensions)
{
    std::vector<unsigned> result(dimensions);
    if (full)
        for (std::size_t j = 0; j < p.fine.regions.size(); ++j)
            result[j] = p.fine.regions[j];
    else
        std::copy(p.blocks8.begin(), p.blocks8.end(), result.begin());
    return result;
}

double rs_cluster_distance(const rs_cluster_profile& p,
                           const std::vector<unsigned>& center,
                           bm::id64_t center_count, bool full)
{
    bm::id64_t total = 0;
    if (full)
    {
        const auto& row = p.fine.regions;
        std::size_t j = 0;
        for (; j < row.size(); ++j) total += rs_abs(row[j], center[j]);
        for (; j < center.size(); ++j) total += center[j];
    }
    else
    {
        const auto& row = p.blocks8;
        std::size_t j = 0;
        for (; j < row.size(); ++j) total += rs_abs(row[j], center[j]);
        for (; j < center.size(); ++j) total += center[j];
    }
    const auto denom = p.fine.count + center_count;
    return denom ? double(total) / double(denom) : 0;
}

std::vector<std::vector<unsigned> > rs_cluster_seeds(
    const std::vector<rs_cluster_profile>& profiles, bool full,
    std::size_t groups, std::size_t dimensions)
{
    std::mt19937_64 rng(rs_cluster_seed);
    std::uniform_int_distribution<std::size_t> pick(0, profiles.size()-1);
    std::vector<std::vector<unsigned> > centers;
    std::vector<bm::id64_t> counts;
    std::vector<bool> chosen(profiles.size(), false);
    auto add_seed = [&](std::size_t i)
    {
        chosen[i] = true;
        centers.push_back(rs_cluster_features(profiles[i], full, dimensions));
        counts.push_back(profiles[i].fine.count);
    };
    add_seed(pick(rng));
    while (centers.size() < groups)
    {
        std::size_t best = profiles.size(); double farthest = -1;
        for (unsigned trial = 0; trial < 64; ++trial)
        {
            const std::size_t i = pick(rng);
            if (chosen[i]) continue;
            double nearest = (std::numeric_limits<double>::max)();
            for (std::size_t g = 0; g < centers.size(); ++g)
                nearest = (std::min)(nearest, rs_cluster_distance(profiles[i], centers[g], counts[g], full));
            if (nearest > farthest || (nearest == farthest && i < best))
            { best = i; farthest = nearest; }
        }
        if (best == profiles.size())
            for (std::size_t i = 0; i < profiles.size(); ++i)
                if (!chosen[i]) { best = i; break; }
        add_seed(best);
    }
    return centers;
}

bool rs_cluster_assign(const std::vector<rs_cluster_profile>& profiles,
                       const std::vector<pair_index_entry>& entries,
                       const std::vector<std::vector<unsigned> >& centers,
                       bool full, const std::vector<std::size_t>& capacities,
                       std::vector<std::size_t>& labels,
                       std::vector<bm::id64_t>& group_memory, double& objective)
{
    const std::size_t n = profiles.size(), groups = centers.size();
    std::vector<bm::id64_t> center_counts(groups);
    for (std::size_t g = 0; g < groups; ++g)
        for (const auto v : centers[g]) center_counts[g] += v;
    std::vector<std::size_t> order(n), counts(groups);
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b)
    { return entries[a].memory != entries[b].memory ? entries[a].memory > entries[b].memory : entries[a].id < entries[b].id; });
    labels.assign(n, groups); group_memory.assign(groups, 0); objective = 0;
    for (const auto i : order)
    {
        std::size_t best = groups; double best_distance = (std::numeric_limits<double>::max)();
        for (std::size_t g = 0; g < groups; ++g)
        {
            if (counts[g] >= capacities[g] || entries[i].memory > rs_cluster_memory - group_memory[g])
                continue;
            const double d = rs_cluster_distance(profiles[i], centers[g], center_counts[g], full);
            if (d < best_distance || (d == best_distance && g < best))
            { best = g; best_distance = d; }
        }
        if (best == groups) return false;
        labels[i] = best; ++counts[best]; group_memory[best] += entries[i].memory;
        objective += best_distance;
    }
    return true;
}

void build_rs_cluster_plan()
{
    if (bv_in_file.empty() || pair_index_in.empty() || rs_cluster_out.empty())
        throw std::runtime_error("RS cluster planning requires -bvin, -pair-index-in, and -rs-cluster-out");
    if (rs_cluster_feature != "regions3" && rs_cluster_feature != "blocks8")
        throw std::runtime_error("-rs-cluster-feature must be regions3 or blocks8");
    if (!rs_cluster_vectors || !rs_cluster_memory)
        throw std::runtime_error("RS cluster vector and memory limits must be positive");
    const std::string report_name = rs_cluster_out + ".report.txt";
    for (const auto& name : {rs_cluster_out, report_name})
    { std::ifstream f(name, std::ios::binary); if (f) throw std::runtime_error("Output already exists: " + name); }
    const auto start = pair_clock::now(); auto last = start;
    std::ofstream report(report_name); report.exceptions(std::ios::failbit|std::ios::badbit);
    report << "Status: incomplete\nFeature=" << rs_cluster_feature
           << "\nSource=" << bv_in_file << "\nIndex=" << pair_index_in
           << "\nSeed=" << rs_cluster_seed << "\nVector cap=" << rs_cluster_vectors
           << "\nDecoded memory cap=" << rs_cluster_memory
           << "\nMinimum bytes=" << rs_min_bytes << "\nByte coverage=" << rs_byte_coverage
           << "\nSerializer level=" << c_level << '\n' << std::flush;

    bool source_complete = false; bm::id64_t indexed_bytes = 0;
    auto index = pair_read_shard_index(source_complete, indexed_bytes);
    std::ifstream input(bv_in_file, std::ios::binary);
    input.seekg(0, std::ios::end);
    const auto file_size = input.tellg(); input.seekg(0);
    if (!input || file_size < std::streampos(0) || bm::id64_t(file_size) < indexed_bytes ||
        (source_complete && bm::id64_t(file_size) != indexed_bytes))
        throw std::runtime_error("RS cluster index input size does not match collection");
    bm::id64_t total_bytes = 0;
    for (const auto& e : index) total_bytes += e.bytes;
    std::size_t cutoff = index.size(); bm::id64_t covered = 0;
    while (cutoff && double(covered) < double(total_bytes)*rs_byte_coverage/100)
        covered += index[--cutoff].bytes;
    std::vector<pair_index_entry> entries;
    bm::id64_t selected_memory = 0, selected_bytes = 0;
    std::vector<rs_cluster_profile> profiles;
    if (!rs_cluster_profiles_in.empty())
        rs_cluster_read_profiles(index,entries,profiles);
    else
        for (std::size_t i = cutoff; i < index.size(); ++i)
        {
            const auto& e = index[i];
            if (e.bytes < rs_min_bytes || e.super_sparse) continue;
            entries.push_back(e);
        }
    if (entries.size() < 2) throw std::runtime_error("Fewer than two eligible RS cluster vectors");
    if (rs_cluster_profiles_in.empty())
        std::sort(entries.begin(), entries.end(), [](const auto& a, const auto& b){return a.id < b.id;});
    for (const auto& e : entries)
    {
        if (e.memory > rs_cluster_memory)
            throw std::runtime_error("A vector exceeds the RS cluster decoded memory cap");
        selected_memory += e.memory; selected_bytes += e.bytes;
    }
    report << "Indexed vectors=" << index.size() << "\nEligible vectors=" << entries.size()
           << "\nBenchmark profiles=" << (rs_cluster_profiles_in.empty() ? "none" : rs_cluster_profiles_in)
           << "\nEligible standalone bytes=" << selected_bytes
           << "\nEligible decoded memory estimate=" << selected_memory << '\n' << std::flush;

    if (rs_cluster_profiles_in.empty()) profiles.resize(entries.size());
    bm::streams_deserializer<pair_bv> reader;
    bm::serializer<pair_bv> serializer; serializer.set_compression_level(c_level);
    bm::serializer<pair_bv>::buffer buffer;
    std::size_t dimensions = 0; bm::id64_t profile_bytes = 0;
    for (std::size_t i = 0; i < entries.size(); ++i)
    {
        auto& p = profiles[i];
        if (rs_cluster_profiles_in.empty())
        {
            auto item = pair_load_indexed_item(input, entries[i], reader, serializer, buffer);
            p.fine = rs_extract(item->bv);
        }
        if (p.fine.count != entries[i].count)
            throw std::runtime_error("RS cluster population does not match pair index");
        rs_cluster_finish_feature(p,rs_cluster_feature == "regions3",dimensions,profile_bytes);
        if (!is_silent && pair_clock::now()-last >= std::chrono::seconds(10))
        { cout << "RS profiles " << i+1 << '/' << entries.size() << endl; last = pair_clock::now(); }
    }
    const double extraction_seconds = rs_seconds(start);
    const bool full = rs_cluster_feature == "regions3";
    const std::size_t n = entries.size();
    const auto ceil_div = [](bm::id64_t a, bm::id64_t b){return 1+(a-1)/b;};
    const bm::id64_t memory_target = rs_cluster_memory - rs_cluster_memory/4;
    const std::size_t initial_groups = std::size_t((std::min)(bm::id64_t(n),
        (std::max)(ceil_div(n,rs_cluster_vectors),
                   selected_memory ? ceil_div(selected_memory,memory_target) : bm::id64_t(1))));
    std::size_t groups = initial_groups;
    std::vector<std::size_t> capacities;
    std::vector<std::size_t> labels;
    std::vector<bm::id64_t> group_memory;
    double objective = 0;
    for (; groups <= n; ++groups)
    {
        capacities.assign(groups, n/groups);
        for (std::size_t g = 0; g < n%groups; ++g) ++capacities[g];
        auto centers = rs_cluster_seeds(profiles, full, groups, dimensions);
        if (!rs_cluster_assign(profiles, entries, centers, full, capacities,
                               labels, group_memory, objective))
            continue;
        std::vector<std::vector<bm::id64_t> > sums(groups, std::vector<bm::id64_t>(dimensions));
        for (std::size_t i = 0; i < n; ++i)
        {
            auto& dst = sums[labels[i]];
            if (full)
                for (std::size_t j = 0; j < profiles[i].fine.regions.size(); ++j)
                    dst[j] += profiles[i].fine.regions[j];
            else
                for (std::size_t j = 0; j < profiles[i].blocks8.size(); ++j)
                    dst[j] += profiles[i].blocks8[j];
        }
        for (std::size_t g = 0; g < groups; ++g)
            for (std::size_t j = 0; j < dimensions; ++j)
                centers[g][j] = unsigned((sums[g][j]+capacities[g]/2)/capacities[g]);
        if (rs_cluster_assign(profiles, entries, centers, full, capacities,
                              labels, group_memory, objective))
            break;
    }
    if (groups > n)
        throw std::runtime_error("RS cluster assignment cannot meet decoded memory cap");

    std::vector<std::vector<std::size_t> > members(groups);
    for (std::size_t i = 0; i < n; ++i) members[labels[i]].push_back(i);
    std::ofstream manifest(rs_cluster_out); manifest.exceptions(std::ios::failbit|std::ios::badbit);
    manifest << "cluster_id,vector_id,source_offset,source_bytes,indexed_input_bytes\n";
    bm::id64_t planned_pairs = 0, peak_memory = 0;
    for (std::size_t g = 0; g < groups; ++g)
    {
        planned_pairs += pair_choose2(members[g].size());
        peak_memory = (std::max)(peak_memory,group_memory[g]);
        for (auto i : members[g])
        {
            const auto& e = entries[i];
            manifest << g << ',' << e.id << ',' << e.offset << ',' << e.source_bytes
                     << ',' << indexed_bytes << '\n';
        }
    }
    manifest.close();
    report << "Initial clusters=" << initial_groups << "\nClusters=" << groups
           << "\nPlanned internal pairs=" << planned_pairs
           << "\nAssignment=sampled farthest seeds, memory-first nearest center, one centroid update"
           << "\nAll eligible pairs=" << pair_choose2(n)
           << "\nOmitted cross-cluster pairs=" << pair_choose2(n)-planned_pairs
           << "\nMaximum cluster decoded memory=" << peak_memory
           << "\nCached feature bytes=" << profile_bytes
           << "\nAssignment objective=" << objective
           << "\nProfile extraction seconds=" << extraction_seconds
           << "\nTotal planning seconds=" << rs_seconds(start) << '\n';
    for (std::size_t g = 0; g < groups; ++g)
        report << "Cluster " << g << ": vectors=" << members[g].size()
               << "; decoded memory MiB=" << double(group_memory[g])/(1024*1024)
               << "; pairs=" << pair_choose2(members[g].size()) << '\n';
    report << "Status: complete\n";
    cout << "RS cluster plan complete. Eligible vectors=" << n << "; clusters=" << groups
         << "; internal pairs=" << planned_pairs << "; peak estimated decoded MiB="
         << double(peak_memory)/(1024*1024) << "; seconds=" << rs_seconds(start) << endl;
}

} // namespace
#endif
