// Read the raw Gov2 u32 collection and measure optimized bvector storage.
#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "bm.h"

namespace {
using count_t = std::uint64_t;
using scan_clock = std::chrono::steady_clock;
constexpr unsigned levels = bm::gap_levels;
constexpr count_t block_bits = bm::bits_in_block;

struct options {
    std::string input, out;
    count_t max_vectors = 0, progress = 1000;
};

count_t number(const char* value) {
    const std::string s(value);
    if (s.empty() || s.find_first_not_of("0123456789") != std::string::npos)
        throw std::runtime_error("Expected a nonnegative integer: " + s);
    return std::stoull(s);
}

options parse(int argc, char** argv) {
    options o;
    for (int i = 1; i < argc; ++i) {
        const std::string key(argv[i]);
        if (++i == argc) throw std::runtime_error("Missing value for " + key);
        if (key == "--input") o.input = argv[i];
        else if (key == "--out") o.out = argv[i];
        else if (key == "--max-vectors") o.max_vectors = number(argv[i]);
        else if (key == "--progress") o.progress = number(argv[i]);
        else throw std::runtime_error("Unknown option: " + key);
    }
    if (o.input.empty() || o.out.empty())
        throw std::runtime_error("Usage: gov2_block_stats --input gov2.sorted --out DIR [--max-vectors N] [--progress N]");
    return o;
}

struct totals {
    count_t vectors = 0, ids = 0, input_bytes = 0;
    count_t logical = 0, zero = 0, full = 0, bit = 0, gap = 0;
    count_t gap_words = 0, gap_level[levels] = {}, gap_level_words[levels] = {};
    count_t vectors_bit = 0, vectors_gap = 0, vectors_mixed = 0;
    count_t vectors_high_gap = 0, level_table_changes = 0;
    bm::gap_word_t first_levels[levels] = {};
    std::vector<count_t> bit_per_vector, gap_per_vector, high_gap_per_vector;
    std::vector<std::array<count_t, levels>> gap_lengths =
        std::vector<std::array<count_t, levels>>(65536);
};

double pct(count_t part, count_t whole) {
    return whole ? 100.0 * double(part) / double(whole) : 0.0;
}

count_t percentile(std::vector<count_t>& values, unsigned p) {
    if (values.empty()) return 0;
    const size_t k = size_t((values.size() - 1) * p / 100);
    std::nth_element(values.begin(), values.begin() + std::vector<count_t>::difference_type(k), values.end());
    return values[k];
}

void inspect(const bm::bvector<>& bv, const bm::bvector<>::statistics& st,
             count_t last_id, count_t vector_index, count_t ids,
             totals& t, std::ofstream& per_vector) {
    const count_t logical = last_id / block_bits + 1;
    count_t bit = 0, gap = 0, full = 0, gap_words = 0, max_gap = 0;
    count_t per_level[levels] = {}, per_level_words[levels] = {};
    const auto& manager = bv.get_blocks_manager();
    auto root = manager.top_blocks_root();
    for (unsigned i = 0; i < manager.top_block_size(); ++i) {
        const count_t base = count_t(i) * bm::set_sub_array_size;
        if (base >= logical) break;
        const unsigned width = unsigned(std::min<count_t>(bm::set_sub_array_size, logical - base));
        const bm::word_t* const* sub = root ? root[i] : nullptr;
        if (reinterpret_cast<const bm::word_t*>(sub) == FULL_BLOCK_FAKE_ADDR) {
            full += width;
            continue;
        }
        if (!sub) continue;
        for (unsigned j = 0; j < width; ++j) {
            const bm::word_t* block = sub[j];
            if (block == FULL_BLOCK_FAKE_ADDR) { ++full; continue; }
            if (!block) continue;
            if (BM_IS_GAP(block)) {
                const auto* gb = BMGAP_PTR(block);
                const unsigned level = bm::gap_level(gb);
                const unsigned len = bm::gap_length(gb);
                if (level >= levels || len >= t.gap_lengths.size())
                    throw std::runtime_error("Invalid GAP level or length");
                ++gap; ++per_level[level]; ++t.gap_lengths[len][level];
                per_level_words[level] += len;
                gap_words += len;
                max_gap = std::max<count_t>(max_gap, len);
            } else ++bit;
        }
    }
    if (bit != st.bit_blocks || gap != st.gap_blocks)
        throw std::runtime_error("Block traversal disagrees with calc_stat at vector " + std::to_string(vector_index));
    for (unsigned level = 0; level < levels; ++level)
        if (per_level[level] != st.gaps_by_level[level])
            throw std::runtime_error("GAP level histogram disagrees with calc_stat");
    if (bit + gap + full > logical)
        throw std::runtime_error("Block counts exceed logical span");
    const count_t zero = logical - bit - gap - full;
    t.logical += logical; t.zero += zero; t.full += full;
    t.bit += bit; t.gap += gap; t.gap_words += gap_words;
    for (unsigned level = 0; level < levels; ++level) {
        t.gap_level[level] += per_level[level];
        t.gap_level_words[level] += per_level_words[level];
    }
    t.vectors_bit += bit != 0;
    t.vectors_gap += gap != 0;
    t.vectors_mixed += bit != 0 && gap != 0;
    t.vectors_high_gap += per_level[2] + per_level[3] != 0;
    t.bit_per_vector.push_back(bit);
    t.gap_per_vector.push_back(gap);
    t.high_gap_per_vector.push_back(per_level[2] + per_level[3]);
    per_vector << vector_index << ',' << ids << ',' << last_id << ','
               << logical << ',' << zero << ',' << full << ',' << bit << ',' << gap;
    for (unsigned level = 0; level < levels; ++level)
        per_vector << ',' << per_level[level];
    per_vector << ',' << gap_words << ',' << max_gap << '\n';
}

void report(const options& o, totals& t, double seconds) {
    std::ofstream lengths(o.out + "/gap_lengths.csv");
    if (!lengths) throw std::runtime_error("Cannot write gap_lengths.csv");
    lengths << "gap_words,level0,level1,level2,level3,total\n";
    for (unsigned len = 0; len < t.gap_lengths.size(); ++len) {
        const auto& by_level = t.gap_lengths[len];
        const count_t sum = by_level[0] + by_level[1] + by_level[2] + by_level[3];
        if (sum)
            lengths << len << ',' << by_level[0] << ',' << by_level[1] << ','
                    << by_level[2] << ',' << by_level[3] << ',' << sum << '\n';
    }
    lengths.close();
    std::ofstream out(o.out + "/report.txt");
    if (!out) throw std::runtime_error("Cannot write report.txt");
    const count_t materialized = t.bit + t.gap;
    const count_t nonzero = materialized + t.full;
    out << std::fixed << std::setprecision(3);
    out << "Input: " << o.input << '\n'
        << "Measurement: raw u32 lists -> bvector::set(BM_SORTED) -> optimize -> calc_stat\n"
        << "Vectors: " << t.vectors << "  Original IDs: " << t.ids
        << "  Input bytes read: " << t.input_bytes << "  Elapsed s: " << seconds << '\n';
    out << "Logical block positions [0..last ID] per vector: " << t.logical << '\n'
        << "Materialized BIT+GAP: " << materialized << "  Nonzero including implicit full: " << nonzero << '\n';
    out << "Type               Blocks       % materialized   % nonzero   % logical\n";
    const std::array<std::pair<const char*, count_t>, 4> kinds{{
        {"BIT", t.bit}, {"GAP", t.gap}, {"implicit full", t.full}, {"implicit zero", t.zero}}};
    for (const auto& kind : kinds) {
        const bool material = kind.first[0] == 'B' || kind.first[0] == 'G';
        const bool occupied = material || std::string(kind.first) == "implicit full";
        out << std::left << std::setw(16) << kind.first << std::right << std::setw(14)
            << kind.second;
        if (material) out << std::setw(18) << pct(kind.second, materialized);
        else out << std::setw(18) << "n/a";
        if (occupied) out << std::setw(12) << pct(kind.second, nonzero);
        else out << std::setw(12) << "n/a";
        out << std::setw(12) << pct(kind.second, t.logical) << '\n';
    }
    out << "GAP level         Blocks       % GAP       Actual words    Capacity words\n";
    for (unsigned i = 0; i < levels; ++i)
        out << std::setw(9) << i << std::setw(15) << t.gap_level[i]
            << std::setw(12) << pct(t.gap_level[i], t.gap)
            << std::setw(19) << t.gap_level_words[i]
            << std::setw(18) << t.first_levels[i] << '\n';
    out << "GAP words (actual lengths): " << t.gap_words
        << "  GAP payload bytes (2*words): " << t.gap_words * 2 << '\n'
        << "BIT payload bytes (8192*blocks): " << t.bit * 8192 << '\n'
        << "Vectors with BIT: " << t.vectors_bit << "  GAP: " << t.vectors_gap
        << "  both: " << t.vectors_mixed << "  GAP level 2 or 3: " << t.vectors_high_gap << '\n'
        << "GAP level-table differences across vectors: " << t.level_table_changes << '\n';
    out << "Per-vector counts        p50       p90       p99\n";
    for (auto item : {std::pair<const char*, std::vector<count_t>*>("BIT", &t.bit_per_vector),
                      {"GAP", &t.gap_per_vector}, {"GAP L2+L3", &t.high_gap_per_vector}})
        out << std::left << std::setw(19) << item.first << std::right
            << std::setw(10) << percentile(*item.second, 50)
            << std::setw(10) << percentile(*item.second, 90)
            << std::setw(10) << percentile(*item.second, 99) << '\n';
    out << "Exact GAP length x level histogram: gap_lengths.csv\n"
        << "Per-vector block counts: per_vector.csv\n"
        << "Resident payload bytes above are structural, not serialized sizes.\n";
    out.close();
}
} // namespace

int main(int argc, char** argv) {
    try {
        const auto o = parse(argc, argv);
        std::ifstream input(o.input, std::ios::binary | std::ios::ate);
        if (!input) throw std::runtime_error("Cannot open input collection");
        const auto file_size = count_t(input.tellg());
        input.seekg(0);
        std::ofstream per_vector(o.out + "/per_vector.csv");
        if (!per_vector) throw std::runtime_error("Cannot write per_vector.csv");
        per_vector << "vector,ids,last_id,logical_blocks,zero_blocks,full_blocks,bit_blocks,gap_blocks,gap_level0,gap_level1,gap_level2,gap_level3,gap_words,max_gap_words\n";
        totals t;
        std::vector<unsigned> ids;
        BM_DECLARE_TEMP_BLOCK(tb)
    const auto start = scan_clock::now();
        while (t.input_bytes < file_size && (!o.max_vectors || t.vectors < o.max_vectors)) {
            std::uint32_t len = 0;
            if (file_size - t.input_bytes < sizeof(len) ||
                !input.read(reinterpret_cast<char*>(&len), sizeof(len)) ||
                !len || count_t(len) * sizeof(unsigned) > file_size - t.input_bytes - sizeof(len))
                throw std::runtime_error("Truncated or invalid input record at vector " + std::to_string(t.vectors));
            ids.resize(len);
            if (!input.read(reinterpret_cast<char*>(ids.data()), std::streamsize(count_t(len) * sizeof(unsigned))))
                throw std::runtime_error("Cannot read input IDs");
            t.input_bytes += sizeof(len) + count_t(len) * sizeof(unsigned);
            for (size_t i = 1; i < ids.size(); ++i)
                if (ids[i] <= ids[i-1])
                    throw std::runtime_error("Input IDs are not strictly sorted at vector " + std::to_string(t.vectors));
            bm::bvector<> bv;
            bv.set(ids.data(), bm::bvector<>::size_type(ids.size()), bm::BM_SORTED);
            bv.optimize(tb);
            bm::bvector<>::statistics st;
            bv.calc_stat(&st);
            if (!t.vectors)
                std::copy(st.gap_levels, st.gap_levels + levels, t.first_levels);
            else if (!std::equal(st.gap_levels, st.gap_levels + levels, t.first_levels))
                ++t.level_table_changes;
            inspect(bv, st, ids.back(), t.vectors, ids.size(), t, per_vector);
            t.ids += ids.size();
            ++t.vectors;
            if (o.progress && t.vectors % o.progress == 0) {
                const double seconds = std::chrono::duration<double>(scan_clock::now() - start).count();
                std::cerr << "Progress vectors=" << t.vectors << " input_bytes=" << t.input_bytes
                          << " bit=" << t.bit << " gap=" << t.gap
                          << " full=" << t.full << " elapsed_s=" << seconds << '\n';
            }
        }
        if (!t.vectors) throw std::runtime_error("Input contains no vectors");
        if (o.max_vectors && t.vectors != o.max_vectors)
            throw std::runtime_error("Input ended before requested vector limit");
        per_vector.close();
        const double seconds = std::chrono::duration<double>(scan_clock::now() - start).count();
        report(o, t, seconds);
        std::cout << "Block statistics complete: vectors=" << t.vectors
                  << " bit=" << t.bit << " gap=" << t.gap
                  << " full=" << t.full << " zero=" << t.zero
                  << " elapsed_s=" << seconds << '\n';
    } catch (const std::exception& e) {
        std::cerr << "Block statistics error: " << e.what() << '\n';
        return 1;
    }
}
