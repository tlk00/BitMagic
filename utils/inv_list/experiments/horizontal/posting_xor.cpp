// Research-only horizontal XOR survey. All codec trials use disposable memory.
#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <cstring>
#include <functional>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

#include "bm.h"
#include "bmserial.h"
#include "../../inv_list_serial_config.h"
#include "group_profile.h"

namespace
{
using BV = bm::bvector<>;
using Mask = bm::id64_t;
using Clock = std::chrono::steady_clock;
using Counts = std::map<std::string, std::uint64_t>;
constexpr unsigned waves = bm::block_waves;
constexpr unsigned words = bm::set_block_digest_wave_size;
constexpr unsigned width = bm::set_sub_array_size;
constexpr Mask all = ~Mask(0);

/// @internal Exposes the existing block selector only inside this survey.
class TrialEncoder : public bm::serializer<BV>
{
public:
    /// Measure one block in context s at slot j; residual enables BIT-to-GAP trials.
    size_t measure(const bm::word_t* p, const bm::bv_sub_survey& s,
                   unsigned j, bool residual);
private:
    std::array<unsigned char, 66560> bytes_; ///< Reused single-record trial buffer.
};

size_t TrialEncoder::measure(const bm::word_t* p, const bm::bv_sub_survey& s,
                             unsigned j, bool residual)
{
    bm::encoder enc(bytes_.data(), bytes_.size());
    const unsigned type = encode_selected_block(p, s, j, enc, residual);
    return enc.size() + (type == bm::set_block_azero || type == bm::set_block_aone);
}

/// @internal Ranked formula for a single source.
struct Candidate
{
    unsigned score; ///< Heuristic gain, used only within one scorer.
    unsigned slot;  ///< Earlier source block.
    Mask mask;      ///< Selected aligned waves.
};

/// Compare sparse and dense masked residual populations for target a and source b.
Candidate population(const bm::word_t* a, const bm::word_t* b, unsigned slot,
                     const bm::block_waves_xor_descr& d, unsigned ones)
{
    unsigned fewer = 0, more = 0;
    Mask sparse = 0, dense = 0;
    for (unsigned w = 0; w < waves; ++w)
    {
        unsigned x = 0;
        for (unsigned k = 0; k < words; ++k)
        {
            x += bm::word_bitcount(a[w * words + k] ^ b[w * words + k]);
        } // for k
        unsigned n = d.sb_bc[w];
        if (x < n)
        {
            fewer += n - x;
            sparse |= Mask(1) << w;
        }
        if (x > n)
        {
            more += x - n;
            dense |= Mask(1) << w;
        }
    } // for w
    unsigned zeros = bm::bits_in_block - ones;
    unsigned initial = std::min(ones, zeros);
    if (ones - fewer <= zeros - more)
    {
        return {initial - (ones - fewer), slot, sparse};
    }
    return {initial - (zeros - more), slot, dense};
}

/// Emit the four best distinct sources and also trial their full masks.
void shortlist(std::vector<Candidate>& ranked, unsigned scorer,
               std::map<std::pair<unsigned, Mask>, unsigned>& trials)
{
    std::stable_sort(ranked.begin(), ranked.end(), [](const Candidate& a, const Candidate& b)
    {
        return a.score > b.score;
    });
    for (size_t i = 0; i < std::min<size_t>(4, ranked.size()); ++i)
    {
        const auto& c = ranked[i];
        trials[{c.slot, c.mask}] |= scorer;
        trials[{c.slot, all}] |= scorer;
    } // for i
}

/// Count exact nonconstant waves with an earlier source; keys borrow group memory.
void repeats(const bm::word_t* a, Counts& counts,
             std::unordered_map<std::string_view, Mask>& seen)
{
    std::vector<std::pair<std::string_view, unsigned>> pending;
    for (unsigned w = 0; w < waves; ++w)
    {
        const auto* p = a + w * words;
        if ((p[0] == 0 || p[0] == ~bm::word_t(0)) &&
            std::all_of(p + 1, p + words, [p](auto v) { return v == p[0]; }))
        {
            continue;
        }
        ++counts["nonconstant_waves"];
        std::string_view key(reinterpret_cast<const char*>(p), words * sizeof(*p));
        Mask bit = Mask(1) << w;
        auto it = seen.find(key);
        if (it != seen.end() && (it->second & ~bit))
        {
            ++counts["different_offset_repeat_waves"];
            counts["different_offset_only_waves"] += !(it->second & bit);
        }
        pending.emplace_back(key, w);
    } // for w
    for (auto entry : pending)
    {
        seen[entry.first] |= Mask(1) << entry.second;
    } // for entry
}

/// Inspect one rebased, optimized group. Accumulate projected savings against
/// its actual level-10 group stream and retain only numeric results.
void inspect(BV& group, TrialEncoder& direct, bm::serializer<BV>& normal,
             TrialEncoder::buffer& buffer, Counts& counts, bm::word_t* residual)
{
    normal.serialize(group, buffer);
    size_t baseline = buffer.size();
    counts["isolated_level10_bytes"] += baseline;
    ++counts["groups"];
    bm::bv_sub_survey s{};
    s.init(&group, 0);
    group.get_blocks_manager().is_sparse_sblock(0, 0, s);
    std::vector<unsigned> slots;
    std::unordered_map<std::string_view, Mask> seen;
    for (unsigned j = 0; j < width; ++j)
    {
        const auto* p = group.get_blocks_manager().get_block_ptr(0, j);
        if (!IS_VALID_ADDR(p))
        {
            continue;
        }
        if (BM_IS_GAP(p))
        {
            auto g = BMGAP_PTR(p);
            unsigned runs = bm::gap_length(g) - 1;
            std::string key = "gap_L" + std::to_string(bm::gap_level(g)) + "_runs_" +
                (runs <= 3 ? "1_3" : runs <= 8 ? "4_8" : runs <= 32 ? "9_32" : "33_plus");
            ++counts[key + "_count"];
            counts[key + "_standalone_bytes"] += direct.measure(p, s, j, false);
            continue;
        }
        slots.push_back(j);
        repeats(p, counts, seen);
    } // for j
    counts["bit_blocks"] += slots.size();
    if (slots.size() < 2)
    {
        return;
    }
    ++counts["searchable_groups"];
    direct.serialize(group, buffer);
    const size_t constituent = buffer.size();
    counts["searchable_baseline_bytes"] += baseline;
    counts["searchable_constituent_bytes"] += constituent;
    const size_t control = std::min(baseline, constituent);
    counts["ungrouping_only_saved_bytes"] += baseline - control;
    std::array<size_t, 3> savings{};
    bm::xor_scanner<BV> scanner;
    for (size_t t = 1; t < slots.size(); ++t)
    {
        const unsigned j = slots[t];
        const auto* a = s.blocks[j];
        scanner.compute_s_block_stats(a);
        auto& d = scanner.get_descr();
        Mask digest = bm::calc_block_digest0(a);
        counts["pairs"] += t;
        std::vector<Candidate> ranked[2];
        for (unsigned scorer = 0; scorer < 2; ++scorer)
        {
            auto start = Clock::now();
            for (size_t r = 0; r < t; ++r)
            {
                unsigned k = slots[r];
                auto* b = s.blocks[k];
                Candidate c{};
                if (!scorer)
                {
                    c = population(a, b, k, d, s.bc_arr[j]);
                }
                else
                {
                    bm::block_xor_match_descr m{};
                    scanner.compute_xor_complexity_descr(a, digest, b, d, m);
                    c = {m.block_gain, k, m.xor_d64};
                }
                if (c.score && c.mask)
                {
                    ranked[scorer].push_back(c);
                }
            } // for r
            auto elapsed = Clock::now() - start;
            counts[scorer ? "run_scoring_ns" : "pop_scoring_ns"] +=
                std::uint64_t(std::chrono::duration_cast<std::chrono::nanoseconds>(elapsed).count());
        } // for scorer
        std::map<std::pair<unsigned, Mask>, unsigned> trials;
        shortlist(ranked[0], 1, trials);
        shortlist(ranked[1], 2, trials);
        const size_t plain = direct.measure(a, s, j, false);
        std::array<size_t, 3> best{plain, plain, plain};
        size_t metadata = 0;
        unsigned distance = 0;
        unsigned selected_waves = 0;
        const char* layout = nullptr;
        for (const auto& trial : trials)
        {
            unsigned k = trial.first.first;
            Mask mask = trial.first.second;
            bool copy = mask == all && !std::memcmp(a, s.blocks[k], sizeof(*a) * bm::set_block_size);
            size_t header = mask == all ? 2 : 10;
            size_t cost = header;
            if (!copy)
            {
                bm::bit_block_xor(residual, a, s.blocks[k], mask);
                cost += direct.measure(residual, s, j, true);
            }
            ++counts["formula_trials"];
            for (unsigned i = 0; i < 2; ++i)
            {
                if (trial.second & (1u << i))
                {
                    best[i] = std::min(best[i], cost);
                }
            } // for i
            if (cost < best[2])
            {
                best[2] = cost;
                metadata = header;
                distance = j - k;
                selected_waves = bm::word_bitcount64(mask);
                layout = copy ? "copy_wins" : mask == all ? "whole_wins" : "masked_wins";
            }
        } // for trial
        if (layout)
        {
            ++counts[layout];
            ++counts["winning_distance_" + std::to_string(distance)];
            ++counts["winning_waves_" + std::to_string(selected_waves)];
            counts["winning_formula_bytes"] += metadata;
        }
        for (unsigned i = 0; i < 3; ++i)
        {
            savings[i] += plain - best[i];
        } // for i
    } // for t
    const char* names[] = {"pop", "run", "union"};
    for (unsigned i = 0; i < 3; ++i)
    {
        if (savings[i] > constituent)
        {
            throw std::logic_error("Block savings exceed constituent stream size");
        }
        size_t projected = constituent - savings[i];
        std::string prefix(names[i]);
        counts[prefix + "_block_saved_bytes"] += savings[i];
        counts[prefix + "_group_saved_bytes"] += projected < control ? control - projected : 0;
        counts[prefix + "_winning_groups"] += projected < control;
    } // for i
}

/// Read sampled records supplied by the runner: ordinal, count and sorted u32 IDs.
/// Print a two-column metric table on success. Any partial/malformed record fails.
int run(bool compare_groups)
{
    TrialEncoder direct;
    bm::serializer<BV> normal;
    inv_list_serial_options options;
    options.level = 10;
    configure_inv_serializer(normal, options, [](auto&) {});
    configure_inv_serializer(direct, options, [](auto&) {});
    direct.set_sparse_cutoff(0);
    direct.set_sblock_large_options(false);
    TrialEncoder::buffer buffer;
    Counts counts{{"different_offset_repeat_waves", 0}, {"different_offset_only_waves", 0},
                  {"copy_wins", 0}, {"whole_wins", 0}, {"masked_wins", 0}};
    BM_DECLARE_TEMP_BLOCK(residual)
    bm::serialization_opcode_collector<BV> collector;
    if (compare_groups)
    {
        std::cout << "ordinal\tsuperblock\tids\toccupied\tbit\tgap\tone_runs\tmax_slot_ids"
                     "\tgroup_opcode\tbaseline_B\tconstituent_B\tbase_nonrecord_B\talt_nonrecord_B"
                     "\tbase_records\talt_records\n";
    }
    unsigned ordinal, n;
    while (std::cin.peek() != EOF)
    {
        if (!std::cin.read(reinterpret_cast<char*>(&ordinal), 4) ||
            !std::cin.read(reinterpret_cast<char*>(&n), 4) || !n)
        {
            throw std::runtime_error("Invalid sampled record header");
        }
        std::vector<unsigned> ids(n);
        if (!std::cin.read(reinterpret_cast<char*>(ids.data()), std::streamsize(n) * 4) ||
            std::adjacent_find(ids.begin(), ids.end(), std::greater_equal<unsigned>()) != ids.end() ||
            ids.back() == bm::id_max)
        {
            throw std::runtime_error("Invalid sampled IDs at ordinal " + std::to_string(ordinal));
        }
        for (size_t start = 0; start < n;)
        {
            constexpr unsigned span = width * bm::bits_in_block;
            unsigned base = ids[start] / span * span;
            size_t end = start;
            while (end < n && ids[end] / span == base / span)
            {
                ids[end++] -= base;
            }
            BV group;
            group.set(ids.data() + start, unsigned(end - start), bm::BM_SORTED);
            group.optimize(residual);
            if (compare_groups)
            {
                profile_group(group, ordinal, base / span, ids.data() + start,
                              unsigned(end - start), normal, direct, buffer, collector);
            }
            else
            {
                inspect(group, direct, normal, buffer, counts, residual);
            }
            start = end;
        } // for start
        counts["ids"] += n;
        if (++counts["vectors"] % 500 == 0)
        {
            std::cerr << "Survey vectors=" << counts["vectors"];
            if (!compare_groups)
            {
                std::cerr << " pairs=" << counts["pairs"];
            }
            std::cerr << '\n';
        }
    } // while sampled records remain
    if (!counts["vectors"] || std::cin.bad())
    {
        throw std::runtime_error("Empty or unreadable sample");
    }
    if (compare_groups)
    {
        return 0;
    }
    for (const auto& entry : counts)
    {
        std::cout << entry.first << '\t' << entry.second << '\n';
    } // for entry
    return 0;
}
} // namespace

/// Run the standalone survey, propagating input or measurement failures.
int main(int argc, char** argv)
{
    try
    {
        bool compare_groups = argc == 2 && std::string_view(argv[1]) == "--group-profile";
        if (argc != 1 && !compare_groups)
        {
            throw std::runtime_error("Usage: gov2_posting_xor [--group-profile]");
        }
        return run(compare_groups);
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
