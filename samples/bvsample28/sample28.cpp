/*
Copyright(c) 2026 Anatoliy Kuznetsov

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy at http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

/** \example sample28.cpp
    Estimate structural compressibility of bit-vectors and their XOR residuals
    with bm::calc_complexity() and bm::calc_complexity_xor(). Compare every
    pair, choose the cheaper vector as the lead, and rank reference candidates.
    \sa bm::bvector_complexity_statistics
    \sa bm::calc_complexity
    \sa bm::calc_complexity_xor
*/
/*! \file sample28.cpp
    \brief Pairwise XOR-complexity analysis and lead-vector selection.
*/

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "bm.h"
#include "bmcomplexity.h"
#include "bmserial.h"

using bvector_type = bm::bvector<>;
using complexity_type = bm::bvector_complexity_statistics;

namespace
{

const bvector_type::size_type block_size = bm::gap_max_bits;
const bvector_type::size_type demo_last = block_size * 8 - 1;

struct demo_vector
{
    explicit demo_vector(const char* vector_name) : name(vector_name) {}

    std::string name;
    bvector_type bv;
    complexity_type complexity;
    bvector_type::statistics storage;
    bm::id64_t cost = 0;
};

struct pair_result
{
    std::size_t first = 0;
    std::size_t second = 0;
    complexity_type residual;
    bm::id64_t residual_cost = 0;
    unsigned score = 0;
    std::size_t lead = 0;
    bm::id64_t independent_cost = 0;
    bm::id64_t referenced_cost = 0;
};

/*
This is deliberately a structural heuristic, not a serialized-byte estimator.
One dense BIT block contributes its 8192-byte payload. Predicted GAP blocks
contribute their used 16-bit words, including headers, without capacity slack.
Uniform block spans retain a heuristic one-unit overhead. Runs in BIT and
full blocks are not charged again: their representation is already covered.

Keeping the formula in one function makes the policy visible and easy to tune.
The statistics themselves remain exact even if this demonstration policy is
later replaced with measured serialization sizes or a trained predictor.
*/
bm::id64_t complexity_cost(const complexity_type& s)
{
    if (!s.count)
        return 0;
    return s.bit_blocks * 8192u +
           s.gap_words * sizeof(bm::gap_word_t) +
           s.zero_spans + s.full_spans;
}

/* Return symmetric XOR compatibility in [0, 100].

The cheaper operand is the lead, so the residual replaces the more expensive
standalone operand. This remains independent of operand order. An empty XOR
scores 100. A residual which does not improve pair storage scores zero.
*/
unsigned compatibility_score(bm::id64_t cost_a, bm::id64_t cost_b,
                             bm::id64_t xor_cost)
{
    const bm::id64_t baseline = (std::max)(cost_a, cost_b);
    if (!xor_cost)
        return 100;
    if (!baseline || xor_cost >= baseline)
        return 0;
    const double improvement = double(baseline - xor_cost) / double(baseline);
    // A real, but sub-percent, predicted saving should not be reported as the
    // zero (no-improvement) category. Truncation also reserves 100 for the
    // exact empty-residual case handled above.
    const unsigned score = unsigned(improvement * 100.0);
    return score ? score : 1;
}

void add_regular_islands(bvector_type& bv)
{
    bv.set_range(1000, 9000);
    bv.set_range(block_size - 50, block_size + 120); // crosses a block boundary
    bv.set_range(block_size * 2 + 4000, block_size * 2 + 18000);
    bv.set_range(block_size * 4, block_size * 5 - 1); // one full block
    for (unsigned i = 0; i < 700; ++i)
        bv.set(block_size * 6 + i * 67);
}

void build_collection(std::vector<std::unique_ptr<demo_vector> >& collection)
{
    collection.emplace_back(new demo_vector("A-base"));
    add_regular_islands(collection.back()->bv);

    collection.emplace_back(new demo_vector("A-copy"));
    collection.back()->bv = collection[0]->bv;

    collection.emplace_back(new demo_vector("A-near"));
    collection.back()->bv = collection[0]->bv;
    collection.back()->bv.set_range(1200, 1280, false);
    collection.back()->bv.set_range(block_size * 3 + 10,
                                    block_size * 3 + 180);
    collection.back()->bv.set(block_size * 7 + 12345);

    collection.emplace_back(new demo_vector("A-complement"));
    collection.back()->bv.set_range(0, demo_last);
    collection.back()->bv ^= collection[0]->bv; // complement in demo universe

    collection.emplace_back(new demo_vector("B-base"));
    collection.back()->bv.set_range(20000, 40000);
    collection.back()->bv.set_range(block_size * 2 - 200,
                                    block_size * 2 + 200);
    collection.back()->bv.set_range(block_size * 3 + 2000,
                                    block_size * 3 + 50000);
    for (unsigned i = 0; i < 1200; ++i)
        collection.back()->bv.set(block_size * 7 + i * 41);

    collection.emplace_back(new demo_vector("B-near"));
    collection.back()->bv = collection[4]->bv;
    collection.back()->bv.set_range(20500, 20700, false);
    collection.back()->bv.set_range(block_size * 5 + 100,
                                    block_size * 5 + 900);

    collection.emplace_back(new demo_vector("sparse"));
    for (unsigned i = 0; i < 80; ++i)
        collection.back()->bv.set(bvector_type::size_type(i) * 6007 + 17);

    collection.emplace_back(new demo_vector("fragmented"));
    for (bvector_type::size_type i = 0; i <= demo_last; i += 3)
        collection.back()->bv.set(i);

    // Normalize physical representation before collecting allocation-sensitive
    // fields such as bit_to_gap_blocks.
    for (std::size_t i = 0; i < collection.size(); ++i)
        collection[i]->bv.optimize();
}

void measure_vectors(std::vector<std::unique_ptr<demo_vector> >& collection)
{
    std::cout << "\n1. Measuring optimized standalone vectors\n"
              << "   cost is the sample's structural heuristic; memory is the"
                 " actual in-memory allocation.\n\n";
    std::cout << std::left << std::setw(15) << "vector"
              << std::right << std::setw(10) << "count"
              << std::setw(8) << "runs" << std::setw(9) << "blk-runs"
              << std::setw(7) << "zero" << std::setw(7) << "full"
              << std::setw(7) << "GAP" << std::setw(7) << "BIT"
              << std::setw(10) << "GAP-words"
              << std::setw(8) << "zspan" << std::setw(8) << "fspan"
              << std::setw(12) << "cost" << std::setw(12) << "memory\n";

    for (std::size_t i = 0; i < collection.size(); ++i)
    {
        demo_vector& item = *collection[i];
        item.complexity = bm::calc_complexity(item.bv);
        item.bv.calc_stat(&item.storage);
        item.cost = complexity_cost(item.complexity);
        const complexity_type& s = item.complexity;
        std::cout << std::left << std::setw(15) << item.name << std::right
                  << std::setw(10) << s.count << std::setw(8) << s.runs
                  << std::setw(9) << s.block_runs
                  << std::setw(7) << s.zero_blocks
                  << std::setw(7) << s.full_blocks
                  << std::setw(7) << s.gap_blocks
                  << std::setw(7) << s.bit_blocks
                  << std::setw(10) << s.gap_words
                  << std::setw(8) << s.zero_spans
                  << std::setw(8) << s.full_spans
                  << std::setw(12) << item.cost
                  << std::setw(12) << item.storage.memory_used << '\n';
    }
}

std::vector<pair_result>
measure_pairs(const std::vector<std::unique_ptr<demo_vector> >& collection)
{
    std::vector<pair_result> pairs;
    std::cout << "\n2. Computing every unordered XOR pair\n"
              << "   XOR statistics and score are symmetric. The lead is the"
                 " cheaper standalone vector.\n\n";
    std::cout << std::left << std::setw(29) << "pair"
              << std::right << std::setw(10) << "xor-count"
              << std::setw(9) << "runs" << std::setw(10) << "blk-runs"
              << std::setw(7) << "GAP" << std::setw(7) << "BIT"
              << std::setw(10) << "GAP-words"
              << std::setw(9) << "BIT->GAP" << std::setw(11) << "xor-cost"
              << std::setw(8) << "score" << "  lead / decision\n";

    for (std::size_t i = 0; i < collection.size(); ++i)
        for (std::size_t j = i + 1; j < collection.size(); ++j)
        {
            pair_result p;
            p.first = i;
            p.second = j;
            p.residual = bm::calc_complexity_xor(collection[i]->bv,
                                                 collection[j]->bv);
            p.residual_cost = complexity_cost(p.residual);
            p.score = compatibility_score(collection[i]->cost,
                                          collection[j]->cost,
                                          p.residual_cost);
            p.lead = collection[i]->cost <= collection[j]->cost ? i : j;
            p.independent_cost = collection[i]->cost + collection[j]->cost;
            p.referenced_cost = collection[p.lead]->cost + p.residual_cost;
            pairs.push_back(p);

            const std::string label = collection[i]->name + " / " +
                                      collection[j]->name;
            std::cout << std::left << std::setw(29) << label << std::right
                      << std::setw(10) << p.residual.count
                      << std::setw(9) << p.residual.runs
                      << std::setw(10) << p.residual.block_runs
                      << std::setw(7) << p.residual.gap_blocks
                      << std::setw(7) << p.residual.bit_blocks
                      << std::setw(10) << p.residual.gap_words
                      << std::setw(9) << p.residual.bit_to_gap_blocks
                      << std::setw(11) << p.residual_cost
                      << std::setw(7) << p.score << "%  "
                      << collection[p.lead]->name << " / "
                      << (p.referenced_cost < p.independent_cost
                              ? "XOR improves" : "keep independent") << '\n';
        }
    return pairs;
}

std::vector<const pair_result*>
print_best_peers(const std::vector<std::unique_ptr<demo_vector> >& collection,
                      const std::vector<pair_result>& pairs)
{
    std::vector<const pair_result*> selected;
    std::cout << "\n3. Best peer for each vector\n"
              << "   This ranks the symmetric scores; ties prefer the lower"
                 " referenced pair cost.\n\n";
    for (std::size_t i = 0; i < collection.size(); ++i)
    {
        const pair_result* best = 0;
        for (std::size_t k = 0; k < pairs.size(); ++k)
        {
            const pair_result& p = pairs[k];
            if (p.first != i && p.second != i)
                continue;
            if (!best || p.score > best->score ||
                (p.score == best->score &&
                 p.referenced_cost < best->referenced_cost))
                best = &p;
        }
        if (!best || !best->score)
        {
            std::cout << "   " << collection[i]->name
                      << " -> no improving peer\n";
            continue;
        }
        const std::size_t peer = best->first == i ? best->second : best->first;
        if (std::find(selected.begin(), selected.end(), best) == selected.end())
            selected.push_back(best);
        std::cout << "   " << std::setw(15) << std::left << collection[i]->name
                  << " -> " << std::setw(15) << collection[peer]->name
                  << std::right << " score=" << std::setw(3) << best->score
                  << "%  lead=" << std::setw(15) << std::left
                  << collection[best->lead]->name << std::right
                  << " pair-cost " << best->independent_cost << " -> "
                  << best->referenced_cost << '\n';
    }
    return selected;
}

/** Measure actual BLOB sizes for the unique positive-score pairs from step 3.
    @param collection Optimized standalone vectors used in the predictions.
    @param selected Selected pairs; pointers remain valid while pairs lives.
    @param pairs All pairs, used to show zero-score negative controls.
    Each residual is materialized and optimized here. One serializer with
    default settings measures all BLOBs consistently. Pair totals exclude
    application framing and reference identifiers. The smaller serialized
    operand is the lead. Negative controls are selected by score, not outcome.
*/
void measure_serialization(
    const std::vector<std::unique_ptr<demo_vector> >& collection,
    const std::vector<const pair_result*>& selected,
    const std::vector<pair_result>& pairs)
{
    std::cout << "\n4. Actual serialization of selected pairs (unique pairs only)\n"
              << "   Default serializer settings; sizes in bytes. Materializing"
                 " and optimizing each XOR.\n"
              << "   Totals include BLOB headers, but exclude reference IDs and"
                 " collection framing.\n";
    bm::serializer<bvector_type> serializer;
    bm::serializer<bvector_type>::buffer buffer;
    std::vector<std::size_t> sizes(collection.size(), 0);
    std::vector<bool> measured(collection.size(), false);
    auto print_pair = [&](const pair_result* p)
    {
        for (std::size_t idx : {p->first, p->second})
            if (!measured[idx])
            {
                serializer.serialize(collection[idx]->bv, buffer);
                sizes[idx] = buffer.size();
                measured[idx] = true;
            }
        bvector_type residual(collection[p->first]->bv);
        residual ^= collection[p->second]->bv;
        residual.optimize();
        serializer.serialize(residual, buffer);
        const std::size_t xor_bytes = buffer.size();
        const std::size_t independent = sizes[p->first] + sizes[p->second];
        const std::size_t measured_lead = sizes[p->first] <= sizes[p->second]
                                       ? p->first : p->second;
        const std::size_t total = sizes[measured_lead] + xor_bytes;
        const double savings = independent
            ? 100.0 * (double(independent) - double(total)) / double(independent)
            : 0.0;
        const std::string label = collection[p->first]->name + " / " +
                                  collection[p->second]->name;
        std::cout << std::left << std::setw(30) << label << std::right
                  << std::setw(6) << p->score
                  << std::fixed << std::setprecision(2) << std::setw(11)
                  << savings << "%" << std::defaultfloat
                  << std::setw(10) << independent << " -> " << std::setw(8)
                  << total << "  " << (total < independent ? "gain" :
                                     (total > independent ? "loss" : "no gain"))
                  << "; lead=" << collection[measured_lead]->name << '\n';
    };
    auto print_header = []()
    {
        std::cout << std::left << std::setw(30) << "pair" << std::right
                  << std::setw(6) << "score" << std::setw(12) << "pair savings"
                  << "   independent -> lead+XOR bytes / result\n";
    };
    std::cout << "\n   --- Recommended best-peer pairs ---\n";
    print_header();
    for (const pair_result* p : selected)
        print_pair(p);
    if (selected.empty())
        std::cout << "   No positive-score candidates to serialize.\n";
    std::cout << "\n   --- NOT recommended: all zero-score pairs ---\n"
              << "   Negative savings mean expansion. A gain here is a missed"
                 " opportunity for the heuristic.\n";
    print_header();
    for (const pair_result& p : pairs)
        if (!p.score)
            print_pair(&p);
}

} // namespace

int main()
{
    try
    {
        std::cout << "BitMagic XOR complexity and reference-selection demo\n"
                  << "Demo universe: [0, " << demo_last << "] (8 blocks)\n"
                  << "Building related, complementary and unrelated vectors...\n";

        std::vector<std::unique_ptr<demo_vector> > collection;
        collection.reserve(8);
        build_collection(collection);
        measure_vectors(collection);
        const std::vector<pair_result> pairs = measure_pairs(collection);
        const auto selected = print_best_peers(collection, pairs);
        measure_serialization(collection, selected, pairs);

        std::cout << "\nDone. See readme.md for the formula, interpretation,"
                     " limitations and possible production policies.\n";
        return 0;
    }
    catch (const std::exception& ex)
    {
        std::cerr << "Error: " << ex.what() << '\n';
        return 1;
    }
}
