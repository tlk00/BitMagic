#ifndef BMBLOCKSELECTION__H__INCLUDED__
#define BMBLOCKSELECTION__H__INCLUDED__
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

/*! \file bmblockselection.h
    \brief Block serialization candidate selection diagnostics (internal).
    Include after bm.h or bm64.h, as for the other serialization utilities.
*/

#include "bmbuffer.h"

namespace bm
{

/**
    \brief Complete-record size observations for one physical bit-block.
    \ingroup bvserial
    \internal

    Each independent candidate includes its ordinary-codec fallback. The kind
    identifies the selected representation: 0 ordinary, 1 periodic step 1,
    2 whole-block BWT, 3 BWT followed by periodic step 1, 4 BWT followed by
    greedy periodic steps (or independent step 2 in its separate experiment),
    and 5 BWT with direct BIC encoding (opcode 96). Sizes are bytes unless noted.
*/
struct block_selection_observation
{
    id64_t vector = 0;                    ///< Input vector ordinal in the serializer run.
    id64_t block = 0;                     ///< Physical bit-block index in the input vector.
    unsigned baseline = 0;               ///< Ordinary-codec record size.
    unsigned periodic = 0;               ///< Independent periodic step-1 size, including fallback.
    unsigned bwt = 0;                    ///< Independent BWT size, including fallback.
    unsigned selected = 0;               ///< Final emitted record size.
    unsigned kind = 0;                   ///< Selected representation kind (0..5).
    unsigned trials = 0;                 ///< Independent periodic residual encoding trials.
    unsigned bwt_payload = 0;            ///< BWT-transformed payload size without outer framing.
    unsigned bwt_candidate = 0;          ///< Complete BWT candidate size before ordinary fallback.
    unsigned previous = 0;               ///< Best size before the composed BWT/periodic trial.
    unsigned composed = 0;               ///< Composed BWT/periodic size, including fallback.
    unsigned compose_trials = 0;         ///< Residual encoding trials for the composed candidate.
    unsigned fourway = 0;                ///< Best size before greedy multistep or step-2 trials.
    unsigned multistep = 0;              ///< Greedy multistep BWT candidate size, including fallback.
    unsigned multistep_trials = 0;       ///< Greedy multistep residual encoding trials.
    unsigned multistep_groups = 0;       ///< Selected arithmetic-progression groups in the greedy trial.
    unsigned multistep_chains = 0;       ///< Selected progression chains in the greedy trial.
    unsigned multistep_last_step = 0;    ///< Last progression step visited by greedy search.
    unsigned multistep_stop = 0;         ///< Greedy stop code: 0 end, 1 marginal, 2 budget, 3 chain cap.
    unsigned multistep_formula = 0;      ///< Selected greedy formula size.
    unsigned multistep_residual = 0;     ///< Selected greedy residual size.
    id64_t multistep_mask = 0;           ///< Selected progression steps, with bit (step-1) set.
    unsigned step2 = 0;                  ///< Independent step-2 BWT candidate size, including fallback.
    unsigned step2_trials = 0;           ///< Independent step-2 residual encoding trials.
    unsigned step2_chains = 0;           ///< Selected progression chains in the step-2 trial.
    unsigned step2_formula = 0;          ///< Selected step-2 formula size.
    unsigned step2_residual = 0;         ///< Selected step-2 residual size.
    unsigned prethreshold_selected = 0;  ///< Best size before applying the minimum BWT saving.
    unsigned prethreshold_kind = 0;      ///< Winner kind before applying the minimum BWT saving.
};

/**
    \brief Accumulated block-selection measurements and per-block observations.
    \ingroup bvserial
    \internal

    \tparam BV Application bvector type providing allocator_type.

    Storage is allocated only when observations are appended. Reserve a known
    block count before collection or merging to avoid repeated allocations.
    Unknown-length collection grows geometrically and retains capacity.
    Observation storage uses BV::allocator_type; scalar counters
    and configuration flags are initialized to zero and false, respectively.
*/
template<class BV>
struct block_selection_statistics
{
    typedef typename BV::allocator_type allocator_type; ///< Application bvector allocator.
    typedef bm::heap_vector<block_selection_observation,
                            allocator_type, true> observations_vect_type; ///< BM observation storage.

    observations_vect_type blocks;       ///< Per-block observations in collection order.
    id64_t baseline_bytes = 0;           ///< Sum of ordinary-codec record sizes.
    id64_t selected_bytes = 0;           ///< Sum of emitted record sizes.
    id64_t saved = 0;                    ///< Bytes saved relative to ordinary serialization.
    id64_t chosen[6] = {};               ///< Emitted block counts for representation kinds 0..5.
    id64_t periodic_ns = 0;              ///< Time spent in independent periodic trials, nanoseconds.
    id64_t bwt_ns = 0;                   ///< Time spent in independent BWT trials, nanoseconds.
    id64_t trials = 0;                   ///< Total independent periodic residual encoding trials.
    id64_t compose_ns = 0;               ///< Time spent in composed BWT/periodic trials, nanoseconds.
    id64_t compose_trials = 0;           ///< Total composed residual encoding trials.
    id64_t previous_bytes = 0;           ///< Sum of best sizes before composition.
    id64_t multistep_ns = 0;             ///< Time spent in greedy multistep trials, nanoseconds.
    id64_t multistep_trials = 0;          ///< Total greedy multistep residual encoding trials.
    id64_t fourway_bytes = 0;            ///< Sum of best sizes before multistep or step-2 trials.
    id64_t step2_ns = 0;                 ///< Time spent in independent step-2 trials, nanoseconds.
    id64_t step2_trials = 0;             ///< Total independent step-2 residual encoding trials.
    id64_t formula_trial_bits[4] = {};    ///< Trial formula bits: header, starts, lengths, padding.
    id64_t formula_candidate_bits[4] = {}; ///< Candidate formula bits in the same four categories.
    id64_t formula_trial_modes[2] = {};   ///< Trial formula counts: gamma lengths, fixed interleaving.
    id64_t formula_candidate_modes[2] = {}; ///< Candidate formula counts for the same two modes.
    id64_t threshold_rejected = 0;       ///< BWT wins rejected by the minimum saving requirement.
    id64_t threshold_lost_bytes = 0;     ///< Savings forgone due to those rejected BWT wins.
    unsigned bwt_min_gain = 0;           ///< Required complete-record BWT saving, in bytes.
    bool composed = false;              ///< Composed BWT/periodic measurements are present.
    bool multistep = false;             ///< Greedy multistep measurements are present.
    bool step2 = false;                 ///< Independent step-2 measurements are present.

    /**
        \brief Reserve observation capacity without changing the observation count.
        \param count Minimum number of observations to accommodate.
        \note Existing observations survive allocation; capacity never shrinks.
            Growth is geometric when count is not known in advance.
    */
    void reserve_blocks(size_t count);

    /**
        \brief Append an observation, preserving capacity for subsequent blocks.
        \param observation Fully initialized observation with kind in 0..5 and
            selected size no greater than baseline. May refer to an existing block.
        \note Only appends the record; the serializer updates aggregate counters.
    */
    void add_block(const block_selection_observation& observation);
};

/**
    \brief Append observations and accumulate block-selection measurements.
    \ingroup bvserial
    \internal
    \tparam BV Application bvector type providing the observation allocator.
    \param d Destination accumulator; observations and counters are preserved.
    \param s Source measurements, appended in their existing order. Self-addition
        is supported and duplicates observations and additive counters.
    \note Configuration flags are combined with OR. A nonzero source bwt_min_gain
        replaces the destination setting, retaining the existing merge policy.
        Observation capacity is reserved before copying the batch.
*/
template<class BV>
inline void block_selection_add(block_selection_statistics<BV>& d,
                                const block_selection_statistics<BV>& s);


// -----------------------------------------------------------------------
// Implementation
// -----------------------------------------------------------------------

template<class BV>
inline void block_selection_statistics<BV>::reserve_blocks(size_t count)
{
    const size_t limit = size_t(-1) / sizeof(block_selection_observation);
    BM_ASSERT(count <= limit);
    const size_t capacity = blocks.capacity();
    if (count <= capacity)
        return;
    size_t next = capacity ? capacity : 64;
    if (next < count)
    {
        next = next <= limit / 2 ? next * 2 : limit;
        if (next < count)
            next = count;
    }
    blocks.reserve(next);
    BM_ASSERT(blocks.capacity() >= count);
}

template<class BV>
inline void block_selection_statistics<BV>::add_block(
                            const block_selection_observation& observation)
{
    BM_ASSERT(observation.kind < 6 && observation.prethreshold_kind < 6);
    BM_ASSERT(observation.selected <= observation.baseline);
    BM_ASSERT(blocks.size() < size_t(-1) / sizeof(block_selection_observation));
    // Preserve references into our own storage if reserving relocates it.
    const block_selection_observation value = observation;
    reserve_blocks(blocks.size() + 1);
    blocks.push_back_no_check(value);
}

template<class BV>
inline void block_selection_add(block_selection_statistics<BV>& d,
                                const block_selection_statistics<BV>& s)
{
    BM_ASSERT(d.selected_bytes <= d.baseline_bytes);
    BM_ASSERT(s.selected_bytes <= s.baseline_bytes);
    const size_t count = s.blocks.size();
    const size_t offset = d.blocks.size();
    BM_ASSERT(count <= size_t(-1) / sizeof(block_selection_observation) - offset);
    d.reserve_blocks(offset + count);
    d.blocks.resize(offset + count);
    for (size_t i = 0; i < count; ++i)
    {
        BM_ASSERT(s.blocks[i].kind < 6 && s.blocks[i].prethreshold_kind < 6);
        BM_ASSERT(s.blocks[i].selected <= s.blocks[i].baseline);
        d.blocks[offset+i] = s.blocks[i];
    } // for i
    d.baseline_bytes += s.baseline_bytes;
    d.selected_bytes += s.selected_bytes;
    d.saved += s.saved;
    d.trials += s.trials;
    d.periodic_ns += s.periodic_ns;
    d.bwt_ns += s.bwt_ns;
    d.compose_ns += s.compose_ns;
    d.compose_trials += s.compose_trials;
    d.previous_bytes += s.previous_bytes;
    d.composed |= s.composed;
    d.multistep_ns += s.multistep_ns;
    d.multistep_trials += s.multistep_trials;
    d.fourway_bytes += s.fourway_bytes;
    d.multistep |= s.multistep;
    d.step2_ns += s.step2_ns;
    d.step2_trials += s.step2_trials;
    d.step2 |= s.step2;
    for (unsigned i = 0; i < 4; ++i)
    {
        d.formula_trial_bits[i] += s.formula_trial_bits[i];
        d.formula_candidate_bits[i] += s.formula_candidate_bits[i];
    } // for i
    for (unsigned i = 0; i < 2; ++i)
    {
        d.formula_trial_modes[i] += s.formula_trial_modes[i];
        d.formula_candidate_modes[i] += s.formula_candidate_modes[i];
    } // for i
    d.threshold_rejected += s.threshold_rejected;
    d.threshold_lost_bytes += s.threshold_lost_bytes;
    if (s.bwt_min_gain)
        d.bwt_min_gain = s.bwt_min_gain;
    for (unsigned i = 0; i < 6; ++i)
    {
        d.chosen[i] += s.chosen[i];
    } // for i
}


} // namespace bm
#endif
