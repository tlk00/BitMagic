#ifndef BMCOMPLEXITY__H__INCLUDED__
#define BMCOMPLEXITY__H__INCLUDED__
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

/*! \file bmcomplexity.h
    \brief Content statistics for bit-vectors and unmaterialized XOR results.
    Include bm.h (or bm64.h) before this header.
*/
#ifndef BM__H__INCLUDED__
# error missing include (bm.h or bm64.h)
#endif

#include "bmxor.h"

namespace bm
{

/*! @defgroup bvcomplexity Bit-vector complexity statistics
    @brief Run-based complexity and block compression profiles for bit-vectors
    and their XOR differences.

    Algorithms in this group compute exact population and one-run counts,
    together with content-based block compression predictors. XOR profiles
    support reference selection without constructing a result bit-vector.
    Predictions complement actual storage statistics from bvector::calc_stat();
    they do not estimate the final serialized byte size.
    @ingroup bvector
*/

/*! Content-based compression predictors, independent of current allocation.

    All counts use the common universe [0, bm::id_max). The reserved bit at
    bm::id_max is treated as zero. Block classifications use physical 64K-bit
    blocks, including this zero padding: a completely set universe therefore
    ends in a mixed block, not a full block.

    gap_blocks and bit_blocks predict the default opt_compress policy for
    mixed blocks: GAP when the number of physical zero/one runs is strictly
    below gap_max_buff_len-4, BIT otherwise. They are not the stored block
    counts and do not predict serialized bytes. For actual allocation counts,
    memory consumption and custom GAP levels, retain bvector::calc_stat().

    zero_blocks + full_blocks + gap_blocks + bit_blocks equals the number of
    physical blocks in the universe, including implicit trailing zero blocks.
    All counters are 64-bit, including in a 32-bit-address build.
    @ingroup bvcomplexity
*/
struct bvector_complexity_statistics
{
    bm::id64_t count;       ///< Population count (XOR population for a pair).
    bm::id64_t runs;        ///< Maximal one-runs; zero for an empty vector.
    bm::id64_t block_runs;  ///< One-runs counted separately inside each block.
    bm::id64_t zero_blocks; ///< Entirely zero physical blocks, allocated or not.
    bm::id64_t full_blocks; ///< Entirely one physical blocks.
    bm::id64_t gap_blocks;  ///< Mixed blocks predicted to fit default GAP policy.
    bm::id64_t gap_words;   ///< Used 16-bit GAP words including headers, summed over predicted GAP blocks; excludes allocation slack.
    bm::id64_t bit_blocks;  ///< Remaining mixed blocks.
    bm::id64_t zero_spans;  ///< Maximal sequences of consecutive zero blocks.
    bm::id64_t full_spans;  ///< Maximal sequences of consecutive full blocks.
    bm::id64_t bit_to_gap_blocks; ///< XOR only: predicted GAP with either input BIT.

    bvector_complexity_statistics() BMNOEXCEPT
        : count(0), runs(0), block_runs(0), zero_blocks(0), full_blocks(0),
          gap_blocks(0), gap_words(0), bit_blocks(0), zero_spans(0), full_spans(0),
          bit_to_gap_blocks(0)
    {}
};

/*! @brief Implementation helpers for complexity statistics.
    @internal
*/
namespace complexity_detail
{
/*! Streaming accumulator. Carries the last bit and block class across skips.
    @internal
*/
class accumulator
{
public:
    /// @internal Accumulated statistics for the consumed block prefix.
    bvector_complexity_statistics stat;
    /// @internal Initialize an empty prefix with a zero predecessor bit.
    accumulator() BMNOEXCEPT : last_bit_(0), last_kind_(2) {}

    /*! Consume n complete uniform blocks without scanning or expanding them.
        @internal
        @param n Number of consecutive physical blocks to consume; may be zero.
        @param one True for all-one blocks, false for all-zero blocks.
    */
    void uniform(bm::id64_t n, bool one) BMNOEXCEPT
    {
        if (!n) return;
        if (one)
        {
            stat.count += n * bm::gap_max_bits;
            stat.block_runs += n;
            stat.runs += !last_bit_;
            stat.full_blocks += n;
            stat.full_spans += (last_kind_ != 1);
        }
        else
        {
            stat.zero_blocks += n;
            stat.zero_spans += (last_kind_ != 0);
        }
        last_bit_ = unsigned(one);
        last_kind_ = unsigned(one);
    }

    /*! Analyze one physical block. The caller has cleared the reserved bit.
        Run starts are ones whose predecessor is zero. Local counting starts
        with a zero predecessor; the global correction joins adjacent runs.
        The physical zero/one run count follows from the two endpoint bits.
        @internal
        @param p Non-null array of bm::set_block_size words to analyze.
        @param input_bit True for an XOR residual with at least one physically
                         BIT input; false for single-vector analysis.
    */
    void block(const bm::word_t* p, bool input_bit) BMNOEXCEPT
    {
        unsigned pop = 0, runs = 0, prev = 0;
        for (unsigned i = 0; i < bm::set_block_size; ++i)
        {
            bm::word_t w = p[i];
            pop += bm::word_bitcount(w);
            runs += bm::word_bitcount(w & ~((w << 1) | prev));
            prev = w >> 31;
        }
        unsigned first = p[0] & 1;
        unsigned all_runs = 2 * runs + 1 - first - prev;
        block_summary(pop, all_runs, first, prev, input_bit);
    }

    /*! Consume precomputed statistics for one physical block.
        @internal
        @param pop Population count of the block result.
        @param all_runs Number of alternating zero and one runs in the result.
        @param first Value of the first bit in the block.
        @param last Value of the last bit in the block.
        @param input_bit True when either XOR input is physically BIT.
    */
    void block_summary(unsigned pop, unsigned all_runs,
                       unsigned first, unsigned last,
                       bool input_bit) BMNOEXCEPT
    {
        if (!pop || pop == bm::gap_max_bits)
        {
            uniform(1, pop != 0);
            return;
        }
        BM_ASSERT(all_runs > 1);
        unsigned runs = (all_runs - 1 + first + last) >> 1;
        stat.count += pop;
        stat.block_runs += runs;
        stat.runs += runs - (last_bit_ & first);
        if (all_runs < bm::gap_max_buff_len - 4)
        {
            ++stat.gap_blocks;
            // Each physical zero/one run has one endpoint word; one more
            // word holds the GAP header. This also predicts GAP length for
            // BIT inputs without converting their result to GAP storage.
            stat.gap_words += all_runs + 1;
            stat.bit_to_gap_blocks += input_bit;
        }
        else
            ++stat.bit_blocks;
        last_bit_ = last;
        last_kind_ = 2;
    }
private:
    /// @internal Last consumed bit, or zero before the first block.
    unsigned last_bit_;
    /// @internal Previous block class: zero (0), full (1), or mixed/initial (2).
    unsigned last_kind_;
};

/*! Expand an operand block into reusable storage; never mutate the operand.
    @internal
    @param p Operand block: null, full-block marker, tagged GAP, or dense BIT.
    @param[out] dst Writable storage for bm::set_block_size words, suitably
                    aligned for block operations and not aliasing the operand.
*/
inline void expand(const bm::word_t* p, bm::word_t* dst) BMNOEXCEPT
{
    if (!p || IS_FULL_BLOCK(p))
    {
        bm::word_t value = p ? ~bm::word_t(0) : 0;
        for (unsigned k = 0; k < bm::set_block_size; ++k) dst[k] = value;
    }
    else if (BM_IS_GAP(p))
        bm::gap_convert_to_bitset(dst, BMGAP_PTR(p));
    else
        ::memcpy(dst, p, sizeof(bm::word_t) * bm::set_block_size);
}

/*! Shared single/pair traversal. A null second operand means single-vector
    analysis. Two stack blocks (16 KiB, suitably aligned) are reused for all
    mixed blocks. No heap allocation and no result bvector are required.
    Uniform/equal subtrees and blocks take constant work per visited entry.
    Mixed blocks are expanded and scanned; this deliberately simple first
    implementation does not yet merge GAP boundary streams directly.
    @internal
    @tparam BV Bit-vector type exposing the BitMagic block manager interface.
    @param a First operand, or the vector to analyze when b is null.
    @param b Second XOR operand in the common universe; null for single-vector
             analysis.
    @return Content statistics for a or for a XOR *b.
*/
template<class BV>
bvector_complexity_statistics scan(const BV& a, const BV* b) BMNOEXCEPT
{
    const typename BV::blocks_manager_type& am = a.get_blocks_manager();
    const typename BV::blocks_manager_type* bmgr = b ? &b->get_blocks_manager() : 0;
    const bm::id64_t blocks = (bm::id64_t(bm::id_max) >> bm::set_block_shift) + 1;
    unsigned tops = unsigned(am.top_block_size());
    if (bmgr && bmgr->top_block_size() > tops)
        tops = unsigned(bmgr->top_block_size());
    const unsigned max_tops = unsigned(blocks / bm::set_sub_array_size);
    if (tops > max_tops) tops = max_tops;

    bm::bit_block_t work_a, work_b;
    bm::word_t* wa = work_a;
    bm::word_t* wb = work_b;
    bm::gap_word_t* gap_tmp = static_cast<bm::gap_word_t*>(work_b);
    accumulator acc;
    for (unsigned i = 0; i < tops; ++i)
    {
        const bm::word_t* const* at = am.get_topblock(i);
        const bm::word_t* const* bt = bmgr ? bmgr->get_topblock(i) : 0;
        bool af = ((const bm::word_t*)at == FULL_BLOCK_FAKE_ADDR);
        bool bf = ((const bm::word_t*)bt == FULL_BLOCK_FAKE_ADDR);
        // The final superblock includes the reserved zero bit. Equal inputs
        // still cancel, but a uniform-one shortcut must not include that bit.
        if (at == bt || ((!at || af) && (!bt || bf) && i + 1 < max_tops))
        {
            acc.uniform(bm::set_sub_array_size, at != bt && (af != bf));
            continue;
        }
        for (unsigned j = 0; j < bm::set_sub_array_size; ++j)
        {
            const bm::word_t* ap = am.get_block_ptr(i, j);
            const bm::word_t* bp = bmgr ? bmgr->get_block_ptr(i, j) : 0;
            bool final = (i + 1 == max_tops && j + 1 == bm::set_sub_array_size);
            if (ap == bp)
            {
                acc.uniform(1, false);
                continue;
            }
            if ((!ap || IS_FULL_BLOCK(ap)) && (!bp || IS_FULL_BLOCK(bp)) && !final)
            {
                acc.uniform(1, bool(ap) != bool(bp));
                continue;
            }
            bool ag = ap && BM_IS_GAP(ap);
            bool bg = bp && BM_IS_GAP(bp);
            if (!b) // single-vector statistics
            {
                if (ag)
                {
                    const bm::gap_word_t* g = BMGAP_PTR(ap);
                    unsigned all_runs = bm::gap_length(g) - 1;
                    unsigned first = *g & 1u;
                    unsigned last = first ^ ((all_runs - 1) & 1u);
                    acc.block_summary(bm::gap_bit_count_unr(g), all_runs,
                                      first, last, false);
                    continue;
                }
                if (ap && !IS_FULL_BLOCK(ap))
                {
                    unsigned all_runs, pop;
                    bm::bit_block_change_bc(ap, &all_runs, &pop);
                    unsigned first = ap[0] & 1u;
                    unsigned last = ap[bm::set_block_size-1] >> 31;
                    acc.block_summary(pop, all_runs, first, last, false);
                    continue;
                }
            }
            bool input_bit = b &&
                ((ap && !IS_FULL_BLOCK(ap) && !BM_IS_GAP(ap)) ||
                 (bp && !IS_FULL_BLOCK(bp) && !BM_IS_GAP(bp)));

            if (ag && bg)
            {
                unsigned dsize;
                bm::gap_operation_xor(BMGAP_PTR(ap), BMGAP_PTR(bp),
                                      gap_tmp, dsize);
                (void)dsize;
                unsigned all_runs = bm::gap_length(gap_tmp) - 1;
                unsigned first = *gap_tmp & 1u;
                unsigned last = first ^ ((all_runs - 1) & 1u);
                acc.block_summary(bm::gap_bit_count_unr(gap_tmp), all_runs,
                                  first, last, false);
                continue;
            }
            if (ag || bg)
            {
                const bm::word_t* gap_block = ag ? ap : bp;
                const bm::gap_word_t* gap = BMGAP_PTR(gap_block);
                const bm::word_t* other = ag ? bp : ap;
                expand(other, wa);
                bm::gap_xor_to_bitset(wa, gap);
                if (final)
                    wa[bm::set_block_size-1] &= ~(bm::word_t(1) << 31);
                // TODO: Add a fused BIT/GAP XOR complexity kernel which
                // computes population and transitions without materializing
                // the dense XOR product in a temporary block.
                acc.block(wa, input_bit);
                continue;
            }
            if (ap && bp && !IS_FULL_BLOCK(ap) && !IS_FULL_BLOCK(bp))
            {
                unsigned all_runs, pop;
                bm::bit_block_xor_change(ap, bp, bm::set_block_size,
                                         &all_runs, &pop);
                unsigned first = (ap[0] ^ bp[0]) & 1u;
                unsigned last = (ap[bm::set_block_size-1] ^
                                 bp[bm::set_block_size-1]) >> 31;
                acc.block_summary(pop, all_runs, first, last, true);
                continue;
            }

            // Mixed dense/uniform blocks need one materialized block. The
            // final physical block also needs its reserved bit cleared.
            expand(ap, wa);
            expand(bp, wb);
            bm::bit_block_xor(wa, wb);
            if (final)
                wa[bm::set_block_size-1] &= ~(bm::word_t(1) << 31);
            // TODO: Add fused BIT/zero-full XOR complexity kernels (including
            // complement statistics) to avoid this remaining materialization.
            acc.block(wa, input_bit);
        }
    }
    acc.uniform(blocks - bm::id64_t(tops) * bm::set_sub_array_size, false);
    return acc.stat;
}
} // namespace complexity_detail

/*! Compute exact population/runs and a normalized block compression profile.
    Ignores current BIT/GAP allocation choices and bvector::size(): the common
    universe is [0, bm::id_max). Use calc_stat() alongside this function when
    caching actual storage statistics. bit_to_gap_blocks is zero for this API.
    Uses two reusable stack blocks; does not allocate or modify bv.
    @ingroup bvcomplexity
    @tparam BV Bit-vector type exposing the BitMagic block manager interface.
    @param bv Source bit-vector in the common universe [0, bm::id_max).
    @return Exact population/run counts and predicted block compression profile.
    @sa calc_complexity_xor, bvector_complexity_statistics
*/
template<class BV>
bvector_complexity_statistics calc_complexity(const BV& bv) BMNOEXCEPT
{
    return complexity_detail::scan(bv, static_cast<const BV*>(0));
}

/*! Compute the same profile for A XOR B without constructing its bvector.
    Every field is symmetric in A and B. runs is the exact run distance:
    half the Hamming distance between zero-padded boundary representations.
    block_runs omits cross-block merging and is a separate compression hint.
    bit_to_gap_blocks counts each predicted GAP residual once when either
    input is physically BIT; it depends on input storage, unlike other fields.
    GAP predictions always use default capacities, not either input's custom
    levels. This score is a reference-selection hint, not a byte-size promise.
    @ingroup bvcomplexity
    @tparam BV Bit-vector type exposing the BitMagic block manager interface.
    @param a First XOR operand in the common universe [0, bm::id_max).
    @param b Second XOR operand in the same universe as a.
    @return Exact XOR population/run counts and predicted residual block profile,
            including the input-storage-dependent BIT-to-GAP count.
    @sa calc_complexity, bvector_complexity_statistics
*/
template<class BV>
bvector_complexity_statistics calc_complexity_xor(const BV& a, const BV& b) BMNOEXCEPT
{
    return complexity_detail::scan(a, &b);
}

} // namespace bm
#endif
