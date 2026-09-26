#ifndef BMNEON__H__INCLUDED__
#define BMNEON__H__INCLUDED__
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

/** @file bmneon.h
    Native AArch64 Advanced SIMD kernels. No SSE types or translation layer.
    Pointers address 32-bit words; vector loads need no 16-byte alignment.
    Whole-block operations consume set_block_size words, digest operations
    consume set_block_digest_wave_size words. Range endpoints are exclusive.
*/
#include <arm_neon.h>
#include <cstring>
#include "bmdef.h"
#include "bmconst.h"
#include "bmutil.h"

namespace bm
{
/** @defgroup NEON Native AArch64 NEON kernels
    Native Advanced SIMD implementations of BitMagic block, digest, GAP,
    and bit-count operations.
    @internal
    @ingroup bvector
*/
/**
    @brief Reduce the four 32-bit lanes with bitwise OR.
    @param v Input four-lane vector.
    @return Bitwise OR of the four 32-bit lanes.
    @ingroup NEON
*/
inline
unsigned neon_or_reduce(uint32x4_t v) BMNOEXCEPT
{
    uint32x2_t w = vorr_u32(vget_low_u32(v), vget_high_u32(v));
    return vget_lane_u32(w, 0) | vget_lane_u32(w, 1);
}

// Per-word counts for transition kernels. Widen before accumulating: byte
// accumulators would overflow after 31 all-one vectors.
/**
    @brief Count set bits independently in each 32-bit lane.
    @param v Input four-lane vector.
    @return Four lane counts, one per input lane.
    @ingroup NEON
*/
inline
uint32x4_t neon_count_lanes(uint32x4_t v) BMNOEXCEPT
{
    return vpaddlq_u16(vpaddlq_u8(vcntq_u8(vreinterpretq_u8_u32(v))));
}

/**
    @brief Bitwise operations available to NEON block and count helpers.
    @ingroup NEON
*/
enum neon_bit_op
{
    neon_identity, neon_and, neon_or, neon_xor, neon_sub
};

/**
    @brief Apply the selected bitwise operation to two vectors.
    @tparam Op Bitwise operation to apply.
    @param a First input vector or block.
    @param b Second input vector or block.
    @return Result of the operation selected by Op.
    @ingroup NEON
*/
template<neon_bit_op Op>
inline
uint32x4_t neon_apply(uint32x4_t a, uint32x4_t b) BMNOEXCEPT
{
    if constexpr (Op == neon_and)
        return vandq_u32(a, b);
    else if constexpr (Op == neon_or)
        return vorrq_u32(a, b);
    else if constexpr (Op == neon_xor)
        return veorq_u32(a, b);
    else if constexpr (Op == neon_sub)
        return vbicq_u32(a, b);
    else
        return a;
}

/**
    @brief Count set bits in a word range, optionally combining each word with a mask range.
    @tparam Op Bitwise operation applied before counting.
    @param first Inclusive first word in the input range.
    @param last Exclusive end of the input range.
    @param mask Optional word range for the operation selected by Op; required unless Op is neon_identity.
    @return Number of set bits in the resulting words.
    @ingroup NEON
*/
template<neon_bit_op Op>
bm::id_t neon_bit_count(const bm::word_t* first,
                        const bm::word_t* last,
                        const bm::word_t* mask = 0) BMNOEXCEPT
{
    bm::id_t count = 0;
    const bm::word_t* p = first;
    // Bound every batch to one block: even summing all eight accumulators
    // leaves each 16-bit lane <= 8192. Widen the final horizontal sum.
    while (last-p >= 4)
    {
        unsigned words = unsigned((last-p > bm::set_block_size)
                                  ? bm::set_block_size : last-p);
        words &= ~3u;
        uint16x8_t acc[8];
        for (unsigned j = 0; j < 8; ++j)
        {
            acc[j] = vdupq_n_u16(0);
        } // for j
        unsigned i = 0;
        for (; i + 32 <= words; i += 32)
        {
            for (unsigned j = 0; j < 8; ++j)
            {
                uint32x4_t v = vld1q_u32(p+i+j*4);
                if constexpr (Op != neon_identity)
                    v = neon_apply<Op>(v, vld1q_u32(mask+i+j*4));
                acc[j] = vpadalq_u8(acc[j], vcntq_u8(vreinterpretq_u8_u32(v)));
            } // for j
        } // for i
        for (; i < words; i += 4)
        {
            uint32x4_t v = vld1q_u32(p+i);
            if constexpr (Op != neon_identity)
                v = neon_apply<Op>(v, vld1q_u32(mask+i));
            acc[0] = vpadalq_u8(acc[0], vcntq_u8(vreinterpretq_u8_u32(v)));
        } // for i
        for (unsigned j = 1; j < 8; ++j)
        {
            acc[0] = vaddq_u16(acc[0], acc[j]);
        } // for j
        count += vaddlvq_u16(acc[0]);
        p += words;
        if constexpr (Op != neon_identity)
            mask += words;
    } // while p
    for (; p < last; ++p)
    {
        unsigned v = *p;
        if constexpr (Op != neon_identity)
        {
            unsigned m = *mask++;
            if constexpr (Op == neon_and)
                v &= m;
            else if constexpr (Op == neon_or)
                v |= m;
            else if constexpr (Op == neon_xor)
                v ^= m;
            else if constexpr (Op == neon_sub)
                v &= ~m;
        }
        count += unsigned(__builtin_popcount(v));
    } // for p
    return count;
}

/**
    @brief Count set bits in the block waves selected by a digest.
    @param block Base address of the full bit block; selected waves are addressed relative to this pointer.
    @param digest Digest bit mask selecting block waves.
    @return Total number of set bits in the selected waves.
    @ingroup NEON
*/
inline
bm::id_t neon_bit_count_digest(const bm::word_t* block,
                               bm::id64_t digest) BMNOEXCEPT
{
    bm::id_t count = 0;
    while (digest)
    {
        unsigned wave = unsigned(__builtin_ctzll(digest));
        const bm::word_t* p = block + wave * bm::set_block_digest_wave_size;
        count += neon_bit_count<neon_identity>(p, p + bm::set_block_digest_wave_size);
        digest &= digest - 1;
    } // while digest
    return count;
}

/**
    @brief AND a full block into dst and report whether the result is nonzero.
    @param dst Destination digest wave; ORed with the intersection src1 AND src2.
    @param src1 First input digest wave.
    @return Nonzero OR reduction of the updated block.
    @ingroup NEON
*/
inline
unsigned neon_and_block(bm::word_t* dst, const bm::word_t* src1) BMNOEXCEPT
{
    static_assert(bm::set_block_size % 16 == 0,
                  "NEON kernel requires complete 16-word batches");
    // 4 independent reduction chains, one per vector in each batch.
    uint32x4_t acc0 = vdupq_n_u32(0);
    uint32x4_t acc1 = vdupq_n_u32(0);
    uint32x4_t acc2 = vdupq_n_u32(0);
    uint32x4_t acc3 = vdupq_n_u32(0);
    for (unsigned i = 0; i < bm::set_block_size; i += 16)
    {
        uint32x4_t v0 = vandq_u32(vld1q_u32(dst + i), vld1q_u32(src1 + i));
        vst1q_u32(dst + i, v0);
        acc0 = vorrq_u32(acc0, v0);

        uint32x4_t v1 = vandq_u32(vld1q_u32(dst + i + 4), vld1q_u32(src1 + i + 4));
        vst1q_u32(dst + i + 4, v1);
        acc1 = vorrq_u32(acc1, v1);

        uint32x4_t v2 = vandq_u32(vld1q_u32(dst + i + 8), vld1q_u32(src1 + i + 8));
        vst1q_u32(dst + i + 8, v2);
        acc2 = vorrq_u32(acc2, v2);

        uint32x4_t v3 = vandq_u32(vld1q_u32(dst + i + 12), vld1q_u32(src1 + i + 12));
        vst1q_u32(dst + i + 12, v3);
        acc3 = vorrq_u32(acc3, v3);
    } // for i
    uint32x4_t acc = vorrq_u32(vorrq_u32(acc0, acc1),
                                vorrq_u32(acc2, acc3));
    return neon_or_reduce(acc);
}

/**
    @brief AND one digest wave into dst.
    @param dst Destination word/block range, updated in place or written by the operation.
    @param src1 First source block or wave.
    @return True when the updated wave is empty.
    @ingroup NEON
*/
inline
bool neon_and_digest(bm::word_t* dst, const bm::word_t* src1) BMNOEXCEPT
{
    static_assert(bm::set_block_digest_wave_size % 16 == 0,
                  "NEON kernel requires complete 16-word batches");
    // 4 independent reduction chains, one per vector in each batch.
    uint32x4_t acc0 = vdupq_n_u32(0);
    uint32x4_t acc1 = vdupq_n_u32(0);
    uint32x4_t acc2 = vdupq_n_u32(0);
    uint32x4_t acc3 = vdupq_n_u32(0);
    for (unsigned i = 0; i < bm::set_block_digest_wave_size; i += 16)
    {
        uint32x4_t v0 = vandq_u32(vld1q_u32(dst + i), vld1q_u32(src1 + i));
        vst1q_u32(dst + i, v0);
        acc0 = vorrq_u32(acc0, v0);

        uint32x4_t v1 = vandq_u32(vld1q_u32(dst + i + 4), vld1q_u32(src1 + i + 4));
        vst1q_u32(dst + i + 4, v1);
        acc1 = vorrq_u32(acc1, v1);

        uint32x4_t v2 = vandq_u32(vld1q_u32(dst + i + 8), vld1q_u32(src1 + i + 8));
        vst1q_u32(dst + i + 8, v2);
        acc2 = vorrq_u32(acc2, v2);

        uint32x4_t v3 = vandq_u32(vld1q_u32(dst + i + 12), vld1q_u32(src1 + i + 12));
        vst1q_u32(dst + i + 12, v3);
        acc3 = vorrq_u32(acc3, v3);
    } // for i
    uint32x4_t acc = vorrq_u32(vorrq_u32(acc0, acc1),
                                vorrq_u32(acc2, acc3));
    return neon_or_reduce(acc) == 0;
}

/**
    @brief AND two input waves and write the result to dst.
    @param dst Destination word/block range, updated in place or written by the operation.
    @param src1 First source block or wave.
    @param src2 Second input digest wave.
    @return True when the result wave is empty.
    @ingroup NEON
*/
inline
bool neon_and_digest_2way(bm::word_t* dst,
                          const bm::word_t* src1, const bm::word_t* src2) BMNOEXCEPT
{
    static_assert(bm::set_block_digest_wave_size % 16 == 0,
                  "NEON kernel requires complete 16-word batches");
    // 4 independent reduction chains, one per vector in each batch.
    uint32x4_t acc0 = vdupq_n_u32(0);
    uint32x4_t acc1 = vdupq_n_u32(0);
    uint32x4_t acc2 = vdupq_n_u32(0);
    uint32x4_t acc3 = vdupq_n_u32(0);
    for (unsigned i = 0; i < bm::set_block_digest_wave_size; i += 16)
    {
        uint32x4_t v0 = vandq_u32(vld1q_u32(src1 + i), vld1q_u32(src2 + i));
        vst1q_u32(dst + i, v0);
        acc0 = vorrq_u32(acc0, v0);

        uint32x4_t v1 = vandq_u32(vld1q_u32(src1 + i + 4), vld1q_u32(src2 + i + 4));
        vst1q_u32(dst + i + 4, v1);
        acc1 = vorrq_u32(acc1, v1);

        uint32x4_t v2 = vandq_u32(vld1q_u32(src1 + i + 8), vld1q_u32(src2 + i + 8));
        vst1q_u32(dst + i + 8, v2);
        acc2 = vorrq_u32(acc2, v2);

        uint32x4_t v3 = vandq_u32(vld1q_u32(src1 + i + 12), vld1q_u32(src2 + i + 12));
        vst1q_u32(dst + i + 12, v3);
        acc3 = vorrq_u32(acc3, v3);
    } // for i
    uint32x4_t acc = vorrq_u32(vorrq_u32(acc0, acc1),
                                vorrq_u32(acc2, acc3));
    return neon_or_reduce(acc) == 0;
}

/**
    @brief AND dst with two input waves in place.
    @param dst Destination word/block range, updated in place or written by the operation.
    @param src1 First source block or wave.
    @param src2 Second source block or wave.
    @return True when the updated wave is empty.
    @ingroup NEON
*/
inline
bool neon_and_digest_3way(bm::word_t* dst,
                          const bm::word_t* src1, const bm::word_t* src2) BMNOEXCEPT
{
    static_assert(bm::set_block_digest_wave_size % 8 == 0,
                  "NEON kernel requires complete 8-word batches");
    // 2 independent reduction chains, one per vector in each batch.
    uint32x4_t acc0 = vdupq_n_u32(0);
    uint32x4_t acc1 = vdupq_n_u32(0);
    for (unsigned i = 0; i < bm::set_block_digest_wave_size; i += 8)
    {
        uint32x4_t v0 = vandq_u32(vandq_u32(vld1q_u32(dst + i), vld1q_u32(src1 + i)), vld1q_u32(src2 + i));
        vst1q_u32(dst + i, v0);
        acc0 = vorrq_u32(acc0, v0);

        uint32x4_t v1 = vandq_u32(vandq_u32(vld1q_u32(dst + i + 4), vld1q_u32(src1 + i + 4)), vld1q_u32(src2 + i + 4));
        vst1q_u32(dst + i + 4, v1);
        acc1 = vorrq_u32(acc1, v1);
    } // for i
    uint32x4_t acc = vorrq_u32(acc0, acc1);
    return neon_or_reduce(acc) == 0;
}

/**
    @brief AND dst with four input waves in place.
    @param dst Destination word/block range, updated in place or written by the operation.
    @param src1 First source block or wave.
    @param src2 Second source block or wave.
    @param src3 Third wave to remove.
    @param src4 Fourth wave to remove.
    @return True when the updated wave is empty.
    @ingroup NEON
*/
inline
bool neon_and_digest_5way(bm::word_t* dst,
    const bm::word_t* src1, const bm::word_t* src2, const bm::word_t* src3, const bm::word_t* src4) BMNOEXCEPT
{
    uint32x4_t acc = vdupq_n_u32(0);
    for (unsigned i = 0; i < bm::set_block_digest_wave_size; i += 4)
    {
        uint32x4_t v = vandq_u32(vandq_u32(vandq_u32(vandq_u32(vld1q_u32(dst + i), vld1q_u32(src1 + i)), vld1q_u32(src2 + i)), vld1q_u32(src3 + i)), vld1q_u32(src4 + i));
        vst1q_u32(dst + i, v);
        acc = vorrq_u32(acc, v);
    } // for i
    return neon_or_reduce(acc) == 0;
}

/**
    @brief OR a full block into dst and test whether every result bit is set.
    @param dst Destination full block; updated with dst OR src1.
    @param src1 Input full block.
    @return True when the updated block is all ones.
    @ingroup NEON
*/
inline
bool neon_or_block(bm::word_t* dst, const bm::word_t* src1) BMNOEXCEPT
{
    static_assert(bm::set_block_size % 16 == 0,
                  "NEON kernel requires complete 16-word batches");
    // 4 independent reduction chains, one per vector in each batch.
    uint32x4_t acc0 = vdupq_n_u32(~0u);
    uint32x4_t acc1 = vdupq_n_u32(~0u);
    uint32x4_t acc2 = vdupq_n_u32(~0u);
    uint32x4_t acc3 = vdupq_n_u32(~0u);
    for (unsigned i = 0; i < bm::set_block_size; i += 16)
    {
        uint32x4_t v0 = vorrq_u32(vld1q_u32(dst + i), vld1q_u32(src1 + i));
        vst1q_u32(dst + i, v0);
        acc0 = vandq_u32(acc0, v0);

        uint32x4_t v1 = vorrq_u32(vld1q_u32(dst + i + 4), vld1q_u32(src1 + i + 4));
        vst1q_u32(dst + i + 4, v1);
        acc1 = vandq_u32(acc1, v1);

        uint32x4_t v2 = vorrq_u32(vld1q_u32(dst + i + 8), vld1q_u32(src1 + i + 8));
        vst1q_u32(dst + i + 8, v2);
        acc2 = vandq_u32(acc2, v2);

        uint32x4_t v3 = vorrq_u32(vld1q_u32(dst + i + 12), vld1q_u32(src1 + i + 12));
        vst1q_u32(dst + i + 12, v3);
        acc3 = vandq_u32(acc3, v3);
    } // for i
    uint32x4_t acc = vandq_u32(vandq_u32(acc0, acc1),
                                vandq_u32(acc2, acc3));
    return vminvq_u32(acc) == ~0u;
}

/**
    @brief OR two input blocks and write the result to dst.
    @param dst Destination word/block range, updated in place or written by the operation.
    @param src1 First source block or wave.
    @param src2 Second source block or wave.
    @return True when the result block is all ones.
    @ingroup NEON
*/
inline
bool neon_or_block_2way(bm::word_t* dst, const bm::word_t* src1, const bm::word_t* src2) BMNOEXCEPT
{
    static_assert(bm::set_block_size % 16 == 0,
                  "NEON kernel requires complete 16-word batches");
    // 4 independent reduction chains, one per vector in each batch.
    uint32x4_t acc0 = vdupq_n_u32(~0u);
    uint32x4_t acc1 = vdupq_n_u32(~0u);
    uint32x4_t acc2 = vdupq_n_u32(~0u);
    uint32x4_t acc3 = vdupq_n_u32(~0u);
    for (unsigned i = 0; i < bm::set_block_size; i += 16)
    {
        uint32x4_t v0 = vorrq_u32(vld1q_u32(src1 + i), vld1q_u32(src2 + i));
        vst1q_u32(dst + i, v0);
        acc0 = vandq_u32(acc0, v0);

        uint32x4_t v1 = vorrq_u32(vld1q_u32(src1 + i + 4), vld1q_u32(src2 + i + 4));
        vst1q_u32(dst + i + 4, v1);
        acc1 = vandq_u32(acc1, v1);

        uint32x4_t v2 = vorrq_u32(vld1q_u32(src1 + i + 8), vld1q_u32(src2 + i + 8));
        vst1q_u32(dst + i + 8, v2);
        acc2 = vandq_u32(acc2, v2);

        uint32x4_t v3 = vorrq_u32(vld1q_u32(src1 + i + 12), vld1q_u32(src2 + i + 12));
        vst1q_u32(dst + i + 12, v3);
        acc3 = vandq_u32(acc3, v3);
    } // for i
    uint32x4_t acc = vandq_u32(vandq_u32(acc0, acc1),
                                vandq_u32(acc2, acc3));
    return vminvq_u32(acc) == ~0u;
}

/**
    @brief OR dst with two input blocks in place.
    @param dst Destination word/block range, updated in place or written by the operation.
    @param src1 First source block or wave.
    @param src2 Second source block or wave.
    @return True when the updated block is all ones.
    @ingroup NEON
*/
inline
bool neon_or_block_3way(bm::word_t* dst, const bm::word_t* src1, const bm::word_t* src2) BMNOEXCEPT
{
    static_assert(bm::set_block_size % 8 == 0,
                  "NEON kernel requires complete 8-word batches");
    // 2 independent reduction chains, one per vector in each batch.
    uint32x4_t acc0 = vdupq_n_u32(~0u);
    uint32x4_t acc1 = vdupq_n_u32(~0u);
    for (unsigned i = 0; i < bm::set_block_size; i += 8)
    {
        uint32x4_t v0 = vorrq_u32(vorrq_u32(vld1q_u32(dst + i), vld1q_u32(src1 + i)), vld1q_u32(src2 + i));
        vst1q_u32(dst + i, v0);
        acc0 = vandq_u32(acc0, v0);

        uint32x4_t v1 = vorrq_u32(vorrq_u32(vld1q_u32(dst + i + 4), vld1q_u32(src1 + i + 4)), vld1q_u32(src2 + i + 4));
        vst1q_u32(dst + i + 4, v1);
        acc1 = vandq_u32(acc1, v1);
    } // for i
    uint32x4_t acc = vandq_u32(acc0, acc1);
    return vminvq_u32(acc) == ~0u;
}

/**
    @brief OR dst with four input blocks in place.
    @param dst Destination word/block range, updated in place or written by the operation.
    @param src1 First source block or wave.
    @param src2 Second source block or wave.
    @param src3 Third source block or wave.
    @param src4 Fourth source block or wave.
    @return True when the updated block is all ones.
    @ingroup NEON
*/
inline
bool neon_or_block_5way(bm::word_t* dst,
        const bm::word_t* src1, const bm::word_t* src2, const bm::word_t* src3, const bm::word_t* src4) BMNOEXCEPT
{
    uint32x4_t acc = vdupq_n_u32(~0u);
    for (unsigned i = 0; i < bm::set_block_size; i += 4)
    {
        uint32x4_t v = vorrq_u32(vorrq_u32(vorrq_u32(vorrq_u32(vld1q_u32(dst + i), vld1q_u32(src1 + i)), vld1q_u32(src2 + i)), vld1q_u32(src3 + i)), vld1q_u32(src4 + i));
        vst1q_u32(dst + i, v);
        acc = vandq_u32(acc, v);
    } // for i
    return vminvq_u32(acc) == ~0u;
}

/**
    @brief XOR a full block into dst.
    @param dst Destination full block; updated with dst XOR src1.
    @param src1 Input full block.
    @return Nonzero OR reduction of the updated block.
    @ingroup NEON
*/
inline
unsigned neon_xor_block(bm::word_t* dst, const bm::word_t* src1) BMNOEXCEPT
{
    static_assert(bm::set_block_size % 16 == 0,
                  "NEON kernel requires complete 16-word batches");
    // 4 independent reduction chains, one per vector in each batch.
    uint32x4_t acc0 = vdupq_n_u32(0);
    uint32x4_t acc1 = vdupq_n_u32(0);
    uint32x4_t acc2 = vdupq_n_u32(0);
    uint32x4_t acc3 = vdupq_n_u32(0);
    for (unsigned i = 0; i < bm::set_block_size; i += 16)
    {
        uint32x4_t v0 = veorq_u32(vld1q_u32(dst + i), vld1q_u32(src1 + i));
        vst1q_u32(dst + i, v0);
        acc0 = vorrq_u32(acc0, v0);

        uint32x4_t v1 = veorq_u32(vld1q_u32(dst + i + 4), vld1q_u32(src1 + i + 4));
        vst1q_u32(dst + i + 4, v1);
        acc1 = vorrq_u32(acc1, v1);

        uint32x4_t v2 = veorq_u32(vld1q_u32(dst + i + 8), vld1q_u32(src1 + i + 8));
        vst1q_u32(dst + i + 8, v2);
        acc2 = vorrq_u32(acc2, v2);

        uint32x4_t v3 = veorq_u32(vld1q_u32(dst + i + 12), vld1q_u32(src1 + i + 12));
        vst1q_u32(dst + i + 12, v3);
        acc3 = vorrq_u32(acc3, v3);
    } // for i
    uint32x4_t acc = vorrq_u32(vorrq_u32(acc0, acc1),
                                vorrq_u32(acc2, acc3));
    return neon_or_reduce(acc);
}

/**
    @brief XOR two input blocks and write the result to dst.
    @param dst Destination word/block range, updated in place or written by the operation.
    @param src1 First source block or wave.
    @param src2 Second source block or wave.
    @return Nonzero OR reduction of the result block.
    @ingroup NEON
*/
inline
unsigned neon_xor_block_2way(bm::word_t* dst, const bm::word_t* src1, const bm::word_t* src2) BMNOEXCEPT
{
    static_assert(bm::set_block_size % 16 == 0,
                  "NEON kernel requires complete 16-word batches");
    // 4 independent reduction chains, one per vector in each batch.
    uint32x4_t acc0 = vdupq_n_u32(0);
    uint32x4_t acc1 = vdupq_n_u32(0);
    uint32x4_t acc2 = vdupq_n_u32(0);
    uint32x4_t acc3 = vdupq_n_u32(0);
    for (unsigned i = 0; i < bm::set_block_size; i += 16)
    {
        uint32x4_t v0 = veorq_u32(vld1q_u32(src1 + i), vld1q_u32(src2 + i));
        vst1q_u32(dst + i, v0);
        acc0 = vorrq_u32(acc0, v0);

        uint32x4_t v1 = veorq_u32(vld1q_u32(src1 + i + 4), vld1q_u32(src2 + i + 4));
        vst1q_u32(dst + i + 4, v1);
        acc1 = vorrq_u32(acc1, v1);

        uint32x4_t v2 = veorq_u32(vld1q_u32(src1 + i + 8), vld1q_u32(src2 + i + 8));
        vst1q_u32(dst + i + 8, v2);
        acc2 = vorrq_u32(acc2, v2);

        uint32x4_t v3 = veorq_u32(vld1q_u32(src1 + i + 12), vld1q_u32(src2 + i + 12));
        vst1q_u32(dst + i + 12, v3);
        acc3 = vorrq_u32(acc3, v3);
    } // for i
    uint32x4_t acc = vorrq_u32(vorrq_u32(acc0, acc1),
                                vorrq_u32(acc2, acc3));
    return neon_or_reduce(acc);
}

/**
    @brief Subtract src1 from dst over a full block (dst AND NOT src1).
    @param dst Destination full block; updated with dst AND NOT src1.
    @param src1 Block to remove from dst.
    @return Nonzero OR reduction of the updated block.
    @ingroup NEON
*/
inline
unsigned neon_sub_block(bm::word_t* dst, const bm::word_t* src1) BMNOEXCEPT
{
    static_assert(bm::set_block_size % 16 == 0,
                  "NEON kernel requires complete 16-word batches");
    // 4 independent reduction chains, one per vector in each batch.
    uint32x4_t acc0 = vdupq_n_u32(0);
    uint32x4_t acc1 = vdupq_n_u32(0);
    uint32x4_t acc2 = vdupq_n_u32(0);
    uint32x4_t acc3 = vdupq_n_u32(0);
    for (unsigned i = 0; i < bm::set_block_size; i += 16)
    {
        uint32x4_t v0 = vbicq_u32(vld1q_u32(dst + i), vld1q_u32(src1 + i));
        vst1q_u32(dst + i, v0);
        acc0 = vorrq_u32(acc0, v0);

        uint32x4_t v1 = vbicq_u32(vld1q_u32(dst + i + 4), vld1q_u32(src1 + i + 4));
        vst1q_u32(dst + i + 4, v1);
        acc1 = vorrq_u32(acc1, v1);

        uint32x4_t v2 = vbicq_u32(vld1q_u32(dst + i + 8), vld1q_u32(src1 + i + 8));
        vst1q_u32(dst + i + 8, v2);
        acc2 = vorrq_u32(acc2, v2);

        uint32x4_t v3 = vbicq_u32(vld1q_u32(dst + i + 12), vld1q_u32(src1 + i + 12));
        vst1q_u32(dst + i + 12, v3);
        acc3 = vorrq_u32(acc3, v3);
    } // for i
    uint32x4_t acc = vorrq_u32(vorrq_u32(acc0, acc1),
                                vorrq_u32(acc2, acc3));
    return neon_or_reduce(acc);
}

/**
    @brief Subtract one digest wave from dst (dst AND NOT src1).
    @param dst Destination word/block range, updated in place or written by the operation.
    @param src1 First source block or wave.
    @return True when the updated wave is empty.
    @ingroup NEON
*/
inline
bool neon_sub_digest(bm::word_t* dst, const bm::word_t* src1) BMNOEXCEPT
{
    static_assert(bm::set_block_digest_wave_size % 16 == 0,
                  "NEON kernel requires complete 16-word batches");
    // 4 independent reduction chains, one per vector in each batch.
    uint32x4_t acc0 = vdupq_n_u32(0);
    uint32x4_t acc1 = vdupq_n_u32(0);
    uint32x4_t acc2 = vdupq_n_u32(0);
    uint32x4_t acc3 = vdupq_n_u32(0);
    for (unsigned i = 0; i < bm::set_block_digest_wave_size; i += 16)
    {
        uint32x4_t v0 = vbicq_u32(vld1q_u32(dst + i), vld1q_u32(src1 + i));
        vst1q_u32(dst + i, v0);
        acc0 = vorrq_u32(acc0, v0);

        uint32x4_t v1 = vbicq_u32(vld1q_u32(dst + i + 4), vld1q_u32(src1 + i + 4));
        vst1q_u32(dst + i + 4, v1);
        acc1 = vorrq_u32(acc1, v1);

        uint32x4_t v2 = vbicq_u32(vld1q_u32(dst + i + 8), vld1q_u32(src1 + i + 8));
        vst1q_u32(dst + i + 8, v2);
        acc2 = vorrq_u32(acc2, v2);

        uint32x4_t v3 = vbicq_u32(vld1q_u32(dst + i + 12), vld1q_u32(src1 + i + 12));
        vst1q_u32(dst + i + 12, v3);
        acc3 = vorrq_u32(acc3, v3);
    } // for i
    uint32x4_t acc = vorrq_u32(vorrq_u32(acc0, acc1),
                                vorrq_u32(acc2, acc3));
    return neon_or_reduce(acc) == 0;
}

/**
    @brief Subtract the second input wave from the first and write to dst.
    @param dst Destination word/block range, updated in place or written by the operation.
    @param src1 First source block or wave.
    @param src2 Second source block or wave.
    @return True when the result wave is empty.
    @ingroup NEON
*/
inline
bool neon_sub_digest_2way(bm::word_t* dst, const bm::word_t* src1, const bm::word_t* src2) BMNOEXCEPT
{
    static_assert(bm::set_block_digest_wave_size % 16 == 0,
                  "NEON kernel requires complete 16-word batches");
    // 4 independent reduction chains, one per vector in each batch.
    uint32x4_t acc0 = vdupq_n_u32(0);
    uint32x4_t acc1 = vdupq_n_u32(0);
    uint32x4_t acc2 = vdupq_n_u32(0);
    uint32x4_t acc3 = vdupq_n_u32(0);
    for (unsigned i = 0; i < bm::set_block_digest_wave_size; i += 16)
    {
        uint32x4_t v0 = vbicq_u32(vld1q_u32(src1 + i), vld1q_u32(src2 + i));
        vst1q_u32(dst + i, v0);
        acc0 = vorrq_u32(acc0, v0);

        uint32x4_t v1 = vbicq_u32(vld1q_u32(src1 + i + 4), vld1q_u32(src2 + i + 4));
        vst1q_u32(dst + i + 4, v1);
        acc1 = vorrq_u32(acc1, v1);

        uint32x4_t v2 = vbicq_u32(vld1q_u32(src1 + i + 8), vld1q_u32(src2 + i + 8));
        vst1q_u32(dst + i + 8, v2);
        acc2 = vorrq_u32(acc2, v2);

        uint32x4_t v3 = vbicq_u32(vld1q_u32(src1 + i + 12), vld1q_u32(src2 + i + 12));
        vst1q_u32(dst + i + 12, v3);
        acc3 = vorrq_u32(acc3, v3);
    } // for i
    uint32x4_t acc = vorrq_u32(vorrq_u32(acc0, acc1),
                                vorrq_u32(acc2, acc3));
    return neon_or_reduce(acc) == 0;
}

/**
    @brief Subtract two input waves from dst in place.
    @param dst Destination word/block range, updated in place or written by the operation.
    @param src1 First source block or wave.
    @param src2 Second source block or wave.
    @return True when the updated wave is empty.
    @ingroup NEON
*/
inline
bool neon_sub_digest_3way(bm::word_t* dst, const bm::word_t* src1, const bm::word_t* src2) BMNOEXCEPT
{
    static_assert(bm::set_block_digest_wave_size % 8 == 0,
                  "NEON kernel requires complete 8-word batches");
    // 2 independent reduction chains, one per vector in each batch.
    uint32x4_t acc0 = vdupq_n_u32(0);
    uint32x4_t acc1 = vdupq_n_u32(0);
    for (unsigned i = 0; i < bm::set_block_digest_wave_size; i += 8)
    {
        uint32x4_t v0 = vbicq_u32(vbicq_u32(vld1q_u32(dst + i), vld1q_u32(src1 + i)), vld1q_u32(src2 + i));
        vst1q_u32(dst + i, v0);
        acc0 = vorrq_u32(acc0, v0);

        uint32x4_t v1 = vbicq_u32(vbicq_u32(vld1q_u32(dst + i + 4), vld1q_u32(src1 + i + 4)), vld1q_u32(src2 + i + 4));
        vst1q_u32(dst + i + 4, v1);
        acc1 = vorrq_u32(acc1, v1);
    } // for i
    uint32x4_t acc = vorrq_u32(acc0, acc1);
    return neon_or_reduce(acc) == 0;
}

/**
    @brief Subtract four input waves from dst in place.
    @param dst Destination word/block range, updated in place or written by the operation.
    @param src1 First source block or wave.
    @param src2 Second source block or wave.
    @param src3 Third source block or wave.
    @param src4 Fourth source block or wave.
    @return True when the updated wave is empty.
    @ingroup NEON
*/
inline
bool neon_sub_digest_5way(bm::word_t* dst,
    const bm::word_t* src1, const bm::word_t* src2, const bm::word_t* src3, const bm::word_t* src4) BMNOEXCEPT
{
    uint32x4_t acc = vdupq_n_u32(0);
    for (unsigned i = 0; i < bm::set_block_digest_wave_size; i += 4)
    {
        uint32x4_t v = vbicq_u32(vbicq_u32(vbicq_u32(vbicq_u32(vld1q_u32(dst + i), vld1q_u32(src1 + i)), vld1q_u32(src2 + i)), vld1q_u32(src3 + i)), vld1q_u32(src4 + i));
        vst1q_u32(dst + i, v);
        acc = vorrq_u32(acc, v);
    } // for i
    return neon_or_reduce(acc) == 0;
}

// Return the emptiness of src1 & src2, not the updated destination.
/**
    @brief OR the intersection of two input waves into dst.
    @param dst Destination word/block range, updated in place or written by the operation.
    @param src1 First source block or wave.
    @param src2 Second source block or wave.
    @return True when the intersection is empty.
    @ingroup NEON
*/
inline
bool neon_and_or_digest_2way(bm::word_t* dst, const bm::word_t* src1,
                                    const bm::word_t* src2) BMNOEXCEPT
{
    static_assert(bm::set_block_digest_wave_size % 8 == 0,
                  "NEON kernel requires complete 8-word batches");
    // 2 independent reduction chains, one per vector in each batch.
    uint32x4_t acc0 = vdupq_n_u32(0);
    uint32x4_t acc1 = vdupq_n_u32(0);
    for (unsigned i = 0; i < bm::set_block_digest_wave_size; i += 8)
    {
        uint32x4_t v0 = vandq_u32(vld1q_u32(src1+i), vld1q_u32(src2+i));
        acc0 = vorrq_u32(acc0, v0);
        vst1q_u32(dst+i, vorrq_u32(vld1q_u32(dst+i), v0));

        uint32x4_t v1 = vandq_u32(vld1q_u32(src1+i + 4), vld1q_u32(src2+i + 4));
        acc1 = vorrq_u32(acc1, v1);
        vst1q_u32(dst+i + 4, vorrq_u32(vld1q_u32(dst+i + 4), v1));
    } // for i
    uint32x4_t acc = vorrq_u32(acc0, acc1);
    return neon_or_reduce(acc) == 0;
}

/**
    @brief Apply a repeated mask to a word range, optionally XORing instead of AND-NOT.
    @tparam Xor Selects XOR when true or mask AND NOT source when false.
    @param dst Output range, with one word written per source word.
    @param src Inclusive first source word.
    @param end Exclusive end of the source range.
    @param mask Word mask applied to each source word.
    @ingroup NEON
*/
template<bool Xor>
inline
void neon_mask(bm::word_t* dst, const bm::word_t* src,
                      const bm::word_t* end, bm::word_t mask) BMNOEXCEPT
{
    uint32x4_t m = vdupq_n_u32(mask);
    for (; end-src >= 4; src += 4, dst += 4)
    {
        uint32x4_t v = vld1q_u32(src);
        vst1q_u32(dst, Xor ? veorq_u32(v, m) : vbicq_u32(m, v));
    } // for src/dst
    for (; src < end; ++src, ++dst)
    {
        *dst = Xor ? (*src ^ mask) : (~*src & mask);
    } // for src/dst
}

/**
    @brief Fill a word range with one repeated value.
    @param dst Destination word/block range, updated in place or written by the operation.
    @param value Value written to each destination word.
    @param size Number of sorted endpoints available.
    @ingroup NEON
*/
inline
void neon_set(bm::word_t* dst, unsigned value, unsigned size) BMNOEXCEPT
{
    uint32x4_t v = vdupq_n_u32(value);
    for (unsigned i = 0; i < size; i += 4)
    {
        vst1q_u32(dst+i, v);
    } // for i
}

// Byte loads also cover unaligned serialized input. Ordinary stores are used
// for stream hooks; no cache-policy promise is part of the VECT interface.
/**
    @brief Copy one full bit block from byte-addressable input.
    @param dst Destination word/block range, updated in place or written by the operation.
    @param src Source data; must contain the range required by the operation.
    @ingroup NEON
*/
inline
void neon_copy(bm::word_t* dst, const void* src) BMNOEXCEPT
{
    const unsigned char* p = static_cast<const unsigned char*>(src);
    for (unsigned i = 0; i < bm::set_block_size; i += 4)
    {
        vst1q_u32(dst+i, vreinterpretq_u32_u8(vld1q_u8(p+i*4)));
    } // for i
}

/**
    @brief Invert every word in a full bit block in place.
    @param block Input or mutable bit block, as indicated by the operation.
    @ingroup NEON
*/
inline
void neon_invert(void* block) BMNOEXCEPT
{
    bm::word_t* dst = static_cast<bm::word_t*>(block);
    for (unsigned i = 0; i < bm::set_block_size; i += 4)
    {
        vst1q_u32(dst+i, vmvnq_u32(vld1q_u32(dst+i)));
    } // for i
}

/**
    @brief Test whether every word in a range is uniformly zero or one.
    @tparam Ones Selects all-ones test when true, all-zeros test when false.
    @param block Input or mutable bit block, as indicated by the operation.
    @param size Number of 32-bit words; must be divisible by four.
    @return True if every word equals the value selected by Ones.
    @ingroup NEON
*/
template<bool Ones>
inline
bool neon_uniform(const void* block, unsigned size) BMNOEXCEPT
{
    const bm::word_t* src = static_cast<const bm::word_t*>(block);
    uint32x4_t acc = vdupq_n_u32(0);
    for (unsigned i = 0; i < size; i += 4)
    {
        uint32x4_t v = vld1q_u32(src+i);
        acc = vorrq_u32(acc, Ones ? vmvnq_u32(v) : v);
    } // for i
    return neon_or_reduce(acc) == 0;
}

/**
    @brief Combine serialized words into dst using OR or AND.
    @tparam Or Selects OR when true and AND when false.
    @param dst Destination word/block range, updated in place or written by the operation.
    @param src Source data; must contain the range required by the operation.
    @param size Number of words, indices, or endpoints to process.
    @return For OR, true if all updated words are all ones; for AND, true if any updated bit remains set.
    @ingroup NEON
*/
template<bool Or>
inline
bool neon_decode_arr(bm::word_t* dst,
                     const unsigned char* src,
                     unsigned size) BMNOEXCEPT
{
    uint32x4_t acc = vdupq_n_u32(Or ? ~0u : 0u);
    unsigned i = 0;
    for (; i + 4 <= size; i += 4)
    {
        uint32x4_t v = vreinterpretq_u32_u8(vld1q_u8(src+i*4));
        v = Or ? vorrq_u32(vld1q_u32(dst+i), v) : vandq_u32(vld1q_u32(dst+i), v);
        vst1q_u32(dst+i, v);
        acc = Or ? vandq_u32(acc, v) : vorrq_u32(acc, v);
    } // for i
    bool result = Or ? vminvq_u32(acc) == ~0u : neon_or_reduce(acc) != 0;
    for (; i < size; ++i)
    {
        unsigned v;
        std::memcpy(&v, src+i*4, sizeof(v));
        if (Or) { dst[i] |= v; result &= dst[i] == ~0u; }
        else { dst[i] &= v; result |= dst[i] != 0; }
    } // for i
    return result;
}

// Count run boundaries by comparing each bit with its predecessor. The first
// bit is its own predecessor: all-zero and all-one blocks both have one run.
/**
    @brief Count bit transitions in a word range and optionally count set bits.
    @tparam Xor XOR the ranges before counting transitions.
    @tparam CountBits Also calculate and store the set-bit count.
    @param block Input or mutable bit block, as indicated by the operation.
    @param other Second range XORed with block when Xor is true; otherwise may be null.
    @param size Number of words in each input range.
    @param gc Output transition/run count.
    @param bc Output set-bit count when CountBits is true; otherwise may be null.
    @ingroup NEON
*/
template<bool Xor, bool CountBits>
inline
void neon_change(const bm::word_t* block, const bm::word_t* other,
                 unsigned size, unsigned* gc, unsigned* bc) BMNOEXCEPT
{
    BM_ASSERT(size);
    unsigned carry = (block[0] ^ (Xor ? other[0] : 0u)) & 1u;
    uint32x4_t previous = vdupq_n_u32(carry);
    uint32x4_t runs = vdupq_n_u32(0), bits = runs;
    unsigned i = 0;
    for (; i + 4 <= size; i += 4)
    {
        uint32x4_t v = vld1q_u32(block+i);
        if (Xor) v = veorq_u32(v, vld1q_u32(other+i));
        uint32x4_t hi = vshrq_n_u32(v, 31);
        uint32x4_t prev = vextq_u32(previous, hi, 3);
        uint32x4_t change = veorq_u32(v, vorrq_u32(vshlq_n_u32(v, 1), prev));
        previous = hi;
        runs = vaddq_u32(runs, neon_count_lanes(change));
        if (CountBits) bits = vaddq_u32(bits, neon_count_lanes(v));
    } // for i
    carry = vgetq_lane_u32(previous, 3);
    unsigned r = 1 + vaddvq_u32(runs), b = vaddvq_u32(bits);
    for (; i < size; ++i)
    {
        unsigned v = block[i] ^ (Xor ? other[i] : 0u);
        r += unsigned(__builtin_popcount(v ^ ((v << 1) | carry)));
        if (CountBits) b += unsigned(__builtin_popcount(v));
        carry = v >> 31;
    } // for i
    *gc = r;
    if (CountBits) *bc = b;
}

/**
    @brief Count runs in a word range.
    @param block Input or mutable bit block, as indicated by the operation.
    @param size Number of words, indices, or endpoints to process.
    @return Number of runs in the range.
    @ingroup NEON
*/
inline
unsigned neon_block_change(const bm::word_t* block, unsigned size) BMNOEXCEPT
{
    unsigned gc;
    neon_change<false, false>(block, 0, size, &gc, 0);
    return gc;
}

// In BM's naming R1 moves bits toward larger indices (word << 1).
/**
    @brief Shift a full block by one bit and return the carry and nonzero status.
    @tparam Right Shift toward larger bit indices when true, otherwise smaller indices.
    @param block Input or mutable bit block, as indicated by the operation.
    @param nonzero Output nonzero status for the shifted block.
    @param carry Incoming carry bit at the low edge for right shifts or high edge for left shifts.
    @return Outgoing carry bit from the opposite edge of the block.
    @ingroup NEON
*/
template<bool Right>
inline
bool neon_shift(bm::word_t* block, unsigned* nonzero, unsigned carry) BMNOEXCEPT
{
    uint32x4_t acc = vdupq_n_u32(0);
    for (unsigned k = 0; k < bm::set_block_size; k += 4)
    {
        unsigned i = Right ? k : bm::set_block_size - 4 - k;
        uint32x4_t v = vld1q_u32(block+i), r;
        if (Right)
        {
            uint32x4_t hi = vshrq_n_u32(v, 31);
            r = vorrq_u32(vshlq_n_u32(v, 1), vextq_u32(vdupq_n_u32(carry), hi, 3));
            carry = vgetq_lane_u32(hi, 3);
        }
        else
        {
            uint32x4_t lo = vshlq_n_u32(v, 31);
            r = vorrq_u32(vshrq_n_u32(v, 1), vextq_u32(lo, vdupq_n_u32(carry << 31), 1));
            carry = vgetq_lane_u32(v, 0) & 1u;
        }
        vst1q_u32(block+i, r);
        acc = vorrq_u32(acc, r);
    } // for k
    *nonzero = neon_or_reduce(acc) != 0;
    return carry != 0;
}

/**
    @brief Shift dst toward larger indices, AND it with mask, and update digest.
    @param block Input or mutable bit block, as indicated by the operation.
    @param carry Incoming carry bit at the low edge of the block.
    @param mask Full block mask applied after shifting.
    @param digest Digest output updated for the resulting block.
    @return True when the resulting block is empty.
    @ingroup NEON
*/
inline
bool neon_shift_r1_and(bm::word_t* block, unsigned carry,
                              const bm::word_t* mask, bm::id64_t* digest) BMNOEXCEPT
{
    bm::id64_t d = *digest, result = 0;
    for (unsigned wave = 0; wave < bm::block_waves; ++wave)
    {
        unsigned off = wave * bm::set_block_digest_wave_size;
        if (!(d & (bm::id64_t(1) << wave)))
        {
            if (carry)
            {
                block[off] = carry & mask[off];
                if (block[off]) result |= bm::id64_t(1) << wave;
                carry = 0;
            }
            continue;
        }
        uint32x4_t acc = vdupq_n_u32(0);
        for (unsigned j = 0; j < bm::set_block_digest_wave_size; j += 4)
        {
            unsigned i = off+j;
            uint32x4_t v = vld1q_u32(block+i), hi = vshrq_n_u32(v, 31);
            uint32x4_t r = vorrq_u32(vshlq_n_u32(v, 1), vextq_u32(vdupq_n_u32(carry), hi, 3));
            carry = vgetq_lane_u32(hi, 3); // carry precedes masking
            r = vandq_u32(r, vld1q_u32(mask+i));
            vst1q_u32(block+i, r);
            acc = vorrq_u32(acc, r);
        } // for j
        if (neon_or_reduce(acc)) result |= bm::id64_t(1) << wave;
    } // for wave
    *digest = result;
    return carry != 0;
}

/**
    @brief Find the first set bit at or after a word offset.
    @param src Full bit block to search.
    @param off Starting word offset.
    @param pos Output absolute bit position when found.
    @return True if a set bit is found.
    @ingroup NEON
*/
inline
bool neon_find_first(const bm::word_t* src, unsigned off, unsigned* pos) BMNOEXCEPT
{
    unsigned i = off;
    for (; i + 4 <= bm::set_block_size; i += 4)
    {
        if (neon_or_reduce(vld1q_u32(src+i))) break;
    } // for i
    for (; i < bm::set_block_size; ++i)
    {
        if (src[i])
        {
            *pos = i*32 + unsigned(__builtin_ctz(src[i]));
            return true;
        }
    } // for i
    return false;
}

/**
    @brief Find the first differing bit in two full blocks.
    @param a First full bit block.
    @param b Second full bit block.
    @param pos Output absolute bit position when a difference is found.
    @return True if a differing bit is found.
    @ingroup NEON
*/
inline bool neon_find_diff(const bm::word_t* a, const bm::word_t* b, unsigned* pos) BMNOEXCEPT
{
    for (unsigned i = 0; i < bm::set_block_size; i += 4)
    {
        uint32x4_t v = veorq_u32(vld1q_u32(a+i), vld1q_u32(b+i));
        if (!neon_or_reduce(v)) continue;
        for (unsigned j = 0; j < 4; ++j)
        {
            if (unsigned w = a[i+j] ^ b[i+j])
            {
                *pos = (i+j)*32 + unsigned(__builtin_ctz(w));
                return true;
            }
        } // for j
    } // for i
    return false;
}

/**
    @brief Copy selected waves or XOR them into dst according to a digest.
    @tparam InPlace When true, leave unselected destination waves unchanged.
    @param dst Destination full bit block.
    @param src Source full bit block.
    @param other Second full block used for selected XOR waves.
    @param digest Digest bits selecting waves to XOR; unselected waves are copied unless InPlace is true.
    @ingroup NEON
*/
template<bool InPlace>
void neon_block_xor(bm::word_t* dst, const bm::word_t* src,
                    const bm::word_t* other, bm::id64_t digest) BMNOEXCEPT
{
    for (unsigned wave = 0; wave < bm::block_waves; ++wave)
    {
        bool use_xor = (digest & (bm::id64_t(1) << wave)) != 0;
        if (InPlace && !use_xor) continue;
        unsigned off = wave * bm::set_block_digest_wave_size;
        for (unsigned j = 0; j < bm::set_block_digest_wave_size; j += 4)
        {
            unsigned i = off+j;
            uint32x4_t v = vld1q_u32(src+i);
            if (use_xor) v = veorq_u32(v, vld1q_u32(other+i));
            vst1q_u32(dst+i, v);
        } // for j
    } // for wave
}

/**
    @brief Find the first value not less than target in an inclusive array range.
    @param arr Ascending array to search.
    @param target Value to locate.
    @param from First array index to inspect.
    @param to Inclusive last array index; the result can be to plus one.
    @return Index of the first value >= target, or to plus one.
    @ingroup NEON
*/
inline
unsigned neon_lower_bound(const unsigned* arr, unsigned target,
                          unsigned from, unsigned to) BMNOEXCEPT
{
    unsigned i = from;
    uint32x4_t t = vdupq_n_u32(target);
    for (; i <= to && to-i >= 3; i += 4)
    {
        if (neon_or_reduce(vcgeq_u32(vld1q_u32(arr+i), t))) break;
    } // for i
    for (; i <= to; ++i)
    {
        if (arr[i] >= target) break;
    } // for i
    return i;
}

/**
    @brief Find the first GAP endpoint not less than target.
    @param arr Sorted GAP endpoint array.
    @param target Value to locate.
    @param size Number of words, indices, or endpoints to process.
    @return Index of the first endpoint >= target, or size.
    @ingroup NEON
*/
inline
unsigned neon_gap_find(const bm::gap_word_t* arr, unsigned target,
                               unsigned size) BMNOEXCEPT
{
    unsigned i = 0;
    uint16x8_t t = vdupq_n_u16(bm::gap_word_t(target));
    for (; i+8 <= size; i += 8)
    {
        if (vmaxvq_u16(vcgeq_u16(vld1q_u16(arr+i), t))) break;
    } // for i
    for (; i < size; ++i)
    {
        if (arr[i] >= target) break;
    } // for i
    return i;
}

/**
    @brief Find a position in a GAP buffer and determine its membership state.
    @param buf Encoded GAP buffer.
    @param pos Bit position to locate.
    @param is_set Output membership state (zero or one).
    @return GAP endpoint index associated with pos.
    @ingroup NEON
*/
inline
unsigned neon_gap_bfind(const bm::gap_word_t* buf, unsigned pos,
                                unsigned* is_set) BMNOEXCEPT
{
    unsigned lo = 1, hi = 1 + (buf[0] >> 3);
    while (hi-lo > 16)
    {
        unsigned mid = lo + (hi-lo)/2;
        if (buf[mid] < pos) lo = mid+1; else hi = mid;
    } // while hi/lo
    unsigned idx = lo + neon_gap_find(buf+lo, pos, hi-lo);
    *is_set = (buf[0] & 1u) ^ ((idx-1) & 1u);
    return idx;
}

/**
    @brief Test whether a bit position is set in an encoded GAP buffer.
    @param buf Encoded GAP buffer.
    @param pos Bit position to test.
    @return One if pos is set, otherwise zero.
    @ingroup NEON
*/
inline
unsigned neon_gap_test(const bm::gap_word_t* buf, unsigned pos) BMNOEXCEPT
{
    unsigned set;
    neon_gap_bfind(buf, pos, &set);
    return set;
}

/**
    @brief Accumulate interval lengths from consecutive GAP endpoint waves.
 
    p points at an end of a one-interval; p[-1] is its start boundary.
    One wave consumes 16 endpoints (eight lengths), matching the SSE hook.

    @param p Pointer to the end endpoint of the first one-interval; p[-1] is its start boundary.
    @param waves Number of 16-endpoint waves to consume.
    @param sum Accumulator incremented by the total length of the one-intervals.
    @return Pointer just past the consumed endpoints.
    @ingroup NEON
*/
inline
const bm::gap_word_t* neon_gap_sum_arr(const bm::gap_word_t* p,
                                             unsigned waves, unsigned* sum) BMNOEXCEPT
{
    uint32x4_t acc = vdupq_n_u32(0);
    for (unsigned i = 0; i < waves; ++i, p += 16)
    {
        uint16x8x2_t v = vld2q_u16(p-1);
        acc = vpadalq_u16(acc, vsubq_u16(v.val[1], v.val[0]));
    } // for i
    *sum += vaddvq_u32(acc);
    return p;
}

/**
    @brief Find the end of the index range belonging to one block number.
    @param idx Global bit indices to set; repeated indices are accepted.
    @param size Number of words, indices, or endpoints to process.
    @param nb Block number to match.
    @param start First index position, inclusive.
    @return First index after the matching block-number run, or size.
    @ingroup NEON
*/
inline
unsigned neon_idx_lookup(const unsigned* idx, unsigned size,
                         unsigned nb, unsigned start) BMNOEXCEPT
{
    unsigned i = start;
    uint32x4_t n = vdupq_n_u32(nb);
    for (; i+4 <= size; i += 4)
    {
        if (vminvq_u32(vceqq_u32(vshrq_n_u32(vld1q_u32(idx+i), bm::set_block_shift), n)) != ~0u)
            break;
    } // for i
    for (; i < size && (idx[i] >> bm::set_block_shift) == nb; ++i) {} // for i
    return i;
}

/**
    @brief Set indexed bits in a bit block.
    @param block Input or mutable bit block, as indicated by the operation.
    @param idx Global bit indices to set.
    @param start First index position, inclusive.
    @param stop End index position, exclusive.
    @ingroup NEON
*/
inline
void neon_set_bits(bm::word_t* block, const unsigned* idx,
                   unsigned start, unsigned stop) BMNOEXCEPT
{
    const uint32x4_t mask = vdupq_n_u32(bm::set_block_mask);
    for (; start+4 <= stop; start += 4)
    {
        uint32x4_t n = vandq_u32(vld1q_u32(idx+start), mask);
        uint32x4_t words = vshrq_n_u32(n, bm::set_word_shift);
        uint32x4_t bits = vshlq_u32(vdupq_n_u32(1),
            vreinterpretq_s32_u32(vandq_u32(n, vdupq_n_u32(bm::set_word_mask))));
        // Sequential scatter preserves repeated indices / shared words.
        unsigned w[4], b[4];
        vst1q_u32(w, words); vst1q_u32(b, bits);
        for (unsigned j = 0; j < 4; ++j)
        {
            block[w[j]] |= b[j];
        } // for j
    } // for start
    for (; start < stop; ++start)
    {
        unsigned n = idx[start] & bm::set_block_mask;
        block[n >> bm::set_word_shift] |= 1u << (n & bm::set_word_mask);
    } // for start
}

/**
    @brief Convert a full bit block into GAP endpoint form.
    @param dest Destination GAP buffer; caller must provide sufficient capacity.
    @param src Full source bit block.
    @param capacity Destination buffer capacity in GAP endpoint words; must be sufficient for all endpoints and the sentinel.
    @return Number of GAP endpoints written, excluding the terminal sentinel.
    @ingroup NEON
*/
inline
unsigned neon_bit_to_gap(bm::gap_word_t* dest, const bm::word_t* src,
                         unsigned capacity) BMNOEXCEPT
{
    // Like the scalar/SSE callers, the caller supplies sufficient capacity.
    unsigned len = 1, carry = src[0] & 1u;
    dest[0] = bm::gap_word_t(carry);
    for (unsigned i = 0; i < bm::set_block_size; i += 4)
    {
        uint32x4_t v = vld1q_u32(src+i), hi = vshrq_n_u32(v, 31);
        uint32x4_t changes = veorq_u32(v,
            vorrq_u32(vshlq_n_u32(v, 1), vextq_u32(vdupq_n_u32(carry), hi, 3)));
        carry = vgetq_lane_u32(hi, 3);
        if (!neon_or_reduce(changes)) continue;
        unsigned words[4]; vst1q_u32(words, changes);
        for (unsigned j = 0; j < 4; ++j)
        {
            for (unsigned w = words[j]; w; w &= w-1)
            {
                BM_ASSERT(len < capacity);
                dest[len++] = bm::gap_word_t((i+j)*32 + unsigned(__builtin_ctz(w)) - 1);
            } // for w
        } // for j
    } // for i
    BM_ASSERT(len < capacity); (void)capacity;
    dest[len] = bm::gap_word_t(bm::gap_max_bits-1);
    dest[0] = bm::gap_word_t((len << 3) | (dest[0] & 1u));
    return len;
}

} // namespace bm

#define VECT_AND_BLOCK(dst, src1) \
    neon_and_block((dst), (src1))

#define VECT_AND_DIGEST(dst, src1) \
    neon_and_digest((dst), (src1))

#define VECT_AND_DIGEST_2WAY(dst, src1, src2) \
    neon_and_digest_2way((dst), (src1), (src2))

#define VECT_AND_DIGEST_3WAY(dst, src1, src2) \
    neon_and_digest_3way((dst), (src1), (src2))

#define VECT_AND_DIGEST_5WAY(dst, src1, src2, src3, src4) \
    neon_and_digest_5way((dst), (src1), (src2), (src3), (src4))

#define VECT_OR_BLOCK(dst, src1) \
    neon_or_block((dst), (src1))

#define VECT_OR_BLOCK_2WAY(dst, src1, src2) \
    neon_or_block_2way((dst), (src1), (src2))

#define VECT_OR_BLOCK_3WAY(dst, src1, src2) \
    neon_or_block_3way((dst), (src1), (src2))

#define VECT_OR_BLOCK_5WAY(dst, src1, src2, src3, src4) \
    neon_or_block_5way((dst), (src1), (src2), (src3), (src4))

#define VECT_XOR_BLOCK(dst, src1) \
    neon_xor_block((dst), (src1))

#define VECT_XOR_BLOCK_2WAY(dst, src1, src2) \
    neon_xor_block_2way((dst), (src1), (src2))

#define VECT_SUB_BLOCK(dst, src1) \
    neon_sub_block((dst), (src1))

#define VECT_SUB_DIGEST(dst, src1) \
    neon_sub_digest((dst), (src1))

#define VECT_SUB_DIGEST_2WAY(dst, src1, src2) \
    neon_sub_digest_2way((dst), (src1), (src2))

#define VECT_SUB_DIGEST_3WAY(dst, src1, src2) \
    neon_sub_digest_3way((dst), (src1), (src2))

#define VECT_SUB_DIGEST_5WAY(dst, src1, src2, src3, src4) \
    neon_sub_digest_5way((dst), (src1), (src2), (src3), (src4))

#define VECT_XOR_ARR_2_MASK(dst, src, end, mask) \
    neon_mask<true>((dst), (src), (end), (mask))

#define VECT_ANDNOT_ARR_2_MASK(dst, src, end, mask) \
    neon_mask<false>((dst), (src), (end), (mask))

#define VECT_BITCOUNT(first, last) \
    neon_bit_count<bm::neon_identity>((first), (last))

#define VECT_BIT_COUNT_DIGEST(src, digest) \
    neon_bit_count_digest((src), (digest))

#define VECT_AND_OR_DIGEST_2WAY(dst, a, b) \
    neon_and_or_digest_2way((dst), (a), (b))

#define VECT_INVERT_BLOCK(dst) \
    neon_invert((dst))

#define VECT_SET_BLOCK(dst, value) \
    neon_set((dst), (value), bm::set_block_size)

#define VECT_BLOCK_SET_DIGEST(dst, value) \
    neon_set((dst), (value), bm::set_block_digest_wave_size)

#define VECT_IS_ZERO_BLOCK(src) \
    neon_uniform<false>((src), bm::set_block_size)

#define VECT_IS_ONE_BLOCK(src) \
    neon_uniform<true>((src), bm::set_block_size)

#define VECT_IS_DIGEST_ZERO(src) \
    neon_uniform<false>((src), bm::set_block_digest_wave_size)

#define VECT_LOWER_BOUND_SCAN_U32(arr, target, from, to) \
    neon_lower_bound((arr), (target), (from), (to))

#define VECT_SHIFT_R1(b, acc, co) \
    neon_shift<true>((b), (acc), (co))

#define VECT_SHIFT_L1(b, acc, co) \
    neon_shift<false>((b), (acc), (co))

#define VECT_SHIFT_R1_AND(b, co, m, d) \
    neon_shift_r1_and((b), (co), (m), (d))

#define VECT_BIT_FIND_FIRST(src, off, pos) \
    neon_find_first((src), (off), (pos))

#define VECT_BIT_FIND_DIFF(a, b, pos) \
    neon_find_diff((a), (b), (pos))

#define VECT_BIT_BLOCK_XOR(dst, src, other, d) \
    neon_block_xor<false>((dst), (src), (other), (d))

#define VECT_BIT_BLOCK_XOR_2WAY(dst, other, d) \
    neon_block_xor<true>((dst), (dst), (other), (d))

#define VECT_GAP_BFIND(buf, pos, is_set) \
    neon_gap_bfind((buf), (pos), (is_set))

#define VECT_GAP_TEST(buf, pos) \
    neon_gap_test((buf), (pos))

#define VECT_ARR_BLOCK_LOOKUP(idx, size, nb, start) \
    neon_idx_lookup((idx), (size), (nb), (start))

#define VECT_SET_BLOCK_BITS(block, idx, start, stop) \
    neon_set_bits((block), (idx), (start), (stop))

#define VECT_BLOCK_CHANGE(block, size) \
    neon_block_change((block), (size))

#define VECT_BLOCK_CHANGE_BC(block, gc, bc) \
    neon_change<false, true>((block), 0, bm::set_block_size, (gc), (bc))

#define VECT_BLOCK_XOR_CHANGE(block, other, size, gc, bc) \
    neon_change<true, true>((block), (other), (size), (gc), (bc))

#define VECT_BIT_TO_GAP(dst, src, len) \
    neon_bit_to_gap((dst), (src), (len))

#define VECT_BITCOUNT_AND(first, last, mask) \
    neon_bit_count<bm::neon_and>((first), (last), (mask))

#define VECT_BITCOUNT_OR(first, last, mask) \
    neon_bit_count<bm::neon_or>((first), (last), (mask))

#define VECT_BITCOUNT_XOR(first, last, mask) \
    neon_bit_count<bm::neon_xor>((first), (last), (mask))

#define VECT_BITCOUNT_SUB(first, last, mask) \
    neon_bit_count<bm::neon_sub>((first), (last), (mask))

#define VECT_COPY_BLOCK(dst, src) \
    neon_copy((dst), (src))

#define VECT_COPY_BLOCK_UNALIGN(dst, src) \
    neon_copy((dst), (src))

#define VECT_STREAM_BLOCK(dst, src) \
    neon_copy((dst), (src))

#define VECT_STREAM_BLOCK_UNALIGN(dst, src) \
    neon_copy((dst), (src))

#endif // BMNEON__H__INCLUDED__
