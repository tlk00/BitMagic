#ifndef BMXORFUNC__H__INCLUDED__
#define BMXORFUNC__H__INCLUDED__
/*
Copyright(c) 2002-2019 Anatoliy Kuznetsov(anatoliy_kuznetsov at yahoo.com)

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

/*! \file bmxor.h
    \brief Functions and utilities for XOR filters (internal)
*/

#include "bmdef.h"
#include "bmutil.h"


namespace bm
{

/*!
    Function calculates number of times when bit value changed
    @internal
*/
inline
unsigned bit_block_xor_change32(const bm::word_t* BMRESTRICT block,
                                const bm::word_t* BMRESTRICT xor_block,
                                unsigned size)
{
    unsigned gap_count = 1;

    bm::word_t  w, w0, w_prev, w_l;
    w = w0 = *block ^ *xor_block;

    const int w_shift = int(sizeof(w) * 8 - 1);
    w ^= (w >> 1);
    BM_INCWORD_BITCOUNT(gap_count, w);
    gap_count -= (w_prev = (w0 >> w_shift)); // negative value correction

    const bm::word_t* block_end = block + size; // bm::set_block_size;
    for (++block, ++xor_block; block < block_end; ++block, ++xor_block)
    {
        w = w0 = *block ^ *xor_block;
        ++gap_count;
        if (!w)
        {
            gap_count -= !w_prev;
            w_prev = 0;
        }
        else
        {
            w ^= (w >> 1);
            BM_INCWORD_BITCOUNT(gap_count, w);

            w_l = w0 & 1;
            gap_count -= (w0 >> w_shift);  // negative value correction
            gap_count -= !(w_prev ^ w_l);  // word border correction

            w_prev = (w0 >> w_shift);
        }
    } // for
    return gap_count;
}


/**
    Structure to compute XOR gap-count profile by sub-block waves

    @ingroup bitfunc
    @internal
*/
struct block_waves_xor_descr
{
    unsigned short sb_change[bm::block_waves];
    unsigned short sb_xor_change[bm::block_waves];
};

/**
    Compute reference (non-XOR) 64-dim complexity descriptor for the
    target block.

    Phase 1 of the XOR filtering process

    @internal
*/
inline
void compute_complexity_descr(
                        const bm::word_t* BMRESTRICT block,
                        block_waves_xor_descr& BMRESTRICT x_descr)
{
    for (unsigned i = 0; i < bm::block_waves; ++i)
    {
        unsigned off = (i * bm::set_block_digest_wave_size);
        const bm::word_t* sub_block = block + off;
        #if defined(VECT_BLOCK_CHANGE)
            unsigned change VECT_BLOCK_CHANGE(sub_block,
                                              bm::set_block_digest_wave_size);
        #else
            unsigned change = bm::bit_block_change32(sub_block,
                                                bm::set_block_digest_wave_size);
        #endif
        x_descr.sb_change[i] = (unsigned short) change;
    } // for i
}


/**
    Compute reference complexity descriptor based on XOR vector.
    Returns the digest of sub-blocks where XOR filtering improved the metric
    (function needs reference to estimate the improvement).

    part of Phase 2 of the XOR filtering process

    @sa compute_sub_block_complexity_descr

    @internal
*/
inline
bm::id64_t compute_xor_complexity_descr(
                        const bm::word_t* BMRESTRICT block,
                        const bm::word_t* BMRESTRICT xor_block,
                        block_waves_xor_descr& BMRESTRICT x_descr)
{
    // TODO: add vectorizations
    //
    bm::id64_t digest = 0;
    for (unsigned i = 0; i < bm::block_waves; ++i)
    {
        unsigned off = (i * bm::set_block_digest_wave_size);
        const bm::word_t* sub_block = block + off;
        const bm::word_t* xor_sub_block = xor_block + off;

        // TODO: SIMD
        unsigned xor_change =
                bm::bit_block_xor_change32(sub_block, xor_sub_block,
                                        bm::set_block_digest_wave_size);

        x_descr.sb_xor_change[i] = (unsigned short)xor_change;
        if (xor_change < x_descr.sb_change[i]) // detected improvement
            digest |= (1ull << i);
    } // for i

    return digest;
}

/**
    Build partial XOR product of 2 bit-blocks using digest mask
*/
inline
void bit_block_xor_product(bm::word_t* BMRESTRICT target_block,
                           const bm::word_t* BMRESTRICT block,
                           const bm::word_t* BMRESTRICT xor_block,
                           bm::id64_t digest)
{
    BM_ASSERT(target_block);
    BM_ASSERT(block);
    BM_ASSERT(xor_block);
    BM_ASSERT(digest);

    for (unsigned i = 0; i < bm::block_waves; ++i)
    {
        unsigned off = (i * bm::set_block_digest_wave_size);
        const bm::word_t* sub_block = block + off;
        bm::word_t* t_sub_block = target_block + off;

        const bm::word_t* sub_block_end = sub_block + bm::set_block_digest_wave_size;

        bm::id64_t mask = (1ull << i);
        if (digest & mask) // XOR filtered sub-block
        {
            const bm::word_t* xor_sub_block = xor_block + off;
            for (; sub_block < sub_block_end; ++sub_block, ++xor_sub_block)
            {
                *t_sub_block = *sub_block ^ *xor_sub_block;
                ++t_sub_block;
            }
        }
        else // copy source
        {
            for (; sub_block < sub_block_end; ++sub_block)
            {
                *t_sub_block = *sub_block;
                ++t_sub_block;
            }
        }
    } // for i
}



} // namespace bm

#endif
