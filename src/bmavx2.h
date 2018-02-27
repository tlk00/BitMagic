#ifndef BMAVX2__H__INCLUDED__
#define BMAVX2__H__INCLUDED__
/*
Copyright(c) 2002-2017 Anatoliy Kuznetsov(anatoliy_kuznetsov at yahoo.com)

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

// some of the algorithms here is based on modified libpopcnt library by Kim Walisch
// https://github.com/kimwalisch/libpopcnt/
//
/*
 * libpopcnt.h - C/C++ library for counting the number of 1 bits (bit
 * population count) in an array as quickly as possible using
 * specialized CPU instructions i.e. POPCNT, AVX2, AVX512, NEON.
 *
 * Copyright (c) 2016 - 2017, Kim Walisch
 * Copyright (c) 2016 - 2017, Wojciech Muła
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this
 *    list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */


/** @defgroup AVX2 AVX2 functions
    Processor specific optimizations for AVX2 instructions (internals)
    @ingroup bvector
    @internal
 */


// Header implements processor specific intrinsics declarations for AVX2
// instruction set
//
#include<emmintrin.h>
#include<immintrin.h>

#include "bmdef.h"


namespace bm
{

#define BM_CSA256(h, l, a, b, c) \
{ \
    __m256i u = _mm256_xor_si256(a, b); \
    h = _mm256_or_si256(_mm256_and_si256(a, b), _mm256_and_si256(u, c)); \
    l = _mm256_xor_si256(u, c); \
}

#define BM_AVX2_BIT_COUNT(ret, v) \
{ \
    __m256i lo = _mm256_and_si256(v, low_mask); \
    __m256i hi = _mm256_and_si256(_mm256_srli_epi16(v, 4), low_mask); \
    __m256i cnt1 = _mm256_shuffle_epi8(lookup1, lo); \
    __m256i cnt2 = _mm256_shuffle_epi8(lookup2, hi); \
    ret = _mm256_sad_epu8(cnt1, cnt2); \
} 

#define BM_AVX2_DECL_LOOKUP1 \
  __m256i lookup1 = _mm256_setr_epi8(4, 5, 5, 6, 5, 6, 6, 7, 5, 6, 6, 7, 6, 7, 7, 8, \
                                     4, 5, 5, 6, 5, 6, 6, 7, 5, 6, 6, 7, 6, 7, 7, 8);
#define BM_AVX2_DECL_LOOKUP2 \
__m256i lookup2 = _mm256_setr_epi8(4, 3, 3, 2, 3, 2, 2, 1, 3, 2, 2, 1, 2, 1, 1, 0, \
                                   4, 3, 3, 2, 3, 2, 2, 1, 3, 2, 2, 1, 2, 1, 1, 0);

#define BM_AVX2_POPCNT_PROLOG \
  BM_AVX2_DECL_LOOKUP1 \
  BM_AVX2_DECL_LOOKUP2 \
  __m256i low_mask = _mm256_set1_epi8(0x0f); \
  __m256i bc;

/*!
    @brief AVX2 Harley-Seal popcount
  The algorithm is based on the paper "Faster Population Counts
  using AVX2 Instructions" by Daniel Lemire, Nathan Kurz and
  Wojciech Mula (23 Nov 2016).
  @see https://arxiv.org/abs/1611.07612

  @ingroup AVX2
*/
inline
bm::id_t avx2_bit_count(const __m256i* BMRESTRICT block,
                        const __m256i* BMRESTRICT block_end)
{
  __m256i cnt = _mm256_setzero_si256();
  __m256i ones = _mm256_setzero_si256();
  __m256i twos = _mm256_setzero_si256();
  __m256i fours = _mm256_setzero_si256();
  __m256i eights = _mm256_setzero_si256();
  __m256i sixteens = _mm256_setzero_si256();
  __m256i twosA, twosB, foursA, foursB, eightsA, eightsB;
  __m256i b, c;

  BM_AVX2_POPCNT_PROLOG
  uint64_t* cnt64;

  do
  {
        b = _mm256_load_si256(block+0); c = _mm256_load_si256(block+1);
        BM_CSA256(twosA, ones, ones, b, c);
        b = _mm256_load_si256(block+2); c = _mm256_load_si256(block+3);
        BM_CSA256(twosB, ones, ones, b, c);
        BM_CSA256(foursA, twos, twos, twosA, twosB);
        b = _mm256_load_si256(block+4); c = _mm256_load_si256(block+5);
        BM_CSA256(twosA, ones, ones, b, c);
        b = _mm256_load_si256(block+6); c = _mm256_load_si256(block+7);
        BM_CSA256(twosB, ones, ones, b, c);
        BM_CSA256(foursB, twos, twos, twosA, twosB);
        BM_CSA256(eightsA, fours, fours, foursA, foursB);
        b = _mm256_load_si256(block+8); c = _mm256_load_si256(block+9);
        BM_CSA256(twosA, ones, ones, b, c);
        b = _mm256_load_si256(block+10); c = _mm256_load_si256(block+11);
        BM_CSA256(twosB, ones, ones, b, c);
        BM_CSA256(foursA, twos, twos, twosA, twosB);
        b = _mm256_load_si256(block+12); c = _mm256_load_si256(block+13);
        BM_CSA256(twosA, ones, ones, b, c);
        b = _mm256_load_si256(block+14); c = _mm256_load_si256(block+15);
        BM_CSA256(twosB, ones, ones, b, c);
        BM_CSA256(foursB, twos, twos, twosA, twosB);
        BM_CSA256(eightsB, fours, fours, foursA, foursB);
        BM_CSA256(sixteens, eights, eights, eightsA, eightsB);
        
        BM_AVX2_BIT_COUNT(bc, sixteens);
        cnt = _mm256_add_epi64(cnt, bc);

        block += 16;
  } while (block < block_end);
  
  cnt = _mm256_slli_epi64(cnt, 4);
  BM_AVX2_BIT_COUNT(bc, eights)
  cnt = _mm256_add_epi64(cnt, _mm256_slli_epi64(bc, 3));
  BM_AVX2_BIT_COUNT(bc, fours);
  cnt = _mm256_add_epi64(cnt, _mm256_slli_epi64(bc, 2));
  BM_AVX2_BIT_COUNT(bc, twos); 
  cnt = _mm256_add_epi64(cnt, _mm256_slli_epi64(bc, 1));
  BM_AVX2_BIT_COUNT(bc, ones);
  cnt = _mm256_add_epi64(cnt, bc);

  cnt64 = (uint64_t*) &cnt;

  return (unsigned)(cnt64[0] + cnt64[1] + cnt64[2] + cnt64[3]);
}

/*!
  @brief AND bit count for two aligned bit-blocks
  @ingroup AVX2
*/
inline
bm::id_t avx2_bit_count_and(const __m256i* BMRESTRICT block,
                            const __m256i* BMRESTRICT block_end,
                            const __m256i* BMRESTRICT mask_block)
{
    uint64_t* cnt64;
    BM_AVX2_POPCNT_PROLOG;
    __m256i cnt = _mm256_setzero_si256();
    __m256i ymm0, ymm1;

    
    do
    {
        ymm0 = _mm256_load_si256(block);
        ymm1 = _mm256_load_si256(mask_block);
        ymm0 = _mm256_and_si256(ymm0, ymm1);
        ++block; ++mask_block;
        BM_AVX2_BIT_COUNT(bc, ymm0)
        cnt = _mm256_add_epi64(cnt, bc);

        ymm0 = _mm256_load_si256(block);
        ymm1 = _mm256_load_si256(mask_block);
        ymm0 = _mm256_and_si256(ymm0, ymm1);
        ++block; ++mask_block;
        BM_AVX2_BIT_COUNT(bc, ymm0)
        cnt = _mm256_add_epi64(cnt, bc);

        ymm0 = _mm256_load_si256(block);
        ymm1 = _mm256_load_si256(mask_block);
        ymm0 = _mm256_and_si256(ymm0, ymm1);
        ++block; ++mask_block;
        BM_AVX2_BIT_COUNT(bc, ymm0)
        cnt = _mm256_add_epi64(cnt, bc);

        ymm0 = _mm256_load_si256(block);
        ymm1 = _mm256_load_si256(mask_block);
        ymm0 = _mm256_and_si256(ymm0, ymm1);
        ++block; ++mask_block;
        BM_AVX2_BIT_COUNT(bc, ymm0)
        cnt = _mm256_add_epi64(cnt, bc);

    } while (block < block_end);

    cnt64 = (uint64_t*)&cnt;
    return (unsigned)(cnt64[0] + cnt64[1] + cnt64[2] + cnt64[3]);
}

inline
bm::id_t avx2_bit_count_or(const __m256i* BMRESTRICT block,
    const __m256i* BMRESTRICT block_end,
    const __m256i* BMRESTRICT mask_block)
{
    uint64_t* cnt64;
    BM_AVX2_POPCNT_PROLOG;
    __m256i cnt = _mm256_setzero_si256();
    do
    {
        __m256i tmp0 = _mm256_load_si256(block);
        __m256i tmp1 = _mm256_load_si256(mask_block);

        tmp0 = _mm256_or_si256(tmp0, tmp1);

        BM_AVX2_BIT_COUNT(bc, tmp0)
        cnt = _mm256_add_epi64(cnt, bc);

        ++block; ++mask_block;

    } while (block < block_end);

    cnt64 = (uint64_t*)&cnt;
    return (unsigned)(cnt64[0] + cnt64[1] + cnt64[2] + cnt64[3]);
}
// experimental code for Harley-Seal Hamming
// (needs more testing)
/*
#if !defined(_MSC_VER)
__attribute__((target("avx2")))
#endif
inline bm::id_t avx2_bit_count_xor(const __m256i* BMRESTRICT block,
                                   const __m256i* BMRESTRICT block_end,
                                   const __m256i* BMRESTRICT m_block
                                   )
{
  __m256i cnt = _mm256_setzero_si256();
  __m256i ones = _mm256_setzero_si256();
  __m256i twos = _mm256_setzero_si256();
  __m256i fours = _mm256_setzero_si256();
  __m256i eights = _mm256_setzero_si256();
  __m256i sixteens = _mm256_setzero_si256();
  __m256i twosA, twosB, foursA, foursB, eightsA, eightsB;
  __m256i b, c;
  __m256i m, n;

  BM_AVX2_POPCNT_PROLOG
  uint64_t* cnt64;

  do
  {
        b = _mm256_load_si256(block+0);   c = _mm256_load_si256(block+1);
        m = _mm256_load_si256(m_block+0); n = _mm256_load_si256(m_block+1);
        b = _mm256_xor_si256(b, m);
        c = _mm256_xor_si256(c, n);
        BM_CSA256(twosA, ones, ones, b, c);
      
        b = _mm256_load_si256(block+2); c = _mm256_load_si256(block+3);
        m = _mm256_load_si256(m_block+2); n = _mm256_load_si256(m_block+3);
        b = _mm256_xor_si256(b, m);
        c = _mm256_xor_si256(c, n);
        BM_CSA256(twosB, ones, ones, b, c);
      
        BM_CSA256(foursA, twos, twos, twosA, twosB);
      
        b = _mm256_load_si256(block+4);   c = _mm256_load_si256(block+5);
        m = _mm256_load_si256(m_block+4); n = _mm256_load_si256(m_block+5);
        b = _mm256_xor_si256(b, m);
        c = _mm256_xor_si256(c, n);
        BM_CSA256(twosA, ones, ones, b, c);
      
        b = _mm256_load_si256(block+6);   c = _mm256_load_si256(block+7);
        m = _mm256_load_si256(m_block+6); n = _mm256_load_si256(m_block+7);
        b = _mm256_xor_si256(b, m);
        c = _mm256_xor_si256(c, n);
        BM_CSA256(twosB, ones, ones, b, c);
      
        BM_CSA256(foursB, twos, twos, twosA, twosB);
        BM_CSA256(eightsA, fours, fours, foursA, foursB);
      
        b = _mm256_load_si256(block+8);   c = _mm256_load_si256(block+9);
        m = _mm256_load_si256(m_block+8); n = _mm256_load_si256(m_block+9);
        b = _mm256_xor_si256(b, m);
        c = _mm256_xor_si256(c, n);
        BM_CSA256(twosA, ones, ones, b, c);
      
        b = _mm256_load_si256(block+10);   c = _mm256_load_si256(block+11);
        m = _mm256_load_si256(m_block+10); n = _mm256_load_si256(m_block+11);
        b = _mm256_xor_si256(b, m);
        c = _mm256_xor_si256(c, n);
        BM_CSA256(twosB, ones, ones, b, c);
        BM_CSA256(foursA, twos, twos, twosA, twosB);
      
        b = _mm256_load_si256(block+12);   c = _mm256_load_si256(block+13);
        m = _mm256_load_si256(m_block+12); n = _mm256_load_si256(m_block+13);
        b = _mm256_xor_si256(b, m);
        c = _mm256_xor_si256(c, n);
        BM_CSA256(twosA, ones, ones, b, c);
      
        b = _mm256_load_si256(block+14);   c = _mm256_load_si256(block+15);
        m = _mm256_load_si256(m_block+15); n = _mm256_load_si256(m_block+15);
        b = _mm256_xor_si256(b, m);
        c = _mm256_xor_si256(c, n);
        BM_CSA256(twosB, ones, ones, b, c);
      
        BM_CSA256(foursB, twos, twos, twosA, twosB);
        BM_CSA256(eightsB, fours, fours, foursA, foursB);
        BM_CSA256(sixteens, eights, eights, eightsA, eightsB);
      
        BM_AVX2_BIT_COUNT(bc, sixteens);
        cnt = _mm256_add_epi64(cnt, bc);

        block += 16; m_block += 16;
      
  } while (block < block_end);
  
  cnt = _mm256_slli_epi64(cnt, 4);
  BM_AVX2_BIT_COUNT(bc, eights)
  cnt = _mm256_add_epi64(cnt, _mm256_slli_epi64(bc, 3));
  BM_AVX2_BIT_COUNT(bc, fours);
  cnt = _mm256_add_epi64(cnt, _mm256_slli_epi64(bc, 2));
  BM_AVX2_BIT_COUNT(bc, twos);
  cnt = _mm256_add_epi64(cnt, _mm256_slli_epi64(bc, 1));
  BM_AVX2_BIT_COUNT(bc, ones);
  cnt = _mm256_add_epi64(cnt, bc);

  cnt64 = (uint64_t*) &cnt;

  return cnt64[0] + cnt64[1] + cnt64[2] + cnt64[3];
}
*/

/*!
  @brief XOR bit count for two aligned bit-blocks
  @ingroup AVX2
*/
inline
bm::id_t avx2_bit_count_xor(const __m256i* BMRESTRICT block,
    const __m256i* BMRESTRICT block_end,
    const __m256i* BMRESTRICT mask_block)
{
    uint64_t* cnt64;
    BM_AVX2_POPCNT_PROLOG
    __m256i cnt = _mm256_setzero_si256();
    __m256i ymm0, ymm1;
    do
    {
        ymm0 = _mm256_load_si256(block);
        ymm1 = _mm256_load_si256(mask_block);
        ymm0 = _mm256_xor_si256(ymm0, ymm1);
        ++block; ++mask_block;
        BM_AVX2_BIT_COUNT(bc, ymm0)
        cnt = _mm256_add_epi64(cnt, bc);

        ymm0 = _mm256_load_si256(block);
        ymm1 = _mm256_load_si256(mask_block);
        ymm0 = _mm256_xor_si256(ymm0, ymm1);
        ++block; ++mask_block;
        BM_AVX2_BIT_COUNT(bc, ymm0);
        cnt = _mm256_add_epi64(cnt, bc);

        ymm0 = _mm256_load_si256(block);
        ymm1 = _mm256_load_si256(mask_block);
        ymm0 = _mm256_xor_si256(ymm0, ymm1);
        ++block; ++mask_block;
        BM_AVX2_BIT_COUNT(bc, ymm0);
        cnt = _mm256_add_epi64(cnt, bc);

        ymm0 = _mm256_load_si256(block);
        ymm1 = _mm256_load_si256(mask_block);
        ymm0 = _mm256_xor_si256(ymm0, ymm1);
        ++block; ++mask_block;
        BM_AVX2_BIT_COUNT(bc, ymm0);
        cnt = _mm256_add_epi64(cnt, bc);

    } while (block < block_end);

    cnt64 = (uint64_t*)&cnt;
    return (unsigned)(cnt64[0] + cnt64[1] + cnt64[2] + cnt64[3]);
}



/*!
  @brief AND NOT bit count for two aligned bit-blocks
  @ingroup AVX2
*/
inline
bm::id_t avx2_bit_count_sub(const __m256i* BMRESTRICT block,
    const __m256i* BMRESTRICT block_end,
    const __m256i* BMRESTRICT mask_block)
{
    uint64_t* cnt64;
    BM_AVX2_POPCNT_PROLOG
    __m256i cnt = _mm256_setzero_si256();
    do
    {
        __m256i tmp0 = _mm256_load_si256(block);
        __m256i tmp1 = _mm256_load_si256(mask_block);

        tmp0 = _mm256_andnot_si256(tmp1, tmp0);

        BM_AVX2_BIT_COUNT(bc, tmp0)
        cnt = _mm256_add_epi64(cnt, bc);

        ++block; ++mask_block;

    } while (block < block_end);

    cnt64 = (uint64_t*)&cnt;
    return (unsigned)(cnt64[0] + cnt64[1] + cnt64[2] + cnt64[3]);
}



/*!
    @brief XOR array elements to specified mask
    *dst = *src ^ mask

    @ingroup AVX2
*/
inline
void avx2_xor_arr_2_mask(__m256i* BMRESTRICT dst,
                         const __m256i* BMRESTRICT src,
                         const __m256i* BMRESTRICT src_end,
                         bm::word_t mask)
{
     __m256i ymm2 = _mm256_set_epi32(mask, mask, mask, mask, mask, mask, mask, mask);
     do
     {
        __m256i ymm1 = _mm256_load_si256(src);

        ymm1 = _mm256_xor_si256(ymm1, ymm2);
        _mm256_store_si256(dst, ymm1);
        ++dst;
        ++src;

     } while (src < src_end);
}


/*!
    @brief Inverts array elements and NOT them to specified mask
    *dst = ~*src & mask

    @ingroup AVX2
*/
inline
void avx2_andnot_arr_2_mask(__m256i* BMRESTRICT dst,
                            const __m256i* BMRESTRICT src,
                            const __m256i* BMRESTRICT src_end,
                            bm::word_t mask)
{
     __m256i ymm2 = _mm256_set_epi32(mask, mask, mask, mask, mask, mask, mask, mask);

     do
     {
        __m256i ymm1 = _mm256_load_si256(src);

        ymm1 = _mm256_andnot_si256(ymm1, ymm2); // ymm1 = (~ymm1) & ymm2
        _mm256_store_si256(dst, ymm1);
        ++dst;
        ++src;

     } while (src < src_end);
}

/*!
    @brief AND array elements against another array
    *dst &= *src

    @ingroup AVX2
*/
inline
void avx2_and_arr(__m256i* BMRESTRICT dst,
                  const __m256i* BMRESTRICT src,
                  const __m256i* BMRESTRICT src_end)
{
    __m256i ymm1, ymm2;
    do
    {
        ymm1 = _mm256_load_si256(src++);
        ymm2 = _mm256_load_si256(dst);
        ymm1 = _mm256_and_si256(ymm1, ymm2);
        _mm256_store_si256(dst++, ymm1);
        
        ymm1 = _mm256_load_si256(src++);
        ymm2 = _mm256_load_si256(dst);
        ymm1 = _mm256_and_si256(ymm1, ymm2);
        _mm256_store_si256(dst++, ymm1);

        ymm1 = _mm256_load_si256(src++);
        ymm2 = _mm256_load_si256(dst);
        ymm1 = _mm256_and_si256(ymm1, ymm2);
        _mm256_store_si256(dst++, ymm1);

        ymm1 = _mm256_load_si256(src++);
        ymm2 = _mm256_load_si256(dst);
        ymm1 = _mm256_and_si256(ymm1, ymm2);
        _mm256_store_si256(dst++, ymm1);

    } while (src < src_end);
}


/*!
    @brief OR array elements against another array
    *dst |= *src

    @ingroup AVX2
*/
inline
void avx2_or_arr(__m256i* BMRESTRICT dst,
                 const __m256i* BMRESTRICT src,
                 const __m256i* BMRESTRICT src_end)
{
    __m256i ymm1, ymm2;
    do
    {
        ymm1 = _mm256_load_si256(src++);
        ymm2 = _mm256_load_si256(dst);
        ymm1 = _mm256_or_si256(ymm1, ymm2);
        _mm256_store_si256(dst++, ymm1);
        
        ymm1 = _mm256_load_si256(src++);
        ymm2 = _mm256_load_si256(dst);
        ymm1 = _mm256_or_si256(ymm1, ymm2);
        _mm256_store_si256(dst++, ymm1);

        ymm1 = _mm256_load_si256(src++);
        ymm2 = _mm256_load_si256(dst);
        ymm1 = _mm256_or_si256(ymm1, ymm2);
        _mm256_store_si256(dst++, ymm1);

        ymm1 = _mm256_load_si256(src++);
        ymm2 = _mm256_load_si256(dst);
        ymm1 = _mm256_or_si256(ymm1, ymm2);
        _mm256_store_si256(dst++, ymm1);
    } while (src < src_end);
}


/*!
    @brief OR array elements against another array
    *dst ^= *src

    @ingroup AVX2
*/
inline
void avx2_xor_arr(__m256i* BMRESTRICT dst,
                  const __m256i* BMRESTRICT src,
                  const __m256i* BMRESTRICT src_end)
{
    __m256i ymm1, ymm2;
    do
    {
        ymm1 = _mm256_load_si256(src++);
        ymm2 = _mm256_load_si256(dst);
        ymm1 = _mm256_xor_si256(ymm1, ymm2);
        _mm256_store_si256(dst++, ymm1);

        ymm1 = _mm256_load_si256(src++);
        ymm2 = _mm256_load_si256(dst);
        ymm1 = _mm256_xor_si256(ymm1, ymm2);
        _mm256_store_si256(dst++, ymm1);

        ymm1 = _mm256_load_si256(src++);
        ymm2 = _mm256_load_si256(dst);
        ymm1 = _mm256_xor_si256(ymm1, ymm2);
        _mm256_store_si256(dst++, ymm1);

        ymm1 = _mm256_load_si256(src++);
        ymm2 = _mm256_load_si256(dst);
        ymm1 = _mm256_xor_si256(ymm1, ymm2);
        _mm256_store_si256(dst++, ymm1);
    } while (src < src_end);
}



/*!
    @brief AND-NOT (SUB) array elements against another array
    *dst &= ~*src

    @ingroup AVX2
*/
inline
void avx2_sub_arr(__m256i* BMRESTRICT dst,
                 const __m256i* BMRESTRICT src,
                 const __m256i* BMRESTRICT src_end)
{
    __m256i ymm1, ymm2;
    do
    {
        ymm1 = _mm256_load_si256(src++);
        ymm2 = _mm256_load_si256(dst);
        ymm1 = _mm256_andnot_si256(ymm1, ymm2);
        _mm256_store_si256(dst++, ymm1);
        
        ymm1 = _mm256_load_si256(src++);
        ymm2 = _mm256_load_si256(dst);
        ymm1 = _mm256_andnot_si256(ymm1, ymm2);
        _mm256_store_si256(dst++, ymm1);

        ymm1 = _mm256_load_si256(src++);
        ymm2 = _mm256_load_si256(dst);
        ymm1 = _mm256_andnot_si256(ymm1, ymm2);
        _mm256_store_si256(dst++, ymm1);

        ymm1 = _mm256_load_si256(src++);
        ymm2 = _mm256_load_si256(dst);
        ymm1 = _mm256_andnot_si256(ymm1, ymm2);
        _mm256_store_si256(dst++, ymm1);

    } while (src < src_end);
}

/*!
    @brief AVX2 block memset
    *dst = value

    @ingroup AVX2
*/

BMFORCEINLINE
void avx2_set_block(__m256i* BMRESTRICT dst,
                    __m256i* BMRESTRICT dst_end,
                    bm::word_t value)
{
    __m256i ymm0 = _mm256_set_epi32 (value, value, value, value, value, value, value, value);
    do
    {
        _mm256_store_si256(dst, ymm0);
    } while (++dst < dst_end);
    
    _mm_sfence();
}



/*!
    @brief AVX2 block copy
    *dst = *src

    @ingroup AVX2
*/
inline
void avx2_copy_block(__m256i* BMRESTRICT dst,
                     const __m256i* BMRESTRICT src,
                     const __m256i* BMRESTRICT src_end)
{
    __m256i ymm0, ymm1, ymm2, ymm3;
    do
    {
        ymm0 = _mm256_load_si256(src+0);
        ymm1 = _mm256_load_si256(src+1);
        ymm2 = _mm256_load_si256(src+2);
        ymm3 = _mm256_load_si256(src+3);
        
        _mm256_store_si256(dst+0, ymm0);
        _mm256_store_si256(dst+1, ymm1);
        _mm256_store_si256(dst+2, ymm2);
        _mm256_store_si256(dst+3, ymm3);
        
        ymm0 = _mm256_load_si256(src+4);
        ymm1 = _mm256_load_si256(src+5);
        ymm2 = _mm256_load_si256(src+6);
        ymm3 = _mm256_load_si256(src+7);
        
        _mm256_store_si256(dst+4, ymm0);
        _mm256_store_si256(dst+5, ymm1);
        _mm256_store_si256(dst+6, ymm2);
        _mm256_store_si256(dst+7, ymm3);
        
        src += 8;
        dst += 8;
        
    } while (src < src_end);
}

/*!
    @brief Invert array elements
    *dst = ~*dst
    or
    *dst ^= *dst

    @ingroup AVX2
*/
inline
void avx2_invert_arr(bm::word_t* BMRESTRICT first,
                     bm::word_t* BMRESTRICT last)
{
    __m256i maskz = _mm256_setzero_si256();
    __m256i ymm1 = _mm256_cmpeq_epi64(maskz, maskz); // set all FF

    __m256i* wrd_ptr = (__m256i*)first;
    __m256i ymm0;
    do
    {
        ymm0 = _mm256_load_si256(wrd_ptr+0);
        ymm0 = _mm256_xor_si256(ymm0, ymm1);
        _mm256_store_si256(wrd_ptr+0, ymm0);

        ymm0 = _mm256_load_si256(wrd_ptr+1);
        ymm0 = _mm256_xor_si256(ymm0, ymm1);
        _mm256_store_si256(wrd_ptr+1, ymm0);

        ymm0 = _mm256_load_si256(wrd_ptr+2);
        ymm0 = _mm256_xor_si256(ymm0, ymm1);
        _mm256_store_si256(wrd_ptr+2, ymm0);

        ymm0 = _mm256_load_si256(wrd_ptr+3);
        ymm0 = _mm256_xor_si256(ymm0, ymm1);
        _mm256_store_si256(wrd_ptr+3, ymm0);
        wrd_ptr += 4;

    } while (wrd_ptr < (__m256i*)last);
}

/*!
    @brief check if block is all zero bits
    @ingroup AVX2
*/
inline
bool avx2_is_all_zero(const __m256i* BMRESTRICT block,
                      const __m256i* BMRESTRICT block_end)
{
    __m256i maskz = _mm256_setzero_si256();

    do
    {
        __m256i w0 = _mm256_load_si256(block+0);
        __m256i w1 = _mm256_load_si256(block+1);
        
        __m256i w = _mm256_or_si256(w0, w1);
        __m256i wcmp= _mm256_cmpeq_epi8(w, maskz); // (w0 | w1) == maskz
        unsigned mask = _mm256_movemask_epi8(wcmp);
        if (mask != ~0u)
        {
            return false;
        }

        block += 2;
    
    } while (block < block_end);
    return true;
}

/*!
    @brief check if block is all one bits
    @ingroup AVX2
*/
inline
bool avx2_is_all_one(const __m256i* BMRESTRICT block,
                     const __m256i* BMRESTRICT block_end)
{
    __m256i maskz = _mm256_setzero_si256();
    __m256i maskF = _mm256_cmpeq_epi8(maskz, maskz); // set FF
    
    do
    {
        __m256i w0 = _mm256_load_si256(block);
        
        __m256i wcmp= _mm256_cmpeq_epi8(w0, maskF); // (w0 == maskF)
        unsigned mask = _mm256_movemask_epi8(wcmp);
        if (mask != ~0u)
        {
            return false;
        }
        ++block;
    } while (block < block_end);
    return true;
}



#define VECT_XOR_ARR_2_MASK(dst, src, src_end, mask)\
    avx2_xor_arr_2_mask((__m256i*)(dst), (__m256i*)(src), (__m256i*)(src_end), (bm::word_t)mask)

#define VECT_ANDNOT_ARR_2_MASK(dst, src, src_end, mask)\
    avx2_andnot_arr_2_mask((__m256i*)(dst), (__m256i*)(src), (__m256i*)(src_end), (bm::word_t)mask)

#define VECT_BITCOUNT(first, last) \
    avx2_bit_count((__m256i*) (first), (__m256i*) (last))

#define VECT_BITCOUNT_AND(first, last, mask) \
    avx2_bit_count_and((__m256i*) (first), (__m256i*) (last), (__m256i*) (mask))

#define VECT_BITCOUNT_OR(first, last, mask) \
    avx2_bit_count_or((__m256i*) (first), (__m256i*) (last), (__m256i*) (mask))

#define VECT_BITCOUNT_XOR(first, last, mask) \
    avx2_bit_count_xor((__m256i*) (first), (__m256i*) (last), (__m256i*) (mask))

#define VECT_BITCOUNT_SUB(first, last, mask) \
    avx2_bit_count_sub((__m256i*) (first), (__m256i*) (last), (__m256i*) (mask))

#define VECT_INVERT_ARR(first, last) \
    avx2_invert_arr((bm::word_t*)first, (bm::word_t*)last);

#define VECT_AND_ARR(dst, src, src_end) \
    avx2_and_arr((__m256i*) dst, (__m256i*) (src), (__m256i*) (src_end))

#define VECT_OR_ARR(dst, src, src_end) \
    avx2_or_arr((__m256i*) dst, (__m256i*) (src), (__m256i*) (src_end))

#define VECT_SUB_ARR(dst, src, src_end) \
    avx2_sub_arr((__m256i*) dst, (__m256i*) (src), (__m256i*) (src_end))

#define VECT_XOR_ARR(dst, src, src_end) \
    avx2_xor_arr((__m256i*) dst, (__m256i*) (src), (__m256i*) (src_end))

#define VECT_COPY_BLOCK(dst, src, src_end) \
    avx2_copy_block((__m256i*) dst, (__m256i*) (src), (__m256i*) (src_end))

#define VECT_SET_BLOCK(dst, dst_end, value) \
    avx2_set_block((__m256i*) dst, (__m256i*) (dst_end), (value))

#define VECT_IS_ZERO_BLOCK(dst, dst_end) \
    avx2_is_all_zero((__m256i*) dst, (__m256i*) (dst_end))

#define VECT_IS_ONE_BLOCK(dst, dst_end) \
    avx2_is_all_one((__m256i*) dst, (__m256i*) (dst_end))


// TODO: write better pipelined AVX2 implementation
/*!
    @brief AVX2 optimized bitcounting and number of GAPs
    @ingroup AVX2
*/
inline
bm::id_t avx2_bit_block_calc_count_change(const __m256i* BMRESTRICT block,
                                          const __m256i* BMRESTRICT block_end,
                                               unsigned* BMRESTRICT bit_count)
{
   int count = (unsigned)(block_end - block)*8;

   bm::word_t  w0, w_prev;
   const int w_shift = sizeof(w0) * 8 - 1;
   bool first_word = true;
   *bit_count = 0;
 
   // first word
   {
       bm::word_t  w;
       const bm::word_t* blk = (const bm::word_t*) block;
       w = w0 = blk[0];
       *bit_count += _mm_popcnt_u32(w);
       w ^= (w >> 1);
       count += _mm_popcnt_u32(w);
       count -= (w_prev = (w0 >> w_shift));
   }

   do
   {
       __m256i b = _mm256_load_si256(block);
       __m256i tmp2 = _mm256_xor_si256(b, _mm256_srli_epi32(b, 1)); // tmp2=(b >> 1) ^ b;
       __m256i tmp3 = _mm256_srli_epi32(b, w_shift); // tmp3 = w0 >> w_shift

       // ---------------------------------------------------------------------
       {
           if (first_word)
           {
               first_word = false;
           }
           else
           {
               w0 = _mm256_extract_epi32(b, 0);
               if (w0)
               {
                   *bit_count += _mm_popcnt_u32(w0);
                   count += _mm_popcnt_u32(_mm256_extract_epi32(tmp2, 0));
                   count -= !(w_prev ^ (w0 & 1));
                   count -= w_prev = _mm256_extract_epi32(tmp3, 0);
               }
               else
               {
                   count -= !w_prev; w_prev ^= w_prev;
               }
           }
           w0 = _mm256_extract_epi32(b, 1);
           if (w0)
           {
               *bit_count += _mm_popcnt_u32(w0);
               count += _mm_popcnt_u32(_mm256_extract_epi32(tmp2, 1));
               count -= !(w_prev ^ (w0 & 1));
               count -= w_prev = _mm256_extract_epi32(tmp3, 1);
           }
           else
           {
               count -= !w_prev; w_prev ^= w_prev;
           }
           w0 = _mm256_extract_epi32(b, 2);
           if (w0)
           {
               *bit_count += _mm_popcnt_u32(w0);
               count += _mm_popcnt_u32(_mm256_extract_epi32(tmp2, 2));
               count -= !(w_prev ^ (w0 & 1));
               count -= w_prev = _mm256_extract_epi32(tmp3, 2);
           }
           else
           {
               count -= !w_prev; w_prev ^= w_prev;
           }
           w0 = _mm256_extract_epi32(b, 3);
           if (w0)
           {
               *bit_count += _mm_popcnt_u32(w0);
               count += _mm_popcnt_u32(_mm256_extract_epi32(tmp2, 3));
               count -= !(w_prev ^ (w0 & 1));
               count -= w_prev = _mm256_extract_epi32(tmp3, 3);
           }
           else
           {
               count -= !w_prev; w_prev ^= w_prev;
           }
           w0 = _mm256_extract_epi32(b, 4);
           if (w0)
           {
               *bit_count += _mm_popcnt_u32(w0);
               count += _mm_popcnt_u32(_mm256_extract_epi32(tmp2, 4));
               count -= !(w_prev ^ (w0 & 1));
               count -= w_prev = _mm256_extract_epi32(tmp3, 4);
           }
           else
           {
               count -= !w_prev; w_prev ^= w_prev;
           }
           w0 = _mm256_extract_epi32(b, 5);
           if (w0)
           {
               *bit_count += _mm_popcnt_u32(w0);
               count += _mm_popcnt_u32(_mm256_extract_epi32(tmp2, 5));
               count -= !(w_prev ^ (w0 & 1));
               count -= w_prev = _mm256_extract_epi32(tmp3, 5);
           }
           else
           {
               count -= !w_prev; w_prev ^= w_prev;
           }
           w0 = _mm256_extract_epi32(b, 6);
           if (w0)
           {
               *bit_count += _mm_popcnt_u32(w0);
               count += _mm_popcnt_u32(_mm256_extract_epi32(tmp2, 6));
               count -= !(w_prev ^ (w0 & 1));
               count -= w_prev = _mm256_extract_epi32(tmp3, 6);
           }
           else
           {
               count -= !w_prev; w_prev ^= w_prev;
           }
           w0 = _mm256_extract_epi32(b, 7);
           if (w0)
           {
               *bit_count += _mm_popcnt_u32(w0);
               count += _mm_popcnt_u32(_mm256_extract_epi32(tmp2, 7));
               count -= !(w_prev ^ (w0 & 1));
               count -= w_prev = _mm256_extract_epi32(tmp3, 7);
           }
           else
           {
               count -= !w_prev; w_prev ^= w_prev;
           }
       }
   } while (++block < block_end);

   return count;
}

/*!
    SSE4.2 optimized bitcounting and number of GAPs
    @ingroup SSE4
*/
inline
bm::id_t sse42_bit_block_calc_count_change(const __m128i* BMRESTRICT block,
                                          const __m128i* BMRESTRICT block_end,
                                               unsigned* BMRESTRICT bit_count)
{
//   __m128i mask1 = _mm_set_epi32(0x1, 0x1, 0x1, 0x1);
   BMREGISTER int count = (unsigned)(block_end - block)*4;

   BMREGISTER bm::word_t  w0, w_prev;
   const int w_shift = sizeof(w0) * 8 - 1;
   bool first_word = true;
   *bit_count = 0;
 
   // first word
   {
       bm::word_t  w;
       const bm::word_t* blk = (const bm::word_t*) block;
       w = w0 = blk[0];
       *bit_count += _mm_popcnt_u32(w);
       w ^= (w >> 1);
       count += _mm_popcnt_u32(w);
       count -= (w_prev = (w0 >> w_shift));
   }

   do
   {
       __m128i b = _mm_load_si128(block);
       __m128i tmp2 = _mm_xor_si128(b, _mm_srli_epi32(b, 1)); // tmp2=(b >> 1) ^ b;
       __m128i tmp3 = _mm_srli_epi32(b, w_shift); // tmp3 = w0 >> w_shift
//       __m128i tmp4 = _mm_and_si128(b, mask1);    // tmp4 = w0 & 1

       // ---------------------------------------------------------------------
       {
           if (first_word)
           {
               first_word = false;
           }
           else
           {
               w0 = _mm_extract_epi32(b, 0);
               if (w0)
               {
                   *bit_count += _mm_popcnt_u32(w0);
                   count += _mm_popcnt_u32(_mm_extract_epi32(tmp2, 0));
                   count -= !(w_prev ^ (w0 & 1));
                   count -= w_prev = _mm_extract_epi32(tmp3, 0);
               }
               else
               {
                   count -= !w_prev; w_prev ^= w_prev;
               }
           }
           w0 = _mm_extract_epi32(b, 1);
           if (w0)
           {
               *bit_count += _mm_popcnt_u32(w0);
               count += _mm_popcnt_u32(_mm_extract_epi32(tmp2, 1));
               count -= !(w_prev ^ (w0 & 1));
               count -= w_prev = _mm_extract_epi32(tmp3, 1);
           }
           else
           {
               count -= !w_prev; w_prev ^= w_prev;
           }
           w0 = _mm_extract_epi32(b, 2);
           if (w0)
           {
               *bit_count += _mm_popcnt_u32(w0);
               count += _mm_popcnt_u32(_mm_extract_epi32(tmp2, 2));
               count -= !(w_prev ^ (w0 & 1));
               count -= w_prev = _mm_extract_epi32(tmp3, 2);
           }
           else
           {
               count -= !w_prev; w_prev ^= w_prev;
           }
           w0 = _mm_extract_epi32(b, 3);
           if (w0)
           {
               *bit_count += _mm_popcnt_u32(w0);
               count += _mm_popcnt_u32(_mm_extract_epi32(tmp2, 3));
               count -= !(w_prev ^ (w0 & 1));
               count -= w_prev = _mm_extract_epi32(tmp3, 3);
           }
           else
           {
               count -= !w_prev; w_prev ^= w_prev;
           }
       }
   } while (++block < block_end);

   return count;
}



/* @brief Gap block population count (array sum) utility
   @param pbuf - unrolled, aligned to 1-start GAP buffer
   @param avx_vect_waves - number of AVX vector lines to process
   @param sum - result acumulator
   @return tail pointer

   @internal
*/
inline
const bm::gap_word_t* avx2_gap_sum_arr(const bm::gap_word_t*  pbuf,
                                       unsigned               avx_vect_waves,
                                       unsigned*              sum)
{
    __m256i xcnt = _mm256_setzero_si256();

    // accumulate odd and even elements of the vector the result is 
    // correct based on modulus 16 (max element value in gap blocks is 65535
    // overflow is not an issue here
    for (unsigned i = 0; i < avx_vect_waves; ++i)
    {
        __m256i ymm0 = _mm256_loadu_si256((__m256i*)(pbuf - 1));
        __m256i ymm1 = _mm256_loadu_si256((__m256i*)(pbuf + 16 - 1));
        __m256i ymm_s2 = _mm256_add_epi16(ymm1, ymm0);
        xcnt = _mm256_add_epi16(xcnt, ymm_s2);
        pbuf += 32;
    }
    // odd minus even vector elements clears the result for 1111 blocks
    // bsrli - byte shifts the vector element by 2 bytes (1 short int)
    xcnt = _mm256_sub_epi16(_mm256_bsrli_epi128(xcnt, 2), xcnt);

    // horizontal sum of vector elements
    // cnt16[0] + cnt16[2] + cnt16[4] + cnt16[6] + cnt16[8] + cnt16[10] + cnt16[12] + cnt16[14];
    //
    xcnt = _mm256_add_epi16(_mm256_bsrli_epi128(xcnt, 4), xcnt);
    xcnt = _mm256_add_epi16(_mm256_bsrli_epi128(xcnt, 8), xcnt);
    __m128i xcnt2 = _mm_add_epi16(_mm256_extracti128_si256(xcnt, 1), _mm256_extracti128_si256(xcnt, 0));

    // extract 32-bit word and mask to take first 16 bits
    *sum += _mm_cvtsi128_si32(xcnt2) & 0xffff;
    return pbuf;
}



} // namespace




#endif
