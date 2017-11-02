#ifndef BMAVX2__H__INCLUDED__
#define BMAVX2__H__INCLUDED__
/*
Copyright(c) 2017 Anatoliy Kuznetsov(anatoliy_kuznetsov at yahoo.com)

Permission is hereby granted, free of charge, to any person
obtaining a copy of this software and associated documentation
files (the "Software"), to deal in the Software without restriction,
including without limitation the rights to use, copy, modify, merge,
publish, distribute, sublicense, and/or sell copies of the Software,
and to permit persons to whom the Software is furnished to do so,
subject to the following conditions:

The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
OTHER DEALINGS IN THE SOFTWARE.

You have to explicitly mention BitMagic project in any derivative product,
its WEB Site, published materials, articles or any other work derived from this
project or based on our code or know-how.

For more information please visit:  http://bitmagic.io

*/

// code here is based on modified libpopcnt library by Kim Walisch
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




// Header implements processor specific intrinsics declarations for AVX2
// instruction set
#include<emmintrin.h>
#include<immintrin.h>

#include "bmdef.h"
#include "bmavx2_util.h"


namespace bm
{


inline void CSA256(__m256i* h, __m256i* l, __m256i a, __m256i b, __m256i c)
{
  __m256i u = _mm256_xor_si256(a, b);
  *h = _mm256_or_si256(_mm256_and_si256(a, b), _mm256_and_si256(u, c));
  *l = _mm256_xor_si256(u, c);
}

inline __m256i avx2_bit_count(__m256i v)
{
  __m256i lookup1 = _mm256_setr_epi8(
      4, 5, 5, 6, 5, 6, 6, 7,
      5, 6, 6, 7, 6, 7, 7, 8,
      4, 5, 5, 6, 5, 6, 6, 7,
      5, 6, 6, 7, 6, 7, 7, 8
  );

  __m256i lookup2 = _mm256_setr_epi8(
      4, 3, 3, 2, 3, 2, 2, 1,
      3, 2, 2, 1, 2, 1, 1, 0,
      4, 3, 3, 2, 3, 2, 2, 1,
      3, 2, 2, 1, 2, 1, 1, 0
  );

  __m256i low_mask = _mm256_set1_epi8(0x0f);
  __m256i lo = _mm256_and_si256(v, low_mask);
  __m256i hi = _mm256_and_si256(_mm256_srli_epi16(v, 4), low_mask);
  __m256i cnt1 = _mm256_shuffle_epi8(lookup1, lo);
  __m256i cnt2 = _mm256_shuffle_epi8(lookup2, hi);

  return _mm256_sad_epu8(cnt1, cnt2);
}

/*
 * AVX2 Harley-Seal popcount (4th iteration).
 * The algorithm is based on the paper "Faster Population Counts
 * using AVX2 Instructions" by Daniel Lemire, Nathan Kurz and
 * Wojciech Mula (23 Nov 2016).
 * @see https://arxiv.org/abs/1611.07612
 */
inline bm::id_t avx2_bit_count(const __m256i* block, const __m256i* block_end)
{
  unsigned i = 0;
  unsigned size = block_end - block;
  uint64_t* cnt64;

  __m256i cnt = _mm256_setzero_si256();
  __m256i ones = _mm256_setzero_si256();
  __m256i twos = _mm256_setzero_si256();
  __m256i fours = _mm256_setzero_si256();
  __m256i eights = _mm256_setzero_si256();
  __m256i sixteens = _mm256_setzero_si256();
  __m256i twosA, twosB, foursA, foursB, eightsA, eightsB;

  for(; i < size; i += 16)
  {
        CSA256(&twosA, &ones, ones, block[i+0], block[i+1]);
        CSA256(&twosB, &ones, ones, block[i+2], block[i+3]);
        CSA256(&foursA, &twos, twos, twosA, twosB);
        CSA256(&twosA, &ones, ones, block[i+4], block[i+5]);
        CSA256(&twosB, &ones, ones, block[i+6], block[i+7]);
        CSA256(&foursB, &twos, twos, twosA, twosB);
        CSA256(&eightsA, &fours, fours, foursA, foursB);
        CSA256(&twosA, &ones, ones, block[i+8], block[i+9]);
        CSA256(&twosB, &ones, ones, block[i+10], block[i+11]);
        CSA256(&foursA, &twos, twos, twosA, twosB);
        CSA256(&twosA, &ones, ones, block[i+12], block[i+13]);
        CSA256(&twosB, &ones, ones, block[i+14], block[i+15]);
        CSA256(&foursB, &twos, twos, twosA, twosB);
        CSA256(&eightsB, &fours, fours, foursA, foursB);
        CSA256(&sixteens, &eights, eights, eightsA, eightsB);

        cnt = _mm256_add_epi64(cnt, avx2_bit_count(sixteens));
  }
  cnt = _mm256_slli_epi64(cnt, 4);
  cnt = _mm256_add_epi64(cnt, _mm256_slli_epi64(avx2_bit_count(eights), 3));
  cnt = _mm256_add_epi64(cnt, _mm256_slli_epi64(avx2_bit_count(fours), 2));
  cnt = _mm256_add_epi64(cnt, _mm256_slli_epi64(avx2_bit_count(twos), 1));
  cnt = _mm256_add_epi64(cnt, avx2_bit_count(ones));

  cnt64 = (uint64_t*) &cnt;

  return cnt64[0] + cnt64[1] + cnt64[2] + cnt64[3];
}


template<class Func>
bm::id_t avx2_bit_count_op(const __m256i* BMRESTRICT block,
                           const __m256i* BMRESTRICT block_end,
                           const __m256i* BMRESTRICT mask_block,
                           Func avx2_func)
{
//    bm::id_t count = 0;
//    id64_t cnt2 = 0;
    uint64_t* cnt64;
    __m256i cnt = _mm256_setzero_si256();
    do
    {
        __m256i tmp0 = _mm256_load_si256(block);
        __m256i tmp1 = _mm256_load_si256(mask_block);

        __m256i b = avx2_func(tmp0, tmp1);
        
        cnt = _mm256_add_epi64(cnt, avx2_bit_count(b));
/*
        count += (unsigned)_mm_popcnt_u64(_mm256_extract_epi64(b, 0));
        count += (unsigned)_mm_popcnt_u64(_mm256_extract_epi64(b, 1));
        count += (unsigned)_mm_popcnt_u64(_mm256_extract_epi64(b, 2));
        count += (unsigned)_mm_popcnt_u64(_mm256_extract_epi64(b, 3));
*/
        ++block; ++mask_block;

    } while (block < block_end);
    
    cnt64 = (uint64_t*) &cnt;
    return cnt64[0] + cnt64[1] + cnt64[2] + cnt64[3];
/*
    cnt2 = cnt64[0] + cnt64[1] + cnt64[2] + cnt64[3];

    if (count != cnt2)
    {
        printf("COUNTING ERR %i %i\n", (int)cnt2, count);
        exit(1);
    }

    return count;
*/
}




#define VECT_XOR_ARR_2_MASK(dst, src, src_end, mask)\
    avx2_xor_arr_2_mask((__m256i*)(dst), (__m256i*)(src), (__m256i*)(src_end), (bm::word_t)mask)

#define VECT_ANDNOT_ARR_2_MASK(dst, src, src_end, mask)\
    avx2_andnot_arr_2_mask((__m256i*)(dst), (__m256i*)(src), (__m256i*)(src_end), (bm::word_t)mask)

#define VECT_BITCOUNT(first, last) \
    avx2_bit_count((__m256i*) (first), (__m256i*) (last))

#define VECT_BITCOUNT_AND(first, last, mask) \
    avx2_bit_count_op((__m256i*) (first), (__m256i*) (last), (__m256i*) (mask), avx2_and)

#define VECT_BITCOUNT_OR(first, last, mask) \
    avx2_bit_count_op((__m256i*) (first), (__m256i*) (last), (__m256i*) (mask), avx2_or)

#define VECT_BITCOUNT_XOR(first, last, mask) \
    avx2_bit_count_op((__m256i*) (first), (__m256i*) (last), (__m256i*) (mask), avx2_xor)

#define VECT_BITCOUNT_SUB(first, last, mask) \
    avx2_bit_count_op((__m256i*) (first), (__m256i*) (last), (__m256i*) (mask), avx2_sub)

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


// TODO: write correct AVX implementation
/*

inline
bm::id_t sse2_bit_block_calc_count_change(const __m128i* BMRESTRICT block,
                                          const __m128i* BMRESTRICT block_end,
                                               unsigned* BMRESTRICT bit_count)
{
   const unsigned mu1 = 0x55555555;
   const unsigned mu2 = 0x33333333;
   const unsigned mu3 = 0x0F0F0F0F;
   const unsigned mu4 = 0x0000003F;

   // Loading masks
   __m128i m1 = _mm_set_epi32 (mu1, mu1, mu1, mu1);
   __m128i m2 = _mm_set_epi32 (mu2, mu2, mu2, mu2);
   __m128i m3 = _mm_set_epi32 (mu3, mu3, mu3, mu3);
   __m128i m4 = _mm_set_epi32 (mu4, mu4, mu4, mu4);
   __m128i mcnt;//, ccnt;
   mcnt = _mm_xor_si128(m1, m1); // bit_cnt = 0
   //ccnt = _mm_xor_si128(m1, m1); // change_cnt = 0

   __m128i tmp1, tmp2;

   int count = (block_end - block)*4; //0;//1;

   bm::word_t  w, w0, w_prev;//, w_l;
   const int w_shift = sizeof(w) * 8 - 1;
   bool first_word = true;
 
   // first word
   {
       const bm::word_t* blk = (const bm::word_t*) block;
       w = w0 = blk[0];
       w ^= (w >> 1);
       BM_INCWORD_BITCOUNT(count, w);
       count -= (w_prev = (w0 >> w_shift)); // negative value correction
   }

   bm::id_t BM_VECT_ALIGN tcnt[4] BM_VECT_ALIGN_ATTR;

   do
   {
       // compute bit-count
       // ---------------------------------------------------------------------
       {
       __m128i b = _mm_load_si128(block);

       // w ^(w >> 1)
       tmp1 = _mm_srli_epi32(b, 1);       // tmp1 = b >> 1
       tmp2 = _mm_xor_si128(b, tmp1);     // tmp2 = tmp1 ^ b;
       _mm_store_si128((__m128i*)tcnt, tmp2);
       

       // compare with zero
       // SSE4: _mm_test_all_zero()
       {
           // b = (b & 0x55555555) + (b >> 1 & 0x55555555);
           //tmp1 = _mm_srli_epi32(b, 1);                    // tmp1 = (b >> 1 & 0x55555555)
           tmp1 = _mm_and_si128(tmp1, m1);
           tmp2 = _mm_and_si128(b, m1);                    // tmp2 = (b & 0x55555555)
           b    = _mm_add_epi32(tmp1, tmp2);               //  b = tmp1 + tmp2

           // b = (b & 0x33333333) + (b >> 2 & 0x33333333);
           tmp1 = _mm_srli_epi32(b, 2);                    // (b >> 2 & 0x33333333)
           tmp1 = _mm_and_si128(tmp1, m2);
           tmp2 = _mm_and_si128(b, m2);                    // (b & 0x33333333)
           b    = _mm_add_epi32(tmp1, tmp2);               // b = tmp1 + tmp2

           // b = (b + (b >> 4)) & 0x0F0F0F0F;
           tmp1 = _mm_srli_epi32(b, 4);                    // tmp1 = b >> 4
           b = _mm_add_epi32(b, tmp1);                     // b = b + (b >> 4)
           b = _mm_and_si128(b, m3);                       //& 0x0F0F0F0F

           // b = b + (b >> 8);
           tmp1 = _mm_srli_epi32 (b, 8);                   // tmp1 = b >> 8
           b = _mm_add_epi32(b, tmp1);                     // b = b + (b >> 8)

           // b = (b + (b >> 16)) & 0x0000003F;
           tmp1 = _mm_srli_epi32 (b, 16);                  // b >> 16
           b = _mm_add_epi32(b, tmp1);                     // b + (b >> 16)
           b = _mm_and_si128(b, m4);                       // (b >> 16) & 0x0000003F;

           mcnt = _mm_add_epi32(mcnt, b);                  // mcnt += b
       }

       }
       // ---------------------------------------------------------------------
       {
           //__m128i b = _mm_load_si128(block);
           // TODO: SSE4...
           //w = _mm_extract_epi32(b, i);               

           const bm::word_t* BMRESTRICT blk = (const bm::word_t*) block;

           if (first_word)
           {
               first_word = false;
           }
           else
           {
               if ((w0=blk[0]))
               {
                   BM_INCWORD_BITCOUNT(count, tcnt[0]);
                   count -= !(w_prev ^ (w0 & 1));
                   count -= w_prev = (w0 >> w_shift);
               }
               else
               {
                   count -= !w_prev; w_prev ^= w_prev;
               }  
           }
           if ((w0=blk[1]))
           {
               BM_INCWORD_BITCOUNT(count, tcnt[1]);
               count -= !(w_prev ^ (w0 & 1));
               count -= w_prev = (w0 >> w_shift);                    
           }
           else
           {
               count -= !w_prev; w_prev ^= w_prev;
           }    
           if ((w0=blk[2]))
           {
               BM_INCWORD_BITCOUNT(count, tcnt[2]);
               count -= !(w_prev ^ (w0 & 1));
               count -= w_prev = (w0 >> w_shift);                    
           }
           else
           {
               count -= !w_prev; w_prev ^= w_prev;
           }      
           if ((w0=blk[3]))
           {
               BM_INCWORD_BITCOUNT(count, tcnt[3]);
               count -= !(w_prev ^ (w0 & 1));
               count -= w_prev = (w0 >> w_shift);                    
           }
           else
           {
               count -= !w_prev; w_prev ^= w_prev;
           }               
       }
   } while (++block < block_end);

   _mm_store_si128((__m128i*)tcnt, mcnt);
   *bit_count = tcnt[0] + tcnt[1] + tcnt[2] + tcnt[3];

   return count;
}

*/




} // namespace




#endif
