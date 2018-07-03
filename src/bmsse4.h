#ifndef BMSSE4__H__INCLUDED__
#define BMSSE4__H__INCLUDED__
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

/*! \file bmsse4.h
    \brief Compute functions for SSE4.2 SIMD instruction set (internal)
*/

#include<mmintrin.h>
#include<emmintrin.h>
#include<smmintrin.h>
#include<nmmintrin.h>

#include "bmdef.h"
#include "bmsse_util.h"
#include "bmutil.h"

namespace bm
{

/** @defgroup SSE4 SSE4.2 funcions (internal)
    Processor specific optimizations for SSE4.2 instructions (internals)
    @internal
    @ingroup bvector
 */

#ifdef __GNUG__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
#endif



/*!
    SSE4.2 optimized bitcounting .
    @ingroup SSE4
*/
inline 
bm::id_t sse4_bit_count(const __m128i* block, const __m128i* block_end)
{
    bm::id_t count = 0;
#ifdef BM64_SSE4
    const bm::id64_t* b = (bm::id64_t*) block;
    const bm::id64_t* b_end = (bm::id64_t*) block_end;
    do
    {
        count += unsigned( _mm_popcnt_u64(b[0]) +
                           _mm_popcnt_u64(b[1]));
        b += 2;
    } while (b < b_end);
#else
    do
    {
        const unsigned* b = (unsigned*) block;
        count += _mm_popcnt_u32(b[0]) +
                 _mm_popcnt_u32(b[1]) +
                 _mm_popcnt_u32(b[2]) +
                 _mm_popcnt_u32(b[3]);
    } while (++block < block_end);
#endif    
    return count;
}

/*!
\internal
*/
BMFORCEINLINE 
unsigned op_xor(unsigned a, unsigned b)
{
    unsigned ret = (a ^ b);
    return ret;
}

/*!
\internal
*/
BMFORCEINLINE 
unsigned op_or(unsigned a, unsigned b)
{
    return (a | b);
}

/*!
\internal
*/
BMFORCEINLINE 
unsigned op_and(unsigned a, unsigned b)
{
    return (a & b);
}


template<class Func>
bm::id_t sse4_bit_count_op(const __m128i* BMRESTRICT block, 
                           const __m128i* BMRESTRICT block_end,
                           const __m128i* BMRESTRICT mask_block,
                           Func sse2_func)
{
    bm::id_t count = 0;
#ifdef BM64_SSE4
    do
    {
        __m128i tmp0 = _mm_load_si128(block);
        __m128i tmp1 = _mm_load_si128(mask_block);        
        __m128i b = sse2_func(tmp0, tmp1);

        count += (unsigned)_mm_popcnt_u64(_mm_extract_epi64(b, 0));
        count += (unsigned)_mm_popcnt_u64(_mm_extract_epi64(b, 1));

        ++block; ++mask_block;
    } while (block < block_end);
#else    
    do
    {
        __m128i tmp0 = _mm_load_si128(block);
        __m128i tmp1 = _mm_load_si128(mask_block);        
        __m128i b = sse2_func(tmp0, tmp1);

        count += _mm_popcnt_u32(_mm_extract_epi32(b, 0));
        count += _mm_popcnt_u32(_mm_extract_epi32(b, 1));
        count += _mm_popcnt_u32(_mm_extract_epi32(b, 2));
        count += _mm_popcnt_u32(_mm_extract_epi32(b, 3));

        ++block; ++mask_block;
    } while (block < block_end);
#endif
    
    return count;
}

/*!
    @brief check if block is all zero bits
    @ingroup SSE4
*/
inline
bool sse4_is_all_zero(const __m128i* BMRESTRICT block,
                      const __m128i* BMRESTRICT block_end)
{
    __m128i w0, w1, w;
    __m128i maskz = _mm_setzero_si128();

    do
    {
        w0 = _mm_load_si128(block+0);
        w1 = _mm_load_si128(block+1);
        
        w = _mm_or_si128(w0, w1);
        if (!_mm_test_all_ones(_mm_cmpeq_epi8(w, maskz))) // (w0 | w1) != maskz
            return false;
        
        w0 = _mm_load_si128(block+2);
        w1 = _mm_load_si128(block+3);
        
        w = _mm_or_si128(w0, w1);
        if (!_mm_test_all_ones(_mm_cmpeq_epi8(w, maskz))) // (w0 | w1) != maskz
            return false;

        block += 4;
    
    } while (block < block_end);
    return true;
}


/*!
    @brief check if block is all zero bits
    @ingroup SSE4
*/
inline
bool sse4_is_all_one(const __m128i* BMRESTRICT block,
                     const __m128i* BMRESTRICT block_end)
{
    do
    {
        __m128i w0 = _mm_load_si128(block);
        if (!_mm_test_all_ones(w0))
        {
            return false;
        }
        ++block;
    } while (block < block_end);
    return true;
}

/*!
    @brief check if wave of pointers is all NULL
    @ingroup AVX2
*/
BMFORCEINLINE
bool sse42_test_all_zero_wave(void* ptr)
{
    __m128i w0 = _mm_loadu_si128((__m128i*)ptr);
    return _mm_testz_si128(w0, w0);
}


#define VECT_XOR_ARR_2_MASK(dst, src, src_end, mask)\
    sse2_xor_arr_2_mask((__m128i*)(dst), (__m128i*)(src), (__m128i*)(src_end), (bm::word_t)mask)

#define VECT_ANDNOT_ARR_2_MASK(dst, src, src_end, mask)\
    sse2_andnot_arr_2_mask((__m128i*)(dst), (__m128i*)(src), (__m128i*)(src_end), (bm::word_t)mask)

#define VECT_BITCOUNT(first, last) \
    sse4_bit_count((__m128i*) (first), (__m128i*) (last)) 

#define VECT_BITCOUNT_AND(first, last, mask) \
    sse4_bit_count_op((__m128i*) (first), (__m128i*) (last), (__m128i*) (mask), sse2_and) 

#define VECT_BITCOUNT_OR(first, last, mask) \
    sse4_bit_count_op((__m128i*) (first), (__m128i*) (last), (__m128i*) (mask), sse2_or) 

#define VECT_BITCOUNT_XOR(first, last, mask) \
    sse4_bit_count_op((__m128i*) (first), (__m128i*) (last), (__m128i*) (mask), sse2_xor) 

#define VECT_BITCOUNT_SUB(first, last, mask) \
    sse4_bit_count_op((__m128i*) (first), (__m128i*) (last), (__m128i*) (mask), sse2_sub) 

#define VECT_INVERT_ARR(first, last) \
    sse2_invert_arr((bm::word_t*)first, (bm::word_t*)last);

#define VECT_AND_ARR(dst, src, src_end) \
    sse2_and_arr((__m128i*) dst, (__m128i*) (src), (__m128i*) (src_end))

#define VECT_OR_ARR(dst, src, src_end) \
    sse2_or_arr((__m128i*) dst, (__m128i*) (src), (__m128i*) (src_end))

#define VECT_OR_ARR_3WAY(dst, src1, src2, src1_end) \
    sse2_or_arr_3way((__m128i*) (dst), (__m128i*) (src1), (__m128i*) (src2), (__m128i*) (src1_end))

#define VECT_OR_ARR_5WAY(dst, src1, src2, src3, src4, src1_end) \
    sse2_or_arr_5way((__m128i*) (dst), (__m128i*) (src1), (__m128i*) (src2), (__m128i*) (src3), (__m128i*) (src4), (__m128i*) (src1_end))

#define VECT_SUB_ARR(dst, src, src_end) \
    sse2_sub_arr((__m128i*) dst, (__m128i*) (src), (__m128i*) (src_end))

#define VECT_XOR_ARR(dst, src, src_end) \
    sse2_xor_arr((__m128i*) dst, (__m128i*) (src), (__m128i*) (src_end))

#define VECT_COPY_BLOCK(dst, src, src_end) \
    sse2_copy_block((__m128i*) dst, (__m128i*) (src), (__m128i*) (src_end))

#define VECT_SET_BLOCK(dst, dst_end, value) \
    sse2_set_block((__m128i*) dst, (__m128i*) (dst_end), (value))

#define VECT_IS_ZERO_BLOCK(dst, dst_end) \
    sse4_is_all_zero((__m128i*) dst, (__m128i*) (dst_end))

#define VECT_IS_ONE_BLOCK(dst, dst_end) \
    sse4_is_all_one((__m128i*) dst, (__m128i*) (dst_end))



/*!
    SSE4.2 optimized bitcounting and number of GAPs
    @ingroup SSE4
*/
inline
bm::id_t sse4_bit_block_calc_count_change(const __m128i* BMRESTRICT block,
                                          const __m128i* BMRESTRICT block_end,
                                               unsigned* BMRESTRICT bit_count)
{
   int count = (unsigned)(block_end - block)*4;

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



#ifdef __GNUG__
// necessary measure to silence false warning from GCC about negative pointer arithmetics
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Warray-bounds"
#endif

/*!
     SSE4.2 check for one to two (variable len) 128 bit SSE lines for gap search results (8 elements)
     \internal
*/
inline
unsigned sse4_gap_find(const bm::gap_word_t* BMRESTRICT pbuf, const bm::gap_word_t pos, const unsigned size)
{
    BM_ASSERT(size <= 16);
    BM_ASSERT(size);

    const unsigned unroll_factor = 8;
    if (size < 4) // for very short vector use conventional scan
    {
        unsigned j;
        for (j = 0; j < size; ++j)
        {
            if (pbuf[j] >= pos)
                break;
        }
        return j;
    }

    __m128i m1, mz, maskF, maskFL;

    mz = _mm_setzero_si128();
    m1 = _mm_loadu_si128((__m128i*)(pbuf)); // load first 8 elements

    maskF = _mm_cmpeq_epi64(mz, mz); // set all FF
    maskFL = _mm_slli_si128(maskF, 4 * 2); // byle shift to make [0000 FFFF] 
    int shiftL= (64 - (unroll_factor - size) * 16);
    maskFL = _mm_slli_epi64(maskFL, shiftL); // additional bit shift to  [0000 00FF]

    m1 = _mm_andnot_si128(maskFL, m1); // m1 = (~mask) & m1
    m1 = _mm_or_si128(m1, maskFL);

    __m128i mp = _mm_set1_epi16(pos);  // broadcast pos into all elements of a SIMD vector
    __m128i  mge_mask = _mm_cmpeq_epi16(_mm_subs_epu16(mp, m1), mz); // unsigned m1 >= mp
    __m128i  c_mask = _mm_slli_epi16(mge_mask, 15); // clear not needed flag bits by shift
    int mi = _mm_movemask_epi8(c_mask);  // collect flag bits
    if (mi)
    {
        // alternative: int bsr_i= bm::bit_scan_fwd(mi) >> 1;
        unsigned bc = _mm_popcnt_u32(mi); // gives us number of elements >= pos
        return unroll_factor - bc;   // address of first one element (target)
    }
    // inspect the next lane with possible step back (to avoid over-read the block boundaries)
    //   GCC gives a false warning for "- unroll_factor" here
    const bm::gap_word_t* BMRESTRICT pbuf2 = pbuf + size - unroll_factor;
    BM_ASSERT(pbuf2 > pbuf || size == 8); // assert in place to make sure GCC warning is indeed false

    m1 = _mm_loadu_si128((__m128i*)(pbuf2)); // load next elements (with possible overlap)
    mge_mask = _mm_cmpeq_epi16(_mm_subs_epu16(mp, m1), mz); // m1 >= mp        
    mi = _mm_movemask_epi8(_mm_slli_epi16(mge_mask, 15));
    unsigned bc = _mm_popcnt_u32(mi); 

    return size - bc;
}

/*!
     SSE4.2 index lookup to check what belongs to the same block (8 elements)
     \internal
*/
inline
unsigned sse4_idx_arr_block_lookup(const unsigned* idx, unsigned size,
                                   unsigned nb, unsigned start)
{
    const unsigned unroll_factor = 8;
    const unsigned len = (size - start);
    const unsigned len_unr = len - (len % unroll_factor);
    unsigned k;

    idx += start;

    __m128i nbM = _mm_set1_epi32(nb);

    for (k = 0; k < len_unr; k+=unroll_factor)
    {
        __m128i idxA = _mm_loadu_si128((__m128i*)(idx+k));
        __m128i idxB = _mm_loadu_si128((__m128i*)(idx+k+4));
        __m128i nbA =  _mm_srli_epi32(idxA, bm::set_block_shift); // idx[k] >> bm::set_block_shift
        __m128i nbB =  _mm_srli_epi32(idxB, bm::set_block_shift);
        
        if (!_mm_test_all_ones(_mm_cmpeq_epi32(nbM, nbA)) |
            !_mm_test_all_ones(_mm_cmpeq_epi32 (nbM, nbB)))
            break;

    } // for k
    for (; k < len; ++k)
    {
        if (nb != unsigned(idx[k] >> bm::set_block_shift))
            break;
    }
    return start + k;
}

/*!
     SSE4.2 bit block gather-scatter
 
     @param arr - destination array to set bits
     @param blk - source bit-block
     @param idx - gather index array
     @param size - gather array size
     @param start - gaher start index
     @param bit_idx - bit to set in the target array
 
     \internal

    C algorithm:
 
    for (unsigned k = start; k < size; ++k)
    {
        nbit = unsigned(idx[k] & bm::set_block_mask);
        nword  = unsigned(nbit >> bm::set_word_shift);
        mask0 = 1u << (nbit & bm::set_word_mask);
        arr[k] |= TRGW(bool(blk[nword] & mask0) << bit_idx);
    }

*/
inline
void sse4_bit_block_gather_scatter(unsigned* BMRESTRICT arr,
                                   const unsigned* BMRESTRICT blk,
                                   const unsigned* BMRESTRICT idx,
                                   unsigned                   size,
                                   unsigned                   start,
                                   unsigned                   bit_idx)
{
    const unsigned unroll_factor = 4;
    const unsigned len = (size - start);
    const unsigned len_unr = len - (len % unroll_factor);
    
    __m128i sb_mask = _mm_set1_epi32(bm::set_block_mask);
    __m128i sw_mask = _mm_set1_epi32(bm::set_word_mask);
    __m128i maskFF  = _mm_set1_epi32(~0u);
    __m128i maskZ   = _mm_xor_si128(maskFF, maskFF);

    __m128i mask_tmp, mask_0;
    
    unsigned BM_ALIGN16 mshift_v[4] BM_ALIGN16ATTR;
    unsigned BM_ALIGN16 mword_v[4] BM_ALIGN16ATTR;

    unsigned k = 0;
    unsigned base = start + k;
    __m128i* idx_ptr = (__m128i*)(idx + base);   // idx[base]
    __m128i* target_ptr = (__m128i*)(arr + base); // arr[base]
    for (; k < len_unr; k+=unroll_factor)
    {
        __m128i nbitA = _mm_and_si128 (_mm_loadu_si128(idx_ptr), sb_mask); // nbit = idx[base] & bm::set_block_mask
        __m128i nwordA = _mm_srli_epi32 (nbitA, bm::set_word_shift); // nword  = nbit >> bm::set_word_shift
        // (nbit & bm::set_word_mask)
        _mm_store_si128 ((__m128i*)mshift_v, _mm_and_si128 (nbitA, sw_mask));
        _mm_store_si128 ((__m128i*)mword_v, nwordA);
        
        // mask0 = 1u << (nbit & bm::set_word_mask);
        //
#if 0
        // ifdefed an alternative SHIFT implementation using SSE and masks
        // (it is not faster than just doing scalar operations)
        {
        __m128i am_0    = _mm_set_epi32(0, 0, 0, ~0u);
        __m128i mask1   = _mm_srli_epi32 (maskFF, 31);
        mask_0 = _mm_and_si128 (_mm_slli_epi32 (mask1, mshift_v[0]), am_0);
        mask_tmp = _mm_and_si128 (_mm_slli_epi32(mask1, mshift_v[1]), _mm_slli_si128 (am_0, 4));
        mask_0 = _mm_or_si128 (mask_0, mask_tmp);
    
        __m128i mask_2 = _mm_and_si128 (_mm_slli_epi32 (mask1, mshift_v[2]),
                                        _mm_slli_si128 (am_0, 8));
        mask_tmp = _mm_and_si128 (
                      _mm_slli_epi32(mask1, mshift_v[3]),
                      _mm_slli_si128 (am_0, 12)
                      );
    
        mask_0 = _mm_or_si128 (mask_0,
                               _mm_or_si128 (mask_2, mask_tmp)); // assemble bit-test mask
        }
#endif
        mask_0 = _mm_set_epi32(1 << mshift_v[3], 1 << mshift_v[2], 1 << mshift_v[1], 1 << mshift_v[0]);


        //  gather for: blk[nword]  (.. & mask0 )
        //
        mask_tmp = _mm_and_si128(_mm_set_epi32(blk[mword_v[3]], blk[mword_v[2]],
                                               blk[mword_v[1]], blk[mword_v[0]]),
                                 mask_0);
        
        // bool(blk[nword]  ...)
        //maskFF  = _mm_set1_epi32(~0u);
        mask_tmp = _mm_cmpeq_epi32 (mask_tmp, maskZ); // set 0xFF where == 0
        mask_tmp = _mm_xor_si128 (mask_tmp, maskFF); // invert
        mask_tmp = _mm_srli_epi32 (mask_tmp, 31); // (bool) 1 only to the 0 pos
        
        mask_tmp = _mm_slli_epi32(mask_tmp, bit_idx); // << bit_idx
        
        _mm_storeu_si128 (target_ptr,                  // arr[base] |= MASK_EXPR
                          _mm_or_si128 (mask_tmp, _mm_loadu_si128(target_ptr)));
        
        ++idx_ptr; ++target_ptr;
        _mm_prefetch((const char*)target_ptr, _MM_HINT_T0);
    }

    for (; k < len; ++k)
    {
        base = start + k;
        unsigned nbit = unsigned(idx[base] & bm::set_block_mask);
        arr[base] |= unsigned(bool(blk[nbit >> bm::set_word_shift] & (1u << (nbit & bm::set_word_mask))) << bit_idx);
    }

}


#ifdef __GNUG__
#pragma GCC diagnostic pop
#endif


#ifdef __GNUG__
#pragma GCC diagnostic pop
#endif


} // namespace




#endif
