#ifndef BMSSE4__H__INCLUDED__
#define BMSSE4__H__INCLUDED__
/*
Copyright(c) 2002-2017 Anatoliy Kuznetsov(anatoliy_kuznetsov at yahoo.com)

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



//    Header implements processor specific intrinsics declarations for SSE2
//    instruction set
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

#define VECT_SUB_ARR(dst, src, src_end) \
    sse2_sub_arr((__m128i*) dst, (__m128i*) (src), (__m128i*) (src_end))

#define VECT_XOR_ARR(dst, src, src_end) \
    sse2_xor_arr((__m128i*) dst, (__m128i*) (src), (__m128i*) (src_end))

#define VECT_COPY_BLOCK(dst, src, src_end) \
    sse2_copy_block((__m128i*) dst, (__m128i*) (src), (__m128i*) (src_end))

#define VECT_SET_BLOCK(dst, dst_end, value) \
    sse2_set_block((__m128i*) dst, (__m128i*) (dst_end), (value))





/*!
    SSE4.2 optimized bitcounting and number of GAPs
    @ingroup SSE4
*/
inline
bm::id_t sse4_bit_block_calc_count_change(const __m128i* BMRESTRICT block,
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

/*! 
    @brief Gap block population count (array sum) utility 
    @param pbuf - unrolled, aligned to 1-start GAP buffer
    @param sse_vect_waves - number of SSE vector lines to process
    @param sum - result acumulator
    @return tail pointer

    @internal
*/
inline
const bm::gap_word_t* sse4_gap_sum_arr(
    const bm::gap_word_t* BMRESTRICT pbuf,
    unsigned sse_vect_waves,
    unsigned* sum)
{
    __m128i xcnt = _mm_setzero_si128();

    for (unsigned i = 0; i < sse_vect_waves; ++i)
    {
        __m128i mm0 = _mm_loadu_si128((__m128i*)(pbuf - 1));
        __m128i mm1 = _mm_loadu_si128((__m128i*)(pbuf + 8 - 1));
        __m128i mm_s2 = _mm_add_epi16(mm1, mm0);
        xcnt = _mm_add_epi16(xcnt, mm_s2);
        pbuf += 16;
    }
    xcnt = _mm_sub_epi16(_mm_srli_epi32(xcnt, 16), xcnt);

    unsigned short* cnt8 = (unsigned short*)&xcnt;
    *sum += (cnt8[0]) + (cnt8[2]) + (cnt8[4]) + (cnt8[6]);
    return pbuf;
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

    const unsigned unroll_factor = 8;
    __m128i m1, maskF, maskF1, maskF2;

    switch (size)
    {
    case 1:
    case 2:
    case 3:
    case 4:
        if (pbuf[0] >= pos) return 0;
        if (pbuf[1] >= pos) return 1;
        if (pbuf[2] >= pos) return 2;
        if (pbuf[3] >= pos) return 3;
        if (pbuf[4] >= pos) return 4;
        BM_ASSERT(0);
        return 0;
    case 5:
        m1 = _mm_loadu_si128((__m128i*)(pbuf)); // load first 8 elements
        maskF = _mm_set1_epi16(-1); // set all to FF
        maskF1 = _mm_srli_si128(maskF, (unroll_factor - 5) * 2);
        maskF2 = _mm_slli_si128(maskF, 5 * 2);
        m1 = _mm_and_si128(m1, maskF1);
        m1 = _mm_or_si128(m1, maskF2);
        break;
    case 6:
        m1 = _mm_loadu_si128((__m128i*)(pbuf)); // load first 8 elements
        maskF = _mm_set1_epi16(-1); // set all to FF
        maskF1 = _mm_srli_si128(maskF, (unroll_factor - 6) * 2);
        maskF2 = _mm_slli_si128(maskF, 6 * 2);
        m1 = _mm_and_si128(m1, maskF1);
        m1 = _mm_or_si128(m1, maskF2);
        break;
    case 7:
        m1 = _mm_loadu_si128((__m128i*)(pbuf)); // load first 8 elements
        maskF = _mm_set1_epi16(-1); // set all to FF
        maskF1 = _mm_srli_si128(maskF, (unroll_factor - 7) * 2);
        maskF2 = _mm_slli_si128(maskF, 7 * 2);
        m1 = _mm_and_si128(m1, maskF1);
        m1 = _mm_or_si128(m1, maskF2);
        break;
    default:
        m1 = _mm_loadu_si128((__m128i*)(pbuf)); // load first 8 elements
    };
    
    __m128i mz = _mm_setzero_si128();
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

    BM_ASSERT(pbuf2 > pbuf); // assert in place to make sure GCC warning is indeed false

    m1 = _mm_loadu_si128((__m128i*)(pbuf2)); // load next elements (with possible overlap)
    mge_mask = _mm_cmpeq_epi16(_mm_subs_epu16(mp, m1), mz); // m1 >= mp
        
    BM_ASSERT(!_mm_test_all_zeros(mge_mask, mge_mask)); // something found

    mi = _mm_movemask_epi8(_mm_slli_epi16(mge_mask, 15));
    unsigned bc = _mm_popcnt_u32(mi); 

    return size - bc;
}
#ifdef __GNUG__
#pragma GCC diagnostic pop
#endif

} // namespace




#endif
