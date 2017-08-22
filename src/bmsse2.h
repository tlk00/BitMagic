#ifndef BMSSE2__H__INCLUDED__
#define BMSSE2__H__INCLUDED__
/*
Copyright(c) 2002-2009 Anatoliy Kuznetsov(anatoliy_kuznetsov at yahoo.com)

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

For more information please visit:  http://bmagic.sourceforge.net

*/



//    Header implements processor specific intrinsics declarations for SSE2
//    instruction set
#include<mmintrin.h>
#include<emmintrin.h>

#include "bmdef.h"
#include "bmsse_util.h"


namespace bm
{


/*!
    SSE2 optimized bitcounting function implements parallel bitcounting
    algorithm for SSE2 instruction set.

<pre>
unsigned CalcBitCount32(unsigned b)
{
    b = (b & 0x55555555) + (b >> 1 & 0x55555555);
    b = (b & 0x33333333) + (b >> 2 & 0x33333333);
    b = (b + (b >> 4)) & 0x0F0F0F0F;
    b = b + (b >> 8);
    b = (b + (b >> 16)) & 0x0000003F;
    return b;
}
</pre>

    @ingroup SSE2

*/
inline 
bm::id_t sse2_bit_count(const __m128i* block, const __m128i* block_end)
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
    __m128i mcnt;
    mcnt = _mm_xor_si128(m1, m1); // cnt = 0

    __m128i tmp1, tmp2;
    do
    {        
        __m128i b = _mm_load_si128(block);
        ++block;

        // b = (b & 0x55555555) + (b >> 1 & 0x55555555);
        tmp1 = _mm_srli_epi32(b, 1);                    // tmp1 = (b >> 1 & 0x55555555)
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
        b = _mm_and_si128(b, m3);                       //           & 0x0F0F0F0F

        // b = b + (b >> 8);
        tmp1 = _mm_srli_epi32 (b, 8);                   // tmp1 = b >> 8
        b = _mm_add_epi32(b, tmp1);                     // b = b + (b >> 8)

        // b = (b + (b >> 16)) & 0x0000003F;
        tmp1 = _mm_srli_epi32 (b, 16);                  // b >> 16
        b = _mm_add_epi32(b, tmp1);                     // b + (b >> 16)
        b = _mm_and_si128(b, m4);                       // (b >> 16) & 0x0000003F;

        mcnt = _mm_add_epi32(mcnt, b);                  // mcnt += b

    } while (block < block_end);


    bm::id_t BM_ALIGN16 tcnt[4] BM_ALIGN16ATTR;
    _mm_store_si128((__m128i*)tcnt, mcnt);

    return tcnt[0] + tcnt[1] + tcnt[2] + tcnt[3];
}



template<class Func>
bm::id_t sse2_bit_count_op(const __m128i* BMRESTRICT block, 
                           const __m128i* BMRESTRICT block_end,
                           const __m128i* BMRESTRICT mask_block,
                           Func sse2_func)
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
    __m128i mcnt;
    mcnt = _mm_xor_si128(m1, m1); // cnt = 0
    do
    {
        __m128i tmp1, tmp2;
        __m128i b = _mm_load_si128(block++);

        tmp1 = _mm_load_si128(mask_block++);
        
        b = sse2_func(b, tmp1);
                        
        // b = (b & 0x55555555) + (b >> 1 & 0x55555555);
        tmp1 = _mm_srli_epi32(b, 1);                    // tmp1 = (b >> 1 & 0x55555555)
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
        b = _mm_and_si128(b, m3);                       //           & 0x0F0F0F0F

        // b = b + (b >> 8);
        tmp1 = _mm_srli_epi32 (b, 8);                   // tmp1 = b >> 8
        b = _mm_add_epi32(b, tmp1);                     // b = b + (b >> 8)
        
        // b = (b + (b >> 16)) & 0x0000003F;
        tmp1 = _mm_srli_epi32 (b, 16);                  // b >> 16
        b = _mm_add_epi32(b, tmp1);                     // b + (b >> 16)
        b = _mm_and_si128(b, m4);                       // (b >> 16) & 0x0000003F;

        mcnt = _mm_add_epi32(mcnt, b);                  // mcnt += b

    } while (block < block_end);

    bm::id_t BM_ALIGN16 tcnt[4] BM_ALIGN16ATTR;
    _mm_store_si128((__m128i*)tcnt, mcnt);

    return tcnt[0] + tcnt[1] + tcnt[2] + tcnt[3];
}




#define VECT_XOR_ARR_2_MASK(dst, src, src_end, mask)\
    sse2_xor_arr_2_mask((__m128i*)(dst), (__m128i*)(src), (__m128i*)(src_end), mask)

#define VECT_ANDNOT_ARR_2_MASK(dst, src, src_end, mask)\
    sse2_andnot_arr_2_mask((__m128i*)(dst), (__m128i*)(src), (__m128i*)(src_end), mask)

#define VECT_BITCOUNT(first, last) \
    sse2_bit_count((__m128i*) (first), (__m128i*) (last)) 

#define VECT_BITCOUNT_AND(first, last, mask) \
    sse2_bit_count_op((__m128i*) (first), (__m128i*) (last), (__m128i*) (mask), sse2_and) 

#define VECT_BITCOUNT_OR(first, last, mask) \
    sse2_bit_count_op((__m128i*) (first), (__m128i*) (last), (__m128i*) (mask), sse2_or) 

#define VECT_BITCOUNT_XOR(first, last, mask) \
    sse2_bit_count_op((__m128i*) (first), (__m128i*) (last), (__m128i*) (mask), sse2_xor) 

#define VECT_BITCOUNT_SUB(first, last, mask) \
    sse2_bit_count_op((__m128i*) (first), (__m128i*) (last), (__m128i*) (mask), sse2_sub) 

#define VECT_INVERT_ARR(first, last) \
    sse2_invert_arr(first, last);

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
   __m128i mcnt, ccnt;
   mcnt = _mm_xor_si128(m1, m1); // bit_cnt = 0
   ccnt = _mm_xor_si128(m1, m1); // change_cnt = 0

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

   bm::id_t BM_ALIGN16 tcnt[4] BM_ALIGN16ATTR;

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

} // namespace




#endif
