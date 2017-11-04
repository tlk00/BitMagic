#ifndef BMAVX2_UTIL__H__INCLUDED__
#define BMAVX2_UTIL__H__INCLUDED__
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



namespace bm
{

/** @defgroup AVX2 AVX2 functions
    Processor specific optimizations for AVX2 instructions (internals)
    @ingroup bvector
 */



/*! 
    @brief XOR array elements to specified mask
    *dst = *src ^ mask

    @ingroup AVX2
*/
BMFORCEINLINE 
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
BMFORCEINLINE 
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
BMFORCEINLINE 
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

    @ingroup SSE2
*/
BMFORCEINLINE 
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
BMFORCEINLINE 
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
BMFORCEINLINE 
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
BMFORCEINLINE 
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
BMFORCEINLINE 
void avx2_invert_arr(bm::word_t* first, bm::word_t* last)
{
    __m256i ymm1 = _mm256_set_epi32(0xFFFFFFFF, 0xFFFFFFFF,
                                    0xFFFFFFFF, 0xFFFFFFFF,
                                    0xFFFFFFFF, 0xFFFFFFFF,
                                    0xFFFFFFFF, 0xFFFFFFFF
                                    );
    __m256i* wrd_ptr = (__m256i*)first;
    __m256i ymm0;
    do 
    {
        ymm0 = _mm256_load_si256(wrd_ptr);
        ymm0 = _mm256_xor_si256(ymm0, ymm1);
        _mm256_store_si256(wrd_ptr, ymm0);
        ++wrd_ptr;

        ymm0 = _mm256_load_si256(wrd_ptr);
        ymm0 = _mm256_xor_si256(ymm0, ymm1);
        _mm256_store_si256(wrd_ptr, ymm0);
        ++wrd_ptr;

        ymm0 = _mm256_load_si256(wrd_ptr);
        ymm0 = _mm256_xor_si256(ymm0, ymm1);
        _mm256_store_si256(wrd_ptr, ymm0);
        ++wrd_ptr;

        ymm0 = _mm256_load_si256(wrd_ptr);
        ymm0 = _mm256_xor_si256(ymm0, ymm1);
        _mm256_store_si256(wrd_ptr, ymm0);
        ++wrd_ptr;

    } while (wrd_ptr < (__m256i*)last);
}



} // namespace



#endif
