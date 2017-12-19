#ifndef BMSSE_UTIL__H__INCLUDED__
#define BMSSE_UTIL__H__INCLUDED__
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

/** @defgroup SSE2 SSE2 functions
    Processor specific optimizations for SSE2 instructions (internals)
    @internal
    @ingroup bvector
 */


/*! 
  @brief SSE2 reinitialization guard class

  SSE2 requires to call _mm_empty() if we are intermixing
  MMX integer commands with floating point arithmetics.
  This class guards critical code fragments where SSE2 integer
  is used.

  As of 2015 _mm_empty() is considered deprecated, and not even recognised
  by some compilers (like MSVC) in 64-bit mode.
  As MMX instructions gets old, we here deprecate and comment out 
  use of _mm_empty()

  @ingroup SSE2
*/
class sse_empty_guard
{
public:
    BMFORCEINLINE sse_empty_guard() 
    {
        //_mm_empty();
    }

    BMFORCEINLINE ~sse_empty_guard() 
    {
        //_mm_empty();
    }
};



/*! 
    @brief XOR array elements to specified mask
    *dst = *src ^ mask

    @ingroup SSE2
*/
BMFORCEINLINE 
void sse2_xor_arr_2_mask(__m128i* BMRESTRICT dst, 
                         const __m128i* BMRESTRICT src, 
                         const __m128i* BMRESTRICT src_end,
                         bm::word_t mask)
{
     __m128i xmm2 = _mm_set_epi32(mask, mask, mask, mask);
     do
     {
        __m128i xmm1 = _mm_load_si128(src);

        xmm1 = _mm_xor_si128(xmm1, xmm2);
        _mm_store_si128(dst, xmm1);
        ++dst;
        ++src;

     } while (src < src_end);
}


/*! 
    @brief Inverts array elements and NOT them to specified mask
    *dst = ~*src & mask

    @ingroup SSE2
*/
BMFORCEINLINE 
void sse2_andnot_arr_2_mask(__m128i* BMRESTRICT dst, 
                            const __m128i* BMRESTRICT src, 
                            const __m128i* BMRESTRICT src_end,
                            bm::word_t mask)
{
     __m128i xmm2 = _mm_set_epi32(mask, mask, mask, mask);
     do
     {
        //_mm_prefetch((const char*)(src)+1024, _MM_HINT_NTA);
        //_mm_prefetch((const char*)(src)+1088, _MM_HINT_NTA);

        __m128i xmm1 = _mm_load_si128(src);

        xmm1 = _mm_andnot_si128(xmm1, xmm2); // xmm1 = (~xmm1) & xmm2 
        _mm_store_si128(dst, xmm1);
        ++dst;
        ++src;

     } while (src < src_end);
}

/*! 
    @brief AND array elements against another array
    *dst &= *src

    @ingroup SSE2
*/
BMFORCEINLINE 
void sse2_and_arr(__m128i* BMRESTRICT dst, 
                  const __m128i* BMRESTRICT src, 
                  const __m128i* BMRESTRICT src_end)
{
    __m128i xmm1, xmm2;
    do
    {
        _mm_prefetch((const char*)(src)+512,  _MM_HINT_NTA);
    
        xmm1 = _mm_load_si128(src++);
        xmm2 = _mm_load_si128(dst);
        xmm1 = _mm_and_si128(xmm1, xmm2);
        _mm_store_si128(dst++, xmm1);
        
        xmm1 = _mm_load_si128(src++);
        xmm2 = _mm_load_si128(dst);
        xmm1 = _mm_and_si128(xmm1, xmm2);
        _mm_store_si128(dst++, xmm1);

        xmm1 = _mm_load_si128(src++);
        xmm2 = _mm_load_si128(dst);
        xmm1 = _mm_and_si128(xmm1, xmm2);
        _mm_store_si128(dst++, xmm1);

        xmm1 = _mm_load_si128(src++);
        xmm2 = _mm_load_si128(dst);
        xmm1 = _mm_and_si128(xmm1, xmm2);
        _mm_store_si128(dst++, xmm1);

    } while (src < src_end);

}


/*! 
    @brief OR array elements against another array
    *dst |= *src

    @ingroup SSE2
*/
BMFORCEINLINE 
void sse2_or_arr(__m128i* BMRESTRICT dst, 
                 const __m128i* BMRESTRICT src, 
                 const __m128i* BMRESTRICT src_end)
{
    __m128i xmm1, xmm2;
    do
    {
        _mm_prefetch((const char*)(src)+512,  _MM_HINT_NTA);
    
        xmm1 = _mm_load_si128(src++);
        xmm2 = _mm_load_si128(dst);
        xmm1 = _mm_or_si128(xmm1, xmm2);
        _mm_store_si128(dst++, xmm1);
        
        xmm1 = _mm_load_si128(src++);
        xmm2 = _mm_load_si128(dst);
        xmm1 = _mm_or_si128(xmm1, xmm2);
        _mm_store_si128(dst++, xmm1);

        xmm1 = _mm_load_si128(src++);
        xmm2 = _mm_load_si128(dst);
        xmm1 = _mm_or_si128(xmm1, xmm2);
        _mm_store_si128(dst++, xmm1);

        xmm1 = _mm_load_si128(src++);
        xmm2 = _mm_load_si128(dst);
        xmm1 = _mm_or_si128(xmm1, xmm2);
        _mm_store_si128(dst++, xmm1);

    } while (src < src_end);
}


/*! 
    @brief OR array elements against another array
    *dst ^= *src

    @ingroup SSE2
*/
BMFORCEINLINE 
void sse2_xor_arr(__m128i* BMRESTRICT dst, 
                  const __m128i* BMRESTRICT src, 
                  const __m128i* BMRESTRICT src_end)
{
    __m128i xmm1, xmm2;
    do
    {
        _mm_prefetch((const char*)(src)+512,  _MM_HINT_NTA);
    
        xmm1 = _mm_load_si128(src++);
        xmm2 = _mm_load_si128(dst);
        xmm1 = _mm_xor_si128(xmm1, xmm2);
        _mm_store_si128(dst++, xmm1);
        
        xmm1 = _mm_load_si128(src++);
        xmm2 = _mm_load_si128(dst);
        xmm1 = _mm_xor_si128(xmm1, xmm2);
        _mm_store_si128(dst++, xmm1);

        xmm1 = _mm_load_si128(src++);
        xmm2 = _mm_load_si128(dst);
        xmm1 = _mm_xor_si128(xmm1, xmm2);
        _mm_store_si128(dst++, xmm1);

        xmm1 = _mm_load_si128(src++);
        xmm2 = _mm_load_si128(dst);
        xmm1 = _mm_xor_si128(xmm1, xmm2);
        _mm_store_si128(dst++, xmm1);

    } while (src < src_end);
}



/*! 
    @brief AND-NOT (SUB) array elements against another array
    *dst &= ~*src

    @ingroup SSE2
*/
BMFORCEINLINE 
void sse2_sub_arr(__m128i* BMRESTRICT dst, 
                 const __m128i* BMRESTRICT src, 
                 const __m128i* BMRESTRICT src_end)
{
    __m128i xmm1, xmm2;
    do
    {
        _mm_prefetch((const char*)(src)+512,  _MM_HINT_NTA);
    
        xmm1 = _mm_load_si128(src++);
        xmm2 = _mm_load_si128(dst);
        xmm1 = _mm_andnot_si128(xmm1, xmm2);
        _mm_store_si128(dst++, xmm1);
        
        xmm1 = _mm_load_si128(src++);
        xmm2 = _mm_load_si128(dst);
        xmm1 = _mm_andnot_si128(xmm1, xmm2);
        _mm_store_si128(dst++, xmm1);

        xmm1 = _mm_load_si128(src++);
        xmm2 = _mm_load_si128(dst);
        xmm1 = _mm_andnot_si128(xmm1, xmm2);
        _mm_store_si128(dst++, xmm1);

        xmm1 = _mm_load_si128(src++);
        xmm2 = _mm_load_si128(dst);
        xmm1 = _mm_andnot_si128(xmm1, xmm2);
        _mm_store_si128(dst++, xmm1);

    } while (src < src_end);    
}

/*! 
    @brief SSE2 block memset
    *dst = value

    @ingroup SSE2
*/

BMFORCEINLINE 
void sse2_set_block(__m128i* BMRESTRICT dst, 
                    __m128i* BMRESTRICT dst_end, 
                    bm::word_t value)
{
    __m128i xmm0 = _mm_set_epi32 (value, value, value, value);
    do
    {            
        _mm_store_si128(dst, xmm0);
/*        
        _mm_store_si128(dst+1, xmm0);
        _mm_store_si128(dst+2, xmm0);
        _mm_store_si128(dst+3, xmm0);

        _mm_store_si128(dst+4, xmm0);
        _mm_store_si128(dst+5, xmm0);
        _mm_store_si128(dst+6, xmm0);
        _mm_store_si128(dst+7, xmm0);

        dst += 8;
*/        
    } while (++dst < dst_end);
    
    _mm_sfence();
}



/*! 
    @brief SSE2 block copy
    *dst = *src

    @ingroup SSE2
*/
BMFORCEINLINE 
void sse2_copy_block(__m128i* BMRESTRICT dst, 
                     const __m128i* BMRESTRICT src, 
                     const __m128i* BMRESTRICT src_end)
{
    __m128i xmm0, xmm1, xmm2, xmm3;
    do
    {
        _mm_prefetch((const char*)(src)+512,  _MM_HINT_NTA);
    
        xmm0 = _mm_load_si128(src+0);
        xmm1 = _mm_load_si128(src+1);
        xmm2 = _mm_load_si128(src+2);
        xmm3 = _mm_load_si128(src+3);
        
        _mm_store_si128(dst+0, xmm0);
        _mm_store_si128(dst+1, xmm1);
        _mm_store_si128(dst+2, xmm2);
        _mm_store_si128(dst+3, xmm3);
        
        xmm0 = _mm_load_si128(src+4);
        xmm1 = _mm_load_si128(src+5);
        xmm2 = _mm_load_si128(src+6);
        xmm3 = _mm_load_si128(src+7);
        
        _mm_store_si128(dst+4, xmm0);
        _mm_store_si128(dst+5, xmm1);
        _mm_store_si128(dst+6, xmm2);
        _mm_store_si128(dst+7, xmm3);
        
        src += 8;
        dst += 8;
        
    } while (src < src_end);    
}

/*! 
    @brief Invert array elements
    *dst = ~*dst
    or
    *dst ^= *dst 

    @ingroup SSE2
*/
BMFORCEINLINE 
void sse2_invert_arr(bm::word_t* first, bm::word_t* last)
{
    __m128i xmm1 = _mm_set_epi32(0xFFFFFFFF, 0xFFFFFFFF, 
                                 0xFFFFFFFF, 0xFFFFFFFF);
    __m128i* wrd_ptr = (__m128i*)first;

    do 
    {
        _mm_prefetch((const char*)(wrd_ptr)+512,  _MM_HINT_NTA);
        
        __m128i xmm0 = _mm_load_si128(wrd_ptr);
        xmm0 = _mm_xor_si128(xmm0, xmm1);
        _mm_store_si128(wrd_ptr, xmm0);
        ++wrd_ptr;
    } while (wrd_ptr < (__m128i*)last);
}

BMFORCEINLINE 
__m128i sse2_and(__m128i a, __m128i b)
{
    return _mm_and_si128(a, b);
}

BMFORCEINLINE 
__m128i sse2_or(__m128i a, __m128i b)
{
    return _mm_or_si128(a, b);
}


BMFORCEINLINE 
__m128i sse2_xor(__m128i a, __m128i b)
{
    return _mm_xor_si128(a, b);
}

BMFORCEINLINE 
__m128i sse2_sub(__m128i a, __m128i b)
{
    return _mm_andnot_si128(b, a);
}



} // namespace



#endif
