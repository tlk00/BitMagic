#ifndef BMUTIL__H__INCLUDED__
#define BMUTIL__H__INCLUDED__
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

#include "bmdef.h"
#include "bmconst.h"

#if defined(_M_AMD64) || defined(_M_X64) 
#include <intrin.h>
#endif

namespace bm
{


/*
// SSE2 vector equality comparison
int _mm_any_eq( vFloat a, vFloat b )
{

    //test a==b for each float in a & b
    vFloat mask = _mm_cmpeq_ps( a, b );

    //copy top bit of each result to maskbits
    int maskBits = _mm_movemask_ps( mask );
    return maskBits != 0;
}
*/

// From:
// http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.37.8562
//
template<typename T>
T bit_scan_fwd(T v)
{
    return 
        DeBruijn_bit_position<true>::_multiply[((word_t)((v & -v) * 0x077CB531U)) >> 27];
}

/**
    Get minimum of 2 values
*/
template<typename T>
T min_value(T v1, T v2)
{
    return v1 < v2 ? v1 : v2;
}


/**
    Fast loop-less function to find LOG2
*/
template<typename T>
T ilog2(T x)
{
    unsigned int l = 0;
    
    if (x >= 1<<16) { x = (T)(x >> 16); l |= 16; }
    if (x >= 1<<8)  { x = (T)(x >> 8);  l |= 8; }
    if (x >= 1<<4)  { x = (T)(x >> 4);  l |= 4; }
    if (x >= 1<<2)  { x = (T)(x >> 2);  l |= 2; }
    if (x >= 1<<1)  l |=1;
    return (T)l;
}

template<>
inline bm::gap_word_t ilog2(gap_word_t x)
{
    unsigned int l = 0;
    if (x >= 1<<8)  { x = (bm::gap_word_t)(x >> 8); l |= 8; }
    if (x >= 1<<4)  { x = (bm::gap_word_t)(x >> 4); l |= 4; }
    if (x >= 1<<2)  { x = (bm::gap_word_t)(x >> 2); l |= 2; }
    if (x >= 1<<1)  l |=1;
    return (bm::gap_word_t)l;
}

/**
     Mini auto-pointer for internal memory management
     @internal
*/
template<class T>
class ptr_guard
{
public:
    ptr_guard(T* p) : ptr_(p) {}
    ~ptr_guard() { delete ptr_; }
private:
    ptr_guard(const ptr_guard<T>& p);
    ptr_guard& operator=(const ptr_guard<T>& p);
private:
    T* ptr_;
};

/**
    Lookup table based integer LOG2
*/
template<typename T>
T ilog2_LUT(T x)
{
    unsigned l = 0;
    if (x & 0xffff0000) 
    {
        l += 16; x >>= 16;
    }
    
    if (x & 0xff00) 
    {
        l += 8; x >>= 8;
    }
    return l + first_bit_table<true>::_idx[x];
}

/**
    Lookup table based short integer LOG2
*/
template<>
inline bm::gap_word_t ilog2_LUT<bm::gap_word_t>(bm::gap_word_t x)
{
    unsigned l = 0;    
    if (x & 0xff00) 
    {
        l += 8;
        x = (bm::gap_word_t)(x >> 8);
    }
    return (bm::gap_word_t)(l + first_bit_table<true>::_idx[x]);
}


// if we are running on x86 CPU we can use inline ASM 

#ifdef BM_x86
#ifdef __GNUG__

inline 
unsigned bsf_asm32(unsigned int v)
{
    unsigned r;
    asm volatile(" bsfl %1, %0": "=r"(r): "rm"(v) );
    return r;
}
 
BMFORCEINLINE unsigned bsr_asm32(register unsigned int v)
{
    unsigned r;
    asm volatile(" bsrl %1, %0": "=r"(r): "rm"(v) );
    return r;
}

#endif  // __GNUG__

#ifdef _MSC_VER

#if defined(_M_AMD64) || defined(_M_X64) // inline assembly not supported

BMFORCEINLINE unsigned int bsr_asm32(unsigned int value)
{
    unsigned long r;
    _BitScanReverse(&r, value);
    return r;
}

BMFORCEINLINE unsigned int bsf_asm32(unsigned int value)
{
    unsigned long r;
    _BitScanForward(&r, value);
    return r;
}

#else

BMFORCEINLINE unsigned int bsr_asm32(register unsigned int value)
{   
  __asm  bsr  eax, value
}

BMFORCEINLINE unsigned int bsf_asm32(register unsigned int value)
{   
  __asm  bsf  eax, value
}

#endif

#endif // _MSC_VER

#endif // BM_x86



} // bm

#endif
