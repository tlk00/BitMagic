#ifndef BMFUNC__H__INCLUDED__
#define BMFUNC__H__INCLUDED__
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

/*! \file bmfunc.h
    \brief Bit manipulation primitives (internal)
*/

#include <memory.h>

#include "bmdef.h"
#include "bmutil.h"


#ifdef _MSC_VER
# pragma warning( disable: 4146 )
#endif

namespace bm
{


inline 
bm::id_t bit_block_calc_count_range(const bm::word_t* block,
                                    bm::word_t left,
                                    bm::word_t right);

inline 
bm::id_t bit_block_any_range(const bm::word_t* block,
                             bm::word_t left,
                             bm::word_t right);

/*!
    @brief Structure with statistical information about bitset's memory 
            allocation details. 
    @ingroup bvector
*/
struct bv_statistics
{
    /// Number of bit blocks.
    unsigned bit_blocks; 
    /// Number of GAP blocks.
    unsigned gap_blocks;  
    /// Estimated maximum of memory required for serialization.
    size_t  max_serialize_mem;
    /// Memory used by bitvector including temp and service blocks
    size_t  memory_used;
    /// Array of all GAP block lengths in the bvector.
    gap_word_t   gap_length[bm::set_total_blocks];
    /// GAP lengths used by bvector
    gap_word_t  gap_levels[bm::gap_levels];


    /// cound bit block
    void add_bit_block()
    {
        ++bit_blocks;
        unsigned mem_used = (unsigned)(sizeof(bm::word_t) * bm::set_block_size);
        memory_used += mem_used;
        max_serialize_mem += mem_used;
    }

    /// count gap block
    void add_gap_block(unsigned capacity, unsigned length)
    {
        (gap_blocks < bm::set_total_blocks) ? gap_length[gap_blocks] = (gap_word_t)length : 0;
        ++gap_blocks;
        unsigned mem_used = (unsigned)(capacity * sizeof(gap_word_t));
        memory_used += mem_used;
        max_serialize_mem += (unsigned)(length * sizeof(gap_word_t));
    }
    
    /// Reset statisctics
    void reset()
    {
        bit_blocks = gap_blocks = 0;
        max_serialize_mem = memory_used = 0;
    }
};

/**
    \brief ad-hoc conditional expressions 
    \internal
*/
template <bool b> struct conditional
{
    static bool test() { return true; }
};
template <> struct conditional<false>
{
    static bool test() { return false; }
};

/*! 
    @defgroup gapfunc GAP functions
    GAP functions implement different opereations on GAP compressed blocks (internals)
    and serve as a minimal building blocks.
    \internal
    @ingroup bvector
 */

/*! 
   @defgroup bitfunc BIT functions
   Bit functions implement different opereations on bit blocks (internals)
   and serve as a minimal building blocks.
   \internal
   @ingroup bvector
 */




/*!
    Returns BSR value
    @ingroup bitfunc
*/
template <class T>
unsigned bit_scan_reverse(T value)
{
    BM_ASSERT(value);
    
    if (bm::conditional<sizeof(T)==8>::test())
    {
        bm::id64_t v8 = value;
        v8 >>= 32;
        unsigned v = (unsigned)v8;
        if (v)
        {
            return bit_scan_reverse32(v) + 32;
        }
    }
    return bit_scan_reverse32((unsigned)value);
}


/*!
    Returns bit count
    @ingroup bitfunc 
*/
BMFORCEINLINE
bm::id_t word_bitcount(bm::id_t w)
{
#if defined(BMSSE42OPT) || defined(BMAVX2OPT)
    return bm::id_t(_mm_popcnt_u32(w));
#else
    return
    bm::bit_count_table<true>::_count[(unsigned char)(w)] + 
    bm::bit_count_table<true>::_count[(unsigned char)((w) >> 8)] + 
    bm::bit_count_table<true>::_count[(unsigned char)((w) >> 16)] + 
    bm::bit_count_table<true>::_count[(unsigned char)((w) >> 24)];
#endif
}

inline
int parallel_popcnt_32(unsigned int n) 
{
   unsigned int tmp;

   tmp = n - ((n >> 1) & 033333333333)
           - ((n >> 2) & 011111111111);
   return ((tmp + (tmp >> 3)) & 030707070707) % 63;
}


/*! 
    Function calculates number of 1 bits in 64-bit word.
    @ingroup bitfunc 
*/
BMFORCEINLINE
unsigned word_bitcount64(bm::id64_t x)
{
#if defined(BMSSE42OPT) || defined(BMAVX2OPT)
#if defined(BM64_SSE4) || defined(BM64_AVX2)
    return unsigned(_mm_popcnt_u64(x));
#else
    return _mm_popcnt_u32(x >> 32) + _mm_popcnt_u32((unsigned)x);
#endif
#else
    x = x - ((x >> 1) & 0x5555555555555555);
    x = (x & 0x3333333333333333) + ((x >> 2) & 0x3333333333333333);
    x = (x + (x >> 4)) & 0x0F0F0F0F0F0F0F0F;
    x = x + (x >> 8);
    x = x + (x >> 16);
    x = x + (x >> 32); 
    return x & 0xFF;
#endif
}

inline 
unsigned bitcount64_4way(bm::id64_t x, bm::id64_t y, 
                    bm::id64_t u, bm::id64_t v)
{
    const bm::id64_t m1 = 0x5555555555555555;
    const bm::id64_t m2 = 0x3333333333333333; 
    const bm::id64_t m3 = 0x0F0F0F0F0F0F0F0F; 
    const bm::id64_t m4 = 0x000000FF000000FF;

    x = x - ((x >> 1) & m1);
    y = y - ((y >> 1) & m1);
    u = u - ((u >> 1) & m1);
    v = v - ((v >> 1) & m1);
    x = (x & m2) + ((x >> 2) & m2);
    y = (y & m2) + ((y >> 2) & m2);
    u = (u & m2) + ((u >> 2) & m2);
    v = (v & m2) + ((v >> 2) & m2);
    x = x + y; 
    u = u + v; 
    x = (x & m3) + ((x >> 4) & m3);
    u = (u & m3) + ((u >> 4) & m3);
    x = x + u; 
    x = x + (x >> 8);
    x = x + (x >> 16);
    x = x & m4; 
    x = x + (x >> 32);
    return x & 0x000001FF;
}



/*!
    Returns number of trailing zeros
    @ingroup bitfunc 
*/
BMFORCEINLINE
bm::id_t word_trailing_zeros(bm::id_t w)
{
    // TODO: find a better variant for MSVC 
#if (defined(BMSSE42OPT) || defined(BMAVX2OPT)) && defined(__GNUC__)
        return bm::id_t(__builtin_ctzl(w));
#else
    // implementation from
    // https://gist.github.com/andrewrk/1883543
    static const int mod37_pos[] =
    {
        -1, 0, 1, 26, 2, 23, 27, 0, 3, 16, 24, 30, 28, 11, 0, 13, 4,
        7, 17, 0, 25, 22, 31, 15, 29, 10, 12, 6, 0, 21, 14, 9, 5,
        20, 8, 19, 18
    };
    return unsigned(mod37_pos[(-w & w) % 37]);
#endif
}




//---------------------------------------------------------------------


/// Returns true if set operation is constant (bitcount)
inline
bool is_const_set_operation(set_operation op)
{
    return (int(op) >= int(set_COUNT));
}

/**
    Bit operations enumeration.
*/
enum operation
{
    BM_AND = set_AND,
    BM_OR  = set_OR,
    BM_SUB = set_SUB,
    BM_XOR = set_XOR
};

/**
    Convert set operation to operation
*/
inline
bm::operation setop2op(bm::set_operation op)
{
    BM_ASSERT(op == set_AND || 
              op == set_OR  || 
              op == set_SUB || 
              op == set_XOR);
    return (bm::operation) op;
}

//---------------------------------------------------------------------

/** 
    Structure carries pointer on bit block with all bits 1
    @ingroup bitfunc 
*/
template<bool T> struct all_set
{
    struct BM_VECT_ALIGN all_set_block
    {
        bm::word_t BM_VECT_ALIGN _p[bm::set_block_size] BM_VECT_ALIGN_ATTR;
        bm::word_t* _p_fullp;

        all_set_block()
        {
            ::memset(_p, 0xFF, sizeof(_p)); // set FULL BLOCK content (all 1s)
            if (bm::conditional<sizeof(void*) == 8>::test())
            {
                const unsigned long long magic_mask = 0xFFFFfffeFFFFfffe;
                ::memcpy(&_p_fullp, &magic_mask, sizeof(magic_mask));
            }
            else 
            {
                const unsigned magic_mask = 0xFFFFfffe;
                ::memcpy(&_p_fullp, &magic_mask, sizeof(magic_mask));
            }
        }
    };

    // version with minimal branching, super-scalar friendly
    //

//#if defined(BM64OPT) || defined(BM64_SSE4) || defined(BM64_AVX2)

    inline
    static bm::id64_t block_type(const bm::word_t* bp)
    {
        bm::id64_t type;
        if (bm::conditional<sizeof(void*) == 8>::test())
        {
            bm::id64_t w = reinterpret_cast<unsigned long long>(bp);
            type = (w & 3) | // FULL BLOCK or GAP
                ((bp == _block._p) << 1);
            type = type ? type : w;
        }
        else
        {
            unsigned w = reinterpret_cast<unsigned long>(bp);
            type = (w & 3) | // FULL BLOCK or GAP
                ((bp == _block._p) << 1);
            type = type ? type : w;
        }
        return type;
    }
    /*
    inline
    static bool is_invalid_addr(const unsigned* bp)
    {
        return 
            (!bp) | ((reinterpret_cast<unsigned long long>(bp) >> 1u) & 1u)
                  | (bp == _block._p);
    }

    inline
    static bool is_full_block(const bm::word_t* bp)
    {
        return ((reinterpret_cast<unsigned long long>(bp) >> 1u) & 1u) | 
               (bp == _block._p);
    }
    inline
    static bool is_valid_block_addr(const bm::word_t* bp)
    {
        return !is_invalid_addr(bp);
    }
    */
//#endif

    BMFORCEINLINE 
    static bool is_full_block(const bm::word_t* bp) 
        { return (bp == _block._p || bp == _block._p_fullp); }

    BMFORCEINLINE 
    static bool is_valid_block_addr(const bm::word_t* bp) 
        { return (bp && !(bp == _block._p || bp == _block._p_fullp)); }


    static all_set_block  _block;
};


template<bool T> typename all_set<T>::all_set_block all_set<T>::_block;

/// XOR swap two scalar variables
template<typename W> 
void xor_swap(W& x, W& y) 
{
    BM_ASSERT(&x != &y);
    x ^= y;
    y ^= x;
    x ^= y;
}


//---------------------------------------------------------------------

/*! 
   \brief Lexicographical comparison of two words as bit strings (reference)
   Auxiliary implementation for testing and reference purposes.
   \param w1 - First word.
   \param w2 - Second word.
   \return  <0 - less, =0 - equal,  >0 - greater.

   @ingroup bitfunc 
*/
template<typename T> int wordcmp0(T w1, T w2)
{
    while (w1 != w2)
    {
        int res = (w1 & 1) - (w2 & 1);
        if (res != 0) return res;
        w1 >>= 1;
        w2 >>= 1;
    }
    return 0;
}


/*
template<typename T> int wordcmp(T w1, T w2)
{
    T diff = w1 ^ w2;
    return diff ? ((w1 & diff & (diff ^ (diff - 1)))? 1 : -1) : 0; 
}
*/
/*! 
   \brief Lexicographical comparison of two words as bit strings.
   Auxiliary implementation for testing and reference purposes.
   \param a - First word.
   \param b - Second word.
   \return  <0 - less, =0 - equal,  >0 - greater.

   @ingroup bitfunc 
*/

template<typename T> int wordcmp(T a, T b)
{
    T diff = a ^ b;
    return diff? ( (a & diff & -diff)? 1 : -1 ) : 0;
}


// Low bit extraction
// x & (x ^ (x-1))



/*! 
   \brief Byte orders recognized by the library.
*/
enum ByteOrder
{
    BigEndian    = 0,
    LittleEndian = 1
};


/**
    Internal structure. Different global settings.
*/
template<bool T> struct globals
{
    struct bo
    {
        ByteOrder  _byte_order;

        bo()
        {
            unsigned x;
            unsigned char *s = (unsigned char *)&x;
            s[0] = 1;
            s[1] = 2;
            s[2] = 3;
            s[3] = 4;

            if(x == 0x04030201) 
            {
                _byte_order = LittleEndian;
                return;
            }

            if(x == 0x01020304) 
            {
                _byte_order = BigEndian;
                return;
            }

            BM_ASSERT(0); // "Invalid Byte Order\n"
            _byte_order = LittleEndian;
        }
    };

    static bo  _bo;

    static ByteOrder byte_order() { return _bo._byte_order; }

};

template<bool T> typename globals<T>::bo globals<T>::_bo;

// ----------------------------------------------------------------------


/*! @brief Returns "true" if all bits in the block are 0
    @ingroup bitfunc
*/
inline
bool bit_is_all_zero(const bm::word_t* BMRESTRICT start)
{
#if defined(BMSSE42OPT) || defined(BMAVX2OPT)
    return VECT_IS_ZERO_BLOCK(start, start + bm::set_block_size);
#else
   const bm::wordop_t* BMRESTRICT blk = (bm::wordop_t*) (start);
   const bm::wordop_t* BMRESTRICT blk_end = (bm::wordop_t*)(start + bm::set_block_size);
   do
   {
        if (blk[0] | blk[1] | blk[2] | blk[3])
            return false;
        blk += 4;
   } while (blk < blk_end);
   return true;
#endif
}

// ----------------------------------------------------------------------

/*!
   \brief Checks if GAP block is all-zero.
   \param buf - GAP buffer pointer.
   \returns true if all-zero.

   @ingroup gapfunc
*/
BMFORCEINLINE
bool gap_is_all_zero(const bm::gap_word_t* buf)
{
    // (almost) branchless variant:
    return (!(*buf & 1u)) & (!(gap_max_bits - 1 - buf[1]));
    //return ((*buf & 1u) == 0) && (buf[1] == bm::gap_max_bits - 1);
}


/*!
   \brief Checks if GAP block is all-one.
   \param buf - GAP buffer pointer.
   \param set_max - max possible bitset length
   \returns true if all-one.

   @ingroup gapfunc
*/
template<typename T>
BMFORCEINLINE
bool gap_is_all_one(const T* buf, unsigned set_max)
{
    return ((*buf & 1u) && (buf[1] == set_max - 1));
}

/*!
   \brief Returs GAP block length.
   \param buf - GAP buffer pointer.
   \returns GAP block length.

   @ingroup gapfunc
*/
template<typename T> T gap_length(const T* buf)
{
    return (T)((*buf >> 3) + 1);
}


/*!
   \brief Returs GAP block capacity
   \param buf - GAP buffer pointer
   \param glevel_len - pointer on GAP header word
   \returns GAP block capacity.

   @ingroup gapfunc
*/
template<typename T>
unsigned gap_capacity(const T* buf, const T* glevel_len)
{
    return glevel_len[(*buf >> 1) & 3];
}


/*!
   \brief Returs GAP block capacity limit.
   \param buf - GAP buffer pointer.
   \param glevel_len - GAP lengths table (gap_len_table)
   \returns GAP block limit.

   @ingroup gapfunc
*/
template<typename T>
unsigned gap_limit(const T* buf, const T* glevel_len)
{
    return glevel_len[(*buf >> 1) & 3]-4;
}


/*!
   \brief Returs GAP blocks capacity level.
   \param buf - GAP buffer pointer.
   \returns GAP block capacity level.

   @ingroup gapfunc
*/
template<typename T>
T gap_level(const T* buf)
{
    return T((*buf >> 1) & 3u);
}


/*!
    \brief GAP block find the last set bit

    \param buf - GAP buffer pointer.
    \param last - index of the last 1 bit
 
    \return 0 if 1 bit was NOT found

    @ingroup gapfunc
*/
template<typename T>
unsigned gap_find_last(const T* buf, unsigned* last)
{
    BM_ASSERT(last);

    T is_set = (*buf) & 1u;
    T end = T((*buf) >> 3u);

    BM_ASSERT(buf[end] == bm::gap_max_bits - 1);

    is_set ^= T((end-1) & 1u);
    if (is_set)
    {
        *last = buf[end];
        return is_set;
    }
    *last = buf[--end];
    return end;
}

/*!
    \brief GAP block find the first set bit

    \param buf - GAP buffer pointer.
    \param last - index of the first 1 bit
 
    \return 0 if 1 bit was NOT found

    @ingroup gapfunc
*/
template<typename T>
unsigned gap_find_first(const T* buf, unsigned* first)
{
    BM_ASSERT(first);

    T is_set = (*buf) & 1u;
    if (is_set)
    {
        *first = 0;
        return is_set;
    }
    if (buf[1] == bm::gap_max_bits - 1)
        return 0;
    *first = buf[1] + 1;
    return 1;
}



/*
   \brief Binary search for the block where bit = pos located.
   \param buf - GAP buffer pointer.
   \param pos - index of the element.
   \param is_set - output. GAP value (0 or 1). 
   \return GAP index.
   @ingroup gapfunc
*/
template<typename T> 
unsigned gap_bfind(const T* buf, unsigned pos, unsigned* is_set)
{
    BM_ASSERT(pos < bm::gap_max_bits);
    *is_set = (*buf) & 1;

    unsigned start = 1;
    unsigned end = 1 + ((*buf) >> 3);

    while ( start != end )
    {
        unsigned curr = (start + end) >> 1;
        if ( buf[curr] < pos )
            start = curr + 1;
        else
            end = curr;
    }
    *is_set ^= ((start-1) & 1);
    return start; 
}


/*!
   \brief Tests if bit = pos is true.
   \param buf - GAP buffer pointer.
   \param pos - index of the element.
   \return true if position is in "1" gap
   @ingroup gapfunc
*/
template<typename T> unsigned gap_test(const T* buf, unsigned pos)
{
    BM_ASSERT(pos < bm::gap_max_bits);

    unsigned start = 1;
    unsigned end = 1 + ((*buf) >> 3);

    if (end - start < 10)
    {
        unsigned sv = *buf & 1;
        unsigned sv1= sv ^ 1;
        if (buf[1] >= pos) return sv;
        if (buf[2] >= pos) return sv1;
        if (buf[3] >= pos) return sv;
        if (buf[4] >= pos) return sv1;
        if (buf[5] >= pos) return sv;
        if (buf[6] >= pos) return sv1;
        if (buf[7] >= pos) return sv;
        if (buf[8] >= pos) return sv1;
        if (buf[9] >= pos) return sv;
        BM_ASSERT(0);
    }
    else
    {
        while (start != end)
        {
            unsigned curr = (start + end) >> 1;
            if (buf[curr] < pos)
                start = curr + 1;
            else
                end = curr;
        }
    }
    return ((*buf) & 1) ^ ((--start) & 1); 
}

/*!
    \brief Tests if bit = pos is true. Analog of bm::gap_test with SIMD unrolling.
    \param buf - GAP buffer pointer.
    \param pos - index of the element.
    \return true if position is in "1" gap
    @ingroup gapfunc
*/
template<typename T> 
unsigned gap_test_unr(const T* buf, const unsigned pos)
{
    BM_ASSERT(pos < bm::gap_max_bits);
#if defined(BMSSE2OPT)
    unsigned start = 1;
    unsigned end = 1 + ((*buf) >> 3);
    unsigned dsize = end - start;

    if (dsize < 17)
    {
        start = bm::sse2_gap_find(buf + 1, (bm::gap_word_t)pos, dsize);
        unsigned res = ((*buf) & 1) ^ ((start) & 1);
        BM_ASSERT(buf[start + 1] >= pos);
        BM_ASSERT(buf[start] < pos || (start == 0));
        BM_ASSERT(res == bm::gap_test(buf, pos));
        return res;
    }
    unsigned arr_end = end;
    while (start != end)
    {
        unsigned curr = (start + end) >> 1;
        if (buf[curr] < pos)
            start = curr + 1;
        else
            end = curr;

        unsigned size = end - start;
        if (size < 16)
        {
            size += (end != arr_end);
            unsigned idx = bm::sse2_gap_find(buf + start, (bm::gap_word_t)pos, size);
            start += idx;

            BM_ASSERT(buf[start] >= pos);
            BM_ASSERT(buf[start - 1] < pos || (start == 1));
            break;
        }
    }

    unsigned res = ((*buf) & 1) ^ ((--start) & 1);

    BM_ASSERT(res == bm::gap_test(buf, pos));
    return res;
//#endif
#elif defined(BMSSE42OPT)
    unsigned start = 1;
    unsigned end = 1 + ((*buf) >> 3);
    unsigned dsize = end - start;

    if (dsize < 17)
    {
        start = bm::sse4_gap_find(buf+1, (bm::gap_word_t)pos, dsize);
        unsigned res = ((*buf) & 1) ^ ((start) & 1);
        BM_ASSERT(buf[start+1] >= pos);
        BM_ASSERT(buf[start] < pos || (start==0));
        BM_ASSERT(res == bm::gap_test(buf, pos));
        return res;
    }
    unsigned arr_end = end;
    while (start != end)
    {
        unsigned curr = (start + end) >> 1;
        if (buf[curr] < pos)
            start = curr + 1;
        else
            end = curr;

        unsigned size = end - start;
        if (size < 16)
        {
            size += (end != arr_end);
            unsigned idx = bm::sse4_gap_find(buf + start, (bm::gap_word_t)pos, size);
            start += idx;

            BM_ASSERT(buf[start] >= pos);
            BM_ASSERT(buf[start - 1] < pos || (start == 1));
            break;
        }
    }
    
    unsigned res = ((*buf) & 1) ^ ((--start) & 1);

    BM_ASSERT(res == bm::gap_test(buf, pos));
    return res;
#else
    return bm::gap_test(buf, pos);
#endif
}


/*! For each non-zero block executes supplied function.
    \internal
*/
template<class T, class F> 
void for_each_nzblock(T*** root, unsigned size1, F& f)
{
    for (unsigned i = 0; i < size1; ++i)
    {
        T** blk_blk = root[i];
        if (!blk_blk) 
        {
            f.on_empty_top(i);
            continue;
        }

        unsigned non_empty_top = 0;
        unsigned r = i * bm::set_array_size;
        unsigned j = 0;
        do
        {
#ifdef BM64_AVX2
            if (!avx2_test_all_zero_wave(blk_blk + j))
            {
                non_empty_top = 1;
                T* blk0 = blk_blk[j + 0];
                T* blk1 = blk_blk[j + 1];
                T* blk2 = blk_blk[j + 2];
                T* blk3 = blk_blk[j + 3];

                unsigned block_idx = r + j + 0;
                if (blk0)
                    f(blk0, block_idx);
                else
                    f.on_empty_block(block_idx);

                if (blk1)
                    f(blk1, block_idx + 1);
                else
                    f.on_empty_block(block_idx + 1);

                if (blk2)
                    f(blk2, block_idx + 2);
                else
                    f.on_empty_block(block_idx + 2);

                if (blk3)
                    f(blk3, block_idx + 3);
                else
                    f.on_empty_block(block_idx + 3);
            }
            else
            {
                f.on_empty_block(r + j + 0); f.on_empty_block(r + j + 1);
                f.on_empty_block(r + j + 2); f.on_empty_block(r + j + 3);
            }
            j += 4;
#elif defined(BM64_SSE4)
            if (!sse42_test_all_zero_wave((blk_blk + j)))
            {
                non_empty_top = 1;
                T* blk0 = blk_blk[j + 0];
                T* blk1 = blk_blk[j + 1];

                unsigned block_idx = r + j + 0;
                if (blk0)
                    f(blk0, block_idx);
                else
                    f.on_empty_block(block_idx);

                ++block_idx;
                if (blk1)
                    f(blk1, block_idx);
                else
                    f.on_empty_block(block_idx);
            }
            else
            {
                f.on_empty_block(r + j + 0);
                f.on_empty_block(r + j + 1);
            }
            j += 2;
#else
            if (blk_blk[j])
            {
                f(blk_blk[j], r + j);
                non_empty_top = 1;
            }
            else
                f.on_empty_block(r + j);
            ++j;
#endif
        } while (j < bm::set_array_size);

        if (non_empty_top == 0)
            f.on_empty_top(i);
    }  // for i
}

/*! For each non-zero block executes supplied function.
*/
template<class T, class F> 
void for_each_nzblock2(T*** root, unsigned size1, F& f)
{
#ifdef BM64_SSE4
    for (unsigned i = 0; i < size1; ++i)
    {
        T** blk_blk;
        if ((blk_blk = root[i])!=0)
        {
            {
                __m128i w0;
                for (unsigned j = 0; j < bm::set_array_size; j+=4)
                {
                    w0 = _mm_loadu_si128((__m128i*)(blk_blk + j));
                    if (!_mm_testz_si128(w0, w0))
                    {
                        T* blk0 = blk_blk[j + 0];
                        T* blk1 = blk_blk[j + 1];

                        if (blk0)
                            f(blk0);
                        if (blk1)
                            f(blk1);
                    }
                    w0 = _mm_loadu_si128((__m128i*)(blk_blk + j + 2));
                    if (!_mm_testz_si128(w0, w0))
                    {
                        T* blk0 = blk_blk[j + 2];
                        T* blk1 = blk_blk[j + 3];
                        if (blk0)
                            f(blk0);
                        if (blk1)
                            f(blk1);
                    }
                } // for j
            }
        }
    }  // for i
#elif defined(BM64_AVX2)
    for (unsigned i = 0; i < size1; ++i)
    {
        T** blk_blk;
        if ((blk_blk = root[i]) != 0)
        {
            {
                for (unsigned j = 0; j < bm::set_array_size; j += 4)
                {
                    __m256i w0 = _mm256_loadu_si256((__m256i*)(blk_blk + j));
                    if (!_mm256_testz_si256(w0, w0))
                    {
                        // as a variant could use: blk0 = (T*)_mm256_extract_epi64(w0, 0);
                        // but it measures marginally slower
                        T* blk0 = blk_blk[j + 0];
                        T* blk1 = blk_blk[j + 1];
                        T* blk2 = blk_blk[j + 2];
                        T* blk3 = blk_blk[j + 3];
                        if (blk0)
                            f(blk0);
                        if (blk1)
                            f(blk1);
                        if (blk2)
                            f(blk2);
                        if (blk3)
                            f(blk3);
                    }
                } // for j
            }
        }
    }  // for i

#else // non-SIMD mode
    for (unsigned i = 0; i < size1; ++i)
    {
        T** blk_blk;
        if ((blk_blk = root[i])!=0) 
        {
            unsigned j = 0;
            do
            {
                if (blk_blk[j])
                    f(blk_blk[j]);
                if (blk_blk[j+1])
                    f(blk_blk[j+1]);
                if (blk_blk[j+2])
                    f(blk_blk[j+2]);
                if (blk_blk[j+3])
                    f(blk_blk[j+3]);
                j += 4;
            } while (j < bm::set_array_size);
        }
    }  // for i
#endif
}


/*! For each non-zero block executes supplied function-predicate.
    Function returns if function-predicate returns true
*/
template<class T, class F> 
bool for_each_nzblock_if(T*** root, unsigned size1, F& f)
{
    unsigned block_idx = 0;
    for (unsigned i = 0; i < size1; ++i)
    {
        T** blk_blk = root[i];
        if (!blk_blk) 
        {
            block_idx += bm::set_array_size;
            continue;
        }

        for (unsigned j = 0;j < bm::set_array_size; ++j, ++block_idx)
        {
            if (blk_blk[j]) 
                if (f(blk_blk[j], block_idx)) return true;
        }
    }
    return false;
}

/*! For each block executes supplied function.
*/
template<class T, class F> 
void for_each_block(T*** root, unsigned size1, F& f)
{
    unsigned block_idx = 0;

    for (unsigned i = 0; i < size1; ++i)
    {
        T** blk_blk = root[i];

        if (blk_blk)
        {
            for (unsigned j = 0;j < bm::set_array_size; ++j, ++block_idx)
            {
                f(blk_blk[j], block_idx);
            }
        }
        else
        {
            for (unsigned j = 0;j < bm::set_array_size; ++j, ++block_idx)
            {
                f(0, block_idx);
            }
        }
    }  
}



/*! Special BM optimized analog of STL for_each
*/
template<class T, class F> F bmfor_each(T first, T last, F f)
{
    do
    {
        f(*first);
        ++first;
    } while (first < last);
    return f;
}

/*! Computes SUM of all elements of the sequence
*/
template<class T> T sum_arr(T* first, T* last)
{
    T sum = 0;
    while (first < last)
    {
        sum += *first;
        ++first;
    }
    return sum;
}




/*!
   \brief Calculates number of bits ON in GAP buffer.
   \param buf - GAP buffer pointer.
   \param dsize - buffer size
   \return Number of non-zero bits.
   @ingroup gapfunc
*/
template<typename T> unsigned gap_bit_count(const T* buf, unsigned dsize=0) 
{
    const T* pcurr = buf;
    if (dsize == 0)
        dsize = (*pcurr >> 3);

    const T* pend = pcurr + dsize;

    unsigned bits_counter = 0;
    ++pcurr;

    if (*buf & 1)
    {
        bits_counter += *pcurr + 1;
        ++pcurr;
    }
    ++pcurr;  // set GAP to 1

    while (pcurr <= pend)
    {
        bits_counter += *pcurr - *(pcurr-1);
        pcurr += 2; // jump to the next positive GAP
    } 

    return bits_counter;
}

/*!
    \brief Calculates number of bits ON in GAP buffer. Loop unrolled version.
    \param buf - GAP buffer pointer.
    \param dsize - buffer size
    \return Number of non-zero bits.
    @ingroup gapfunc
*/
template<typename T> unsigned gap_bit_count_unr(const T* buf)
{
    const T* pcurr = buf;
    unsigned dsize = (*pcurr >> 3);

    unsigned cnt = 0;
    pcurr = buf + 1; // set up start position
    T first_one = *buf & 1;
    if (first_one)
    {
        cnt += *pcurr + 1;
        ++pcurr;
    }
    ++pcurr;  // set GAP to 1

    #if defined(BMAVX2OPT)
    if (dsize > 34)
    {
        const unsigned unr_factor = 32;
        unsigned waves = (dsize-2) / unr_factor;
        pcurr = avx2_gap_sum_arr(pcurr, waves, &cnt);
    }
    #elif defined(BMSSE42OPT) || defined(BMSSE2OPT)
    if (dsize > 18)
    {
        const unsigned unr_factor = 16;
        unsigned waves = (dsize - 2) / unr_factor;
        pcurr = sse2_gap_sum_arr(pcurr, waves, &cnt);
    }
    #else
    if (dsize > 10)
    {
        const unsigned unr_factor = 8;
        unsigned waves = (dsize - 2) / unr_factor;
        for (unsigned i = 0; i < waves; i += unr_factor)
        {
            cnt += pcurr[0] - pcurr[0 - 1];
            cnt += pcurr[2] - pcurr[2 - 1];
            cnt += pcurr[4] - pcurr[4 - 1];
            cnt += pcurr[6] - pcurr[6 - 1];

            pcurr += unr_factor;
        } // for
    }
    #endif
    
    const T* pend = buf + dsize;
    for ( ; pcurr <= pend ; pcurr+=2)
    {
        cnt += *pcurr - *(pcurr - 1);
    }
    BM_ASSERT(cnt == gap_bit_count(buf));
    return cnt;
}



/*!
   \brief Counts 1 bits in GAP buffer in the closed [left, right] range.
   \param buf - GAP buffer pointer.
   \param left - leftmost bit index to start from
   \param right- rightmost bit index
   \return Number of non-zero bits.
   @ingroup gapfunc
*/
template<typename T>
unsigned gap_bit_count_range(const T* const buf, T left, T right)
{
    BM_ASSERT(left <= right);
    
    const T* pcurr = buf;
    const T* pend = pcurr + (*pcurr >> 3);
    
    unsigned bits_counter = 0;
    unsigned is_set;
    unsigned start_pos = bm::gap_bfind(buf, left, &is_set);
    is_set = ~(is_set - 1u); // 0xFFF.. if true (mask for branchless code)

    pcurr = buf + start_pos;
    if (right <= *pcurr) // we are in the target gap right now
    {
        bits_counter = unsigned(right - left + 1u) & is_set; // & is_set == if(is_set)
        return bits_counter;
    }
    bits_counter += unsigned(*pcurr - left + 1u) & is_set;

    unsigned prev_gap = *pcurr++;
    for (is_set ^= ~0u; right > *pcurr; is_set ^= ~0u)
    {
        bits_counter += (*pcurr - prev_gap) & is_set;
        if (pcurr == pend) 
            return bits_counter;
        prev_gap = *pcurr++;
    }
    bits_counter += unsigned(right - prev_gap) & is_set;
    return bits_counter;
}

/*!
    \brief GAP block find position for the rank

    \param block - bit block buffer pointer
    \param rank - rank to find (must be > 0)
    \param nbit_from - start bit position in block
    \param nbit_pos - found position
 
    \return 0 if position with rank was found, or
              the remaining rank (rank - population count)

    @ingroup gapfunc
*/
template<typename T>
bm::id_t gap_find_rank(const T* const block,
                       bm::id_t rank,
                       bm::gap_word_t  nbit_from,
                       unsigned&       nbit_pos)
{
    BM_ASSERT(block);
    BM_ASSERT(rank);

    const T* pcurr = block;
    const T* pend = pcurr + (*pcurr >> 3);

    unsigned bits_counter = 0;
    unsigned is_set;
    unsigned start_pos = bm::gap_bfind(block, nbit_from, &is_set);
    is_set = ~(is_set - 1u); // 0xFFF.. if true (mask for branchless code)

    pcurr = block + start_pos;
    bits_counter += unsigned(*pcurr - nbit_from + 1u) & is_set;
    if (bits_counter >= rank) // found!
    {
        nbit_pos = nbit_from + rank - 1u;
        return 0;
    }
    rank -= bits_counter;
    unsigned prev_gap = *pcurr++;
    for (is_set ^= ~0u; pcurr <= pend; is_set ^= ~0u)
    {
        bits_counter = (*pcurr - prev_gap) & is_set;
        if (bits_counter >= rank) // found!
        {
            nbit_pos = prev_gap + rank;
            return 0;
        }
        rank -= bits_counter;
        prev_gap = *pcurr++;
    } // for
    
    return rank;
}
                       


/*!
    \brief Counts 1 bits in GAP buffer in the closed [0, right] range.
    \param buf - GAP buffer pointer.
    \param right- rightmost bit index
    \return Number of non-zero bits.
    @ingroup gapfunc
*/
template<typename T>
unsigned gap_bit_count_to(const T* const buf, T right)
{
    const T* pcurr = buf;
    const T* pend = pcurr + (*pcurr >> 3);

    unsigned bits_counter = 0;
    unsigned is_set = ~((unsigned(*buf) & 1u) - 1u); // 0xFFF.. if true (mask for branchless code)
    BM_ASSERT(is_set == 0u || is_set == ~0u);
    pcurr = buf + 1;

    if (right <= *pcurr) // we are in the target block right now
    {
        bits_counter = (right + 1u) & is_set; // & is_set == if (is_set)
        return bits_counter;
    }
    bits_counter += (*pcurr + 1u) & is_set;

    unsigned prev_gap = *pcurr++;
    for (is_set ^= ~0u; right > *pcurr; is_set ^= ~0u)
    {
        bits_counter += (*pcurr - prev_gap) & is_set;
        if (pcurr == pend)
            return bits_counter;
        prev_gap = *pcurr++;
    }
    bits_counter += (right - prev_gap) & is_set;
    return bits_counter;
}


/*! 
    D-GAP block for_each algorithm
    
    D-Gap Functor is called for each element but last one.
    
   \param gap_buf - GAP buffer 
   \param func - functor object
    
*/
template<class T, class Func> 
void for_each_dgap(const T* gap_buf, Func& func)
{
    const T* pcurr = gap_buf;
    const T* pend = pcurr + (*pcurr >> 3);
    ++pcurr;
    
    T prev = *pcurr;
    func((T)(prev + 1)); // first element incremented to avoid 0
    ++pcurr;
    do
    {
        func((T)(*pcurr - prev)); // all others are [N] - [N-1]
        prev = *pcurr;
    } while (++pcurr < pend);
}

/** d-Gap copy functor
    @internal
*/
template<typename T> struct d_copy_func
{
    d_copy_func(T* dg_buf) : dgap_buf_(dg_buf) {}
    void operator()(T dgap) { *dgap_buf_++ = dgap; }

    T* dgap_buf_;
};

/*! 
   \brief Convert GAP buffer into D-GAP buffer
   
   Delta GAP representation is DGAP[N] = GAP[N] - GAP[N-1]    
   
   \param gap_buf - GAP buffer 
   \param dgap_buf - Delta-GAP buffer
   \param copy_head - flag to copy GAP header
   
   \internal
   
   @ingroup gapfunc
*/
template<typename T>
T* gap_2_dgap(const T* gap_buf, T* dgap_buf, bool copy_head=true)
{
    if (copy_head) // copy GAP header
    {
        *dgap_buf++ = *gap_buf;
    }

    d_copy_func<T> copy_func(dgap_buf);
    for_each_dgap<T, d_copy_func<T> >(gap_buf, copy_func);
    return copy_func.dgap_buf_;
}

/*! 
   \brief Convert D-GAP buffer into GAP buffer
   
   GAP representation is GAP[N] = DGAP[N] + DGAP[N-1]    
   
   \param dgap_buf - Delta-GAP buffer
   \param gap_header - GAP header word
   \param gap_buf  - GAP buffer

   \internal
   @ingroup gapfunc
*/
template<typename T>
void dgap_2_gap(const T* dgap_buf, T* gap_buf, T gap_header=0)
{
    const T* pcurr = dgap_buf;
    unsigned len;    
    if (!gap_header) // GAP header is already part of the stream
    {
        len = *pcurr >> 3;
        *gap_buf++ = *pcurr++; // copy GAP header
    }
    else // GAP header passed as a parameter
    {
        len = gap_header >> 3;
        *gap_buf++ = gap_header; // assign GAP header
    }    
    --len; // last element is actually not encoded
    const T* pend = pcurr + len;

    *gap_buf = *pcurr++; // copy first element
    if (*gap_buf == 0) 
        *gap_buf = 65535; // fix +1 overflow
    else
        *gap_buf = T(*gap_buf - 1);
    
    for (++gap_buf; pcurr < pend; ++pcurr)
    {
        T prev = *(gap_buf-1); // don't remove temp(undef expression!)           
        *gap_buf++ = T(*pcurr + prev);
    }
    *gap_buf = 65535; // add missing last element  
}


/*! 
   \brief Lexicographical comparison of GAP buffers.
   \param buf1 - First GAP buffer pointer.
   \param buf2 - Second GAP buffer pointer.
   \return  <0 - less, =0 - equal,  >0 - greater.

   @ingroup gapfunc
*/
template<typename T> int gapcmp(const T* buf1, const T* buf2)
{
    const T* pcurr1 = buf1;
    const T* pend1 = pcurr1 + (*pcurr1 >> 3);
    unsigned bitval1 = *buf1 & 1;
    ++pcurr1;

    const T* pcurr2 = buf2;
    unsigned bitval2 = *buf2 & 1;
    ++pcurr2;

    while (pcurr1 <= pend1)
    {
        if (*pcurr1 == *pcurr2)
        {
            if (bitval1 != bitval2)
            {
                return (bitval1) ? 1 : -1;
            }
        }
        else
        {
            if (bitval1 == bitval2)
            {
                if (bitval1)
                {
                    return (*pcurr1 < *pcurr2) ? -1 : 1;
                }
                else
                {
                    return (*pcurr1 < *pcurr2) ? 1 : -1;
                }
            }
            else
            {
                return (bitval1) ? 1 : -1;
            }
        }

        ++pcurr1; ++pcurr2;

        bitval1 ^= 1;
        bitval2 ^= 1;
    }

    return 0;
}


/*!
   \brief Abstract operation for GAP buffers. 
          Receives functor F as a template argument
   \param dest - destination memory buffer.
   \param vect1 - operand 1 GAP encoded buffer.
   \param vect1_mask - XOR mask for starting bitflag for vector1 
   can be 0 or 1 (1 inverts the vector)
   \param vect2 - operand 2 GAP encoded buffer.
   \param vect2_mask - same as vect1_mask
   \param f - operation functor.
   \param dlen - destination length after the operation

   \note Internal function.
   @internal

   @ingroup gapfunc
*/
template<typename T, class F> 
void gap_buff_op(T*         BMRESTRICT dest, 
                 const T*   BMRESTRICT vect1,
                 unsigned   vect1_mask, 
                 const T*   BMRESTRICT vect2,
                 unsigned   vect2_mask, 
                 F&         f,
                 unsigned&  dlen)
{
    const T*  cur1 = vect1;
    const T*  cur2 = vect2;

    T bitval1 = (T)((*cur1++ & 1) ^ vect1_mask);
    T bitval2 = (T)((*cur2++ & 1) ^ vect2_mask);
    
    T bitval = (T) f(bitval1, bitval2);
    T bitval_prev = bitval;

    T* res = dest;
    *res = bitval;
    ++res;

    while (1)
    {
        bitval = (T) f(bitval1, bitval2);

        // Check if GAP value changes and we need to 
        // start the next one.
        if (bitval != bitval_prev)
        {
            ++res;
            bitval_prev = bitval;
        }

        if (*cur1 < *cur2)
        {
            *res = *cur1;
            ++cur1;
            bitval1 ^= 1;
        }
        else // >=
        {
            *res = *cur2;
            if (*cur2 < *cur1)
            {
                bitval2 ^= 1;                
            }
            else  // equal
            {
                if (*cur2 == (bm::gap_max_bits - 1))
                {
                    break;
                }

                ++cur1;
                bitval1 ^= 1;
                bitval2 ^= 1;
            }
            ++cur2;
        }

    } // while

    dlen = (unsigned)(res - dest);
    *dest = (T)((*dest & 7) + (dlen << 3));
}

/*!
   \brief Abstract distance test operation for GAP buffers. 
          Receives functor F as a template argument
   \param vect1 - operand 1 GAP encoded buffer.
   \param vect1_mask - XOR mask for starting bitflag for vector1 
                       can be 0 or 1 (1 inverts the vector)
   \param vect2 - operand 2 GAP encoded buffer.
   \param vect2_mask - same as vect1_mask
   \param f - operation functor.
   \note Internal function.
   \return non zero value if operation result returns any 1 bit 

   @ingroup gapfunc
*/
template<typename T, class F> 
unsigned gap_buff_any_op(const T*   BMRESTRICT vect1,
                         unsigned              vect1_mask, 
                         const T*   BMRESTRICT vect2,
                         unsigned              vect2_mask, 
                         F                     f)
{
    const T*  cur1 = vect1;
    const T*  cur2 = vect2;

    unsigned bitval1 = (*cur1++ & 1) ^ vect1_mask;
    unsigned bitval2 = (*cur2++ & 1) ^ vect2_mask;
    
    unsigned bitval = f(bitval1, bitval2);
    if (bitval)
        return bitval;
    unsigned bitval_prev = bitval;

    while (1)
    {
        bitval = f(bitval1, bitval2);
        if (bitval)
            return bitval;

        if (bitval != bitval_prev)
            bitval_prev = bitval;

        if (*cur1 < *cur2)
        {
            ++cur1;
            bitval1 ^= 1;
        }
        else // >=
        {
            if (*cur2 < *cur1)
            {
                bitval2 ^= 1;                
            }
            else  // equal
            {
                if (*cur2 == (bm::gap_max_bits - 1))
                {
                    break;
                }

                ++cur1;
                bitval1 ^= 1;
                bitval2 ^= 1;
            }
            ++cur2;
        }

    } // while

    return 0;
}



/*!
   \brief Abstract distance(similarity) operation for GAP buffers. 
          Receives functor F as a template argument
   \param vect1 - operand 1 GAP encoded buffer.
   \param vect2 - operand 2 GAP encoded buffer.
   \param f - operation functor.
   \note Internal function.

   @ingroup gapfunc
*/
template<typename T, class F> 
unsigned gap_buff_count_op(const T*  vect1, const T*  vect2, F f)
{
    unsigned count;// = 0;
    const T* cur1 = vect1;
    const T* cur2 = vect2;

    unsigned bitval1 = (*cur1++ & 1);
    unsigned bitval2 = (*cur2++ & 1);
    unsigned bitval = count = f(bitval1, bitval2);
    unsigned bitval_prev = bitval;

    //if (bitval) ++count;
    
    T res, res_prev;
    res = res_prev = 0;

    while (1)
    {
        bitval = f(bitval1, bitval2);

        // Check if GAP value changes and we need to 
        // start the next one.
        if (bitval != bitval_prev)
        {
            bitval_prev = bitval;
            res_prev = res;
        }

        if (*cur1 < *cur2)
        {
            res = *cur1;
            if (bitval)
            {
                count += res - res_prev; 
                res_prev = res;
            }
            ++cur1;
            bitval1 ^= 1;
        }
        else // >=
        {
            res = *cur2;
            if (bitval)
            {
                count += res - res_prev; 
                res_prev = res;
            }
            if (*cur2 < *cur1)
            {
                bitval2 ^= 1;                
            }
            else  // equal
            {
                if (*cur2 == (bm::gap_max_bits - 1))
                {
                    break;
                }

                ++cur1;
                bitval1 ^= 1;
                bitval2 ^= 1;
            }
            ++cur2;
        }

    } // while

    return count;
}



/*!
   \brief Sets or clears bit in the GAP buffer.

   \param val - new bit value
   \param buf - GAP buffer.
   \param pos - Index of bit to set.
   \param is_set - (OUT) flag if bit was actually set.

   \return New GAP buffer length. 

   @ingroup gapfunc
*/
template<typename T> unsigned gap_set_value(unsigned val, 
                                            T* BMRESTRICT buf, 
                                            unsigned pos, 
                                            unsigned* BMRESTRICT is_set)
{
    BM_ASSERT(pos < bm::gap_max_bits);
    unsigned curr = gap_bfind(buf, pos, is_set);

    T end = (T)(*buf >> 3);
    if (*is_set == val)
    {
        *is_set = 0;
        return end;
    }
    *is_set = 1;

    T* pcurr = buf + curr;
    T* pprev = pcurr - 1;
    T* pend = buf + end;

    // Special case, first bit GAP operation. There is no platform beside it.
    // initial flag must be inverted.
    if (pos == 0)
    {
        *buf ^= 1;
        if ( buf[1] ) // We need to insert a 1 bit platform here.
        {
            ::memmove(&buf[2], &buf[1], (end - 1) * sizeof(gap_word_t));
            buf[1] = 0;
            ++end;
        }
        else // Only 1 bit in the GAP. We need to delete the first GAP.
        {
            pprev = buf + 1;
            pcurr = pprev + 1;
            do
            {
                *pprev++ = *pcurr++;
            } while (pcurr < pend);
            --end;
        }
    }
    else if (curr > 1 && ((unsigned)(*pprev))+1 == pos) // Left border bit
    {
       ++(*pprev);
       if (*pprev == *pcurr)  // Curr. GAP to be merged with prev.GAP.
       {
            --end;
            if (pcurr != pend) // GAP merge: 2 GAPS to be deleted 
            {
                --end;
                ++pcurr;
                do
                {
                    *pprev++ = *pcurr++;
                } while (pcurr < pend);
            }
       }    
    }
    else if (*pcurr == pos) // Rightmost bit in the GAP. Border goes left.
    {
        --(*pcurr);       
        if (pcurr == pend)
        {
           ++end;
        }
    }
    else  // Worst case we need to split current block.
    {
        ::memmove(pcurr+2, pcurr,(end - curr + 1)*sizeof(T));
        *pcurr++ = (T)(pos - 1);
        *pcurr = (T)pos;
        end = (T)(end + 2);
    }

    // Set correct length word.
    *buf = (T)((*buf & 7) + (end << 3));

    buf[end] = bm::gap_max_bits - 1;
    return end;
}

/*!
   \brief Add new value to the end of GAP buffer.

   \param buf - GAP buffer.
   \param pos - Index of bit to set.

   \return New GAP buffer length. 

   @ingroup gapfunc
*/
template<typename T> 
unsigned gap_add_value(T* buf, unsigned pos)
{
    BM_ASSERT(pos < bm::gap_max_bits);

    T end = (T)(*buf >> 3);
    T curr = end;
    T* pcurr = buf + end;
    T* pend  = pcurr;
    T* pprev = pcurr - 1;

    // Special case, first bit GAP operation. There is no platform beside it.
    // initial flag must be inverted.
    if (pos == 0)
    {
        *buf ^= 1;
        if ( buf[1] ) // We need to insert a 1 bit platform here.
        {
            ::memmove(&buf[2], &buf[1], (end - 1) * sizeof(gap_word_t));
            buf[1] = 0;
            ++end;
        }
        else // Only 1 bit in the GAP. We need to delete the first GAP.
        {
            pprev = buf + 1;
            pcurr = pprev + 1;
            do
            {
                *pprev++ = *pcurr++;
            } while (pcurr < pend);
            --end;
        }
    }
    else if (((unsigned)(*pprev))+1 == pos && (curr > 1) ) // Left border bit
    {
       ++(*pprev);
       if (*pprev == *pcurr)  // Curr. GAP to be merged with prev.GAP.
       {
            --end;
            if (pcurr != pend) // GAP merge: 2 GAPS to be deleted 
            {
                // TODO: should never get here...
                --end;
                ++pcurr;
                do
                {
                    *pprev++ = *pcurr++;
                } while (pcurr < pend);
            }
       } 
    }
    else if (*pcurr == pos) // Rightmost bit in the GAP. Border goes left.
    {
        --(*pcurr);       
        if (pcurr == pend)
        {
           ++end;
        }
    }
    else  // Worst case we need to split current block.
    {
        *pcurr++ = (T)(pos - 1);
        *pcurr = (T)pos;
        end = (T)(end+2);
    }

    // Set correct length word.
    *buf = (T)((*buf & 7) + (end << 3));

    buf[end] = bm::gap_max_bits - 1;
    return end;
}

/*!
   \brief Convert array to GAP buffer.

   \param buf - GAP buffer.
   \param arr - array of values to set
   \param len - length of the array

   \return New GAP buffer length. 

   @ingroup gapfunc
*/

template<typename T> 
unsigned gap_set_array(T* buf, const T* arr, unsigned len)
{
    *buf = (T)((*buf & 6u) + (1u << 3)); // gap header setup

    T* pcurr = buf + 1;

    unsigned i = 0;
    T curr = arr[i];
    if (curr != 0) // need to add the first gap: (0 to arr[0]-1)
    {
        *pcurr = (T)(curr - 1);
        ++pcurr;
    }
    else
    {
        ++(*buf); // GAP starts with 1
    }
    T prev = curr; 
    T acc = prev;

    for (i = 1; i < len; ++i)
    {
        curr = arr[i];
        if (curr == prev + 1)
        {
            ++acc;
            prev = curr;
        }
        else
        {
            *pcurr++ = acc;
            acc = curr;
            *pcurr++ = (T)(curr-1);
        }
        prev = curr;
    }
    *pcurr = acc;
    if (acc != bm::gap_max_bits - 1)
    {
        ++pcurr;
        *pcurr = bm::gap_max_bits - 1;
    }

    unsigned end = unsigned(pcurr - buf);

    *buf = (T)((*buf & 7) + (end << 3));
    return end+1;
}


//------------------------------------------------------------------------

/**
    \brief Compute number of GAPs in bit-array
    \param arr - array of BITs
    \param len - array length

    @ingroup gapfunc
*/
template<typename T> 
unsigned bit_array_compute_gaps(const T* arr, 
                                unsigned len)
{
    unsigned gap_count = 1;
    T prev = arr[0];
    if (prev > 0)
        ++gap_count;
    for (unsigned i = 1; i < len; ++i)
    {
        T curr = arr[i];
        if (curr != prev + 1)
        {
            gap_count += 2;
        }
        prev = curr;
    }
    return gap_count;
}


//------------------------------------------------------------------------

/**
    \brief Searches for the next 1 bit in the GAP block
    \param buf - GAP buffer
    \param nbit - bit index to start checking from.
    \param prev - returns previously checked value
 
    \return 0 if not found

    @ingroup gapfunc
*/
template<typename T>
unsigned gap_find_in_block(const T* buf,
                           unsigned nbit,
                           bm::id_t* prev)
{
    BM_ASSERT(nbit < bm::gap_max_bits);

    unsigned bitval;
    unsigned gap_idx = bm::gap_bfind(buf, nbit, &bitval);

    if (bitval) // positive block.
    {
       return 1u;
    }

    unsigned val = buf[gap_idx] + 1;
    *prev += val - nbit;
 
    return (val != bm::gap_max_bits);  // no bug here.
}

/*! 
    \brief Set 1 bit in a block
    
    @ingroup bitfunc
*/
BMFORCEINLINE
void set_bit(unsigned* dest, unsigned  bitpos)
{
    unsigned nbit  = unsigned(bitpos & bm::set_block_mask); 
    unsigned nword = unsigned(nbit >> bm::set_word_shift); 
    nbit &= bm::set_word_mask;
    dest[nword] |= unsigned(1u << nbit);
}

/*! 
    \brief Test 1 bit in a block
    
    @ingroup bitfunc
*/
BMFORCEINLINE
unsigned test_bit(const unsigned* block, unsigned  bitpos)
{
    unsigned nbit  = unsigned(bitpos & bm::set_block_mask); 
    unsigned nword = unsigned(nbit >> bm::set_word_shift); 
    nbit &= bm::set_word_mask;
    return (block[nword] >> nbit) & 1u;
}


/*! 
   \brief Sets bits to 1 in the bitblock.
   \param dest - Bitset buffer.
   \param bitpos - Offset of the start bit.
   \param bitcount - number of bits to set.

   @ingroup bitfunc
*/  
inline
void or_bit_block(unsigned* dest, unsigned bitpos, unsigned bitcount)
{
    unsigned nbit  = unsigned(bitpos & bm::set_block_mask); 
    unsigned nword = unsigned(nbit >> bm::set_word_shift); 
    nbit &= bm::set_word_mask;

    bm::word_t* word = dest + nword;

    if (bitcount == 1)  // special case (only 1 bit to set)
    {
        *word |= bitcount << nbit;
        return;
    }

    if (nbit) // starting position is not aligned
    {
        unsigned right_margin = nbit + bitcount;

        // here we checking if we setting bits only in the current
        // word. Example: 00111000000000000000000000000000 (32 bits word)

        if (right_margin < 32) 
        {
            unsigned mask = 
                block_set_table<true>::_right[nbit] & 
                block_set_table<true>::_left[right_margin-1];
            *word |= mask;
            return;
        }
        *word |= block_set_table<true>::_right[nbit];
        bitcount -= 32 - nbit;
        ++word;
    }
    // now word aligned, use unrolled loops
    //
    for ( ;bitcount >= 64; bitcount-=64, word+=2)
        word[0] = word[1] = ~0u;
    if (bitcount >= 32)
    {
        *word++ = ~0u; bitcount -= 32;
    }
    BM_ASSERT(bitcount < 32);

    if (bitcount) 
        *word |= block_set_table<true>::_left[bitcount-1];
}


/*! 
   \brief SUB (AND NOT) bit interval to 1 in the bitblock.
   \param dest - Bitset buffer.
   \param bitpos - Offset of the start bit.
   \param bitcount - number of bits to set.

   @ingroup bitfunc
*/  
inline void sub_bit_block(unsigned* dest, 
                          unsigned bitpos, 
                          unsigned bitcount)
{
    unsigned nbit  = unsigned(bitpos & bm::set_block_mask); 
    unsigned nword = unsigned(nbit >> bm::set_word_shift); 
    nbit &= bm::set_word_mask;

    bm::word_t* word = dest + nword;

    if (bitcount == 1)  // special case (only 1 bit to set)
    {
        *word &= ~unsigned(1 << nbit);
        return;
    }

    if (nbit) // starting position is not aligned
    {
        unsigned right_margin = nbit + bitcount;

        // here we checking if we setting bits only in the current
        // word. Example: 00111000000000000000000000000000 (32 bits word)

        if (right_margin < 32) 
        {
            unsigned mask = 
                block_set_table<true>::_right[nbit] & 
                block_set_table<true>::_left[right_margin-1];
            *word &= ~mask;
            return;
        }
        *word &= ~block_set_table<true>::_right[nbit];
        bitcount -= 32 - nbit;
        ++word;
    }
    // now word aligned, use unrolled loops
    //
    for ( ;bitcount >= 64; bitcount-=64, word+=2)
        word[0] = word[1] = 0u;
    if (bitcount >= 32)
    {
        *word++ = 0u; bitcount -= 32;
    }
    if (bitcount)
    {
        *word &= ~block_set_table<true>::_left[bitcount-1];
    }
}


/*! 
   \brief XOR bit interval to 1 in the bitblock.
   \param dest - Bitset buffer.
   \param bitpos - Offset of the start bit.
   \param bitcount - number of bits to set.

   @ingroup bitfunc
*/  
inline void xor_bit_block(unsigned* dest, 
                          unsigned bitpos, 
                          unsigned bitcount)
{
    unsigned nbit  = unsigned(bitpos & bm::set_block_mask); 
    unsigned nword = unsigned(nbit >> bm::set_word_shift); 
    nbit &= bm::set_word_mask;

    bm::word_t* word = dest + nword;

    if (bitcount == 1)  // special case (only 1 bit to set)
    {
        *word ^= unsigned(1 << nbit);
        return;                             
    }

    if (nbit) // starting position is not aligned
    {
        unsigned right_margin = nbit + bitcount;

        // here we checking if we setting bits only in the current
        // word. Example: 00111000000000000000000000000000 (32 bits word)

        if (right_margin < 32) 
        {
            unsigned mask = 
                block_set_table<true>::_right[nbit] & 
                block_set_table<true>::_left[right_margin-1];
            *word ^= mask;
            return;
        }
        *word ^= block_set_table<true>::_right[nbit];
        bitcount -= 32 - nbit;
        ++word;
    }
    for ( ;bitcount >= 64; bitcount-=64, word+=2)
    {
        word[0] ^= ~0u; word[1] ^= ~0u;
    }
    if (bitcount >= 32)
    {
        *word++ ^= ~0u; bitcount -= 32;
    }
    if (bitcount)
        *word ^= block_set_table<true>::_left[bitcount-1];
}


/*!
   \brief SUB (AND NOT) GAP block to bitblock.
   \param dest - bitblock buffer pointer.
   \param pcurr  - GAP buffer pointer.

   @ingroup gapfunc
*/
template<typename T> 
void gap_sub_to_bitset(unsigned* dest, const T*  pcurr)
{
    BM_ASSERT(dest && pcurr);
    
    const T* pend = pcurr + (*pcurr >> 3);
    if (*pcurr & 1)  // Starts with 1
    {
        bm::sub_bit_block(dest, 0, 1 + pcurr[1]);
        ++pcurr;
    }
    for (pcurr += 2; pcurr <= pend; pcurr += 2)
    {
        BM_ASSERT(*pcurr > pcurr[-1]);
        bm::sub_bit_block(dest, 1 + pcurr[-1], *pcurr - pcurr[-1]);
    }
}



/*!
   \brief XOR GAP block to bitblock.
   \param dest - bitblock buffer pointer.
   \param pcurr  - GAP buffer pointer.

   @ingroup gapfunc
*/
template<typename T> 
void gap_xor_to_bitset(unsigned* dest, const T*  pcurr)
{
    BM_ASSERT(dest && pcurr);

    const T* pend = pcurr + (*pcurr >> 3);
    if (*pcurr & 1)  // Starts with 1
    {
        bm::xor_bit_block(dest, 0, 1 + pcurr[1]);
        ++pcurr;
    }
    for (pcurr += 2; pcurr <= pend; pcurr += 2)
    {
        BM_ASSERT(*pcurr > pcurr[-1]);
        bm::xor_bit_block(dest, 1 + pcurr[-1], *pcurr - pcurr[-1]);
    }
}


/*!
   \brief Adds(OR) GAP block to bitblock.
   \param dest - bitblock buffer pointer.
   \param pcurr  - GAP buffer pointer.

   @ingroup gapfunc
*/
template<typename T>
void gap_add_to_bitset(unsigned* dest, const T*  pcurr)
{
    BM_ASSERT(dest && pcurr);
    
    const T* pend = pcurr + (*pcurr >> 3);
    if (*pcurr & 1)  // Starts with 1
    {
        bm::or_bit_block(dest, 0, 1 + pcurr[1]);
        ++pcurr;
    }
    for (pcurr += 2; pcurr <= pend; pcurr += 2)
    {
        BM_ASSERT(*pcurr > pcurr[-1]);
        bm::or_bit_block(dest, 1 + pcurr[-1], *pcurr - pcurr[-1]);
    }
}

/*!
   \brief ANDs GAP block to bitblock.
   \param dest - bitblock buffer pointer.
   \param pcurr  - GAP buffer pointer.

   @ingroup gapfunc
*/
template<typename T> 
void gap_and_to_bitset(unsigned* dest, const T*  pcurr)
{
    BM_ASSERT(dest && pcurr);
    
    const T* pend = pcurr + (*pcurr >> 3);
    if (!(*pcurr & 1) )  // Starts with 0
    {
        bm::sub_bit_block(dest, 0, pcurr[1] + 1); // (not AND) - SUB [0] gaps
        ++pcurr;
    }
    for (pcurr += 2; pcurr <= pend; pcurr += 2) // now we are in GAP "0" again
    {
        BM_ASSERT(*pcurr > *(pcurr-1));
        bm::sub_bit_block(dest, pcurr[-1]+1, *pcurr - pcurr[-1]);
    }
}


/*!
   \brief Compute bitcount of bit block AND masked by GAP block
   \param block - bitblock buffer pointer
   \param pcurr  - GAP buffer pointer
   \return bitcount - cardinality of the AND product

   @ingroup gapfunc
*/
template<typename T> 
bm::id_t gap_bitset_and_count(const unsigned* block, const T*  pcurr)
{
    BM_ASSERT(block);
    const T* pend = pcurr + (*pcurr >> 3);
    bm::id_t count = 0;
    if (*pcurr & 1)  // Starts with 1
    {
        count += bm::bit_block_calc_count_range(block, 0, pcurr[1]);
        ++pcurr;
    }
    for (pcurr +=2 ;pcurr <= pend; pcurr += 2)
    {
        count += bm::bit_block_calc_count_range(block, pcurr[-1]+1, *pcurr);
    }
    return count;
}


/*!
   \brief Bitcount test of bit block AND masked by GAP block.
   \param block - bitblock buffer pointer
   \param pcurr  - GAP buffer pointer
   \return non-zero value if AND produces any result

   @ingroup gapfunc
*/
template<typename T> 
bm::id_t gap_bitset_and_any(const unsigned* block, const T* pcurr)
{
    BM_ASSERT(block);

    const T* pend = pcurr + (*pcurr >> 3);
    bm::id_t count = 0;
    if (*pcurr & 1)  // Starts with 1
    {
        count = bm::bit_block_any_range(block, 0, pcurr[1]);
        ++pcurr;
    }
    for (pcurr +=2 ;!count && pcurr <= pend; pcurr += 2)
    {
        count = bm::bit_block_any_range(block, pcurr[-1]+1, *pcurr);
    }
    return count;
}



/*!
   \brief Compute bitcount of bit block SUB masked by GAP block.
   \param block - bitblock buffer pointer.
   \param buf  - GAP buffer pointer.
   \return bit-count result of AND NOT operation

   @ingroup gapfunc
*/
template<typename T> 
bm::id_t gap_bitset_sub_count(const unsigned* block, const T*  buf)
{
    BM_ASSERT(block);

    const T* pcurr = buf;
    const T* pend = pcurr + (*pcurr >> 3);
    ++pcurr;

    bm::id_t count = 0;

    if (!(*buf & 1))  // Starts with 0
    {
        count += bit_block_calc_count_range(block, 0, *pcurr);
        ++pcurr;
    }
    ++pcurr; // now we are in GAP "0" again

    for (;pcurr <= pend; pcurr+=2)
    {
        count += bm::bit_block_calc_count_range(block, *(pcurr-1)+1, *pcurr);
    }
    return count;
}


/*!
   \brief Compute bitcount test of bit block SUB masked by GAP block
   \param block - bitblock buffer pointer
   \param buf  - GAP buffer pointer
   \return non-zero value if AND NOT produces any 1 bits

   @ingroup gapfunc
*/
template<typename T> 
bm::id_t gap_bitset_sub_any(const unsigned* block, const T*  buf)
{
    BM_ASSERT(block);

    const T* pcurr = buf;
    const T* pend = pcurr + (*pcurr >> 3);
    ++pcurr;

    bm::id_t count = 0;

    if (!(*buf & 1))  // Starts with 0
    {
        count += bit_block_any_range(block, 0, *pcurr);
        if (count)
            return count;
        ++pcurr;
    }
    ++pcurr; // now we are in GAP "0" again

    for (; !count && pcurr <= pend; pcurr+=2)
    {
        count += bm::bit_block_any_range(block, *(pcurr-1)+1, *pcurr);
    }
    return count;
}



/*!
   \brief Compute bitcount of bit block XOR masked by GAP block
   \param block - bitblock buffer pointer
   \param buf  - GAP buffer pointer
   \return bit count value of XOR operation

   @ingroup gapfunc
*/
template<typename T> 
bm::id_t gap_bitset_xor_count(const unsigned* block, const T*  buf)
{
    BM_ASSERT(block);

    const T* pcurr = buf;
    const T* pend = pcurr + (*pcurr >> 3);
    ++pcurr;

    unsigned bitval = *buf & 1;
    
    bm::id_t count = bm::bit_block_calc_count_range(block, 0, *pcurr);
    if (bitval)
    {
        count = *pcurr + 1 - count;
    }
    
    for (bitval^=1, ++pcurr; pcurr <= pend; bitval^=1, ++pcurr)
    {
        T prev = (T)(*(pcurr-1)+1);
        bm::id_t c = bit_block_calc_count_range(block, prev, *pcurr);
        
        if (bitval) // 1 gap; means Result = Total_Bits - BitCount;
            c = (*pcurr - prev + 1) - c;
        count += c;
    }
    return count;
}

/*!
   \brief Compute bitcount test of bit block XOR masked by GAP block.
   \param block - bitblock buffer pointer
   \param buf  - GAP buffer pointer
   \return non-zero value if XOR returns anything

   @ingroup gapfunc
*/
template<typename T> 
bm::id_t gap_bitset_xor_any(const unsigned* block, const T*  buf)
{
    BM_ASSERT(block);

    const T* pcurr = buf;
    const T* pend = pcurr + (*pcurr >> 3);
    ++pcurr;

    unsigned bitval = *buf & 1;
    
    bm::id_t count = bit_block_any_range(block, 0, *pcurr);
    if (bitval)
        count = *pcurr + 1 - count;
    
    for (bitval^=1, ++pcurr; !count && pcurr <= pend; bitval^=1, ++pcurr)
    {
        T prev = (T)(*(pcurr-1)+1);
        bm::id_t c = bit_block_any_range(block, prev, *pcurr);
        
        if (bitval) // 1 gap; means Result = Total_Bits - BitCount;
            c = (*pcurr - prev + 1) - c;
        count += c;
    }
    return count;
}



/*!
   \brief Compute bitcount of bit block OR masked by GAP block.
   \param block - bitblock buffer pointer.
   \param buf  - GAP buffer pointer.
   \return bit count of OR

   @ingroup gapfunc
*/
template<typename T> 
bm::id_t gap_bitset_or_count(const unsigned* block, const T*  buf)
{
    BM_ASSERT(block);

    const T* pcurr = buf;
    const T* pend = pcurr + (*pcurr >> 3);
    ++pcurr;

    unsigned bitval = *buf & 1;
    
    bm::id_t count = bitval ? *pcurr + 1
                            : bm::bit_block_calc_count_range(block, 0, *pcurr);
    for (bitval^=1, ++pcurr; pcurr <= pend; bitval^=1, ++pcurr)
    {
        T prev = (T)(*(pcurr-1)+1);
        bm::id_t c =
            bitval ? (*pcurr - prev + 1)
                   : bm::bit_block_calc_count_range(block, prev, *pcurr);
        count += c;
    }
    return count;
}

/*!
   \brief Compute bitcount test of bit block OR masked by GAP block
   \param block - bitblock buffer pointer
   \param buf  - GAP buffer pointer
   \return non zero value if union (OR) returns anything

   @ingroup gapfunc
*/
template<typename T> 
bm::id_t gap_bitset_or_any(const unsigned* block, const T*  buf)
{
    bool b = !bm::gap_is_all_zero(buf) ||
             !bm::bit_is_all_zero(block);
    return b;
}



/*!
   \brief Bitblock memset operation. 

   \param dst - destination block.
   \param value - value to set.

   @ingroup bitfunc
*/
inline 
void bit_block_set(bm::word_t* BMRESTRICT dst, bm::word_t value)
{
//#ifdef BMVECTOPT
//    VECT_SET_BLOCK(dst, dst + bm::set_block_size, value);
//#else
    ::memset(dst, int(value), bm::set_block_size * sizeof(bm::word_t));
//#endif
}


/*!
   \brief GAP block to bitblock conversion.
   \param dest - bitblock buffer pointer.
   \param buf  - GAP buffer pointer.

   @ingroup gapfunc
*/
template<typename T> 
void gap_convert_to_bitset(unsigned* dest, const T*  buf)
{
    bm::bit_block_set(dest, 0);
    bm::gap_add_to_bitset(dest, buf);
}


/*!
   \brief GAP block to bitblock conversion.
   \param dest - bitblock buffer pointer.
   \param buf  - GAP buffer pointer.
   \param dest_len - length/size of the destination buffer.

   @ingroup gapfunc
*/
template<typename T> 
void gap_convert_to_bitset(unsigned* dest, const T*  buf,  unsigned  dest_len)
{
   ::memset(dest, 0, dest_len * sizeof(unsigned));
    bm::gap_add_to_bitset(dest, buf);
}


/*!
   \brief Smart GAP block to bitblock conversion.

    Checks if GAP block is ALL-ZERO or ALL-ON. In those cases returns 
    pointer on special static bitblocks.

   \param dest - bitblock buffer pointer.
   \param buf  - GAP buffer pointer.
   \param set_max - max possible bitset length

   @ingroup gapfunc
*/
template<typename T> 
unsigned* gap_convert_to_bitset_smart(unsigned* dest,
                                      const T* buf,
                                      id_t set_max)
{
    if (buf[1] == set_max - 1)
        return (buf[0] & 1) ? FULL_BLOCK_REAL_ADDR : 0;

    bm::gap_convert_to_bitset(dest, buf);
    return dest;
}


/*!
   \brief Calculates sum of all words in GAP block. (For debugging purposes)
   \note For debugging and testing ONLY.
   \param buf - GAP buffer pointer.
   \return Sum of all words.

   @ingroup gapfunc
   @internal
*/
template<typename T> unsigned gap_control_sum(const T* buf)
{
    unsigned end = *buf >> 3;

    const T* pcurr = buf;
    const T* pend = pcurr + (*pcurr >> 3);
    ++pcurr;

    if (*buf & 1)  // Starts with 1
    {
        ++pcurr;
    }
    ++pcurr; // now we are in GAP "1" again

    while (pcurr <= pend)
    {
        BM_ASSERT(*pcurr > *(pcurr-1));
        pcurr += 2;
    }
    return buf[end];
}


/*! 
   \brief Sets all bits to 0 or 1 (GAP)
   \param buf - GAP buffer pointer.
   \param set_max - max possible bitset length
   \param value - value to set

   @ingroup gapfunc
*/
template<class T> void gap_set_all(T* buf, 
                                   unsigned set_max,
                                   unsigned value)
{
    BM_ASSERT(value == 0 || value == 1);
    *buf = (T)((*buf & 6u) + (1u << 3) + value);
    *(++buf) = (T)(set_max - 1);
}


/*!
    \brief Init gap block so it has block in it (can be whole block)
    \param buf  - GAP buffer pointer.
    \param from - one block start
    \param to   - one block end
    \param value - (block value)1 or 0
 
   @ingroup gapfunc
*/
template<class T> 
void gap_init_range_block(T* buf,
                          T  from,
                          T  to,
                          T  value)
                          //unsigned set_max)
{
    BM_ASSERT(value == 0 || value == 1);
    const unsigned set_max = bm::bits_in_block;

    unsigned gap_len;
    if (from == 0)
    {
        if (to == set_max - 1)
        {
            bm::gap_set_all(buf, set_max, value);
        }
        else
        {
            gap_len = 2;
            buf[1] = (T)to;
            buf[2] = (T)(set_max - 1);
            buf[0] = (T)((*buf & 6u) + (gap_len << 3) + value);
        }
        return;
    }
    // from != 0

    value = !value;
    if (to == set_max - 1)
    {
        gap_len = 2;
        buf[1] = (T)(from - 1);
        buf[2] = (T)(set_max - 1);
    }
    else
    {
        gap_len = 3;
        buf[1] = (T) (from - 1);
        buf[2] = (T) to;
        buf[3] = (T)(set_max - 1);
    }
    buf[0] =  (T)((*buf & 6u) + (gap_len << 3) + value);
}


/*! 
   \brief Inverts all bits in the GAP buffer.
   \param buf - GAP buffer pointer.

   @ingroup gapfunc
*/
template<typename T> void gap_invert(T* buf)
{ 
    *buf ^= 1;
}


#ifdef __GNUG__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
#endif

/*!
   \brief Sets GAP block capacity level.
   \param buf - GAP buffer pointer.
   \param level new GAP block capacity level.

   @ingroup gapfunc
*/
template<typename T> 
void set_gap_level(T* buf, int level)
{
    BM_ASSERT(level >= 0);
    BM_ASSERT(unsigned(level) < bm::gap_levels);
    
    *buf = (T)(((level & 3) << 1) | (*buf & 1) | (*buf & ~7));
}
#ifdef __GNUG__
#pragma GCC diagnostic pop
#endif



/*!
   \brief Calculates GAP block capacity level.
   \param len - GAP buffer length.
   \param glevel_len - GAP lengths table
   \return GAP block capacity level. 
            -1 if block does not fit any level.
   @ingroup gapfunc
*/
template<typename T>
inline int gap_calc_level(unsigned len, const T* glevel_len)
{
    if (len <= unsigned(glevel_len[0]-4)) return 0;
    if (len <= unsigned(glevel_len[1]-4)) return 1;
    if (len <= unsigned(glevel_len[2]-4)) return 2;
    if (len <= unsigned(glevel_len[3]-4)) return 3;

    BM_ASSERT(bm::gap_levels == 4);
    return -1;
}

/*! @brief Returns number of free elements in GAP block array. 
    Difference between GAP block capacity on this level and actual GAP length.
    
    @param buf - GAP buffer pointer
    @param glevel_len - GAP lengths table
    
    @return Number of free GAP elements
    @ingroup gapfunc
*/
template<typename T>
inline unsigned gap_free_elements(const T* buf, const T* glevel_len)
{
    unsigned len = gap_length(buf);
    unsigned capacity = gap_capacity(buf, glevel_len);
    return capacity - len;
}


/*! 
   \brief Lexicographical comparison of BIT buffers.
   \param buf1 - First buffer pointer.
   \param buf2 - Second buffer pointer.
   \param len - Buffer length in elements (T).
   \return  <0 - less, =0 - equal,  >0 - greater.

   @ingroup bitfunc 
*/
template<typename T> 
int bitcmp(const T* buf1, const T* buf2, unsigned len)
{
    BM_ASSERT(len);
    const T* pend1 = buf1 + len; 
    do
    {
        T w1 = *buf1++;
        T w2 = *buf2++;
        T diff = w1 ^ w2;
    
        if (diff)
        { 
            return (w1 & diff & -diff) ? 1 : -1;
        }

    } while (buf1 < pend1);

    return 0;
}


/*! 
   \brief Converts bit block to GAP. 
   \param dest - Destinatio GAP buffer.
   \param src - Source bitblock buffer.
   \param bits - Number of bits to convert.
   \param dest_len - length of the dest. buffer.
   \return  New length of GAP block or 0 if conversion failed 
   (insufficicent space).

   @ingroup gapfunc
*/
template<typename T>
unsigned bit_convert_to_gap(T* BMRESTRICT dest,
                            const unsigned* BMRESTRICT src,
                            bm::id_t bits,
                            unsigned dest_len)
{
    T* BMRESTRICT pcurr = dest;
    T* BMRESTRICT end = dest + dest_len; 
    unsigned bitval = (*src) & 1u;
    *pcurr = (T)bitval;

    ++pcurr;
    *pcurr = 0;
    unsigned bit_idx = 0;
    unsigned bitval_next;

    unsigned val = *src;

    do
    {
        // We can fast pace if *src == 0 or *src = 0xffffffff

        while (val == 0 || val == 0xffffffff)
        {
           bitval_next = val ? 1 : 0;
           if (bitval != bitval_next)
           {
               *pcurr++ = (T)(bit_idx-1); 
               BM_ASSERT((pcurr-1) == (dest+1) || *(pcurr-1) > *(pcurr-2));
               if (pcurr >= end)
               {
                   return 0; // OUT of memory
               }
               bitval = bitval_next;
           }
           bit_idx += unsigned(sizeof(*src) * 8);
           if (bit_idx >= bits)
           {
               goto complete;
           }
           ++src;
           val = *src;
        }

        unsigned mask = 1;
        while (mask)
        {
            // Now plain bitshifting. TODO: Optimization wanted.

            bitval_next = val & mask ? 1 : 0;
            if (bitval != bitval_next)
            {
                *pcurr++ = (T)(bit_idx-1);
                BM_ASSERT((pcurr-1) == (dest+1) || *(pcurr-1) > *(pcurr-2));
                bitval = bitval_next;
                if (pcurr >= end)
                    return 0; // OUT of memory
            }
            mask <<= 1;
            ++bit_idx;
        } // while mask

        if (bit_idx >= bits)
            goto complete;

        ++src;
        val = *src;

    } while(1);

complete:
    *pcurr = (T)(bit_idx-1);
    unsigned len = (unsigned)(pcurr - dest);
    *dest = (T)((*dest & 7) + (len << 3));
    return len;
}


/*!
   \brief Iterate gap block as delta-bits with a functor 
   @ingroup gapfunc
*/
template<class T, class F>
void for_each_gap_dbit(const T* buf, F& func)
{
    const T* pcurr = buf;
    const T* pend = pcurr + (*pcurr >> 3);

    ++pcurr;

    unsigned prev = 0;
    unsigned first_inc;

    if (*buf & 1)
    {
        first_inc = 0;
        unsigned to = *pcurr;
        for (unsigned i = 0; i <= to; ++i) 
        {
            func(1);
        }
        prev = to;
        ++pcurr;
    }
    else
    {
        first_inc = 1;
    }
    ++pcurr;  // set GAP to 1

    while (pcurr <= pend)
    {
        unsigned from = *(pcurr-1)+1;
        unsigned to = *pcurr;
        if (first_inc)
        {
            func(from - prev + first_inc);
            first_inc = 0;
        }
        else
        {
            func(from - prev);
        }

        for (unsigned i = from+1; i <= to; ++i) 
        {
            func(1);
        }
        prev = to;
        pcurr += 2; // jump to the next positive GAP
    }
}

/*!
   \brief Convert gap block into array of ints corresponding to 1 bits
   @ingroup gapfunc
*/
template<typename D, typename T>
D gap_convert_to_arr(D* BMRESTRICT       dest, 
                     const T* BMRESTRICT buf,
                     unsigned            dest_len,
                     bool                invert = false)
{
    BMREGISTER const T* BMRESTRICT pcurr = buf;
    BMREGISTER const T* pend = pcurr + (*pcurr >> 3);

    D* BMRESTRICT dest_curr = dest;
    ++pcurr;

    int bitval = (*buf) & 1;
    if (invert) 
        bitval = !bitval; // invert the GAP buffer

    if (bitval)
    {
        if (unsigned(*pcurr + 1) >= dest_len)
            return 0; // insufficient space
        dest_len -= *pcurr;
        T to = *pcurr;
        for (T i = 0; ;++i) 
        {
            *dest_curr++ = i;
            if (i == to) break;
        }
        ++pcurr;
    }
    ++pcurr;  // set GAP to 1

    while (pcurr <= pend)
    {
        unsigned pending = *pcurr - *(pcurr-1);
        if (pending >= dest_len)
            return 0;
        dest_len -= pending;
        T from = (T)(*(pcurr-1)+1);
        T to = *pcurr;
        for (T i = from; ;++i) 
        {
            *dest_curr++ = i;
            if (i == to) break;
        }
        pcurr += 2; // jump to the next positive GAP
    }
    return (D) (dest_curr - dest);
}



/*! 
    @brief Bitcount for bit string
    
    Function calculates number of 1 bits in the given array of words.
    Make sure the addresses are aligned.

    @ingroup bitfunc 
*/
inline 
bm::id_t bit_block_count(const bm::word_t* block)
{
    const bm::word_t* block_end = block + bm::set_block_size;
    bm::id_t count = 0;

#ifdef BMVECTOPT
    count = VECT_BITCOUNT(block, block_end);
#else  
#ifdef BM64OPT
    // 64-bit optimized algorithm. No sparse vect opt.
    // instead it uses 4-way parallel pipelined version

    const bm::id64_t* b1 = (bm::id64_t*) block;
    const bm::id64_t* b2 = (bm::id64_t*) block_end;
    do 
    {
        count += bitcount64_4way(b1[0], b1[1], b1[2], b1[3]);
        b1 += 4;
    } while (b1 < b2);
#else
    // For 32 bit code the fastest method is
    // to use bitcount table for each byte in the block.
    // As optimization for sparse bitsets used bits accumulator
    // to collect ON bits using bitwise OR. 
    bm::word_t  acc = *block++;
    do
    {
        bm::word_t in = *block++;
        bm::word_t acc_prev = acc;
        acc |= in;
        if (acc_prev &= in)  // accumulator miss: counting bits
        {
            BM_INCWORD_BITCOUNT(count, acc);
            acc = acc_prev;
        }
    } while (block < block_end);

    BM_INCWORD_BITCOUNT(count, acc); // count-in remaining accumulator 

#endif
#endif	
    return count;
}

/*!
    @brief Bitcount for bit string

    Added for back-compatibility purposes, not block aligned,
    not SIMD accelerated

    @ingroup bitfunc
*/
inline
bm::id_t bit_block_calc_count(const bm::word_t* block,
                              const bm::word_t* block_end)
{
    bm::id_t count = 0;
    bm::word_t  acc = *block++;
    do
    {
        bm::word_t in = *block++;
        bm::word_t acc_prev = acc;
        acc |= in;
        if (acc_prev &= in)  // accumulator miss: counting bits
        {
            BM_INCWORD_BITCOUNT(count, acc);
            acc = acc_prev;
        }
    } while (block < block_end);

    BM_INCWORD_BITCOUNT(count, acc); // count-in remaining accumulator
    return count;
}



/*!
    Function calculates number of times when bit value changed 
    (1-0 or 0-1).
    
    For 001 result is 2
        010 - 3
        011 - 2
        111 - 1
    
    @ingroup bitfunc 
*/

inline 
bm::id_t bit_count_change(bm::word_t w)
{
    unsigned count = 1;
    w ^= (w >> 1);

    count += bm::word_bitcount(w);
    count -= (w >> ((sizeof(w) * 8) - 1));
    return count;
}


/*!
    Function calculates number of times when bit value changed 
    @internal
*/
inline
void bit_count_change32(const bm::word_t* block, 
                        const bm::word_t* block_end,
                        unsigned*         bit_count,
                        unsigned*         gap_count)
{
    BM_ASSERT(block < block_end);
    BM_ASSERT(bit_count);
    BM_ASSERT(gap_count);
    
    *gap_count = 1;
    *bit_count = 0;

    bm::word_t  w, w0, w_prev, w_l;     
    w = w0 = *block;
    
    BM_INCWORD_BITCOUNT(*bit_count, w);
    
    const int w_shift = int(sizeof(w) * 8 - 1);
    w ^= (w >> 1);
    BM_INCWORD_BITCOUNT(*gap_count, w);
    *gap_count -= (w_prev = (w0 >> w_shift)); // negative value correction

    for (++block ;block < block_end; ++block)
    {
        w = w0 = *block;
        ++(*gap_count);

        if (!w)
        {       
            *gap_count -= !w_prev;
            w_prev = 0;
        }
        else
        {
            BM_INCWORD_BITCOUNT(*bit_count, w);
        
            w ^= (w >> 1);
            BM_INCWORD_BITCOUNT(*gap_count, w);
            
            w_l = w0 & 1;
            *gap_count -= (w0 >> w_shift);  // negative value correction
            *gap_count -= !(w_prev ^ w_l);  // word border correction
            
            w_prev = (w0 >> w_shift);
        }
    } // for
}


/*!
    Function calculates number of times when bit value changed 
    (1-0 or 0-1) in the bit block.
    Also calulates number of bits ON.
    
    @param bit_count - OUT total number of bits ON
    @param block - bit-block start pointer
    @param block_end - bit-block end pointer
    
    @return number of 1-0, 0-1 transitions
        
    @ingroup bitfunc 
*/
inline 
bm::id_t bit_block_calc_count_change(const bm::word_t* block, 
                                     const bm::word_t* block_end,
                                     unsigned*         bit_count)
{
#if defined(BMSSE2OPT) || defined(BMSSE42OPT) || defined (BMAVX2OPT)

#ifdef BMAVX2OPT
    // TODO: debug true avx2 function (stress test failure)
    // temp use SSE4.2 variant
    return sse42_bit_block_calc_count_change(
        (const __m128i*)block, (const __m128i*)block_end, bit_count);
/*
    return avx2_bit_block_calc_count_change(
        (const __m256i*)block, (const __m256i*)block_end, bit_count);
*/
#else
#ifdef BMSSE42OPT
    return sse4_bit_block_calc_count_change(
        (const __m128i*)block, (const __m128i*)block_end, bit_count);
#else
# ifdef BMSSE2OPT
    return sse2_bit_block_calc_count_change(
        (const __m128i*)block, (const __m128i*)block_end, bit_count);
# endif
#endif
#endif

#else // non-SIMD code

    BM_ASSERT(block < block_end);
    BM_ASSERT(bit_count);
    
    
#ifdef BM64OPT
    bm::id64_t count = 1;
    *bit_count = 0;

    // 64-bit optimized algorithm.

    const bm::id64_t* b1 = (bm::id64_t*) block;
    const bm::id64_t* b2 = (bm::id64_t*) block_end;

    bm::id64_t w, w0, w_prev, w_l;
    w = w0 = *b1;
    
    *bit_count = word_bitcount64(w);
    
    const int w_shift = sizeof(w) * 8 - 1;
    w ^= (w >> 1);
    count += word_bitcount64(w);
    count -= (w_prev = (w0 >> w_shift)); // negative value correction

    
    for (++b1 ;b1 < b2; ++b1)
    {
        w = w0 = *b1;
        
        ++count;
        
        if (!w)
        {
            count -= !w_prev;
            w_prev = 0;
        }
        else
        {
            *bit_count += word_bitcount64(w);
            w ^= (w >> 1);
            count += word_bitcount64(w);
            
            w_l = w0 & 1;
            count -= (w0 >> w_shift);  // negative value correction
            count -= !(w_prev ^ w_l);  // word border correction
            
            w_prev = (w0 >> w_shift);
        }
    } // for
    return (bm::id_t) count;

#else
    unsigned gap_count;
    bit_count_change32(block, block_end, bit_count, &gap_count);
    return gap_count;
#endif

#endif
}


/*!
    Function calculates number of 1 bits in the given array of words in
    the range between left anf right bits (borders included)
    Make sure the addr is aligned.

    @ingroup bitfunc
*/
inline 
bm::id_t bit_block_calc_count_range(const bm::word_t* block,
                                    bm::word_t left,
                                    bm::word_t right)
{
    BM_ASSERT(left <= right);
    
    unsigned nword, nbit, bitcount, count;
    nbit = left & bm::set_word_mask;
    const bm::word_t* word = 
        block + (nword = unsigned(left >> bm::set_word_shift));
    if (left == right)  // special case (only 1 bit to check)
    {
        return (*word >> nbit) & 1u;
    }
    
    count = 0;
    bitcount = right - left + 1u;
    if (nbit) // starting position is not aligned
    {
        unsigned right_margin = nbit + right - left;
        if (right_margin < 32)
        {
            unsigned mask =
                block_set_table<true>::_right[nbit] &
                block_set_table<true>::_left[right_margin];
            return bm::word_bitcount(*word & mask);
        }
        count = bm::word_bitcount(*word & block_set_table<true>::_right[nbit]);
        bitcount -= 32 - nbit;
        ++word;
    }
    
    // now when we are word aligned, we can count bits the usual way
    //
    #if defined(BM64_SSE4) || defined(BM64_AVX2)
        for ( ;bitcount >= 64; bitcount-=64)
        {
            bm::id64_t w64 = *((bm::id64_t*)word); // x86 unaligned(!) read
            count += unsigned(_mm_popcnt_u64(w64));
            word += 2;
        }
        if (bitcount >= 32)
        {
            count += bm::word_bitcount(*word++);
            bitcount-=32;
        }
    #else
        for ( ;bitcount >= 32; bitcount-=32, ++word)
            count += bm::word_bitcount(*word);
    #endif
    BM_ASSERT(bitcount < 32);
    
    if (bitcount)  // we have a tail to count
    {
        count += bm::word_bitcount(
                        *word & block_set_table<true>::_left[bitcount-1]);
    }
    return count;
}

/*!
    Function calculates number of 1 bits in the given array of words in
    the range between 0 anf right bits (borders included)
    Make sure the addr is aligned.

    @ingroup bitfunc
*/
inline
bm::id_t bit_block_calc_count_to(const bm::word_t*  block,
                                 bm::word_t         right)
{
    BM_ASSERT(block);
    if (!right) // special case, first bit check
        return *block & 1u;
    bm::id_t count = 0;

    unsigned bitcount = right + 1;

    // AVX2 or 64-bit loop unroll
    #if defined(BMAVX2OPT)
        BM_AVX2_POPCNT_PROLOG
    
        __m256i cnt = _mm256_setzero_si256();
        bm::id64_t* cnt64;
    
        for ( ;bitcount >= 256; bitcount -= 256)
        {
            const __m256i* src = (__m256i*)block;
            __m256i xmm256 = _mm256_load_si256(src);
            BM_AVX2_BIT_COUNT(bc, xmm256);
            cnt = _mm256_add_epi64(cnt, bc);

            block += 8;
        }
        cnt64 = (bm::id64_t*)&cnt;
        count += (unsigned)(cnt64[0] + cnt64[1] + cnt64[2] + cnt64[3]);
    #endif
    
    for ( ;bitcount >= 64; bitcount -= 64)
    {
        bm::id64_t* p = (bm::id64_t*)block;
        count += bm::word_bitcount64(*p);
        block += 2;
    }
    if (bitcount >= 32)
    {
        count += bm::word_bitcount(*block++);
        bitcount-=32;
    }

    if (bitcount)  // we have a tail to count
    {
        count +=
        bm::word_bitcount(*block & block_set_table<true>::_left[bitcount-1]);
    }
    return count;
}



/*!
    Cyclic rotation of bit-block left by 1 bit
    @ingroup bitfunc
*/
inline
void bit_block_rotate_left_1(bm::word_t* block)
{
    bm::word_t co_flag = (block[0] >> 31) & 1; // carry over bit
    for (unsigned i = 0; i < bm::set_block_size-1; ++i)
    {
        block[i] = (block[i] << 1) | (block[i + 1] >> 31);
    } 
    block[set_block_size - 1] = (block[set_block_size - 1] << 1) | co_flag;
}

/*!
    Unrolled cyclic rotation of bit-block left by 1 bit
    @ingroup bitfunc
*/
inline
void bit_block_rotate_left_1_unr(bm::word_t* block)
{
    bm::word_t co_flag = (block[0] >> 31) & 1; // carry over bit
    const unsigned unroll_factor = 4;
    bm::word_t w0, w1, w2, w3;

    unsigned i;
    for (i = 0; i < bm::set_block_size - unroll_factor; i += unroll_factor)
    {
        w0 = block[i + 1] >> 31;
        w1 = block[i + 2] >> 31;
        w2 = block[i + 3] >> 31;
        w3 = block[i + 4] >> 31;

        block[0 + i] = (block[0 + i] << 1) | w0;
        block[1 + i] = (block[1 + i] << 1) | w1;
        block[2 + i] = (block[2 + i] << 1) | w2;
        block[3 + i] = (block[3 + i] << 1) | w3;
    }
    block[i] = (block[i] << 1) | (block[i + 1] >> 31);
    block[i + 1] = (block[i + 1] << 1) | (block[i + 2] >> 31);
    block[i + 2] = (block[i + 2] << 1) | (block[i + 3] >> 31);
    block[set_block_size - 1] = (block[set_block_size - 1] << 1) | co_flag;
}



/*!
    Function calculates if there is any number of 1 bits 
    in the given array of words in the range between left anf right bits 
    (borders included). Make sure the addresses are aligned.

    @ingroup bitfunc
*/
inline 
bm::id_t bit_block_any_range(const bm::word_t* block,
                             bm::word_t left,
                             bm::word_t right)
{
    BM_ASSERT(left <= right);
    
    unsigned nbit  = left; // unsigned(left & bm::set_block_mask);
    unsigned nword = unsigned(nbit >> bm::set_word_shift);
    nbit &= bm::set_word_mask;

    const bm::word_t* word = block + nword;

    if (left == right)  // special case (only 1 bit to check)
    {
        return (*word >> nbit) & 1;
    }
    unsigned acc;
    unsigned bitcount = right - left + 1;

    if (nbit) // starting position is not aligned
    {
        unsigned right_margin = nbit + (right - left);
        if (right_margin < 32)
        {
            unsigned mask =
                block_set_table<true>::_right[nbit] &
                block_set_table<true>::_left[right_margin];
            acc = *word & mask;
            return acc;
        }
        else
        {
            acc = *word & block_set_table<true>::_right[nbit];
            if (acc) 
                return acc;
            bitcount -= 32 - nbit;
        }
        ++word;
    }

    // now when we are word aligned, we can check bits the usual way
    for ( ;bitcount >= 32; bitcount -= 32)
    {
        acc = *word++;
        if (acc) 
            return acc;
    }

    if (bitcount)  // we have a tail to count
    {
        acc = (*word) & block_set_table<true>::_left[bitcount-1];
        if (acc) 
            return acc;
    }

    return 0;
}

// ----------------------------------------------------------------------

/*! Function inverts block of bits 
    @ingroup bitfunc 
*/
template<typename T> void bit_invert(T* start, T* end)
{
#ifdef BMVECTOPT
    VECT_INVERT_ARR(start, end);
#else
    do
    {
        start[0] = ~start[0];
        start[1] = ~start[1];
        start[2] = ~start[2];
        start[3] = ~start[3];
        start+=4;
    } while (start < end);
#endif
}

// ----------------------------------------------------------------------

/*! @brief Returns "true" if all bits in the block are 1
    @ingroup bitfunc 
*/
inline
bool is_bits_one(const bm::wordop_t* start,
                 const bm::wordop_t* end)
{
#if defined(BMSSE42OPT) || defined(BMAVX2OPT)
    return VECT_IS_ONE_BLOCK(start, end);
#else
   do
   {
        bm::wordop_t tmp = 
            start[0] & start[1] & start[2] & start[3];
        if (tmp != bm::all_bits_mask) 
            return false;
        start += 4;
   } while (start < end);
   return true;
#endif
}

// ----------------------------------------------------------------------

// GAP blocks manipulation functions:

/*! \brief GAP and functor */
BMFORCEINLINE unsigned and_op(unsigned v1, unsigned v2)
{
    return v1 & v2;
}


/*! \brief GAP xor functor */
BMFORCEINLINE unsigned xor_op(unsigned v1, unsigned v2)
{
    return v1 ^ v2;
}


/*! \brief GAP or functor */
BMFORCEINLINE unsigned or_op(unsigned v1, unsigned v2)
{
    return v1 | v2;
}

/*! \brief GAP or functor */
BMFORCEINLINE unsigned sub_op(unsigned v1, unsigned v2)
{
    return v1 & ~v2;
}


/*!
   \brief GAP AND operation.
   
   Function performs AND logical operation on gap vectors.
   If possible function put the result into vect1 and returns this
   pointer.  Otherwise result is put into tmp_buf, which should be 
   twice of the vector size.

   \param vect1   - operand 1
   \param vect2   - operand 2
   \param tmp_buf - pointer on temporary buffer
   \param dsize   - out size of the destination
   \return Result pointer (tmp_buf OR vect1)

   @ingroup gapfunc
*/
BMFORCEINLINE 
gap_word_t* gap_operation_and(const gap_word_t* BMRESTRICT vect1,
                              const gap_word_t* BMRESTRICT vect2,
                              gap_word_t*       BMRESTRICT tmp_buf,
                              unsigned&         dsize)
{
    gap_buff_op(tmp_buf, vect1, 0, vect2, 0, and_op, dsize);
    return tmp_buf;
}

/*!
   \brief GAP AND operation test.
   
   Function performs AND logical operation on gap vectors.
   If possible function put the result into vect1 and returns this
   pointer.  Otherwise result is put into tmp_buf, which should be 
   twice of the vector size.

   \param vect1   - operand 1
   \param vect2   - operand 2
   \return non zero value if operation returns any 1 bit

   @ingroup gapfunc
*/
BMFORCEINLINE 
unsigned gap_operation_any_and(const gap_word_t* BMRESTRICT vect1,
                                      const gap_word_t* BMRESTRICT vect2)
{
    return gap_buff_any_op(vect1, 0, vect2, 0, and_op);
}


/*!
   \brief GAP bitcount AND operation test.
   
   \param vect1   - operand 1
   \param vect2   - operand 2
   \return bitcount of vect1 AND vect2

   @ingroup gapfunc
*/
BMFORCEINLINE 
unsigned gap_count_and(const gap_word_t* BMRESTRICT vect1,
                       const gap_word_t* BMRESTRICT vect2)
{
    return gap_buff_count_op(vect1, vect2, and_op);
}



/*!
   \brief GAP XOR operation.
   
   Function performs XOR logical operation on gap vectors.
   If possible function put the result into vect1 and returns this
   pointer.  Otherwise result is put into tmp_buf, which should be 
   twice of the vector size.

   \param vect1   - operand 1
   \param vect2   - operand 2
   \param tmp_buf - pointer on temporary buffer
   \param dsize   - out destination size
   \return Result pointer (tmp_buf)

   @ingroup gapfunc
*/
BMFORCEINLINE 
gap_word_t* gap_operation_xor(const gap_word_t*  BMRESTRICT vect1,
                              const gap_word_t*  BMRESTRICT vect2,
                              gap_word_t*        BMRESTRICT tmp_buf,
                              unsigned&                     dsize)
{
    gap_buff_op(tmp_buf, vect1, 0, vect2, 0, bm::xor_op, dsize);
    return tmp_buf;
}


/*!
   \brief GAP XOR operation test.
   
   Function performs AND logical operation on gap vectors.
   If possible function put the result into vect1 and returns this
   pointer.  Otherwise result is put into tmp_buf, which should be 
   twice of the vector size.

   \param vect1   - operand 1
   \param vect2   - operand 2
   \return non zero value if operation returns any 1 bit

   @ingroup gapfunc
*/
BMFORCEINLINE 
unsigned gap_operation_any_xor(const gap_word_t* BMRESTRICT vect1,
                               const gap_word_t* BMRESTRICT vect2)
{
    return gap_buff_any_op(vect1, 0, vect2, 0, bm::xor_op);
}

/*!
   \brief GAP bitcount XOR operation test.
   
   \param vect1   - operand 1
   \param vect2   - operand 2
   \return bitcount of vect1 XOR vect2

   @ingroup gapfunc
*/
BMFORCEINLINE 
unsigned gap_count_xor(const gap_word_t* BMRESTRICT vect1,
                       const gap_word_t* BMRESTRICT vect2)
{
    return gap_buff_count_op(vect1, vect2, bm::xor_op);
}


/*!
   \brief GAP OR operation.
   
   Function performs OR logical oparation on gap vectors.
   If possible function put the result into vect1 and returns this
   pointer.  Otherwise result is put into tmp_buf, which should be 
   twice of the vector size.

   \param vect1   - operand 1
   \param vect2   - operand 2
   \param tmp_buf - pointer on temporary buffer
   \param dsize   - out destination size

   \return Result pointer (tmp_buf)

   @ingroup gapfunc
*/
inline 
gap_word_t* gap_operation_or(const gap_word_t*  BMRESTRICT vect1,
                             const gap_word_t*  BMRESTRICT vect2,
                             gap_word_t*        BMRESTRICT tmp_buf,
                             unsigned&                     dsize)
{
    gap_buff_op(tmp_buf, vect1, 1, vect2, 1, bm::and_op, dsize);
    gap_invert(tmp_buf);
    return tmp_buf;
}

/*!
   \brief GAP bitcount OR operation test.
   
   \param vect1   - operand 1
   \param vect2   - operand 2
   \return bitcount of vect1 OR vect2

   @ingroup gapfunc
*/
BMFORCEINLINE 
unsigned gap_count_or(const gap_word_t* BMRESTRICT vect1,
                      const gap_word_t* BMRESTRICT vect2)
{
    return gap_buff_count_op(vect1, vect2, bm::or_op);
}



/*!
   \brief GAP SUB (AND NOT) operation.
   
   Function performs SUB logical operation on gap vectors.
   If possible function put the result into vect1 and returns this
   pointer.  Otherwise result is put into tmp_buf, which should be 
   twice of the vector size.

   \param vect1   - operand 1
   \param vect2   - operand 2
   \param tmp_buf - pointer on temporary buffer
   \param dsize   - out destination size

   \return Result pointer (tmp_buf)

   @ingroup gapfunc
*/
inline gap_word_t* gap_operation_sub(const gap_word_t*  BMRESTRICT vect1,
                                     const gap_word_t*  BMRESTRICT vect2,
                                     gap_word_t*        BMRESTRICT tmp_buf,
                                     unsigned&                     dsize)
{
    gap_buff_op(tmp_buf, vect1, 0, vect2, 1, and_op, dsize);    
    return tmp_buf;
}


/*!
   \brief GAP SUB operation test.
   
   Function performs AND logical operation on gap vectors.
   If possible function put the result into vect1 and returns this
   pointer.  Otherwise result is put into tmp_buf, which should be 
   twice of the vector size.

   \param vect1   - operand 1
   \param vect2   - operand 2
   \return non zero value if operation returns any 1 bit

   @ingroup gapfunc
*/
BMFORCEINLINE 
unsigned gap_operation_any_sub(const gap_word_t* BMRESTRICT vect1,
                               const gap_word_t* BMRESTRICT vect2)
{
    return gap_buff_any_op(vect1, 0, vect2, 1, bm::and_op);    
}


/*!
\brief GAP bitcount SUB (AND NOT) operation test.

\param vect1   - operand 1
\param vect2   - operand 2
\return bitcount of vect1 SUB (AND NOT) vect2

@ingroup gapfunc
*/
BMFORCEINLINE 
unsigned gap_count_sub(const gap_word_t* BMRESTRICT vect1,
                       const gap_word_t* BMRESTRICT vect2)
{
    return gap_buff_count_op(vect1, vect2, bm::sub_op);
}


// ----------------------------------------------------------------------

// BIT blocks manipulation functions:


/*!
   \brief Bitblock copy operation. 

   \param dst - destination block.
   \param src - source block.

   @ingroup bitfunc
*/
inline 
void bit_block_copy(bm::word_t* BMRESTRICT dst, const bm::word_t* BMRESTRICT src)
{
#ifdef BMVECTOPT
    VECT_COPY_BLOCK(dst, src, src + bm::set_block_size);
#else
    ::memcpy(dst, src, bm::set_block_size * sizeof(bm::word_t));
#endif
}


/*!
   \brief Plain bitblock AND operation. 
   Function does not analyse availability of source and destination blocks.

   \param dst - destination block.
   \param src - source block.
 
   \return 0 if AND operation did not produce anything (no 1s in the output)

   @ingroup bitfunc
*/
inline 
bm::id64_t bit_block_and(bm::word_t* BMRESTRICT dst, const bm::word_t* BMRESTRICT src)
{
    BM_ASSERT(dst);
    BM_ASSERT(src);
    BM_ASSERT(dst != src);

#ifdef BMVECTOPT
    bm::id64_t acc = VECT_AND_ARR(dst, src, src + bm::set_block_size);
    return acc;
#else
    unsigned arr_sz = bm::set_block_size / 2;

    const bm::bit_block_t::bunion_t* BMRESTRICT src_u = (const bm::bit_block_t::bunion_t*)src;
    bm::bit_block_t::bunion_t* BMRESTRICT dst_u = (bm::bit_block_t::bunion_t*)dst;

    bm::id64_t acc = 0;
    for (unsigned i = 0; i < arr_sz; i+=4)
    {
        acc |= dst_u->w64[i] &= src_u->w64[i];
        acc |= dst_u->w64[i+1] &= src_u->w64[i+1];
        acc |= dst_u->w64[i+2] &= src_u->w64[i+2];
        acc |= dst_u->w64[i+3] &= src_u->w64[i+3];
    }
    return acc;
#endif
}


/*!
   \brief Function ANDs two bitblocks and computes the bitcount. 
   Function does not analyse availability of source blocks.

   \param src1     - first bit block
   \param src2     - second bit block

   @ingroup bitfunc
*/
inline 
unsigned bit_block_and_count(const bm::word_t* BMRESTRICT src1,
                             const bm::word_t* BMRESTRICT src2)
{
    unsigned count;
    const bm::word_t* src1_end = src1 + bm::set_block_size;
#ifdef BMVECTOPT
    count = VECT_BITCOUNT_AND(src1, src1_end, src2);
#else  
    count = 0;  
# ifdef BM64OPT
    const bm::id64_t* b1 = (bm::id64_t*) src1;
    const bm::id64_t* b1_end = (bm::id64_t*) src1_end;
    const bm::id64_t* b2 = (bm::id64_t*) src2;
    do
    {
        count += bitcount64_4way(b1[0] & b2[0], 
                                 b1[1] & b2[1], 
                                 b1[2] & b2[2], 
                                 b1[3] & b2[3]);
        b1 += 4;
        b2 += 4;
    } while (b1 < b1_end);
# else
    do
    {
        BM_INCWORD_BITCOUNT(count, src1[0] & src2[0]);
        BM_INCWORD_BITCOUNT(count, src1[1] & src2[1]);
        BM_INCWORD_BITCOUNT(count, src1[2] & src2[2]);
        BM_INCWORD_BITCOUNT(count, src1[3] & src2[3]);

        src1+=4;
        src2+=4;
    } while (src1 < src1_end);
# endif
#endif    
    return count;
}


/*!
   \brief Function ANDs two bitblocks and tests for any bit. 
   Function does not analyse availability of source blocks.

   \param src1     - first bit block
   \param src2     - second bit block

   @ingroup bitfunc
*/
inline 
unsigned bit_block_and_any(const bm::word_t* src1, 
                           const bm::word_t* src2)
{
    unsigned count = 0;
    const bm::word_t* src1_end = src1 + bm::set_block_size;
    do
    {
        count = (src1[0] & src2[0]) |
                (src1[1] & src2[1]) |
                (src1[2] & src2[2]) |
                (src1[3] & src2[3]);

        src1+=4; src2+=4;
    } while ((src1 < src1_end) && !count);
    return count;
}




/*!
   \brief Function XORs two bitblocks and computes the bitcount. 
   Function does not analyse availability of source blocks.

   \param src1     - first bit block
   \param src2     - second bit block

   @ingroup bitfunc
*/
inline 
unsigned bit_block_xor_count(const bm::word_t* BMRESTRICT src1,
                             const bm::word_t* BMRESTRICT src2)
{
    unsigned count;
    const bm::word_t* BMRESTRICT src1_end = src1 + bm::set_block_size;
#ifdef BMVECTOPT
    count = VECT_BITCOUNT_XOR(src1, src1_end, src2);
#else  
    count = 0;  
# ifdef BM64OPT
    const bm::id64_t* b1 = (bm::id64_t*) src1;
    const bm::id64_t* b1_end = (bm::id64_t*) src1_end;
    const bm::id64_t* b2 = (bm::id64_t*) src2;
    do
    {
        count += bitcount64_4way(b1[0] ^ b2[0], 
                                 b1[1] ^ b2[1], 
                                 b1[2] ^ b2[2], 
                                 b1[3] ^ b2[3]);
        b1 += 4;
        b2 += 4;
    } while (b1 < b1_end);
# else
    do
    {
        BM_INCWORD_BITCOUNT(count, src1[0] ^ src2[0]);
        BM_INCWORD_BITCOUNT(count, src1[1] ^ src2[1]);
        BM_INCWORD_BITCOUNT(count, src1[2] ^ src2[2]);
        BM_INCWORD_BITCOUNT(count, src1[3] ^ src2[3]);

        src1+=4;
        src2+=4;
    } while (src1 < src1_end);
# endif
#endif
    return count;
}


/*!
   \brief Function XORs two bitblocks and and tests for any bit.
   Function does not analyse availability of source blocks.

   \param src1     - first bit block.
   \param src2     - second bit block.

   @ingroup bitfunc
*/
inline 
unsigned bit_block_xor_any(const bm::word_t* BMRESTRICT src1,
                           const bm::word_t* BMRESTRICT src2)
{
    unsigned count = 0;
    const bm::word_t* BMRESTRICT src1_end = src1 + bm::set_block_size;
    do
    {
        count = (src1[0] ^ src2[0]) |
                (src1[1] ^ src2[1]) |
                (src1[2] ^ src2[2]) |
                (src1[3] ^ src2[3]);

        src1+=4; src2+=4;
    } while (!count && (src1 < src1_end));
    return count;
}




/*!
   \brief Function SUBs two bitblocks and computes the bitcount. 
   Function does not analyse availability of source blocks.

   \param src1     - first bit block.
   \param src2     - second bit block.

   @ingroup bitfunc
*/
inline 
unsigned bit_block_sub_count(const bm::word_t* BMRESTRICT src1,
                             const bm::word_t* BMRESTRICT src2)
{
    unsigned count;
    const bm::word_t* BMRESTRICT src1_end = src1 + bm::set_block_size;
#ifdef BMVECTOPT
    count = VECT_BITCOUNT_SUB(src1, src1_end, src2);
#else  
    count = 0;  
# ifdef BM64OPT
    const bm::id64_t* b1 = (bm::id64_t*) src1;
    const bm::id64_t* b1_end = (bm::id64_t*) src1_end;
    const bm::id64_t* b2 = (bm::id64_t*) src2;
    do
    {
        count += bitcount64_4way(b1[0] & ~b2[0], 
                                 b1[1] & ~b2[1], 
                                 b1[2] & ~b2[2], 
                                 b1[3] & ~b2[3]);
        b1 += 4;
        b2 += 4;
    } while (b1 < b1_end);
# else
    do
    {
        BM_INCWORD_BITCOUNT(count, src1[0] & ~src2[0]);
        BM_INCWORD_BITCOUNT(count, src1[1] & ~src2[1]);
        BM_INCWORD_BITCOUNT(count, src1[2] & ~src2[2]);
        BM_INCWORD_BITCOUNT(count, src1[3] & ~src2[3]);

        src1+=4;
        src2+=4;
    } while (src1 < src1_end);
# endif
#endif
    return count;
}

/*!
   \brief Function SUBs two bitblocks and and tests for any bit.
   Function does not analyse availability of source blocks.

   \param src1     - first bit block.
   \param src2     - second bit block.

   @ingroup bitfunc
*/
inline 
unsigned bit_block_sub_any(const bm::word_t* BMRESTRICT src1,
                           const bm::word_t* BMRESTRICT src2)
{
    unsigned count = 0;
    const bm::word_t* BMRESTRICT src1_end = src1 + bm::set_block_size;

    do
    {
        count = (src1[0] & ~src2[0]) |
                (src1[1] & ~src2[1]) |
                (src1[2] & ~src2[2]) |
                (src1[3] & ~src2[3]);

        src1+=4; src2+=4;
    } while ((src1 < src1_end) && (count == 0));
    return count;
}



/*!
   \brief Function ORs two bitblocks and computes the bitcount. 
   Function does not analyse availability of source blocks.

   \param src1     - first bit block
   \param src2     - second bit block.

   @ingroup bitfunc
*/
inline 
unsigned bit_block_or_count(const bm::word_t* src1, 
                            const bm::word_t* src2)
{
    unsigned count;
    const bm::word_t* src1_end = src1 + bm::set_block_size;
#ifdef BMVECTOPT
    count = VECT_BITCOUNT_OR(src1, src1_end, src2);
#else  
    count = 0;  
# ifdef BM64OPT
    const bm::id64_t* b1 = (bm::id64_t*) src1;
    const bm::id64_t* b1_end = (bm::id64_t*) src1_end;
    const bm::id64_t* b2 = (bm::id64_t*) src2;
    do
    {
        count += bitcount64_4way(b1[0] | b2[0], 
                                 b1[1] | b2[1], 
                                 b1[2] | b2[2], 
                                 b1[3] | b2[3]);
        b1 += 4;
        b2 += 4;
    } while (b1 < b1_end);
# else
    do
    {
        BM_INCWORD_BITCOUNT(count, src1[0] | src2[0]);
        BM_INCWORD_BITCOUNT(count, src1[1] | src2[1]);
        BM_INCWORD_BITCOUNT(count, src1[2] | src2[2]);
        BM_INCWORD_BITCOUNT(count, src1[3] | src2[3]);

        src1+=4;
        src2+=4;
    } while (src1 < src1_end);
# endif
#endif
    return count;
}

/*!
   \brief Function ORs two bitblocks and and tests for any bit.
   Function does not analyse availability of source blocks.

   \param src1     - first bit block.
   \param src2     - second bit block.

   @ingroup bitfunc
*/
inline 
unsigned bit_block_or_any(const bm::word_t* BMRESTRICT src1,
                          const bm::word_t* BMRESTRICT src2)
{
    unsigned count = 0;
    const bm::word_t* BMRESTRICT src1_end = src1 + bm::set_block_size;
    do
    {
        count = (src1[0] | src2[0]) |
                (src1[1] | src2[1]) |
                (src1[2] | src2[2]) |
                (src1[3] | src2[3]);

        src1+=4; src2+=4;
    } while (!count && (src1 < src1_end));
    return count;
}




/*!
   \brief bitblock AND operation. 

   \param dst - destination block.
   \param src - source block.

   \returns pointer on destination block. 
    If returned value  equal to src means that block mutation requested. 
    NULL is valid return value.

   @ingroup bitfunc
*/
inline bm::word_t* bit_operation_and(bm::word_t* BMRESTRICT dst, 
                                     const bm::word_t* BMRESTRICT src)
{
    BM_ASSERT(dst || src);

    bm::word_t* ret = dst;

    if (IS_VALID_ADDR(dst))  // The destination block already exists
    {
        if (!IS_VALID_ADDR(src))
        {
            if (IS_EMPTY_BLOCK(src))
            {
                //If the source block is zero 
                //just clean the destination block
                return 0;
            }
        }
        else
        {
            // Regular operation AND on the whole block.
            //
            auto any = bm::bit_block_and(dst, src);
            if (!any)
                ret = 0;
        }
    }
    else // The destination block does not exist yet
    {
        if(!IS_VALID_ADDR(src))
        {
            if(IS_EMPTY_BLOCK(src)) 
            {
                // The source block is empty.
                // One argument empty - all result is empty.
                return 0;
            }
            // Here we have nothing to do.
            // Src block is all ON, dst block remains as it is
        }
        else // destination block does not exists, src - valid block
        {
            if (IS_FULL_BLOCK(dst))
                return const_cast<bm::word_t*>(src);
            // Nothng to do.
            // Dst block is all ZERO no combination required.
        }
    }

    return ret;
}


/*!
   \brief Performs bitblock AND operation and calculates bitcount of the result. 

   \param src1     - first bit block.
   \param src2     - second bit block.

   \returns bitcount value 

   @ingroup bitfunc
*/
inline 
bm::id_t bit_operation_and_count(const bm::word_t* BMRESTRICT src1,
                                 const bm::word_t* BMRESTRICT src2)
{
    if (IS_EMPTY_BLOCK(src1) || IS_EMPTY_BLOCK(src2))
        return 0;
    return bit_block_and_count(src1, src2);
}

/*!
   \brief Performs bitblock AND operation test. 

   \param src1     - first bit block.
   \param src2     - second bit block.

   \returns non zero if there is any value 

   @ingroup bitfunc
*/
inline 
bm::id_t bit_operation_and_any(const bm::word_t* BMRESTRICT src1,
                               const bm::word_t* BMRESTRICT src2)
{
    if (IS_EMPTY_BLOCK(src1) || IS_EMPTY_BLOCK(src2))
        return 0;
    return bit_block_and_any(src1, src2);
}



/*!
   \brief Performs bitblock SUB operation and calculates bitcount of the result. 

   \param src1      - first bit block.
   \param src2      - second bit block

   \returns bitcount value 

   @ingroup bitfunc
*/
inline 
bm::id_t bit_operation_sub_count(const bm::word_t* BMRESTRICT src1, 
                                 const bm::word_t* BMRESTRICT src2)
{
    if (IS_EMPTY_BLOCK(src1))
        return 0;
    
    if (IS_EMPTY_BLOCK(src2)) // nothing to diff
    {
        return bit_block_count(src1);
    }
    return bit_block_sub_count(src1, src2);
}


/*!
   \brief Performs inverted bitblock SUB operation and calculates 
          bitcount of the result. 

   \param src1      - first bit block.
   \param src2      - second bit block

   \returns bitcount value 

   @ingroup bitfunc
*/
inline 
bm::id_t bit_operation_sub_count_inv(const bm::word_t* BMRESTRICT src1, 
                                     const bm::word_t* BMRESTRICT src2)
{
    return bit_operation_sub_count(src2, src1);
}


/*!
   \brief Performs bitblock test of SUB operation. 

   \param src1      - first bit block.
   \param src2      - second bit block

   \returns non zero value if there are any bits

   @ingroup bitfunc
*/
inline 
bm::id_t bit_operation_sub_any(const bm::word_t* BMRESTRICT src1, 
                               const bm::word_t* BMRESTRICT src2)
{
    if (IS_EMPTY_BLOCK(src1))
        return 0;
    
    if (IS_EMPTY_BLOCK(src2)) // nothing to diff
        return !bit_is_all_zero(src1);
    return bit_block_sub_any(src1, src2);
}



/*!
   \brief Performs bitblock OR operation and calculates bitcount of the result. 

   \param src1     - first bit block.
   \param src1_end - first bit block end
   \param src2     - second bit block.

   \returns bitcount value 

   @ingroup bitfunc
*/
inline 
bm::id_t bit_operation_or_count(const bm::word_t* BMRESTRICT src1,
                                const bm::word_t* BMRESTRICT src2)
{
    if (IS_EMPTY_BLOCK(src1))
    {
        if (!IS_EMPTY_BLOCK(src2))
            return bit_block_count(src2);
        else
            return 0; // both blocks are empty        
    }
    else
    {
        if (IS_EMPTY_BLOCK(src2))
            return bit_block_count(src1);
    }

    return bit_block_or_count(src1, src2);
}

/*!
   \brief Performs bitblock OR operation test. 

   \param src1     - first bit block.
   \param src2     - second bit block.

   \returns non zero value if there are any bits

   @ingroup bitfunc
*/
inline 
bm::id_t bit_operation_or_any(const bm::word_t* BMRESTRICT src1,
                              const bm::word_t* BMRESTRICT src2)
{
    if (IS_EMPTY_BLOCK(src1))
    {
        if (!IS_EMPTY_BLOCK(src2))
            return !bit_is_all_zero(src2);
        else
            return 0; // both blocks are empty        
    }
    else
    {
        if (IS_EMPTY_BLOCK(src2))
            return !bit_is_all_zero(src1);
    }

    return bit_block_or_any(src1, src2);
}



/*!
   \brief Plain bitblock OR operation. 
   Function does not analyse availability of source and destination blocks.

   \param dst - destination block.
   \param src - source block.

   @ingroup bitfunc
*/
inline 
bool bit_block_or(bm::word_t* BMRESTRICT dst, 
                  const bm::word_t* BMRESTRICT src)
{
#ifdef BMVECTOPT
    return VECT_OR_ARR(dst, src, src + bm::set_block_size);
#else
    const bm::wordop_t* BMRESTRICT wrd_ptr = (wordop_t*)src;
    const bm::wordop_t* BMRESTRICT wrd_end = (wordop_t*)(src + bm::set_block_size);
    bm::wordop_t* BMRESTRICT dst_ptr = (wordop_t*)dst;

    bm::wordop_t acc = 0;
    const bm::wordop_t not_acc = acc = ~acc;

    do
    {
        acc &= (dst_ptr[0] |= wrd_ptr[0]);
        acc &= (dst_ptr[1] |= wrd_ptr[1]);
        acc &= (dst_ptr[2] |= wrd_ptr[2]);
        acc &= (dst_ptr[3] |= wrd_ptr[3]);

        dst_ptr+=4;wrd_ptr+=4;

    } while (wrd_ptr < wrd_end);
    return acc == not_acc;
#endif
}

/*!
   \brief 3 way (target, source1, source2) bitblock OR operation.
   Function does not analyse availability of source and destination blocks.

   \param dst - destination block.
   \param src - source block.
 
   @return 1 if produced block of ALL ones

   @ingroup bitfunc
*/
inline
bool bit_block_or_3way(bm::word_t* BMRESTRICT dst,
                        const bm::word_t* BMRESTRICT src1,
                        const bm::word_t* BMRESTRICT src2)
{
#ifdef BMVECTOPT
    return VECT_OR_ARR_3WAY(dst, src1, src2, src1 + bm::set_block_size);
#else
    const bm::wordop_t* BMRESTRICT wrd_ptr1 = (wordop_t*)src1;
    const bm::wordop_t* BMRESTRICT wrd_end1 = (wordop_t*)(src1 + set_block_size);
    const bm::wordop_t* BMRESTRICT wrd_ptr2 = (wordop_t*)src2;
    bm::wordop_t* BMRESTRICT dst_ptr = (wordop_t*)dst;

    bm::wordop_t acc = 0;
    const bm::wordop_t not_acc = acc = ~acc;
    do
    {
        acc &= (dst_ptr[0] |= wrd_ptr1[0] | wrd_ptr2[0]);
        acc &= (dst_ptr[1] |= wrd_ptr1[1] | wrd_ptr2[1]);
        acc &= (dst_ptr[2] |= wrd_ptr1[2] | wrd_ptr2[2]);
        acc &= (dst_ptr[3] |= wrd_ptr1[3] | wrd_ptr2[3]);
        
        dst_ptr+=4; wrd_ptr1+=4;wrd_ptr2+=4;
        
    } while (wrd_ptr1 < wrd_end1);
    return acc == not_acc;
#endif
}


/*!
   \brief 5 way (target, source1, source2) bitblock OR operation.
   Function does not analyse availability of source and destination blocks.

   \param dst - destination block.
   \param src - source block.
 
   @return 1 if produced block of ALL ones

   @ingroup bitfunc
*/
inline
bool bit_block_or_5way(bm::word_t* BMRESTRICT dst,
                        const bm::word_t* BMRESTRICT src1,
                        const bm::word_t* BMRESTRICT src2,
                        const bm::word_t* BMRESTRICT src3,
                        const bm::word_t* BMRESTRICT src4)
{
#ifdef BMVECTOPT
    return VECT_OR_ARR_5WAY(dst, src1, src2, src3, src4, src1 + bm::set_block_size);
#else
    const bm::wordop_t* BMRESTRICT wrd_ptr1 = (wordop_t*)src1;
    const bm::wordop_t* BMRESTRICT wrd_end1 = (wordop_t*)(src1 + set_block_size);
    const bm::wordop_t* BMRESTRICT wrd_ptr2 = (wordop_t*)src2;
    const bm::wordop_t* BMRESTRICT wrd_ptr3 = (wordop_t*)src3;
    const bm::wordop_t* BMRESTRICT wrd_ptr4 = (wordop_t*)src4;
    bm::wordop_t* BMRESTRICT dst_ptr = (wordop_t*)dst;

    bm::wordop_t acc = 0;
    const bm::wordop_t not_acc = acc = ~acc;
    do
    {
        acc &= (dst_ptr[0] |= wrd_ptr1[0] | wrd_ptr2[0] | wrd_ptr3[0] | wrd_ptr4[0]);
        acc &= (dst_ptr[1] |= wrd_ptr1[1] | wrd_ptr2[1] | wrd_ptr3[1] | wrd_ptr4[1]);
        acc &= (dst_ptr[2] |= wrd_ptr1[2] | wrd_ptr2[2] | wrd_ptr3[2] | wrd_ptr4[2]);
        acc &= (dst_ptr[3] |= wrd_ptr1[3] | wrd_ptr2[3] | wrd_ptr3[3] | wrd_ptr4[3]);
        
        dst_ptr+=4;
        wrd_ptr1+=4;wrd_ptr2+=4;wrd_ptr3+=4;wrd_ptr4+=4;
        
    } while (wrd_ptr1 < wrd_end1);
    return acc == not_acc;
#endif
}




/*!
   \brief Block OR operation. Makes analysis if block is 0 or FULL. 

   \param dst - destination block.
   \param src - source block.

   \returns pointer on destination block. 
    If returned value  equal to src means that block mutation requested. 
    NULL is valid return value.

   @ingroup bitfunc
*/
inline 
bm::word_t* bit_operation_or(bm::word_t* BMRESTRICT dst, 
                             const bm::word_t* BMRESTRICT src)
{
    BM_ASSERT(dst || src);

    bm::word_t* ret = dst;

    if (IS_VALID_ADDR(dst)) // The destination block already exists
    {
        if (!IS_VALID_ADDR(src))
        {
            if (IS_FULL_BLOCK(src))
            {
                // if the source block is all set 
                // just set the destination block
                ::memset(dst, 0xFF, bm::set_block_size * sizeof(bm::word_t));
            }
        }
        else
        {
            // Regular operation OR on the whole block
            bm::bit_block_or(dst, src);
        }
    }
    else // The destination block does not exist yet
    {
        if (!IS_VALID_ADDR(src))
        {
            if (IS_FULL_BLOCK(src)) 
            {
                // The source block is all set, because dst does not exist
                // we can simply replace it.
                return const_cast<bm::word_t*>(FULL_BLOCK_FAKE_ADDR);
            }
        }
        else
        {
            if (dst == 0)
            {
                // The only case when we have to allocate the new block:
                // Src is all zero and Dst does not exist
                return const_cast<bm::word_t*>(src);
            }
        }
    }
    return ret;
}

/*!
   \brief Plain bitblock SUB (AND NOT) operation. 
   Function does not analyse availability of source and destination blocks.

   \param dst - destination block.
   \param src - source block.

   \return 0 if SUB operation did not produce anything (no 1s in the output)

   @ingroup bitfunc
*/
inline
bm::id64_t bit_block_sub(bm::word_t* BMRESTRICT dst, const bm::word_t* BMRESTRICT src)
{
#ifdef BMVECTOPT
    bm::id64_t acc = VECT_SUB_ARR(dst, src, src + bm::set_block_size);
    return acc;
#else
    unsigned arr_sz = bm::set_block_size / 2;

    const bm::bit_block_t::bunion_t* BMRESTRICT src_u = (const bm::bit_block_t::bunion_t*)src;
    bm::bit_block_t::bunion_t* BMRESTRICT dst_u = (bm::bit_block_t::bunion_t*)dst;

    bm::id64_t acc = 0;
    for (unsigned i = 0; i < arr_sz; i+=4)
    {
        acc |= dst_u->w64[i] &= ~src_u->w64[i];
        acc |= dst_u->w64[i+1] &= ~src_u->w64[i+1];
        acc |= dst_u->w64[i+2] &= ~src_u->w64[i+2];
        acc |= dst_u->w64[i+3] &= ~src_u->w64[i+3];
    }
    return acc;
#endif
}



/*!
   \brief bitblock SUB operation. 

   \param dst - destination block.
   \param src - source block.

   \returns pointer on destination block. 
    If returned value  equal to src means that block mutation requested. 
    NULL is valid return value.

   @ingroup bitfunc
*/
inline 
bm::word_t* bit_operation_sub(bm::word_t* BMRESTRICT dst, 
                              const bm::word_t* BMRESTRICT src)
{
    BM_ASSERT(dst || src);

    bm::word_t* ret = dst;
    if (IS_VALID_ADDR(dst))  //  The destination block already exists
    {
        if (!IS_VALID_ADDR(src))
        {
            if (IS_FULL_BLOCK(src))
            {
                // If the source block is all set
                // just clean the destination block
                return 0;
            }
        }
        else
        {
            auto any = bm::bit_block_sub(dst, src);
            if (!any)
                ret = 0;
        }
    }
    else // The destination block does not exist yet
    {
        if (!IS_VALID_ADDR(src))
        {
            if (IS_FULL_BLOCK(src)) 
            {
                // The source block is full
                return 0;
            }
        }
        else
        {
            if (IS_FULL_BLOCK(dst))
            {
                // The only case when we have to allocate the new block:
                // dst is all set and src exists
                return const_cast<bm::word_t*>(src);                  
            }
        }
    }
    return ret;
}


/*!
   \brief Plain bitblock XOR operation. 
   Function does not analyse availability of source and destination blocks.

   \param dst - destination block.
   \param src - source block.

   @ingroup bitfunc
*/
inline 
void bit_block_xor(bm::word_t* BMRESTRICT dst, 
                   const bm::word_t* BMRESTRICT src)
{
#ifdef BMVECTOPT
    VECT_XOR_ARR(dst, src, src + bm::set_block_size);
#else
    const bm::wordop_t* BMRESTRICT wrd_ptr = (wordop_t*) src;
    const bm::wordop_t* BMRESTRICT wrd_end = 
                            (wordop_t*) (src + bm::set_block_size);
    bm::wordop_t* BMRESTRICT dst_ptr = (wordop_t*)dst;

    // Regular XOR operation on the whole block.
    do
    {
        dst_ptr[0] ^= wrd_ptr[0];
        dst_ptr[1] ^= wrd_ptr[1];
        dst_ptr[2] ^= wrd_ptr[2];
        dst_ptr[3] ^= wrd_ptr[3];

        dst_ptr+=4;
        wrd_ptr+=4;
    } while (wrd_ptr < wrd_end);
#endif
}

/*!
   \brief bitblock AND NOT with constant ~0 mask operation.

   \param dst - destination block.
   \param src - source block.

   @ingroup bitfunc
*/
inline
void bit_andnot_arr_ffmask(bm::word_t* BMRESTRICT dst,
                           const bm::word_t* BMRESTRICT src)
{
    const bm::word_t* BMRESTRICT src_end = src + bm::set_block_size;
#ifdef BMVECTOPT
    VECT_ANDNOT_ARR_2_MASK(dst, src, src_end, ~0u);
#else
    bm::wordop_t* dst_ptr = (wordop_t*)dst;
    const bm::wordop_t* wrd_ptr = (wordop_t*) src;
    const bm::wordop_t* wrd_end = (wordop_t*) src_end;

    do
    {
        dst_ptr[0] = bm::all_bits_mask & ~wrd_ptr[0];
        dst_ptr[1] = bm::all_bits_mask & ~wrd_ptr[1];
        dst_ptr[2] = bm::all_bits_mask & ~wrd_ptr[2];
        dst_ptr[3] = bm::all_bits_mask & ~wrd_ptr[3];
        dst_ptr+=4; wrd_ptr+=4;
    } while (wrd_ptr < wrd_end);
#endif
}

/*!
   \brief bitblock XOR operation. 

   \param dst - destination block.
   \param src - source block.

   \returns pointer on destination block. 
    If returned value  equal to src means that block mutation requested. 
    NULL is valid return value.

   @ingroup bitfunc
*/
inline 
bm::word_t* bit_operation_xor(bm::word_t* BMRESTRICT dst, 
                              const bm::word_t* BMRESTRICT src)
{
    BM_ASSERT(dst || src);
    if (src == dst) return 0;  // XOR rule  

    bm::word_t* ret = dst;

    if (IS_VALID_ADDR(dst))  //  The destination block already exists
    {           
        if (!src) return dst;
        
        bit_block_xor(dst, src);
    }
    else // The destination block does not exist yet
    {
        if (!src) return dst;      // 1 xor 0 = 1

        // Here we have two cases:
        // if dest block is full or zero if zero we need to copy the source
        // otherwise XOR loop against 0xFF...
        //BM_ASSERT(dst == 0);
        return const_cast<bm::word_t*>(src);  // src is the final result               
    }
    return ret;
}

/*!
   \brief Performs bitblock XOR operation and calculates bitcount of the result. 

   \param src1 - bit block start ptr
   \param src1_end - bit block end ptr
   \param src2 - second bit block

   \returns bitcount value 

   @ingroup bitfunc
*/
inline 
bm::id_t bit_operation_xor_count(const bm::word_t* BMRESTRICT src1,
                                 const bm::word_t* BMRESTRICT src2)
{
    if (IS_EMPTY_BLOCK(src1) || IS_EMPTY_BLOCK(src2))
    {
        if (IS_EMPTY_BLOCK(src1) && IS_EMPTY_BLOCK(src2))
            return 0;
        const bm::word_t* block = IS_EMPTY_BLOCK(src1) ? src2 : src1;
        return bit_block_count(block);
    }
    return bit_block_xor_count(src1, src2);
}

/*!
   \brief Performs bitblock XOR operation test. 

   \param src1 - bit block start ptr
   \param src2 - second bit block ptr

   \returns non zero value if there are bits

   @ingroup bitfunc
*/
inline 
bm::id_t bit_operation_xor_any(const bm::word_t* BMRESTRICT src1,
                               const bm::word_t* BMRESTRICT src2)
{
    if (IS_EMPTY_BLOCK(src1) || IS_EMPTY_BLOCK(src2))
    {
        if (IS_EMPTY_BLOCK(src1) && IS_EMPTY_BLOCK(src2))
            return 0;
        const bm::word_t* block = IS_EMPTY_BLOCK(src1) ? src2 : src1;
        return !bit_is_all_zero(block);
    }
    return bit_block_xor_any(src1, src2);
}



/**
    \brief Inspects block for full zero words 

    \param blk - bit block pointer
    \param data_size - data size

    \return size of all non-zero words

    @ingroup bitfunc
*/

template<class T>
unsigned bit_count_nonzero_size(const T*     blk, 
                                unsigned     data_size)
{
    BM_ASSERT(blk && data_size);
    unsigned count = 0;
    const T* blk_end = blk + data_size - 2;

    do
    {
        if (*blk == 0) 
        {
            // scan fwd to find 0 island length
            const T* blk_j = blk + 1;
            for (; blk_j < blk_end; ++blk_j)
            {
                if (*blk_j != 0)
                    break;
            }
            blk = blk_j-1;
            count += (unsigned)sizeof(gap_word_t);
        }
        else
        {
            // scan fwd to find non-0 island length
            const T* blk_j = blk + 1;
            for ( ; blk_j < blk_end; ++blk_j)
            {
                if (*blk_j == 0)
                {
                    // look ahead to identify and ignore short 0-run
                    if (blk_j[1] | blk_j[2])
                    {
                        // skip zero word
                        ++blk_j;
                        continue;
                    }
                    break;
                }
            }
            count += unsigned(sizeof(gap_word_t));
            // count all bit-words now
            count += unsigned(blk_j - blk) * unsigned(sizeof(T));
            blk = blk_j;
        }
        ++blk;
    }
    while(blk < blk_end); 

    return count + unsigned(2 * sizeof(T));
}


/**
    \brief Searches for the next 1 bit in the BIT block
    \param data - BIT buffer
    \param nbit - bit index to start checking from
    \param prev - returns previously checked value

    @ingroup bitfunc
*/
inline 
int bit_find_in_block(const bm::word_t* data, 
                      unsigned          nbit, 
                      bm::id_t*         prev)
{
    bm::id_t p = *prev;
    int found = 0;

    for(;;)
    {
        unsigned nword  = nbit >> bm::set_word_shift;
        if (nword >= bm::set_block_size)
            break;

        bm::word_t val = data[nword] >> (p & bm::set_word_mask);
        if (val)
        {
            unsigned trail_z = bm::word_trailing_zeros(val);
            if (trail_z)
            {
                val >>= trail_z;
                p += trail_z;
                BM_ASSERT(val & 1);
            }
            ++found;
            break;
        }
        else
        {
           p    += (bm::set_word_mask + 1) - (nbit & bm::set_word_mask);
           nbit += (bm::set_word_mask + 1) - (nbit & bm::set_word_mask);
        }
    }
    *prev = p;
    return found;
}

/*!
    \brief BIT block find the last set bit (backward search)

    \param block - bit block buffer pointer
    \param last - index of the last 1 bit (out)
    \return 0 is not found

    @ingroup bitfunc
*/
inline
unsigned bit_find_last(const bm::word_t* block, unsigned* last)
{
    BM_ASSERT(block);
    BM_ASSERT(last);

    for (unsigned i = bm::set_block_size-1; true; --i)
    {
        bm::word_t w = block[i];
        if (w)
        {
            unsigned idx = bm::bit_scan_reverse(w);
            *last = unsigned(idx + (i * 8u * sizeof(bm::word_t)));
            return w;
        }
        if (i == 0)
            break;
    } // for i
    return 0u;
}

/*!
    \brief BIT block find the first set bit

    \param block - bit block buffer pointer
    \param first - index of the first 1 bit (out)
    \return 0 if not foud

    @ingroup bitfunc
*/
inline
unsigned bit_find_first(const bm::word_t* block, unsigned* first)
{
    BM_ASSERT(block);
    BM_ASSERT(first);

    for (unsigned i = 0; i < bm::set_block_size; ++i)
    {
        bm::word_t w = block[i];
        if (w)
        {
            unsigned idx = bm::word_trailing_zeros(w);
            *first = unsigned(idx + (i * 8u * sizeof(bm::word_t)));
            return w;
        }
    } // for i
    return 0u;
}

/*!
    \brief BIT block find position for the rank

    \param block - bit block buffer pointer
    \param rank - rank to find (must be > 0)
    \param nbit_from - start bit position in block
    \param nbit_pos - found position
 
    \return 0 if position with rank was found, or
              the remaining rank (rank - population count)

    @ingroup bitfunc
*/
inline
bm::id_t bit_find_rank(const bm::word_t* const block,
                       bm::id_t          rank,
                       unsigned          nbit_from,
                       unsigned&         nbit_pos)
{
    BM_ASSERT(block);
    BM_ASSERT(rank);
    
    unsigned nword  = nbit_from >> bm::set_word_shift;
    BM_ASSERT(nword < bm::set_block_size);

    bm::id_t nbit = (nbit_from & bm::set_word_mask);
    unsigned pos = nbit_from;

    if (nbit)
    {
        bm::id_t w = block[nword];
        w >>= nbit;
        for (; w; w >>= 1u)
        {
            rank -= (w & 1u);
            if (!rank)
            {
                nbit_pos = pos;
                return rank;
            }
            ++nbit; ++pos;
        } // for
        pos += unsigned(32u - nbit);
        ++nword;
    }
    
    for (; nword < bm::set_block_size; ++nword)
    {
        bm::id_t w = block[nword];
        bm::id_t bc = bm::word_bitcount(w);
        if (rank > bc)
        {
            rank -= bc;
            pos += 32u;
            continue;
        }
        for (; w; w >>= 1u)
        {
            rank -= (w & 1u);
            if (!rank)
            {
                nbit_pos = pos;
                return rank;
            }
            ++pos;
        } // for
    } // for nword
    return rank;
}


/*!
   \brief Templated algorithm to unpacks octet based word into list of ON bit indexes
   \param w - value
   \param func - bit functor 

   @ingroup bitfunc
*/
template<typename T, typename F> 
void bit_for_each_4(T w, F& func)
{
    for (unsigned sub_octet = 0; w != 0; w >>= 4, sub_octet += 4)
    {
        switch (w & 15) // 1111
        {
        case 0: // 0000
            break;
        case 1: // 0001
            func(sub_octet);
            break;
        case 2: // 0010
            func(sub_octet + 1);
            break;
        case 3:	// 0011
            func(sub_octet, sub_octet + 1);
            break;
        case 4: // 0100
            func(sub_octet + 2);
            break;
        case 5: // 0101
            func(sub_octet, sub_octet + 2);
            break;
        case 6: // 0110
            func(sub_octet + 1, sub_octet + 2);
            break;
        case 7: // 0111
            func(sub_octet, sub_octet + 1, sub_octet + 2);
            break;
        case 8: // 1000
            func(sub_octet + 3);
            break;
        case 9: // 1001
            func(sub_octet, sub_octet + 3);
            break;
        case 10: // 1010
            func(sub_octet + 1, sub_octet + 3);
            break;
        case 11: // 1011
            func(sub_octet, sub_octet + 1, sub_octet + 3);
            break;
        case 12: // 1100
            func(sub_octet + 2, sub_octet + 3);
            break;
        case 13: // 1101
            func(sub_octet, sub_octet + 2, sub_octet + 3);
            break;
        case 14: // 1110
            func(sub_octet + 1, sub_octet + 2, sub_octet + 3);
            break;
        case 15: // 1111
            func(sub_octet, sub_octet + 1, sub_octet + 2, sub_octet + 3);
            break;
        default:
            BM_ASSERT(0);
            break;
        }
        
    } // for
}


/*!
   \brief Templated algorithm to unpacks word into list of ON bit indexes
   \param w - value
   \param func - bit functor 

   @ingroup bitfunc
*/
template<typename T, typename F> 
void bit_for_each(T w, F& func)
{
    // Note: 4-bit table method works slower than plain check approach
    for (unsigned octet = 0; w != 0; w >>= 8, octet += 8)
    {
        if (w & 1)   func(octet + 0);
        if (w & 2)   func(octet + 1);
        if (w & 4)   func(octet + 2);
        if (w & 8)   func(octet + 3);
        if (w & 16)  func(octet + 4);
        if (w & 32)  func(octet + 5);
        if (w & 64)  func(octet + 6);
        if (w & 128) func(octet + 7);
        
    } // for
}

/*! @brief Adaptor to copy 1 bits to array
    @internal
*/
template<typename B> class copy_to_array_functor
{
public:
    copy_to_array_functor(B* bits): bp_(bits)
    {}

    B* ptr() { return bp_; }
    
    void operator()(unsigned bit_idx) { *bp_++ = (B)bit_idx; }
    
    void operator()(unsigned bit_idx0, 
                    unsigned bit_idx1) 
    { 
        bp_[0] = (B)bit_idx0; bp_[1] = (B)bit_idx1;
        bp_+=2;
    }
    
    void operator()(unsigned bit_idx0, 
                    unsigned bit_idx1, 
                    unsigned bit_idx2) 
    { 
        bp_[0] = (B)bit_idx0; bp_[1] = (B)bit_idx1; bp_[2] = (B)bit_idx2;
        bp_+=3;
    }
    
    void operator()(unsigned bit_idx0, 
                    unsigned bit_idx1, 
                    unsigned bit_idx2, 
                    unsigned bit_idx3) 
    { 
        bp_[0] = (B)bit_idx0; bp_[1] = (B)bit_idx1;
        bp_[2] = (B)bit_idx2; bp_[3] = (B)bit_idx3;
        bp_+=4;
    }

private:
    copy_to_array_functor(const copy_to_array_functor&);
    copy_to_array_functor& operator=(const copy_to_array_functor&);
private:
    B* bp_;
};


/*!
   \brief Unpacks word into list of ON bit indexes (quad-bit based)
   \param w - value
   \param bits - pointer on the result array 
   \return number of bits in the list

   @ingroup bitfunc
*/
template<typename T,typename B> unsigned bit_list_4(T w, B* bits)
{
    copy_to_array_functor<B> func(bits);
    bit_for_each_4(w, func);
    return (unsigned)(func.ptr() - bits);
}

/*!
    \brief Unpacks word into list of ON bit indexes using popcnt method
    \param w - value
    \param bits - pointer on the result array
    \param offset - value to add to bit position (programmed shift)
    \return number of bits in the list

    @ingroup bitfunc
    @internal
*/
template<typename B>
unsigned short bitscan_popcnt(bm::id_t w, B* bits, unsigned short offs)
{
    unsigned pos = 0;
    while (w)
    {
        bm::id_t t = w & -w;
        bits[pos++] = (B)(bm::word_bitcount(t - 1) + offs);
        w &= w - 1;
    }
    return (unsigned short)pos;
}

/*!
    \brief Unpacks word into list of ON bit indexes using popcnt method
    \param w - value
    \param bits - pointer on the result array
    \return number of bits in the list
 
    @ingroup bitfunc
    @internal
*/
template<typename B>
unsigned short bitscan_popcnt(bm::id_t w, B* bits)
{
    unsigned pos = 0;
    while (w)
    {
        bm::id_t t = w & -w;
        bits[pos++] = (B)(bm::word_bitcount(t - 1));
        w &= w - 1;
    }
    return (unsigned short)pos;
}


/*!
    \brief Unpacks 64-bit word into list of ON bit indexes using popcnt method
    \param w - value
    \param bits - pointer on the result array
    \param offs - bit address offset to add (0 - default)
    \return number of bits in the list

    @ingroup bitfunc
*/
template<typename B>
unsigned short bitscan_popcnt64(bm::id64_t w, B* bits)
{
    unsigned short pos = 0;
    while (w)
    {
        bm::id64_t t = w & -w;
        bits[pos++] = (B) bm::word_bitcount64(t - 1);
        w &= w - 1;
    }
    return pos;
}

template<typename V, typename B>
unsigned short bitscan(V w, B* bits)
{
    if (bm::conditional<sizeof(V) == 8>::test())
    {
        return bm::bitscan_popcnt64(w, bits);
    }
    else
    {
        return bm::bitscan_popcnt((bm::word_t)w, bits);
    }
}




/*!
   \brief Unpacks word into list of ON bit indexes
   \param w - value
   \param bits - pointer on the result array 
   \return number of bits in the list

   @ingroup bitfunc
*/
template<typename T,typename B> unsigned bit_list(T w, B* bits)
{
    copy_to_array_functor<B> func(bits);
    bit_for_each(w, func);
    return (unsigned)(func.ptr() - bits);
}


/*!
    @brief Choose best representation for a bit-block
    @ingroup bitfunc 
*/
inline
bm::set_representation best_representation(unsigned bit_count,
                                           unsigned total_possible_bitcount,
                                           unsigned gap_count,
                                           unsigned block_size)
{
    unsigned arr_size = unsigned(sizeof(bm::gap_word_t) * bit_count + sizeof(bm::gap_word_t));
    unsigned gap_size = unsigned(sizeof(bm::gap_word_t) * gap_count + sizeof(bm::gap_word_t));
    unsigned inv_arr_size = unsigned(sizeof(bm::gap_word_t) * (total_possible_bitcount - bit_count) + sizeof(bm::gap_word_t));

    if ((gap_size < block_size) && (gap_size < arr_size) && (gap_size < inv_arr_size))
    {
        return bm::set_gap;
    }

    if (arr_size < inv_arr_size)
    {
        if ((arr_size < block_size) && (arr_size < gap_size))
        {
            return bm::set_array1;
        }
    }
    else
    {
        if ((inv_arr_size < block_size) && (inv_arr_size < gap_size))
        {
            return bm::set_array0;
        }
    }
    return bm::set_bitset;
}

/*!
    @brief Convert bit block into an array of ints corresponding to 1 bits
    @ingroup bitfunc 
*/
template<typename T> T bit_convert_to_arr(T* BMRESTRICT dest, 
                                          const unsigned* BMRESTRICT src, 
                                          bm::id_t bits, 
                                          unsigned dest_len,
                                          unsigned mask = 0)
{
    T* BMRESTRICT pcurr = dest;
    for (unsigned bit_idx=0; bit_idx < bits; ++src,bit_idx += unsigned(sizeof(*src) * 8))
    {
        unsigned val = *src ^ mask; // possible to invert value by XOR 0xFF..
        if (val == 0) 
        {
            continue;
        }
        if (pcurr + sizeof(val)*8 >= dest + dest_len) // insufficient space
        {
            return 0;
        }
        unsigned char b_list[64];
        unsigned word_bit_cnt  = bm::bitscan_popcnt(val, b_list);
        for (unsigned j = 0; j < word_bit_cnt; ++j)
        {
            *pcurr++ = (T)(b_list[j] + bit_idx);
        }
    } // for
    return (T)(pcurr - dest);
}

/**
    \brief Checks all conditions and returns true if block consists of only 0 bits
    \param blk - Blocks's pointer
    \param deep_scan - flag to do full bit block verification (scan)
                       when deep scan is not requested result can be approximate
    \returns true if all bits are in the block are 0

    @internal
*/
inline
bool check_block_zero(const bm::word_t* blk, bool  deep_scan)
{
    if (!blk) return true;

    bool ret;
    if (BM_IS_GAP(blk))
        ret = gap_is_all_zero(BMGAP_PTR(blk));
    else
        ret = deep_scan ? bm::bit_is_all_zero(blk) : 0;
    return ret;
}


/**
    \brief Checks if block has only 1 bits
    \param blk - Block's pointer
    \param deep_scan - flag to do full bit block verification (scan)
                       when deep scan is not requested result can be approximate
    \return true if block consists of 1 bits.

    @internal
*/
inline
bool check_block_one(const bm::word_t* blk, bool deep_scan)
{
    if (blk == 0) return false;

    if (BM_IS_GAP(blk))
        return bm::gap_is_all_one(BMGAP_PTR(blk), bm::gap_max_bits);

    if (IS_FULL_BLOCK(blk))
        return true;
    
    if (!deep_scan)
        return false; // block exists - presume it has 0 bits

    return bm::is_bits_one((wordop_t*)blk,
                           (wordop_t*)(blk + bm::set_block_size));
}



/*! @brief Calculates memory overhead for number of gap blocks sharing 
           the same memory allocation table (level lengths table).
    @ingroup gapfunc
*/
template<typename T> 
unsigned gap_overhead(const T* length, 
                      const T* length_end, 
                      const T* glevel_len)
{
    BM_ASSERT(length && length_end && glevel_len);

    unsigned overhead = 0;
    for (;length < length_end; ++length)
    {
        unsigned len = *length;
        int level = gap_calc_level(len, glevel_len);
        BM_ASSERT(level >= 0 && level < (int)bm::gap_levels);
        unsigned capacity = glevel_len[level];
        BM_ASSERT(capacity >= len);
        overhead += capacity - len;
    }
    return overhead;
}


/*! @brief Finds optimal gap blocks lengths.
    @param length - first element of GAP lengths array
    @param length_end - end of the GAP lengths array
    @param glevel_len - destination GAP lengths array
    @ingroup gapfunc
*/
template<typename T>
bool improve_gap_levels(const T* length,
                        const T* length_end,
                        T*       glevel_len)
{
    BM_ASSERT(length && length_end && glevel_len);

    size_t lsize = size_t(length_end - length);

    BM_ASSERT(lsize);
    
    gap_word_t max_len = 0;
    unsigned i;
    for (i = 0; i < lsize; ++i)
    {
        if (length[i] > max_len)
            max_len = length[i];
    }
    if (max_len < 5 || lsize <= bm::gap_levels)
    {
        glevel_len[0] = T(max_len + 4);
        for (i = 1; i < bm::gap_levels; ++i)
        {
            glevel_len[i] = bm::gap_max_buff_len;
        }
        return true;
    }

    glevel_len[bm::gap_levels-1] = T(max_len + 5);

    unsigned min_overhead = gap_overhead(length, length_end, glevel_len);
    bool is_improved = false;

    // main problem solving loop
    //
    for (i = bm::gap_levels-2; ; --i)
    {
        unsigned opt_len = 0;
        unsigned j;
        bool imp_flag = false;
        gap_word_t gap_saved_value = glevel_len[i];
        for (j = 0; j < lsize; ++j)
        {
            glevel_len[i] = T(length[j] + 4);
            unsigned ov = gap_overhead(length, length_end, glevel_len);
            if (ov <= min_overhead)
            {
                min_overhead = ov;                
                opt_len = length[j]+4;
                imp_flag = true;
            }
        }
        if (imp_flag) 
        {
            glevel_len[i] = (T)opt_len; // length[opt_idx]+4;
            is_improved = true;
        }
        else 
        {
            glevel_len[i] = gap_saved_value;
        }
        if (i == 0) 
            break;
    }
    
    // Remove duplicates
    //
    T val = *glevel_len;
    T* gp = glevel_len;
    T* res = glevel_len;
    for (i = 0; i < bm::gap_levels; ++i)
    {
        if (val != *gp)
        {
            val = *gp;
            *++res = val;
        }
        ++gp;
    }

    // Filling the "unused" part with max. possible value
    while (++res < (glevel_len + bm::gap_levels)) 
    {
        *res = bm::gap_max_buff_len;
    }

    return is_improved;

}



/**
    Bit-block get adapter, takes bitblock and represents it as a 
    get_32() accessor function
    \internal
*/
class bitblock_get_adapter
{
public:
    bitblock_get_adapter(const bm::word_t* bit_block) : b_(bit_block) {}
    
    BMFORCEINLINE
    bm::word_t get_32() { return *b_++; }
private:
    const bm::word_t*  b_;
};


/**
    Bit-block store adapter, takes bitblock and saves results into it
    \internal
*/
class bitblock_store_adapter
{
public:
    bitblock_store_adapter(bm::word_t* bit_block) : b_(bit_block) {}
    BMFORCEINLINE
    void push_back(bm::word_t w) { *b_++ = w; }
private:
    bm::word_t* b_;
};

/**
    Bit-block sum adapter, takes values and sums it
    /internal
*/
class bitblock_sum_adapter
{
public:
    bitblock_sum_adapter() : sum_(0) {}
    BMFORCEINLINE
    void push_back(bm::word_t w) { this->sum_+= w; }
    /// Get accumulated sum
    bm::word_t sum() const { return this->sum_; }
private:
    bm::word_t sum_;
};

/**
    Adapter to get words from a range stream 
    (see range serialized bit-block)
    \internal
*/
template<class DEC> class decoder_range_adapter
{
public: 
    decoder_range_adapter(DEC& dec, unsigned from_idx, unsigned to_idx)
    : decoder_(dec),
      from_(from_idx),
      to_(to_idx),
      cnt_(0)
    {}

    bm::word_t get_32()
    {
        if (cnt_ < from_ || cnt_ > to_)
        {    
            ++cnt_; return 0;
        }
        ++cnt_;
        return decoder_.get_32();
    }

private:
    DEC&     decoder_;
    unsigned from_;
    unsigned to_;
    unsigned cnt_;
};


/*!
    Abstract recombination algorithm for two bit-blocks
    Bit blocks can come as dserialization decoders or bit-streams
*/
template<class It1, class It2, class BinaryOp, class Encoder>
void bit_recomb(It1& it1, It2& it2, 
                BinaryOp& op, 
                Encoder& enc, 
                unsigned block_size = bm::set_block_size)
{
    for (unsigned i = 0; i < block_size; ++i)
    {
        bm::word_t w1 = it1.get_32();
        bm::word_t w2 = it2.get_32();
        bm::word_t w = op(w1, w2);
        enc.push_back( w );
    } // for
}

/// Bit AND functor
template<typename W> struct bit_AND
{
    W operator()(W w1, W w2) { return w1 & w2; }
};

/// Bit OR functor
template<typename W> struct bit_OR
{
    W operator()(W w1, W w2) { return w1 | w2; }
};

/// Bit SUB functor
template<typename W> struct bit_SUB
{
    W operator()(W w1, W w2) { return w1 & ~w2; }
};

/// Bit XOR functor
template<typename W> struct bit_XOR
{
    W operator()(W w1, W w2) { return w1 ^ w2; }
};

/// Bit ASSIGN functor
template<typename W> struct bit_ASSIGN
{
    W operator()(W, W w2) { return w2; }
};

/// Bit COUNT functor
template<typename W> struct bit_COUNT
{
    W operator()(W w1, W w2) 
    {
        w1 = 0;
        BM_INCWORD_BITCOUNT(w1, w2);
        return w1;
    }
};

/// Bit COUNT AND functor
template<typename W> struct bit_COUNT_AND
{
    W operator()(W w1, W w2) 
    {
        W r = 0;
        BM_INCWORD_BITCOUNT(r, w1 & w2);
        return r;
    }
};

/// Bit COUNT XOR functor
template<typename W> struct bit_COUNT_XOR
{
    W operator()(W w1, W w2) 
    {
        W r = 0;
        BM_INCWORD_BITCOUNT(r, w1 ^ w2);
        return r;
    }
};

/// Bit COUNT OR functor
template<typename W> struct bit_COUNT_OR
{
    W operator()(W w1, W w2) 
    {
        W r = 0;
        BM_INCWORD_BITCOUNT(r, w1 | w2);
        return r;
    }
};


/// Bit COUNT SUB AB functor
template<typename W> struct bit_COUNT_SUB_AB
{
    W operator()(W w1, W w2) 
    {
        W r = 0;
        BM_INCWORD_BITCOUNT(r, w1 & (~w2));
        return r;
    }
};

/// Bit SUB BA functor
template<typename W> struct bit_COUNT_SUB_BA
{
    W operator()(W w1, W w2) 
    {
        W r = 0;
        BM_INCWORD_BITCOUNT(r, w2 & (~w1));
        return r;
    }
};

/// Bit COUNT A functor
template<typename W> struct bit_COUNT_A
{
    W operator()(W w1, W )
    {
        W r = 0;
        BM_INCWORD_BITCOUNT(r, w1);
        return r;
    }
};

/// Bit COUNT B functor
template<typename W> struct bit_COUNT_B
{
    W operator()(W, W w2)
    {
        W r = 0;
        BM_INCWORD_BITCOUNT(r, w2);
        return r;
    }
};

typedef 
void (*gap_operation_to_bitset_func_type)(unsigned*, 
                                          const gap_word_t*);

typedef 
gap_word_t* (*gap_operation_func_type)(const gap_word_t* BMRESTRICT,
                                       const gap_word_t* BMRESTRICT,
                                       gap_word_t*       BMRESTRICT,
                                       unsigned& );

typedef
bm::id_t (*bit_operation_count_func_type)(const bm::word_t* BMRESTRICT,
                                          const bm::word_t* BMRESTRICT);


template<bool T> 
struct operation_functions
{
    static 
        gap_operation_to_bitset_func_type gap2bit_table_[bm::set_END];
    static 
        gap_operation_func_type gapop_table_[bm::set_END];
    static
        bit_operation_count_func_type bit_op_count_table_[bm::set_END];

    static
    gap_operation_to_bitset_func_type gap_op_to_bit(unsigned i)
    {
        return gap2bit_table_[i];
    }

    static
    gap_operation_func_type gap_operation(unsigned i)
    {
        return gapop_table_[i];
    }

    static
    bit_operation_count_func_type bit_operation_count(unsigned i)
    {
        return bit_op_count_table_[i];
    }
};

template<bool T>
gap_operation_to_bitset_func_type 
operation_functions<T>::gap2bit_table_[bm::set_END] = {
    &gap_and_to_bitset<bm::gap_word_t>,    // set_AND
    &gap_add_to_bitset<bm::gap_word_t>,    // set_OR
    &gap_sub_to_bitset<bm::gap_word_t>,    // set_SUB
    &gap_xor_to_bitset<bm::gap_word_t>,    // set_XOR
    0
};

template<bool T>
gap_operation_func_type 
operation_functions<T>::gapop_table_[bm::set_END] = {
    &gap_operation_and,    // set_AND
    &gap_operation_or,     // set_OR
    &gap_operation_sub,    // set_SUB
    &gap_operation_xor,    // set_XOR
    0
};


template<bool T>
bit_operation_count_func_type 
operation_functions<T>::bit_op_count_table_[bm::set_END] = {
    0,                            // set_AND
    0,                            // set_OR
    0,                            // set_SUB
    0,                            // set_XOR
    0,                            // set_ASSIGN
    0,                            // set_COUNT
    &bit_operation_and_count,     // set_COUNT_AND
    &bit_operation_xor_count,     // set_COUNT_XOR
    &bit_operation_or_count,      // set_COUNT_OR
    &bit_operation_sub_count,     // set_COUNT_SUB_AB
    &bit_operation_sub_count_inv, // set_COUNT_SUB_BA
    0,                            // set_COUNT_A
    0,                            // set_COUNT_B
};


const unsigned short set_bitscan_wave_size = 2;
/*!
    \brief Unpacks word wave (Nx 32-bit words)
    \param w_ptr - pointer on wave start
    \param bits - pointer on the result array
    \return number of bits in the list

    @ingroup bitfunc
    @internal
*/
inline
unsigned short bitscan_wave(const bm::word_t* w_ptr, unsigned char* bits)
{
    bm::word_t w0, w1;
    unsigned short cnt0;

    w0 = w_ptr[0];
    w1 = w_ptr[1];
    
#if defined(BMAVX2OPT) || defined(BMSSE42OPT)
    // combine into 64-bit word and scan (when HW popcnt64 is available)
    bm::id64_t w = (bm::id64_t(w1) << 32) | w0;
    cnt0 = (unsigned short) bm::bitscan_popcnt64(w, bits);
#else
    unsigned short cnt1;
    // decode wave as two 32-bit bitscan decodes
    cnt0 = w0 ? bm::bitscan_popcnt(w0, bits) : 0;
    cnt1 = w1 ? bm::bitscan_popcnt(w1, bits + cnt0, 32) : 0;
    cnt0 = (unsigned short)(cnt0 + cnt1);
#endif
    return cnt0;
}


/**
    bit index to word gather-scatter algorithm
    @ingroup bitfunc
    @internal
*/
template<typename TRGW, typename IDX, typename SZ>
void bit_block_gather_scatter(TRGW* arr, const bm::word_t* blk,
                              const IDX* idx, SZ size, unsigned start, unsigned bit_idx)
{
#if defined(BM64_SSE4)
    // optimized for unsigned
    if (bm::conditional<sizeof(TRGW)==4 && sizeof(IDX)==4>::test())
    {
        sse4_bit_block_gather_scatter(arr, blk, idx, size, start, bit_idx);
        return;
    }
#elif defined(BM64_AVX2)
    if (bm::conditional<sizeof(TRGW)==4 && sizeof(IDX)==4>::test())
    {
        avx2_bit_block_gather_scatter(arr, blk, idx, size, start, bit_idx);
        return;
    }
#endif

    const unsigned len = (size - start);
    const unsigned len_unr = len - (len % 2);
    unsigned k;
    for (k = 0; k < len_unr; k+=2)
    {
        const unsigned base = start + k;

        const unsigned nbitA = unsigned(idx[base] & bm::set_block_mask);
        arr[base]   |= TRGW(bool(blk[nbitA >> bm::set_word_shift] & (1u << (nbitA & bm::set_word_mask))) << bit_idx); 
        const unsigned nbitB = unsigned(idx[base + 1] & bm::set_block_mask);
        arr[base+1] |= TRGW(bool(blk[nbitB >> bm::set_word_shift] & (1u << (nbitB & bm::set_word_mask))) << bit_idx);
    }
    for (; k < len; ++k)
    {
        unsigned nbit = unsigned(idx[start + k] & bm::set_block_mask);
        arr[start + k] |= TRGW(bool(blk[nbit >> bm::set_word_shift] & (1u << (nbit & bm::set_word_mask))) << bit_idx);
    }
}

/**
    block boundaries look ahead
    @internal
*/
inline
unsigned idx_arr_block_lookup(const unsigned* idx, unsigned size, unsigned nb, unsigned start)
{
    BM_ASSERT(idx);
    
    if (nb == unsigned(idx[size-1] >> bm::set_block_shift))
        return size;
    
#if defined(BMAVX2OPT)
    return avx2_idx_arr_block_lookup(idx, size, nb, start);
#elif defined(BMSSE42OPT)
    return sse4_idx_arr_block_lookup(idx, size, nb, start);
#else
    // use 64-bit variable for parallel compare (SIMD style)
    if (bm::conditional< sizeof(void*)==8 >::test())
    {
        const unsigned len = (size - start);
        const unsigned len_unr = len - (len % 2);
        
        idx += start;
        bm::id64_t nb64 = (nb << 31) | nb;
        
        unsigned k;
        for (k = 0; k < len_unr; k+=2)
        {
            bm::id64_t i64 =
                ((idx[k] >> bm::set_block_shift) << 31) | (idx[k+1] >> bm::set_block_shift);
            if (nb64 != i64)
                break;
        }
        for (; k < len; ++k)
        {
            if (nb != unsigned(idx[k] >> bm::set_block_shift))
                break;
        }
        return start + k;
    }
    
    for (;(start < size) &&
          (nb == unsigned(idx[start] >> bm::set_block_shift)); ++start)
    {}
    return start;
#endif
}


// --------------------------------------------------------------
// Functions to work with int values stored in 64-bit pointers
// --------------------------------------------------------------

/*!
    \brief helper union to interpret pointer as integers 
    @internal
*/
union ptr_payload_t
{
    bm::word_t* blk;
    bm::id64_t  i64;
    unsigned short i16[4];
};

/*!
    Test presense of value in payload pointer
    @internal
*/
inline
bm::id64_t ptrp_test(ptr_payload_t ptr, bm::gap_word_t v)
{
    if (v == 0)
    {
        return (ptr.i16[1] == 0);
    }
    bm::id64_t r = (ptr.i16[1] == v) | (ptr.i16[2] == v) | (ptr.i16[3] == v);
    return r;
}

// --------------------------------------------------------------------------

} // namespace bm

#endif
