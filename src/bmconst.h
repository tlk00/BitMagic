#ifndef BMCONST__H__INCLUDED__
#define BMCONST__H__INCLUDED__
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

/*! \file bmconst.h
    \brief Various constants and tables
*/

namespace bm
{

#if defined(_WIN32) || defined (_WIN64)
typedef unsigned __int64 id64_t;
#else
typedef unsigned long long int id64_t;
#endif

typedef unsigned int   id_t;
typedef unsigned int   word_t;
typedef unsigned short short_t;

#ifndef BM_DEFAULT_POOL_SIZE
# define BM_DEFAULT_POOL_SIZE 4096
#endif


const unsigned id_max = 0xFFFFFFFF;

// Data Block parameters

const unsigned set_block_size  = 2048u;
const unsigned set_block_shift = 16u;
const unsigned set_block_mask  = 0xFFFFu;
const unsigned set_blkblk_mask = 0xFFFFFFu;

const unsigned set_block_plain_size = set_block_size / 32u;
const unsigned set_block_plain_cnt = (unsigned)(sizeof(bm::word_t) * 8u);

// Word parameters

const unsigned set_word_shift = 5u;
const unsigned set_word_mask  = 0x1Fu;


// GAP related parameters.

typedef unsigned short gap_word_t;

const unsigned gap_max_buff_len = 1280;
const unsigned gap_max_bits = 65536;
const unsigned gap_equiv_len = (unsigned)
   ((sizeof(bm::word_t) * bm::set_block_size) / sizeof(gap_word_t));
const unsigned gap_levels = 4;
const unsigned gap_max_level = bm::gap_levels - 1;


// Block Array parameters

const unsigned set_array_size = 256u;
const unsigned set_array_shift = 8u;
const unsigned set_array_mask  = 0xFFu;
const unsigned set_total_blocks = (bm::set_array_size * bm::set_array_size);

const unsigned bits_in_block = bm::set_block_size * (unsigned)(sizeof(bm::word_t) * 8);
const unsigned bits_in_array = bm::bits_in_block * bm::set_array_size;


#if defined(BM64OPT) || defined(BM64_SSE4)

typedef id64_t  wordop_t;
const id64_t    all_bits_mask = 0xffffffffffffffff;
const unsigned set_block_size_op  = bm::set_block_size / 2;

#else

typedef word_t wordop_t;
const word_t all_bits_mask = 0xffffffff;
const unsigned set_block_size_op  = bm::set_block_size;

#endif


/*!
   @brief Block allocation strategies.
   @ingroup bvector
*/
enum strategy
{
    BM_BIT = 0, //!< No GAP compression strategy. All new blocks are bit blocks.
    BM_GAP = 1  //!< GAP compression is ON.
};

/**
    Codes of set operations
    @ingroup bvector
*/
enum set_operation
{
    set_AND         = 0,
    set_OR          = 1,
    set_SUB         = 2,
    set_XOR         = 3,
    set_ASSIGN      = 4,
    set_COUNT       = 5,
    set_COUNT_AND   = 6,
    set_COUNT_XOR   = 7,
    set_COUNT_OR    = 8,
    set_COUNT_SUB_AB= 9,
    set_COUNT_SUB_BA= 10,
    set_COUNT_A     = 11,
    set_COUNT_B     = 12,

    set_END
};


/*!
   @brief Sort order declaration
   @ingroup bvector
*/
enum sort_order
{
    BM_UNSORTED = 0,      //!< input set is NOT sorted
    BM_SORTED = 1,        //!< input set is sorted (ascending order)
    BM_SORTED_UNIFORM = 2,//!< input set is sorted and belongs to one address block
    BM_UNKNOWN = 3        //!< sort order unknown
};


/*!
    @brief set representation variants
    @internal
*/
enum set_representation
{
    set_bitset  = 0,  //!< Simple bitset
    set_gap     = 1,  //!< GAP-RLE compression
    set_array1  = 2,  //!< array of set 1 values
    set_array0  = 3   //!< array of 0 values
};

/*!
   @brief NULL-able value support
   @ingroup bvector
*/
enum null_support
{
    use_null = 0, //!< support "non-assigned" or "NULL" logic
    no_null  = 1   //!< do not support NULL values
};


/**
    Internal structure. Copyright information.
    @internal
*/
template<bool T> struct _copyright
{
    static const char _p[];
    static const unsigned _v[3];
};

template<bool T> const char _copyright<T>::_p[] = 
    "BitMagic C++ Library. v.3.12.5 (c) 2002-2018 Anatoliy Kuznetsov.";
template<bool T> const unsigned _copyright<T>::_v[3] = {3, 12, 5};


template<bool T> struct DeBruijn_bit_position
{
    static const unsigned _multiply[32];
};

/**
    DeBruijn majic table
    @internal
*/
template<bool T>
const unsigned DeBruijn_bit_position<T>::_multiply[32] = { 
  0, 1, 28, 2, 29, 14, 24, 3, 30, 22, 20, 15, 25, 17, 4, 8, 
  31, 27, 13, 23, 21, 19, 16, 7, 26, 12, 18, 6, 11, 5, 10, 9
};

/**
    Structure keeps index of first right 1 bit for every byte.
    @ingroup bitfunc 
*/
template<bool T> struct first_bit_table
{
    static const signed char _idx[256];
};

template<bool T>
const signed char first_bit_table<T>::_idx[256] = {
  -1, 0, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
    6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
    6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
    6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
    6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
};

//---------------------------------------------------------------------

/** Structure to aid in counting bits
    table contains count of bits in 0-255 diapason of numbers

   @ingroup bitfunc
*/
template<bool T> struct bit_count_table 
{
  static const unsigned char _count[256];
};

template<bool T>
const unsigned char bit_count_table<T>::_count[256] = {
    0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4,1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
    1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
    1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
    2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
    1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
    2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
    2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
    3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8
};

//---------------------------------------------------------------------

/** Structure keeps all-left/right ON bits masks. 
    @ingroup bitfunc 
*/
template<bool T> struct block_set_table
{
    static const unsigned _left[32];
    static const unsigned _right[32];
};

template<bool T>
const unsigned block_set_table<T>::_left[32] = {
    0x1, 0x3, 0x7, 0xf, 0x1f, 0x3f, 0x7f, 0xff, 0x1ff, 0x3ff, 0x7ff,
    0xfff, 0x1fff, 0x3fff, 0x7fff, 0xffff, 0x1ffff, 0x3ffff, 0x7ffff,
    0xfffff, 0x1fffff, 0x3fffff, 0x7fffff, 0xffffff, 0x1ffffff, 0x3ffffff,
    0x7ffffff, 0xfffffff, 0x1fffffff, 0x3fffffff, 0x7fffffff, 0xffffffff
};

template<bool T>
const unsigned block_set_table<T>::_right[32] = {
    0xffffffff, 0xfffffffe, 0xfffffffc, 0xfffffff8, 0xfffffff0,
    0xffffffe0, 0xffffffc0, 0xffffff80, 0xffffff00, 0xfffffe00,
    0xfffffc00, 0xfffff800, 0xfffff000, 0xffffe000, 0xffffc000,
    0xffff8000, 0xffff0000, 0xfffe0000, 0xfffc0000, 0xfff80000,
    0xfff00000, 0xffe00000, 0xffc00000, 0xff800000, 0xff000000,
    0xfe000000, 0xfc000000, 0xf8000000, 0xf0000000, 0xe0000000,
    0xc0000000, 0x80000000
};

//---------------------------------------------------------------------



/*! @brief Default GAP lengths table.
    @ingroup gapfunc
*/
template<bool T> struct gap_len_table
{
    static const gap_word_t _len[bm::gap_levels];
};

template<bool T>
const gap_word_t gap_len_table<T>::_len[bm::gap_levels] = 
                { 128, 256, 512, bm::gap_max_buff_len }; 


/*! @brief Alternative GAP lengths table. 
    Good for for memory saver mode and very sparse bitsets.

    @ingroup gapfunc
*/
template<bool T> struct gap_len_table_min
{
    static const gap_word_t _len[bm::gap_levels];
};

template<bool T>
const gap_word_t gap_len_table_min<T>::_len[bm::gap_levels] = 
                                { 32, 96, 128, 512 }; 


/*! @brief Non-linear size growth GAP lengths table.
    @ingroup gapfunc
*/
template<bool T> struct gap_len_table_nl
{
    static const gap_word_t _len[bm::gap_levels];
};

template<bool T>
const gap_word_t gap_len_table_nl<T>::_len[bm::gap_levels] =
                { 32, 128, 512, bm::gap_max_buff_len };

/*!
    @brief codes for supported SIMD optimizations
*/
enum simd_codes
{
    simd_none  = 0,   ///!< No SIMD or any other optimization
    simd_sse2  = 1,   ///!< Intel SSE2
    simd_sse42 = 2,   ///!< Intel SSE4.2
    simd_avx2  = 5    ///!< Intel AVX2
};



} // namespace

#endif

