#ifndef LIBBM_INCLUDED_H__
#define LIBBM_INCLUDED_H__
/*
BitMagic Library License

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

#include <stddef.h>

/* Error codes */

/* General purpose codes */
#define BM_OK (0)
#define BM_ERR_BADALLOC (1)
#define BM_ERR_BADARG (2)
#define BM_ERR_RANGE (3)
#define BM_ERR_CPU   (4)
#define BM_ERR_SERIALFORMAT (5)
#define BM_ERR_BAD_VALUE (6)
#define BM_ERR_RANK_SELECT_IDX_MISSING (7)

/* Error codes for Java/JNI incapsulation */
#define BM_ERR_DETACHED (101)
#define BM_ERR_JVM_NOT_SUPPORTED (102)
#define BM_ERR_JVM_OUT_OF_MEMORY (103)

/*
    error codes and messages
*/
#define BM_OK_MSG           "BM-00: All correct"
#define BM_ERR_BADALLOC_MSG "BM-01: Allocation error"
#define BM_ERR_BADARG_MSG   "BM-02: Invalid or missing function argument"
#define BM_ERR_RANGE_MSG    "BM-03: Incorrect range or index"
#define BM_ERR_CPU_MSG      "BM-04: Incorrect CPU vectorization (SIMD) version"
#define BM_ERR_SERIALFORMAT_MSG "BM-05: Serialization format error"
#define BM_ERR_BAD_VALUE_MSG "BM-06: Bad value"
#define BM_ERR_RANK_SELECT_IDX_MISSING_MSG "BM-07: Rank-Select index not constructed, call sync() first"



#define BM_ERR_DETACHED_MSG    "BM-101: Current thread no attached to JVM"
#define BM_ERR_JVM_NOT_SUPPORTED_MSG    "BM-102: JVM version not supported"
#define BM_ERR_JVM_OUT_OF_MEMORY_MSG    "BM-103: Out of memory error"

#define BM_UNK_MSG          "BM-XX: Unknown error"


/*
    List of supported SIMD versions
*/
#define BM_SIMD_NO    0
#define BM_SIMD_SSE2  1
#define BM_SIMD_SSE42 2
#define BM_SIMD_AVX2  5



/* bit-vector handle */
#define BM_BVHANDLE void*
/* bit-vector enumerator handle */
#define BM_BVEHANDLE void*


/* arguments codes and values */
#define BM_TRUE 1
#define BM_FALSE 0

// -------------------------------------------------------
// Windows DLL import/export definitions
// to enable DLL configuration #define "BM_USE_DLL" 
// otherwise it assumes static linkage
//
#if defined(_WIN32) 
#ifdef BMDLLEXPORTS
#    define BM_API_EXPORT __declspec(dllexport)
#else
#if defined(BM_USE_DLL)
#    define BM_API_EXPORT __declspec(dllimport)
#else 
#   define BM_API_EXPORT
#endif
#endif
#else
#   define BM_API_EXPORT
#endif


#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/*  bit vector statistics used for serialization and memory management
*/
struct BM_bvector_statistics
{
    size_t bit_blocks;
    size_t gap_blocks;
    size_t max_serialize_mem;
    size_t memory_used;
};


/* -------------------------------------------- */
/* General purpose functions                    */
/* -------------------------------------------- */

/* Initialize libbm runtime before use*/
BM_API_EXPORT int BM_init(void*);

/**
    return copyright info string and version information.
*/
BM_API_EXPORT const char* BM_version(int* major, int* minor, int* patch);

/**
    return SIMD version used to build binaries
    one of BM_SIMD_* defines
*/
BM_API_EXPORT int BM_simd_version(void);


/**
    return error message by code
*/
BM_API_EXPORT const char* BM_error_msg(int errcode);

/* ------------------------------------------------- */
/* bvector construction, swap and sizing methods     */
/* ------------------------------------------------- */

/* construction and setters                          */

/* construct bvector handle 
   bv_max - maximum number of allowed bits (if 0 - allows maximum)
   (it is recommened to set size to maximum always and do not use size params)
*/
BM_API_EXPORT int BM_bvector_construct(BM_BVHANDLE* h, unsigned int bv_max);

/* init bvector handle to finilize construction
    This step is optional, unless we need to call "_no_check()" functions, which bypass
    certain checks
*/
BM_API_EXPORT int BM_bvector_init(BM_BVHANDLE h);


/* construct bvector handle as a copy
   hfrom - another handle to copy from
*/
BM_API_EXPORT int BM_bvector_construct_copy(BM_BVHANDLE* h, BM_BVHANDLE hfrom);

/* construct bvector handle as a read-only copy
   hfrom - another handle to copy from
*/
BM_API_EXPORT int BM_bvector_construct_copy_ro(BM_BVHANDLE* h, BM_BVHANDLE hfrom);

/* construct bvector handle as a read-write copy
   hfrom - another handle to copy from
*/
BM_API_EXPORT int BM_bvector_construct_copy_rw(BM_BVHANDLE* h, BM_BVHANDLE hfrom);


/* destroy bvector handle */
BM_API_EXPORT int BM_bvector_free(BM_BVHANDLE h);

/* get bit vector size
   psize - return size the bit vector
*/
BM_API_EXPORT int BM_bvector_get_size(BM_BVHANDLE h, unsigned int* psize);

/* get bit vector capacity
   pcap - return caapacity the bit vector
*/
/*int BM_bvector_get_capacity(BM_BVHANDLE h, unsigned int* pcap);*/


/* resize bit vector
   new_size - new requested size
*/
BM_API_EXPORT int BM_bvector_set_size(BM_BVHANDLE h, unsigned int new_size);

/* swap two bit-vectors
*/
BM_API_EXPORT int BM_bvector_swap(BM_BVHANDLE h1, BM_BVHANDLE h2);


/* -------------------------------------------- */
/* bvector functions to set and clear bits      */
/* -------------------------------------------- */


/* set bit to 1 or 0 
   i - index of a bit to set
   val - value (0 | 1)
*/
BM_API_EXPORT int BM_bvector_set_bit(BM_BVHANDLE h, unsigned int i, int val);

/* set list of bits to 1
   idx - index of bits to set
   idx_size - size of the array to set
*/ 
BM_API_EXPORT int BM_bvector_set_bits(BM_BVHANDLE h, unsigned int* idx, unsigned int idx_size);


/* set bit to 1 without extra checks (faster).
   Use of this function requires full initialization by BM_bvector_init(); 
   i - index of a bit to set
*/
BM_API_EXPORT int BM_bvector_set_bit_no_check(BM_BVHANDLE h, unsigned int i);


/* set bit only if current value equals the condition
   i - index of a bit to set
   val - value (0 | 1)
   condition - expected current value of bit i
   pchanged - optional return value, if bit was actually changed
*/
BM_API_EXPORT
int BM_bvector_set_bit_conditional(BM_BVHANDLE  h,
                                   unsigned int i,
                                   int          val,
                                   int          condition,
                                   int*         pchanged);
    
/* flip bit
   i - index of a bit to flip
*/
BM_API_EXPORT 
int BM_bvector_flip_bit(BM_BVHANDLE h, unsigned int i);
    
/* inc bit at position
   i           - index of a bit to flip
   1 + 0 = 1 (no carry over)
   1 + 1 = 0 (1 carry over)
   carry_over  - carry over bit
*/
BM_API_EXPORT 
int BM_bvector_inc_bit(BM_BVHANDLE h, unsigned int i, int* carry_over);


/* set all bits to 1
*/
BM_API_EXPORT int BM_bvector_set(BM_BVHANDLE h);



/* Sets all bits in the specified closed interval [left,right]
   This method DOES NOT resize vector.

   left  - interval start
   right - interval end (closed interval)
   value - value to set interval in
 
*/
BM_API_EXPORT int BM_bvector_set_range(BM_BVHANDLE h,
                         unsigned int left,
                         unsigned int right,
                         int          value);
    
/* invert all bits in the bit vector
*/
BM_API_EXPORT int BM_bvector_invert(BM_BVHANDLE h);
    


/* set all bits to 0 and (optionally) free unused memory
    free_mem - flag to release unused memory
*/
BM_API_EXPORT int BM_bvector_clear(BM_BVHANDLE h, int free_mem);

/* find 1 bit index in the vector and set it to 0

   i - index of bit to search from
   pnext - return index of the next 1 bit. 0 - means no more 1 bits.
*/
BM_API_EXPORT int BM_bvector_extract_next(BM_BVHANDLE h, unsigned int i, unsigned int* pnext);



/* -------------------------------------------- */
/* bvector functions to read bits               */
/* -------------------------------------------- */


/* get bit value */
BM_API_EXPORT int BM_bvector_get_bit(BM_BVHANDLE h, unsigned int i, int* pval);


/* bitcount
   pcount - return number of ON bits in the vector
*/
BM_API_EXPORT int BM_bvector_count(BM_BVHANDLE h, unsigned int* pcount);

/* range bitcount
   left  - interval start
   right - interval end (closed interval)
   pcount - return number of ON bits in the vector
*/
BM_API_EXPORT 
int BM_bvector_count_range(BM_BVHANDLE   h,
                           unsigned int  left,
                           unsigned int  right,
                           unsigned int* pcount);

/* check if there are any bits set 
   pval - return non-zero value if any bits are ON
*/
BM_API_EXPORT int BM_bvector_any(BM_BVHANDLE h, int* pval);

/* Finds index of 1 bit starting from position
   from - initial search position
   ppos - found position of 1 bit (>= from)
   pfound - 0 if nothing found
*/
BM_API_EXPORT int BM_bvector_find(BM_BVHANDLE h,
                    unsigned int from, unsigned int* ppos, int* pfound);

/* Finds index of 1 bit starting from the end of the vector
    ppos - found position of 1 bit (from the end)
    pfound - 0 if nothing found
*/
BM_API_EXPORT int BM_bvector_find_reverse(BM_BVHANDLE h,
                     unsigned int* ppos, int* pfound);

/* Finds bit-vector position (index) with specified rank
   starting from another poistion (if 0 rank is absolute).
 
   rank - search rank value
   from - bit-vector position to find from
   pidx - (out) position with relative rank
   pfound - 0 if nothing found
*/
BM_API_EXPORT int BM_bvector_find_rank(BM_BVHANDLE h,
                                       unsigned int rank,
                                       unsigned int from,
                                       unsigned int* pidx,
                                       int* pfound);


/* find first 1 bit index in the vector

   pi - return index of first bit 
   found - return 0 if first bit not found (empty vector)
*/
BM_API_EXPORT int BM_bvector_get_first(BM_BVHANDLE h, unsigned int* pi,
                                       int* pfound);

/* find 1 bit index in the vector

   i - index of bit to search from
   pnext - return index of the next 1 bit. 0 - means no more 1 bits.
*/
BM_API_EXPORT int BM_bvector_get_next(BM_BVHANDLE h, unsigned int i,
                                      unsigned int* pnext);



/* -------------------------------------------- */
/* bvector operations                           */
/* -------------------------------------------- */

/* Lexicographical comparison of two bit vectors
   pres - returns -1 if h1 less than h2, 1 - greater, 0 - equal.
*/
BM_API_EXPORT int BM_bvector_compare(BM_BVHANDLE h1, BM_BVHANDLE h2, int* pres);

/* Find first mismatch between two vectors
   pi - return found mismatch index
   pfound - return 0 if not found (vectors are identical)
*/
BM_API_EXPORT int BM_bvector_find_first_mismatch(BM_BVHANDLE h1, BM_BVHANDLE h2,
                 unsigned int* pi,
                 int* pfound);



/* Perform bit vector memory optimization
   opt_mode - optimization level:
    (0 - default, 1 - free empty blocks, 2 - free empty and full blocks, 3 - GAP compress)
   pstat - optional post optimization statistics
*/
BM_API_EXPORT
int BM_bvector_optimize(BM_BVHANDLE h,
                        int opt_mode,
                        struct BM_bvector_statistics* pstat);

/* freeze bit-vector into read-only mode. bit-vector becomes immutable.
*/
BM_API_EXPORT int BM_bvector_freeze(BM_BVHANDLE h);

/* Check if bit-vector is read-only.
   bit-vector becomes read-only after BM_bvector_freeze()
*/
BM_API_EXPORT int BM_bvector_is_ro(BM_BVHANDLE h, int* pval);


/* Perform calculate bit vector statistics
   pstat - bit vector statistics
*/
BM_API_EXPORT
int BM_bvector_calc_stat(BM_BVHANDLE h,
                         struct BM_bvector_statistics* pstat);

/* perform logical operation on two bit vectors
   hdst = hdst {OR/AND/XOR/SUB} hsrc
   opcode - operation code 
       AND - 0
       OR  - 1
       SUB - 2
       XOR - 3
*/
BM_API_EXPORT int BM_bvector_combine_operation(BM_BVHANDLE hdst, BM_BVHANDLE hsrc, int opcode);

/* perform logical AND operation on two bit vectors
   hdst = hdst AND hsrc
*/
BM_API_EXPORT int BM_bvector_combine_AND(BM_BVHANDLE hdst, BM_BVHANDLE hsrc);

/* perform 3-operand logical AND (set intersect) on two source bit vectors
   hdst = hsrc1 AND hsrc2
   compress - 0 (do not compress, 1 - compress)
*/
BM_API_EXPORT int BM_bvector_combine_AND_2sc(BM_BVHANDLE hdst, BM_BVHANDLE hsrc1, BM_BVHANDLE hsrc2, int compress);

/* perform logical OR operation on two bit vectors
   hdst = hdst OR hsrc
*/
BM_API_EXPORT int BM_bvector_combine_OR(BM_BVHANDLE hdst, BM_BVHANDLE hsrc);

/* perform 3-operand logical OR (set union) on two source bit vectors
   hdst = hsrc1 AND hsrc2
   compress - 0 (do not compress, 1 - compress)
*/
BM_API_EXPORT int BM_bvector_combine_OR_2sc(BM_BVHANDLE hdst, BM_BVHANDLE hsrc1, BM_BVHANDLE hsrc2, int compress);

/* perform logical SUB operation on two bit vectors
   hdst = hdst SUB hsrc
*/
BM_API_EXPORT int BM_bvector_combine_SUB(BM_BVHANDLE hdst, BM_BVHANDLE hsrc);

/* perform 3-operand logical AND NOT (set subtraction) between two source bit vectors
   hdst = hsrc1 SUB hsrc2
   compress - 0 (do not compress, 1 - compress)
*/
BM_API_EXPORT int BM_bvector_combine_SUB_2sc(BM_BVHANDLE hdst, BM_BVHANDLE hsrc1, BM_BVHANDLE hsrc2, int compress);

/* perform logical XOR operation on two bit vectors
   hdst = hdst XOR hsrc
*/
BM_API_EXPORT int BM_bvector_combine_XOR(BM_BVHANDLE hdst, BM_BVHANDLE hsrc);

/* perform 3-operand logical XOR (exclusive OR) between two source bit vectors
   hdst = hsrc1 XOR hsrc2
   compress - 0 (do not compress, 1 - compress)
*/
BM_API_EXPORT int BM_bvector_combine_XOR_2sc(BM_BVHANDLE hdst, BM_BVHANDLE hsrc1, BM_BVHANDLE hsrc2, int compress);


/* perform logical OR operation on two bit vectors
   hdst = hdst OR hsrc
   merge operation is faster than just OR, but it destroys the
   source vector (by borrowing its memory)
*/
BM_API_EXPORT int BM_bvector_merge(BM_BVHANDLE hdst, BM_BVHANDLE hsrc);


/* Logical shift right by 1
*/
BM_API_EXPORT int BM_bvector_rshift1(BM_BVHANDLE hdst);

/* -------------------------------------------- */
/* bvector operations with arrays               */
/* -------------------------------------------- */


/* perform logical OR operation over bit vector and an array
   hdst - destination bit vector handle
   arr_begin - array start
   arr_end   - array end (defined as start + size)
*/
BM_API_EXPORT
int BM_bvector_combine_OR_arr(BM_BVHANDLE hdst,
                               const unsigned int* arr_begin,
                               const unsigned int* arr_end);

/* perform logical XOR operation over bit vector and an array
   hdst - destination bit vector handle
   arr_begin - array start
   arr_end   - array end (defined as start + size)
*/
BM_API_EXPORT
int BM_bvector_combine_XOR_arr(BM_BVHANDLE hdst,
                               const unsigned int* arr_begin,
                               const unsigned int* arr_end);


/* perform logical SUB(MINUS) operation over bit vector and an array
   hdst - destination bit vector handle
   arr_begin - array start
   arr_end   - array end (defined as start + size)
*/
BM_API_EXPORT
int BM_bvector_combine_SUB_arr(BM_BVHANDLE hdst,
                               const unsigned int* arr_begin,
                               const unsigned int* arr_end);


/* perform logical AND operation over bit vector and an array
   hdst - destination bit vector handle
   arr_begin - array start
   arr_end   - array end (defined as start + size)
*/
BM_API_EXPORT
int BM_bvector_combine_AND_arr(BM_BVHANDLE hdst,
                               const unsigned int* arr_begin,
                               const unsigned int* arr_end);

/* perform logical AND operation over bit vector and a sorted array
   sort order has to ascending
   hdst - destination bit vector handle
   arr_begin - array start
   arr_end   - array end (defined as start + size)
*/
BM_API_EXPORT
int BM_bvector_combine_AND_arr_sorted(BM_BVHANDLE hdst,
                                      const unsigned int* arr_begin,
                                      const unsigned int* arr_end);


/* -------------------------------------------- */
/* bvector traversal/enumerator                 */
/* -------------------------------------------- */

/* construct bvector enumerator for ON bit index traversal
   (starting from first ON bit)
   h  - handle of source bvector
   peh - pointer on enumerator to be created
*/
BM_API_EXPORT int BM_bvector_enumerator_construct(BM_BVHANDLE h, BM_BVEHANDLE* peh);

/* construct bvector enumerator for ON bit index traversal
   starting from position
   h  - handle of source bvector
   peh - pointer on enumerator to be created
   pos - start position, if 0 - it starts from the first available bit
*/
BM_API_EXPORT
int BM_bvector_enumerator_construct_from(BM_BVHANDLE h,
                                         BM_BVEHANDLE* peh,
                                         unsigned int pos);


/* destroy bvector enumerator handle */
BM_API_EXPORT int BM_bvector_enumerator_free(BM_BVEHANDLE eh);

/* Check if enumerator is valid or reached the end of traversal
   pvalid - returns 0 if enumerator is no longer valid
   (empty vector or end of traversal)
*/
BM_API_EXPORT int BM_bvector_enumerator_is_valid(BM_BVEHANDLE eh, int* pvalid);

/* Return current enumerator value (traversal position)
   pvalue - returns bit traversal position
*/
BM_API_EXPORT int BM_bvector_enumerator_get_value(BM_BVEHANDLE eh, unsigned int* pvalue);

/* Advance enumerator to next traversal position.
   pvalid - (optional) returns 0 if traversal ended
   pvalue - (optional) current value
*/
BM_API_EXPORT
int BM_bvector_enumerator_next(BM_BVEHANDLE eh,
                               int* pvalid, unsigned int* pvalue);


/* Advance enumerator to the specieid traversal position.
   If bit-vector position is not 1 it will automatically go to the next available
   1 bit.
   pos    - requested position
   pvalid - (optional) returns 0 if traversal ended
   pvalue - (optional) current value
*/
BM_API_EXPORT
int BM_bvector_enumerator_goto(BM_BVEHANDLE eh, unsigned int pos,
                               int* pvalid, unsigned int* pvalue);


/* -------------------------------------------- */
/* bvector serialization                        */
/* -------------------------------------------- */

/*  serialize bit vector
    buf - buffer pointer 
      (should be allocated using BM_bvector_statistics.max_serialize_mem)
    buf_size - size of the buffer in bytes
    pblob_size - size of the serialized BLOB
*/
BM_API_EXPORT
int BM_bvector_serialize(BM_BVHANDLE h,
                         char*       buf,
                         size_t      buf_size,
                         size_t*     pblob_size);
    
/*  deserialize bit vector
    buf - buffer pointer 
      (should be allocated using BM_bvector_statistics.max_serialize_mem)
    buf_size - size of the buffer in bytes
*/
BM_API_EXPORT
int BM_bvector_deserialize(BM_BVHANDLE   h,
                           const char*   buf,
                           size_t        buf_size);
    
    
/* -------------------------------------------- */
/* bvector algorithms                           */
/* -------------------------------------------- */

/* compute population count of AND of two const bit vectors
   pcount - bit count of AND of two vectors
*/
BM_API_EXPORT int BM_bvector_count_AND(BM_BVHANDLE h1, BM_BVHANDLE h2, unsigned int* pcount);

/* return true if AND operation of two vectors produce any result
   (this is faster than count_AND)
   pany - non-zero if any bits were found
*/
BM_API_EXPORT int BM_bvector_any_AND(BM_BVHANDLE h1, BM_BVHANDLE h2, unsigned int* pany);

/* compute population count of XOR of two const bit vectors
   pcount - bit count of XOR of two vectors
*/
BM_API_EXPORT int BM_bvector_count_XOR(BM_BVHANDLE h1, BM_BVHANDLE h2, unsigned int* pcount);

/* return true if XOR operation of two vectors produce any result
   (this is faster than count_XOR)
   pany - non-zero if any bits were found
*/
BM_API_EXPORT int BM_bvector_any_XOR(BM_BVHANDLE h1, BM_BVHANDLE h2, unsigned int* pany);

/* compute population count of SUB of two const bit vectors
   pcount - bit count of SUB of two vectors
*/
BM_API_EXPORT int BM_bvector_count_SUB(BM_BVHANDLE h1, BM_BVHANDLE h2, unsigned int* pcount);

/* return true if SUB operation of two vectors produce any result
   (this is faster than count_SUB)
   pany - non-zero if any bits were found
*/
BM_API_EXPORT int BM_bvector_any_SUB(BM_BVHANDLE h1, BM_BVHANDLE h2, unsigned int* pany);

/* compute population count of OR of two const bit vectors
   pcount - bit count of OR of two vectors
*/
BM_API_EXPORT int BM_bvector_count_OR(BM_BVHANDLE h1, BM_BVHANDLE h2, unsigned int* pcount);

/* return true if SUB operation of two vectors produce any result
   (this is faster than count_SUB)
   pany - non-zero if any bits were found
*/
BM_API_EXPORT int BM_bvector_any_OR(BM_BVHANDLE h1, BM_BVHANDLE h2, unsigned int* pany);


#ifdef __cplusplus
}
#endif

#endif
