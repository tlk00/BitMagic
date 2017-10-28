#ifndef LIBBM_INCLUDED_H__
#define LIBBM_INCLUDED_H__
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


#include <stddef.h>

/* Error codes */

#define BM_OK (0)
#define BM_ERR_BADALLOC (1)
#define BM_ERR_BADARG (2)
#define BM_ERR_RANGE (3)

/*
    error codes and messages
*/
#define BM_OK_MSG           "BM-00: All correct"
#define BM_ERR_BADALLOC_MSG "BM-01: Allocation error"
#define BM_ERR_BADARG_MSG   "BM-02: Invalid or missing function argument"
#define BM_ERR_RANGE_MSG    "BM-03: Incorrect range or index"

#define BM_UNK_MSG          "BM-XX: Unknown error"


/* bit-vector handle */
#define BM_BVHANDLE void*
/* bit-vector enumerator handle */
#define BM_BVEHANDLE void*


/* arguments codes and values */
#define BM_TRUE 1
#define BM_FALSE 0


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
int BM_init(void*);

/**
    return copyright info string and version information.
*/
const char* BM_version(int* major, int* minor, int* patch);

/**
    return error message by code
*/
const char* BM_error_msg(int errcode);

/* ------------------------------------------------- */
/* bvector construction, swap and sizing methods     */
/* ------------------------------------------------- */

/* construction and setters                          */

/* construct bvector handle 
   bv_max - maximum number of allowed bits (if 0 - allows maximum)
   (it is recommened to set size to maximum always and do not use size params)
*/
int BM_bvector_construct(BM_BVHANDLE* h, unsigned int bv_max);

/* construct bvector handle as a copy
   hfrom - another handle to copy from
*/
int BM_bvector_construct_copy(BM_BVHANDLE* h, BM_BVHANDLE hfrom);


/* destroy bvector handle */
int BM_bvector_free(BM_BVHANDLE h);

/* get bit vector size
   psize - return size the bit vector
*/
int BM_bvector_get_size(BM_BVHANDLE h, unsigned int* psize);

/* get bit vector capacity
   pcap - return caapacity the bit vector
*/
/*int BM_bvector_get_capacity(BM_BVHANDLE h, unsigned int* pcap);*/


/* resize bit vector
   new_size - new requested size
*/
int BM_bvector_set_size(BM_BVHANDLE h, unsigned int new_size);

/* swap two bit-vectors
*/
int BM_bvector_swap(BM_BVHANDLE h1, BM_BVHANDLE h2);


/* -------------------------------------------- */
/* bvector functions to set and clear bits      */
/* -------------------------------------------- */


/* set bit 
   i - index of a bit to set
   val - value (0 | 1)
*/
int BM_bvector_set_bit(BM_BVHANDLE h, unsigned int i, int val);


/* set bit only if current value equals the condition
   i - index of a bit to set
   val - value (0 | 1)
   condition - expected current value of bit i
   pchanged - optional return value, if bit was actually changed
*/
int BM_bvector_set_bit_conditional(BM_BVHANDLE  h,
                                   unsigned int i,
                                   int          val,
                                   int          condition,
                                   int*         pchanged);
    
/* flip bit
   i - index of a bit to flip
*/
int BM_bvector_flip_bit(BM_BVHANDLE h, unsigned int i);
    

/* set all bits to 1
*/
int BM_bvector_set(BM_BVHANDLE h);



/* Sets all bits in the specified closed interval [left,right]
   This method DOES NOT resize vector.

   left  - interval start
   right - interval end (closed interval)
   value - value to set interval in
 
*/
int BM_bvector_set_range(BM_BVHANDLE h,
                         unsigned int left,
                         unsigned int right,
                         int          value);
    
/* invert all bits in the bit vector
*/
int BM_bvector_invert(BM_BVHANDLE h);
    


/* set all bits to 0 and (optionally) free unused memory
    free_mem - flag to release unused memory
*/
int BM_bvector_clear(BM_BVHANDLE h, int free_mem);

/* find 1 bit index in the vector and set it to 0

   i - index of bit to search from
   pnext - return index of the next 1 bit. 0 - means no more 1 bits.
*/
int BM_bvector_extract_next(BM_BVHANDLE h, unsigned int i, unsigned int* pnext);



/* -------------------------------------------- */
/* bvector functions to read bits               */
/* -------------------------------------------- */


/* get bit value */
int BM_bvector_get_bit(BM_BVHANDLE h, unsigned int i, int* pval);


/* bitcount
   pcount - return number of ON bits in the vector
*/
int BM_bvector_count(BM_BVHANDLE h, unsigned int* pcount);

/* range bitcount
   left  - interval start
   right - interval end (closed interval)
   pcount - return number of ON bits in the vector
*/
int BM_bvector_count_range(BM_BVHANDLE   h,
                           unsigned int  left,
                           unsigned int  right,
                           unsigned int* pcount);

/* check if there are any bits set 
   pval - return non-zero value if any bits are ON
*/
int BM_bvector_any(BM_BVHANDLE h, int* pval);

/* Finds index of 1 bit starting from position
   from - initial search position
   ppos - found position of 1 bit (>= from)
   pfound - 0 if nothing found
*/
int BM_bvector_find(BM_BVHANDLE h,
                    unsigned int from, unsigned int* ppos, int* pfound);

/* find first 1 bit index in the vector

   pi - return index of first bit 
   found - return 0 if first bit not found (empty vector)
*/
int BM_bvector_get_first(BM_BVHANDLE h, unsigned int* pi, int* pfound);

/* find 1 bit index in the vector

   i - index of bit to search from
   pnext - return index of the next 1 bit. 0 - means no more 1 bits.
*/
int BM_bvector_get_next(BM_BVHANDLE h, unsigned int i, unsigned int* pnext);



/* -------------------------------------------- */
/* bvector operations                           */
/* -------------------------------------------- */

/* Lexicographical comparison of two bit vectors
   pres - returns -1 if h1 less than h2, 1 - greater, 0 - equal.
*/
int BM_bvector_compare(BM_BVHANDLE h1, BM_BVHANDLE h2, int* pres);



/* Perform bit vector memory optimization
   opt_mode - optimization level:
    (0 - default, 1 - free empty blocks, 2 - free empty and full blocks, 3 - GAP compress)
   pstat - optional post optimization statistics
*/
int BM_bvector_optimize(BM_BVHANDLE h,
                        int opt_mode,
                        struct BM_bvector_statistics* pstat);
    
/* Perform calculate bit vector statistics
   pstat - bit vector statistics
*/
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
int BM_bvector_combine_operation(BM_BVHANDLE hdst, BM_BVHANDLE hsrc, int opcode);
    
int BM_bvector_combine_AND(BM_BVHANDLE hdst, BM_BVHANDLE hsrc);
int BM_bvector_combine_OR(BM_BVHANDLE hdst, BM_BVHANDLE hsrc);
int BM_bvector_combine_SUB(BM_BVHANDLE hdst, BM_BVHANDLE hsrc);
int BM_bvector_combine_XOR(BM_BVHANDLE hdst, BM_BVHANDLE hsrc);

/* -------------------------------------------- */
/* bvector traversal/enumerator                 */
/* -------------------------------------------- */

/* construct bvector enumerator for bit index traversal
   h  - handle of source bvector
   peh - pointer on enumerator to be created
*/
int BM_bvector_enumerator_construct(BM_BVHANDLE h, BM_BVEHANDLE* peh);

/* destroy bvector enumerator handle */
int BM_bvector_enumerator_free(BM_BVEHANDLE eh);

/* Check if enumerator is valid or reached the end of traversal
   pvalid - returns 0 if enumerator is no longer valid
   (empty vector or end of traversal)
*/
int BM_bvector_enumerator_is_valid(BM_BVEHANDLE eh, int* pvalid);

/* Return current enumerator value (traversal position)
   pvalue - returns bit traversal position
*/
int BM_bvector_enumerator_get_value(BM_BVEHANDLE eh, unsigned int* pvalue);

/* Advance enumerator to next traversal position.
   pvalid - (optional) returns 0 if traversal ended
   pvalue - (optional) current value
*/
int BM_bvector_enumerator_next(BM_BVEHANDLE eh,
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
int BM_bvector_serialize(BM_BVHANDLE h,
                         char*       buf,
                         size_t      buf_size,
                         size_t*     pblob_size);
    
/*  deserialize bit vector
    buf - buffer pointer 
      (should be allocated using BM_bvector_statistics.max_serialize_mem)
    buf_size - size of the buffer in bytes
*/
int BM_bvector_deserialize(BM_BVHANDLE   h,
                           const char*   buf,
                           size_t        buf_size);
    



#ifdef __cplusplus
}
#endif

#endif
