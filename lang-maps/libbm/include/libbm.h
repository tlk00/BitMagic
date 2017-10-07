#ifndef LIBBM_INCLUDED_H__
#define LIBBM_INCLUDED_H__

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

/* arguments codes and values */
#define BM_TRUE 1
#define BM_FALSE 0


#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

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

/* -------------------------------------------- */
/* bvector construction and sizing methods      */
/* -------------------------------------------- */

/* construction and setters                     */

/* construct bvector handle 
   bv_max - maximum number of allowed bits (if 0 - allows maximum)
   (it is recommened to set size to maximum always and do not use size params)
*/
int BM_bvector_construct(BM_BVHANDLE* h, unsigned int bv_max);

/* destroy bvector handle */
int BM_bvector_free(BM_BVHANDLE h);

/* get bit vector size
   psize - return size the bit vector
*/
int BM_bvector_get_size(BM_BVHANDLE h, unsigned int* psize);

/* get bit vector capacity
   pcap - return caapacity the bit vector
*/
int BM_bvector_get_capacity(BM_BVHANDLE h, unsigned int* pcap);


/* resize bit vector
   new_size - new requested size 
*/
int BM_bvector_set_size(BM_BVHANDLE h, unsigned int new_size);



/* -------------------------------------------- */
/* bvector functions to set bits                */
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


/* set all bits to 0 and (optionally) free unused memory
    free_mem - flag to release unused memory
*/
int BM_bvector_clear(BM_BVHANDLE h, int free_mem);


/* -------------------------------------------- */
/* bvector functions to read bits               */
/* -------------------------------------------- */


/* get bit value */
int BM_bvector_get_bit(BM_BVHANDLE h, unsigned int i, int* pval);


/* bitcount
   pcount - return number of ON bits in the vector
*/
int BM_bvector_count(BM_BVHANDLE h, unsigned int* pcount);



#ifdef __cplusplus
}
#endif

#endif
