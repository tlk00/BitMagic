#ifndef LIBBM_INCLUDED_H__
#define LIBBM_INCLUDED_H__

/* Error codes */

#define BM_OK (0)
#define BM_ERR_BADALLOC (1)

#define BM_BVHANDLE void*

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/* Initialize libbm runtime before use*/
int BM_init(void*);

/* -------------------------------------------- */
/* bvector functions                            */
/* -------------------------------------------- */

/* construction and setters                     */

/* construct bvector handle 
   bv_max - maximum number of allowed bits (if 0 - allows maximum)
*/
int BM_bvector_construct(BM_BVHANDLE* h, unsigned int bv_max);

/* destroy bvector handle */
int BM_bvector_free(BM_BVHANDLE h);

/* set bit 
   i - index of a bit to set
*/
int BM_bvector_set_bit(BM_BVHANDLE h, unsigned int i);



/* -------------------------------------------- */
/* read-only getters                            */


/* get bit value */
int BM_bvector_get_bit(BM_BVHANDLE h, unsigned int i, unsigned int* ret);



#ifdef __cplusplus
}
#endif

#endif
