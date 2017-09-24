#ifndef LIBBM_INCLUDED_H__
#define LIBBM_INCLUDED_H__

/* Error codes */

#define BM_OK (0)
#define BM_ERR_BADALLOC (1)

#define BM_BVHANDLE void*

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */


int BM_init(void*);

int BM_construct(BM_BVHANDLE* h);
int BM_destruct(BM_BVHANDLE h);
int BM_set_bit(BM_BVHANDLE h, unsigned int i);
int BM_get_bit(BM_BVHANDLE h, unsigned int i, unsigned int* ret);



#ifdef __cplusplus
}
#endif

#endif
