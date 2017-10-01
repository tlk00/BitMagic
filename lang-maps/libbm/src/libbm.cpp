#include "libbm.h"
#include "try_throw_catch.h"
static jmp_buf ex_buf__;

#define BM_NO_STL
#define BMALLOC__H__INCLUDED__

#include "bmconst.h"
#include "bmcalloc.h"
#include "bm.h"

#include <new>

typedef libbm::standard_allocator              TBM_Alloc;
typedef bm::bvector<libbm::standard_allocator> TBM_bvector;


int BM_init(void*)
{
	return BM_OK;
}


int BM_bvector_construct(BM_BVHANDLE* h, unsigned int bv_max)
{
	TRY
	{
		void* mem = ::malloc(sizeof(TBM_bvector));
		if (mem == 0)
		{
			*h = 0;
			return BM_ERR_BADALLOC;
		}
        if (bv_max == 0)
        {
            bv_max = bm::id_max;
        }
		// placement new just to call the constructor
		TBM_bvector* bv = new(mem) TBM_bvector(bm::BM_GAP,
                                               bm::gap_len_table<true>::_len,
                                               bv_max,
                                               TBM_Alloc());
		*h = bv;
	}
	CATCH (BM_ERR_BADALLOC)
	{
		*h = 0;
		return BM_ERR_BADALLOC;
	}
	ETRY;

	return BM_OK;
}


int BM_bvector_free(BM_BVHANDLE h)
{
	if (!h)
		return BM_ERR_BADARG;
	TBM_bvector* bv = (TBM_bvector*)h;
	bv->~TBM_bvector();
	::free(h);

	return BM_OK;
}

int BM_bvector_set_bit(BM_BVHANDLE h, unsigned int i, unsigned int val)
{
	if (!h)
		return BM_ERR_BADARG;
	TRY
	{
        TBM_bvector* bv = (TBM_bvector*)h;
        bv->set(i, val);
	}
	CATCH (BM_ERR_BADALLOC)
	{
		return BM_ERR_BADALLOC;
	}
	ETRY;
	return BM_OK;
}

int BM_bvector_get_bit(BM_BVHANDLE h, unsigned int i, unsigned int* pval)
{
	if (!h || !pval)
		return BM_ERR_BADARG;
	TRY
	{
        const TBM_bvector* bv = (TBM_bvector*)h;
        *pval = bv->test(i);
	}
	CATCH (BM_ERR_BADALLOC)
	{
		return BM_ERR_BADALLOC;
	}
	ETRY;
	return BM_OK;
}


int BM_bvector_count(BM_BVHANDLE h, unsigned int* pcount)
{
	if (!h || !pcount)
		return BM_ERR_BADARG;
	TRY
	{
        const TBM_bvector* bv = (TBM_bvector*)h;
        *pcount = bv->count();
	}
	CATCH (BM_ERR_BADALLOC)
	{
		return BM_ERR_BADALLOC;
	}
	ETRY;
	return BM_OK;
}


