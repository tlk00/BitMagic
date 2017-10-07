
// -----------------------------------------------------------------

int BM_init(void*)
{
	return BM_OK;
}

// -----------------------------------------------------------------

const char* BM_version(int* major, int* minor, int* patch)
{
    if (major)
        *major = bm::_copyright<true>::_v[0];
    if (minor)
        *minor = bm::_copyright<true>::_v[1];
    if (patch)
        *patch = bm::_copyright<true>::_v[2];
    
    return bm::_copyright<true>::_p;
}

// -----------------------------------------------------------------

const char* BM_error_msg(int errcode)
{
    switch (errcode)
    {
    case BM_OK:
        return BM_OK_MSG;
    case BM_ERR_BADALLOC:
        return BM_ERR_BADALLOC_MSG;
    case BM_ERR_BADARG:
        return BM_ERR_BADARG_MSG;
    case BM_ERR_RANGE:
        return BM_ERR_RANGE_MSG;
    }
    return BM_UNK_MSG;
}

// -----------------------------------------------------------------


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
		TBM_bvector* bv = new(mem) TBM_bvector(bm::BM_BIT,
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

// -----------------------------------------------------------------

int BM_bvector_free(BM_BVHANDLE h)
{
	if (!h)
		return BM_ERR_BADARG;
	TBM_bvector* bv = (TBM_bvector*)h;
	bv->~TBM_bvector();
	::free(h);

	return BM_OK;
}

// -----------------------------------------------------------------

int BM_bvector_swap(BM_BVHANDLE h1, BM_BVHANDLE h2)
{
	if (!h1 || !h2)
		return BM_ERR_BADARG;
	TRY
	{
        TBM_bvector* bv1 = (TBM_bvector*)h1;
        TBM_bvector* bv2 = (TBM_bvector*)h2;
        bv1->swap(*bv2);
	}
	CATCH (BM_ERR_BADALLOC)
	{
		return BM_ERR_BADALLOC;
	}
	ETRY;
	return BM_OK;

}

// -----------------------------------------------------------------


int BM_bvector_get_size(BM_BVHANDLE h, unsigned int* psize)
{
	if (!h || !psize)
		return BM_ERR_BADARG;
	TRY
	{
        const TBM_bvector* bv = (TBM_bvector*)h;
        *psize = bv->size();
	}
	CATCH (BM_ERR_BADALLOC)
	{
		return BM_ERR_BADALLOC;
	}
	ETRY;
	return BM_OK;
}

// -----------------------------------------------------------------

int BM_bvector_get_capacity(BM_BVHANDLE h, unsigned int* pcap)
{
	if (!h)
		return BM_ERR_BADARG;
	TRY
	{
        const TBM_bvector* bv = (TBM_bvector*)h;
        if (pcap)
        {
            *pcap = bv->capacity();
        }
	}
	CATCH (BM_ERR_BADALLOC)
	{
		return BM_ERR_BADALLOC;
	}
	ETRY;
	return BM_OK;
}

// -----------------------------------------------------------------


int BM_bvector_set_size(BM_BVHANDLE h, unsigned int new_size)
{
	if (!h)
		return BM_ERR_BADARG;
	TRY
	{
        TBM_bvector* bv = (TBM_bvector*)h;
        bv->resize(new_size);
	}
	CATCH (BM_ERR_BADALLOC)
	{
		return BM_ERR_BADALLOC;
	}
	ETRY;
	return BM_OK;
}


// -----------------------------------------------------------------

int BM_bvector_set_bit(BM_BVHANDLE h, unsigned int i, int val)
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

// -----------------------------------------------------------------

int BM_bvector_flip_bit(BM_BVHANDLE h, unsigned int i)
{
	if (!h)
		return BM_ERR_BADARG;
	TRY
	{
        TBM_bvector* bv = (TBM_bvector*)h;
        bv->flip(i);
	}
	CATCH (BM_ERR_BADALLOC)
	{
		return BM_ERR_BADALLOC;
	}
	ETRY;
	return BM_OK;
}



// -----------------------------------------------------------------

int BM_bvector_set_bit_conditional(BM_BVHANDLE  h,
                                   unsigned int i,
                                   int          val,
                                   int          condition,
                                   int*         pchanged)
{
    unsigned int sz;
	if (!h)
		return BM_ERR_BADARG;
    
	TRY
	{
        TBM_bvector* bv = (TBM_bvector*)h;
        sz = bv->size();
        if (i >= sz)
            return BM_ERR_RANGE;
        
        bool b_changed = bv->set_bit_conditional(i, val, condition);
        if (pchanged)
        {
            *pchanged = b_changed;
        }
	}
	CATCH (BM_ERR_BADALLOC)
	{
		return BM_ERR_BADALLOC;
	}
	ETRY;
	return BM_OK;
}

// -----------------------------------------------------------------

int BM_bvector_set(BM_BVHANDLE h)
{
	if (!h)
		return BM_ERR_BADARG;
	TRY
	{
        TBM_bvector* bv = (TBM_bvector*)h;
        bv->set();
	}
	CATCH (BM_ERR_BADALLOC)
	{
		return BM_ERR_BADALLOC;
	}
	ETRY;
	return BM_OK;
}

// -----------------------------------------------------------------

int BM_bvector_set_range(BM_BVHANDLE h,
                         unsigned int left,
                         unsigned int right,
                         int          value)
{
	if (!h)
		return BM_ERR_BADARG;
    if (left > right)
		return BM_ERR_BADARG;
    
	TRY
	{
        TBM_bvector* bv = (TBM_bvector*)h;
        bv->set_range(left, right, value);
	}
	CATCH (BM_ERR_BADALLOC)
	{
		return BM_ERR_BADALLOC;
	}
	ETRY;
	return BM_OK;
}


// -----------------------------------------------------------------

int BM_bvector_clear(BM_BVHANDLE h, int free_mem)
{
	if (!h)
		return BM_ERR_BADARG;
	TRY
	{
        TBM_bvector* bv = (TBM_bvector*)h;
        bv->clear(free_mem);
	}
	CATCH (BM_ERR_BADALLOC)
	{
		return BM_ERR_BADALLOC;
	}
	ETRY;
	return BM_OK;
}

// -----------------------------------------------------------------

int BM_bvector_invert(BM_BVHANDLE h)
{
	if (!h)
		return BM_ERR_BADARG;
	TRY
	{
        TBM_bvector* bv = (TBM_bvector*)h;
        bv->invert();
	}
	CATCH (BM_ERR_BADALLOC)
	{
		return BM_ERR_BADALLOC;
	}
	ETRY;
	return BM_OK;
}

// -----------------------------------------------------------------

int BM_bvector_any(BM_BVHANDLE h, int* pval)
{
	if (!h || !pval)
		return BM_ERR_BADARG;
	TRY
	{
        const TBM_bvector* bv = (TBM_bvector*)h;
        *pval = bv->any();
	}
	CATCH (BM_ERR_BADALLOC)
	{
		return BM_ERR_BADALLOC;
	}
	ETRY;
	return BM_OK;
}



// -----------------------------------------------------------------

int BM_bvector_get_bit(BM_BVHANDLE h, unsigned int i,  int* pval)
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

// -----------------------------------------------------------------

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

// -----------------------------------------------------------------


int BM_bvector_count_range(BM_BVHANDLE h,
                           unsigned int  left,
                           unsigned int  right,
                           unsigned int* pcount)
{
	if (!h || !pcount)
		return BM_ERR_BADARG;
    if (left > right)
		return BM_ERR_BADARG;
    
	TRY
	{
        const TBM_bvector* bv = (TBM_bvector*)h;
        *pcount = bv->count_range(left, right);
	}
	CATCH (BM_ERR_BADALLOC)
	{
		return BM_ERR_BADALLOC;
	}
	ETRY;
	return BM_OK;
}

// -----------------------------------------------------------------


