#include "bmserial.h"
#include "bmalgo.h"


typedef bm::bvector<libbm::standard_allocator>::enumerator TBM_bvector_enumerator;


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
    if (h == 0)
        return BM_ERR_BADARG;
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

int BM_bvector_construct_copy(BM_BVHANDLE* h, BM_BVHANDLE hfrom)
{
    if (h == 0 || !hfrom)
        return BM_ERR_BADARG;
    TRY
    {
        void* mem = ::malloc(sizeof(TBM_bvector));
        if (mem == 0)
        {
            *h = 0;
            return BM_ERR_BADALLOC;
        }
        
        const TBM_bvector* bv_from = (TBM_bvector*)hfrom;
        
        // placement new just to call the copy constructor
        TBM_bvector* bv = new(mem) TBM_bvector(*bv_from);
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

int BM_bvector_find(BM_BVHANDLE h,
                    unsigned int from, unsigned int* ppos, int* pfound)
{
    if (!h || !ppos || !pfound)
        return BM_ERR_BADARG;
    
    TRY
    {
        const TBM_bvector* bv = (TBM_bvector*)h;
        *pfound = bv->find(from, *ppos);
    }
    CATCH (BM_ERR_BADALLOC)
    {
        return BM_ERR_BADALLOC;
    }
    ETRY;
    return BM_OK;
}


// -----------------------------------------------------------------

int BM_bvector_get_first(BM_BVHANDLE h, unsigned int* pi, int* pfound)
{
    if (!h || !pi || !pfound)
        return BM_ERR_BADARG;
    
    TRY
    {
        const TBM_bvector* bv = (TBM_bvector*)h;
        *pi = bv->get_first();
        if (*pi == 0)
        {
            if (!(bv->test(0)))
            {
                *pfound = 0;
                return BM_OK;
            }
        }
        *pfound = 1;
    }
    CATCH (BM_ERR_BADALLOC)
    {
        return BM_ERR_BADALLOC;
    }
    ETRY;
    return BM_OK;
}

// -----------------------------------------------------------------


int BM_bvector_get_next(BM_BVHANDLE h, unsigned int i, unsigned int* pnext)
{
    if (!h || !pnext)
        return BM_ERR_BADARG;
    
    TRY
    {
        const TBM_bvector* bv = (TBM_bvector*)h;
        *pnext = bv->get_next(i);
    }
    CATCH (BM_ERR_BADALLOC)
    {
        return BM_ERR_BADALLOC;
    }
    ETRY;
    return BM_OK;
}

// -----------------------------------------------------------------

int BM_bvector_extract_next(BM_BVHANDLE h, unsigned int i, unsigned int* pnext)
{
    if (!h || !pnext)
        return BM_ERR_BADARG;
    
    TRY
    {
        TBM_bvector* bv = (TBM_bvector*)h;
        *pnext = bv->extract_next(i);
    }
    CATCH (BM_ERR_BADALLOC)
    {
        return BM_ERR_BADALLOC;
    }
    ETRY;
    return BM_OK;
}

// -----------------------------------------------------------------


int BM_bvector_compare(BM_BVHANDLE h1, BM_BVHANDLE h2, int* pres)
{
    if (!h1 || !h2 || !pres)
        return BM_ERR_BADARG;
    
    TRY
    {
        const TBM_bvector* bv1 = (TBM_bvector*)h1;
        const TBM_bvector* bv2 = (TBM_bvector*)h2;
        
        *pres = bv1->compare(*bv2);
    }
    CATCH (BM_ERR_BADALLOC)
    {
        return BM_ERR_BADALLOC;
    }
    ETRY;
    return BM_OK;
    
}

// -----------------------------------------------------------------


int BM_bvector_optimize(BM_BVHANDLE            h,
                        int                    opt_mode,
                        struct BM_bvector_statistics* pstat)
{
    if (!h)
        return BM_ERR_BADARG;
    TBM_bvector::optmode omode = TBM_bvector::opt_compress;
    TBM_bvector::statistics stat;
    
    switch (opt_mode)
    {
    case 1: omode = TBM_bvector::opt_free_0; break;
    case 2: omode = TBM_bvector::opt_free_01; break;
    }
    
    TRY
    {
        BM_DECLARE_TEMP_BLOCK(tb)
    
        TBM_bvector* bv = (TBM_bvector*)h;
        bv->optimize(tb, omode, &stat);
        
        if (pstat)
        {
            pstat->bit_blocks = stat.bit_blocks;
            pstat->gap_blocks = stat.gap_blocks;
            pstat->max_serialize_mem = stat.max_serialize_mem;
            pstat->memory_used = stat.memory_used;
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


int BM_bvector_calc_stat(BM_BVHANDLE h,
                         struct BM_bvector_statistics* pstat)
{
    if (!h || !pstat)
        return BM_ERR_BADARG;
    TBM_bvector::statistics stat;
    
    TRY
    {
        const TBM_bvector* bv = (TBM_bvector*)h;
        bv->calc_stat(&stat);
        
        pstat->bit_blocks = stat.bit_blocks;
        pstat->gap_blocks = stat.gap_blocks;
        pstat->max_serialize_mem = stat.max_serialize_mem;
        pstat->memory_used = stat.memory_used;
    }
    CATCH (BM_ERR_BADALLOC)
    {
        return BM_ERR_BADALLOC;
    }
    ETRY;
    return BM_OK;
}

// -----------------------------------------------------------------


int BM_bvector_combine_operation(BM_BVHANDLE hdst, BM_BVHANDLE hsrc, int opcode)
{
    if (!hdst || !hsrc)
        return BM_ERR_BADARG;

    bm::operation opc;
    switch (opcode)
    {
    case 0: opc = bm::BM_AND; break;
    case 1: opc = bm::BM_OR;  break;
    case 2: opc = bm::BM_SUB; break;
    case 3: opc = bm::BM_XOR; break;
    default:
        return BM_ERR_BADARG;
    }
    
    TRY
    {
        TBM_bvector* bv1 = (TBM_bvector*)hdst;
        const TBM_bvector* bv2 = (TBM_bvector*)hsrc;
        
        bv1->combine_operation(*bv2, opc);
    }
    CATCH (BM_ERR_BADALLOC)
    {
        return BM_ERR_BADALLOC;
    }
    ETRY;
    return BM_OK;
}


int BM_bvector_combine_AND(BM_BVHANDLE hdst, BM_BVHANDLE hsrc)
{
    return BM_bvector_combine_operation(hdst, hsrc, 0);
}

int BM_bvector_combine_OR(BM_BVHANDLE hdst, BM_BVHANDLE hsrc)
{
    return BM_bvector_combine_operation(hdst, hsrc, 1);
}

int BM_bvector_combine_SUB(BM_BVHANDLE hdst, BM_BVHANDLE hsrc)
{
    return BM_bvector_combine_operation(hdst, hsrc, 2);
}

int BM_bvector_combine_XOR(BM_BVHANDLE hdst, BM_BVHANDLE hsrc)
{
    return BM_bvector_combine_operation(hdst, hsrc, 3);
}

// -----------------------------------------------------------------


int BM_bvector_combine_AND_arr(BM_BVHANDLE hdst,
                               const unsigned int* arr_begin,
                               const unsigned int* arr_end)
{
    if (!hdst)
        return BM_ERR_BADARG;
    
    TRY
    {
        TBM_bvector* bv = (TBM_bvector*)hdst;
        bm::combine_and(*bv, arr_begin, arr_end);
    }
    CATCH (BM_ERR_BADALLOC)
    {
        return BM_ERR_BADALLOC;
    }
    ETRY;
    return BM_OK;
}

// -----------------------------------------------------------------


int BM_bvector_combine_AND_arr_sorted(BM_BVHANDLE hdst,
                                      const unsigned int* arr_begin,
                                      const unsigned int* arr_end)
{
    if (!hdst)
        return BM_ERR_BADARG;
    
    TRY
    {
        TBM_bvector* bv = (TBM_bvector*)hdst;
        bm::combine_and_sorted(*bv, arr_begin, arr_end);
    }
    CATCH (BM_ERR_BADALLOC)
    {
        return BM_ERR_BADALLOC;
    }
    ETRY;
    return BM_OK;
}

// -----------------------------------------------------------------


int BM_bvector_combine_OR_arr(BM_BVHANDLE hdst,
                               const unsigned int* arr_begin,
                               const unsigned int* arr_end)
{
    if (!hdst)
        return BM_ERR_BADARG;
    
    TRY
    {
        TBM_bvector* bv = (TBM_bvector*)hdst;
        bm::combine_or(*bv, arr_begin, arr_end);
    }
    CATCH (BM_ERR_BADALLOC)
    {
        return BM_ERR_BADALLOC;
    }
    ETRY;
    return BM_OK;
}

// -----------------------------------------------------------------

int BM_bvector_combine_XOR_arr(BM_BVHANDLE hdst,
                               const unsigned int* arr_begin,
                               const unsigned int* arr_end)
{
    if (!hdst)
        return BM_ERR_BADARG;
    
    TRY
    {
        TBM_bvector* bv = (TBM_bvector*)hdst;
        bm::combine_xor(*bv, arr_begin, arr_end);
    }
    CATCH (BM_ERR_BADALLOC)
    {
        return BM_ERR_BADALLOC;
    }
    ETRY;
    return BM_OK;
}


// -----------------------------------------------------------------

int BM_bvector_combine_SUB_arr(BM_BVHANDLE hdst,
                               const unsigned int* arr_begin,
                               const unsigned int* arr_end)
{
    if (!hdst)
        return BM_ERR_BADARG;
    
    TRY
    {
        TBM_bvector* bv = (TBM_bvector*)hdst;
        bm::combine_sub(*bv, arr_begin, arr_end);
    }
    CATCH (BM_ERR_BADALLOC)
    {
        return BM_ERR_BADALLOC;
    }
    ETRY;
    return BM_OK;
}


// -----------------------------------------------------------------


int BM_bvector_serialize(BM_BVHANDLE h,
                         char*       buf,
                         size_t      buf_size,
                         size_t*     pblob_size)
{
    if (!h || !pblob_size)
        return BM_ERR_BADARG;
    
    TRY
    {
        BM_DECLARE_TEMP_BLOCK(tb)
    
        const TBM_bvector* bv = (TBM_bvector*)h;
        
        bm::serializer<TBM_bvector> bvs(TBM_bvector::allocator_type(), tb);
        bvs.set_compression_level(4);
        
        *pblob_size = bvs.serialize(*bv, (unsigned char*)buf, (unsigned)buf_size);
    }
    CATCH (BM_ERR_BADALLOC)
    {
        return BM_ERR_BADALLOC;
    }
    ETRY;
    return BM_OK;
}


// -----------------------------------------------------------------

int BM_bvector_deserialize(BM_BVHANDLE   h,
                           const char*   buf,
                           size_t        /*buf_size*/)
{
    if (!h)
        return BM_ERR_BADARG;
    
    TRY
    {
        TBM_bvector* bv = (TBM_bvector*)h;
        bm::deserialize(*bv, (const unsigned char*)buf);
    }
    CATCH (BM_ERR_BADALLOC)
    {
        return BM_ERR_BADALLOC;
    }
    ETRY;
    return BM_OK;
}

// -----------------------------------------------------------------

int BM_bvector_enumerator_construct(BM_BVHANDLE h, BM_BVEHANDLE* peh)
{
    return BM_bvector_enumerator_construct_from(h, peh, 0);
}

// -----------------------------------------------------------------


int BM_bvector_enumerator_construct_from(BM_BVHANDLE h,
                                         BM_BVEHANDLE* peh,
                                         unsigned int  pos)
{
    if (h == 0 || peh == 0)
        return BM_ERR_BADARG;

    TRY
    {
        TBM_bvector* bv = (TBM_bvector*)h;

        void* mem = ::malloc(sizeof(TBM_bvector_enumerator));
        if (mem == 0)
        {
            *peh = 0;
            return BM_ERR_BADALLOC;
        }
        // placement new just to call the constructor
        TBM_bvector_enumerator* bvenum = new(mem) TBM_bvector_enumerator(bv, pos);
        *peh = bvenum;
        
    }
    CATCH (BM_ERR_BADALLOC)
    {
        *peh = 0;
        return BM_ERR_BADALLOC;
    }
    ETRY;

    return BM_OK;
}

// -----------------------------------------------------------------

int BM_bvector_enumerator_free(BM_BVEHANDLE eh)
{
    if (!eh)
        return BM_ERR_BADARG;
    TBM_bvector_enumerator* bvenum = (TBM_bvector_enumerator*)eh;
    bvenum->~TBM_bvector_enumerator();
    ::free(eh);

    return BM_OK;
}

// -----------------------------------------------------------------

int BM_bvector_enumerator_is_valid(BM_BVEHANDLE eh, int* valid)
{
    if (!eh || !valid)
        return BM_ERR_BADARG;
    
    TRY
    {
        TBM_bvector_enumerator* bvenum = (TBM_bvector_enumerator*)eh;
        *valid = bvenum->valid();
    }
    CATCH (BM_ERR_BADALLOC)
    {
        return BM_ERR_BADALLOC;
    }
    ETRY;
    return BM_OK;
}

// -----------------------------------------------------------------

int BM_bvector_enumerator_get_value(BM_BVEHANDLE eh, unsigned int* pvalue)
{
    if (!eh || !pvalue)
        return BM_ERR_BADARG;
    
    TRY
    {
        TBM_bvector_enumerator* bvenum = (TBM_bvector_enumerator*)eh;
        *pvalue = bvenum->value();
    }
    CATCH (BM_ERR_BADALLOC)
    {
        return BM_ERR_BADALLOC;
    }
    ETRY;
    return BM_OK;
}

// -----------------------------------------------------------------

int BM_bvector_enumerator_next(BM_BVEHANDLE eh,
                               int* pvalid, unsigned int* pvalue)
{
    if (!eh)
        return BM_ERR_BADARG;
    
    TRY
    {
        TBM_bvector_enumerator* bvenum = (TBM_bvector_enumerator*)eh;
        bvenum->go_up();
        if (pvalid)
        {
            *pvalid = bvenum->valid();
        }
        if (pvalue)
        {
            *pvalue = bvenum->value();
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

int BM_bvector_enumerator_goto(BM_BVEHANDLE eh, unsigned int pos,
                               int* pvalid, unsigned int* pvalue)
{
    if (!eh)
        return BM_ERR_BADARG;
    
    TRY
    {
        TBM_bvector_enumerator* bvenum = (TBM_bvector_enumerator*)eh;
        bvenum->go_to(pos);
        if (pvalid)
        {
            *pvalid = bvenum->valid();
        }
        if (pvalue)
        {
            *pvalue = bvenum->value();
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

int BM_bvector_count_AND(BM_BVHANDLE h1, BM_BVHANDLE h2, unsigned int* pcount)
{
    if (!h1 || !h2 || !pcount)
        return BM_ERR_BADARG;
    TRY
    {
        const TBM_bvector* bv1 = (TBM_bvector*)h1;
        const TBM_bvector* bv2 = (TBM_bvector*)h2;
        *pcount = bm::count_and(*bv1, *bv2);
    }
    CATCH (BM_ERR_BADALLOC)
    {
        return BM_ERR_BADALLOC;
    }
    ETRY;
    return BM_OK;
}

// -----------------------------------------------------------------

int BM_bvector_any_AND(BM_BVHANDLE h1, BM_BVHANDLE h2, unsigned int* pany)
{
    if (!h1 || !h2 || !pany)
        return BM_ERR_BADARG;
    TRY
    {
        const TBM_bvector* bv1 = (TBM_bvector*)h1;
        const TBM_bvector* bv2 = (TBM_bvector*)h2;
        *pany = bm::any_and(*bv1, *bv2);
    }
    CATCH (BM_ERR_BADALLOC)
    {
        return BM_ERR_BADALLOC;
    }
    ETRY;
    return BM_OK;
}

// -----------------------------------------------------------------

int BM_bvector_count_XOR(BM_BVHANDLE h1, BM_BVHANDLE h2, unsigned int* pcount)
{
    if (!h1 || !h2 || !pcount)
        return BM_ERR_BADARG;
    TRY
    {
        const TBM_bvector* bv1 = (TBM_bvector*)h1;
        const TBM_bvector* bv2 = (TBM_bvector*)h2;
        *pcount = bm::count_xor(*bv1, *bv2);
    }
    CATCH (BM_ERR_BADALLOC)
    {
        return BM_ERR_BADALLOC;
    }
    ETRY;
    return BM_OK;
}

// -----------------------------------------------------------------

int BM_bvector_any_XOR(BM_BVHANDLE h1, BM_BVHANDLE h2, unsigned int* pany)
{
    if (!h1 || !h2 || !pany)
        return BM_ERR_BADARG;
    TRY
    {
        const TBM_bvector* bv1 = (TBM_bvector*)h1;
        const TBM_bvector* bv2 = (TBM_bvector*)h2;
        *pany = bm::any_xor(*bv1, *bv2);
    }
    CATCH (BM_ERR_BADALLOC)
    {
        return BM_ERR_BADALLOC;
    }
    ETRY;
    return BM_OK;
}

// -----------------------------------------------------------------

int BM_bvector_count_SUB(BM_BVHANDLE h1, BM_BVHANDLE h2, unsigned int* pcount)
{
    if (!h1 || !h2 || !pcount)
        return BM_ERR_BADARG;
    TRY
    {
        const TBM_bvector* bv1 = (TBM_bvector*)h1;
        const TBM_bvector* bv2 = (TBM_bvector*)h2;
        *pcount = bm::count_sub(*bv1, *bv2);
    }
    CATCH (BM_ERR_BADALLOC)
    {
        return BM_ERR_BADALLOC;
    }
    ETRY;
    return BM_OK;
}

// -----------------------------------------------------------------

int BM_bvector_any_SUB(BM_BVHANDLE h1, BM_BVHANDLE h2, unsigned int* pany)
{
    if (!h1 || !h2 || !pany)
        return BM_ERR_BADARG;
    TRY
    {
        const TBM_bvector* bv1 = (TBM_bvector*)h1;
        const TBM_bvector* bv2 = (TBM_bvector*)h2;
        *pany = bm::any_sub(*bv1, *bv2);
    }
    CATCH (BM_ERR_BADALLOC)
    {
        return BM_ERR_BADALLOC;
    }
    ETRY;
    return BM_OK;
}

// -----------------------------------------------------------------

int BM_bvector_count_OR(BM_BVHANDLE h1, BM_BVHANDLE h2, unsigned int* pcount)
{
    if (!h1 || !h2 || !pcount)
        return BM_ERR_BADARG;
    TRY
    {
        const TBM_bvector* bv1 = (TBM_bvector*)h1;
        const TBM_bvector* bv2 = (TBM_bvector*)h2;
        *pcount = bm::count_or(*bv1, *bv2);
    }
    CATCH (BM_ERR_BADALLOC)
    {
        return BM_ERR_BADALLOC;
    }
    ETRY;
    return BM_OK;
}

// -----------------------------------------------------------------

int BM_bvector_any_OR(BM_BVHANDLE h1, BM_BVHANDLE h2, unsigned int* pany)
{
    if (!h1 || !h2 || !pany)
        return BM_ERR_BADARG;
    TRY
    {
        const TBM_bvector* bv1 = (TBM_bvector*)h1;
        const TBM_bvector* bv2 = (TBM_bvector*)h2;
        *pany = bm::any_or(*bv1, *bv2);
    }
    CATCH (BM_ERR_BADALLOC)
    {
        return BM_ERR_BADALLOC;
    }
    ETRY;
    return BM_OK;
}

// -----------------------------------------------------------------

