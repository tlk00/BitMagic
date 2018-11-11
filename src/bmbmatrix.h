#ifndef BMBMATRIX__H__INCLUDED__
#define BMBMATRIX__H__INCLUDED__
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

/*! \file bmbmatirx.h
    \brief basic bit-matrix classes
*/

#include <stddef.h>
#include "bmconst.h"

#ifndef BM_NO_STL
#include <stdexcept>
#endif

#include "bm.h"
#include "bmtrans.h"
#include "bmdef.h"



namespace bm
{

/**
    Basic dense bit-matrix class.
 
    Container of row-major bit-vectors, forming a bit-matrix.
    This class uses dense form of row storage.
    It is applicable as a build block for other sparse containers and
    succinct data structures, implementing high level abstractions.

    @ingroup bmagic
    @internal
*/
template<typename BV>
class basic_bmatrix
{
public:
    typedef BV                                       bvector_type;
    typedef bvector_type*                            bvector_type_ptr;
    typedef const bvector_type*                      bvector_type_const_ptr;
    typedef typename BV::allocator_type              allocator_type;
    typedef typename bvector_type::allocation_policy allocation_policy_type;
    typedef typename allocator_type::allocator_pool_type allocator_pool_type;
    typedef bm::id_t                                 size_type;

public:
    // ------------------------------------------------------------
    /*! @name Construction, assignment                          */
    ///@{

    basic_bmatrix(size_type rsize,
                  allocation_policy_type ap = allocation_policy_type(),
                  size_type bv_max_size = bm::id_max,
                  const allocator_type&   alloc  = allocator_type());
    ~basic_bmatrix() BMNOEXEPT;
    
    /*! copy-ctor */
    basic_bmatrix(const basic_bmatrix<BV>& bbm);
    basic_bmatrix<BV>& operator = (const basic_bmatrix<BV>& bbm)
    {
        copy_from(bbm);
        return *this;
    }
    
#ifndef BM_NO_CXX11
    /*! move-ctor */
    basic_bmatrix(basic_bmatrix<BV>&& bbm) BMNOEXEPT;

    /*! move assignmment operator */
    basic_bmatrix<BV>& operator = (basic_bmatrix<BV>&& bbm) BMNOEXEPT
    {
        if (this != &bbm)
        {
            free_rows();
            swap(bbm);
        }
        return *this;
    }
#endif

    void set_allocator_pool(allocator_pool_type* pool_ptr) { pool_ = pool_ptr; }

    ///@}
    
    // ------------------------------------------------------------
    /*! @name content manipulatio                                */
    ///@{

    /*! Swap content */
    void swap(basic_bmatrix<BV>& bbm) BMNOEXEPT;
    
    /*! Copy content */
    void copy_from(const basic_bmatrix<BV>& bbm);
    
    ///@}

    // ------------------------------------------------------------
    /*! @name row access                                         */
    ///@{

    /*! Get row bit-vector */
    const bvector_type* row(size_type i) const;

    /*! Get row bit-vector */
    bvector_type_const_ptr get_row(size_type i) const;

    /*! Get row bit-vector */
    bvector_type* get_row(size_type i);
    
    size_type rows() const { return rsize_; }
    
    /*! Make sure row is constructed, return bit-vector */
    bvector_type_ptr construct_row(size_type row);

    /*! Make sure row is copy-constructed, return bit-vector */
    bvector_type_ptr construct_row(size_type row, const bvector_type& bv);

    void destruct_row(size_type row);
    ///@}

public:
    // ------------------------------------------------------------
    /*! @name Utility function                                   */
    ///@{
    
    /// Test if 4 rows from i are not NULL
    bool test_4rows(unsigned i) const;
    ///@}


protected:
    void allocate_rows(size_type rsize);
    void free_rows() BMNOEXEPT;

    bvector_type* construct_bvector(const bvector_type* bv) const;
    void destruct_bvector(bvector_type* bv) const;

    static
    void throw_bad_alloc() { BV::throw_bad_alloc(); }

protected:
    size_type                bv_size_;
    allocator_type           alloc_;
    allocation_policy_type   ap_;
    allocator_pool_type*     pool_;
    
    bvector_type_ptr*        bv_rows_;
    size_type                rsize_;
};


//---------------------------------------------------------------------
//
//---------------------------------------------------------------------

template<typename BV>
basic_bmatrix<BV>::basic_bmatrix(size_type rsize,
              allocation_policy_type ap,
              size_type bv_max_size,
              const allocator_type&   alloc)
: bv_size_(bv_max_size),
  alloc_(alloc),
  ap_(ap),
  pool_(0),
  bv_rows_(0),
  rsize_(0)
{
    allocate_rows(rsize);
}

//---------------------------------------------------------------------

template<typename BV>
basic_bmatrix<BV>::~basic_bmatrix() BMNOEXEPT
{
    free_rows();
}

//---------------------------------------------------------------------

template<typename BV>
basic_bmatrix<BV>::basic_bmatrix(const basic_bmatrix<BV>& bbm)
: bv_size_(bbm.bv_size_),
  alloc_(bbm.alloc_),
  ap_(bbm.ap_),
  pool_(0),
  bv_rows_(0),
  rsize_(0)
{
    copy_from(bbm);
}

//---------------------------------------------------------------------

template<typename BV>
basic_bmatrix<BV>::basic_bmatrix(basic_bmatrix<BV>&& bbm) BMNOEXEPT
: bv_size_(bbm.bv_size_),
  alloc_(bbm.alloc_),
  ap_(bbm.ap_),
  pool_(0),
  bv_rows_(0),
  rsize_(0)
{
    swap(bbm);
}

//---------------------------------------------------------------------

template<typename BV>
const typename basic_bmatrix<BV>::bvector_type*
basic_bmatrix<BV>::row(size_type i) const
{
    BM_ASSERT(i < rsize_);
    return bv_rows_[i];
}

//---------------------------------------------------------------------

template<typename BV>
const typename basic_bmatrix<BV>::bvector_type*
basic_bmatrix<BV>::get_row(size_type i) const
{
    BM_ASSERT(i < rsize_);
    return bv_rows_[i];
}

//---------------------------------------------------------------------

template<typename BV>
typename basic_bmatrix<BV>::bvector_type*
basic_bmatrix<BV>::get_row(size_type i)
{
    BM_ASSERT(i < rsize_);
    return bv_rows_[i];
}

//---------------------------------------------------------------------

template<typename BV>
bool basic_bmatrix<BV>::test_4rows(unsigned j) const
{
    BM_ASSERT((j + 4) <= rsize_);
#if defined(BM64_SSE4)
        __m128i w0 = _mm_loadu_si128((__m128i*)(bv_rows_ + j));
        __m128i w1 = _mm_loadu_si128((__m128i*)(bv_rows_ + j + 2));
        w0 = _mm_or_si128(w0, w1);
        return !_mm_testz_si128(w0, w0);
#elif defined(BM64_AVX2) || defined(BM64_AVX512)
        __m256i w0 = _mm256_loadu_si256((__m256i*)(bv_rows_ + j));
        return !_mm256_testz_si256(w0, w0);
#else
        bool b = bv_rows_[j + 0] || bv_rows_[j + 1] || bv_rows_[j + 2] || bv_rows_[j + 3];
        return b;
#endif
}

//---------------------------------------------------------------------

template<typename BV>
void basic_bmatrix<BV>::copy_from(const basic_bmatrix<BV>& bbm)
{
    if (this == &bbm) // nothing to do
        return;
    free_rows();

    bv_size_ = bbm.bv_size_;
    alloc_ = bbm.alloc_;
    ap_ = bbm.ap_;

    size_type rsize = bbm.rsize_;
    if (rsize)
    {
        bv_rows_ = (bvector_type_ptr*)alloc_.alloc_ptr(rsize);
        if (!bv_rows_)
            throw_bad_alloc();
        else
        {
            rsize_ = rsize;
            for (size_type i = 0; i < rsize_; ++i)
            {
                const bvector_type_ptr bv = bbm.bv_rows_[i];
                bv_rows_[i] = bv ?  construct_bvector(bv) : 0;
            }
        }
    }

}


//---------------------------------------------------------------------

template<typename BV>
void basic_bmatrix<BV>::allocate_rows(size_type rsize)
{
    BM_ASSERT(!bv_rows_);
    
    if (rsize)
    {
        bv_rows_ = (bvector_type_ptr*)alloc_.alloc_ptr(rsize);
        if (!bv_rows_)
            throw_bad_alloc();
        else
        {
            rsize_ = rsize;
            for (size_type i = 0; i < rsize; ++i)
                bv_rows_[i] = 0;
        }
    }
}

//---------------------------------------------------------------------

template<typename BV>
void basic_bmatrix<BV>::free_rows() BMNOEXEPT
{
    for (size_type i = 0; i < rsize_; ++i)
    {
        bvector_type_ptr bv = bv_rows_[i];
        if (bv)
        {
            destruct_bvector(bv);
            bv_rows_[i] = 0;
        }
    } // for i
    if (bv_rows_)
        alloc_.free_ptr(bv_rows_, rsize_);
    bv_rows_ = 0;
}

//---------------------------------------------------------------------

template<typename BV>
void basic_bmatrix<BV>::swap(basic_bmatrix<BV>& bbm) BMNOEXEPT
{
    if (this == &bbm)
        return;
    
    bm::xor_swap(bv_size_, bbm.bv_size_);

    allocator_type alloc_tmp = alloc_;
    alloc_ = bbm.alloc_;
    bbm.alloc_ = alloc_tmp;

    allocation_policy_type ap_tmp = ap_;
    ap_ = bbm.ap_;
    bbm.ap_ = ap_tmp;
    
    allocator_pool_type*     pool_tmp = pool_;
    pool_ = bbm.pool_;
    bbm.pool_ = pool_tmp;

    bm::xor_swap(rsize_, bbm.rsize_);
    
    bvector_type_ptr* rtmp = bv_rows_;
    bv_rows_ = bbm.bv_rows_;
    bbm.bv_rows_ = rtmp;
}

//---------------------------------------------------------------------

template<typename BV>
typename basic_bmatrix<BV>::bvector_type_ptr
basic_bmatrix<BV>::construct_row(size_type row)
{
    BM_ASSERT(row < rsize_);
    bvector_type_ptr bv = bv_rows_[row];
    if (!bv)
    {
        bv = bv_rows_[row] = construct_bvector(0);
    }
    return bv;
}

//---------------------------------------------------------------------

template<typename BV>
typename basic_bmatrix<BV>::bvector_type_ptr
basic_bmatrix<BV>::construct_row(size_type row, const bvector_type& bv_src)
{
    BM_ASSERT(row < rsize_);
    bvector_type_ptr bv = bv_rows_[row];
    if (bv)
    {
        *bv = bv_src;
    }
    else
    {
        bv = bv_rows_[row] = construct_bvector(&bv_src);
    }
    return bv;
}


//---------------------------------------------------------------------

template<typename BV>
void basic_bmatrix<BV>::destruct_row(size_type row)
{
    BM_ASSERT(row < rsize_);
    bvector_type_ptr bv = bv_rows_[row];
    if (bv)
    {
        destruct_bvector(bv);
        bv_rows_[row] = 0;
    }
}


//---------------------------------------------------------------------

template<typename BV>
typename basic_bmatrix<BV>::bvector_type*
basic_bmatrix<BV>::construct_bvector(const bvector_type* bv) const
{
    bvector_type* rbv = 0;
#ifdef BM_NO_STL   // C compatibility mode
    void* mem = ::malloc(sizeof(bvector_type));
    if (mem == 0)
    {
        BM_THROW(false, BM_ERR_BADALLOC);
    }
    rbv = bv ? new(mem) bvector_type(*bv)
             : new(mem) bvector_type(ap_.strat, ap_.glevel_len,
                                     bv_size_,
                                     alloc_);
#else
    rbv = bv ? new bvector_type(*bv)
             : new bvector_type(ap_.strat, ap_.glevel_len,
                                bv_size_,
                                alloc_);
#endif
    return rbv;
}

//---------------------------------------------------------------------

template<typename BV>
void basic_bmatrix<BV>::destruct_bvector(bvector_type* bv) const
{
#ifdef BM_NO_STL   // C compatibility mode
    bv->~TBM_bvector();
    ::free((void*)bv);
#else
    delete bv;
#endif
}

//---------------------------------------------------------------------


} // namespace

#endif
