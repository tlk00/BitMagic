#ifndef BMSTRSPARSEVEC__H__INCLUDED__
#define BMSTRSPARSEVEC__H__INCLUDED__
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

/*! \file bmstrsparsevec.h
    \brief string sparse vector based on bit-transposed matrix
*/

#include <stddef.h>
#include "bmconst.h"

#ifndef BM_NO_STL
#include <stdexcept>
#endif

#include "bm.h"
#include "bmtrans.h"
#include "bmalgo.h"
#include "bmbuffer.h"
#include "bmbmatrix.h"
#include "bmdef.h"

namespace bm
{

template<typename CharType, typename BV, unsigned MAX_STR_SIZE>
class str_sparse_vector : public base_sparse_vector<CharType, BV, MAX_STR_SIZE>
{
public:
    typedef CharType                                 value_type;
    typedef bm::id_t                                 size_type;
    typedef BV                                       bvector_type;
    typedef bvector_type*                            bvector_type_ptr;
    typedef const bvector_type*                      bvector_type_const_ptr;
    typedef const value_type&                        const_reference;
    typedef typename BV::allocator_type              allocator_type;
    typedef typename bvector_type::allocation_policy allocation_policy_type;
    typedef typename bvector_type::enumerator        bvector_enumerator_type;
    typedef typename allocator_type::allocator_pool_type allocator_pool_type;
    typedef bm::basic_bmatrix<BV>                    bmatrix_type;
    typedef base_sparse_vector<CharType, BV, MAX_STR_SIZE> parent_type;

    /*! Statistical information about  memory allocation details. */
    struct statistics : public bv_statistics
    {};

public:

    /*!
        \brief Sparse vector constructor
     
        \param null_able - defines if vector supports NULL values flag
            by default it is OFF, use bm::use_null to enable it
        \param ap - allocation strategy for underlying bit-vectors
        Default allocation policy uses BM_BIT setting (fastest access)
        \param bv_max_size - maximum possible size of underlying bit-vectors
        Please note, this is NOT size of svector itself, it is dynamic upper limit
        which should be used very carefully if we surely know the ultimate size
        \param alloc - allocator for bit-vectors
     
        \sa bvector<>
        \sa bm::bvector<>::allocation_policy
        \sa bm::startegy
    */
    str_sparse_vector(bm::null_support null_able = bm::no_null,
                      allocation_policy_type ap = allocation_policy_type(),
                      size_type bv_max_size = bm::id_max,
                      const allocator_type&   alloc  = allocator_type());

    /*! copy-ctor */
    str_sparse_vector(const str_sparse_vector& str_sv);
    
    /*! copy assignmment operator */
    str_sparse_vector<CharType, BV, MAX_STR_SIZE>& operator = (
                const str_sparse_vector<CharType, BV, MAX_STR_SIZE>& str_sv)
    {
        if (this != &str_sv)
            parent_type::copy_from(str_sv);
        return *this;
    }
#ifndef BM_NO_CXX11
    /*! move-ctor */
    str_sparse_vector(str_sparse_vector<CharType, BV, MAX_STR_SIZE>&& str_sv) BMNOEXEPT
    {
        parent_type::swap(str_sv);
    }

    /*! move assignmment operator */
    str_sparse_vector<CharType, BV, MAX_STR_SIZE>& operator =
            (str_sparse_vector<CharType, BV, MAX_STR_SIZE>&& str_sv) BMNOEXEPT
    {
        if (this != &str_sv)
        {
            swap(str_sv);
        }
        return *this;
    }
#endif

public:

    // ------------------------------------------------------------
    /*! @name Element access */
    ///@{

    /*!
        \brief set specified element with bounds checking and automatic resize
        \param idx  - element index
        \param str  - string to set (zero terminated)
    */
    void set(size_type idx, const value_type* str);

    void get(size_type idx, value_type* str, size_type buf_size) const;
    
    ///@}

    // ------------------------------------------------------------
    /*! @name Memory optimization                                */
    ///@{

    /*!
        \brief run memory optimization for all vector plains
        \param temp_block - pre-allocated memory block to avoid unnecessary re-allocs
        \param opt_mode - requested compression depth
        \param stat - memory allocation statistics after optimization
    */
    void optimize(bm::word_t* temp_block = 0,
                  typename bvector_type::optmode opt_mode = bvector_type::opt_compress,
                  typename str_sparse_vector<CharType, BV, MAX_STR_SIZE>::statistics* stat = 0);

    ///@}

protected:

    /*! \brief set value without checking boundaries */
    void set_value(size_type idx, const value_type* str);

    /*! \brief set value without checking boundaries or support of NULL */
    void set_value_no_null(size_type idx, const value_type* str);

};

//---------------------------------------------------------------------
//---------------------------------------------------------------------


template<class CharType, class BV, unsigned MAX_STR_SIZE>
str_sparse_vector<CharType, BV, MAX_STR_SIZE>::str_sparse_vector(
        bm::null_support null_able,
        allocation_policy_type  ap,
        size_type               bv_max_size,
        const allocator_type&   alloc)
: parent_type(null_able, ap, bv_max_size, alloc)
{}


//---------------------------------------------------------------------

template<class CharType, class BV, unsigned MAX_STR_SIZE>
str_sparse_vector<CharType, BV, MAX_STR_SIZE>::str_sparse_vector(
                                        const str_sparse_vector& str_sv)
: parent_type(str_sv)
{}

//---------------------------------------------------------------------


template<class CharType, class BV, unsigned MAX_STR_SIZE>
void str_sparse_vector<CharType, BV, MAX_STR_SIZE>::set(
                                size_type idx, const value_type* str)
{
    if (idx >= this->size())
    {
        this->size_ = idx+1;
    }
    set_value(idx, str);
}

//---------------------------------------------------------------------

template<class CharType, class BV, unsigned MAX_STR_SIZE>
void str_sparse_vector<CharType, BV, MAX_STR_SIZE>::set_value(
                                size_type idx, const value_type* str)
{
    set_value_no_null(idx, str);
    bvector_type* bv_null = this->get_null_bvect();
    if (bv_null)
        bv_null->set_bit_no_check(idx);
}

//---------------------------------------------------------------------

template<class CharType, class BV, unsigned MAX_STR_SIZE>
void str_sparse_vector<CharType, BV, MAX_STR_SIZE>::set_value_no_null(
                                size_type idx, const value_type* str)
{
    // calculate logical block coordinates and masks
    //
    unsigned i = 0;
    for (; i < MAX_STR_SIZE; ++i)
    {
        CharType ch = str[i];
        if (!ch)
            break;
        this->bmatr_.set_octet(idx, i, ch);
    } // for i
    
    // clear all the extra plains above string size
    i = i * 8 + 1;
    for (; i < parent_type::sv_value_plains; ++i)
    {
        bvector_type* bv = this->bmatr_.get_row(i);
        if (bv)
            bv->clear_bit_no_check(idx);
    }
}

//---------------------------------------------------------------------

template<class CharType, class BV, unsigned MAX_STR_SIZE>
void str_sparse_vector<CharType, BV, MAX_STR_SIZE>::get(
            size_type idx, value_type* str, size_type buf_size) const
{
    unsigned i = 0;
    for (; i < MAX_STR_SIZE; ++i)
    {
        if (i < buf_size)
            str[i] = 0;
        else
            break;
        CharType ch = CharType(this->bmatr_.get_octet(idx, i));
        if (ch == 0)
            break;
        str[i] = ch;
    } // for i
}

//---------------------------------------------------------------------

template<class CharType, class BV, unsigned MAX_STR_SIZE>
void str_sparse_vector<CharType, BV, MAX_STR_SIZE>::optimize(
      bm::word_t* temp_block,
      typename bvector_type::optmode opt_mode,
      typename str_sparse_vector<CharType, BV, MAX_STR_SIZE>::statistics* st)
{
    typename bvector_type::statistics stbv;
    parent_type::optimize(temp_block, opt_mode, &stbv);
    
    if (st)
    {
        st->bit_blocks += stbv.bit_blocks;
        st->gap_blocks += stbv.gap_blocks;
        st->max_serialize_mem += stbv.max_serialize_mem + 8;
        st->memory_used += stbv.memory_used;
    }
}

//---------------------------------------------------------------------


} // namespace

#endif
