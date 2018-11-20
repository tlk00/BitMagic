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
//    typedef const value_type&                        const_reference;
    typedef typename BV::allocator_type              allocator_type;
    typedef typename bvector_type::allocation_policy allocation_policy_type;
    typedef typename bvector_type::enumerator        bvector_enumerator_type;
    typedef typename allocator_type::allocator_pool_type allocator_pool_type;
    typedef bm::basic_bmatrix<BV>                    bmatrix_type;
    typedef base_sparse_vector<CharType, BV, MAX_STR_SIZE> parent_type;

    /*! Statistical information about  memory allocation details. */
    struct statistics : public bv_statistics
    {};


    /**
         Reference class to access elements via common [] operator
         @ingroup sv
    */
    class const_reference
    {
    public:
        const_reference(const str_sparse_vector<CharType, BV, MAX_STR_SIZE>& str_sv,
                  size_type idx) BMNOEXEPT
        : str_sv_(str_sv), idx_(idx)
        {}
        
        operator const value_type*() const
        {
            str_sv_.get(idx_, buf_, MAX_STR_SIZE);
            return &(buf_[0]);
        }

        bool operator==(const const_reference& ref) const
                                { return bool(*this) == bool(ref); }
        bool is_null() const { return str_sv_.is_null(idx_); }
    private:
        const str_sparse_vector<CharType, BV, MAX_STR_SIZE>& str_sv_;
        size_type                                            idx_;
        mutable CharType                                    buf_[MAX_STR_SIZE];
    };

    /**
         Reference class to access elements via common [] operator
         @ingroup sv
    */
    class reference
    {
    public:
        reference(str_sparse_vector<CharType, BV, MAX_STR_SIZE>& str_sv,
                  size_type idx) BMNOEXEPT
        : str_sv_(str_sv), idx_(idx)
        {}
        
        operator const value_type*() const
        {
            str_sv_.get(idx_, buf_, MAX_STR_SIZE);
            return &(buf_[0]);
        }

        reference& operator=(const reference& ref)
        {
            // TO DO: implement element copy bit by bit
            str_sv_.set(idx_, (const value_type*)ref);
            return *this;
        }

        reference& operator=(const value_type* str)
        {
            str_sv_.set(idx_, str);
            return *this;
        }
        bool operator==(const reference& ref) const
                                { return bool(*this) == bool(ref); }
        bool is_null() const { return str_sv_.is_null(idx_); }
    private:
        str_sparse_vector<CharType, BV, MAX_STR_SIZE>& str_sv_;
        size_type                                      idx_;
        mutable CharType                               buf_[MAX_STR_SIZE];
    };


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
    /*! @name String element access */
    ///@{

    /** \brief Operator to get write access to an element  */
    reference operator[](size_type idx) { return reference(*this, idx); }

    /** \brief Operator to get read access to an element  */
    const_reference operator[](size_type idx) const
                                    { return const_reference(*this, idx); }

    /*!
        \brief set specified element with bounds checking and automatic resize
        \param idx  - element index (vector auto-resized if needs to)
        \param str  - string to set (zero terminated)
    */
    void set(size_type idx, const value_type* str);

    /*!
        \brief get specified element
     
        \param idx  - element index (vector auto-resized if needs to)
        \param str  - string buffer
        \param buf_size - string buffer size
    */
    void get(size_type idx, value_type* str, size_type buf_size) const;
    
    /*!
        \brief set specified element with bounds checking and automatic resize
        \param idx  - element index (vector auto-resized if needs to)
        \param str  - string to set
                    (STL class with size() support, like basic_string)
    */
    template<typename StrType>
    void assign(size_type idx, const StrType& str)
    {
        if (idx >= this->size())
            this->size_ = idx+1;

        size_type sz = (str.size() < MAX_STR_SIZE) ? str.size() : MAX_STR_SIZE;
        if(!sz)
        {
            this->clear_value_plains_from(0, idx);
            return;
        }
        unsigned i = 0;
        for (; i < sz; ++i)
        {
            CharType ch = str[i];
            this->bmatr_.set_octet(idx, i, ch);
        } // for i
        this->bmatr_.set_octet(idx, sz, 0);
        this->clear_value_plains_from(sz*8+1, idx);
    }

    /*!
        \brief get specified string element
     
        \param idx  - element index (vector auto-resized if needs to)
        \param str  - string to get [out]
    */
    template<typename StrType>
    void get(size_type idx, StrType& str) const
    {
        str.clear();
        for (unsigned i = 0; i < MAX_STR_SIZE; ++i)
        {
            CharType ch = CharType(this->bmatr_.get_octet(idx, i));
            if (ch == 0)
                break;
            str.push_back(ch);
        } // for i
    }


    ///@}

    // ------------------------------------------------------------
    /*! @name Memory optimization/compression                    */
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
        this->size_ = idx+1;
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
    unsigned i = 0;
    for (; i < MAX_STR_SIZE; ++i)
    {
        CharType ch = str[i];
        if (!ch)
        {
            this->clear_value_plains_from(i*8, idx);
            return;
        }
        this->bmatr_.set_octet(idx, i, ch);
    } // for i
}

//---------------------------------------------------------------------

template<class CharType, class BV, unsigned MAX_STR_SIZE>
void str_sparse_vector<CharType, BV, MAX_STR_SIZE>::get(
            size_type idx, value_type* str, size_type buf_size) const
{
    for (unsigned i = 0; i < MAX_STR_SIZE; ++i)
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
