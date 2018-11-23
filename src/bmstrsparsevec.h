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

/*!
   \brief sparse vector for strings with compression using bit transposition method
 
   Initial string is bit-transposed into bit-planes so collection may use less
   memory due to prefix sum compression in bit-plains.
 
   @ingroup sv
*/
template<typename CharType, typename BV, unsigned MAX_STR_SIZE>
class str_sparse_vector : public base_sparse_vector<CharType, BV, MAX_STR_SIZE>
{
public:
    typedef CharType                                 value_type;
    typedef bm::id_t                                 size_type;
    typedef BV                                       bvector_type;
    typedef bvector_type*                            bvector_type_ptr;
    typedef const bvector_type*                      bvector_type_const_ptr;
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

        size_type sz = size_type((str.size() < MAX_STR_SIZE) ? str.size() : MAX_STR_SIZE);
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
        \brief push back a string
        \param str  - string to set
                    (STL class with size() support, like basic_string)
    */
    template<typename StrType>
    void push_back(const StrType& str)
    {
        assign(this->size_, str);
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
    /*! @name Element comparison functions       */
    ///@{

    int compare(size_type idx, const value_type* str) const;

    ///@}

    
    
    // ------------------------------------------------------------
    /*! @name Size, etc       */
    ///@{

    /*! \brief return size of the vector
        \return size of sparse vector
    */
    size_type size() const { return this->size_; }
    
    /*! \brief return true if vector is empty
        \return true if empty
    */
    bool empty() const { return (size() == 0); }
    
    /*! \brief resize vector
        \param sz - new size
    */
    void resize(size_type sz) { parent_type::resize(sz); }
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

    /*!
        @brief Calculates memory statistics.

        Function fills statistics structure containing information about how
        this vector uses memory and estimation of max. amount of memory
        bvector needs to serialize itself.

        @param st - pointer on statistics structure to be filled in.

        @sa statistics
    */
    void calc_stat(struct str_sparse_vector<CharType, BV, MAX_STR_SIZE>::statistics* st) const;
    
    ///@}


    // ------------------------------------------------------------
    /*! @name Various traits                                     */
    //@{
    
    /** \brief trait if sparse vector is "compressed" (false)
    */
    static
    bool is_compressed() { return false; }

    ///@}

    /*! \brief syncronize internal structures */
    void sync(bool /*force*/) {}

    /*!
        \brief check if another sparse vector has the same content and size
     
        \param sv        - sparse vector for comparison
        \param null_able - flag to consider NULL vector in comparison (default)
                           or compare only value content plains
     
        \return true, if it is the same
    */
    bool equal(const str_sparse_vector<CharType, BV, MAX_STR_SIZE>& sv,
               bm::null_support null_able = bm::use_null) const
    {
        return parent_type::equal(sv, null_able);
    }

    /**
        \brief find position of compressed element by its rank
    */
    static
    bool find_rank(bm::id_t rank, bm::id_t& pos);
    
    /**
        \brief size of sparse vector (may be different for RSC)
    */
    size_type effective_size() const { return size(); }


protected:

    /*! \brief set value without checking boundaries */
    void set_value(size_type idx, const value_type* str);

    /*! \brief set value without checking boundaries or support of NULL */
    void set_value_no_null(size_type idx, const value_type* str);

    size_type size_internal() const { return size(); }
    void resize_internal(size_type sz) { resize(sz); }

protected:
    template<class SVect> friend class sparse_vector_serializer;
    template<class SVect> friend class sparse_vector_deserializer;
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
    for (unsigned i = 0; i < MAX_STR_SIZE; ++i)
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

template<class CharType, class BV, unsigned MAX_STR_SIZE>
void str_sparse_vector<CharType, BV, MAX_STR_SIZE>::calc_stat(
    struct str_sparse_vector<CharType, BV, MAX_STR_SIZE>::statistics* st) const
{
    BM_ASSERT(st);
    typename bvector_type::statistics stbv;
    parent_type::calc_stat(&stbv);
    
    st->reset();
    
    st->bit_blocks += stbv.bit_blocks;
    st->gap_blocks += stbv.gap_blocks;
    st->max_serialize_mem += stbv.max_serialize_mem + 8;
    st->memory_used += stbv.memory_used;
}


//---------------------------------------------------------------------

template<class CharType, class BV, unsigned MAX_STR_SIZE>
int str_sparse_vector<CharType, BV, MAX_STR_SIZE>::compare(
                     size_type idx,
                     const value_type* str) const
{
    BM_ASSERT(str);
    int res = 0;
    
    for (unsigned i = 0; i < MAX_STR_SIZE; ++i)
    {
        CharType ch = str[i];
        res = this->bmatr_.compare_octet(idx, i, ch);
        if (res || !ch)
            return res;
    } // for
    return res;
}

//---------------------------------------------------------------------

template<class CharType, class BV, unsigned MAX_STR_SIZE>
bool str_sparse_vector<CharType, BV, MAX_STR_SIZE>::find_rank(bm::id_t rank,
                                                              bm::id_t& pos)
{
    BM_ASSERT(rank);
    pos = rank - 1;
    return true;
}

//---------------------------------------------------------------------

} // namespace

#endif
