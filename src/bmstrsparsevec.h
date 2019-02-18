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
   memory due to prefix sum (GAP) compression in bit-plains.
 
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
    
    enum octet_plains
    {
        sv_octet_plains = MAX_STR_SIZE
    };
    
    typedef
    bm::heap_matrix<unsigned char,
                    MAX_STR_SIZE, // ROWS
                    256,          // COLS
                    typename bvector_type::allocator_type>
                                    plain_octet_matrix_type;

    struct is_remap_support { enum trait { value = true }; };
    struct is_rsc_support { enum trait { value = false }; };

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
        remap_flags_ = str_sv.remap_flags_;
        remap_matrix1_ = str_sv.remap_matrix1_;
        remap_matrix2_ = str_sv.remap_matrix2_;
        return *this;
    }
#ifndef BM_NO_CXX11
    /*! move-ctor */
    str_sparse_vector(str_sparse_vector<CharType, BV, MAX_STR_SIZE>&& str_sv) BMNOEXEPT
    {
        parent_type::swap(str_sv);
        remap_flags_ = str_sv.remap_flags_;
        remap_matrix1_.swap(str_sv.remap_matrix1_);
        remap_matrix2_.swap(str_sv.remap_matrix2_);
    }

    /*! move assignmment operator */
    str_sparse_vector<CharType, BV, MAX_STR_SIZE>& operator =
            (str_sparse_vector<CharType, BV, MAX_STR_SIZE>&& str_sv) BMNOEXEPT
    {
        if (this != &str_sv)
        {
            this->swap(str_sv);
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
        \brief insert the specified element
        \param idx  - element index (vector auto-resized if needs to)
        \param str  - string to set (zero terminated)
    */
    void insert(size_type idx, const value_type* str);

    /*!
        \brief get specified element
     
        \param idx  - element index
        \param str  - string buffer
        \param buf_size - string buffer size
     
        @return string length
    */
    size_type get(size_type idx, value_type* str, size_type buf_size) const;
    
    /*!
        \brief set specified element with bounds checking and automatic resize
     
        This is an equivalent of set() method, but templetized to be
        more compatible with the STL std::string and the likes
     
        \param idx  - element index (vector auto-resized if needs to)
        \param str  - input string
                      expected an STL class with size() support,
                      like basic_string<> or vector<char>
    */
    template<typename StrType>
    void assign(size_type idx, const StrType& str)
    {
        if (idx >= this->size())
            this->size_ = idx+1;

        size_type str_size = size_type(str.size());
        size_type sz = size_type((str_size < MAX_STR_SIZE) ? str_size : MAX_STR_SIZE-1);
        if (!sz)
        {
            this->clear_value_plains_from(0, idx);
            return;
        }
        unsigned i = 0;
        for (; i < sz; ++i)
        {
            CharType ch = str[i];
            if (remap_flags_) // compressional re-mapping is in effect
            {
                unsigned char remap_value = remap_matrix2_.get(i, unsigned(ch));
                BM_ASSERT(remap_value);
                ch = CharType(remap_value);
            }
            this->bmatr_.set_octet(idx, i, (unsigned char)ch);
            if (!ch)
                break;
        } // for i
        if (idx > sz)
            return;
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
        \brief push back a string (zero terminated)
        \param str  - string to set
    */
    void push_back(const value_type* str) { set(this->size_, str); }


    /*!
        \brief get specified string element
     
        Template method expects an STL-compatible type basic_string<>
     
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
            if (remap_flags_)
            {
                const unsigned char* remap_row = remap_matrix1_.row(i);
                unsigned char remap_value = remap_row[unsigned(ch)];
                BM_ASSERT(remap_value);
                if (!remap_value) // unknown dictionary element
                {
                    // TODO: add an exception throw here
                    break;
                }
                ch = CharType(remap_value);
            }
            str.push_back(ch);
        } // for i
    }

    /*! Swap content */
    void swap(str_sparse_vector& str_sv) BMNOEXEPT;

    ///@}
    
    // ------------------------------------------------------------
    /*! @name Element comparison functions       */
    ///@{

    /**
        \brief Compare vector element with argument lexicographically
     
        NOTE: for a re-mapped container, input string may have no correct
        remapping, in this case we have an ambiguity
        (we know it is not equal (0) but LT or GT?).
        Behavior is undefined.
     
        \param idx - vactor element index
        \param str - argument to compare with
     
        \return 0 - equal, < 0 - vect[i] < str, >0 otherwise
    */
    int compare(size_type idx, const value_type* str) const;
    
    
    /**
        \brief Find size of common prefix between two vector elements in octets
        \return size of common prefix
    */
    unsigned common_prefix_length(size_type idx1, size_type idx2) const;

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
    
    /*! \brief get maximum string length capacity
        \return maximum string length sparse vector can take
    */
    static size_type max_str() { return sv_octet_plains; }
    
    /*! \brief get effective string length used in vector
    
        Method returns efficiency, how close are we
        to reserved maximum.
    
        \return current string length maximum
    */
    size_type effective_max_str() const;
    
    /*! \brief get effective string length used in vector
        \return current string length maximum
    */
    size_type effective_vector_max() const { return effective_max_str(); }
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
    ///@{
    
    /** \brief trait if sparse vector is "compressed" (false)
    */
    static
    bool is_compressed() { return false; }

    ///@}

    // ------------------------------------------------------------
    /*! @name remapping, succinct utilities
        Remapping implements reduction of dit-depth thus improves
        search performance. Remapping limits farther modifications
        of sparse vector.
    */
    ///@{
    
    /**
        Get remapping status (true|false)
    */
    bool is_remap() const { return remap_flags_ != 0; }
    
    /**
        Build remapping profile and load content from another sparse vector
        \param str_sv - source sparse vector (assumed it is not remapped)
    */
    void remap_from(const str_sparse_vector& str_sv);

    /*!
        Calculate flags which octets are present on each byte-plain.
        @internal
    */
    void calc_octet_stat(plain_octet_matrix_type& octet_matrix) const;

    static
    void build_octet_remap(
                plain_octet_matrix_type& octet_remap_matrix1,
                plain_octet_matrix_type& octet_remap_matrix2,
                const plain_octet_matrix_type& octet_occupancy_matrix);
    /*!
        remap string from external (ASCII) system to matrix internal code
        @return true if remapping was ok, false if found incorrect value
                for the plain
        @internal
    */
    static
    bool remap_tosv(value_type*       sv_str,
                    size_type         buf_size,
                    const value_type* str,
                    const plain_octet_matrix_type& octet_remap_matrix2);
    
    /*!
        remap string from external (ASCII) system to matrix internal code
        @internal
    */
    bool remap_tosv(value_type*       sv_str,
                    size_type         buf_size,
                    const value_type* str) const
    {
        return remap_tosv(sv_str, buf_size, str, remap_matrix2_);
    }
    /*!
        remap string from internal code to external (ASCII) system
        @return true if remapping was ok, false if found incorrect value
                for the plain
        @internal
    */
    static
    bool remap_fromsv(value_type*       str,
                      size_type         buf_size,
                      const value_type* sv_str,
                      const plain_octet_matrix_type& octet_remap_matrix1);
    /*!
        re-calculate remap matrix2 based on matrix1
        @internal
    */
    void recalc_remap_matrix2();

    ///@}


    // ------------------------------------------------------------

    /*! \brief syncronize internal structures */
    void sync(bool force);

    /*!
        \brief check if another sparse vector has the same content and size
     
        \param sv        - sparse vector for comparison
        \param null_able - flag to consider NULL vector in comparison (default)
                           or compare only value content plains
     
        \return true, if it is the same
    */
    bool equal(const str_sparse_vector<CharType, BV, MAX_STR_SIZE>& sv,
               bm::null_support null_able = bm::use_null) const;

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

    /*! \brief insert value without checking boundaries */
    void insert_value(size_type idx, const value_type* str);

    /*! \brief insert value without checking boundaries or support of NULL */
    void insert_value_no_null(size_type idx, const value_type* str);


    size_type size_internal() const { return size(); }
    void resize_internal(size_type sz) { resize(sz); }

    size_t remap_size() const { return remap_matrix1_.get_buffer().size(); }
    const unsigned char* get_remap_buffer() const
                { return remap_matrix1_.get_buffer().buf(); }
    unsigned char* init_remap_buffer()
    {
        remap_matrix1_.init();
        return remap_matrix1_.get_buffer().data();
    }
    void set_remap() { remap_flags_ = 1; }
    
protected:
    template<class SVect> friend class sparse_vector_serializer;
    template<class SVect> friend class sparse_vector_deserializer;
    
protected:
    unsigned                 remap_flags_;   ///< remapping status
    plain_octet_matrix_type  remap_matrix1_; ///< octet remap table 1
    plain_octet_matrix_type  remap_matrix2_; ///< octet remap table 2
};

//---------------------------------------------------------------------
//---------------------------------------------------------------------


template<class CharType, class BV, unsigned MAX_STR_SIZE>
str_sparse_vector<CharType, BV, MAX_STR_SIZE>::str_sparse_vector(
        bm::null_support null_able,
        allocation_policy_type  ap,
        size_type               bv_max_size,
        const allocator_type&   alloc)
: parent_type(null_able, ap, bv_max_size, alloc),
  remap_flags_(0)
{}


//---------------------------------------------------------------------

template<class CharType, class BV, unsigned MAX_STR_SIZE>
str_sparse_vector<CharType, BV, MAX_STR_SIZE>::str_sparse_vector(
                                        const str_sparse_vector& str_sv)
: parent_type(str_sv),
  remap_flags_(str_sv.remap_flags_),
  remap_matrix1_(str_sv.remap_matrix1_),
  remap_matrix2_(str_sv.remap_matrix2_)
{}

//---------------------------------------------------------------------

template<class CharType, class BV, unsigned MAX_STR_SIZE>
void str_sparse_vector<CharType, BV, MAX_STR_SIZE>::swap(str_sparse_vector& str_sv) BMNOEXEPT
{
    parent_type::swap(str_sv);
    bm::xor_swap(remap_flags_, str_sv.remap_flags_);
    remap_matrix1_.swap(str_sv.remap_matrix1_);
    remap_matrix2_.swap(str_sv.remap_matrix2_);
}

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
void str_sparse_vector<CharType, BV, MAX_STR_SIZE>::insert(
                                size_type idx, const value_type* str)
{
    if (idx >= this->size())
    {
        this->size_ = idx+1;
        set_value(idx, str);
        return;
    }
    insert_value(idx, str);
    this->size_++;
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
        
        if (remap_flags_) // compressional re-mapping is in effect
        {
            unsigned char remap_value = remap_matrix2_.get(i, unsigned(ch));
            BM_ASSERT(remap_value);
            if (!remap_value) // unknown dictionary element
            {
                this->clear_value_plains_from(i*8, idx);
                return;
            }
            ch = CharType(remap_value);
        }
        this->bmatr_.set_octet(idx, i, (unsigned char)ch);
    } // for i
}

//---------------------------------------------------------------------

template<class CharType, class BV, unsigned MAX_STR_SIZE>
void str_sparse_vector<CharType, BV, MAX_STR_SIZE>::insert_value(
                                    size_type idx, const value_type* str)
{
    insert_value_no_null(idx, str);
    bvector_type* bv_null = this->get_null_bvect();
    if (bv_null)
        bv_null->insert(idx, true);
}

//---------------------------------------------------------------------

template<class CharType, class BV, unsigned MAX_STR_SIZE>
void str_sparse_vector<CharType, BV, MAX_STR_SIZE>::insert_value_no_null(
                                        size_type idx, const value_type* str)
{
    for (unsigned i = 0; i < MAX_STR_SIZE; ++i)
    {
        CharType ch = str[i];
        if (!ch)
        {
            this->insert_clear_value_plains_from(i*8, idx);
            return;
        }
        
        if (remap_flags_) // compressional re-mapping is in effect
        {
            unsigned char remap_value = remap_matrix2_.get(i, unsigned(ch));
            BM_ASSERT(remap_value);
            if (!remap_value) // unknown dictionary element
            {
                this->insert_clear_value_plains_from(i*8, idx);
                return;
            }
            ch = CharType(remap_value);
        }
        this->bmatr_.insert_octet(idx, i, (unsigned char)ch);
    } // for i
}


//---------------------------------------------------------------------

template<class CharType, class BV, unsigned MAX_STR_SIZE>
typename str_sparse_vector<CharType, BV, MAX_STR_SIZE>::size_type
str_sparse_vector<CharType, BV, MAX_STR_SIZE>::get(
            size_type idx, value_type* str, size_type buf_size) const
{
    size_type i = 0;
    for (; i < MAX_STR_SIZE; ++i)
    {
        if (i < buf_size)
            str[i] = 0;
        else
            break;
        CharType ch = CharType(this->bmatr_.get_octet(idx, i));
        if (ch == 0)
        {
            str[i] = 0;
            break;
        }
        str[i] = ch;
    } // for i
    if (remap_flags_)
    {
        remap_matrix1_.remap(str, i);
    }
    return i;
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
    st->ptr_sub_blocks += stbv.ptr_sub_blocks;
    
    st->max_serialize_mem += stbv.max_serialize_mem + 8;
    st->memory_used += stbv.memory_used;
    
    size_t remap_mem_usage = sizeof(remap_flags_);
    remap_mem_usage += remap_matrix1_.get_buffer().mem_usage();
    remap_mem_usage += remap_matrix2_.get_buffer().mem_usage();

    st->memory_used += remap_mem_usage;
    if (remap_flags_) // use of remapping requires some extra storage
    {
        st->max_serialize_mem += remap_mem_usage;
    }
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
        if (remap_flags_ && ch)
        {
            unsigned char remap_value = remap_matrix2_.get(i, unsigned(ch));
            if (!remap_value) // unknown dictionary element
            {
                return -1; // TODO: what would be the best return value here...
            }
            ch = CharType(remap_value);
        }

        res = this->bmatr_.compare_octet(idx, i, ch);
        if (res || !ch)
            break;
    } // for
    return res;
}

//---------------------------------------------------------------------

template<class CharType, class BV, unsigned MAX_STR_SIZE>
unsigned str_sparse_vector<CharType, BV, MAX_STR_SIZE>::common_prefix_length(
                                          size_type idx1, size_type idx2) const
{
    unsigned i = 0;
    for (; i < MAX_STR_SIZE; ++i)
    {
        CharType ch1 = CharType(this->bmatr_.get_octet(idx1, i));
        CharType ch2 = CharType(this->bmatr_.get_octet(idx2, i));
        if (!ch1 || !ch2)
        {
            if (i) 
                --i;
            break;
        }
        if (ch1 != ch2)
        {
            break;
        }
    } // for

    return i;
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

template<class CharType, class BV, unsigned MAX_STR_SIZE>
typename str_sparse_vector<CharType, BV, MAX_STR_SIZE>::size_type
str_sparse_vector<CharType, BV, MAX_STR_SIZE>::effective_max_str() const
{
    for (int i = MAX_STR_SIZE-1; i >= 0; --i)
    {
        unsigned octet_plain = unsigned(i) * unsigned(sizeof(CharType)) * 8;
        for (unsigned j = 0; j < sizeof(CharType) * 8; ++j)
        {
            if (this->bmatr_.row(octet_plain+j))
                return unsigned(i)+1;
        } // for j
    } // for i
    return 0;
}

//---------------------------------------------------------------------

template<class CharType, class BV, unsigned MAX_STR_SIZE>
void str_sparse_vector<CharType, BV, MAX_STR_SIZE>::calc_octet_stat(
                    plain_octet_matrix_type& octet_matrix) const
{
    octet_matrix.init();
    octet_matrix.set_zero();
    
    size_type size = this->size();
    
    for (unsigned i = 0; i < MAX_STR_SIZE; ++i)
    {
        unsigned char* row = octet_matrix.row(i);
        
        // TODO: optimize partial transposition
        for (size_type j = 0; j < size; ++j)
        {
            unsigned char ch = this->bmatr_.get_octet(j, i);
            unsigned k = ch;
            if (k)
                row[k] = 1;
        } // for j
    } // for i
}

//---------------------------------------------------------------------

template<class CharType, class BV, unsigned MAX_STR_SIZE>
void str_sparse_vector<CharType, BV, MAX_STR_SIZE>::build_octet_remap(
            plain_octet_matrix_type& octet_remap_matrix1,
            plain_octet_matrix_type& octet_remap_matrix2,
            const plain_octet_matrix_type& octet_occupancy_matrix)
{
    octet_remap_matrix1.init();
    octet_remap_matrix1.set_zero();
    octet_remap_matrix2.init();
    octet_remap_matrix2.set_zero();

    for (unsigned i = 0; i < octet_occupancy_matrix.rows(); ++i)
    {
        const unsigned char* row = octet_occupancy_matrix.row(i);
        unsigned char* remap_row1 = octet_remap_matrix1.row(i);
        unsigned char* remap_row2 = octet_remap_matrix2.row(i);
        unsigned count = 1;
        for (unsigned j = 1; j < octet_occupancy_matrix.cols(); ++j)
        {
            if (row[j]) // octet is present
            {
                // set two remapping table values
                remap_row1[count] = (unsigned char)j;
                remap_row2[j] = (unsigned char)count;
                ++count;
                BM_ASSERT(count < 256);
            }
        } // for j
    } // for i
}

//---------------------------------------------------------------------

template<class CharType, class BV, unsigned MAX_STR_SIZE>
void str_sparse_vector<CharType, BV, MAX_STR_SIZE>::recalc_remap_matrix2()
{
    BM_ASSERT(remap_flags_);
    
    remap_matrix2_.init();
    remap_matrix2_.set_zero();
    
    for (unsigned i = 0; i < remap_matrix1_.rows(); ++i)
    {
        const unsigned char* remap_row1 = remap_matrix1_.row(i);
              unsigned char* remap_row2 = remap_matrix2_.row(i);
        for (unsigned j = 1; j < remap_matrix1_.cols(); ++j)
        {
            if (remap_row1[j])
            {
                unsigned count = remap_row1[j];
                remap_row2[count] = (unsigned char)j;
                BM_ASSERT(count < 256);
            }
        } // for j
    } // for i
}

//---------------------------------------------------------------------

template<class CharType, class BV, unsigned MAX_STR_SIZE>
bool str_sparse_vector<CharType, BV, MAX_STR_SIZE>::remap_tosv(
                   value_type*       sv_str,
                   size_type         buf_size,
                   const value_type* str,
                   const plain_octet_matrix_type& octet_remap_matrix2)
{
    for (unsigned i = 0; i < buf_size; ++i)
    {
        CharType ch = str[i];
        if (ch == 0)
        {
            sv_str[i] = ch;
            break;
        }
        const unsigned char* remap_row = octet_remap_matrix2.row(i);
        unsigned char remap_value = remap_row[unsigned(ch)];
        if (!remap_value) // unknown dictionary element
        {
            return false;
        }
        sv_str[i] = CharType(remap_value);
    } // for i
    return true;
}

//---------------------------------------------------------------------

template<class CharType, class BV, unsigned MAX_STR_SIZE>
bool str_sparse_vector<CharType, BV, MAX_STR_SIZE>::remap_fromsv(
                         value_type* str,
                         size_type         buf_size,
                         const value_type* sv_str,
                         const plain_octet_matrix_type& octet_remap_matrix1)
{
    for (unsigned i = 0; i < buf_size; ++i)
    {
        CharType ch = sv_str[i];
        if (ch == 0)
        {
            str[i] = ch;
            break;
        }
        const unsigned char* remap_row = octet_remap_matrix1.row(i);
        unsigned char remap_value = remap_row[unsigned(ch)];
        if (!remap_value) // unknown dictionary element
        {
            return false;
        }
        str[i] = CharType(remap_value);
    } // for i
    return true;
}

//---------------------------------------------------------------------

template<class CharType, class BV, unsigned MAX_STR_SIZE>
void str_sparse_vector<CharType, BV, MAX_STR_SIZE>::remap_from(const str_sparse_vector& str_sv)
{
    if (str_sv.is_remap())
    {
        *this = str_sv;
        return;
    }
    this->clear();
    if (str_sv.empty()) // no content to remap
    {
        return;
    }
    
    plain_octet_matrix_type omatrix; // occupancy map
    str_sv.calc_octet_stat(omatrix);
    
    str_sv.build_octet_remap(remap_matrix1_, remap_matrix2_, omatrix);
    remap_flags_ = 1; // turn ON remapped mode
    
    // load content
    // TODO: optimization (current implementation is a naive "get-set")
    value_type temp_str[MAX_STR_SIZE+1];
    for (size_type i = 0; i < str_sv.size(); ++i)
    {
        str_sv.get(i, &temp_str[0], MAX_STR_SIZE);
        this->set(i, &temp_str[0]);
    } // for i
}

//---------------------------------------------------------------------

template<class CharType, class BV, unsigned MAX_STR_SIZE>
void str_sparse_vector<CharType, BV, MAX_STR_SIZE>::sync(bool /*force*/)
{
    if (remap_flags_)
    {
        recalc_remap_matrix2();
    }
}

//---------------------------------------------------------------------

template<class CharType, class BV, unsigned MAX_STR_SIZE>
bool str_sparse_vector<CharType, BV, MAX_STR_SIZE>::equal(
                const str_sparse_vector<CharType, BV, MAX_STR_SIZE>& sv,
                bm::null_support null_able) const
{
    // at this point both vectors should have the same remap settings
    // to be considered "equal".
    // Strictly speaking this is incorrect, because re-map and non-remap
    // vectors may have the same content

    if (remap_flags_ != sv.remap_flags_)
        return false;
    if (remap_flags_)
    {
        bool b;
        b = remap_matrix1_.get_buffer().equal(sv.remap_matrix1_.get_buffer());
        if (!b)
            return b;
        b = remap_matrix2_.get_buffer().equal(sv.remap_matrix2_.get_buffer());
        if (!b)
            return b;
    }
    return parent_type::equal(sv, null_able);
}

//---------------------------------------------------------------------

} // namespace

#endif
