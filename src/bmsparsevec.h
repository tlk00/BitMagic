#ifndef BMSPARSEVEC_H__INCLUDED__
#define BMSPARSEVEC_H__INCLUDED__
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

/*! \file bmsparsevec.h
    \brief Sparse constainer sparse_vector<> for integer types using bit-transposition transform
*/

#include <memory.h>

#ifndef BM_NO_STL
#include <stdexcept>
#endif

#include "bm.h"
#include "bmtrans.h"
#include "bmalgo.h"
#include "bmbuffer.h"
#include "bmdef.h"

namespace bm
{

/** \defgroup svector Sparse and compressed vectors
    Sparse vector for integer types using bit transposition transform
 
    @ingroup bmagic
 */


/** \defgroup sv sparse_vector<>
    Sparse vector for integer types using bit transposition transform
 
    @ingroup bmagic
 */

// forward declarations
//
template<class Val, class SV> class compressed_sparse_vector;


/*!
   \brief sparse vector with runtime compression using bit transposition method
 
   Sparse vector implements variable bit-depth storage model.
   Initial data is bit-transposed into bit-planes so initial each element
   may use less memory than the original native data type prescribes.
   For example, 32-bit integer may only use 20 bits.
 
   Another level of compression is provided by bit-vector (BV template parameter)
   used for storing bit planes. bvector<> implements varians of on the fly block
   compression, so if a significant area of a sparse vector uses less bits - it
   will save memory.
 
   Overall it provides variable bit-depth compression, sparse compression in
   bit-plains.
 
   @ingroup sv
*/
template<class Val, class BV>
class sparse_vector
{
public:

    enum bit_plains
    {
        sv_plains = (sizeof(Val) * 8 + 1),
        sv_value_plains = (sizeof(Val) * 8)
    };

    typedef Val                                      value_type;
    typedef bm::id_t                                 size_type;
    typedef BV                                       bvector_type;
    typedef bvector_type*                            bvector_type_ptr;
	typedef const value_type&                        const_reference;
    typedef typename BV::allocator_type              allocator_type;
    typedef typename bvector_type::allocation_policy allocation_policy_type;
    typedef typename bvector_type::enumerator        bvector_enumerator_type;
    typedef typename allocator_type::allocator_pool_type allocator_pool_type;
    

    /*! Statistical information about  memory allocation details. */
    struct statistics : public bv_statistics
    {};

    /**
         Reference class to access elements via common [] operator
         @ingroup sv
    */
    class reference
    {
    public:
        reference(sparse_vector<Val, BV>& sv, size_type idx) BMNOEXEPT
        : sv_(sv), idx_(idx)
        {}
        operator value_type() const { return sv_.get(idx_); }
        reference& operator=(const reference& ref)
        {
            sv_.set(idx_, (value_type)ref);
            return *this;
        }
        reference& operator=(value_type val)
        {
            sv_.set(idx_, val);
            return *this;
        }
        bool operator==(const reference& ref) const
                                { return bool(*this) == bool(ref); }
        bool is_null() const { return sv_.is_null(idx_); }
    private:
        sparse_vector<Val, BV>& sv_;
        size_type               idx_;
    };
    
    
    /**
        Const iterator to traverse the sparse vector.
     
        Implementation uses buffer for decoding so, competing changes
        to the original vector may not match the iterator returned values.
     
        This iterator keeps an operational buffer for 8K elements,
        so memory footprint is not negligable (about 64K for unsigned int)
     
        @ingroup sv
    */
    class const_iterator
    {
    public:
    friend class sparse_vector;

#ifndef BM_NO_STL
        typedef std::input_iterator_tag  iterator_category;
#endif
        typedef sparse_vector<Val, BV>                     sparse_vector_type;
        typedef sparse_vector_type*                        sparse_vector_type_ptr;
        typedef typename sparse_vector_type::value_type    value_type;
        typedef typename sparse_vector_type::size_type     size_type;
        typedef typename sparse_vector_type::bvector_type  bvector_type;
        typedef typename bvector_type::allocator_type      allocator_type;
        typedef typename bvector_type::allocator_type::allocator_pool_type allocator_pool_type;
        typedef bm::byte_buffer<allocator_type>            buffer_type;

        typedef unsigned                    difference_type;
        typedef unsigned*                   pointer;
        typedef value_type&                 reference;
        

    public:
        const_iterator();
        const_iterator(const sparse_vector_type* sv);
        const_iterator(const sparse_vector_type* sv, bm::id_t pos);
        const_iterator(const const_iterator& it);
        

        bool operator==(const const_iterator& it) const
                                { return (pos_ == it.pos_) && (sv_ == it.sv_); }
        bool operator!=(const const_iterator& it) const
                                { return ! operator==(it); }
        bool operator < (const const_iterator& it) const
                                { return pos_ < it.pos_; }
        bool operator <= (const const_iterator& it) const
                                { return pos_ <= it.pos_; }
        bool operator > (const const_iterator& it) const
                                { return pos_ > it.pos_; }
        bool operator >= (const const_iterator& it) const
                                { return pos_ >= it.pos_; }

        /*! \brief Get current position (value) */
        value_type operator*() const { return this->value(); }
        
        
        /*! \brief Advance to the next available value */
        const_iterator& operator++() { this->advance(); return *this; }

        /*! \brief Advance to the next available value */
        const_iterator& operator++(int)
            { const_iterator tmp(*this);this->advance(); return tmp; }


        /*! \brief Get current position (value) */
        value_type value() const;
        
        /*! \brief Get NULL status */
        bool is_null() const;
        
        /// Returns true if iterator is at a valid position
        bool valid() const { return pos_ != bm::id_max; }
        
        /// Invalidate current iterator
        void invalidate() { pos_ = bm::id_max; }
        
        /// Current position (index) in the vector
        bm::id_t pos() const { return pos_; }
        
        /// re-position to a specified position
        void go_to(bm::id_t pos);
        
        /// advance iterator forward by one
        void advance();
        
        void skip_zero_values();
    private:
        enum buf_size_e
        {
            n_buf_size = 1024 * 8
            //n_buf_size = 65535 * 10
        };
        
    private:
        const bm::sparse_vector<Val, BV>* sv_;
        bm::id_t                          pos_;     ///!< Position
        mutable buffer_type               buffer_;  ///!< value buffer
        mutable value_type*               buf_ptr_; ///!< position in the buffer
        mutable allocator_pool_type       pool_;
    };
    
    /**
        Back insert iterator implements buffered insert, faster than generic
        access assignment.
     
        Limitations for bufferen inserter:
        1. Do not use more than one inserter competitively
        2. Use method flush() at the end to send the rest of accumulated buffer
        flush is happening automatically on destruction, but if flush produces an
        exception (for whatever reason) it will be an exception in destructor.
        As such, explicit flush() is safer way to finilize the sparse vector load.

        @ingroup sv
    */
    class back_insert_iterator
    {
    public:
#ifndef BM_NO_STL
        typedef std::output_iterator_tag  iterator_category;
#endif
        typedef sparse_vector<Val, BV>                     sparse_vector_type;
        typedef sparse_vector_type*                        sparse_vector_type_ptr;
        typedef typename sparse_vector_type::value_type    value_type;
        typedef typename sparse_vector_type::size_type     size_type;
        typedef typename sparse_vector_type::bvector_type  bvector_type;
        typedef typename bvector_type::allocator_type      allocator_type;
        typedef typename bvector_type::allocator_type::allocator_pool_type allocator_pool_type;
        typedef bm::byte_buffer<allocator_type>            buffer_type;

        typedef void difference_type;
        typedef void pointer;
        typedef void reference;
        
    public:
        back_insert_iterator();
        back_insert_iterator(sparse_vector_type* sv);
        back_insert_iterator(const back_insert_iterator& bi);
        
        back_insert_iterator& operator=(const back_insert_iterator& bi)
        {
            BM_ASSERT(bi.empty());
            this->flush(); sv_ = bi.sv_;
            return *this;
        }

        ~back_insert_iterator();
        
        /** push value to the vector */
        back_insert_iterator& operator=(value_type v) { this->add(v); return *this; }
        /** noop */
        back_insert_iterator& operator*() { return *this; }
        /** noop */
        back_insert_iterator& operator++() { return *this; }
        /** noop */
        back_insert_iterator& operator++( int ) { return *this; }
        
        /** add value to the container*/
        void add(value_type v);
        
        /** add NULL (no-value) to the container */
        void add_null();
        
        /** add a series of consequitve NULLs (no-value) to the container */
        void add_null(size_type count);

        /** return true if insertion buffer is empty */
        bool empty() const;
        
        /** flush the accumulated buffer */
        void flush();
    protected:
    
        /** add value to the buffer without touyching the NULL vector
            @param v - value to push back
            @return index of added value in the internal buffer
            @internal
        */
        unsigned add_value(value_type v);
        
    private:
        enum buf_size_e
        {
            n_buf_size = 1024 * 8
        };
    private:
        bm::sparse_vector<Val, BV>* sv_;      ///!< pointer on the parent vector
        bvector_type*               bv_null_; ///!< not NULL vector pointer
        buffer_type                 buffer_;  ///!< value buffer
        value_type*                 buf_ptr_; ///!< position in the buffer
    };
    
    friend const_iterator;
    friend back_insert_iterator;

public:
    // ------------------------------------------------------------
    /*! @name Construction and assignment  */
    ///@{

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
    sparse_vector(bm::null_support null_able = bm::no_null,
                  allocation_policy_type ap = allocation_policy_type(),
                  size_type bv_max_size = bm::id_max,
                  const allocator_type&   alloc  = allocator_type());
    
    /*! copy-ctor */
    sparse_vector(const sparse_vector<Val, BV>& sv);

    /*! copy assignmment operator */
    sparse_vector<Val,BV>& operator = (const sparse_vector<Val, BV>& sv)
    {
        if (this != &sv)
        {
            clear();
            resize(sv.size());
            bv_size_ = sv.bv_size_;
            alloc_ = sv.alloc_;
            effective_plains_ = sv.effective_plains_;

            for (size_type i = 0; i < sv.stored_plains(); ++i)
            {
                const bvector_type* bv = sv.plains_[i];
                plains_[i] = bv ? construct_bvector(bv) : 0;
            } // for i
        }
        return *this;
    }

#ifndef BM_NO_CXX11
    /*! move-ctor */
    sparse_vector(sparse_vector<Val, BV>&& sv) BMNOEXEPT;


    /*! move assignmment operator */
    sparse_vector<Val,BV>& operator = (sparse_vector<Val, BV>&& sv) BMNOEXEPT
    {
        if (this != &sv)
        {
            clear();
            swap(sv);
        }
        return *this;
    }
#endif

    ~sparse_vector() BMNOEXEPT;
    ///@}

    
    // ------------------------------------------------------------
    /*! @name Element access */
    ///@{

    /** \brief Operator to get write access to an element  */
    reference operator[](size_type idx) { return reference(*this, idx); }

    /*!
        \brief get specified element without bounds checking
        \param idx - element index
        \return value of the element
    */
    value_type operator[](size_type idx) const { return this->get(idx); }

    /*!
        \brief access specified element with bounds checking
        \param idx - element index
        \return value of the element
    */
    value_type at(size_type idx) const;
    /*!
        \brief get specified element without bounds checking
        \param idx - element index
        \return value of the element
    */
    value_type get(bm::id_t idx) const;

    /*!
        \brief set specified element with bounds checking and automatic resize
        \param idx - element index
        \param v   - element value
    */
    void set(size_type idx, value_type v);

    /*!
        \brief increment specified element by one
        \param idx - element index
    */
    void inc(size_type idx);

    /*!
        \brief push value back into vector
        \param v   - element value
    */
    void push_back(value_type v);

    /*!
        \brief clear specified element with bounds checking and automatic resize
        \param idx - element index
    */
    void clear(size_type idx, bool set_null = false);

    ///@}

    // ------------------------------------------------------------
    /*! @name Iterator access */
    //@{

    /** Provide const iterator access to container content  */
    const_iterator begin() const;

    /** Provide const iterator access to the end    */
    const_iterator end() const { return const_iterator(this, bm::id_max); }

    /** Get const_itertor re-positioned to specific element
    @param idx - position in the sparse vector
    */
    const_iterator get_const_iterator(size_type idx) const { return const_iterator(this, idx); }
 
    /** Provide back insert iterator
    Back insert iterator implements buffered insertion, which is faster, than random access
    or push_back
    */
    back_insert_iterator get_back_inserter() { return back_insert_iterator(this); }
    ///@}


    // ------------------------------------------------------------
    /*! @name Various traits                                     */
    //@{
    /**
        \brief check if container supports NULL(unassigned) values
    */
    bool is_nullable() const;
    
    /**
        \brief Get bit-vector of assigned values or NULL
        (if not constructed that way)
    */
    const bvector_type* get_null_bvector() const;
    
    /** \brief test if specified element is NULL
        \param idx - element index
        \return true if it is NULL false if it was assigned or container
        is not configured to support assignment flags
    */
    bool is_null(size_type idx) const;
    
    /** \brief set specified element to unassigned value (NULL)
        \param idx - element index
    */
    void set_null(size_type idx);
    
    /** \brief trait if sparse vector is "compressed" (false)
    */
    static
    bool is_compressed() { return false; }
    
    ///@}


    // ------------------------------------------------------------
    /*! @name Loading of sparse vector from C-stype array       */
    //@{
    /*!
        \brief Import list of elements from a C-style array
        \param arr  - source array
        \param size - source size
        \param offset - target index in the sparse vector
    */
    void import(const value_type* arr, size_type size, size_type offset = 0);
    
    /*!
        \brief Import list of elements from a C-style array (pushed back)
        \param arr  - source array
        \param size - source size
    */
    void import_back(const value_type* arr, size_type size);
    ///@}

    // ------------------------------------------------------------
    /*! @name Export content to C-stype array       */
    ///@{

    /*!
        \brief Bulk export list of elements to a C-style array
     
        For efficiency, this is left as a low level function,
        it does not do any bounds checking on the target array, it will
        override memory and crash if you are not careful with allocation
        and request size.
     
        \param arr  - dest array
        \param idx_from - index in the sparse vector to export from
        \param size - decoding size (array allocation should match)
        \param zero_mem - set to false if target array is pre-initialized
                          with 0s to avoid performance penalty
     
        \return number of actually exported elements (can be less than requested)
     
        \sa gather
    */
    size_type decode(value_type* arr,
                     size_type   idx_from,
                     size_type   size,
                     bool        zero_mem = true) const;
    
    
    /*!
        \brief Gather elements to a C-style array
     
        Gather collects values from different locations, for best
        performance feed it with sorted list of indexes.
     
        Faster than one-by-one random access.
     
        For efficiency, this is left as a low level function,
        it does not do any bounds checking on the target array, it will
        override memory and crash if you are not careful with allocation
        and request size.
     
        \param arr  - dest array
        \param idx - index list to gather elements
        \param size - decoding index list size (array allocation should match)
        \param sorted_idx - sort order directive for the idx array
                            (BM_UNSORTED, BM_SORTED, BM_UNKNOWN)
        Sort order affects both performance and correctness(!), use BM_UNKNOWN
        if not sure.
     
        \return number of actually exported elements (can be less than requested)
     
        \sa decode
    */
    size_type gather(value_type* arr,
                     const size_type* idx,
                     size_type   size,
                     bm::sort_order sorted_idx) const;
    ///@}

    /*! \brief content exchange
    */
    void swap(sparse_vector<Val, BV>& sv) BMNOEXEPT;

    // ------------------------------------------------------------
    /*! @name Clear                                              */
    ///@{

    /*! \brief resize to zero, free memory */
    void clear() BMNOEXEPT;

    /*!
        \brief clear range (assign bit 0 for all plains)
        \param left  - interval start
        \param right - interval end (closed interval)
        \param set_null - set cleared values to unassigned (NULL)
    */
    sparse_vector<Val, BV>& clear_range(size_type left, size_type right, bool set_null = false);

    ///@}

    // ------------------------------------------------------------
    /*! @name Size, etc       */
    ///@{

    /*! \brief return size of the vector
        \return size of sparse vector
    */
    size_type size() const;
    
    /*! \brief return true if vector is empty
        \return true if empty
    */
    bool empty() const;
    
    /*! \brief resize vector
        \param sz - new size
    */
    void resize(size_type sz);
    ///@}
        
    // ------------------------------------------------------------
    /*! @name Comparison       */
    ///@{

    /*!
        \brief check if another sparse vector has the same content and size
     
        \param sv        - sparse vector for comparison
        \param null_able - flag to consider NULL vector in comparison (default)
                           or compare only value content plains
     
        \return true, if it is the same
    */
    bool equal(const sparse_vector<Val, BV>& sv,
               bm::null_support null_able = bm::use_null) const;
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
                  typename sparse_vector<Val, BV>::statistics* stat = 0);
    /*!
       \brief Optimize sizes of GAP blocks

       This method runs an analysis to find optimal GAP levels for all bit plains
       of the vector.
    */
    void optimize_gap_size();

    /*!
        @brief Calculates memory statistics.

        Function fills statistics structure containing information about how
        this vector uses memory and estimation of max. amount of memory
        bvector needs to serialize itself.

        @param st - pointer on statistics structure to be filled in.

        @sa statistics
    */
    void calc_stat(struct sparse_vector<Val, BV>::statistics* st) const;
    ///@}

    // ------------------------------------------------------------
    /*! @name Merge, split, partition data                        */
    ///@{

    /*!
        \brief join all with another sparse vector using OR operation
        \param sv - argument vector to join with
        \return slf reference
    */
    sparse_vector<Val, BV>& join(const sparse_vector<Val, BV>& sv);

    /**
        @brief copy range of values from another sparse vector
     
        Copy [left..right] values from the source vector,
        clear everything outside the range.
     
        \param sv - source vector
        \param left  - index from in losed diapason of [left..right]
        \param right - index to in losed diapason of [left..right]
    */
    void copy_range(const sparse_vector<Val, BV>& sv,
                    size_type left, size_type right);
    
    /**
        @brief Apply value filter, defined by mask vector
     
        All bit-plains are ANDed against the filter mask.
    */
    void filter(const bvector_type& bv_mask);

    ///@}
    

    // ------------------------------------------------------------
    /*! @name Access to internals                                */
    ///@{

    /*!
        \brief get access to bit-plain, function checks and creates a plain
        \return bit-vector for the bit plain
    */
    bvector_type_ptr get_plain(unsigned i);

    /*!
        \brief get read-only access to bit-plain
        \return bit-vector for the bit plain or NULL
    */
    const bvector_type_ptr get_plain(unsigned i) const { return plains_[i]; }

    /*!
        \brief get total number of bit-plains in the vector
    */
    static unsigned plains() { return value_bits(); }

    /** Number of stored bit-plains (value plains + extra */
    static unsigned stored_plains() { return value_bits()+1; }

    /*!
        \brief get access to bit-plain as is (can return NULL)
    */
    bvector_type_ptr plain(unsigned i) { return plains_[i]; }
    const bvector_type_ptr plain(unsigned i) const { return plains_[i]; }
    
    /*!
        \brief free memory in bit-plain
    */
    void free_plain(unsigned i);
               
    /*!
        \brief Bulk export list of elements to a C-style array
     
        Use of all extract() methods is restricted.
        Please consider decode() for the same purpose.
     
        \param arr  - dest array
        \param size - dest size
        \param offset - target index in the sparse vector to export from
        \param zero_mem - set to false if target array is pre-initialized
                          with 0s to avoid performance penalty
     
        \return number of exported elements
     
        \sa decode
     
        @internal
    */
    size_type extract(value_type* arr,
                      size_type size,
                      size_type offset = 0,
                      bool      zero_mem = true,
                      allocator_pool_type* pool_ptr = 0) const;

    /** \brief extract small window without use of masking vector
        \sa decode
        @internal
    */
    size_type extract_range(value_type* arr,
                            size_type size,
                            size_type offset,
                            bool      zero_mem = true) const;
    
    /** \brief extract medium window without use of masking vector
        \sa decode
        @internal
    */
    size_type extract_plains(value_type* arr,
                             size_type size,
                             size_type offset,
                             bool      zero_mem = true) const;
    
    /** \brief address translation for this type of container
        \internal
    */
    static
    size_type translate_address(size_type i) { return i; }
    
    /**
        \brief throw range error
        \internal
    */
    static
    void throw_range_error(const char* err_msg);

    /**
        \brief throw bad alloc
        \internal
    */
    static
    void throw_bad_alloc();

    /** Number of effective bit-plains in the value type */
    unsigned effective_plains() const { return effective_plains_ + 1; }

    /**
    \brief find position of compressed element by its rank
    */
    bool find_rank(bm::id_t rank, bm::id_t& pos) const;

    ///@}

    /// Set allocator pool for local (non-threaded)
    /// memory cyclic(lots of alloc-free ops) opertations
    ///
    void set_allocator_pool(allocator_pool_type* pool_ptr);


private:

    /*! \brief free all internal vectors
    */
    void free_vectors() BMNOEXEPT;
    
    /** Number of total bit-plains in the value type*/
    static unsigned value_bits() { return sv_value_plains; }
    
    /** plain index for the "NOT NULL" flags plain */
    static unsigned null_plain() { return value_bits(); }

protected:
    /*! \brief set value without checking boundaries */
    void set_value(size_type idx, value_type v);
    
    /*! \brief set value without checking boundaries or support of NULL */
    void set_value_no_null(size_type idx, value_type v);

    /*! \brief push value back into vector without NULL semantics */
    void push_back_no_null(value_type v);


    const bm::word_t* get_block(unsigned p, unsigned i, unsigned j) const;

    bvector_type* construct_bvector(const bvector_type* bv) const;
    void destruct_bvector(bvector_type* bv) const;
    bvector_type* get_null_bvect() { return plains_[this->null_plain()]; }


    template<class V, class SVect> friend class compressed_sparse_vector;
    template<class SVect> friend class sparse_vector_scanner;

private:
    size_type                bv_size_;
    allocator_type           alloc_;
    allocation_policy_type   ap_;
    allocator_pool_type*     pool_;
    
    bvector_type_ptr         plains_[sv_plains];
    size_type                size_;
    unsigned                 effective_plains_;
};


//---------------------------------------------------------------------
//---------------------------------------------------------------------


template<class Val, class BV>
sparse_vector<Val, BV>::sparse_vector(
        bm::null_support null_able,
        allocation_policy_type  ap,
        size_type               bv_max_size,
        const allocator_type&   alloc)
: bv_size_(bv_max_size),
  alloc_(alloc),
  ap_(ap),
  size_(0),
  effective_plains_(0)
{
    ::memset(plains_, 0, sizeof(plains_));
    if (null_able == bm::use_null)
    {
        unsigned i = null_plain();
        plains_[i] = construct_bvector(0);
        plains_[i]->init();
    }
}

//---------------------------------------------------------------------

template<class Val, class BV>
sparse_vector<Val, BV>::sparse_vector(const sparse_vector<Val, BV>& sv)
: bv_size_ (sv.bv_size_),
  alloc_(sv.alloc_),
  ap_(sv.ap_),
  size_(sv.size_),
  effective_plains_(sv.effective_plains_)
{
    if (this != &sv)
    {
        for (size_type i = 0; i < stored_plains(); ++i)
        {
            const bvector_type* bv = sv.plains_[i];
            plains_[i] = bv ? construct_bvector(bv) : 0;
        }
    }
}

//---------------------------------------------------------------------
#ifndef BM_NO_CXX11

template<class Val, class BV>
sparse_vector<Val, BV>::sparse_vector(sparse_vector<Val, BV>&& sv) BMNOEXEPT
{
    bv_size_ = 0;
    alloc_ = sv.alloc_;
    ap_ = sv.ap_;
    size_ = sv.size_;
    effective_plains_ = sv.effective_plains_;
        
    for (size_type i = 0; i < stored_plains(); ++i)
    {
        plains_[i] = sv.plains_[i];
        sv.plains_[i] = 0;
    }
    sv.size_ = 0;
}

#endif


//---------------------------------------------------------------------

template<class Val, class BV>
sparse_vector<Val, BV>::~sparse_vector() BMNOEXEPT
{
    free_vectors();
}

//---------------------------------------------------------------------

template<class Val, class BV>
void sparse_vector<Val, BV>::swap(sparse_vector<Val, BV>& sv) BMNOEXEPT
{
    if (this != &sv)
    {
        bm::xor_swap(bv_size_, sv.bv_size_);
        
        allocator_type alloc_tmp = alloc_;
        alloc_ = sv.alloc_;
        sv.alloc_ = alloc_tmp;
        
        allocation_policy_type ap_tmp = ap_;
        ap_ = sv.ap_;
        sv.ap_ = ap_tmp;
        
        for (size_type i = 0; i < sv_plains; ++i)
        {
            if (plains_[i] != sv.plains_[i])
            {
                bvector_type* bv_tmp = plains_[i];
                plains_[i] = sv.plains_[i];
                sv.plains_[i] = bv_tmp;
            }
        } // for i
        
        bm::xor_swap(size_, sv.size_);
        bm::xor_swap(effective_plains_, sv.effective_plains_);
    }
}

//---------------------------------------------------------------------

template<class Val, class BV>
void sparse_vector<Val, BV>::throw_range_error(const char* err_msg)
{
#ifndef BM_NO_STL
    throw std::range_error(err_msg);
#else
    BM_ASSERT_THROW(false, BM_ERR_RANGE);
#endif
}

//---------------------------------------------------------------------

template<class Val, class BV>
void sparse_vector<Val, BV>::throw_bad_alloc()
{
#ifndef BM_NO_STL
    throw std::bad_alloc();
#else
    BM_ASSERT_THROW(false, BM_BAD_ALLOC);
#endif
}


//---------------------------------------------------------------------

template<class Val, class BV>
typename sparse_vector<Val, BV>::bvector_type*
sparse_vector<Val, BV>::construct_bvector(const bvector_type* bv) const
{
    bvector_type* rbv = 0;
#ifdef BM_NO_STL   // C compatibility mode
    void* mem = ::malloc(sizeof(bvector_type));
    if (mem == 0)
    {
        BM_ASSERT_THROW(false, BM_ERR_BADALLOC);
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

template<class Val, class BV>
bool sparse_vector<Val, BV>::is_nullable() const
{
    return (plains_[this->null_plain()] != 0);
}

//---------------------------------------------------------------------

template<class Val, class BV>
const typename sparse_vector<Val, BV>::bvector_type*
sparse_vector<Val, BV>::get_null_bvector() const
{
    return plains_[this->null_plain()];
}

//---------------------------------------------------------------------

template<class Val, class BV>
bool sparse_vector<Val, BV>::is_null(size_type idx) const
{
    const bvector_type* bv_null = get_null_bvector();
    return (bv_null) ? (!bv_null->test(idx)) : false;
}

//---------------------------------------------------------------------

template<class Val, class BV>
void sparse_vector<Val, BV>::set_null(size_type idx)
{
    clear(idx, true);
}

//---------------------------------------------------------------------

template<class Val, class BV>
void sparse_vector<Val, BV>::import(const value_type* arr,
                                    size_type         size,
                                    size_type         offset)
{
    unsigned char b_list[sizeof(Val)*8];
    unsigned row_len[sizeof(Val)*8] = {0, };
    
    const unsigned transpose_window = 256;
    bm::tmatrix<bm::id_t, sizeof(Val)*8, transpose_window> tm; // matrix accumulator
    
    if (size == 0)
        throw_range_error("sparse_vector range error (import size 0)");
    
    // clear all plains in the range to provide corrrect import of 0 values
    this->clear_range(offset, offset + size - 1);
    
    // transposition algorithm uses bitscan to find index bits and store it
    // in temporary matrix (list for each bit plain), matrix here works
    // when array gets to big - the list gets loaded into bit-vector using
    // bulk load algorithm, which is faster than single bit access
    //
    
    size_type i;
    for (i = 0; i < size; ++i)
    {
        unsigned bcnt = bm::bitscan(arr[i], b_list);
        const unsigned bit_idx = i + offset;
        
        for (unsigned j = 0; j < bcnt; ++j)
        {
            unsigned p = b_list[j];
            unsigned rl = row_len[p];
            tm.row(p)[rl] = bit_idx;
            row_len[p] = ++rl;
            
            if (rl == transpose_window)
            {
                bvector_type* bv = get_plain(p);
                const bm::id_t* r = tm.row(p);
                bm::combine_or(*bv, r, r + rl);
                row_len[p] = 0;
                tm.row(p)[0] = 0;
            }
        } // for j
    } // for i
    
    // process incomplete transposition lines
    //
    for (unsigned k = 0; k < tm.rows(); ++k)
    {
        unsigned rl = row_len[k];
        if (rl)
        {
            bvector_type* bv = get_plain(k);
            const bm::id_t* r = tm.row(k);
            bm::combine_or(*bv, r, r + rl);
        }
    } // for k
    
    
    if (i + offset > size_)
        size_ = i + offset;
    
    bvector_type* bv_null = get_null_bvect();
    if (bv_null) // configured to support NULL assignments
    {
        bv_null->set_range(offset, offset + size - 1);
    }
}

//---------------------------------------------------------------------

template<class Val, class BV>
void sparse_vector<Val, BV>::import_back(const value_type* arr,
                                         size_type         size)
{
    this->import(arr, size, this->size());
}

//---------------------------------------------------------------------

template<class Val, class BV>
typename sparse_vector<Val, BV>::size_type
sparse_vector<Val, BV>::decode(value_type* arr,
                               size_type   idx_from,
                               size_type   size,
                               bool        zero_mem) const
{
    if (size < 32)
    {
        return extract_range(arr, size, idx_from, zero_mem);
    }
    if (size < 1024)
    {
        return extract_plains(arr, size, idx_from, zero_mem);
    }
    return extract(arr, size, idx_from, zero_mem);
}

//---------------------------------------------------------------------

template<class Val, class BV>
typename sparse_vector<Val, BV>::size_type
sparse_vector<Val, BV>::gather(value_type*       arr,
                               const size_type*  idx,
                               size_type         size,
                               bm::sort_order    sorted_idx) const
{
    BM_ASSERT(arr);
    BM_ASSERT(idx);
    BM_ASSERT(size);

    if (size == 1) // corner case: get 1 value
    {
        arr[0] = this->get(idx[0]);
        return size;
    }
    ::memset(arr, 0, sizeof(value_type)*size);
    
    for (unsigned i = 0; i < size;)
    {
        bool sorted_block = true;
        
        // look ahead for the depth of the same block
        //          (speculate more than one index lookup per block)
        //
        unsigned nb = unsigned(idx[i] >> bm::set_block_shift);
        unsigned r = i;
        
        switch (sorted_idx)
        {
        case BM_UNKNOWN:
            {
                size_type idx_prev = idx[r];
                for (; (r < size) && (nb == unsigned(idx[r] >> bm::set_block_shift)); ++r)
                {
                    sorted_block = !(idx[r] < idx_prev); // sorted check
                    idx_prev = idx[r];
                }
            }
            break;
        case BM_UNSORTED:
            sorted_block = false;
            
            for (; r < size; ++r)
                if (nb != unsigned(idx[r] >> bm::set_block_shift))
                    break;
            break;            
            // no break(!) intentional fall through
        case BM_SORTED:
            r = bm::idx_arr_block_lookup(idx, size, nb, r);
            break;
        case BM_SORTED_UNIFORM:
            r = size;
            break;
        default:
            BM_ASSERT(0);
        } // switch
        
        // single element hit, use plain random access
        if (r == i+1)
        {
            arr[i] = this->get(idx[i]);
            ++i;
            continue;
        }

        // process block co-located elements at ones for best (CPU cache opt)
        //
        unsigned i0 = nb >> bm::set_array_shift; // top block address
        unsigned j0 = nb &  bm::set_array_mask;  // address in sub-block
        
        unsigned eff_plains = effective_plains();
        for (unsigned j = 0; j < eff_plains; ++j)
        {
            const bm::word_t* blk = get_block(j, i0, j0);
            if (!blk)
                continue;
            value_type vm;
            if (blk == FULL_BLOCK_FAKE_ADDR)
            {
                vm = (1u << j);
                for (unsigned k = i; k < r; ++k)
                    arr[k] |= vm;
                continue;
            }
            if (BM_IS_GAP(blk))
            {
                const bm::gap_word_t* gap_blk = BMGAP_PTR(blk);
                unsigned is_set;
                
                if (sorted_block) // b-search hybrid with scan lookup
                {
                    for (unsigned k = i; k < r; )
                    {
                        unsigned nbit = unsigned(idx[k] & bm::set_block_mask);
                        
                        unsigned gidx = bm::gap_bfind(gap_blk, nbit, &is_set);
                        unsigned gap_value = gap_blk[gidx];
                        if (is_set)
                        {
                            arr[k] |= vm = (1u << j);
                            for (++k; k < r; ++k) // speculative look-up
                            {
                                if (unsigned(idx[k] & bm::set_block_mask) <= gap_value)
                                    arr[k] |= vm;
                                else
                                    break;
                            }
                        }
                        else // 0 GAP - skip. not set
                        {
                            for (++k;
                                 (k < r) &&
                                 (unsigned(idx[k] & bm::set_block_mask) <= gap_value);
                                 ++k) {}
                        }
                    } // for k
                }
                else // unsorted block gather request: b-search lookup
                {
                    for (unsigned k = i; k < r; ++k)
                    {
                        unsigned nbit = unsigned(idx[k] & bm::set_block_mask);
                        is_set = bm::gap_test_unr(gap_blk, nbit);
                        arr[k] |= value_type(bool(is_set) << j);
                    } // for k
                }
                continue;
            }
            bm::bit_block_gather_scatter(arr, blk, idx, r, i, j);
        } // for (each plain)
        
        i = r;

    } // for i

    return size;
}

//---------------------------------------------------------------------

template<class Val, class BV>
typename sparse_vector<Val, BV>::size_type
sparse_vector<Val, BV>::extract_range(value_type* arr,
                                      size_type size,
                                      size_type offset,
                                      bool      zero_mem) const
{
    if (size == 0)
        return 0;
    if (zero_mem)
        ::memset(arr, 0, sizeof(value_type)*size);

    size_type start = offset;
    size_type end = start + size;
    if (end > size_)
    {
        end = size_;
    }
    
    // calculate logical block coordinates and masks
    //
    unsigned nb = unsigned(start >>  bm::set_block_shift);
    unsigned i0 = nb >> bm::set_array_shift; // top block address
    unsigned j0 = nb &  bm::set_array_mask;  // address in sub-block
    unsigned nbit = unsigned(start & bm::set_block_mask);
    unsigned nword  = unsigned(nbit >> bm::set_word_shift);
    unsigned mask0 = 1u << (nbit & bm::set_word_mask);
    const bm::word_t* blk = 0;
    unsigned is_set;
    
    for (unsigned j = 0; j < sizeof(Val)*8; ++j)
    {
        blk = get_block(j, i0, j0);
        bool is_gap = BM_IS_GAP(blk);
        
        for (unsigned k = start; k < end; ++k)
        {
            unsigned nb1 = unsigned(k >>  bm::set_block_shift);
            if (nb1 != nb) // block switch boundaries
            {
                nb = nb1;
                i0 = nb >> bm::set_array_shift;
                j0 = nb &  bm::set_array_mask;
                blk = get_block(j, i0, j0);
                is_gap = BM_IS_GAP(blk);
            }
        
            if (!blk)
                continue;
            nbit = unsigned(k & bm::set_block_mask);
            if (is_gap)
            {
                is_set = bm::gap_test_unr(BMGAP_PTR(blk), nbit);
            }
            else
            {
                if (blk == FULL_BLOCK_FAKE_ADDR)
                {
                    is_set = 1;
                }
                else
                {
                    nword  = unsigned(nbit >> bm::set_word_shift);
                    mask0 = 1u << (nbit & bm::set_word_mask);
                    is_set = (blk[nword] & mask0);
                }
            }
            size_type idx = k - offset;
            value_type vm = (bool) is_set;
            vm <<= j;
            arr[idx] |= vm;
            
        } // for k

    } // for j
    return 0;
}

//---------------------------------------------------------------------

template<class Val, class BV>
typename sparse_vector<Val, BV>::size_type
sparse_vector<Val, BV>::extract_plains(value_type* arr,
                                       size_type   size,
                                       size_type   offset,
                                       bool        zero_mem) const
{
    if (size == 0)
        return 0;

    if (zero_mem)
        ::memset(arr, 0, sizeof(value_type)*size);
    
    size_type start = offset;
    size_type end = start + size;
    if (end > size_)
    {
        end = size_;
    }
    
    for (size_type i = 0; i < value_bits(); ++i)
    {
        const bvector_type* bv = plains_[i];
        if (!bv)
            continue;
       
        value_type mask = 1;
        mask <<= i;
        typename BV::enumerator en(bv, offset);
        for (;en.valid(); ++en)
        {
            size_type idx = *en - offset;
            if (idx >= size)
                break;
            arr[idx] |= mask;
        } // for
        
    } // for i

    return 0;
}


//---------------------------------------------------------------------

template<class Val, class BV>
typename sparse_vector<Val, BV>::size_type
sparse_vector<Val, BV>::extract(value_type* arr,
                                size_type   size,
                                size_type   offset,
                                bool        zero_mem,
                                allocator_pool_type* pool_ptr) const
{
    /// Decoder functor
    /// @internal
    ///
    struct sv_decode_visitor_func
    {
        sv_decode_visitor_func(value_type* varr,
                               value_type  mask,
                               size_type   off)
        : arr_(varr), mask_(mask), off_(off)
        {}
        
        void add_bits(bm::id_t arr_offset, const unsigned char* bits, unsigned bits_size)
        {
            size_type idx_base = arr_offset - off_;
            
            const value_type m = mask_;
            unsigned i = 0;
            for (; i < bits_size; ++i)
            {
                arr_[idx_base + bits[i]] |= m;
            }
        }
        
        void add_range(bm::id_t arr_offset, unsigned sz)
        {
            size_type idx_base = arr_offset - off_;
            const value_type m = mask_;
            for (unsigned i = 0;i < sz; ++i)
                arr_[i + idx_base] |= m;
        }
        value_type*  arr_;
        value_type   mask_;
        size_type    off_;
    };


    if (size == 0)
        return 0;

    if (zero_mem)
        ::memset(arr, 0, sizeof(value_type)*size);
    
    size_type start = offset;
    size_type end = start + size;
    if (end > size_)
    {
        end = size_;
    }
    
	bool masked_scan = !(offset == 0 && size == this->size());

    if (masked_scan) // use temp vector to decompress the area
    {
        // for large array extraction use logical opertions
        // (faster due to vectorization)
        bvector_type bv_mask;
        bv_mask.set_allocator_pool(pool_ptr);
        
        for (size_type i = 0; i < value_bits(); ++i)
        {
            const bvector_type* bv = plains_[i];
            if (bv)
            {
                bv_mask.set_range(offset, end - 1);
                bv_mask.bit_and(*bv); // TODO: make copy-and fusion
            
                sv_decode_visitor_func func(arr, (value_type(1) << i), offset);
                bm::for_each_bit(bv_mask, func);
                bv_mask.set_range(offset, end - 1, false);
            }
        } // for i
    }
    else
    {
        for (size_type i = 0; i < value_bits(); ++i)
        {
            const bvector_type* bv = plains_[i];
            if (bv)
            {
                sv_decode_visitor_func func(arr, (value_type(1) << i), 0);
                bm::for_each_bit(*bv, func);
            }
        } // for i
    }

    return end - start;
}

//---------------------------------------------------------------------

template<class Val, class BV>
typename sparse_vector<Val, BV>::size_type
sparse_vector<Val, BV>::size() const
{
    return size_;
}

//---------------------------------------------------------------------

template<class Val, class BV>
bool sparse_vector<Val, BV>::empty() const
{
    return (size_ == 0);
}


//---------------------------------------------------------------------

template<class Val, class BV>
void sparse_vector<Val, BV>::resize(typename sparse_vector<Val, BV>::size_type sz)
{
    if (sz == size_)  // nothing to do
        return;
    
    if (!sz) // resize to zero is an equivalent of non-destructive deallocation
    {
        clear();
        return;
    }
    
    if (sz < size_) // vector shrink
    {
        this->clear_range(sz, size_-1, true);   // clear the tails and NULL vect
    }
    
    size_ = sz;
}

//---------------------------------------------------------------------


template<class Val, class BV>
typename sparse_vector<Val, BV>::bvector_type_ptr
   sparse_vector<Val, BV>::get_plain(unsigned i)
{
    bvector_type_ptr bv = plains_[i];
    if (!bv)
    {
        bv = construct_bvector(0);
        bv->init();
        plains_[i] = bv;
        if (i > effective_plains_ && i < value_bits())
            effective_plains_ = i;
    }
    return bv;
}

//---------------------------------------------------------------------

template<class Val, class BV>
typename sparse_vector<Val, BV>::value_type
sparse_vector<Val, BV>::at(typename sparse_vector<Val, BV>::size_type idx) const
{
    if (idx >= size_)
        throw_range_error("sparse vector range error");
    return this->get(idx);
}

//---------------------------------------------------------------------

template<class Val, class BV>
const bm::word_t* sparse_vector<Val, BV>::get_block(unsigned p, unsigned i, unsigned j) const
{
    const bvector_type* bv = this->plains_[p];
    if (bv)
    {
        const typename bvector_type::blocks_manager_type& bman = bv->get_blocks_manager();
        return bman.get_block_ptr(i, j);
    }
    return 0;
}

//---------------------------------------------------------------------

template<class Val, class BV>
typename sparse_vector<Val, BV>::value_type
sparse_vector<Val, BV>::get(bm::id_t i) const
{
    BM_ASSERT(i < size_);
    
    value_type v = 0;
    
    // calculate logical block coordinates and masks
    //
    unsigned nb = unsigned(i >>  bm::set_block_shift);
    unsigned i0 = nb >> bm::set_array_shift; // top block address
    unsigned j0 = nb &  bm::set_array_mask;  // address in sub-block
    unsigned nbit = unsigned(i & bm::set_block_mask);
    unsigned nword  = unsigned(nbit >> bm::set_word_shift);
    unsigned mask0 = 1u << (nbit & bm::set_word_mask);
    const bm::word_t* blk;
    const bm::word_t* blka[4];
    
    unsigned eff_plains = effective_plains();
    for (unsigned j = 0; j < eff_plains; j+=4)
    {
#if defined(BM64_SSE4)
        __m128i w0 = _mm_loadu_si128((__m128i*)(plains_ + j));
        __m128i w1 = _mm_loadu_si128((__m128i*)(plains_ + j + 2));
        w0 = _mm_or_si128(w0, w1);
        if (_mm_testz_si128(w0, w0))
        {
            BM_ASSERT(!(plains_[j + 0] || plains_[j + 1] || plains_[j + 2] || plains_[j + 3]));
            continue;
        }
#elif defined(BM64_AVX2)
        __m256i w0 = _mm256_loadu_si256((__m256i*)(plains_ + j));
        if (_mm256_testz_si256(w0, w0))
        {
            BM_ASSERT(!(plains_[j + 0] || plains_[j + 1] || plains_[j + 2] || plains_[j + 3]));
            continue;
        }
#else
        bool b = plains_[j + 0] || plains_[j + 1] || plains_[j + 2] || plains_[j + 3];
        if (!b)
            continue;
#endif
        blka[0] = get_block(j+0, i0, j0);
        blka[1] = get_block(j+1, i0, j0);
        blka[2] = get_block(j+2, i0, j0);
        blka[3] = get_block(j+3, i0, j0);
        unsigned is_set;

        if ((blk = blka[0+0])!=0)
        {
            if (blk == FULL_BLOCK_FAKE_ADDR)
                is_set = 1;
            else
                is_set = (BM_IS_GAP(blk)) ? bm::gap_test_unr(BMGAP_PTR(blk), nbit) : (blk[nword] & mask0);
            value_type vm = (bool) is_set;
            vm <<= (j+0);
            v |= vm;
        }
        if ((blk = blka[0+1])!=0)
        {
            if (blk == FULL_BLOCK_FAKE_ADDR)
                is_set = 1;
            else
                is_set = (BM_IS_GAP(blk)) ? bm::gap_test_unr(BMGAP_PTR(blk), nbit) : (blk[nword] & mask0);
            value_type vm = (bool) is_set;
            vm <<= (j+1);
            v |= vm;
        }
        if ((blk = blka[0+2])!=0)
        {
            if (blk == FULL_BLOCK_FAKE_ADDR)
                is_set = 1;
            else
                is_set = (BM_IS_GAP(blk)) ? bm::gap_test_unr(BMGAP_PTR(blk), nbit) : (blk[nword] & mask0);
            value_type vm = (bool) is_set;
            vm <<= (j+2);
            v |= vm;
        }
        if ((blk = blka[0+3])!=0)
        {
            if (blk == FULL_BLOCK_FAKE_ADDR)
                is_set = 1;
            else
                is_set = (BM_IS_GAP(blk)) ? bm::gap_test_unr(BMGAP_PTR(blk), nbit) : (blk[nword] & mask0);
            value_type vm = (bool) is_set;
            vm <<= (j+3);
            v |= vm;
        }

    } // for j
    
    return v;
}


//---------------------------------------------------------------------

template<class Val, class BV>
void sparse_vector<Val, BV>::set(size_type idx, value_type v)
{ 
    if (idx >= size_)
    {
        size_ = idx+1;
    }
    set_value(idx, v);
}

//---------------------------------------------------------------------

template<class Val, class BV>
void sparse_vector<Val, BV>::clear(size_type idx, bool set_null)
{
    if (idx >= size_)
        size_ = idx+1;

    set_value(idx, (value_type)0);
    if (set_null)
    {
        bvector_type* bv_null = get_null_bvect();
        if (bv_null)
            bv_null->set(idx, false);
    }
}

//---------------------------------------------------------------------

template<class Val, class BV>
void sparse_vector<Val, BV>::push_back(value_type v)
{
    set_value(size_, v);
    ++size_;
}

//---------------------------------------------------------------------

template<class Val, class BV>
void sparse_vector<Val, BV>::push_back_no_null(value_type v)
{
    set_value_no_null(size_, v);
    ++size_;
}

//---------------------------------------------------------------------

template<class Val, class BV>
void sparse_vector<Val, BV>::set_value(size_type idx, value_type v)
{
    set_value_no_null(idx, v);
    bvector_type* bv_null = get_null_bvect();
    if (bv_null)
        bv_null->set_bit_no_check(idx);
}

//---------------------------------------------------------------------

template<class Val, class BV>
void sparse_vector<Val, BV>::set_value_no_null(size_type idx, value_type v)
{
    // calculate logical block coordinates and masks
    //
    unsigned nb = unsigned(idx >>  bm::set_block_shift);
    unsigned i0 = nb >> bm::set_array_shift; // top block address
    unsigned j0 = nb &  bm::set_array_mask;  // address in sub-block

    // clear the plains where needed
    unsigned eff_plains = effective_plains();
    unsigned bsr = v ? bm::bit_scan_reverse(v) : 0u;
        
    for (unsigned i = bsr; i < eff_plains; ++i)
    {
        const bm::word_t* blk = get_block(i, i0, j0);
        if (blk)
        {
            BM_ASSERT(plains_[i]);
            bvector_type* bv = plains_[i];
            bv->clear_bit_no_check(idx);
        }
    }
    if (v)
    {
        value_type mask = 1u;
        for (unsigned j = 0; j <= bsr; ++j)
        {
            if (v & mask)
            {
                bvector_type* bv = get_plain(j);
                bv->set_bit_no_check(idx);
            }
            else
            {
                const bm::word_t* blk = get_block(j, i0, j0);
                if (blk)
                {
                    BM_ASSERT(plains_[j]);
                    bvector_type* bv = plains_[j];
                    bv->clear_bit_no_check(idx);
                }
            }
            mask <<=  1;
        }
    }

}

//---------------------------------------------------------------------

template<class Val, class BV>
void sparse_vector<Val, BV>::inc(size_type idx)
{
    if (idx >= size_)
        size_ = idx+1;

    for (unsigned i = 0; i < sv_value_plains; ++i)
    {
        bvector_type* bv = get_plain(i);
        bool carry_over = bv->inc(idx);
        if (!carry_over)
            break;
    }
    bvector_type* bv_null = get_null_bvect();
    if (bv_null)
        bv_null->set_bit_no_check(idx);
}

//---------------------------------------------------------------------

template<class Val, class BV>
void sparse_vector<Val, BV>::clear() BMNOEXEPT
{
    for (size_type i = 0; i < value_bits(); ++i)
    {
        bvector_type* bv = plains_[i];
        if (bv)
        {
            destruct_bvector(bv);
            plains_[i] = 0;
        }
    }
    size_ = 0;
    bvector_type* bv_null = get_null_bvect();
    if (bv_null)
    {
        bv_null->clear(true);
        bv_null->init();
    }
}

//---------------------------------------------------------------------

template<class Val, class BV>
void sparse_vector<Val, BV>::destruct_bvector(bvector_type* bv) const
{
#ifdef BM_NO_STL   // C compatibility mode
    bv->~TBM_bvector();
    ::free((void*)bv);
#else
    delete bv;
#endif
}

//---------------------------------------------------------------------

template<class Val, class BV>
void sparse_vector<Val, BV>::free_vectors() BMNOEXEPT
{
    for (size_type i = 0; i < stored_plains(); ++i)
    {
        bvector_type* bv = plains_[i];
        destruct_bvector(bv);
    }
}

//---------------------------------------------------------------------

template<class Val, class BV>
bool sparse_vector<Val, BV>::find_rank(bm::id_t rank, bm::id_t& pos) const
{
    BM_ASSERT(rank);
    pos = rank - 1; 
    return true;
}

//---------------------------------------------------------------------


template<class Val, class BV>
void sparse_vector<Val, BV>::free_plain(unsigned i)
{
    BM_ASSERT(i < stored_plains());
    bvector_type* bv = plains_[i];
    destruct_bvector(bv);
    plains_[i] = 0;
}

//---------------------------------------------------------------------

template<class Val, class BV>
sparse_vector<Val, BV>&
sparse_vector<Val, BV>::clear_range(
    typename sparse_vector<Val, BV>::size_type left,
    typename sparse_vector<Val, BV>::size_type right,
    bool set_null)
{
    if (right < left)
    {
        return clear_range(right, left);
    }
    unsigned eff_plains = effective_plains();
    for (unsigned i = 0; i < eff_plains; ++i)
    {
        bvector_type* bv = plains_[i];
        if (bv)
        {
            bv->set_range(left, right, false);
        }
    } // for i
    
    if (set_null)
    {
        bvector_type* bv_null = get_null_bvect();
        if (bv_null)
            bv_null->set_range(left, right, false);
    }
    
    return *this;
}

//---------------------------------------------------------------------

template<class Val, class BV>
void sparse_vector<Val, BV>::calc_stat(
     struct sparse_vector<Val, BV>::statistics* st) const
{
    BM_ASSERT(st);
    
    st->reset();

    unsigned stored_plains = this->stored_plains();
    for (unsigned j = 0; j < stored_plains; ++j)
    {
        const bvector_type* bv = this->plains_[j];
        if (bv)
        {
            typename bvector_type::statistics stbv;
            bv->calc_stat(&stbv);
            
            st->bit_blocks += stbv.bit_blocks;
            st->gap_blocks += stbv.gap_blocks;
            st->max_serialize_mem += stbv.max_serialize_mem + 8;
            st->memory_used += stbv.memory_used;
        }
    } // for j
    // header accounting
    st->max_serialize_mem += 1 + 1 + 1 + 1 + 8 + (8 * this->stored_plains());

}

//---------------------------------------------------------------------

template<class Val, class BV>
void sparse_vector<Val, BV>::optimize(
    bm::word_t*                                  temp_block, 
    typename bvector_type::optmode               opt_mode,
    typename sparse_vector<Val, BV>::statistics* st)
{
    if (st)
    {
        st->bit_blocks = st->gap_blocks = 0;
        st->max_serialize_mem = st->memory_used = 0;
    }
    bvector_type* bv_null = this->get_null_bvect();
    
    unsigned stored_plains = this->stored_plains();
    for (unsigned j = 0; j < stored_plains; ++j)
    {
        bvector_type* bv = this->plains_[j];
        if (bv)
        {
            if (bv != bv_null) // protect the NULL vector from de-allocation
            {
                if (!bv->any())  // empty vector?
                {
                    destruct_bvector(this->plains_[j]);
                    this->plains_[j] = 0;
                    continue;
                }
            }
            
            typename bvector_type::statistics stbv;
            bv->optimize(temp_block, opt_mode, &stbv);
            
            if (st)
            {
                st->bit_blocks += stbv.bit_blocks;
                st->gap_blocks += stbv.gap_blocks;
                st->max_serialize_mem += stbv.max_serialize_mem + 8;
                st->memory_used += stbv.memory_used;
            }
        }
    } // for j

}

//---------------------------------------------------------------------


template<class Val, class BV>
void sparse_vector<Val, BV>::optimize_gap_size()
{
    unsigned stored_plains = this->stored_plains();
    for (unsigned j = 0; j < stored_plains; ++j)
    {
        bvector_type* bv = this->plains_[j];
        if (bv)
        {
            bv->optimize_gap_size();
        }
    }
}

//---------------------------------------------------------------------

template<class Val, class BV>
sparse_vector<Val, BV>&
sparse_vector<Val, BV>::join(const sparse_vector<Val, BV>& sv)
{
    size_type arg_size = sv.size();
    if (size_ < arg_size)
    {
        resize(arg_size);
    }
    bvector_type* bv_null = this->get_null_bvect();
    unsigned plains;
    if (bv_null)
        plains = this->stored_plains();
    else
        plains = this->plains();
    
    for (unsigned j = 0; j < plains; ++j)
    {
        const bvector_type* arg_bv = sv.plains_[j];
        if (arg_bv)
        {
            bvector_type* bv = this->plains_[j];
            if (!bv)
            {
                bv = get_plain(j);
            }
            *bv |= *arg_bv;
        }
    } // for j
    
    // our vector is NULL-able but argument is not (assumed all values are real)
    if (bv_null && !sv.is_nullable())
    {
        bv_null->set_range(0, arg_size-1);
    }
    
    return *this;
}

//---------------------------------------------------------------------

template<class Val, class BV>
void sparse_vector<Val, BV>::copy_range(const sparse_vector<Val, BV>& sv,
                                        typename sparse_vector<Val, BV>::size_type left,
                                        typename sparse_vector<Val, BV>::size_type right)
{
    if (left > right)
        bm::xor_swap(left, right);
    
    bvector_type* bv_null = this->get_null_bvect();
    unsigned plains;
    if (bv_null)
    {
        plains = this->stored_plains();
        const bvector_type* bv_null_arg = sv.get_null_bvector();
        if (bv_null_arg)
            bv_null->copy_range(*bv_null_arg, left, right);
    }
    else
        plains = this->plains();
    
    for (unsigned j = 0; j < plains; ++j)
    {
        const bvector_type* arg_bv = sv.plains_[j];
        if (arg_bv)
        {
            bvector_type* bv = this->plains_[j];
            if (!bv)
                bv = get_plain(j);
            bv->copy_range(*arg_bv, left, right);
        }
    } // for j
    this->resize(sv.size());
}
//---------------------------------------------------------------------

template<class Val, class BV>
void sparse_vector<Val, BV>::filter(
                    const typename sparse_vector<Val, BV>::bvector_type& bv_mask)
{
    bvector_type* bv_null = this->get_null_bvect();
    unsigned plains;
    if (bv_null)
    {
        plains = this->stored_plains();
        bv_null->bit_and(bv_mask);
    }
    else
        plains = this->plains();
    
    for (unsigned j = 0; j < plains; ++j)
    {
        bvector_type* bv = this->plains_[j];
        if (bv)
            bv->bit_and(bv_mask);
    }
}


//---------------------------------------------------------------------

template<class Val, class BV>
bool sparse_vector<Val, BV>::equal(const sparse_vector<Val, BV>& sv,
                                   bm::null_support null_able) const
{
    size_type arg_size = sv.size();
    if (size_ != arg_size)
    {
        return false;
    }
    unsigned plains = this->plains();
    for (unsigned j = 0; j < plains; ++j)
    {
        const bvector_type* bv = plains_[j];
        const bvector_type* arg_bv = sv.plains_[j];
        if (bv == arg_bv) // same NULL
            continue;
        // check if any not NULL and not empty
        if (!bv && arg_bv)
            return !arg_bv->any();
        if (bv && !arg_bv)
            return !bv->any();
        // both not NULL
        int cmp = bv->compare(*arg_bv);
        if (cmp != 0)
            return false;
    } // for j
    
    if (null_able == bm::use_null)
    {
        const bvector_type* bv_null = this->get_null_bvector();
        const bvector_type* bv_null_arg = sv.get_null_bvector();
        
        // check the NULL vectors
        if (bv_null == bv_null_arg)
            return true;
        if (!bv_null || !bv_null_arg)
            return false;
        BM_ASSERT(bv_null);
        BM_ASSERT(bv_null_arg);
        int cmp = bv_null->compare(*bv_null);
        if (cmp != 0)
            return false;
    }

    return true;
}

//---------------------------------------------------------------------

template<class Val, class BV>
typename sparse_vector<Val, BV>::const_iterator
sparse_vector<Val, BV>::begin() const
{
    typedef typename sparse_vector<Val, BV>::const_iterator it_type;
    return it_type(this);
}

//---------------------------------------------------------------------

template<class Val, class BV>
void sparse_vector<Val, BV>::set_allocator_pool(
    typename sparse_vector<Val, BV>::allocator_pool_type* pool_ptr)
{
    pool_ = pool_ptr;
}


//---------------------------------------------------------------------
//
//---------------------------------------------------------------------


template<class Val, class BV>
sparse_vector<Val, BV>::const_iterator::const_iterator()
: sv_(0), pos_(bm::id_max), buf_ptr_(0)
{}

//---------------------------------------------------------------------

template<class Val, class BV>
sparse_vector<Val, BV>::const_iterator::const_iterator(
                        const typename sparse_vector<Val, BV>::const_iterator& it)
: sv_(it.sv_), pos_(it.pos_), buf_ptr_(0)
{}

//---------------------------------------------------------------------

template<class Val, class BV>
sparse_vector<Val, BV>::const_iterator::const_iterator(
  const typename sparse_vector<Val, BV>::const_iterator::sparse_vector_type* sv)
: sv_(sv), buf_ptr_(0)
{
    BM_ASSERT(sv_);
    pos_ = sv_->empty() ? bm::id_max : 0u;
}

//---------------------------------------------------------------------

template<class Val, class BV>
sparse_vector<Val, BV>::const_iterator::const_iterator(
 const typename sparse_vector<Val, BV>::const_iterator::sparse_vector_type* sv,
 bm::id_t pos)
: sv_(sv), buf_ptr_(0)
{
    BM_ASSERT(sv_);
    this->go_to(pos);
}

//---------------------------------------------------------------------

template<class Val, class BV>
void sparse_vector<Val, BV>::const_iterator::go_to(bm::id_t pos)
{
    pos_ = (!sv_ || pos >= sv_->size()) ? bm::id_max : pos;
    buf_ptr_ = 0;
}

//---------------------------------------------------------------------

template<class Val, class BV>
void sparse_vector<Val, BV>::const_iterator::advance()
{
    if (pos_ == bm::id_max) // nothing to do, we are at the end
        return;
    ++pos_;
    if (pos_ >= sv_->size())
        pos_ = bm::id_max;
    else
    {
        if (buf_ptr_)
        {
            ++buf_ptr_;
            if (buf_ptr_ - ((value_type*)buffer_.data()) >= n_buf_size)
                buf_ptr_ = 0;
        }
    }
}

//---------------------------------------------------------------------

template<class Val, class BV>
typename sparse_vector<Val, BV>::const_iterator::value_type
sparse_vector<Val, BV>::const_iterator::value() const
{
    BM_ASSERT(this->valid());
    value_type v;
    
    if (!buf_ptr_)
    {
        buffer_.reserve(n_buf_size * sizeof(value_type));
        buf_ptr_ = (value_type*)(buffer_.data());
        sv_->extract(buf_ptr_, n_buf_size, pos_, true, &pool_);
    }
    v = *buf_ptr_;
    return v;
}

//---------------------------------------------------------------------

template<class Val, class BV>
void sparse_vector<Val, BV>::const_iterator::skip_zero_values()
{
    value_type v = value();
    if (buf_ptr_)
    {
        v = *buf_ptr_;
        value_type* buf_end = ((value_type*)buffer_.data()) + n_buf_size;
        while(!v)
        {
            ++pos_;
            if (++buf_ptr_ < buf_end)
                v = *buf_ptr_;
            else
                break;
        }
        if (pos_ >= sv_->size())
        {
            pos_ = bm::id_max;
            return;
        }
        if (buf_ptr_ >= buf_end)
            buf_ptr_ = 0;
    }
}

//---------------------------------------------------------------------

template<class Val, class BV>
bool sparse_vector<Val, BV>::const_iterator::is_null() const
{
    return sv_->is_null(pos_);
}


//---------------------------------------------------------------------
//
//---------------------------------------------------------------------


template<class Val, class BV>
sparse_vector<Val, BV>::back_insert_iterator::back_insert_iterator()
: sv_(0), bv_null_(0), buf_ptr_(0)
{}

//---------------------------------------------------------------------

template<class Val, class BV>
sparse_vector<Val, BV>::back_insert_iterator::back_insert_iterator(
   typename sparse_vector<Val, BV>::back_insert_iterator::sparse_vector_type* sv)
: sv_(sv), buf_ptr_(0)
{
    bv_null_ = sv_? sv_->get_null_bvect() : 0;
}

//---------------------------------------------------------------------

template<class Val, class BV>
sparse_vector<Val, BV>::back_insert_iterator::back_insert_iterator(
    const typename sparse_vector<Val, BV>::back_insert_iterator::back_insert_iterator& bi)
: sv_(bi.sv_), buf_ptr_(0)
{
    BM_ASSERT(bi.empty());
}

//---------------------------------------------------------------------

template<class Val, class BV>
sparse_vector<Val, BV>::back_insert_iterator::~back_insert_iterator()
{
    this->flush();
}

//---------------------------------------------------------------------

template<class Val, class BV>
void sparse_vector<Val, BV>::back_insert_iterator::add(
         typename sparse_vector<Val, BV>::back_insert_iterator::value_type v)
{
    unsigned buf_idx = this->add_value(v);
    if (bv_null_)
    {
        typename sparse_vector<Val, BV>::size_type sz = sv_->size();
        bv_null_->set_bit_no_check(sz + buf_idx);
    }
}

//---------------------------------------------------------------------

template<class Val, class BV>
unsigned sparse_vector<Val, BV>::back_insert_iterator::add_value(
         typename sparse_vector<Val, BV>::back_insert_iterator::value_type v)
{
    BM_ASSERT(sv_);
    if (!buf_ptr_) // not allocated (yet)
    {
        buffer_.reserve(n_buf_size * sizeof(value_type));
        buf_ptr_ = (value_type*)(buffer_.data());
        *buf_ptr_ = v;
        ++buf_ptr_;
        return 0;
    }
    if (buf_ptr_ - ((value_type*)buffer_.data()) >= n_buf_size)
    {
        this->flush();
        buf_ptr_ = (value_type*)(buffer_.data());
    }
    
    *buf_ptr_ = v;
    size_type buf_idx = size_type(buf_ptr_ - ((value_type*)buffer_.data()));
    ++buf_ptr_;
    return buf_idx;
}

//---------------------------------------------------------------------

template<class Val, class BV>
void sparse_vector<Val, BV>::back_insert_iterator::add_null()
{
    this->add_value(value_type(0));
}

//---------------------------------------------------------------------

template<class Val, class BV>
void sparse_vector<Val, BV>::back_insert_iterator::add_null(
    typename sparse_vector<Val, BV>::back_insert_iterator::size_type count)
{
    for (size_type i = 0; i < count; ++i) // TODO: optimization
        this->add_value(value_type(0));
}

//---------------------------------------------------------------------

template<class Val, class BV>
bool sparse_vector<Val, BV>::back_insert_iterator::empty() const
{
    return (!buf_ptr_ || !sv_);
}

//---------------------------------------------------------------------

template<class Val, class BV>
void sparse_vector<Val, BV>::back_insert_iterator::flush()
{
    if (this->empty())
        return;
    value_type* d = (value_type*)buffer_.data();
    sv_->import_back(d, size_type(buf_ptr_ - d));
    buf_ptr_ = 0;
}

//---------------------------------------------------------------------


} // namespace bm

#include "bmundef.h"


#endif
