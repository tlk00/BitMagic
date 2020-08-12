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

/*! \file bmbmatrix.h
    \brief basic bit-matrix class and utilities
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
    typedef typename bvector_type::size_type             size_type;
    typedef typename bvector_type::block_idx_type        block_idx_type;
    typedef unsigned char                                octet_type;

public:
    // ------------------------------------------------------------
    /*! @name Construction, assignment                          */
    ///@{

    basic_bmatrix(size_type rsize,
                  allocation_policy_type ap = allocation_policy_type(),
                  size_type bv_max_size = bm::id_max,
                  const allocator_type&   alloc  = allocator_type());
    ~basic_bmatrix() BMNOEXCEPT;
    
    /*! copy-ctor */
    basic_bmatrix(const basic_bmatrix<BV>& bbm);
    basic_bmatrix<BV>& operator=(const basic_bmatrix<BV>& bbm)
    {
        copy_from(bbm);
        return *this;
    }
    
#ifndef BM_NO_CXX11
    /*! move-ctor */
    basic_bmatrix(basic_bmatrix<BV>&& bbm) BMNOEXCEPT;

    /*! move assignmment operator */
    basic_bmatrix<BV>& operator = (basic_bmatrix<BV>&& bbm) BMNOEXCEPT
    {
        if (this != &bbm)
        {
            free_rows();
            swap(bbm);
        }
        return *this;
    }
#endif

    void set_allocator_pool(allocator_pool_type* pool_ptr) BMNOEXCEPT
    { pool_ = pool_ptr; }

    ///@}
    
    // ------------------------------------------------------------
    /*! @name content manipulation                                */
    ///@{

    /*! Swap content */
    void swap(basic_bmatrix<BV>& bbm) BMNOEXCEPT;
    
    /*! Copy content */
    void copy_from(const basic_bmatrix<BV>& bbm);
    
    ///@}

    // ------------------------------------------------------------
    /*! @name row access                                         */
    ///@{

    /*! Get row bit-vector. Can return NULL */
    const bvector_type* row(size_type i) const BMNOEXCEPT;

    /*! Get row bit-vector. Can return NULL */
    bvector_type_const_ptr get_row(size_type i) const BMNOEXCEPT;

    /*! Get row bit-vector. Can return NULL */
    bvector_type* get_row(size_type i) BMNOEXCEPT;
    
    /*! get number of value rows */
    size_type rows() const BMNOEXCEPT { return rsize_; }
    
    /*! Make sure row is constructed, return bit-vector */
    bvector_type_ptr construct_row(size_type row);

    /*! Make sure row is copy-constructed, return bit-vector */
    bvector_type_ptr construct_row(size_type row, const bvector_type& bv);

    /*! destruct/deallocate row */
    void destruct_row(size_type row);

    /*! clear row bit-vector */
    void clear_row(size_type row, bool free_mem);

    ///@}
    
    
    // ------------------------------------------------------------
    /*! @name octet access and transposition                     */
    ///@{

    /*!
        Bit-transpose an octet and assign it to a bit-matrix
     
        @param pos - column position in the matrix
        @param octet_idx - octet based row position (1 octet - 8 rows)
        @param octet - value to assign
    */
    void set_octet(size_type pos, size_type octet_idx, unsigned char octet);

    /*!
        Bit-transpose and insert an octet and assign it to a bit-matrix
     
        @param pos - column position in the matrix
        @param octet_idx - octet based row position (1 octet - 8 rows)
        @param octet - value to assign
    */
    void insert_octet(size_type pos, size_type octet_idx, unsigned char octet);

    /*!
        return octet from the matrix
     
        @param pos - column position in the matrix
        @param octet_idx - octet based row position (1 octet - 8 rows)
    */
    unsigned char get_octet(size_type pos, size_type octet_idx) const BMNOEXCEPT;
    
    /*!
        Compare vector[pos] with octet
     
        It uses regulat comparison of chars to comply with the (signed)
        char sort order.
     
        @param pos - column position in the matrix
        @param octet_idx - octet based row position (1 octet - 8 rows)
        @param octet - octet value to compare
     
        @return 0 - equal, -1 - less(vect[pos] < octet), 1 - greater
    */
    int compare_octet(size_type pos,
                      size_type octet_idx, char octet) const BMNOEXCEPT;
    
    ///@}

public:

    // ------------------------------------------------------------
    /*! @name Utility function                                   */
    ///@{
    
    /// Test if 4 rows from i are not NULL
    bool test_4rows(unsigned i) const BMNOEXCEPT;

    /// Get low level internal access to
    const bm::word_t* get_block(size_type p,
                                unsigned i, unsigned j) const BMNOEXCEPT;
    
    unsigned get_half_octet(size_type pos, size_type row_idx) const BMNOEXCEPT;

    /*!
        \brief run memory optimization for all bit-vector rows
        \param temp_block - pre-allocated memory block to avoid re-allocs
        \param opt_mode - requested compression depth
        \param stat - memory allocation statistics after optimization
    */
    void optimize(
        bm::word_t* temp_block = 0,
        typename bvector_type::optmode opt_mode = bvector_type::opt_compress,
        typename bvector_type::statistics* stat = 0);
    
    /*! Optimize block in all planes
        @internal
    */
    void optimize_block(block_idx_type nb);

    ///@}


protected:
    void allocate_rows(size_type rsize);
    void free_rows() BMNOEXCEPT;

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

/**
    Base class for bit-transposed sparse vector construction
 
    @ingroup bmagic
    @internal
*/
template<typename Val, typename BV, unsigned MAX_SIZE>
class base_sparse_vector
{
public:
    enum bit_plains
    {
        sv_plains = (sizeof(Val) * 8 * MAX_SIZE + 1),
        sv_value_plains = (sizeof(Val) * 8 * MAX_SIZE)
    };

    enum vector_capacity
    {
        max_vector_size = MAX_SIZE
    };

    typedef Val                                      value_type;
    typedef BV                                       bvector_type;
    typedef typename BV::size_type                   size_type;
    typedef bvector_type*                            bvector_type_ptr;
    typedef const bvector_type*                      bvector_type_const_ptr;
    typedef const value_type&                        const_reference;
    typedef typename BV::allocator_type              allocator_type;
    typedef typename bvector_type::allocation_policy allocation_policy_type;
    typedef typename bvector_type::enumerator        bvector_enumerator_type;
    typedef typename allocator_type::allocator_pool_type allocator_pool_type;
    typedef bm::basic_bmatrix<BV>                        bmatrix_type;

public:
    base_sparse_vector();

    base_sparse_vector(bm::null_support        null_able,
                       allocation_policy_type  ap,
                       size_type               bv_max_size,
                       const allocator_type&   alloc);
    
    base_sparse_vector(const base_sparse_vector<Val, BV, MAX_SIZE>& bsv);

#ifndef BM_NO_CXX11
    /*! move-ctor */
    base_sparse_vector(base_sparse_vector<Val, BV, MAX_SIZE>&& bsv) BMNOEXCEPT
    {
        bmatr_.swap(bsv.bmatr_);
        size_ = bsv.size_;
        effective_plains_ = bsv.effective_plains_;
        bsv.size_ = 0;
    }
#endif

    void swap(base_sparse_vector<Val, BV, MAX_SIZE>& bsv) BMNOEXCEPT;

    size_type size() const BMNOEXCEPT { return size_; }
    
    void resize(size_type new_size);
    
    void clear_range(size_type left, size_type right, bool set_null);

    /*! \brief resize to zero, free memory
        @param free_mem - fully destroys the plane vectors if true
    */
    void clear_all(bool free_mem = true) BMNOEXCEPT;
    
    /*! return true if empty */
    bool empty() const BMNOEXCEPT { return size() == 0; }

public:

    // ------------------------------------------------------------
    /*! @name Various traits                                     */
    //@{
    /**
        \brief check if container supports NULL(unassigned) values
    */
    bool is_nullable() const BMNOEXCEPT
        { return bmatr_.get_row(this->null_plain()) != 0; }

    /**
        \brief Get bit-vector of assigned values or NULL
        (if not constructed that way)
    */
    const bvector_type* get_null_bvector() const BMNOEXCEPT
        { return bmatr_.get_row(this->null_plain()); }
    
    /** \brief test if specified element is NULL
        \param idx - element index
        \return true if it is NULL false if it was assigned or container
        is not configured to support assignment flags
    */
    bool is_null(size_type idx) const BMNOEXCEPT;
    

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
    bvector_type_const_ptr
    get_plain(unsigned i) const BMNOEXCEPT { return bmatr_.row(i); }

    /*!
        \brief get total number of bit-plains in the vector
    */
    static unsigned plains() BMNOEXCEPT { return value_bits(); }

    /** Number of stored bit-plains (value plains + extra */
    static unsigned stored_plains() BMNOEXCEPT { return value_bits()+1; }


    /** Number of effective bit-plains in the value type */
    unsigned effective_plains() const BMNOEXCEPT
                                { return effective_plains_ + 1; }

    /*!
        \brief get access to bit-plain as is (can return NULL)
    */
    bvector_type_ptr plain(unsigned i) BMNOEXCEPT { return bmatr_.get_row(i); }
    bvector_type_const_ptr plain(unsigned i) const BMNOEXCEPT
                                    { return bmatr_.get_row(i); }

    bvector_type* get_null_bvect() { return bmatr_.get_row(this->null_plain());}

    /*!
        \brief free memory in bit-plain
    */
    void free_plain(unsigned i) { bmatr_.destruct_row(i); }
    
    /*!
        return mask of allocated bit-plains
        1 in the mask - means bit-plain N is present
        returns 64-bit unsigned mask for sub 64-bit types (like int)
        unallocated mask bits will be zero extended
     
        @return 64-bit mask
        @internal
    */
    bm::id64_t get_plains_mask(unsigned element_idx) const BMNOEXCEPT;

    /*!
        get read-only access to inetrnal bit-matrix
    */
    const bmatrix_type& get_bmatrix() const BMNOEXCEPT { return bmatr_; }
    ///@}
    
    /*!
        \brief run memory optimization for all bit-vector rows
        \param temp_block - pre-allocated memory block to avoid unnecessary re-allocs
        \param opt_mode - requested compression depth
        \param stat - memory allocation statistics after optimization
    */
    void optimize(bm::word_t* temp_block = 0,
                  typename bvector_type::optmode opt_mode = bvector_type::opt_compress,
                  typename bvector_type::statistics* stat = 0);

    /*!
        @brief Calculates memory statistics.

        Function fills statistics structure containing information about how
        this vector uses memory and estimation of max. amount of memory
        bvector needs to serialize itself.

        @param st - pointer on statistics structure to be filled in.

        @sa statistics
    */
    void calc_stat(typename bvector_type::statistics* st) const BMNOEXCEPT;

    /*!
        \brief check if another sparse vector has the same content and size
     
        \param sv        - sparse vector for comparison
        \param null_able - flag to consider NULL vector in comparison (default)
                           or compare only value content plains
     
        \return true, if it is the same
    */
    bool equal(const base_sparse_vector<Val, BV, MAX_SIZE>& sv,
               bm::null_support null_able = bm::use_null) const BMNOEXCEPT;

protected:
    void copy_from(const base_sparse_vector<Val, BV, MAX_SIZE>& bsv);

    /*!
        clear column in all value plains
        \param plain_idx - row (plain index to start from)
        \param idx       - bit (column) to clear
    */
    void clear_value_plains_from(unsigned plain_idx, size_type idx);

    /*!
        insert false (clear) column in all value plains
        \param plain_idx - row (plain index to start from)
        \param idx       - bit (column) to clear insert
    */
    void insert_clear_value_plains_from(unsigned plain_idx, size_type idx);
    
    /*!
        erase bit (column) from all plains
        \param idx - bit (column) to erase
    */
    void erase_column(size_type idx);
    
    /*!
        insert (NOT) NULL value
    */
    void insert_null(size_type idx, bool not_null);

protected:
    typedef typename bvector_type::block_idx_type block_idx_type;

    /** Number of total bit-plains in the value type*/
    static unsigned value_bits() BMNOEXCEPT
    {
        return base_sparse_vector<Val, BV, MAX_SIZE>::sv_value_plains;
    }
    
    /** plain index for the "NOT NULL" flags plain */
    static unsigned null_plain() BMNOEXCEPT { return value_bits(); }
    
    /** optimize block in all matrix plains */
    void optimize_block(block_idx_type nb)
    {
        bmatr_.optimize_block(nb);
    }

    /**
        Perform copy_range() on a set of plains
    */
    void copy_range_plains(
        const base_sparse_vector<Val, BV, MAX_SIZE>& bsv,
        typename base_sparse_vector<Val, BV, MAX_SIZE>::size_type left,
        typename base_sparse_vector<Val, BV, MAX_SIZE>::size_type right,
        bm::null_support splice_null);

protected:
    bmatrix_type             bmatr_;              ///< bit-transposed matrix
    size_type                size_;               ///< array size
    unsigned                 effective_plains_;

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
basic_bmatrix<BV>::~basic_bmatrix() BMNOEXCEPT
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
basic_bmatrix<BV>::basic_bmatrix(basic_bmatrix<BV>&& bbm) BMNOEXCEPT
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
basic_bmatrix<BV>::row(size_type i) const BMNOEXCEPT
{
    BM_ASSERT(i < rsize_);
    return bv_rows_[i];
}

//---------------------------------------------------------------------

template<typename BV>
const typename basic_bmatrix<BV>::bvector_type*
basic_bmatrix<BV>::get_row(size_type i) const BMNOEXCEPT
{
    BM_ASSERT(i < rsize_);
    return bv_rows_[i];
}

//---------------------------------------------------------------------

template<typename BV>
typename basic_bmatrix<BV>::bvector_type*
basic_bmatrix<BV>::get_row(size_type i) BMNOEXCEPT
{
    BM_ASSERT(i < rsize_);
    return bv_rows_[i];
}

//---------------------------------------------------------------------

template<typename BV>
bool basic_bmatrix<BV>::test_4rows(unsigned j) const BMNOEXCEPT
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
        bool b = bv_rows_[j + 0] || bv_rows_[j + 1] ||
                 bv_rows_[j + 2] || bv_rows_[j + 3];
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
        bv_rows_ = (bvector_type_ptr*)alloc_.alloc_ptr(unsigned(rsize));
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
void basic_bmatrix<BV>::free_rows() BMNOEXCEPT
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
    {
        alloc_.free_ptr(bv_rows_, unsigned(rsize_));
    }
    bv_rows_ = 0;
}

//---------------------------------------------------------------------

template<typename BV>
void basic_bmatrix<BV>::swap(basic_bmatrix<BV>& bbm) BMNOEXCEPT
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
void basic_bmatrix<BV>::clear_row(size_type row, bool free_mem)
{
    BM_ASSERT(row < rsize_);
    bvector_type_ptr bv = bv_rows_[row];
    if (bv)
    {
        if (free_mem)
        {
            destruct_bvector(bv);
            bv_rows_[row] = 0;
        }
        else
        {
            bv->clear(true);
        }
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

template<typename BV>
const bm::word_t*
basic_bmatrix<BV>::get_block(size_type p,
                             unsigned i, unsigned j) const BMNOEXCEPT
{
    bvector_type_const_ptr bv = this->row(p);
    if (bv)
    {
        const typename bvector_type::blocks_manager_type& bman =
                                            bv->get_blocks_manager();
        return bman.get_block_ptr(i, j);
    }
    return 0;
}


//---------------------------------------------------------------------

template<typename BV>
void basic_bmatrix<BV>::set_octet(size_type pos,
                                  size_type octet_idx,
                                  unsigned char octet)
{
    BM_ASSERT(octet_idx * 8u < rsize_);
    
    size_type oct = octet;
    size_type row = octet_idx * 8;
    size_type row_end = row + 8;
    for (; row < row_end; ++row)
    {
        bvector_type* bv = this->get_row(row);
        if (oct & 1u)
        {
            if (!bv)
            {
                bv = this->construct_row(row);
                bv->init();
            }
            bv->set_bit_no_check(pos);
        }
        else
        {
            if (bv)
                bv->clear_bit_no_check(pos);
        }
        oct >>= 1;
        if (!oct)
            break;
    } // for
    
    // clear the tail
    for (++row; row < row_end; ++row)
    {
        bvector_type* bv = this->get_row(row);
        if (bv)
            bv->clear_bit_no_check(pos);
    } // for
}

//---------------------------------------------------------------------

template<typename BV>
void basic_bmatrix<BV>::insert_octet(size_type pos,
                                     size_type octet_idx,
                                     unsigned char octet)
{
    BM_ASSERT(octet_idx * 8u < rsize_);
    
    size_type oct = octet;
    size_type row = octet_idx * 8;
    size_type row_end = row + 8;
    for (; row < row_end; ++row)
    {
        bvector_type* bv = this->get_row(row);
        if (oct & 1u)
        {
            if (!bv)
            {
                bv = this->construct_row(row);
                bv->init();
                bv->set_bit_no_check(pos);
            }
            else
            {
                bv->insert(pos, true);
            }
        }
        else
        {
            if (bv)
                bv->insert(pos, false);
        }
        oct >>= 1;
        if (!oct)
            break;
    } // for
    
    // clear the tail
    for (++row; row < row_end; ++row)
    {
        bvector_type* bv = this->get_row(row);
        if (bv)
            bv->insert(pos, false);
    } // for
}


//---------------------------------------------------------------------

template<typename BV>
unsigned char
basic_bmatrix<BV>::get_octet(size_type pos, size_type octet_idx) const BMNOEXCEPT
{
    unsigned v = 0;

    block_idx_type nb = (pos >>  bm::set_block_shift);
    unsigned i0 = unsigned(nb >> bm::set_array_shift); // top block address
    unsigned j0 = unsigned(nb &  bm::set_array_mask);  // address in sub-block

    const bm::word_t* blk;
    const bm::word_t* blka[8];
    unsigned nbit = unsigned(pos & bm::set_block_mask);
    unsigned nword  = unsigned(nbit >> bm::set_word_shift);
    unsigned mask0 = 1u << (nbit & bm::set_word_mask);
    
    unsigned row_idx = unsigned(octet_idx * 8);

    blka[0] = get_block(row_idx+0, i0, j0);
    blka[1] = get_block(row_idx+1, i0, j0);
    blka[2] = get_block(row_idx+2, i0, j0);
    blka[3] = get_block(row_idx+3, i0, j0);
    blka[4] = get_block(row_idx+4, i0, j0);
    blka[5] = get_block(row_idx+5, i0, j0);
    blka[6] = get_block(row_idx+6, i0, j0);
    blka[7] = get_block(row_idx+7, i0, j0);
    unsigned is_set;
    
    if ((blk = blka[0])!=0)
    {
        if (blk == FULL_BLOCK_FAKE_ADDR)
            is_set = 1;
        else
            is_set = (BM_IS_GAP(blk)) ? bm::gap_test_unr(BMGAP_PTR(blk), nbit) : (blk[nword] & mask0);
        v |= (unsigned)bool(is_set);
    }
    if ((blk = blka[1])!=0)
    {
        if (blk == FULL_BLOCK_FAKE_ADDR)
            is_set = 1;
        else
            is_set = (BM_IS_GAP(blk)) ? bm::gap_test_unr(BMGAP_PTR(blk), nbit) : (blk[nword] & mask0);
        v |= unsigned(bool(is_set)) << 1u;
    }
    if ((blk = blka[2])!=0)
    {
        if (blk == FULL_BLOCK_FAKE_ADDR)
            is_set = 1;
        else
            is_set = (BM_IS_GAP(blk)) ? bm::gap_test_unr(BMGAP_PTR(blk), nbit) : (blk[nword] & mask0);
        v |= unsigned(bool(is_set)) << 2u;
    }
    if ((blk = blka[3])!=0)
    {
        if (blk == FULL_BLOCK_FAKE_ADDR)
            is_set = 1;
        else
            is_set = (BM_IS_GAP(blk)) ? bm::gap_test_unr(BMGAP_PTR(blk), nbit) : (blk[nword] & mask0);
        v |= unsigned(bool(is_set)) << 3u;
    }
    
    
    if ((blk = blka[4])!=0)
    {
        if (blk == FULL_BLOCK_FAKE_ADDR)
            is_set = 1;
        else
            is_set = (BM_IS_GAP(blk)) ? bm::gap_test_unr(BMGAP_PTR(blk), nbit) : (blk[nword] & mask0);
        v |= unsigned(bool(is_set)) << 4u;
    }
    if ((blk = blka[5])!=0)
    {
        if (blk == FULL_BLOCK_FAKE_ADDR)
            is_set = 1;
        else
            is_set = (BM_IS_GAP(blk)) ? bm::gap_test_unr(BMGAP_PTR(blk), nbit) : (blk[nword] & mask0);
        v |= unsigned(bool(is_set)) << 5u;
    }
    if ((blk = blka[6])!=0)
    {
        if (blk == FULL_BLOCK_FAKE_ADDR)
            is_set = 1;
        else
            is_set = (BM_IS_GAP(blk)) ? bm::gap_test_unr(BMGAP_PTR(blk), nbit) : (blk[nword] & mask0);
        v |= unsigned(bool(is_set)) << 6u;
    }
    if ((blk = blka[7])!=0)
    {
        if (blk == FULL_BLOCK_FAKE_ADDR)
            is_set = 1;
        else
            is_set = (BM_IS_GAP(blk)) ? bm::gap_test_unr(BMGAP_PTR(blk), nbit) : (blk[nword] & mask0);
        v |= unsigned(bool(is_set)) << 7u;
    }
   
    return (unsigned char)v;
}

//---------------------------------------------------------------------

template<typename BV>
int basic_bmatrix<BV>::compare_octet(size_type pos,
                                     size_type octet_idx,
                                     char      octet) const BMNOEXCEPT
{
    char value = char(get_octet(pos, octet_idx));
    return (value > octet) - (value < octet);
}

//---------------------------------------------------------------------

template<typename BV>
unsigned
basic_bmatrix<BV>::get_half_octet(size_type pos, size_type row_idx) const BMNOEXCEPT
{
    unsigned v = 0;

    block_idx_type nb = (pos >>  bm::set_block_shift);
    unsigned i0 = unsigned(nb >> bm::set_array_shift); // top block address
    unsigned j0 = unsigned(nb &  bm::set_array_mask);  // address in sub-block

    const bm::word_t* blk;
    const bm::word_t* blka[4];
    unsigned nbit = unsigned(pos & bm::set_block_mask);
    unsigned nword  = unsigned(nbit >> bm::set_word_shift);
    unsigned mask0 = 1u << (nbit & bm::set_word_mask);

    blka[0] = get_block(row_idx+0, i0, j0);
    blka[1] = get_block(row_idx+1, i0, j0);
    blka[2] = get_block(row_idx+2, i0, j0);
    blka[3] = get_block(row_idx+3, i0, j0);
    unsigned is_set;
    
    if ((blk = blka[0])!=0)
    {
        if (blk == FULL_BLOCK_FAKE_ADDR)
            is_set = 1;
        else
            is_set = (BM_IS_GAP(blk)) ? bm::gap_test_unr(BMGAP_PTR(blk), nbit) : (blk[nword] & mask0);
        v |= unsigned(bool(is_set));
    }
    if ((blk = blka[1])!=0)
    {
        if (blk == FULL_BLOCK_FAKE_ADDR)
            is_set = 1;
        else
            is_set = (BM_IS_GAP(blk)) ? bm::gap_test_unr(BMGAP_PTR(blk), nbit) : (blk[nword] & mask0);
        v |= unsigned(bool(is_set)) << 1u;
    }
    if ((blk = blka[2])!=0)
    {
        if (blk == FULL_BLOCK_FAKE_ADDR)
            is_set = 1;
        else
            is_set = (BM_IS_GAP(blk)) ? bm::gap_test_unr(BMGAP_PTR(blk), nbit) : (blk[nword] & mask0);
        v |= unsigned(bool(is_set)) << 2u;
    }
    if ((blk = blka[3])!=0)
    {
        if (blk == FULL_BLOCK_FAKE_ADDR)
            is_set = 1;
        else
            is_set = (BM_IS_GAP(blk)) ? bm::gap_test_unr(BMGAP_PTR(blk), nbit) : (blk[nword] & mask0);
        v |= unsigned(bool(is_set)) << 3u;
    }
    return v;
}

//---------------------------------------------------------------------

template<typename BV>
void basic_bmatrix<BV>::optimize(bm::word_t* temp_block,
                  typename bvector_type::optmode opt_mode,
                  typename bvector_type::statistics* st)
{
    if (st)
        st->reset();

    BM_DECLARE_TEMP_BLOCK(tb);
    if (!temp_block)
        temp_block = tb;

    for (unsigned k = 0; k < rsize_; ++k)
    {
        bvector_type* bv = get_row(k);
        if (bv)
        {
            typename bvector_type::statistics stbv;
            stbv.reset();
            bv->optimize(temp_block, opt_mode, st ? &stbv : 0);
            if (st)
            {
                st->add(stbv);
            }
        }
    } // for k
}

//---------------------------------------------------------------------

template<typename BV>
void basic_bmatrix<BV>::optimize_block(block_idx_type nb)
{
    for (unsigned k = 0; k < rsize_; ++k)
    {
        bvector_type* bv = get_row(k);
        if (bv)
        {
            unsigned i, j;
            bm::get_block_coord(nb, i, j);
            typename bvector_type::blocks_manager_type& bman =
                                                bv->get_blocks_manager();
            bman.optimize_bit_block(i, j);
        }
    } // for k

}

//---------------------------------------------------------------------
//---------------------------------------------------------------------



template<class Val, class BV, unsigned MAX_SIZE>
base_sparse_vector<Val, BV, MAX_SIZE>::base_sparse_vector()
: bmatr_(sv_plains, allocation_policy_type(), bm::id_max, allocator_type()),
  size_(0),
  effective_plains_(0)
{
}

//---------------------------------------------------------------------

template<class Val, class BV, unsigned MAX_SIZE>
base_sparse_vector<Val, BV, MAX_SIZE>::base_sparse_vector(
        bm::null_support        null_able,
        allocation_policy_type  ap,
        size_type               bv_max_size,
        const allocator_type&       alloc)
: bmatr_(sv_plains, ap, bv_max_size, alloc),
  size_(0),
  effective_plains_(0)
{
    if (null_able == bm::use_null)
    {
        unsigned i = null_plain();
        bmatr_.construct_row(i)->init();
    }
}

//---------------------------------------------------------------------

template<class Val, class BV, unsigned MAX_SIZE>
base_sparse_vector<Val, BV, MAX_SIZE>::base_sparse_vector(
                        const base_sparse_vector<Val, BV, MAX_SIZE>& bsv)
: bmatr_(bsv.bmatr_),
  size_(bsv.size_),
  effective_plains_(bsv.effective_plains_)
{
}

//---------------------------------------------------------------------

template<class Val, class BV, unsigned MAX_SIZE>
void base_sparse_vector<Val, BV, MAX_SIZE>::copy_from(
                const base_sparse_vector<Val, BV, MAX_SIZE>& bsv)
{
    resize(bsv.size());
    effective_plains_ = bsv.effective_plains_;

    unsigned ni = this->null_plain();
    unsigned plains = bsv.stored_plains();
    for (size_type i = 0; i < plains; ++i)
    {
        bvector_type* bv = bmatr_.get_row(i);
        const bvector_type* bv_src = bsv.bmatr_.row(i);
        
        if (i == ni) // NULL plain copy
        {
            if (bv && !bv_src) // special case (copy from not NULL)
            {
                if (size_)
                    bv->set_range(0, size_-1);
                continue;
            }
        }

        if (bv)
            bmatr_.destruct_row(i);
        if (bv_src)
            bmatr_.construct_row(i, *bv_src);
    } // for i
}

//---------------------------------------------------------------------

template<class Val, class BV, unsigned MAX_SIZE>
void base_sparse_vector<Val, BV, MAX_SIZE>::swap(
                 base_sparse_vector<Val, BV, MAX_SIZE>& bsv) BMNOEXCEPT
{
    if (this != &bsv)
    {
        bmatr_.swap(bsv.bmatr_);

        bm::xor_swap(size_, bsv.size_);
        bm::xor_swap(effective_plains_, bsv.effective_plains_);
    }
}

//---------------------------------------------------------------------

template<class Val, class BV, unsigned MAX_SIZE>
void base_sparse_vector<Val, BV, MAX_SIZE>::clear_all(bool free_mem) BMNOEXCEPT
{
    unsigned plains = value_bits();
    for (size_type i = 0; i < plains; ++i)
        bmatr_.clear_row(i, free_mem);
    size_ = 0;
    bvector_type* bv_null = get_null_bvect();
    if (bv_null)
        bv_null->clear(true);
}

//---------------------------------------------------------------------

template<class Val, class BV, unsigned MAX_SIZE>
void base_sparse_vector<Val, BV, MAX_SIZE>::clear_range(
        typename base_sparse_vector<Val, BV, MAX_SIZE>::size_type left,
        typename base_sparse_vector<Val, BV, MAX_SIZE>::size_type right,
        bool set_null)
{
    if (right < left)
    {
        return clear_range(right, left, set_null);
    }
    unsigned plains = value_bits();
    for (unsigned i = 0; i < plains; ++i)
    {
        bvector_type* bv = this->bmatr_.get_row(i);
        if (bv)
            bv->set_range(left, right, false);
    } // for i
    
    if (set_null)
    {
        bvector_type* bv_null = this->get_null_bvect();
        if (bv_null)
            bv_null->set_range(left, right, false);
    }
}

//---------------------------------------------------------------------

template<class Val, class BV, unsigned MAX_SIZE>
void base_sparse_vector<Val, BV, MAX_SIZE>::resize(size_type sz)
{
    if (sz == size())  // nothing to do
        return;
    if (!sz) // resize to zero is an equivalent of non-destructive deallocation
    {
        clear_all();
        return;
    }
    if (sz < size()) // vector shrink
        clear_range(sz, this->size_-1, true);   // clear the tails and NULL vect
    size_ = sz;
}


//---------------------------------------------------------------------

template<class Val, class BV, unsigned MAX_SIZE>
bool base_sparse_vector<Val, BV, MAX_SIZE>::is_null(
                                    size_type idx) const BMNOEXCEPT
{
    const bvector_type* bv_null = get_null_bvector();
    return (bv_null) ? (!bv_null->test(idx)) : false;
}

//---------------------------------------------------------------------

template<class Val, class BV, unsigned MAX_SIZE>
void base_sparse_vector<Val, BV, MAX_SIZE>::insert_null(size_type idx,
                                                        bool      not_null)
{
    bvector_type* bv_null = this->get_null_bvect();
    if (bv_null)
        bv_null->insert(idx, not_null);
}

//---------------------------------------------------------------------

template<class Val, class BV, unsigned MAX_SIZE>
typename base_sparse_vector<Val, BV, MAX_SIZE>::bvector_type_ptr
                base_sparse_vector<Val, BV, MAX_SIZE>::get_plain(unsigned i)
{
    bvector_type_ptr bv = bmatr_.get_row(i);
    if (!bv)
    {
        bv = bmatr_.construct_row(i);
        bv->init();
        if (i > effective_plains_ && i < value_bits())
            effective_plains_ = i;
    }
    return bv;
}

//---------------------------------------------------------------------

template<class Val, class BV, unsigned MAX_SIZE>
bm::id64_t base_sparse_vector<Val, BV, MAX_SIZE>::get_plains_mask(
                                        unsigned element_idx) const BMNOEXCEPT
{
    BM_ASSERT(element_idx < MAX_SIZE);
    bm::id64_t mask = 0;
    bm::id64_t mask1 = 1;
    const unsigned plains = sizeof(value_type) * 8;
    unsigned bidx = 0;
    for (unsigned i = element_idx * plains; i < (element_idx+1) * plains; i+=4)
    {
        mask |= get_plain(i+0) ? (mask1 << (bidx+0)) : 0ull;
        mask |= get_plain(i+1) ? (mask1 << (bidx+1)) : 0ull;
        mask |= get_plain(i+2) ? (mask1 << (bidx+2)) : 0ull;
        mask |= get_plain(i+3) ? (mask1 << (bidx+3)) : 0ull;
        bidx += 4;
    } // for i
    return mask;
}

//---------------------------------------------------------------------

template<class Val, class BV, unsigned MAX_SIZE>
void base_sparse_vector<Val, BV, MAX_SIZE>::optimize(bm::word_t* temp_block,
                                    typename bvector_type::optmode opt_mode,
                                    typename bvector_type::statistics* st)
{
    typename bvector_type::statistics stbv;
    bmatr_.optimize(temp_block, opt_mode, &stbv);
    if (st)
        st->add(stbv);
    
    bvector_type* bv_null = this->get_null_bvect();
    
    unsigned stored_plains = this->stored_plains();
    for (unsigned j = 0; j < stored_plains; ++j)
    {
        bvector_type* bv = this->bmatr_.get_row(j);
        if (bv && (bv != bv_null)) // protect the NULL vector from de-allocation
        {
            // TODO: check if this can be done within optimize loop
            if (!bv->any())  // empty vector?
            {
                this->bmatr_.destruct_row(j);
                continue;
            }
        }
    } // for j
}

//---------------------------------------------------------------------

template<class Val, class BV, unsigned MAX_SIZE>
void base_sparse_vector<Val, BV, MAX_SIZE>::calc_stat(
                    typename bvector_type::statistics* st) const BMNOEXCEPT
{
    BM_ASSERT(st);
    
    st->reset();

    unsigned stored_plains = this->stored_plains();
    for (unsigned j = 0; j < stored_plains; ++j)
    {
        const bvector_type* bv = this->bmatr_.row(j);
        if (bv)
        {
            typename bvector_type::statistics stbv;
            bv->calc_stat(&stbv);
            st->add(stbv);
        }
        else
        {
            st->max_serialize_mem += 8;
        }
    } // for j
    
    // header accounting
    st->max_serialize_mem += 1 + 1 + 1 + 1 + 8 + (8 * this->stored_plains());
    st->max_serialize_mem += 1 + 8; // extra header fields for large bit-matrixes
}

//---------------------------------------------------------------------

template<class Val, class BV, unsigned MAX_SIZE>
void base_sparse_vector<Val, BV, MAX_SIZE>::clear_value_plains_from(
                                    unsigned plain_idx, size_type idx)
{
    for (unsigned i = plain_idx; i < sv_value_plains; ++i)
    {
        bvector_type* bv = this->bmatr_.get_row(i);
        if (bv)
            bv->clear_bit_no_check(idx);
    }
}

//---------------------------------------------------------------------

template<class Val, class BV, unsigned MAX_SIZE>
void base_sparse_vector<Val, BV, MAX_SIZE>::insert_clear_value_plains_from(
                                        unsigned plain_idx, size_type idx)
{
    for (unsigned i = plain_idx; i < sv_value_plains; ++i)
    {
        bvector_type* bv = this->bmatr_.get_row(i);
        if (bv)
            bv->insert(idx, false);
    }
}

//---------------------------------------------------------------------

template<class Val, class BV, unsigned MAX_SIZE>
void base_sparse_vector<Val, BV, MAX_SIZE>::erase_column(size_type idx)
{
    for (unsigned i = 0; i < sv_value_plains; ++i)
    {
        bvector_type* bv = this->bmatr_.get_row(i);
        if (bv)
            bv->erase(idx);
    }
}

//---------------------------------------------------------------------

template<class Val, class BV, unsigned MAX_SIZE>
bool base_sparse_vector<Val, BV, MAX_SIZE>::equal(
            const base_sparse_vector<Val, BV, MAX_SIZE>& sv,
             bm::null_support null_able) const BMNOEXCEPT
{
    size_type arg_size = sv.size();
    if (this->size_ != arg_size)
    {
        return false;
    }
    unsigned plains = this->plains();
    for (unsigned j = 0; j < plains; ++j)
    {
        const bvector_type* bv = this->bmatr_.get_row(j);
        const bvector_type* arg_bv = sv.bmatr_.get_row(j);
        if (bv == arg_bv) // same NULL
            continue;
        // check if any not NULL and not empty
        if (!bv && arg_bv)
        {
            if (arg_bv->any())
                return false;
            continue;
        }
        if (bv && !arg_bv)
        {
            if (bv->any())
                return false;
            continue;
        }
        // both not NULL
        bool eq = bv->equal(*arg_bv);
        if (!eq)
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
        bool eq = bv_null->equal(*bv_null_arg);
        if (!eq)
            return false;
    }
    return true;
}

//---------------------------------------------------------------------

template<class Val, class BV, unsigned MAX_SIZE>
void base_sparse_vector<Val, BV, MAX_SIZE>::copy_range_plains(
        const base_sparse_vector<Val, BV, MAX_SIZE>& bsv,
        typename base_sparse_vector<Val, BV, MAX_SIZE>::size_type left,
        typename base_sparse_vector<Val, BV, MAX_SIZE>::size_type right,
        bm::null_support splice_null)
{
    bvector_type* bv_null = get_null_bvect();
    const bvector_type* bv_null_arg = bsv.get_null_bvector();
    unsigned plains;
    if (bv_null)
    {
        plains = this->stored_plains();
        if (bv_null_arg && (splice_null == bm::use_null))
            bv_null->copy_range(*bv_null_arg, left, right);
        --plains;
    }
    else
        plains = this->plains();

    for (unsigned j = 0; j < plains; ++j)
    {
        const bvector_type* arg_bv = bsv.bmatr_.row(j);
        if (arg_bv)
        {
            bvector_type* bv = this->bmatr_.get_row(j);
            if (!bv)
                bv = this->get_plain(j);
            bv->copy_range(*arg_bv, left, right);
        }
    } // for j

}

//---------------------------------------------------------------------

} // namespace

#endif
