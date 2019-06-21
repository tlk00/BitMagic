#ifndef BMSPARSEVEC_COMPR_H__INCLUDED__
#define BMSPARSEVEC_COMPR_H__INCLUDED__
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

/*! \file bmsparsevec_compr.h
    \brief Compressed sparse container rsc_sparse_vector<> for integer types
*/

#include <memory.h>

#ifndef BM_NO_STL
#include <stdexcept>
#endif

#include "bmsparsevec.h"
#include "bmdef.h"

namespace bm
{


/*!
   \brief Rank-Select compressed sparse vector
 
   Container uses Rank-Select method of compression, where
   all NULL columns gets dropped, effective address of columns is computed
   using address bit-vector.

   Use for cases, where memory efficiency is preferable over
   single element access latency.
 
   \ingroup sv
*/
template<class Val, class SV>
class rsc_sparse_vector
{
public:
    enum bit_plains
    {
        sv_plains = (sizeof(Val) * 8 + 1),
        sv_value_plains = (sizeof(Val) * 8)
    };

    typedef Val                                      value_type;
    typedef const value_type&                        const_reference;
    typedef typename SV::size_type                   size_type;
    typedef SV                                       sparse_vector_type;
    typedef typename SV::const_iterator              sparse_vector_const_iterator;
    typedef typename SV::bvector_type                bvector_type;
    typedef bvector_type*                            bvector_type_ptr;
    typedef const bvector_type*                      bvector_type_const_ptr;
    typedef typename bvector_type::allocator_type    allocator_type;
    typedef typename bvector_type::allocation_policy allocation_policy_type;
    typedef typename bvector_type::rs_index_type     rs_index_type;
    typedef typename bvector_type::enumerator        bvector_enumerator_type;

    enum vector_capacity
    {
        max_vector_size = 1
    };

    struct is_remap_support { enum trait { value = false }; };
    struct is_rsc_support { enum trait { value = true }; };

public:
    /*! Statistical information about  memory allocation details. */
    struct statistics : public bv_statistics
    {};
    
public:
    /**
         Reference class to access elements via common [] operator
    */
    class reference
    {
    public:
        reference(rsc_sparse_vector<Val, SV>& csv, size_type idx) BMNOEXEPT
        : csv_(csv), idx_(idx)
        {}
        operator value_type() const { return csv_.get(idx_); }
        bool operator==(const reference& ref) const
                                { return bool(*this) == bool(ref); }
        bool is_null() const { return csv_.is_null(idx_); }
    private:
        rsc_sparse_vector<Val, SV>& csv_;
        size_type                          idx_;
    };

public:
    // ------------------------------------------------------------
    /*! @name Construction and assignment  */
    //@{

    rsc_sparse_vector(bm::null_support null_able = bm::use_null,
                      allocation_policy_type ap = allocation_policy_type(),
                      size_type bv_max_size = bm::id_max,
                      const allocator_type&   alloc  = allocator_type());
    ~rsc_sparse_vector();
    
    /*! copy-ctor */
    rsc_sparse_vector(const rsc_sparse_vector<Val, SV>& csv);
    
    
    /*! copy assignmment operator */
    rsc_sparse_vector<Val,SV>& operator = (const rsc_sparse_vector<Val, SV>& csv)
    {
        if (this != &csv)
        {
            sv_ = csv.sv_;
            max_id_ = csv.max_id_;
            in_sync_ = csv.in_sync_;
            if (in_sync_)
            {
                bv_blocks_ptr_->copy_from(*(csv.bv_blocks_ptr_));
            }
        }
        return *this;
    }
    
#ifndef BM_NO_CXX11
    /*! move-ctor */
    rsc_sparse_vector(rsc_sparse_vector<Val,SV>&& csv) BMNOEXEPT;

    /*! move assignmment operator */
    rsc_sparse_vector<Val,SV>& operator=(rsc_sparse_vector<Val,SV>&& csv) BMNOEXEPT
    {
        if (this != &csv)
        {
            clear();
            sv_.swap(csv.sv_);
            max_id_ = csv.max_id_; in_sync_ = csv.in_sync_;
            if (in_sync_)
            {
                bv_blocks_ptr_->copy_from(*(csv.bv_blocks_ptr_));
            }
        }
        return *this;
    }
#endif

    //@}
    // ------------------------------------------------------------
    /*! @name Size, etc                                          */
    ///@{

    /*! \brief return size of the vector
        \return size of sparse vector
    */
    size_type size() const;
    
    /*! \brief return true if vector is empty
        \return true if empty
    */
    bool empty() const { return sv_.empty(); }
    
    ///@}

    // ------------------------------------------------------------
    /*! @name Element access */
    //@{

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
    value_type get(size_type idx) const;
    
    /*!
        \brief set specified element with bounds checking and automatic resize
     
        Method cannot insert elements, so every new idx has to be greater or equal
        than what it used before. Elements must be loaded in a sorted order.
     
        \param idx - element index
        \param v   - element value
    */
    void push_back(size_type idx, value_type v);
    
    /*!
        \brief set specified element with bounds checking and automatic resize
        \param idx - element index
        \param v   - element value
    */
    void set(size_type idx, value_type v);
    
    /*!
        \brief set specified element to NULL
        RSC vector actually erases element when it is set to NULL (expensive).
        \param idx - element index
    */
    void set_null(size_type idx);


    
    /** \brief test if specified element is NULL
        \param idx - element index
        \return true if it is NULL false if it was assigned or container
        is not configured to support assignment flags
    */
    bool is_null(size_type idx) const;
    
    /**
        \brief Get bit-vector of assigned values (or NULL)
    */
    const bvector_type* get_null_bvector() const;

    /**
        \brief find position of compressed element by its rank
        \param rank - rank  (virtual index in sparse vector)
        \param idx  - index (true position)
    */
    bool find_rank(size_type rank, size_type& idx) const;

    //@}
    
    // ------------------------------------------------------------
    /*! @name Export content to C-stype array                    */
    ///@{
    
    size_type decode(value_type* arr,
                     size_type   idx_from,
                     size_type   size,
                     bool        zero_mem = true) const;

    ///@}

    
    // ------------------------------------------------------------
    /*! @name Various traits                                     */
    //@{
    /**
        \brief if container supports NULL(unassigned) values (true)
    */
    bool is_nullable() const { return true; }
    
    /** \brief trait if sparse vector is "compressed" (true)
    */
    static
    bool is_compressed() { return true; }

    ///@}

    
    // ------------------------------------------------------------
    /*! @name Comparison */
    //@{

    /*!
        \brief check if another vector has the same content
        \return true, if it is the same
    */
    bool equal(const rsc_sparse_vector<Val, SV>& csv) const;
    //@}


    // ------------------------------------------------------------
    /*! @name Load-Export compressed vector with data */
    
    //@{

    
    /*!
        \brief Load compressed vector from a sparse vector (with NULLs)
        \param sv_src - source sparse vector
    */
    void load_from(const sparse_vector_type& sv_src);
    
    /*!
        \brief Exort compressed vector to a sparse vector (with NULLs)
        \param sv - target sparse vector
    */
    void load_to(sparse_vector_type& sv) const;
    
    //@}

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
                  statistics* stat = 0);
    
    /*! \brief resize to zero, free memory
    */
    void clear() BMNOEXEPT;
    
    /*!
        @brief Calculates memory statistics.

        Function fills statistics structure containing information about how
        this vector uses memory and estimation of max. amount of memory
        bvector needs to serialize itself.

        @param st - pointer on statistics structure to be filled in.

        @sa statistics
    */
    void calc_stat(struct rsc_sparse_vector<Val, SV>::statistics* st) const;

    ///@}


    // ------------------------------------------------------------
    /*! @name Fast access structues sync                         */
    //@{
    /*!
        \brief Re-calculate prefix sum table used for rank search
        \param force - force recalculation even if it is already recalculated
    */
    void sync(bool force);

    /*!
        \brief Re-calculate prefix sum table used for rank search (if necessary)
    */
    void sync() { sync(false); }

    /*!
        \brief returns true if prefix sum table is in sync with the vector
    */
    bool in_sync() const { return in_sync_; }
    
    /*!
        \brief Unsync the prefix sum table
    */
    void unsync() { in_sync_ = false; }
    ///@}

    // ------------------------------------------------------------
    /*! @name Access to internals                                */
    ///@{

    /*!
        \brief get access to bit-plain, function checks and creates a plain
        \return bit-vector for the bit plain
    */
    bvector_type_const_ptr get_plain(unsigned i) const { return sv_.get_plain(i); }

    bvector_type_ptr get_plain(unsigned i)  { return sv_.get_plain(i); }
    
    /*!
        Number of effective bit-plains in the value type
    */
    unsigned effective_plains() const { return sv_.effective_plains(); }
    
    /*!
        \brief get total number of bit-plains in the vector
    */
    static unsigned plains() { return sparse_vector_type::plains(); }

    /** Number of stored bit-plains (value plains + extra */
    static unsigned stored_plains() { return sparse_vector_type::stored_plains(); }

    /*!
        \brief access dense vector
    */
    const sparse_vector_type& get_sv() const { return sv_; }

    /*!
        \brief size of internal dense vector
    */
    size_type effective_size() const { return sv_.size(); }

    /**
        \brief Always 1 (non-matrix type)
    */
    size_type effective_vector_max() const { return 1; }

    ///@}
    
protected:
    enum octet_plains
    {
        sv_octet_plains = sizeof(value_type)
    };

    /*!
        \brief Resolve logical address to access via rank compressed address
     
        \param idx    - input id to resolve
        \param idx_to - output id
     
        \return true if id is known and resolved successfully
    */
    bool resolve(size_type idx, size_type* idx_to) const;
    
    void resize_internal(size_type sz) { sv_.resize_internal(sz); }
    size_type size_internal() const { return sv_.size(); }

    bool is_remap() const { return false; }
    size_t remap_size() const { return 0; }
    const unsigned char* get_remap_buffer() const { return 0; }
    unsigned char* init_remap_buffer() { return 0; }
    void set_remap() { }
    
    void push_back_no_check(size_type idx, value_type v);


private:
    void construct_bv_blocks();
    void free_bv_blocks();

protected:
    template<class SVect> friend class sparse_vector_scanner;
    template<class SVect> friend class sparse_vector_serializer;
    template<class SVect> friend class sparse_vector_deserializer;

private:
    sparse_vector_type            sv_;       ///< transpose-sparse vector for "dense" packing
    size_type                     max_id_;   ///< control variable for sorted load
    bool                          in_sync_;  ///< flag if prefix sum is in-sync with vector
    rs_index_type*                bv_blocks_ptr_ = 0; ///< prefix sum for rank translation
};

//---------------------------------------------------------------------
//---------------------------------------------------------------------

template<class Val, class SV>
rsc_sparse_vector<Val, SV>::rsc_sparse_vector(bm::null_support null_able,
                                              allocation_policy_type ap,
                                              size_type bv_max_size,
                                              const allocator_type&   alloc)
: sv_(null_able, ap, bv_max_size, alloc),
  max_id_(0), in_sync_(false)
{
    BM_ASSERT(null_able == bm::use_null);
    BM_ASSERT(int(sv_value_plains) == int(SV::sv_value_plains));
    
    construct_bv_blocks();
}

//---------------------------------------------------------------------

template<class Val, class SV>
rsc_sparse_vector<Val, SV>::~rsc_sparse_vector()
{
    free_bv_blocks();
}

//---------------------------------------------------------------------

template<class Val, class SV>
rsc_sparse_vector<Val, SV>::rsc_sparse_vector(
                          const rsc_sparse_vector<Val, SV>& csv)
: sv_(csv.sv_),
  max_id_(csv.max_id_),
  in_sync_(csv.in_sync_)
{
    BM_ASSERT(int(sv_value_plains) == int(SV::sv_value_plains));
    
    construct_bv_blocks();
    if (in_sync_)
    {
        bv_blocks_ptr_->copy_from(*(csv.bv_blocks_ptr_));
    }
}

//---------------------------------------------------------------------

template<class Val, class SV>
rsc_sparse_vector<Val, SV>::rsc_sparse_vector(rsc_sparse_vector<Val,SV>&& csv) BMNOEXEPT
: sv_(bm::use_null),
  max_id_(0), in_sync_(false)
{
    if (this != &csv)
    {
        sv_.swap(csv.sv_);
        max_id_ = csv.max_id_; in_sync_ = csv.in_sync_;
        
        bv_blocks_ptr_ = csv.bv_blocks_ptr_; csv.bv_blocks_ptr_ = 0;
    }
}

//---------------------------------------------------------------------

template<class Val, class SV>
typename rsc_sparse_vector<Val, SV>::size_type
rsc_sparse_vector<Val, SV>::size() const
{
    return max_id_+1;
}

//---------------------------------------------------------------------

template<class Val, class SV>
void rsc_sparse_vector<Val, SV>::push_back(size_type idx, value_type v)
{
    if (sv_.empty())
    {}
    else
    if (idx <= max_id_)
    {
        sv_.throw_range_error("compressed sparse vector push_back() range error");
    }
    push_back_no_check(idx, v);
}

//---------------------------------------------------------------------

template<class Val, class SV>
void rsc_sparse_vector<Val, SV>::push_back_no_check(size_type idx, value_type v)
{
    bvector_type* bv_null = sv_.get_null_bvect();
    BM_ASSERT(bv_null);
    
    bv_null->set_bit_no_check(idx);
    sv_.push_back_no_null(v);
    
    max_id_ = idx;
    in_sync_ = false;
}

//---------------------------------------------------------------------

template<class Val, class SV>
void rsc_sparse_vector<Val, SV>::set_null(size_type idx)
{
    bvector_type* bv_null = sv_.get_null_bvect();
    BM_ASSERT(bv_null);
    
    bool found = bv_null->test(idx);
    if (found)
    {
        size_type sv_idx = bv_null->count_range(0, idx);
        bv_null->clear_bit_no_check(idx);
        sv_.erase(--sv_idx);
    }
}

//---------------------------------------------------------------------

template<class Val, class SV>
void rsc_sparse_vector<Val, SV>::set(size_type idx, value_type v)
{
    bvector_type* bv_null = sv_.get_null_bvect();
    BM_ASSERT(bv_null);
    
    bool found = bv_null->test(idx);
    size_type sv_idx = bv_null->count_range(0, idx); // TODO: make test'n'count
    
    if (found)
    {
        sv_.set(--sv_idx, v);
    }
    else
    {
        sv_.insert_value_no_null(sv_idx, v);
        bv_null->set_bit_no_check(idx);

        if (idx > max_id_)
            max_id_ = idx;
        in_sync_ = false;
    }
}

//---------------------------------------------------------------------

template<class Val, class SV>
bool rsc_sparse_vector<Val, SV>::equal(
                    const rsc_sparse_vector<Val, SV>& csv) const
{
    if (this == &csv)
        return true;
    if (max_id_ != csv.max_id_)
        return false;
    bool same_sv = sv_.equal(csv.sv_);
    return same_sv;
}

//---------------------------------------------------------------------

template<class Val, class SV>
void rsc_sparse_vector<Val, SV>::load_from(
                                    const sparse_vector_type& sv_src)
{
    max_id_ = 0;

    bvector_type* bv_null = sv_.get_null_bvect();
    BM_ASSERT(bv_null);

    const bvector_type* bv_null_src = sv_src.get_null_bvector();
    if (!bv_null_src) // dense vector (no NULL columns)
    {
        sv_ = sv_src;
        BM_ASSERT(sv_.get_null_bvect());
    }
    else
    {
        sv_.clear();
        *bv_null = *bv_null_src;
        
        bm::rank_compressor<bvector_type> rank_compr; // re-used for plains
        
        unsigned src_plains = sv_src.plains();
        for (unsigned i = 0; i < src_plains; ++i)
        {
            const bvector_type* bv_src_plain = sv_src.get_plain(i);
            if (bv_src_plain)
            {
                bvector_type* bv_plain = sv_.get_plain(i);
                rank_compr.compress(*bv_plain, *bv_null, *bv_src_plain);
            }
        } // for
        size_type count = bv_null->count(); // set correct sizes
        sv_.resize(count);
    }
    
    sync(true);
}

//---------------------------------------------------------------------

template<class Val, class SV>
void rsc_sparse_vector<Val, SV>::load_to(sparse_vector_type& sv) const
{
    sv.clear();
    
    const bvector_type* bv_null_src = this->get_null_bvector();
    if (!bv_null_src)
    {
        BM_ASSERT(bv_null_src);
        return;
    }
    
    bvector_type* bv_null = sv.get_null_bvect();
    BM_ASSERT(bv_null);
    *bv_null = *bv_null_src;
    
    bm::rank_compressor<bvector_type> rank_compr; // re-used for plains

    unsigned src_plains = sv_.plains();
    for (unsigned i = 0; i < src_plains; ++i)
    {
        const bvector_type* bv_src_plain = sv_.get_plain(i);
        if (bv_src_plain)
        {
            bvector_type* bv_plain = sv.get_plain(i);
            rank_compr.decompress(*bv_plain, *bv_null, *bv_src_plain);
        }
    } // for
    sv.resize(this->size());
}

//---------------------------------------------------------------------

template<class Val, class SV>
void rsc_sparse_vector<Val, SV>::sync(bool force)
{
    if (in_sync_ && force == false)
        return;  // nothing to do
    const bvector_type* bv_null = sv_.get_null_bvector();
    BM_ASSERT(bv_null);
    bv_null->build_rs_index(bv_blocks_ptr_); // compute popcount prefix list
    
    // sync the max-id
    bool found = bv_null->find_reverse(max_id_);
    if (!found)
    {
        BM_ASSERT(!bv_null->any());
        max_id_ = 0;
    }
    in_sync_ = true;
}

//---------------------------------------------------------------------

template<class Val, class SV>
bool rsc_sparse_vector<Val, SV>::resolve(size_type idx, size_type* idx_to) const
{
    BM_ASSERT(idx_to);
    
    const bvector_type* bv_null = sv_.get_null_bvector();
    if (in_sync_)
    {
        *idx_to = bv_null->count_to_test(idx, *bv_blocks_ptr_);
    }
    else  // slow access
    {
        bool found = bv_null->test(idx);
        if (!found)
        {
            *idx_to = 0;
        }
        else
        {
            *idx_to = bv_null->count_range(0, idx);
        }
    }
    return bool(*idx_to);
}

//---------------------------------------------------------------------

template<class Val, class SV>
typename rsc_sparse_vector<Val, SV>::value_type
rsc_sparse_vector<Val, SV>::at(size_type idx) const
{
    size_type sv_idx;
    bool found = resolve(idx, &sv_idx);
    if (!found)
    {
        sv_.throw_range_error("compressed collection item not found");
    }
    return sv_.at(--sv_idx);
}

//---------------------------------------------------------------------

template<class Val, class SV>
typename rsc_sparse_vector<Val, SV>::value_type
rsc_sparse_vector<Val, SV>::get(size_type idx) const
{
    size_type sv_idx;
    bool found = resolve(idx, &sv_idx);
    if (!found)
        return value_type();
    
    return sv_.get(--sv_idx);
}

//---------------------------------------------------------------------

template<class Val, class SV>
bool rsc_sparse_vector<Val, SV>::is_null(size_type idx) const
{
    const bvector_type* bv_null = sv_.get_null_bvector();
    BM_ASSERT(bv_null);
    return !(bv_null->test(idx));
}

//---------------------------------------------------------------------

template<class Val, class SV>
void rsc_sparse_vector<Val, SV>::optimize(bm::word_t*  temp_block,
                    typename bvector_type::optmode opt_mode,
                    statistics* st)
{
    sv_.optimize(temp_block, opt_mode, (typename sparse_vector_type::statistics*)st);
    if (st)
    {
        st->memory_used += sizeof(rs_index_type);
        st->memory_used += bv_blocks_ptr_->get_total() *
             (sizeof(typename rs_index_type::size_type)
             + sizeof(typename rs_index_type::sb_pair_type));
    }
}

//---------------------------------------------------------------------

template<class Val, class SV>
void rsc_sparse_vector<Val, SV>::clear() BMNOEXEPT
{
    sv_.clear();
    in_sync_ = false;  max_id_ = 0;
}

//---------------------------------------------------------------------

template<class Val, class SV>
void rsc_sparse_vector<Val, SV>::calc_stat(
            struct rsc_sparse_vector<Val, SV>::statistics* st) const
{
    BM_ASSERT(st);
    sv_.calc_stat((typename sparse_vector_type::statistics*)st);
    if (st)
    {
        st->memory_used += sizeof(rs_index_type);
        st->memory_used += bv_blocks_ptr_->get_total() *
                   (sizeof(typename rs_index_type::size_type)
                   + sizeof(typename rs_index_type::sb_pair_type));
    }
}

//---------------------------------------------------------------------

template<class Val, class SV>
const typename rsc_sparse_vector<Val, SV>::bvector_type*
rsc_sparse_vector<Val, SV>::get_null_bvector() const
{
    return sv_.get_null_bvector();
}

//---------------------------------------------------------------------

template<class Val, class SV>
bool
rsc_sparse_vector<Val, SV>::find_rank(size_type rank, size_type& idx) const
{
    BM_ASSERT(rank);
    bool b;
    const bvector_type* bv_null = get_null_bvector();
    if (in_sync())
        b = bv_null->select(rank, idx, *bv_blocks_ptr_);
    else
        b = bv_null->find_rank(rank, 0, idx);
    return b;
}

//---------------------------------------------------------------------


template<class Val, class SV>
typename rsc_sparse_vector<Val, SV>::size_type
rsc_sparse_vector<Val, SV>::decode(value_type* arr,
                                   size_type   idx_from,
                                   size_type   size,
                                   bool        /*zero_mem*/) const
{
    BM_ASSERT(arr);
    BM_ASSERT(in_sync_);  // call sync() before decoding
    BM_ASSERT(bv_blocks_ptr_);
    
    if (size == 0)
        return 0;
    if (idx_from >= this->size())
        return 0;
    
    if ((bm::id_max - size) <= idx_from)
        size = bm::id_max - idx_from;

    const bvector_type* bv_null = sv_.get_null_bvector();

    size_type rank = bv_null->count_to(idx_from, *bv_blocks_ptr_);
    bool b = bv_null->test(idx_from);
    
    bvector_enumerator_type en_i = bv_null->get_enumerator(idx_from);
    size_type i = *en_i;
    if (idx_from + size <= i)  // empty space (all zeros)
    {
        ::memset(arr, 0, sizeof(value_type)*size);
        return size;
    }
    rank -= b;
    sparse_vector_const_iterator it = sv_.get_const_iterator(rank);
    i = 0;
    while (it.valid())
    {
        if (!en_i.valid())
            break;
        size_type en_idx = *en_i;
        while (idx_from < en_idx) // zero the empty prefix
        {
            arr[i] ^= arr[i];
            ++i; ++idx_from;
            if (i == size)
                return i;
        }
        BM_ASSERT(idx_from == en_idx);
        arr[i] = *it;
        ++i; ++idx_from;
        if (i == size)
            return i;
        
        en_i.advance();
        it.advance();
    } // while
    
    return i;
}

//---------------------------------------------------------------------

template<class Val, class SV>
void rsc_sparse_vector<Val, SV>::construct_bv_blocks()
{
    if (bv_blocks_ptr_)
        return;
    bv_blocks_ptr_ =
        (rs_index_type*) bm::aligned_new_malloc(sizeof(rs_index_type));
    bv_blocks_ptr_ = new(bv_blocks_ptr_) rs_index_type(); // placement new
}

//---------------------------------------------------------------------

template<class Val, class SV>
void rsc_sparse_vector<Val, SV>::free_bv_blocks()
{
    if (bv_blocks_ptr_)
    {
        bv_blocks_ptr_->~rs_index_type();
        bm::aligned_free(bv_blocks_ptr_);
    }
}

//---------------------------------------------------------------------


} // namespace bm

#include "bmundef.h"


#endif
