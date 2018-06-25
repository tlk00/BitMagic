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
    \brief Compressed sparse constainer compressed_sparse_vector<> for integer types
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
   \brief compressed sparse vector for NULL-able sparse vectors
 
   \ingroup sv
*/
template<class Val, class SV>
class compressed_sparse_vector
{
public:

    enum bit_plains
    {
        sv_plains = (sizeof(Val) * 8 + 1),
        sv_value_plains = (sizeof(Val) * 8)
    };

    typedef Val                                      value_type;
    typedef const value_type&                        const_reference;
    typedef bm::id_t                                 size_type;
    typedef SV                                       sparse_vector_type;
    typedef typename SV::bvector_type                bvector_type;
    typedef bvector_type*                            bvector_type_ptr;
    typedef typename bvector_type::allocator_type    allocator_type;
    typedef typename bvector_type::allocation_policy allocation_policy_type;
    typedef typename bvector_type::blocks_count      bvector_blocks_psum_type;
    
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
        reference(compressed_sparse_vector<Val, SV>& csv, size_type idx) BMNOEXEPT
        : csv_(csv), idx_(idx)
        {}
        operator value_type() const { return csv_.get(idx_); }
        bool operator==(const reference& ref) const
                                { return bool(*this) == bool(ref); }
        bool is_null() const { return csv_.is_null(idx_); }
    private:
        compressed_sparse_vector<Val, SV>& csv_;
        size_type                          idx_;
    };

public:
    // ------------------------------------------------------------
    /*! @name Construction and assignment  */
    //@{

    compressed_sparse_vector(allocation_policy_type ap = allocation_policy_type(),
                             size_type bv_max_size = bm::id_max,
                             const allocator_type&   alloc  = allocator_type());
    /*! copy-ctor */
    compressed_sparse_vector(const compressed_sparse_vector<Val, SV>& csv);
    
    /*! Compression constructor, loads data from non-compressed sparse vector */
    compressed_sparse_vector(const sparse_vector_type& sv);
    
    /*! copy assignmment operator */
    compressed_sparse_vector<Val,SV>& operator = (const compressed_sparse_vector<Val, SV>& csv)
    {
        if (this != &csv)
            sv_ = csv.sv_;
        return *this;
    }
    
#ifndef BM_NO_CXX11
    /*! move-ctor */
    compressed_sparse_vector(compressed_sparse_vector<Val,SV>&& csv) BMNOEXEPT;

    /*! move assignmment operator */
    compressed_sparse_vector<Val,SV>& operator=(compressed_sparse_vector<Val,SV>&& csv) BMNOEXEPT
    {
        if (this != &csv)
        {
            clear();
            sv_.swap(csv.sv_);
            max_id_ = csv.max_id_; in_sync_ = csv.in_sync_;
            if (in_sync_)
                bv_blocks_.copy_from(csv.bv_blocks_);
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
    value_type at(bm::id_t idx) const;
    
    /*!
        \brief get specified element without bounds checking
        \param idx - element index
        \return value of the element
    */
    value_type get(bm::id_t idx) const;

    
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
    */
    bool find_rank(bm::id_t rank, bm::id_t& pos) const;

    //@}
    
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
    bool equal(const compressed_sparse_vector<Val, SV>& csv) const;
    //@}

    // ------------------------------------------------------------
    /*! @name Load-Export compressed vector with data */
    
    //@{

    /*!
        \brief set specified element with bounds checking and automatic resize
     
        Method cannot insert elements, so every new idx has to be greater or equal
        than what it used before. Elements must be loaded in a sorted order.
     
        \param idx - element index
        \param v   - element value
    */
    void push_back(size_type idx, value_type v);
    
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


    // ------------------------------------------------------------
    /*! @name Read-acess fast sync */
    //@{
    /*!
        \brief Re-calculate prefix sum table used for rank search
        \param force - force recalculation even if it is already recalculated
    */
    void sync(bool force = false);
    
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
    const bvector_type_ptr get_plain(unsigned i) const { return sv_.get_plain(i); }

    /*!
        Number of effective bit-plains in the value type
    */
    unsigned effective_plains() const { return sv_.effective_plains(); }
    
    /*!
        \brief get total number of bit-plains in the vector
    */
    static unsigned plains() { return sparse_vector_type::plains(); }

    /*!
        \brief access dense vector
    */
    const sparse_vector_type& get_sv() const { return sv_; }
    
    ///@}
    
protected:
    /*!
        \brief Resolve logical address to access via rank compressed address
     
        \param idx    - input id to resolve
        \param idx_to - output id
     
        \return true if id is known and resolved successfully
    */
    bool resolve(bm::id_t idx, bm::id_t* idx_to) const;

private:
    sparse_vector_type            sv_;       ///< transpose-sparse vector for "dense" packing
    bm::id_t                      max_id_;   ///< control variable for sorted load
    bool                          in_sync_;  ///< flag if prefix sum is in-sync with vector
    bvector_blocks_psum_type      bv_blocks_; ///< prefix sum for rank translation
};

//---------------------------------------------------------------------
//---------------------------------------------------------------------

template<class Val, class SV>
compressed_sparse_vector<Val, SV>::compressed_sparse_vector(allocation_policy_type ap,
                                                            size_type bv_max_size,
                                                            const allocator_type&   alloc)
    : sv_(bm::use_null, ap, bv_max_size, alloc),
      max_id_(0), in_sync_(false)
{}

//---------------------------------------------------------------------

template<class Val, class SV>
compressed_sparse_vector<Val, SV>::compressed_sparse_vector(
                          const compressed_sparse_vector<Val, SV>& csv)
: sv_(csv.sv_),
  max_id_(csv.max_id_),
  in_sync_(csv.in_sync_)
{
    if (in_sync_)
    {
        bv_blocks_.copy_from(csv.bv_blocks_);
    }
}

//---------------------------------------------------------------------

template<class Val, class SV>
compressed_sparse_vector<Val, SV>::compressed_sparse_vector(compressed_sparse_vector<Val,SV>&& csv) BMNOEXEPT
: sv_(bm::use_null),
  max_id_(0), in_sync_(false)
{
    if (this != &csv)
    {
        sv_.swap(csv.sv_);
        max_id_ = csv.max_id_; in_sync_ = csv.in_sync_;
        if (in_sync_)
            bv_blocks_.copy_from(csv.bv_blocks_);
    }
}

//---------------------------------------------------------------------

template<class Val, class SV>
compressed_sparse_vector<Val, SV>::compressed_sparse_vector(
    const typename compressed_sparse_vector<Val, SV>::sparse_vector_type& sv)
: sv_(bm::use_null),
  max_id_(0), in_sync_(false)
{
    load_from(sv);
    sync();
}


//---------------------------------------------------------------------

template<class Val, class SV>
typename compressed_sparse_vector<Val, SV>::size_type
compressed_sparse_vector<Val, SV>::size() const
{
    return max_id_+1;
}

//---------------------------------------------------------------------

template<class Val, class SV>
void compressed_sparse_vector<Val, SV>::push_back(size_type idx, value_type v)
{
    if (sv_.empty())
    {
    }
    else
    if (idx <= max_id_)
    {
        sv_.throw_range_error("compressed sparse vector push_back() range error");
    }
    
    bvector_type* bv_null = sv_.get_null_bvect();
    BM_ASSERT(bv_null);
    
    bv_null->set(idx);
    sv_.push_back_no_null(v);
    
    max_id_ = idx;
    in_sync_ = false;
}

//---------------------------------------------------------------------

template<class Val, class SV>
bool compressed_sparse_vector<Val, SV>::equal(const compressed_sparse_vector<Val, SV>& csv) const
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
void compressed_sparse_vector<Val, SV>::load_from(const sparse_vector_type& sv_src)
{
    sv_.clear();
    max_id_ = 0;
    
    const bvector_type* bv_null_src = sv_src.get_null_bvector();
    if (!bv_null_src)
    {
        BM_ASSERT(bv_null_src);
        return;
    }
    
    bvector_type* bv_null = sv_.get_null_bvect();
    BM_ASSERT(bv_null);
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
    
    unsigned count = bv_null->count(); // set correct sizes
    sv_.resize(count);
    
    in_sync_ = false;
    
    bool found = bv_null->find_reverse(max_id_);
    if (!found)
    {
        BM_ASSERT(!bv_null->any());
        max_id_ = 0;
    }
}

//---------------------------------------------------------------------

template<class Val, class SV>
void compressed_sparse_vector<Val, SV>::load_to(sparse_vector_type& sv) const
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
void compressed_sparse_vector<Val, SV>::sync(bool force)
{
    if (in_sync_ && force == false)
        return;  // nothing to do
    const bvector_type* bv_null = sv_.get_null_bvector();
    if (bv_null)
        bv_null->running_count_blocks(&bv_blocks_); // compute popcount prefix list
    in_sync_ = true;
}

//---------------------------------------------------------------------

template<class Val, class SV>
bool compressed_sparse_vector<Val, SV>::resolve(bm::id_t idx, bm::id_t* idx_to) const
{
    BM_ASSERT(idx_to);
    
    const bvector_type* bv_null = sv_.get_null_bvector();
    if (in_sync_)
    {
        *idx_to = bv_null->count_to_test(idx, bv_blocks_);
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
typename compressed_sparse_vector<Val, SV>::value_type
compressed_sparse_vector<Val, SV>::at(bm::id_t idx) const
{
    bm::id_t sv_idx;
    bool found = resolve(idx, &sv_idx);
    if (!found)
    {
        sv_.throw_range_error("compressed collection item not found");
    }
    return sv_.at(--sv_idx);
}

//---------------------------------------------------------------------

template<class Val, class SV>
typename compressed_sparse_vector<Val, SV>::value_type
compressed_sparse_vector<Val, SV>::get(bm::id_t idx) const
{
    bm::id_t sv_idx;
    bool found = resolve(idx, &sv_idx);
    if (!found)
        return value_type();
    
    return sv_.get(--sv_idx);
}

//---------------------------------------------------------------------

template<class Val, class SV>
bool compressed_sparse_vector<Val, SV>::is_null(size_type idx) const
{
    const bvector_type* bv_null = sv_.get_null_bvector();
    BM_ASSERT(bv_null);
    return !(bv_null->test(idx));
}

//---------------------------------------------------------------------

template<class Val, class SV>
void compressed_sparse_vector<Val, SV>::optimize(bm::word_t*  temp_block,
                    typename bvector_type::optmode opt_mode,
                    statistics* stat)
{
    sv_.optimize(temp_block, opt_mode, (typename sparse_vector_type::statistics*)stat);
}

//---------------------------------------------------------------------

template<class Val, class SV>
void compressed_sparse_vector<Val, SV>::clear() BMNOEXEPT
{
    sv_.clear();
    in_sync_ = false;  max_id_ = 0;
}

//---------------------------------------------------------------------

template<class Val, class SV>
const typename compressed_sparse_vector<Val, SV>::bvector_type*
compressed_sparse_vector<Val, SV>::get_null_bvector() const
{
    return sv_.get_null_bvector();
}

//---------------------------------------------------------------------

template<class Val, class SV>
bool compressed_sparse_vector<Val, SV>::find_rank(bm::id_t rank, bm::id_t& pos) const
{
    BM_ASSERT(rank);

    const bvector_type* bv_null = get_null_bvector();
    bool b = bv_null->find_rank(rank, 0, pos);
    return b;
}

//---------------------------------------------------------------------

} // namespace bm

#include "bmundef.h"


#endif
