#ifndef BMSPARSEVEC_ALGO__H__INCLUDED__
#define BMSPARSEVEC_ALGO__H__INCLUDED__
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
/*! \file bmsparsevec_algo.h
    \brief Algorithms for sparse_vector<>
*/

#include "bmdef.h"
#include "bmsparsevec.h"
#include "bmaggregator.h"
#include "bmbuffer.h"
#include "bmdef.h"

#ifdef _MSC_VER
# pragma warning( disable: 4146 )
#endif


/** \defgroup svalgo Sparse vector algorithms
    Sparse vector algorithms
    \ingroup svector
 */


namespace bm
{


/*!
    \brief Clip dynamic range for signal higher than specified
    
    \param  svect - sparse vector to do clipping
    \param  high_bit - max bit (inclusive) allowed for this signal vector
    
    
    \ingroup svalgo
    \sa dynamic_range_clip_low
*/
template<typename SV>
void dynamic_range_clip_high(SV& svect, unsigned high_bit)
{
    unsigned sv_plains = svect.plains();
    
    BM_ASSERT(sv_plains > high_bit && high_bit > 0);
    
    typename SV::bvector_type bv_acc;
    unsigned i;
    
    // combine all the high bits into accumulator vector
    for (i = high_bit+1; i < sv_plains; ++i)
    {
        typename SV::bvector_type* bv_plain = svect.plain(i);
        if (bv_plain)
        {
            bv_acc.bit_or(*bv_plain);
            svect.free_plain(i);
        }
    } // for i
    
    // set all bits ON for all low vectors, which happen to be clipped
    for (i = high_bit; true; --i)
    {
        typename SV::bvector_type* bv_plain = svect.get_plain(i);
        bv_plain->bit_or(bv_acc);
        if (i == 0)
            break;
    } // for i
}


/*!
    \brief Clip dynamic range for signal lower than specified (boost)
    
    \param  svect - sparse vector to do clipping
    \param  low_bit - low bit (inclusive) allowed for this signal vector
    
    \ingroup svalgo
    \sa dynamic_range_clip_high 
*/
template<typename SV>
void dynamic_range_clip_low(SV& svect, unsigned low_bit)
{
    if (low_bit == 0) return; // nothing to do
    BM_ASSERT(svect.plains() > low_bit);
    
    unsigned sv_plains = svect.plains();
    typename SV::bvector_type bv_acc1;
    unsigned i;
    
    // combine all the high bits into accumulator vector
    for (i = low_bit+1; i < sv_plains; ++i)
    {
        typename SV::bvector_type* bv_plain = svect.plain(i);
        if (bv_plain)
            bv_acc1.bit_or(*bv_plain);
    } // for i
    
    // accumulate all vectors below the clipping point
    typename SV::bvector_type bv_acc2;
    typename SV::bvector_type* bv_low_plain = svect.get_plain(low_bit);
    
    for (i = low_bit-1; true; --i)
    {
        typename SV::bvector_type* bv_plain = svect.plain(i);
        if (bv_plain)
        {
            bv_acc2.bit_or(*bv_plain);
            svect.free_plain(i);
            if (i == 0)
                break;
        }
    } // for i
    
    // now we want to set 1 in the clipping low plain based on
    // exclusive or (XOR) between upper and lower parts)
    // as a result high signal (any bits in the upper plains) gets
    // slightly lower value (due to clipping) low signal gets amplified
    // (lower contrast algorithm)
    
    bv_acc1.bit_xor(bv_acc2);
    bv_low_plain->bit_or(bv_acc1);
}


/**
    \brief algorithms for sparse_vector scan/seach
 
    Scanner uses properties of bit-vector plains to answer questions
    like "find all sparse vector elements equivalent to XYZ".

    Class uses fast algorithms based on properties of bit-plains.
    This is NOT a brute force, direct scan.
 
    @ingroup svalgo
    @ingroup setalgo
*/
template<typename SV>
class sparse_vector_scanner
{
public:
    typedef typename SV::bvector_type       bvector_type;
    typedef const bvector_type*             bvector_type_const_ptr;
    typedef bvector_type*                   bvector_type_ptr;
    typedef typename SV::value_type         value_type;
    typedef typename SV::size_type          size_type;
    typedef typename bvector_type::allocator_type::allocator_pool_type allocator_pool_type;
    
public:
    sparse_vector_scanner();

    /**
        \brief bind sparse vector for all searches
        
        \param sv - input sparse vector to bind for searches
        \param sorted - source index is sorted, build index for binary search
    */
    void bind(const SV& sv, bool sorted);

    /**
        \brief reset sparse vector binding
    */
    void reset_binding();

    /**
        \brief find all sparse vector elements EQ to search value

        Find all sparse vector elements equivalent to specified value

        \param sv - input sparse vector
        \param value - value to search for
        \param bv_out - output bit-vector (search result masks 1 elements)
    */
    void find_eq(const SV&                  sv,
                 typename SV::value_type    value,
                 typename SV::bvector_type& bv_out);

    /**
        \brief find first sparse vector element

        Find all sparse vector elements equivalent to specified value.
        Works well if sperse vector represents unordered set

        \param sv - input sparse vector
        \param value - value to search for
        \param pos - output found sparse vector element index
     
        \return true if found
    */
    bool find_eq(const SV&                  sv,
                 typename SV::value_type    value,
                 typename SV::size_type&    pos);
    
    /**
        \brief find first sparse vector element (string)
    */
    bool find_eq_str(const SV&                      sv,
                     const typename SV::value_type* str,
                     typename SV::size_type&        pos);

    /**
        \brief binary find first sparse vector element (string)
        Sparse vector must be attached (bind())
        @sa bind
    */
    bool find_eq_str(const typename SV::value_type* str,
                     typename SV::size_type&        pos);

    /**
        \brief binary find first sparse vector element (string)     
        Sparse vector must be sorted.
    */
    bool bfind_eq_str(const SV&                      sv,
                     const typename SV::value_type* str,
                     typename SV::size_type&        pos);

    /**
        \brief binary find first sparse vector element (string)
        Sparse vector must be sorted and attached
        @sa bind
    */
    bool bfind_eq_str(const typename SV::value_type* str,
                      typename SV::size_type&        pos);

    /**
        \brief find all sparse vector elements EQ to 0
        \param sv - input sparse vector
        \param bv_out - output bit-vector (search result masks 1 elements)
    */
    void find_zero(const SV&                  sv,
                   typename SV::bvector_type& bv_out);

    /*!
        \brief Find non-zero elements
        Output vector is computed as a logical OR (join) of all plains

        \param  sv - input sparse vector
        \param  bv_out - output bit-bector of non-zero elements
    */
    void find_nonzero(const SV& sv, typename SV::bvector_type& bv_out);

    /**
        \brief invert search result ("EQ" to "not EQ")

        \param  sv - input sparse vector
        \param  bv_out - output bit-bector of non-zero elements
    */
    void invert(const SV& sv, typename SV::bvector_type& bv_out);

    /**
        \brief find all values A IN (C, D, E, F)
        \param  sv - input sparse vector
        \param  start - start iterator (set to search) 
        \param  end   - end iterator (set to search)
        \param  bv_out - output bit-bector of non-zero elements
     */
    template<typename It>
    void find_eq(const SV&  sv,
                 It    start, 
                 It    end,
                 typename SV::bvector_type& bv_out)
    {
        typename bvector_type::mem_pool_guard mp_guard;
        mp_guard.assign_if_not_set(pool_, bv_out); // set algorithm-local memory pool to avoid heap contention

        bvector_type bv1;
        typename bvector_type::mem_pool_guard mp_guard1(pool_, bv1);
        bool any_zero = false;
        for (; start < end; ++start)
        {
            value_type v = *start;
            any_zero |= (v == 0);
            bool found = find_eq_with_nulls(sv, v, bv1);
            if (found)
                bv_out.bit_or(bv1);
            else
            {
                BM_ASSERT(!bv1.any());
            }
        } // for
        if (any_zero)
            correct_nulls(sv, bv_out);
    }

    /// For testing purposes only
    ///
    /// @internal
    void find_eq_with_nulls_horizontal(const SV&   sv,
                 typename SV::value_type           value,
                 typename SV::bvector_type&        bv_out);

    /** Exclude possible NULL values from the result vector
        \param  sv - input sparse vector
        \param  bv_out - output bit-bector of non-zero elements
    */
    void correct_nulls(const SV&   sv, typename SV::bvector_type& bv_out);
    
protected:

    /// set search boundaries (hint for the aggregator)
    void set_search_range(size_type from, size_type to);
    
    /// reset (disable) search range
    void reset_search_range();

    /// find value (may include NULL indexes)
    bool find_eq_with_nulls(const SV&   sv,
                            typename SV::value_type         value,
                            typename SV::bvector_type&      bv_out,
                            typename SV::size_type          search_limit = 0);

    /// find first value (may include NULL indexes)
    bool find_first_eq(const SV&   sv,
                       typename SV::value_type         value,
                       bm::id_t&                       idx);
    
    /// find first string value (may include NULL indexes)
    bool find_first_eq(const SV&                       sv,
                       const typename SV::value_type*  str,
                       bm::id_t&                       idx,
                       bool                            remaped);

    
    /// Prepare aggregator for AND-SUB (EQ) search
    bool prepare_and_sub_aggregator(const SV&   sv,
                                    typename SV::value_type   value);

    /// Prepare aggregator for AND-SUB (EQ) search (string)
    bool prepare_and_sub_aggregator(const SV&  sv,
                                    const typename SV::value_type*  str,
                                    unsigned octet_start);

    /// Rank-Select decompression for RSC vectors
    void decompress(const SV&   sv, typename SV::bvector_type& bv_out);

    int compare_str(const SV& sv, size_type idx, const value_type* str);
    
protected:
    sparse_vector_scanner(const sparse_vector_scanner&) = delete;
    void operator=(const sparse_vector_scanner&) = delete;

protected:

    enum vector_capacity
    {
        max_columns = SV::max_vector_size
    };
    
    typedef bm::heap_matrix<value_type,
        bm::set_total_blocks,
        max_columns,
        typename bvector_type::allocator_type>  heap_matrix_type0;
    typedef bm::heap_matrix<value_type,
        bm::set_total_blocks*3,
        max_columns,
        typename bvector_type::allocator_type>  heap_matrix_type3;

private:
    allocator_pool_type                pool_;
    bvector_type                       bv_tmp_;
    bm::aggregator<bvector_type>       agg_;
    bm::rank_compressor<bvector_type>  rank_compr_;
    
    size_type                          mask_from_;
    size_type                          mask_to_;
    bool                               mask_set_;
    
    const SV*                          bound_sv_;
    heap_matrix_type0                  block0_elements_cache_; ///< cache for elements[0] of each block
    heap_matrix_type3                  block3_elements_cache_; ///< cache for elements[16384x] of each block
    size_type                          effective_str_max_;
    
    value_type                         remap_value_vect_[SV::max_vector_size];
    /// masks of allocated bit-plains (1 - means there is a bit-plain)
    bm::id64_t                         vector_plain_masks_[SV::max_vector_size];
};


/*!
    \brief Integer set to set transformation (functional image in groups theory)
    https://en.wikipedia.org/wiki/Image_(mathematics)
 
    Input sets gets translated through the function, which is defined as
    "one to one (or NULL)" binary relation object (sparse_vector).
    Also works for M:1 relationships.
 
    \ingroup svalgo
    \ingroup setalgo
*/
template<typename SV>
class set2set_11_transform
{
public:
    typedef typename SV::bvector_type       bvector_type;
    typedef typename SV::value_type         value_type;
    typedef typename SV::size_type          size_type;
    typedef typename bvector_type::allocator_type::allocator_pool_type allocator_pool_type;
public:

    set2set_11_transform();
    ~set2set_11_transform();


    /**  Get read access to zero-elements vector
         Zero vector gets populated after attach_sv() is called
         or as a side-effect of remap() with immediate sv param
    */
    const bvector_type& get_bv_zero() const { return bv_zero_; }

    /** Perform remapping (Image function) (based on attached translation table)

    \param bv_in  - input set, defined as a bit-vector
    \param bv_out - output set as a bit-vector

    @sa attach_sv
    */
    void remap(const bvector_type&        bv_in,
               bvector_type&              bv_out);

    /** Perform remapping (Image function)
   
     \param bv_in  - input set, defined as a bit-vector
     \param sv_brel   - binary relation (translation table) sparse vector
     \param bv_out - output set as a bit-vector
    */
    void remap(const bvector_type&        bv_in,
               const    SV&               sv_brel,
               bvector_type&              bv_out);
    
    /** Remap single set element
   
     \param id_from  - input value
     \param sv_brel  - translation table sparse vector
     \param id_to    - out value
     
     \return - true if value was found and remapped
    */
    bool remap(size_type id_from, const SV& sv_brel, size_type& id_to);


    /** Run remap transformation
   
     \param bv_in  - input set, defined as a bit-vector
     \param sv_brel   - translation table sparse vector
     \param bv_out - output set as a bit-vector
     
     @sa remap
    */
    void run(const bvector_type&        bv_in,
             const    SV&               sv_brel,
             bvector_type&              bv_out)
    {
        remap(bv_in, sv_brel, bv_out);
    }
    
    /** Attach a translation table vector for remapping (Image function)

        \param sv_brel   - binary relation sparse vector pointer
                          (pass NULL to detach)
        \param compute_stats - flag to perform computation of some statistics
                               later used in remapping. This only make sense
                               for series of remappings on the same translation
                               vector.
        @sa remap
    */
    void attach_sv(const SV* sv_brel, bool compute_stats = false);

    
protected:
    void one_pass_run(const bvector_type&        bv_in,
                      const    SV&               sv_brel,
                      bvector_type&              bv_out);
    
    /// @internal
    template<unsigned BSIZE>
    struct gather_buffer
    {
        size_type   BM_VECT_ALIGN gather_idx_[BSIZE] BM_VECT_ALIGN_ATTR;
        value_type  BM_VECT_ALIGN buffer_[BSIZE] BM_VECT_ALIGN_ATTR;
    };
    
    enum gather_window_size
    {
        sv_g_size = 1024 * 8
    };
    typedef gather_buffer<sv_g_size>  gather_buffer_type;
    

protected:
    set2set_11_transform(const set2set_11_transform&) = delete;
    void operator=(const set2set_11_transform&) = delete;
    
protected:
    const SV*              sv_ptr_;    ///< current translation table vector
    gather_buffer_type*    gb_;        ///< intermediate buffers
    bvector_type           bv_product_;///< temp vector
    
    bool                   have_stats_; ///< flag of statistics presense
    bvector_type           bv_zero_;   ///< bit-vector for zero elements
    
    allocator_pool_type    pool_;
};



//----------------------------------------------------------------------------
//
//----------------------------------------------------------------------------

template<typename SV>
set2set_11_transform<SV>::set2set_11_transform()
: sv_ptr_(0), gb_(0), have_stats_(false)
{
    gb_ = (gather_buffer_type*)::malloc(sizeof(gather_buffer_type));
    if (!gb_)
    {
        SV::throw_bad_alloc();
    }
}

//----------------------------------------------------------------------------

template<typename SV>
set2set_11_transform<SV>::~set2set_11_transform()
{
    if (gb_)
        ::free(gb_);
}


//----------------------------------------------------------------------------

template<typename SV>
void set2set_11_transform<SV>::attach_sv(const SV*  sv_brel, bool compute_stats)
{
    sv_ptr_ = sv_brel;
    if (!sv_ptr_)
    {
        have_stats_ = false;
    }
    else
    {
        if (sv_brel->empty() || !compute_stats)
            return; // nothing to do
        const bvector_type* bv_non_null = sv_brel->get_null_bvector();
        if (bv_non_null)
            return; // already have 0 elements vector
        

        typename bvector_type::mem_pool_guard mp_g_z;
        mp_g_z.assign_if_not_set(pool_, bv_zero_);

        bm::sparse_vector_scanner<SV> scanner;
        scanner.find_zero(*sv_brel, bv_zero_);
        have_stats_ = true;
    }
}

//----------------------------------------------------------------------------

template<typename SV>
bool set2set_11_transform<SV>::remap(size_type  id_from,
                                     const SV&  sv_brel,
                                     size_type& id_to)
{
    if (sv_brel.empty())
        return false; // nothing to do

    const bvector_type* bv_non_null = sv_brel.get_null_bvector();
    if (bv_non_null)
    {
        if (!bv_non_null->test(id_from))
            return false;
    }
    size_type idx = sv_brel.translate_address(id_from);
    id_to = sv_brel.get(idx);
    return true;
}

//----------------------------------------------------------------------------

template<typename SV>
void set2set_11_transform<SV>::remap(const bvector_type&        bv_in,
                                     const    SV&               sv_brel,
                                     bvector_type&              bv_out)
{
    typename bvector_type::mem_pool_guard mp_g_out, mp_g_p, mp_g_z;
    mp_g_out.assign_if_not_set(pool_, bv_out);
    mp_g_p.assign_if_not_set(pool_, bv_product_);
    mp_g_z.assign_if_not_set(pool_, bv_zero_);

    attach_sv(&sv_brel);
    
    remap(bv_in, bv_out);

    attach_sv(0); // detach translation table
}

template<typename SV>
void set2set_11_transform<SV>::remap(const bvector_type&        bv_in,
                                     bvector_type&              bv_out)
{
    BM_ASSERT(sv_ptr_);

    bv_out.clear();

    if (sv_ptr_ == 0 || sv_ptr_->empty())
        return; // nothing to do

    bv_out.init(); // just in case to "fast set" later

    typename bvector_type::mem_pool_guard mp_g_out, mp_g_p, mp_g_z;
    mp_g_out.assign_if_not_set(pool_, bv_out);
    mp_g_p.assign_if_not_set(pool_, bv_product_);
    mp_g_z.assign_if_not_set(pool_, bv_zero_);


    const bvector_type* enum_bv;

    const bvector_type * bv_non_null = sv_ptr_->get_null_bvector();
    if (bv_non_null)
    {
        // TODO: optimize with 2-way ops
        bv_product_ = bv_in;
        bv_product_.bit_and(*bv_non_null);
        enum_bv = &bv_product_;
    }
    else // non-NULL vector is not available
    {
        enum_bv = &bv_in;
        // if we have any elements mapping into "0" on the other end
        // we map it once (chances are there are many duplicates)
        //
        
        if (have_stats_) // pre-attached translation statistics
        {
            bv_product_ = bv_in;
            unsigned cnt1 = bv_product_.count();
            bv_product_.bit_sub(bv_zero_);
            unsigned cnt2 = bv_product_.count();
            
            BM_ASSERT(cnt2 <= cnt1);
            
            if (cnt1 != cnt2) // mapping included 0 elements
                bv_out.set_bit_no_check(0);
            
            enum_bv = &bv_product_;
        }
    }

    

    unsigned buf_cnt, nb_old, nb;
    buf_cnt = nb_old = 0;
    
    typename bvector_type::enumerator en(enum_bv->first());
    for (; en.valid(); ++en)
    {
        typename SV::size_type idx = *en;
        idx = sv_ptr_->translate_address(idx);
        
        nb = unsigned(idx >> bm::set_block_shift);
        if (nb != nb_old) // new blocks
        {
            if (buf_cnt)
            {
                sv_ptr_->gather(&gb_->buffer_[0], &gb_->gather_idx_[0], buf_cnt, BM_SORTED_UNIFORM);
                bv_out.set(&gb_->buffer_[0], buf_cnt, BM_SORTED);
                buf_cnt ^= buf_cnt;
            }
            nb_old = nb;
            gb_->gather_idx_[buf_cnt++] = idx;
        }
        else
        {
            gb_->gather_idx_[buf_cnt++] = idx;
        }
        
        if (buf_cnt == sv_g_size)
        {
            sv_ptr_->gather(&gb_->buffer_[0], &gb_->gather_idx_[0], buf_cnt, BM_SORTED_UNIFORM);
            bv_out.set(&gb_->buffer_[0], buf_cnt, bm::BM_SORTED);
            buf_cnt ^= buf_cnt;
        }
    } // for en
    if (buf_cnt)
    {
        sv_ptr_->gather(&gb_->buffer_[0], &gb_->gather_idx_[0], buf_cnt, BM_SORTED_UNIFORM);
        bv_out.set(&gb_->buffer_[0], buf_cnt, bm::BM_SORTED);
    }
}


//----------------------------------------------------------------------------

template<typename SV>
void set2set_11_transform<SV>::one_pass_run(const bvector_type&        bv_in,
                                            const    SV&               sv_brel,
                                            bvector_type&              bv_out)
{
    if (sv_brel.empty())
        return; // nothing to do

    bv_out.init();

    typename SV::const_iterator it = sv_brel.begin();
    for (; it.valid(); ++it)
    {
        typename SV::value_type t_id = *it;
        bm::id_t idx = it.pos();
        if (bv_in.test(idx))
        {
            bv_out.set_bit_no_check(t_id);
        }
    } // for
}


//----------------------------------------------------------------------------
//
//----------------------------------------------------------------------------

template<typename SV>
sparse_vector_scanner<SV>::sparse_vector_scanner()
{
    mask_set_ = false;
    mask_from_ = mask_to_ = bm::id_max;

    bound_sv_ = 0;
    effective_str_max_ = 0;
}

//----------------------------------------------------------------------------

template<typename SV>
void sparse_vector_scanner<SV>::bind(const SV&  sv, bool sorted)
{
    bound_sv_ = &sv;
    if (sorted)
    {
        block0_elements_cache_.init();
        block0_elements_cache_.set_zero();
        block3_elements_cache_.init();
        block3_elements_cache_.set_zero();
        
        effective_str_max_ = sv.effective_vector_max();
        // fill in elements cache
        for (size_type i = 0; i < sv.size(); i+= bm::gap_max_bits)
        {
            unsigned nb = unsigned(i >> bm::set_block_shift);
            value_type* s0 = block0_elements_cache_.row(nb);
            sv.get(i, s0, block0_elements_cache_.cols());
            
            for (size_type k = 0; k < 3; ++k)
            {
                value_type* s1 = block3_elements_cache_.row(nb * 3 + k);
                size_type idx = i + (k+1) * bm::sub_block3_size;
                sv.get(idx, s1, block3_elements_cache_.cols());
            } // for k
        } // for i
    }
    // pre-calculate vector plain masks
    //
    for (unsigned i = 0; i < SV::max_vector_size; ++i)
    {
        vector_plain_masks_[i] = sv.get_plains_mask(i);
    } // for i
    
}

//----------------------------------------------------------------------------

template<typename SV>
void sparse_vector_scanner<SV>::reset_binding()
{
    bound_sv_ = 0;
    effective_str_max_ = 0;
}

//----------------------------------------------------------------------------

template<typename SV>
void sparse_vector_scanner<SV>::find_zero(const SV&                  sv,
                                          typename SV::bvector_type& bv_out)
{
    if (sv.size() == 0)
    {
        bv_out.clear();
        return;
    }

    find_nonzero(sv, bv_out);
    if (sv.is_compressed())
    {
        bv_out.invert();
        bv_out.set_range(sv.effective_size(), bm::id_max - 1, false);
        decompress(sv, bv_out);
    }
    else
    {
        invert(sv, bv_out);
    }
    correct_nulls(sv, bv_out);
}

//----------------------------------------------------------------------------

template<typename SV>
void sparse_vector_scanner<SV>::invert(const SV& sv, typename SV::bvector_type& bv_out)
{
    if (sv.size() == 0)
    {
        bv_out.clear();
        return;
    }
    bv_out.invert();
    const bvector_type* bv_null = sv.get_null_bvector();
    if (bv_null) // correct result to only use not NULL elements
        bv_out &= *bv_null;
    else
        bv_out.set_range(sv.size(), bm::id_max - 1, false);
}

//----------------------------------------------------------------------------

template<typename SV>
void sparse_vector_scanner<SV>::correct_nulls(const SV&   sv,
                           typename SV::bvector_type& bv_out)
{
    const bvector_type* bv_null = sv.get_null_bvector();
    if (bv_null) // correct result to only use not NULL elements
        bv_out.bit_and(*bv_null);
}

//----------------------------------------------------------------------------

template<typename SV>
bool sparse_vector_scanner<SV>::find_eq_with_nulls(const SV&    sv,
                                    typename SV::value_type     value,
                                    typename SV::bvector_type&  bv_out,
                                    typename SV::size_type      search_limit)
{
    if (sv.empty())
        return false; // nothing to do

    if (!value)
    {
        find_zero(sv, bv_out);
        return bv_out.any();
    }
    agg_.reset();
    
    bool found = prepare_and_sub_aggregator(sv, value);
    if (!found)
    {
        bv_out.clear();
        return found;
    }

    bool any = (search_limit == 1);
    found = agg_.combine_and_sub(bv_out, any);
    agg_.reset();
    return found;
}

//----------------------------------------------------------------------------

template<typename SV>
bool sparse_vector_scanner<SV>::find_first_eq(const SV&   sv,
                               typename SV::value_type    value,
                               bm::id_t&                  idx)
{
    if (sv.empty())
        return false; // nothing to do
    
    BM_ASSERT(value); // cannot handle zero values
    if (!value)
        return false;

    agg_.reset();
    bool found = prepare_and_sub_aggregator(sv, value);
    if (!found)
        return found;
    found = agg_.find_first_and_sub(idx);
    agg_.reset();
    return found;
}


//----------------------------------------------------------------------------

template<typename SV>
bool sparse_vector_scanner<SV>::find_first_eq(const SV&                   sv,
                                          const typename SV::value_type*  str,
                                          bm::id_t&                       idx,
                                          bool                            remaped)
{
    if (sv.empty())
        return false; // nothing to do
    BM_ASSERT(*str);

    if (!*str)
        return false;

    agg_.reset();
    unsigned common_prefix_len = 0;
    
    if (mask_set_)
    {
        agg_.set_range_hint(mask_from_, mask_to_);
        common_prefix_len = sv.common_prefix_length(mask_from_, mask_to_);
    }
    
    if (remaped)
    {
        str = remap_value_vect_;
    }
    else
    {
        if (sv.is_remap() && str != remap_value_vect_)
        {
            bool r = sv.remap_tosv(remap_value_vect_, SV::max_vector_size, str);
            if (!r)
                return r;
            str = remap_value_vect_;
        }
    }
    
    bool found = prepare_and_sub_aggregator(sv, str, common_prefix_len);
    if (!found)
        return found;
    
    found = agg_.find_first_and_sub(idx);
    agg_.reset();
    return found;
}

//----------------------------------------------------------------------------

template<typename SV>
bool sparse_vector_scanner<SV>::prepare_and_sub_aggregator(const SV&  sv,
                                      const typename SV::value_type*  str,
                                      unsigned octet_start)
{
    unsigned char bits[64];

    int len = 0;
    for (; str[len] != 0; ++len)
    {}
    BM_ASSERT(len);

    // use reverse order (faster for sorted arrays)
    for (int octet_idx = len-1; octet_idx >= 0; --octet_idx)
    {
        if (unsigned(octet_idx) < octet_start) // common prefix
            break;

        unsigned value = unsigned(str[octet_idx]) & 0xFF;
        BM_ASSERT(value != 0);
        
        bm::id64_t plains_mask;
        if (&sv == bound_sv_)
            plains_mask = vector_plain_masks_[octet_idx];
        else
            plains_mask = sv.get_plains_mask(unsigned(octet_idx));

        if ((value & plains_mask) != value) // pre-screen for impossible values
            return false; // found non-existing plain

        // prep the lists for combined AND-SUB aggregator
        //
        unsigned short bit_count_v = bm::bitscan(value, bits);
        for (unsigned i = 0; i < bit_count_v; ++i)
        {
            unsigned bit_idx = bits[i];
            BM_ASSERT(value & (value_type(1) << bit_idx));
            unsigned plain_idx = (unsigned(octet_idx) * 8) + bit_idx;
            const bvector_type* bv = sv.get_plain(plain_idx);
            agg_.add(bv);
        } // for i
        
        unsigned vinv = unsigned(value);
        if (bm::conditional<sizeof(value_type) == 1>::test())
        {
            vinv ^= 0xFF;
        }
        else // 2-byte char
        {
            BM_ASSERT(sizeof(value_type) == 2);
            vinv ^= 0xFFFF;
        }
        // exclude the octet bits which are all not set (no vectors)
        vinv &= unsigned(plains_mask);
        for (unsigned octet_plain = (unsigned(octet_idx) * 8); vinv; vinv &= vinv-1)
        {
            bm::id_t t = vinv & -vinv;
            unsigned bit_idx = bm::word_bitcount(t - 1);
            unsigned plain_idx = octet_plain + bit_idx;
            const bvector_type* bv = sv.get_plain(plain_idx);
            BM_ASSERT(bv);
            agg_.add(bv, 1); // agg to SUB group
        } // for
     } // for octet_idx
    
    // add all vectors above string len to the SUB operation group
    //
    typename SV::size_type plain_idx = unsigned(len * 8) + 1;
    typename SV::size_type plains;
    if (&sv == bound_sv_)
        plains = effective_str_max_ * unsigned(sizeof(value_type)) * 8;
    else
        plains = sv.plains();
    for (; plain_idx < plains; ++plain_idx)
    {
        bvector_type_const_ptr bv = sv.get_plain(plain_idx);
        if (bv)
            agg_.add(bv, 1); // agg to SUB group
    } // for
    return true;
}

//----------------------------------------------------------------------------

template<typename SV>
bool sparse_vector_scanner<SV>::prepare_and_sub_aggregator(const SV&   sv,
                                           typename SV::value_type   value)
{
    unsigned char bits[sizeof(value) * 8];
    unsigned short bit_count_v = bm::bitscan(value, bits);
    BM_ASSERT(bit_count_v);

    // prep the lists for combined AND-SUB aggregator
    //   (backward order has better chance for bit reduction on AND)
    //
    for (unsigned i = bit_count_v; i > 0; --i)
    {
        unsigned bit_idx = bits[i-1];
        BM_ASSERT(value & (value_type(1) << bit_idx));
        const bvector_type* bv = sv.get_plain(bit_idx);
        if (bv)
            agg_.add(bv);
        else
            return false;
    }
    
    unsigned sv_plains = sv.effective_plains();
    for (unsigned i = 0; (i < sv_plains) && value; ++i)
    {
        bvector_type_const_ptr bv = sv.get_plain(i);
        if (bv && !(value & (value_type(1) << i)))
            agg_.add(bv, 1); // agg to SUB group
    } // for i
    return true;
}


//----------------------------------------------------------------------------

template<typename SV>
void sparse_vector_scanner<SV>::find_eq_with_nulls_horizontal(const SV&  sv,
    typename SV::value_type    value,
    typename SV::bvector_type& bv_out)
{
    if (sv.empty())
        return; // nothing to do

    if (!value)
    {
        find_zero(sv, bv_out);
        return;
    }

    unsigned char bits[sizeof(value) * 8];
    unsigned short bit_count_v = bm::bitscan(value, bits);
    BM_ASSERT(bit_count_v);

    // aggregate AND all matching vectors
    {
        const bvector_type* bv_plain = sv.get_plain(bits[--bit_count_v]);
        if (bv_plain)
            bv_out = *bv_plain;
        else // plain not found
        {
            bv_out.clear();
            return;
        }
    }
    for (unsigned i = 0; i < bit_count_v; ++i)
    {
        const bvector_type* bv_plain = sv.get_plain(bits[i]);
        if (bv_plain)
            bv_out &= *bv_plain;
        else // mandatory plain not found - empty result!
        {
            bv_out.clear();
            return;
        }
    } // for i

    // SUB all other plains
    unsigned sv_plains = sv.effective_plains();
    for (unsigned i = 0; (i < sv_plains) && value; ++i)
    {
        const bvector_type* bv_plain = sv.get_plain(i);
        if (bv_plain && !(value & (value_type(1) << i)))
            bv_out -= *bv_plain;
    }
}

//----------------------------------------------------------------------------

template<typename SV>
bool sparse_vector_scanner<SV>::find_eq_str(const typename SV::value_type* str,
                                            typename SV::size_type&        pos)
{
    BM_ASSERT(bound_sv_);
    return find_eq_str(*bound_sv_, str, pos);
}

//----------------------------------------------------------------------------

template<typename SV>
bool sparse_vector_scanner<SV>::find_eq_str(const SV&                      sv,
                                            const typename SV::value_type* str,
                                            typename SV::size_type&        pos)
{
    bool found = false;
    if (sv.empty())
        return found;
    if (*str)
    {
        bool remaped = false;
        if (bm::conditional<SV::is_remap_support::value>::test()) // test remapping trait
        {
            if (sv.is_remap() && str != remap_value_vect_)
            {
                remaped = sv.remap_tosv(remap_value_vect_, SV::max_vector_size, str);
                if (!remaped)
                    return remaped;
                str = remap_value_vect_;
            }
        }
    
        bm::id_t found_pos;
        found = find_first_eq(sv, str, found_pos, remaped);
        if (found)
        {
            pos = found_pos;
            if (bm::conditional<SV::is_rsc_support::value>::test()) // test rank/select trait
            {
                if (sv.is_compressed()) // if compressed vector - need rank translation
                    found = sv.find_rank(found_pos + 1, pos);
            }
        }
    }
    else // search for zero value
    {
        // TODO: implement optimized version which would work witout temp vector
        typename SV::bvector_type bv_zero;
        find_zero(sv, bv_zero);
        found = bv_zero.find(pos);
    }
    return found;
}

//----------------------------------------------------------------------------

template<typename SV>
bool sparse_vector_scanner<SV>::bfind_eq_str(const SV&                      sv,
                                             const typename SV::value_type* str,
                                             typename SV::size_type&        pos)
{
    bool found = false;
    if (sv.empty())
        return found;

    if (*str)
    {
        bool remaped = false;
        // test search pre-condition based on remap tables
        if (bm::conditional<SV::is_remap_support::value>::test())
        {
            if (sv.is_remap() && str != remap_value_vect_)
            {
                remaped = sv.remap_tosv(
                                remap_value_vect_, SV::max_vector_size, str);
                if (!remaped)
                    return remaped;
            }
        }
        
        reset_search_range();
        
        // narrow down the search
        const unsigned min_distance_cutoff = bm::gap_max_bits + bm::gap_max_bits / 2;
        size_type l, r, dist;
        l = 0; r = sv.size()-1;
        bm::id_t found_pos;
        
        // binary search to narrow down the search window
        while (l <= r)
        {
            dist = r - l;
            if (dist < min_distance_cutoff)
            {
                // we are in an narrow <2 blocks window, but still may be in two
                // different neighboring blocks, lets try to narrow
                // to exactly one block
                
                unsigned nb_l = unsigned(l >> bm::set_block_shift);
                unsigned nb_r = unsigned(r >> bm::set_block_shift);
                if (nb_l != nb_r)
                {
                    size_type mid = nb_r * bm::gap_max_bits;
                    if (mid < r)
                    {
                        int cmp = this->compare_str(sv, mid, str);
                        if (cmp < 0) // mid < str
                            l = mid;
                        else
                            r = mid-(cmp!=0); // branchless if (cmp==0) r=mid;
                        BM_ASSERT(l < r);
                    }
                    nb_l = unsigned(l >> bm::set_block_shift);
                    nb_r = unsigned(r >> bm::set_block_shift);
                }
                
                if (nb_l == nb_r)
                {
                    unsigned max_nb = sv.size() >> bm::set_block_shift;
                    if (nb_l != max_nb)
                    {
                        // linear in-place fixed depth scan to identify the sub-range
                        size_type mid = nb_r * bm::gap_max_bits + bm::sub_block3_size;
                        int cmp = this->compare_str(sv, mid, str);
                        if (cmp < 0)
                        {
                            l = mid;
                            mid = nb_r * bm::gap_max_bits + bm::sub_block3_size * 2;
                            cmp = this->compare_str(sv, mid, str);
                            if (cmp < 0)
                            {
                                l = mid;
                                mid = nb_r * bm::gap_max_bits + bm::sub_block3_size * 3;
                                cmp = this->compare_str(sv, mid, str);
                                if (cmp < 0)
                                    l = mid;
                                else
                                    r = mid;
                            }
                            else
                            {
                                r = mid;
                            }
                        }
                        else
                        {
                            r = mid;
                        }
                    }
                }
                
                set_search_range(l, r);
                break;
            }

            typename SV::size_type mid = dist/2+l;
            unsigned nb = unsigned(mid >> bm::set_block_shift);
            mid = nb * bm::gap_max_bits;
            if (mid <= l)
            {
                if (nb == 0 && r > bm::gap_max_bits)
                    mid = bm::gap_max_bits;
                else
                    mid = dist / 2 + l;
            }
            BM_ASSERT(mid > l);

            int cmp = this->compare_str(sv, mid, str);
            if (cmp == 0)
            {
                found_pos = mid;
                found = true;
                set_search_range(l, mid);
                break;
            }
            if (cmp < 0)
                l = mid+1;
            else
                r = mid-1;
        } // while

        // use linear search (range is set)
        found = find_first_eq(sv, str, found_pos, remaped);
        if (found)
        {
            pos = found_pos;
            if (bm::conditional<SV::is_rsc_support::value>::test()) // test rank/select trait
            {
                if (sv.is_compressed()) // if compressed vector - need rank translation
                    found = sv.find_rank(found_pos + 1, pos);
            }
        }
        reset_search_range();
    }
    else // search for zero value
    {
        // TODO: implement optimized version which would work without temp vector
        typename SV::bvector_type bv_zero;
        find_zero(sv, bv_zero);
        found = bv_zero.find(pos);
    }
    return found;
}

//----------------------------------------------------------------------------

template<typename SV>
bool sparse_vector_scanner<SV>::bfind_eq_str(const typename SV::value_type* str,
                                             typename SV::size_type&        pos)
{
    BM_ASSERT(bound_sv_);
    return bfind_eq_str(*bound_sv_, str, pos);
}

//----------------------------------------------------------------------------

template<typename SV>
int sparse_vector_scanner<SV>::compare_str(const SV& sv,
                                           size_type idx,
                                           const value_type* str)
{
    if (bound_sv_ == &sv)
    {
        unsigned nb = unsigned(idx >> bm::set_block_shift);
        unsigned nbit = unsigned(idx & bm::set_block_mask);
        if (nbit == 0) // access to sentinel, first block element
        {
            value_type* s0 = block0_elements_cache_.row(nb);
            if (*s0 == 0) // uninitialized element
            {
                sv.get(idx, s0, block0_elements_cache_.cols());
            }
            int res = 0;
            for (unsigned i = 0; i < block0_elements_cache_.cols(); ++i)
            {
                char octet = str[i]; char value = s0[i];
                res = (value > octet) - (value < octet);
                if (res || !octet)
                    break;
            } // for i
            BM_ASSERT(res == sv.compare(idx, str));
            return res;
        }
        else
        {
            if (nbit % bm::sub_block3_size == 0) // TODO: use AND mask here
            {
                size_type k = nbit / bm::sub_block3_size - 1;
                value_type* s1 = block3_elements_cache_.row(nb * 3 + k);
                int res = 0;
                for (unsigned i = 0; i < block3_elements_cache_.cols(); ++i)
                {
                    char octet = str[i]; char value = s1[i];
                    res = (value > octet) - (value < octet);
                    if (res || !octet)
                        break;
                } // for i
                BM_ASSERT(res == sv.compare(idx, str));
                return res;
            }
        }
    }
    return sv.compare(idx, str);
}

//----------------------------------------------------------------------------

template<typename SV>
void sparse_vector_scanner<SV>::find_eq(const SV&                  sv,
                                        typename SV::value_type    value,
                                        typename SV::bvector_type& bv_out)
{
    if (sv.empty())
    {
        bv_out.clear();
        return; // nothing to do
    }
    if (!value)
    {
        find_zero(sv, bv_out);
        return;
    }

    find_eq_with_nulls(sv, value, bv_out, 0);
    
    decompress(sv, bv_out);
    correct_nulls(sv, bv_out);
}

//----------------------------------------------------------------------------

template<typename SV>
bool sparse_vector_scanner<SV>::find_eq(const SV&                  sv,
                                        typename SV::value_type    value,
                                        typename SV::size_type&    pos)
{
    if (!value) // zero value - special case
    {
        bvector_type bv_zero;
        find_eq(sv, value, bv_zero);
        bool found = bv_zero.find(pos);
        return found;
    }

    bm::id_t found_pos;
    bool found = find_first_eq(sv, value, found_pos);
    if (found)
    {
        pos = found_pos;
        if (bm::conditional<SV::is_rsc_support::value>::test()) // test rank/select trait
        {
            if (sv.is_compressed()) // if compressed vector - need rank translation
                found = sv.find_rank(found_pos + 1, pos);
        }
    }
    return found;
}

//----------------------------------------------------------------------------

template<typename SV>
void sparse_vector_scanner<SV>::find_nonzero(const SV& sv, 
                                             typename SV::bvector_type& bv_out)
{
    agg_.reset(); // in case if previous scan was interrupted
    for (unsigned i = 0; i < sv.plains(); ++i)
        agg_.add(sv.get_plain(i));
    agg_.combine_or(bv_out);
    agg_.reset();
}

//----------------------------------------------------------------------------

template<typename SV>
void sparse_vector_scanner<SV>::decompress(const SV&   sv,
                                           typename SV::bvector_type& bv_out)
{
    if (!sv.is_compressed())
        return; // nothing to do
    const bvector_type* bv_non_null = sv.get_null_bvector();
    BM_ASSERT(bv_non_null);

    // TODO: implement faster decompressor for small result-sets
    rank_compr_.decompress(bv_tmp_, *bv_non_null, bv_out);
    bv_out.swap(bv_tmp_);
}

//----------------------------------------------------------------------------

template<typename SV>
void sparse_vector_scanner<SV>::set_search_range(size_type from, size_type to)
{
    BM_ASSERT(from < to);
    mask_from_ = from;
    mask_to_ = to;
    mask_set_ = true;
}

//----------------------------------------------------------------------------

template<typename SV>
void sparse_vector_scanner<SV>::reset_search_range()
{
    mask_set_ = false;
    mask_from_ = mask_to_ = bm::id_max;
}


//----------------------------------------------------------------------------
//
//----------------------------------------------------------------------------


} // namespace bm

#include "bmundef.h"

#endif
