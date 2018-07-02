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
    typedef typename SV::value_type         value_type;
    typedef typename SV::size_type          size_type;
    typedef typename bvector_type::allocator_type::allocator_pool_type allocator_pool_type;
    
public:
    sparse_vector_scanner() {}

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
        \param bv_out - output bit-vector (search result masks 1 elements)
    */
    bool find_eq(const SV&                  sv,
                 typename SV::value_type    value,
                 typename SV::size_type&    pos);


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
            find_eq_with_nulls(sv, v, bv1);
            bv_out.bit_or(bv1);
        } // for
        if (any_zero)
            correct_nulls(sv, bv_out);
    }

    void find_eq_with_nulls(const SV&   sv,
        typename SV::value_type           value,
        typename SV::bvector_type&        bv_out);

    void correct_nulls(const SV&   sv,
        typename SV::bvector_type& bv_out);


protected:
    sparse_vector_scanner(const sparse_vector_scanner&) = delete;
    void operator=(const sparse_vector_scanner&) = delete;
private:
    allocator_pool_type  pool_;
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
     \param it_to    - out value
     
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
                bm::combine_or(bv_out, &gb_->buffer_[0], &gb_->buffer_[buf_cnt]);
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
            bm::combine_or(bv_out, &gb_->buffer_[0], &gb_->buffer_[buf_cnt]);
            buf_cnt ^= buf_cnt;
        }
    } // for en
    if (buf_cnt)
    {
        sv_ptr_->gather(&gb_->buffer_[0], &gb_->gather_idx_[0], buf_cnt, BM_SORTED_UNIFORM);
        bm::combine_or(bv_out, &gb_->buffer_[0], &gb_->buffer_[buf_cnt]);
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

    //const typename SV::bvector_type * bv_non_null = sv_brel.get_null_bvector();
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
void sparse_vector_scanner<SV>::find_zero(const SV&                  sv,
                                          typename SV::bvector_type& bv_out)
{
    find_nonzero(sv, bv_out);
    invert(sv, bv_out);
}

//----------------------------------------------------------------------------

template<typename SV>
void sparse_vector_scanner<SV>::invert(const SV& sv, typename SV::bvector_type& bv_out)
{
    if (sv.size() == 0)
        return;
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
void sparse_vector_scanner<SV>::find_eq_with_nulls(const SV&  sv,
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
    //
    {
        const bvector_type* bv_plain = sv.get_plain(bits[--bit_count_v]);
        if (bv_plain)
        {
            bv_out = *bv_plain;
        }
        else // plain not found
        {
            bv_out.clear(true);
            return;
        }
    }
    for (unsigned i = 0; i < bit_count_v; ++i)
    {
        const bvector_type* bv_plain = sv.get_plain(bits[i]);
        if (bv_plain)
        {
            bv_out &= *bv_plain;
            // TODO: better detect when accumulator is empty to break early
        }
        else // mandatory plain not found - empty result!
        {
            bv_out.clear(true);
            return;
        }
    } // for i

    // SUB all other plains
    //
    unsigned sv_plains = sv.effective_plains();
    for (unsigned i = 0; (i < sv_plains) && value; ++i)
    {
        const bvector_type* bv_plain = sv.get_plain(i);
        if (bv_plain && !(value & (value_type(1) << i)))
        {
            // TODO: better detect when result is empty to break early
            bv_out -= *bv_plain;
        }
    } // for i
}

//----------------------------------------------------------------------------

template<typename SV>
void sparse_vector_scanner<SV>::find_eq(const SV&                  sv,
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

    find_eq_with_nulls(sv, value, bv_out);
    correct_nulls(sv, bv_out);
}

//----------------------------------------------------------------------------

template<typename SV>
bool sparse_vector_scanner<SV>::find_eq(const SV&                  sv,
                                        typename SV::value_type    value,
                                        typename SV::size_type&    pos)
{
    bvector_type bv_tmp;
    bv_tmp.set_allocator_pool(&pool_);
    
    find_eq_with_nulls(sv, value, bv_tmp);
    bm::id_t found_pos;
    bool found = bv_tmp.find(found_pos);
    
    if (found)
    {
        if (sv.is_compressed()) // if compressed vector - need rank translation
            found = sv.find_rank(found_pos + 1, pos);
        else
            pos = found_pos;
        
        if (!value && found)
        {
            const bvector_type* bv_null = sv.get_null_bvector();
            if (bv_null) // correct result to only use not NULL elements
                found = bv_null->test(pos);
        }
    }

    return found;
}

//----------------------------------------------------------------------------

template<typename SV>
void sparse_vector_scanner<SV>::find_nonzero(const SV& sv, 
                                             typename SV::bvector_type& bv_out)
{
    bvector_type* scan_list[SV::sv_plains];

    unsigned cnt = 0;
    for (unsigned i = 0; i < sv.plains(); ++i)
    {
        const bvector_type* bv = sv.get_plain(i);
        if (bv)
        {
            const typename bvector_type::blocks_manager_type& bman = bv->get_blocks_manager();
            if (bman.is_init())
                scan_list[cnt++] = sv.get_plain(i);
        }
    }
    if (cnt)
    {
        bm::aggregator<typename SV::bvector_type> agg;
        agg.combine_or(bv_out, scan_list, cnt);
    }
    else
    {
        bv_out.clear();
    }
}

//----------------------------------------------------------------------------
//
//----------------------------------------------------------------------------


} // namespace bm

#include "bmundef.h"

#endif
