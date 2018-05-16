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
public:
    /** Perform transformation
   
     \param bv_in  - input set, defined as a bit-vector
     \param sv_brel   - binary relation sparse vector
     \param bv_out - output set as a bit-vector
    */
    void run(typename SV::bvector_type& bv_in,
             const    SV&               sv_brel,
             typename SV::bvector_type& bv_out)
    {
        if (sv_brel.empty())
            return; // nothing to do
        
        bv_out.init(); // just in case to "fast set" later
        
        const typename SV::bvector_type * bv_non_null = sv_brel.get_null_bvector(); 
        if (bv_non_null) // NULL-able association vector
        {
            bv_product_ = bv_in;
            bv_product_.bit_and(*bv_non_null);
        }
        else
        {
            bv_product_.clear(true);
            bv_product_.set_range(0, sv_brel.size()-1);
            bv_product_.bit_and(bv_in);
        }
        
        typename SV::bvector_type::enumerator en(bv_product_.first());
        for (; en.valid(); ++en)
        {
            auto idx = *en;
            idx = sv_brel.translate_address(idx);
            typename SV::value_type translated_id = sv_brel.get(idx);
            bv_out.set_bit_no_check(translated_id);
        } // for en
    }
    
protected:
    bvector_type   bv_product_;
};

/**
    \brief algorithms for sparse_vector scan/seach
 
    Scanner uses properties of bit-vector plains to answer questions
    like "find all sparse vector elements equivalent to XYZ".

    Class uses fast algorithms based on properties of bit-plains.
    This is NOT a brute force, direct scan.
 
    @ingroup svalgo
*/
template<typename SV>
class sparse_vector_scanner
{
public:
    typedef typename SV::bvector_type       bvector_type;
    typedef typename SV::value_type         value_type;
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

        for (; start < end; ++start)
        {
            value_type v = *start;
            find_eq_with_nulls(sv, v, bv1);
            bv_out.bit_or(bv1);
        } // for
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

template<typename SV>
void sparse_vector_scanner<SV>::find_nonzero(const SV& sv, 
                                             typename SV::bvector_type& bv_out)
{
    bool first = true;
    for (unsigned i = 0; i < sv.plains(); ++i)
    {
        const typename SV::bvector_type* bv_plain = sv.plain(i);
        if (bv_plain)
        {
            if (first) // first found plain - use simple assignment
            {
                bv_out = *bv_plain;
                first = false;
            }
            else // for everything else - use OR
            {
                bv_out.bit_or(*bv_plain);
            }
        }
    } // for i
    if (first) // no plains were found, just clear the result
        bv_out.clear(true);
}

//----------------------------------------------------------------------------
//
//----------------------------------------------------------------------------


} // namespace bm

#include "bmundef.h"

#endif
