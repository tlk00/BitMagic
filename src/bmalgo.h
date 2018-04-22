#ifndef BMALGO__H__INCLUDED__
#define BMALGO__H__INCLUDED__
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

/*! \file bmalgo.h
    \brief Algorithms for bvector<> (main include)
*/

#include "bm.h"
#include "bmfunc.h"
#include "bmdef.h"

#include "bmalgo_impl.h"



namespace bm
{

/**
    @brief bit-vector visitor scanner to traverse each 1 bit using C++ visitor
 
    @param bv - bit vector to scan
    @param bit_functor (should support add_bits() and add_range() methods
 
    \ingroup setalgo
*/
template<class BV, class Func>
void for_each_bit(const BV&    bv,
                  Func&        bit_functor)
{
    const typename BV::blocks_manager_type& bman = bv.get_blocks_manager();
    bm::word_t*** blk_root = bman.top_blocks_root();
    
    unsigned tsize = bman.top_block_size();
    for (unsigned i = 0; i < tsize; ++i)
    {
        bm::word_t** blk_blk = blk_root[i];
        if (!blk_blk)
        {
            continue;
        }
        unsigned r = i * bm::set_array_size;
        for (unsigned j = 0; j < bm::set_array_size; ++j)
        {
            const bm::word_t* block = blk_blk[j];
            if (block)
            {
                unsigned nb = r + j;
                if (BM_IS_GAP(block))
                {
                    bm::for_each_gap_blk(BMGAP_PTR(block),
                                         nb * bm::bits_in_block,
                                         bit_functor);
                }
                else
                {
                    // TODO: optimize FULL BLOCK ADDRESS scan
                    block = BLOCK_ADDR_SAN(block);
                    
                    bm::for_each_bit_blk(block, nb * bm::bits_in_block,
                                         bit_functor);
                }
            }
        } // for j
    }  // for i

}


/**
    @brief bit-vector visitor scanner to traverse each 1 bit using C callback
 
    @param bv - bit vector to scan
    @param handle_ptr - handle to private memory used by callback
    @param callback_ptr - callback function
 
    \ingroup setalgo
 
    @sa bit_visitor_callback_type
*/
template<class BV>
void visit_each_bit(const BV&                 bv,
                    void*                     handle_ptr,
                    bit_visitor_callback_type callback_ptr)
{
    // private adaptor for C-style callbacks
    struct callback_adaptor
    {
        callback_adaptor(void* h, bit_visitor_callback_type cb_func)
        : handle_(h), func_(cb_func)
        {}
        
        void add_bits(bm::id_t offset, const unsigned char* bits, unsigned size)
        {
            for (unsigned i = 0; i < size; ++i)
                func_(handle_, offset + bits[i]);
        }
        void add_range(bm::id_t offset, unsigned size)
        {
            for (unsigned i = 0; i < size; ++i)
                func_(handle_, offset + i);
        }
        
        void* handle_;
        bit_visitor_callback_type func_;
    };
    
    callback_adaptor func(handle_ptr, callback_ptr);
    bm::for_each_bit(bv, func);
}


/**
    Algorithms for rank compression of bit-vector

    1. Source vector (bv_src) is a subset of index vector (bv_idx)
    2. As a subset it can be collapsed using bit-rank method, where each position
    in the source vector is defined by population count (range) [0..index_position] (count_range())
    As a result all integer set of source vector gets re-mapped in
    accord with the index vector.
 
    \ingroup setalgo
*/
template<class BV>
class bvector_rank_compressor
{
public:
    typedef BV                         bvector_type;
    typedef typename BV::blocks_count  block_count_type;
public:
    /**
    Basic algorithm based on two palallel iterators/enumerators set of source
    vector gets re-mapped in accord with the index/rank vector.

    \param bv_target - target bit-vector
    \param bv_idx    - index (rank) vector used for address recalculation
    \param bv_src    - source vector for re-mapping
    */
    void compress(BV& bv_target, const BV& bv_idx, const BV& bv_src);
    
    /**
    \brief Source vector priority + index based rank
    
    @sa compress
    */
    void compress_by_source(BV& bv_target,
                                const BV& bv_idx,
                                const block_count_type& bc_idx,
                                const BV& bv_src);
};

template<class BV>
void bvector_rank_compressor<BV>::compress(BV& bv_target,
                                           const BV& bv_idx,
                                           const BV& bv_src)
{
    bv_target.clear();
    bv_target.init();

    if (&bv_idx == &bv_src)
    {
        bv_target = bv_src;
        return;
    }
    
    typedef typename BV::enumerator enumerator_t;
    enumerator_t en_s = bv_src.first();
    enumerator_t en_i = bv_idx.first();

    bm::id_t r_idx = 0;
    bm::id_t i, s;
    for (; en_i.valid(); )
    {
        if (!en_s.valid())
            return;
        i = *en_i; s = *en_s;

        BM_ASSERT(s >= i);
        BM_ASSERT(bv_idx.test(i));

        if (s < i)
            return;
        
        if (i == s)
        {
            bv_target.set_bit_no_check(r_idx++);
            ++en_i; ++en_s;
        }
        else
        {
            if (s > i)
            {
                if ((s - i) >= 256) // sufficiently far away, jump
                {
                    bm::id_t r_dist = bv_idx.count_range(i + 1, s);
                    en_i.go_to(s);
                    BM_ASSERT(en_i.valid());
                    r_idx += r_dist;
                }
                else  // small distance, iterate to close the gap
                {
                    for (; s > i; ++r_idx)
                    {
                        ++en_i;
                        i = *en_i;
                        BM_ASSERT(en_i.valid());
                        if (!en_i.valid())
                            return;
                    } // for
                }
            }
        }
    } // for
}

template<class BV>
void bvector_rank_compressor<BV>::compress_by_source(BV& bv_target,
                                           const BV& bv_idx,
                                           const block_count_type& bc_idx,
                                           const BV& bv_src)
{
    /// Rand compressor visitor functor
    /// @internal
    ///
    struct visitor_func
    {
        visitor_func(bvector_type&       bv_out,
                     const bvector_type& bv_index,
                     const block_count_type& bc_index)
        : bv_target_(bv_out),
          bv_index_(bv_index),
          bc_index_(bc_index)
        {}
        
        void add_bits(bm::id_t arr_offset, const unsigned char* bits, unsigned bits_size)
        {
            for (unsigned i = 0; i < bits_size; ++i)
            {
                bm::id_t idx = arr_offset + bits[i];
                BM_ASSERT(bv_index_.test(idx));

                bm::id_t r_idx = bv_index_.count_to(idx, bc_index_) - 1;
                bv_target_.set_bit_no_check(r_idx);
            }
        }
        void add_range(bm::id_t arr_offset, unsigned sz)
        {
            for (unsigned i = 0; i < sz; ++i)
            {
                bm::id_t idx = i + arr_offset;
                BM_ASSERT(bv_index_.test(idx));

                bm::id_t r_idx = bv_index_.count_to(idx, bc_index_) - 1;
                bv_target_.set_bit_no_check(r_idx);
            }
        }
        
        bvector_type&           bv_target_;
        const bvector_type&     bv_index_;
        const block_count_type& bc_index_;
    };
    // ------------------------------------


    bv_target.clear();
    bv_target.init();

    if (&bv_idx == &bv_src)
    {
        bv_target = bv_src;
        return;
    }
    visitor_func func(bv_target, bv_idx, bc_idx);
    bm::for_each_bit(bv_src, func);
}


} // bm

#include "bmundef.h"

#endif
