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

#ifndef BM__H__INCLUDED__
// BitMagic utility headers do not include main "bm.h" declaration 
// #include "bm.h" or "bm64.h" explicitly 
# error missing include (bm.h or bm64.h)
#endif

#include "bmfunc.h"
#include "bmdef.h"

#include "bmalgo_impl.h"



namespace bm
{

#define BM_SCANNER_OP(x) \
    if (0 != (block = blk_blk[j+x])) \
    { \
        if (BM_IS_GAP(block)) \
        { \
            bm::for_each_gap_blk(BMGAP_PTR(block), (r+j+x)*bm::bits_in_block,\
                                 bit_functor); \
        } \
        else \
        { \
            bm::for_each_bit_blk(block, (r+j+x)*bm::bits_in_block,bit_functor); \
        } \
    }
    

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
    typedef typename BV::size_type size_type;

    const typename BV::blocks_manager_type& bman = bv.get_blocks_manager();
    bm::word_t*** blk_root = bman.top_blocks_root();
    
    if (!blk_root)
        return;
    
    unsigned tsize = bman.top_block_size();
    for (unsigned i = 0; i < tsize; ++i)
    {
        bm::word_t** blk_blk = blk_root[i];
        if (!blk_blk)
            continue;

        if ((bm::word_t*)blk_blk == FULL_BLOCK_FAKE_ADDR)
            blk_blk = FULL_SUB_BLOCK_REAL_ADDR;
        
        const bm::word_t* block;
        size_type r = size_type(i) * bm::set_sub_array_size;
        unsigned j = 0;
        do
        {
        #ifdef BM64_AVX2
            if (!avx2_test_all_zero_wave(blk_blk + j))
            {
                BM_SCANNER_OP(0)
                BM_SCANNER_OP(1)
                BM_SCANNER_OP(2)
                BM_SCANNER_OP(3)
            }
            j += 4;
        #elif defined(BM64_SSE4)
            if (!sse42_test_all_zero_wave(blk_blk + j))
            {
                BM_SCANNER_OP(0)
                BM_SCANNER_OP(1)
            }
            j += 2;
        #else
            BM_SCANNER_OP(0)
            ++j;
        #endif
        
        } while (j < bm::set_sub_array_size);
        
    }  // for i
}

#undef BM_SCANNER_OP

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
    typedef typename BV::size_type size_type;
    // private adaptor for C-style callbacks
    struct callback_adaptor
    {
        callback_adaptor(void* h, bit_visitor_callback_type cb_func)
        : handle_(h), func_(cb_func)
        {}
        
        void add_bits(size_type offset, const unsigned char* bits, unsigned size)
        {
            for (unsigned i = 0; i < size; ++i)
                func_(handle_, offset + bits[i]);
        }
        void add_range(size_type offset, unsigned size)
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
template<typename BV>
class rank_compressor
{
public:
    typedef BV                         bvector_type;
    typedef typename BV::size_type     size_type;
    typedef typename BV::rs_index_type rs_index_type;
    enum buffer_cap
    {
        n_buffer_cap = 1024
    };
public:

    /**
    Rank decompression
    */
    void decompress(BV& bv_target, const BV& bv_idx, const BV& bv_src);

    /**
    Rank compression algorithm based on two palallel iterators/enumerators set of source
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
                                const rs_index_type& bc_idx,
                                const BV& bv_src);
};


// ------------------------------------------------------------------------
//
// ------------------------------------------------------------------------


template<class BV>
void rank_compressor<BV>::compress(BV& bv_target,
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
    size_type ibuffer[n_buffer_cap];
    size_type b_size;
    
    typedef typename BV::enumerator enumerator_t;
    enumerator_t en_s = bv_src.first();
    enumerator_t en_i = bv_idx.first();

    size_type r_idx = b_size = 0;
    size_type i, s;
    
    for (; en_i.valid(); )
    {
        if (!en_s.valid())
            break;
        i = *en_i; s = *en_s;

        BM_ASSERT(s >= i);
        BM_ASSERT(bv_idx.test(i));

        if (i == s)
        {
            ibuffer[b_size++] = r_idx++;
            if (b_size == n_buffer_cap)
            {
                bm::combine_or(bv_target, ibuffer+0, ibuffer+b_size);
                b_size ^= b_size; // = 0
            }
            ++en_i; ++en_s;
            continue;
        }
        BM_ASSERT(s > i);
        
        size_type dist = s - i;
        if (dist >= 64) // sufficiently far away, jump
        {
            size_type r_dist = bv_idx.count_range(i + 1, s);
            r_idx += r_dist;
            en_i.go_to(s);
            BM_ASSERT(en_i.valid());
        }
        else  // small distance, iterate to close the gap
        {
            for (; s > i; ++r_idx)
            {
                ++en_i;
                i = *en_i;
            } // for
            BM_ASSERT(en_i.valid());
        }
    } // for
    
    if (b_size)
    {
        bm::combine_or(bv_target, ibuffer+0, ibuffer+b_size);
    }

}

// ------------------------------------------------------------------------


template<class BV>
void rank_compressor<BV>::decompress(BV& bv_target,
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
    
    size_type r_idx, i, s, b_size;
    size_type ibuffer[n_buffer_cap];
    
    b_size = r_idx = 0;

    typedef typename BV::enumerator enumerator_t;
    enumerator_t en_s = bv_src.first();
    enumerator_t en_i = bv_idx.first();
    for (; en_i.valid(); )
    {
        if (!en_s.valid())
            break;
        s = *en_s;
        i = *en_i;
        if (s == r_idx)
        {
            ibuffer[b_size++] = i;
            if (b_size == n_buffer_cap)
            {
                bm::combine_or(bv_target, ibuffer+0, ibuffer+b_size);
                b_size ^= b_size; // = 0
            }
            ++en_i; ++en_s; ++r_idx;
            continue;
        }
        // source is "faster" than index, need to re-align
        BM_ASSERT(s > r_idx);
        size_type rank = s - r_idx + 1u;
        size_type new_pos = 0;
        
        if (rank < 256)
        {
            en_i.skip(s - r_idx);
            BM_ASSERT(en_i.valid());
            new_pos = *en_i;
        }
        else
        {
            bv_idx.find_rank(rank, i, new_pos);
            BM_ASSERT(new_pos);
            en_i.go_to(new_pos);
            BM_ASSERT(en_i.valid());
        }
        
        r_idx = s;
        ibuffer[b_size++] = new_pos;
        if (b_size == n_buffer_cap)
        {
            bm::combine_or(bv_target, ibuffer+0, ibuffer+b_size);
            b_size ^= b_size; // = 0
        }
        ++en_i; ++en_s; ++r_idx;
        
    } // for en
    
    if (b_size)
    {
        bm::combine_or(bv_target, ibuffer+0, ibuffer+b_size);
    }
}

// ------------------------------------------------------------------------

template<class BV>
void rank_compressor<BV>::compress_by_source(BV& bv_target,
                                             const BV& bv_idx,
                                             const rs_index_type& bc_idx,
                                             const BV& bv_src)
{
    /// Rank compressor visitor (functor)
    /// @internal
    struct visitor_func
    {
        visitor_func(bvector_type&        bv_out,
                     const bvector_type&  bv_index,
                     const rs_index_type& bc_index)
        : bv_target_(bv_out),
          bv_index_(bv_index),
          bc_index_(bc_index)
        {}
        
        void add_bits(size_type arr_offset, const unsigned char* bits, unsigned bits_size)
        {
            for (unsigned i = 0; i < bits_size; ++i)
            {
                size_type idx = arr_offset + bits[i];
                BM_ASSERT(bv_index_.test(idx));

                size_type r_idx = bv_index_.count_to(idx, bc_index_) - 1;
                bv_target_.set_bit_no_check(r_idx);
            }
        }
        void add_range(size_type arr_offset, size_type sz)
        {
            for (size_type i = 0; i < sz; ++i)
            {
                size_type idx = i + arr_offset;
                BM_ASSERT(bv_index_.test(idx));

                size_type r_idx = bv_index_.count_to(idx, bc_index_) - 1;
                bv_target_.set_bit_no_check(r_idx);
            }
        }
        
        bvector_type&           bv_target_;
        const bvector_type&     bv_index_;
        const rs_index_type&    bc_index_;
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
