#ifndef BMAGGREGATOR__H__INCLUDED__
#define BMAGGREGATOR__H__INCLUDED__
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

/*! \file bmaggregator.h
    \brief Algorithms for fast aggregation of N bvectors
*/

#include "bm.h"
#include "bmfunc.h"
#include "bmdef.h"

#include "bmalgo_impl.h"


namespace bm
{


/**
    Algorithms for fast aggregation of a group of bit-vectors
 
    Algorithms of this class use cache locality optimizations and efficient
    on cases, wehen we need to apply the same logical operation (aggregate)
    more than 2x vectors.
 
    \ingroup setalgo
*/
template<typename BV>
class aggregator
{
public:
    typedef BV                         bvector_type;
    typedef bvector_type*              bvector_type_ptr;
    
public:
    aggregator();
    ~aggregator();
    
    
    /**
        Aggregate group of vectors using logical OR
        \param bv_target - target vector
        \param bv_src    - array of pointers on bit-vector aggregate arguments
        \param src_size  - size of bv_src (how many vectors to aggregate)
    */
    void combine_or(bvector_type& bv_target,
                    const bvector_type_ptr* bv_src, unsigned src_size);

    /**
        Aggregate group of vectors using logical AND
        \param bv_target - target vector
        \param bv_src    - array of pointers on bit-vector aggregate arguments
        \param src_size  - size of bv_src (how many vectors to aggregate)
    */
    void combine_and(bvector_type& bv_target,
                    const bvector_type_ptr* bv_src, unsigned src_size);


    /// Horizontal OR aggregation (potentially slower) method
    void combine_or_horizontal(bvector_type& bv_target,
                               const bvector_type_ptr* bv_src, unsigned src_size);
    /// Horizontal AND aggregation (potentially slower) method
    void combine_and_horizontal(bvector_type& bv_target,
                                const bvector_type_ptr* bv_src, unsigned src_size);

protected:
    typedef typename bvector_type::blocks_manager_type blocks_manager_type;
    
    void combine_or(unsigned i, unsigned j,
                    bvector_type& bv_target,
                    const bvector_type_ptr* bv_src, unsigned src_size);

    void combine_and(unsigned i, unsigned j,
                    bvector_type& bv_target,
                    const bvector_type_ptr* bv_src, unsigned src_size);

    static
    unsigned resize_target(bvector_type& bv_target,
                           const bvector_type_ptr* bv_src,
                           unsigned src_size);
    
    bm::word_t* sort_input_blocks_or(const bvector_type_ptr* bv_src,
                                     unsigned src_size,
                                     unsigned i, unsigned j,
                                     unsigned* arg_blk_count,
                                     unsigned* arg_blk_gap_count);
    
    bm::word_t* sort_input_blocks_and(const bvector_type_ptr* bv_src,
                                      unsigned src_size,
                                      unsigned i, unsigned j,
                                      unsigned* arg_blk_count,
                                      unsigned* arg_blk_gap_count);

    bool process_bit_blocks_or(blocks_manager_type& bman_target,
                               unsigned i, unsigned j, unsigned block_count);

    bool process_gap_blocks_or(blocks_manager_type& bman_target,
                               unsigned i, unsigned j, unsigned block_count);
    static
    unsigned find_effective_sub_block_size(unsigned i,
                                           const bvector_type_ptr* bv_src,
                                           unsigned src_size);
private:
    struct data_buffer
    {
        BM_DECLARE_TEMP_BLOCK(tb);
        bm::gap_word_t        gap_res_buf1[bm::gap_equiv_len * 3]; ///< temp 1
        bm::gap_word_t        gap_res_buf2[bm::gap_equiv_len * 3]; ///< temp 2
        bm::gap_word_t        gap_res_buf3[bm::gap_equiv_len * 6]; ///< temp 3
        const bm::word_t*     v_arg_blk[256];     ///< source blocks list
        const bm::gap_word_t* v_arg_blk_gap[256]; ///< source GAP blocks list
    };
    
    aggregator(const aggregator&) = delete;
    aggregator& operator=(const aggregator&) = delete;
    
private:
    data_buffer*          db_; ///< data buffer pointer
};


// ------------------------------------------------------------------------
//
// ------------------------------------------------------------------------

#if defined(BMSSE2OPT) || defined(BMSSE42OPT)
#define BM_ALLOC_ALIGN 16
#endif
#if defined(BMAVX2OPT)
#define BM_ALLOC_ALIGN 32
#endif

template<typename BV>
aggregator<BV>::aggregator()
{
#if defined(BM_ALLOC_ALIGN)
#ifdef _MSC_VER
    db_ = (data_buffer*) ::_aligned_malloc(sizeof(data_buffer), BM_ALLOC_ALIGN);
#else
    db_ = (data_buffer*) ::_mm_malloc(sizeof(data_buffer), BM_ALLOC_ALIGN);
#endif
#else
    db_ = (data_buffer*) ::malloc(sizeof(data_buffer));
#endif
    if (!db_)
    {
        BV::throw_bad_alloc();
    }
}

// ------------------------------------------------------------------------

template<typename BV>
aggregator<BV>::~aggregator()
{
#ifdef BM_ALLOC_ALIGN
# ifdef _MSC_VER
    ::_aligned_free(db_);
#else
    ::_mm_free(db_);
# endif
#else
    ::free(db_);
#endif
}

#undef BM_ALLOC_ALIGN

// ------------------------------------------------------------------------

template<typename BV>
void aggregator<BV>::combine_or(bvector_type& bv_target,
                        const bvector_type_ptr* bv_src, unsigned src_size)
{
    BM_ASSERT(src_size);

    unsigned top_blocks = resize_target(bv_target, bv_src, src_size);
    for (unsigned i = 0; i < top_blocks; ++i)
    {
        unsigned set_array_max = find_effective_sub_block_size(i, bv_src, src_size);
        unsigned j = 0;
        do
        {
            combine_or(i, j, bv_target, bv_src, src_size);
            ++j;
        } while (j < set_array_max);
    } // for i
}

// ------------------------------------------------------------------------

template<typename BV>
void aggregator<BV>::combine_and(bvector_type& bv_target,
                        const bvector_type_ptr* bv_src, unsigned src_size)
{
    BM_ASSERT(src_size);

    unsigned top_blocks = resize_target(bv_target, bv_src, src_size);
    for (unsigned i = 0; i < top_blocks; ++i)
    {
        unsigned set_array_max = find_effective_sub_block_size(i, bv_src, src_size);
        unsigned j = 0;
        do
        {
            combine_and(i, j, bv_target, bv_src, src_size);
            ++j;
        } while (j < set_array_max);
    } // for i
}

// ------------------------------------------------------------------------

template<typename BV>
unsigned
aggregator<BV>::find_effective_sub_block_size(unsigned i,
                                              const bvector_type_ptr* bv_src,
                                              unsigned src_size) 
{
    unsigned max_size = 1;
/*    
    #ifdef BM64_AVX2
        max_size = 4;
    #elif defined(BM64_SSE4)
        max_size = 2;
    #else
        max_size = 2;
    #endif
*/
    
    for (unsigned k = 0; k < src_size; ++k)
    {
        const bvector_type* bv = bv_src[k];
        BM_ASSERT(bv);
        const typename bvector_type::blocks_manager_type& bman_arg = bv->get_blocks_manager();

        const bm::word_t* const* blk_blk_arg = bman_arg.get_topblock(i);
        if (!blk_blk_arg)
            continue;

        for (unsigned j = bm::set_array_size - 1; j > max_size; --j)
        {
            if (blk_blk_arg[j])
            {
                max_size = j; break;
            }
        }
/*
        // backward scan to find last non-NULL block
        for (unsigned j = bm::set_array_size-1; j > max_size; )
        {
        #ifdef BM64_AVX2
            BM_ASSERT(j - 4);
            const void* addr = &(blk_blk_arg[j - 4]);
            if (!avx2_test_all_zero_wave(addr))
            {
                max_size = j; break;
            }
            j -= 4;
        #elif defined(BM64_SSE4)
            const void* addr = blk_blk_arg + j - 2;
            if (!sse42_test_all_zero_wave(addr))
            {
                max_size = j; break;
            }
            j -= 2;
        #else
            BM_ASSERT(j - 2);
            if (blk_blk_arg[j] || blk_blk_arg[j-1])
            {
                max_size = j; break;
            }
            j -= 2;
        #endif
        } // for j
*/
    } // for k
    ++max_size;
    BM_ASSERT(max_size <= bm::set_array_size);
    
    return max_size;
}

// ------------------------------------------------------------------------

template<typename BV>
void aggregator<BV>::combine_or(unsigned i, unsigned j,
                                bvector_type& bv_target,
                                const bvector_type_ptr* bv_src,
                                unsigned src_size)
{
    typename bvector_type::blocks_manager_type& bman_target = bv_target.get_blocks_manager();

    unsigned arg_blk_count = 0;
    unsigned arg_blk_gap_count = 0;
    bm::word_t* blk =
        sort_input_blocks_or(bv_src, src_size, i, j,
                             &arg_blk_count, &arg_blk_gap_count);

    BM_ASSERT(blk == 0 || blk == FULL_BLOCK_FAKE_ADDR);

    if (blk == FULL_BLOCK_FAKE_ADDR) // nothing to do - golden block(!)
    {
        bman_target.check_alloc_top_subblock(i);
        bman_target.set_block_ptr(i, j, blk);
    }
    else
    {
        blk = db_->tb;
        if (arg_blk_count || arg_blk_gap_count)
        {
            bool all_one =
                process_bit_blocks_or(bman_target, i, j, arg_blk_count);
            if (!all_one)
            {
                if (arg_blk_gap_count)
                {
                    all_one =
                        process_gap_blocks_or(bman_target, i, j, arg_blk_gap_count);
                }
                if (!all_one)
                {
                    // we have some results, allocate block and copy from temp
                    bman_target.check_alloc_top_subblock(i);
                    blk = bman_target.get_allocator().alloc_bit_block();
                    bman_target.set_block_ptr(i, j, blk);
                    bm::bit_block_copy(blk, db_->tb);
                }
            }
        }
    }

}

// ------------------------------------------------------------------------

template<typename BV>
void aggregator<BV>::combine_and(unsigned i, unsigned j,
                                 bvector_type& bv_target,
                                 const bvector_type_ptr* bv_src,
                                 unsigned src_size)
{
    typename bvector_type::blocks_manager_type& bman_target = bv_target.get_blocks_manager();

    unsigned arg_blk_count = 0;
    unsigned arg_blk_gap_count = 0;
    bm::word_t* blk =
        sort_input_blocks_and(bv_src, src_size,
                              i, j,
                              &arg_blk_count, &arg_blk_gap_count);

    BM_ASSERT(blk == 0 || blk == FULL_BLOCK_FAKE_ADDR);

    if (!blk) // nothing to do - golden block(!)
    {
        bman_target.check_alloc_top_subblock(i);
        bman_target.set_block_ptr(i, j, blk);
    }
    else
    {
        blk = db_->tb;
        if (arg_blk_count || arg_blk_gap_count)
        {
            bool all_one =
                process_bit_blocks_or(bman_target, i, j, arg_blk_count);
            if (!all_one)
            {
                if (arg_blk_gap_count)
                {
                    all_one =
                        process_gap_blocks_or(bman_target, i, j, arg_blk_gap_count);
                }
                if (!all_one)
                {
                    // we have some results, allocate block and copy from temp
                    bman_target.check_alloc_top_subblock(i);
                    blk = bman_target.get_allocator().alloc_bit_block();
                    bman_target.set_block_ptr(i, j, blk);
                    bm::bit_block_copy(blk, db_->tb);
                }
            }
        }
    }

}


// ------------------------------------------------------------------------

template<typename BV>
bool aggregator<BV>::process_gap_blocks_or(blocks_manager_type& bman_target,
                                           unsigned i, unsigned j,
                                           unsigned arg_blk_gap_count)
{
    bm::word_t* blk = db_->tb;
    bool all_one;

    unsigned k = 0;
#if 0
    unsigned unroll_factor = 4;
    unsigned len = arg_blk_gap_count - k;
    unsigned len_unr = len - (len % unroll_factor);
    
    const bm::gap_word_t* res;
    unsigned res_len;

    for (; k < len_unr; k+=unroll_factor)
    {
        {
            const bm::gap_word_t* gap1 = db_->v_arg_blk_gap[k];
            const bm::gap_word_t* gap2 = db_->v_arg_blk_gap[k+1];
            res = bm::gap_operation_or(gap1, gap2, db_->gap_res_buf1, res_len);
            BM_ASSERT(res == db_->gap_res_buf1);
            if (bm::gap_is_all_one(res, bm::gap_max_bits))
            {
                bman_target.set_block(i, j, FULL_BLOCK_FAKE_ADDR, false);
                return true; // golden block found (all 1s)!
            }
        }
        {
            const bm::gap_word_t* gap1 = db_->v_arg_blk_gap[k+2];
            const bm::gap_word_t* gap2 = db_->v_arg_blk_gap[k+3];
            res = bm::gap_operation_or(gap1, gap2, db_->gap_res_buf2, res_len);
            BM_ASSERT(res == db_->gap_res_buf2);
            if (bm::gap_is_all_one(res, bm::gap_max_bits))
            {
                bman_target.set_block(i, j, FULL_BLOCK_FAKE_ADDR, false);
                return true; // golden block found (all 1s)!
            }
        }
        
        res = bm::gap_operation_or(db_->gap_res_buf1, db_->gap_res_buf2, db_->gap_res_buf3, res_len);
        BM_ASSERT(res == db_->gap_res_buf3);
        if (bm::gap_is_all_one(res, bm::gap_max_bits))
        {
            bman_target.set_block(i, j, FULL_BLOCK_FAKE_ADDR, false);
            return true; // golden block found (all 1s)!
        }
        
        // accumulate the result of 4x gap OR
        bm::gap_add_to_bitset(blk, res);
    }
    unroll_factor = 2;
    len = arg_blk_gap_count - k;
    len_unr = len - (len % unroll_factor);

    for (; k < len_unr; k+=unroll_factor)
    {
        {
            const bm::gap_word_t* gap1 = db_->v_arg_blk_gap[k];
            const bm::gap_word_t* gap2 = db_->v_arg_blk_gap[k+1];
            res = bm::gap_operation_or(gap1, gap2, db_->gap_res_buf1, res_len);
            BM_ASSERT(res == db_->gap_res_buf1);
            if (bm::gap_is_all_one(res, bm::gap_max_bits))
            {
                bman_target.set_block(i, j, FULL_BLOCK_FAKE_ADDR, false);
                return true; // golden block found (all 1s)!
            }
        }
        // accumulate the result of 2x gap OR
        bm::gap_add_to_bitset(blk, res);
    }
#endif

    for (; k < arg_blk_gap_count; ++k)
    {
        bm::gap_add_to_bitset(blk, db_->v_arg_blk_gap[k]);
    }
    
    all_one = bm::is_bits_one((bm::wordop_t*) blk,
              (bm::wordop_t*) (blk + bm::set_block_size));
    if (all_one)
    {
        bman_target.set_block(i, j, FULL_BLOCK_FAKE_ADDR, false);
        return true;
    }
    
    return false;
}

// ------------------------------------------------------------------------

template<typename BV>
bool aggregator<BV>::process_bit_blocks_or(blocks_manager_type& bman_target,
                                           unsigned i, unsigned j,
                                           unsigned arg_blk_count)
{
    bm::word_t* blk = db_->tb;
    bool all_one;

    unsigned k = 0;

    if (arg_blk_count)  // copy the first block
        bm::bit_block_copy(blk, db_->v_arg_blk[k++]);
    else
        bm::bit_block_set(blk, 0);

    // process all BIT blocks
    //
    unsigned unroll_factor, len, len_unr;
    
    unroll_factor = 4;
    len = arg_blk_count - k;
    len_unr = len - (len % unroll_factor);
    for( ;k < len_unr; k+=unroll_factor)
    {
        all_one = bm::bit_block_or_5way(blk,
                                        db_->v_arg_blk[k], db_->v_arg_blk[k+1],
                                        db_->v_arg_blk[k+2], db_->v_arg_blk[k+3]);
        if (all_one)
        {
            BM_ASSERT(blk == db_->tb);
            BM_ASSERT(bm::is_bits_one((bm::wordop_t*) blk, (bm::wordop_t*) (blk + bm::set_block_size)));
            bman_target.set_block(i, j, FULL_BLOCK_FAKE_ADDR, false);
            return true;
        }
    } // for k

    unroll_factor = 2;
    len = arg_blk_count - k;
    len_unr = len - (len % unroll_factor);
    for( ;k < len_unr; k+=unroll_factor)
    {
        all_one = bm::bit_block_or_3way(blk, db_->v_arg_blk[k], db_->v_arg_blk[k+1]);
        if (all_one)
        {
            BM_ASSERT(blk == db_->tb);
            BM_ASSERT(bm::is_bits_one((bm::wordop_t*) blk, (bm::wordop_t*) (blk + bm::set_block_size)));
            bman_target.set_block(i, j, FULL_BLOCK_FAKE_ADDR, false);
            return true;
        }
    } // for k

    for (; k < arg_blk_count; ++k)
    {
        all_one = bm::bit_block_or(blk, db_->v_arg_blk[k]);
        if (all_one)
        {
            BM_ASSERT(blk == db_->tb);
            BM_ASSERT(bm::is_bits_one((bm::wordop_t*) blk, (bm::wordop_t*) (blk + bm::set_block_size)));
            bman_target.set_block(i, j, FULL_BLOCK_FAKE_ADDR, false);
            return true;
        }
    } // for k

    return false;
}

// ------------------------------------------------------------------------

template<typename BV>
unsigned aggregator<BV>::resize_target(bvector_type& bv_target,
                            const bvector_type_ptr* bv_src, unsigned src_size)
{
    typename bvector_type::blocks_manager_type& bman_target = bv_target.get_blocks_manager();
    if (!bman_target.is_init())
        bman_target.init_tree();
    else
        bv_target.clear(true);
    
    unsigned top_blocks = bman_target.top_block_size();
    auto size = bv_target.size();

    // pre-scan to do target size harmonization
    for (unsigned i = 0; i < src_size; ++i)
    {
        const bvector_type* bv = bv_src[i];
        BM_ASSERT(bv);
        const typename bvector_type::blocks_manager_type& bman_arg = bv->get_blocks_manager();
        unsigned arg_top_blocks = bman_arg.top_block_size();
        if (arg_top_blocks > top_blocks)
            top_blocks = bman_target.reserve_top_blocks(arg_top_blocks);
        auto arg_size = bv->size();
        if (arg_size > size)
        {
            bv_target.resize(arg_size);
            size = arg_size;
        }
    } // for i
    return top_blocks;
}

// ------------------------------------------------------------------------

template<typename BV>
bm::word_t* aggregator<BV>::sort_input_blocks_or(const bvector_type_ptr* bv_src,
                                                 unsigned src_size,
                                                 unsigned i, unsigned j,
                                                 unsigned* arg_blk_count,
                                                 unsigned* arg_blk_gap_count)
{
    bm::word_t* blk = 0;
    for (unsigned k = 0; k < src_size; ++k)
    {
        const bvector_type* bv = bv_src[k];
        BM_ASSERT(bv);
        const typename bvector_type::blocks_manager_type& bman_arg = bv->get_blocks_manager();
        const bm::word_t* arg_blk = bman_arg.get_block_ptr(i, j);
        if (!arg_blk)
            continue;
        if (BM_IS_GAP(arg_blk))
        {
            db_->v_arg_blk_gap[*arg_blk_gap_count] = BMGAP_PTR(arg_blk);
            (*arg_blk_gap_count)++;
        }
        else // FULL or bit block
        {
            if (IS_FULL_BLOCK(arg_blk))
            {
                blk = FULL_BLOCK_FAKE_ADDR;
                *arg_blk_gap_count = *arg_blk_count = 0; // nothing to do
                break;
            }
            db_->v_arg_blk[*arg_blk_count] = arg_blk;
            (*arg_blk_count)++;
        }
    } // for k
    return blk;
}

// ------------------------------------------------------------------------

template<typename BV>
bm::word_t* aggregator<BV>::sort_input_blocks_and(const bvector_type_ptr* bv_src,
                                                  unsigned src_size,
                                                  unsigned i, unsigned j,
                                                  unsigned* arg_blk_count,
                                                  unsigned* arg_blk_gap_count)
{
    bm::word_t* blk = 0;
    for (unsigned k = 0; k < src_size; ++k)
    {
        const bvector_type* bv = bv_src[k];
        BM_ASSERT(bv);
        const typename bvector_type::blocks_manager_type& bman_arg = bv->get_blocks_manager();
        const bm::word_t* arg_blk = bman_arg.get_block_ptr(i, j);
        if (!arg_blk)
        {
            blk = 0;
            *arg_blk_gap_count = *arg_blk_count = 0; // nothing to do
            break;
        }
        if (BM_IS_GAP(arg_blk))
        {
            db_->v_arg_blk_gap[*arg_blk_gap_count] = BMGAP_PTR(arg_blk);
            (*arg_blk_gap_count)++;
        }
        else // FULL or bit block
        {
            db_->v_arg_blk[*arg_blk_count] = arg_blk;
            (*arg_blk_count)++;
        }
    } // for k
    return blk;
}



// ------------------------------------------------------------------------

template<typename BV>
void aggregator<BV>::combine_or_horizontal(bvector_type& bv_target,
                     const bvector_type_ptr* bv_src, unsigned src_size)
{
    BM_ASSERT(src_size);
    
    const bvector_type* bv = bv_src[0];
    bv_target = *bv;
    
    for (unsigned i = 1; i < src_size; ++i)
    {
        bv = bv_src[i];
        BM_ASSERT(bv);
        bv_target.bit_or(*bv);
    }
}

// ------------------------------------------------------------------------

template<typename BV>
void aggregator<BV>::combine_and_horizontal(bvector_type& bv_target,
                     const bvector_type_ptr* bv_src, unsigned src_size)
{
    BM_ASSERT(src_size);
    
    const bvector_type* bv = bv_src[0];
    bv_target = *bv;
    
    for (unsigned i = 1; i < src_size; ++i)
    {
        bv = bv_src[i];
        BM_ASSERT(bv);
        bv_target.bit_and(*bv);
    }
}

// ------------------------------------------------------------------------

} // bm

#include "bmundef.h"

#endif
