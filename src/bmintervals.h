#ifndef BMINTERVALS__H__INCLUDED__
#define BMINTERVALS__H__INCLUDED__

/*
Copyright(c) 2002-2020 Anatoliy Kuznetsov(anatoliy_kuznetsov at yahoo.com)

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
/*! \file bmintervals.h
    \brief Algorithms for bit ranges and intervals
*/

#ifndef BM__H__INCLUDED__
// BitMagic utility headers do not include main "bm.h" declaration
// #include "bm.h" or "bm64.h" explicitly
# error missing include (bm.h or bm64.h)
#endif

#include "bmdef.h"

/** \defgroup bvintervals Algorithms for bit ranges
    Algorithms for bit ranges and intervals
    @ingroup bvector
 */


namespace bm
{

//----------------------------------------------------------------------------

/*!
    \brief Returns true if range is all 1s flanked with 0s
    Function performs the test on a closed range [left, right]
    true interval is all 1s AND test(left-1)==false AND test(right+1)==false
    Examples:
        01110 [1,3] - true
        11110 [0,3] - true
        11110 [1,3] - false
    \param bv   - bit-vector for check
   \param left - index of first bit start checking
   \param right - index of last bit
   \return true/false
   @sa is_all_one_range
*/
template<class BV>
bool is_interval(const BV& bv,
                 typename BV::size_type left,
                 typename BV::size_type right) BMNOEXCEPT
{
    typedef typename BV::block_idx_type block_idx_type;

    const typename BV::blocks_manager_type& bman = bv.get_blocks_manager();

    if (!bman.is_init())
        return false; // nothing to do

    if (right < left)
        bm::xor_swap(left, right);
    if (left == bm::id_max) // out of range
        return false;
    if (right == bm::id_max)
        --right;

    block_idx_type nblock_left = (left >> bm::set_block_shift);
    block_idx_type nblock_right = (right >> bm::set_block_shift);

    if (nblock_left == nblock_right) // same block (fast case)
    {
        unsigned nbit_left = unsigned(left  & bm::set_block_mask);
        unsigned nbit_right = unsigned(right  & bm::set_block_mask);
        if ((nbit_left > 0) && (nbit_right < bm::gap_max_bits-1))
        {
            unsigned i0, j0;
            bm::get_block_coord(nblock_left, i0, j0);
            const bm::word_t* block = bman.get_block_ptr(i0, j0);
            bool b = bm::block_is_interval(block, nbit_left, nbit_right);
            return b;
        }
    }
    bool is_left, is_right, is_all_one;
    is_left = left > 0 ? bv.test(left-1) : false;
    if (is_left == false)
    {
        is_right = (right < (bm::id_max - 1)) ? bv.test(right + 1) : false;
        if (is_left == false && is_right == false)
        {
            is_all_one = bv.is_all_one_range(left, right);
            return is_all_one;
        }
    }
    return false;
}


//----------------------------------------------------------------------------

/*!

    \brief Reverse find index of first 1 bit gap (01110) starting from position
    Reverse scan for the first 1 in a block of continious 1s.
    Method employs closed interval semantics: 0[pos..from]

    \param bv   - bit-vector for search
    \param from - position to start reverse search from
    \param pos - [out] index of the found first 1 bit in a gap of bits
    \return true if search returned result, false if not found
           (start point is zero)

    \sa is_interval, find_interval_end
    \ingroup bvintervals
*/
template<class BV>
bool find_interval_start(const BV& bv,
                         typename BV::size_type from,
                         typename BV::size_type& pos)  BMNOEXCEPT
{
    typedef typename BV::size_type size_type;
    typedef typename BV::block_idx_type block_idx_type;

    const typename BV::blocks_manager_type& bman = bv.get_blocks_manager();

    if (!bman.is_init())
        return false; // nothing to do
    if (!from)
    {
        pos = from;
        return bv.test(from);
    }

    block_idx_type nb = (from >> bm::set_block_shift);
    unsigned i0, j0;
    bm::get_block_coord(nb, i0, j0);

    size_type base_idx;
    unsigned found_nbit;

    const bm::word_t* block = bman.get_block_ptr(i0, j0);
    if (!block)
        return false;
    unsigned nbit = unsigned(from & bm::set_block_mask);
    unsigned res = bm::block_find_interval_start(block, nbit, &found_nbit);

    switch (res)
    {
    case 0: // not interval
        return false;
    case 1: // interval found
        pos = found_nbit + (nb * bm::gap_max_bits);
        return true;
    case 2: // keep scanning
        base_idx = bm::get_block_start<size_type>(i0, j0);
        pos = base_idx + found_nbit;
        if (!nb)
            return true;
        break;
    default:
        BM_ASSERT(0);
    } // switch

    --nb;
    bm::get_block_coord(nb, i0, j0);
    bm::word_t*** blk_root = bman.top_blocks_root();

    for (unsigned i = i0; true; --i)
    {
        bm::word_t** blk_blk = blk_root[i];
        if (!blk_blk)
            return true;
        if ((bm::word_t*)blk_blk == FULL_BLOCK_FAKE_ADDR)
        {
            pos = bm::get_super_block_start<size_type>(i);
            if (!i)
                break;
            continue;
        }
        unsigned j = (i == i0) ? j0 : 255;
        for (; true; --j)
        {
            if ((bm::word_t*)blk_blk == FULL_BLOCK_FAKE_ADDR)
            {
                pos = bm::get_block_start<size_type>(i, j);
                goto loop_j_end; // continue
            }

            block = blk_blk[j];
            if (!block)
                return true;

            res = bm::block_find_interval_start(block,
                                            bm::gap_max_bits-1, &found_nbit);
            switch (res)
            {
            case 0: // not interval (but it was the interval, so last result
                return true;
            case 1: // interval found
                base_idx = bm::get_block_start<size_type>(i, j);
                pos = base_idx + found_nbit;
                return true;
            case 2: // keep scanning
                pos = bm::get_block_start<size_type>(i, j);
                break;
            default:
                BM_ASSERT(0);
            } // switch

            loop_j_end: // continue point
            if (!j)
                break;
        } // for j

        if (!i)
            break;
    } // for i

    return true;
}


//----------------------------------------------------------------------------

/*!
   \brief Reverse find index of first 1 bit gap (01110) starting from position
   Reverse scan for the first 1 in a block of continious 1s.
   Method employs closed interval semantics: 0[pos..from]

   \param bv   - bit-vector for search
   \param from - position to start reverse search from
   \param pos - [out] index of the found first 1 bit in a gap of bits
   \return true if search returned result, false if not found
           (start point is zero)

   \sa is_interval, find_interval_end
    \ingroup bvintervals
*/
template <typename BV>
bool find_interval_end(const BV& bv,
                       typename BV::size_type from,
                       typename BV::size_type & pos)  BMNOEXCEPT
{
    typedef typename BV::block_idx_type block_idx_type;

    if (from == bm::id_max)
        return false;
    const typename BV::blocks_manager_type& bman = bv.get_blocks_manager();

    if (!bman.is_init())
        return false; // nothing to do
    if (from == bm::id_max-1)
    {
        pos = from;
        return bv.test(from);
    }

    block_idx_type nb = (from >> bm::set_block_shift);
    unsigned i0, j0;
    bm::get_block_coord(nb, i0, j0);

    unsigned found_nbit;

    const bm::word_t* block = bman.get_block_ptr(i0, j0);
    if (!block)
        return false;
    unsigned nbit = unsigned(from & bm::set_block_mask);
    unsigned res = bm::block_find_interval_end(block, nbit, &found_nbit);
    switch (res)
    {
    case 0: // not interval
        return false;
    case 1: // interval found
        pos = found_nbit + (nb * bm::gap_max_bits);
        return true;
    case 2: // keep scanning
        pos = found_nbit + (nb * bm::gap_max_bits);
        break;
    default:
        BM_ASSERT(0);
    } // switch

    block_idx_type nblock_right = (bm::id_max >> bm::set_block_shift);
    unsigned i_from, j_from, i_to, j_to;
    bm::get_block_coord(nblock_right, i_to, j_to);
    block_idx_type top_size = bman.top_block_size();
    if (i_to >= top_size)
        i_to = unsigned(top_size-1);

    ++nb;
    bm::word_t*** blk_root = bman.top_blocks_root();
    bm::get_block_coord(nb, i_from, j_from);

    for (unsigned i = i_from; i <= i_to; ++i)
    {
        bm::word_t** blk_blk = blk_root[i];
        if (!blk_blk)
            return true;
        if ((bm::word_t*)blk_blk == FULL_BLOCK_FAKE_ADDR)
        {
            if (i > i_from)
            {
                pos += bm::gap_max_bits * bm::set_sub_array_size;
                continue;
            }
            else
            {
                // TODO: optimization to avoid scanning rest of the super block
            }
        }

        unsigned j = (i == i_from) ? j_from : 0;
        do
        {
            if ((bm::word_t*)blk_blk == FULL_BLOCK_FAKE_ADDR)
            {
                pos += bm::gap_max_bits;
                continue;
            }

            block = blk_blk[j];
            if (!block)
                return true;

            res = bm::block_find_interval_end(block, 0, &found_nbit);
            switch (res)
            {
            case 0: // not interval (but it was the interval, so last result
                return true;
            case 1: // interval found
                pos += found_nbit+1;
                return true;
            case 2: // keep scanning
                pos += bm::gap_max_bits;
                break;
            default:
                BM_ASSERT(0);
            } // switch
        } while (++j < bm::set_sub_array_size);
    } // for i

    return true;
}

//----------------------------------------------------------------------------

/// @internal
///
template<typename BV>
class interval_enumerator
{
public:
#ifndef BM_NO_STL
        typedef std::input_iterator_tag  iterator_category;
#endif
        typedef BV                                         bvector_type;
        typedef typename bvector_type::size_type           size_type;
        typedef typename bvector_type::allocator_type      allocator_type;
        typedef bm::byte_buffer<allocator_type>            buffer_type;

public:

    interval_enumerator()
        : bv_(0), pos_(bm::id_max), end_pos_(bm::id_max), gap_ptr_(0)
    {}

    interval_enumerator(const BV& bv)
        : bv_(&bv), pos_(bm::id_max), end_pos_(bm::id_max), gap_ptr_(0)
    {
        go_to_impl(0, true);
    }

    size_type start() const BMNOEXCEPT;
    size_type end() const BMNOEXCEPT;

    bool valid() const BMNOEXCEPT;

    void go_to(size_type pos);

    bool advance();

protected:
    typedef typename bvector_type::block_idx_type       block_idx_type;

    bool go_to_impl(size_type pos, bool extend_start);

    /// Turn FSM into invalid state (out of range)
    void invalidate() BMNOEXCEPT;


private:
    const BV*                  bv_;      ///!< bit-vector for traversal
    mutable buffer_type        buffer_;  ///!< buffer for decoded block
    size_type                  pos_;     ///!< current position
    size_type                  end_pos_; ///!< current interval end position
    const bm::gap_word_t*      gap_ptr_; ///!< current pointer in GAP block
};


//----------------------------------------------------------------------------
//
//----------------------------------------------------------------------------

template<typename BV>
typename interval_enumerator<BV>::size_type
interval_enumerator<BV>::start() const BMNOEXCEPT
{
    return pos_;
}

//----------------------------------------------------------------------------

template<typename BV>
typename interval_enumerator<BV>::size_type
interval_enumerator<BV>::end() const BMNOEXCEPT
{
    return end_pos_;
}

//----------------------------------------------------------------------------

template<typename BV>
bool interval_enumerator<BV>::valid() const BMNOEXCEPT
{
    return (pos_ != bm::id_max);
}

//----------------------------------------------------------------------------

template<typename BV>
void interval_enumerator<BV>::invalidate() BMNOEXCEPT
{
    pos_ = end_pos_ = bm::id_max;
}

//----------------------------------------------------------------------------

template<typename BV>
void interval_enumerator<BV>::go_to(size_type pos)
{
    go_to_impl(pos, true);
}

//----------------------------------------------------------------------------

template<typename BV>
bool interval_enumerator<BV>::go_to_impl(size_type pos, bool extend_start)
{
    if (!bv_ || !bv_->is_init() || pos == pos < bm::id_max)
    {
        invalidate();
        return false;
    }

    bool found;
    size_type start_pos;

    // go to prolog: identify the true interval start position
    //
    if (extend_start)
    {
        found = bm::find_interval_start(*bv_, pos, start_pos);
        if (!found)
        {
            found = bv_->find(pos, start_pos);
            if (!found)
            {
                invalidate();
                return false;
            }
        }
    }
    else
    {
        found = bv_->find(pos, start_pos);
        if (!found)
        {
            invalidate();
            return false;
        }
    }

    // start position established, start decoding from it
    pos_ = pos = start_pos;

    block_idx_type nb = (pos >> bm::set_block_shift);
    const typename BV::blocks_manager_type& bman = bv_->get_blocks_manager();
    {
        unsigned i0, j0;
        bm::get_block_coord(nb, i0, j0);
        const bm::word_t* block = bman.get_block_ptr(i0, j0);
        BM_ASSERT(block);

        if (block == FULL_BLOCK_FAKE_ADDR)
        {
            // super-long interval, find the end of it
            found = bm::find_interval_end(*bv_, pos, end_pos_);
            BM_ASSERT(found);
            return true;
        }

        if (BM_IS_GAP(block))
        {
            const bm::gap_word_t* BMRESTRICT gap_block = BMGAP_PTR(block);
            unsigned nbit = unsigned(pos  & bm::set_block_mask);

            unsigned is_set;
            unsigned gap_pos = bm::gap_bfind(gap_block, nbit, &is_set);
            BM_ASSERT(is_set);

            end_pos_ = (nb * bm::gap_max_bits) + gap_block[gap_pos];
            if (gap_block[gap_pos] == bm::gap_max_bits-1)
            {
                // it is the end of the GAP block - run search
                //
                if (end_pos_ == bm::id_max-1)
                {
                    gap_ptr_ = 0;
                    return true;
                }
                found = bm::find_interval_end(*bv_, end_pos_+1, start_pos);
                if (found)
                {
                    end_pos_ = start_pos;
                    gap_ptr_ = 0;
                    return true;
                }
            }
            gap_ptr_ = gap_block + gap_pos;
            return true;
        }
    }
    BM_ASSERT(0);
    return false;
}

//----------------------------------------------------------------------------

template<typename BV>
bool interval_enumerator<BV>::advance()
{
    if (end_pos_ == bm::id_max-1)
    {
        invalidate();
        return false;
    }
    bool found;
    if (gap_ptr_) // in GAP block
    {
        ++gap_ptr_; // 0 - GAP
        if (*gap_ptr_ == bm::gap_max_bits-1) // GAP block end
        {
            return go_to_impl(end_pos_+ 1, false);
        }
        unsigned prev = *gap_ptr_;

        ++gap_ptr_; // 1 - GAP
        block_idx_type nb = (pos_ >> bm::set_block_shift);
        pos_ = (nb * bm::gap_max_bits) + prev + 1;
        if (*gap_ptr_ == bm::gap_max_bits-1) // GAP block end
        {
            found = bm::find_interval_end(*bv_, pos_, end_pos_);
            BM_ASSERT(found);
            gap_ptr_ = 0;
            return true;
        }
        end_pos_ = (nb * bm::gap_max_bits) + *gap_ptr_;
        return true;
    }

    return go_to_impl(end_pos_+ 1, false);
}

//----------------------------------------------------------------------------
//
//----------------------------------------------------------------------------


} // namespace bm

#include "bmundef.h"

#endif
