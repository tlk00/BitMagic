#ifndef BMRS__H__INCLUDED__
#define BMRS__H__INCLUDED__
/*
Copyright(c) 2002-2019 Anatoliy Kuznetsov(anatoliy_kuznetsov at yahoo.com)

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

/*! \file bmrs.h
    \brief Rank-Select index structures
*/

namespace bm
{

/**
    @brief Rank-Select acceleration index
 
    Index uses two-level acceleration structure:
    bcount - running total popcount for all (possible) blocks
    (missing blocks give duplicate counts as POPCNT(N-1) + 0).
    subcount - sub-count inside blocks
 
    @ingroup bvector
*/
template<typename BVAlloc>
class rs_index
{
public:
    typedef BVAlloc        bv_allocator_type;
    
#ifdef BM64ADDR
    typedef bm::id64_t                                   size_type;
    typedef bm::id64_t                                   block_idx_type;
#else
    typedef bm::id_t                                     size_type;
    typedef bm::id_t                                     block_idx_type;
#endif

    typedef bm::pair<bm::gap_word_t, bm::gap_word_t> sb_pair_type;

public:
    rs_index() : total_blocks_(0) {}
    rs_index(const rs_index& rsi) BMNOEXEPT;
    
    /// init arrays to zeros
    void init() BMNOEXEPT;
    
    /// reserve the capacity for block count
    void resize(block_idx_type new_size);
    
    /// copy rs index
    void copy_from(const rs_index& rsi) BMNOEXEPT;
    
    /// set total blocks
    void set_total(size_type t) { total_blocks_ = t; }
    
    /// get total blocks
    size_type get_total() const { return total_blocks_; }
    
    /// return bit-count for specified block
    size_type count(block_idx_type nb) const;
    
    /// return total bit-count for the index
    size_type count() const;

    /// return running bit-count for specified block
    size_type bcount(block_idx_type nb) const;
    
    /// return reference of sub-index pair
    const sb_pair_type& sub_count(block_idx_type nb) const;

    /// set sub-block index pair
    void set_sub_count(block_idx_type nb,
                       bm::gap_word_t first, bm::gap_word_t second);

    /// determine the sub-range within a bit-block
    unsigned find_sub_range(block_idx_type block_bit_pos) const;
    
    /// determine block sub-range for rank search
    bm::gap_word_t select_sub_range(block_idx_type nb, unsigned& rank) const;
    
    /// assign bcount value for block nb
    void set_bcount(block_idx_type nb, unsigned v);

    /// C-style direct accessor to top-level blocks vector
    const unsigned* bcount_begin() const { return block_count_.begin(); }
    
private:
    typedef bm::heap_vector<size_type, bv_allocator_type>     uvector_type;
    typedef bm::heap_vector<sb_pair_type, bv_allocator_type>  sbvector_type;

    uvector_type   block_count_;  ///< top level accumulated blocks count
    sbvector_type  sub_count_;    ///< second level (subblock) counts
    size_type      total_blocks_; ///< total bit-blocks in the index
};

//---------------------------------------------------------------------
//
//---------------------------------------------------------------------

template<typename BVAlloc>
rs_index<BVAlloc>::rs_index(const rs_index<BVAlloc>& rsi) BMNOEXEPT
{
    copy_from(rsi);
}

//---------------------------------------------------------------------


template<typename BVAlloc>
void rs_index<BVAlloc>::init() BMNOEXEPT
{
    block_count_.resize(0);
    sub_count_.resize(0);
    total_blocks_ = 0;
}

//---------------------------------------------------------------------

template<typename BVAlloc>
void rs_index<BVAlloc>::resize(block_idx_type new_size)
{
    block_count_.resize(new_size);
    sub_count_.resize(new_size);
}

//---------------------------------------------------------------------

template<typename BVAlloc>
void rs_index<BVAlloc>::copy_from(const rs_index& rsi) BMNOEXEPT
{
    block_count_ = rsi.block_count_;
    sub_count_ = rsi.sub_count_;
    this->total_blocks_ = rsi.total_blocks_;
}

//---------------------------------------------------------------------

template<typename BVAlloc>
typename rs_index<BVAlloc>::size_type
rs_index<BVAlloc>::count(block_idx_type nb) const
{
    if (nb >= total_blocks_)
        return count();
    return (nb == 0) ? block_count_[nb]
                     : block_count_[nb] - block_count_[nb-1];
}

//---------------------------------------------------------------------

template<typename BVAlloc>
typename rs_index<BVAlloc>::size_type
rs_index<BVAlloc>::count() const
{
    if (!total_blocks_)
        return 0;
    return block_count_[total_blocks_ - 1];
}

//---------------------------------------------------------------------

template<typename BVAlloc>
typename rs_index<BVAlloc>::size_type
rs_index<BVAlloc>::bcount(block_idx_type nb) const
{
    if (nb >= total_blocks_)
        return count();
    return block_count_[nb];
}

//---------------------------------------------------------------------

template<typename BVAlloc>
const typename rs_index<BVAlloc>::sb_pair_type&
rs_index<BVAlloc>::sub_count(block_idx_type nb) const
{
//    BM_ASSERT(nb < total_blocks_);
    if (nb > total_blocks_)
        nb = total_blocks_;
    const sb_pair_type& sbp = sub_count_[nb];
    return sbp;
}

//---------------------------------------------------------------------

template<typename BVAlloc>
void rs_index<BVAlloc>::set_sub_count(block_idx_type nb,
                                      bm::gap_word_t first,
                                      bm::gap_word_t second)
{
    sb_pair_type& sbp = sub_count_[nb];
    sbp.first = first;
    sbp.second = second;
}

//---------------------------------------------------------------------

template<typename BVAlloc>
unsigned rs_index<BVAlloc>::find_sub_range(block_idx_type block_bit_pos) const
{
    return (block_bit_pos <= rs3_border0) ? 0 :
            (block_bit_pos > rs3_border1) ? 2 : 1;
}

//---------------------------------------------------------------------

template<typename BVAlloc>
bm::gap_word_t rs_index<BVAlloc>::select_sub_range(block_idx_type nb,
                                                   unsigned& rank) const
{
    const sb_pair_type& sbp = sub_count(nb);
    if (rank > sbp.first)
    {
        rank -= sbp.first;
        if (rank > sbp.second)
        {
            rank -= sbp.second;
            return rs3_border1 + 1;
        }
        else
            return rs3_border0 + 1;
    }
    return 0;
}
//---------------------------------------------------------------------

template<typename BVAlloc>
void rs_index<BVAlloc>::set_bcount(block_idx_type nb, unsigned v)
{
    block_count_[nb] = v;
}

//---------------------------------------------------------------------


}
#endif
