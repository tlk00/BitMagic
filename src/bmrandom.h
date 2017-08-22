#ifndef BMRANDOM__H__INCLUDED__
#define BMRANDOM__H__INCLUDED__
/*
Copyright(c) 2009 Anatoliy Kuznetsov(anatoliy_kuznetsov at yahoo.com)

Permission is hereby granted, free of charge, to any person 
obtaining a copy of this software and associated documentation 
files (the "Software"), to deal in the Software without restriction, 
including without limitation the rights to use, copy, modify, merge, 
publish, distribute, sublicense, and/or sell copies of the Software, 
and to permit persons to whom the Software is furnished to do so, 
subject to the following conditions:

The above copyright notice and this permission notice shall be included 
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, 
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, 
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR 
OTHER DEALINGS IN THE SOFTWARE.

For more information please visit:  http://bmagic.sourceforge.net

*/

#include "bm.h"
#include "bmfunc.h"
#include "bmdef.h"

#include <stdlib.h>
#include <algorithm>


namespace bm
{

/*!
    Class implements algorithm for random subset generation.

    Implemented method tries to be fair, but doesn't guarantee
    true randomeness.

    Performace note: 
    Class holds temporary buffers and variables, so it is recommended to
    re-use instances over multiple calls.

    \ingroup setalgo
*/
template<class BV>
class random_subset
{
public:
    random_subset();
    ~random_subset();

    /// Get random subset of input vector
    ///
    /// @param bv_out - destination vector
    /// @param bv_in  - input vector
    /// @param count  - number of bits to pick
    ///
    void sample(BV& bv_out, const BV& bv_in, unsigned count);
    

private:
    typedef 
        typename BV::blocks_manager_type  blocks_manager_type;

private:
    void get_subset(BV& bv_out, 
                    const BV& bv_in, 
                    unsigned  bv_in_count,
                    unsigned count);

    unsigned find_max_block();
    
    void get_random_subset(bm::word_t*       blk_out, 
                           const bm::word_t* blk_src,
                           unsigned          count);
    static 
    unsigned process_word(bm::word_t*       blk_out, 
                          const bm::word_t* blk_src,
                          unsigned          i,
                          unsigned          count);

    static
    void get_random_array(bm::word_t*       blk_out, 
                          bm::gap_word_t*   bit_list,
                          unsigned          bit_list_size,
                          unsigned          count);


private:
    random_subset(const random_subset&);
    random_subset& operator=(const random_subset&);
private:
    unsigned*         block_counts_; 
    unsigned short*   block_bits_take_;
    unsigned          blocks_;
    bm::gap_word_t    bit_list_[bm::gap_max_bits];
    unsigned*         block_candidates_;
    unsigned          candidates_count_;
    bm::word_t*       sub_block_;
};


///////////////////////////////////////////////////////////////////



template<class BV>
random_subset<BV>::random_subset()
: block_counts_(new unsigned[bm::set_total_blocks]),
  block_bits_take_(new bm::gap_word_t[bm::set_total_blocks]),
  block_candidates_(new unsigned[bm::set_total_blocks]),
  candidates_count_(0),
  sub_block_(new bm::word_t[bm::set_block_size])
{
}

template<class BV>
random_subset<BV>::~random_subset()
{
    delete [] block_counts_;
    delete [] block_bits_take_;
    delete [] block_candidates_;
    delete [] sub_block_;
}

template<class BV>
void random_subset<BV>::sample(BV&       bv_out, 
                               const BV& bv_in, 
                               unsigned  count)
{
    if (count == 0)
    {
        bv_out.clear(true);
        return;
    }

    unsigned bcnt = bv_in.count();
    if (count >= bcnt)
    {
        bv_out = bv_in;
        return;
    }
    if (count > bcnt/2) 
    {
        // build the complement vector and subtract it
        BV tmp_bv;
        unsigned delta_count = bcnt - count;

        get_subset(tmp_bv, bv_in, bcnt, delta_count);
        bv_out = bv_in;
        bv_out -= tmp_bv;
        return;
    }

    get_subset(bv_out, bv_in, bcnt, count);
    bv_out.forget_count();
}

template<class BV>
void random_subset<BV>::get_subset(BV&       bv_out, 
                                   const BV& bv_in, 
                                   unsigned  bcnt,
                                   unsigned  count)
{
    bv_out.clear(true);
    bv_out.resize(bv_in.size());

    const blocks_manager_type& bman_in = bv_in.get_blocks_manager();
    blocks_manager_type& bman_out = bv_out.get_blocks_manager();

    bm::word_t* tmp_block = bman_out.check_allocate_tempblock();
    candidates_count_ = 0;


    // compute bit-counts in all blocks
    //
    unsigned block_count = blocks_ = bv_in.count_blocks(block_counts_);
    for (unsigned i = 0; i <= block_count; ++i)
    {
        if (block_counts_[i])
        {
            float block_percent = 
                ((float)(block_counts_[i]+1)) / (float)bcnt;
            float bits_to_take = ((float)count) * block_percent; 
            bits_to_take += 0.99f;

            unsigned t = (unsigned)bits_to_take;
            block_bits_take_[i] = (bm::gap_word_t)t;
            
            if (block_bits_take_[i] == 0)
            {
                block_bits_take_[i] = 1;
            }
            else
            if (block_bits_take_[i] > block_counts_[i])
                block_bits_take_[i] = (gap_word_t)block_counts_[i];

            BM_ASSERT(bman_in.get_block(i));
        }
        else
        {
            block_bits_take_[i] = 0;
        }
    } // for i

    for (unsigned take_count = 0; count; count -= take_count) 
    {
        unsigned i = find_max_block();
        take_count = block_bits_take_[i];
        if (take_count > count)
            take_count = count;
        if (take_count == 0)
            continue;
        const bm::word_t* blk_src = bman_in.get_block(i);
        BM_ASSERT(blk_src);

        // allocate target block
        bm::word_t* blk_out = bman_out.get_block(i);
        if (blk_out != 0)
        {
            blk_out = bman_out.deoptimize_block(i);
        } 
        else
        {            
            blk_out = bman_out.get_allocator().alloc_bit_block();
            bman_out.set_block(i, blk_out);
        }
        if (take_count == block_counts_[i])
        {
            // copy the whole src block
            if (BM_IS_GAP(blk_src))
            {
                gap_convert_to_bitset(blk_out, BMGAP_PTR(blk_src));
            }
            else
            {
                bm::bit_block_copy(blk_out, blk_src);                
            }
            block_bits_take_[i] = 0; // exclude block from any farther picking
            continue;
        }
        bit_block_set(blk_out, 0);

        if (block_counts_[i] < 4096) // use array shuffle
        {
            unsigned arr_len;
            // convert source block to bit-block
            if (BM_IS_GAP(blk_src))
            {
                arr_len = gap_convert_to_arr(bit_list_,
                                             BMGAP_PTR(blk_src),
                                             bm::gap_max_bits);
            }
            else // bit-block
            {
                arr_len = bit_convert_to_arr(bit_list_, 
                                             blk_src, 
                                             bm::gap_max_bits, 
                                             bm::gap_max_bits,
                                             0);
            }
            BM_ASSERT(arr_len);
            get_random_array(blk_out, bit_list_, arr_len, take_count);
        }
        else // dense block
        {
            // convert source block to bit-block
            if (BM_IS_GAP(blk_src))
            {
                gap_convert_to_bitset(tmp_block, BMGAP_PTR(blk_src));
                blk_src = tmp_block;
            }

            // pick random bits source block to target
            get_random_subset(blk_out, blk_src, take_count);
        }
        
        block_bits_take_[i] = 0; // exclude block from any farther picking
    }
}

template<class BV>
void random_subset<BV>::get_random_subset(bm::word_t*       blk_out, 
                                          const bm::word_t* blk_src,
                                          unsigned          take_count)
{
    for (unsigned rounds = 0; take_count && rounds < 10; ++rounds)
    {
        // pick random scan start and scan direction
        //
        unsigned i = rand() % bm::set_block_size;
        unsigned new_count;

        for (; i < bm::set_block_size && take_count; ++i)
        {
            if (blk_src[i] && (blk_out[i] != blk_src[i]))
            {
                take_count -= new_count = 
                    process_word(blk_out, blk_src, i, take_count);
            }
        } // for i

    } // for
    // if masked scan did not produce enough results..
    //
    if (take_count)
    {
        // Find all vacant bits: do logical (src SUB out)
        for (unsigned i = 0; i < bm::set_block_size; ++i)
        {
            sub_block_[i] = blk_src[i] & ~blk_out[i];
        }
        // now transform vacant bits to array, then pick random elements
        //
        unsigned arr_len = bit_convert_to_arr(bit_list_, 
                                              sub_block_, 
                                              bm::gap_max_bits, 
                                              bm::gap_max_bits,
                                              0);
        BM_ASSERT(arr_len);
        get_random_array(blk_out, bit_list_, arr_len, take_count);        
    }
}

template<class BV>
unsigned random_subset<BV>::process_word(bm::word_t*       blk_out, 
                                         const bm::word_t* blk_src,
                                         unsigned          i,
                                         unsigned          count)
{
    unsigned new_bits, mask;

    do 
    {    
        mask = rand();
        mask ^= mask << 16;
    } while (mask == 0);

    unsigned src_rand = blk_src[i] & mask;    
    new_bits = src_rand & ~blk_out[i];

    if (new_bits)
    {
        unsigned new_count = bm::word_bitcount(new_bits);

        // check if we accidentally picked more bits than needed
        if (new_count > count)
        {
            BM_ASSERT(count);

            unsigned char blist[64];
            unsigned arr_size = bm::bit_list_4(new_bits, blist);
            BM_ASSERT(arr_size == new_count);
            std::random_shuffle(blist, blist + arr_size);
            unsigned value = 0;
            for (unsigned j = 0; j < count; ++j)
            {
                value |= (1 << blist[j]);
            }
            new_bits = value;
            new_count = count;

            BM_ASSERT(bm::word_bitcount(new_bits) == count);
            BM_ASSERT((new_bits & ~blk_src[i]) == 0);
        }

        blk_out[i] |= new_bits;
        return new_count;
    }

    return 0;    
}


template<class BV>
void random_subset<BV>::get_random_array(bm::word_t*       blk_out, 
                                         bm::gap_word_t*   bit_list,
                                         unsigned          bit_list_size,
                                         unsigned          count)
{
    std::random_shuffle(bit_list, bit_list + bit_list_size);
    for (unsigned i = 0; i < count; ++i)
    {
        bm::set_bit(blk_out, bit_list[i]);
    }
}

template<class BV>
unsigned random_subset<BV>::find_max_block() 
{
    if (candidates_count_)
    {
        return block_candidates_[--candidates_count_];
    }

    unsigned candidate = 0;
    unsigned max_i = 0;
    for (unsigned i = 0; i <= blocks_; ++i)
    {
        if (block_bits_take_[i] == 0) continue;
        if (block_bits_take_[i] == candidate)
        {
            block_candidates_[candidates_count_++] = i;
        }
        else
        {
            unsigned diff = abs((int)block_bits_take_[i] - (int)candidate);
            double d = (double)diff / (double)candidate;

            if (d < 0.20f) // delta is statistically insignificant
            {
                block_candidates_[candidates_count_++] = i;
            }
            else
            if (block_bits_take_[i] > candidate)
            {
                 candidate = block_bits_take_[i];
                 max_i = i;
                 candidates_count_ = 0;
                 block_candidates_[candidates_count_++] = i;
            }
        }
    }

    if (candidates_count_ > 1)
    {
        std::random_shuffle(block_candidates_, block_candidates_ + candidates_count_);
        return find_max_block();
    }
    return max_i;
}


} // namespace

#include "bmundef.h"

#endif
