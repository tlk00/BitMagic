#ifndef BM_BLOCKS__H__INCLUDED__
#define BM_BLOCKS__H__INCLUDED__
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


#include "bmfwd.h"

#ifdef _MSC_VER
#pragma warning( push )
#pragma warning( disable : 4100)
#endif


namespace bm
{

/*!
   @brief bitvector blocks manager
        Embedded class managing bit-blocks on very low level.
        Includes number of functor classes used in different bitset algorithms. 
   @ingroup bvector
   @internal
*/

template<class Alloc>
class blocks_manager
{
public:
    template<typename TAlloc> friend class bvector;

    typedef Alloc allocator_type;

    /** Base functor class (block visitor)*/
    class bm_func_base
    {
    public:
        bm_func_base(blocks_manager& bman) : bm_(bman) {}

        void on_empty_top(unsigned /* top_block_idx*/ ) {}
        void on_empty_block(unsigned /* block_idx*/ ) {}
    private:
        bm_func_base(const bm_func_base&);
        bm_func_base& operator=(const bm_func_base&);
    protected:
        blocks_manager&  bm_;
    };


    /** Base functor class connected for "constant" functors*/
    class bm_func_base_const
    {
    public:
        bm_func_base_const(const blocks_manager& bman) : bm_(bman) {}

        void on_empty_top(unsigned /* top_block_idx*/ ) {}
        void on_empty_block(unsigned /* block_idx*/ ) {}
    private:
        bm_func_base_const(const bm_func_base_const&);
        bm_func_base_const& operator=(const bm_func_base_const&);
    protected:
        const blocks_manager&  bm_;
    };


    /** Base class for bitcounting functors */
    class block_count_base : public bm_func_base_const
    {
    protected:
        block_count_base(const blocks_manager& bm) 
            : bm_func_base_const(bm) {}

        bm::id_t block_count(const bm::word_t* block) const
        {
            return this->bm_.block_bitcount(block);
        }
    };


    /** Bitcounting functor */
    class block_count_func : public block_count_base
    {
    public:
        block_count_func(const blocks_manager& bm) 
            : block_count_base(bm), count_(0) {}

        bm::id_t count() const { return count_; }

        void operator()(const bm::word_t* block)
        {
            count_ += this->block_count(block);
        }

    private:
        bm::id_t count_;
    };


    /** Bitcounting functor filling the block counts array*/
    class block_count_arr_func : public block_count_base
    {
    public:
        block_count_arr_func(const blocks_manager& bm, unsigned* arr) 
            : block_count_base(bm), arr_(arr), last_idx_(0) 
        {
            arr_[0] = 0;
        }

        void operator()(const bm::word_t* block, unsigned idx)
        {
            while (++last_idx_ < idx)
            {
                arr_[last_idx_] = 0;
            }
            arr_[idx] = this->block_count(block);
            last_idx_ = idx;
        }

        unsigned last_block() const { return last_idx_; }

    private:
        unsigned*  arr_;
        unsigned   last_idx_;
    };

    /** bit value change counting functor */
    class block_count_change_func : public bm_func_base_const
    {
    public:
        block_count_change_func(const blocks_manager& bm) 
            : bm_func_base_const(bm),
                count_(0),
                prev_block_border_bit_(0)
        {}

        bm::id_t block_count(const bm::word_t* block, unsigned idx)
        {
            bm::id_t cnt = 0;
            bm::id_t first_bit;
            
            if (IS_FULL_BLOCK(block) || (block == 0))
            {
                cnt = 1;
                if (idx)
                {
                    first_bit = block ? 1 : 0;
                    cnt -= !(prev_block_border_bit_ ^ first_bit);
                }
                prev_block_border_bit_ = block ? 1 : 0;
            }
            else
            {
                if (BM_IS_GAP(block))
                {
                    gap_word_t* gap_block = BMGAP_PTR(block);
                    cnt = gap_length(gap_block) - 1;
                    if (idx)
                    {
                        first_bit = bm::gap_test_unr(gap_block, 0);
                        cnt -= !(prev_block_border_bit_ ^ first_bit);
                    }
                        
                    prev_block_border_bit_ = 
                        bm::gap_test_unr(gap_block, gap_max_bits-1);
                }
                else // bitset
                {
                    unsigned bit_count;
                    cnt = bit_block_calc_count_change(block, block + bm::set_block_size, &bit_count);
                    if (idx)
                    {
                        first_bit = block[0] & 1;
                        cnt -= !(prev_block_border_bit_ ^ first_bit);
                    }
                    prev_block_border_bit_ = 
                        block[set_block_size-1] >> ((sizeof(block[0]) * 8) - 1);
                    
                }
            }
            return cnt;
        }
        
        bm::id_t count() const { return count_; }

        void operator()(const bm::word_t* block, unsigned idx)
        {
            count_ += block_count(block, idx);
        }

    private:
        bm::id_t   count_;
        bm::id_t   prev_block_border_bit_;
    };


    /** Functor detects if any bit set*/
    class block_any_func : public bm_func_base_const
    {
    public:
        block_any_func(const blocks_manager& bm) 
            : bm_func_base_const(bm) 
        {}

        bool operator()(const bm::word_t* block, unsigned /*idx*/)
        {
            if (BM_IS_GAP(block)) // gap block
                return (!gap_is_all_zero(BMGAP_PTR(block)));
            if (IS_FULL_BLOCK(block)) 
                return true;
            return (!bit_is_all_zero(block));
        }
    };

    /*! Change GAP level lengths functor */
    class gap_level_func : public bm_func_base
    {
    public:
        gap_level_func(blocks_manager& bm, const gap_word_t* glevel_len)
            : bm_func_base(bm),
                glevel_len_(glevel_len)
        {
            BM_ASSERT(glevel_len);
        }

        void operator()(bm::word_t* block, unsigned idx)
        {
            blocks_manager& bman = this->bm_;
            
            if (!BM_IS_GAP(block))
                return;

            gap_word_t* gap_blk = BMGAP_PTR(block);

            // TODO: Use the same code as in the optimize functor
            if (gap_is_all_zero(gap_blk))
            {
                bman.set_block_ptr(idx, 0);
                bman.get_allocator().free_gap_block(gap_blk,
                                                    bman.glen());
            }
            else 
            if (gap_is_all_one(gap_blk, bm::gap_max_bits))
            {
                bman.set_block_ptr(idx, FULL_BLOCK_FAKE_ADDR);
                bman.get_allocator().free_gap_block(gap_blk,
                                                    bman.glen());
                return;
            }

            unsigned len = gap_length(gap_blk);
            int level = gap_calc_level(len, glevel_len_);
            if (level == -1)
            {
                bm::word_t* blk = 
                    bman.get_allocator().alloc_bit_block();
                bman.set_block_ptr(idx, blk);
                bm::gap_convert_to_bitset(blk, gap_blk);
            }
            else
            {
                gap_word_t* gap_blk_new = 
                    bman.allocate_gap_block(unsigned(level), gap_blk, glevel_len_);

                bm::word_t* p = (bm::word_t*) gap_blk_new;
                BMSET_PTRGAP(p);
                bman.set_block_ptr(idx, p);
            }
            bman.get_allocator().free_gap_block(gap_blk, bman.glen());
        }

    private:
        const gap_word_t* glevel_len_;
    };


    /*! Bitblock optimization functor */
    class block_opt_func : public bm_func_base
    {
    public:
        block_opt_func(blocks_manager& bm, 
                       bm::word_t*     temp_block,
                       int             opt_mode,
                        bv_statistics*  bv_stat=0) 
            : bm_func_base(bm),
              temp_block_(temp_block),
              opt_mode_(opt_mode),
              stat_(bv_stat),
              empty_(0)
        {
            BM_ASSERT(temp_block);
        }

        void on_empty_top(unsigned i)
        {
            BM_ASSERT(this->bm_.is_init());
            bm::word_t*** blk_root = this->bm_.top_blocks_root();
            bm::word_t** blk_blk = blk_root[i];
            if (blk_blk) 
            {
                this->bm_.get_allocator().free_ptr(blk_blk);
                blk_root[i] = 0;
            }
            if (stat_)
            {
                stat_->max_serialize_mem += (unsigned)(sizeof(unsigned) + 1);
            }
        }
        void on_empty_block(unsigned /* block_idx*/ ) { ++empty_; }

        void operator()(bm::word_t* block, unsigned idx)
        {
            blocks_manager& bman = this->bm_;
            if (IS_FULL_BLOCK(block)) 
            {
                ++empty_;
                return;
            }

            if (stat_) 
            {
                stat_->max_serialize_mem += empty_ << 2;
                empty_ = 0;
            }

            gap_word_t* gap_blk;
            if (BM_IS_GAP(block)) // gap block
            {
                gap_blk = BMGAP_PTR(block);
                // check if block is empty (or all 1)
                if (gap_is_all_zero(gap_blk))
                {
                    bman.set_block_ptr(idx, 0);
                    this->free_block(gap_blk);
                    ++empty_;
                }
                else 
                if (gap_is_all_one(gap_blk, bm::gap_max_bits))
                {
                    bman.set_block_ptr(idx, FULL_BLOCK_FAKE_ADDR);
                    this->free_block(gap_blk);
                    ++empty_;
                }
                else
                {
                    // regular gap block - compute statistics
                    if (stat_)
                    {
                        stat_->add_gap_block(
                                    bm::gap_capacity(gap_blk, bman.glen()),
                                    bm::gap_length(gap_blk));
                    }
                }
            }
            else // bit block
            {
                if (opt_mode_ < 3) // free_01 optimization
                {  
                    bm::wordop_t* blk1 = (wordop_t*)block;
                    bm::wordop_t* blk2 = (wordop_t*)(block + bm::set_block_size);

                    bool b = bm::bit_is_all_zero((bm::word_t*)blk1);
                    if (b)
                    {
                        bman.get_allocator().free_bit_block(block);
                        bman.set_block_ptr(idx, 0);
                        ++empty_;
                    } 
                    else
                    if (opt_mode_ == 2) // check if it is all 1 block
                    {
                        b = bm::is_bits_one(blk1, blk2);
                        if (b) 
                        {
                            bman.get_allocator().free_bit_block(block);
                            bman.set_block_ptr(idx, FULL_BLOCK_FAKE_ADDR);
                            ++empty_;
                        }
                    }

                    if (!b && stat_)
                        stat_->add_bit_block();

                    return;
                }
            
                // try to compress
            
                gap_word_t* tmp_gap_blk = (gap_word_t*)temp_block_;
                *tmp_gap_blk = bm::gap_max_level << 1;

                unsigned threashold = bman.glen(bm::gap_max_level)-4;

                unsigned len = bit_convert_to_gap(tmp_gap_blk, 
                                                    block, 
                                                    bm::gap_max_bits, 
                                                    threashold);
                if (len)    // compression successful                
                {                
                    bman.get_allocator().free_bit_block(block);

                    // check if new gap block can be eliminated
                    if (gap_is_all_zero(tmp_gap_blk))
                    {
                        bman.set_block_ptr(idx, 0);
                        ++empty_;
                    }
                    else if (gap_is_all_one(tmp_gap_blk, bm::gap_max_bits))
                    {
                        bman.set_block_ptr(idx, FULL_BLOCK_FAKE_ADDR);
                        ++empty_;
                    }
                    else
                    {
                        int level = bm::gap_calc_level(len, bman.glen());
                        BM_ASSERT(level >= 0);
                        gap_blk = 
                            bman.allocate_gap_block(unsigned(level), tmp_gap_blk);
                        bman.set_block_gap_ptr(idx, gap_blk);
                        if (stat_)
                        {
                            stat_->add_gap_block(
                                    bm::gap_capacity(gap_blk, bman.glen()),
                                    bm::gap_length(gap_blk));
                        }
                    }
                }  
                else  // non-compressable bit-block
                {
                    if (stat_)
                        stat_->add_bit_block();
                }
            }
        }
    private:
        void free_block(gap_word_t* gap_blk)
        {
            this->bm_.get_allocator().free_gap_block(gap_blk,
                                                     this->bm_.glen());
        }

    private:
        bm::word_t*         temp_block_;
        int                 opt_mode_;
        bv_statistics*      stat_;
        unsigned            empty_;
    };

    /** Bitblock invert functor */
    class block_invert_func : public bm_func_base
    {
    public:
        block_invert_func(blocks_manager& bm) 
            : bm_func_base(bm) {}

        void operator()(bm::word_t* block, unsigned idx)
        {
            if (!block)
                this->bm_.set_block(idx, FULL_BLOCK_FAKE_ADDR);
            else
            if (IS_FULL_BLOCK(block))
                this->bm_.set_block_ptr(idx, 0);
            else
            {
                if (BM_IS_GAP(block)) // gap block
                    gap_invert(BMGAP_PTR(block));
                else  // bit block
                {
                    bm::wordop_t* wrd_ptr = (wordop_t*) block;
                    bm::wordop_t* wrd_end = 
                            (wordop_t*) (block + bm::set_block_size);
                    bm::bit_invert(wrd_ptr, wrd_end);
                }
            }

        }
    };

    /** Set block zero functor */
    class block_zero_func : public bm_func_base
    {
    public:
        block_zero_func(blocks_manager& bm) 
        : bm_func_base(bm), alloc_(bm.get_allocator())
        {}

        void operator()(bm::word_t* block, unsigned idx)
        {
            if (BM_IS_GAP(block))
                bm::gap_set_all(BMGAP_PTR(block), bm::gap_max_bits, 0);
            else  // BIT block
            {
                if (IS_FULL_BLOCK(block))
                    this->bm_.set_block_ptr(idx, 0);
                else
                {
                    BM_ASSERT(block);
                    {
                        // faster to free and re-aquire than set to zero
                        alloc_.free_bit_block(block);
                        this->bm_.set_block_ptr(idx, 0);
                    }
                }
            }
        }
    private:
        allocator_type& alloc_;
    };

    /** Fill block with all-one bits functor */
    class block_one_func : public bm_func_base
    {
    public:
        block_one_func(blocks_manager& bm) : bm_func_base(bm) {}

        void operator()(bm::word_t* block, unsigned idx)
        {
            if (!IS_FULL_BLOCK(block))
                this->bm_.set_block_all_set(idx);
        }
    };

public:
    blocks_manager()
    : max_bits_(bm::id_max),
      top_blocks_(0),
      temp_block_(0),
      alloc_(Alloc())
    {
        ::memcpy(glevel_len_, bm::gap_len_table<true>::_len, sizeof(glevel_len_));
        top_block_size_ = 1;
    }

    blocks_manager(const gap_word_t* glevel_len, 
                    bm::id_t          max_bits,
                    const Alloc&      alloc = Alloc())
        : max_bits_(max_bits),
          top_blocks_(0),
          temp_block_(0),
          alloc_(alloc)
    {
        ::memcpy(glevel_len_, glevel_len, sizeof(glevel_len_));
        top_block_size_ = 1;
    }

    blocks_manager(const blocks_manager& blockman)
        : max_bits_(blockman.max_bits_),
          top_blocks_(0),
          top_block_size_(blockman.top_block_size_),
        #ifdef BM_DISBALE_BIT_IN_PTR
            gap_flags_(blockman.gap_flags_),
        #endif
            temp_block_(0),
            alloc_(blockman.alloc_)
    {
        ::memcpy(glevel_len_, blockman.glevel_len_, sizeof(glevel_len_));

        if (blockman.is_init())
        {
            reserve_top_blocks(blockman.top_block_size());
            this->copy(blockman);
        }
    }
    
#ifndef BM_NO_CXX11
    blocks_manager(blocks_manager&& blockman) BMNOEXEPT
        : max_bits_(blockman.max_bits_),
          top_blocks_(0),
          top_block_size_(blockman.top_block_size_),
          temp_block_(0),
          alloc_(blockman.alloc_)
    {
        ::memcpy(glevel_len_, blockman.glevel_len_, sizeof(glevel_len_));
        move_from(blockman);
    }
#endif

    ~blocks_manager() BMNOEXEPT
    {
        if (temp_block_)
            alloc_.free_bit_block(temp_block_);
        destroy_tree();
    }
    
    /*! \brief Swaps content 
        \param bm  another blocks manager
    */
    void swap(blocks_manager& bm) BMNOEXEPT
    {
        BM_ASSERT(this != &bm);

        word_t*** btmp = top_blocks_;
        top_blocks_ = bm.top_blocks_;
        bm.top_blocks_ = btmp;

        bm::xor_swap(this->max_bits_, bm.max_bits_);
        bm::xor_swap(this->top_block_size_, bm.top_block_size_);

        BM_ASSERT(sizeof(glevel_len_) / sizeof(glevel_len_[0]) == bm::gap_levels); // paranoiya check
        for (unsigned i = 0; i < bm::gap_levels; ++i)
        {
            bm::xor_swap(glevel_len_[i], bm.glevel_len_[i]);
        }
    }
    
    /*! \brief implementation of moving semantics
    */
    void move_from(blocks_manager& bm) BMNOEXEPT
    {
        deinit_tree();
        swap(bm);
        alloc_ = bm.alloc_;
        if (!temp_block_)  // this does not have temp_block, borrow it from the donor
        {
            temp_block_ = bm.temp_block_;
            bm.temp_block_ = 0;
        }
    }
    

    void free_ptr(bm::word_t** ptr)
    {
        if (ptr) alloc_.free_ptr(ptr);
    }

    /**
        \brief Compute size of the block array needed to store bits
        \param bits_to_store - supposed capacity (number of bits)
        \return size of the top level block
    */
    unsigned compute_top_block_size(bm::id_t bits_to_store)
    {
        if (bits_to_store == bm::id_max)  // working in full-range mode
            return bm::set_array_size;

        unsigned top_block_sz = (unsigned)
            (bits_to_store / (bm::set_block_size * sizeof(bm::word_t) *
                                                bm::set_array_size * 8));
        if (top_block_sz < bm::set_array_size) ++top_block_sz;
        return top_block_sz;
    }

    /**
        Returns current capacity (bits)
    */
    bm::id_t capacity() const
    {
        // arithmetic overflow protection...
        return top_block_size_ == bm::set_array_size ? bm::id_max :
            top_block_size_ * bm::set_array_size * bm::bits_in_block;
    }

    /**
        \brief Finds block in 2-level blocks array  
        \param nb - Index of block (logical linear number)
        \return block adress or NULL if not yet allocated
    */
    bm::word_t* get_block(unsigned nb) const
    {
        if (!top_blocks_)
            return 0;
        unsigned block_idx = nb >> bm::set_array_shift;
        if (block_idx >= top_block_size_)
        {
            return 0;
        }
        bm::word_t** blk_blk = top_blocks_[block_idx];
        bm::word_t* ret = blk_blk ? blk_blk[nb & bm::set_array_mask] : 0;
        if (ret == FULL_BLOCK_FAKE_ADDR)
            ret = FULL_BLOCK_REAL_ADDR;
        return ret;
    }

    /**
    \brief Finds const block in 2-level blocks array (returns unsanitized address!)
    \param nb - Index of block (logical linear number)
    \return block adress or NULL if not yet allocated or FULL_BLOCK_FAKE_ADDR
    */
    bm::word_t* get_block_ptr(unsigned nb) const
    {
        unsigned block_idx = nb >> bm::set_array_shift;
        if (!top_blocks_ || (block_idx >= top_block_size_))
            return 0;
        bm::word_t** blk_blk = top_blocks_[block_idx];
        return blk_blk ? blk_blk[nb & bm::set_array_mask] : 0;
    }

    /**
        \brief Finds block in 2-level blocks array  
        Specilized version of get_block(unsigned), returns an additional
        condition when there are no more blocks

        \param nb - Index of block (logical linear number)
        \param no_more_blocks - 1 if there are no more blocks at all
        \return block adress or NULL if not yet allocated
    */
    bm::word_t* get_block(unsigned nb, int* no_more_blocks) const
    {
        BM_ASSERT(top_blocks_);
        unsigned block_idx = nb >> bm::set_array_shift;
        if (block_idx >= top_block_size_)
        {
            *no_more_blocks = 1;
            return 0;
        }
        *no_more_blocks = 0;
        bm::word_t** blk_blk = top_blocks_[block_idx];
        bm::word_t* ret = blk_blk ? blk_blk[nb & bm::set_array_mask] : 0;
        if (ret == FULL_BLOCK_FAKE_ADDR)
            ret = FULL_BLOCK_REAL_ADDR;
        return ret;
    }

    /** 
    Recalculate absolute block address into coordinates
    */
    static
    BMFORCEINLINE
    void get_block_coord(unsigned nb, unsigned& i, unsigned& j)
    {
        i = nb >> bm::set_array_shift; // top block address
        j = nb &  bm::set_array_mask;  // address in sub-block
    }

    /**
    Find the next non-zero block starting from nb
    \param nb - block index
    \param deep_scan - flag to perform detailed bit-block analysis
    @return bm::set_total_blocks - no more blocks
    */
    unsigned find_next_nz_block(unsigned nb, bool deep_scan = true) const
    {
        if (is_init())
        {
            unsigned i,j;
            get_block_coord(nb, i, j);
            for (;i < top_block_size_; ++i)
            { 
                bm::word_t** blk_blk = top_blocks_[i];
                if (!blk_blk)
                { 
                    nb += bm::set_array_size - j;
                }
                else
                   for (;j < bm::set_array_size; ++j, ++nb)
                   {
                       bm::word_t* blk = blk_blk[j];
                       if (blk && !bm::check_block_zero(blk, deep_scan))
                           return nb;
                   } // for j
                j = 0;
            } // for i
        } // is_init()

        return bm::set_total_blocks;
    }

    /**
        \brief Finds block in 2-level blocks array
        \param i - top level block index
        \param j - second level block index
        \return block adress or NULL if not yet allocated
    */
    const bm::word_t* get_block(unsigned i, unsigned j) const
    {
        if (!top_blocks_ || i >= top_block_size_) return 0;

        const bm::word_t* const* blk_blk = top_blocks_[i];
        const bm::word_t* ret = (blk_blk == 0) ? 0 : blk_blk[j];
        return (ret == FULL_BLOCK_FAKE_ADDR) ? FULL_BLOCK_REAL_ADDR : ret;
    }

    /**
        \brief Finds block in 2-level blocks array (unsinitized)
        \param i - top level block index
        \param j - second level block index
        \return block adress or NULL if not yet allocated
    */
    const bm::word_t* get_block_ptr(unsigned i, unsigned j) const
    {
        if (!top_blocks_ || i >= top_block_size_) return 0;

        const bm::word_t* const* blk_blk = top_blocks_[i];
        const bm::word_t* ret = (blk_blk == 0) ? 0 : blk_blk[j];
        return ret;
    }
    /**
        \brief Finds block in 2-level blocks array (unsinitized)
        \param i - top level block index
        \param j - second level block index
        \return block adress or NULL if not yet allocated
    */
    bm::word_t* get_block_ptr(unsigned i, unsigned j)
    {
        if (!top_blocks_ || i >= top_block_size_) return 0;

        bm::word_t* const* blk_blk = top_blocks_[i];
        bm::word_t* ret = (blk_blk == 0) ? 0 : blk_blk[j];
        return ret;
    }


    /**
        \brief Function returns top-level block in 2-level blocks array
        \param i - top level block index
        \return block adress or NULL if not yet allocated
    */
    const bm::word_t* const * get_topblock(unsigned i) const
    {
        return (!top_blocks_ || i >= top_block_size_) ? 0 : top_blocks_[i];
    }

    /** 
        \brief Returns root block in the tree.
    */
    bm::word_t*** top_blocks_root() const
    {
        blocks_manager* bm = 
            const_cast<blocks_manager*>(this);
        return bm->top_blocks_root();
    }

    void set_block_all_set(unsigned nb)
    {
        bm::word_t* block = this->get_block(nb);
        set_block(nb, const_cast<bm::word_t*>(FULL_BLOCK_FAKE_ADDR));

        if (BM_IS_GAP(block))
            alloc_.free_gap_block(BMGAP_PTR(block), glevel_len_);
        else
        {
            if (IS_VALID_ADDR(block))
                alloc_.free_bit_block(block);
        }
    }

    /**
        Create(allocate) bit block. Old block (if exists) gets deleted.
    */
    bm::word_t* alloc_bit_block(unsigned nb)
    {
        bm::word_t* block = this->get_allocator().alloc_bit_block();
        bm::word_t* old_block = set_block(nb, block);
        if (IS_VALID_ADDR(old_block))
        {
            if (BM_IS_GAP(old_block))
                alloc_.free_gap_block(BMGAP_PTR(old_block), glen());
            else
                alloc_.free_bit_block(old_block);
        }
        return block;
    }

    /**
        Create all-zeros bit block. Old block (if exists) gets deleted.
    */
    bm::word_t* make_bit_block(unsigned nb)
    {
        bm::word_t* block = this->alloc_bit_block(nb);
        bit_block_set(block, 0);
        return block;
    }

    /**
        Create bit block as a copy of source block (bit or gap).
        Old block (if exists) gets deleted.
    */
    bm::word_t* copy_bit_block(unsigned          nb, 
                               const bm::word_t* block_src, int is_src_gap)
    {
        if (block_src == 0)
        {
            zero_block(nb);
            return 0;
        }
        bm::word_t* block = alloc_bit_block(nb);
        if (is_src_gap)
        {
            gap_word_t* gap_block = BMGAP_PTR(block_src);
            gap_convert_to_bitset(block, gap_block);
        }
        else
        {            
            bit_block_copy(block, block_src);
        }
        return block;
    }
    

    /**
        Clone GAP block from another GAP
        It can mutate into a bit-block if does not fit
    */
    bm::word_t* clone_gap_block(const bm::gap_word_t* gap_block, bool& gap_res)
    {
        BM_ASSERT(gap_block);
        
        bm::word_t* new_block;
        unsigned len = bm::gap_length(gap_block);
        int new_level = bm::gap_calc_level(len, glen());
        if (new_level < 0) // mutation
        {
            gap_res = false;
            new_block = alloc_.alloc_bit_block();
            bm::gap_convert_to_bitset(new_block, gap_block);
        }
        else
        {
            gap_res = true;
            new_block = (bm::word_t*)
                    get_allocator().alloc_gap_block(unsigned(new_level), glen());
            ::memcpy(new_block, gap_block, len * sizeof(bm::gap_word_t));
            bm::set_gap_level(new_block, new_level);
        }
        return new_block;
    }
    
    
    /**
        Copy block from another vector.
        Note:Target block is always replaced through re-allocation.
    */
    bm::word_t* copy_block(unsigned idx, const blocks_manager& bm_src)
    {
        const bm::word_t* block = bm_src.get_block(idx);
        if (block == 0)
        {
            zero_block(idx);
            return 0;
        }
        bm::word_t* new_blk = 0;
        bool is_gap = BM_IS_GAP(block);

        if (is_gap)
        {
            bm::gap_word_t* gap_block = BMGAP_PTR(block);
            new_blk = clone_gap_block(gap_block, is_gap); // is_gap - output
        }
        else
        {
            if (IS_FULL_BLOCK(block))
            {
                new_blk = FULL_BLOCK_FAKE_ADDR;
            }
            else
            {
                new_blk = alloc_.alloc_bit_block();
                bm::bit_block_copy(new_blk, block);
            }
        }
        set_block(idx, new_blk, is_gap);
        return new_blk;
    }


    /** 
        Function checks if block is not yet allocated, allocates it and sets to
        all-zero or all-one bits. 

        If content_flag == 1 (ALLSET block) requested and taken block is already ALLSET,
        function will return NULL

        initial_block_type and actual_block_type : 0 - bitset, 1 - gap
    */
    bm::word_t* check_allocate_block(unsigned nb,
                                     unsigned content_flag,
                                     int      initial_block_type,
                                     int*     actual_block_type,
                                     bool     allow_null_ret=true)
    {
        bm::word_t* block = this->get_block_ptr(nb);

        if (!IS_VALID_ADDR(block)) // NULL block or ALLSET
        {
            // if we wanted ALLSET and requested block is ALLSET return NULL
            unsigned block_flag = IS_FULL_BLOCK(block);
            *actual_block_type = initial_block_type;
            if (block_flag == content_flag && allow_null_ret)
            {
                return 0; // it means nothing to do for the caller
            }

            if (initial_block_type == 0) // bitset requested
            {
                block = alloc_.alloc_bit_block();
                // initialize block depending on its previous status
                bit_block_set(block, block_flag ? 0xFF : 0);
                set_block(nb, block);
            }
            else // gap block requested
            {
                bm::gap_word_t* gap_block = allocate_gap_block(0);
                gap_set_all(gap_block, bm::gap_max_bits, block_flag);
                set_block(nb, (bm::word_t*)gap_block, true/*gap*/);
                return (bm::word_t*)gap_block;
            }
        }
        else // block already exists
        {
            *actual_block_type = BM_IS_GAP(block);
        }

        return block;
    }

    /**
        Function checks if block is not yet allocated, allocates and returns
    */
    bm::word_t* check_allocate_block(unsigned nb, int initial_block_type)
    {
        bm::word_t* block = this->get_block_ptr(nb);

        if (!IS_VALID_ADDR(block)) // NULL block or ALLSET
        {
            // if we wanted ALLSET and requested block is ALLSET return NULL
            unsigned block_flag = IS_FULL_BLOCK(block);
            if (initial_block_type == 0) // bitset requested
            {
                block = alloc_.alloc_bit_block();
                // initialize block depending on its previous status
                bit_block_set(block, block_flag ? 0xFF : 0);
                set_block(nb, block);
            }
            else // gap block requested
            {
                bm::gap_word_t* gap_block = allocate_gap_block(0);
                gap_set_all(gap_block, bm::gap_max_bits, block_flag);
                block = (bm::word_t*)gap_block;
                set_block(nb, block, true/*gap*/);
                BMSET_PTRGAP(block);
            }
        }
        return block;
    }


    /*! @brief Fills all blocks with 0.
        @param free_mem - if true function frees the resources
    */
    void set_all_zero(bool free_mem)
    {
        if (!is_init()) return;
        
        unsigned top_size = this->top_block_size();
        if (free_mem)
        {
            deinit_tree(); // TODO: optimization of top-level realloc
        }
        else
        {
            block_zero_func zero_func(*this);
            for_each_nzblock(top_blocks_, top_size,  zero_func);
        }
    }

    /*! Replaces all blocks with ALL_ONE block.
    */
    void set_all_one()
    {
        if (!is_init())
            init_tree();
        block_one_func func(*this);
        for_each_block(top_blocks_, top_block_size_,
                                bm::set_array_size, func);
    }
    
    
    bm::word_t** alloc_top_subblock(unsigned nblk_blk)
    {
        BM_ASSERT(top_blocks_[nblk_blk] == 0);

        bm::word_t** p = (bm::word_t**)alloc_.alloc_ptr();
        ::memset(top_blocks_[nblk_blk] = p, 0,
            bm::set_array_size * sizeof(bm::word_t*));
        return p;
    }
    
    bm::word_t** check_alloc_top_subblock(unsigned nblk_blk)
    {
        if(top_blocks_[nblk_blk] == 0)
        {
            return alloc_top_subblock(nblk_blk);
        }
        return top_blocks_[nblk_blk];
    }


    /**
        Places new block into descriptors table, returns old block's address.
        Old block is NOT deleted.
    */
    bm::word_t* set_block(unsigned nb, bm::word_t* block)
    {
        bm::word_t* old_block;
        
        if (!is_init())
            init_tree();

        // never use real full block adress, it may be instantiated in another DLL 
        // (unsafe static instantiation on windows)
        if (block == FULL_BLOCK_REAL_ADDR)
            block = FULL_BLOCK_FAKE_ADDR;

        // top block index
        unsigned nblk_blk = nb >> bm::set_array_shift;
        reserve_top_blocks(nblk_blk+1);
        
        // If first level array not yet allocated, allocate it and
        // assign block to it
        if (top_blocks_[nblk_blk] == 0)
        {
            alloc_top_subblock(nblk_blk);
            /*
            top_blocks_[nblk_blk] = (bm::word_t**)alloc_.alloc_ptr();
            ::memset(top_blocks_[nblk_blk], 0,
                bm::set_array_size * sizeof(bm::word_t*));
            */
            old_block = 0;
        }
        else
        {
            old_block = top_blocks_[nblk_blk][nb & bm::set_array_mask];
        }

        // NOTE: block will be replaced without freeing, potential memory leak?
        top_blocks_[nblk_blk][nb & bm::set_array_mask] = block;

        return old_block;
    }

    /**
    Allocate an place new GAP block (copy of provided block)
    */
    bm::word_t* set_gap_block(unsigned      nb,
                          const gap_word_t* gap_block_src,
                          int               level)
    {
        BM_ASSERT(top_blocks_);
        if (level < 0)
        {
            bm::word_t* blk = alloc_.alloc_bit_block();
            set_block(nb, blk);
            gap_convert_to_bitset(blk, gap_block_src);
            return blk;
        }
        else
        {
            gap_word_t* gap_blk = alloc_.alloc_gap_block(
                                                unsigned(level), this->glen());
            gap_word_t* gap_blk_ptr = BMGAP_PTR(gap_blk);
            ::memcpy(gap_blk_ptr, gap_block_src, 
                                  gap_length(gap_block_src) * sizeof(gap_word_t));
            set_gap_level(gap_blk_ptr, level);
            set_block(nb, (bm::word_t*)gap_blk, true /*GAP*/);
            return (bm::word_t*)gap_blk;
        }
    }

    /**
        Places new block into descriptors table, returns old block's address.
        Old block is not deleted.
    */
    bm::word_t* set_block(unsigned nb, bm::word_t* block, bool gap)
    {
        unsigned i, j;
        get_block_coord(nb, i, j);
        reserve_top_blocks(i + 1);
        
        return set_block(i, j, block, gap);
    }

    /**
        Places new block into descriptors table, returns old block's address.
        Old block is not deleted.
    */
    bm::word_t* set_block(unsigned i, unsigned j, bm::word_t* block, bool gap)
    {
        BM_ASSERT(i < top_block_size_);
     
        bm::word_t* old_block;
        if (block)
        {
            if (block == FULL_BLOCK_REAL_ADDR)
                block = FULL_BLOCK_FAKE_ADDR;
            else
                block =
                (bm::word_t*) (gap ? BMPTR_SETBIT0(block) : BMPTR_CLEARBIT0(block));
        }

        // If first level array not yet allocated, allocate it and
        // assign block to it
        if (!top_blocks_[i])
        {
            top_blocks_[i] = (bm::word_t**)alloc_.alloc_ptr();
            ::memset(top_blocks_[i], 0, bm::set_array_size * sizeof(void*));
            old_block = 0;
        }
        else
            old_block = top_blocks_[i][j];

        // NOTE: block will be replaced without freeing, potential memory leak?
        top_blocks_[i][j] = block;
        
        return old_block;
    }


    /**
        Places new block into blocks table.
    */
    BMFORCEINLINE
    void set_block_ptr(unsigned nb, bm::word_t* block)
    {
        BM_ASSERT((nb >> bm::set_array_shift) < top_block_size_);
        BM_ASSERT(is_init());
        BM_ASSERT(top_blocks_[nb >> bm::set_array_shift]);
        
        top_blocks_[nb >> bm::set_array_shift][nb & bm::set_array_mask] =
            (block == FULL_BLOCK_REAL_ADDR) ? FULL_BLOCK_FAKE_ADDR : block;
    }

    /**
        Places new block into blocks table.
    */
    BMFORCEINLINE
    void set_block_ptr(unsigned i, unsigned j, bm::word_t* block)
    {
        BM_ASSERT(is_init());
        BM_ASSERT(i < top_block_size_);
        BM_ASSERT(top_blocks_[i]);
        
        top_blocks_[i][j] =
            (block == FULL_BLOCK_REAL_ADDR) ? FULL_BLOCK_FAKE_ADDR : block;
    }

    /** 
        \brief Converts block from type gap to conventional bitset block.
        \param nb - Block's index. 
        \param gap_block - Pointer to the gap block, if NULL block nb is taken
        \return new bit block's memory
    */
    bm::word_t* convert_gap2bitset(unsigned nb, const gap_word_t* gap_block=0)
    {
        BM_ASSERT(is_init());
        
        unsigned i, j;
        get_block_coord(nb, i, j);
        reserve_top_blocks(i);
        
        return convert_gap2bitset(i, j, gap_block);
    }

    /**
        \brief Converts block from type gap to conventional bitset block.
        \param i - top index.
        \param j - secondary index.
        \param gap_block - Pointer to the gap block, if NULL block [i,j] is taken
        \return new bit block's memory
    */
    bm::word_t* convert_gap2bitset(unsigned i, unsigned j,
                                   const gap_word_t* gap_block=0)
    {
        BM_ASSERT(is_init());
        
        if (!top_blocks_[i])
        {
            alloc_top_subblock(i);
            /*
            top_blocks_[i] = (bm::word_t**)alloc_.alloc_ptr();
            ::memset(top_blocks_[i], 0, bm::set_array_size * sizeof(void*));
            */
        }
        bm::word_t* block = top_blocks_[i][j];
        gap_block = gap_block ? gap_block : BMGAP_PTR(block);

        BM_ASSERT(IS_VALID_ADDR((bm::word_t*)gap_block));

        bm::word_t* new_block = alloc_.alloc_bit_block();
        bm::gap_convert_to_bitset(new_block, gap_block);
        
        top_blocks_[i][j] = new_block;

        // new block will replace the old one(no deletion)
        if (block)
            alloc_.free_gap_block(BMGAP_PTR(block), glen());

        return new_block;
    }


    /**
        Make sure block turns into true bit-block if it is GAP or a full block
    */
    bm::word_t* deoptimize_block(unsigned nb)
    {
        bm::word_t* block = this->get_block(nb);
        if (BM_IS_GAP(block))
        {
            gap_word_t* gap_block = BMGAP_PTR(block);
            
            bm::word_t* new_block = alloc_.alloc_bit_block();
            gap_convert_to_bitset(new_block, gap_block);
            alloc_.free_gap_block(gap_block, this->glen());
            
            set_block_ptr(nb, new_block);
            return new_block;
        }
        if (IS_FULL_BLOCK(block)) 
        {
            bm::word_t* new_block = alloc_.alloc_bit_block();
            bm::bit_block_copy(new_block, FULL_BLOCK_REAL_ADDR);
            
            set_block_ptr(nb, new_block);
            return new_block;
        }
        return block;
    }

    /**
        Free block, make it zero pointer in the tree
    */
    void zero_block(unsigned nb)
    {
        unsigned i, j;
        get_block_coord(nb, i, j);
        if (!top_blocks_ || i >= top_block_size_)
            return;
        zero_block(i, j);
    }

    /**
    Free block, make it zero pointer in the tree
    */
    void zero_block(unsigned i, unsigned j)
    {
        BM_ASSERT(top_blocks_ && i < top_block_size_);
        
        bm::word_t** blk_blk = top_blocks_[i];
        if (blk_blk)
        {
            bm::word_t* block = blk_blk[j];
            blk_blk[j] = 0;

            if (BM_IS_GAP(block))
                alloc_.free_gap_block(BMGAP_PTR(block), glen());
            else
                if (IS_VALID_ADDR(block))
                    alloc_.free_bit_block(block);
        }
    }
    
    /**
    Free block, make it zero pointer in the tree
    */
    void zero_gap_block_ptr(unsigned i, unsigned j)
    {
        BM_ASSERT(top_blocks_ && i < top_block_size_);
        
        bm::word_t** blk_blk = top_blocks_[i];
        bm::word_t* block = blk_blk[j];

        BM_ASSERT(blk_blk);
        BM_ASSERT(BM_IS_GAP(block));

        blk_blk[j] = 0;
        alloc_.free_gap_block(BMGAP_PTR(block), glen());
    }


    /**
        Count number of bits ON in the block
    */
    static
    bm::id_t block_bitcount(const bm::word_t* block)
    {
        BM_ASSERT(block);
        id_t count;
        if (BM_IS_GAP(block))
        {
            count = gap_bit_count_unr(BMGAP_PTR(block));
        }
        else // bitset
        {
            count = (IS_FULL_BLOCK(block)) ? bm::bits_in_block
                                           : bit_block_count(block);
        }
        return count;
    }

    /**
        \brief Extends GAP block to the next level or converts it to bit block.
        \param nb - Block's linear index.
        \param blk - Blocks's pointer 

        \return new GAP block pointer or NULL if block type mutated
    */
    bm::gap_word_t* extend_gap_block(unsigned nb, gap_word_t* blk)
    {
        unsigned level = bm::gap_level(blk);
        unsigned len = bm::gap_length(blk);
        if (level == bm::gap_max_level || len >= gap_max_buff_len)
        {
            deoptimize_block(nb);
        }
        else
        {
            bm::gap_word_t* new_gap_blk = allocate_gap_block(++level, blk);
            bm::word_t* new_blk = (bm::word_t*)new_gap_blk;

            BMSET_PTRGAP(new_blk);

            set_block_ptr(nb, new_blk);
            alloc_.free_gap_block(blk, glen());

            return new_gap_blk;
        }
        return 0;
    }
    /**
        Mark pointer as GAP and assign to the blocks tree
    */
    void set_block_gap_ptr(unsigned nb, gap_word_t* gap_blk)
    {
        bm::word_t* block = (bm::word_t*)BMPTR_SETBIT0(gap_blk);
        set_block_ptr(nb, block);
    }


    /*! Returns temporary block, allocates if needed. */
    bm::word_t* check_allocate_tempblock()
    {
        return temp_block_ ? temp_block_ 
                            : (temp_block_ = alloc_.alloc_bit_block());
    }
    
    /*! deallocate temp block */
    void free_temp_block()
    {
        if (temp_block_)
        {
            alloc_.free_bit_block(temp_block_);
            temp_block_ = 0;
        }
    }

    /*! Assigns new GAP lengths vector */
    void set_glen(const gap_word_t* glevel_len)
    {
        ::memcpy(glevel_len_, glevel_len, sizeof(glevel_len_));
    }


    bm::gap_word_t* allocate_gap_block(unsigned level, 
                                       const gap_word_t* src = 0,
                                       const gap_word_t* glevel_len = 0)
    {
        if (!glevel_len)
            glevel_len = glevel_len_;
        gap_word_t* ptr = alloc_.alloc_gap_block(level, glevel_len);
        if (src)
        {
            unsigned len = gap_length(src);
            ::memcpy(ptr, src, len * sizeof(gap_word_t));
            // Reconstruct the mask word using the new level code.
            *ptr = (gap_word_t)(((len-1) << 3) | (level << 1) | (*src & 1));
        }
        else
        {
            *ptr = (gap_word_t)(level << 1);
        }
        return ptr;
    }

    unsigned mem_used() const
    {
        unsigned m_used = (unsigned)sizeof(*this);
        m_used += (unsigned)(temp_block_ ? sizeof(word_t) * bm::set_block_size : 0);
        m_used += (unsigned)(sizeof(bm::word_t**) * top_block_size_);

        #ifdef BM_DISBALE_BIT_IN_PTR
        m_used += (unsigned)(gap_flags_.mem_used() - sizeof(gap_flags_));
        #endif
        
        if (is_init())
        {
            for (unsigned i = 0; i < top_block_size_; ++i)
            {
                m_used += (unsigned)
                    (top_blocks_[i] ? sizeof(void*) * bm::set_array_size : 0);
            }
        }

        return m_used;
    }

    /** Returns true if second level block pointer is 0.
    */
    bool is_subblock_null(unsigned nsub) const
    {
        BM_ASSERT(top_blocks_);
        if (nsub >= top_block_size_)
            return true;
        return top_blocks_[nsub] == NULL;
    }

    bm::word_t*** top_blocks_root()
    {
        return top_blocks_;
    }

    /*! \brief Returns current GAP level vector
    */
    const gap_word_t* glen() const
    {
        return glevel_len_;
    }

    /*! \brief Returns GAP level length for specified level
        \param level - level number
    */
    unsigned glen(unsigned level) const
    {
        return glevel_len_[level];
    }
    
    /*! \brief Returns size of the top block array in the tree 
    */
    unsigned top_block_size() const
    {
        return top_block_size_;
    }

    /**
        \brief reserve capacity for specified number of bits
    */
    void reserve(unsigned max_bits)
    {
        if (max_bits) 
        {
            unsigned b = compute_top_block_size(max_bits);
            reserve_top_blocks(b);
        }
    }

    /*!
        \brief Make sure blocks manager has enough blocks capacity
    */
    unsigned reserve_top_blocks(unsigned top_blocks)
    {
        BM_ASSERT(top_blocks <= bm::set_array_size);
        //BM_ASSERT(is_init());

        if (top_blocks_ && top_blocks <= top_block_size_)
            return top_block_size_; // nothing to do
        
        bm::word_t*** new_blocks = 
            (bm::word_t***)alloc_.alloc_ptr(top_blocks);

        unsigned i = 0;
        if (top_blocks_)
        {
            for (; i < top_block_size_; ++i)
                new_blocks[i] = top_blocks_[i];
            alloc_.free_ptr(top_blocks_, top_block_size_);
        }
        for (; i < top_blocks; ++i)
            new_blocks[i] = 0;
        
        top_blocks_ = new_blocks;
        top_block_size_ = top_blocks;
        return top_block_size_;
    }
    
    /** \brief Returns reference on the allocator
    */
    allocator_type& get_allocator() { return alloc_; }

    /** \brief Returns allocator
    */
    allocator_type get_allocator() const { return alloc_; }
    
    
    /// if tree of blocks already up
    bool is_init() const { return top_blocks_ != 0; }
    
    /// allocate first level of descr. of blocks 
    void init_tree()
    {
        BM_ASSERT(top_blocks_ == 0);
        
        if (top_block_size_)
        {
            top_blocks_ = (bm::word_t***) alloc_.alloc_ptr(top_block_size_);
            ::memset(top_blocks_, 0, top_block_size_ * sizeof(bm::word_t**));
        }
        else
        {
            top_blocks_ = 0;
        }
    }

private:

    void operator =(const blocks_manager&);
    
    // ----------------------------------------------------------------
    #define BM_FREE_OP(x) blk = blk_blk[j + x]; \
        if (IS_VALID_ADDR(blk)) \
        { \
            if (BM_IS_GAP(blk)) \
                alloc_.free_gap_block(BMGAP_PTR(blk), glen()); \
            else \
                alloc_.free_bit_block(blk); \
        } 

    /** destroy tree, free memory in all blocks and control structures
        Note: pointers are NOT assigned to zero(!)
    */
    void destroy_tree() BMNOEXEPT
    {
        if (!top_blocks_) 
            return;

        unsigned top_blocks = top_block_size();
        for (unsigned i = 0; i < top_blocks; ++i)
        {
            bm::word_t** blk_blk = top_blocks_[i];
            if (!blk_blk) 
                continue;
            unsigned j = 0; bm::word_t* blk;
            do
            {
            #ifdef BM64_AVX2
                if (!avx2_test_all_zero_wave(blk_blk + j))
                {
                    BM_FREE_OP(0)
                    BM_FREE_OP(1)
                    BM_FREE_OP(2)
                    BM_FREE_OP(3)
                }
                j += 4;
            #elif defined(BM64_SSE4)
                if (!sse42_test_all_zero_wave(blk_blk + j))
                {
                    BM_FREE_OP(0)
                    BM_FREE_OP(1)
                }
                j += 2;
            #else
                BM_FREE_OP(0)
                ++j;
            #endif
            } while (j < bm::set_array_size);

            alloc_.free_ptr(top_blocks_[i]); // free second level
        } // for i

        alloc_.free_ptr(top_blocks_, top_block_size_); // free the top
    }
    #undef BM_FREE_OP 

    void deinit_tree() BMNOEXEPT
    {
        destroy_tree();
        top_blocks_ = 0; top_block_size_ = 0;
    }
    
    // ----------------------------------------------------------------
    
    void copy(const blocks_manager& blockman,
              unsigned block_from = 0, unsigned block_to = 65535)
    {
        unsigned arg_top_blocks = blockman.top_block_size();
        this->reserve_top_blocks(arg_top_blocks);
        
        bm::word_t*** blk_root = top_blocks_root();
        bm::word_t*** blk_root_arg = blockman.top_blocks_root();
        
        if (!blk_root_arg)
            return;
        
        unsigned i_from, j_from, i_to, j_to;
        get_block_coord(block_from, i_from, j_from);
        get_block_coord(block_to, i_to, j_to);
        
        if (i_to >= arg_top_blocks-1)
        {
            i_to = arg_top_blocks-1;
            j_to = bm::set_array_size-1;
        }

        for (unsigned i = i_from; i <= i_to; ++i)
        {
            bm::word_t** blk_blk_arg = blk_root_arg[i];
            if (!blk_blk_arg)
                continue;
            
            BM_ASSERT(blk_root[i] == 0);

            bm::word_t** blk_blk = blk_root[i] = (bm::word_t**)alloc_.alloc_ptr();
            ::memset(blk_blk, 0, bm::set_array_size * sizeof(bm::word_t*));
            
            unsigned j = (i == i_from) ? j_from : 0;
            unsigned j_limit = (i == i_to) ? j_to+1 : bm::set_array_size;
            bm::word_t* blk;
            const bm::word_t* blk_arg;
            do
            {
                blk = blk_blk[j]; blk_arg = blk_blk_arg[j];
                if (blk_arg)
                {
                    bool is_gap = BM_IS_GAP(blk_arg);
                    if (is_gap)
                    {
                        blk = clone_gap_block(BMGAP_PTR(blk_arg), is_gap);
                        if (is_gap)
                            BMSET_PTRGAP(blk);
                    }
                    else
                    {
                        if (blk_arg == FULL_BLOCK_FAKE_ADDR /*IS_FULL_BLOCK(blk_arg)*/)
                            blk = FULL_BLOCK_FAKE_ADDR;
                        else
                        {
                            BM_ASSERT(!IS_FULL_BLOCK(blk_arg));
                            blk = alloc_.alloc_bit_block();
                            bm::bit_block_copy(blk, blk_arg);
                        }
                    }
                    blk_blk[j] = blk;
                }
                ++j;
            } while (j < j_limit);
        } // for i
    }


private:
    /// maximum addresable bits
    bm::id_t                               max_bits_;
    /// Tree of blocks.
    bm::word_t***                          top_blocks_;
    /// Size of the top level block array in blocks_ tree
    unsigned                               top_block_size_;
    /// Temp block.
    bm::word_t*                            temp_block_; 
    /// vector defines gap block lengths for different levels 
    gap_word_t                             glevel_len_[bm::gap_levels];
    /// allocator
    allocator_type                         alloc_;
};

/**
    Bit block buffer guard
    \internal
*/
template<class BlocksManager>
class bit_block_guard
{
public:
    bit_block_guard(BlocksManager& bman, bm::word_t* blk=0) 
        : bman_(bman), 
          block_(blk)
    {}
    ~bit_block_guard()
    {
        if (IS_VALID_ADDR(block_))
            bman_.get_allocator().free_bit_block(block_, 3);
    }
    void attach(bm::word_t* blk)
    {
        if (IS_VALID_ADDR(block_))
            bman_.get_allocator().free_bit_block(block_);
        block_ = blk;
    }
    bm::word_t* allocate()
    {
        attach(bman_.get_allocator().alloc_bit_block(3));
        return block_;
    }
    bm::word_t* get() { return block_; }

private:
    bit_block_guard(const bit_block_guard&);
    bit_block_guard& operator=(const bit_block_guard&);
private:
    BlocksManager& bman_;
    bm::word_t*    block_;
};


}

#ifdef _MSC_VER
#pragma warning( pop )
#endif

#endif

