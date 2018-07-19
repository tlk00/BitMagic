#ifndef BM__H__INCLUDED__
#define BM__H__INCLUDED__
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

/*! \file bm.h
    \brief Compressed bit-vector bvector<> container, set algebraic methods, traversal iterators
*/


// define BM_NO_STL if you use BM in "STL free" environment and want
// to disable any references to STL headers
#ifndef BM_NO_STL
# include <iterator>
# include <initializer_list>
#endif

#include <limits.h>

#ifdef _MSC_VER
#pragma warning( push )
#pragma warning( disable : 4311 4312 4127)
#endif


#include "bmdef.h"
#include "bmconst.h"
#include "bmsimd.h"
#include "bmfwd.h"

# define BM_DECLARE_TEMP_BLOCK(x) bm::bit_block_t x;

#include "bmfunc.h"
#include "encoding.h"
#include "bmalloc.h"
#include "bmblocks.h"


extern "C"
{
    /**
    Callback type to visit (callback style) bits in bit-vector(s)
    
    @param handle_ptr - custom pointer to callback specific data
    @param bit_idx - number/index of visited bit
    
    @ingroup bvector
    */
    typedef int (*bit_visitor_callback_type)(void* handle_ptr, bm::id_t bit_idx);
}


namespace bm
{

/** @defgroup bmagic BitMagic Library
    BitMagic C++ Library
    For more information please visit: http://bitmagic.io
*/


/** @defgroup bvector bvector<> container
    The Main bvector<> Group
    bvector<> template: front end of the BitMagic library.
 
    @ingroup bmagic
*/

/** @defgroup bvit bvector<> iterators
    Iterators for compressed bit-vector traversal
    @ingroup bvector
*/


/*!
   @brief bitvector 
   Bit-vector container with runtime compression of bits
 
   @ingroup bvector
*/
template<class Alloc> 
class bvector
{
public:
    typedef Alloc                                        allocator_type;
    typedef typename allocator_type::allocator_pool_type allocator_pool_type;
    typedef blocks_manager<Alloc>                        blocks_manager_type;
    typedef bm::id_t                                     size_type; 

    /** Statistical information about bitset's memory allocation details. */
    struct statistics : public bv_statistics
    {};

    /**
        @brief Class reference implements an object for bit assignment.
        Since C++ does not provide with build-in bit type supporting l-value 
        operations we have to emulate it.

        @ingroup bvector
    */
    class reference
    {
    public:
        reference(bvector<Alloc>& bv, bm::id_t position) 
        : bv_(bv),
          position_(position)
        {}

        reference(const reference& ref)
        : bv_(ref.bv_), 
          position_(ref.position_)
        {
            bv_.set(position_, ref.bv_.get_bit(position_));
        }
        
        operator bool() const
        {
            return bv_.get_bit(position_);
        }

        const reference& operator=(const reference& ref) const
        {
            bv_.set(position_, (bool)ref);
            return *this;
        }

        const reference& operator=(bool value) const
        {
            bv_.set(position_, value);
            return *this;
        }

        bool operator==(const reference& ref) const
        {
            return bool(*this) == bool(ref);
        }

        /*! Bitwise AND. Performs operation: bit = bit AND value */
        const reference& operator&=(bool value) const
        {
            bv_.set_bit_and(position_, value);
            return *this;
        }

        /*! Bitwise OR. Performs operation: bit = bit OR value */
        const reference& operator|=(bool value) const
        {
            if (value != bv_.get_bit(position_))
            {
                bv_.set_bit(position_);
            }
            return *this;
        }

        /*! Bitwise exclusive-OR (XOR). Performs operation: bit = bit XOR value */
        const reference& operator^=(bool value) const
        {
            bv_.set(position_, value != bv_.get_bit(position_));
            return *this;
        }

        /*! Logical Not operator */
        bool operator!() const
        {
            return !bv_.get_bit(position_);
        }

        /*! Bit Not operator */
        bool operator~() const
        {
            return !bv_.get_bit(position_);
        }

        /*! Negates the bit value */
        reference& flip()
        {
            bv_.flip(position_);
            return *this;
        }

    private:
        bvector<Alloc>&   bv_;       //!< Reference variable on the parent.
        bm::id_t          position_; //!< Position in the parent bitvector.
    };

    typedef bool const_reference;

    /*!
        @brief Base class for all iterators.
        @ingroup bvit
    */
    class iterator_base
    {
    friend class bvector;
    public:
        iterator_base() : bv_(0), position_(bm::id_max), block_(0) {}

        bool operator==(const iterator_base& it) const
        {
            return (position_ == it.position_) && (bv_ == it.bv_);
        }

        bool operator!=(const iterator_base& it) const
        {
            return ! operator==(it);
        }

        bool operator < (const iterator_base& it) const
        {
            return position_ < it.position_;
        }

        bool operator <= (const iterator_base& it) const
        {
            return position_ <= it.position_;
        }

        bool operator > (const iterator_base& it) const
        {
            return position_ > it.position_;
        }

        bool operator >= (const iterator_base& it) const
        {
            return position_ >= it.position_;
        }

        /**
           \fn bool bm::bvector::iterator_base::valid() const
           \brief Checks if iterator is still valid. Analog of != 0 comparison for pointers.
           \returns true if iterator is valid.
        */
        bool valid() const
        {
            return position_ != bm::id_max;
        }

        /**
           \fn bool bm::bvector::iterator_base::invalidate() 
           \brief Turns iterator into an invalid state.
        */
        void invalidate()
        {
            position_ = bm::id_max;
        }
        
        /** \brief Compare FSMs for testing purposes
            \internal
        */
        bool compare_state(const iterator_base& ib) const
        {
            if (this->bv_ != ib.bv_)                 return false;
            if (this->position_ != ib.position_)     return false;
            if (this->block_ != ib.block_)           return false;
            if (this->block_type_ != ib.block_type_) return false;
            if (this->block_idx_ != ib.block_idx_)   return false;
            
            const block_descr& bd = this->bdescr_;
            const block_descr& ib_db = ib.bdescr_;

            if (this->block_type_ == 0) // bit block
            {
                if (bd.bit_.ptr != ib_db.bit_.ptr) return false;
                if (bd.bit_.idx != ib_db.bit_.idx) return false;
                if (bd.bit_.cnt != ib_db.bit_.cnt) return false;
                if (bd.bit_.pos != ib_db.bit_.pos) return false;
                for (unsigned i = 0; i < bd.bit_.cnt; ++i)
                {
                    if (bd.bit_.bits[i] != bd.bit_.bits[i]) return false;
                }
            }
            else // GAP block
            {
                if (bd.gap_.ptr != ib_db.gap_.ptr) return false;
                if (bd.gap_.gap_len != ib_db.gap_.gap_len) return false;
            }
            return true;
        }

    public:

        /** Information about current bitblock. */
        struct bitblock_descr
        {
            const bm::word_t*   ptr;      //!< Word pointer.
            unsigned char       bits[set_bitscan_wave_size*32]; //!< bit list
            unsigned short      idx;      //!< Current position in the bit list
            unsigned short      cnt;      //!< Number of ON bits
            bm::id_t            pos;      //!< Last bit position decode before
        };

        /** Information about current DGAP block. */
        struct dgap_descr
        {
            const gap_word_t*   ptr;       //!< Word pointer.
            gap_word_t          gap_len;   //!< Current dgap length.
        };

    protected:
        bm::bvector<Alloc>*     bv_;         //!< Pointer on parent bitvector
        bm::id_t                position_;   //!< Bit position (bit idx)
        const bm::word_t*       block_;      //!< Block pointer.(NULL-invalid)
        unsigned                block_type_; //!< Type of block. 0-Bit, 1-GAP
        unsigned                block_idx_;  //!< Block index

        /*! Block type dependent information for current block. */
        union block_descr
        {
            bitblock_descr   bit_;  //!< BitBlock related info.
            dgap_descr       gap_;  //!< DGAP block related info.
        } bdescr_;
    };

    /*!
        @brief Output iterator iterator designed to set "ON" bits based on
        input sequence of integers (bit indeces).

        STL container can be converted to bvector using this iterator
        Insert iterator guarantees the vector will be dynamically resized
        (set_bit does not do that).

        @note
        If you have many bits to set it is a good idea to use output iterator
        instead of explicitly calling set, because iterator may implement
        some performance specific tricks to make sure bulk insert is fast.

        @ingroup bvit
    */
    class insert_iterator
    {
    public:
#ifndef BM_NO_STL
        typedef std::output_iterator_tag  iterator_category;
#endif
        typedef unsigned value_type;
        typedef void difference_type;
        typedef void pointer;
        typedef void reference;

        insert_iterator() : bvect_(0), max_bit_(0) {}

        insert_iterator(bvector<Alloc>& bvect)
            : bvect_(&bvect), 
              max_bit_(bvect.size())
        {
            bvect_->init();
        }
        
        insert_iterator(const insert_iterator& iit)
        : bvect_(iit.bvect_),
          max_bit_(iit.max_bit_)
        {
            bvect_->init();
        }
        
        insert_iterator& operator=(const insert_iterator& ii)
        {
            bvect_ = ii.bvect_; 
            max_bit_ = ii.max_bit_;
            return *this;
        }

        insert_iterator& operator=(bm::id_t n)
        {
            BM_ASSERT(n < bm::id_max);
            BM_ASSERT_THROW(n < bm::id_max, BM_ERR_RANGE);

            if (n >= max_bit_) 
            {
                max_bit_ = n;
                if (n >= bvect_->size()) 
                {
                    bm::id_t new_size = (n == bm::id_max) ? bm::id_max : n + 1;
                    bvect_->resize(new_size);
                }
            }

            bvect_->set_bit_no_check(n);
            return *this;
        }
        
        /*! Returns *this without doing anything (no-op) */
        insert_iterator& operator*() { return *this; }
        /*! Returns *this. This iterator does not move (no-op) */
        insert_iterator& operator++() { return *this; }
        /*! Returns *this. This iterator does not move (no-op)*/
        insert_iterator& operator++(int) { return *this; }
        
    protected:
        bm::bvector<Alloc>*   bvect_;
        bm::id_t              max_bit_;
    };


    /*! @brief Constant iterator designed to enumerate "ON" bits
        @ingroup bvit
    */
    class enumerator : public iterator_base
    {
    public:
#ifndef BM_NO_STL
        typedef std::input_iterator_tag  iterator_category;
#endif
        typedef unsigned   value_type;
        typedef unsigned   difference_type;
        typedef unsigned*  pointer;
        typedef unsigned&  reference;

    public:
        enumerator() : iterator_base()
        {}
        
        /*! @brief Construct enumerator associated with a vector.
            This construction creates unpositioned iterator with status
            valid() == false. It can be re-positioned using go_first() or go_to()
        */
        enumerator(const bvector<Alloc>* bv)
            : iterator_base()
        {
            this->bv_ = const_cast<bvector<Alloc>*>(bv);
        }

        /*! @brief Construct enumerator for bit vector
            @param bv  bit-vector pointer
            @param pos bit position in the vector
                       if position is 0, it finds the next 1 or becomes not valid
                       (en.valid() == false)
        */
        enumerator(const bvector<Alloc>* bv, bm::id_t pos)
            : iterator_base()
        { 
            this->bv_ = const_cast<bvector<Alloc>*>(bv);
            this->go_to(pos);
        }

        /*! \brief Get current position (value) */
        bm::id_t operator*() const
        { 
            return this->position_; 
        }

        /*! \brief Get current position (value) */
        bm::id_t value() const
        {
            return this->position_;
        }
        /*! \brief Advance enumerator forward to the next available bit */
        enumerator& operator++()
        {
            return this->go_up();
        }

        /*! \brief Advance enumerator forward to the next available bit.
             Possibly do NOT use this operator it is slower than the pre-fix increment.
         */
        enumerator operator++(int)
        {
            enumerator tmp = *this;
            this->go_up();
            return tmp;
        }


        /*! \brief Position enumerator to the first available bit */
        void go_first()
        {
            BM_ASSERT(this->bv_);
            
            blocks_manager_type* bman = &(this->bv_->blockman_);
            if (!bman->is_init())
            {
                this->invalidate();
                return;
            }
            
            bm::word_t*** blk_root = bman->top_blocks_root();

            this->block_idx_ = this->position_= 0;
            unsigned i, j;

            for (i = 0; i < bman->top_block_size(); ++i)
            {
                bm::word_t** blk_blk = blk_root[i];

                if (blk_blk == 0) // not allocated
                {
                    this->block_idx_ += bm::set_array_size;
                    this->position_ += bm::bits_in_array;
                    continue;
                }

                for (j = 0; j < bm::set_array_size; ++j,++(this->block_idx_))
                {
                    this->block_ = blk_blk[j];

                    if (this->block_ == 0)
                    {
                        this->position_ += bits_in_block;
                        continue;
                    }

                    if (BM_IS_GAP(this->block_))
                    {
                        this->block_type_ = 1;
                        if (search_in_gapblock())
                        {
                            return;
                        }
                    }
                    else
                    {
                        if (this->block_ == FULL_BLOCK_FAKE_ADDR)
                            this->block_ = FULL_BLOCK_REAL_ADDR;

                        this->block_type_ = 0;
                        if (search_in_bitblock())
                        {
                            return;
                        }
                    }
            
                } // for j

            } // for i

            this->invalidate();
        }

        /*! \brief Advance enumerator to the next available bit */
        enumerator& go_up()
        {
            BM_ASSERT(this->valid());
            BM_ASSERT_THROW(this->valid(), BM_ERR_RANGE);

            // Current block search.
            //
            
            block_descr_type* bdescr = &(this->bdescr_);
            switch (this->block_type_)
            {
            case 0:   //  BitBlock
                {
                    // check if we can get the value from the bits traversal cache
                    unsigned short idx = ++(bdescr->bit_.idx);
                    if (idx < bdescr->bit_.cnt)
                    {
                        this->position_ = bdescr->bit_.pos + bdescr->bit_.bits[idx];
                        return *this;
                    }
                    this->position_ +=
                        (bm::set_bitscan_wave_size * 32) - bdescr->bit_.bits[--idx];
                    
                    bdescr->bit_.ptr += bm::set_bitscan_wave_size;
                    if (decode_bit_group(bdescr))
                    {
                        return *this;
                    }
                }
                break;
            case 1:   // DGAP Block
                {
                    ++this->position_;
                    if (--(bdescr->gap_.gap_len))
                    {
                        return *this;
                    }

                    // next gap is "OFF" by definition.
                    if (*(bdescr->gap_.ptr) == bm::gap_max_bits - 1)
                    {
                        break;
                    }
                    gap_word_t prev = *(bdescr->gap_.ptr);
                    unsigned int val = *(++(bdescr->gap_.ptr));
                    
                    this->position_ += val - prev;
                    // next gap is now "ON"
                    if (*(bdescr->gap_.ptr) == bm::gap_max_bits - 1)
                    {
                        break;
                    }
                    prev = *(bdescr->gap_.ptr);
                    val = *(++(bdescr->gap_.ptr));
                    bdescr->gap_.gap_len = (gap_word_t)(val - prev);
                    return *this;  // next "ON" found;
                }
            default:
                BM_ASSERT(0);

            } // switch

            if (search_in_blocks())
                return *this;
            
            this->invalidate();
            return *this;
        }
        
        /*!
            @brief Skip to specified relative rank
            @param rank - number of ON bits to go for
        */
        enumerator& skip_to_rank(bm::id_t rank)
        {
            --rank;
            if (!rank)
                return *this;
            return skip(rank);
        }
        
        /*!
            @brief Skip specified number of bits from enumeration
            @param rank - number of ON bits to skip
        */
        enumerator& skip(bm::id_t rank)
        {
            if (!this->valid() || !rank)
                return *this;
            for (; rank; --rank)
            {
                block_descr_type* bdescr = &(this->bdescr_);
                switch (this->block_type_)
                {
                case 0:   //  BitBlock
                    for (; rank; --rank)
                    {
                        unsigned short idx = ++(bdescr->bit_.idx);
                        if (idx < bdescr->bit_.cnt)
                        {
                            this->position_ = bdescr->bit_.pos + bdescr->bit_.bits[idx];
                            continue;
                        }
                        this->position_ +=
                            (bm::set_bitscan_wave_size * 32) - bdescr->bit_.bits[--idx];
                        bdescr->bit_.ptr += bm::set_bitscan_wave_size;
                        
                        if (!decode_bit_group(bdescr, rank))
                            break;
                    } // for rank
                    break;
                case 1:   // DGAP Block
                    for (; rank; --rank) // TODO: better skip logic
                    {
                        ++this->position_;
                        if (--(bdescr->gap_.gap_len))
                        {
                            continue;
                        }

                        // next gap is "OFF" by definition.
                        if (*(bdescr->gap_.ptr) == bm::gap_max_bits - 1)
                        {
                            break;
                        }
                        gap_word_t prev = *(bdescr->gap_.ptr);
                        unsigned int val = *(++(bdescr->gap_.ptr));
                        
                        this->position_ += val - prev;
                        // next gap is now "ON"
                        if (*(bdescr->gap_.ptr) == bm::gap_max_bits - 1)
                        {
                            break;
                        }
                        prev = *(bdescr->gap_.ptr);
                        val = *(++(bdescr->gap_.ptr));
                        bdescr->gap_.gap_len = (gap_word_t)(val - prev);
                    } // for rank
                    break;
                default:
                    BM_ASSERT(0);
                } // switch
                
                if (!rank)
                    return *this;

                if (!search_in_blocks())
                {
                    this->invalidate();
                    return *this;
                }
            } // for rank
            return *this;
        }
        
        /*!
            @brief go to a specific position in the bit-vector (or next)
        */
        enumerator& go_to(bm::id_t pos)
        {
            if (pos == 0)
            {
                go_first();
                return *this;
            }

            pos = this->bv_->check_or_next(pos); // find the true pos
            if (pos == 0) // no bits available
            {
                this->invalidate();
                return *this;
            }
            
            this->position_ = pos;
            unsigned nb = this->block_idx_ = unsigned(pos >>  bm::set_block_shift);
            bm::bvector<Alloc>::blocks_manager_type& bman =
                                                 this->bv_->get_blocks_manager();
            unsigned i0 = nb >> bm::set_array_shift; // top block address
            unsigned j0 = nb &  bm::set_array_mask;  // address in sub-block
            this->block_ = bman.get_block(i0, j0);

            BM_ASSERT(this->block_);
            
            this->block_type_ = (bool)BM_IS_GAP(this->block_);
            typedef typename iterator_base::block_descr block_descr_type;

            block_descr_type* bdescr = &(this->bdescr_);
            unsigned nbit = unsigned(pos & bm::set_block_mask);

            if (this->block_type_) // gap
            {
                this->position_ = nb * bm::set_block_size * 32;
                search_in_gapblock();
                
                if (this->position_ == pos)
                {
                    return *this;
                }
                this->position_ = pos;

                gap_word_t* gptr = BMGAP_PTR(this->block_);
                unsigned is_set;
                unsigned gpos = bm::gap_bfind(gptr, nbit, &is_set);
                BM_ASSERT(is_set);
                
                bdescr->gap_.ptr = gptr + gpos;
                if (gpos == 1)
                {
                    bdescr->gap_.gap_len = bm::gap_word_t(gptr[gpos] - (nbit - 1));
                }
                else
                {
                    bm::gap_word_t interval = bm::gap_word_t(gptr[gpos] - gptr[gpos - 1]);
                    bm::gap_word_t interval2 = bm::gap_word_t(nbit - gptr[gpos - 1]);
                    bdescr->gap_.gap_len = bm::gap_word_t(interval - interval2 + 1);
                }
            }
            else // bit
            {
                if (nbit == 0)
                {
                    search_in_bitblock();
                    return *this;
                }

                unsigned nword  = unsigned(nbit >> bm::set_word_shift);
                
                // check if we need to step back to match the wave
                unsigned parity = nword % bm::set_bitscan_wave_size;
                bdescr->bit_.ptr = this->block_ + (nword - parity);
                bdescr->bit_.cnt = bm::bitscan_wave(bdescr->bit_.ptr, bdescr->bit_.bits);
                BM_ASSERT(bdescr->bit_.cnt);
                bdescr->bit_.pos = (nb * bm::set_block_size * 32) + ((nword - parity) * 32);
                bdescr->bit_.idx = 0;
                nbit &= bm::set_word_mask;
                nbit += 32 * parity;
                for (unsigned i = 0; i < bdescr->bit_.cnt; ++i)
                {
                    if (bdescr->bit_.bits[i] == nbit)
                        return *this;
                    bdescr->bit_.idx++;
                } // for
                BM_ASSERT(0);
            }
            return *this;
        }


    private:
        typedef typename iterator_base::block_descr block_descr_type;
        
        
        bool decode_wave(block_descr_type* bdescr)
        {
            bdescr->bit_.cnt = bm::bitscan_wave(bdescr->bit_.ptr, bdescr->bit_.bits);
            if (bdescr->bit_.cnt) // found
            {
                bdescr->bit_.idx ^= bdescr->bit_.idx; // = 0;
                bdescr->bit_.pos = this->position_;
                this->position_ += bdescr->bit_.bits[0];
                return true;
            }
            return false;
        }
        
        bool decode_bit_group(block_descr_type* bdescr)
        {
            const word_t* block_end = this->block_ + bm::set_block_size;
            
            for (; bdescr->bit_.ptr < block_end;)
            {
                if (decode_wave(bdescr))
                    return true;
                this->position_ += bm::set_bitscan_wave_size * 32; // wave size
                bdescr->bit_.ptr += bm::set_bitscan_wave_size;
            } // for
            return false;
        }
        
        bool decode_bit_group(block_descr_type* bdescr, bm::id_t& rank)
        {
            const word_t* block_end = this->block_ + bm::set_block_size;
            
            for (; bdescr->bit_.ptr < block_end;)
            {
                const bm::id64_t* w64_p = (bm::id64_t*)bdescr->bit_.ptr;
                bm::id64_t w64 = *w64_p;
                unsigned cnt = bm::word_bitcount64(w64);
                if (rank > cnt)
                {
                    rank -= cnt;
                }
                else
                {
                    if (decode_wave(bdescr))
                        return true;
                }
                this->position_ += bm::set_bitscan_wave_size * 32; // wave size
                bdescr->bit_.ptr += bm::set_bitscan_wave_size;
            } // for
            return false;
        }

        bool search_in_bitblock()
        {
            BM_ASSERT(this->block_type_ == 0);
            
            block_descr_type* bdescr = &(this->bdescr_);
            bdescr->bit_.ptr = this->block_;
            
            return decode_bit_group(bdescr);
        }

        bool search_in_gapblock()
        {
            BM_ASSERT(this->block_type_ == 1);

            typedef typename iterator_base::block_descr block_descr_type;
            block_descr_type* bdescr = &(this->bdescr_);            

            bdescr->gap_.ptr = BMGAP_PTR(this->block_);
            unsigned bitval = *(bdescr->gap_.ptr) & 1;

            ++(bdescr->gap_.ptr);

            for (;true;)
            {
                BMREGISTER unsigned val = *(bdescr->gap_.ptr);

                if (bitval)
                {
                    gap_word_t* first = BMGAP_PTR(this->block_) + 1;
                    if (bdescr->gap_.ptr == first)
                    {
                        bdescr->gap_.gap_len = (gap_word_t)(val + 1);
                    }
                    else
                    {
                        bdescr->gap_.gap_len = 
                             (gap_word_t)(val - *(bdescr->gap_.ptr-1));
                    }
           
                    return true;
                }
                this->position_ += val + 1;
                if (val == bm::gap_max_bits - 1)
                {
                    break;
                }
                bitval ^= 1;
                ++(bdescr->gap_.ptr);
            }
            return false;
        }
        
        bool search_in_blocks()
        {
            ++(this->block_idx_);
            unsigned i = this->block_idx_ >> bm::set_array_shift;
            unsigned top_block_size = this->bv_->blockman_.top_block_size();
            for (; i < top_block_size; ++i)
            {
                bm::word_t** blk_blk = this->bv_->blockman_.top_blocks_root()[i];
                if (blk_blk == 0)
                {
                    this->block_idx_ += bm::set_array_size;
                    this->position_ += bm::bits_in_array;
                    continue;
                }

                unsigned j = this->block_idx_ & bm::set_array_mask;

                for(; j < bm::set_array_size; ++j, ++(this->block_idx_))
                {
                    this->block_ = blk_blk[j];

                    if (this->block_ == 0)
                    {
                        this->position_ += bm::bits_in_block;
                        continue;
                    }

                    this->block_type_ = BM_IS_GAP(this->block_);
                    if (this->block_type_)
                    {
                        if (search_in_gapblock())
                            return true;
                    }
                    else
                    {
                        if (this->block_ == FULL_BLOCK_FAKE_ADDR)
                            this->block_ = FULL_BLOCK_REAL_ADDR;
                        if (search_in_bitblock())
                            return true;
                    }
                } // for j
            } // for i
            return false;
        }
    };
    
    /*!
        @brief Constant iterator designed to enumerate "ON" bits
        counted_enumerator keeps bitcount, ie number of ON bits starting
        from the position 0 in the bit string up to the currently enumerated bit
        
        When increment operator called current position is increased by 1.
        
        @ingroup bvit
    */
    class counted_enumerator : public enumerator
    {
    public:
#ifndef BM_NO_STL
        typedef std::input_iterator_tag  iterator_category;
#endif
        counted_enumerator() : bit_count_(0){}
        
        counted_enumerator(const enumerator& en)
        : enumerator(en)
        {
            if (this->valid())
                bit_count_ = 1;
        }
        
        counted_enumerator& operator=(const enumerator& en)
        {
            enumerator* me = this;
            *me = en;
            if (this->valid())
                this->bit_count_ = 1;
            return *this;
        }
        
        counted_enumerator& operator++()
        {
            this->go_up();
            if (this->valid())
                ++(this->bit_count_);
            return *this;
        }

        counted_enumerator operator++(int)
        {
            counted_enumerator tmp(*this);
            this->go_up();
            if (this->valid())
                ++bit_count_;
            return tmp;
        }
        
        /*! @brief Number of bits ON starting from the .
        
            Method returns number of ON bits fromn the bit 0 to the current bit 
            For the first bit in bitvector it is 1, for the second 2 
        */
        bm::id_t count() const { return bit_count_; }
    private:
        /*! Function closed for usage */
        counted_enumerator& go_to(bm::id_t pos);

    private:
        bm::id_t   bit_count_;
    };

    /*! 
        Resource guard for bvector<>::set_allocator_pool()
        @ingroup bvector
        @internal
    */
    class mem_pool_guard
    {
    public:
        mem_pool_guard() : bv_(0)
        {}

        mem_pool_guard(allocator_pool_type& pool, bvector<Alloc>& bv)
            : bv_(&bv)
        {
            bv.set_allocator_pool(&pool);
        }
        ~mem_pool_guard()
        {
            if (bv_)
                bv_->set_allocator_pool(0);
        }

        void assign_if_not_set(allocator_pool_type& pool, bvector<Alloc>& bv)
        {
            if (bv.get_allocator_pool() == 0) // alloc pool not set yet
            {
                BM_ASSERT(!bv_);
                bv_ = &bv;
                bv.set_allocator_pool(&pool);
            }
        }

    private:
        mem_pool_guard(const mem_pool_guard&) = delete;
        void operator=(const mem_pool_guard&) = delete;
    private:
        bvector<Alloc>* bv_; ///< garded object
    };


    friend class iterator_base;
    friend class enumerator;
    template<class BV> friend class aggregator;

public:
    /*! @brief memory allocation policy
    
        Defualt memory allocation policy uses BM_BIT, and standard 
        GAP levels tune-ups
    */
    struct allocation_policy
    {
        bm::strategy  strat;
        const         gap_word_t* glevel_len;
        
        allocation_policy(bm::strategy s=BM_BIT,
                          const gap_word_t* glevels = bm::gap_len_table<true>::_len)
        : strat(s), glevel_len(glevels)
        {}
    };
    
    /*! @brief structure for bit counts per block (prefix sum)
        Structure is used to accelerate bit range scans
    */
    struct blocks_count
    {
        unsigned  BM_VECT_ALIGN cnt[bm::set_total_blocks] BM_VECT_ALIGN_ATTR;
        
        blocks_count() {}
        blocks_count(const blocks_count& bc) BMNOEXEPT
        {
            copy_from(bc);
        }
        void copy_from(const blocks_count& bc) BMNOEXEPT
        {
            ::memcpy(this->cnt, bc.cnt, sizeof(this->cnt));
        }
    };

public:
    /*! @name Construction, initialization, assignment */
    //@{
    
    /*!
        \brief Constructs bvector class
        \param strat - operation mode strategy, 
                       BM_BIT - default strategy, bvector use plain bitset 
                       blocks, (performance oriented strategy).
                       BM_GAP - memory effitent strategy, bvector allocates 
                       blocks as array of intervals(gaps) and convert blocks 
                       into plain bitsets only when enthropy grows.
        \param glevel_len 
           - pointer on C-style array keeping GAP block sizes. 
            bm::gap_len_table<true>::_len - default value set
            (use bm::gap_len_table_min<true>::_len for very sparse vectors)
            (use bm::gap_len_table_nl<true>::_len non-linear GAP growth)
        \param bv_size
          - bvector size (number of bits addressable by bvector), bm::id_max means 
          "no limits" (recommended). 
          bit vector allocates this space dynamically on demand.
        \param alloc - alllocator for this instance
         
        \sa bm::gap_len_table bm::gap_len_table_min set_new_blocks_strat
    */
    bvector(strategy          strat      = BM_BIT,
            const gap_word_t* glevel_len = bm::gap_len_table<true>::_len,
            size_type         bv_size    = bm::id_max,
            const Alloc&      alloc      = Alloc()) 
    : blockman_(glevel_len, bv_size, alloc),
      new_blocks_strat_(strat),
      size_(bv_size)
    {}

    /*!
        \brief Constructs bvector class
    */
    bvector(size_type         bv_size,
            strategy          strat      = BM_BIT,
            const gap_word_t* glevel_len = bm::gap_len_table<true>::_len,
            const Alloc&      alloc      = Alloc()) 
    : blockman_(glevel_len, bv_size, alloc),
      new_blocks_strat_(strat),
      size_(bv_size)
    {}

    /*!
        \brief Copy constructor
    */
    bvector(const bvector<Alloc>& bvect)
        :  blockman_(bvect.blockman_),
           new_blocks_strat_(bvect.new_blocks_strat_),
           size_(bvect.size_)
    {}
    
    /*!
        \brief Copy constructor for range copy [left..right]
     
        \sa copy_range
    */
    bvector(const bvector<Alloc>& bvect, bm::id_t left, bm::id_t right)
        : blockman_(bvect.blockman_.glevel_len_, bvect.blockman_.max_bits_, bvect.blockman_.alloc_),
          new_blocks_strat_(bvect.new_blocks_strat_),
          size_(bvect.size_)
    {
        if (!bvect.blockman_.is_init())
            return;
        if (left > right)
            bm::xor_swap(left, right);
        copy_range_no_check(bvect, left, right);
    }

    
    ~bvector() BMNOEXEPT
    {}
    /*!
        \brief Explicit post-construction initialization
    */
    void init();

    /*! 
        \brief Copy assignment operator
    */
    bvector& operator=(const bvector<Alloc>& bvect)
    {
        if (this != &bvect)
        {
            blockman_.deinit_tree();
            blockman_.copy(bvect.blockman_);
            resize(bvect.size());
        }
        return *this;
    }

#ifndef BM_NO_CXX11
    /*!
        \brief Move constructor
    */
    bvector(bvector<Alloc>&& bvect) BMNOEXEPT
    {
        blockman_.move_from(bvect.blockman_);
        size_ = bvect.size_;
        new_blocks_strat_ = bvect.new_blocks_strat_;
    }


    /*!
        \brief Brace constructor
    */
    bvector(std::initializer_list<bm::id_t> il) 
        : blockman_(bm::gap_len_table<true>::_len, bm::id_max, Alloc()),
          new_blocks_strat_(BM_BIT),
          size_(bm::id_max)
    {
        init();
        std::initializer_list<bm::id_t>::const_iterator it_start = il.begin();
        std::initializer_list<bm::id_t>::const_iterator it_end = il.end();
        for (; it_start < it_end; ++it_start)
        {
            this->set_bit_no_check(*it_start);
        }
    }
    
    /*! 
        \brief Move assignment operator
    */
    bvector& operator=(bvector<Alloc>&& bvect) BMNOEXEPT
    {
        this->move_from(bvect);
        return *this;
    }
#endif
    /*!
        \brief Move bvector content from another bvector
    */
    void move_from(bvector<Alloc>& bvect) BMNOEXEPT;
    
    /*! \brief Exchanges content of bv and this bvector.
    */
    void swap(bvector<Alloc>& bvect) BMNOEXEPT;
    
    //@}


    reference operator[](bm::id_t n)
    {
        if (n >= size_)
        {
            bm::id_t new_size = (n == bm::id_max) ? bm::id_max : n + 1;
            resize(new_size);
        }
        return reference(*this, n);
    }


    bool operator[](bm::id_t n) const
    {
        BM_ASSERT(n < size_);
        return get_bit(n);
    }

    void operator &= (const bvector<Alloc>& bv)
    {
        bit_and(bv);
    }

    void operator ^= (const bvector<Alloc>& bv)
    {
        bit_xor(bv);
    }

    void operator |= (const bvector<Alloc>& bv)
    {
        bit_or(bv);
    }

    void operator -= (const bvector<Alloc>& bv)
    {
        bit_sub(bv);
    }

    bool operator < (const bvector<Alloc>& bv) const
    {
        return compare(bv) < 0;
    }

    bool operator <= (const bvector<Alloc>& bv) const
    {
        return compare(bv) <= 0;
    }

    bool operator > (const bvector<Alloc>& bv) const
    {
        return compare(bv) > 0;
    }

    bool operator >= (const bvector<Alloc>& bv) const
    {
        return compare(bv) >= 0;
    }

    bool operator == (const bvector<Alloc>& bv) const
    {
        return compare(bv) == 0;
    }

    bool operator != (const bvector<Alloc>& bv) const
    {
        return compare(bv) != 0;
    }

    bvector<Alloc> operator~() const
    {
        return bvector<Alloc>(*this).invert();
    }
    
    Alloc get_allocator() const
    {
        return blockman_.get_allocator();
    }

    /// Set allocator pool for local (non-threaded) 
    /// memory cyclic(lots of alloc-free ops) opertations
    ///
    void set_allocator_pool(allocator_pool_type* pool_ptr)
    {
        blockman_.get_allocator().set_pool(pool_ptr);
    }

    /// Get curent allocator pool (if set)
    /// @return pointer to the current pool or NULL
    allocator_pool_type* get_allocator_pool()
    {
        return blockman_.get_allocator().get_pool();
    }

    // --------------------------------------------------------------------
    /*! @name Bit access/modification methods  */
    //@{

    /*!
       \brief Sets bit n.
       \param n - index of the bit to be set. 
       \param val - new bit value
       \return  TRUE if bit was changed
    */
    bool set_bit(bm::id_t n, bool val = true)
    {
        BM_ASSERT_THROW(n < bm::id_max, BM_ERR_RANGE);

        if (!blockman_.is_init())
            blockman_.init_tree();
        if (n >= size_)
        {
            bm::id_t new_size = (n == bm::id_max) ? bm::id_max : n + 1;
            resize(new_size);
        }
        
        return set_bit_no_check(n, val);
    }

    /*!
       \brief Sets bit n using bit AND with the provided value.
       \param n - index of the bit to be set. 
       \param val - new bit value
       \return  TRUE if bit was changed
    */
    bool set_bit_and(bm::id_t n, bool val = true)
    {
        BM_ASSERT(n < size_);
        BM_ASSERT_THROW(n < size_, BM_ERR_RANGE);
        
        if (!blockman_.is_init())
        {
            blockman_.init_tree();
        }
        return and_bit_no_check(n, val);
    }
    
    /*!
       \brief Increment the specified element
     
       Bit increment rules:
        0 + 1 = 1 (no carry over)
        1 + 1 = 0 (with carry over returned)
     
       \param n - index of the bit to be set
       \return  TRUE if carry over created (1+1)
    */
    bool inc(bm::id_t n);
    

    /*!
       \brief Sets bit n only if current value equals the condition
       \param n - index of the bit to be set. 
       \param val - new bit value
       \param condition - expected current value
       \return TRUE if bit was changed
    */
    bool set_bit_conditional(bm::id_t n, bool val, bool condition)
    {
        if (val == condition) return false;
        if (!blockman_.is_init())
            blockman_.init_tree();
        if (n >= size_)
        {
            bm::id_t new_size = (n == bm::id_max) ? bm::id_max : n + 1;
            resize(new_size);
        }

        return set_bit_conditional_impl(n, val, condition);
    }


    /*!
        \brief Sets bit n if val is true, clears bit n if val is false
        \param n - index of the bit to be set
        \param val - new bit value
        \return *this
    */
    bvector<Alloc>& set(bm::id_t n, bool val = true)
    {
        set_bit(n, val);
        return *this;
    }

    /*!
       \brief Sets every bit in this bitset to 1.
       \return *this
    */
    bvector<Alloc>& set()
    {
        set_range(0, size_ - 1, true);
        return *this;
    }
    
    /*!
        \brief Set bit without checking preconditions (size, etc)
     
        Fast set bit method, without safety net.
        Make sure you call bvector<>::init() before setting bits with this
        function.
     
        \param n - bit number
    */
    void set_bit_no_check(bm::id_t n);


    /*!
        \brief Sets all bits in the specified closed interval [left,right]
        Interval must be inside the bvector's size. 
        This method DOES NOT resize vector.
        
        \param left  - interval start
        \param right - interval end (closed interval)
        \param value - value to set interval in
        
        \return *this
    */
    bvector<Alloc>& set_range(bm::id_t left,
                              bm::id_t right,
                              bool     value = true);
    
    /*!
        \brief Copy all bits in the specified closed interval [left,right]

        \param bvect - source bit-vector
        \param left  - interval start
        \param right - interval end (closed interval)
    */
    void copy_range(const bvector<Alloc>& bvect,
                    bm::id_t left,
                    bm::id_t right);

    /*!
       \brief Clears bit n.
       \param n - bit's index to be cleaned.
       \return true if bit was cleared
    */
    bool clear_bit(bm::id_t n) { return set_bit(n, false); }
    
    /*!
       \brief Clears bit n without precondiion checks
       \param n - bit's index to be cleaned.
    */
    void clear_bit_no_check(bm::id_t n) { set_bit_no_check(n, false); }

    
    
    /*!
       \brief Clears every bit in the bitvector.

       \param free_mem if "true" (default) bvector frees the memory,
       otherwise sets blocks to 0.
    */
    void clear(bool free_mem = false)
    {
        blockman_.set_all_zero(free_mem);
    }

    /*!
       \brief Clears every bit in the bitvector.
       \return *this;
    */
    bvector<Alloc>& reset()
    {
        clear(true);
        return *this;
    }
    
    /*!
       \brief Flips bit n
       \return *this
    */
    bvector<Alloc>& flip(bm::id_t n)
    {
        //set(n, !get_bit(n));
        this->inc(n);
        return *this;
    }

    /*!
       \brief Flips all bits
       \return *this
       @sa invert
    */
    bvector<Alloc>& flip()
    {
        return invert();
    }

    //@}
    // --------------------------------------------------------------------

    
    /*! Function erturns insert iterator for this bitvector */
    insert_iterator inserter()
    {
        return insert_iterator(*this);
    }

    // --------------------------------------------------------------------
    /*! @name Size and capacity
        By default bvector is dynamically sized, manual control methods
        available
    */
    //@{

    /** \brief Returns bvector's capacity (number of bits it can store) */
    size_type capacity() const 
    {
        return blockman_.capacity();
    }

    /*! \brief return current size of the vector (bits) */
    size_type size() const 
    {
        return size_;
    }

    /*!
        \brief Change size of the bvector
        \param new_size - new size in bits
    */
    void resize(size_type new_size);
    
    //@}
    // --------------------------------------------------------------------

    /*! @name Population counting and ranking methods
    */
    //@{

    /*!
       \brief population cout (count of ON bits)
       \return Total number of bits ON.
    */
    bm::id_t count() const;

    /*! \brief Computes bitcount values for all bvector blocks
        \param arr - pointer on array of block bit counts
        \return Index of the last block counted. 
        This number +1 gives you number of arr elements initialized during the
        function call.
    */
    unsigned count_blocks(unsigned* arr) const
    {
        bm::word_t*** blk_root = blockman_.top_blocks_root();
        if (blk_root == 0)
            return 0;
        typename blocks_manager_type::block_count_arr_func func(blockman_, &(arr[0]));
        for_each_nzblock(blk_root, blockman_.top_block_size(), func);
        return func.last_block();
    }

    /*!
       \brief Returns count of 1 bits in the given range
     
       \param left - index of first bit start counting from
       \param right - index of last bit 
       \param block_count_arr - optional parameter (bitcount by bvector blocks)
              calculated by count_blocks method. Used to improve performance of
              wide range searches
       \return population count in the diapason
    */
    bm::id_t count_range(bm::id_t left, 
                         bm::id_t right, 
                         const unsigned* block_count_arr=0) const;
    
    /*! \brief compute running total of all blocks in bit vector
        \param blocks_cnt - out pointer to counting structure, holding the array
        Function will fill full array of running totals
        \sa count_to
    */
    void running_count_blocks(blocks_count* blocks_cnt) const;
    
    /*!
       \brief Returns count of 1 bits (population) in [0..right] range.
     
       This operation is also known as rank of bit N.
     
       \param n - index of bit to rank
       \param blocks_cnt - block count structure to accelerate search
                           should be prepared using running_count_blocks
       \return population count in the range
       \sa running_count_blocks
       \sa count_to_test
    */
    bm::id_t count_to(bm::id_t n, const blocks_count&  blocks_cnt) const;

    /*!
        \brief Returns count of 1 bits (population) in [0..right] range if test(right) == true
     
        This is conditional rank operation, which is faster than test()
        plus count_to()
     
        \param n - index of bit to test and rank
        \param blocks_cnt - block count structure to accelerate search
               should be prepared using running_count_blocks

        \return population count in the diapason or 0 if right bit test failed

        \sa running_count_blocks
        \sa count_to
    */
    bm::id_t count_to_test(bm::id_t n, const blocks_count&  blocks_cnt) const;


    /*! Recalculate bitcount (deprecated)
    */
    bm::id_t recalc_count()
    {
        return count();
    }
    
    /*!
        Disables count cache. (deprecated).
    */
    void forget_count()
    {
    }
    //@}
    
    // --------------------------------------------------------------------
    /*! @name Bit access (read-only)  */
    //@{

    /*!
       \brief returns true if bit n is set and false is bit n is 0. 
       \param n - Index of the bit to check.
       \return Bit value (1 or 0)
    */
    bool get_bit(bm::id_t n) const;

    /*!
       \brief returns true if bit n is set and false is bit n is 0. 
       \param n - Index of the bit to check.
       \return Bit value (1 or 0)
    */
    bool test(bm::id_t n) const 
    { 
        return get_bit(n); 
    }
    //@}

    // --------------------------------------------------------------------
    /*! @name Check for empty-ness of container  */
    //@{

    /*!
       \brief Returns true if any bits in this bitset are set, and otherwise returns false.
       \return true if any bit is set
    */
    bool any() const
    {
        word_t*** blk_root = blockman_.top_blocks_root();
        if (!blk_root) 
            return false;
        typename blocks_manager_type::block_any_func func(blockman_);
        return for_each_nzblock_if(blk_root, blockman_.top_block_size(), func);
    }

    /*!
        \brief Returns true if no bits are set, otherwise returns false.
    */
    bool none() const
    {
        return !any();
    }
    //@}
    // --------------------------------------------------------------------


    /*! @name Scan and find bits and indexes */
    //@{
    
    /*!
       \fn bool bvector::find(bm::id_t& pos) const
       \brief Finds index of first 1 bit
       \param pos - index of the found 1 bit
       \return true if search returned result
       \sa get_first, get_next, extract_next, find_reverse
    */
    bool find(bm::id_t& pos) const;

    /*!
       \fn bool bvector::find(bm::id_t from, bm::id_t& pos) const
       \brief Finds index of 1 bit starting from position
       \param from - position to start search from
       \param pos - index of the found 1 bit
       \return true if search returned result
       \sa get_first, get_next, extract_next, find_reverse
    */
    bool find(bm::id_t from, bm::id_t& pos) const;

    /*!
       \fn bm::id_t bvector::get_first() const
       \brief find first 1 bit in vector. 
       Function may return 0 and this requires an extra check if bit 0 is 
       actually set or bit-vector is empty
     
       \return Index of the first 1 bit, may return 0
       \sa get_next, find, extract_next, find_reverse
    */
    bm::id_t get_first() const { return check_or_next(0); }

    /*!
       \fn bm::id_t bvector::get_next(bm::id_t prev) const
       \brief Finds the number of the next bit ON.
       \param prev - Index of the previously found bit. 
       \return Index of the next bit which is ON or 0 if not found.
       \sa get_first, find, extract_next, find_reverse
    */
    bm::id_t get_next(bm::id_t prev) const
    {
        return (++prev == bm::id_max) ? 0 : check_or_next(prev);
    }

    /*!
       \fn bm::id_t bvector::extract_next(bm::id_t prev)
       \brief Finds the number of the next bit ON and sets it to 0.
       \param prev - Index of the previously found bit. 
       \return Index of the next bit which is ON or 0 if not found.
       \sa get_first, get_next, find_reverse
    */
    bm::id_t extract_next(bm::id_t prev)
    {
        return (++prev == bm::id_max) ? 0 : check_or_next_extract(prev);
    }

    /*!
       \brief Finds last index of 1 bit
       \param pos - index of the last found 1 bit
       \return true if search returned result
       \sa get_first, get_next, extract_next, find
    */
    bool find_reverse(bm::id_t& pos) const;
    
    /*!
       \brief Finds dynamic range of bit-vector [first, last]
       \param first - index of the first found 1 bit
       \param last - index of the last found 1 bit
       \return true if search returned result
       \sa get_first, get_next, extract_next, find, find_reverse
    */
    bool find_range(bm::id_t& first, bm::id_t& last) const;
    
    /*!
        \brief Find bit-vector position for the specified rank(bitcount)
     
        Rank based search, counts number of 1s from specified position until
        finds the ranked position relative to start from position.
        In other words: range population count between from and pos == rank.
     
        \param rank - rank to find (bitcount)
        \param from - start positioon for rank search
        \param pos  - position with speciefied rank (relative to from position)
     
        \return true if requested rank was found
    */
    bool find_rank(bm::id_t rank, bm::id_t from, bm::id_t& pos) const;

    
    //@}


    // --------------------------------------------------------------------
    /*! @name Set Algebra operations  */
    //@{

    /*!
       \brief Logical OR operation.
       \param vect - Argument vector.
    */
    bm::bvector<Alloc>& bit_or(const  bm::bvector<Alloc>& vect)
    {
        combine_operation_or(vect);
        return *this;
    }

    /*!
       \brief Logical AND operation.
       \param bv - argument vector.
    */
    bm::bvector<Alloc>& bit_and(const bm::bvector<Alloc>& bv)
    {
        combine_operation_and(bv);
        return *this;
    }

    /*!
       \brief Logical XOR operation.
       \param vect - argument vector.
    */
    bm::bvector<Alloc>& bit_xor(const bm::bvector<Alloc>& vect)
    {
        combine_operation(vect, BM_XOR);
        return *this;
    }

    /*!
       \brief Logical SUB operation.
       \param vect - Argument vector.
    */
    bm::bvector<Alloc>& bit_sub(const bm::bvector<Alloc>& bv)
    {
        combine_operation_sub(bv);
        return *this;
    }
    
    /*!
        \brief Inverts all bits.
    */
    bvector<Alloc>& invert();
    
    /*! \brief perform a set-algebra operation by operation code
    */
    void combine_operation(const bm::bvector<Alloc>& bvect,
                            bm::operation            opcode);

    /*! \brief perform a set-algebra operation OR
    */
    void combine_operation_or(const bm::bvector<Alloc>& bvect);

    /*! \brief perform a set-algebra operation AND
    */
    void combine_operation_and(const bm::bvector<Alloc>& bvect);

    /*! \brief perform a set-algebra operation MINUS (AND NOT)
    */
    void combine_operation_sub(const bm::bvector<Alloc>& bvect);

    // @}

    // --------------------------------------------------------------------
    /*! @name Iterator-traversal methods  */
    //@{

    /**
       \brief Returns enumerator pointing on the first non-zero bit.
    */
    enumerator first() const { return get_enumerator(0); }

    /**
       \fn bvector::enumerator bvector::end() const
       \brief Returns enumerator pointing on the next bit after the last.
    */
    enumerator end() const
    {
        typedef typename bvector<Alloc>::enumerator enumerator_type;
        return enumerator_type(this);
    }
    
    /**
       \brief Returns enumerator pointing on specified or the next available bit.
    */
    enumerator get_enumerator(bm::id_t pos) const
    {
        typedef typename bvector<Alloc>::enumerator enumerator_type;
        return enumerator_type(this, pos);
    }
    
    //@}

    // --------------------------------------------------------------------
    /*! @name Memory management and compression  */
    
    //@{
    
    /*!
       @brief Calculates bitvector statistics.

       @param st - pointer on statistics structure to be filled in.

       Function fills statistics structure containing information about how
       this vector uses memory and estimation of max. amount of memory
       bvector needs to serialize itself.

       @sa statistics
    */
    void calc_stat(struct bm::bvector<Alloc>::statistics* st) const;

    /*!
       \brief Sets new blocks allocation strategy.
       \param strat - Strategy code 0 - bitblocks allocation only.
                      1 - Blocks mutation mode (adaptive algorithm)
    */
    void set_new_blocks_strat(strategy strat) 
    { 
        new_blocks_strat_ = strat; 
    }

    /*!
       \brief Returns blocks allocation strategy.
       \return - Strategy code 0 - bitblocks allocation only.
                 1 - Blocks mutation mode (adaptive algorithm)
       \sa set_new_blocks_strat
    */
    strategy  get_new_blocks_strat() const 
    { 
        return new_blocks_strat_; 
    }

    /*! 
        \brief Optimization mode
        Every next level means additional checks (better compression vs time)
        \sa optimize
    */
    enum optmode
    {
        opt_free_0    = 1, ///< Free unused 0 blocks
        opt_free_01   = 2, ///< Free unused 0 and 1 blocks
        opt_compress  = 3  ///< compress blocks when possible
    };

    /*!
       \brief Optimize memory bitvector's memory allocation.
   
       Function analyze all blocks in the bitvector, compresses blocks 
       with a regular structure, frees some memory. This function is recommended
       after a bulk modification of the bitvector using set_bit, clear_bit or
       logical operations.

       Optionally function can calculate vector post optimization statistics
       
       @sa optmode, optimize_gap_size
    */
    void optimize(bm::word_t* temp_block = 0,
                  optmode opt_mode       = opt_compress,
                  statistics* stat       = 0);

    /*!
       \brief Optimize sizes of GAP blocks

       This method runs an analysis to find optimal GAP levels for the 
       specific vector. Current GAP compression algorithm uses several fixed
       GAP sizes. By default bvector uses some reasonable preset. 
    */
    void optimize_gap_size();


    /*!
        @brief Sets new GAP lengths table. All GAP blocks will be reallocated 
        to match the new scheme.

        @param glevel_len - pointer on C-style array keeping GAP block sizes. 
    */
    void set_gap_levels(const gap_word_t* glevel_len);
    
    //@}
    
    // --------------------------------------------------------------------
    
    /*! @name Comparison  */
    //@{

    /*!
        \brief Lexicographical comparison with a bitvector.

        Function compares current bitvector with the provided argument 
        bit by bit and returns -1 if our bitvector less than the argument, 
        1 - greater, 0 - equal.
    */
    int compare(const bvector<Alloc>& bvect) const;
    
    //@}

    // --------------------------------------------------------------------
    /*! @name Open internals   */
    //@{

    /*! @brief Allocates temporary block of memory. 

        Temp block can be passed to bvector functions requiring some temp memory
        for their operation. (like serialize)
        
        @note method is marked const, but it's not quite true, since
        it can in some cases modify the state of the block allocator
        (if it has a state). (Can be important in MT programs).

        @sa free_tempblock
    */
    bm::word_t* allocate_tempblock() const
    {
        blocks_manager_type* bm = 
            const_cast<blocks_manager_type*>(&blockman_);
        return bm->get_allocator().alloc_bit_block();
    }

    /*! @brief Frees temporary block of memory. 

        @note method is marked const, but it's not quite true, since
        it can in some cases modify the state of the block allocator
        (if it has a state). (Can be important in MT programs).

        @sa allocate_tempblock
    */
    void free_tempblock(bm::word_t* block) const
    {
        blocks_manager_type* bm = 
            const_cast<blocks_manager_type*>(&blockman_);
        bm->get_allocator().free_bit_block(block);
    }


    /*!
        \brief get access to internal block by number
     
        This is an internal method, declared public to add low level
        optimized algorithms outside of the main bvector<> class.
        Use it AT YOUR OWN RISK!
     
        \param nb - block number
     
        \internal
    */
    const bm::word_t* get_block(unsigned nb) const
    {
        return blockman_.get_block(nb);
    }
    
    /*!
    @internal
    */
    void combine_operation_with_block(unsigned nb,
                                      const bm::word_t* arg_blk,
                                      bool arg_gap,
                                      bm::operation opcode)
    {
        bm::word_t* blk = const_cast<bm::word_t*>(get_block(nb));
        bool gap = BM_IS_GAP(blk);
        combine_operation_with_block(nb, gap, blk, arg_blk, arg_gap, opcode);
    }
    
    
    const blocks_manager_type& get_blocks_manager() const
    {
        return blockman_;
    }
    blocks_manager_type& get_blocks_manager()
    {
        return blockman_;
    }

    //@}
    
    
private:

    bm::id_t check_or_next(bm::id_t prev) const;
    
    /// set bit in GAP block withlength extension control
    bool gap_block_set(bm::gap_word_t* gap_blk,
                       bool val, unsigned nblock, unsigned nbit);
    
    /// check if specified bit is 1, and set it to 0
    /// if specified bit is 0, scan for the next 1 and returns it
    /// if no 1 found returns 0
    bm::id_t check_or_next_extract(bm::id_t prev);

    /**
        \brief Set specified bit without checking preconditions (size, etc)
    */
    bool set_bit_no_check(bm::id_t n, bool val);

    /**
        \brief AND specified bit without checking preconditions (size, etc)
    */
    bool and_bit_no_check(bm::id_t n, bool val);

    bool set_bit_conditional_impl(bm::id_t n, bool val, bool condition);


    void combine_operation_with_block(unsigned nb,
                                      bool gap,
                                      bm::word_t* blk,
                                      const bm::word_t* arg_blk,
                                      bool arg_gap,
                                      bm::operation opcode);
    
    void combine_operation_block_or(unsigned i,
                                    unsigned j,
                                    bm::word_t* blk,
                                    const bm::word_t* arg_blk);

    void combine_operation_block_and(unsigned i,
                                     unsigned j,
                                     bm::word_t* blk,
                                     const bm::word_t* arg_blk);
    void combine_operation_block_sub(unsigned i,
                                     unsigned j,
                                     bm::word_t* blk,
                                     const bm::word_t* arg_blk);

    void assign_gap_result(unsigned              nb,
                           const bm::gap_word_t* res,
                           unsigned              res_len,
                           bm::word_t*           blk,
                           gap_word_t*           tmp_buf);

    void assign_gap_result(unsigned              i,
                           unsigned              j,
                           const bm::gap_word_t* res,
                           unsigned              res_len,
                           bm::word_t*           blk,
                           gap_word_t*           tmp_buf);
    
    void copy_range_no_check(const bvector<Alloc>& bvect,
                             bm::id_t left,
                             bm::id_t right);

private:
    /**
       \brief Extends GAP block to the next level or converts it to bit block.
       \param nb - Block's linear index.
       \param blk - Blocks's pointer 
    */
    void extend_gap_block(unsigned nb, gap_word_t* blk)
    {
        blockman_.extend_gap_block(nb, blk);
    }

    /**
       \brief Set range without validity/bouds checking
    */
    void set_range_no_check(bm::id_t left,
                            bm::id_t right);
    /**
        \brief Clear range without validity/bouds checking
    */
    void clear_range_no_check(bm::id_t left,
                              bm::id_t right);

private:
    blocks_manager_type  blockman_;         //!< bitblocks manager
    strategy             new_blocks_strat_; //!< block allocation strategy
    size_type            size_;             //!< size in bits
};





//---------------------------------------------------------------------

template<class Alloc> 
inline bvector<Alloc> operator& (const bvector<Alloc>& v1,
                                 const bvector<Alloc>& v2)
{
#ifdef BM_USE_EXPLICIT_TEMP
    bvector<Alloc> ret(v1);
    ret.bit_and(v2);
    return ret;
#else    
    return bvector<Alloc>(v1).bit_and(v2);
#endif
}

//---------------------------------------------------------------------

template<class Alloc> 
inline bvector<Alloc> operator| (const bvector<Alloc>& v1,
                                 const bvector<Alloc>& v2)
{
#ifdef BM_USE_EXPLICIT_TEMP
    bvector<Alloc> ret(v1);
    ret.bit_or(v2);
    return ret;
#else    
    return bvector<Alloc>(v1).bit_or(v2);
#endif
}

//---------------------------------------------------------------------

template<class Alloc> 
inline bvector<Alloc> operator^ (const bvector<Alloc>& v1,
                                 const bvector<Alloc>& v2)
{
#ifdef BM_USE_EXPLICIT_TEMP
    bvector<Alloc> ret(v1);
    ret.bit_xor(v2);
    return ret;
#else    
    return bvector<Alloc>(v1).bit_xor(v2);
#endif
}

//---------------------------------------------------------------------

template<class Alloc> 
inline bvector<Alloc> operator- (const bvector<Alloc>& v1,
                                 const bvector<Alloc>& v2)
{
#ifdef BM_USE_EXPLICIT_TEMP
    bvector<Alloc> ret(v1);
    ret.bit_sub(v2);
    return ret;
#else    
    return bvector<Alloc>(v1).bit_sub(v2);
#endif
}


// -----------------------------------------------------------------------

template<typename Alloc>
void bvector<Alloc>::init()
{
    if (!blockman_.is_init())
        blockman_.init_tree();
}

// -----------------------------------------------------------------------

template<typename Alloc>
void bvector<Alloc>::move_from(bvector<Alloc>& bvect) BMNOEXEPT
{
    if (this != &bvect)
    {
        blockman_.move_from(bvect.blockman_);
        size_ = bvect.size_;
        new_blocks_strat_ = bvect.new_blocks_strat_;
    }
}



// -----------------------------------------------------------------------

template<typename Alloc> 
bvector<Alloc>& bvector<Alloc>::set_range(bm::id_t left,
                                          bm::id_t right,
                                          bool     value)
{
    if (!blockman_.is_init())
    {
        if (!value)
            return *this; // nothing to do
        blockman_.init_tree();
    }

    if (right < left)
    {
        return set_range(right, left, value);
    }

    BM_ASSERT_THROW(right < bm::id_max, BM_ERR_RANGE);
    if (right >= size_) // this vect shorter than the arg.
    {
        bm::id_t new_size = (right == bm::id_max) ? bm::id_max : right + 1;
        resize(new_size);
    }

    BM_ASSERT(left <= right);
    if (left >= size_)
    {
        std::cout << "size:" << size_ << " left=" << left << std::endl;
    }
    BM_ASSERT(left < size_);

    if (value)
        set_range_no_check(left, right);
    else
        clear_range_no_check(left, right);

    return *this;
}

// -----------------------------------------------------------------------

template<typename Alloc> 
bm::id_t bvector<Alloc>::count() const
{
    if (!blockman_.is_init())
        return 0;
    
    word_t*** blk_root = blockman_.top_blocks_root();
    if (!blk_root) 
    {
        return 0;
    }    
    typename blocks_manager_type::block_count_func func(blockman_);
    for_each_nzblock2(blk_root, blockman_.top_block_size(), func);

    return func.count();
}

// -----------------------------------------------------------------------

template<typename Alloc> 
void bvector<Alloc>::resize(size_type new_size)
{
    if (size_ == new_size) return; // nothing to do
    if (size_ < new_size) // size grows 
    {
        if (!blockman_.is_init())
            blockman_.init_tree();
    
        blockman_.reserve(new_size);
        size_ = new_size;
    }
    else // shrink
    {
        set_range(new_size, size_ - 1, false); // clear the tail
        size_ = new_size;
    }
}

// -----------------------------------------------------------------------

template<typename Alloc>
void bvector<Alloc>::running_count_blocks(blocks_count* blocks_cnt) const
{
    BM_ASSERT(blocks_cnt);
    
    ::memset(blocks_cnt->cnt, 0, sizeof(blocks_cnt->cnt));
    
    this->count_blocks(blocks_cnt->cnt);  // simple count
    // compute running count
    for (unsigned i = 1; i < bm::set_total_blocks; ++i)
    {
        blocks_cnt->cnt[i] += blocks_cnt->cnt[i-1];
    }
}
// -----------------------------------------------------------------------

template<typename Alloc>
bm::id_t bvector<Alloc>::count_to(bm::id_t right,
                                  const blocks_count&  blocks_cnt) const
{
    if (!blockman_.is_init())
        return 0;

    unsigned nblock_right = unsigned(right >>  bm::set_block_shift);    
    unsigned nbit_right = unsigned(right & bm::set_block_mask);

    // running count of all blocks before target
    //
    bm::id_t cnt = nblock_right ? blocks_cnt.cnt[nblock_right-1] : 0;

    const bm::word_t* block = blockman_.get_block_ptr(nblock_right);
    if (!block)
        return cnt;

    bool gap = BM_IS_GAP(block);
    if (gap)
    {
        unsigned c = bm::gap_bit_count_to(BMGAP_PTR(block), (gap_word_t)nbit_right);
        BM_ASSERT(c == bm::gap_bit_count_range(BMGAP_PTR(block), (gap_word_t)0, (gap_word_t)nbit_right));
        cnt += c;
    }
    else
    {
        cnt += 
            (block == FULL_BLOCK_FAKE_ADDR) ? (nbit_right + 1) : bm::bit_block_calc_count_to(block, nbit_right);
    }

    return cnt;
}

// -----------------------------------------------------------------------

template<typename Alloc>
bm::id_t bvector<Alloc>::count_to_test(bm::id_t right,
                                       const blocks_count&  blocks_cnt) const
{
    if (!blockman_.is_init())
        return 0;

    unsigned nblock_right = unsigned(right >> bm::set_block_shift);
    unsigned nbit_right = unsigned(right & bm::set_block_mask);

    // running count of all blocks before target
    //
    bm::id_t cnt = 0;
    
    const bm::word_t* block = blockman_.get_block_ptr(nblock_right);
    if (!block)
        return 0;

    bool gap = BM_IS_GAP(block);
    if (gap)
    {
        bm::gap_word_t *gap_blk = BMGAP_PTR(block);
        if (bm::gap_test_unr(gap_blk, (gap_word_t)nbit_right))
            cnt = bm::gap_bit_count_to(gap_blk, (gap_word_t)nbit_right);
        else
            return 0;
    }
    else
    {
        if (block == FULL_BLOCK_FAKE_ADDR)
        {
            cnt = nbit_right + 1;
        }
        else
        {
            unsigned nword = nbit_right >> bm::set_word_shift;
            unsigned w = block[nword];
            w &= (1u << (nbit_right & bm::set_word_mask));
            if (w)
                cnt = bm::bit_block_calc_count_to(block, nbit_right);
            else
                return 0;
        }
    }
    cnt += nblock_right ? blocks_cnt.cnt[nblock_right - 1] : 0;
    return cnt;
}


// -----------------------------------------------------------------------

template<typename Alloc> 
bm::id_t bvector<Alloc>::count_range(bm::id_t left, 
                                     bm::id_t right,
                                     const unsigned* block_count_arr) const
{
    BM_ASSERT(left <= right);

    BM_ASSERT_THROW(right < bm::id_max, BM_ERR_RANGE);
    BM_ASSERT_THROW(left <= right, BM_ERR_RANGE);

    
    if (!blockman_.is_init())
        return 0;

    unsigned cnt = 0;

    // calculate logical number of start and destination blocks
    unsigned nblock_left  = unsigned(left  >>  bm::set_block_shift);
    unsigned nblock_right = unsigned(right >>  bm::set_block_shift);

    const bm::word_t* block = blockman_.get_block(nblock_left);
    bool left_gap = BM_IS_GAP(block);

    unsigned nbit_left  = unsigned(left  & bm::set_block_mask); 
    unsigned nbit_right = unsigned(right & bm::set_block_mask); 

    unsigned r = 
        (nblock_left == nblock_right) ? nbit_right : (bm::bits_in_block-1);

    typename blocks_manager_type::block_count_func func(blockman_);

    if (block)
    {
        if ((nbit_left == 0) && (r == (bm::bits_in_block-1))) // whole block
        {
            if (block_count_arr)
            {
                cnt += block_count_arr[nblock_left];
            }
            else
            {
                func(block);//, nblock_left);
            }
        }
        else
        {
            if (left_gap)
            {
                cnt += gap_bit_count_range(BMGAP_PTR(block),
                                            (gap_word_t)nbit_left,
                                            (gap_word_t)r);
            }
            else
            {
                cnt += bit_block_calc_count_range(block, nbit_left, r);
            }
        }
    }

    if (nblock_left == nblock_right)  // in one block
    {
        return cnt + func.count();
    }

    for (unsigned nb = nblock_left+1; nb < nblock_right; ++nb)
    {
        block = blockman_.get_block(nb);
        if (block_count_arr)
        {
            cnt += block_count_arr[nb];
        }
        else 
        {
            if (block)
                func(block);//, nb);
        }
    }
    cnt += func.count();

    block = blockman_.get_block(nblock_right);
    bool right_gap = BM_IS_GAP(block);

    if (block)
    {
        if (right_gap)
        {
            cnt += gap_bit_count_range(BMGAP_PTR(block),
                                        (gap_word_t)0,
                                        (gap_word_t)nbit_right);
        }
        else
        {
            cnt += bit_block_calc_count_range(block, 0, nbit_right);
        }
    }

    return cnt;
}

// -----------------------------------------------------------------------

template<typename Alloc>
bvector<Alloc>& bvector<Alloc>::invert()
{
    blockman_.reserve_top_blocks(bm::set_array_size);

    bm::word_t*** blk_root = blockman_.top_blocks_root();
    typename blocks_manager_type::block_invert_func func(blockman_);    
    for_each_block(blk_root, blockman_.top_block_size(), func);
    if (size_ == bm::id_max) 
    {
        set_bit_no_check(bm::id_max, false);
    } 
    else
    {
        clear_range_no_check(size_, bm::id_max);
    }

    return *this;
}

// -----------------------------------------------------------------------

template<typename Alloc> 
bool bvector<Alloc>::get_bit(bm::id_t n) const
{    
    BM_ASSERT(n < size_);
    BM_ASSERT_THROW((n < size_), BM_ERR_RANGE);

    // calculate logical block number
    unsigned nblock = unsigned(n >>  bm::set_block_shift); 
    const bm::word_t* block = blockman_.get_block_ptr(nblock); // get unsanitized block ptr

    if (block)
    {
        // calculate word number in block and bit
        unsigned nbit = unsigned(n & bm::set_block_mask); 
        if (BM_IS_GAP(block))
        {
            return gap_test_unr(BMGAP_PTR(block), nbit);
        }
        else 
        {
            block = BLOCK_ADDR_SAN(block); // FULL BLOCK ADDR check
            unsigned nword  = unsigned(nbit >> bm::set_word_shift);
            unsigned w = block[nword];
            return w & (1u << (nbit & bm::set_word_mask));
        }
    }
    return false;
}

// -----------------------------------------------------------------------

template<typename Alloc> 
void bvector<Alloc>::optimize(bm::word_t* temp_block,
                              optmode     opt_mode,
                              statistics* stat)
{
    if (!blockman_.is_init())
    {
        if (stat)
            calc_stat(stat);
        return;
    }
    word_t*** blk_root = blockman_.top_blocks_root();

    if (!temp_block)
        temp_block = blockman_.check_allocate_tempblock();

    typename 
        blocks_manager_type::block_opt_func  opt_func(blockman_, 
                                                temp_block, 
                                                (int)opt_mode,
                                                stat);
    if (stat)
    {
        stat->reset();
        ::memcpy(stat->gap_levels,
                blockman_.glen(), sizeof(gap_word_t) * bm::gap_levels);
        stat->max_serialize_mem = (unsigned)sizeof(id_t) * 4;
    }

    for_each_nzblock(blk_root, blockman_.top_block_size(), opt_func);

    if (stat)
    {
        size_t safe_inc = stat->max_serialize_mem / 10; // 10% increment
        if (!safe_inc) safe_inc = 256;
        stat->max_serialize_mem += safe_inc;
        stat->memory_used += (unsigned)(sizeof(*this) - sizeof(blockman_));
        stat->memory_used += blockman_.mem_used();
    }
    
    // the expectation is that we don't need to keep temp block if we
    // optimizing memory usage
    blockman_.free_temp_block();
}

// -----------------------------------------------------------------------

template<typename Alloc> 
void bvector<Alloc>::optimize_gap_size()
{
    if (!blockman_.is_init())
        return;

    struct bvector<Alloc>::statistics st;
    calc_stat(&st);

    if (!st.gap_blocks)
        return;

    gap_word_t opt_glen[bm::gap_levels];
    ::memcpy(opt_glen, st.gap_levels, bm::gap_levels * sizeof(*opt_glen));

    improve_gap_levels(st.gap_length, 
                            st.gap_length + st.gap_blocks, 
                            opt_glen);
    
    set_gap_levels(opt_glen);
}

// -----------------------------------------------------------------------

template<typename Alloc> 
void bvector<Alloc>::set_gap_levels(const gap_word_t* glevel_len)
{
    if (blockman_.is_init())
    {
        word_t*** blk_root = blockman_.top_blocks_root();
        typename 
            blocks_manager_type::gap_level_func  gl_func(blockman_, glevel_len);
        for_each_nzblock(blk_root, blockman_.top_block_size(),gl_func);
    }

    blockman_.set_glen(glevel_len);
}

// -----------------------------------------------------------------------

template<typename Alloc> 
int bvector<Alloc>::compare(const bvector<Alloc>& bv) const
{
    int res;
    unsigned bn = 0;
    
    unsigned top_blocks = blockman_.top_block_size();
    unsigned bvect_top_blocks = bv.blockman_.top_block_size();

    if (bvect_top_blocks > top_blocks) top_blocks = bvect_top_blocks;

    for (unsigned i = 0; i < top_blocks; ++i)
    {
        const bm::word_t* const* blk_blk = blockman_.get_topblock(i);
        const bm::word_t* const* arg_blk_blk = bv.blockman_.get_topblock(i);

        if (blk_blk == arg_blk_blk) 
        {
            bn += bm::set_array_size;
            continue;
        }

        for (unsigned j = 0; j < bm::set_array_size; ++j, ++bn)
        {
            const bm::word_t* arg_blk = arg_blk_blk ? arg_blk_blk[j] : 0;
            const bm::word_t* blk = blk_blk ? blk_blk[j] : 0;
            if (arg_blk == FULL_BLOCK_FAKE_ADDR)
                arg_blk = FULL_BLOCK_REAL_ADDR;
            if (blk == FULL_BLOCK_FAKE_ADDR)
                blk = FULL_BLOCK_REAL_ADDR;
            if (blk == arg_blk) continue;

            // If one block is zero we check if the other one has at least 
            // one bit ON

            if (!blk || !arg_blk)  
            {
                const bm::word_t*  pblk;
                bool is_gap;

                if (blk)
                {
                    pblk = blk;
                    res = 1;
                    is_gap = BM_IS_GAP(blk);
                }
                else
                {
                    pblk = arg_blk;
                    res = -1;
                    is_gap = BM_IS_GAP(arg_blk);
                }

                if (is_gap)
                {
                    if (!bm::gap_is_all_zero(BMGAP_PTR(pblk)))
                    {
                        return res;
                    }
                }
                else
                {
                    if (!bit_is_all_zero(pblk))
                    {
                        return res;
                    }
                }

                continue;
            }

            bool arg_gap = BM_IS_GAP(arg_blk);
            bool gap = BM_IS_GAP(blk);

            if (arg_gap != gap)
            {
                bm::wordop_t temp_blk[bm::set_block_size_op]; 
                bm::wordop_t* blk1;
                bm::wordop_t* blk2;

                if (gap)
                {
                    gap_convert_to_bitset((bm::word_t*)temp_blk, 
                                            BMGAP_PTR(blk));

                    blk1 = (bm::wordop_t*)temp_blk;
                    blk2 = (bm::wordop_t*)arg_blk;
                }
                else
                {
                    gap_convert_to_bitset((bm::word_t*)temp_blk, 
                                            BMGAP_PTR(arg_blk));

                    blk1 = (bm::wordop_t*)blk;
                    blk2 = (bm::wordop_t*)temp_blk;

                }                        
                res = bitcmp(blk1, blk2, bm::set_block_size_op);  

            }
            else
            {
                if (gap)
                {
                    res = gapcmp(BMGAP_PTR(blk), BMGAP_PTR(arg_blk));
                }
                else
                {
                    res = bitcmp((bm::wordop_t*)blk, 
                                    (bm::wordop_t*)arg_blk, 
                                    bm::set_block_size_op);
                }
            }

            if (res != 0)
            {
                return res;
            }
        

        } // for j

    } // for i

    return 0;
}

// -----------------------------------------------------------------------

template<typename Alloc>
void bvector<Alloc>::swap(bvector<Alloc>& bvect) BMNOEXEPT
{
    if (this != &bvect)
    {
        blockman_.swap(bvect.blockman_);
        bm::xor_swap(size_,bvect.size_);
    }
}

// -----------------------------------------------------------------------

template<typename Alloc> 
void bvector<Alloc>::calc_stat(struct bvector<Alloc>::statistics* st) const
{
    BM_ASSERT(st);
    
    st->reset();
    ::memcpy(st->gap_levels, 
             blockman_.glen(), sizeof(gap_word_t) * bm::gap_levels);

    unsigned empty_blocks = 0;

    st->max_serialize_mem = unsigned(sizeof(id_t) * 4);

    unsigned top_size = blockman_.top_block_size();
    for (unsigned i = 0; i < top_size; ++i)
    {
        const bm::word_t* const* blk_blk = blockman_.get_topblock(i);

        if (!blk_blk) 
        {
            st->max_serialize_mem += unsigned(sizeof(unsigned) + 1);
            continue;
        }

        for (unsigned j = 0; j < bm::set_array_size; ++j)
        {
            const bm::word_t* blk = blk_blk[j];
            if (IS_VALID_ADDR(blk))
            {
                st->max_serialize_mem += unsigned(empty_blocks << 2);
                empty_blocks = 0;

                if (BM_IS_GAP(blk))
                {
                    bm::gap_word_t* gap_blk = BMGAP_PTR(blk);
                    unsigned cap = bm::gap_capacity(gap_blk, blockman_.glen());
                    unsigned len = gap_length(gap_blk);
                    st->add_gap_block(cap, len);
                }
                else // bit block
                {
                    st->add_bit_block();
                }
            }
            else
            {
                ++empty_blocks;
            }
        }
    }  

    size_t safe_inc = st->max_serialize_mem / 10; // 10% increment
    if (!safe_inc) safe_inc = 256;
    st->max_serialize_mem += safe_inc;

    // Calc size of different odd and temporary things.

    st->memory_used += unsigned(sizeof(*this) - sizeof(blockman_));
    st->memory_used += blockman_.mem_used();
}

// -----------------------------------------------------------------------

template<class Alloc>
void bvector<Alloc>::set_bit_no_check(bm::id_t n)
{
    BM_ASSERT(blockman_.is_init());
    BM_ASSERT_THROW(n < bm::id_max, BM_ERR_RANGE);
    
    bool val = true; // set bit
    
    // calculate logical block number
    unsigned nblock = unsigned(n >>  bm::set_block_shift);
    // calculate word number in block and bit
    unsigned nbit   = unsigned(n & bm::set_block_mask);

    int block_type;
    bm::word_t* blk =
        blockman_.check_allocate_block(nblock,
                                       val,
                                       get_new_blocks_strat(),
                                       &block_type);
    if (!blk) return; // nothing to do

    if (block_type) // gap block
    {
        bm::gap_word_t* gap_blk = BMGAP_PTR(blk);
        gap_block_set(gap_blk, val, nblock, nbit);
    }
    else  // bit block
    {
        unsigned nword  = nbit >> bm::set_word_shift;
        nbit &= bm::set_word_mask;
        blk[nword] |= (1u << nbit); // set bit
    }
}

// -----------------------------------------------------------------------


template<class Alloc> 
bool bvector<Alloc>::set_bit_no_check(bm::id_t n, bool val)
{
    // calculate logical block number
    unsigned nblock = unsigned(n >>  bm::set_block_shift); 

    int block_type;
    bm::word_t* blk = 
        blockman_.check_allocate_block(nblock, 
                                        val,
                                        get_new_blocks_strat(), 
                                        &block_type);

    if (!blk) return false;

    // calculate word number in block and bit
    unsigned nbit   = unsigned(n & bm::set_block_mask); 

    if (block_type) // gap
    {
        bm::gap_word_t* gap_blk = BMGAP_PTR(blk);
        unsigned is_set = gap_block_set(gap_blk, val, nblock, nbit);
        return is_set;
    }
    else  // bit block
    {
        unsigned nword  = unsigned(nbit >> bm::set_word_shift); 
        nbit &= bm::set_word_mask;

        bm::word_t* word = blk + nword;
        bm::word_t  mask = (((bm::word_t)1) << nbit);

        if (val)
        {
            if ( ((*word) & mask) == 0 )
            {
                *word |= mask; // set bit
                return true;
            }
        }
        else
        {
            if ((*word) & mask)
            {
                *word &= ~mask; // clear bit
                return true;
            }
        }
    }
    return false;
}

// -----------------------------------------------------------------------

template<class Alloc>
bool bvector<Alloc>::gap_block_set(bm::gap_word_t* gap_blk,
                                   bool val, unsigned nblock, unsigned nbit)
{
    unsigned is_set, new_block_len;
    new_block_len =
        bm::gap_set_value(val, gap_blk, nbit, &is_set);
    if (is_set)
    {
        unsigned threshold =  bm::gap_limit(gap_blk, blockman_.glen());
        if (new_block_len > threshold)
        {
            extend_gap_block(nblock, gap_blk);
        }
    }
    return is_set;
}

// -----------------------------------------------------------------------

template<class Alloc>
bool bvector<Alloc>::inc(bm::id_t n)
{
    // calculate logical block number
    unsigned nblock = unsigned(n >>  bm::set_block_shift);
    bm::word_t* blk =
        blockman_.check_allocate_block(nblock,
                                       get_new_blocks_strat());
    BM_ASSERT(blk);

    unsigned nbit   = unsigned(n & bm::set_block_mask);

    unsigned is_set;
    if (BM_IS_GAP(blk))
    {
        bm::gap_word_t* gap_blk = BMGAP_PTR(blk);
        is_set = (bm::gap_test_unr(gap_blk, nbit) != 0);
        this->gap_block_set(gap_blk, !is_set, nblock, nbit); // flip
    }
    else // bit block
    {
        unsigned nword  = unsigned(nbit >> bm::set_word_shift);
        nbit &= bm::set_word_mask;

        bm::word_t* word = blk + nword;
        bm::word_t  mask = (((bm::word_t)1) << nbit);
        is_set = ((*word) & mask);
        
        *word = (is_set) ? (*word & ~mask) : (*word | mask);
    }
    return is_set;
}

// -----------------------------------------------------------------------

template<class Alloc> 
bool bvector<Alloc>::set_bit_conditional_impl(bm::id_t n, 
                                              bool     val, 
                                              bool     condition)
{
    // calculate logical block number
    unsigned nblock = unsigned(n >>  bm::set_block_shift); 

    int block_type;
    bm::word_t* blk =
        blockman_.check_allocate_block(nblock, 
                                       val,
                                       get_new_blocks_strat(), 
                                       &block_type);
    if (!blk) 
        return false;

    // calculate word number in block and bit
    unsigned nbit   = unsigned(n & bm::set_block_mask); 

    if (block_type == 1) // gap
    {
        bm::gap_word_t* gap_blk = BMGAP_PTR(blk);
        bool old_val = (bm::gap_test_unr(gap_blk, nbit) != 0);

        if (old_val != condition) 
        {
            return false;
        }

        if (val != old_val)
        {
            unsigned is_set = gap_block_set(gap_blk, val, nblock, nbit);
            BM_ASSERT(is_set);
            return is_set;
        }
    }
    else  // bit block
    {
        unsigned nword  = unsigned(nbit >> bm::set_word_shift); 
        nbit &= bm::set_word_mask;

        bm::word_t* word = blk + nword;
        bm::word_t  mask = (((bm::word_t)1) << nbit);
        bool is_set = ((*word) & mask) != 0;

        if (is_set != condition)
        {
            return false;
        }
        if (is_set != val)    // need to change bit
        {
            if (val)          // set bit
            {
                *word |= mask;
            }
            else               // clear bit
            {
                *word &= ~mask;
            }
            return true;
        }
    }
    return false;

}

// -----------------------------------------------------------------------


template<class Alloc> 
bool bvector<Alloc>::and_bit_no_check(bm::id_t n, bool val)
{
    // calculate logical block number
    unsigned nblock = unsigned(n >>  bm::set_block_shift); 

    int block_type;
    bm::word_t* blk =
        blockman_.check_allocate_block(nblock, 
                                       val,
                                       get_new_blocks_strat(), 
                                       &block_type);
    if (!blk) return false;

    // calculate word number in block and bit
    unsigned nbit   = unsigned(n & bm::set_block_mask); 

    if (block_type == 1) // gap
    {
        bm::gap_word_t* gap_blk = BMGAP_PTR(blk);
        bool old_val = (bm::gap_test_unr(gap_blk, nbit) != 0);

        bool new_val = val & old_val;
        if (new_val != old_val)
        {
            unsigned is_set = gap_block_set(gap_blk, val, nblock, nbit);
            BM_ASSERT(is_set);
            return is_set;
        }
    }
    else  // bit block
    {
        unsigned nword  = unsigned(nbit >> bm::set_word_shift); 
        nbit &= bm::set_word_mask;

        bm::word_t* word = blk + nword;
        bm::word_t  mask = (((bm::word_t)1) << nbit);
        bool is_set = ((*word) & mask) != 0;

        bool new_val = is_set & val;
        if (new_val != val)    // need to change bit
        {
            if (new_val)       // set bit
            {
                *word |= mask;
            }
            else               // clear bit
            {
                *word &= ~mask;
            }
            return true;
        }
    }
    return false;
}

//---------------------------------------------------------------------

template<class Alloc>
bool bvector<Alloc>::find(bm::id_t from, bm::id_t& pos) const
{
    BM_ASSERT_THROW(from < bm::id_max, BM_ERR_RANGE);

    if (from == 0)
    {
        return find(pos);
    }
    pos = check_or_next(from);
    return (pos != 0);
}

//---------------------------------------------------------------------

template<class Alloc>
bool bvector<Alloc>::find_reverse(bm::id_t& pos) const
{
    bool found;
    
    unsigned top_blocks = blockman_.top_block_size();
    for (unsigned i = top_blocks-1; true; --i)
    {
        const bm::word_t* const* blk_blk = blockman_.get_topblock(i);
        if (blk_blk)
        {
            for (unsigned short j = bm::set_array_size-1; true; --j)
            {
                const bm::word_t* blk = blk_blk[j];
                if (blk)
                {
                    if (blk == FULL_BLOCK_FAKE_ADDR)
                        blk = FULL_BLOCK_REAL_ADDR;
                    
                    bool is_gap = BM_IS_GAP(blk);
                    found = is_gap ? bm::gap_find_last(BMGAP_PTR(blk), &pos)
                                   : bm::bit_find_last(blk, &pos);
                    if (found)
                    {
                        unsigned base_idx = i * bm::set_array_size * bm::gap_max_bits;
                        base_idx += j * bm::gap_max_bits;
                        pos += base_idx;
                        return found;
                    }
                }
                
                if (j == 0)
                    break;
            } // for j
        } // if blk_blk
        
        if (i == 0)
            break;
    } // for i
    return false;
}

//---------------------------------------------------------------------

template<class Alloc>
bool bvector<Alloc>::find(bm::id_t& pos) const
{
    bool found;
    
    unsigned top_blocks = blockman_.top_block_size();
    for (unsigned short i = 0; i < top_blocks; ++i)
    {
        const bm::word_t* const* blk_blk = blockman_.get_topblock(i);
        if (blk_blk)
        {
            for (unsigned short j = 0; j < bm::set_array_size; ++j)
            {
                const bm::word_t* blk = blk_blk[j];
                if (blk)
                {
                    if (blk == FULL_BLOCK_FAKE_ADDR)
                    {
                        found = true; pos = 0;
                    }
                    else
                    {
                        bool is_gap = BM_IS_GAP(blk);
                        found = (is_gap) ? bm::gap_find_first(BMGAP_PTR(blk), &pos)
                                         : bm::bit_find_first(blk, &pos);
                    }
                    if (found)
                    {
                        unsigned base_idx = i * bm::set_array_size * bm::gap_max_bits;
                        base_idx += j * bm::gap_max_bits;
                        pos += base_idx;
                        return found;
                    }
                }
            } // for j
        } // if blk_blk
    } // for i
    return false;
}

//---------------------------------------------------------------------

template<class Alloc>
bool bvector<Alloc>::find_range(bm::id_t& in_first, bm::id_t& in_last) const
{
    bool found = find(in_first);
    if (found)
    {
        found = find_reverse(in_last);
        BM_ASSERT(found);
    }
    return found;
}

//---------------------------------------------------------------------

template<class Alloc>
bool bvector<Alloc>::find_rank(bm::id_t rank, bm::id_t from, bm::id_t& pos) const
{
    BM_ASSERT_THROW(from < bm::id_max, BM_ERR_RANGE);

    bool ret = false;
    
    if (!rank || !blockman_.is_init())
        return ret;
    
    unsigned nb  = unsigned(from  >>  bm::set_block_shift);
    unsigned nb_right  = unsigned((bm::id_max-1)  >>  bm::set_block_shift);
    bm::gap_word_t nbit = bm::gap_word_t(from & bm::set_block_mask);
    unsigned bit_pos = 0;

    for (; nb < nb_right; ++nb)
    {
        int no_more_blocks;
        const bm::word_t* block = blockman_.get_block(nb, &no_more_blocks);
        if (block)
        {
            if (BM_IS_GAP(block))
            {
                const bm::gap_word_t* const gap_block = BMGAP_PTR(block);
                rank = bm::gap_find_rank(gap_block, rank, nbit, bit_pos);
            }
            else
            {
                rank = bm::bit_find_rank(block, rank, nbit, bit_pos);
            }
            
            if (!rank) // target found
            {
                bm::id_t prev_pos = nb ? (nb * bm::set_block_size * 32) : 0;
                pos = bit_pos + prev_pos;
                return true;
            }
        }
        else
        {
            if (no_more_blocks)
                break;
        }
        nbit = 0;
    } // for nb
    
    return ret;
}

//---------------------------------------------------------------------

template<class Alloc> 
bm::id_t bvector<Alloc>::check_or_next(bm::id_t prev) const
{
    if (!blockman_.is_init())
        return 0;
    for (;;)
    {
        unsigned nblock = unsigned(prev >> bm::set_block_shift); 
        if (nblock >= bm::set_total_blocks) 
            break;

        if (blockman_.is_subblock_null(nblock >> bm::set_array_shift))
        {
            prev += (bm::set_blkblk_mask + 1) - (prev & bm::set_blkblk_mask);
        }
        else
        {
            unsigned nbit = unsigned(prev & bm::set_block_mask);
            int no_more_blocks;
            const bm::word_t* block = 
                blockman_.get_block(nblock, &no_more_blocks);

            if (no_more_blocks) 
            {
                BM_ASSERT(block == 0);
                break;
            }

            if (block)
            {
                if (BM_IS_GAP(block))
                {
                    if (bm::gap_find_in_block(BMGAP_PTR(block),
                                                nbit,
                                                &prev))
                    {
                        return prev;
                    }
                }
                else
                {
                    if (IS_FULL_BLOCK(block)) return prev;
                    if (bm::bit_find_in_block(block, nbit, &prev))
                    {
                        return prev;
                    }
                }
            }
            else
            {
                prev += (bm::set_block_mask + 1) - (prev & bm::set_block_mask);
            }

        }
        if (!prev) break;
    }

    return 0;
}


//---------------------------------------------------------------------

template<class Alloc> 
bm::id_t bvector<Alloc>::check_or_next_extract(bm::id_t prev)
{
    if (!blockman_.is_init())
        return 0;

    for (;;)
    {
        unsigned nblock = unsigned(prev >> bm::set_block_shift); 
        if (nblock >= bm::set_total_blocks) break;

        if (blockman_.is_subblock_null(nblock >> bm::set_array_shift))
        {
            prev += (bm::set_blkblk_mask + 1) -
                            (prev & bm::set_blkblk_mask);
        }
        else
        {
            unsigned nbit = unsigned(prev & bm::set_block_mask);

            int no_more_blocks;
            bm::word_t* block = 
                blockman_.get_block(nblock, &no_more_blocks);

            if (no_more_blocks) 
            {
                BM_ASSERT(block == 0);
                break;
            }


            if (block)
            {
                if (IS_FULL_BLOCK(block))
                {
                    set(prev, false);
                    return prev;
                }
                if (BM_IS_GAP(block))
                {
                    unsigned is_set;
                    unsigned new_block_len = 
                        gap_set_value(0, BMGAP_PTR(block), nbit, &is_set);
                    if (is_set)
                    {
                        unsigned threshold =
                            bm::gap_limit(BMGAP_PTR(block), blockman_.glen());
                        if (new_block_len > threshold) 
                        {
                            extend_gap_block(nblock, BMGAP_PTR(block));
                        }
                        return prev;
                    }
                    else
                    {
                        if (bm::gap_find_in_block(BMGAP_PTR(block),
                                                    nbit,
                                                    &prev))
                        {
                            set(prev, false);
                            return prev;
                        }
                    }
                }
                else // bit block
                {
                    if (bm::bit_find_in_block(block, nbit, &prev)) 
                    {
                        unsigned nbit1 =
                            unsigned(prev & bm::set_block_mask); 
                        unsigned nword = 
                            unsigned(nbit1 >> bm::set_word_shift);
                        nbit1 &= bm::set_word_mask;
                        bm::word_t* word = block + nword;
                        bm::word_t  mask = ((bm::word_t)1) << nbit1;
                        *word &= ~mask;

                        return prev;
                    }
                }
            }
            else
            {
                prev += (bm::set_block_mask + 1) - 
                            (prev & bm::set_block_mask);
            }

        }
        if (!prev) break;
    }

    return 0;
}

//---------------------------------------------------------------------

#define BM_OR_OP(x)  \
    { \
        blk = blk_blk[j+x]; arg_blk = blk_blk_arg[j+x]; \
        if (blk != arg_blk) \
            combine_operation_block_or(i, j+x, blk, arg_blk); \
    }

template<class Alloc>
void bvector<Alloc>::combine_operation_or(const bm::bvector<Alloc>& bv)
{
    if (!bv.blockman_.is_init())
        return;

    unsigned top_blocks = blockman_.top_block_size();
    if (size_ < bv.size_) // this vect shorter than the arg.
    {
        size_ = bv.size_;
    }
    unsigned arg_top_blocks = bv.blockman_.top_block_size();
    top_blocks = blockman_.reserve_top_blocks(arg_top_blocks);

    
    bm::word_t*** blk_root = blockman_.top_blocks_root();
    bm::word_t*** blk_root_arg = bv.blockman_.top_blocks_root();

    for (unsigned i = 0; i < top_blocks; ++i)
    {
        bm::word_t** blk_blk = blk_root[i];
        bm::word_t** blk_blk_arg = (i < arg_top_blocks) ? blk_root_arg[i] : 0;
        if (blk_blk == blk_blk_arg || !blk_blk_arg) // nothing to do (0 OR 0 == 0)
            continue;
        if (!blk_blk)
            blk_blk = blockman_.alloc_top_subblock(i);

        unsigned j = 0;
        bm::word_t* blk;
        const bm::word_t* arg_blk;
        do
        {
        #if defined(BM64_AVX2)
            BM_OR_OP(0)
            BM_OR_OP(1)
            BM_OR_OP(2)
            BM_OR_OP(3)
            j += 4;
        #elif defined(BM64_SSE4)
            BM_OR_OP(0)
            BM_OR_OP(1)
            j += 2;
        #else
            BM_OR_OP(0)
            ++j;
        #endif
        } while (j < bm::set_array_size);
    } // for i
}

#undef BM_OR_OP


//---------------------------------------------------------------------

#define BM_AND_OP(x)  if (0 != (blk = blk_blk[j+x])) \
    { \
        if (0 != (arg_blk = blk_blk_arg[j+x])) \
            combine_operation_block_and(i, j+x, blk, arg_blk); \
        else \
            blockman_.zero_block(i, j+x); \
    }

template<class Alloc>
void bvector<Alloc>::combine_operation_and(const bm::bvector<Alloc>& bv)
{
    if (!blockman_.is_init())
        return;  // nothing to do, already empty
    if (!bv.blockman_.is_init())
    {
        clear(true);
        return;
    }

    unsigned top_blocks = blockman_.top_block_size();
    if (size_ < bv.size_) // this vect shorter than the arg.
    {
        size_ = bv.size_;
    }
    unsigned arg_top_blocks = bv.blockman_.top_block_size();
    top_blocks = blockman_.reserve_top_blocks(arg_top_blocks);

    
    bm::word_t*** blk_root = blockman_.top_blocks_root();
    bm::word_t*** blk_root_arg = bv.blockman_.top_blocks_root();

    for (unsigned i = 0; i < top_blocks; ++i)
    {
        bm::word_t** blk_blk = blk_root[i];
        if (!blk_blk) // nothing to do (0 AND 1 == 0)
            continue;
        bm::word_t** blk_blk_arg = (i < arg_top_blocks) ? blk_root_arg[i] : 0;
        if (!blk_blk_arg) // free a whole group of blocks
        {
            for (unsigned j = 0; j < bm::set_array_size; ++j)
                blockman_.zero_block(i, j);
            continue;
        }
        unsigned j = 0;
        bm::word_t* blk;
        const bm::word_t* arg_blk;
        do
        {
        #ifdef BM64_AVX2
            if (!avx2_test_all_zero_wave(blk_blk + j))
            {
                BM_AND_OP(0)
                BM_AND_OP(1)
                BM_AND_OP(2)
                BM_AND_OP(3)
            }
            j += 4;
        #elif defined(BM64_SSE4)
            if (!sse42_test_all_zero_wave(blk_blk + j))
            {
                BM_AND_OP(0)
                BM_AND_OP(1)
            }
            j += 2;
        #else
            BM_AND_OP(0)
            ++j;
        #endif
        } while (j < bm::set_array_size);
    } // for i
}

#undef BM_AND_OP


//---------------------------------------------------------------------

#define BM_SUB_OP(x) \
    if ((0 != (blk = blk_blk[j+x])) && (0 != (arg_blk = blk_blk_arg[j+x]))) \
        combine_operation_block_sub(i, j+x, blk, arg_blk);


template<class Alloc>
void bvector<Alloc>::combine_operation_sub(const bm::bvector<Alloc>& bv)
{
    if (!blockman_.is_init())
        return;  
    if (!bv.blockman_.is_init())
        return;

    unsigned top_blocks = blockman_.top_block_size();
    if (size_ < bv.size_) // this vect shorter than the arg.
    {
        size_ = bv.size_;
    }
    unsigned arg_top_blocks = bv.blockman_.top_block_size();
    top_blocks = blockman_.reserve_top_blocks(arg_top_blocks);

    bm::word_t*** blk_root = blockman_.top_blocks_root();
    bm::word_t*** blk_root_arg = bv.blockman_.top_blocks_root();

    for (unsigned i = 0; i < top_blocks; ++i)
    {
        bm::word_t** blk_blk = blk_root[i];
        bm::word_t** blk_blk_arg = (i < arg_top_blocks) ? blk_root_arg[i] : 0;
        if (!blk_blk || !blk_blk_arg) // nothing to do (0 AND NOT 1 == 0)
            continue;
        bm::word_t* blk;
        const bm::word_t* arg_blk;
        unsigned j = 0;
        do
        {
        #ifdef BM64_AVX2
            if (!avx2_test_all_zero_wave(blk_blk + j))
            {
                BM_SUB_OP(0)
                BM_SUB_OP(1)
                BM_SUB_OP(2)
                BM_SUB_OP(3)
            }
            j += 4;
        #elif defined(BM64_SSE4)
            if (!sse42_test_all_zero_wave(blk_blk + j))
            {
                BM_SUB_OP(0)
                BM_SUB_OP(1)
            }
            j += 2;
        #else
            BM_SUB_OP(0)
            ++j;
        #endif
        } while (j < bm::set_array_size);
    } // for i
}

#undef BM_SUB_OP

//---------------------------------------------------------------------

template<class Alloc> 
void bvector<Alloc>::combine_operation(
                                  const bm::bvector<Alloc>& bv,
                                  bm::operation             opcode)
{
    if (!blockman_.is_init())
    {
        if (opcode == BM_AND || opcode == BM_SUB)
        {
            return;
        }
        blockman_.init_tree();
    }

    unsigned top_blocks = blockman_.top_block_size();
    unsigned arg_top_blocks = bv.blockman_.top_block_size();
    
    if (arg_top_blocks > top_blocks)
        top_blocks = blockman_.reserve_top_blocks(arg_top_blocks);

    if (size_ < bv.size_) // this vect shorter than the arg.
    {
        size_ = bv.size_;
        // stretch our capacity
        blockman_.reserve_top_blocks(arg_top_blocks);
        top_blocks = blockman_.top_block_size();
    }
    else 
    if (size_ > bv.size_) // this vector larger
    {
        if (opcode == BM_AND) // clear the tail with zeros
        {
            set_range(bv.size_, size_ - 1, false);
            if (arg_top_blocks < top_blocks)
            {
                // not to scan blocks we already swiped
                top_blocks = arg_top_blocks;
            }
        }
    }
    
    bm::word_t*** blk_root = blockman_.top_blocks_root();
    unsigned block_idx = 0;
    unsigned i, j;

    // calculate effective top size to avoid overscan
    top_blocks = blockman_.top_block_size();
    if (top_blocks < bv.blockman_.top_block_size())
    {
        if (opcode != BM_AND)
        {
            top_blocks = bv.blockman_.top_block_size();
        }
    }

    for (i = 0; i < top_blocks; ++i)
    {
        bm::word_t** blk_blk = blk_root[i];
        if (blk_blk == 0) // not allocated
        {
            if (opcode == BM_AND) // 0 AND anything == 0
            {
                block_idx += bm::set_array_size;
                continue; 
            }
            const bm::word_t* const* bvbb = bv.blockman_.get_topblock(i);
            if (bvbb == 0) // skip it because 0 OP 0 == 0 
            {
                block_idx += bm::set_array_size;
                continue; 
            }
            // 0 - self, non-zero argument
            unsigned r = i * bm::set_array_size;
            for (j = 0; j < bm::set_array_size; ++j)
            {
                const bm::word_t* arg_blk = bv.blockman_.get_block(i, j);
                if (arg_blk )
                    combine_operation_with_block(r + j,
                                                 0, 0, 
                                                 arg_blk, BM_IS_GAP(arg_blk), 
                                                 opcode);
            } // for j
            continue;
        }

        if (opcode == BM_AND)
        {
            unsigned r = i * bm::set_array_size;
            for (j = 0; j < bm::set_array_size; ++j)
            {            
                bm::word_t* blk = blk_blk[j];
                if (blk)
                {
                    const bm::word_t* arg_blk = bv.blockman_.get_block(i, j);
                    if (arg_blk)
                        combine_operation_with_block(r + j,
                                                     BM_IS_GAP(blk), blk, 
                                                     arg_blk, BM_IS_GAP(arg_blk),
                                                     opcode);                    
                    else
                        blockman_.zero_block(i, j);
                }

            } // for j
        }
        else // OR, SUB, XOR
        {
            unsigned r = i * bm::set_array_size;
            for (j = 0; j < bm::set_array_size; ++j)
            {            
                bm::word_t* blk = blk_blk[j];
                const bm::word_t* arg_blk = bv.blockman_.get_block(i, j);
                if (arg_blk || blk)
                    combine_operation_with_block(r + j, BM_IS_GAP(blk), blk, 
                                                 arg_blk, BM_IS_GAP(arg_blk),
                                                 opcode);
            } // for j
        }
    } // for i

}

//---------------------------------------------------------------------

template<class Alloc>
void bvector<Alloc>::combine_operation_block_or(
        unsigned i, unsigned j,
        bm::word_t* blk, const bm::word_t* arg_blk)
{
    bm::gap_word_t tmp_buf[bm::gap_equiv_len * 3]; // temporary result
    if (IS_FULL_BLOCK(blk) || !arg_blk) // all bits are set
        return; // nothing to do
    
    if (IS_FULL_BLOCK(arg_blk))
    {
        if (blk)
            blockman_.zero_block(i, j); // free target block and assign FULL
        blockman_.set_block_ptr(i, j, FULL_BLOCK_FAKE_ADDR);
        return;
    }
    
    if (BM_IS_GAP(blk)) // our block GAP-type
    {
        bm::gap_word_t* gap_blk = BMGAP_PTR(blk);
        if (BM_IS_GAP(arg_blk))  // both blocks GAP-type
        {
            const bm::gap_word_t* res;
            unsigned res_len;
            res = bm::gap_operation_or(gap_blk,
                                        BMGAP_PTR(arg_blk),
                                        tmp_buf,
                                        res_len);
            BM_ASSERT(res == tmp_buf);
            if (bm::gap_is_all_one(res, bm::gap_max_bits))
            {
                blockman_.zero_gap_block_ptr(i, j);
                blockman_.set_block_ptr(i, j, FULL_BLOCK_FAKE_ADDR);
            }
            else
            {
                assign_gap_result(i, j, res, ++res_len, blk, tmp_buf);
            }
            return;
        }
        // GAP or BIT
        //
        bm::word_t* new_blk = blockman_.get_allocator().alloc_bit_block();
        bm::bit_block_copy(new_blk, arg_blk); // TODO: 3 way operation
        bm::gap_add_to_bitset(new_blk, gap_blk);
        
        blockman_.get_allocator().free_gap_block(gap_blk, blockman_.glen());
        blockman_.set_block_ptr(i, j, new_blk);
        
        return;
    }
    else // our block is BITSET-type (but NOT FULL_BLOCK we checked it)
    {
        if (BM_IS_GAP(arg_blk)) // argument block is GAP-type
        {
            const bm::gap_word_t* arg_gap = BMGAP_PTR(arg_blk);
            if (!blk)
            {
                bool gap = true;
                blk = blockman_.clone_gap_block(arg_gap, gap);
                blockman_.set_block(i, j, blk, gap);
                return;
            }

            // BIT & GAP
            bm::gap_add_to_bitset(blk, arg_gap);
            return;
        } // if arg_gap
    }
    
    if (!blk)
    {
        blk = blockman_.alloc_.alloc_bit_block();
        bm::bit_block_copy(blk, arg_blk);
        blockman_.set_block_ptr(i, j, blk);
        return;
    }

    bool all_one = bm::bit_block_or(blk, arg_blk);
    if (all_one)
    {
        BM_ASSERT(bm::is_bits_one((bm::wordop_t*) blk, (bm::wordop_t*) (blk + bm::set_block_size)));
        blockman_.get_allocator().free_bit_block(blk);
        blockman_.set_block_ptr(i, j, FULL_BLOCK_FAKE_ADDR);
    }

}


//---------------------------------------------------------------------


template<class Alloc>
void bvector<Alloc>::combine_operation_block_and(
        unsigned i, unsigned j, bm::word_t* blk, const bm::word_t* arg_blk)
{
    BM_ASSERT(arg_blk && blk);
    
    if (IS_FULL_BLOCK(arg_blk))
        return;  // nothing to do
    
    gap_word_t tmp_buf[bm::gap_equiv_len * 3]; // temporary result
    
    if (BM_IS_GAP(blk)) // our block GAP-type
    {
        bm::gap_word_t* gap_blk = BMGAP_PTR(blk);
        if (BM_IS_GAP(arg_blk))  // both blocks GAP-type
        {
            const bm::gap_word_t* res;
            unsigned res_len;
            res = bm::gap_operation_and(gap_blk,
                                        BMGAP_PTR(arg_blk),
                                        tmp_buf,
                                        res_len);
            BM_ASSERT(res == tmp_buf);

            if (bm::gap_is_all_zero(res))
                blockman_.zero_gap_block_ptr(i, j);
            else
                assign_gap_result(i, j, res, ++res_len, blk, tmp_buf);
            return;
        }
        // GAP & BIT
        //
        bm::word_t* new_blk = blockman_.get_allocator().alloc_bit_block();
        bm::bit_block_copy(new_blk, arg_blk); // TODO: 3 way operation
        bm::gap_and_to_bitset(new_blk, gap_blk);
        
        bool empty = bm::bit_is_all_zero(new_blk);
        if (empty) // operation converged bit-block to empty
        {
            blockman_.get_allocator().free_bit_block(new_blk);
            new_blk = 0;
        }
        blockman_.get_allocator().free_gap_block(gap_blk, blockman_.glen());
        blockman_.set_block_ptr(i, j, new_blk);
        
        return;
    }
    else // our block is BITSET-type or FULL_BLOCK
    {
        if (BM_IS_GAP(arg_blk)) // argument block is GAP-type
        {
            const bm::gap_word_t* arg_gap = BMGAP_PTR(arg_blk);
            if (bm::gap_is_all_zero(arg_gap))
            {
                blockman_.zero_block(i, j);
                return;
            }
            // FULL & GAP is common when AND with set_range() mask
            //
            if (IS_FULL_BLOCK(blk)) // FULL & gap = gap
            {
                bool  is_new_gap;
                bm::word_t* new_blk = blockman_.clone_gap_block(arg_gap, is_new_gap);
                if (is_new_gap)
                    BMSET_PTRGAP(new_blk);
                
                blockman_.set_block_ptr(i, j, new_blk);
                
                return;
            }
            // BIT & GAP
            bm::gap_and_to_bitset(blk, arg_gap);
            bool empty = bm::bit_is_all_zero(blk);
            if (empty) // operation converged bit-block to empty
                blockman_.zero_block(i, j);
            
            return;
        } // if arg_gap
    }
    
    // FULL & bit is common when AND with set_range() mask
    //
    if (IS_FULL_BLOCK(blk)) // FULL & bit = bit
    {
        bm::word_t* new_blk = blockman_.get_allocator().alloc_bit_block();
        bm::bit_block_copy(new_blk, arg_blk);
        blockman_.set_block_ptr(i, j, new_blk);
        
        return;
    }
    auto any_bits = bm::bit_block_and(blk, arg_blk);
    if (!any_bits)
    {
        blockman_.get_allocator().free_bit_block(blk);
        blockman_.set_block_ptr(i, j, 0);
    }
}

//---------------------------------------------------------------------

template<class Alloc>
void bvector<Alloc>::combine_operation_block_sub(
                unsigned i, unsigned j, bm::word_t* blk, const bm::word_t* arg_blk)
{
    BM_ASSERT(arg_blk && blk);

    if (IS_FULL_BLOCK(arg_blk))
    {
        blockman_.zero_block(i, j);
        return;
    }

    gap_word_t tmp_buf[bm::gap_equiv_len * 3]; // temporary result
    
    if (BM_IS_GAP(blk)) // our block GAP-type
    {
        bm::gap_word_t* gap_blk = BMGAP_PTR(blk);
        if (BM_IS_GAP(arg_blk))  // both blocks GAP-type
        {
            const bm::gap_word_t* res;
            unsigned res_len;
            res = bm::gap_operation_sub(gap_blk,
                                        BMGAP_PTR(arg_blk),
                                        tmp_buf,
                                        res_len);

            BM_ASSERT(res == tmp_buf);
            BM_ASSERT(!(res == tmp_buf && res_len == 0));

            if (bm::gap_is_all_zero(res))
                blockman_.zero_gap_block_ptr(i, j);
            else
                assign_gap_result(i, j, res, ++res_len, blk, tmp_buf);
            return;
        }
        // else: argument is BITSET-type (own block is GAP)
        blk = blockman_.convert_gap2bitset(i, j, gap_blk);
        // fall through to bit-block to bit-block operation
    }
    else // our block is BITSET-type or FULL_BLOCK
    {
        if (BM_IS_GAP(arg_blk)) // argument block is GAP-type
        {
            if (!IS_FULL_BLOCK(blk))  // gap combined to bitset
            {
                bm::gap_sub_to_bitset(blk, BMGAP_PTR(arg_blk));
                bool empty = bm::bit_is_all_zero(blk);
                if (empty) // operation converged bit-block to empty
                    blockman_.zero_block(i, j);
                return;
            }
            // the worst case: convert arg block to bitset
            arg_blk =
                gap_convert_to_bitset_smart(blockman_.check_allocate_tempblock(),
                                            BMGAP_PTR(arg_blk),
                                            bm::gap_max_bits);
        }
    }

    // Now here we combine two plain bitblocks using supplied bit function.
    bm::word_t* dst = blk;

    bm::word_t* ret;
    if (!dst || !arg_blk)
        return;

    ret = bm::bit_operation_sub(dst, arg_blk);
    if (ret && ret == arg_blk)
    {
        ret = blockman_.get_allocator().alloc_bit_block();
        bm::bit_andnot_arr_ffmask(ret, arg_blk);
    }

    if (ret != dst) // block mutation
    {
        blockman_.set_block_ptr(i, j, ret);
        if (IS_VALID_ADDR(dst))
            blockman_.get_allocator().free_bit_block(dst);
    }
}

//---------------------------------------------------------------------

template<class Alloc> 
void 
bvector<Alloc>::combine_operation_with_block(unsigned          nb,
                                             bool              gap,
                                             bm::word_t*       blk,
                                             const bm::word_t* arg_blk,
                                             bool              arg_gap,
                                             bm::operation     opcode)
{
    gap_word_t tmp_buf[bm::gap_equiv_len * 3]; // temporary result            
    const bm::gap_word_t* res;
    unsigned res_len;

    if (opcode == BM_OR || opcode == BM_XOR)
    {        
        if (!blk && arg_gap) 
        {
            blk = blockman_.clone_gap_block(BMGAP_PTR(arg_blk), gap);
            blockman_.set_block(nb, blk, gap);
            return;
        }
    }

    if (gap) // our block GAP-type
    {
        if (arg_gap)  // both blocks GAP-type
        {
            {
                gap_operation_func_type gfunc =
                    operation_functions<true>::gap_operation(opcode);
                BM_ASSERT(gfunc);
                res = (*gfunc)(BMGAP_PTR(blk),
                               BMGAP_PTR(arg_blk),
                               tmp_buf,
                               res_len);
            }
            BM_ASSERT(res == tmp_buf);
            BM_ASSERT(!(res == tmp_buf && res_len == 0));

            // if as a result of the operation gap block turned to zero
            // we can now replace it with NULL
            if (gap_is_all_zero(res))
                blockman_.zero_block(nb);
            else
                assign_gap_result(nb, res, ++res_len,  blk, tmp_buf);
            return;
        }
        else // argument is BITSET-type (own block is GAP)
        {
            // since we can not combine blocks of mixed type
            // we need to convert our block to bitset
           
            if (arg_blk == 0)  // Combining against an empty block
            {
                switch (opcode)
                {
                case BM_AND:  // ("Value" AND  0) == 0
                    blockman_.zero_block(nb);
                    return;
                case BM_OR: case BM_SUB: case BM_XOR:
                    return; // nothing to do
                default:
                    return; // nothing to do
                }
            }
            gap_word_t* gap_blk = BMGAP_PTR(blk);

            blk = blockman_.convert_gap2bitset(nb, gap_blk);
        }
    }
    else // our block is BITSET-type
    {
        if (arg_gap) // argument block is GAP-type
        {
            if (IS_VALID_ADDR(blk))
            {
                // special case, maybe we can do the job without
                // converting the GAP argument to bitblock
                gap_operation_to_bitset_func_type gfunc =
                    operation_functions<true>::gap_op_to_bit(opcode);
                BM_ASSERT(gfunc);
                (*gfunc)(blk, BMGAP_PTR(arg_blk));

                if (opcode != BM_OR)
                {
                    bool b = bm::bit_is_all_zero(blk);
                    if (b) // operation converged bit-block to empty
                        blockman_.zero_block(nb);
                }

                return;
            }
            
            // the worst case we need to convert argument block to
            // bitset type.
            gap_word_t* temp_blk = (gap_word_t*) blockman_.check_allocate_tempblock();
            arg_blk =
                gap_convert_to_bitset_smart((bm::word_t*)temp_blk,
                                            BMGAP_PTR(arg_blk),
                                            bm::gap_max_bits);
        
        }
    }

    // Now here we combine two plain bitblocks using supplied bit function.
    bm::word_t* dst = blk;

    bm::word_t* ret;
    if (dst == 0 && arg_blk == 0)
    {
        return;
    }

    switch (opcode)
    {
    case BM_AND:
        ret = bit_operation_and(dst, arg_blk);
        goto copy_block;
    case BM_XOR:
        ret = bit_operation_xor(dst, arg_blk);
        if (ret && (ret == arg_blk) && IS_FULL_BLOCK(dst))
        {
            ret = blockman_.get_allocator().alloc_bit_block();
#ifdef BMVECTOPT
        VECT_XOR_ARR_2_MASK(ret,
                            arg_blk,
                            arg_blk + bm::set_block_size,
                            ~0u);
#else
            bm::wordop_t* dst_ptr = (wordop_t*)ret;
            const bm::wordop_t* wrd_ptr = (wordop_t*) arg_blk;
            const bm::wordop_t* wrd_end =
            (wordop_t*) (arg_blk + bm::set_block_size);

            do
            {
                dst_ptr[0] = bm::all_bits_mask ^ wrd_ptr[0];
                dst_ptr[1] = bm::all_bits_mask ^ wrd_ptr[1];
                dst_ptr[2] = bm::all_bits_mask ^ wrd_ptr[2];
                dst_ptr[3] = bm::all_bits_mask ^ wrd_ptr[3];

                dst_ptr+=4;
                wrd_ptr+=4;

            } while (wrd_ptr < wrd_end);
#endif
            break;
        }
        goto copy_block;
    case BM_OR:
        ret = bit_operation_or(dst, arg_blk);
    copy_block:
        if (ret && (ret == arg_blk) && !IS_FULL_BLOCK(ret))
        {
        ret = blockman_.get_allocator().alloc_bit_block();
        bit_block_copy(ret, arg_blk);
        }
        break;

    case BM_SUB:
        ret = bit_operation_sub(dst, arg_blk);
        if (ret && ret == arg_blk)
        {
            ret = blockman_.get_allocator().alloc_bit_block();
            bm::bit_andnot_arr_ffmask(ret, arg_blk);
        }
        break;
    default:
        BM_ASSERT(0);
        ret = 0;
    }

    if (ret != dst) // block mutation
    {
        blockman_.set_block(nb, ret);
        if (IS_VALID_ADDR(dst))
        {
            blockman_.get_allocator().free_bit_block(dst);
        }
    }
}

//---------------------------------------------------------------------

template<class Alloc>
void bvector<Alloc>::assign_gap_result(
                       unsigned              nb,
                       const bm::gap_word_t* res,
                       unsigned              res_len,
                       bm::word_t*           blk,
                       gap_word_t*           tmp_buf)
{
    int level = bm::gap_level(BMGAP_PTR(blk));
    BM_ASSERT(level >= 0);
    unsigned threshold = unsigned(blockman_.glen(unsigned(level)) - 4u);


    int new_level = bm::gap_calc_level(res_len, blockman_.glen());
    if (new_level < 0)
    {
        blockman_.convert_gap2bitset(nb, res);
        return;
    }

    if (res_len > threshold)
    {
        BM_ASSERT(new_level >= 0);
        gap_word_t* new_blk =
            blockman_.allocate_gap_block(unsigned(new_level), res);
        bm::set_gap_level(new_blk, new_level);

        bm::word_t* p = (bm::word_t*)new_blk;
        BMSET_PTRGAP(p);

        if (blk)
        {
            blockman_.set_block_ptr(nb, p);
            blockman_.get_allocator().free_gap_block(BMGAP_PTR(blk),
                                                     blockman_.glen());
        }
        else
        {
            blockman_.set_block(nb, p, true); // set GAP block
        }
        return;
    }

    // gap operation result is in the temporary buffer
    // we copy it back to the gap_block

    BM_ASSERT(blk);

    bm::set_gap_level(tmp_buf, level);
    ::memcpy(BMGAP_PTR(blk), tmp_buf, res_len * sizeof(gap_word_t));
}


//---------------------------------------------------------------------

template<class Alloc>
void bvector<Alloc>::assign_gap_result(
                       unsigned              i,
                       unsigned              j,
                       const bm::gap_word_t* res,
                       unsigned              res_len,
                       bm::word_t*           blk,
                       gap_word_t*           tmp_buf)
{
    int level = bm::gap_level(BMGAP_PTR(blk));
    BM_ASSERT(level >= 0);
    unsigned threshold = unsigned(blockman_.glen(unsigned(level)) - 4u);

    int new_level = bm::gap_calc_level(res_len, blockman_.glen());
    if (new_level < 0)
    {
        blockman_.convert_gap2bitset(i, j, res);
        return;
    }

    if (res_len > threshold)
    {
        BM_ASSERT(new_level >= 0);
        gap_word_t* new_blk =
            blockman_.allocate_gap_block(unsigned(new_level), res);
        bm::set_gap_level(new_blk, new_level);

        bm::word_t* p = (bm::word_t*)new_blk;
        BMSET_PTRGAP(p);

        if (blk)
        {
            blockman_.set_block_ptr(i, j, p);
            blockman_.get_allocator().free_gap_block(BMGAP_PTR(blk),
                                                     blockman_.glen());
        }
        else
        {
            blockman_.set_block(i, j, p, true); // set GAP block
        }
        return;
    }

    // gap operation result is in the temporary buffer
    // we copy it back to the gap_block

    BM_ASSERT(blk);

    bm::set_gap_level(tmp_buf, level);
    ::memcpy(BMGAP_PTR(blk), tmp_buf, res_len * sizeof(gap_word_t));
}


//---------------------------------------------------------------------

template<class Alloc> 
void bvector<Alloc>::set_range_no_check(bm::id_t left,
                                        bm::id_t right)
{
    unsigned nblock_left  = unsigned(left  >>  bm::set_block_shift);
    unsigned nblock_right = unsigned(right >>  bm::set_block_shift);

    unsigned nbit_right = unsigned(right & bm::set_block_mask);

    unsigned r = 
        (nblock_left == nblock_right) ? nbit_right :(bm::bits_in_block-1);

    bm::gap_word_t tmp_gap_blk[5];
    tmp_gap_blk[0] = 0; // just to silence GCC warning on uninit var

    // Set bits in the starting block

    unsigned nb;
    bm::word_t* block; //= blockman_.get_block(nblock_left);
    unsigned nbit_left  = unsigned(left  & bm::set_block_mask);
    if ((nbit_left == 0) && (r == bm::bits_in_block - 1)) // full block
    {
        nb = nblock_left;
    }
    else
    {
        gap_init_range_block<gap_word_t>(tmp_gap_blk,
                                         (gap_word_t)nbit_left, 
                                         (gap_word_t)r, 
                                         (gap_word_t)1);
        block = blockman_.get_block(nblock_left);
        combine_operation_with_block(nblock_left,
            BM_IS_GAP(block),
            block,
            (bm::word_t*) tmp_gap_blk,
            1, BM_OR);

        if (nblock_left == nblock_right)  // in one block
            return;
        nb = nblock_left+1;
    }

    // Set all full blocks between left and right
    //
    unsigned nb_to = nblock_right + (nbit_right ==(bm::bits_in_block-1));
    for (; nb < nb_to; ++nb)
    {
        block = blockman_.get_block(nb);
        if (IS_FULL_BLOCK(block)) 
            continue;
        blockman_.set_block_all_set(nb);            
    } // for

    if (nb_to > nblock_right)
        return;

    block = blockman_.get_block(nblock_right);

    gap_init_range_block<gap_word_t>(tmp_gap_blk, 
                                     (gap_word_t)0, 
                                     (gap_word_t)nbit_right, 
                                     (gap_word_t)1);

    combine_operation_with_block(nblock_right,
        BM_IS_GAP(block),
        block,
        (bm::word_t*) tmp_gap_blk,
        1, BM_OR);
}

//---------------------------------------------------------------------

template<class Alloc>
void bvector<Alloc>::clear_range_no_check(bm::id_t left,
                                          bm::id_t right)
{
    unsigned nb;

    // calculate logical number of start and destination blocks
    unsigned nblock_left = unsigned(left >> bm::set_block_shift);
    unsigned nblock_right = unsigned(right >> bm::set_block_shift);

    unsigned nbit_right = unsigned(right & bm::set_block_mask);
    unsigned r =
        (nblock_left == nblock_right) ? nbit_right : (bm::bits_in_block - 1);

    bm::gap_word_t tmp_gap_blk[5];
    tmp_gap_blk[0] = 0; // just to silence GCC warning on uninit var
  
    // Set bits in the starting block
    bm::word_t* block;// = blockman_.get_block(nblock_left);
    
    unsigned nbit_left = unsigned(left  & bm::set_block_mask);
    if ((nbit_left == 0) && (r == bm::bits_in_block - 1)) // full block
    {
        nb = nblock_left;
    }
    else
    {
        bm::gap_init_range_block<gap_word_t>(tmp_gap_blk,
            (gap_word_t)nbit_left,
            (gap_word_t)r,
            (gap_word_t)0);
        block = blockman_.get_block(nblock_left);
        combine_operation_with_block(nblock_left,
            BM_IS_GAP(block),
            block,
            (bm::word_t*) tmp_gap_blk,
            1,
            BM_AND);
        
        if (nblock_left == nblock_right)  // in one block
            return;
        nb = nblock_left + 1;
    }

    // Clear all full blocks between left and right

    unsigned nb_to = nblock_right + (nbit_right == (bm::bits_in_block - 1));
    for (; nb < nb_to; ++nb)
    {
        int no_more_blocks;
        block = blockman_.get_block(nb, &no_more_blocks);
        if (no_more_blocks)
        {
            BM_ASSERT(block == 0);
            return;
        }
        if (!block)  // nothing to do
            continue;
        blockman_.zero_block(nb);
    } // for

    if (nb_to > nblock_right)
        return;

    block = blockman_.get_block(nblock_right);
    gap_init_range_block<gap_word_t>(tmp_gap_blk,
        (gap_word_t)0,
        (gap_word_t)nbit_right,
        (gap_word_t)0);
    
    combine_operation_with_block(nblock_right,
        BM_IS_GAP(block),
        block,
        (bm::word_t*) tmp_gap_blk,
        1,
        BM_AND);
}


//---------------------------------------------------------------------

template<class Alloc>
void bvector<Alloc>::copy_range(const bvector<Alloc>& bvect,
                                bm::id_t left,
                                bm::id_t right)
{
    if (!bvect.blockman_.is_init())
    {
        clear();
        return;
    }
    
    if (blockman_.is_init())
    {
        blockman_.deinit_tree();
    }
    if (left > right)
        bm::xor_swap(left, right);

    copy_range_no_check(bvect, left, right);
}


//---------------------------------------------------------------------

template<class Alloc>
void bvector<Alloc>::copy_range_no_check(const bvector<Alloc>& bvect,
                                         bm::id_t left,
                                         bm::id_t right)
{
    BM_ASSERT(left <= right);
    BM_ASSERT_THROW(right < bm::id_max, BM_ERR_RANGE);
    
    // copy all block(s) belonging to our range
    unsigned nblock_left  = unsigned(left  >>  bm::set_block_shift);
    unsigned nblock_right = unsigned(right >>  bm::set_block_shift);
    
    blockman_.copy(bvect.blockman_, nblock_left, nblock_right);
    
    // clear the flanks
    //
    if (left)
    {
        bm::id_t from =
            (left + bm::gap_max_bits >= left) ? 0u : left - bm::gap_max_bits;
        clear_range_no_check(from, left-1);
    }
    if (right+1 < bm::id_max)
    {
        clear_range_no_check(right+1, bm::id_max-1);
    }
}


} // namespace

#include "bmundef.h"

#ifdef _MSC_VER
#pragma warning( pop )
#endif


#endif
