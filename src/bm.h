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

/**
bit-block array wrapped into union for correct interpretation of
32-bit vs 64-bit access vs SIMD
@internal
*/
namespace bm
{
    struct bit_block_t
    {
        union bunion_t
        {
            bm::word_t BM_VECT_ALIGN w32[bm::set_block_size] BM_VECT_ALIGN_ATTR;
            bm::id64_t BM_VECT_ALIGN w64[bm::set_block_size / 2] BM_VECT_ALIGN_ATTR;
#ifdef BMAVX2OPT
            __m256i  BM_VECT_ALIGN w256[bm::set_block_size / 8] BM_VECT_ALIGN_ATTR;
#endif
#if defined(BMSSE2OPT) || defined(BMSSE42OPT)
            __m128i  BM_VECT_ALIGN w128[bm::set_block_size / 4] BM_VECT_ALIGN_ATTR;
#endif
        } b_;

        operator bm::word_t*() { return &(b_.w32[0]); }
        operator const bm::word_t*() const { return &(b_.w32[0]); }
        explicit operator bm::id64_t*() { return &b_.w64[0]; }
        explicit operator const bm::id64_t*() const { return &b_.w64[0]; }
#ifdef BMAVX2OPT
        explicit operator __m256i*() { return &b_.w256[0]; }
        explicit operator const __m256i*() const { return &b_.w256[0]; }
#endif
#if defined(BMSSE2OPT) || defined(BMSSE42OPT)
        explicit operator __m128i*() { return &b_.w128[0]; }
        explicit operator const __m128i*() const { return &b_.w128[0]; }
#endif

        const bm::word_t* begin() const { return (b_.w32 + 0); }
        bm::word_t* begin() { return (b_.w32 + 0); }
        const bm::word_t* end() const { return (b_.w32 + bm::set_block_size); }
        bm::word_t* end() { return (b_.w32 + bm::set_block_size); }
    };
}
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
    */
    typedef int (*bit_visitor_callback_type)(void* handle_ptr, bm::id_t bit_idx);
}


namespace bm
{


#ifdef BMCOUNTOPT

# define BMCOUNT_INC ++count_;
# define BMCOUNT_DEC --count_;
# define BMCOUNT_VALID(x) count_is_valid_ = x;
# define BMCOUNT_SET(x) count_ = x; count_is_valid_ = true;
# define BMCOUNT_ADJ(x) if (x) ++count_; else --count_;

#else

# define BMCOUNT_INC
# define BMCOUNT_DEC
# define BMCOUNT_VALID(x)
# define BMCOUNT_SET(x)
# define BMCOUNT_ADJ(x)

#endif





/** @defgroup bmagic BitMagic Library
    BitMagic C++ Library
    For more information please visit:
    http://bmagic.sourceforge.net
    https://github.com/tlk00/BitMagic
 
*/


/** @defgroup bvector bvector<>
    The Main bvector<> Group
    bvector<> template: front end of the BitMagic library.
    @ingroup bmagic
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

    typedef Alloc  allocator_type;
    typedef blocks_manager<Alloc>      blocks_manager_type;
    /** Type used to count bits in the bit vector */
    typedef bm::id_t                   size_type; 

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
        @ingroup bvector
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

        @ingroup bvector
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
                    bvect_->resize(n + 1);
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
        @ingroup bvector
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

        #ifdef BMCOUNTOPT
            if (this->bv_->count_is_valid_ && 
                this->bv_->count_ == 0)
            {
                this->invalidate();
                return;
            }
        #endif

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


            // next bit not present in the current block
            // keep looking in the next blocks.
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
                        {
                            return *this;
                        }
                    }
                    else
                    {
                        if (this->block_ == FULL_BLOCK_FAKE_ADDR)
                            this->block_ = FULL_BLOCK_REAL_ADDR;
                        if (search_in_bitblock())
                        {
                            return *this;
                        }
                    }

                } // for j

            } // for i

            this->invalidate();
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
                    bdescr->gap_.gap_len = gptr[gpos] - bm::gap_word_t(nbit - 1);
                }
                else
                {
                    bm::gap_word_t interval = gptr[gpos] - gptr[gpos - 1];
                    bm::gap_word_t interval2 = (bm::gap_word_t)(nbit - gptr[gpos - 1]);
                    bdescr->gap_.gap_len = interval - interval2 + 1;
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
        
        
        bool decode_bit_group(block_descr_type* bdescr)
        {
            const word_t* block_end = this->block_ + bm::set_block_size;
            
            for (; bdescr->bit_.ptr < block_end;)
            {
                bdescr->bit_.cnt = bm::bitscan_wave(bdescr->bit_.ptr, bdescr->bit_.bits);
                if (bdescr->bit_.cnt) // found
                {
                    bdescr->bit_.idx = 0;
                    bdescr->bit_.pos = this->position_;
                    this->position_ += bdescr->bit_.bits[0];
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

    };
    
    /*!
        @brief Constant iterator designed to enumerate "ON" bits
        counted_enumerator keeps bitcount, ie number of ON bits starting
        from the position 0 in the bit string up to the currently enumerated bit
        
        When increment operator called current position is increased by 1.
        
        @ingroup bvector
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


    friend class iterator_base;
    friend class enumerator;

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

#ifdef BMCOUNTOPT
    bvector(strategy          strat      = BM_BIT,
            const gap_word_t* glevel_len = bm::gap_len_table<true>::_len,
            size_type         bv_size    = bm::id_max,
            const Alloc&      alloc      = Alloc()) 
    : count_(0),
      count_is_valid_(true),
      blockman_(glevel_len, bv_size, alloc),
      new_blocks_strat_(strat),
      size_(bv_size)
    {}

    bvector(size_type         bv_size,
            bm::strategy      strat      = BM_BIT,
            const gap_word_t* glevel_len = bm::gap_len_table<true>::_len,
            const Alloc&      alloc = Alloc()) 
    : count_(0),
      count_is_valid_(true),
      blockman_(glevel_len, bv_size, alloc),
      new_blocks_strat_(strat),
      size_(bv_size)
    {}


    bvector(const bm::bvector<Alloc>& bvect)
     : count_(bvect.count_),
       count_is_valid_(bvect.count_is_valid_),
       blockman_(bvect.blockman_),
       new_blocks_strat_(bvect.new_blocks_strat_),
       size_(bvect.size_)
    {}

#else
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
    
    ~bvector() BMNOEXEPT
    {}

#endif

    /*!
        \brief Explicit post-construction initialization
    */
    void init()
    {
        if (!blockman_.is_init())
            blockman_.init_tree();
    }

    /*! 
        \brief Copy assignment operator
    */
    bvector& operator=(const bvector<Alloc>& bvect)
    {
        if (this != &bvect)
        {
            clear(true); // memory free cleaning
            resize(bvect.size());
            bit_or(bvect);
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
            
#ifdef BMCOUNTOPT
        count_ = bvect.count_;
        count_is_valid_ = bvect.count_is_valid_;
#endif
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
#ifdef BMCOUNTOPT
        count_ = (bm::id_t)il.size();
        count_is_valid_ = true;
#endif
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
    void move_from(bvector<Alloc>& bvect) BMNOEXEPT
    {
        if (this != &bvect)
        {
            blockman_.move_from(bvect.blockman_);
            size_ = bvect.size_;
            new_blocks_strat_ = bvect.new_blocks_strat_;
            
    #ifdef BMCOUNTOPT
            count_ = bvect.count_;
            count_is_valid_ = bvect.count_is_valid_;
    #endif
        }
    }
    
    /*! \brief Exchanges content of bv and this bvector.
    */
    void swap(bvector<Alloc>& bvect) BMNOEXEPT;
    
    //@}


    reference operator[](bm::id_t n)
    {
        BM_ASSERT(n < size_);
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
        BM_ASSERT(n < size_);
        BM_ASSERT_THROW(n < size_, BM_ERR_RANGE);

        if (!blockman_.is_init())
            blockman_.init_tree();
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
            blockman_.init_tree();
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
        BM_ASSERT(n < size_);
        BM_ASSERT_THROW(n < size_, BM_ERR_RANGE);

        if (val == condition) return false;
        if (!blockman_.is_init())
            blockman_.init_tree();
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
        BMCOUNT_VALID(false)
        set_range(0, size_ - 1, true);
        return *this;
    }
    
    /*!
        \brief Set bit without checking preconditions (size, etc)
     
        Fast set bit method, without safety net.
        Make sure you call bvector<>::init() before setting bits with this
        function.
     
        Side effect: invalidates bit-count optimization (BMCOUNTOPT)
     
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
       \brief Clears bit n.
       \param n - bit's index to be cleaned.
       \return true if bit was cleared
    */
    bool clear_bit(bm::id_t n)
    {
        return set_bit(n, false);
    }
    
    /*!
       \brief Clears every bit in the bitvector.

       \param free_mem if "true" (default) bvector frees the memory,
       otherwise sets blocks to 0.
    */
    void clear(bool free_mem = false)
    {
        blockman_.set_all_zero(free_mem);
        BMCOUNT_SET(0);
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
        set(n, !get_bit(n));
        //this->inc(n);
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

    /*! @name Population counting methods
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
        for_each_nzblock(blk_root, blockman_.effective_top_block_size(), 
                         func);
        return func.last_block();
    }

    /*!
       \brief Returns count of 1 bits in the given diapason.
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
       \param right - index of last bit
       \param blocks_cnt - block count structure to accelerate search
                           should be prepared using running_count_blocks
       \return population count in the diapason
       \sa running_count_blocks
       \sa count_to_test
    */
    bm::id_t count_to(bm::id_t right, const blocks_count&  blocks_cnt) const;

    /*!
        \brief Returns count of 1 bits (population) in [0..right] range if test(right) == true.
        Faster than test() plus count_to()
        \param right - index of last bit
        \param blocks_cnt - block count structure to accelerate search
               should be prepared using running_count_blocks

        \return population count in the diapason or 0 if right bit test failed

        \sa running_count_blocks
        \sa count_to
    */
    bm::id_t count_to_test(bm::id_t right, const blocks_count&  blocks_cnt) const;

    /*! Recalculate bitcount
        this function only make sense when BMCOUNTOPT is defined
        and bvector<> keeps its bitcount. Otherwise, equivalent of cout().
    */
    bm::id_t recalc_count()
    {
        BMCOUNT_VALID(false)
        return count();
    }
    
    /*!
        Disables count cache. Next call to count() or recalc_count()
        restores count caching.
        
        @note Works only if BMCOUNTOPT enabled(defined). 
        Othewise does nothing.
    */
    void forget_count()
    {
        BMCOUNT_VALID(false)    
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
    #ifdef BMCOUNTOPT
        if (count_is_valid_)
            return count_ != 0;
    #endif
        
        word_t*** blk_root = blockman_.top_blocks_root();
        if (!blk_root) 
            return false;
        typename blocks_manager_type::block_any_func func(blockman_);
        return for_each_nzblock_if(blk_root, 
                                   blockman_.effective_top_block_size(),
                                   func);
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


    /*! @name Scan and find bit and indexes */
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
        BMCOUNT_VALID(false);
        combine_operation(vect, BM_OR);
        return *this;
    }

    /*!
       \brief Logical AND operation.
       \param vect - Argument vector.
    */
    bm::bvector<Alloc>& bit_and(const bm::bvector<Alloc>& vect)
    {
        BMCOUNT_VALID(false);
        combine_operation(vect, BM_AND);
        return *this;
    }

    /*!
       \brief Logical XOR operation.
       \param vect - Argument vector.
    */
    bm::bvector<Alloc>& bit_xor(const bm::bvector<Alloc>& vect)
    {
        BMCOUNT_VALID(false);
        combine_operation(vect, BM_XOR);
        return *this;
    }

    /*!
       \brief Logical SUB operation.
       \param vect - Argument vector.
    */
    bm::bvector<Alloc>& bit_sub(const bm::bvector<Alloc>& vect)
    {
        BMCOUNT_VALID(false);
        combine_operation(vect, BM_SUB);
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

    // @}

    // --------------------------------------------------------------------
    /*! @name Iterator-traversal methods  */
    //@{

    /**
       \brief Returns enumerator pointing on the first non-zero bit.
    */
    enumerator first() const
    {
        return get_enumerator(0);
    }

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
private:
#if 0
    void combine_count_operation_with_block(unsigned nb,
                                            const bm::word_t* arg_blk,
                                            bool arg_gap,
                                            bm::operation opcode)
    {
        const bm::word_t* blk = get_block(nb);
        bool gap = BM_IS_GAP(blk);
        combine_count_operation_with_block(nb, gap, blk, arg_blk, arg_gap, opcode);
    }
#endif


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
       \brief Set range without validity checking
    */
    void set_range_no_check(bm::id_t left,
                            bm::id_t right,
                            bool     value);
public:



private:

// This block defines two additional hidden variables used for bitcount
// optimization, in rare cases can make bitvector thread unsafe.
#ifdef BMCOUNTOPT
    mutable id_t      count_;            //!< number of 1 bits in the vector
    mutable bool      count_is_valid_;   //!< actualization flag
#endif

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

    BM_ASSERT(left < size_);
    BM_ASSERT(right < size_);

    BMCOUNT_VALID(false)
    BM_SET_MMX_GUARD

    set_range_no_check(left, right, value);

    return *this;
}

// -----------------------------------------------------------------------

template<typename Alloc> 
bm::id_t bvector<Alloc>::count() const
{
#ifdef BMCOUNTOPT
    if (count_is_valid_) return count_;
#endif
    if (!blockman_.is_init())
        return 0;
    
    word_t*** blk_root = blockman_.top_blocks_root();
    if (!blk_root) 
    {
        BMCOUNT_SET(0);
        return 0;
    }    
    typename blocks_manager_type::block_count_func func(blockman_);
    for_each_nzblock2(blk_root, blockman_.effective_top_block_size(), 
                      func);

    BMCOUNT_SET(func.count());
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
    BMCOUNT_VALID(false)
    BM_SET_MMX_GUARD
    
    if (!blockman_.is_init())
        blockman_.init_tree();

    bm::word_t*** blk_root = blockman_.top_blocks_root();
    typename blocks_manager_type::block_invert_func func(blockman_);    
    for_each_block(blk_root, blockman_.top_block_size(), func);
    if (size_ == bm::id_max) 
    {
        set_bit_no_check(bm::id_max, false);
    } 
    else
    {
        set_range_no_check(size_, bm::id_max, false);
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

    for_each_nzblock(blk_root, blockman_.effective_top_block_size(),
                     opt_func);

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
    
    unsigned top_blocks = blockman_.effective_top_block_size();
    unsigned bvect_top_blocks = bv.blockman_.effective_top_block_size();

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
                    if (!gap_is_all_zero(BMGAP_PTR(pblk), bm::gap_max_bits))
                    {
                        return res;
                    }
                }
                else
                {
                    bm::wordop_t* blk1 = (wordop_t*)pblk;
                    bm::wordop_t* blk2 = 
                        (wordop_t*)(pblk + bm::set_block_size);
                    if (!bit_is_all_zero(blk1, blk2))
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
#ifdef BMCOUNTOPT
        bm::xor_swap(count_, bvect.count_);
        bm::xor_swap(count_is_valid_, bvect.count_is_valid_);
#endif
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

    unsigned top_size = blockman_.effective_top_block_size();
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
                    unsigned capacity = bm::gap_capacity(gap_blk, blockman_.glen());
                    unsigned length = gap_length(gap_blk);
                    
                    st->add_gap_block(capacity, length);
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
    
    BMCOUNT_VALID(false)
    
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
                BMCOUNT_INC;
                return true;
            }
        }
        else
        {
            if ((*word) & mask)
            {
                *word &= ~mask; // clear bit
                BMCOUNT_DEC;
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
        BMCOUNT_ADJ(val)

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
    if (!blockman_.is_init())
        blockman_.init_tree();

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
        bool curr_val = (bm::gap_test_unr(gap_blk, nbit) != 0);
        is_set = this->gap_block_set(gap_blk, !curr_val, nblock, nbit); // flip
    }
    else // bit block
    {
        unsigned nword  = unsigned(nbit >> bm::set_word_shift);
        nbit &= bm::set_word_mask;

        bm::word_t* word = blk + nword;
        bm::word_t  mask = (((bm::word_t)1) << nbit);
        is_set = ((*word) & mask);
        
        if (is_set) // need to clear the bit (because we do flip here)
        {
            *word &= ~mask;
            BMCOUNT_DEC;
        }
        else // set the bit
        {
            *word |= mask;
            BMCOUNT_INC;
        }
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
                BMCOUNT_INC;
            }
            else               // clear bit
            {
                *word &= ~mask;
                BMCOUNT_DEC;
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
                BMCOUNT_INC;
            }
            else               // clear bit
            {
                *word &= ~mask;
                BMCOUNT_DEC;
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
    
    unsigned top_blocks = blockman_.effective_top_block_size();
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
    
    unsigned top_blocks = blockman_.effective_top_block_size();
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
bool bvector<Alloc>::find_range(bm::id_t& first, bm::id_t& last) const
{
    bool found = find(first);
    if (found)
    {
        found = find_reverse(last);
        BM_ASSERT(found);
    }
    return found;
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
                        BMCOUNT_DEC
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
                        BMCOUNT_DEC

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

    if (size_ == bv.size_)
    {
        BM_ASSERT(top_blocks >= arg_top_blocks);
    }
    else
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

    BM_SET_MMX_GUARD

    // calculate effective top size to avoid overscan
    top_blocks = blockman_.effective_top_block_size();
    if (top_blocks < bv.blockman_.effective_top_block_size())
    {
        if (opcode != BM_AND)
        {
            top_blocks = bv.blockman_.effective_top_block_size();
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
    int      level;
    unsigned threshold;


    if (opcode == BM_OR || opcode == BM_XOR)
    {        
        if (!blk && arg_gap) 
        {
            res = BMGAP_PTR(arg_blk);
            res_len = bm::gap_length(res);
            level = -1;
            threshold = 0;
            goto assign_gap_result;
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
                ++res_len;// = bm::gap_length(res);

                BM_ASSERT(!(res == tmp_buf && res_len == 0));

                // if as a result of the operation gap block turned to zero
                // we can now replace it with NULL
                if (gap_is_all_zero(res, bm::gap_max_bits))
                {
                    blockman_.zero_block(nb);
                    return;
                }

                // mutation check

                level = gap_level(BMGAP_PTR(blk));
                threshold = blockman_.glen(level)-4;

            assign_gap_result:
                int new_level = gap_calc_level(res_len, blockman_.glen());
                if (new_level == -1)
                {
                    blockman_.convert_gap2bitset(nb, res);
                    return;
                }

                if (res_len > threshold)
                {
                    gap_word_t* new_blk = 
                        blockman_.allocate_gap_block(new_level, res);
                    set_gap_level(new_blk, new_level);

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

                set_gap_level(tmp_buf, level);
                ::memcpy(BMGAP_PTR(blk), tmp_buf, res_len * sizeof(gap_word_t));
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
                if (opcode == BM_AND)
                {
                    unsigned gap_cnt = gap_bit_count_unr(gap_blk);
                    if (gap_cnt < 128)
                    {
                        
                        gap_word_t arr_len = 
                            gap_convert_to_arr(tmp_buf, gap_blk, 
                                               bm::gap_equiv_len-10);
                        BM_ASSERT(gap_cnt == arr_len);
                        blockman_.zero_block(nb);
                        unsigned arr_i = 0;
                        int block_type;
                        blk =
                            blockman_.check_allocate_block(nb,
                                                           true,
                                                           BM_GAP,
                                                           &block_type,
                                                           false //no null return
                                                           );
                        BM_ASSERT(block_type==1); // GAP
                        gap_blk = BMGAP_PTR(blk);
                        threshold = bm::gap_limit(gap_blk, blockman_.glen());
                        for (; arr_i < arr_len; ++arr_i)
                        {
                            gap_word_t bit_idx = tmp_buf[arr_i];
                            if (bm::test_bit(arg_blk, bit_idx))
                            {
                                unsigned is_set;
                                unsigned new_block_len =
                                    gap_set_value(true, gap_blk, bit_idx, &is_set);
                                BM_ASSERT(is_set);
                                if (new_block_len > threshold) 
                                {
                                    gap_blk = 
                                        blockman_.extend_gap_block(nb, gap_blk);
                                    if (gap_blk == 0) // mutated into bit-block
                                    {
                                        blk = blockman_.check_allocate_block(
                                                         nb,
                                                         true,
                                                         this->get_new_blocks_strat(),
                                                         &block_type,
                                                         false // no null return
                                                         );  
                                        BM_ASSERT(block_type == 0); // BIT
                                        // target block degraded into plain bit-block
                                        for (++arr_i; arr_i < arr_len; ++arr_i)
                                        {
                                            bit_idx = tmp_buf[arr_i];
                                            if (bm::test_bit(arg_blk, bit_idx))
                                            {
                                                or_bit_block(blk, bit_idx, 1);
                                            }
                                        } // for arr_i
                                        return;
                                    } // if gap mutated
                                }
                            } // for arr_i
                        }

                        return;
                    }                    
                } // BM_AND

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
#ifdef BMVECTOPT
                VECT_ANDNOT_ARR_2_MASK(ret, 
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
                    dst_ptr[0] = bm::all_bits_mask & ~wrd_ptr[0];
                    dst_ptr[1] = bm::all_bits_mask & ~wrd_ptr[1];
                    dst_ptr[2] = bm::all_bits_mask & ~wrd_ptr[2];
                    dst_ptr[3] = bm::all_bits_mask & ~wrd_ptr[3];

                    dst_ptr+=4;
                    wrd_ptr+=4;

                } while (wrd_ptr < wrd_end);
#endif
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
void bvector<Alloc>::set_range_no_check(bm::id_t left,
                                        bm::id_t right,
                                        bool     value)
{
    // calculate logical number of start and destination blocks
    unsigned nblock_left  = unsigned(left  >>  bm::set_block_shift);
    unsigned nblock_right = unsigned(right >>  bm::set_block_shift);

    bm::word_t* block = blockman_.get_block(nblock_left);
    bool left_gap = BM_IS_GAP(block);

    unsigned nbit_left  = unsigned(left  & bm::set_block_mask); 
    unsigned nbit_right = unsigned(right & bm::set_block_mask); 

    unsigned r = 
        (nblock_left == nblock_right) ? nbit_right :(bm::bits_in_block-1);

        bm::gap_word_t tmp_gap_blk[5] = {0,};

    // Set bits in the starting block

    unsigned nb;
    if ((nbit_left == 0) && (r == bm::bits_in_block - 1)) // full block
    {
        nb = nblock_left;
    }
    else
    {
        gap_init_range_block<gap_word_t>(tmp_gap_blk,
                                         (gap_word_t)nbit_left, 
                                         (gap_word_t)r, 
                                         (gap_word_t)value, 
                                         bm::bits_in_block);

        combine_operation_with_block(nblock_left, 
                                    left_gap, 
                                    block,
                                    (bm::word_t*) tmp_gap_blk,
                                    1,
                                    value ? BM_OR : BM_AND);

        if (nblock_left == nblock_right)  // in one block
            return;
        nb = nblock_left+1;
    }

    // Set (or clear) all full blocks between left and right
    
    unsigned nb_to = nblock_right + (nbit_right ==(bm::bits_in_block-1));
            
    if (value)
    {
        for (; nb < nb_to; ++nb)
        {
            block = blockman_.get_block(nb);
            if (IS_FULL_BLOCK(block)) 
                continue;

            bool is_gap = BM_IS_GAP(block);

            blockman_.set_block(nb, FULL_BLOCK_FAKE_ADDR);
            blockman_.set_block_bit(nb);
            
            if (is_gap)
            {
                blockman_.get_allocator().free_gap_block(BMGAP_PTR(block), 
                                                            blockman_.glen());
            }
            else
            {
                if (IS_VALID_ADDR(block))
                    blockman_.get_allocator().free_bit_block(block);
            }
            
        } // for
    }
    else // value == 0
    {
        for (; nb < nb_to; ++nb)
        {
            block = blockman_.get_block(nb);
            if (block == 0)  // nothing to do
                continue;
            bool is_gap = BM_IS_GAP(block);
            blockman_.set_block(nb, 0, false /*bit*/);
            //blockman_.set_block_bit(nb);

            if (is_gap) 
            {
                blockman_.get_allocator().free_gap_block(BMGAP_PTR(block),
                                                         blockman_.glen());
            }
            else
            {
                if (IS_VALID_ADDR(block))
                    blockman_.get_allocator().free_bit_block(block);
            }

        } // for
    } // if value else 

    if (nb_to > nblock_right)
        return;

    block = blockman_.get_block(nblock_right);
    bool right_gap = BM_IS_GAP(block);

    gap_init_range_block<gap_word_t>(tmp_gap_blk, 
                                     (gap_word_t)0, 
                                     (gap_word_t)nbit_right, 
                                     (gap_word_t)value, 
                                     bm::bits_in_block);

    combine_operation_with_block(nblock_right, 
                                    right_gap, 
                                    block,
                                    (bm::word_t*) tmp_gap_blk,
                                    1,
                                    value ? BM_OR : BM_AND);

}

//---------------------------------------------------------------------


} // namespace

#include "bmundef.h"

#ifdef _MSC_VER
#pragma warning( pop )
#endif


#endif
