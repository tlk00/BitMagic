
#ifndef BM_SPARSE_VEC_FLOAT_INCLUDED
#define BM_SPARSE_VEC_FLOAT_INCLUDED
/*
Copyright(c) 2026 Anatoliy Kuznetsov(tolikkuznetsov66 at gmail.com)

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

/*! \file sparse_vector_float.h
    \brief Wrapper class for sparse_vector<>, sparse_vector_float that allows storage of floats
*/

#include <memory.h>

#include "bm.h"
#include "bmsparsevec.h"

namespace bm
{

/** \defgroup svector Sparse and compressed vectors
    Sparse vector for float types by splitting floats into sign, exponent and mantissa
    and storing those in sparse vectors
 
    @ingroup bmagic
 */

/*!
   \brief Wrapper for sparse vector that supports floats
 
   @ingroup sv
*/
template<class BV>
class sparse_vector_float
{
    friend class sparse_vector_float_serialized;
public:
    //
    typedef float                                   value_type;
    typedef BV                                      bvector_type;
    typedef typename bvector_type::size_type        size_type;
    typedef bm::sparse_vector<unsigned int, BV>     sparse_vector_u32;
    typedef typename BV::allocator_type             allocator_type;
    typedef typename bvector_type::allocation_policy allocation_policy_type;

    struct statistics : public bv_statistics
    {};
    
    /*
    const_iterator for traversing the sparse_vector_float
    */
    class const_iterator
    {
    public:
    friend class sparse_vector_float;

#ifndef BM_NO_STL
        typedef std::input_iterator_tag  iterator_category;
#endif
        typedef sparse_vector_float                        sparse_vector_type;
        typedef sparse_vector_type*                        sparse_vector_type_ptr;
        typedef typename sparse_vector_type::value_type    value_type;
        typedef typename sparse_vector_type::size_type     size_type;
        typedef typename sparse_vector_type::bvector_type  bvector_type;
        typedef typename bvector_type::allocator_type      allocator_type;
        typedef typename bvector_type::allocator_type::allocator_pool_type allocator_pool_type;
        typedef bm::byte_buffer<allocator_type>            buffer_type;

        typedef unsigned                    difference_type;
        typedef unsigned*                   pointer;
        typedef value_type&                 reference;

    public:
        const_iterator() BMNOEXCEPT;
        const_iterator(const sparse_vector_type* sv) BMNOEXCEPT;
        const_iterator(const sparse_vector_type* sv, size_type pos) BMNOEXCEPT;
        const_iterator(const const_iterator& it) BMNOEXCEPT;

        
        bool operator==(const const_iterator& it) const BMNOEXCEPT
                                { return (pos_ == it.pos_) && (sv_ == it.sv_); }
        bool operator!=(const const_iterator& it) const BMNOEXCEPT
                                { return ! operator==(it); }
        bool operator < (const const_iterator& it) const BMNOEXCEPT
                                { return pos_ < it.pos_; }
        bool operator <= (const const_iterator& it) const BMNOEXCEPT
                                { return pos_ <= it.pos_; }
        bool operator > (const const_iterator& it) const BMNOEXCEPT
                                { return pos_ > it.pos_; }
        bool operator >= (const const_iterator& it) const BMNOEXCEPT
                                { return pos_ >= it.pos_; }

        /// \brief Get current position (value)
        value_type operator*() const  { return this->value(); }
        
        
        /// \brief Advance to the next available value
        const_iterator& operator++() BMNOEXCEPT { this->advance(); return *this; }
        /// \brief Advance to the next available value
        ///
        const_iterator operator++(int)
            { const_iterator tmp(*this);this->advance(); return tmp; }


        /// \brief Get current position (value)
        value_type value() const;
        
        /// \brief Get NULL status
        bool is_null() const BMNOEXCEPT;
        
        /// Returns true if iterator is at a valid position
        bool valid() const BMNOEXCEPT { return pos_ != bm::id_max; }
        
        /// Invalidate current iterator
        void invalidate() BMNOEXCEPT { pos_ = bm::id_max; }
        
        /// Current position (index) in the vector
        size_type pos() const BMNOEXCEPT{ return pos_; }
        
        /// re-position to a specified position
        void go_to(size_type pos) BMNOEXCEPT;
        
        /// advance iterator forward by one
        /// @return true if it is still valid
        bool advance() BMNOEXCEPT;
        
    private:
        const sparse_vector_type*         sv_;      ///!< ptr to parent
        size_type                         pos_;     ///!< Position
        sparse_vector_u32::const_iterator            exp_it_;
        sparse_vector_u32::const_iterator            mant_it_;
    };

    const_iterator begin() const;
    const_iterator end() const { return const_iterator(this, bm::id_max); };


    /*!
        \brief Sparse vector constructor
     
        \param null_able - defines if vector supports NULL values flag
            by default it is OFF, use bm::use_null to enable it
        \param ap - allocation strategy for underlying bit-vectors
        Default allocation policy uses BM_BIT setting (fastest access)
        \param bv_max_size - maximum possible size of underlying bit-vectors
        Please note, this is NOT size of svector itself, it is dynamic upper limit
        which should be used very carefully if we surely know the ultimate size
        \param alloc - allocator for bit-vectors
        
        \sa bvector<>
        \sa bm::bvector<>::allocation_policy
        \sa bm::startegy
    */
    sparse_vector_float(bm::null_support null_able = bm::no_null,
                        allocation_policy_type ap = allocation_policy_type(),
                        size_type bv_max_size = bm::id_max,
                        const allocator_type&   alloc  = allocator_type());

    sparse_vector_float(const sparse_vector_float& svf);
    ~sparse_vector_float();

    sparse_vector_float& operator=(const sparse_vector_float& svf);
    bool operator==(const sparse_vector_float& svf) const;
    bool operator!=(const sparse_vector_float& svf) const;


    void swap(sparse_vector_float& svf);
    

    /*! \brief return size of the vector
        \return size of sparse vector
    */
    size_type size() const BMNOEXCEPT { return this->mantissas_.size(); };

    /*! \brief return true if vector is empty
        \return true if empty
    */
    bool empty() const BMNOEXCEPT { return (size() == 0); }

    /*!
        \brief get specified element without bounds checking
        \param idx - element index
        \return value of the element
    */
    value_type get(size_type idx) const BMNOEXCEPT;

    /*!
        \brief push value back into vector
        \param v   - element value
    */
    void push_back(value_type v);

    /*!
        \brief set specified element with bounds checking and automatic resize
        \param idx - element index
        \param v   - element value
    */
    void set(size_type idx, value_type v);

    /*! \brief resize to zero, free memory */
    void clear() BMNOEXCEPT;

    /*!
        \brief clear range (assign bit 0 for all planes)
        \param left  - interval start
        \param right - interval end (closed interval)
        \param set_null - set cleared values to unassigned (NULL)
    */
    sparse_vector_float<BV>& clear_range(size_type left,
                                        size_type right,
                                        bool set_null = false);

     /*!
        \brief check if another sparse vector float has the same content and size
     
        \param sv        - sparse vector for comparison
        \param null_able - flag to consider NULL vector in comparison (default)
                           or compare only value content planes
     
        \return true, if it is the same
    */
    bool equal(const sparse_vector_float<BV>& sv,
               bm::null_support null_able = bm::use_null) const BMNOEXCEPT;

    /**
        \brief Compare vector element with argument
     
        \param idx - vactor element index
        \param val - argument to compare with
        \param epsilon - amount of precision
     
        \return 0 - equal, -1 - vect[i] < val, 1 otherwise
    */
    int compare(size_type idx, const value_type val, float epsilon = std::numeric_limits<float>::epsilon()) const BMNOEXCEPT;

    /*!
        \brief join all with another sparse vector using OR operation
            NOTE: if you join 2 floats together and one is not 0.0, ie -1.5f and -2.5f
            it is possible for the float to become NaN if the exponent becomes entirely 11111111
        \param sv - argument vector to join with
        \return slf reference
        @sa merge
    */
    sparse_vector_float<BV>& join(const sparse_vector_float<BV>& svf);

    /*!
        \brief merge with another sparse vector using OR operation
        Merge is different from join(), because it borrows data from the source
        vector, so it gets modified.
     
        \param sv - [in, out]argument vector to join with (vector mutates)
     
        \return slf reference
        @sa join
    */
    sparse_vector_float<BV>& merge(sparse_vector_float<BV>& svf);

    /**
        @brief copy range of values from another sparse vector
     
        Copy [left..right] values from the source vector,
        clear everything outside the range.
     
        \param sv - source vector
        \param left  - index from in losed diapason of [left..right]
        \param right - index to in losed diapason of [left..right]
        \param slice_null - "use_null" copy range for NULL vector or
                             do not copy it
    */
    void copy_range(const sparse_vector_float<BV>& svf,
                    size_type left, size_type right,
                    bm::null_support slice_null = bm::use_null);


    /*!
        \brief Import list of elements from a C-style array
        \param arr  - source array
        \param arr_size - source size
        \param offset - target index in the sparse vector
        \param set_not_null - import should register in not null vector
    */
    void import(const value_type* arr,
                size_type arr_size,
                size_type offset = 0,
                bool      set_not_null = true);

    void optimize(bm::word_t* temp_block = 0,
          typename bvector_type::optmode opt_mode = bvector_type::opt_compress);

    /*!
        @brief Calculates memory statistics.

        Function fills statistics structure containing information about how
        this vector uses memory and estimation of max. amount of memory
        bvector needs to serialize itself.

        @param st - pointer on statistics structure to be filled in.

        @sa statistics
    */
    void calc_stat(
        struct sparse_vector_float<BV>::statistics* st) const BMNOEXCEPT;

protected:
    enum buf_size_e{
        n_buf_size = 1024 * 8
    };

private:
    bm::bvector<>       signs_;
    sparse_vector_u32   exponents_;
    sparse_vector_u32   mantissas_;
    
};

//---------------------------------------------------------------------
//sparse_vector_float methods

//---------------------------------------------------------------------

template<class BV>
sparse_vector_float<BV>::sparse_vector_float(bm::null_support null_able,
                                             allocation_policy_type ap,
                                             size_type bv_max_size,
                                             const allocator_type&   alloc)
:signs_(null_able, ap, bv_max_size, alloc),
 exponents_(null_able, ap, bv_max_size, alloc),
 mantissas_(null_able, ap, bv_max_size, alloc)
{}

//---------------------------------------------------------------------

template<class BV>
sparse_vector_float<BV>::sparse_vector_float(const sparse_vector_float& svf)
:signs_(svf.signs_), exponents_(svf.exponents_), mantissas_(svf.mantissas_)
{}

//---------------------------------------------------------------------

template<class BV>
sparse_vector_float<BV>::~sparse_vector_float()
{}

//---------------------------------------------------------------------

template<class BV>
sparse_vector_float<BV>::const_iterator sparse_vector_float<BV>::begin() const
{
    typedef typename sparse_vector_float::const_iterator it_type;
    return it_type(this);
}

//---------------------------------------------------------------------

template<class BV>
sparse_vector_float<BV>& sparse_vector_float<BV>::operator=(const sparse_vector_float& svf)
{
    signs_ = svf.signs_;
    exponents_ = svf.exponents_;
    mantissas_ = svf.mantissas_;
    return *this;
}

//---------------------------------------------------------------------

template<class BV>
bool sparse_vector_float<BV>::operator==(const sparse_vector_float<BV>& svf) const
{
    if (mantissas_.size() == svf.mantissas_.size() && signs_ == svf.signs_)
    {
        for (size_t i = 0; i < mantissas_.size(); i++)
        {
            unsigned int thisSVF = exponents_.get(i);
            unsigned int otherSVF = svf.exponents_.get(i);

            if (thisSVF != otherSVF) return false;

            thisSVF = mantissas_.get(i);
            otherSVF = svf.mantissas_.get(i);

            if (thisSVF != otherSVF) return false;
        }

        return true;
    }
    return false;
}

//---------------------------------------------------------------------

template<class BV>
bool sparse_vector_float<BV>::operator!=(const sparse_vector_float<BV>& svf) const
{
    return !operator==(svf);
}

//---------------------------------------------------------------------

template<class BV>
void sparse_vector_float<BV>::swap(sparse_vector_float<BV> &svf)
{
    signs_.swap(svf.signs_);
    exponents_.swap(svf.exponents_);
    mantissas_.swap(svf.mantissas_);
}

//---------------------------------------------------------------------

template<class BV>
typename sparse_vector_float<BV>::value_type 
sparse_vector_float<BV>::get(typename sparse_vector_float<BV>::size_type idx) const BMNOEXCEPT
{
    unsigned int sign = signs_.test(idx) ? 1 : 0;
    unsigned int exponent = exponents_.get(idx);
    unsigned int mantissa = mantissas_.get(idx);
    
    unsigned int bits = (sign << 31) | (exponent << 23) | mantissa;

    float toReturn;
    memcpy(&toReturn, &bits, sizeof(float));

    return toReturn;
}

//---------------------------------------------------------------------

template<class BV>
void sparse_vector_float<BV>::push_back(value_type v)
{
    unsigned int bits;
    memcpy(&bits, &v, sizeof(float));

    unsigned int sign     = (bits >> 31) & 0x1;
    if (sign == 1) signs_.set(mantissas_.size());
    unsigned int exponent = (bits >> 23) & 0xFF;
    exponents_.push_back(exponent);
    unsigned int mantissa =  bits        & 0x7FFFFF;
    mantissas_.push_back(mantissa);
}

//---------------------------------------------------------------------

template<class BV>
void sparse_vector_float<BV>::set(size_type idx, value_type v)
{

    unsigned int bits;
    memcpy(&bits, &v, sizeof(float));

    unsigned int sign     = (bits >> 31) & 0x1;
    bool neg = false;
    if (sign == 1) neg = true;

    signs_.set(idx, neg);
    
    unsigned int exponent = (bits >> 23) & 0xFF;
    exponents_.set(idx, exponent);
    unsigned int mantissa =  bits        & 0x7FFFFF;
    mantissas_.set(idx, mantissa);
}

//---------------------------------------------------------------------

template<class BV>
void sparse_vector_float<BV>::clear() BMNOEXCEPT
{
    signs_.clear();
    exponents_.clear();
    mantissas_.clear();
}

//---------------------------------------------------------------------

template<class BV>
sparse_vector_float<BV>& sparse_vector_float<BV>::clear_range(size_type left,
                                    size_type right,
                                    bool set_null = false)
{
    signs_.clear_range(left, right);
    exponents_.clear_range(left, right, set_null);
    mantissas_.clear_range(left, right, set_null);

    return *this;
}

//---------------------------------------------------------------------

template<class BV>
bool sparse_vector_float<BV>::equal(const sparse_vector_float<BV>& sv,
                                    bm::null_support null_able = bm::use_null) const BMNOEXCEPT
{
    if (signs_.equal(sv.signs_) && exponents_.equal(sv.exponents_, null_able) && mantissas_.equal(sv.mantissas_, null_able)) return true;
    return false;
}

//---------------------------------------------------------------------

template<class BV>
int sparse_vector_float<BV>::compare(size_type idx, const value_type val, float epsilon) const BMNOEXCEPT
{
    value_type svf_value = get(idx);

    float diff = svf_value - val;

    if (std::fabs(diff) <= epsilon)     return 0;
    else if (diff > 0)                  return 1;
    else                                return -1;
}

//---------------------------------------------------------------------


template<class BV>
sparse_vector_float<BV>& sparse_vector_float<BV>::join(const sparse_vector_float<BV>& svf)
{
    signs_.join(svf.signs_);
    exponents_.join(svf.exponents_);
    mantissas_.join(svf.mantissas_);
    
    return *this;
}

//---------------------------------------------------------------------

template<class BV>
sparse_vector_float<BV>& sparse_vector_float<BV>::merge(sparse_vector_float<BV>& svf)
{
    signs_.merge(svf.signs_);
    exponents_.merge(svf.exponents_);
    mantissas_.merge(svf.mantissas_);
    
    return *this;
}


//---------------------------------------------------------------------

template<class BV>
void sparse_vector_float<BV>::copy_range(const sparse_vector_float<BV>& svf,
                                        size_type left, size_type right,
                                        bm::null_support slice_null)
{
    signs_.copy_range(svf.signs_, left, right, slice_null);
    exponents_.copy_range(svf.exponents_, left, right, slice_null);
    mantissas_.copy_range(svf.mantissas_, left, right, slice_null);
}


//---------------------------------------------------------------------

template<class BV>
void sparse_vector_float<BV>::import(const value_type* arr, 
                                size_type arr_size, 
                                size_type offset, 
                                bool set_not_null)
{
    if (offset > mantissas_.size())
    {
        for (size_type i = mantissas_.size(); i < offset; ++i)
        {
            push_back(0.0f);
        }
    }

    for (size_type i = 0; i < arr_size; ++i) {
        unsigned int bits;
        memcpy(&bits, &arr[i], sizeof(float));

        unsigned int sign     = (bits >> 31) & 0x1;
        unsigned int exponent = (bits >> 23) & 0xFF;
        unsigned int mantissa =  bits        & 0x7FFFFF;

        size_type idx = offset + i;

        bool neg = false;
        if (sign == 1)
            neg = true;

        signs_.set(idx, neg);
        exponents_.set(idx, exponent);
        mantissas_.set(idx, mantissa);
    }

    size_type new_size = offset + arr_size;
}

//---------------------------------------------------------------------

template<class BV>
void sparse_vector_float<BV>::optimize(bm::word_t* temp_block, typename bvector_type::optmode opt_mode)
{
    signs_.optimize(temp_block, opt_mode);
    exponents_.optimize(temp_block, opt_mode);
    mantissas_.optimize(temp_block, opt_mode);
}

//---------------------------------------------------------------------

template<class BV>
void sparse_vector_float<BV>::calc_stat(struct sparse_vector_float<BV>::statistics* st) const BMNOEXCEPT
{
    BM_ASSERT(st);

    bm::bvector<>::statistics signStat;
    signs_.calc_stat(signStat);

    sparse_vector_u32::statistics expStat, mantStat;
    exponents_.calc_stat(expStat);
    mantissas_.calc_stat(mantStat);

    st->bit_blocks        = signStat.bit_blocks   + expStat.bit_blocks   + mantStat.bit_blocks;
    st->gap_blocks        = signStat.gap_blocks   + expStat.gap_blocks   + mantStat.gap_blocks;
    st->ptr_sub_blocks    = signStat.ptr_sub_blocks + expStat.ptr_sub_blocks + mantStat.ptr_sub_blocks;
    st->memory_used       = signStat.memory_used  + expStat.memory_used  + mantStat.memory_used;
    st->max_serialize_mem = signStat.max_serialize_mem + expStat.max_serialize_mem + mantStat.max_serialize_mem;
    st->gap_cap_overhead  = signStat.gap_cap_overhead  + expStat.gap_cap_overhead  + mantStat.gap_cap_overhead;
    st->bv_count          = 1 + expStat.bv_count + mantStat.bv_count; // 1 for signs bvector

    // gap levels - take the max across all components
    for (int i = 0; i < bm::gap_levels; i++) {
        st->gap_levels[i] = std::max({signStat.gap_levels[i], 
                                      expStat.gap_levels[i], 
                                      mantStat.gap_levels[i]});
        st->gaps_by_level[i] = signStat.gaps_by_level[i] 
                              + expStat.gaps_by_level[i] 
                              + mantStat.gaps_by_level[i];
    }
}

//---------------------------------------------------------------------

//---------------------------------------------------------------------
//const_iterator methods

template<class BV>
sparse_vector_float<BV>::const_iterator::const_iterator() BMNOEXCEPT
: sv_(0), pos_(bm::id_max), exp_it_(), mant_it_() 
{}

//---------------------------------------------------------------------

template<class BV>
sparse_vector_float<BV>::const_iterator::const_iterator(const sparse_vector_type* sv) BMNOEXCEPT
: sv_(sv), exp_it_(sv->exponents_.begin()), mant_it_(sv->mantissas_.begin())
{
    BM_ASSERT(sv_);
    pos_ = sv_->empty() ? bm::id_max : 0u;
}

//---------------------------------------------------------------------

template<class BV>
sparse_vector_float<BV>::const_iterator::const_iterator(const sparse_vector_type* sv, size_type pos) BMNOEXCEPT
: sv_(sv), exp_it_(sv->exponents_.begin()), mant_it_(sv->mantissas_.begin()) 
{
    BM_ASSERT(sv_);
    this->go_to(pos);
}

//---------------------------------------------------------------------

template<class BV>
sparse_vector_float<BV>::const_iterator::const_iterator(const const_iterator& it) BMNOEXCEPT
: sv_(it.sv_), pos_(it.pos_), exp_it_(it.exp_it_), mant_it_(it.mant_it_) 
{}

//---------------------------------------------------------------------

template<class BV>
sparse_vector_float<BV>::const_iterator::value_type sparse_vector_float<BV>::const_iterator::value() const
{
    unsigned int sign = sv_->signs_.test(pos_) ? 1 : 0;
    unsigned int exponent = exp_it_.value();
    unsigned int mantissa = mant_it_.value();
    
    unsigned int bits = (sign << 31) | (exponent << 23) | mantissa;

    float toReturn;
    memcpy(&toReturn, &bits, sizeof(float));

    return toReturn;
}

//---------------------------------------------------------------------

template<class BV>
bool sparse_vector_float<BV>::const_iterator::is_null() const BMNOEXCEPT
{
    return false; // no such thing as null
}

//---------------------------------------------------------------------

template<class BV>
void sparse_vector_float<BV>::const_iterator::go_to(size_type pos) BMNOEXCEPT
{
    pos_ = (!sv_ || pos >= sv_->size()) ? bm::id_max : pos;
    exp_it_.go_to(pos_);
    mant_it_.go_to(pos_);
}

//---------------------------------------------------------------------

template<class BV>
bool sparse_vector_float<BV>::const_iterator::advance() BMNOEXCEPT
{
    if (pos_ == bm::id_max) // nothing to do, we are at the end
        return false;
    
    ++pos_;
    if (pos_ >= sv_->size())
    {
        this->invalidate();
        return false;
    }

    ++exp_it_;
    ++mant_it_;

    return true;
}

//---------------------------------------------------------------------

}//namespace bm

#endif
