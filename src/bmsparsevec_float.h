
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
#include <cstring>
#include <cmath>

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
   \brief Wrapper for sparse vector in order to add float support
    
   This wrapper takes in floats and splits them up into their sign, exponent, and mantissa
   It stores the sign in a bvector, and the exponent and mantissa in sparse vectors as unsigned ints
 
   @ingroup sv
*/

template<class T>
struct is_rsc_sparse_vector : std::false_type {};

template<class Val, class SV>
struct is_rsc_sparse_vector<bm::rsc_sparse_vector<Val, SV>> : std::true_type {};

template<class SV>
class sparse_vector_float
{
    template<class SVect, unsigned S_FACTOR> friend class sparse_vector_scanner;
    template<class SVect> friend class sparse_vector_float_serializer;
    template<class SVect> friend class sparse_vector_float_deserializer;
public:
    //
    typedef float                                   value_type;
    typedef typename SV::bvector_type               bvector_type;
    typedef typename bvector_type::size_type        size_type;
    typedef SV                                      sparse_vector_u;
    typedef typename bvector_type::allocator_type   allocator_type;
    typedef typename bvector_type::allocation_policy allocation_policy_type;

    struct statistics : public bv_statistics
    {};

    /**
         Reference class to access elements via common [] operator
         @ingroup sv
    */
    class reference
    {
    public:
        reference(sparse_vector_float<SV>& svf, size_type idx) BMNOEXCEPT
        : svf_(svf), idx_(idx)
        {}

        operator value_type() const BMNOEXCEPT { return svf_.get(idx_); }

        reference& operator=(const reference& ref)
        {
            svf_.set(idx_, (value_type)ref);
            return *this;
        }

        reference& operator=(value_type val)
        {
            svf_.set(idx_, val);
            return *this;
        }

        bool operator==(const reference& ref) const BMNOEXCEPT
                                { return bool(*this) == bool(ref); }

        //bool is_null() const BMNOEXCEPT { return svf_.is_null(idx_); }
    private:
        sparse_vector_float<SV>& svf_;
        size_type               idx_;
    };


    /*
        Const iterator for traversing the sparse_vector_float

        Implementation uses the const_iterators for the exponent and mantissa sparse_vectors
        and a reference to the sparse_vector_float for signs
    */
    class const_iterator
    {
    public:
    friend class sparse_vector_float;

#ifndef BM_NO_STL
        typedef std::input_iterator_tag  iterator_category;
#endif
        typedef sparse_vector_float                        sparse_vector_type;
        typedef typename sparse_vector_type::value_type    value_type;
        typedef typename sparse_vector_type::size_type     size_type;

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
        const sparse_vector_type*         sv_;                ///!< ptr to parent
        size_type                         pos_;               ///!< Position
        typename sparse_vector_u::const_iterator   exp_it_;   ///!< exponent iterator
        typename sparse_vector_u::const_iterator   mant_it_;  ///!< mantissa iterator
    };


    /**
        Back insert iterator for inserting to the end of the sparse_vector_float
    */
    class back_insert_iterator
    {
    public:
        typedef sparse_vector_float<SV>                    sparse_vector_type;
        typedef typename sparse_vector_type::value_type    value_type;

        typedef void difference_type;
        typedef void pointer;
        typedef void reference;
        
    public:
        /*! @name Construction and assignment  */
        back_insert_iterator();
        back_insert_iterator(sparse_vector_type* svf);
        back_insert_iterator(const back_insert_iterator& bi);
        
        /** move constructor */
        back_insert_iterator(back_insert_iterator&& bi) BMNOEXCEPT;

        ~back_insert_iterator();

        /** push value to the vector */
        //back_insert_iterator&
        void operator=(value_type v) { this->add(v); /*return *this;*/ }
        /** noop */
        back_insert_iterator& operator*() { return *this; }
        /** noop */
        back_insert_iterator& operator++() { return *this; }
        /** noop */
        back_insert_iterator& operator++( int ) { return *this; }
        
        /** add value to the container*/
        void add(value_type v);

    private:
        bm::sparse_vector_float<SV>* svf_ = 0;      ///!< pointer on the parent vector
    };

    /*! \brief returns a const_iterator for this sparse_vector_float pointing to the first index */
    const_iterator begin() const;

    /*! \brief returns a const_iterator for this sparse_vector_float pointing to the last index */
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

    /*! \brief copy constructor*/
    sparse_vector_float(const sparse_vector_float& svf);
    ~sparse_vector_float();

    sparse_vector_float& operator=(const sparse_vector_float& svf);
    bool operator==(const sparse_vector_float& svf) const;
    bool operator!=(const sparse_vector_float& svf) const;


    /*! \brief swaps the elements in this sparse_vector_float and the given sparse_vector_float */
    void swap(sparse_vector_float& svf);

    /*! \brief return size of the vector
        \return size of sparse vector float
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
    sparse_vector_float<SV>& clear_range(size_type left,
                                        size_type right,
                                        bool set_null = false);

     /*!
        \brief check if another sparse vector float has the same content and size
     
        \param sv        - sparse vector for comparison
        \param null_able - flag to consider NULL vector in comparison (default)
                           or compare only value content planes
     
        \return true, if it is the same
    */
    bool equal(const sparse_vector_float<SV>& sv,
               bm::null_support null_able = bm::use_null) const BMNOEXCEPT;

    /**
        \brief Compare vector element with argument
     
        \param idx - vector element index
        \param val - argument to compare with
        \param epsilon - amount of precision
     
        \return 0 - equal, -1 - vect[i] < val, 1 otherwise
    */
    int compare(size_type idx, const value_type val, float epsilon = std::numeric_limits<float>::epsilon()) const BMNOEXCEPT;

    /*!
        \brief join all with another sparse vector float using OR operation
            NOTE: if you join 2 floats together and one is not 0.0, ie -1.5f and -2.5f
            it is possible for the float to become NaN if the exponent becomes entirely 11111111
        \param sv - argument vector to join with
        \return slf reference
        @sa merge
    */
    sparse_vector_float<SV>& join(const sparse_vector_float<SV>& svf);

    /*!
        \brief merge with another sparse vector float using OR operation
        Merge is different from join(), because it borrows data from the source
        vector, so it gets modified.
     
        \param sv - [in, out]argument vector to join with (vector mutates)
     
        \return slf reference
        @sa join
    */
    sparse_vector_float<SV>& merge(sparse_vector_float<SV>& svf);

    /**
        @brief copy range of values from another sparse vector float
     
        Copy [left..right] values from the source vector,
        clear everything outside the range.
     
        \param sv - source vector
        \param left  - index from in losed diapason of [left..right]
        \param right - index to in losed diapason of [left..right]
        \param slice_null - "use_null" copy range for NULL vector or
                             do not copy it
    */
    void copy_range(const sparse_vector_float<SV>& svf,
                    size_type left, size_type right,
                    bm::null_support slice_null = bm::use_null);


    /*!
        \brief Import list of elements from a C-style array
        \param arr  - source array
        \param arr_size - source size
        \param offset - target index in the sparse vector float
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
        struct sparse_vector_float<SV>::statistics* st) const BMNOEXCEPT;

    /*! \brief syncronize internal structures, build fast access index
    */
    void sync(bool /*force*/, bool /*sync_size*/);

    /*!
        \brief Bulk export list of elements to a C-style array
     
        Use of all extract() methods is restricted.
        Please consider decode() for the same purpose.
     
        \param arr  - dest array
        \param size - dest size
        \param offset - target index in the sparse vector to export from
        \param zero_mem - set to false if target array is pre-initialized
                          with 0s to avoid performance penalty   
        \return effective size(number) of exported elements
     
        \sa decode
     
        @internal
    */
    size_type extract(value_type* arr,
                      size_type size,
                      size_type offset = 0,
                      bool      zero_mem = true) const BMNOEXCEPT2;

    /** \brief extract small window without use of masking vector
        \sa decode
        @internal
    */
    size_type extract_range(value_type* arr,
                            size_type size,
                            size_type offset,
                            bool      zero_mem = true) const;

    /*!
        \brief Bulk export list of elements to a C-style array
     
        For efficiency, this is left as a low level function,
        it does not do any bounds checking on the target array, it will
        override memory and crash if you are not careful with allocation
        and request size.
     
        \param arr  - dest array
        \param idx_from - index in the sparse vector to export from
        \param dec_size - decoding size (array allocation should match)
        \param zero_mem - set to false if target array is pre-initialized
                          with 0s to avoid performance penalty
     
        \return number of actually exported elements (can be less than requested)
     
        \sa gather
    */
    size_type decode(value_type* arr,
                     size_type   idx_from,
                     size_type   dec_size,
                     bool        zero_mem = true) const;

    /*!
        \brief Gather elements to a C-style array
     
        Gather collects values from different locations, for best
        performance feed it with sorted list of indexes.
     
        Faster than one-by-one random access.
     
        For efficiency, this is left as a low level function,
        it does not do any bounds checking on the target array, it will
        override memory and crash if you are not careful with allocation
        and request size.
     
        \param arr  - dest array
        \param idx - index list to gather elements
        \param size - decoding index list size (array allocation should match)
        \param sorted_idx - sort order directive for the idx array
                            (BM_UNSORTED, BM_SORTED, BM_UNKNOWN)
        Sort order affects both performance and correctness(!), use BM_UNKNOWN
        if not sure.
     
        \return number of actually exported elements (can be less than requested)
     
        \sa decode
    */
    size_type gather(value_type* arr,
                     const size_type* idx,
                     size_type   size,
                     bm::sort_order sorted_idx) const;


    /**
        @brief Turn sparse vector float into immutable mode
        Read-only (immutable) vector uses less memory and allows faster searches.
        Before freezing it is recommenede to call optimize() to get full memory saving effect
        @sa optimize
     */
    void freeze();

    /** Returns true if vector is in read-only mode.
         @sa freeze
    */
    bool is_ro() const BMNOEXCEPT;

protected:
    enum buf_size_e{
        n_buf_size = 1024 * 8
    };
//private:
public:
    bvector_type      signs_;      ///!< sign bit vector
    sparse_vector_u   exponents_;  ///!< exponent sparse vector
    sparse_vector_u   mantissas_;  ///!< mantissa sparse vector
    
};

//---------------------------------------------------------------------
//sparse_vector_float methods

template<class SV>
sparse_vector_float<SV>::sparse_vector_float(bm::null_support null_able,
                                             allocation_policy_type ap,
                                             size_type bv_max_size,
                                             const allocator_type&   alloc)
:signs_(ap.strat, ap.glevel_len, bv_max_size, alloc),
 exponents_(null_able, ap, bv_max_size, alloc),
 mantissas_(null_able, ap, bv_max_size, alloc)
{}

//---------------------------------------------------------------------

template<class SV>
sparse_vector_float<SV>::sparse_vector_float(const sparse_vector_float& svf)
:signs_(svf.signs_), exponents_(svf.exponents_), mantissas_(svf.mantissas_)
{}

//---------------------------------------------------------------------

template<class SV>
sparse_vector_float<SV>::~sparse_vector_float()
{}

//---------------------------------------------------------------------

template<class SV>
typename sparse_vector_float<SV>::const_iterator sparse_vector_float<SV>::begin() const
{
    typedef typename sparse_vector_float::const_iterator it_type;
    return it_type(this);
}

//---------------------------------------------------------------------

template<class SV>
sparse_vector_float<SV>& sparse_vector_float<SV>::operator=(const sparse_vector_float& svf)
{
    signs_ = svf.signs_;
    exponents_ = svf.exponents_;
    mantissas_ = svf.mantissas_;
    return *this;
}

//---------------------------------------------------------------------

template<class SV>
bool sparse_vector_float<SV>::operator==(const sparse_vector_float<SV>& svf) const
{
    if (mantissas_.size() == svf.mantissas_.size() && signs_ == svf.signs_)
    {
        for (size_type i = 0; i < mantissas_.size(); i++)
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

template<class SV>
bool sparse_vector_float<SV>::operator!=(const sparse_vector_float<SV>& svf) const
{
    return !operator==(svf);
}

//---------------------------------------------------------------------

template<class SV>
void sparse_vector_float<SV>::swap(sparse_vector_float<SV> &svf)
{
    if constexpr (!is_rsc_sparse_vector<SV>::value)
    {
        signs_.swap(svf.signs_);
        exponents_.swap(svf.exponents_);
        mantissas_.swap(svf.mantissas_);
    }
}

//---------------------------------------------------------------------

template<class SV>
typename sparse_vector_float<SV>::value_type 
sparse_vector_float<SV>::get(typename sparse_vector_float<SV>::size_type idx) const BMNOEXCEPT
{
    unsigned int sign = signs_.test(idx) ? 1 : 0;
    unsigned int exponent = exponents_.get(idx);
    unsigned int mantissa = mantissas_.get(idx);
    
    unsigned int bits = (sign << 31) | (exponent << 23) | mantissa;

    float toReturn;
    std::memcpy(&toReturn, &bits, sizeof(float));

    return toReturn;
}

//---------------------------------------------------------------------

template<class SV>
void sparse_vector_float<SV>::push_back(value_type v)
{
    unsigned int bits;
    std::memcpy(&bits, &v, sizeof(float));

    unsigned int sign     = (bits >> 31) & 0x1;
    if (sign == 1) signs_.set(mantissas_.size());
    unsigned int exponent = (bits >> 23) & 0xFF;
    exponents_.push_back(exponent);
    unsigned int mantissa =  bits        & 0x7FFFFF;
    mantissas_.push_back(mantissa);
}

//---------------------------------------------------------------------

template<class SV>
void sparse_vector_float<SV>::set(size_type idx, value_type v)
{

    unsigned int bits;
    std::memcpy(&bits, &v, sizeof(float));

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

template<class SV>
void sparse_vector_float<SV>::clear() BMNOEXCEPT
{
    signs_.clear();
    exponents_.clear();
    mantissas_.clear();
}

//---------------------------------------------------------------------

template<class SV>
sparse_vector_float<SV>& sparse_vector_float<SV>::clear_range(size_type left,
                                    size_type right,
                                    bool set_null)
{
    if constexpr (!is_rsc_sparse_vector<SV>::value)
    {
        signs_.clear_range(left, right);
        exponents_.clear_range(left, right, set_null);
        mantissas_.clear_range(left, right, set_null);
    }

    return *this;
}

//---------------------------------------------------------------------

template<class SV>
bool sparse_vector_float<SV>::equal(const sparse_vector_float<SV>& sv,
                                    bm::null_support null_able) const BMNOEXCEPT
{
    if constexpr (!is_rsc_sparse_vector<SV>::value)
    {
        if (signs_.equal(sv.signs_) && exponents_.equal(sv.exponents_, null_able) && mantissas_.equal(sv.mantissas_, null_able)) return true;
        return false;
    }

    if (signs_.equal(sv.signs_) && exponents_.equal(sv.exponents_) && mantissas_.equal(sv.mantissas_)) return true;
    return false;
}

//---------------------------------------------------------------------

template<class SV>
int sparse_vector_float<SV>::compare(size_type idx, const value_type val, float epsilon) const BMNOEXCEPT
{
    value_type svf_value = get(idx);

    float diff = svf_value - val;

    if (std::fabs(diff) <= epsilon)     return 0;
    else if (diff > 0)                  return 1;
    else                                return -1;
}

//---------------------------------------------------------------------


template<class SV>
sparse_vector_float<SV>& sparse_vector_float<SV>::join(const sparse_vector_float<SV>& svf)
{
    if constexpr (!is_rsc_sparse_vector<SV>::value)
    {
        signs_ |= svf.signs_;
        exponents_.join(svf.exponents_);
        mantissas_.join(svf.mantissas_);
    }
    
    return *this;
}

//---------------------------------------------------------------------

template<class SV>
sparse_vector_float<SV>& sparse_vector_float<SV>::merge(sparse_vector_float<SV>& svf)
{
    if constexpr (!is_rsc_sparse_vector<SV>::value)
    {
        signs_.merge(svf.signs_);
        exponents_.merge(svf.exponents_);
        mantissas_.merge(svf.mantissas_);
    }
    return *this;
}


//---------------------------------------------------------------------

template<class SV>
void sparse_vector_float<SV>::copy_range(const sparse_vector_float<SV>& svf,
                                        size_type left, size_type right,
                                        bm::null_support slice_null)
{
    signs_.copy_range(svf.signs_, left, right);
    if constexpr (!is_rsc_sparse_vector<SV>::value)
    {
        exponents_.copy_range(svf.exponents_, left, right, slice_null);
        mantissas_.copy_range(svf.mantissas_, left, right, slice_null);
    }else{
        exponents_.copy_range(svf.exponents_, left, right);
        mantissas_.copy_range(svf.mantissas_, left, right);
    }
}


//---------------------------------------------------------------------

template<class SV>
void sparse_vector_float<SV>::import(const value_type* arr, 
                                size_type arr_size, 
                                size_type offset, 
                                bool set_not_null)
{
    const size_type CHUNK = 256;
    unsigned int exp_buf[CHUNK];
    unsigned int mant_buf[CHUNK];
 
    size_type remaining = arr_size;
    size_type pos = 0;

    if constexpr (is_rsc_sparse_vector<SV>::value)
    {
        typename sparse_vector_u::sparse_vector_type tempExpSV(bm::use_null);
        typename sparse_vector_u::sparse_vector_type tempMantSV(bm::use_null);
        bm::bvector<> nullBits(CHUNK);

        while (remaining > 0)
        {
            size_type chunk = std::min(remaining, CHUNK);
            size_type idx = offset + pos;
    
            for (size_type i = 0; i < chunk; i++)
            {
                if(arr[pos+i] != arr[pos+i]){
                    nullBits.set(i, true);
                    continue;
                }
                unsigned int bits;
                std::memcpy(&bits, &arr[pos + i], sizeof(float));
    
                unsigned int sign = (bits >> 31) & 0x1;
                exp_buf[i]  = (bits >> 23) & 0xFF;
                mant_buf[i] =  bits        & 0x7FFFFF;
    
                signs_.set(idx + i, sign == 1);
            }

            if(nullBits.count() != CHUNK)
            {
                tempExpSV.import(exp_buf, chunk, idx, set_not_null);
                tempMantSV.import(mant_buf, chunk, idx, set_not_null);

                bm::bvector<>::enumerator first = nullBits.first();
                bm::bvector<>::enumerator end = nullBits.end();
                while(first != end)
                {
                    tempExpSV.clear(*first + pos, true);
                    tempMantSV.clear(*first + pos, true);
                    ++first;
                }
            }
            
            pos       += chunk;
            remaining -= chunk;
            nullBits.clear();
        }

        exponents_.load_from(tempExpSV);
        mantissas_.load_from(tempMantSV);
    }
    else
    {
        while (remaining > 0)
        {
            size_type chunk = std::min(remaining, CHUNK);
            size_type idx = offset + pos;
    
            for (size_type i = 0; i < chunk; i++)
            {
                unsigned int bits;
                std::memcpy(&bits, &arr[pos + i], sizeof(float));
    
                unsigned int sign = (bits >> 31) & 0x1;
                exp_buf[i]  = (bits >> 23) & 0xFF;
                mant_buf[i] =  bits        & 0x7FFFFF;
    
                signs_.set(idx + i, sign == 1);
            }

            exponents_.import(exp_buf, chunk, idx, set_not_null);
            mantissas_.import(mant_buf, chunk, idx, set_not_null);

    
            pos       += chunk;
            remaining -= chunk;
        }
    }
}

//---------------------------------------------------------------------

template<class SV>
void sparse_vector_float<SV>::optimize(bm::word_t* temp_block, typename bvector_type::optmode opt_mode)
{
    signs_.optimize(temp_block, opt_mode);
    exponents_.optimize(temp_block, opt_mode);
    mantissas_.optimize(temp_block, opt_mode);
}

//---------------------------------------------------------------------

template<class SV>
void sparse_vector_float<SV>::calc_stat(struct sparse_vector_float<SV>::statistics* st) const BMNOEXCEPT
{
    BM_ASSERT(st);

    bm::bvector<>::statistics signStat;
    signs_.calc_stat(signStat);

    typename sparse_vector_u::statistics expStat, mantStat;
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
    for (int i = 0; i < bm::gap_levels; i++)
    {
        st->gap_levels[i] = std::max({signStat.gap_levels[i], 
                                      expStat.gap_levels[i], 
                                      mantStat.gap_levels[i]});
        st->gaps_by_level[i] = signStat.gaps_by_level[i] 
                              + expStat.gaps_by_level[i] 
                              + mantStat.gaps_by_level[i];
    }
}

//---------------------------------------------------------------------

template<class SV>
void sparse_vector_float<SV>::sync(bool /*force*/, bool /*sync_size*/)
{
    exponents_.sync();
    mantissas_.sync();
}

//---------------------------------------------------------------------

template<class SV>
typename sparse_vector_float<SV>::size_type 
sparse_vector_float<SV>::extract(value_type* arr,
                                 size_type size,
                                 size_type offset,
                                 bool      zero_mem) const BMNOEXCEPT2
{
    if constexpr (!is_rsc_sparse_vector<SV>::value)
    {
        const size_type CHUNK = 256;

        unsigned int exp_buf[CHUNK];
        unsigned int mant_buf[CHUNK];

        size_type remaining = size;
        size_type pos = offset;
        size_type toReturn = 0;

        while (remaining > 0)
        {
            size_type chunk = std::min(remaining, CHUNK);

            toReturn += exponents_.extract(exp_buf, chunk, pos, zero_mem);
            toReturn += mantissas_.extract(mant_buf, chunk, pos, zero_mem);

            for (size_type i = 0; i < chunk; i++)
            {
                unsigned int sign     = signs_.test(pos + i) ? 1 : 0;
                unsigned int exponent = exp_buf[i];
                unsigned int mantissa = mant_buf[i];

                unsigned int bits = (sign << 31) | (exponent << 23) | mantissa;
                std::memcpy(&arr[pos - offset + i], &bits, sizeof(float));
            }

            pos       += chunk;
            remaining -= chunk;
        }
        return toReturn;
    }
}

//---------------------------------------------------------------------

template<class SV>
typename sparse_vector_float<SV>::size_type 
sparse_vector_float<SV>::extract_range(value_type* arr,
                                       size_type size,
                                       size_type offset,
                                       bool      zero_mem) const
{
    if constexpr (!is_rsc_sparse_vector<SV>::value)
    {
        const size_type CHUNK = 256;

        unsigned int exp_buf[CHUNK];
        unsigned int mant_buf[CHUNK];

        size_type remaining = size;
        size_type pos = offset;
        size_type toReturn = 0;

        while (remaining > 0)
        {
            size_type chunk = std::min(remaining, CHUNK);

            toReturn += exponents_.extract_range(exp_buf, chunk, pos, zero_mem);
            toReturn += mantissas_.extract_range(mant_buf, chunk, pos, zero_mem);

            for (size_type i = 0; i < chunk; i++)
            {
                unsigned int sign     = signs_.test(pos + i) ? 1 : 0;
                unsigned int exponent = exp_buf[i];
                unsigned int mantissa = mant_buf[i];

                unsigned int bits = (sign << 31) | (exponent << 23) | mantissa;
                std::memcpy(&arr[pos - offset + i], &bits, sizeof(float));
            }

            pos       += chunk;
            remaining -= chunk;
        }
        return toReturn;
    }
}

//---------------------------------------------------------------------

template<class SV>
typename sparse_vector_float<SV>::size_type
sparse_vector_float<SV>::decode(value_type* arr,
                                size_type   idx_from,
                                size_type   dec_size,
                                bool        zero_mem) const
{
    if constexpr (!is_rsc_sparse_vector<SV>::value)
    {
        return extract(arr, dec_size, idx_from, zero_mem);
    }
    else
    {
        const size_type CHUNK = 256;
        unsigned int exp_buf[CHUNK];
        unsigned int mant_buf[CHUNK];
        size_type remaining = dec_size;
        size_type pos = idx_from;
        size_type toReturn = 0;
 
        while (remaining > 0)
        {
            size_type chunk = std::min(remaining, CHUNK);
 
            toReturn += exponents_.decode(exp_buf, pos, chunk, zero_mem);
            toReturn += mantissas_.decode(mant_buf, pos, chunk, zero_mem);
 
            for (size_type i = 0; i < chunk; i++)
            {
                unsigned int sign     = signs_.test(pos + i) ? 1 : 0;
                unsigned int exponent = exp_buf[i];
                unsigned int mantissa = mant_buf[i];
                unsigned int bits = (sign << 31) | (exponent << 23) | mantissa;
                std::memcpy(&arr[pos - idx_from + i], &bits, sizeof(float));
            }
 
            pos       += chunk;
            remaining -= chunk;
        }
        return toReturn;
    }
}

//---------------------------------------------------------------------

template<class SV>
typename sparse_vector_float<SV>::size_type 
sparse_vector_float<SV>::gather(value_type* arr,
                                const size_type* idx,
                                size_type   size,
                                bm::sort_order sorted_idx) const
{
    const size_type CHUNK = 256;

    unsigned int exp_buf[CHUNK];
    unsigned int mant_buf[CHUNK];

    size_type remaining = size;
    size_type pos = 0;
    size_type toReturn = 0;

    while (remaining > 0)
    {
        size_type chunk = std::min(remaining, CHUNK);

        if constexpr (!is_rsc_sparse_vector<SV>::value)
        {
            toReturn += exponents_.gather(exp_buf, idx + pos, chunk, sorted_idx);
            toReturn += mantissas_.gather(mant_buf, idx + pos, chunk, sorted_idx);
        }
        else
        {
            size_type idx_tmp_buf[CHUNK];

            toReturn += exponents_.gather(exp_buf, idx + pos, idx_tmp_buf, chunk, sorted_idx);
            toReturn += mantissas_.gather(mant_buf, idx + pos, idx_tmp_buf, chunk, sorted_idx);
        }

        for (size_type i = 0; i < chunk; i++)
        {
            unsigned int sign     = signs_.test(idx[pos + i]) ? 1 : 0;
            unsigned int exponent = exp_buf[i];
            unsigned int mantissa = mant_buf[i];

            unsigned int bits = (sign << 31) | (exponent << 23) | mantissa;
            std::memcpy(&arr[pos + i], &bits, sizeof(float));
        }

        remaining -= chunk;
        pos += chunk;
        
    }
    return toReturn;
}

//---------------------------------------------------------------------

template<class SV>
void sparse_vector_float<SV>::freeze()
{
    signs_.freeze();
    exponents_.freeze();
    mantissas_.freeze();
}

//---------------------------------------------------------------------

template<class SV>
bool sparse_vector_float<SV>::is_ro() const BMNOEXCEPT
{
    return mantissas_.is_ro();
}

//---------------------------------------------------------------------
//const_iterator methods

template<class SV>
sparse_vector_float<SV>::const_iterator::const_iterator() BMNOEXCEPT
: sv_(0), pos_(bm::id_max), exp_it_(), mant_it_() 
{}

//---------------------------------------------------------------------

template<class SV>
sparse_vector_float<SV>::const_iterator::const_iterator(const sparse_vector_type* sv) BMNOEXCEPT
: sv_(sv), exp_it_(sv->exponents_.begin()), mant_it_(sv->mantissas_.begin())
{
    BM_ASSERT(sv_);
    pos_ = sv_->empty() ? bm::id_max : 0u;
}

//---------------------------------------------------------------------

template<class SV>
sparse_vector_float<SV>::const_iterator::const_iterator(const sparse_vector_type* sv, size_type pos) BMNOEXCEPT
: sv_(sv), exp_it_(sv->exponents_.begin()), mant_it_(sv->mantissas_.begin()) 
{
    BM_ASSERT(sv_);
    this->go_to(pos);
}

//---------------------------------------------------------------------

template<class SV>
sparse_vector_float<SV>::const_iterator::const_iterator(const const_iterator& it) BMNOEXCEPT
: sv_(it.sv_), pos_(it.pos_), exp_it_(it.exp_it_), mant_it_(it.mant_it_) 
{}

//---------------------------------------------------------------------

template<class SV>
typename sparse_vector_float<SV>::const_iterator::value_type
sparse_vector_float<SV>::const_iterator::value() const
{
    if (is_null())
    {
        return std::numeric_limits<float>::quiet_NaN();
    }
    unsigned int sign = sv_->signs_.test(pos_) ? 1 : 0;
    unsigned int exponent = exp_it_.value();
    unsigned int mantissa = mant_it_.value();
    
    unsigned int bits = (sign << 31) | (exponent << 23) | mantissa;

    float toReturn;
    std::memcpy(&toReturn, &bits, sizeof(float));

    return toReturn;
}

//---------------------------------------------------------------------

template<class SV>
bool sparse_vector_float<SV>::const_iterator::is_null() const BMNOEXCEPT
{
    return mant_it_.is_null();
}

//---------------------------------------------------------------------

template<class SV>
void sparse_vector_float<SV>::const_iterator::go_to(size_type pos) BMNOEXCEPT
{
    pos_ = (!sv_ || pos >= sv_->size()) ? bm::id_max : pos;
    exp_it_.go_to(pos_);
    mant_it_.go_to(pos_);
}

//---------------------------------------------------------------------

template<class SV>
bool sparse_vector_float<SV>::const_iterator::advance() BMNOEXCEPT
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
//back_insert_iterator methods

template<class SV>
sparse_vector_float<SV>::back_insert_iterator::back_insert_iterator()
{}

//---------------------------------------------------------------------

template<class SV>
sparse_vector_float<SV>::back_insert_iterator::back_insert_iterator(sparse_vector_type* svf)
: svf_(svf)
{}

//---------------------------------------------------------------------

template<class SV>
sparse_vector_float<SV>::back_insert_iterator::back_insert_iterator(const back_insert_iterator& bi)
: svf_(bi.svf_)
{}

//---------------------------------------------------------------------

template<class SV>
sparse_vector_float<SV>::back_insert_iterator::~back_insert_iterator()
{}

//---------------------------------------------------------------------

template<class SV>
void sparse_vector_float<SV>::back_insert_iterator::add(value_type v)
{
    svf_->push_back(v);
}

//---------------------------------------------------------------------

template<class SV>
sparse_vector_float<SV>::back_insert_iterator::back_insert_iterator(back_insert_iterator&& bi) BMNOEXCEPT
: svf_(bi.svf_)
{
    bi.svf_ = nullptr;
}

}//namespace bm

#endif
