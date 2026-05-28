
#include <memory.h>

#ifndef BM_SPARSE_VEC_FLOAT
#define BM_SPARSE_VEC_FLOAT

#include "bm.h"
#include "bmsparsevec.h"

typedef bm::sparse_vector<unsigned int, bm::bvector<>> sparse_vector_u32;

namespace bm{

class sparse_vector_float{
    friend class sparse_vector_float_serialized;
public:
    //
    typedef float                                       value_type;
    typedef bm::bvector<>                               bvector_type;
    typedef bvector_type*                               bvector_type_ptr;
    typedef typename bvector_type::size_type            size_type;
    typedef typename bvector_type::block_idx_type       block_idx_type;
    typedef const bvector_type*                         bvector_type_const_ptr;
	typedef const value_type&                           const_reference;
    typedef typename bvector_type::allocator_type       allocator_type;
    typedef typename bvector_type::allocation_policy    allocation_policy_type;
    typedef typename bvector_type::enumerator           bvector_enumerator_type;
    typedef typename allocator_type::allocator_pool_type    allocator_pool_type;
    typedef bm::basic_bmatrix<bvector_type>                 bmatrix_type;
    
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

    sparse_vector_float();
    sparse_vector_float(const sparse_vector_float&);
    ~sparse_vector_float();

    sparse_vector_float& operator=(const sparse_vector_float& svf);
    bool operator==(const sparse_vector_float& svf) const;
    bool operator!=(const sparse_vector_float& svf) const;


    void swap(sparse_vector_float& svf);
    

    /*! \brief return size of the vector
        \return size of sparse vector
    */
    size_type size() const BMNOEXCEPT { return this->size_; };

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

protected:
    enum buf_size_e{
        n_buf_size = 1024 * 8
    };

private:
    size_type size_ = 0;
    bm::bvector<> signs;
    sparse_vector_u32 exponents;
    sparse_vector_u32 mantissas;
    
};

//---------------------------------------------------------------------
//sparse_vector_float methods

sparse_vector_float::sparse_vector_float(){}

//---------------------------------------------------------------------

sparse_vector_float::sparse_vector_float(const sparse_vector_float&){}

//---------------------------------------------------------------------

sparse_vector_float::~sparse_vector_float(){}

//---------------------------------------------------------------------

sparse_vector_float::const_iterator sparse_vector_float::begin() const{
    typedef typename sparse_vector_float::const_iterator it_type;
    return it_type(this);
}

//---------------------------------------------------------------------

sparse_vector_float& sparse_vector_float::operator=(const sparse_vector_float& svf){
    signs = svf.signs;
    exponents = svf.exponents;
    mantissas = svf.mantissas;
    size_ = svf.size_;
    return *this;
}

//---------------------------------------------------------------------

bool sparse_vector_float::operator==(const sparse_vector_float& svf) const{
    if (size_ == svf.size_ && signs == svf.signs){
        //(exponents == svf.exponents) && (mantissas == svf.mantissas)
        for(size_t i = 0; i < size_; i++){
            unsigned int thisSVF = exponents.get(i);
            unsigned int otherSVF = svf.exponents.get(i);

            if(thisSVF != otherSVF) return false;

            thisSVF = mantissas.get(i);
            otherSVF = svf.mantissas.get(i);

            if(thisSVF != otherSVF) return false;
        }

        return true;
    }
    return false;
}

//---------------------------------------------------------------------

bool sparse_vector_float::operator!=(const sparse_vector_float& svf) const{
    return !operator==(svf);
}

//---------------------------------------------------------------------
void sparse_vector_float::swap(sparse_vector_float &svf){
    signs.swap(svf.signs);
    exponents.swap(svf.exponents);
    mantissas.swap(svf.mantissas);
    std::swap(size_, svf.size_);
}

//---------------------------------------------------------------------

typename sparse_vector_float::value_type 
sparse_vector_float::get(typename sparse_vector_float::size_type idx) const BMNOEXCEPT{
    unsigned int sign = signs.test(idx) ? 1 : 0;
    unsigned int exponent = exponents.get(idx);
    unsigned int mantissa = mantissas.get(idx);
    
    unsigned int bits = (sign << 31) | (exponent << 23) | mantissa;

    float toReturn;
    memcpy(&toReturn, &bits, sizeof(float));

    return toReturn;
}

//---------------------------------------------------------------------

void sparse_vector_float::push_back(value_type v)
{
    unsigned int bits;
    memcpy(&bits, &v, sizeof(float));

    unsigned int sign     = (bits >> 31) & 0x1;
    if(sign == 1) signs.set(size_);
    unsigned int exponent = (bits >> 23) & 0xFF;
    exponents.push_back(exponent);
    unsigned int mantissa =  bits        & 0x7FFFFF;
    mantissas.push_back(mantissa);

    ++(this->size_);
}

//---------------------------------------------------------------------

void sparse_vector_float::set(size_type idx, value_type v){
    if (idx >= size()) this->size_ = idx+1;

    unsigned int bits;
    memcpy(&bits, &v, sizeof(float));

    unsigned int sign     = (bits >> 31) & 0x1;
    if(sign == 1) signs.set(idx);
    else signs.clear(idx);
    unsigned int exponent = (bits >> 23) & 0xFF;
    exponents.set(idx, exponent);
    unsigned int mantissa =  bits        & 0x7FFFFF;
    mantissas.set(idx, mantissa);
}

//---------------------------------------------------------------------

void sparse_vector_float::import(const value_type* arr, size_type arr_size, size_type offset, bool set_not_null){
    if (offset > size_) {
        for (size_type i = size_; i < offset; ++i) {
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

        if (sign == 1)
            signs.set(idx);
        else
            signs.clear(idx);

        exponents.set(idx, exponent);
        mantissas.set(idx, mantissa);
    }

    size_type new_size = offset + arr_size;
    if (new_size > size_)
        size_ = new_size;
}

//---------------------------------------------------------------------

void sparse_vector_float::optimize(bm::word_t* temp_block, typename bvector_type::optmode opt_mode){
    signs.optimize(temp_block, opt_mode);
    exponents.optimize(temp_block, opt_mode);
    mantissas.optimize(temp_block, opt_mode);
}

//---------------------------------------------------------------------

//---------------------------------------------------------------------
//const_iterator methods

sparse_vector_float::const_iterator::const_iterator() BMNOEXCEPT
: sv_(0), pos_(bm::id_max), exp_it_(), mant_it_() {}

//---------------------------------------------------------------------

sparse_vector_float::const_iterator::const_iterator(const sparse_vector_type* sv) BMNOEXCEPT
: sv_(sv), exp_it_(sv->exponents.begin()), mant_it_(sv->mantissas.begin()){
    BM_ASSERT(sv_);
    pos_ = sv_->empty() ? bm::id_max : 0u;
}

//---------------------------------------------------------------------

sparse_vector_float::const_iterator::const_iterator(const sparse_vector_type* sv, size_type pos) BMNOEXCEPT
: sv_(sv), exp_it_(sv->exponents.begin()), mant_it_(sv->mantissas.begin()) {
    BM_ASSERT(sv_);
    this->go_to(pos);
}

//---------------------------------------------------------------------

sparse_vector_float::const_iterator::const_iterator(const const_iterator& it) BMNOEXCEPT
: sv_(it.sv_), pos_(it.pos_), exp_it_(it.exp_it_), mant_it_(it.mant_it_) {}

//---------------------------------------------------------------------

sparse_vector_float::const_iterator::value_type sparse_vector_float::const_iterator::value() const{
    unsigned int sign = sv_->signs.test(pos_) ? 1 : 0;
    unsigned int exponent = exp_it_.value();
    unsigned int mantissa = mant_it_.value();
    
    unsigned int bits = (sign << 31) | (exponent << 23) | mantissa;

    float toReturn;
    memcpy(&toReturn, &bits, sizeof(float));

    return toReturn;
}

//---------------------------------------------------------------------

bool sparse_vector_float::const_iterator::is_null() const BMNOEXCEPT{
    return false; // no such thing as null
}

//---------------------------------------------------------------------

void sparse_vector_float::const_iterator::go_to(size_type pos) BMNOEXCEPT{
    pos_ = (!sv_ || pos >= sv_->size()) ? bm::id_max : pos;
    exp_it_.go_to(pos_);
    mant_it_.go_to(pos_);
}

//---------------------------------------------------------------------

bool sparse_vector_float::const_iterator::advance() BMNOEXCEPT{
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
