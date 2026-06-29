#ifndef BM_SPARSE_VEC_FLOAT_SERIAL_INCLUDED
#define BM_SPARSE_VEC_FLOAT_SERIAL_INCLUDED

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

/*! \file sparse_vector_float_serial.h
    \brief Serialization for sparse_vector_float
*/

#include <memory.h>
#include <cstring>

#include "bm.h"

#include "bmsparsevec_float.h"
#include "bmsparsevec.h"
#include "bmsparsevec_serial.h"

namespace bm
{

template<class SV>
class sparse_vector_float_serial_layout
{
    template<class SVT> friend class sparse_vector_float_serializer;
    template<class SVT> friend class sparse_vector_float_deserializer;
public:
    typedef typename SV::value_type         value_type;
    typedef typename SV::bvector_type       bvector_type;
    typedef typename SV::sparse_vector_u    sparse_vector_u;

    /*! \brief constructor for serial_layout s*/
    sparse_vector_float_serial_layout() BMNOEXCEPT {}
    ~sparse_vector_float_serial_layout() {if (buf_)free(buf_);}
    sparse_vector_float_serial_layout(const sparse_vector_float_serial_layout&);

    
    /*!
        \brief resize capacity
        \param capacity - new capacity
        \return new buffer or 0 if failed
    */
    unsigned char* reserve(size_t capacity);

    /// return current serialized size
    size_t size();

    /// Set new serialized size
    void resize(size_t ssize);

    /// Set new serialized size (shrink)
    void shrink(size_t ssize);

    /// return serialization buffer capacity
    size_t  capacity() const BMNOEXCEPT;
    
    /// free memory
    void freemem() BMNOEXCEPT;

    /// Return serialization buffer pointer
    const unsigned char* buf() const BMNOEXCEPT { return buf_;  }

    /// Return serialization buffer pointer
    const unsigned char* data() const BMNOEXCEPT { return buf_;  }
    
private:
    unsigned char* buf_ = nullptr;  ///!< buf storing the serialized sparse_vector_float
    size_t size_ = 0;               ///!< how much space the buf is using
    size_t capacity_ = 0;           ///!< how much space is allocated to the buf
};


template<typename SV>
class sparse_vector_float_serializer
{
public:
    typedef typename SV::bvector_type                bvector_type;
    typedef typename
    bm::serializer<bvector_type>::bv_ref_vector_type bv_ref_vector_type;
    typedef typename
    bm::serializer<bvector_type>::xor_sim_model_type xor_sim_model_type;
    typedef typename SV::sparse_vector_u             sparse_vector_type;

public:
    /// @brief constructor
    sparse_vector_float_serializer() {}

    /**
        Add skip-markers for faster range deserialization

        @param enable - TRUE searilization will add bookmark codes
        @param bm_interval - bookmark interval in (number of blocks)
                            (suggested between 4 and 512)
        smaller interval means more bookmarks added to the skip list thus
        more increasing the BLOB size
    */
    void set_bookmarks(bool enable, unsigned bm_interval = 256) BMNOEXCEPT;

    /**
        Enable XOR compression on vector serialization
        @sa set_xor_ref
        @sa disable_xor_compression
     */
    void enable_xor_compression() BMNOEXCEPT
        { set_xor_ref(true); }

    /**
        Disable XOR compression on serialization
     */
    void disable_xor_compression() BMNOEXCEPT
        { set_xor_ref(false); }

    /** Turn ON and OFF XOR compression of sparse vectors
        Enables XOR reference compression for the sparse vector float.
        Default: disabled
        Reference bit-vectors from the sparse vector itself
    */
    void set_xor_ref(bool is_enabled) BMNOEXCEPT;

    /** Enable external XOR serialization via external reference vectors
       (data frame ref. vector).
       This method is useful when we serialize a group of related
       sparse vectors which benefits from the XOR referencial compression

       @param bv_ref_ptr - external reference vector
       if NULL - resets the use of reference vector
    */
    void set_xor_ref(const bv_ref_vector_type* bv_ref_ptr) BMNOEXCEPT;

    /**
        Calculate XOR similarity model for ref_vector
        refernece vector must be associated before
        @sa set_ref_vectors, set_sim_model
        @internal
     */
    void compute_sim_model(xor_sim_model_type& sim_model,
                           const bv_ref_vector_type& ref_vect,
                           const bm::xor_sim_params& params);

    /**
        Attach serizalizer to a pre-computed similarity model
        @param sim_model - pointer to external computed model
     */
    void set_sim_model(const xor_sim_model_type* sim_model) BMNOEXCEPT;

    /**
        Returns the XOR reference compression status (enabled/disabled)
    */
    bool is_xor_ref() const BMNOEXCEPT;

    /*!
        \brief Serialize sparse vector into a memory buffer(s) structure
     
        \param sv                 - sparse vector to serialize
        \param sv_layout  - buffer structure to keep the result
        as defined in bm::serialization_flags
    */
    void serialize(const SV&                        sv,
                   sparse_vector_float_serial_layout<SV>& sv_layout);


protected:
    bm::serializer<bvector_type>                     signSerializer_;       ///!< serializer for the sign bvector
    bm::sparse_vector_serializer<sparse_vector_type> exponentSerializer_;   ///!< serializer for the exponent sparse_vector
    bm::sparse_vector_serializer<sparse_vector_type> mantissaSerializer_;   ///!< serializer for the mantissa sparse_vector
    
};

/**
    sparse vector de-serializer
*/
template<typename SV>
class sparse_vector_float_deserializer
{
    public:
    typedef typename SV::bvector_type       bvector_type;
    typedef typename SV::size_type          size_type;
    typedef typename bm::serializer<bvector_type>::bv_ref_vector_type bv_ref_vector_type;
    typedef typename SV::sparse_vector_u   sparse_vector_type;

public:
    /// @brief constructor
    sparse_vector_float_deserializer();
    ~sparse_vector_float_deserializer();

    /**
        Set deserialization finalization to force deserialized vectors into READONLY (or READWRITE) mode.
        Performance impact: Turning ON finalization will make deserialization a lit slower,
        because each bit-vector will be re-converted into new mode (READONLY).
        Following (search) operations may perform a bit faster.

        @param is_final - finalization code
                        (use bm::finalization::READONLY to produce an immutable vector)
     */
    void set_finalization(bm::finalization is_final);

    /** Set external XOR reference vectors
        (data frame referenece vectors)

        @param bv_ref_ptr - external reference vector
        if NULL - resets the use of reference
    */
    void set_xor_ref(bv_ref_vector_type* bv_ref_ptr);

    /*!
        Deserialize sparse vector

        @param sv - [out] target sparse vector to populate
        @param buf - input BLOB source memory pointer
        @param clear_sv - if true clears the target vector (sv)

        @sa deserialize_range
    */
    void deserialize(SV& sv,
                     const unsigned char* buf,
                     bool clear_sv = true);

    /*!
        Deserialize sparse vector float for the range [from, to]

        @param sv - [out] target sparse vector float to populate
        @param buf - input BLOB  source memory pointer
        @param from - start vector index for deserialization range
        @param to - end vector index for deserialization range
        @param clear_sv - if true clears the target vector

    */
    void deserialize_range(SV& sv, const unsigned char* buf,
                           size_type from, size_type to,
                           bool clear_sv = true);

    /*!
        Better use deserialize_range()
        @sa deserialize_range
    */
    void deserialize(SV& sv, const unsigned char* buf,
                     size_type from, size_type to)
    {
        deserialize_range(sv, buf, from, to);
    }

    /*!
        Deserialize sparse vector using address mask vector
        Address mask defines (by set bits) which vector elements to be extracted
        from the compressed BLOB

        @param sv - [out] target sparse vector to populate
        @param buf - source memory pointer
        @param mask_bv - AND mask bit-vector (address gather vector)
    */
    void deserialize(SV& sv,
                     const unsigned char* buf,
                     const bvector_type& mask_bv);

private:
    bm::sparse_vector_deserializer<sparse_vector_type> exponentDeserializer_;   ///!< Exponents deserializer
    bm::sparse_vector_deserializer<sparse_vector_type> mantissaDeserializer_;   ///!< Mantissas deserializer
};

//---------------------------------------------------------------------

/*!
    \brief Serialize sparse vector float into a memory buffer(s) structure
 
    \param sv         - sparse vector float to serialize
    \param sv_layout  - buffer structure to keep the result
    \param temp_block - temporary buffer 
                        (allocate with BM_DECLARE_TEMP_BLOCK(x) for speed)
    
    \ingroup svserial
    
    @sa serialization_flags
    @sa sparse_vector_deserializer
*/
template<class SV>
void sparse_vector_float_serialize(
                const SV&                              sv,
                sparse_vector_float_serial_layout<SV>& sv_layout,
                bm::word_t*                            temp_block = 0)
{
    (void)temp_block;
    bm::sparse_vector_float_serializer<SV> sv_serializer;
//  sv_serializer.enable_xor_compression();
    sv_serializer.serialize(sv, sv_layout);
}

//---------------------------------------------------------------------

/*!
    \brief Deserialize sparse vector
    \param sv         - target sparse vector
    \param buf        - source memory buffer
    \param temp_block - temporary block buffer to avoid re-allocations
 
    \return 0  (error processing via std::logic_error)
 
    \ingroup svserial
    @sa sparse_vector_deserializer
*/
template<class SV>
int sparse_vector_float_deserialize(SV& sv,
                              const unsigned char* buf,
                              bm::word_t* temp_block=0)
{
    (void)temp_block;
    bm::sparse_vector_float_deserializer<SV> sv_deserializer;
    sv_deserializer.deserialize(sv, buf);
    return 0;
}

//---------------------------------------------------------------------
//sparse_vec_float_serial_layout methods

template<class SV>
sparse_vector_float_serial_layout<SV>::sparse_vector_float_serial_layout(const sparse_vector_float_serial_layout<SV>&)
{}

//---------------------------------------------------------------------

template<class SV>
unsigned char* sparse_vector_float_serial_layout<SV>::reserve(size_t capacity)
{
    if (capacity > capacity_)
    {
            unsigned char* new_buf = (unsigned char*)realloc(buf_, capacity);
            if (new_buf)
            {
                buf_ = new_buf;
                capacity_ = capacity;
            }
        }
    return buf_;
}

//---------------------------------------------------------------------

template<class SV>
size_t sparse_vector_float_serial_layout<SV>::size()
{
    return size_;
}

//---------------------------------------------------------------------

template<class SV>
void sparse_vector_float_serial_layout<SV>::resize(size_t ssize)
{
    reserve(ssize);
    size_ = ssize;
}

//---------------------------------------------------------------------

template<class SV>
void sparse_vector_float_serial_layout<SV>::shrink(size_t ssize)
{
    if (ssize < size_)
        size_ = ssize;
}

//---------------------------------------------------------------------

template<class SV>
size_t sparse_vector_float_serial_layout<SV>::capacity() const BMNOEXCEPT
{
    return capacity_;
}

//---------------------------------------------------------------------

template<class SV>
void sparse_vector_float_serial_layout<SV>::freemem() BMNOEXCEPT
{
    if (buf_)
    {
        ::free(buf_);
        buf_      = nullptr;
        capacity_ = 0;
        size_     = 0;
    }
}

//---------------------------------------------------------------------
//sparse_vec_float_serializer methods

template<class SV>
void sparse_vector_float_serializer<SV>::set_bookmarks(bool enable, unsigned bm_interval) BMNOEXCEPT
{
    signSerializer_.set_bookmarks(enable, bm_interval);
    exponentSerializer_.set_bookmarks(enable, bm_interval);
    mantissaSerializer_.set_bookmarks(enable, bm_interval);
}

template<class SV>
void sparse_vector_float_serializer<SV>::set_xor_ref(bool is_enabled) BMNOEXCEPT
{
    exponentSerializer_.set_xor_ref(is_enabled);
    mantissaSerializer_.set_xor_ref(is_enabled);
}

template<class SV>
void sparse_vector_float_serializer<SV>::set_xor_ref(const bv_ref_vector_type* bv_ref_ptr) BMNOEXCEPT
{
    exponentSerializer_.set_xor_ref(bv_ref_ptr);
    mantissaSerializer_.set_xor_ref(bv_ref_ptr);
}

template<class SV>
void sparse_vector_float_serializer<SV>::compute_sim_model(xor_sim_model_type& sim_model,
                                                            const bv_ref_vector_type& ref_vect,
                                                            const bm::xor_sim_params& params)
{
    exponentSerializer_.compute_sim_model(sim_model, ref_vect, params);
}

template<class SV>
void sparse_vector_float_serializer<SV>::set_sim_model(const xor_sim_model_type* sim_model) BMNOEXCEPT
{
    signSerializer_.set_sim_model(sim_model);
    exponentSerializer_.set_sim_model(sim_model);
    mantissaSerializer_.set_sim_model(sim_model);
}

template<class SV>
bool sparse_vector_float_serializer<SV>::is_xor_ref() const BMNOEXCEPT
{
    return mantissaSerializer_.is_xor_ref();
}

template<class SV>
void sparse_vector_float_serializer<SV>::serialize(const SV&                        sv,
                                                    sparse_vector_float_serial_layout<SV>& sv_layout)
{
    typedef typename SV::sparse_vector_u sparse_vector_u_type;

    typename serializer<typename SV::bvector_type>::buffer signBufTemp;
    sparse_vector_serial_layout<sparse_vector_u_type> expLayTemp;
    sparse_vector_serial_layout<sparse_vector_u_type> mantLayTemp;

    signSerializer_.serialize(sv.signs_, signBufTemp);
    exponentSerializer_.serialize(sv.exponents_, expLayTemp);
    mantissaSerializer_.serialize(sv.mantissas_, mantLayTemp);

    size_t sign_size_ = signBufTemp.size();
    size_t exp_size_ = expLayTemp.size();
    size_t mant_size_ = mantLayTemp.size();

    size_t header_size = 3 + sizeof(size_t) * 3;
    size_t totalSize = header_size + signBufTemp.size() + expLayTemp.size() + mantLayTemp.size();

    sv_layout.resize(totalSize);
    unsigned char* dest = sv_layout.buf_;

    std::memcpy(dest, "bf0", 3);
    dest += 3;

    std::memcpy(dest, &sign_size_, sizeof(size_t));
    dest += sizeof(size_t);
    std::memcpy(dest, &exp_size_, sizeof(size_t));
    dest += sizeof(size_t);
    std::memcpy(dest, &mant_size_, sizeof(size_t));
    dest += sizeof(size_t);

    std::memcpy(dest, signBufTemp.data(), sign_size_);
    dest += sign_size_;
    std::memcpy(dest, expLayTemp.buf(), exp_size_);
    dest += exp_size_;
    std::memcpy(dest, mantLayTemp.buf(), mant_size_);

}

//---------------------------------------------------------------------
//sparse_vec_float_deserializer methods

template<class SV>
sparse_vector_float_deserializer<SV>::sparse_vector_float_deserializer()
: exponentDeserializer_(), mantissaDeserializer_()
{}

//---------------------------------------------------------------------

template<class SV>
sparse_vector_float_deserializer<SV>::~sparse_vector_float_deserializer()
{}

//---------------------------------------------------------------------

template<class SV>
void sparse_vector_float_deserializer<SV>::set_finalization(bm::finalization is_final)
{
    exponentDeserializer_.set_finalization(is_final);
    mantissaDeserializer_.set_finalization(is_final);
}

//---------------------------------------------------------------------

template<class SV>
void sparse_vector_float_deserializer<SV>::set_xor_ref(bv_ref_vector_type* bv_ref_ptr)
{
    exponentDeserializer_.set_xor_ref(bv_ref_ptr);
    mantissaDeserializer_.set_xor_ref(bv_ref_ptr);
}

//---------------------------------------------------------------------

template<class SV>
void sparse_vector_float_deserializer<SV>::deserialize(SV& sv,
                                                        const unsigned char* buf,
                                                        bool clear_sv )
{
    if (clear_sv)
        sv.clear();

    size_t sign_size, exp_size, mant_size;
    
    buf+=3;
    
    std::memcpy(&sign_size, buf, sizeof(size_t)); buf += sizeof(size_t);
    std::memcpy(&exp_size,  buf, sizeof(size_t)); buf += sizeof(size_t);
    std::memcpy(&mant_size, buf, sizeof(size_t)); buf += sizeof(size_t);

    bm::deserialize(sv.signs_, buf);
    buf += sign_size;

    exponentDeserializer_.deserialize(sv.exponents_, buf, clear_sv);
    buf += exp_size;

    mantissaDeserializer_.deserialize(sv.mantissas_, buf, clear_sv);
}

//---------------------------------------------------------------------

template<class SV>
void sparse_vector_float_deserializer<SV>::deserialize_range(SV& sv, const unsigned char* buf,
                                                            size_type from, size_type to,
                                                            bool clear_sv)
{
    if (clear_sv)
        sv.clear();
    
    const unsigned char* ptr = buf;
    
    ptr+=3;

    size_t sign_size, exp_size, mant_size;
    std::memcpy(&sign_size, ptr, sizeof(size_t)); 
    ptr += sizeof(size_t);
    std::memcpy(&exp_size,  ptr, sizeof(size_t)); 
    ptr += sizeof(size_t);
    std::memcpy(&mant_size, ptr, sizeof(size_t)); 
    ptr += sizeof(size_t);
    
    bm::deserialize_range(sv.signs_,    ptr,              from, to);
    ptr += sign_size;
    exponentDeserializer_.deserialize_range(sv.exponents_, ptr, from, to);
    ptr += exp_size;
    mantissaDeserializer_.deserialize_range(sv.mantissas_, ptr, from, to);
}

//---------------------------------------------------------------------

template<class SV>
void sparse_vector_float_deserializer<SV>::deserialize(SV& sv,
                                                        const unsigned char* buf,
                                                        const bvector_type& mask_bv)
{
    sv.clear();

    const unsigned char* ptr = buf;
    ptr+=3;

    size_t sign_size, exp_size, mant_size;
    std::memcpy(&sign_size, ptr, sizeof(size_t)); 
    ptr += sizeof(size_t);
    std::memcpy(&exp_size,  ptr, sizeof(size_t)); 
    ptr += sizeof(size_t);
    std::memcpy(&mant_size, ptr, sizeof(size_t)); 
    ptr += sizeof(size_t);

    bm::deserialize(sv.signs_, ptr);
    sv.signs_ &= mask_bv;
    ptr += sign_size;
    exponentDeserializer_.deserialize(sv.exponents_, ptr, mask_bv);
    ptr += exp_size;
    mantissaDeserializer_.deserialize(sv.mantissas_, ptr, mask_bv);
}

//---------------------------------------------------------------------


}//namespace bm
#endif
