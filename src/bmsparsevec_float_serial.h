#ifndef BM_SPARSE_VEC_FLOAT_SERIAL_INCLUDED
#define BM_SPARSE_VEC_FLOAT_SERIAL_INCLUDED

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

/*! \file sparse_vector_float_serial.h
    \brief Serialization for sparse_vector_float
*/

#include <memory.h>

#include "bm.h"

#include "bmsparsevec_float.h"
#include "bmsparsevec.h"
#include "bmsparsevec_serial.h"

namespace bm{

template<class SV>
class sparse_vector_float_serial_layout{
 
public:
    typedef typename SV::value_type   value_type;
    typedef typename SV::bvector_type bvector_type;
    typedef bm::sparse_vector<unsigned int, bvector_type> sparse_vector_u32;

    
    sparse_vector_float_serial_layout() BMNOEXCEPT {}
    ~sparse_vector_float_serial_layout() {}
    sparse_vector_float_serial_layout(const sparse_vector_float_serial_layout&);

    
    /*!
        \brief resize capacity
        \param capacity - new capacity
        \return new buffer or 0 if failed
    */
    unsigned char* reserve(size_t capacity)

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
    
private:
    serializer<bm::bvector<>>::buffer sign_buf;
    sparse_vector_serial_layout<sparse_vector_u32> exp_lay;
    sparse_vector_serial_layout<sparse_vector_u32> mant_lay;
    
};


template<typename SV>
class sparse_vector_float_serializer
{
public:
    typedef typename SV::bvector_type       bvector_type;
    typedef const bvector_type*             bvector_type_const_ptr;
    typedef bvector_type*                   bvector_type_ptr;
    typedef typename SV::value_type         value_type;
    typedef typename SV::size_type          size_type;
    typedef typename bvector_type::allocator_type       alloc_type;
    typedef typename alloc_type::allocator_pool_type    allocator_pool_type;
    typedef typename
    bm::serializer<bvector_type>::bv_ref_vector_type bv_ref_vector_type;
    typedef typename
    bm::serializer<bvector_type>::xor_sim_model_type xor_sim_model_type;

public:
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
        Enables XOR reference compression for the sparse vector.
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

    ///@}

    /*! @name Serialization                                     */
    ///@{

    /*!
        \brief Serialize sparse vector into a memory buffer(s) structure
     
        \param sv                 - sparse vector to serialize
        \param sv_layout  - buffer structure to keep the result
        as defined in bm::serialization_flags
    */
    void serialize(const SV&                        sv,
                   sparse_vector_float_serial_layout<SV>& sv_layout);


protected:
    bm::serializer<bvector_type>     signSerializer_;
    bm::sparse_vector_serializer<bvector_type> exponentSerializer_;
    bm::sparse_vector_serializer<bvector_type> mantissaSerializer_;
    
};

/**
    sparse vector de-serializer

*/
template<typename SV>
class sparse_vector_deserializer
{

};


/*!
    \brief Serialize sparse vector into a memory buffer(s) structure
 
    \param sv         - sparse vector to serialize
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
//    sv_serializer.enable_xor_compression();
    sv_serializer.serialize(sv, sv_layout);
}


//---------------------------------------------------------------------
//sparse_vec_float_serial_layout methods

template<class SV>
sparse_vector_float_serial_layout<SV>::sparse_vector_float_serial_layout(const sparse_vector_float_serial_layout<SV>&)
{}

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
    //todo, unsure how this should work currently
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
    signSerializer_.serialize(sv.signs, sv_layout.sign_buf);
    exponentSerializer_.serialize(sv.floats, sv_layout.exp_lay);
    mantissaSerializer_.serialize(sv.mantissas, sv_layout.mant_lay);
}

//---------------------------------------------------------------------
//sparse_vec_float_deserializer methods


//---------------------------------------------------------------------


}//namespace bm
#endif
