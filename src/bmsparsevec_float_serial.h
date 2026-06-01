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

template<class BV>
class sparse_vector_float_serialized{
 
public:
    typedef bm::sparse_vector<unsigned int, BV> sparse_vector_u32;

    
    sparse_vector_float_serialized();
    sparse_vector_float_serialized(const sparse_vector_float_serialized&);
    ~sparse_vector_float_serialized();

    /*!
        \brief Serializes a given sparse_vector_float
        \param svf   sparse_vector_float to serialize
    */
    void serialize(sparse_vector_float<BV>& svf);

    /*!
        \brief Deserializes a sparse_vector_float
        \param svf   sparse_vector_float to deserialize into
    */
    void deserialize(sparse_vector_float<BV>& svf);

    /// return current serialized size
    size_t size();
        
private:
    bm::bvector<>::size_type size_;
    serializer<bm::bvector<>>::buffer sign_buf;
    sparse_vector_serial_layout<sparse_vector_u32> exp_lay;
    sparse_vector_serial_layout<sparse_vector_u32> mant_lay;
    
};
//---------------------------------------------------------------------

template<class BV>
sparse_vector_float_serialized<BV>::sparse_vector_float_serialized(){}

template<class BV>
sparse_vector_float_serialized<BV>::sparse_vector_float_serialized(const sparse_vector_float_serialized<BV>&){}

template<class BV>
sparse_vector_float_serialized<BV>::~sparse_vector_float_serialized(){}

//---------------------------------------------------------------------

template<class BV>
void sparse_vector_float_serialized<BV>::serialize(sparse_vector_float<BV>& svf){
    size_ = svf.size_;
    serializer<bm::bvector<>> bvs;
    bvs.serialize(svf.signs, sign_buf);
    sparse_vector_serialize(svf.exponents, exp_lay);
    sparse_vector_serialize(svf.mantissas, mant_lay);
    
}

//---------------------------------------------------------------------

template<class BV>
void sparse_vector_float_serialized<BV>::deserialize(sparse_vector_float<BV>& svf){
    const unsigned char* buf = sign_buf.buf();
    bm::deserialize(svf.signs, buf);

    bm::sparse_vector_deserialize(svf.exponents, exp_lay.buf());

    bm::sparse_vector_deserialize(svf.mantissas, mant_lay.buf());
    svf.size_ = size_;
}

//---------------------------------------------------------------------

template<class BV>
size_t sparse_vector_float_serialized<BV>::size(){
    return sign_buf.size() + exp_lay.size() + mant_lay.size();
}

}//namespace bm
#endif
