
#include <memory.h>

#include "bm.h"
#include "bmsparsevec.h"
#include "bmsparsevec_serial.h"
#include "bmsparsevec_float.h"

typedef bm::sparse_vector<unsigned int, bm::bvector<>> sparse_vector_u32;

namespace bm{

class sparse_vector_float_serialized{
 
public:
    sparse_vector_float_serialized();
    sparse_vector_float_serialized(const sparse_vector_float_serialized&);
    ~sparse_vector_float_serialized();

    void serialize(sparse_vector_float& svf);
    void deserialize(sparse_vector_float& svf);
        
private:
    serializer<bm::bvector<>>::buffer sign_buf;
    sparse_vector_serial_layout<sparse_vector_u32> exp_lay;
    sparse_vector_serial_layout<sparse_vector_u32> mant_lay;
    
};
//---------------------------------------------------------------------

sparse_vector_float_serialized::sparse_vector_float_serialized(){}
sparse_vector_float_serialized::sparse_vector_float_serialized(const sparse_vector_float_serialized&){}
sparse_vector_float_serialized::~sparse_vector_float_serialized(){}

//---------------------------------------------------------------------

void sparse_vector_float_serialized::serialize(sparse_vector_float& svf){
    serializer<bm::bvector<>> bvs;
    bvs.serialize(svf.signs, sign_buf);
    sparse_vector_serialize(svf.exponents, exp_lay);
    sparse_vector_serialize(svf.mantissas, mant_lay);
}

//---------------------------------------------------------------------

void sparse_vector_float_serialized::deserialize(sparse_vector_float& svf){
    const unsigned char* buf = sign_buf.buf();
    bm::deserialize(svf.signs, buf);

    bm::sparse_vector_deserialize(svf.exponents, exp_lay.buf());

    bm::sparse_vector_deserialize(svf.mantissas, mant_lay.buf());
}

//---------------------------------------------------------------------

}//namespace bm

