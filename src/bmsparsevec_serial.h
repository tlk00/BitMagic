#ifndef BMSPARSEVEC_SERIAL__H__INCLUDED__
#define BMSPARSEVEC_SERIAL__H__INCLUDED__
/*
Copyright(c) 2017 Anatoliy Kuznetsov(anatoliy_kuznetsov at yahoo.com)

Permission is hereby granted, free of charge, to any person 
obtaining a copy of this software and associated documentation 
files (the "Software"), to deal in the Software without restriction, 
including without limitation the rights to use, copy, modify, merge, 
publish, distribute, sublicense, and/or sell copies of the Software, 
and to permit persons to whom the Software is furnished to do so, 
subject to the following conditions:

The above copyright notice and this permission notice shall be included 
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, 
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, 
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR 
OTHER DEALINGS IN THE SOFTWARE.

For more information please visit:  http://bmagic.sourceforge.net

*/

#include "bmdef.h"
#include "bmsparsevec.h"

namespace bm
{

/*!
    \brief layout class for serialization buffer structure
    
    Class keeps a memory block sized for the target vector
 
*/
template<class SV>
struct sparse_vector_serial_layout
{
    sparse_vector_serial_layout()
    : buffer_(0), capacity_(0), serialized_size_(0)
    {
    }
    
    ~sparse_vector_serial_layout()
    {
        if (buffer_)
            free(buffer_);
    }
    
    /// return current serialized size
    size_t  size() const { return serialized_size_; }
    
    /// return serialization buffer capacity
    size_t  capacity() const { return capacity_; }
    
private:
    sparse_vector_serial_layout(const sparse_vector_serial_layout&);
    void operator=(const sparse_vector_serial_layout&);
public:
    unsigned char* buffer_;                       ///< serialization buffer
    size_t         capacity_;                      ///< buffer capacity
    size_t         serialized_size_;               ///< serialized size
    
    unsigned char* plain_ptrs_[SV::value_type*8]; ///< pointers on serialized bit-palins
    unsigned plane_size_[SV::value_type*8];       ///< serialized plain size
};

// -------------------------------------------------------------------------

template<class SV>
class sparse_vec_serializer
{
public:
    sparse_vec_serializer();
    
    
private:
    sparse_vec_serializer(const sparse_vec_serializer&);
    sparse_vec_serializer& operator=(const sparse_vec_serializer&);

};

// -------------------------------------------------------------------------

template<class SV>
sparse_vec_serializer<SV>::sparse_vec_serializer()
{
}

// -------------------------------------------------------------------------


template<class SV>
void sparse_vector_serialize(
                const SV&                        sv,
                sparse_vector_serial_layout<SV>& sv_layout,
                unsigned                         bv_serialization_flags = 0)
{
}

// -------------------------------------------------------------------------

    
} // namespace bm

#include "bmundef.h"

#endif
