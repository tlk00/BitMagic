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

template<class SV>
struct sparse_vector_serial_layout
{
    sparse_vector_serial_layout()
    : buffer(0), serialized_size(0)
    {
    }
    
    ~sparse_vector_serial_layout()
    {
        delete [] buffer;
    }
private:
    sparse_vector_serial_layout(const sparse_vector_serial_layout&);
    void operator=(const sparse_vector_serial_layout&);
public:
    unsigned char* buffer;                       ///< serialization buffer
    unsigned       serialized_size;              ///< serialized size
    
    unsigned char* plain_ptrs[SV::value_type*8]; ///< pointers on serialized bit-palins
    unsigned plane_size[[SV::value_type*8];      ///< serialized plain size
};


template<class SV>
void sparse_vector_serialize(
                const SV&                        sv,
                sparse_vector_serial_layout<SV>& sv_layout,
                unsigned                         bv_serialization_flags = 0)
             
{
}

    
} // namespace bm

#include "bmundef.h"

#endif
