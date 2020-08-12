/*
Copyright(c) 2002-2019 Anatoliy Kuznetsov(anatoliy_kuznetsov at yahoo.com)

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

/** \example svsample08.cpp
  Example of how to serialize bm::sparse_vector<> container
 
  \sa bm::sparse_vector
  \sa bm::sparse_vector_deserializer
  \sa bm::sparse_vector_serial_layout
  \sa bm::sparse_vector_serializer

*/

/*! \file svsample08.cpp
    \brief Example: sparse_vector<> selective de-serialization (gather)
    and range deserialization

    @sa strsvsample05.cpp
*/

#include <iostream>
#include <vector>

#include "bm.h"
#include "bmsparsevec.h"
#include "bmsparsevec_serial.h"
#include "bmundef.h" /* clear the pre-proc defines from BM */

using namespace std;

/// Print sparse vector content
template<typename SV> void PrintSV(const SV& sv)
{
    cout << "size() = " << sv.size() << " : ";
    auto it = sv.begin(); auto it_end = sv.end();
    for (; it != it_end; ++it)
    {
        if (it.is_null())
            cout << "NULL";
        else
            cout << *it;
        cout << ", ";
    }
    cout << endl;
}

typedef bm::sparse_vector<unsigned, bm::bvector<> > svector_u32;

int main(void)
{
    try
    {
        svector_u32 sv1(bm::use_null);
        svector_u32 sv2(bm::use_null);

        // add sample data using back insert iterator
        //
        {
            auto bit = sv1.get_back_inserter();
            bit = 10;
            bit = 11;
            bit.add_null();
            bit = 13;
            bit = 14;
            bit.add_null(2);
            bit = 256;
            bit.flush();
        }
        PrintSV(sv1); // size() = 8 : 10, 11, NULL, 13, 14, NULL, NULL, 256,

        // serialize sparse vector
        bm::sparse_vector_serial_layout<svector_u32> sv_lay;
        {
            BM_DECLARE_TEMP_BLOCK(tb)
            // optimize memory allocation of sparse vector
            sv1.optimize(tb);


            // construct a serializer utility class, setup serialization parameters
            //
            // please note, use of "set_bookmarks()" to enable fast range
            // deserialization. Bookmarks somewhat increase the BLOB size but allow
            // more effeiciently skip parts which we would not need (paging) and
            // avoid decompression of blocks we would never need
            //
            // This example sets "64" as a bookmarks parameter, but you have to
            // experiment with what works for you, between 4 and 512
            //
            // Each block corresponds to 64K vector element
            // making bookmarks after each block does not make much sense
            // because decode is reasonably fast and some residual throw away
            // is usually ok.
            //

            bm::sparse_vector_serializer<svector_u32> sv_serializer;
            sv_serializer.set_bookmarks(true, 64);

            // serialization, creates a compressed BLOB
            sv_serializer.serialize(sv1, sv_lay);
        }

        // get serialization pointer
        const unsigned char* buf = sv_lay.buf();

        // instantiate the deserializer object to do all multiple
        // deserialization operations (faster)
        //
        bm::sparse_vector_deserializer<svector_u32> sv_deserial;

        // Case 1:
        // Here we define de-serialization indexes as a bit-vector (mask)
        // so deserializer gathers only elements selected by the bit-vector
        // all vector elements not set in the gather vector will be 0 (or NULL)
        //
        {
            svector_u32::bvector_type bv_mask;
            bv_mask.set(0);
            bv_mask.set(1);
            sv_deserial.deserialize(sv2, buf, bv_mask);
        }
        PrintSV(sv2); // size() = 8 : 10, 11, NULL, 0, 0, NULL, NULL, 0,

        // Case 2:
        // Here we define de-serialization indexes as a closed range [1..4]
        // It is an equivalent of setting selection bits in the mask vector
        // but defines the intent more clearly and works faster
        {
            sv_deserial.deserialize_range(sv2, buf, 1, 4);
        }
        PrintSV(sv2); // size() = 8 : 0, 11, NULL, 13, 14, NULL, NULL, 0,

    }
    catch(std::exception& ex)
    {
        std::cerr << ex.what() << std::endl;
        return 1;
    }
    
    

    return 0;
}

