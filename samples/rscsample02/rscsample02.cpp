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

/** \example rscsample02.cpp
  Example of how to partially de-serialize bm::rsc_sparse_vector<> container
 
  \sa bm::rsc_sparse_vector
  \sa bm::rsc_sparse_vector<>::back_insert_iterator

*/

/*! \file rscsample02.cpp
    \brief Example: rsc_sparse_vector<> selective and range de-serialization

    @sa strsvsample05.cpp
*/

#include <iostream>
#include <vector>

#include "bm.h"
#include "bmsparsevec.h"
#include "bmsparsevec_compr.h"
#include "bmsparsevec_serial.h"
#include "bmundef.h" /* clear the pre-proc defines from BM */

using namespace std;

/// Print sparse vector content
template<typename SV> void PrintSV(const SV& sv)
{
    auto sz = sv.size();
    cout << "size() = " << sz << " : ";

    for (typename SV::size_type i = 0; i < sz; ++i)
    {
        if (sv.is_null(i))
            cout << "NULL";
        else
            cout << sv.get(i);
        cout << ", ";
    }
    cout << endl;
}

typedef bm::sparse_vector<unsigned, bm::bvector<> >         sparse_vector_u32;
typedef bm::rsc_sparse_vector<unsigned, sparse_vector_u32 > rsc_sparse_vector_u32;


int main(void)
{
    try
    {
        rsc_sparse_vector_u32 csv1;
        rsc_sparse_vector_u32 csv2;

        // add sample data using back insert iterator
        //
        {
            auto bit = csv1.get_back_inserter();
            bit = 10;
            bit = 11;
            bit.add_null();
            bit = 13;
            bit = 14;
            bit.add_null(2);
            bit = 256;
            bit.flush();

            csv1.sync();
        }
        PrintSV(csv1); // size() = 8 : 10, 11, NULL, 13, 14, NULL, NULL, 256,

        // serialize sparse vector
        bm::sparse_vector_serial_layout<rsc_sparse_vector_u32> sv_lay;
        {
            BM_DECLARE_TEMP_BLOCK(tb)
            // optimize memory allocation of sparse vector
            csv1.optimize(tb);

            // configure serializer, set bookmarking option for faster
            // range deserialization
            //
            bm::sparse_vector_serializer<rsc_sparse_vector_u32> sv_serializer;
            sv_serializer.set_bookmarks(true, 64);

            sv_serializer.serialize(csv1, sv_lay); // serialization
        }

        // get serialization pointer
        const unsigned char* buf = sv_lay.buf();

        // instantiate the deserializer object to do all multiple
        // deserialization operations (faster)
        //
        bm::sparse_vector_deserializer<rsc_sparse_vector_u32> sv_deserial;

        // Case 1:
        // Here we define de-serialization indexes as a bit-vector (mask)
        // so deserializer gathers only elements selected by the bit-vector
        // all vector elements not set in the gather vector will be 0 (or NULL)
        //
        {
            rsc_sparse_vector_u32::bvector_type bv_mask;
            bv_mask.set(0);
            bv_mask.set(1);
            sv_deserial.deserialize(csv2, buf, bv_mask);
        }
        PrintSV(csv2); // size() = 8 : 10, 11, NULL, 0, 0, NULL, NULL, 0,

        // Case 2:
        // Here we define de-serialization indexes as a closed range [1..4]
        // It is an equivalent of setting selection bits in the mask vector
        // but defines the intent more clearly and works faster
        {
            sv_deserial.deserialize_range(csv2, buf, 1, 4);
        }
        PrintSV(csv2); // size() = 8 : 0, 11, NULL, 13, 14, NULL, NULL, 0,

    }
    catch(std::exception& ex)
    {
        std::cerr << ex.what() << std::endl;
        return 1;
    }
    
    

    return 0;
}

