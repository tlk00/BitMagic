/*
Copyright(c) 2002-2023 Anatoliy Kuznetsov(anatoliy_kuznetsov at yahoo.com)

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

/** \example bmbsample01.cpp
  Example of how to serialize bm::basic_bmatrix<> template class
 
  \sa bm::basic_bmatrix
  \sa bm::sparse_vector_serializer
  \sa bm::sparse_vector_deserializer
  \sa sparse_vector_serial_layout

*/

/*! \file bmbsample01.cpp
    \brief Example: basic_bmatrix<> serialization

*/

#include <iostream>
#include <vector>
#include <assert.h>

#include "bm.h"
#include "bmbmatrix.h"

#include "bmsparsevec.h"
#include "bmsparsevec_serial.h"
#include "bmundef.h" /* clear the pre-proc defines from BM */

using namespace std;

typedef bm::basic_bmatrix< bm::bvector<> >       bmatr_32;
typedef bm::sparse_vector_serializer<bmatr_32>   bmb_serializer_type;
typedef bm::sparse_vector_deserializer<bmatr_32> bmb_deserializer_type;


int main(void)
{
    try
    {
        bmatr_32 bmatr0(0); // create sample row major bit-matrix

        // fill in some random content (add random rows, set random bits)
        //
        for (unsigned i = 0; i < 1000; i+= (unsigned) rand()%25)
        {
            bmatr_32::bvector_type* bv = bmatr0.construct_row(i);
            for (unsigned j = 0; j < 150000; j += (unsigned)rand()%3)
                bv->set(j);
        } // for i
        bmatr0.optimize(); // optimize all bit-bectors


        // Serialization
        //

        // memory layout object for bmatr0
        //
        bm::sparse_vector_serial_layout<bmatr_32> bmatr_lay;

        {
            bmb_serializer_type bmbser;     // serializer class
            // set bookmarks for faster range deserialization
            bmbser.set_bookmarks(true, 16);
            // enable XOR filter
            //   in this case likely useless, since matrix is
            //   populated with the random data
            //
            bmbser.enable_xor_compression();

            // run serialization, layout now contains BLOB
            bmbser.serialize(bmatr0, bmatr_lay);

            cout << "Serialized size = " << bmatr_lay.size() << endl;
        }

        // De-serialization 1 (regular case)
        //
        {
        bmb_deserializer_type bmatr_deserial; // deserializer instance

            const unsigned char* buf = bmatr_lay.buf();
            bmatr_32 bmatr1(0);

            bmatr_deserial.deserialize(bmatr1, buf);

            bool eq = bmatr0.equal(bmatr1);
            assert(eq); (void)eq;

            assert(!bmatr0.is_ro()); // restored matrix is writable
        }


        // De-serialization 2 (read-only)
        //
        {
        bmb_deserializer_type bmatr_deserial_ro; // deserializer instance
        // configure for read-only finalization
        // (read-only uses less memroty and less heap fragmentation)
        //
        bmatr_deserial_ro.set_finalization(bm::finalization::READONLY);

            const unsigned char* buf = bmatr_lay.buf();
            bmatr_32 bmatr1(0);

            bmatr_deserial_ro.deserialize(bmatr1, buf);

            bool eq = bmatr0.equal(bmatr1);
            assert(eq); (void)eq;

            assert(bmatr1.is_ro()); // restored matrix is NOT writable
        }

        // De-serialization 3 (read-only, get only specific range)
        //
        {
        bmb_deserializer_type bmatr_deserial_ro; // deserializer instance
        // configure for read-only finalization
        // (read-only uses less memroty and less heap fragmentation)
        //
        bmatr_deserial_ro.set_finalization(bm::finalization::READONLY);

            const unsigned char* buf = bmatr_lay.buf();
            bmatr_32 bmatr1(0);

            bmatr_deserial_ro.deserialize_range(bmatr1, buf, 10, 100);

            bool eq = bmatr0.equal(bmatr1); (void)eq;
            assert(!eq); // not the same matrix
            assert(bmatr1.is_ro()); // restored matrix is NOT writable
        }


        // De-serialization 3 (read-only, get only specific range)
        //
        {
        bmb_deserializer_type bmatr_deserial_ro; // deserializer instance
        bmatr_deserial_ro.set_finalization(bm::finalization::READONLY);

        bmatr_32::bvector_type bv_mask; // AND vector [10..100]...[200..210]
        bv_mask.set_range(10, 100);
        bv_mask.set_range(200, 210);

        const unsigned char* buf = bmatr_lay.buf();
        bmatr_32 bmatr1(0);
        bmatr_deserial_ro.deserialize(bmatr1, buf, bv_mask);

        assert(bmatr1.is_ro());


        // modify the original matrix to compare
        //
        for (unsigned i = 0; i < bmatr0.rows(); ++i)
        {
            bmatr_32::bvector_type* bv = bmatr0.get_row(i);
            if (bv)
                bv->bit_and(bv_mask);
        } // for i
        bool eq = bmatr0.equal(bmatr1); (void)eq;
        assert(eq); // not the same matrix


        }

    }
    catch(std::exception& ex)
    {
        std::cerr << ex.what() << std::endl;
        return 1;
    }

    return 0;
}

