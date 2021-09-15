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

/** \example strsvsample05.cpp

  Example of how to use bm::str_sparse_vector<> - succinct container for
  bit-transposed string collections for deserialization of only select elements
  from the serialized BLOB
 
  \sa bm::str_sparse_vector
  \sa bm::sparse_vector_deserializer
  \sa bm::sparse_vector_serializer

*/

/*! \file strsvsample05.cpp
    \brief Example: str_sparse_vector<> gather deserialization example
 
    This example loads a range of a sparse vector from an STL container to save
    memory and improve deserialization performance
*/

#include <iostream>
#include <string>
#include <vector>
#include <assert.h>

#include "bm.h"
#include "bmstrsparsevec.h"
#include "bmsparsevec_serial.h"
#include "bmundef.h" /* clear the pre-proc defines from BM */

using namespace std;

typedef bm::bvector<> bvector_type;
typedef bm::str_sparse_vector<char, bvector_type, 5> str_sv_type;




int main(void)
{
    try
    {
       str_sv_type str_sv1;
       str_sv_type str_sv2;
       str_sv_type str_sv3;

       {
           str_sv_type str_sv0;
           // here we generate collection of k-mer (4-mer) strings
           // imitating a DNA sequence
           {
               auto bi = str_sv0.get_back_inserter();
               for (unsigned i = 0; i < 100000; ++i)
               {
                    bi = "ATGC";
                    bi = "GCTA";
                    bi = "GCAA";
                    bi = "TATA";
               } // for
           }
           str_sv1.remap_from(str_sv0); // SV1 now contains a remapped(smaller) copy of SV0
       }
       BM_DECLARE_TEMP_BLOCK(tb)
       str_sv1.optimize(tb);


        // calculate memory footprint
        //
        str_sv_type::statistics st;
        str_sv1.calc_stat(&st);

        cout << "Used memory: " << st.memory_used << std::endl;

        bm::sparse_vector_serial_layout<str_sv_type> sv_lay;

        // construct a serializer utility class, setup serialization parameters
        //
        // please note, use of "set_bookmarks()" to enable fast range
        // deserialization. Bookmarks somewhat increase the BLOB size but allow
        // more effeiciently skip parts which we would not need (paging) and
        // avoid decompression of blocks we would never need
        //
        // This example sets "128" as a bookmarks parameter, but you have to
        // experiment with what works for you, between 4 and 512
        //
        // Each block corresponds to 64K vector element
        // making bookmarks after each block does not make much sense
        // because decode is reasonably fast and some residual throw away
        // is usually ok.
        //
        bm::sparse_vector_serializer<str_sv_type> sv_serializer;
        sv_serializer.set_bookmarks(true, 128);


        // run str-vector serialization with compression
        //
        sv_serializer.serialize(str_sv1, sv_lay);

        const unsigned char* buf = sv_lay.buf();
        cout << "Serialized size = " << sv_lay.size() << endl;

        // instantiate deserializer utility class
        //
        bm::sparse_vector_deserializer<str_sv_type> sv_deserial;


        bvector_type::size_type from = 100000;
        bvector_type::size_type to = from + 65536;
        {
            // 1.
            // one way to deserialize is to provide a mask vector
            // specifying which sparse vector elements needs to be
            // decompressed from the BLOB
            // mask vector does not necessarily has to be just one range
            //
            bvector_type bv_mask;
            bv_mask.set_range(from, to);
            sv_deserial.deserialize(str_sv2, buf, bv_mask);


            // 2.
            // If it is just one range (common use case for paging)
            // it is faster and cleaner to use deserialize_range().
            // It will produce the same result as with (1) just faster.
            //
            sv_deserial.deserialize_range(str_sv3, buf, from, to);

            // run a quick comparison, that selected range matches values in
            // the container str_sv2, str_sv3
            //
            char s1[16]; char s2[16]; char s3[16];
            for (bvector_type::size_type j = from; j < to; ++j)
            {
                str_sv1.get(j, s1, sizeof(s1));
                str_sv2.get(j, s2, sizeof(s2));
                str_sv3.get(j, s3, sizeof(s3));

                int cmp;
                cmp = ::strcmp(s1, s2);
                assert(cmp==0);
                cmp = ::strcmp(s1, s3);
                assert(cmp==0); (void)cmp;

            } // for j
            cout << "Gather deserialization check OK" << endl;
        }

    }
    catch(std::exception& ex)
    {
        std::cerr << ex.what() << std::endl;
        return 1;
    }
    

    return 0;
}

