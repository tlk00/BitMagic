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

/** \example sample14.cpp
 Exmaple demonstrates bitvector serialization/deserialization and set-operations on
 searialized BLOBs
 
  \sa bm::serializer
  \sa bm::deserialize
*/

/*! \file sample14.cpp
    \brief Example: bvector<> set operations on serialized/compressed BLOBs
*/


#include <stdlib.h>
#include <iostream>
#include <vector>

#include "bm.h"
#include "bmserial.h"
#include "bmundef.h" /* clear the pre-proc defines from BM */

using namespace std;


const unsigned MAX_VALUE = 1000000;

static
void fill_bvector(bm::bvector<>* bv)
{
    for (unsigned i = 0; i < MAX_VALUE; ++i)
    {
        if ((rand() % 10))
        {
            bv->set(i);
        }
    }
}


int main(void)
{
    bm::operation_deserializer<bm::bvector<> > od;
    try
    {
        bm::bvector<>   bv1;
        bm::bvector<>   bv2;

        fill_bvector(&bv1);
        fill_bvector(&bv2);
        
        cout << "bv1 count = " << bv1.count() << endl;
        cout << "bv2 count = " << bv2.count() << endl;


        // Prepare a serializer class
        //  for best performance - create serilizer once and reuse it
        //
        BM_DECLARE_TEMP_BLOCK(tb)
        bm::serializer<bm::bvector<> > bvs(tb);
        bvs.set_compression_level(4);

        bm::bvector<>::statistics st1;
        bm::bvector<>::statistics st2;
        
        // compress bit-vectors and compute statistics
        // (later re-used in serialization)
        //
        bv1.optimize(tb, bm::bvector<>::opt_compress, &st1);
        bv1.optimize(tb, bm::bvector<>::opt_compress, &st2);

        // declare serialization buffers
        bm::serializer<bm::bvector<> >::buffer sbuf1;
        bm::serializer<bm::bvector<> >::buffer sbuf2;

        // perform serialization
        //
        bvs.serialize(bv1, sbuf1, &st1);
        bvs.serialize(bv2, sbuf2, &st2);


        // Serialized bvectors (sbuf1 and sbuf2) now ready to be
        // saved to a database, file or send over a network.
        // to simulate this we just copy content to std::vector<>
        //
        
        std::vector<unsigned char> vect1;
        std::vector<unsigned char> vect2;
        
        vect1.resize(sbuf1.size());
        vect2.resize(sbuf2.size());
        
        ::memcpy(vect1.data(), sbuf1.buf(), sbuf1.size());
        ::memcpy(vect2.data(), sbuf2.buf(), sbuf2.size());


        // Simple deserialization.
        //
        bm::bvector<>  bv3;

        // As a result of desrialization bv3 will contain all bits from
        // bv1 and bv3:
        //   bv3 = bv1 OR bv2

        bm::deserialize(bv3, sbuf1.buf());
        bm::deserialize(bv3, sbuf2.buf());

        bv3.optimize(tb);
        
        // A few examples of operation deserializations
        //  set algebraic operation over bit-vector and a BLOB
        //
        bm::bvector<>  bv4(bv3);
        
        // bv4 = (bv1 OR bv2) AND bv1
        // this must be equal to bv1 ?
        //
        od.deserialize(bv4, vect1.data(), bm::set_AND);
        
        cout << "bv4 count = " << bv4.count() << endl;
        bool eq = bv1.equal(bv4);
        if (!eq)
        {
            cerr << "Logical error detected!" << endl;
            return 1;
        }
        else
            cout << "bv4 is equal to bv1" << endl;
        
        bm::bvector<>  bv5(bv3);


        // if we need just set count, we can get it faster
        // via set_COUNT_ operations
        // use of COUNT operations does not materialize a target vector
        //
        // POPCNT((bv1 OR bv2) MINUS bv1)
        auto cnt_sub =
        od.deserialize(bv3, sbuf1.buf(), bm::set_COUNT_SUB_AB);
        cout << "minus count = " << cnt_sub << endl;

        
        // or we can actually perform the operation and get the full vector
        // bv5 = (bv1 OR bv2) MINUS bv1
        //
        od.deserialize(bv5, sbuf1.buf(), tb, bm::set_SUB);
        auto bv5_cnt = bv5.count();
        cout << "bv5 count = " << bv5_cnt << endl;

        if (cnt_sub != bv5_cnt)
        {
            cerr << "Logical error!" << endl;
            return 1;
        }
    }
    catch(std::exception& ex)
    {
        std::cerr << ex.what() << std::endl;
        return 1;
    }
    

    return 0;
}

