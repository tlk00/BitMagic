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

/** \example sample4.cpp
 Exmaple demonstrates bitvector serialization/deserialization.
 
For more information please visit:  http://bmagic.sourceforge.net

  \sa bm::serializer
  \sa bm::deserialize
*/

/*! \file sample4.cpp
    \brief Example: bvector<> serialization/deserialization.
*/
#include <stdlib.h>
#include <iostream>

#include "bm.h"
#include "bmserial.h"
#include "bmundef.h" /* clear the pre-proc defines from BM */

using namespace std;


// This exmaple demonstrates bitvector serialization/deserialization.


const unsigned MAX_VALUE = 1000000;

// This procedure creates very dense bitvector.
// The resulting set will consists mostly from ON (1) bits
// interrupted with small gaps of 0 bits.

static
void fill_bvector(bm::bvector<>* bv)
{
    for (unsigned i = 0; i < MAX_VALUE; ++i)
    {
        if (rand() % 2500)
        {
            bv->set_bit(i);
        }
    }
}

static
void print_statistics(const bm::bvector<>& bv)
{
    bm::bvector<>::statistics st;
    bv.calc_stat(&st);

    cout << "Bits count:" << bv.count() << endl;
    cout << "Bit blocks:" << st.bit_blocks << endl;
    cout << "GAP blocks:" << st.gap_blocks << endl;
    cout << "Memory used:"<< st.memory_used << endl;
    cout << "Max.serialize mem.:" << st.max_serialize_mem << endl << endl;;
}

// This is fairly low level serialization method, saves into a buffer
//
static
unsigned char* serialize_bvector(bm::serializer<bm::bvector<> >& bvs, 
                                 bm::bvector<>& bv)
{
    // It is good to optimize vector before serialization.
    
	// scratch memory block
    BM_DECLARE_TEMP_BLOCK(tb)
    bm::bvector<>::statistics st;
    bv.optimize(tb, bm::bvector<>::opt_compress, &st);

    cout << "Bits count:" << bv.count() << endl;
    cout << "Bit blocks:" << st.bit_blocks << endl;
    cout << "GAP blocks:" << st.gap_blocks << endl;
    cout << "Memory used:"<< st.memory_used << endl;
    cout << "Max.serialize mem.:" << st.max_serialize_mem << endl;

    // Allocate serialization buffer.
    unsigned char*  buf = new unsigned char[st.max_serialize_mem];

    // Serialization to memory.
    size_t len = bvs.serialize(bv, buf, st.max_serialize_mem);

    cout << "Serialized size:" << len << endl << endl;
    return buf;
}



int main(void)
{
    unsigned char* buf1 = 0;
    unsigned char* buf2 = 0;

    try
    {
        bm::bvector<>   bv1;
        bm::bvector<>   bv2;

        bv2.set_new_blocks_strat(bm::BM_GAP);  //  set DGAP compression mode ON

        fill_bvector(&bv1);
        fill_bvector(&bv2);

        // Prepare a serializer class
        //  for best performance it is best to create serilizer once and reuse it
        //  (saves a memory allocations). Do NOT use it concurrently.
        //
        bm::serializer<bm::bvector<> > bvs;

        // next settings provide lowest size
        bvs.byte_order_serialization(false);
        bvs.gap_length_serialization(false);

        // 1: serializes into a new allocated buffer (needs to be deallocated)
        buf1 = serialize_bvector(bvs, bv1);
        
        // 2: use serialization buffer class (automatic RAI, freed on destruction)
        bm::serializer<bm::bvector<> >::buffer sbuf;
        {
            bvs.serialize(bv2, sbuf);
            buf2 = sbuf.data();
            auto sz = sbuf.size();
            cout << "BV2 Serialized size:" << sz << endl;
        }

        // Serialized bvectors (buf1 and buf2) now ready to be
        // saved to a database, file or send over a network.

        // ...

        // Deserialization.

        bm::bvector<>  bv3;

        // As a result of desrialization bv3 will contain all bits from
        // bv1 and bv3:
        //   bv3 = bv1 OR bv2

        bm::deserialize(bv3, buf1);
        bm::deserialize(bv3, buf2);

        print_statistics(bv3);

        // After a complex operation we can try to optimize bv3.

        bv3.optimize();

        print_statistics(bv3);
        
        // control check using OR on the original vectors
        {
            bm::bvector<> bv_c;
            bv_c.bit_or(bv1, bv2, bm::bvector<>::opt_compress);
            int cmp = bv_c.compare(bv3);
            if (cmp != 0)
            {
                std::cerr << "Error: bug in serialization" << std::endl;
            }
            bm::serializer<bm::bvector<> >::buffer sbuf2;
            
            // destructive serialization, "bv_c" should not be used after
            bvs.optimize_serialize_destroy(bv_c, sbuf2);
            
            cout << "BV_C Serialized size:" << sbuf2.size() << endl;
        }
        
    }
    catch(std::exception& ex)
    {
        std::cerr << ex.what() << std::endl;
        delete [] buf1;
        return 1;
    }
    
    delete [] buf1;

    return 0;
}

