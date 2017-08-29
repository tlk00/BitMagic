/*
Copyright(c) 2002-2009 Anatoliy Kuznetsov(anatoliy_kuznetsov at yahoo.com)

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
*/

/** \example sample4.cpp
 Exmaple demonstrates bitvector serialization/deserialization.
 
For more information please visit:  http://bmagic.sourceforge.net

  \sa bm::serializer
  \sa bm::deserialize 

*/

#include <stdlib.h>
#include <iostream>
#include "bm.h"
#include "bmserial.h"

using namespace std;


// This exmaple demonstrates bitvector serialization/deserialization.



const unsigned MAX_VALUE = 1000000;

// This procedure creates very dense bitvector.
// The resulting set will consists mostly from ON (1) bits
// interrupted with small gaps of 0 bits.

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


unsigned char* serialize_bvector(bm::serializer<bm::bvector<> >& bvs, 
                                 bm::bvector<>& bv)
{
    // It is reccomended to optimize vector before serialization.
    bv.optimize();  

    bm::bvector<>::statistics st;
    bv.calc_stat(&st);

    cout << "Bits count:" << bv.count() << endl;
    cout << "Bit blocks:" << st.bit_blocks << endl;
    cout << "GAP blocks:" << st.gap_blocks << endl;
    cout << "Memory used:"<< st.memory_used << endl;
    cout << "Max.serialize mem.:" << st.max_serialize_mem << endl;

    // Allocate serialization buffer.
    unsigned char*  buf = new unsigned char[st.max_serialize_mem];

    // Serialization to memory.
    unsigned len = bvs.serialize(bv, buf, 0);


    cout << "Serialized size:" << len << endl << endl;

    return buf;
}


int main(void)
{
    bm::bvector<>   bv1;    
    bm::bvector<>   bv2;

    bv2.set_new_blocks_strat(bm::BM_GAP);  //  set DGAP compression mode ON

    fill_bvector(&bv1);
    fill_bvector(&bv2);

    // Prepare a serializer class 
    //  for best performance it is best to create serilizer once and reuse it
    //  (saves a lot of memory allocations)
    //
    bm::serializer<bm::bvector<> > bvs;

    // next settings provide lowest serilized size 
    bvs.byte_order_serialization(false);
    bvs.gap_length_serialization(false);
    bvs.set_compression_level(4);


    unsigned char* buf1 = serialize_bvector(bvs, bv1);
    unsigned char* buf2 = serialize_bvector(bvs, bv2);

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

    delete [] buf1;
    delete [] buf2;

    return 0;
}

