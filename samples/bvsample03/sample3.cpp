/*
Copyright(c) 2002-2005 Anatoliy Kuznetsov(anatoliy_kuznetsov at yahoo.com)

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

/** \example sample3.cpp
 Exmaple demonstrates using bitvectors with different initial
 block allocation strategy. 
 Bitvector 1 (bv1) by default working without RLE compression option
 (best performance, maximum memory consumption). 
 Bitvector 2 (bv2) will be working in compression mode and use less memory.
 
  \sa bm::bvector<>::set_new_blocks_strat() 

  For more information please visit: http://bmagic.sourceforge.net

*/

#include <stdlib.h>
#include <iostream>
#include "bm.h"

using namespace std;

const unsigned MAX_VALUE = 1000000;

// This procedure creates very dense bitvectors.
// The resulting set will consists mostly from ON (1) bits
// interrupted with small gaps of 0 bits.

void fill_bvector(bm::bvector<>* bv1, bm::bvector<>* bv2)
{
    for (unsigned i = 0; i < MAX_VALUE; ++i)
    {
        if (rand() % 2500)
        {
            bv1->set_bit(i);
            bv2->set_bit(i);
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


int main(void)
{
    bm::bvector<>   bv1;    
    bm::bvector<>   bv2;

    bv2.set_new_blocks_strat(bm::BM_GAP);  //  set DGAP compression mode ON

    fill_bvector(&bv1, &bv2);  // Fill both bvectors with the same values

    // For a given distrubution statistics should demonstrate
    // lower memory consumption for the vector with compression

    print_statistics(bv1);    
    print_statistics(bv2);

    // Now run optimization procedure for bv1 and see statistics.

    bv1.optimize();

    print_statistics(bv1);    

    return 0;
}

