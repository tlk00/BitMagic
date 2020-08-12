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

/** \example sample3.cpp
 Exmaple demonstrates using bitvectors with different initial
 block allocation strategy. 
 Bitvector 1 (bv1) by default working without RLE compression option
 (best performance, maximum memory consumption). 
 Bitvector 2 (bv2) will be working in compression mode and use less memory.
 
  \sa bm::bvector<>::set_new_blocks_strat() 

  For more information please visit: http://bmagic.sourceforge.net

*/
/*! \file sample3.cpp
    \brief Example: bvector<> with different allocation/compression strategies
*/
#include <stdlib.h>
#include <iostream>

#include "bm.h"
#include "bmundef.h" /* clear the pre-proc defines from BM */

using namespace std;

const unsigned MAX_VALUE = 1000000;

// This procedure creates very dense bitvectors.
// The resulting set will consists mostly from ON (1) bits
// interrupted with small gaps of 0 bits.
static
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


int main(void)
{
    try
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
    }
    catch(std::exception& ex)
    {
        std::cerr << ex.what() << std::endl;
    }

    return 0;
}

