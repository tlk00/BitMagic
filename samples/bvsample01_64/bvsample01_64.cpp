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

/** \example bvsample01_64.cpp
  Example how to use bvector<> in 64-bit mode
 */

/*! \file bvsample01_64.cpp
    \brief Example: how to use 64-bit mode

    By default BitMagic uses 32-bit address mode even when if it is
    compiled for 64-bit version. This way it supports 2^32-1 address space.


    There are two ways to run BitMagic in 64-bit address mode:

    1. define BM64ADDR in your project
    Or
    2. include "bm64.h"

    Limitations: 
    - you CANNOT use 32-bit and 64-bit in the same compilation unit.
    - Current implementation internally is 48-bit (which is a lot),
      so your range will be [0..2^48-1]
    - 32-bit vectors can be serialized and read as 64-bit, but 
      not vice versa.
*/
#include <iostream>
#include <assert.h>

#include "bm64.h"
#include "bmundef.h" /* clear the pre-proc defines from BM */

using namespace std;

int main(void)
{
    try
    {
        bm::bvector<>   bv { 1, 2, 3, bm::id_max-1 };    // Bitvector variable declaration with init list
        bm::bvector<>::size_type count = bv.count();

        cout << "BitCount = " << count << endl;
        cout << "Max possible ID = " << bm::id_max-1 << endl;

        bm::bvector<>::size_type first, last;
        bool range_found = bv.find_range(first, last);
        assert(range_found);
        (void)range_found; // to silence the warning on unused var

        cout << "[" << first << ", " << last << "]" << endl;

        bm::bvector<>   bv_full;
        bv_full.set(); // set the whole vector to 1111111....11111

        auto full_count = bv_full.count();
        cout << "Full vector bitcount = " << full_count << endl;

        bv_full.set(bm::id_max - 1, false);

        bv &= bv_full;
        bm::bvector<>::enumerator en = bv.first();
        for (; en.valid(); ++en)
        {
            bm::id64_t idx = *en;
            cout << idx << ", ";
        }
        cout << endl;
    }
    catch(std::exception& ex)
    {
        std::cerr << ex.what() << std::endl;
        return 1;
    }
    

    return 0;
}

