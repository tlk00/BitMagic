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

/** \example sample18a.cpp
  Example of bit-vector import from an external bit-stream
  \sa bm::bvector::bulk_insert_iterator
*/

/*! \file sample18a.cpp
    \brief Example: import from an external bit-stream
*/

#include <stdlib.h>
#include <assert.h> 
#include <iostream>

#include "bm.h"
#include "bmbvimport.h"
#include "bmundef.h" /* clear the pre-proc defines from BM */

using namespace std;

int main(void)
{
    try
    {
        bm::bvector<>   bv;

        // specify bits in an array of unsigned ints (32-bit)
        // we expect that each set bit can be efficiently imported into a
        // bvector<>. Import function breaks the bit stream into blocks,
        // each block can be optimized/compressed on the fly, while data is
        // still in the CPU cache

        unsigned int arr[2058] = {0, };
        arr[0] = 1 << 16;
        arr[2047] = 1u << 31;
        arr[2048] = 1u << 7;

        bm::bit_import_u32(bv, arr, sizeof(arr)/sizeof(arr[0]), true/* optimize on the fly*/);
        auto cnt = bv.count();
        cout << "Imported " << cnt << " bits." << endl;

        assert(cnt == 3);
        assert(bv.test(16));
        assert(bv.test(65535));
        assert(bv.test(65536 + 7));

    }
    catch(std::exception& ex)
    {
        std::cerr << ex.what() << std::endl;
        return 1;
    }

    return 0;
}


