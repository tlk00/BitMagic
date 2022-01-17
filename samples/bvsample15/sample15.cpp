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

/** \example sample15.cpp
Example for finding first and last bits in bit-vector (dynamic range).

Ranges of bit-vectors can be used to find probability of intersection.
For instance, in some corner cases AND product can be predicted empty if vectors
belong to different ranges.

    @sa bm::bvector::find()
    @sa bm::bvector::find_reverse()
    @sa bm::bvector::find_range()
*/

/*! \file sample15.cpp
    \brief Example: bvector<> methods to find last bit and bit-vectors effective range
*/

#include <stdlib.h>
#include <iostream>
#include <vector>
#include <cassert>

#include "bm.h"
#include "bmundef.h" /* clear the pre-proc defines from BM */

using namespace std;


const unsigned MAX_VALUE = 10000000;

static
void fill_bvector(bm::bvector<>* bv)
{
    unsigned start = MAX_VALUE / (unsigned)(rand()%10);
    for (unsigned i = start; i < MAX_VALUE; ++i)
    {
        if ((rand() % 10))
        {
            bv->set(i);
        }
    }
}


int main(void)
{
    try
    {
        bm::bvector<>   bv1;
        bm::bvector<>   bv2;

        fill_bvector(&bv1);
        fill_bvector(&bv2);
        
        cout << "bv1 count = " << bv1.count() << endl;
        cout << "bv2 count = " << bv2.count() << endl;

        bool found;
        bm::bvector<>::size_type first, last, pos, second;
        
        found = bv1.find(first);
        if (found)
            cout << "bv1 first = " << first << endl;

        // make a repeat find starting on a discovered position
        //
        // find will return the same position (if it is set)
        found = bv1.find(first, pos);
        assert (found);
        {
            cout << "bv1 pos = " << pos << endl; // will be same as first
            assert(pos == first);
        }

        // Q: what if we need a second set position?
        // A: increament previously found index and find() again
        //
        found = bv1.find(first+1, second); // use first + 1 
        assert (found);
        {
            cout << "bv1 second = " << second << endl;
            assert(second > first);
            assert(bv1.test(second));
        }

        found = bv1.find_reverse(last);
        if (found)
            cout << "bv1 last = " << last << endl;

        found = bv2.find(first);
        if (found)
            cout << "bv2 first = " << first << endl;
        
        found = bv2.find_reverse(last);
        if (found)
            cout << "bv2 last = " << last << endl;
        
        found = bv1.find_range(first, last);
        if (found)
            cout << "bv1 range = [" << first << ", " << last << "]" << endl;
        
        found = bv2.find_range(first, last);
        if (found)
            cout << "bv2 range = [" << first << ", " << last << "]" << endl;

    }
    catch(std::exception& ex)
    {
        std::cerr << ex.what() << std::endl;
        return 1;
    }
    

    return 0;
}

