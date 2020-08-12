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

/** \example sample18.cpp
  Example of bulk insert iterator
  \sa bm::bvector::insert_iterator
  \sa bm::bvector::bulk_insert_iterator 

*/

/*! \file sample18.cpp
    \brief Example: bulk insert iterator

    Bulk insert iterator uses internal buffer to accumulate a batch 
    of bit indexes and flush it faster. 
    This delay makes its behavior different from classic 
    deterministic STL iterators.
*/

#include <stdlib.h>
#include <iostream>

#include "bm.h"
#include "bmundef.h" /* clear the pre-proc defines from BM */

using namespace std;

int main(void)
{
    try
    {
        bm::bvector<>   bv;

        // use regular insert iterator to add bits to vector
        bm::bvector<>::insert_iterator iit = bv.inserter();
        for (unsigned i = 5; i != 0; --i)
        {
            iit = i;
            cout << bv.count() << ", ";  // note that bits are added immediately
        }
        cout << endl;

        bv.clear();  // clear the vector

        // bulk insert adds bits faster (on sorted data), but provides no guarantee 
        // that it happens immediately. Instead it uses an internal buffer to accumulate it first
        // and periodically flushes it into the bit-vector.
        // 
        {
            bm::bvector<>::bulk_insert_iterator bulk_iit(bv);
            for (bm::bvector<>::size_type i = 5; i != 0; --i)
            {
                bulk_iit = i;
                cout << bv.count() << ", ";  // note that bits are NOT added immediately
            }
            cout << endl;
        } // <= scope is a point when data flush happens (on destruction)
        cout << bv.count() << endl; // now all bits are added

        bv.clear();

        // bulk content flash can be forced, before going out-of-scope
        // which may be problematic due to possible exceptions from the destructor.
        // 
        {
            bm::bvector<>::bulk_insert_iterator bulk_iit(bv);
            for (bm::bvector<>::size_type i = 0; i < 5; ++i)
                bulk_iit = i;
            cout << bv.count() << endl; // not added yet
            bulk_iit.flush();           // force bulk iterator to submit the data
            cout << bv.count() << endl; // now added!
        }

    }
    catch(std::exception& ex)
    {
        std::cerr << ex.what() << std::endl;
        return 1;
    }

    return 0;
}


