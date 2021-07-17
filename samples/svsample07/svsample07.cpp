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

/** \example svsample07.cpp
  Example of how to use bm::str_sparse_vector<> - succinct container for
  sorted bit-transposed collection of integers

  \sa bm::sparse_vector
  \sa bm::sparse_vector_scanner<>::lower_bound
*/

/*! \file svsample07.cpp
    \brief Example: sparse_vector<> lower bound search
*/

#include <iostream>
#include <vector>
#include <chrono>
#include <algorithm>
#include <random>
#include <algorithm>
#include <stdexcept>

#include "bm.h"
#include "bmsparsevec.h"
#include "bmsparsevec_algo.h"
#include "bmundef.h" /* clear the pre-proc defines from BM */

using namespace std;

typedef bm::sparse_vector<bm::id_t, bm::bvector<> > sparse_vector_u32;

// generate collection of strings from integers and shuffle it
//
static
void generate_set(vector<unsigned>& vec)
{
    const unsigned max_coll = 50000;
   
    vec.resize(0);
    for (unsigned i = 10; i < max_coll; i += rand() % 3)
    {
        vec.emplace_back(i);
    } // for i
    
    // shuffle the data set
    //
    std::random_device rd;
    std::mt19937       g(rd());
    std::shuffle(vec.begin(), vec.end(), g);
}

// insertion sort takes data from unsorted vector places it into sparse vector
// maintaining correct sorted order (for fast search)
//
static
void insertion_sort(sparse_vector_u32& sv, const vector<unsigned>& vec)
{
    // scanner object is re-used throught the processing
    //
    bm::sparse_vector_scanner<sparse_vector_u32> scanner;
    sparse_vector_u32::size_type pos;

    for (const unsigned u : vec)
    {
        bool found = scanner.bfind(sv, u, pos);
        (void)found; // just to silence the unused variable warning
        
        sv.insert(pos, u);
        
    } // for u
}


int main(void)
{
    try
    {
        sparse_vector_u32 sv;
        vector<unsigned> vec;
        generate_set(vec);
        
        insertion_sort(sv, vec);
        
        sv.optimize(); // mmemory optimization after active editing
        
        // validate the results to match STL sort
        std::sort(vec.begin(), vec.end());
        {
            vector<unsigned>::const_iterator vit = vec.begin();
            sparse_vector_u32::const_iterator it = sv.begin();
            sparse_vector_u32::const_iterator it_end = sv.end();
            for (; it != it_end; ++it, ++vit)
            {
                unsigned u = *it;
                if (*vit != u)
                {
                    cerr << "Mismatch at:" << u << "!=" << *vit << endl;
                    return 1;
                }
            } // for
        }
        cout << "Sort validation Ok." << endl;
    }
    catch(std::exception& ex)
    {
        std::cerr << ex.what() << std::endl;
        return 1;
    }

    return 0;
}
