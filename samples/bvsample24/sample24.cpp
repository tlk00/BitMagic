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

/** \example sample24.cpp
@brief Example for finding bit-vector ranges with the specifed number of ON bits

The use case here:
   - search algorithm runs a query, the result set represented as a bit-vector

   - second stage needs to calculate statistics (or run a sub-query) on each
     element of a result set and we want to do it in-parallel so we need to
     create batches (or jobs) with approximately equal complexity

   - we assume that each found element in the bitset has approximately same
     complexity, so to model the situation we need to find range boundaries with
     the specifed population count

    @sa bm::rank_range_split
    @sa bm::bvector::get_enumerator
    @sa bm::bvector::enumerator

*/

/*! \file sample24.cpp
    \brief Example: demo for bm::rank_range_split
*/

#include <stdlib.h>
#include <iostream>
#include <vector>
#include <utility>

#include "bm.h"
#include "bmalgo.h"
#include "bmundef.h" /* clear the pre-proc defines from BM */

using namespace std;

typedef bm::bvector<>::size_type bv_size_type;
typedef std::vector<std::pair<bv_size_type, bv_size_type> > bv_ranges_vector;

int main(void)
{
    try
    {
        const unsigned ranges = 3; // split into 3 ranges
        bm::bvector<>   bv { 10, 20, 100, 200, 300, 655000, bm::id_max-1 };


        bv_size_type cnt = bv.count();
        bv_size_type split_rank = cnt / ranges; // target population count

        bv_ranges_vector pair_vect;

        // run the search for splits, traget bin size is split_rank

        cout << "Target split rank:" << split_rank << endl;

        bm::rank_range_split(bv, split_rank, pair_vect);

        // go through each range and print all bit values
        for (size_t k = 0; k < pair_vect.size(); ++k)
        {
            const auto& p = pair_vect[k];
            cout << k << ": [" << p.first << ".." << p.second << "] ";

            // find and print bits in the target range with bvector<>::enumerator
            bm::bvector<>::enumerator en = bv.get_enumerator(p.first);
            for (; en.valid(); ++en)
            {
                auto v = *en;
                if (v > p.second)
                    break;
                cout << v << ", ";
            } // for en
            cout << endl;
        } // for k

    }
    catch(std::exception& ex)
    {
        std::cerr << ex.what() << std::endl;
        return 1;
    }
    

    return 0;
}

