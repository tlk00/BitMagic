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

/** \example svsample07a.cpp
  Example of how to search bm::sparse_vector_scanner<> with sparse vector scanner
  and inspect search results using bm::interval_enumerator<>
  to find special values: unique, co-locvated together in the results or dispersed.

  \sa bm::sparse_vector
  \sa bm::sparse_vector_scanner
  \sa bm::interval_enumerator
  \sa bvintervals
  \sa sample23

*/

/*! \file svsample07a.cpp
    \brief Example: sparse_vector<>  search
*/

#include <assert.h>
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
#include "bmintervals.h"
#include "bmundef.h" /* clear the pre-proc defines from BM */

using namespace std;

typedef bm::bvector<> bvector_type;
typedef bm::sparse_vector<unsigned, bvector_type > sparse_vector_u32;
typedef bm::interval_enumerator<bvector_type > interval_enumerator_type;



int main(void)
{
    try
    {
        sparse_vector_u32 sv;

        // use back inserter to add values to the vector
        {
            sparse_vector_u32::back_insert_iterator bit(&sv);
            bit = 0;
            bit = bm::id_max;
            bit = 17;
            bit = 17;
            bit = 5;
            bit = 18;
            bit = 178;
            bit = 178;
            bit = 17;
            bit = 0;
            bit = bm::id_max;
            bit.flush();
        }
        sv.optimize(); // mmemory optimization after active editing

        // just a print-out
        {
            cout << "Sparse vector:" << endl;
            std::for_each(sv.begin(), sv.end(), [](unsigned v) {
                    cout << v << ",";
                });
            cout << endl << endl;
        }

        bm::sparse_vector_scanner<sparse_vector_u32> scanner;

        scanner.bind(sv, false);

        bvector_type bv_res;
        // connect bv_res with allocation pool scanner for better
        // this is optional but known to improve performance on repeated searches
        typename bvector_type::mem_pool_guard mp_guard;
        mp_guard.assign_if_not_set(scanner.get_bvector_alloc_pool(), bv_res);

        sparse_vector_u32::const_iterator it = sv.begin();
        sparse_vector_u32::const_iterator it_end = sv.end();

        // setup two bit-vectors for positive and negative values we seen
        // if values are signed - maintain two separatevectors
        // for positives and negatives
        bvector_type bv_seen_values(bm::BM_GAP);
        bool id_max_seen(false); // bit-vector cannot take bm::id_max value
        for (;it != it_end; ++it)
        {
            auto v = *it;
            if (bvector_type::size_type(v) == bm::id_max)
            {
                if (id_max_seen)
                    continue;
                id_max_seen = true;
            }
            else
            {
                // scanner search is an expensive operation do avoid duplicate lookups
                // maintain bit-vectors of already seen values
                //
                if (bv_seen_values.test(bvector_type::size_type(v)))
                    continue;
                bv_seen_values.set(bvector_type::size_type(v));
            }

            scanner.find_eq(sv, v, bv_res);

            // we have our results, we can now inspect the result bit-vector
            // using interval enumerator to find which values are co-located
            // in one bunch as 000011000...0000
            // or unique as 000010000...00
            // or form non-unique, multi-occurence patterns as 000111000..011..
            //

            interval_enumerator_type ien(bv_res);
            assert(ien.valid()); // as we took value from the search vector
            bvector_type::size_type range_cnt = ien.end() - ien.start() + 1;

            // try to advance
            if (!ien.advance()) // only one interval, cannot advance more
            {
                if (range_cnt == 1)
                {
                    cout << "Value = " << v << " is unique" << endl;
                    continue;
                }
                cout << "Value = " << v << " is colocated" << endl;
                continue;
            }
            cout << "Value = " << v << " is not colocated" << endl;
        } // for it
    }
    catch(std::exception& ex)
    {
        std::cerr << ex.what() << std::endl;
        return 1;
    }

    return 0;
}
