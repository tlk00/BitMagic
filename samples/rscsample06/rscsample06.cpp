/*
Copyright(c) 2002-2020 Anatoliy Kuznetsov(anatoliy_kuznetsov at yahoo.com)

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

/** \example rscsample06.cpp
  Example of how to use bm::rsc_sparse_vector<>::gather to extract values in random order


  \sa bm::rsc_sparse_vector
  \sa bm::rsc_sparse_vector::sync
  \sa bm::rsc_sparse_vector::gather
*/

/*! \file rscsample06.cpp
    \brief Example: bm::rsc_sparse_vector<>::gather

*/

#include <iostream>
#include <vector>
#include <cassert>

#include "bm.h"
#include "bmsparsevec.h"
#include "bmsparsevec_compr.h"
#include "bmundef.h" /* clear the pre-proc defines from BM */

using namespace std;

/// Print sparse vector not NULL elements
template<typename SV> void PrintSV(const SV& sv)
{
    typename SV::const_iterator it = sv.begin();
    typename SV::const_iterator it_end = sv.end();

    for (size_t i = 0; it != it_end; ++it, ++i)
    {
        if (!it.is_null())
            cout << i << ": " << *it << ", ";
    }
    cout << endl;
}

typedef bm::sparse_vector<unsigned, bm::bvector<> >         sparse_vector_u32;
typedef bm::rsc_sparse_vector<unsigned, sparse_vector_u32 > rsc_sparse_vector_u32;


int main(void)
{
    try
    {
        // vector of not NULL elements
        bm::bvector<>         bv1 { 10, 20, 100, 200, 65536 };
        rsc_sparse_vector_u32 csv1(bv1);

        csv1.sync(); // build Rank-Select index for faster access
        bv1.clear(); // we don't need bv1 anymore

        // modify vector elements but only touch the not NULL elements
        // for fast access
        //
        csv1.set(100, 1);
        csv1.set(20, 1);
        csv1.inc(10);
        csv1.inc(65536, 1);

        assert(csv1.in_sync()); // make sure Rank-Select index is still alive

        // gather some of the values in random order
        //

        // idx vector to request RSC vector elements
        bm::bvector<>::size_type idx[] = { 20, 0, 10, 65536 };

        // scratch buffer for Rank-Select coordinates recalculation
        // MUST be the same size as index buffer
        //
        bm::bvector<>::size_type idx_buf[sizeof(idx)/sizeof(idx[0])];


        rsc_sparse_vector_u32::value_type values[sizeof(idx)/sizeof(idx[0])];


        bm::bvector<>::size_type req_size = 4;

        // gather using bm::BM_UNKNOWM sort order (most universal).
        // if sort order is ascending use BM_SORTED for faster operation
        //
        csv1.gather(&values[0], &idx[0], &idx_buf[0], req_size, bm::BM_UNKNOWN);


        // it is possible to use gather() to request NULL values, in this case
        // the secondary array buffer will contain the magic value "bm::id_max"
        //
        for (size_t i = 0; i < req_size; ++i)
        {
            auto ix = idx[i];
            if (idx_buf[i] == bm::id_max)
                cout << "[" << ix << "] = " << "NULL" << endl;
            else
                cout << "[" << ix << "] = " << values[i] << endl;
        } // for

        PrintSV(csv1); // 10: 1, 20: 1, 100: 1, 200: 0, 65536: 1
    }
    catch(std::exception& ex)
    {
        std::cerr << ex.what() << std::endl;
        return 1;
    }
    
    

    return 0;
}

