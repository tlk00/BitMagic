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

/** \example svsample10.cpp
    Scanner searches: GT, GE, LT, LE, RANGE[from..to]
 
  \sa bm::sparse_vector
  \sa bm::sparse_vector_scanner<>::find_gt
  \sa bm::sparse_vector_scanner<>::find_ge
  \sa bm::sparse_vector_scanner<>::find_lt
  \sa bm::sparse_vector_scanner<>::find_le
  \sa bm::sparse_vector_scanner<>::find_range

*/

/*! \file svsample10.cpp
    \brief Example: Succinct vector searches: GT, GE, LT, LE, RANGE[from..to]
*/

#include <assert.h>
#include <stdlib.h>

#include <iostream>
#include <vector>
#include <utility>

#include "bm.h"
#include "bmsparsevec.h"
#include "bmsparsevec_algo.h"
#include "bmundef.h" /* clear the pre-proc defines from BM */

using namespace std;


typedef bm::sparse_vector<unsigned, bm::bvector<> > svector_u32;
typedef bm::sparse_vector<int, bm::bvector<> >      svector_i32;

template<typename SV, typename BV>
void PrintResults(const SV& sv, const BV& bv)
{
    auto en = bv.get_enumerator(0);
    if (!en.valid())
    {
        cout << "<EMPTY>" << endl;
        return;
    }
    auto cnt = bv.count();
    cout << "size=" << cnt << "  ";
    for (; en.valid(); ++en)
    {
        auto idx = *en;
        auto v = sv.get(idx);
        cout << idx << ":" << v << ", ";
    } // for
    cout << endl;
}


int main(void)
{
    try
    {
        svector_u32 sv1;
        svector_i32 sv2(bm::use_null);

        // fill in vectors with some test data
        //
        {
            auto bit = sv1.get_back_inserter();
            bit = 1;
            bit = 2;
            bit = 20;
            bit = 0;

            bit.flush();
        }
        {
            auto bit = sv2.get_back_inserter();
            bit = 1;
            bit = -2;
            bit = 0;
            bit.add_null();
            bit = 30;

            bit.flush();
        }

        bm::sparse_vector_scanner<svector_u32> scanner_u32;
        bm::sparse_vector_scanner<svector_i32> scanner_i32;

        // Here we perform scanner searches of various types
        // Please note that NULL values are never included into
        // the search results.
        // Any comparison with NULL results as false.

        // bm::sparse_vector_scanner performs searches without reverse
        // transformation, using logical ops, returning a result-set
        // bvector<>
        //

        // GT search
        {
            bm::bvector<> bv_res;
            scanner_u32.find_gt(sv1, 1, bv_res);
            PrintResults(sv1, bv_res); // size=2  1:2, 2:20,
        }

        // LT search
        {
            bm::bvector<> bv_res;
            scanner_u32.find_lt(sv1, 0, bv_res);
            PrintResults(sv1, bv_res); // <EMPTY>
        }

        // GE search
        {
            bm::bvector<> bv_res;
            scanner_i32.find_ge(sv2, -1, bv_res);
            PrintResults(sv2, bv_res); // size=3  0:1, 2:0, 4:30,
        }

        // LE search
        {
            bm::bvector<> bv_res;
            scanner_i32.find_le(sv2, 30, bv_res);
            PrintResults(sv2, bv_res); // size=4  0:1, 1:-2, 2:0, 4:30,
        }

        // range search [from..to] closed range
        {
            bm::bvector<> bv_res;
            scanner_i32.find_range(sv2, -1, 30, bv_res);
            PrintResults(sv2, bv_res); // size=3  0:1, 2:0, 4:30,
        }

        // a bit more complex search expression
        // (> 10) OR (<= -1)
        //
        // (> 10) AND (<= -1)
        //
        {
            bm::bvector<> bv_res_gt;
            bm::bvector<> bv_res_le;

            scanner_i32.find_gt(sv2, 10, bv_res_gt); // > 10
            scanner_i32.find_le(sv2, -1, bv_res_le); // <= -1

            bm::bvector<> bv_res;
            bv_res.bit_or(bv_res_gt, bv_res_le); // OR two together

            PrintResults(sv2, bv_res); // size=2  1:-2, 4:30,

            bv_res.bit_and(bv_res_gt, bv_res_le); // AND two together - empty set
            PrintResults(sv2, bv_res); // <EMPTY>

        }

    }
    catch(std::exception& ex)
    {
        std::cerr << ex.what() << std::endl;
        return 1;
    }
    
    

    return 0;
}

