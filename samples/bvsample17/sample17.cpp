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

/** \example sample17.cpp
Example rank and select operations.

*/

/*! \file sample17.cpp
    \brief Example: rank and select operations using rank-select index
*/

#include <stdlib.h>
#include <iostream>
#include <vector>
#include <memory>

#include "bm.h"
#include "bmundef.h" /* clear the pre-proc defines from BM */

int main(void)
{
    try
    {
        bm::bvector<>   bv { 1, 20, 30, 31 }; // init a test bit-vector
        
        // construct rank-select index
        std::unique_ptr<bm::bvector<>::rs_index_type>
                                  rs_idx(new bm::bvector<>::rs_index_type());
        bv.build_rs_index(rs_idx.get());

        // lets find a few ranks here
        //
        auto r1 = bv.rank(20, *rs_idx);
        std::cout << r1 << std::endl;  // 2
        
        r1 = bv.rank(21, *rs_idx);
        std::cout << r1 << std::endl;  // still 2

        r1 = bv.rank(30, *rs_idx);
        std::cout << r1 << std::endl;  // 3

        // position value corrected rank
        // one special case of rank function returns rank-1
        // if position bit is set or just rank, otherwise
        //
        // this is an equivalent of
        // bv.count_range(0, n) - bv.text(n)
        // (just faster, because of the fused rank-test)
        //

        auto r1c = bv.rank_corrected(31, *rs_idx); // 3
        std::cout << r1c << std::endl;
        r1c = bv.rank_corrected(32, *rs_idx); // 4
        std::cout << r1c << std::endl;
        r1c = bv.rank_corrected(33, *rs_idx); // 4
        std::cout << r1c << std::endl;


        // now perform a search for a position for a rank
        //
        bm::bvector<>::size_type pos;
        bool found = bv.select(2, pos, *rs_idx);
        if (found)
            std::cout << pos << std::endl;  // 20
        else
            std::cout << "Rank not found." << std::endl;
        
        found = bv.select(2, pos, *rs_idx);
        if (found)
            std::cout << pos << std::endl;  // 30
        else
            std::cout << "Rank not found." << std::endl;

        found = bv.select(5, pos, *rs_idx);
        if (found)
            std::cout << pos << std::endl;
        else
            std::cout << "Rank not found." << std::endl;  // this!

    }
    catch(std::exception& ex)
    {
        std::cerr << ex.what() << std::endl;
        return 1;
    }
    

    return 0;
}


