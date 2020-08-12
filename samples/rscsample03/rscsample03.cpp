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

/** \example rscsample03.cpp
  Example of how to use bm::rsc_sparse_vector<>::const_iterator
 
  \sa bm::rsc_sparse_vector
  \sa bm::rsc_sparse_vector::const_iterator
  \sa bm::rsc_sparse_vector::const_iterator::go_to

*/

/*! \file rscsample03.cpp
    \brief Example: bm::rsc_sparse_vector<>::const_iterator

    @sa rscsample02.cpp
*/

#include <iostream>
#include <vector>

#include "bm.h"
#include "bmsparsevec.h"
#include "bmsparsevec_compr.h"
#include "bmundef.h" /* clear the pre-proc defines from BM */

using namespace std;

/// Print sparse vector content
template<typename SV> void PrintSV(const SV& sv)
{
    typename SV::const_iterator it = sv.begin();
    typename SV::const_iterator it_end = sv.end();

    for (; it != it_end; ++it)
    {
        if (it.is_null())
            cout << "NULL";
        else
            cout << *it;
        cout << ", ";
    }
    cout << endl;
}

typedef bm::sparse_vector<unsigned, bm::bvector<> >         sparse_vector_u32;
typedef bm::rsc_sparse_vector<unsigned, sparse_vector_u32 > rsc_sparse_vector_u32;


int main(void)
{
    try
    {
        rsc_sparse_vector_u32 csv1;

        // add sample data using back insert iterator
        //
        {
            auto bit = csv1.get_back_inserter();
            bit = 10;
            bit = 11;
            bit.add_null();
            bit = 13;
            bit = 14;
            bit.add_null(2);
            bit = 256;
            bit.flush();
        }

        // sync call updates rank-select index of the succinct container and
        // it is necessary for the correct work of the const_iterator
        // and bm::rsc_sparse_vector<>::decode()
        // and bm::rsc_sparse_vector<>::decode_buf() methods
        //
        csv1.sync();

        PrintSV(csv1); // 10, 11, NULL, 13, 14, NULL, NULL, 256,


        // another way to use iterator at a random position
        // is to use rsc_sparse_vector<>::get_const_iterator()
        {
            rsc_sparse_vector_u32::const_iterator it =
                                                csv1.get_const_iterator(2);
            if (it.valid())  // check if iterator is avlid and not at the end
            {
                do
                {
                    auto v = it.value();
                    if (it.is_null())
                        cout << "NULL";
                    else
                        cout << v;
                    cout << ", ";
                } while (it.advance()); // while it is able to advance
                cout << endl;

                // iterator can be re-positioned
                //
                it.go_to(3); // position iterator at 3
                cout << it.value() << endl; // 13
            }
        }

    }
    catch(std::exception& ex)
    {
        std::cerr << ex.what() << std::endl;
        return 1;
    }
    
    

    return 0;
}

