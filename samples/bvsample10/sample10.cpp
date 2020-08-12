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

/** \example sample10.cpp
  Example of how to get random subset of a bit-vector
 
  \sa bm::random_subset  
 */

/*! \file sample10.cpp
    \brief Example: bvector<> generation of random sub-set
*/


#include <iostream>

#include "bm.h"
#include "bmrandom.h"
#include "bmundef.h" /* clear the pre-proc defines from BM */

using namespace std;

template<class T> void PrintContainer(T first, T last)
{
    if (first == last)
        cout << "<EMPTY SET>";
    else
        for(;first != last; ++first)
            cout << *first << ";";
    cout << endl;
}


int main(void)
{
    try
    {
        bm::bvector<>   bv;
        // -----------------------------------------------
        // set some bits
        //
        bm::bvector<>::size_type i;
        for (i = 0; i < 30000; i+=10)
        {
            bv.set(i);
        }

        for (i = 300000; i < 400000; i+=100)
        {
            bv.set(i);
        }

        bm::bvector<>   bvsubset1;
        bm::bvector<>   bvsubset2;

        // random sampler instance can be shared between calls
        //
        bm::random_subset<bm::bvector<> > rand_sampler;
        rand_sampler.sample(bvsubset1, bv, 20);
        rand_sampler.sample(bvsubset2, bv, 20);
     
        PrintContainer(bvsubset1.first(), bvsubset1.end());
        cout << endl;
        PrintContainer(bvsubset2.first(), bvsubset2.end());
    }
    catch(std::exception& ex)
    {
        std::cerr << ex.what() << std::endl;
        return 1;
    }

    return 0;
}

