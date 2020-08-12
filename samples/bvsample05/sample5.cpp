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

/** \example sample5.cpp
 Example demonstrates using enumerators - the fastest way to retrieve 
 indexes of 1 bits from the bitvector. This approach works faster than
 get_first()/get_next() functions.
 
  \sa bm::bvector<>::enumerator 
  \sa bm::bvector<>::first()
  \sa bm::bvector<>::end()
  \sa bm::bvector<>::get_enumerator()
*/

/*! \file sample5.cpp
    \brief Example: bvector<>::enumerator use
*/

#include <iostream>
#include <algorithm>

#include "bm.h"
#include "bmundef.h" /* clear the pre-proc defines from BM */

using namespace std;

inline
void Print(bm::bvector<>::size_type n)
{
    cout << n << endl;;
}

int main(void)
{
    try
    {
        bm::bvector<>   bv;

        bv[10] = true;
        bv[100] = true;
        bv[10000] = true;
        bv[65536] = true;
        bv[65537] = true;
        bv[65538] = true;
        bv[65540] = true;

        bm::bvector<>::enumerator en = bv.first();
        bm::bvector<>::enumerator en_end = bv.end();

        while (en < en_end)
        {
            cout << *en << ", ";
            ++en;  // Fastest way to increment enumerator
        }
        cout << endl;

        en = bv.first();

        // This is not the fastest way to do the job, because for_each
        // often will try to calculate difference between iterators,
        // which is expensive for enumerators.
        // But it can be useful for some STL loyal applications.

        std::for_each(en, en_end, Print);
        cout << endl;

        // example to illustrate random positioning of enumerator 
        // go to a random bit number, enumerator automatically finds the available bit
        //
        en.go_to(65537);
        for (; en.valid(); ++en)
        {
            cout << *en << ", ";
        }
        cout << endl;
    }
    catch(std::exception& ex)
    {
        std::cerr << ex.what() << std::endl;
        return 1;
    }
    return 0;
}
