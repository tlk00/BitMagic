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

/** \example sample1.cpp
  Example how to use bvector<> to set bits and then retrieve indexes of ON bits
 

  \sa bm::bvector<>::get_next() 
  \sa bm::bvector<>::get_first() 
  \sa bm::bvector<>::set()
  \sa bm::bvector<>::count() 
  \sa bm::bvector<>::clear()
 */

/*! \file sample1.cpp
    \brief Example: bvector<> set bits and then retrieve indexes of ON bits
*/
#include <iostream>

#include "bm.h"
#include "bmundef.h" /* clear the pre-proc defines from BM */

using namespace std;

int main(void)
{
    try
    {
        bm::bvector<>   bv { 1, 2, 3 };    // Bitvector variable declaration with init list

        cout << "1. bitcount: " << bv.count() << endl;

        // Set some bits.

        bv.set(10);
        bv.set(100);
        bv.set(1000000);

        // New bitvector's count.

        cout << "2. bitcount: " << bv.count() << endl;


        // Print the bitvector.

        auto value = bv.get_first();
        do
        {
            cout << value;
            value = bv.get_next(value);
            if (value)
            {
                cout << ",";
            }
            else
            {
                break;
            }
        } while(1);

        cout << endl;

        bv.clear();   // Clean up.

        cout << "3. bitcount: " << bv.count() << endl;

        // We also can use operators to set-clear bits;

        bv[10] = true;
        bv[100] = true;
        bv[10000] = true;

        cout << "4. bitcount: " << bv.count() << endl;

        if (bv[10])
        {
            bv[10] = false;
        }

        cout << "5. bitcount: " << bv.count() << endl;
    }
    catch(std::exception& ex)
    {
        std::cerr << ex.what() << std::endl;
        return 1;
    }
    

    return 0;
}

