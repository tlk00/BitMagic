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

/** \example sample2.cpp
  Example demonstrates using set operations AND, OR, XOR, etc.
  \sa bvsetalgebra.cpp
*/

/*! \file sample2.cpp
    \brief Example: bvector<> set algebra operations AND, OR, XOR, etc.
*/

#include <iostream>

#include "bm.h"
#include "bmundef.h" /* clear the pre-proc defines from BM */

using namespace std;

static
void print_bvector(const bm::bvector<>& bv)
{
    bm::bvector<>::size_type value = bv.get_first();
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
}

int main(void)
{
    try
    {
        bm::bvector<>   bv1;
        bm::bvector<>   bv2;
        bm::bvector<>   bv3;

        bv1.set(10);
        bv1.set(100);
        bv1.set(1000000);


        bv2.set(10);
        bv2.set(100);

        // Logical AND operation on bv2 (bv1 is the argument)
        // bv2 = bv2 AND bv1

        bv3 = bv1 & bv2;
        print_bvector(bv3);

        bv2 &= bv1;  // You also can use: bv2.bit_and(bv1);
        print_bvector(bv2);
        
        // bv2 = bv2 OR bv1

        bv3 = bv1 | bv2;
        print_bvector(bv3);

        bv2 |= bv1;  //  You can also use: bv2.bit_or(bv1);
        print_bvector(bv2);

        
        bv1.set(1000000, false);
        
        // bv2 = bv2 SUB bv1

        bv3 = bv2 - bv1;
        print_bvector(bv3);

        bv2 -= bv1;   // You can also use: bv2.bit_sub(bv1);
        print_bvector(bv2);

        // bv2 XOR bv1

        bv3 = bv2 ^ bv1;
        print_bvector(bv3);

        // product of XOR is a mismatch vector (definition of XOR)
        //
        {
            bm::bvector<>::size_type pos;
            bool f = bv3.find(pos);
            if (f)
            {
                cout << "XOR mismatch position = " << pos << endl;
            }

            // if we need to find only first mismatch we don't need a full
            // XOR product, we can use bvector<>::find_first_mismatch()
            //
            f = bv2.find_first_mismatch(bv1, pos);
            if (f)
            {
                cout << "search mismatch position = " << pos << endl;
            }
        }

        bv2 ^= bv1;  // You can also use: bv2.bit_xor(bv1);
        print_bvector(bv2);

        // For lexicographical comparison there is set of overloaded
        // operators and function compare (see also bvector<>::equal() )

        if (bv2 == bv3)
        {
            cerr << "Equivalent. Comparison result = "
                 << bv2.compare(bv3) << endl;
        }
        else
        {
            cout << "Error." << endl;
            return 1;
        }
    }
    catch(std::exception& ex)
    {
        std::cerr << ex.what() << std::endl;
    }


    return 0;
}
