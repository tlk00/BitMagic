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

/** \example sample21.cpp

  Example for shifting - erase of bits.
 
  \sa bm::bvector::erase
  \sa bm::bvector::shift_left
*/

/*! \file sample21.cpp
    \brief Example: bvector<> - bit-shifts
*/

#include <iostream>

#include "bm.h"
#include "bmundef.h" /* clear the pre-proc defines from BM */

using namespace std;

// Utility template function used to print container
template<class T> void PrintContainer(T first, T last)
{
    if (first == last)
        cout << "<EMPTY SET>";
    else
        for(;first != last; ++first)
            cout << *first << ", ";
    cout << endl;
}

int main(void)
{
    try
    {
        bm::bvector<>   bv { 1, 10, 100 };
        bv.optimize(); // run memory compression

        cout << "Source set:";
        PrintContainer(bv.first(), bv.end());
        
        // shift the vector right by 1
        //
        bool carry_over = bv.shift_left();
        cout << "CO=" << carry_over << endl;
        PrintContainer(bv.first(), bv.end()); // 0, 9, 99

        // erase 0 bit (equivalent of shift_left)
        //
        bv.erase(0);
        PrintContainer(bv.first(), bv.end()); // 8, 98
        
        // erase/insert manipulations with bit-vector may de-optimize
        // the bit-vector so it needs re-optimization at some point
        //
        bv.optimize();

    }
    catch(std::exception& ex)
    {
        std::cerr << ex.what() << std::endl;
        return 1;
    }
    
    
    return 0;
}

