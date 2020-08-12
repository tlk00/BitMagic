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

/** \example sample20.cpp

  Example for shifting - insertion of bits.
 
  \sa bm::bvector::insert
  \sa bm::bvector::shift_right
*/

/*! \file sample20.cpp
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

        cout << "Source set:";
        PrintContainer(bv.first(), bv.end());
        
        // shift the vector right by 1
        //
        bv.shift_right();
        PrintContainer(bv.first(), bv.end()); // 2, 11, 101

        // insert 0 bit into position 0 (equivalent of shift_right)
        //
        bv.insert(0, false);
        PrintContainer(bv.first(), bv.end()); // 3, 12, 102
        
        // insert at a random position
        bv.insert(4, true);
        PrintContainer(bv.first(), bv.end()); // 3, 4, 13, 103

    }
    catch(std::exception& ex)
    {
        std::cerr << ex.what() << std::endl;
        return 1;
    }
    
    
    return 0;
}

