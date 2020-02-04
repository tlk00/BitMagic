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

/** \example sample22.cpp

  Example on ranges and intervals.

  BitMagic bvector<> implements high performance operation on RLE
  coded bit-vectors, transparently supporting all logical operations
  like intersections. Serialization uses compressive encoding
  (binary interpolative codes) to efficiently store collections 
  of intervals.

  This creates a use case to use compressed bit-vector as an engine
  for ranges and intervals.
 
  \sa bm::bvector::erase
  \sa bm::bvector::shift_left
*/

/*! \file sample22.cpp
    \brief Example: bvector<> - range and interval functions
*/

#include <iostream>

#include "bm.h"

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
        bm::bvector<>   bv(bm::BM_GAP); // use RLE compressed vector from the start
        
        bv.set_range(100, 110); // sets a range of 1s [100, 110] .....11111....

        bv.optimize(); // RLE memory compression

        cout << "Source set:";
        PrintContainer(bv.first(), bv.end());

        cout << "bvector<>::is_all_one_range() demo" << endl;
        // check is range has no 0s in it "....1111...."
        bool all_one = bv.is_all_one_range(100, 110); // true
        cout << all_one << endl;
        all_one = bv.is_all_one_range(100, 111); // false (last bit is 0)
        cout << all_one << endl;

        cout << "bvector<>::is_interval() demo" << endl;
        // verify if the range is all 1s AND flanked by 0s "...011110..."
        bool is_int = bv.is_interval(100, 110); // true
        cout << is_int << endl;
        is_int = bv.is_interval(99, 110); // false (internal range is not all 1s)
        cout << is_int << endl;

        cout << "bvector<>::any_range() demo" << endl;
        // Check is specified inetrval contains at least one 1
        bool any_one = bv.any_range(0, 99); // false
        cout << any_one << endl;
        any_one = bv.any_range(0, 100); // true 
        cout << any_one << endl;

        bv.clear_range(99, 100); // clear bit in [99..100]
        PrintContainer(bv.first(), bv.end());

        bv.keep_range(99, 105); // clear everything outside [99..105]
        PrintContainer(bv.first(), bv.end());
    }
    catch(std::exception& ex)
    {
        std::cerr << ex.what() << std::endl;
        return 1;
    }
    
    
    return 0;
}

