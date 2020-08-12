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
 
  \sa bm::bvector::set_range
  \sa bm::bvector::clear_range
  \sa bm::bvector::keep_range
  \sa bm::bvector::is_all_one_range
  \sa bm::bvector::any_range
  \sa bm::bvector::find_reverse()
  \sa bm::is_interval
  \sa bm::find_interval_end
  \sa bm::find_interval_start
  \sa bm::deserialize_range

  \sa sample23.cpp
  \sa bvintervals
*/

/*! \file sample22.cpp
    \brief Example: bvector<> - ranges and intervals functions
*/

#include <iostream>
#include <assert.h>

#include "bm.h"
#include "bmserial.h"
#include "bmintervals.h"
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
        bm::bvector<> bv(bm::BM_GAP); // use RLE compressed vector from the start
        
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
        bool is_int = bm::is_interval(bv, 100, 110); // true
        cout << is_int << endl;
        is_int = bm::is_interval(bv, 99, 110); // false (internal range is not all 1s)
        cout << is_int << endl;

        cout << "bvector<>::any_range() demo" << endl;
        // Check is specified interval contains at least one 1
        bool any_one = bv.any_range(0, 99); // false
        cout << any_one << endl;
        any_one = bv.any_range(0, 100); // true 
        cout << any_one << endl;

        cout << "bvector<>::find_reverse() demo" << endl;
        bm::bvector<>::size_type pos;

        bool found = bv.find_reverse(256, pos);
        assert(found);
        cout << pos << endl; // 110


        // interval boundaries detection
        //
        //
        cout << "bvector<>::find_interval demo" << endl;

        // interval end search from interval start
        found = bm::find_interval_end(bv, 100, pos);
        if (found)
            cout << pos << endl; // 110
        else
            cout << "Not found." << endl;

        // interval end start from a non-interval location
        // - it will not find anything
        found = bm::find_interval_end(bv, 99, pos);
        if (found)
            cout << pos << endl; 
        else
            cout << "Not found." << endl; // This ! 

        // start search from a position within the interval
        found = bm::find_interval_start(bv, 105, pos);
        if (found)
            cout << pos << endl; // 100
        else
            cout << "Not found." << endl;

        bm::bvector<> bv3(bv);
        

        // range editing 
        //
        cout << endl;

        bm::bvector<> bv2;
        bv2.copy_range(bv, 101, 105); // make a copy of [101..105]

        bv.clear_range(99, 100); // clear bits in [99..100]
        PrintContainer(bv.first(), bv.end());

        bv.keep_range(99, 105); // clear everything outside [99..105]
        PrintContainer(bv.first(), bv.end()); // print edited vector

        PrintContainer(bv2.first(), bv2.end()); // print range copy vector

        bool eq = bv.equal(bv2); // make sure both vectors are the same
        cout << eq << endl; // true

        // demo (de)serialization range
        // BM can serialize a bit-vector and selectively restore only 
        // part of it (range). Adding block bookmarks makes target range
        // extraction faster. 
        //
        {
            // configure serializer to use bookmarks every 64 blocks
            // for faster extraction. Higher number gives you smaller BLOB
            // but slower speed. Tweak it as needed.
            // 
            // (range deserialization would work even without bookmarks)
            //
            bm::serializer<bm::bvector<> > ser;
            ser.set_bookmarks(true, 64); 
            
            bm::serializer<bm::bvector<> >::buffer buf;
            ser.serialize(bv3, buf);
            cout << "BLOB size=" << buf.size() << endl;

            // equivalent of a copy_range right from a compressed BLOB
            bm::bvector<> bv4;
            bm::deserialize_range(bv4, buf.data(), 101, 105);

            PrintContainer(bv4.first(), bv4.end()); // print deserialized range vector

            eq = bv4.equal(bv2); // make sure both vectors are the same
            cout << eq << endl; // true
        }

    }
    catch(std::exception& ex)
    {
        std::cerr << ex.what() << std::endl;
        return 1;
    }
        
    return 0;
}

