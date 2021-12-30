/*
Copyright(c) 2002-2021 Anatoliy Kuznetsov(anatoliy_kuznetsov at yahoo.com)

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

/*! \example sample25.cpp

This example illustrates verious traversal methods based on iterator/enumerator or algorithms.

    \sa bm::bvector
    \sa bm::bvector::get_enumerator
    \sa bm::bvector::enumerator
    \sa bm::visit_each_bit
    \sa bm::visit_each_bit_range
    \sa bm::for_each_bit
    \sa bm::for_each_bit_range
*/

/*! \file sample25.cpp
    \brief Example: demo for bit-vector traversal techniques
*/

#include <stdlib.h>
#include <iostream>
#include <vector>
#include <utility>
#include <cassert>

#include "bm.h"
#include "bmalgo.h"
#include "bmundef.h" /* clear the pre-proc defines from BM */

using namespace std;

extern "C" {

static
int bit_visitor_callback(void* handle_ptr, bm::id_t bit_idx) noexcept
{
    assert(handle_ptr);
    unsigned* cnt = (unsigned*)handle_ptr;
    if (*cnt >= 5)
        return -1; // negative code indicates we requested to stop iteration
    *cnt += 1;
    cout << bit_idx << ", ";
    return 0;
}

} // extern C

/**
* Visitor calss for bm::for_each_bit() algorithm
* It is NOT a classic functor (with function call operator() overload)
* For efficiency it needs to support two methods: add_bits and add_range
* which corresponds to different internal methods of representation of 
* sets in the bm::bvector<> (bit-stream and D-GAP/RLE representation) 
* 
* This method somewhat exposes internals and requires writing a custom
* object, but it is also the fastest method, used inside BitMagic 
* library a lot.
*/
struct bit_visitor_functor
{
    using size_type = bm::bvector<>::size_type;

    /**
    * Decoded bits callback
    * @param offset - index offset in the bit-vector (decoding address context)
    * @param bits - array of set bit indexes within the address context
                    offset + bits[i] gives us global index in the bit_vector
      @param size - number of decoded bits (size of bits array)
    */
    int add_bits(size_type offset, 
                 const unsigned char* bits, 
                 unsigned size)
    {
        for (unsigned i = 0; i < size; ++i)
        {
            if (cnt >= 5)
                return -1; // negative code indicates we requested to stop iteration
            cout << offset + bits[i] << ", ";
            ++cnt;
        }
        return 0;
    }


    /**
    * Decoded range callback
    * @param offset - index offset in the bit-vector (decoding address context)
      @param size - number of consequitive ON bits (size of bits array)
    */
    int add_range(size_type offset, size_type size)
    {
        for (size_type i = 0; i < size; ++i)
        {
            if (cnt >= 5)
                return -1; // negative code indicates we requested to stop iteration
            cout << offset + i << ", ";
            ++cnt;
        }
        return 0;
    }
    unsigned cnt = 0; // counter variable to limit the traversal
};



int main(void)
{
    try
    {
        bm::bvector<>   bv { 0, 1, 10, 20, 100, 200, 300, 655000, bm::id_max-1 };
        bv.optimize(); // compress the bit-vector

        // CASE 1:
        // print first 5 set elements using bm::bvector::enumerator
        {
            auto en = bv.get_enumerator(0); // 0 is a starting position in bit-vector
            for (unsigned cnt = 0; en.valid() && cnt < 5; ++en, ++cnt)
                cout << *en << ", ";
        }
        cout << endl;

        // CASE 2:
        // Use C-style callback with bm::visit_each_bit() algorithm
        // callback function receives the execution context as a void ptr, free for
        // reinterpretation, and can signal to interrupt the cycle by negative return code
        // return code is passed back from the visit algorithm
        {
            unsigned cnt = 0;
            void* ctx_ptr = &cnt;
            int res = bm::visit_each_bit(bv, ctx_ptr, bit_visitor_callback);
            cout << "\nvisitor return code = " << res << endl; // -1
        }

        // CASE 3:
        // use bm::for_each_bit() and a visitor class to traverse the bit-vector
        // This method is the fastest, but it exposes some details on internals 
        // of the bit-vector
        {
            bit_visitor_functor func; // declare the visitor object
            int res = bm::for_each_bit(bv, func);
            cout << "\nvisitor return code = " << res << endl; // -1
        }

        // If we want to traverse a range in bit-vector
        // BitMagic offers range traversal algorithms with [from..to] 
        // closed range semantics
        //

        bm::bvector<>::size_type from(15), to(550);
        cout << "\nRange traversal: [" << from << ".." << to << "]" << endl;

        // CASE 2a: bm::visit_each_bit_range
        {
            unsigned cnt = 0;
            int res = bm::visit_each_bit_range(bv, from, to, (void*)&cnt, bit_visitor_callback);
            // note that return code would be 0 in this case,
            // because range [15,500] has less than 5 set elements
            cout << "\nvisitor return code = " << res << endl; // -1
        }


        // CASE 3a: bm::for_each_bit_range
        {
            bit_visitor_functor func; // declare the visitor object
            int res = bm::for_each_bit_range(bv, from, to, func);
            cout << "\nvisitor return code = " << res << endl; // 0
        }


    }
    catch(std::exception& ex)
    {
        std::cerr << ex.what() << std::endl;
        return 1;
    }
    

    return 0;
}

