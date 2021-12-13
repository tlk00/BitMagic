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

/** \example sample12.cpp
  Example of how to use various bit setting techniques. 
  Several techniques benchmarked to better illustrate relative performance.

  \sa bm::bvector<>::set() 
  \sa bm::bvector<>::set_bit()
  \sa bm::bvector<>::set_bit_conditional()
  \sa bm::bvector<>::set_range()
  \sa bm::bvector<>::clear_bit()
  \sa bm::bvector<>::reset()
  \sa bm::bvector<>::flip()
  \sa bm::bvector<>::swap()
  \sa bm::bvector<>::extract_next()
  \sa bm::bvector<>::set_bit_no_check()
  \sa bm::combine_or()
 */

/*! \file sample12.cpp
    \brief Example: bvector<> analysis of bit setting methods
*/

#include <iostream>
#include <vector>

#define BM64ADDR

#include "bm.h"
#include "bmalgo.h"
#include "bmtimer.h"
#include "bmundef.h" /* clear the pre-proc defines from BM */

using namespace std;

// timing storage for benchmarking
bm::chrono_taker<>::duration_map_type  timing_map;

// depending of the build it may be "unsigned int" or 64-bit "unsigned long long"
typedef bm::bvector<>::size_type bm_size_type;

const unsigned benchmark_count = 1000;
bm_size_type       vector_max = 4000000;


// Utility template function used to print container
template<class T> void PrintContainer(T first, T last)
{
    if (first == last)
        std::cout << "<EMPTY SET>";
    else
        for (; first != last; ++first)
            std::cout << *first << ";";
    std::cout << std::endl;
}


static
void generate_test_vectors(std::vector<bm_size_type> &v1,
                           std::vector<bm_size_type> &v2,
                           std::vector<bm_size_type> &v3)
{
    bm_size_type j;
    for (j = 0; j < vector_max; j += 2)
    {
        v1.push_back(j);
    }
    for (j = 0; j < vector_max; j += 5)
    {
        v2.push_back(j);
    }
    for (j = 0; j < vector_max; j += 120)
    {
        v3.push_back(j);
    }
}


// stress test for bm::bvector<>::set_bit() 
//
static
void bv_set_bit_test()
{
    bm::chrono_taker tt1(cout, "1. bvector<>::set_bit()", benchmark_count, &timing_map);
    bm::bvector<> bv1, bv2, bv3;

    for (unsigned i = 0; i < benchmark_count; ++i)
    {
        bm_size_type j;
        for (j = 0; j < vector_max; j += 2)
        {
            bv1.set_bit(j, true);
        }
        for (j = 0; j < vector_max; j += 10)
        {
            bv2.set_bit(j, true);
        }
        for (j = 0; j < vector_max; j += 120)
        {
            bv3.set_bit(j, true);
        }
        bv1.reset();
        bv2.reset();
        bv3.reset();
    } // for
}


// stress test for bm::bvector<>::set_bit() 
//
static
void bv_set_bit_no_check_test()
{
    bm::chrono_taker tt1(cout, "2. bvector<>::set_bit_no_check()", benchmark_count, &timing_map);
    bm::bvector<> bv1, bv2, bv3;

    bv1.init();
    bv2.init();
    bv3.init();

    for (unsigned i = 0; i < benchmark_count; ++i)
    {
        bm::id_t j;
        for (j = 0; j < vector_max; j += 2)
        {
            bv1.set_bit_no_check(j);
        }
        for (j = 0; j < vector_max; j += 10)
        {
            bv2.set_bit_no_check(j);
        }
        for (j = 0; j < vector_max; j += 120)
        {
            bv3.set_bit_no_check(j);
        }
    }
}

static
void combine_or_test(std::vector<bm_size_type> &v1,
                     std::vector<bm_size_type> &v2,
                     std::vector<bm_size_type> &v3)
{
    bm::chrono_taker tt1(cout, "3. combine_or()", benchmark_count, &timing_map);
    bm::bvector<> bv1, bv2, bv3;

    for (unsigned i = 0; i < benchmark_count; ++i)
    {
        bm::combine_or(bv1, v1.begin(), v1.end());
        bm::combine_or(bv2, v2.begin(), v2.end());
        bm::combine_or(bv3, v3.begin(), v3.end());
    }
}

static
void bvector_bulk_set_test(std::vector<bm_size_type> &v1,
                           std::vector<bm_size_type> &v2,
                           std::vector<bm_size_type> &v3)
{
    bm::chrono_taker tt1(cout, "3. bvector<>::set() array", benchmark_count, &timing_map);
    bm::bvector<> bv1, bv2, bv3;

    for (unsigned i = 0; i < benchmark_count; ++i)
    {
        bv1.set(&v1[0], bm_size_type(v1.size()));
        bv2.set(&v2[0], bm_size_type(v2.size()));
        bv3.set(&v3[0], bm_size_type(v3.size()));
    }
}



int main(void)
{
    try
    {
        // 0. create bvector, use brace initialization to add some initial data
        //
        bm::bvector<>   bv1 { 2, 3, 4 };
        
        PrintContainer(bv1.first(), bv1.end()); // 2, 3, 4
        bv1.clear();

        // 1. Set some bits using regular bvector<>::set() method
        //
        bv1.set(10);
        bv1.set(256);
        bv1.set(1000000);

        PrintContainer(bv1.first(), bv1.end()); // 10, 256, 1000000

        // 2. now use bvector<>::set_bit()
        // it returns a report if target actually changed
        //

        bm::id_t bits[] = { 256, 512, 10 };
        unsigned cnt = 0;

        for (unsigned i = 0; i < sizeof(bits) / sizeof(bits[0]); ++i)
        {
            bool b = bv1.set_bit(bits[i], true);
            cnt += b;
        }
        std::cout << "Number of bits changed:" << cnt << std::endl;
        PrintContainer(bv1.first(), bv1.end());

        // 3. set and clear some bits using bvector<>::set_bit_conditional()
        // method sets bit n only if current value equals the condition
        //
        bool b;
        b = bv1.set_bit_conditional(5, true, false); // set bit 5 to true if it is false (yes)
        std::cout << "Bit 5 set:" << (b ? " yes " : " no ") << std::endl; // (yes)
        
        b = bv1.set_bit_conditional(256, true, false); // set bit 256 to true if it is false
        std::cout << "Bit 256 set:" << (b ? " yes " : " no ") << std::endl; // (no)
        
        b = bv1.set_bit_conditional(256, true, true); // set bit 256 to true if it is true 
        std::cout << "Bit 256 set:" << (b ? " yes " : " no ") << std::endl; // (no)

        PrintContainer(bv1.first(), bv1.end());

        // 4. set or clear multiple bits using bvector::set_range()
        // This method is faster than calling bvector::set() many times in a row
        //
        bv1.set_range(10, 15, true); // set all bits in [10..15] closed interval 
        PrintContainer(bv1.first(), bv1.end());

        bv1.set_range(10, 12, false); // clear all bits in [10..15] closed interval 
        PrintContainer(bv1.first(), bv1.end());

        // 5. bvector::clear_bit() - same as set_bit() just a syntax sugar
        //
        b = bv1.clear_bit(13);
        std::cout << "Bit 13 set:" << (b ? " yes " : " no ") << std::endl; // (yes)

        // 6. bvector<>::reset() - clears all the bits, frees the blocks memory
        // same as bvector<>::clear(true);
        //
        bv1.reset();
        PrintContainer(bv1.first(), bv1.end());  // <EMPTY>

        // 7. use bm::combine_or() to set bits
        //
        bm::combine_or(bv1, &bits[0], &bits[0] + (sizeof(bits) / sizeof(bits[0])));
        PrintContainer(bv1.first(), bv1.end()); // 10, 256, 512

        // 8. use bvector<>::flip( ) to flip a bit
        //
        bv1.flip(256); 
        bv1.flip(257);
        PrintContainer(bv1.first(), bv1.end()); // 10, 257, 512

        // 9. bvector<>::swap() to flip content of two bit-vectors
        //
        bm::bvector<> bv2;
        bv1.swap(bv2);
        PrintContainer(bv1.first(), bv1.end());  // <EMPTY>
        PrintContainer(bv2.first(), bv2.end()); // 10, 257, 512

        // 10. use bvector<>::extract_next() to find ON bit and turn it to 0
        // this function is useful for building FIFO queues on bvector
        //
        bm_size_type p = 1;
        for (p = bv2.extract_next(p); p != 0; p = bv2.extract_next(p))
        {
            std::cout << "Extracted p = " << p << std::endl;
        } 
        PrintContainer(bv2.first(), bv2.end()); // <EMPTY>

        // 11. use bvector<>::set_bit_no_check() to set bits faster
        //
        bm::bvector<> bv3;
        bv3.init();   // the key here you MUST call init() before setting bits this way
                      
        bv3.set_bit_no_check(10);
        bv3.set_bit_no_check(100);
        bv3.set_bit_no_check(1000);

        PrintContainer(bv3.first(), bv3.end()); // 10, 100, 1000


        std::vector<bm_size_type> v1, v2, v3;
        generate_test_vectors(v1, v2, v3);

        for (unsigned k = 0; k < 1000; ++k)
        {
            // run a few CPU "pre-heat" loops to get consistent results
            bm_size_type s = 0;
            bm_size_type i;
            for (i = 0; i < v1.size(); ++i)
                s = s + s* v1[i] + i;
            for (i = 0; i < v2.size(); ++i)
                s = s + s* v2[i] + i;
            std::cout << s << "\r";
        }

        std::cout << std::endl << "Running benchmarks..." << std::endl;

        bv_set_bit_test();
        bv_set_bit_no_check_test();
        combine_or_test(v1, v2, v3);
        bvector_bulk_set_test(v1, v2, v3);

        // print all test timing results
        //
        std::cout << std::endl;
        bm::chrono_taker<>::print_duration_map(cout, timing_map, bm::chrono_taker<>::ct_all);
    }
    catch(std::exception& ex)
    {
        std::cerr << ex.what() << std::endl;
        return 1;
    }

    return 0;
}

