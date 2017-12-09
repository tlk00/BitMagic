/*
Copyright(c) 2002-2017 Anatoliy Kuznetsov(anatoliy_kuznetsov at yahoo.com)

Permission is hereby granted, free of charge, to any person
obtaining a copy of this software and associated documentation
files (the "Software"), to deal in the Software without restriction,
including without limitation the rights to use, copy, modify, merge,
publish, distribute, sublicense, and/or sell copies of the Software,
and to permit persons to whom the Software is furnished to do so,
subject to the following conditions:

The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
OTHER DEALINGS IN THE SOFTWARE.

You have to explicitly mention BitMagic project in any derivative product,
its WEB Site, published materials, articles or any other work derived from this
project or based on our code or know-how.

For more information please visit:  http://bitmagic.io

*/

/** \example sample11.cpp
  Example of how to use various bit counting techniques

  \sa bm::bvector<>::count() 
  \sa bm::bvector<>::count_range()
  \sa bm::bvector<>::count_to() 
  \sa bm::count_and() 
 */


#include <iostream>

#include "bm.h"
#include "bmalgo.h"
#include "bmtimer.h"

using namespace std;

// timing storage for benchmarking
bm::chrono_taker::duration_map_type  timing_map;

const unsigned benchmark_count = 10000;
unsigned       vector_max = 40000000;

/// generate pseudo-random bit-vector for testing
///
void generate_bvector(bm::bvector<>& bv)
{
    unsigned i;
    for (i = 0; i < 30000; i += 10)
    {
        bv.set(i);
    }
    for (i = 300000; i < vector_max; i += 100)
    {
        bv.set(i);
    }
}


/// simple population count for the whole vector
///
void bv_count_test(const bm::bvector<>& bv)
{
    bm::chrono_taker tt1("1. bvector<>::count()", benchmark_count, &timing_map);

    unsigned cnt = 0;
    for (unsigned i = 0; i < benchmark_count; ++i)
    {
        cnt += bv.count();
    }
    // this is mostly to prevent compiler to optimize loop away
    std::cout << "Count test finished." << cnt << std::endl;
}

/// count_range() test
///
void bv_count_range(const bm::bvector<>& bv)
{
    bm::chrono_taker tt1("2. bvector<>::count_range()", benchmark_count, &timing_map);

    unsigned cnt = 0;
    for (unsigned i = 0; i < benchmark_count; ++i)
    {
        unsigned from = rand() % vector_max;
        unsigned to = rand() % vector_max;
        if (from > to)
            swap(from, to);
        cnt += bv.count_range(from, to);
    }
    // this is mostly to prevent compiler to optimize loop away
    std::cout << "Count range test finished." << cnt << std::endl;
}

/// count_range() test using pre-calculated blocks bit count
///
void bv_count_range_acc(const bm::bvector<>& bv)
{
    // build a block population count list, used for count_range() acceleration
    unsigned  blocks_cnt[bm::set_total_blocks];
    bv.count_blocks(blocks_cnt);

    bm::chrono_taker tt1("3. bvector<>::count_range() with blocks list", benchmark_count, &timing_map);

    unsigned cnt = 0;
    for (unsigned i = 0; i < benchmark_count; ++i)
    {
        unsigned from = rand() % vector_max;
        unsigned to = rand() % vector_max;
        if (from > to)
            swap(from, to);
        cnt += bv.count_range(from, to, blocks_cnt); // use blocks count for acceleration
    }
    // this is mostly to prevent compiler to optimize loop away
    std::cout << "Count range with blocks test finished." << cnt << std::endl;
}

/// count_to() test using pre-calculated blocks bit count
///
void bv_count_to_acc(const bm::bvector<>& bv)
{
    // build a block population count list, used for count_to() acceleration
    bm::bvector<>::blocks_count bc;
    bv.running_count_blocks(&bc);

    bm::chrono_taker tt1("4. bvector<>::count_to() with blocks list", benchmark_count, &timing_map);

    unsigned cnt = 0;
    for (unsigned i = 0; i < benchmark_count; ++i)
    {
        unsigned to = rand() % vector_max;
        cnt += bv.count_to(to, bc); // use blocks count for acceleration
    }
    // this is mostly to prevent compiler to optimize loop away
    std::cout << "Count to with blocks test finished." << cnt << std::endl;
}


/// count_range implemented via two count_to() calls using pre-calculated running count
///
void bv_count_to_range_acc(const bm::bvector<>& bv)
{
    // build a block population count list, used for count_to() acceleration
    bm::bvector<>::blocks_count bc;
    bv.running_count_blocks(&bc);

    bm::chrono_taker tt1("5. bvector<>::count_to to simulate count_range()", benchmark_count, &timing_map);

    unsigned cnt = 0;
    for (unsigned i = 0; i < benchmark_count; ++i)
    {
        unsigned from = 1 + rand() % vector_max;
        unsigned to = rand() % vector_max;
        if (from > to)
            swap(from, to);
        
        unsigned cnt_to = bv.count_to(to, bc);
        unsigned cnt_from = bv.count_to(from - 1, bc);
        unsigned cnt_r = cnt_to - cnt_from;
        cnt += cnt_r;
    }
    // this is mostly to prevent compiler to optimize loop away
    std::cout << "Count range via count_to test finished." << cnt << std::endl;
}

/// count_range implemented via bm::count_and
///
/// this method can be used, when we need co compute multiple ranges in one call
///
void bv_count_and(const bm::bvector<>& bv)
{
    bm::chrono_taker tt1("6. bm::count_and with mask vector", benchmark_count, &timing_map);

    unsigned cnt = 0;
    bm::bvector<> mask_bv(bm::BM_GAP);

    for (unsigned i = 0; i < benchmark_count; ++i)
    {
        unsigned from = 1 + rand() % vector_max;
        unsigned to = rand() % vector_max;
        if (from > to)
            swap(from, to);

        mask_bv.set_range(from, to, true); // set mask vector

        cnt += bm::count_and(bv, mask_bv);
        mask_bv.clear(true); // clear and free memory (faster)
    }
    // this is mostly to prevent compiler to optimize loop away
    std::cout << "count AND finished." << cnt << std::endl;
}



int main(void)
{
    try
    {
        bm::bvector<>   bv;
        generate_bvector(bv);

        // Test 1.
        // Uses plain bvector<>::count() to compute global population count
        // This function would benefit from SIMD (SSE42 / AVX2) acceleration
        //
        bv_count_test(bv);

        // Test 2.
        // Uses bvector<>::count_range() to compute population count in a randomly generated
        // region of a bit-vector.
        // This is should be naturally faster than Test 1, because it range is less than the whole
        //
        bv_count_range(bv);

        // Test 3.
        // Uses bvector<>::count_range() together with bvector<>::count_blocks() 
        // (pre-calculated bit-count for each block).
        // It make sense to use this method if bit-vector is constant (or chnages infrequently)
        // and we need to do many range counting calculations
        //
        bv_count_range_acc(bv);

        // Test 4.
        // Uses bvector<>::count_to() to compute population count to a specified element.
        // Equivalent of count_range(0, to);
        // This method uses acceleration structure using bvector<>::running_count_blocks()
        // It is similar to count_range acceleration, but uses a different (faster) algorithm
        //
        bv_count_to_acc(bv);

        // Test 5.
        // Uses bvector<>::count_to() twice to simulate count_range()
        // using counting difference:
        // count_r = count_to(0, from) - count_to(0, to-1)
        // This method can actually be faster than count_range()
        //
        bv_count_to_range_acc(bv);

        // Test 6.
        // Compute range population count via a mask vector and logical AND operation.
        // Not the fastest method, but can be useful, when multiple ranges needs to be computed
        //
        bv_count_and(bv);

        // print all test timing results
        bm::chrono_taker::print_duration_map(timing_map, bm::chrono_taker::ct_ops_per_sec);

     }
    catch(std::exception& ex)
    {
        std::cerr << ex.what() << std::endl;
        return 1;
    }

    return 0;
}

