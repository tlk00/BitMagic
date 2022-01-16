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

/** \example sample11.cpp
  Example of how to use various bit counting techniques

  \sa bm::bvector::count()
  \sa bm::bvector::count_range()
  \sa bm::bvector::count_to()
  \sa bm::bvector::count_to_test()
  \sa bm::count_and()
  \sa bm::bvector::counted_enumerator
 */

/*! \file sample11.cpp
    \brief Example: bvector<> bit-counting techniques analysis
*/

#include <iostream>
#include <random>
#include <memory>

#include "bm.h"
#include "bmalgo.h"
#include "bmtimer.h"
#include "bmundef.h" /* clear the pre-proc defines from BM */

using namespace std;

// timing storage for benchmarking
bm::chrono_taker<>::duration_map_type  timing_map;

const unsigned benchmark_count = 10000;
unsigned            vector_max = 400000000;

std::random_device rand_dev;  
std::mt19937 gen(rand_dev()); // mersenne_twister_engine 
std::uniform_int_distribution<> rand_dis(1, int(vector_max)); // generate uniform numebrs for [1, vector_max]


/// generate pseudo-random bit-vector, mix of blocks
/// 
static
void generate_bvector(bm::bvector<>& bv)
{
    bm::bvector<>::size_type i, j;
    for (i = 0; i < vector_max;)
    {
        // generate bit-blocks
        for (j = 0; j < 65535*8; i += 10, j++)
        {
            bv.set(i);
        }
        if (i > vector_max)
            break;
        // generate GAP (compressed) blocks
        for (j = 0; j < 65535; i += 120, j++)
        {
            unsigned len = rand() % 64;
            bv.set_range(i, i + len);
            i += len;
            if (i > vector_max)
                break;
        }
    }

    // compress vector
    BM_DECLARE_TEMP_BLOCK(tb)
    bv.optimize(tb);

    // compute bit-vector statistics
    bm::bvector<>::statistics st;
    bv.calc_stat(&st);

    std::cout << "Bit-vector statistics: GAP (compressed blocks)=" << st.gap_blocks
              << ", BIT (uncompressed blocks)=" << st.bit_blocks
              << std::endl << std::endl;
}

/// "pre-heat" CPU to minimize dynamic overclocking effects
///
static
bm::bvector<>::size_type pre_heat(const bm::bvector<>& bv)
{
    bm::bvector<>::size_type cnt = 0;
    bm::bvector<>::size_type m = 1;
    for (unsigned i = 0; i < benchmark_count; ++i)
    {
        cnt += bv.count();
        m+=cnt*cnt;
    }
    return m;
}



/// simple population count for the whole vector
///
static
void bv_count_test(const bm::bvector<>& bv)
{
    bm::bvector<>::size_type cnt = 0;

    {
        bm::chrono_taker tt1(cout, "1. bvector<>::count()", benchmark_count / 2, &timing_map);
        for (unsigned i = 0; i < benchmark_count / 2; ++i)
        {
            cnt += bv.count();
        }
    }
    // this is mostly to prevent compiler to optimize loop away
    std::cout << "Count test finished." << cnt << "\r";
}

/// count_range() test
///
static
void bv_count_range(const bm::bvector<>& bv)
{
    bm::bvector<>::size_type cnt = 0;
    {
        bm::chrono_taker tt1(cout, "2. bvector<>::count_range()", benchmark_count, &timing_map);
        for (unsigned i = 0; i < benchmark_count; ++i)
        {
            unsigned from = unsigned(rand_dis(gen));
            unsigned to = unsigned(rand_dis(gen));
            if (from > to)
                swap(from, to);
            cnt += bv.count_range(from, to);
        }
    }
    // this is mostly to prevent compiler to optimize loop away
    std::cout << "Count range test finished." << cnt << "\r";
}

/// count_range() test using pre-calculated blocks bit count
///
static
void bv_count_range_acc(const bm::bvector<>& bv)
{
    bm::bvector<>::size_type cnt = 0;
    
    std::unique_ptr<bm::bvector<>::rs_index_type> rs(new bm::bvector<>::rs_index_type());
    bv.build_rs_index(rs.get());

    {
        bm::chrono_taker tt1(cout, "3. bvector<>::count_range() with rs_index", benchmark_count, &timing_map);
        cnt = 0;
        for (unsigned i = 0; i < benchmark_count; ++i)
        {
            unsigned from = unsigned(rand_dis(gen));
            unsigned to = unsigned(rand_dis(gen));
            if (from > to)
                swap(from, to);
            cnt += bv.count_range(from, to, *rs); // use rs index for acceleration
        } // for i
    }
    // this is mostly to prevent compiler to optimize loop away
    std::cout << "Count range with blocks test finished." << cnt << "\r";
}

/// count_to() test using pre-calculated rank-select index
///
static
void bv_count_to_acc(const bm::bvector<>& bv)
{
    bm::bvector<>::size_type cnt = 0;
    
    // build a block population count list, used for count_to() acceleration
    std::unique_ptr<bm::bvector<>::rs_index_type> rs(new bm::bvector<>::rs_index_type());
    bv.build_rs_index(rs.get());

    {
        bm::chrono_taker tt1(cout, "4. bvector<>::count_to() with rs_index", benchmark_count, &timing_map);

        for (unsigned i = 0; i < benchmark_count; ++i)
        {
            unsigned to = unsigned(rand_dis(gen));
            cnt += bv.count_to(to, *rs); // use rank-select index for acceleration
        }
    }
    // this is mostly to prevent compiler to optimize loop away
    std::cout << "Count to with blocks test finished." << cnt << "\r";
}


/// count_range implemented via two count_to() calls using pre-calculated
/// rank-select index
///
static
void bv_count_to_range_acc(const bm::bvector<>& bv)
{
    bm::bvector<>::size_type cnt = 0;
    
    // build a block population count list, used for count_to() acceleration
    std::unique_ptr<bm::bvector<>::rs_index_type> rs(new bm::bvector<>::rs_index_type());
    bv.build_rs_index(rs.get());

    {
        bm::chrono_taker tt1(cout, "5. bvector<>::count_to to simulate count_range()", benchmark_count, &timing_map);

        for (unsigned i = 0; i < benchmark_count; ++i)
        {
            unsigned from = unsigned(rand_dis(gen));
            unsigned to = unsigned(rand_dis(gen));
            if (from > to)
                swap(from, to);
            
            bm::bvector<>::size_type cnt_to = bv.count_to(to, *rs);
            bm::bvector<>::size_type cnt_from = bv.count_to(from - 1, *rs);
            bm::bvector<>::size_type cnt_r = cnt_to - cnt_from;
            cnt += cnt_r;
        }
    }
    // this is mostly to prevent compiler to optimize loop away
    std::cout << "Count range via count_to test finished." << cnt << "\r";
}

/// count_range implemented via bm::count_and
///
/// this method can be used, when we need co compute multiple ranges in one call
///
static
void bv_count_and(const bm::bvector<>& bv)
{
    bm::bvector<>::size_type cnt = 0;
    {
        bm::chrono_taker tt1(cout, "6. bm::count_and with mask vector", benchmark_count, &timing_map);

        bm::bvector<> mask_bv(bm::BM_GAP); // use compressed mask, better seluts on long ranges
        for (unsigned i = 0; i < benchmark_count; ++i)
        {
            unsigned from = unsigned(rand_dis(gen));
            unsigned to = unsigned(rand_dis(gen));
            if (from > to)
                swap(from, to);

            mask_bv.set_range(from, to, true); // set mask vector

            cnt += bm::count_and(bv, mask_bv);
            mask_bv.clear(true); // clear and free memory (faster)
        }
    }
    // this is mostly to prevent compiler to optimize loop away
    std::cout << "count AND finished." << cnt << "\r";
}

/// count_to implemented via bm::bvector<>::counted_enumerator
///
/// Counted enumerator is an iterator automata, which counts the running population count 
/// along the iteration sequence
///
static
void bv_counted_enumerator(const bm::bvector<>& bv)
{
    bm::bvector<>::size_type cnt = 0;
    {
        // This is a slow method so we use less iterators
        bm::chrono_taker tt1(cout, "7. bm::bvector<>::counted_enumerator", benchmark_count/20, &timing_map);

        for (unsigned i = 0; i < benchmark_count/20; ++i)
        {
            unsigned to = unsigned(rand_dis(gen));
            bm::bvector<>::counted_enumerator en = bv.first();
            for (; en.valid(); ++en)
            {
                if (*en > to)
                    break;
            }
            cnt += en.count();
        }
    }
    std::cout << "counted_enumerator finished." << cnt << "\r";
}




int main(void)
{
    try
    {
        bm::bvector<>   bv;
        generate_bvector(bv);
        
        /// pre-heat CPU to minimize dynamic overclocking
        unsigned s = pre_heat(bv);
        std::cout << s << "\r";


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

        // Test 7.
        // Compute cout using counted_enumerator iterator
        // method combines iteratrion over bit vector and sliding population count
        bv_counted_enumerator(bv);


        // print all test timing results
        //
        std::cout << "                                                        "
                  << std::endl;
                  
        bm::chrono_taker<>::print_duration_map(cout, timing_map, bm::chrono_taker<>::ct_ops_per_sec);
    }
    catch(std::exception& ex)
    {
        std::cerr << ex.what() << std::endl;
        return 1;
    }

    return 0;
}

