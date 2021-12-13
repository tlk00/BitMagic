/*
Copyright(c) 2018 Anatoliy Kuznetsov(anatoliy_kuznetsov at yahoo.com)

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

/** \example xsample02.cpp
  Counting sort using bit-vector and sparse vector to build histogram of unsigned ints.
  Benchmark compares different histogram buiding techniques using BitMagic and std::sort()
 
  Histogram construction, based on integer events is a common problem,
  this demo studies different approaches, potential for parallelization and other
  aspects.
*/

/*! \file xsample02.cpp
    \brief Example: sparse_vector<> used for counting sort / historgam construction
*/


#include <iostream>
#include <memory>
#include <map>
#include <vector>
#include <chrono>
#include <algorithm>
#include <random>
#include <stdexcept>

#include <future>
#include <thread>

using namespace std;

#include "bm.h"
#include "bmtimer.h"
#include "bmsparsevec.h"
#include "bmundef.h" /* clear the pre-proc defines from BM */

// ----------------------------------------------------
// Global parameters and types
// ----------------------------------------------------

const unsigned  value_max = 1250000;    // range of variants of events [0..max]
const unsigned  test_size = 250000000;  // number of events (ints) to generate

// -------------------------------------------
// Random generator
// -------------------------------------------

std::random_device rand_dev;
std::mt19937 gen(rand_dev());
std::uniform_int_distribution<> rand_dis(1, value_max); // generate uniform numebrs for [1, vector_max]


typedef bm::sparse_vector<unsigned, bm::bvector<> > sparse_vector_u32;
typedef std::map<unsigned, unsigned>                map_u32;


// timing storage for benchmarking
bm::chrono_taker<>::duration_map_type  timing_map;


// -------------------------------------------
// Counting sort / histogram construction (std::map)
// -------------------------------------------

static
void sort_map(map_u32& hmap, const std::vector<unsigned>& vin)
{
    for (auto v : vin)
    {
        hmap[v]++;
    }
}


// -------------------------------------------
// Counting sort / histogram construction (naive)
// -------------------------------------------

// This sorting method uses sparse_vector<> as a storage but implements increment
// as an get-inc-put operations (decoding-encoding every value in the sum)
//
static
void counting_sort_naive(sparse_vector_u32& sv_out, const std::vector<unsigned>& vin)
{
    for (auto v : vin)
    {
        auto count = sv_out.get(v);
        sv_out.set(v, count + 1);
    }
}

// -------------------------------------------
// Counting sort / histogram construction
// -------------------------------------------

// This sorting method uses sparse_vector<> as a storage but implements increment
// but increment was implemented as a bm::sparse_vector::inc() method
// which in turn is based on uses bvector<>::inc()
//
// This approach is faster than decoding-encoding used in naive counting sort
//
static
void counting_sort(sparse_vector_u32& sv_out, const std::vector<unsigned>& vin)
{
    for(auto v : vin)
        sv_out.inc(v);
}

// --------------------------------------------------
// Counting sort / histogram construction (parallel)
// --------------------------------------------------

// parallel subproblem for all even numbers: (v & 1) == 0
inline 
unsigned counting_sort_subbatch(sparse_vector_u32* sv_out, const std::vector<unsigned>* vin)
{
    for (size_t i = 0; i < vin->size(); i++)
    {
        auto v = (*vin)[i];
        if ((v & 1) == 0)
            sv_out->inc(v);
    }
    return 0;
}

// Parallel histogram construction uses a very simple divide and conquer technique
// splitting by even/odd numbers, uses std::async() for parallelization
//
// (should be possible to do a lot better than that)
//
static
void counting_sort_parallel(sparse_vector_u32& sv_out, const std::vector<unsigned>& vin)
{
    sparse_vector_u32 sv_out2(bm::use_null);
    // process evens in parallel
    std::future<unsigned> f1 = std::async(std::launch::async, counting_sort_subbatch, &sv_out2, &vin);

    // process all odd elements
    for (size_t i = 0; i < vin.size(); i++)
    {
        auto v = vin[i];
        if (v & 1)
            sv_out.inc(v);
    }
    f1.wait();
    
    // merge effectively performs logical OR on all plains, only it
    // borrows memory blocks from the argument vector, so it gets changed
    // (which is ok here, since sv_out2 is a temporary)
    //
    sv_out.merge(sv_out2);
}

// Test utility. It also illustrates histogram access method
//
static
void print_sorted(const sparse_vector_u32& sv)
{
    const sparse_vector_u32::bvector_type* bv_null = sv.get_null_bvector();
    sparse_vector_u32::bvector_type::enumerator en = bv_null->first();
    
    for (; en.valid(); ++en)
    {
        unsigned v = *en;
        unsigned cnt = sv.get(v);
        for (unsigned j = 0; j < cnt; ++j)
        {
            std::cout << v << ", ";
        } // for
    } // for en
    std::cout << std::endl;
}

// Test utility for std::map
//
static
void print_sorted(const map_u32& hmap)
{
    map_u32::const_iterator it = hmap.begin();
    map_u32::const_iterator it_end = hmap.end();
    
    for (; it != it_end; ++it)
    {
        unsigned v = it->first;
        unsigned cnt = it->second;
        for (unsigned j = 0; j < cnt; ++j)
        {
            std::cout << v << ", ";
        } // for
    } // for en
    std::cout << std::endl;
}


// build histogram using sorted vector
//
static
void build_histogram(sparse_vector_u32& sv_out, const std::vector<unsigned>& vin)
{
    if (vin.empty())
        return;
    unsigned start = vin[0];
    unsigned count = 0; // histogram counter
    for (auto v : vin)
    {
        if (v == start)
        {
            ++count;
        }
        else
        {
            sv_out.set(start, count);
            start = v; count = 1; 
        }
    }
    if (count)
    {
        sv_out.set(start, count);
    }
}




int main(void)
{
    try
    {
        // try simple input vector as a model
        //
        {
            std::vector<unsigned> v {10, 1, 5, 4, 8, 8, 8} ;
            sparse_vector_u32 r_sv(bm::use_null);  // result vector
            
            counting_sort(r_sv, v);

            print_sorted(r_sv); // 1, 4, 5, 8, 8, 8, 10,
            
            sparse_vector_u32 p_sv(bm::use_null);
            counting_sort_parallel(p_sv, v);
            print_sorted(r_sv);
            
            map_u32  h_map;
            sort_map(h_map, v);
            print_sorted(h_map);


            std::sort(v.begin(), v.end());
            sparse_vector_u32 h_sv(bm::use_null);  // histogram vector
            build_histogram(h_sv, v);
            if (!r_sv.equal(h_sv))
            {
                std::cerr << "Error: Histogram comparison failed!" << std::endl;
                print_sorted(h_sv);
                return 1;
            }

        }
        
        // run benchmarks
        //
        std::vector<unsigned> v;
        
        // generate vector of random numbers
        for (unsigned i = 0; i < test_size; ++i)
        {
            v.push_back(unsigned(rand_dis(gen)));
        }
        std::cout << "test vector generation ok" << std::endl;

        sparse_vector_u32 r_sv(bm::use_null);
        sparse_vector_u32 h_sv(bm::use_null);
        sparse_vector_u32 n_sv(bm::use_null);
        sparse_vector_u32 p_sv(bm::use_null);
        map_u32  h_map;

        {
            bm::chrono_taker tt1(cout, "1. counting sort ", 1, &timing_map);
            counting_sort(r_sv, v);
        }

        {
            bm::chrono_taker tt1(cout, "3. counting sort (naive) ", 1, &timing_map);
            counting_sort_naive(n_sv, v);
        }

        {
            bm::chrono_taker tt1(cout, "4. counting sort (parallel) ", 1, &timing_map);
            counting_sort_parallel(p_sv, v);
        }

        {
            bm::chrono_taker tt1(cout, "5. counting sort (map) ", 1, &timing_map);
            sort_map(h_map, v);
        }
        
        {
            bm::chrono_taker tt1(cout, "2. std::sort() + histogram", 1, &timing_map);
            std::sort(v.begin(), v.end());
            build_histogram(h_sv, v);
        }


        // quality assurance checks
        //
        if (!r_sv.equal(h_sv) || !n_sv.equal(h_sv)) 
        {
            std::cerr << "Error: Histogram comparison failed!" << std::endl;
            return 1;
        }
        if (!r_sv.equal(p_sv))
        {
            std::cerr << "Error: Histogram comparison failed for parallel sort!" << std::endl;
            return 1;
        }

        // compute memory consumption of sparse vector
        {
            std::cout << std::endl;

            BM_DECLARE_TEMP_BLOCK(tb);
            sparse_vector_u32::statistics st;
            r_sv.optimize(tb, sparse_vector_u32::bvector_type::opt_compress, &st);
            
            std::cout << "Sparse vector memory usage:" << st.memory_used / (1024*1024)<< "MB" << std::endl;
            std::cout << "vector<unsigned> usage:" << v.size() * sizeof(v[0]) / (1024 * 1024) << "MB" << std::endl << std::endl;
        }


        bm::chrono_taker<>::print_duration_map(cout, timing_map, bm::chrono_taker<>::ct_ops_per_sec);

    }
    catch(std::exception& ex)
    {
        std::cerr << ex.what() << std::endl;
        return 1;
    }

    return 0;
}

