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
  Bucket sort using bit-vector and sparse vecctor
 
*/


#include <iostream>
#include <memory>
#include <map>
#include <vector>
#include <chrono>
#include <algorithm>
#include <random>
#include <stdexcept>

#include "bm.h"
#include "bmtimer.h"
#include "bmsparsevec.h"


// ----------------------------------------------------
// Global parameters and types
// ----------------------------------------------------

const unsigned  value_max = 250000;
const unsigned  test_size = 100000000;

std::random_device rand_dev;
std::mt19937 gen(rand_dev());
std::uniform_int_distribution<> rand_dis(1, value_max); // generate uniform numebrs for [1, vector_max]


typedef bm::sparse_vector<bm::id_t, bm::bvector<> > sparse_vector_u32;


// timing storage for benchmarking
bm::chrono_taker::duration_map_type  timing_map;


static
void bucket_sort(sparse_vector_u32& sv_out, const std::vector<unsigned>& vin)
{
    sparse_vector_u32::bvector_type* bv_null =
        const_cast<sparse_vector_u32::bvector_type*>(sv_out.get_null_bvector());
    
    for(auto v : vin)
    {
        if (!bv_null->test(v)) // first occurence of v
        {
            bv_null->set_bit_no_check(v); // fast set
        }
        else
        {
            unsigned v_cnt_prev = sv_out.get(v); // get previous occurence count
            ++v_cnt_prev;
            sv_out.set(v, v_cnt_prev);
        }
    } // for v
}

void print_sorted(const sparse_vector_u32& sv)
{
    const sparse_vector_u32::bvector_type* bv_null = sv.get_null_bvector();
    sparse_vector_u32::bvector_type::enumerator en = bv_null->first();
    
    for (; en.valid(); ++en)
    {
        unsigned v = *en;
        std::cout << v << ", ";
        unsigned cnt = sv.get(v);
        for (unsigned j = 0; j < cnt; ++j)
        {
            std::cout << v << ", ";
        } // for
    } // for en
    std::cout << std::endl;
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
            
            bucket_sort(r_sv, v);

            print_sorted(r_sv); // 1, 4, 5, 8, 8, 8, 10,
        }
        
        // run benchmark
        //
        std::vector<unsigned> v;
        
        // generate vector of random numbers
        for (unsigned i = 0; i < test_size; ++i)
        {
            v.push_back(rand_dis(gen));
        }
        std::cout << "test vector generation ok" << std::endl;

        
        {
            sparse_vector_u32 r_sv(bm::use_null);
            {
            bm::chrono_taker tt1("1. bucket sort", 1, &timing_map);

            bucket_sort(r_sv, v);
            }
        }
        
        {
            bm::chrono_taker tt1("2. std::sort", 1, &timing_map);
            std::sort(v.begin(), v.end());
        }
        
        std::cout << "                                                        "
                  << std::endl;
        
        bm::chrono_taker::print_duration_map(timing_map, bm::chrono_taker::ct_ops_per_sec);

    }
    catch(std::exception& ex)
    {
        std::cerr << ex.what() << std::endl;
        return 1;
    }

    return 0;
}

