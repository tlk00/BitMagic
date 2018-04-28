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

/** \example svsample06.cpp
  Example how to search for an element.
 
  \sa bm::sparse_vector<>
*/

/*! \file svsample06.cpp
    \brief Example: sparse_vector<> scan search
*/


#include <iostream>
#include <vector>
#include <chrono>
#include <algorithm>
#include <random>
#include <stdexcept>

#include "bmsparsevec.h"
#include "bmsparsevec_algo.h"
#include "bmtimer.h"


using namespace std;

typedef bm::sparse_vector<bm::id_t, bm::bvector<> > sparse_vector_u32;


// ----------------------------------------------------
// Global parameters and types
// ----------------------------------------------------

const unsigned  value_max = 1250000;    // range of variants of events [0..max]
const unsigned  test_size = 250000000;  // vector size to generate

// -------------------------------------------
// Random generator
// -------------------------------------------

std::random_device rand_dev;
std::mt19937 gen(rand_dev());
std::uniform_int_distribution<> rand_dis(1, value_max); // generate uniform numebrs for [1, vector_max]

// timing storage for benchmarking
bm::chrono_taker::duration_map_type  timing_map;


// Function to generate test vector set with some NULL values stored as a
// separate bit-bector
//
static
void generate_test_set(std::vector<unsigned>& vect,
                       bm::bvector<>&         bv_null,
                       sparse_vector_u32&     sv
                       )
{
    vect.resize(test_size);
    bv_null.reset();
    
    for (unsigned i = 0; i < test_size; ++i)
    {
        unsigned v = rand_dis(gen);
        vect[i] = v;
        bv_null.set(i);
        
        sv.set(i, v);
        
        if (i % 64 == 0)
        {
            i += 5;  // insert a small NULL plate (unassigned values)
        }
    } // for
}


static
void vector_search(const std::vector<unsigned>& vect,
                   const bm::bvector<>&         bv_null,
                   unsigned                     value,
                   bm::bvector<>&               bv_res)
{
    bv_res.clear(true);
    for (size_t i = 0; i < vect.size(); ++i)
    {
        if (vect[i] == value)
        {
            bv_res.set_bit_no_check(i);
        }
    } // for
    bv_res &= bv_null;
}




static
void print_svector(const sparse_vector_u32& sv)
{
    if (sv.size() == 0)
    {
        cout << sv.size() << ": [ EMPTY ]" << endl;
        return;
    }
    cout << sv.size() << ": [ ";
    for (unsigned i = 0; i < sv.size(); ++i)
    {
        unsigned v = sv.at(i);
        bool is_null = sv.is_null(i);
        
        if (is_null)
            cout << "NULL";
        else
            cout << v << "";
        
        if (i == sv.size()-1)
            cout << " ]";
        else
            cout << ", ";
    }
    cout << endl;
}

inline
void print_bvector(const bm::bvector<>& bv)
{
    cout << "( count = " << bv.count() << ")" << ": [";
    
    bm::bvector<>::enumerator en = bv.first();
    for (; en.valid(); ++en)
    {
        cout << *en << ", ";
    }
    cout << "]" << endl;
}


int main(void)
{
    try
    {
        {
            sparse_vector_u32 sv(bm::use_null);
            
            // 2 -> 25
            // 3 -> 35
            // 7 -> 75
            // 1000 -> 2000
            // 256 -> 2001
            
            sv.set(2, 25);
            sv.set(3, 35);
            sv.set(7, 75);
            sv.set(1000, 2000);
            sv.set(256, 2001);

            bm::bvector<> bv_found;
            
            bm::sparse_vector_scan<sparse_vector_u32> scanner;
            scanner.find_eq(sv, 25, bv_found);
            
            print_bvector(bv_found);
        }
        
        std::vector<unsigned> vect;
        bm::bvector<> bv_null;
        bm::bvector<> bv_res;
        sparse_vector_u32 sv(bm::use_null);
        
        
        unsigned seach_repeats = 100;
        
        generate_test_set(vect, bv_null, sv);
        
        
        {
            bm::chrono_taker tt1("1. std::vector<> scan ", seach_repeats, &timing_map);
            
            for (unsigned i = 0; i < seach_repeats; ++i)
            {
                unsigned vs = rand_dis(gen);
                vector_search(vect, bv_null, vs, bv_res);
            } // for
        }

        {
            bm::sparse_vector_scan<sparse_vector_u32> scanner;

            bm::chrono_taker tt1("2. sparse_vector<> scan ", seach_repeats, &timing_map);
            
            for (unsigned i = 0; i < seach_repeats; ++i)
            {
                unsigned vs = rand_dis(gen);
                scanner.find_eq(sv, vs, bv_res);
            } // for
        }

        
        bm::chrono_taker::print_duration_map(timing_map, bm::chrono_taker::ct_ops_per_sec);

    }
    catch(std::exception& ex)
    {
        std::cerr << ex.what() << std::endl;
        return 1;
    }

    return 0;
}

