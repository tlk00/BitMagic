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
   Search/scan for elements in unordered, non-unique sparse vector
 
  \sa bm::sparse_vector<>::const_iterator
  \sa bm::sparse_vector<>::back_insert_iterator
  \sa bm::sparse_vector_scanner
*/

/*! \file svsample06.cpp
    \brief Example: sparse_vector<> scan search (non-ordered set functionality)
*/

#include <iostream>
#include <vector>
#include <chrono>
#include <algorithm>
#include <random>
#include <stdexcept>

#include "bm.h"
#include "bmsparsevec.h"
#include "bmsparsevec_algo.h"
#include "bmtimer.h"
#include "bmundef.h" /* clear the pre-proc defines from BM */

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
bm::chrono_taker<>::duration_map_type  timing_map;



// Function to generate test vector set with some NULL values stored as a
// separate bit-bector
//
static
void generate_test_set(std::vector<unsigned>& vect,
                        bm::bvector<>&         bv_null,
                        sparse_vector_u32&     sv)
{
    // back insert iterator is faster than random element access for sparse vector
    //
    sparse_vector_u32::back_insert_iterator bi(sv.get_back_inserter());

    vect.resize(test_size);
    bv_null.reset();

    for (unsigned i = 0; i < test_size; ++i)
    {
        unsigned v = unsigned(rand_dis(gen));

        vect[i] = v;
        bv_null[i] = true; // not NULL(assigned) element

        *bi = v; // push back an element to sparse vector

        if (i % 64 == 0)
        {
            bi.add_null(5);  // add 5 unassigned elements using back inserter
            i += 5;  // insert a small NULL plate (unassigned values)
        }
    } // for
}


// plain scan in std::vector<>, matching values are indexed 
// in result bit-vector (subset projection)
// values are added, so multiple calls result in subset addition
static
void vector_search(const std::vector<unsigned>& vect,
                   const bm::bvector<>&         bv_null,
                   unsigned                     value,
                   bm::bvector<>&               bv_res)
{
    bv_res.init(); // always use init() if set_bit_no_check()
    for (size_t i = 0; i < vect.size(); ++i)
    {
        if (vect[i] == value)
            bv_res.set_bit_no_check((bm::id_t)i);
    } // for
    bv_res &= bv_null; // correct results to only include non-NULL values
}


inline
void print_bvector(const bm::bvector<>& bv)
{
    cout << "( count = " << bv.count() << ")" << ": [";
    
    bm::bvector<>::enumerator en = bv.first();
    for (; en.valid(); ++en)
        cout << *en << ", ";
    cout << "]" << endl;
}


int main(void)
{
    try
    {
        // First, lets run, simple (printable) search case
        //
        {
            sparse_vector_u32 sv(bm::use_null);
            
            sv.set(2, 25);
            sv.set(3, 35);
            sv.set(7, 75);
            sv.set(1000, 2000);
            sv.set(256, 2001);
            sv.set(77, 25);

            bm::bvector<> bv_found;  // search results vector
            
            bm::sparse_vector_scanner<sparse_vector_u32> scanner; // scanner class
            scanner.find_eq(sv, 25, bv_found); // seach for all values == 25
            
            print_bvector(bv_found); // print results 

            scanner.invert(sv, bv_found); // invert search results to NOT EQ
            print_bvector(bv_found);  // print all != 25
        }

        
        std::vector<unsigned> vect;
        bm::bvector<> bv_null;
        sparse_vector_u32 sv(bm::use_null);

        {
            bm::chrono_taker<> tt1(cout, "0. test set generate ", 1, &timing_map);
            generate_test_set(vect, bv_null, sv);
        }

        unsigned search_repeats = 500;

        // generate a search vector for benchmarking
        //
        std::vector<unsigned> search_vect;
        {
            bm::bvector<> bv_tmp;
            search_vect.reserve(search_repeats);
            for (unsigned i = 0; i < search_repeats;)
            {
                bm::id_t idx = bm::id_t(rand_dis(gen));
                if (!bv_tmp.test(idx)) // check if number is unique
                {
                    search_vect.push_back(idx);
                    bv_tmp[idx] = 1;
                    ++i;
                }
            }
        }
        
        // run benchmarks
        //
        bm::bvector<> bv_res1;
        bm::bvector<> bv_res2;
        bm::bvector<> bv_res3;

        {
            bm::chrono_taker tt1(cout, "1. std::vector<> scan ", search_repeats, &timing_map);
            
            for (unsigned i = 0; i < search_repeats; ++i)
            {
                unsigned vs = search_vect[i];
                vector_search(vect, bv_null, vs, bv_res1);
            } // for
        }

        {
            bm::chrono_taker tt1(cout, "2. sparse_vector<> scan ", search_repeats, &timing_map);

            bm::sparse_vector_scanner<sparse_vector_u32> scanner;
            scanner.find_eq(sv, search_vect.begin(), search_vect.end(), bv_res2);
        }

        // check jus in case if results look correct
        if (bv_res1.compare(bv_res2) != 0)
        {
            std::cerr << "2. Search result mismatch!" << std::endl;
        }

        {
            bv_res3.init(); // always init before calling "set_bit_no_check()"
            
            bm::chrono_taker tt1(cout, "3. sparse_vector<>::const_iterator search ", search_repeats, &timing_map);

            // prepare a unique search set
            bm::bvector<> bv_search(bm::BM_GAP);
            bm::combine_or(bv_search, search_vect.begin(), search_vect.end());

            sparse_vector_u32::const_iterator it = sv.begin();
            sparse_vector_u32::const_iterator it_end = sv.end();
            for (; it != it_end; ++it)
            {
                unsigned v = *it;
                if (bv_search.test(v))
                {
                    bv_res3.set_bit_no_check(it.pos());
                }
            } // for
        }

        // paranoiya check
        if (bv_res1.compare(bv_res3) != 0)
        {
            std::cerr << "3. Search result mismatch!" << std::endl;
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

