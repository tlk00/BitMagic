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

/** \example strsvsample02a.cpp

   Example shows how to use optimized bm::str_sparse_vector::compare() functions with std::sort
   And then use bm::str_sparse_vector::remap() to reduce memory footprint.
   Succinct methods are very responsive to sort order and other regularities in the data.

  \sa bm::str_sparse_vector
  \sa bm::str_sparse_vector::const_iterator
  \sa bm::str_sparse_vector::back_insert_iterator
  \sa bm::str_sparse_vector::compare
  \sa bm::str_sparse_vector::remap
  \sa bm::str_sparse_vector::optimize

*/

/*! \file strsvsample02a.cpp
    \brief Example: str_sparse_vector<>  sort example
*/

#include <iostream>
#include <string>
#include <vector>
#include <random>
#include <algorithm>

#include "bm.h"
#include "bmstrsparsevec.h"
#include "bmsparsevec_algo.h"

#include "bmtimer.h"
#include "bmdbg.h"

#include "bmundef.h" /* clear the pre-proc defines from BM */

using namespace std;

typedef bm::bvector<> bvector_type;
typedef bm::str_sparse_vector<char, bvector_type, 3> str_sv_type;


/// generate collection of strings from integers and shuffle it
///
static
void generate_string_set(vector<string>& str_vec,
                         const unsigned max_coll = 250000)
{
    str_vec.resize(0);
    string str;
    for (unsigned i = 10; i < max_coll; i += unsigned(rand() % 3))
    {
        str = to_string(i);
        str_vec.emplace_back(str);
    } // for i
    
    // shuffle the data set
    //
    std::random_device rd;
    std::mt19937       g(rd());
    std::shuffle(str_vec.begin(), str_vec.end(), g);
}

int main(void)
{
    try
    {
        str_sv_type str_sv;
        std::vector<string> str_vec;

        generate_string_set(str_vec);
        {
            auto bit = str_sv.get_back_inserter();
            for (auto it = str_vec.begin(); it != str_vec.end(); ++it)
                bit = *it;
            bit.flush(); // important to avoid descructor exceptions
            str_sv.optimize();
        }
        cout << endl;

        // This approach to sort uses an index vector to perform the sort on
        // plain uncompressed vector of indexes offers better sort performance
        // because swap is faster.
        //
        // In this example we are swapping index elements, leaving the original
        // vector in place with an option to copy it later
        //

        std::vector<uint32_t> index(str_sv.size());
        std::generate(index.begin(), index.end(),
                                         [n = 0] () mutable { return n++; });
        {{
        bm::chrono_taker tt(cout, "1.std::sort() of index: ", 0); // timing
#if (1)
        // fastest variant uses local cache to keep one of the comparison vars
        // to reduce access to the succinct vector, right variable gets more hits
        // due to specifics of sort implementation
        // (improves comparison performance, alternative variant provided)
        //

        std::sort(index.begin(), index.end(),
            [&str_sv](const uint32_t l, const uint32_t r)
            {
                static thread_local string last_right_str; // caching variable
                static thread_local uint32_t last_right = uint32_t(-1);
                if (last_right != r)
                {
                    last_right = r;
                    str_sv.get(last_right, last_right_str);
                }
                return str_sv.compare(l, last_right_str.c_str()) < 0;
            } // lambda
        ); // sort
#else

        // possible alternative is to use sort with compare (l, r) it is
        // measurably slower (due to hight latency more access to succint vector)
        // (code is somewhat simpler)

        std::sort(index.begin(), index.end(),
            [&str_sv](const uint32_t l, const uint32_t r)
                    { return str_sv.compare(l, r) < 0; } // lambda
        ); // sort
#endif
        }}

        // copy the sorted vector
        //
        str_sv_type str_sv_sorted;
        {
            char buf[1024]; // value buffer sized for the max.possible len
            auto bit = str_sv_sorted.get_back_inserter();
            for (auto it = index.begin(); it != index.end(); ++it)
            {
                auto i = *it;                     // get index
                str_sv.get(i, buf, sizeof(buf));  // read value by index
                bit = (const char*)buf;           // to back inserter
            }
            bit.flush();
        }

        assert(str_sv_sorted.size()==str_sv.size());


        // validate the results to match STL sort
        //
        std::sort(str_vec.begin(), str_vec.end());
        {
            std::vector<string>::const_iterator sit = str_vec.begin();
            auto it = str_sv_sorted.begin();
            auto it_end = str_sv_sorted.end();
            for (; it != it_end; ++it, ++sit)
            {
                string s = *it;
                if (*sit != s)
                {
                    cerr << "Mismatch at:" << s << "!=" << *sit << endl;
                    return 1;
                }
            } // for
        }
        cout << "Sort validation Ok." << endl << endl;
        cout << "Memory footprint statistics:\n" << endl;
        //
        //
        // compute succinct container stats before remapping
        str_sv_type::statistics st_sorted;
        str_sv_sorted.calc_stat(&st_sorted);
        cout << "SV sorted(before remap) memory_used : "
             << st_sorted.memory_used << endl;

        // memory remapping turns vector into basically read-only
        // (some careful modifications are allowed), but often saves memory
        // due to frequency based succinct codes
        //
        str_sv_sorted.remap();
        str_sv_sorted.optimize();


        str_sv_type::statistics st;
        str_sv.calc_stat(&st);
        str_sv_sorted.calc_stat(&st_sorted);

        cout << "SV unsorted memory_used             : "
             << st.memory_used << endl;

        cout << "SV sorted(after remap) memory_used  : "
             << st_sorted.memory_used << endl;

    }
    catch(std::exception& ex)
    {
        std::cerr << ex.what() << std::endl;
        return 1;
    }
    

    return 0;
}

