/*
Copyright(c) 2002-2022 Anatoliy Kuznetsov(anatoliy_kuznetsov at yahoo.com)

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

/** \example strsvsample08.cpp
    Succinct container binary search

  \sa bm::str_sparse_vector
  \sa bm::sparse_vector_scanner
  \sa bm::sparse_vector_scanner::bind
  \sa bm::sparse_vector_scanner::bfind_eq_str
  \sa bm::str_sparse_vector::freeze
  \sa bm::str_sparse_vector::remap
  \sa bm::str_sparse_vector::optimize
*/

/*! \file strsvsample06.cpp
    \brief Example: Succinct container binary search
*/

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <random>
#include <cassert>

#include "bm.h"
#include "bmstrsparsevec.h"
#include "bmsparsevec_algo.h"
#include "bmtimer.h"
#include "bmundef.h" /* clear the pre-proc defines from BM */

using namespace std;

typedef bm::bvector<> bvector_type;
typedef bm::str_sparse_vector<char, bvector_type, 16> str_sv_type;


static
void GenerateTestStrCollection(std::vector<string>& str_coll, unsigned max_coll)
{
    string prefix = "az";
    string str;
    for (unsigned i = 0; i < max_coll; ++i)
    {
        str = prefix;
        str.append(to_string(i));
        str_coll.emplace_back(str);

        if (i % 1024 == 0) // generate new prefix
        {
            prefix.clear();
            unsigned prefix_len = (unsigned)rand() % 5;
            for (unsigned j = 0; j < prefix_len; ++j)
            {
                char cch = char('a' + (unsigned)rand() % 26);
                prefix.push_back(cch);
            } // for j
        }
    } // for i
}


int main(void)
{
    const unsigned max_coll = 30000000;

    try
    {
        std::vector<string> str_coll;
        str_sv_type str_sv;

        cout << "Prepare test data" << endl;
        cout << " generation ..." << endl;
        GenerateTestStrCollection(str_coll, max_coll);

        cout << " sort..." << endl;
        std::sort(str_coll.begin(), str_coll.end());

        str_sv_type str_sv_srt; // sorted succinct vector
        // fill in using back-insert iterator (fastest mathod to fill in)
        {
            auto bi = str_sv_srt.get_back_inserter();
            for (auto str : str_coll)
               bi = str;
            // important to explicitly flush
            //   ( because destruction flush is not exception safe )
            bi.flush();
        }
        cout << " remap-optimize-freeze..." << endl;
        str_sv_srt.remap();     // apply remapping compression
        str_sv_srt.optimize();  // RLE compression
        str_sv_srt.freeze();    // turn into read-only mode

        // pick the test sets
        std::vector<string> str_coll_test1;
        {
            const unsigned pick_factor = 5;
            for (size_t i = 0; i < size_t(str_coll.size()); i+=pick_factor)
            {
                const string& s = str_coll[i];
                str_coll_test1.push_back(s);
            } // for i
            std::random_device rd;
            std::mt19937 g(rd());
            std::shuffle(str_coll_test1.begin(), str_coll_test1.end(), g);
        }

        cout << "Run test" << endl;


        size_t sum_pos4=0, sum_pos32=0;

        // code below illustrates bm::sparse_vector_scanner<> search sampling
        // parameters: 4 and 32 where 32 offers better performance.
        //
        // the use case here is to create a vector-scanner pair, bind then
        // and repeat multiple searches to vector using its scanner
        // (scanner maintains the sampled search index (for sorted searches).
        //
        // Sampling approach works well, because BM uses hybrid binary search
        // first narrowing down the seach using uncomprssed samples O(Nlog(N))
        // then at the end running a vector search using logical ops with
        // search space prunning.
        //
        // all of the above implements fast search in sorted-compressed array
        // without decompression
        //

        {
            str_sv_type::size_type pos;
            bm::sparse_vector_scanner<str_sv_type, 4> scanner;
            scanner.bind(str_sv_srt, true); // bind sorted vector

            bm::chrono_taker tt(cout, "bm::sparse_vector_scanner<>::bfind_eq_str() [4] ", 1);
            for (unsigned i = 0; i < unsigned(str_coll_test1.size()); ++i)
            {
                const string& s = str_coll_test1[i];
                bool found = scanner.bfind_eq_str(s.c_str(), pos);
                if (!found)
                {
                    cerr << "String bfind_eq_str() failure!" << endl;
                    assert(0); exit(1);
                }
                sum_pos4 += pos;
            } // for
        }

        {
            str_sv_type::size_type pos;
            bm::sparse_vector_scanner<str_sv_type, 32> scanner;
            scanner.bind(str_sv_srt, true); // bind sorted vector

            bm::chrono_taker tt(cout, "bm::sparse_vector_scanner<>::bfind_eq_str() [32] ", 1);
            for (unsigned i = 0; i < unsigned(str_coll_test1.size()); ++i)
            {
                const string& s = str_coll_test1[i];
                bool found = scanner.bfind_eq_str(s.c_str(), pos);
                if (!found)
                {
                    cerr << "String bfind_eq_str() failure!" << endl;
                    assert(0); exit(1);
                }
                sum_pos32 += pos;
            } // for
        }

        assert(sum_pos4 == sum_pos32);

    }
    catch(std::exception& ex)
    {
        std::cerr << ex.what() << std::endl;
        return 1;
    }


    return 0;
}

