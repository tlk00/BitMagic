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

/** \example strsvsample09.cpp
  Example of how to use bm::str_sparse_vector<> - succinct container for sorting in compressive memory

 
  \sa bm::str_sparse_vector
  \sa bm::str_sparse_vector::swap
  \sa bm::str_sparse_vector::compare
*/

/*! \file strsvsample09.cpp
    \brief Example: str_sparse_vector<> sorting example
*/

#include <iostream>
#include <string>
#include <vector>
#include <random>
#include <algorithm>
//#define BMAVX2OPT
//#define BMSSE42OPT
#include "bm.h"
#include "bmstrsparsevec.h"
#include "bmsparsevec_algo.h"

#include "bmdbg.h"
#include "bmtimer.h"

#include "bmundef.h" /* clear the pre-proc defines from BM */

using namespace std;

typedef bm::bvector<> bvector_type;

// define the sparse vector type for 'char' type using bvector as
// a container of bits for bit-transposed planes
typedef bm::str_sparse_vector<char, bvector_type, 16> str_sv_type;


/// generate collection of strings from integers with common prefixes
/// ... and shuffle it
static
void generate_string_set(vector<string>& str_vec,
                         const unsigned max_coll = 850000,
                         unsigned repeat = 220)
{
    str_vec.resize(0);
    string str;
    for (unsigned i = 10; i < max_coll; i += unsigned(rand() % 3))
    {
        switch (rand()%8)
        {
        case 0: str = "xnssv"; break;
        default: str = "xrs";  break;
        }
        str.append(to_string(i));

        for (unsigned k = 0; k < repeat; ++k, ++i) // add more of the same string
            str_vec.emplace_back(str);

    } // for i

    std::random_device rd;
    std::mt19937       g(rd());
    std::shuffle(str_vec.begin() + str_vec.size()/2, str_vec.end(), g);

}


/// quick-sort
///
void quicksort(str_sv_type& strsv, int first, int last)
{
    using stype = str_sv_type::size_type;
    int i, j, pivot;
    if (first < last)
    {
        pivot = i= first;
        j = last;
        while (i <j)
        {
            while((i < last) && (strsv.compare(stype(i), stype(pivot)) <= 0))
                i++;
            while(strsv.compare(stype(j), stype(pivot)) > 0) // number[j]>number[pivot])
                j--;
            if (i < j)
                strsv.swap(stype(i), stype(j));
        } // while
        strsv.swap(stype(pivot), stype(j));

        quicksort(strsv, first, j-1);
        quicksort(strsv, j+1, last);
    }
}


/// Faster variant of quicksort, uses different variant of pivot compare, with decompressed argument
///
/// optimizations:
/// 1. use of bm::str_sparse_vector<>::compare() function friendly for re-use of pivot element
/// 2.  tail call recursion eleimination
///
void quicksort2(str_sv_type& strsv, int first, int last)
{
    using stype = str_sv_type::size_type;
    int i, j, pivot;

    // fixed size for simplicity (in prod code needs dynamic buffer handling)
    static str_sv_type::value_type pivot_buf[128];
    while (first < last)
    {
        pivot = i = first;
        j = last;

        // save the pivor to re-use it in strsv.compare(..)
        strsv.get(stype(pivot), pivot_buf, sizeof(pivot_buf));

        while (i < j)
        {
            while((i < last) && (strsv.compare(stype(i), pivot_buf) <= 0))
                ++i;
            while(strsv.compare(stype(j), pivot_buf) > 0)
                --j;
            if (i < j)
                strsv.swap(stype(i), stype(j));
        } // while
        strsv.swap(stype(pivot), stype(j));

        quicksort2(strsv, first, j-1);
        first = j+1; // tail recursion
    } // while
}


/// insertion sort for performance comnparison
///
static
void insertion_sort(str_sv_type& str_sv, const vector<string>& str_vec)
{
    // scanner object is re-used throught the processing
    //
    bm::sparse_vector_scanner<str_sv_type> scanner;

    for (const string& s : str_vec)
    {
        const char* cs = s.c_str();
        str_sv_type::size_type pos;
        bool found = scanner.lower_bound_str(str_sv, cs, pos);
        (void)found; // just to silence the unused variable warning

        str_sv.insert(pos, cs);

    } // for s
}



bm::chrono_taker<>::duration_map_type  timing_map; // timing stats

int main(void)
{
    try
    {
        str_sv_type str_sv, str_sv2, str_sv3;

        vector<string> str_vec;
        generate_string_set(str_vec);

        // load compact vector
        cout << "Loading " << str_vec.size() << " elements..." << endl;
        {
            auto bi = str_sv.get_back_inserter();
            for (const string& term : str_vec)
                bi = term;
            bi.flush();
        }
        // remap succinct vector into optimal codes
        // (after final load of content)
        //
        str_sv.remap();
        {
            BM_DECLARE_TEMP_BLOCK(tb)
            str_sv.optimize(tb);
        }

       //  print_svector_stat(cout, str_sv);

        str_sv2 = str_sv;

        // calculate and print memory usage statistics
        {
            str_sv_type::statistics st;
            str_sv.calc_stat(&st);
            size_t std_vect_mem = sizeof(str_vec);
            for (const string& term : str_vec)
                std_vect_mem += term.size() + sizeof(term);

            cout << "std::vector<string> mem    = " << std_vect_mem << endl;
            cout << "Succinct vector vector mem = " << st.memory_used << endl;
        }


        cout << "Quick Sort... 1" << endl;
        {
        bm::chrono_taker tt1(cout, "1. quick sort (succint)", 1, &timing_map);
        quicksort(str_sv, 0, (int)str_sv.size()-1);
        }

        cout << "Quick Sort... 2" << endl;
        {
        bm::chrono_taker tt1(cout, "2. quick sort 2 (succint)", 1, &timing_map);
        quicksort2(str_sv2, 0, (int)str_sv2.size()-1);
        }
/*
        str_sv2.optimize();
        print_svector_stat(cout, str_sv);

        cout << "Quick Sort... (on sorted)" << endl;
        {
            bm::chrono_taker tt1(cout, "3. quick sort 2 (succint) (pr-sorted)", 1, &timing_map);
            quicksort2(str_sv2, 0, (int)str_sv2.size() - 1);
        }
*/

        cout << "Insertion Sort... " << endl;
        {
        bm::chrono_taker tt1(cout, "3. insertion sort (succint)", 1, &timing_map);
        insertion_sort(str_sv3, str_vec);
        }

        cout << "std::sort..." << endl;
        {
        bm::chrono_taker tt1(cout, "4. std::sort()", 1, &timing_map);
        std::sort(str_vec.begin(), str_vec.end());
        }

        // validation of different sort methods
        {
            bool eq = str_sv.equal(str_sv2);
            if (!eq)
            {
                cerr << "post-sort vector mismatch! (qsort 1-2)" << endl;
                assert(0); exit(1);
            }

            // validate the results to match STL sort
            {
                string s, s3;
                vector<string>::const_iterator sit = str_vec.begin();
                str_sv_type::const_iterator it = str_sv.begin();
                str_sv_type::const_iterator it_end = str_sv.end();
                str_sv_type::const_iterator it3 = str_sv3.begin();
                for (; it != it_end; ++it, ++sit, ++it3)
                {
                    s = *it;
                    if (*sit != s)
                    {
                        cerr << "Mismatch at:" << s << "!=" << *sit << endl;
                        assert(0);
                        return 1;
                    }
                    s3 = *it3;
                    if (s != s3)
                    {
                        cerr << "Mismatch at:" << s << "!=" << s3 << endl;
                        return 1;
                    }

                } // for
            }
            cout << "Sort validation Ok." << endl;
        }

        bm::chrono_taker<>::print_duration_map(cout, timing_map, bm::chrono_taker<>::ct_time);

    }
    catch(std::exception& ex)
    {
        std::cerr << ex.what() << std::endl;
        return 1;
    }
    

    return 0;
}

