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

/** \example strsvsample02.cpp
  Example of how to use bm::str_sparse_vector<> - succinct container for
  bit-transposed string collections
 
  \sa bm::str_sparse_vector
  \sa bm::sparse_vector_scanner
*/

/*! \file strsvsample02.cpp
    \brief Example: str_sparse_vector<> insertion sort example
*/

#include <iostream>
#include <string>
#include <vector>
#include <random>
#include <algorithm>

#include "bm.h"
#include "bmstrsparsevec.h"
#include "bmsparsevec_algo.h"
#include "bmundef.h" /* clear the pre-proc defines from BM */

using namespace std;

typedef bm::bvector<> bvector_type;

// define the sparse vector type for 'char' type using bvector as
// a container of bits for bit-transposed planes
// 32 - is maximum string length for this container.
//      Memory allocation is dynamic using sparse techniques, so this number
//      just defines the max capacity.
//
typedef bm::str_sparse_vector<char, bvector_type, 32> str_sv_type;


// generate collection of strings from integers and shuffle it
//
static
void generate_string_set(vector<string>& str_vec)
{
    const unsigned max_coll = 50000;
   
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

// insertion sort takes data from unsorted vector places it into sparse vector
// maintaining correct sorted order (for fast search)
//
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


int main(void)
{
    try
    {
        str_sv_type str_sv;

        vector<string> str_vec;
        generate_string_set(str_vec);
        
        insertion_sort(str_sv, str_vec);
        
        {
            BM_DECLARE_TEMP_BLOCK(tb)
            str_sv.optimize(tb);
        }
        
        // validate the results to match STL sort
        std::sort(str_vec.begin(), str_vec.end());
        {
            vector<string>::const_iterator sit = str_vec.begin();
            str_sv_type::const_iterator it = str_sv.begin();
            str_sv_type::const_iterator it_end = str_sv.end();
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
        cout << "Sort validation Ok." << endl;
    }
    catch(std::exception& ex)
    {
        std::cerr << ex.what() << std::endl;
        return 1;
    }
    

    return 0;
}

