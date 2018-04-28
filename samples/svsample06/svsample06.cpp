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

/*! \file svsample05.cpp
    \brief Example: sparse_vector<> used for set 2 set remapping (theory of groups Image)
*/


#include <iostream>
#include <vector>
#include "bmsparsevec.h"
#include "bmsparsevec_algo.h"

using namespace std;

typedef bm::sparse_vector<bm::id_t, bm::bvector<> > sparse_vector_u32;

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
        
        print_bvector(bv_found); // ( count = 4): [25, 35, 2000, 2001, ]
    }
    catch(std::exception& ex)
    {
        std::cerr << ex.what() << std::endl;
        return 1;
    }

    return 0;
}

