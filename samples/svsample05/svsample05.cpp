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

/** \example svsample05.cpp
  Example how to use for bvector re-mapping / transformation
  based on sparse_vector as a translation table.

  Example discusses how use of NULL (unassigned) sparse_vector<> semantics affects
  behavior of image transformation algorithm.
 
  \sa bm::sparse_vector<>
  \sa bm::set2set_11_transform
*/

/*! \file svsample05.cpp
    \brief Example: sparse_vector<> used for set 2 set remapping (theory of groups Image)
*/


#include <iostream>
#include <vector>

#include "bm.h"
#include "bmsparsevec.h"
#include "bmsparsevec_algo.h"
#include "bmundef.h" /* clear the pre-proc defines from BM */

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
        // initialize sparse_vector as a binary relation (transform function)
        // to translate from one set to another (please note it uses NULL values)
        //
        sparse_vector_u32 sv_transform(bm::use_null);
        
        // transformation function defined as a table:
        // 2 -> 25
        // 3 -> 35
        // 7 -> 75
        // 1000 -> 2000
        // 256 -> 2001
        
        sv_transform.set(2, 25);
        sv_transform.set(3, 35);
        sv_transform.set(7, 75);
        sv_transform.set(1000, 2000);
        sv_transform.set(256, 2001);
        
        print_svector(sv_transform);
        
        // initialize input bit-vector
        //
        bm::bvector<> bv_in { 1, 2, 3, 1000, 256 };
        
        bm::bvector<> bv_out;
       
        bm::set2set_11_transform<sparse_vector_u32> func;
        func.run(bv_in, sv_transform, bv_out);
        
        std::cout << "re-mapped vector:" << endl;
        print_bvector(bv_out); // ( count = 4): [25, 35, 2000, 2001, ]

        // Please note, that remapping request for "1" cannot be satisfied
        // as sv_transform vector did not assign this element so it has
        // special unassigned, NULL value

        bm::id_t from = 1;
        bm::id_t to;
        bool found = func.remap(1, sv_transform, to);
        if (!found)
            std::cout << from << " not mapped" << std::endl;  // THIS works!
        else
            std::cout << from << " mapped to " << to << std::endl;


        // now lets try another experiment here, construct image transform 
        // function as default NOT using NULL values

        sparse_vector_u32 sv_transform2;

        sv_transform2.set(2, 25);
        sv_transform2.set(3, 35);
        sv_transform2.set(7, 75);
        sv_transform2.set(1000, 2000);
        sv_transform2.set(256, 2001);

        found = func.remap(1, sv_transform2, to);
        if (!found)
            std::cout << from << " not mapped" << std::endl;
        else
            std::cout << from << " mapped to " << to << std::endl; // Now THIS works!

        // Please note that remapping for unassigned value - works and maps to 0
        // because transformation sparse vector does not keep information about
        // unassigned values, all unassigned are now implicitly assumed to be "0"

        func.run(bv_in, sv_transform2, bv_out);

        std::cout << "re-mapped vector2:" << endl;
        print_bvector(bv_out); // ( count = 5): [0, 25, 35, 2000, 2001, ]

    }
    catch(std::exception& ex)
    {
        std::cerr << ex.what() << std::endl;
        return 1;
    }

    return 0;
}

