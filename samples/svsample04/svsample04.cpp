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

/** \example svsample04.cpp
  Example how to use NULL-enabled sparse vector
 
  \sa bm::sparse_vector
  \sa bm::sparse_vector::clear
  \sa bm::sparse_vector::set_null
  \sa bm::sparse_vector::push_back
  \sa bm::sparse_vector::decode
  \sa bm::sparse_vector::resize
  \sa bm::sparse_vector::join
  \sa bm::bvector::enumerator

*/

/*! \file svsample04.cpp
    \brief Example: sparse_vector<> with NULL (unset) values
*/

#include <iostream>
#include <vector>

#include "bm.h"
#include "bmsparsevec.h"
#include "bmundef.h" /* clear the pre-proc defines from BM */

using namespace std;

static
void print_svector(const bm::sparse_vector<unsigned, bm::bvector<> >& sv)
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

int main(void)
{
    try
    {
        // sparse vector can support unassigned value semantics
        // (also known as NULL value in databases)
        //
        // to enable this function sparse_vector<> needs to be constructed using
        // special flag ( bm::use_null )
        //
        bm::sparse_vector<unsigned, bm::bvector<> > sv1(bm::use_null);
        
        sv1.resize(10);
        sv1.set(2, 25);
        sv1.set(3, 35);
        sv1.set(7, 75);
        
        print_svector(sv1); // 10: [ NULL, NULL, 25, 35, NULL, NULL, NULL, 75, NULL, NULL ]
        
        sv1.set_null(7);

        print_svector(sv1); // 10: [ NULL, NULL, 25, 35, NULL, NULL, NULL, NULL, NULL, NULL ]
        
        unsigned arr[3] = {1,2,3};
        sv1.import(arr, 3, 9); // import from a C-style array (by index 9)
        print_svector(sv1); // 12: [ NULL, NULL, 25, 35, NULL, NULL, NULL, NULL, NULL, 1, 2, 3 ]
        
        sv1.clear(2, true); // clear element and set to NULL
        print_svector(sv1); // 12: [ NULL, NULL, NULL, 35, NULL, NULL, NULL, NULL, NULL, 1, 2, 3 ]

        sv1.clear(2, false); // clear element (not NULL anymore)
        print_svector(sv1); // 12: [ NULL, NULL, 0, 35, NULL, NULL, NULL, NULL, NULL, 1, 2, 3 ]

        sv1.clear(3, false); // clear element (it stays NOT NULL)
        print_svector(sv1); // 12: [ NULL, NULL, NULL, 0, NULL, NULL, NULL, NULL, NULL, 1, 2, 3 ]

        sv1.clear();
        print_svector(sv1); // 0: [ EMPTY ]
        
        sv1.resize(3);
        print_svector(sv1); // 3: [ NULL, NULL, NULL ]
        
        sv1.push_back(10);
        sv1.push_back(20);
        sv1.push_back(30);

        print_svector(sv1); // 5: [ NULL, NULL, NULL, 10, 20, 30 ]

        // bulk clear the vector (can be a lot faster than random element access)
        {
            bm::bvector<> bv_clear { 3, 5 }; // clear sv1[3]=0 and sv1[5]=0
            sv1.clear(bv_clear); // clear multiple elemets at once
            print_svector(sv1); // 5: [ NULL, NULL, NULL, 0, 20, 0 ]
            // please note, that clear does not set elements to NULL (but to 0)
        }

        sv1.set(1, 10);
        sv1.set(2, 20);
        sv1.push_back(40);

        print_svector(sv1); // 5: [ NULL, 10, 20, 0, 20, 0, 40 ]
        // bulk set_null the vector 
        {
            bm::bvector<> bv_set_null { 0, 1, 3, 5 }; // NULLs sv1[3] and sv1[5]
            sv1.set_null(bv_set_null); // clear multiple elemets at once
            print_svector(sv1); // 5: [ NULL, NULL, 20, NULL, 20, NULL, 40]
            // please note, that clear both sets to zero and sets it to NULL
        }



        
        // construct not-null-able sparse vector
        //
        bm::sparse_vector<unsigned, bm::bvector<> > sv2;
        sv2.push_back(100);
        sv2.push_back(200);
        
        sv1.join(sv2);  // merge two vectors together
        
        print_svector(sv1); // 5: [ 100, 200, NULL, 10, 20 ]
        
        // construct another NULL-able vector for another join
        bm::sparse_vector<unsigned, bm::bvector<> > sv3(bm::use_null);
        sv3.resize(9);
        sv3.push_back(300);
        
        // this join will fully respect the NULL flags assigned in both sparse vectors
        sv1.join(sv3);
        
        print_svector(sv1); // 10: [ 100, 200, NULL, 10, 20, NULL, NULL, NULL, NULL, 300 ]

        // traverse and print the non-NULL values
        // using bit-vector of assigned values
        //
        //          [0] = 100, [1] = 200, [3] = 10, [4] = 20, [9] = 300
        //
        {
            const bm::bvector<>* bv_non_null = sv1.get_null_bvector();
            if (bv_non_null)
            {
                bm::bvector<>::enumerator en = bv_non_null->first();

                for (; en.valid(); ++en)
                {
                    unsigned idx = *en;
                    unsigned v = sv1[idx];
                    std::cout << "[" << idx << "] = " << v << ", ";
                }
                std::cout << endl;
            }
        }


        // decode sparse_vector with NULL values
        // decode ignores NULLs (in our case replaces values with '0s'
        //
        {
            std::vector<unsigned> v1(sv1.size());
            sv1.decode(&v1[0], 0, sv1.size()); // extract elements starting from 0
            
            // 100,200,0,10,20,0,0,0,0,300,
            for (unsigned i = 0; i < v1.size(); ++i)
            {
                cout << v1[i] << ",";
            }
            cout << endl;
        }
    }
    catch(std::exception& ex)
    {
        std::cerr << ex.what() << std::endl;
        return 1;
    }

    return 0;
}

