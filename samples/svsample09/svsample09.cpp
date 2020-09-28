/*
Copyright(c) 2002-2019 Anatoliy Kuznetsov(anatoliy_kuznetsov at yahoo.com)

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

/** \example svsample09.cpp
    Mismatch search. bm::sparse_vector_find_first_mismatch
 
  \sa bm::sparse_vector
  \sa bm::sparse_vector_find_first_mismatch

  \sa bm::bvector<>::find_first_mismatch()

*/

/*! \file svsample09.cpp
    \brief Example: Use of sparse vector mismatch search
*/

#include <assert.h>
#include <stdlib.h>

#include <iostream>
#include <vector>
#include <utility>

#include "bm.h"
#include "bmsparsevec.h"
#include "bmsparsevec_algo.h"
#include "bmundef.h" /* clear the pre-proc defines from BM */

using namespace std;


typedef bm::sparse_vector<unsigned, bm::bvector<> > svector_u32;


int main(void)
{
    try
    {
        svector_u32::size_type pos;

        cout << "sparse vectors without NULL:" << endl;
        {
            svector_u32 sv1, sv2;

            sv1.push_back(1);
            sv1.push_back(2);
            sv1.push_back(2);
            sv1.push_back(3);

            sv2 = sv1;

            bool found = bm::sparse_vector_find_first_mismatch(sv1, sv2, pos);
            if (found)
                std::cout << "Mismatch found." << endl; // this would be an error
            else
                std::cout << "Mismatch not found" << endl; // identical vectors

            // modify sv2 to introduce a mismatch

            sv2[2] = 0;
            found = bm::sparse_vector_find_first_mismatch(sv1, sv2, pos);
            if (found)
            {
                std::cout << "Mismatch found at: " << pos << endl; // at 2
            }
        }

        // bm::sparse_vector_find_first_mismatch works with NULL-able vectors
        // first NULL (unsassigned) between vectors is a mimatch
        //
        cout << endl << "sparse vectors with NULL:" << endl;
        {
            svector_u32 sv1(bm::use_null);
            svector_u32 sv2(bm::use_null);

            sv1[1] = 1;
            sv1[2] = 2;
            sv1.set_null(3); // set element 3 to NULL
            sv1[4] = 0;

            sv2 = sv1;

            bool found = bm::sparse_vector_find_first_mismatch(sv1, sv2, pos);
            if (found)
                std::cout << "Mismatch found." << endl; // this would be an error
            else
                std::cout << "Mismatch not found" << endl; // identical vectors

            sv2[4] = 10;
            found = bm::sparse_vector_find_first_mismatch(sv1, sv2, pos);
            if (found)
            {
                std::cout << "Mismatch found at: " << pos << endl; // at 4
            }


            // set element 3 to 0, now we have a situation when element 3 in both
            // vectors has value 0, except it is NULL value in vector 1
            // mismatch search should detect it
            //
            sv2[3] = 0;
            found = bm::sparse_vector_find_first_mismatch(sv1, sv2, pos);
            if (found)
            {
                std::cout << "Mismatch found at: " << pos << endl; // at 3
            }
        }

    }
    catch(std::exception& ex)
    {
        std::cerr << ex.what() << std::endl;
        return 1;
    }
    
    

    return 0;
}

