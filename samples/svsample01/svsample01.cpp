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

/** \example svsample01.cpp
  Example of how to use bm::sparse_vector<> template class to set values
 
  \sa bm::sparse_vector<>::import
  \sa bm::sparse_vector<>::at
  \sa bm::sparse_vector<>::optimize
  \sa bm::sparse_vector<>::size
*/

/*! \file svsample01.cpp
    \brief Example: sparse_vector<> container set values
*/


#include <iostream>

#include "bm.h"
#include "bmsparsevec.h"

using namespace std;

int main(void)
{
    try
    {
        bm::sparse_vector<unsigned, bm::bvector<> > sv1;
        
        unsigned arr[3] = {1,2,3};
        sv1.import(arr, 3); // import from a C-style array (fastest way to populate)

        // optimize memory allocation of sparse vector
        {
            BM_DECLARE_TEMP_BLOCK(tb)
            sv1.optimize(tb);
        }
        
        cout << "sv1.size() = " << sv1.size() << endl;
        cout << "sv[]:";
        
        // print the vector elements using direct access operator
        for (unsigned i = 0; i < sv1.size(); ++i)
        {
            cout << sv1.at(i) << ",";
        }
        cout << endl;

        // add more at the end
        unsigned arr2[5] = {10, 20, 30, 40, 50};
        sv1.import(arr2, 5, sv1.size());
        
        cout << "sv1.size() = " << sv1.size() << endl;
        cout << "sv[]:";
        for (unsigned i = 0; i < sv1.size(); ++i)
        {
            cout << sv1.at(i) << ",";
        }
        cout << endl;

        
    }
    catch(std::exception& ex)
    {
        std::cerr << ex.what() << std::endl;
        return 1;
    }
    

    return 0;
}

