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

/** \example svsample03.cpp
  Example of how to merge sparse vectors and extract values fast
 
  \sa bm::sparse_vector<>::set
  \sa bm::sparse_vector<>::import
  \sa bm::sparse_vector<>::decode
  \sa bm::sparse_vector<>::join
*/

/*! \file svsample03.cpp
    \brief Example: sparse_vector<> merge and fast extraction of content
*/


#include <iostream>
#include <vector>
#include "bm.h"
#include "bmsparsevec.h"
#include "bmundef.h" /* clear the pre-proc defines from BM */

using namespace std;

int main(void)
{
    bm::sparse_vector<unsigned, bm::bvector<> > sv1;
    bm::sparse_vector<unsigned, bm::bvector<> > sv2;
    unsigned i;
    unsigned arr[3] = {1,2,3};
    
    sv1.import(arr, 3); // import from a C-style array (fastest way to populate)

    // optimize memory allocation of sparse vector
    {
        BM_DECLARE_TEMP_BLOCK(tb)
        sv1.optimize(tb);
    }
    cout << "sv1.size() = " << sv1.size() << ": ";
    for (i = 0; i < sv1.size(); ++i)
    {
        cout << sv1.at(i) << ",";
    }
    cout << endl;

    // populate second sparse vector using dunamic resize assignment
    //
    for (i = 65536; i < 65536+10; ++i)
    {
        sv2.set(i, 256+i);
    }
    cout << "sv2.size() = " << sv2.size() << endl;

    cout << "Perform sparse_vector<>::join()" << endl;
    // join sv2 into sv1
    // operation implemented using bitset OR on all bit-planes.
    sv1.join(sv2);
    sv1.optimize(); // please note, we use default optimize here (it is a bit slower)
    
    cout << "Now sv1.size() = " << sv1.size() << endl;
    
    
    std::vector<unsigned> v1(16);
    sv1.decode(&v1[0], 65530, 16); // extract 16 elements starting from 65530
    for (i = 0; i < 16; ++i)
    {
        cout << v1[i] << ",";
    }
    cout << endl;
    
    

    return 0;
}

