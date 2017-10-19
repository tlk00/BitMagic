/*
Copyright(c) 2002-2017 Anatoliy Kuznetsov(anatoliy_kuznetsov at yahoo.com)

Permission is hereby granted, free of charge, to any person
obtaining a copy of this software and associated documentation
files (the "Software"), to deal in the Software without restriction,
including without limitation the rights to use, copy, modify, merge,
publish, distribute, sublicense, and/or sell copies of the Software,
and to permit persons to whom the Software is furnished to do so,
subject to the following conditions:

The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
OTHER DEALINGS IN THE SOFTWARE.

You have to explicitly mention BitMagic project in any derivative product,
its WEB Site, published materials, articles or any other work derived from this
project or based on our code or know-how.

For more information please visit:  http://bitmagic.io

*/

/** \example svsample03.cpp
  Example of how to merge sparse vectors and extract values fast
 
  \sa bm::sparse_vector<>::set
  \sa bm::sparse_vector<>::import
  \sa bm::sparse_vector<>::export
  \sa bm::sparse_vector<>::join
 */

#include <iostream>
#include <vector>
#include "bmsparsevec.h"

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
    sv1.extract(&v1[0], 16, 65530); // extract 16 elements starting from 65530
    for (i = 0; i < 16; ++i)
    {
        cout << v1[i] << ",";
    }
    cout << endl;
    
    

    return 0;
}

