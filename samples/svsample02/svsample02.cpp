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

/** \example svsample02.cpp
  Example of how to serialize bm::sparse_vector<> template class
 
  \sa bm::sparse_vector<>
  \sa bm::sparse_vector<>::push_back
  \sa bm::sparse_vector<>::equal
  \sa bm::sparse_vector_serialize
  \sa bm::sparse_vector_deserialize
  
 */

#include <iostream>
#include <vector>

#include "bmsparsevec.h"
#include "bmsparsevec_serial.h"

using namespace std;

typedef bm::sparse_vector<unsigned, bm::bvector<> > svector;

int main(void)
{
    try
    {
        svector sv1;
        svector sv2;

        // temp buffer to avoid unnecessary re-allocations
        BM_DECLARE_TEMP_BLOCK(tb)


        for (unsigned i = 0; i < 128000; ++i)
        {
            sv1.push_back(8);
        }
        
        // optimize memory allocation of sparse vector
        sv1.optimize(tb);
        
        bm::sparse_vector_serial_layout<svector> sv_lay;
        int res = bm::sparse_vector_serialize(sv1, sv_lay, tb);
        if (res != 0)
        {
            cerr << "Sparse vector serialization error!" << endl;
            exit(1);
        }
        
        // copy serialization buffer to some other location
        // to simulate data-base storage or network transaction
        //
        const unsigned char* buf = sv_lay.buf();
        size_t buf_size = sv_lay.size();
        
        vector<unsigned char> tmp_buf(buf_size);
        ::memcpy(&tmp_buf[0], buf, buf_size);
        
        res = bm::sparse_vector_deserialize(sv2, &tmp_buf[0], tb);
        if (res != 0)
        {
            cerr << "De-Serialization error in TestEqualSparseVectors()" << endl;
            exit(1);
        }
        
        
        if (!sv1.equal(sv2) )
        {
            cerr << "Error! Please report a bug to BitMagic project support." << endl;
        }
        
        cout << sv2.size() << endl;
    }
    catch(std::exception& ex)
    {
        std::cerr << ex.what() << std::endl;
        return 1;
    }
    
    

    return 0;
}

