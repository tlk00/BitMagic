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

/** \example svsample02.cpp
  Example of how to serialize bm::sparse_vector<> template class
 
  \sa bm::sparse_vector<>
  \sa bm::sparse_vector<>::push_back
  \sa bm::sparse_vector<>::equal
  \sa bm::sparse_vector_serialize
  \sa bm::sparse_vector_deserialize
  
*/

/*! \file svsample02.cpp
    \brief Example: sparse_vector<> serialization
*/

#include <iostream>
#include <vector>

#include "bm.h"
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
        bm::sparse_vector_serialize(sv1, sv_lay, tb);
        
        // copy serialization buffer to some other location
        // to simulate data-base storage or network transaction
        //
        const unsigned char* buf = sv_lay.buf();
        size_t buf_size = sv_lay.size();
        
        vector<unsigned char> tmp_buf(buf_size);
        ::memcpy(&tmp_buf[0], buf, buf_size);
        
        int res = bm::sparse_vector_deserialize(sv2, &tmp_buf[0], tb);
        if (res != 0)
        {
            cerr << "De-Serialization error!" << endl;
            return 1;
        }
        
        
        if (!sv1.equal(sv2) )
        {
            cerr << "Error! Please report a bug to BitMagic project support." << endl;
            return 1;
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

