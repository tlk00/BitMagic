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

/** \example rscsample01.cpp
  Example of how to use bm::rsc_sparse_vector<> template class
 
  \sa bm::sparse_vector
  \sa bm::rsc_sparse_vector
  \sa bm::sparse_vector::push_back
  \sa bm::sparse_vector_serialize
  \sa bm::sparse_vector_deserialize
  
*/

/*! \file rscsample01.cpp
    \brief Example: rsc_sparse_vector<> usage
 
    rsc_sparse_vector<> is a sparse vector which uses bit-transposition and
    rank-select succinct method of compression of NULL (unassigned) values.
    Unassigned values are dropped (as transposed columns) from the bit-matrix.
 
    rsc_sparse_vector<> is basically a read-only structure, which can be used
    for compact data storage and search. Random access to elements is possible
    with a penalty of bit-vector Rank or Select operations.
*/

#include <iostream>
#include <vector>

#include "bm.h"
#include "bmsparsevec.h"
#include "bmsparsevec_compr.h"
#include "bmsparsevec_serial.h"
#include "bmundef.h" /* clear the pre-proc defines from BM */

using namespace std;

typedef bm::sparse_vector<unsigned, bm::bvector<> >         sparse_vector_u32;
typedef bm::rsc_sparse_vector<unsigned, sparse_vector_u32 > rsc_sparse_vector_u32;

/// Prints the vector using is_null() and get()
/// Please note, that random access is not the fastest, const_iterator 
/// is preferred in many cases (better performance)
/// 
template<class SV>
void print_svector(const SV& sv, bool show_nulls = false)
{
    std::cout << sv.size() << ": ";
    for (unsigned i = 0; i < sv.size(); ++i)
    {
        if (show_nulls)
        {
            if (sv.is_null(i))
                std::cout << "NULL, ";
            else
            {
                typename SV::value_type v = sv.get(i);
                std::cout << v << ", ";
            }
        }
        else
        {
            typename SV::value_type v = sv.get(i);
            std::cout << v << ", ";
        }
    }
    std::cout << std::endl;
}

/// Prints succinct vector using try_get() random access method
/// 
///
template<class SV>
void print_svector2(const SV& sv, bool show_nulls = false)
{
    using value_type = typename SV::value_type;
    value_type v;
    std::cout << sv.size() << ": ";
    for (unsigned i = 0; i < sv.size(); ++i)
    {
        if (show_nulls)
        {
            if (!sv.try_get(i, v))
                std::cout << "NULL, ";
            else
                std::cout << v << ", ";
        }
        else
        {
            v = sv.get(i);
            std::cout << v << ", ";
        }
    }
    std::cout << std::endl;
}

int main(void)
{
    // temp buffer to avoid unnecessary re-allocations
    BM_DECLARE_TEMP_BLOCK(tb)

    try
    {
        sparse_vector_u32 sv1(bm::use_null); // use null is needed to build rsc vector
        rsc_sparse_vector_u32 csv2;
        
        // fill in sparse vector leaving some unassigned gaps
        for (unsigned i = 0; i < 15; i+=3)
        {
            sv1[i] = i;
        }
        
        print_svector(sv1); // print sparse vector disregard NULL values (printed as 0)
        print_svector2(sv1, true); // print sparse vector show NULLs (unassigned)

        csv2.load_from(sv1); // load rank-select-compacted (rsc) sparse vector
        
        // print results - it should look the same
        print_svector(csv2);
        print_svector2(csv2, true);

        // serialize rsc vector
        
        // optimize memory allocation of sparse vector
        csv2.optimize(tb);
        
        bm::sparse_vector_serial_layout<rsc_sparse_vector_u32> sv_lay;
        bm::sparse_vector_serialize(csv2, sv_lay, tb);

        // memory copy just simulates network or DB transaction
        const unsigned char* buf = sv_lay.buf();
        size_t buf_size = sv_lay.size();

        vector<unsigned char> tmp_buf(buf_size);
        ::memcpy(&tmp_buf[0], buf, buf_size);


        rsc_sparse_vector_u32 csv3;  // target vector

        int res = bm::sparse_vector_deserialize(csv3, &tmp_buf[0], tb);
        if (res != 0)
        {
            std::cerr << "De-Serialization error!" << std::endl;
            return 1;
        }
        if (!csv3.equal(csv2) )
        {
            cerr << "Error! Please report a bug to BitMagic project support." << endl;
            return 1;
        }

        // unload rsc to sv, this makes a round-trip to an editable form
        sparse_vector_u32 sv3(bm::use_null);
        csv3.load_to(sv3);
        
        if (!sv3.equal(sv1) )
        {
            std::cerr << "Error! Please report a bug to BitMagic project support." << std::endl;
            return 1;
        }
        
        print_svector(sv1, true); // print sparse vector again

    }
    catch(std::exception& ex)
    {
        std::cerr << ex.what() << std::endl;
        return 1;
    }
    
    

    return 0;
}

