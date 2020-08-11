/*
Copyright(c) 2002-2020 Anatoliy Kuznetsov(anatoliy_kuznetsov at yahoo.com)

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

/** \example rscsample05.cpp
  Example of how to use collaborative compression (XOR compression)
  for a group of vectors (data-frame)

  \sa bm::rsc_sparse_vector
*/

/*! \file rscsample05.cpp
    \brief Example: collaborative compression (XOR compression)

    @sa rscsample02.cpp
*/

#include <iostream>
#include <vector>
#include <cassert>

#include "bm.h"
#include "bmsparsevec.h"
#include "bmsparsevec_compr.h"
#include "bmsparsevec_serial.h"

using namespace std;


typedef bm::bvector<>                                     bvector_type;
typedef bm::sparse_vector<unsigned, bvector_type>         sparse_vector_u32;
typedef bm::rsc_sparse_vector<unsigned, sparse_vector_u32 > rsc_sparse_vector_u32;

typedef bm::sparse_vector_serializer<rsc_sparse_vector_u32> csv_serializer_type;
typedef bm::sparse_vector_deserializer<rsc_sparse_vector_u32> csv_deserializer_type;

/// sample data-frame structure with multiple bm vectors
///
struct sample_data_frame
{
    rsc_sparse_vector_u32 csv1;
    rsc_sparse_vector_u32 csv2;
    rsc_sparse_vector_u32 csv3;
};

/// generate some data just to illustrate the case
///
static
void fill_test_data(sample_data_frame& df)
{
    for (unsigned i = 0; i < 65536; i+=7)
    {
        df.csv1.set(i, 4);
        df.csv2.set(i, 8);
        df.csv3.set(i, 17);
    }

    // rebuild Rank-Select index once data are loaded
    df.csv1.sync();
    df.csv2.sync();
    df.csv3.sync();
}

/// paranoiya check to make sure data frame matches pre-generated values
static
void test_data(sample_data_frame& df)
{
    for (unsigned i = 0; i < 65536; i+=7)
    {
        auto v1 = df.csv1.get(i);
        assert(v1 == 4);
        auto v2 = df.csv2.get(i);
        assert(v2 == 8);
        auto v3 = df.csv3.get(i);
        assert(v3 == 17);
    }
}

/**
    Copy buffer content into the buffer
    @internal
 */
template<typename SVLay>
unsigned char* copy_buffer(unsigned char* buf_ptr, const SVLay& sv_lay)
{
    auto s = sv_lay.size();
    ::memcpy(buf_ptr, sv_lay.buf(), s);
    return buf_ptr + s;
}

/// serialize with disabled XOR compression
///
static
void serialize_df0(const sample_data_frame& df,
                  std::vector<unsigned char>& buf,
                  csv_serializer_type& csv_ser)
{
    csv_ser.set_xor_ref(false); // disable XOR compression

    // buffers for serialization
    bm::sparse_vector_serial_layout<rsc_sparse_vector_u32> sv_lay1, sv_lay2, sv_lay3;

    // serialize all data-frame vectors in their natural order (1, 2, 3 ... N)
    //
    csv_ser.serialize(df.csv1, sv_lay1);
    csv_ser.serialize(df.csv2, sv_lay2);
    csv_ser.serialize(df.csv3, sv_lay3);

    size_t sz = (sizeof(size_t) * 3) +
                sv_lay1.size() + sv_lay2.size() + sv_lay3.size();

    buf.resize(sz);

    unsigned char* buf_ptr = buf.data(); // low level access to vector memory

    // serialize data-frame header with sizes of containers
    {
        auto s = sv_lay1.size();
        ::memcpy(buf_ptr, &s, sizeof(s));
        buf_ptr += sizeof(s);
    }
    {
        auto s = sv_lay2.size();
        ::memcpy(buf_ptr, &s, sizeof(s));
        buf_ptr += sizeof(s);
    }
    {
        auto s = sv_lay3.size();
        ::memcpy(buf_ptr, &s, sizeof(s));
        buf_ptr += sizeof(s);
    }
    // save all serialization buffers to one BLOB
    buf_ptr = copy_buffer(buf_ptr, sv_lay1);
    buf_ptr = copy_buffer(buf_ptr, sv_lay2);
    buf_ptr = copy_buffer(buf_ptr, sv_lay3);

}

/// Simple (individual) de-serialization of vectors in the data-frame
///
static
void deserialize_df0(sample_data_frame& df,
                     const std::vector<unsigned char>& buf,
                     csv_deserializer_type& csv_dser)
{
    assert(buf.size() > sizeof(size_t)*3);

    size_t sz1, sz2, sz3;
    const unsigned char* data_ptr = buf.data();

    // read the header to be able to calculate deserialization
    // offsets within BLOB
    //
    ::memcpy(&sz1, data_ptr, sizeof(size_t));
    data_ptr += sizeof(size_t);
    ::memcpy(&sz2, data_ptr, sizeof(size_t));
    data_ptr += sizeof(size_t);
    ::memcpy(&sz3, data_ptr, sizeof(size_t));
    data_ptr += sizeof(size_t);

    csv_dser.deserialize(df.csv1, data_ptr);
    data_ptr += sz1;

    csv_dser.deserialize(df.csv2, data_ptr);
    data_ptr += sz2;

    csv_dser.deserialize(df.csv3, data_ptr);
}




/// Data frame serialization using collaborative method (XOR compression)
///
static
void serialize_df2(const sample_data_frame& df,
                  std::vector<unsigned char>& buf,
                  csv_serializer_type& csv_ser)
{
    // reference vector(s) to keep all bit-planes for collaborative compresson
    csv_serializer_type::bv_ref_vector_type bv_ref;

    // Build the list of reference vectors participating in our data-frame
    // Important: add references in reverse (sic!) order (N ... 3, 2, 1)
    // Important: add all data frame vectors
    bv_ref.add_vectors(df.csv3.get_bmatrix());
    bv_ref.add_vectors(df.csv2.get_bmatrix());
    bv_ref.add_vectors(df.csv1.get_bmatrix());

    csv_ser.set_xor_ref(&bv_ref); // connect reference vector to serializer

    // buffers for serialization
    bm::sparse_vector_serial_layout<rsc_sparse_vector_u32> sv_lay1, sv_lay2, sv_lay3;

    // serialize all data-frame vectors in their natural order (1, 2, 3 ... N)
    //
    csv_ser.serialize(df.csv1, sv_lay1);
    csv_ser.serialize(df.csv2, sv_lay2);
    csv_ser.serialize(df.csv3, sv_lay3);

    size_t sz = (sizeof(size_t) * 3) +
                sv_lay1.size() + sv_lay2.size() + sv_lay3.size();

    buf.resize(sz);

    unsigned char* buf_ptr = buf.data(); // low level access to vector memory

    // serialize data-frame header with sizes of containers
    {
        auto s = sv_lay1.size();
        ::memcpy(buf_ptr, &s, sizeof(s));
        buf_ptr += sizeof(s);
    }
    {
        auto s = sv_lay2.size();
        ::memcpy(buf_ptr, &s, sizeof(s));
        buf_ptr += sizeof(s);
    }
    {
        auto s = sv_lay3.size();
        ::memcpy(buf_ptr, &s, sizeof(s));
        buf_ptr += sizeof(s);
    }
    // save all serialization buffers to one BLOB
    //
    buf_ptr = copy_buffer(buf_ptr, sv_lay1);
    buf_ptr = copy_buffer(buf_ptr, sv_lay2);
    buf_ptr = copy_buffer(buf_ptr, sv_lay3);


    // if serializer is re-used we need to disconnect it from the
    // current frame reference vectors

    csv_ser.set_xor_ref(nullptr);
}

/// Collaborative de-serialization of vectors in the data-frame
///
static
void deserialize_df2(sample_data_frame& df,
                     const std::vector<unsigned char>& buf,
                     csv_deserializer_type& csv_dser)
{
    assert(buf.size() > sizeof(size_t)*3);

    size_t sz1, sz2, sz3;
    const unsigned char* data_ptr = buf.data();

    // read the header to be able to calculate deserialization
    // offsets within BLOB
    //
    ::memcpy(&sz1, data_ptr, sizeof(size_t));
    data_ptr += sizeof(size_t);
    ::memcpy(&sz2, data_ptr, sizeof(size_t));
    data_ptr += sizeof(size_t);
    ::memcpy(&sz3, data_ptr, sizeof(size_t));
    data_ptr += sizeof(size_t);

    // reference vector for XOR deserialization
    //
    bm::sparse_vector_deserializer<rsc_sparse_vector_u32>::bv_ref_vector_type bv_ref;

    // ----------------------------------------------------------
    // first pass: build reference vectors, pre-construct vectors
    //
    csv_dser.deserialize_structure(df.csv1, data_ptr);
    data_ptr += sz1;
    csv_dser.deserialize_structure(df.csv2, data_ptr);
    data_ptr += sz2;
    csv_dser.deserialize_structure(df.csv3, data_ptr);

    // construct the reference vector
    // Important: add references in reverse (sic!) order (N ... 3, 2, 1)
    // Important: add all data frame vectors
    bv_ref.add_vectors(df.csv3.get_bmatrix());
    bv_ref.add_vectors(df.csv2.get_bmatrix());
    bv_ref.add_vectors(df.csv1.get_bmatrix());

    csv_dser.set_xor_ref(&bv_ref);


    // -------------------------------------------------------------
    // second pass: data deserialization
    //

    // get deserialization start-pointer again
    data_ptr = buf.data() + 3 * sizeof(size_t);


    // it is important that second pass uses 'false' as a third arument
    // to keep pre-created vectors structure, which is otherwise destroyed
    //

    csv_dser.deserialize(df.csv1, data_ptr, false);
    data_ptr += sz1;

    csv_dser.deserialize(df.csv2, data_ptr, false);
    data_ptr += sz2;

    csv_dser.deserialize(df.csv3, data_ptr, false);


    // disconnect deserialized from the reference vector
    //
    csv_dser.set_xor_ref(nullptr);

}


int main(void)
{
    try
    {
        std::vector<unsigned char> buf0, buf2;

        {
            sample_data_frame    df1;
            csv_serializer_type csv_ser;

            fill_test_data(df1);
            test_data(df1);

            serialize_df0(df1, buf0, csv_ser);
            cout << buf0.size() << endl;

            serialize_df2(df1, buf2, csv_ser);
            cout << buf2.size() << endl;
        }

        // run simple deserialization here, test results
        {
            sample_data_frame    df0;
            csv_deserializer_type csv_dser;

            deserialize_df0(df0, buf0, csv_dser);
            test_data(df0);
        }

        // collaborative deserialization (XOR decode)
        {
            sample_data_frame    df2;
            csv_deserializer_type csv_dser;

            deserialize_df2(df2, buf2, csv_dser);
            test_data(df2);
        }


    }
    catch(std::exception& ex)
    {
        std::cerr << ex.what() << std::endl;
        return 1;
    }
    
    

    return 0;
}

