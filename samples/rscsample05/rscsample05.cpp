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
  for a group of bit-transposed sparse vectors in a data-frame

  \sa bm::sparse_vector
  \sa bm::rsc_sparse_vector
  \sa bm::sparse_vector_serializer
  \sa bm::sparse_vector_deserializer

  \sa rscsample02.cpp
  \sa svsample02.cpp
*/

/*! \file rscsample05.cpp
    \brief Example: collaborative compression (XOR compression)
*/

#include <iostream>
#include <vector>
#include <cassert>

#include "bm.h"
#include "bmsparsevec.h"
#include "bmsparsevec_compr.h"
#include "bmsparsevec_serial.h"

using namespace std;

// type definitions
//

typedef bm::bvector<>                                     bvector_type;
typedef bm::sparse_vector<unsigned short, bvector_type>   sparse_vector_u16;
typedef bm::sparse_vector<unsigned, bvector_type>         sparse_vector_u32;
typedef bm::rsc_sparse_vector<unsigned, sparse_vector_u32 > rsc_sparse_vector_u32;

typedef bm::sparse_vector_serializer<rsc_sparse_vector_u32> csv_serializer_type;
typedef bm::sparse_vector_serializer<sparse_vector_u16> sv16_serializer_type;
typedef bm::sparse_vector_deserializer<rsc_sparse_vector_u32> csv_deserializer_type;
typedef bm::sparse_vector_deserializer<sparse_vector_u16> sv16_deserializer_type;

/// sample data-frame structure with multiple bm vectors
///
struct sample_data_frame
{
    rsc_sparse_vector_u32 csv1;
    rsc_sparse_vector_u32 csv2;
    rsc_sparse_vector_u32 csv3;

    sparse_vector_u16     sv0; // non-compressed sparse vector

    sample_data_frame()
        : sv0(bm::use_null)
    {}
};

/// Estimate raw data size
///
inline
size_t raw_data_size(const sample_data_frame& df)
{
    size_t sz = 0;
    if (df.csv1.size())
        sz += df.csv1.size() * sizeof(df.csv1.get(0));
    if (df.csv2.size())
        sz += df.csv2.size() * sizeof(df.csv2.get(0));
    if (df.csv3.size())
        sz += df.csv3.size() * sizeof(df.csv3.get(0));
    if (df.sv0.size())
        sz += df.sv0.size() * sizeof(df.sv0.get(0));
    return sz;
}

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
        df.sv0.set(i, i % 256);
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
        assert(v1 == 4); (void)v1;
        auto v2 = df.csv2.get(i);
        assert(v2 == 8); (void)v2;
        auto v3 = df.csv3.get(i);
        assert(v3 == 17); (void)v3;

        auto v0 = df.sv0.get(i);
        assert(v0 == i % 256); (void)v0;
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
                  csv_serializer_type& csv_ser,
                  sv16_serializer_type& sv16_ser)
{
    csv_ser.set_xor_ref(false); // disable XOR compression

    // buffers for serialization
    bm::sparse_vector_serial_layout<rsc_sparse_vector_u32> sv_lay1, sv_lay2, sv_lay3;
    bm::sparse_vector_serial_layout<sparse_vector_u16> sv_lay0;

    // serialize all data-frame vectors in their natural order (1, 2, 3 ... N)
    //
    csv_ser.serialize(df.csv1, sv_lay1);
    csv_ser.serialize(df.csv2, sv_lay2);
    csv_ser.serialize(df.csv3, sv_lay3);
    sv16_ser.serialize(df.sv0, sv_lay0);

    size_t sz = (sizeof(size_t) * 3) +
                sv_lay1.size() + sv_lay2.size() + sv_lay3.size() + sv_lay0.size();

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
    {
        auto s = sv_lay0.size();
        ::memcpy(buf_ptr, &s, sizeof(s));
        buf_ptr += sizeof(s);
    }

    // save all serialization buffers to one BLOB
    buf_ptr = copy_buffer(buf_ptr, sv_lay1);
    buf_ptr = copy_buffer(buf_ptr, sv_lay2);
    buf_ptr = copy_buffer(buf_ptr, sv_lay3);
    buf_ptr = copy_buffer(buf_ptr, sv_lay0);

}

/// Simple (individual) de-serialization of vectors in the data-frame
///
static
void deserialize_df0(sample_data_frame& df,
                     const std::vector<unsigned char>& buf,
                     csv_deserializer_type& csv_dser,
                     sv16_deserializer_type& sv16_dser)
{
    assert(buf.size() > sizeof(size_t)*3);

    size_t sz1, sz2, sz3, sz0;
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
    ::memcpy(&sz0, data_ptr, sizeof(size_t));
    data_ptr += sizeof(size_t);

    csv_dser.deserialize(df.csv1, data_ptr);
    data_ptr += sz1;

    csv_dser.deserialize(df.csv2, data_ptr);
    data_ptr += sz2;

    csv_dser.deserialize(df.csv3, data_ptr);
    data_ptr += sz3;

    sv16_dser.deserialize(df.sv0, data_ptr);
    data_ptr += sz0;

}




/// Data frame serialization using collaborative method (XOR compression)
///
static
void serialize_df2(const sample_data_frame& df,
                  std::vector<unsigned char>& buf,
                  csv_serializer_type& csv_ser,
                  sv16_serializer_type& sv16_ser)
{
    try
    {
        // reference vector(s) to keep all bit-planes for collaborative compresson
        csv_serializer_type::bv_ref_vector_type bv_ref;

        // Build the list of reference vectors participating in our data-frame
        // Important: add references in reverse (sic!) order (N ... 3, 2, 1)
        // Important: add all data frame vectors

        bv_ref.add_sparse_vector(df.sv0); // (!) last vectors is added first
        bv_ref.add_sparse_vector(df.csv3);
        bv_ref.add_sparse_vector(df.csv2);
        bv_ref.add_sparse_vector(df.csv1);

        csv_ser.set_xor_ref(&bv_ref); // connect reference vector to serializer
        sv16_ser.set_xor_ref(&bv_ref); // connect reference vector to sv16 serializer

        // compute XOR similarity model - it is common for all serializers
        // and must be added after set_xor_ref()
        //
        bm::xor_sim_params x_params;
        csv_serializer_type::xor_sim_model_type sim_model;
        csv_ser.compute_sim_model(sim_model, bv_ref, x_params);

        // add similarity model to each serializer
        //
        csv_ser.set_sim_model(&sim_model);
        sv16_ser.set_sim_model(&sim_model);

        // buffers for serialization
        bm::sparse_vector_serial_layout<rsc_sparse_vector_u32> sv_lay1, sv_lay2, sv_lay3;
        bm::sparse_vector_serial_layout<sparse_vector_u16> sv_lay0;

        // serialize all data-frame vectors in their natural order (1, 2, 3 ... N)
        //
        csv_ser.serialize(df.csv1, sv_lay1);
        csv_ser.serialize(df.csv2, sv_lay2);
        csv_ser.serialize(df.csv3, sv_lay3);
        sv16_ser.serialize(df.sv0,  sv_lay0);

        size_t sz = (sizeof(size_t) * 4) +
             sv_lay1.size() + sv_lay2.size() + sv_lay3.size() + sv_lay0.size();

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
        {
            auto s = sv_lay0.size();
            ::memcpy(buf_ptr, &s, sizeof(s));
            buf_ptr += sizeof(s);
        }

        // save all serialization buffers to one BLOB
        //
        buf_ptr = copy_buffer(buf_ptr, sv_lay1);
        buf_ptr = copy_buffer(buf_ptr, sv_lay2);
        buf_ptr = copy_buffer(buf_ptr, sv_lay3);
        buf_ptr = copy_buffer(buf_ptr, sv_lay0);


        // if serializer is re-used we need to disconnect it from the
        // current frame reference vectors
        csv_ser.set_xor_ref(nullptr);
        sv16_ser.set_xor_ref(nullptr);
    }
    catch (...)
    {
        // catch block is used to guarantee that serialiers
        // are not associated with a dead reference vector
        //
        csv_ser.set_xor_ref(nullptr);
        sv16_ser.set_xor_ref(nullptr);
        throw;
    }
}

/// Collaborative de-serialization of vectors in the data-frame
///
static
void deserialize_df2(sample_data_frame& df,
                     const std::vector<unsigned char>& buf,
                     csv_deserializer_type& csv_dser,
                     sv16_deserializer_type& sv16_dser)
{
    assert(buf.size() > sizeof(size_t)*3);

    try
    {
        size_t sz1, sz2, sz3, sz0;
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
        ::memcpy(&sz0, data_ptr, sizeof(size_t));
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
        data_ptr += sz3;
        sv16_dser.deserialize_structure(df.sv0, data_ptr);
        data_ptr += sz0;

        // construct the reference vector
        // Important: add references in reverse (sic!) order (N ... 3, 2, 1)
        // Important: add all data frame vectors
        //
        bv_ref.add_vectors(df.sv0.get_bmatrix()); // (!) last comes first
        bv_ref.add_vectors(df.csv3.get_bmatrix());
        bv_ref.add_vectors(df.csv2.get_bmatrix());
        bv_ref.add_vectors(df.csv1.get_bmatrix());

        // all de-serializers in the data-frame connect to the same set of refs
        csv_dser.set_xor_ref(&bv_ref);
        sv16_dser.set_xor_ref(&bv_ref);


        // -------------------------------------------------------------
        // second pass: data deserialization
        //

        // get deserialization start-pointer again
        data_ptr = buf.data() + (4 * sizeof(size_t));


        // it is important that second pass uses 'false' as a third arument
        // to keep pre-created vectors structure, which is otherwise destroyed
        //

        csv_dser.deserialize(df.csv1, data_ptr, false);
        data_ptr += sz1;

        csv_dser.deserialize(df.csv2, data_ptr, false);
        data_ptr += sz2;

        csv_dser.deserialize(df.csv3, data_ptr, false);
        data_ptr += sz3;

        sv16_dser.deserialize(df.sv0, data_ptr, false);
        data_ptr += sz0;


        // disconnect deserialized from the reference vector
        // before leaving the scope
        //
        csv_dser.set_xor_ref(nullptr);
        sv16_dser.set_xor_ref(nullptr);
    }
    catch (...)
    {
        // catch block is used to guarantee that de-serialiers
        // are not associated with a dead reference vector
        //
        csv_dser.set_xor_ref(nullptr);
        sv16_dser.set_xor_ref(nullptr);
        throw;
    }


}


int main(void)
{
    try
    {
        std::vector<unsigned char> buf0, buf2;

        {
            sample_data_frame    df1; // sample data-frame
            fill_test_data(df1);      // add some test data
            test_data(df1);

            size_t raw_size = raw_data_size(df1);
            cout << "raw size: " << raw_size << endl;

            // declare serializers for sparse vectors of different types
            csv_serializer_type    csv_ser;
            sv16_serializer_type   sv16_ser;

            serialize_df0(df1, buf0, csv_ser, sv16_ser);
            cout << "Plain serializarion: "  << buf0.size() << endl;

            serialize_df2(df1, buf2, csv_ser, sv16_ser);
            cout << "XOR data frame serialization: " << buf2.size() << endl;
        }

        // de-serialiers for sparse vectors of different types
        // please note that de-serializers are reusable
        //
        csv_deserializer_type csv_dser;
        sv16_deserializer_type sv16_dser;

        // run simple deserialization here, test results
        {
            sample_data_frame    df0; // empty data frame to read into

            deserialize_df0(df0, buf0, csv_dser, sv16_dser);

            test_data(df0); // check to make sure we are OK
        }

        // collaborative deserialization (XOR decode)
        {
            sample_data_frame    df2;

            deserialize_df2(df2, buf2, csv_dser, sv16_dser);

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

