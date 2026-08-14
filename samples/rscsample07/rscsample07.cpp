/*
Copyright(c) 2002-2026 Anatoliy Kuznetsov(anatoliy_kuznetsov at yahoo.com)

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

/** \example rscsample07.cpp
  Example of a data-frame style group where one RSC sparse vector owns the
  NOT NULL plane and other sparse vectors attach to it.

  Shared NULL planes are useful when several columns are populated together and
  therefore have the same NULL/not-NULL shape. The leading vector owns the
  bit-vector lifetime. Secondary vectors reference the same bit-vector and skip
  serializing it by default.

  \sa bm::sparse_vector
  \sa bm::rsc_sparse_vector
  \sa bm::sparse_vector_serializer
  \sa bm::sparse_vector_deserializer
  \sa bm::sparse_vector_serializer::set_serialize_external_null

  \sa rscsample05.cpp
  \sa svsample04.cpp
*/

/*! \file rscsample07.cpp
    \brief Example: shared NULL plane for RSC sparse-vector data frames
*/

#include <iostream>
#include <cassert>
#include <stdexcept>

#include "bm.h"
#include "bmsparsevec.h"
#include "bmsparsevec_compr.h"
#include "bmsparsevec_serial.h"
#include "bmundef.h" /* clear the pre-proc defines from BM */

using namespace std;

typedef bm::bvector<> bvector_type;
typedef bm::sparse_vector<unsigned, bvector_type> sparse_vector_u32;
typedef bm::rsc_sparse_vector<unsigned, sparse_vector_u32> rsc_sparse_vector_u32;

/// Small data-frame with one RSC owner and two attached secondary columns.
struct sample_data_frame
{
    rsc_sparse_vector_u32 master_id;   ///< owns the NOT NULL plane
    rsc_sparse_vector_u32 rsc_value;   ///< RSC follower, may share RS index
    sparse_vector_u32     plain_value; ///< plain sparse-vector follower

    sample_data_frame()
        : plain_value(bm::use_null)
    {}
};

/// Attach secondary columns to the master's NOT NULL plane.
///
/// This is a post-construction assembly step. The master owns the bvector and
/// must outlive every attached follower. RSC-to-RSC attachment can also reuse
/// the master's rank-select index when the master is already synced.
static
void attach_followers(sample_data_frame& df)
{
    df.master_id.sync();
    df.rsc_value.attach_null_bvector(df.master_id);

    bvector_type* bv_null = df.master_id.get_null_bvector();
    assert(bv_null);
    df.plain_value.attach_null_bvector(*bv_null);
}

/// Generate columns with identical NULL shape and independent values.
static
void fill_test_data(sample_data_frame& df)
{
    rsc_sparse_vector_u32::back_insert_iterator master_bi =
        df.master_id.get_back_inserter();
    rsc_sparse_vector_u32::back_insert_iterator rsc_bi =
        df.rsc_value.get_back_inserter();

    for (unsigned i = 0; i < 65536 + 257; ++i)
    {
        const bool is_not_null = (i % 5u) != 2u;
        if (is_not_null)
        {
            master_bi.add(i + 11);
            rsc_bi.add((i * 3u) + 7u);
            df.plain_value.push_back((i * 5u) + 17u);
        }
        else
        {
            master_bi.add_null();
            rsc_bi.add_null();
            df.plain_value.push_back_null();
        }
    }
    master_bi.flush();
    rsc_bi.flush();

    attach_followers(df);
}

/// Validate values, NULL shape and attachment identity.
static
void validate_data(const sample_data_frame& df)
{
    assert(!df.master_id.is_null_external());
    assert(df.rsc_value.is_null_external());
    assert(df.plain_value.is_null_external());
    assert(df.master_id.get_null_bvector() == df.rsc_value.get_null_bvector());
    assert(df.master_id.get_null_bvector() == df.plain_value.get_null_bvector());

    rsc_sparse_vector_u32::const_iterator master_it = df.master_id.begin();
    for (unsigned i = 0; i < df.master_id.size(); ++i, ++master_it)
    {
        const bool is_null = ((i % 5u) == 2u);
        assert(df.master_id.is_null(i) == is_null);
        assert(df.rsc_value.is_null(i) == is_null);
        assert(df.plain_value.is_null(i) == is_null);
        assert(master_it.valid());
        assert(master_it.is_null() == is_null);

        sparse_vector_u32::const_iterator sv_it =
            df.plain_value.get_const_iterator(i);
        assert(sv_it.valid());
        assert(sv_it.is_null() == is_null);

        if (!is_null)
        {
            assert(df.master_id.get(i) == i + 11);
            assert(df.rsc_value.get(i) == (i * 3u) + 7u);
            assert(df.plain_value.get(i) == (i * 5u) + 17u);
        }
    }
}

template<typename SV>
static
void serialize_vector(const SV& sv,
                      bm::sparse_vector_serial_layout<SV>& layout,
                      bool serialize_external_null,
                      bool disable_xor = false)
{
    bm::sparse_vector_serializer<SV> ser;
    ser.set_serialize_external_null(serialize_external_null);
    if (disable_xor)
        ser.disable_xor_compression();
    ser.serialize(sv, layout);
}

template<typename SV>
static
void deserialize_vector(SV& sv,
                        const bm::sparse_vector_serial_layout<SV>& layout)
{
    bm::sparse_vector_deserializer<SV> deser;
    deser.deserialize(sv, layout.buf());
}

/// Demonstrate skipped external NULL serialization and required load order.
static
void serialization_demo(const sample_data_frame& df)
{
    bm::sparse_vector_serial_layout<rsc_sparse_vector_u32> master_layout;
    bm::sparse_vector_serial_layout<rsc_sparse_vector_u32> rsc_keep_layout;
    bm::sparse_vector_serial_layout<rsc_sparse_vector_u32> rsc_skip_layout;
    bm::sparse_vector_serial_layout<sparse_vector_u32> sv_keep_layout;
    bm::sparse_vector_serial_layout<sparse_vector_u32> sv_skip_layout;

    // The master is self-contained. Followers skip the external NULL plane.
    serialize_vector(df.master_id, master_layout, true);
    serialize_vector(df.rsc_value, rsc_keep_layout, true, true);
    serialize_vector(df.rsc_value, rsc_skip_layout, false, true);
    serialize_vector(df.plain_value, sv_keep_layout, true, true);
    serialize_vector(df.plain_value, sv_skip_layout, false, true);

    cout << "RSC follower keep NULL: " << rsc_keep_layout.size() << " bytes" << endl;
    cout << "RSC follower skip NULL: " << rsc_skip_layout.size() << " bytes" << endl;
    cout << "SV  follower keep NULL: " << sv_keep_layout.size() << " bytes" << endl;
    cout << "SV  follower skip NULL: " << sv_skip_layout.size() << " bytes" << endl;

    assert(rsc_skip_layout.size() < rsc_keep_layout.size());
    assert(sv_skip_layout.size() < sv_keep_layout.size());

    // A skipped-NULL follower cannot be deserialized as a standalone vector.
    try
    {
        sparse_vector_u32 bad_sv(bm::use_null);
        deserialize_vector(bad_sv, sv_skip_layout);
        assert(0);
    }
    catch (const std::logic_error&)
    {
        cout << "Skipped NULL follower rejected standalone deserialization" << endl;
    }

    // Correct restore order: master first, attach followers, then load followers.
    sample_data_frame df2;
    deserialize_vector(df2.master_id, master_layout);
    attach_followers(df2);

    bm::sparse_vector_serial_layout<rsc_sparse_vector_u32> rsc_skip_xor_layout;
    bm::sparse_vector_serial_layout<sparse_vector_u32> sv_skip_xor_layout;
    serialize_vector(df.rsc_value, rsc_skip_xor_layout, false);
    serialize_vector(df.plain_value, sv_skip_xor_layout, false);

    deserialize_vector(df2.rsc_value, rsc_skip_xor_layout);
    deserialize_vector(df2.plain_value, sv_skip_xor_layout);
    df2.master_id.sync();
    df2.rsc_value.attach_null_bvector(df2.master_id);

    validate_data(df2);

    sparse_vector_u32 owned_sv(bm::use_null);
    bm::sparse_vector_serial_layout<sparse_vector_u32> sv_full_layout;
    serialize_vector(df.plain_value, sv_full_layout, true);
    deserialize_vector(owned_sv, sv_full_layout);

    sparse_vector_u32::statistics owned_stat, shared_stat;
    owned_sv.calc_stat(&owned_stat);
    df2.plain_value.calc_stat(&shared_stat);
    assert(shared_stat.memory_used < owned_stat.memory_used);
    assert(shared_stat.bv_count + 1 == owned_stat.bv_count);
}

int main(void)
{
    try
    {
        sample_data_frame df;
        fill_test_data(df);
        validate_data(df);
        serialization_demo(df);

        // Shape-changing edits through followers update the shared master plane.
        // RSC containers must invalidate/rebuild RS state after such edits.
        df.plain_value.set(2, 777);
        assert(!df.master_id.is_null(2));
        df.master_id.invalidate_rs_index();
        df.rsc_value.invalidate_rs_index();
        df.master_id.sync();
        df.rsc_value.attach_null_bvector(df.master_id);

        cout << "rscsample07 OK" << endl;
    }
    catch (std::exception& ex)
    {
        cerr << ex.what() << endl;
        return 1;
    }
    return 0;
}
