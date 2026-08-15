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
  serializing it by default. The sample also shows range restore with
  serializer bookmarks, where the same requested range is restored into the
  master and all attached followers.

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
                      bool disable_xor = false,
                      bool use_bookmarks = false)
{
    bm::sparse_vector_serializer<SV> ser;
    ser.set_serialize_external_null(serialize_external_null);
    if (disable_xor)
        ser.disable_xor_compression();
    if (use_bookmarks)
        ser.set_bookmarks(true, 64);
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

template<typename SV>
static
void deserialize_vector_range(SV& sv,
                              const bm::sparse_vector_serial_layout<SV>& layout,
                              unsigned from,
                              unsigned to)
{
    bm::sparse_vector_deserializer<SV> deser;
    deser.deserialize_range(sv, layout.buf(), from, to);
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

    // Correct restore order: assemble the connected group, load the master,
    // then load followers whose external NULL plane was skipped.
    sample_data_frame df2;
    attach_followers(df2);
    const bvector_type* bv_null_df2 = df2.master_id.get_null_bvector();
    deserialize_vector(df2.master_id, master_layout);
    assert(bv_null_df2 == df2.master_id.get_null_bvector());
    assert(bv_null_df2 == df2.rsc_value.get_null_bvector());
    assert(bv_null_df2 == df2.plain_value.get_null_bvector());

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

/// Compare one restored range against the original data-frame columns.
static
void validate_range(const sample_data_frame& src,
                    const sample_data_frame& restored,
                    unsigned from,
                    unsigned to)
{
    assert(restored.rsc_value.is_null_external());
    assert(restored.plain_value.is_null_external());
    assert(restored.master_id.get_null_bvector() ==
           restored.rsc_value.get_null_bvector());
    assert(restored.master_id.get_null_bvector() ==
           restored.plain_value.get_null_bvector());

    const unsigned check_to = (to < src.master_id.size()) ?
                                to : unsigned(src.master_id.size() - 1);
    for (unsigned i = from; i <= check_to; ++i)
    {
        assert(src.master_id.is_null(i) == restored.master_id.is_null(i));
        assert(src.rsc_value.is_null(i) == restored.rsc_value.is_null(i));
        assert(src.plain_value.is_null(i) == restored.plain_value.is_null(i));

        unsigned src_v = 0, restored_v = 0;
        bool src_found = src.master_id.try_get(i, src_v);
        bool restored_found = restored.master_id.try_get(i, restored_v);
        assert(src_found == restored_found);
        if (src_found)
            assert(src_v == restored_v);

        src_found = src.rsc_value.try_get(i, src_v);
        restored_found = restored.rsc_value.try_get(i, restored_v);
        assert(src_found == restored_found);
        if (src_found)
            assert(src_v == restored_v);

        src_found = src.plain_value.try_get(i, src_v);
        restored_found = restored.plain_value.try_get(i, restored_v);
        assert(src_found == restored_found);
        if (src_found)
            assert(src_v == restored_v);
    }
}

/// Demonstrate follower edit and explicit RS rebuild for follow-up use.
///
/// This changes a value at an already NOT NULL position, so the shared NULL
/// shape is preserved. Shape-changing edits must be coordinated across all RSC
/// columns; otherwise their compressed value ranks no longer match the master.
static
void rs_rebuild_demo(sample_data_frame& df)
{
    const unsigned pos = 0;
    assert(!df.master_id.is_null(pos));

    df.plain_value.set(pos, 777);
    assert(!df.master_id.is_null(pos));
    assert(df.plain_value.get(pos) == 777);

    // Explicitly rebuild the master RS index and reattach the RSC follower so
    // subsequent rank/select operations can use the shared index again.
    df.master_id.invalidate_rs_index();
    df.rsc_value.invalidate_rs_index();
    df.master_id.sync();
    df.rsc_value.attach_null_bvector(df.master_id);

    assert(df.master_id.in_sync());
    assert(df.rsc_value.in_sync());
    assert(df.rsc_value.is_rs_index_external());
    assert(df.rsc_value.get_RS() == df.master_id.get_RS());
}

/// Demonstrate range restore for a connected group.
///
/// Bookmarks make range deserialization faster by adding skip points to the
/// serialized bit-vector streams. All connected columns must be restored with
/// the same range so the shared master NULL plane and value planes agree.
static
void range_restore_demo(const sample_data_frame& df)
{
    const unsigned from = 1024;
    const unsigned to = 8192;

    bm::sparse_vector_serial_layout<rsc_sparse_vector_u32> master_layout;
    bm::sparse_vector_serial_layout<rsc_sparse_vector_u32> rsc_layout;
    bm::sparse_vector_serial_layout<sparse_vector_u32> sv_layout;

    // Use bookmarks for fast range restore. Followers still skip external NULL.
    serialize_vector(df.master_id, master_layout, true, false, true);
    serialize_vector(df.rsc_value, rsc_layout, false, false, true);
    serialize_vector(df.plain_value, sv_layout, false, false, true);

    // Assemble the connected group before deserialization. The master owns the
    // shared NULL plane; followers require it when their serialized NULL is skipped.
    sample_data_frame df_range;
    attach_followers(df_range);
    const bvector_type* bv_null = df_range.master_id.get_null_bvector();
    deserialize_vector_range(df_range.master_id, master_layout, from, to);
    assert(bv_null == df_range.master_id.get_null_bvector());
    assert(bv_null == df_range.rsc_value.get_null_bvector());

    deserialize_vector_range(df_range.rsc_value, rsc_layout, from, to);
    deserialize_vector_range(df_range.plain_value, sv_layout, from, to);

    // Reattach the RSC follower after master sync so it can share the RS index.
    df_range.master_id.sync();
    df_range.rsc_value.attach_null_bvector(df_range.master_id);

    validate_range(df, df_range, from, to);
    cout << "Range restore [" << from << ", " << to << "] OK" << endl;
}

int main(void)
{
    try
    {
        sample_data_frame df;
        fill_test_data(df);
        validate_data(df);
        serialization_demo(df);
        range_restore_demo(df);

        rs_rebuild_demo(df);

        cout << "rscsample07 OK" << endl;
    }
    catch (std::exception& ex)
    {
        cerr << ex.what() << endl;
        return 1;
    }
    return 0;
}
