/*
Copyright(c) 2026 Anatoliy Kuznetsov(tolikkuznetsov66 at gmail.com)

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

/** \example svfsample05.cpp
  Example of how to search bm::sparse_vector_float<> with
  bm::sparse_vector_scanner<>, serialize the search result sets, and later use
  the restored result sets for deserialization-index assisted gather
  deserialization from a compressed sparse-vector BLOB.

  The sample models a two-stage application. Stage 1 has the full vector in
  memory, searches for interesting values, and serializes both the result sets
  and the vector. Stage 2 keeps only serialized BLOBs, restores the result sets,
  and gathers the selected float values without restoring the whole vector.

  \sa bm::sparse_vector_float
  \sa bm::sparse_vector_scanner<>::find_gt_float
  \sa bm::sparse_vector_scanner<>::find_lt_float
  \sa bm::sparse_vector_float_serializer
  \sa bm::sparse_vector_float_deserializer
  \sa bm::sparse_vector_float_deserialization_index
  \sa bm::serializer<bm::bvector<> >
*/

/*! \file svfsample05.cpp
    \brief Example: sparse_vector_float search result sets and indexed gather
           deserialization.
*/

#include <iostream>
#include <iomanip>
#include <limits>

#include <bm.h>
#include <bmserial.h>
#include <bmsparsevec_algo.h>
#include <bmsparsevec_float.h>
#include <bmsparsevec_float_serial.h>

typedef bm::sparse_vector_float<
            bm::sparse_vector<unsigned int, bm::bvector<> > > sparseVecFloat;
typedef sparseVecFloat::bvector_type bvector_type;
typedef bm::serializer<bvector_type>::buffer bv_blob_type;
typedef bm::sparse_vector_float_serial_layout<sparseVecFloat> svf_blob_type;
typedef bm::sparse_vector_float_deserializer<sparseVecFloat> svf_deserializer_type;
typedef svf_deserializer_type::deserialization_index_type svf_dindex_type;

const sparseVecFloat::size_type vector_size = 512 * 1024;
const float positive_threshold = 8.0f;
const float negative_threshold = -8.0f;

struct stage1_data
{
    svf_blob_type svf_blob;
    bv_blob_type positive_blob;
    bv_blob_type negative_blob;
    bv_blob_type anomaly_blob;
};

struct aggregate_stats
{
    bm::id64_t count = 0;
    double sum = 0.0;
    float min_value = std::numeric_limits<float>::max();
    float max_value = -std::numeric_limits<float>::max();
};

static float synthetic_value(sparseVecFloat::size_type i)
{
    // Most values stay close to zero, which is typical for sensor deltas,
    // residuals, normalized scores, or other anomaly-search inputs.
    float v = float(int(i % 97) - 48) * 0.015f;

    // Periodic positive and negative spikes are deterministic stand-ins for
    // anomalies. Keeping the distribution deterministic makes the sample output
    // stable and easier to compare when experimenting with the code.
    if (i && (i % 4096) == 0)
        v += 16.0f + float(i % 11) * 0.25f;

    if (i && (i % 6143) == 0)
        v -= 18.0f + float(i % 7) * 0.35f;

    return v;
}

static void serialize_bvector(const bvector_type& bv, bv_blob_type& blob)
{
    // Result sets are ordinary bit-vectors. Serialization turns them into
    // compact BLOBs which can be saved and re-used as gather masks later.
    bm::serializer<bvector_type> bvs;
    bvs.serialize(bv, blob);
}

static void deserialize_bvector(bvector_type& bv, const bv_blob_type& blob)
{
    bv.clear();
    bm::deserialize(bv, blob.buf());
}

static void print_result_set(const char* name,
                             const bvector_type& bv,
                             const bv_blob_type& blob)
{
    std::cout << "  " << name << ": " << bv.count()
              << " positions, serialized result-set BLOB = "
              << blob.size() << " bytes" << std::endl;
}

static double bytes_to_mb(size_t bytes)
{
    return double(bytes) / (1024.0 * 1024.0);
}

static void build_search_and_serialize(stage1_data& data)
{
    std::cout << "Stage 1: build float sparse vector" << std::endl;

    // A plain vector is the simplest baseline: one contiguous float per
    // logical row. This estimates only the payload, not allocator overhead.
    const size_t plain_vector_bytes = vector_size * sizeof(float);

    sparseVecFloat svf;
    for (sparseVecFloat::size_type i = 0; i < vector_size; ++i)
        svf.push_back(synthetic_value(i));

    BM_DECLARE_TEMP_BLOCK(tb)
    svf.optimize(tb);

    std::cout << "  vector size = " << svf.size() << " values" << std::endl;
    std::cout << "  plain std::vector<float> payload estimate = "
              << bytes_to_mb(plain_vector_bytes) << " MB" << std::endl;

    std::cout << "Stage 1: search for synthetic anomalies" << std::endl;

    bm::sparse_vector_scanner<sparseVecFloat> scanner;
    bvector_type positive_bv;
    bvector_type negative_bv;

    scanner.find_gt_float(svf, positive_threshold, positive_bv);
    scanner.find_lt_float(svf, negative_threshold, negative_bv);

    bvector_type anomaly_bv(positive_bv);
    anomaly_bv.bit_or(negative_bv);
    anomaly_bv.optimize(tb);

    serialize_bvector(positive_bv, data.positive_blob);
    serialize_bvector(negative_bv, data.negative_blob);
    serialize_bvector(anomaly_bv, data.anomaly_blob);

    print_result_set("positive spikes", positive_bv, data.positive_blob);
    print_result_set("negative spikes", negative_bv, data.negative_blob);
    print_result_set("all anomalies", anomaly_bv, data.anomaly_blob);

    std::cout << "Stage 1: serialize float sparse vector with bookmarks"
              << std::endl;

    bm::sparse_vector_float_serializer<sparseVecFloat> svf_serializer;
    svf_serializer.set_bookmarks(true, 64);
    svf_serializer.serialize(svf, data.svf_blob);

    std::cout << "  serialized float-vector BLOB = "
              << data.svf_blob.size() << " bytes" << std::endl;

    // svf leaves scope here. Stage 2 demonstrates the later retrieval phase
    // where only serialized BLOBs remain resident.
}

static aggregate_stats aggregate_selected_values(const sparseVecFloat& svf,
                                                 const bvector_type& mask_bv)
{
    aggregate_stats st;

    // for_each_sparse() accepts the same bit-vector as a filter, so the
    // aggregation visits exactly the values requested by gather deserialization.
    // The visitor receives value, NULL flag, and original sparse-vector index.
    auto accumulate = [&st](float v,
                            bool is_null,
                            sparseVecFloat::size_type /* idx */) -> int
    {
        if (is_null)
            return 0;

        st.sum += v;
        if (v < st.min_value)
            st.min_value = v;
        if (v > st.max_value)
            st.max_value = v;
        ++st.count;
        return 0;
    };

    bm::for_each_sparse(svf, mask_bv, accumulate);

    return st;
}

static void gather_and_report(const char* name,
                              svf_deserializer_type& deserializer,
                              const unsigned char* svf_blob,
                              const bv_blob_type& result_blob)
{
    bvector_type mask_bv;
    deserialize_bvector(mask_bv, result_blob);

    // The mask bit-vector is the gather layout: every set bit names one
    // original vector position to restore from the compressed float-vector BLOB.
    sparseVecFloat gathered_svf;
    deserializer.deserialize(gathered_svf, svf_blob, mask_bv);

    aggregate_stats st = aggregate_selected_values(gathered_svf, mask_bv);

    std::cout << "  " << name
              << ": count=" << st.count
              << ", sum=" << st.sum
              << ", avg=" << (st.count ? st.sum / double(st.count) : 0.0)
              << ", min=" << (st.count ? st.min_value : 0.0f)
              << ", max=" << (st.count ? st.max_value : 0.0f)
              << std::endl;
}

static void restore_and_compute(const stage1_data& data)
{
    const unsigned char* svf_blob = data.svf_blob.buf();

    std::cout << "Stage 2: build deserialization index" << std::endl;

    svf_dindex_type dindex;
    {
        // The index is built once for the serialized float-vector BLOB. It is
        // kept alive while the deserializer uses it for repeated gather reads.
        svf_deserializer_type index_builder;
        index_builder.construct_deserialization_index(dindex, svf_blob);
        dindex.optimize();
    }

    std::cout << "Stage 2: gather values from serialized BLOB" << std::endl;

    svf_deserializer_type deserializer;
    deserializer.set_deserialization_index(&dindex);
    deserializer.set_deserialization_index_use(true);

    gather_and_report("positive spikes", deserializer,
                      svf_blob, data.positive_blob);
    gather_and_report("negative spikes", deserializer,
                      svf_blob, data.negative_blob);
    gather_and_report("all anomalies", deserializer,
                      svf_blob, data.anomaly_blob);
}

int main(void)
{
    try
    {
        std::cout << std::fixed << std::setprecision(3);

        stage1_data data;
        build_search_and_serialize(data);

        std::cout << std::endl;

        restore_and_compute(data);
    }
    catch (std::exception& ex)
    {
        std::cerr << ex.what() << std::endl;
        return 1;
    }
    return 0;
}
