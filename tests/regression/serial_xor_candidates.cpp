/** \file serial_xor_candidates.cpp
    Focused reference-XOR regression for explicit options and level defaults.
    Checks allocator-backed model reuse, disabled-mode bytes and round trips.
*/
#include <cstring>
#include <iostream>
#include <random>
#include <vector>
#include <sstream>
#include <stdexcept>
#include "bm.h"
#include "bmsparsevec.h"
#include "bmsparsevec_serial.h"
#include "bmdbg.h"

/** Verify compact masked references and indexed reads across superblocks. */
static void test_h64_gather()
{
    using BV = bm::bvector<>;
    using SV = bm::sparse_vector<unsigned, BV>;
    SV source;
    std::mt19937 random(8041);
    BV selection;
    unsigned selected_blocks[] = {0, 1, 256, 257};
    for (unsigned block : selected_blocks)
    {
        const unsigned start = block * bm::gap_max_bits;
        for (unsigned i = 0; i < bm::gap_max_bits; ++i)
        {
            const unsigned ref = random() & 1u;
            const unsigned low = i < 8192 ? ref : unsigned(random() % 128 == 0);
            source.set(start + i, (ref << 1) | low);
        } // for i
        selection.set_range(start + 1000, start + 1060);
    } // for block
    source.optimize();
    bm::sparse_vector_serializer<SV> writer;
    writer.set_xor_ref(true);
    writer.set_bookmarks(true, 4);
    writer.get_bv_serializer().set_xor_ref8_h64(true);
    bm::sparse_vector_serial_layout<SV> layout;
    writer.serialize(source, layout);
    bm::sparse_vector_deserializer<SV> reader;
    bm::sparse_vector_deserializer<SV>::deserialization_index_type index;
    reader.construct_deserialization_index(index, layout.data());
    bm::serialization_opcode_profile profile;
    bm::serialization_opcode_collector<BV> collector;
    collector.add_blob(layout.get_plane(0), index.get_row(0)->serialized_size(), profile);
    if (!profile.reconciles() || !profile.opcode[bm::set_block_xor_ref8_h64].count)
    {
        throw std::runtime_error("Synthetic H64 references missing");
    }
    // OR into empty, sparse GAP, dense BIT and full target blocks.
    BV target;
    target.set(bm::gap_max_bits + 9000);
    const unsigned dense = 256 * bm::gap_max_bits;
    for (unsigned i = 0; i < bm::gap_max_bits; i += 2)
    {
        target.set(dense + i);
    } // for i
    target.set_range(dense + bm::gap_max_bits, dense + 2 * bm::gap_max_bits - 1);
    target.optimize();
    BV expected_or(target);
    expected_or |= *source.get_bmatrix().get_row(0);
    bm::bv_ref_vector<BV> references;
    references.add(source.get_bmatrix().get_row(1), 1);
    bm::deserialize(target, layout.get_plane(0), nullptr, &references);
    if (target.compare(expected_or))
    {
        throw std::runtime_error("H64 decode into populated blocks failed");
    }
    reader.set_deserialization_index(&index);
    reader.set_deserialization_index_use();
    SV actual, expected;
    reader.deserialize(actual, layout.data(), selection);
    expected.resize(source.size());
    for (auto en = selection.first(); en.valid(); ++en)
    {
        expected.set(*en, source.get(*en));
    } // for en
    if (!actual.equal(expected))
    {
        throw std::runtime_error("H64 indexed gather failed");
    }
}

/** Check the four selected blocks use the gather layout for this level.
    @param profile Wire report for the completed target plane.
    @param level Serializer compression level.
*/
static void check_gather_opcodes(const bm::serialization_opcode_profile& profile,
                                 unsigned level)
{
    const unsigned legacy = bm::set_block_xor_gather;
    const unsigned shared = bm::set_block_xor_parallel_gather;
    const auto expected = level == 6 ? legacy : shared;
    const auto other = level == 6 ? shared : legacy;
    std::cout << "Gather level=" << level << " opcode 100=" << profile.opcode[legacy].count
              << " opcode 102=" << profile.opcode[shared].count << std::endl;
    if (!profile.reconciles() || profile.opcode[expected].count != 4 ||
        profile.opcode[other].count != 0)
    {
        throw std::runtime_error("Incorrect gather opcode or record count");
    }
}

/** Verify parallel gather references and indexed reads across superblocks.
    @param level Compression level for residual and ordinary block trials.
*/
static void test_parallel_gather(unsigned level)
{
    using BV = bm::bvector<>;
    using SV = bm::sparse_vector<unsigned, BV>;
    SV source;
    std::mt19937 random(8041);
    BV selection;
    unsigned selected_blocks[] = {0, 1, 256, 257};
    for (unsigned block : selected_blocks)
    {
        const unsigned start = block * bm::gap_max_bits;
        for (unsigned i = 0; i < bm::gap_max_bits; ++i)
        {
            const unsigned ref = random() & 1u;
            const unsigned other = random() & 1u;
            const unsigned low = ((i / 1024) % 2 ? other : ref) ^ unsigned(random() % 256 == 0);
            source.set(start + i, (other << 2) | (ref << 1) | low);
        } // for i
        selection.set_range(start + 1000, start + 1060);
    } // for block
    source.optimize();
    bm::sparse_vector_serializer<SV> writer;
    writer.set_xor_ref(true);
    writer.set_bookmarks(true, 4);
    writer.set_compression_level(level);
    if (level < 8)
    {
        writer.get_bv_serializer().set_xor_ref8_h64(true);
        writer.get_bv_serializer().set_xor_gather(true, true);
    }
    bm::sparse_vector_serial_layout<SV> layout;
    writer.serialize(source, layout);
    bm::sparse_vector_deserializer<SV> reader;
    bm::sparse_vector_deserializer<SV>::deserialization_index_type index;
    reader.construct_deserialization_index(index, layout.data());
    bm::serialization_opcode_profile profile;
    bm::serialization_opcode_collector<BV> collector;
    collector.add_blob(layout.get_plane(0), index.get_row(0)->serialized_size(), profile);
    check_gather_opcodes(profile, level);
    // External collection IDs can be offset or deliberately nonmonotonic.
    for (unsigned base : {50u, 500u})
    {
        bm::bv_ref_vector<BV> external;
        external.add(source.get_bmatrix().get_row(0), base);
        external.add(source.get_bmatrix().get_row(1), 100);
        external.add(source.get_bmatrix().get_row(2), 200);
        bm::sparse_vector_serializer<SV> external_writer;
        bm::sparse_vector_serializer<SV>::xor_sim_model_type model;
        external_writer.set_xor_ref(&external);
        external_writer.set_compression_level(level);
        external_writer.get_bv_serializer().set_xor_gather(true, true);
        external_writer.compute_sim_model(model, external, bm::xor_sim_params());
        external_writer.set_sim_model(&model);
        bm::sparse_vector_serial_layout<SV> external_layout;
        external_writer.serialize(source, external_layout);
        decltype(index) ext_index;
        reader.construct_deserialization_index(ext_index, external_layout.data());
        bm::serialization_opcode_profile observed;
        collector.set_ref_vectors(&external);
        collector.set_gather_target_id(base);
        collector.set_collect_ids();
        collector.add_blob(external_layout.get_plane(0), ext_index.get_row(0)->serialized_size(),
            observed, [base](const unsigned char* rec, size_t)
            {
                if ((*rec == bm::set_block_xor_gather || *rec == bm::set_block_xor_parallel_gather) && bool(rec[1] & 2) != (base == 50))
                    throw std::runtime_error("Incorrect gather coordinate mode");
            });
        check_gather_opcodes(observed, level);
        if (observed.ids != source.get_bmatrix().get_row(0)->count())
            throw std::runtime_error("External gather reconstruction failed");
        BV recovered;
        bm::deserializer<BV, bm::decoder> ext_reader;
        ext_reader.set_ref_vectors(&external);
        ext_reader.set_gather_target_id(base);
        ext_reader.deserialize(recovered, external_layout.get_plane(0));
        if (recovered.compare(*source.get_bmatrix().get_row(0)))
            throw std::runtime_error("External gather plane mismatch");
    } // for base
    // OR into empty, sparse GAP, dense BIT and full target blocks.
    BV target;
    target.set(bm::gap_max_bits + 9000);
    const unsigned dense = 256 * bm::gap_max_bits;
    for (unsigned i = 0; i < bm::gap_max_bits; i += 2)
    {
        target.set(dense + i);
    } // for i
    target.set_range(dense + bm::gap_max_bits, dense + 2 * bm::gap_max_bits - 1);
    target.optimize();
    BV expected_or(target);
    expected_or |= *source.get_bmatrix().get_row(0);
    bm::bv_ref_vector<BV> references;
    references.add(source.get_bmatrix().get_row(1), 1);
    references.add(source.get_bmatrix().get_row(2), 2);
    bm::deserializer<BV, bm::decoder> standalone;
    standalone.set_ref_vectors(&references);
    standalone.set_gather_target_id(0);
    standalone.deserialize(target, layout.get_plane(0));

    if (target.compare(expected_or))
    {
        throw std::runtime_error("H64 decode into populated blocks failed");
    }
    // Stream readers use the same record parser and must preserve boundaries.
    std::string bytes(reinterpret_cast<const char*>(layout.data()), layout.size());
    std::stringstream stream(bytes, std::ios::in | std::ios::binary);
    bm::streams_sparse_vector_deserializer<SV> stream_reader;
    SV stream_result;
    if (!stream_reader.deserialize(stream_result, stream) || !source.equal(stream_result))
    {
        throw std::runtime_error("Gather stream round trip failed");
    }
    reader.set_deserialization_index(&index);
    reader.set_deserialization_index_use();
    SV actual, expected;
    reader.deserialize(actual, layout.data(), selection);
    expected.resize(source.size());
    for (auto en = selection.first(); en.valid(); ++en)
    {
        expected.set(*en, source.get(*en));
    } // for en
    if (!actual.equal(expected))
    {
        throw std::runtime_error("H64 indexed gather failed");
    }
}

/** Exercise correlated planes at levels 6..8 and alternate the trial policy.
    @return Zero after all invariants pass; throws on failure.
*/
int main()
{
    test_h64_gather();
    for (unsigned level = 6; level <= 10; ++level)
    {
        test_parallel_gather(level);
    } // for level
    using SV = bm::sparse_vector<unsigned, bm::bvector<> >;
    SV source;
    std::mt19937 random(4591);
    for (unsigned i = 0; i < 131073; ++i)
    {
        unsigned value = 0;
        const unsigned base = random() & 1u;
        for (unsigned p = 0; p < 8; ++p)
        {
            value |= (base ^ unsigned(random() % (20 + p) == 0)) << p;
        } // for p
        source.set(i, value);
    } // for i
    source.optimize();
    bm::sparse_vector_serializer<SV> writer;
    writer.set_xor_ref(true);
    bm::id64_t alternatives = 0;
    std::vector<unsigned char> level6_bytes;
    for (unsigned level = 6; level <= 8; ++level)
    {
        writer.set_compression_level(level);
        bm::sparse_vector_serial_layout<SV> baseline, trial, restored_policy;
        writer.get_bv_serializer().set_xor_candidate_trials(false);
        writer.serialize(source, baseline);
        if (level == 6)
        {
            level6_bytes.assign(baseline.data(), baseline.data() + baseline.size());
        }
        writer.get_bv_serializer().set_xor_candidate_trials(true);
        writer.serialize(source, trial);
        alternatives += writer.get_bv_serializer().get_xor_trial_statistics().alternatives;
        SV decoded;
        if (bm::sparse_vector_deserialize(decoded, trial.data()) || !source.equal(decoded))
        {
            throw std::runtime_error("Four-candidate round trip failed");
        }
        if (trial.size() > baseline.size())
        {
            throw std::runtime_error("Four-candidate output grew");
        }
        writer.get_bv_serializer().set_xor_candidate_trials(false);
        writer.serialize(source, restored_policy);
        if (baseline.size() != restored_policy.size() ||
            std::memcmp(baseline.data(), restored_policy.data(), baseline.size()))
        {
            throw std::runtime_error("Disabled policy changed after reuse");
        }
    } // for level
    writer.set_compression_level(6);
    bm::sparse_vector_serial_layout<SV> reset_layout;
    writer.serialize(source, reset_layout);
    if (reset_layout.size() != level6_bytes.size() ||
        std::memcmp(reset_layout.data(), level6_bytes.data(), level6_bytes.size()))
    {
        throw std::runtime_error("Returning to level 6 changed its stream");
    }
    if (!alternatives)
    {
        throw std::runtime_error("Alternative trials were not exercised");
    }
    std::cout << "XOR candidate regression passed; alternatives=" << alternatives << '\n';
}
