/** @file serial_sblock_hierarchical_full.cpp
    @brief Exercise opcode 109 (bm::set_sblock_bienc_wave_hierarchical) with
    a FULL-block map, residual counts above 65,535 and indexed decoding.
*/
#include <cassert>
#include "bm.h"
#include "bmserial.h"
#include "bmserial_profile.h"

/** Verify one FULL-aware level-8 candidate and its ordinary fallback.
    @param residual_count Number of set positions outside FULL blocks.
    @param two_full Include a second FULL block when true.
    @return True if opcode 109 (bm::set_sblock_bienc_wave_hierarchical) won.
*/
static bool check_full_case(unsigned residual_count, bool two_full)
{
    using bvector_type = bm::bvector<>;
    using serializer_type = bm::serializer<bvector_type>;
    bvector_type source;
    source.set_range(0, bm::gap_max_bits - 1);
    if (two_full)
    {
        source.set_range(7 * bm::gap_max_bits, 8 * bm::gap_max_bits - 1);
    }
    for (unsigned i = 0; i < residual_count; ++i)
    {
        source.set(2 * bm::gap_max_bits + 2 * i);
    } // for i
    const unsigned distant = 2 * bm::set_sub_total_bits + 17;
    source.set(distant);
    source.optimize();

    serializer_type candidate;
    candidate.set_compression_level(8);
    serializer_type::buffer full_blob;
    candidate.serialize(source, full_blob);

    bvector_type restored;
    bm::deserialize(restored, full_blob.data());
    assert(restored == source);
    assert(IS_FULL_BLOCK(restored.get_blocks_manager().get_block_ptr(0, 0)));

    bm::serialization_opcode_profile profile;
    bm::serialization_opcode_collector<bvector_type> collector;
    unsigned flagged = 0;
    collector.add_blob(full_blob.data(), full_blob.size(), profile,
        [&](const unsigned char* record, size_t)
        {
            if (record[0] == bm::set_sblock_bienc_wave_hierarchical)
            {
                bm::decoder decoder(record + 1);
                bm::bit_in<bm::decoder> bits(decoder);
                bits.get_bits(8); // Width, extension, index and payload fields.
                flagged += bits.get_bits(1);
            }
        });
    assert(profile.reconciles());

    bvector_type merged;
    merged.set(7 * bm::gap_max_bits + 9);
    merged.set(distant + 1);
    bvector_type expected(merged);
    expected |= source;
    bm::deserialize(merged, full_blob.data());
    assert(merged == expected);

    bm::deserializer<bvector_type, bm::decoder> indexed;
    bm::deserialization_index<bvector_type> index;
    bvector_type scan;
    indexed.set_deserialization_index_construct(&index);
    indexed.deserialize(scan, full_blob.data());
    // Index construction scans the records without materializing the vector.
    assert(!scan.any());
    indexed.unset_deserialization_index();
    indexed.set_deserialization_index_use(&index);
    indexed.set_range(2 * bm::gap_max_bits,
                      4 * bm::gap_max_bits - 1);
    bvector_type range_result, wanted, clipped;
    indexed.deserialize(range_result, full_blob.data());
    wanted.copy_range(source, 2 * bm::gap_max_bits,
                      4 * bm::gap_max_bits - 1);
    clipped.copy_range(range_result, 2 * bm::gap_max_bits,
                       4 * bm::gap_max_bits - 1);
    assert(clipped == wanted);
    bvector_type unexpected(range_result);
    unexpected -= source;
    assert(!unexpected.any());
    return flagged != 0;
}

/** Run focused FULL-map wire, OR and indexed-range checks. */
int main()
{
    assert(check_full_case(70000, false));
    check_full_case(60000, true);
    check_full_case(130000, false);
}
