// Focused opcode 107 regression. Build with -std=c++17 -Isrc.
#include "bm.h"
#include "bmserial.h"
#include "bmserial_profile.h"
#include <cstring>
#include <stdexcept>
#include <vector>

/** Verify copied-word reconstruction and the virtual-zero shared boundary. */
int main()
{
    typedef bm::heap_vector<bm::word_t, bm::standard_allocator, true> words_type;
    words_type copied, restored, expected;
    copied.resize(1870);
    restored.resize(bm::set_block_size);
    expected.resize(bm::set_block_size);
    bm::bit_block_set(expected.data(), 0);
    unsigned state = 107;
    for (unsigned i = 0; i < copied.size(); ++i)
    {
        state = state * 1664525u + 1013904223u;
        copied[i] = state;
        expected[31 + i] = state;
    } // for i
    std::vector<unsigned char> storage(32768), joint(32768);
    bm::heap_vector<unsigned, bm::standard_allocator, true> positions;
    for (unsigned mode = 0; mode < 2; ++mode)
    {
        bm::encoder enc(storage.data(), storage.size());
        enc.put_8(bm::set_block_bit_0runs_arithmetic);
        enc.put_8(0);
        enc.put_16(31);
        enc.put_16(1870);
        enc.put_16(147);
        const unsigned payload = bm::arithmetic_words_encoder::encode(
            enc, copied.data(), unsigned(copied.size()), bool(mode), false);
        bm::decoder source(storage.data());
        if (source.get_8() != bm::set_block_bit_0runs_arithmetic)
            throw std::runtime_error("Wrong opcode");
        bm::arithmetic_0runs_decoder::decode(source, restored.data(), positions);
        if (source.size() != enc.size() ||
            ::memcmp(restored.data(), expected.data(), 8192))
            throw std::runtime_error("Copied-word inverse failed");
        bm::decoder skipped(storage.data() + 1);
        bm::arithmetic_0runs_decoder::decode(skipped, nullptr, positions);
        if (skipped.size() + 1 != enc.size())
            throw std::runtime_error("Copied-word skip crossed the boundary");
        for (unsigned phase : {0u, 1u, 17u, 31u})
        {
            bm::encoder dest(joint.data(), joint.size());
            {
                bm::bit_out<bm::encoder> output(dest);
                for (unsigned i = 0; i < phase; ++i)
                    output.put_bit(i & 1u);
                bm::decoder body(storage.data() + 1);
                bm::bit_in<bm::decoder> input(body);
                for (unsigned i = 0; i < 7 * 8 + 19 + payload; ++i)
                    output.put_bit(input.get_bit());
                output.gamma(7);
            }
            bm::decoder input(joint.data());
            bm::bit_in<bm::decoder> bits(input);
            for (unsigned i = 0; i < phase; ++i)
            {
                if (bits.get_bit() != (i & 1u))
                    throw std::runtime_error("Shared prefix changed");
            } // for i
            bm::arithmetic_0runs_decoder::decode(input, restored.data(), positions, &bits);
            if (::memcmp(restored.data(), expected.data(), 8192) || bits.gamma() != 7)
                throw std::runtime_error("Shared coordinate boundary changed");
        } // for phase
        storage[2] = storage[3] = 0;
        bool rejected = false;
        try
        {
            bm::decoder invalid(storage.data() + 1);
            bm::arithmetic_0runs_decoder::decode(invalid, restored.data(), positions);
        }
        catch (const std::exception&)
        {
            rejected = true;
        }
        if (!rejected)
            throw std::runtime_error("Zero-length run accepted");
    } // for mode

    typedef bm::bvector<> BV;
    BV source;
    for (unsigned block_number : {0u, 257u})
    {
        bm::word_t* block = source.get_blocks_manager().make_bit_block(block_number);
        bm::bit_block_set(block, 0);
        unsigned value = 107 + block_number;
        for (unsigned i = 31; i < 1901; ++i)
        {
            value = value * 1664525u + 1013904223u;
            block[i] = (value & 0x0f0f0f0fu) | 0x30303030u;
        } // for i
    } // for block_number
    bm::serializer<BV> writer;
    writer.set_compression_level(7);
    writer.set_sparse_cutoff(0);
    writer.set_bookmarks(true, 4);
    BV::statistics stats;
    source.calc_stat(&stats);
    std::vector<unsigned char> blob(stats.max_serialize_mem);
    blob.resize(writer.serialize(source, blob.data(), blob.size()));
    bm::serialization_opcode_collector<BV> collector;
    bm::serialization_opcode_profile profile;
    collector.add_blob(blob.data(), blob.size(), profile);
    if (profile.opcode[bm::set_block_bit_0runs_arithmetic].count != 2)
        throw std::runtime_error("Opcode 107 was not selected twice");
    if (profile.arithmetic_0runs_modes[0][1])
    {
        throw std::runtime_error("L7 selected byte-mode copied-word arithmetic");
    }
    bm::deserializer<BV, bm::decoder> reader;
    BV restored_vector;
    reader.deserialize(restored_vector, blob.data());
    if (restored_vector != source)
        throw std::runtime_error("Vector inverse failed");
    bm::deserialization_index<BV> index;
    reader.set_deserialization_index_construct(&index);
    BV traversed;
    reader.deserialize(traversed, blob.data());
    reader.unset_deserialization_index();
    reader.set_deserialization_index_use(&index);
    reader.set_range(257u * 65536u, 258u * 65536u - 1u);
    BV ranged;
    reader.deserialize(ranged, blob.data());
    BV wanted_range;
    wanted_range.copy_range(source, 257u * 65536u, 258u * 65536u - 1u);
    BV actual_range;
    actual_range.copy_range(ranged, 257u * 65536u, 258u * 65536u - 1u);
    if (actual_range != wanted_range)
        throw std::runtime_error("Indexed range inverse failed");
    for (bm::set_operation operation : {bm::set_OR, bm::set_AND,
                                        bm::set_XOR, bm::set_SUB})
    {
        BV target;
        target.set_range(500, 80000);
        BV expected_target(target);
        switch (operation)
        {
        case bm::set_OR: expected_target |= source; break;
        case bm::set_AND: expected_target &= source; break;
        case bm::set_XOR: expected_target ^= source; break;
        case bm::set_SUB: expected_target -= source; break;
        default: break;
        } // switch operation
        bm::operation_deserializer<BV> set_reader;
        set_reader.deserialize(target, blob.data(), operation);
        if (target != expected_target)
            throw std::runtime_error("Set operation inverse failed");
    } // for operation
}
