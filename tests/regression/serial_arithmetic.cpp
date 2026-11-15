/** Focused opcode-106 regression: physical payload limits, shared bit boundaries,
    both alphabets, stream/index traversal, OR semantics and level reset policy.
    Build with -Isrc; all checks remain active in release configurations.
*/
#include <cstring>
#include <iostream>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include "bm.h"
#include "bmserial.h"
#include "bmsparsevec.h"
#include "bmsparsevec_serial.h"
#include "bmserial_profile.h"
#include "bmdbg.h"

using BV = bm::bvector<>;
using Bytes = std::vector<unsigned char>;

/** Require a test invariant without depending on NDEBUG. */
static void check_at(bool valid, unsigned line)
{
    if (!valid)
    {
        throw std::runtime_error("Opcode 106 regression failed at line " + std::to_string(line));
    }
}

#define check(value) check_at((value), __LINE__)

/** Test forced alphabets and every shared-accumulator phase, including dry skip. */
static void test_kernel()
{
    bm::heap_vector<bm::word_t, BV::allocator_type, true> block, decoded;
    block.resize(bm::set_block_size);
    decoded.resize(bm::set_block_size);
    std::mt19937 rng(10681);
    Bytes storage(32768), joint(65536);
    for (unsigned pattern = 0; pattern < 5; ++pattern)
    {
        for (unsigned i = 0; i < bm::set_block_size; ++i)
        {
            block[i] = pattern == 0 ? 0 : pattern == 1 ? ~0u :
                       pattern == 2 ? 0xAAAAAAAAu : pattern == 3 ? rng() : rng() & rng();
        } // for i
        for (unsigned mode = 0; mode < 2; ++mode)
        {
            bm::encoder enc(storage.data(), storage.size());
            const unsigned count = bm::arithmetic_block_encoder::encode(enc, block.data(), bool(mode));
            bm::decoder source(storage.data());
            bm::arithmetic_block_decoder::decode(source, decoded.data());
            check(!::memcmp(block.data(), decoded.data(), 8192));
            check(source.size() == enc.size());
            for (unsigned phase = 0; phase < 32; ++phase)
            {
                bm::encoder dest(joint.data(), joint.size());
                {
                    bm::bit_out<bm::encoder> bits(dest);
                    for (unsigned bit = 0; bit < phase; ++bit)
                    {
                        bits.put_bit(bit & 1);
                    } // for bit
                    bm::decoder in(storage.data());
                    bm::bit_in<bm::decoder> body(in);
                    for (unsigned bit = 0; bit < 19 + count; ++bit)
                    {
                        bits.put_bit(body.get_bit());
                    } // for bit
                    bits.gamma(5);
                    bits.flush_if_full();
                    bits.put_bit(1);
                    bits.gamma(31);
                }
                for (unsigned dry = 0; dry < 2; ++dry)
                {
                    bm::decoder in(joint.data());
                    bm::bit_in<bm::decoder> bits(in);
                    for (unsigned bit = 0; bit < phase; ++bit)
                    {
                        check(bits.get_bit() == (bit & 1));
                    } // for bit
                    bm::arithmetic_block_decoder::decode(in, dry ? nullptr : decoded.data(), &bits);
                    if (!dry)
                    {
                        check(!::memcmp(block.data(), decoded.data(), 8192));
                    }
                    check(bits.gamma() == 5);
                    check(bits.get_bit() == 1);
                    check(bits.gamma() == 31);
                } // for dry
                // Both guard bits must still be consumed and validated when
                // their checks are combined. Reject every nonzero pair at
                // every enclosing bit-accumulator phase.
                for (unsigned guard = 1; guard < 4; ++guard)
                {
                    bm::decoder original(storage.data());
                    bm::bit_in<bm::decoder> input(original);
                    bm::encoder dest(joint.data(), joint.size());
                    {
                        bm::bit_out<bm::encoder> bits(dest);
                        for (unsigned bit = 0; bit < phase; ++bit)
                        {
                            bits.put_bit(bit & 1);
                        } // for bit
                        for (unsigned bit = 0; bit < 19 + count; ++bit)
                        {
                            unsigned value = input.get_bit();
                            if (bit >= 19 + count - 2)
                            {
                                value = (guard >> (bit - (19 + count - 2))) & 1;
                            }
                            bits.put_bit(value);
                        } // for bit
                    }
                    bool rejected = false;
                    try
                    {
                        bm::decoder source(joint.data());
                        bm::bit_in<bm::decoder> bits(source);
                        for (unsigned bit = 0; bit < phase; ++bit)
                        {
                            bits.get_bit();
                        } // for bit
                        bm::arithmetic_block_decoder::decode(source, decoded.data(), &bits);
                    }
                    catch (const std::exception&)
                    {
                        rejected = true;
                    }
                    check(rejected);
                } // for guard
            } // for phase
            if (pattern == 3)
            {
                // A short declared payload cannot borrow following metadata.
                bm::encoder corrupt(storage.data(), storage.size());
                {
                    bm::bit_out<bm::encoder> bits(corrupt);
                    bits.put_bits(34 | (mode << 18), 19);
                }
                bool rejected = false;
                try
                {
                    bm::decoder bad(storage.data());
                    bm::arithmetic_block_decoder::decode(bad, decoded.data());
                }
                catch (const std::exception&)
                {
                    rejected = true;
                }
                check(rejected);
            }
        } // for mode
    } // for pattern
}

/** Make a vector with nibble- and byte-biased blocks spanning superblocks. */
static BV make_vector()
{
    BV bv;
    std::mt19937 rng(10682);
    unsigned blocks[] = {0, 1, 256, 257};
    for (unsigned block : blocks)
    {
        for (unsigned byte = 0; byte < 8192; ++byte)
        {
            unsigned value = (rng() & 15) * 17;
            if (!(block & 1))
            {
                value = 0;
                for (unsigned bit = 0; bit < 8; ++bit)
                {
                    value |= unsigned(rng() % 10 < 4) << bit;
                } // for bit
            }
            for (unsigned bit = 0; bit < 8; ++bit)
            {
                if ((value >> bit) & 1)
                {
                    bv.set(BV::size_type(block) * 65536 + byte * 8 + bit);
                }
            } // for bit
        } // for byte
    } // for block
    bv.optimize();
    return bv;
}

/** Exercise ordinary/indexed/stream/set-operation readers and level identity. */
static void test_serialization()
{
    const BV source = make_vector();
    for (unsigned level : {6u, 7u, 8u, 9u, 10u})
    {
        bm::serializer<BV> writer;
        writer.set_compression_level(level);
        unsigned attempts[2] = {};
        writer.set_arithmetic_trial_observer(
            [](void* context, const bm::serializer<BV>::arithmetic_trial_observation& trial)
            {
                ++static_cast<unsigned*>(context)[trial.mode];
            }, attempts);
        check(writer.get_arithmetic_coding() == (level >= 7));
        writer.set_bookmarks(true, 4);
        writer.set_sparse_cutoff(0); // Exercise block codecs rather than sparse superblock routing.
        BV::statistics stats;
        source.calc_stat(&stats);
        Bytes off(stats.max_serialize_mem), on(stats.max_serialize_mem);
        writer.set_arithmetic_coding(false);
        off.resize(writer.serialize(source, off.data(), off.size()));
        writer.set_arithmetic_coding(true);
        on.resize(writer.serialize(source, on.data(), on.size()));
        if (level < 7)
        {
            check(on == off);
            continue;
        }
        bm::serialization_opcode_collector<BV> collector;
        bm::serialization_opcode_profile profile;
        collector.add_blob(on.data(), on.size(), profile);
        check(profile.reconciles());
        if (level == 7 || level == 8)
        {
            check(attempts[0] != 0 && attempts[1] == 0);
            for (unsigned context = 0; context < 4; ++context)
            {
                check(profile.arithmetic_modes[context][1] == 0);
                check(profile.arithmetic_0runs_modes[context][1] == 0);
            } // for context
        }
        if (level >= 9)
        {
            check(attempts[0] != 0 && attempts[1] != 0);
        }
        if (level == 8)
        {
            check(profile.opcode[bm::set_block_arithmetic].count != 0);
        }
        bm::deserializer<BV, bm::decoder> reader;
        BV actual;
        reader.deserialize(actual, on.data());
        check(actual == source);
        for (unsigned target_kind = 0; target_kind < 3; ++target_kind)
        {
            BV target;
            target.set(65000);
            if (target_kind == 1)
            {
                target.get_blocks_manager().deoptimize_block(0, 0, true);
            }
            else if (target_kind == 2)
            {
                target.set_range(0, 65535);
            }
            BV expected(target);
            expected |= source;
            reader.deserialize(target, on.data());
            check(target == expected);
        } // for target_kind
        bm::deserialization_index<BV> index;
        BV dry;
        reader.set_deserialization_index_construct(&index);
        check(reader.deserialize(dry, on.data()) == on.size());
        reader.unset_deserialization_index();
        reader.set_deserialization_index_use(&index);
        reader.set_range(256 * 65536 + 100, 256 * 65536 + 900);
        BV ranged, clipped, expected;
        reader.deserialize(ranged, on.data());
        clipped.copy_range(ranged, 256 * 65536 + 100, 256 * 65536 + 900);
        expected.copy_range(source, 256 * 65536 + 100, 256 * 65536 + 900);
        check(clipped == expected);
        std::stringstream stream(std::string(reinterpret_cast<const char*>(on.data()), on.size()));
        bm::streams_decoder input(stream);
        bm::streams_deserializer<BV> channel;
        BV streamed;
        check(channel.deserialize(streamed, input) && streamed == source);
        for (bm::set_operation op : {bm::set_OR, bm::set_AND, bm::set_XOR, bm::set_SUB})
        {
            BV target;
            target.set_range(1000, 80000);
            BV wanted(target);
            switch (op)
            {
            case bm::set_OR: wanted |= source; break;
            case bm::set_AND: wanted &= source; break;
            case bm::set_XOR: wanted ^= source; break;
            case bm::set_SUB: wanted -= source; break;
            default: break;
            } // switch op
            bm::operation_deserializer<BV> operation;
            operation.deserialize(target, on.data(), op);
            check(target == wanted);
        } // for op
        if (level == 8)
        {
            for (bm::set_operation op : {bm::set_COUNT_AND, bm::set_COUNT_XOR,
                 bm::set_COUNT_OR, bm::set_COUNT_SUB_AB, bm::set_COUNT_SUB_BA,
                 bm::set_COUNT_A, bm::set_COUNT_B})
            {
                BV target;
                target.set_range(1000, 80000);
                const BV before(target);
                BV expected(target);
                switch (op)
                {
                case bm::set_COUNT_AND: expected &= source; break;
                case bm::set_COUNT_XOR: expected ^= source; break;
                case bm::set_COUNT_OR: expected |= source; break;
                case bm::set_COUNT_SUB_AB: expected -= source; break;
                case bm::set_COUNT_SUB_BA: expected = source; expected -= target; break;
                case bm::set_COUNT_B: expected = source; break;
                default: break;
                } // switch op
                bm::operation_deserializer<BV> operation;
                check(operation.deserialize(target, on.data(), op) == expected.count());
                check(target == before);
            } // for op
        }
    } // for level
}

/** Check shared reference/gather arithmetic and indexed sparse-vector reads.
    @param gather False uses one source plane; true alternates two per wave.
*/
static void test_xor(bool gather)
{
    using SV = bm::sparse_vector<unsigned, BV>;
    SV source;
    BV selection;
    std::mt19937 rng(10683);
    for (unsigned block : {0u, 1u, 256u, 257u})
    {
        const unsigned start = block * 65536;
        for (unsigned bit = 0; bit < 65536; ++bit)
        {
            const unsigned r1 = bit < 49152 ? rng() & 1 : 0;
            const unsigned r2 = gather && bit < 49152 ? rng() & 1 : 0;
            const unsigned ref = gather && ((bit / 1024) & 1) ? r2 : r1;
            const unsigned low = bit < 49152 ? ref ^ unsigned(rng() % 10 < 4) : rng() & 1;
            source.set(start + bit, low | (r1 << 1) | (r2 << 2));
        } // for bit
        selection.set_range(start + 600, start + 700);
    } // for block
    source.optimize();
    for (unsigned level : {8u, 9u, 10u})
    {
        bm::sparse_vector_serializer<SV> writer;
        writer.set_compression_level(level);
        writer.set_xor_ref(true);
        writer.set_bookmarks(true, 4);
        writer.get_bv_serializer().set_xor_gather(gather);
        bm::sparse_vector_serial_layout<SV> layout, baseline;
        writer.get_bv_serializer().set_arithmetic_coding(false);
        writer.serialize(source, baseline);
        writer.get_bv_serializer().set_arithmetic_coding(true);
        writer.serialize(source, layout);
        check(layout.size() <= baseline.size());
        bm::sparse_vector_deserializer<SV> reader;
        SV actual;
        reader.deserialize(actual, layout.data());
        check(actual.equal(source));
        typename bm::sparse_vector_deserializer<SV>::deserialization_index_type index;
        reader.construct_deserialization_index(index, layout.data());
        bm::serialization_opcode_collector<BV> collector;
        bm::serialization_opcode_profile profile;
        for (unsigned row = 0; row < index.rows(); ++row)
        {
            const auto* plane = index.get_row(row);
            if (plane && plane->serialized_size())
            {
                collector.add_blob(layout.data() + index.plane_offset(row), plane->serialized_size(), profile);
            }
        } // for row
        check(profile.reconciles());
        const unsigned context = gather ? 3 : 2;
        if (!profile.arithmetic_modes[context][0] && !profile.arithmetic_modes[context][1])
        {
            std::cerr << "No shared arithmetic at level=" << level << " gather=" << gather << '\n';
            bm::print_serialization_opcode_profile(std::cerr, profile);
        }
        check(profile.arithmetic_modes[context][0] + profile.arithmetic_modes[context][1] != 0);
        reader.set_deserialization_index(&index);
        reader.set_deserialization_index_use();
        SV gathered, expected;
        reader.deserialize(gathered, layout.data(), selection);
        expected.resize(source.size());
        for (auto en = selection.first(); en.valid(); ++en)
        {
            expected.set(*en, source.get(*en));
        } // for en
        check(gathered.equal(expected));
    } // for level
}

/** Verify alphabet overrides, preset resets and both-mode decoding. */
static void test_mode_selection()
{
    using Mode = bm::arithmetic_coding_mode;
    bm::serializer<BV> writer;
    check(writer.get_arithmetic_coding_mode() == Mode::all);
    writer.set_compression_level(7);
    check(writer.get_arithmetic_coding_mode() == Mode::nibble);
    writer.set_arithmetic_coding_mode(Mode::byte);
    writer.set_compression_level(8);
    check(writer.get_arithmetic_coding_mode() == Mode::nibble);
    bool rejected = false;
    try
    {
        writer.set_arithmetic_coding_mode(static_cast<Mode>(99));
    }
    catch (const std::exception&)
    {
        rejected = true;
    }
    check(rejected && writer.get_arithmetic_coding_mode() == Mode::nibble);
    for (unsigned level : {9u, 10u})
    {
        writer.set_arithmetic_coding_mode(Mode::byte);
        writer.set_compression_level(level);
        check(writer.get_arithmetic_coding_mode() == Mode::all);
    } // for level
    const BV source = make_vector();
    BV::statistics stats;
    source.calc_stat(&stats);
    for (unsigned level : {7u, 8u})
    {
        for (Mode mode : {Mode::nibble, Mode::byte, Mode::all})
        {
            writer.set_compression_level(level);
            writer.set_sparse_cutoff(0);
            writer.set_arithmetic_coding_mode(mode);
            check(writer.get_arithmetic_coding_mode() == mode);
            unsigned attempts[2] = {};
            writer.set_arithmetic_trial_observer(
                [](void* context, const bm::serializer<BV>::arithmetic_trial_observation& trial)
                {
                    ++static_cast<unsigned*>(context)[trial.mode];
                }, attempts);
            Bytes data(stats.max_serialize_mem);
            data.resize(writer.serialize(source, data.data(), data.size()));
            writer.set_arithmetic_trial_observer(nullptr, nullptr);
            if (mode != Mode::byte)
            {
                check(attempts[0] != 0);
            }
            if (mode != Mode::nibble)
            {
                check(attempts[1] != 0);
            }
            if (mode != Mode::all)
            {
                const unsigned forbidden = mode == Mode::nibble ? 1 : 0;
                check(attempts[forbidden] == 0);
                bm::serialization_opcode_collector<BV> collector;
                bm::serialization_opcode_profile profile;
                collector.add_blob(data.data(), data.size(), profile);
                check(profile.reconciles());
                for (unsigned context = 0; context < 4; ++context)
                {
                    check(profile.arithmetic_modes[context][forbidden] == 0);
                    check(profile.arithmetic_0runs_modes[context][forbidden] == 0);
                } // for context
            }
            BV restored;
            bm::deserialize(restored, data.data());
            check(restored == source);
        } // for mode
    } // for level
}

/** Run focused shared-bit and serialization coverage. */
int main()
{
    test_mode_selection();
    test_kernel();
    test_serialization();
    test_xor(false);
    test_xor(true);
    std::cout << "Arithmetic opcode 106: alphabets, payload limits, all shared phases, indexes, set/count operations and L6 identity PASS\n";
}
