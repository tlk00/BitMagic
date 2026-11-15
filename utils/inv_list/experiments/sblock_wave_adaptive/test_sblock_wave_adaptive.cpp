/** Focused opcode-108 scoring, reconstruction and reserved-bit regression.
    This executable is development-only and is not part of inv_list.
*/
#include "bm.h"
#include "bmserial.h"
#include "bmserial_profile.h"
#include <algorithm>
#include <cassert>
#include <iostream>
#include <random>
#include <vector>

/** Compare dry selection with full records and round-trip every valid layout.
    @param positions Source superblock positions, sorted and deduplicated locally.
    @param sb Absolute superblock index used in each candidate header.
*/
static void check(std::vector<unsigned> positions, unsigned sb)
{
    std::sort(positions.begin(), positions.end());
    positions.erase(std::unique(positions.begin(), positions.end()), positions.end());
    bm::sblock_wave_adaptive_model<> model;
    const unsigned n = unsigned(positions.size());
    const bm::id64_t incumbent = 10000000;
    const auto selected = bm::sblock_wave_adaptive_select(
        model, positions.data(), n, sb, incumbent);
    bm::sblock_wave_adaptive_choice best;
    best.bytes = incumbent;
    std::vector<unsigned char> bytes(2000000);
    std::vector<unsigned> restored(n);
    model.build(positions.data(), n);
    for (unsigned width = 0; width < 6; ++width)
    {
        model.wave_size = 128u << width;
        model.wave_count = (1u << 24) / model.wave_size;
        bm::sblock_wave_record record;
        if (!model.prepare(positions.data(), n, record))
        {
            continue;
        }
        for (unsigned mode = 0; mode < 16; ++mode)
        {
            if ((mode & 8) && (mode & 1))
            {
                continue;
            }
            bm::encoder enc(bytes.data(), bytes.size());
            enc.put_8(bm::set_sblock_bienc_wave_runs);
            {
                bm::bit_out<bm::encoder> out(enc);
                out.put_bits(width, 3);
                bm::sblock_wave_cost cost;
                auto position = [&]()
                {
                    return bm::id64_t(enc.size()) * 8 + out.pending_bits();
                };
                model.encode(out, position, sb, mode, cost);
            }
            if (enc.size() < best.bytes)
            {
                best.bytes = enc.size();
                best.width = width;
                best.mode = mode;
            }
            assert(!(bytes[1] & 8));
            bm::decoder dec(bytes.data() + 1);
            bm::bit_in<bm::decoder> in(dec);
            unsigned decoded_sb = 0;
            const unsigned count = bm::sblock_wave_adaptive_model<>::decode(
                in, restored.data(), decoded_sb);
            assert(count == n && restored == positions && decoded_sb == sb);
            assert(size_t(dec.get_pos() - bytes.data()) == enc.size());
            bytes[1] |= 8;
            bm::decoder invalid(bytes.data() + 1);
            bm::bit_in<bm::decoder> invalid_in(invalid);
            bool rejected = false;
            try
            {
                bm::sblock_wave_adaptive_model<>::decode(
                    invalid_in, restored.data(), decoded_sb);
            }
            catch (const std::logic_error&)
            {
                rejected = true;
            }
            assert(rejected);
        }
    }
    assert(selected.bytes == best.bytes && selected.width == best.width &&
           selected.mode == best.mode);
    const auto tie = bm::sblock_wave_adaptive_select(
        model, positions.data(), n, sb, best.bytes);
    assert(tie.bytes == best.bytes);
}

/** Exercise boundary, periodic and random inputs and actual serializer routing. */
int main()
{
    check({0, 1}, 0);
    check({0, 16777215}, 7);
    check({127, 128, 4095, 4096}, 255);
    std::vector<unsigned> positions;
    for (unsigned i = 0; i < 65535; ++i)
    {
        positions.push_back(i);
    }
    check(positions, 65535);
    positions.clear();
    for (unsigned i = 0; i < 65535; ++i)
    {
        positions.push_back(i * 256 + 3);
    }
    check(positions, 70000);
    std::mt19937 random(34129);
    for (unsigned trial = 0; trial < 30; ++trial)
    {
        positions.clear();
        for (unsigned i = 0; i < 2 + trial * 131; ++i)
        {
            positions.push_back(random() % (1u << 24));
        }
        check(positions, trial);
    }
    bm::bvector<> source;
    for (unsigned i = 0; i < 100; ++i)
    {
        source.set(1 + i * 163841);
    }
    source.optimize();
    bm::serializer<bm::bvector<>> writer;
    writer.set_compression_level(8);
    bm::serializer<bm::bvector<>>::buffer bytes;
    writer.serialize(source, bytes);
    bm::bvector<> restored;
    bm::deserialize(restored, bytes.data());
    assert(restored == source);
    bm::serialization_opcode_collector<bm::bvector<>> collector;
    bm::serialization_opcode_profile profile;
    collector.add_blob(bytes.data(), bytes.size(), profile);
    assert(profile.opcode[108].count && !profile.opcode[90].count &&
           !profile.opcode[97].count && profile.reconciles());
    std::cout << "Opcode 108: dry-score parity, all layouts, ties and extension rejection PASS\n";
}
