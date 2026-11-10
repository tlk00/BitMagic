#include "bm.h"
#include "bmserial.h"
#include "bmserial_profile.h"
#include <algorithm>
#include <cassert>
#include <cstring>
#include <iostream>
#include <random>
#include <vector>

// Independent decoder for the wire model, checked against the library reader.
static std::vector<unsigned> decode(const unsigned char* data, unsigned expected_sb)
{
    bm::decoder decoder(data);
    assert(decoder.get_8() == bm::set_sblock_bienc_wave1024);
    bm::bit_in<bm::decoder> in(decoder);
    unsigned mode = in.get_bits(2), sb_mode = in.get_bits(2), sb = 0;
    switch (sb_mode)
    {
    case 0: sb = in.gamma() - 1; break;
    case 1: sb = in.get_bits(8); break;
    case 2: sb = in.get_16_no(); break;
    default: sb = in.get_32_no(); break;
    }
    assert(sb == expected_sb);
    unsigned n = in.gamma(), count = in.gamma(), d = in.gamma() - 1;
    std::vector<unsigned> map(count), occupied, values(n);
    in.bic_decode_u32_cm(map.data(), count, 0, 16383);
    if (mode & 2)
    {
        for (unsigned wave = 0, j = 0; wave < 16384; ++wave)
            if (j < count && map[j] == wave) ++j;
            else occupied.push_back(wave);
    }
    else occupied = map;
    unsigned upper = unsigned(occupied.size()) * 1024 - 1 - (n - 1) * d;
    if (mode & 1)
    {
        values.front() = in.gamma() - 1;
        values.back() = upper - (in.gamma() - 1);
        if (n > 2) in.bic_decode_u32_cm(values.data() + 1, n - 2,
                                       values.front() + 1, values.back() - 1);
    }
    else in.bic_decode_u32_cm(values.data(), n, 0, upper);
    for (unsigned i = 0; i < n; ++i)
    {
        unsigned packed = values[i] + i*d;
        values[i] = occupied.at(packed / 1024) * 1024 + packed % 1024;
    }
    return values;
}

static void trial(const std::vector<unsigned>& input, unsigned sb)
{
    bm::sblock_wave1024_model<> model;
    bm::sblock_wave_record r;
    if (!model.prepare(input.data(), unsigned(input.size()), r))
    { assert(model.occupied.size() == 16384); return; }
    std::vector<unsigned char> buffer((input.size() + 16384)*4 + 4096);
    for (unsigned mode = 0; mode < 4; ++mode)
    {
        bm::encoder encoder(buffer.data(), buffer.size());
        encoder.put_8(bm::set_sblock_bienc_wave1024);
        auto& c = r.costs[mode];
        {
            bm::bit_out<bm::encoder> out(encoder);
            auto position = [&]() -> bm::id64_t
                { return bm::id64_t(encoder.size())*8 + out.pending_bits(); };
            model.encode(out, position, sb, mode, c);
        }
        c.bytes = encoder.size(); c.padding = c.bytes*8 - c.header - c.map - c.bic;
        assert(c.reconciles());
        assert(decode(buffer.data(), sb) == input);
        bm::decoder library_decoder(buffer.data() + 1);
        struct reader : bm::deseriaizer_base<bm::decoder, unsigned>
        { using bm::deseriaizer_base<bm::decoder, unsigned>::read_bic_sb_arr; };
        reader lib_reader;
        std::vector<unsigned> library_values(input.size());
        unsigned library_sb = 0;
        unsigned library_count = lib_reader.read_bic_sb_arr(
            library_decoder, bm::set_sblock_bienc_wave1024, library_values.data(), &library_sb);
        assert(library_count == input.size() && library_sb == sb &&
               library_values == input);
        if (c.bytes < r.costs[r.best].bytes) r.best = mode;
    }

}

int main()
{
    trial({0, 1}, 0); // d=0, trailing-only, n=2
    trial({0, 1023}, 7); // upper bound collapses to 1 after correction
    trial({16776192, 16777215}, 255); // leading only, final wave
    trial({0, 16777215}, 65535); // extreme interior gap
    trial({1023, 1024, 3072, 4095}, 70000); // cross-wave adjacency forces d=0
    std::vector<unsigned> dense;
    for (unsigned i = 0; i < 5000; ++i) dense.push_back(i);
    trial(dense, 8);
    std::vector<unsigned> waves;
    for (unsigned i = 0; i < 16384; ++i) waves.push_back(i*1024);
    trial(waves, 256); // no empty waves: excluded
    waves.erase(waves.begin()+8000); trial(waves, 256); // exactly one empty wave
    waves.erase(waves.begin()+8000); trial(waves, 256); // exactly two empty waves
    std::mt19937 rng(42);
    for (unsigned round = 0; round < 200; ++round)
    {
        std::vector<unsigned> a;
        for (unsigned i = 0, n = 2 + rng()%2000; i < n; ++i)
            a.push_back(rng()%16777216);
        std::sort(a.begin(), a.end()); a.erase(std::unique(a.begin(), a.end()), a.end());
        trial(a, round);
    }
    // Both the ordinary and experimental streams must decode to the same bits.
    unsigned selected = 0;
    for (unsigned pattern = 0; pattern < 5; ++pattern)
    {
        bm::bvector<> source;
        const unsigned counts[] = {2000, 2000, 100, 500, 1000};
        const unsigned steps[] = {5003, 5003, 163841, 17011, 16273};
        const unsigned bases[] = {0, 100, 1, 2, 100};
        for (unsigned i = 0; i < counts[pattern]; ++i)
            source.set(bases[pattern] + i*steps[pattern]);
        source.optimize();
        bm::serializer<bm::bvector<> > old_writer, wave_writer;
        old_writer.set_compression_level(6); wave_writer.set_compression_level(6);
        wave_writer.set_sblock_wave1024_options(true);
        bm::serializer<bm::bvector<> >::buffer old_bytes, wave_bytes;
        old_writer.serialize(source, old_bytes); wave_writer.serialize(source, wave_bytes);
        bm::serialization_opcode_collector<bm::bvector<>> collector;
        bm::serialization_opcode_profile profile;
        collector.add_blob(wave_bytes.data(),wave_bytes.size(),profile);
        assert(profile.reconciles());
        selected += unsigned(profile.opcode[90].count);
        assert(wave_bytes.size() <= old_bytes.size());
        if (!profile.opcode[90].count)
            assert(old_bytes.size() == wave_bytes.size() &&
                   !std::memcmp(old_bytes.data(), wave_bytes.data(), old_bytes.size()));
        bm::bvector<> decoded_old, decoded_wave;
        bm::deserialize(decoded_old, old_bytes.data());
        bm::deserialize(decoded_wave, wave_bytes.data());
        assert(decoded_old == source && decoded_wave == source);
    }
    assert(selected > 0);
    // Adaptive opcode 92 remains a rival. Emitting 90 may only improve its result.
    bm::bvector<> varying;
    unsigned at = 10;
    varying.set(at);
    for (unsigned i = 1; i < 32 * 96; ++i)
    {
        at += ((i / 32) & 1) ? 64u : 2u;
        varying.set(at);
    }
    varying.optimize();
    bm::serializer<bm::bvector<> > ordinary, adaptive, wave, combined;
    ordinary.set_compression_level(6); adaptive.set_compression_level(6);
    wave.set_compression_level(6); combined.set_compression_level(6);
    adaptive.set_sblock_adaptive_options(true);
    wave.set_sblock_wave1024_options(true);
    combined.set_sblock_adaptive_options(true); combined.set_sblock_wave1024_options(true);
    bm::serializer<bm::bvector<> >::buffer ordinary_blob, adaptive_blob,
                                           wave_blob, combined_blob;
    ordinary.serialize(varying, ordinary_blob);
    adaptive.serialize(varying, adaptive_blob); combined.serialize(varying, combined_blob);
    wave.serialize(varying, wave_blob);
    assert(combined_blob.size() <= ordinary_blob.size());
    assert(combined_blob.size() <= adaptive_blob.size());
    assert(combined_blob.size() <= wave_blob.size());
    bm::bvector<> checked;
    bm::deserialize(checked, ordinary_blob.data()); assert(checked == varying);
    checked.clear();
    bm::deserialize(checked, adaptive_blob.data()); assert(checked == varying);
    checked.clear();
    bm::deserialize(checked, wave_blob.data()); assert(checked == varying);
    checked.clear();
    bm::deserialize(checked, combined_blob.data()); assert(checked == varying);
    std::cout << "wave model: all four modes round-trip; serializer byte identity PASS\n";
}
