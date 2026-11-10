#include "bm.h"
#include "bmserial.h"
#include "bmserial_profile.h"
#include "bmfio.h"
#include <cassert>
#include <iostream>
#include <sstream>
#include <thread>
#include <vector>

using BV = bm::bvector<>;

struct record_reader : bm::deseriaizer_base<bm::decoder, unsigned>
{
    using bm::deseriaizer_base<bm::decoder, unsigned>::read_bic_sb_arr;
};

static void malformed(unsigned mode, bool negative, unsigned code)
{
    unsigned char storage[512] = {};
    bm::encoder enc(storage, sizeof(storage));
    {
        bm::bit_out<bm::encoder> out(enc);
        out.put_32_no(0); out.gamma(1); out.put_bits(0,2); out.put_bits(mode,2);
        out.gamma(64); out.put_24_no(0); out.put_24_no(100); out.gamma(1);
        out.put_bit(1); out.put_bit(0); out.put_bit(0); // three raw tail flags
        if (mode == 2) out.put_bit(negative);
        out.gamma(code);
    }
    bm::decoder dec(storage);
    record_reader reader;
    unsigned values[64] = {}, sb = 0;
    bool rejected = false;
    try { reader.read_bic_sb_arr(dec, bm::set_sblock_bienc_adaptive, values, &sb); }
    catch (const std::logic_error&) { rejected = true; }
    assert(rejected);
}

static void check(unsigned len, unsigned offset, unsigned low)
{
    BV original;
    unsigned pos = offset + 10;
    original.set(pos);
    for (unsigned i = 1; i < len; ++i)
    {
        // Repeated and increasing/decreasing corrections, with partial tails.
        const unsigned gaps[] = {low, 96, 128, 96, low, 192};
        pos += gaps[(i / 128) % 6] + (i % 11 == 0 ? 3u : 0u);
        original.set(pos);
    }
    original.optimize();
    bm::serializer<BV> ordinary;
    ordinary.set_compression_level(6);
    bm::serializer<BV>::buffer plain;
    ordinary.serialize(original, plain);
    bm::serializer<BV> writer;
    writer.set_compression_level(6);
    writer.set_sblock_adaptive_options(true);
    bm::serializer<BV>::buffer bytes;
    writer.serialize(original, bytes);
    assert(bytes.size() <= plain.size());
    BV decoded;
    bm::deserialize(decoded, bytes.data());
    assert(decoded == original);
    bm::serialization_opcode_collector<BV> collector;
    bm::serialization_opcode_profile profile;
    collector.add_blob(bytes.data(),bytes.size(),profile);
    assert(profile.reconciles());
    if (len >= 4093) assert(profile.opcode[92].count == 1);
    std::stringstream stream(std::ios::in | std::ios::out | std::ios::binary);
    bm::streams_encoder output(stream);
    assert(writer.serialize_and_finish(original, output));
    stream.seekg(0);
    bm::streams_decoder input(stream);
    bm::streams_deserializer<BV> reader;
    BV streamed;
    assert(reader.deserialize(streamed, input));
    assert(streamed == original);
}

int main()
{
    for (unsigned prefix = 0; prefix < 32; ++prefix)
    {
        unsigned char data[256] = {};
        bm::encoder enc(data, sizeof(data));
        const unsigned codes[] = {1, 1, 2, 3, 17, 65536, 0xffffff};
        {
            bm::bit_out<bm::encoder> out(enc);
            if (prefix) out.put_bits(0, prefix);
            unsigned bits = prefix;
            for (unsigned code : codes) { out.gamma(code); bits += bm::gamma_size(code); }
            assert(enc.size()*8 + out.pending_bits() == bits);
        }
        bm::decoder dec(data);
        bm::bit_in<bm::decoder> in(dec);
        if (prefix) assert(in.get_bits(prefix) == 0);
        for (unsigned code : codes) assert(in.gamma() == code);
    }
    malformed(3, false, 1); // reserved predictor
    malformed(2, true, 1); // previous-relative underflow
    malformed(2, false, bm::set_sub_total_bits); // overflow
    malformed(1, false, 1); // retired global-relative predictor
    {
        unsigned char storage[32] = {}; bm::encoder enc(storage, sizeof(storage));
        { bm::bit_out<bm::encoder> out(enc); out.put_32_no(0); out.gamma(64); }
        bm::decoder dec(storage); record_reader reader; unsigned values[64] = {}, sb = 0;
        bool rejected = false;
        try { reader.read_bic_sb_arr(dec, bm::set_sblock_bienc_adaptive, values, &sb); }
        catch (const std::logic_error&) { rejected = true; }
        assert(rejected); // retired fixed-step header
    }
    {
        // No complete eight-flag group: the old gate rejected this profitable
        // adaptive transform even though its complete record beats BIC.
        BV short_array;
        unsigned pos = 10;
        short_array.set(pos);
        for (unsigned i = 1; i < 96; ++i)
        { pos += i >= 32 && i < 64 ? 2048u : 2u; short_array.set(pos); }
        short_array.optimize();
        bm::serializer<BV> writer;
        writer.set_compression_level(6);
        writer.set_sblock_adaptive_options(true);
        bm::serializer<BV>::buffer bytes;
        writer.serialize(short_array, bytes);
        assert(writer.get_compression_stat()[bm::set_sblock_bienc_adaptive] == 1);
        BV restored;
        bm::deserialize(restored, bytes.data());
        assert(restored == short_array);

    }
    for (unsigned n : {64u, 96u, 257u, 519u, 4093u, 65000u}) check(n, 0, 2);
    std::vector<std::thread> workers;
    for (unsigned t = 0; t < 10; ++t)
        workers.emplace_back([t]() { check(4093 + t, 2u * bm::set_sub_total_bits, 1 + t % 3); });
    for (auto& t : workers) t.join();
    std::cout << "delta steps/predictors, exact padded costs, streaming and 10 workers PASS\n";
}
