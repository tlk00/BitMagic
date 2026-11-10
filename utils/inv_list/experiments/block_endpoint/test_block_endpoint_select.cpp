#include "bm.h"
#include "bmserial.h"
#include "bmserial_profile.h"
#include <cassert>
#include <iostream>
#include <random>
#include <vector>

int main()
{
    // Every delta16s branch, including accumulator padding boundaries.
    for (unsigned v = 0; v < 65536; ++v)
    {
        unsigned char data[32];
        bm::encoder enc(data, sizeof(data));
        bm::bit_out<bm::encoder> out(enc);
        out.delta16s(v);
        assert(out.pending_bits() == bm::delta16s_bits(v));
    }
    std::mt19937 sizes(19);
    for (unsigned round = 0; round < 2000; ++round)
    {
        std::vector<bm::gap_word_t> positions;
        unsigned stride = 1 + sizes() % 64;
        unsigned first = sizes() % 1024;
        unsigned count = 2 + sizes() % ((65535 - first) / stride);
        for (unsigned i = 0; i < count; ++i)
            positions.push_back(bm::gap_word_t(first + i * stride));
        for (unsigned layout = 0; layout < 2; ++layout)
        {
            unsigned char data[131072];
            bm::encoder enc(data, sizeof(data));
            bm::bit_out<bm::encoder> out(enc);
            unsigned tail = 65536 - positions.back();
            out.delta16s(count); out.put_bits(layout, 1);
            unsigned header = bm::delta16s_bits(count) + 1;
            if (layout) {
                out.delta16s(first); out.delta16s(tail);
                header += bm::delta16s_bits(first) + bm::delta16s_bits(tail);
            }
            const auto* arr = positions.data() + layout;
            unsigned n = count - 2 * layout;
            auto lo = bm::gap_word_t(layout ? first + 1 : 0);
            auto hi = bm::gap_word_t(layout ? positions.back() - 1 : 65535);
            out.flush_if_full(); out.bic_encode_u16(arr, n, lo, hi); out.flush();
            bm::bic_u16_bit_counter counter{header, ~0u};
            assert(bm::bic_encode_u16_cm_sink(counter, arr, n, lo, hi));
            assert(enc.size() == ((counter.bits + 31) / 32) * 4);
            // Exact tie and one-word-smaller budgets exercise early rejection.
            unsigned padded = unsigned(enc.size()) * 8;
            for (unsigned budget : {padded, padded - 32}) {
                bm::bic_u16_bit_counter bounded{header, budget};
                bool fits = header <= budget &&
                    bm::bic_encode_u16_cm_sink(bounded, arr, n, lo, hi);
                assert(fits == (counter.bits <= budget));
            }
        }
    }
    bm::serialization_opcode_collector<bm::bvector<> > collector;
    bm::serialization_opcode_profile profile;
    std::mt19937 rng(42);
    for (unsigned round = 0; round < 20; ++round)
    {
        bm::bvector<> source;
        for (unsigned i = 0; i < 1800; ++i)
            source.set(rng() % 65536);
        source.optimize();
        bm::serializer<bm::bvector<> > ordinary, selected;
        ordinary.set_compression_level(6);
        selected.set_compression_level(6);
        selected.set_block_endpoint_selection(true);
        bm::serializer<bm::bvector<> >::buffer old_blob, new_blob;
        ordinary.serialize(source, old_blob);
        selected.serialize(source, new_blob);
        assert(new_blob.size() <= old_blob.size());
        bm::bvector<> decoded;
        bm::deserialize(decoded, old_blob.data());
        assert(decoded == source);
        decoded.clear();
        bm::deserialize(decoded, new_blob.data());
        assert(decoded == source);
        collector.add_blob(new_blob.data(), new_blob.size(), profile);
    }
    std::cout << "Endpoint selection: whole-stream round trip and exact layout choice PASS\n";
}
