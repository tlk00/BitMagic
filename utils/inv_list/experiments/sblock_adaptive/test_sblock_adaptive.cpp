#include "bm.h"
#include "bmserial.h"
#include "bmfio.h"
#include <cassert>
#include <iostream>
#include <sstream>
#include <thread>
#include <vector>

int main()
{
    typedef bm::bvector<> bv_type;
    bv_type original;
    unsigned pos = 10;
    original.set(pos);
    for (unsigned i = 1; i < 96; ++i)
    {
        pos += (i >= 32 && i < 64) ? 2048u : 2u;
        original.set(pos);
    }
    original.optimize();

    bm::serializer<bv_type> writer;
    writer.set_compression_level(6);
    writer.set_sblock_adaptive_options(true);
    bm::serializer<bv_type>::buffer encoded;
    writer.serialize(original, encoded);
    assert(writer.get_compression_stat()[bm::set_sblock_bienc_adaptive] == 1);

    bv_type decoded;
    bm::deserialize(decoded, encoded.data());
    assert(decoded == original);

    std::stringstream stream(std::ios::in | std::ios::out | std::ios::binary);
    bm::streams_encoder output(stream);
    assert(writer.serialize_and_finish(original, output));
    stream.seekg(0);
    bm::streams_decoder input(stream);
    bm::streams_deserializer<bv_type> reader;
    bv_type streamed;
    assert(reader.deserialize(streamed, input));
    assert(streamed == original);
    std::cout << "adaptive BIC short-window round trip PASS\n";

    // Every 32-position window changes its minimum gap. This gives several
    // all-one flag groups and exercises the new grouped opcode.
    bv_type alternating;
    pos = 10;
    alternating.set(pos);
    for (unsigned i = 1; i < 32 * 96; ++i)
    {
        pos += ((i / 32) & 1) ? 64u : 2u;
        alternating.set(pos);
    }
    alternating.optimize();
    bm::serializer<bv_type> ordinary_writer;
    ordinary_writer.set_compression_level(6);
    bm::serializer<bv_type>::buffer ordinary_bytes;
    ordinary_writer.serialize(alternating, ordinary_bytes);
    bm::serializer<bv_type> grouped_writer;
    grouped_writer.set_compression_level(6);
    grouped_writer.set_sblock_adaptive_options(true);
    bm::serializer<bv_type>::buffer grouped_bytes;
    grouped_writer.serialize(alternating, grouped_bytes);
    assert(grouped_bytes.size() < ordinary_bytes.size());
    assert(grouped_writer.get_compression_stat()[bm::set_sblock_bienc_adaptive] == 1);
    bv_type grouped_decoded;
    bm::deserialize(grouped_decoded, grouped_bytes.data());
    assert(grouped_decoded == alternating);
    std::stringstream grouped_stream(std::ios::in | std::ios::out | std::ios::binary);
    bm::streams_encoder grouped_output(grouped_stream);
    assert(grouped_writer.serialize_and_finish(alternating, grouped_output));
    grouped_stream.seekg(0);
    bm::streams_decoder grouped_input(grouped_stream);
    bm::streams_deserializer<bv_type> grouped_reader;
    bv_type grouped_streamed;
    assert(grouped_reader.deserialize(grouped_streamed, grouped_input));
    assert(grouped_streamed == alternating);
    std::cout << "grouped flags round trip PASS\n";

    std::vector<std::thread> workers;
    for (unsigned t = 0; t < 10; ++t)
        workers.emplace_back([&alternating]() {
            bm::serializer<bv_type> local_writer;
            local_writer.set_compression_level(6);
            local_writer.set_sblock_adaptive_options(true);
            for (unsigned i = 0; i < 20; ++i)
            {
                bm::serializer<bv_type>::buffer bytes;
                local_writer.serialize(alternating, bytes);
                assert(local_writer.get_compression_stat()[bm::set_sblock_bienc_adaptive] > 0);
                bv_type restored;
                bm::deserialize(restored, bytes.data());
                assert(restored == alternating);
            }
        });
    for (auto& worker : workers) worker.join();
    std::cout << "grouped flags 10-worker round trips PASS\n";
}
