// Focused library decoder regression. Build with -std=c++17 -Isrc.
#include "bm.h"
#include "bmserial.h"
#include "bmserial_profile.h"
#include <cassert>
#include <iostream>
#include <sstream>
#include <vector>

using BV = bm::bvector<>;
using Bytes = std::vector<unsigned char>;

static void normalized(const BV& actual)
{
    // Optimizing a copy must not find further BIT -> GAP/null/full conversions.
    BV expected(actual);
    expected.optimize();
    BV::statistics a, b;
    actual.calc_stat(&a); expected.calc_stat(&b);
    if (a.bit_blocks != b.bit_blocks || a.gap_blocks != b.gap_blocks)
        std::cerr << "representation actual=" << a.bit_blocks << "/" << a.gap_blocks
                  << " optimized=" << b.bit_blocks << "/" << b.gap_blocks << std::endl;
    assert(a.bit_blocks == b.bit_blocks && a.gap_blocks == b.gap_blocks);
}

static void check(const Bytes& bytes, const BV& expected)
{
    bm::deserializer<BV, bm::decoder> reader;
    BV actual;
    reader.deserialize(actual, bytes.data());
    assert(actual == expected); normalized(actual);
    // Reuse the reader with a nonempty BIT destination: preserve merge behavior.
    BV merged; merged.set(65000);
    auto& manager = merged.get_blocks_manager();
    manager.deoptimize_block(0, 0, true);
    BV want(expected); want.set(65000);
    reader.deserialize(merged, bytes.data());
    assert(merged == want);
    assert(!BM_IS_GAP(manager.get_block_ptr(0, 0)));
    // Reuse it again on an empty destination, including preallocated empty BIT.
    actual.clear(); actual.get_blocks_manager().deoptimize_block(0, 0, true);
    reader.deserialize(actual, bytes.data());
    assert(actual == expected); normalized(actual);
    std::string data(reinterpret_cast<const char*>(bytes.data()), bytes.size());
    std::stringstream stream(data);
    bm::streams_decoder input(stream); bm::streams_deserializer<BV> channel;
    BV streamed; assert(channel.deserialize(streamed, input));
    assert(streamed == expected); normalized(streamed);
    bm::deserialization_index<BV> index; BV dry;
    reader.set_deserialization_index_construct(&index);
    reader.deserialize(dry, bytes.data()); assert(dry.none());
    reader.unset_deserialization_index(); reader.set_deserialization_index_use(&index);
    reader.set_range(100, 100000);
    BV range, wanted; wanted.copy_range(expected, 100, 100000);
    reader.deserialize(range, bytes.data());
    BV clipped; clipped.copy_range(range, 100, 100000); assert(clipped == wanted);
    BV extra(range); extra -= expected; assert(extra.none()); normalized(range);
}

static void block_records()
{
    bm::heap_vector<bm::gap_word_t, BV::allocator_type, true> positions;
    bm::bwt_full_workspace<> workspace;
    for (unsigned pattern = 0; pattern < 5; ++pattern)
    {
        BM_DECLARE_TEMP_BLOCK(block)
        BM_DECLARE_TEMP_BLOCK(transformed)
        bm::bit_block_set(block, 0);
        for (unsigned k = 0; k < 65536; ++k)
        {
            bool one = pattern == 0 ? (k >= 100 && k < 900) :
                       pattern == 1 ? !(k >= 100 && k < 900) :
                       pattern == 2 ? true : pattern == 3 ? false : bool(k & 1);
            if (one) bm::set_bit(block, k);
        }
        BV expected;
        for (unsigned k = 0; k < 65536; ++k)
            if (block[k >> 5] & (1u << (k & 31))) expected.set(k);
        for (unsigned code : {unsigned(bm::set_block_bit), unsigned(bm::set_block_bwt), unsigned(bm::set_block_bwt_unpacked)})
        {
            Bytes bytes(65536);
            bm::encoder enc(bytes.data(), bytes.size());
            enc.put_8(bm::BM_HM_NO_BO | bm::BM_HM_NO_GAPL);
            enc.put_8(static_cast<unsigned char>(code));
            if (code == bm::set_block_bit) enc.put_32(block, bm::set_block_size);
            else
            {
                unsigned short primary = bm::bwt_transform_full(block, transformed, workspace);
                if (code == bm::set_block_bwt_unpacked)
                    bm::bwt_bic_encoder<>::encode(enc, transformed, primary, positions);
                else
                {
                    enc.put_16(8193); enc.put_8(bm::set_block_bit);
                    enc.put_32(transformed, bm::set_block_size);
                    enc.put_8(bm::set_block_bwt_full_post); enc.put_16(primary);
                }
            }
            enc.put_8(bm::set_block_end); bytes.resize(enc.size());
            check(bytes, expected);
        }
    }
    // GAP-XOR record: singleton base and two corrected exception positions.
    Bytes bytes(1024); bm::encoder enc(bytes.data(), bytes.size());
    enc.put_8(bm::BM_HM_NO_BO | bm::BM_HM_NO_GAPL);
    enc.put_8(bm::set_block_gap_xor_framed);
    enc.put_8(bm::set_block_bit_1bit); enc.put_16(5);
    { bm::bit_out<bm::encoder> out(enc); out.gamma(2); out.gamma(3);
      bm::gap_word_t ids[] = {7, 8}; out.bic_encode_u16(ids, 2, 0, 65533); out.flush(); }
    enc.put_8(bm::set_block_end); bytes.resize(enc.size());
    BV expected; expected.set(5); expected.set(7); expected.set(10); check(bytes, expected);
}

static void large_records()
{
    for (bool inverse : {false, true})
    for (unsigned count : {65536u, 131072u})
    {
        std::vector<unsigned> ids(count), scratch;
        BV expected;
        if (inverse) expected.set_range(0, bm::set_sub_total_bits-1);
        for (unsigned k = 0; k < count; ++k)
        {
            ids[k] = 17 + k;
            if (inverse) expected.clear_bit(ids[k]); else expected.set(ids[k]);
        }
        Bytes bytes(count * 4 + 1024); bm::encoder enc(bytes.data(), bytes.size());
        enc.put_8(bm::BM_HM_NO_BO | bm::BM_HM_NO_GAPL | bm::BM_HM_SPARSE);
        enc.put_8(bm::set_sblock_bienc_large);
        { bm::bit_out<bm::encoder> out(enc);
          bm::sblock_large_encoder<decltype(out)>::encode(out, ids.data(), count, 0, 0, scratch, inverse); }
        enc.put_8(bm::set_block_end); bytes.resize(enc.size()); check(bytes, expected);
    }
}

static void wave_records()
{
    std::vector<unsigned> ids(96);
    BV expected;
    for (unsigned k = 0; k < ids.size(); ++k)
    { ids[k] = 17 + k * 1031; expected.set(ids[k]); }
    for (unsigned width : {0u, 128u, 256u, 512u, 1024u, 2048u, 4096u})
    for (unsigned mode = 0; mode < 4; ++mode)
    {
        Bytes bytes(65536); bm::encoder enc(bytes.data(), bytes.size());
        enc.put_8(bm::BM_HM_NO_BO | bm::BM_HM_NO_GAPL | bm::BM_HM_SPARSE);
        enc.put_8(width ? bm::set_sblock_bienc_wave_adaptive : bm::set_sblock_bienc_wave1024);
        bm::sblock_wave_record record; bm::sblock_wave_cost cost;
        { bm::bit_out<bm::encoder> out(enc);
          auto position = [&]() { return bm::id64_t(enc.size()) * 8 + out.pending_bits(); };
          if (width)
          {
              bm::sblock_wave_adaptive_model<> model;
              model.build(ids.data(), unsigned(ids.size()));
              model.wave_size = width; model.wave_count = (1u << 24) / width;
              assert(model.prepare(ids.data(), unsigned(ids.size()), record));
              unsigned code = 0; while ((128u << code) != width) ++code;
              out.put_bits(code, 3); model.encode(out, position, 0, mode, cost);
          }
          else
          {
              bm::sblock_wave1024_model<> model;
              assert(model.prepare(ids.data(), unsigned(ids.size()), record));
              model.encode(out, position, 0, mode, cost);
          }
        }
        enc.put_8(bm::set_block_end); bytes.resize(enc.size()); check(bytes, expected);
    }
    // Exercise legacy 68 and adaptive correction 92 through the actual selector.
    BV source; unsigned at = 10; source.set(at);
    for (unsigned k = 1; k < 96; ++k)
    { at += k >= 32 && k < 64 ? 2048u : 2u; source.set(at); }
    source.optimize();
    for (bool adaptive : {false, true})
    {
        bm::serializer<BV> writer; writer.set_compression_level(6);
        writer.set_sblock_adaptive_options(adaptive);
        bm::serializer<BV>::buffer blob; writer.serialize(source, blob);
        assert(writer.get_compression_stat()[adaptive ? 92 : 68] == 1);
        check(Bytes(blob.data(), blob.data() + blob.size()), source);
    }
}

static void reference_xor()
{
    BV target, reference;
    for (unsigned k = 0; k < 500; ++k) target.set(k * 131);
    reference = target;
    for (unsigned k = 0; k < 20; ++k) reference.flip(k * 197 + 1);
    target.optimize(); reference.optimize();
    using Writer = bm::serializer<BV>;
    Writer::bv_ref_vector_type refs; refs.add(&target, 0); refs.add(&reference, 1);
    Writer::xor_sim_model_type model; Writer writer;
    writer.set_ref_vectors(&refs); writer.compute_sim_model(model, refs, bm::xor_sim_params());
    writer.set_sim_model(&model); writer.set_curr_ref_idx(0);
    Writer::buffer blob; writer.serialize(target, blob);
    bm::deserializer<BV, bm::decoder> reader; reader.set_ref_vectors(&refs);
    BV actual; reader.deserialize(actual, blob.data());
    assert(actual == target); normalized(actual);
    reader.set_range(100, 20000); actual.clear(); reader.deserialize(actual, blob.data());
    BV wanted, clipped; wanted.copy_range(target, 100, 20000); clipped.copy_range(actual, 100, 20000);
    assert(clipped == wanted); normalized(actual);
}

int main()
{
    block_records(); large_records(); wave_records(); reference_xor();
    // The small-vector plain ID-list header bypasses the block dispatch.
    unsigned char bytes[1024] = {}; bm::encoder enc(bytes, sizeof(bytes));
    enc.put_8(bm::BM_HM_ID_LIST | bm::BM_HM_NO_BO); enc.put_32(3);
    BV expected; for (unsigned id : {5u, 100u, 65000u}) { enc.put_32(id); expected.set(id); }
    BV decoded; bm::deserialize(decoded, bytes); assert(decoded == expected); normalized(decoded);
    std::cout << "GAP normalization: sparse/dense/solid/BIT, 68/70/72/90/92/94/96/97/98, reference XOR, empty-only, stream/index/range PASS\n";
}
