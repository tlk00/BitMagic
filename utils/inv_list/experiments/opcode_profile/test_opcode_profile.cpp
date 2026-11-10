#include "bm.h"
#include "bmserial.h"
#include "bmserial_profile.h"
#include "bmdbg.h"
#include <cassert>
#include <sstream>
#include <string>

using BV = bm::bvector<>;

int main()
{
    bm::serialization_opcode_collector<BV> collector;
    collector.set_collect_ids();
    bm::serialization_opcode_profile profile;
    bm::serializer<BV> writer;
    writer.set_compression_level(6);
    BV v;
    v.set(5); v.set(65536+19); v.set(4*65536+111);
    v.optimize();
    bm::serializer<BV>::buffer blob;
    writer.serialize(v, blob);
    collector.add_blob(blob.data(), blob.size(), profile);
    assert(profile.reconciles() && profile.blobs == 1 && profile.bytes == blob.size());
    assert(profile.records >= 2 && profile.header_bytes && (profile.opcode[bm::set_block_end].count || profile.opcode[bm::set_block_azero].count));
    assert(profile.cardinalities_reconcile() && profile.ids == v.count());
    auto before = profile;
    bool rejected = false;
    try { collector.add_blob(blob.data(), blob.size()-1, profile); }
    catch (const std::exception&) { rejected = true; }
    assert(rejected && profile.blobs == before.blobs && profile.bytes == before.bytes);
    collector.add_blob(blob.data(), blob.size(), profile);
    assert(profile.reconciles() && profile.bytes == 2*blob.size());

    // Packed zero-run is a marker byte with no separate payload.
    const unsigned char packed[] = {bm::BM_HM_NO_BO | bm::BM_HM_NO_GAPL, 0x83, bm::set_block_end};
    collector.add_blob(packed, sizeof(packed), profile);
    assert(profile.opcode[0x83].count == 1 && profile.opcode[0x83].bytes == 1);
    assert(profile.reconciles());

    // Legacy plain ID-list format has no block markers.
    const unsigned char ids[] = {bm::BM_HM_NO_BO | bm::BM_HM_NO_GAPL | bm::BM_HM_ID_LIST, 0, 0, 0, 0};
    collector.add_blob(ids, sizeof(ids), profile);
    assert(profile.id_list_bytes == sizeof(ids) && profile.reconciles());
    assert(profile.cardinalities_reconcile() && profile.ids == 2*v.count());
    // Multi-block ones, holes, sparse superblocks, explicit size and bookmarks.
    for (unsigned mode = 0; mode < 5; ++mode) {
        BV mixed;
        mixed.set_range(65536, 4*65536-1);
        mixed.set(2*65536+17, false);
        for (unsigned i=0; i<5000; ++i) mixed.set(8*65536+i*503);
        mixed.resize(8*1024*1024);
        mixed.optimize();
        writer.set_bookmarks(true, 4);
        writer.set_sblock_adaptive_options(mode == 1);
        writer.set_sblock_wave1024_options(mode == 2);
        writer.set_sblock_wave_adaptive_options(mode == 3);
        writer.set_gap_xor_options(mode == 4);
        writer.serialize(mixed, blob);
        bm::serialization_opcode_profile p;
        collector.add_blob(blob.data(),blob.size(),p);
        assert(p.reconciles() && p.cardinalities_reconcile() && p.ids == mixed.count());
        // Enabling profiling never changes serialization bytes.
        bm::serializer<BV>::buffer again;
        writer.serialize(mixed, again);
        assert(blob.size() == again.size() && !memcmp(blob.data(),again.data(),blob.size()));
    }
    // BWT/chain records count restored positions, not transformed residual IDs.
    writer.set_block_selection_step2_options(); writer.set_bwt_shared_options();
    for (unsigned inverse=0; inverse<2; ++inverse) {
        BV bwt;
        for (unsigned i=0;i<65536;++i)
            if ((i%101<13) != bool(inverse)) bwt.set(i);
        writer.serialize(bwt,blob);
        bm::serialization_opcode_profile bp;
        collector.add_blob(blob.data(),blob.size(),bp);
        assert(bp.cardinalities_reconcile() && bp.ids == bwt.count());
    }
    // Direct fixtures give known per-opcode denominators.
    const unsigned char known[] = {bm::BM_HM_NO_BO | bm::BM_HM_NO_GAPL,
        bm::set_block_1one, bm::set_block_bit_1bit, 7, 0,
        bm::set_block_1zero, bm::set_block_end};
    bm::serialization_opcode_profile p;
    collector.add_blob(known,sizeof(known),p);
    assert(p.cardinalities_reconcile() && p.ids == 65537);
    assert(p.opcode[bm::set_block_1one].ids == 65536);
    assert(p.opcode[bm::set_block_bit_1bit].ids == 1);
    assert(p.opcode[bm::set_block_1zero].ids == 0);
    const unsigned char plain_ids[] = {
        bm::BM_HM_NO_BO | bm::BM_HM_NO_GAPL | bm::BM_HM_ID_LIST,
        2,0,0,0, 3,0,0,0, 7,0,0,0};
    bm::serialization_opcode_profile plain;
    collector.add_blob(plain_ids,sizeof(plain_ids),plain);
    assert(plain.cardinalities_reconcile() && plain.ids == 2 && plain.id_list_ids == 2);
    std::ostringstream out;
    bm::print_serialization_opcode_profile(out, profile);
    assert(out.str().find("reconciled=1") != std::string::npos);
    assert(out.str().find("packed_zero_run") != std::string::npos);
}
