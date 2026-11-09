// Focused preset regression: build with -std=c++17 -Isrc.
#include "bm.h"
#include "bmserial.h"
#include <cassert>
#include <cstring>
#include <iostream>
using BV = bm::bvector<>;
using Serializer = bm::serializer<BV>;

static void explicit_options(Serializer& s, unsigned level)
{
    s.set_compression_level(6);
    if (level >= 7) {
        s.set_block_endpoint_selection(true);
        s.set_sblock_correction_selection(true);
    }
    if (level >= 8) {
        s.set_gap_xor_options(true);
        s.set_sblock_adaptive_options(true);
        s.set_sblock_wave_adaptive_options(true);
        s.set_sblock_large_options(true);
    }
    if (level >= 9) {
        s.set_block_selection_step2_options();
        s.set_periodic_step1_interleaved(4);
        s.set_bwt_shared_options();
        s.set_block_selection_bwt_min_gain(level == 9 ? 16 : 1);
    }
}
static Serializer::buffer encode(Serializer& s, const BV& bv)
{
    Serializer::buffer b; s.serialize(bv, b);
    BV restored; bm::deserialize(restored, b.data()); assert(restored == bv);
    return b;
}
static void same(const Serializer::buffer& a, const Serializer::buffer& b)
{ assert(a.size() == b.size()); assert(!std::memcmp(a.data(), b.data(), a.size())); }
int main()
{
    assert(bm::set_compression_max == 10 && bm::set_compression_default == 6);
    for (unsigned pattern = 0; pattern < 6; ++pattern) {
        BV bv;
        if (pattern == 0) bv.set(17);
        if (pattern == 1) for (unsigned k=0; k<65536; ++k) if (k%101<13) bv.set(k);
        if (pattern == 2) for (unsigned k=0; k<65536; ++k) if (k%113<104) bv.set(k);
        if (pattern == 3) for (unsigned k=0; k<65536; ++k) bv.set(17+k*67);
        if (pattern == 4) for (unsigned k=0; k<131072; ++k) bv.set(17+k*3);
        if (pattern == 5) for (unsigned k=0; k<65536; ++k) if (k%19!=0) bv.set(k);
        bv.optimize();
        size_t previous = 0;
        for (unsigned level=6; level<=10; ++level) {
            Serializer preset, manual;
            preset.set_compression_level(level); explicit_options(manual,level);
            auto bytes=encode(preset,bv); same(bytes,encode(manual,bv));
            if (level>6) assert(bytes.size()<=previous);
            previous=bytes.size();
            preset.set_compression_level(10);
            preset.set_block_selection_bwt_min_gain(234);
            preset.set_sblock_large_inverse_limit(300000);
            preset.set_compression_level(level);
            same(bytes,encode(preset,bv));
            preset.set_compression_level(11);
            assert(preset.get_compression_level()==level);
            if (level>=9) assert(preset.get_block_selection_statistics().bwt_min_gain==(level==9?16:1));
        }
        Serializer reused;
        for (unsigned level=0; level<=6; ++level) {
            reused.set_compression_level(10); reused.set_compression_level(level);
            Serializer fresh; fresh.set_compression_level(level);
            same(encode(fresh,bv),encode(reused,bv));
        }
    }
    std::cout << "Levels 0..10: explicit-policy identity, round trips, preset reset, invalid levels, monotonic sizes PASS\n";
}
