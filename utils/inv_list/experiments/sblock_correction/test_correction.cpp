#include "bm.h"
#include "bmserial.h"
#include "bmserial_profile.h"
#include <cassert>
#include <cstring>
#include <iostream>
#include <random>
using BV=bm::bvector<>;
struct writer : bm::serializer<BV> {
    std::vector<unsigned char> record(const BV& bv,int mode) {
        std::vector<unsigned char> out(1024*1024);
        bm::encoder enc(out.data(),out.size());
        bienc_arr_sblock(bv,0,enc,mode);
        out.resize(enc.size());return out;
    }
};
int main() {
    std::mt19937 rng(68);
    unsigned wins=0,ties=0;
    for (unsigned k=0;k<240;++k) {
        BV bv; unsigned n=3+rng()%600, pos=rng()%200000;
        unsigned gap=k<8 ? (1u<<k) : 1+rng()%500;
        for(unsigned i=0;i<n;++i) { bv.set(pos);pos+=gap+(k%3 ? rng()%2000 : 0); }
        bv.optimize();
        writer s; s.set_sblock_correction_selection(true);
        auto old=s.record(bv,2), zero=s.record(bv,0), corrected=s.record(bv,1), selected=s.record(bv,-1);
        const auto& expected=zero.size()<corrected.size()?zero:corrected;
        assert(selected.size()==expected.size());
        if(selected.size()==old.size()) { assert(selected==old);++ties; } else ++wins;
        std::vector<unsigned char> blob={bm::BM_HM_NO_BO|bm::BM_HM_NO_GAPL};
        blob.insert(blob.end(),selected.begin(),selected.end());blob.push_back(bm::set_block_end);
        BV decoded;bm::deserialize(decoded,blob.data());assert(decoded==bv);
        bm::serialization_opcode_collector<BV> direct;direct.set_collect_ids();
        bm::serialization_opcode_profile dp;direct.add_blob(blob.data(),blob.size(),dp);
        assert(dp.ids==bv.count());
        // Real serializer selection, including the retained rival codecs.
        for(unsigned rivals=0;rivals<2;++rivals) {
            bm::serializer<BV> a,b;
            a.set_compression_level(6);b.set_compression_level(6);
            a.set_sblock_adaptive_options(rivals);b.set_sblock_adaptive_options(rivals);
            a.set_sblock_wave_adaptive_options(rivals);b.set_sblock_wave_adaptive_options(rivals);
            b.set_sblock_correction_selection(true);
            bm::serializer<BV>::buffer ab,bb;
            a.serialize(bv,ab);b.serialize(bv,bb);assert(bb.size()<=ab.size());
            BV restored;bm::deserialize(restored,bb.data());assert(restored==bv);
            bm::serialization_opcode_collector<BV> collector;collector.set_collect_ids();
            bm::serialization_opcode_profile p;collector.add_blob(bb.data(),bb.size(),p);
            assert(p.cardinalities_reconcile() && p.ids==bv.count());
            for(unsigned c: {68u,90u,92u,97u}) assert(p.opcode[c].count==b.get_compression_stat()[c]);
        }
    }
    assert(wins && ties);
    std::cout<<"Correction min-cost, legacy decoding, ties, rival selection and profiling PASS; wins="<<wins<<" ties="<<ties<<'\n';
}
