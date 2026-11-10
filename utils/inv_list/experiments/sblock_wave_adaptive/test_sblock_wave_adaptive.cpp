#include "bm.h"
#include "bmserial.h"
#include "bmserial_profile.h"
#include "test_model.h"
#include <cassert>
#include <random>
#include <iostream>

// The pre-existing exploration model is deliberately independent of the new
// occupancy merging/counting selector: exhaustively materialize all layouts.
static void check(std::vector<unsigned> positions,unsigned sb)
{
    std::sort(positions.begin(),positions.end());
    positions.erase(std::unique(positions.begin(),positions.end()),positions.end());
    bm::sblock_wave_adaptive_model<> model;
    const bm::id64_t incumbent=10000000;
    auto selected=bm::sblock_wave_adaptive_select(model,positions.data(),unsigned(positions.size()),sb,incumbent);
    control_wave_test::trial_model control;
    bm::id64_t best=incumbent;unsigned bw=0,bm_=0;
    std::vector<unsigned char> bytes(2000000);
    for(unsigned w=0;w<6;++w) {
        control.wave_size=128u<<w;control.wave_count=(1u<<24)/control.wave_size;
        bm::sblock_wave_record r;if(!control.prepare(positions.data(),unsigned(positions.size()),r))continue;
        for(unsigned mode=0;mode<4;++mode) {
            bm::encoder enc(bytes.data(),bytes.size());enc.put_8(97);
            {bm::bit_out<bm::encoder> out(enc);out.put_bits(w,3);bm::sblock_wave_cost c;
             auto pos=[&](){return bm::id64_t(enc.size())*8+out.pending_bits();};
             control.encode(out,pos,sb,mode,c);}
            if(enc.size()<best) {best=enc.size();bw=w;bm_=mode;}
            // Every width/layout round trips through the new production reader.
            bm::decoder dec(bytes.data()+1);bm::bit_in<bm::decoder> in(dec);
            std::vector<unsigned> restored(positions.size());unsigned decoded_sb;
            unsigned n=bm::sblock_wave_adaptive_model<>::decode(in,restored.data(),decoded_sb);
            assert(n==positions.size() && restored==positions && decoded_sb==sb);
            assert(size_t(dec.get_pos()-bytes.data())==enc.size());
        }
    }
    assert(selected.bytes==best && selected.width==bw && selected.mode==bm_);
    auto tie=bm::sblock_wave_adaptive_select(model,positions.data(),unsigned(positions.size()),sb,best);
    assert(tie.bytes==best); // strict smaller, no emitted candidate on ties
}
int main()
{
    check({0,1},0);check({0,16777215},7);check({127,128,4095,4096},255);
    std::vector<unsigned> positions;
    for(unsigned i=0;i<65535;++i)positions.push_back(i);check(positions,65535);
    positions.clear();for(unsigned i=0;i<131072;i+=2)positions.push_back(i*128+3);
    positions.pop_back();check(positions,70000); // 32-bit map IDs > 65535
    std::mt19937 random(34129);
    for(unsigned trial=0;trial<30;++trial) {
        positions.clear();for(unsigned i=0;i<2+trial*131;++i)positions.push_back(random()%(1u<<24));
        check(positions,trial);
    }
    bm::bvector<> bv;
    for(unsigned i=0;i<100;++i)bv.set(1+i*163841);bv.optimize();
    bm::serializer<bm::bvector<>> ser;ser.set_compression_level(6);
    ser.set_sblock_wave_adaptive_options(true);
    bm::serializer<bm::bvector<>>::buffer blob;ser.serialize(bv,blob);
    bm::bvector<> restored;bm::deserialize(restored,blob.data());assert(restored==bv);
    bm::serialization_opcode_collector<bm::bvector<>> collector;bm::serialization_opcode_profile profile;
    collector.add_blob(blob.data(),blob.size(),profile);
    assert(profile.opcode[97].count>0 && profile.reconciles());
    // Invalid selector codes are rejected before using them as shifts.
    unsigned char invalid[4]={6,0,0,0};bm::decoder dec(invalid);bm::bit_in<bm::decoder> in(dec);
    bool rejected=false;unsigned sb;unsigned output[2];
    try {bm::sblock_wave_adaptive_model<>::decode(in,output,sb);}catch(const std::logic_error&){rejected=true;}
    assert(rejected);
    std::cout<<"Adaptive waves: exhaustive score/choice parity, all-width/layout round trips, collector, invalid width and tie fallback PASS\n";
}
