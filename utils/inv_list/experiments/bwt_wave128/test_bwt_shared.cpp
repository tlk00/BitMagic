#include "bm.h"
#include "bmserial.h"
#include "bmserial_profile.h"
#include <cassert>
#include <cstring>
#include <fstream>
#include <iostream>
#include <iterator>
#include <random>
int main(int argc, char** argv)
{
    using BV=bm::bvector<>;
    using S=bm::serializer<BV>;
    for (unsigned pattern=0;pattern<4;++pattern) {
        BV bv;
        std::mt19937 rng(219+pattern);
        for(unsigned i=0;i<65536;++i) {
            bool bit=pattern==0 ? i%101<13 : pattern==1 ? i%113<104 :
                     pattern==2 ? rng()%100<24 : rng()%100<78;
            if(bit) bv.set(i);
        }
        S ser;ser.set_block_selection_step2_options();ser.set_bwt_shared_options();
        S::buffer blob;ser.serialize(bv,blob);
        BV restored;bm::deserialize(restored,blob.data());assert(restored==bv);
        bm::serialization_opcode_collector<BV> c;
        bm::serialization_opcode_profile p;c.add_blob(blob.data(),blob.size(),p);
        assert(p.reconciles());
        assert(!p.opcode[95].count);
        if(argc==2) { // Optional pre-cleanup fixtures, produced by the old binary.
            std::ifstream f(std::string(argv[1])+std::to_string(pattern),std::ios::binary);
            assert(f.good());
            std::vector<unsigned char> old((std::istreambuf_iterator<char>(f)),{});
            assert(old.size()==blob.size() && !std::memcmp(old.data(),blob.data(),old.size()));
            BV previous;bm::deserialize(previous,old.data());assert(previous==bv);
        }
    }
    bm::bwt_full_workspace<> decode_workspace;
    // Retired decoder-only opcode-96 polarity/formula variants must fail in
    // both materializing and skip/index traversal, without consuming a tail.
    for (unsigned variant = 0; variant < 4; ++variant) {
        unsigned char bytes[1024] = {}; bm::encoder enc(bytes, sizeof(bytes));
        {
            bm::bit_out<bm::encoder> out(enc);
            out.put_bit(variant == 0); out.put_bits(0, 16);
            out.put_bits(variant == 1, 9); out.put_bit(variant == 2);
            out.gamma(1); // empty residual
            out.gamma(variant == 3 ? 2 : 1); // retired chain formula
            out.flush();
        }
        BM_DECLARE_TEMP_BLOCK(actual)
        for (bool skip : {false, true}) {
            bm::decoder dec(bytes); bool rejected = false;
            try { bm::bwt_bic_decoder<>::decode(dec, skip ? nullptr : actual, decode_workspace); }
            catch (const std::logic_error&) { rejected = true; }
            assert(rejected);
        }
    }
    std::cout<<"Shared 72/96: round trips, profiling, one BWT, retired variant rejection, optional old-byte identity PASS\n";
}
