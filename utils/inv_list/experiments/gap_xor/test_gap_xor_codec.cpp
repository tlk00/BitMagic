#include "bm.h"
#include "bmserial.h"
#include "bmfio.h"
#include <cassert>
#include <iostream>
#include <sstream>
#include <thread>

struct reader : bm::deseriaizer_base<bm::decoder, bm::bvector<>::block_idx_type>
{
    bm::word_t ex0[bm::set_block_size*4] = {}, ex1[bm::set_block_size*4] = {};
    bm::gap_word_t ids[65536] = {};
    reader() { ex0_arr_=ex0; ex1_arr_=ex1; id_array_=ids; }
    void decode(const unsigned char* data, bm::gap_word_t* output, unsigned type=bm::set_block_gap_xor_framed) {
        bm::decoder dec(data); bm::gap_word_t head=0;
        read_gap_block(dec,type,output,head);
    }
};

static void roundtrip(unsigned worker)
{
    for(unsigned n : {32u,63u,96u,160u,320u,500u}) {
        bm::bvector<> bv;
        const unsigned offset=worker*65536u;
        bv.set_range(offset+200,offset+65000);
        unsigned pos=1000;
        for(unsigned i=0;i<n;++i) { pos += (i/32)%2 ? 90u : 3u; bv.clear_bit(offset+pos); }
        bv.optimize();
        bm::serializer<bm::bvector<> > ordinary, framed;
        ordinary.set_compression_level(6);
        framed.set_compression_level(6); framed.set_gap_xor_options(true);
        bm::serializer<bm::bvector<> >::buffer before, middle;
        ordinary.serialize(bv,before); framed.serialize(bv,middle);
        assert(middle.size()<=before.size());
        bm::bvector<> framed_decoded; bm::deserialize(framed_decoded,middle.data());
        assert(framed_decoded==bv);
        bm::deserializer<bm::bvector<>, bm::decoder> partial;
        bm::deserializer<bm::bvector<>, bm::decoder>::deserialization_index_type index;
        bm::bvector<> indexed;
        partial.set_deserialization_index_construct(&index);
        assert(partial.deserialize(indexed,middle.data()) == middle.size());
        partial.unset_deserialization_index();
        partial.set_deserialization_index_use(&index);
        bm::bvector<> requested, skipped;
        requested.set(worker+1); partial.set_block_digest_vector_use(&requested);
        partial.deserialize(skipped,middle.data()); assert(skipped.none());
        bm::operation_deserializer<bm::bvector<> > op;
        bm::bvector<> combined; op.deserialize(combined,middle.data(),bm::set_OR); assert(combined==bv);
        std::stringstream stream(std::ios::in|std::ios::out|std::ios::binary);
        bm::streams_encoder output(stream); assert(framed.serialize_and_finish(bv,output));
        stream.seekg(0); bm::streams_decoder input(stream); bm::streams_deserializer<bm::bvector<> > sr;
        bm::bvector<> streamed; assert(sr.deserialize(streamed,input)); assert(streamed==bv);
    }
}
int main() {
    assert(!bm::superblock_delta_model::eligible(32,16));
    {
        unsigned char data[128]={}; bm::encoder enc(data,sizeof(data));
        enc.put_8(bm::set_block_bit_1bit); enc.put_16(5);
        { bm::bit_out<bm::encoder> out(enc);
          out.gamma(2); out.gamma(3);
          bm::gap_word_t corrected[]={7,8};
          out.bic_encode_u16(corrected,2,0,65533); out.flush(); }
        reader r; bm::gap_word_t output[8192+8]={};
        r.decode(data,output,bm::set_block_gap_xor_framed);
        assert(bm::gap_test(output,5) && bm::gap_test(output,7) && bm::gap_test(output,10));
    }
    roundtrip(0);
    std::vector<std::thread> workers;
    for(unsigned i=0;i<10;++i) workers.emplace_back([i]{roundtrip(i);});
    for(auto& t:workers)t.join();
    std::cout<<"GAP-XOR BIC codec, streaming, operations and 10 workers PASS\n";
}
