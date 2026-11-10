#include "bm.h"
#include "bmserial.h"
#include "bmserial_profile.h"
#include <cassert>
#include <iostream>
#include <pthread.h>
#include <sstream>
using BV=bm::bvector<>;
static void* test(void*) {
    unsigned wins=0;
    for(unsigned n:{65535u,65536u,131071u,131072u,196607u,196608u})
    for(unsigned sb:{0u,1u,8u,255u}) {
        BV bv;unsigned base=sb*bm::set_sub_total_bits;
        for(unsigned i=0;i<n;++i)bv.set(base+17+i*67);
        bv.optimize();bm::serializer<BV> ordinary,trial;
        ordinary.set_compression_level(6);trial.set_compression_level(6);trial.set_sblock_large_options(true);
        bm::serializer<BV>::buffer a,b,c;ordinary.serialize(bv,a);trial.serialize(bv,b);
        assert(b.size()<=a.size());BV restored;bm::deserialize(restored,b.data());assert(restored==bv);
        bm::serialization_opcode_collector<BV> collector;collector.set_collect_ids();bm::serialization_opcode_profile p;
        collector.add_blob(b.data(),b.size(),p);assert(p.reconciles()&&p.cardinalities_reconcile()&&p.ids==n);
        bool eligible=n>=65536&&n<=196607;
        assert(bool(p.opcode[98].count)==eligible);wins+=eligible;
        if(!eligible)assert(a.size()==b.size()&&!memcmp(a.data(),b.data(),a.size()));
        trial.set_sblock_large_options(false);trial.serialize(bv,c);
        assert(c.size()==a.size()&&!memcmp(c.data(),a.data(),a.size()));
        // Full-stream dry index and bounded streaming decode are exercised above.
        if(eligible){
            trial.set_sblock_large_options(true);
            std::stringstream stream;bm::streams_encoder channel(stream);
            assert(trial.serialize(bv,channel)&&channel.finish());
            assert(stream.str().size()==b.size()&&!memcmp(stream.str().data(),b.data(),b.size()));
            bm::deserialization_index<BV> index;BV dry;
            bm::deserializer<BV,bm::decoder> reader;
            reader.set_deserialization_index_construct(&index);reader.deserialize(dry,b.data());
            reader.unset_deserialization_index();reader.set_deserialization_index_use(&index);
            const unsigned from=base+1000,to=base+20000;
            reader.set_range(from,to);BV range;reader.deserialize(range,b.data());
            assert(range.count_range(from,to)==bv.count_range(from,to));
            for(auto e=range.first();e.valid();++e)assert(bv.test(*e));
            for(size_t cut:{size_t(1),b.size()/2,b.size()-1}) {
                bool rejected=false;try{bm::serialization_opcode_profile q;collector.add_blob(b.data(),cut,q);}catch(...){rejected=true;}
                assert(rejected);
            }
        }
    }
    // Inverse wire counts include empty/singleton/two-element lists and the
    // old positive-count boundaries. Test the codec without selection bias.
    for(unsigned n:{0u,1u,2u,3u,65535u,65536u,196607u,262144u}) {
        std::vector<unsigned> ids(n),scratch,decoded;
        for(unsigned k=0;k<n;++k)ids[k]=17+k*31;
        std::vector<unsigned char> bytes(size_t(n)*4+1024);
        bm::encoder enc(bytes.data(),bytes.size());
        {bm::bit_out<bm::encoder> bits(enc);bm::sblock_large_encoder<decltype(bits)>::encode(bits,ids.data(),n,8,0,scratch,true);}
        {bm::decoder probe(bytes.data());bm::bit_in<bm::decoder> header(probe);
         assert(header.get_bits(1)==1);unsigned flags=header.get_bits(8);
         assert(bool(flags&bm::sblock_large_format::range3)==(n>=65536));
         assert(((flags&bm::sblock_large_format::range3)?header.get_24_no():header.get_16_no())==n);}
        bm::decoder dec(bytes.data());unsigned sb=0;bool inverse=false;
        bm::bit_in<bm::decoder> bits(dec);
        unsigned count=bm::sblock_large_decoder<decltype(bits)>::decode(bits,decoded,sb,inverse);
        assert(inverse&&count==n&&sb==8&&std::equal(ids.begin(),ids.end(),decoded.begin()));
    }
    // Dense groups, mixed small inverse / legacy / positive records, last
    // address-space group, indexed ranges, streaming, OR restoration and caps.
    for(unsigned zero_count:{1024u,65536u})
    for(unsigned group:{0u,1u,255u}) {
        BV bv;BV::size_type first=BV::size_type(group)*bm::set_sub_total_bits;
        BV::size_type last=std::min(first+bm::set_sub_total_bits-1,BV::size_type(bm::id_max-1));
        bv.set_range(first,last);
        for(unsigned k=0;k<zero_count;++k)bv.clear_bit(first+17+k*((bm::set_sub_total_bits-64)/zero_count));
        if(group==0){for(unsigned k=0;k<65536;++k)bv.set(2u*bm::set_sub_total_bits+k*67);bv.set(3u*bm::set_sub_total_bits+123);}
        bv.optimize();bm::serializer<BV> s;s.set_compression_level(6);s.set_sblock_large_options(true);
        bm::serializer<BV>::buffer blob;s.serialize(bv,blob);
        assert(s.get_compression_stat()[98]>=1);
        BV restored;bm::deserialize(restored,blob.data());assert(restored==bv);
        // Deserialization combines with an existing vector; zero exceptions
        // must not clear bits already present in that vector.
        restored.clear();restored.set(first+17);BV expected=bv;expected.set(first+17);
        bm::deserialize(restored,blob.data());assert(restored==expected);
        bm::serialization_opcode_collector<BV> collector;collector.set_collect_ids();bm::serialization_opcode_profile profile;
        collector.add_blob(blob.data(),blob.size(),profile);assert(profile.ids==bv.count()&&profile.cardinalities_reconcile());
        bm::deserialization_index<BV> index;bm::deserializer<BV,bm::decoder> reader;BV dry;
        reader.set_deserialization_index_construct(&index);reader.deserialize(dry,blob.data());
        reader.unset_deserialization_index();reader.set_deserialization_index_use(&index);reader.set_range(first+10,first+100000);
        BV range;reader.deserialize(range,blob.data());BV want;want.copy_range(bv,first+10,first+100000);assert(range==want);
        std::stringstream stream;bm::streams_encoder channel(stream);assert(s.serialize(bv,channel)&&channel.finish());
        assert(stream.str().size()==blob.size()&&!memcmp(stream.str().data(),blob.data(),blob.size()));
        for(size_t cut:{size_t(1),blob.size()/2,blob.size()-1}){bool rejected=false;
            try{bm::serialization_opcode_profile q;collector.add_blob(blob.data(),cut,q);}catch(...){rejected=true;}assert(rejected);}
        s.set_sblock_large_inverse_limit(512);s.serialize(bv,blob);
        BV capped;bm::deserialize(capped,blob.data());assert(capped==bv);
    }
    // Large record, legacy sparse group, and a zero run crossing a boundary.
    {
        BV bv;
        for(unsigned i=0;i<65536;++i)bv.set(i*67);
        for(unsigned i=0;i<1000;++i)bv.set(2*bm::set_sub_total_bits+i*73);
        for(unsigned i=0;i<65535;++i)bv.set(255u*bm::set_sub_total_bits+i*67);
        bv.optimize();bm::serializer<BV> s;s.set_compression_level(6);s.set_sblock_large_options(true);
        bm::serializer<BV>::buffer blob;s.serialize(bv,blob);BV decoded;bm::deserialize(decoded,blob.data());assert(decoded==bv);
        bm::serialization_opcode_collector<BV> collector;collector.set_collect_ids();bm::serialization_opcode_profile p;
        collector.add_blob(blob.data(),blob.size(),p);assert(p.ids==bv.count()&&p.cardinalities_reconcile());
    }
    // Dense runs must retain the exact ordinary bytes, including cross-group runs.
#ifdef BM64ADDR
    {
        BV bv;BV::size_type base=BV::size_type(65536)*bm::set_sub_total_bits;
        for(unsigned i=0;i<196607;++i)bv.set(base+17+i*67);
        bv.optimize();bm::serializer<BV> s;s.set_compression_level(6);s.set_sblock_large_options(true);
        bm::serializer<BV>::buffer blob;s.serialize(bv,blob);BV decoded;bm::deserialize(decoded,blob.data());assert(decoded==bv);
        bm::serialization_opcode_collector<BV> collector;collector.set_collect_ids();bm::serialization_opcode_profile p;
        collector.add_blob(blob.data(),blob.size(),p);assert(p.opcode[98].count==1&&p.ids==bv.count());
    }
#endif
    for(unsigned start:{0u,bm::set_sub_total_bits-65536u}) {
        BV bv;bv.set_range(start,start+131071);bv.optimize();bm::serializer<BV> s;
        bm::serializer<BV>::buffer a,b;s.set_compression_level(6);s.serialize(bv,a);
        s.set_sblock_large_options(true);s.serialize(bv,b);
        assert(a.size()==b.size()&&!memcmp(a.data(),b.data(),a.size()));
    }
    std::cout<<"Large count boundary/fallback/profile/truncation PASS; wins="<<wins<<'\n';return nullptr;
}
int main(int argc,char**) {
    if(argc>1){pthread_attr_t attr;pthread_attr_init(&attr);pthread_attr_setstacksize(&attr,256*1024);
        pthread_t thread;assert(!pthread_create(&thread,&attr,test,nullptr));pthread_attr_destroy(&attr);pthread_join(thread,nullptr);}
    else test(nullptr);
}
