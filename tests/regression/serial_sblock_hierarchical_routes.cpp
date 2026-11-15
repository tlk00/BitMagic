/** Experimental 109 ordinary, OR, range and indexed skip checks across superblocks. */
#include "bm.h"
#include "bmserial.h"
#include "bmserial_profile.h"
#include <iostream>
int main()
{
 bm::bvector<> a;
 for(unsigned base:{0u,2u*(1u<<24)})for(unsigned b:{0u,3u,97u})for(unsigned i=1;i<6000;i+=7)a.set(base+b*65536+i);
 for(unsigned policy:{0u,1u,2u})
 {
 bm::serializer<bm::bvector<>> s;s.set_compression_level(10);s.set_sblock_wave_experiment(policy);
 bm::serializer<bm::bvector<>>::buffer buf;s.serialize(a,buf);
 bm::bvector<> b;bm::deserialize(b,buf.data());if(b!=a)return 1;
 bm::serialization_opcode_collector<bm::bvector<>> col;bm::serialization_opcode_profile p;
 unsigned seen=0;col.add_blob(buf.data(),buf.size(),p,[&](const unsigned char* ptr,size_t){seen+=*ptr==109;});
 if(policy==2 && !seen)return 2;
 if(!p.reconciles())return 3;
 bm::deserializer<bm::bvector<>,bm::decoder> indexed;
 bm::deserialization_index<bm::bvector<>> index;
 bm::bvector<> dry;
 indexed.set_deserialization_index_construct(&index);
 indexed.deserialize(dry,buf.data());
 indexed.unset_deserialization_index();
 indexed.set_deserialization_index_use(&index);
 indexed.set_range(2u*(1u<<24)+65536,2u*(1u<<24)+100*65536);
 bm::bvector<> indexed_result;
 indexed.deserialize(indexed_result,buf.data());
 bm::bvector<> wanted,clipped;
 wanted.copy_range(a,2u*(1u<<24)+65536,2u*(1u<<24)+100*65536);
 clipped.copy_range(indexed_result,2u*(1u<<24)+65536,2u*(1u<<24)+100*65536);
 if(clipped!=wanted)return 6;
 bm::bvector<> unexpected(indexed_result);unexpected-=a;if(unexpected.any())return 7;
 bm::bvector<> slice;
 bm::operation_deserializer<bm::bvector<>> reader;
 reader.deserialize_range(slice,buf.data(),65536,100*65536);
 bm::bvector<> expected(a);expected.keep_range(65536,100*65536);
 if(slice!=expected)return 5;
 bm::bvector<> merged;merged.set(7);bm::deserialize(merged,buf.data());b.set(7);if(merged!=b)return 4;
 std::cout<<policy<<" bytes="<<buf.size()<<" opcode109="<<seen<<"\n";
 }
}
