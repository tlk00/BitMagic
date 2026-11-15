/** Experimental opcode-109 transform round trips for all six widths and four payloads. */
#include "bm.h"
#include "bmserial.h"
#include "bmserial_profile.h"
#include <vector>
#include <iostream>
int main()
{
 bm::sblock_hier_model<bm::standard_allocator> h;
 std::vector<unsigned> ids;
 // More than eight distinct pairs exercises dictionary misses as literals.
 unsigned state = 17, position = 0;
 for (unsigned group = 0; group < 2500; ++group)
 {
     state = state * 1664525u + 1013904223u;
     const unsigned run = 1 + state % 7;
     for (unsigned j = 0; j < run; ++j)
     {
         ids.push_back(position++);
     } // for j
     state = state * 1664525u + 1013904223u;
     position += 1 + state % 23;
     if (group % 43 == 0)
     {
         position += 16384;
     }
 } // for group
 h.build(ids.data(),unsigned(ids.size()));
 std::vector<unsigned char> buf(1000000);
 std::vector<unsigned> dst(ids.size());
 for (unsigned superblock : {5u,255u,256u,65536u,0x123456u})
 for(unsigned w=0;w<6;++w)for(unsigned m=0;m<4;++m)
 {
 h.wave.wave_size=128u<<w;h.wave.wave_count=unsigned(h.blocks.size())*65536/h.wave.wave_size;
 bm::sblock_wave_record r;h.wave.prepare(h.source.data(),unsigned(ids.size()),r,true);h.pairs();
 bm::encoder enc(buf.data(),buf.size());
 {bm::bit_out<bm::encoder> out(enc);h.encode(out,superblock,w,m);}
 bm::decoder dec(buf.data());bm::bit_in<bm::decoder> in(dec);
 bm::heap_vector<unsigned,bm::standard_allocator,true> blocks,waves,full;
 unsigned sb;unsigned n=h.decode(in,dst.data(),sb,blocks,waves,full);
 if(!full.empty())return 2;
 if(n!=ids.size()||sb!=superblock||dst!=ids)return 1;
 }
 std::cout<<"120 header/width/payload round trips passed\n";
}
