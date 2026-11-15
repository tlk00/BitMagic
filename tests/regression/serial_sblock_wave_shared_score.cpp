/** Compare consolidated exact choices against the original independent scorers, including ties and count limits. */
#include "bm.h"
#include "bmserial.h"
#include <vector>
#include <iostream>
using Alloc=bm::standard_allocator;
int main()
{
 unsigned rng=17;
 for(unsigned trial=0;trial<184;++trial)
 {
  std::vector<unsigned> ids;
  unsigned pos=0;
  for(unsigned i=0;i<100+trial*7;++i)
  {
   rng=rng*1664525u+1013904223u;
   pos+=1+rng%(trial%3==0?20000:trial%3==1?23:400);
   if(pos>=(1u<<24))break;
   ids.push_back(pos);
  }
  if (trial>=180)
  {
      ids.clear();
      unsigned count=trial==180?2:trial==181?65535:trial==182?5000:256;
      for(unsigned i=0;i<count;++i)
          ids.push_back(trial==180?i*6553600u:unsigned((bm::id64_t(i)*(1u<<24))/count));
  }
  for(unsigned policy:{1u,2u,3u})
  {
   unsigned occupied=0, prior=~0u;
   for(unsigned id:ids)
   {
       if((id>>16)!=prior) ++occupied;
       prior=id>>16;
   }
   const unsigned expected_policy=policy==3?(occupied<=16?2u:1u):policy;
   bm::sblock_wave_adaptive_model<Alloc> oldflat,newflat;
   bm::sblock_hier_model<Alloc> oldhier,newhier;
   bm::id64_t rival=trial%5==0?20:1000000;
   bm::sblock_wave_adaptive_choice expected;expected.bytes=rival;
   bool expected_h=false;
   if(expected_policy!=2)expected=bm::sblock_wave_adaptive_select(oldflat,ids.data(),unsigned(ids.size()),trial,rival);
   oldhier.build(ids.data(),unsigned(ids.size()));
   for(unsigned w=0;w<6;++w)
   {
    oldhier.wave.wave_size=128u<<w;
    oldhier.wave.wave_count=unsigned(oldhier.blocks.size())*65536/oldhier.wave.wave_size;
    bm::sblock_wave_record r;
    oldhier.wave.prepare(oldhier.source.data(),unsigned(ids.size()),r,true);
    oldhier.pairs();
    for(unsigned mode=0;mode<4;++mode)
    {
     if(mode==3&&!oldhier.entries)continue;
     bm::sblock_wave_bit_counter counter;
     oldhier.encode(counter,trial,w,mode);
     bm::id64_t bytes=1+4*((counter.bits+31)/32);
     if(bytes<expected.bytes){expected.bytes=bytes;expected.width=w;expected.mode=mode;expected_h=true;}
    }
   }
   bool actual_h;
   auto actual=bm::sblock_wave_shared_select(newflat,newhier,ids.data(),unsigned(ids.size()),trial,rival,policy,actual_h,nullptr,nullptr);
   if(actual.bytes!=expected.bytes||actual_h!=expected_h||
      (actual.bytes<rival&&(actual.width!=expected.width||actual.mode!=expected.mode)))
   {std::cerr<<"mismatch "<<trial<<' '<<policy<<' '<<expected.bytes<<' '<<actual.bytes<<'\n';return 1;}
  }
 }
 std::cout<<"552 exact selector decisions match old scorers\n";
}
