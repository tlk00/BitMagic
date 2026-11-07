// Build once against a released src/ with -DBM_WRITE_RELEASE_FIXTURES, then
// against current src/ without that define. Both binaries use the same fixtures.
#include "bm.h"
#include "bmserial.h"
#ifndef BM_WRITE_RELEASE_FIXTURES
#include "bmfio.h"
#include <sstream>
#endif
#include <cassert>
#include <fstream>
#include <iostream>
#include <iterator>
#include <string>
#include <vector>

using BV = bm::bvector<>;
using Writer = bm::serializer<BV>;

static BV fixture(unsigned pattern)
{
    BV bv;
    switch (pattern)
    {
    case 0: break;
    case 1: bv.set(17); break;
    case 2: for (unsigned i=0;i<500;++i) bv.set(i*131); break;
    case 3:
        bv.set_range(200,65000);
        for (unsigned i=0;i<320;++i) bv.clear_bit(1000+i*3);
        break;
    case 4: for (unsigned i=0;i<65536;++i) if (i%97<37) bv.set(i); break;
    case 5:
        { unsigned x=219;
          for (unsigned i=0;i<65536;++i) {x^=x<<13;x^=x>>17;x^=x<<5;if(x%100<41)bv.set(i);} }
        break;
    case 6: for (unsigned i=0;i<6000;++i) bv.set(17+i*1031); break;
    case 7: bv.set_range(0,65535);bv.set_range(300*65536,306*65536-1);break;
    case 8: for (unsigned i=0;i<320;++i) bv.set_range(i*167,i*167+i%99+1);break;
    case 9:
        bv.set_range(65536,131071);
        for (unsigned i=0;i<16;++i) bv.clear_bit(65536+17+i*131);
        break;
    case 10: for (unsigned i=0;i<65000;++i) bv.set(17+i*67);break;
    case 11:
        for (unsigned i=0;i<6000;++i) {bv.set(bm::set_sub_total_bits+17+i*101);bv.set(512*65536+13+i*7);}
        break;
    case 12: bv.set(0);bv.set(bm::id_max-1);break;
    default: assert(false);
    }
    bv.optimize();return bv;
}

static void check(const std::string& path, Writer& writer, const BV& source,
                  const Writer::bv_ref_vector_type* refs=nullptr)
{
#ifdef BM_WRITE_RELEASE_FIXTURES
    Writer::buffer blob;writer.serialize(source,blob);
    std::ofstream output(path,std::ios::binary);
    output.write(reinterpret_cast<const char*>(blob.data()),std::streamsize(blob.size()));
    assert(output.good());
#else
    (void)writer;
    std::ifstream input(path,std::ios::binary);assert(input.good());
    std::vector<unsigned char> blob((std::istreambuf_iterator<char>(input)),{});
    assert(!blob.empty());
    BV restored;bm::deserialize(restored,blob.data(),0,refs);assert(restored==source);
    BV merged;merged.set(42);BV expected(source);expected.set(42);
    bm::deserialize(merged,blob.data(),0,refs);assert(merged==expected);
    std::stringstream stream(std::string(reinterpret_cast<const char*>(blob.data()),blob.size()));
    bm::streams_decoder channel(stream);bm::streams_deserializer<BV> reader;
    reader.set_ref_vectors(refs);BV streamed;
    assert(reader.deserialize(streamed,channel) && streamed==source && channel.tell()==blob.size());
#endif
    (void)refs;
}

int main(int argc,char** argv)
{
    assert(argc==2);unsigned files=0;
    for (unsigned pattern=0;pattern<13;++pattern)
    for (unsigned level=0;level<=6;++level)
    for (bool bookmarks:{false,true})
    {
        BV bv=fixture(pattern);Writer writer;writer.set_compression_level(level);
        writer.set_bookmarks(bookmarks,4);
        check(std::string(argv[1])+"/"+std::to_string(pattern)+"-"+std::to_string(level)+"-"+std::to_string(bookmarks)+".bv",writer,bv);
        ++files;
    }
    BV target=fixture(2),reference(target);
    for(unsigned k=0;k<20;++k) reference.flip(k*197+1);
    reference.optimize();Writer::bv_ref_vector_type refs;refs.add(&target,0);refs.add(&reference,1);
    Writer::xor_sim_model_type model;Writer writer;
    writer.set_ref_vectors(&refs);writer.compute_sim_model(model,refs,bm::xor_sim_params());
    writer.set_sim_model(&model);writer.set_curr_ref_idx(0);
    for (unsigned level:{0u,6u})
    {
        writer.set_compression_level(level);
        check(std::string(argv[1])+"/ref-"+std::to_string(level)+".bv",writer,target,&refs);++files;
    }
#ifdef BM_WRITE_RELEASE_FIXTURES
    std::cout<<"Wrote "<<files<<" released fixtures PASS\n";
#else
    std::cout<<"Decoded "<<files<<" released fixtures with current library: empty/OR/stream, levels 0..6, bookmarks and reference XOR PASS\n";
#endif
}
