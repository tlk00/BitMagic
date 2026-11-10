// Single-thread serializer benchmark. Input preparation and validation are untimed.
#include "bm.h"
#include "bmserial.h"
#include <chrono>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <stdexcept>
int main(int argc,char**argv) try {
    if(argc!=4) throw std::runtime_error("bench_endpoint INPUT OUTPUT MAX_VECTORS");
    std::ifstream in(argv[1],std::ios::binary);
    std::ofstream out(argv[2],std::ios::binary);
    if(!in || !out) throw std::runtime_error("Cannot open input/output");
    const unsigned limit=unsigned(std::stoul(argv[3]));
    using clock=std::chrono::steady_clock;
    double encode=0,decode=0; unsigned vectors=0;
    bm::id64_t ids=0,bytes=0; std::vector<unsigned> input;
    BM_DECLARE_TEMP_BLOCK(tb)
    for(;vectors<limit;++vectors) {
        unsigned n=0;
        if(!in.read(reinterpret_cast<char*>(&n),4) || !n) throw std::runtime_error("Invalid record length");
        input.resize(n);
        if(!in.read(reinterpret_cast<char*>(input.data()),std::streamsize(n)*4)) throw std::runtime_error("Truncated record");
        bm::bvector<> bv;
        bv.set(input.data(),n,bm::BM_SORTED);bv.optimize(tb);
        bm::serializer<bm::bvector<>> serializer;
        serializer.byte_order_serialization(false);serializer.gap_length_serialization(false);
        serializer.set_compression_level(6);serializer.set_block_endpoint_selection(true);
        bm::serializer<bm::bvector<>>::buffer blob;
        auto begin=clock::now();serializer.serialize(bv,blob);
        encode+=std::chrono::duration<double>(clock::now()-begin).count();
        bm::bvector<> restored;begin=clock::now();bm::deserialize(restored,blob.data());
        decode+=std::chrono::duration<double>(clock::now()-begin).count();
        if(restored!=bv) throw std::runtime_error("Round trip mismatch");
        out.write(reinterpret_cast<const char*>(blob.data()),std::streamsize(blob.size()));
        if(!out) throw std::runtime_error("Output write failed");
        ids+=n;bytes+=blob.size();
        if((vectors+1)%1000==0) std::cerr<<"vectors="<<vectors+1<<" encode_s="<<encode<<'\n';
    }
    out.close();if(!out) throw std::runtime_error("Output close failed");
    std::cout<<std::setprecision(12)<<"{\"vectors\":"<<vectors<<",\"ids\":"<<ids<<",\"bytes\":"<<bytes
             <<",\"encode_s\":"<<encode<<",\"decode_s\":"<<decode<<",\"verified\":true}\n";
} catch(const std::exception&e) {std::cerr<<e.what()<<'\n';return 1;}
