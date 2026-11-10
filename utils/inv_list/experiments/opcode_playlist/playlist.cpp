#include "bm.h"
#include "bmserial.h"
#include "bmfio.h"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
using BV=bm::bvector<>;
static BV read(const char* path) {
    std::ifstream in(path,std::ios::binary);
    if(!in)throw std::runtime_error("Cannot open playlist");
    bm::streams_decoder dec(in);bm::streams_deserializer<BV> reader;BV bv;
    if(!reader.deserialize(bv,dec) || in.peek()!=std::char_traits<char>::eof())
        throw std::runtime_error("Playlist is truncated or has trailing data");
    return bv;
}
int main(int argc,char** argv) {
    try {
        if(argc!=4)throw std::runtime_error("Usage: gov2_playlist pack INDICES.txt OUT.bv | unpack IN.bv OUT.txt");
        if(std::string(argv[1])=="pack") {
            std::ifstream input(argv[2]);if(!input)throw std::runtime_error("Cannot read indices");
            BV bv;std::uint64_t n;std::uint64_t count=0;
            while(input>>n) {
                if(n>=bm::id_max || bv.test(unsigned(n)))throw std::runtime_error("Duplicate/invalid playlist index");
                bv.set(unsigned(n));++count;
            }
            if(!input.eof() || !count)throw std::runtime_error("Empty/malformed index file");
            bv.optimize();bm::serializer<BV> serializer;bm::serializer<BV>::buffer blob;
            serializer.serialize(bv,blob);
            std::ofstream output(argv[3],std::ios::binary);
            output.write(reinterpret_cast<const char*>(blob.data()),std::streamsize(blob.size()));
            output.close();if(!output)throw std::runtime_error("Cannot write playlist");
            if(read(argv[3])!=bv)throw std::runtime_error("Playlist round trip failed");
            std::cout<<"Playlist verified indices="<<count<<" bytes="<<blob.size()<<'\n';
        } else if(std::string(argv[1])=="unpack") {
            BV bv=read(argv[2]);std::ofstream out(argv[3]);
            for(auto e=bv.first();e.valid();++e)out<<*e<<'\n';
            out.close();if(!out)throw std::runtime_error("Cannot write decoded indices");
        } else throw std::runtime_error("Unknown playlist operation");
    } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
