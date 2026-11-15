/** @file profile_examples.cpp
    Research-only opcode attribution for retained zstd survey BLOBs.
    Read .bv files in one examples directory and validate each through the
    production opcode collector. Emit tab-separated per-vector opcode counts,
    complete bytes and logical ID counts. Caller owns output redirection.
*/
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <stdexcept>
#include <vector>
#include "bm.h"
#include "bmserial.h"
#include "bmserial_profile.h"
#include "bmdbg.h"

/// Profile borrowed filenames from one directory. Return nonzero on any failure.
int main(int argc, char** argv)
{
    try
    {
        if (argc != 2)
        {
            throw std::runtime_error("Expected one examples directory");
        }
        bm::serialization_opcode_collector<bm::bvector<> > collector;
        std::cout << "file\topcode\tname\trecords\tbytes\tids\n";
        for (const auto& file : std::filesystem::directory_iterator(argv[1]))
        {
            if (file.path().extension() != ".bv" ||
                file.path().filename().string().compare(0, 2, "._") == 0)
            {
                continue;
            }
            std::ifstream input(file.path(), std::ios::binary);
            std::vector<unsigned char> blob((std::istreambuf_iterator<char>(input)),
                                            std::istreambuf_iterator<char>());
            if (input.bad() || blob.empty())
            {
                throw std::runtime_error("Unreadable example");
            }
            bm::serialization_opcode_profile profile;
            collector.add_blob(blob.data(), blob.size(), profile);
            if (!profile.reconciles())
            {
                throw std::logic_error("Profile mismatch");
            }
            for (unsigned code = 0; code < 256; ++code)
            {
                const auto& r = profile.opcode[code];
                if (r.count)
                {
                    std::cout << file.path().filename().string() << '\t' << code << '\t'
                              << bm::serialization_opcode_name(code) << '\t'
                              << r.count << '\t' << r.bytes << '\t' << r.ids << '\n';
                }
            } // for code
        } // for file
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
