/*
Copyright(c) 2002-2017 Anatoliy Kuznetsov(anatoliy_kuznetsov at yahoo.com)

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

For more information please visit:  http://bitmagic.io
*/

#include <iostream>
#include <chrono>
#include <thread>
#include <time.h>
#include <stdio.h>


#ifdef _MSC_VER
#pragma warning( push )
#pragma warning( disable : 4996)
#endif

#include <vector>
#include <chrono>
#include <map>

#include "bm.h"
#include "bmalgo.h"
#include "bmserial.h"
#include "bmsparsevec.h"
#include "bmsparsevec_compr.h"
#include "bmsparsevec_algo.h"
#include "bmsparsevec_serial.h"
#include "bmalgo_similarity.h"


#include "bmdbg.h"
#include "bmtimer.h"

using namespace std;

static
void show_help()
{
    std::cerr
      << "BitMagic Sparse Vector Analysis Utility (c) 2017"            << std::endl
      << "-bvin  bv-file              -- bv file to load"              << std::endl
      << "-svin  sv-input-file        -- 32-bit sparse_vector file to load" << std::endl
      << "-rsc64in rsc-64-bit-file    -- 64-bit rsc sparse vector to load"  << std::endl
      << "-u32in u32-input-file       -- raw 32-bit unsigned int file" << std::endl
      << "-svout sv-output-file       -- sv output file to produce"    << std::endl
      << "-u32out u32-output-file     -- raw 32-bit output file to produce" << std::endl
      << "-diag (-d)                  -- print statistics/diagnostics info" << std::endl
      << "-timing (-t)                -- evaluate timing/duration of operations" << std::endl
      ;
}




// Arguments
//
std::string  bv_in_file;
std::string  sv_in_file;
std::string  rsc64_in_file;
std::string  u32_in_file;
std::string  sv_out_file;
std::string  u32_out_file;
bool         is_diag = false;
bool         is_timing = false;

static
int parse_args(int argc, char *argv[])
{
    for (int i = 1; i < argc; ++i)
    {
        std::string arg = argv[i];
        if ((arg == "-h") || (arg == "--help"))
        {
            show_help();
            return 0;
        }
        
        if (arg == "-svin" || arg == "--svin")
        {
            if (i + 1 < argc)
            {
                sv_in_file = argv[++i];
            }
            else
            {
                std::cerr << "Error: -svin requires file name" << std::endl;
                return 1;
            }
            continue;
        }

        if (arg == "-rsc64in" || arg == "--rsc64in")
        {
            if (i + 1 < argc)
            {
                rsc64_in_file = argv[++i];
            }
            else
            {
                std::cerr << "Error: -rsc64in requires file name" << std::endl;
                return 1;
            }
            continue;
        }

        if (arg == "-bvin" || arg == "--bvin")
        {
            if (i + 1 < argc)
            {
                bv_in_file = argv[++i];
            }
            else
            {
                return 1;
                std::cerr << "Error: -bvin requires file name" << std::endl;
            }
            continue;
        }

        if (arg == "-u32in" || arg == "--u32in")
        {
            if (i + 1 < argc)
            {
                u32_in_file = argv[++i];
            }
            else
            {
                std::cerr << "Error: -u32in requires file name" << std::endl;
                return 1;
            }
            continue;
        }
        
        if (arg == "-svout" || arg == "--svout")
        {
            if (i + 1 < argc)
            {
                sv_out_file = argv[++i];
            }
            else
            {
                std::cerr << "Error: -svout requires file name" << std::endl;
                return 1;
            }
            continue;
        }

        if (arg == "-u32out" || arg == "--u32out")
        {
            if (i + 1 < argc)
            {
                u32_out_file = argv[++i];
            }
            else
            {
                std::cerr << "Error: -u32out requires file name" << std::endl;
                return 1;
            }
            continue;
        }


        if (arg == "-diag" || arg == "--diag" || arg == "-d" || arg == "--d")
            is_diag = true;
        if (arg == "-timing" || arg == "--timing" || arg == "-t" || arg == "--t")
            is_timing = true;
        
        
    } // for i
    return 0;
}


// Globals
//
typedef bm::sparse_vector<unsigned, bm::bvector<> > sparse_vector_u32;
typedef bm::sparse_vector<unsigned long long, bm::bvector<> > sparse_vector_u64;
typedef bm::rsc_sparse_vector<unsigned, sparse_vector_u32> rsc_sparse_vector_u32;
typedef bm::rsc_sparse_vector<bm::id64_t, sparse_vector_u64> rsc_sparse_vector_u64;


bm::bvector<>          bv_inp;
sparse_vector_u32      sv_u32_in;
sparse_vector_u64      sv_u64_in;
rsc_sparse_vector_u64  rsc_u64_in(bm::use_null);
sparse_vector_u32      sv_u32_out;
bool                   sv_u32_in_flag = false;
bool                   rsc_u64_in_flag = false;
std::vector<unsigned>  vect_u32_in;
std::vector<unsigned>  vect_u32_out;

bm::chrono_taker<>::duration_map_type  timing_map;


// load bvector from a file
//
static
void load_bv(const std::string& fname, bm::bvector<>& bv)
{
    bm::LoadBVector(fname.c_str(), bv);
}

// load sparse_vector from a file
//
static
int load_sv(const std::string& fname, sparse_vector_u32& sv)
{
    std::vector<unsigned char> buffer;

    // read the input buffer, validate errors
    bm::chrono_taker tt1(cout, "serialized sparse vector BLOB read", 1, &timing_map);
    auto ret = bm::read_dump_file(fname, buffer);
    
    tt1.stop(is_timing);
    
    if (ret != 0)
    {
        std::cerr << "Failed to read file:" << fname << std::endl;
        return 2;
    }
    if (buffer.size() == 0)
    {
        std::cerr << "Empty input file:" << fname << std::endl;
        return 3;
    }
    
    // deserialize
    //
    bm::chrono_taker tt2(cout, "sparse vector deserialization", 1, &timing_map);
    const unsigned char* buf = &buffer[0];
    BM_DECLARE_TEMP_BLOCK(tb)
    auto res = bm::sparse_vector_deserialize(sv, buf, tb);
    tt2.stop(is_timing);
    if (res != 0)
    {
        std::cerr << "Sparse vector deserialization failed ("
                  << fname << ")"
                  << std::endl;
        return 4;
    }
    return 0;
}

// load sparse_vector from a file
//
static
int load_rsc64(const std::string& fname, rsc_sparse_vector_u64& csv)
{
    std::vector<unsigned char> buffer;

    // read the input buffer, validate errors
    bm::chrono_taker tt1(cout, "serialized rsc sparse vector BLOB read", 1, &timing_map);
    auto ret = bm::read_dump_file(fname, buffer);
    
    tt1.stop(is_timing);
    
    if (ret != 0)
    {
        std::cerr << "Failed to read file:" << fname << std::endl;
        return 2;
    }
    if (buffer.size() == 0)
    {
        std::cerr << "Empty input file:" << fname << std::endl;
        return 3;
    }
    
    // deserialize
    //
    bm::chrono_taker tt2(cout, "rsc sparse vector deserialization", 1, &timing_map);
    const unsigned char* buf = &buffer[0];
    BM_DECLARE_TEMP_BLOCK(tb)
    auto res = bm::sparse_vector_deserialize(csv, buf, tb);
    tt2.stop(is_timing);
    if (res != 0)
    {
        std::cerr << "rsc sparse vector deserialization failed ("
                  << fname << ")"
                  << std::endl;
        return 4;
    }
    /*
    std::vector<unsigned long long> vect;
    vect.resize(13107001);
    csv.decode(vect.data(), 26214000, 13107000);
    */
    return 0;
}


// load raw unsigned file
//
static
int load_u32(const std::string& fname, std::vector<unsigned>& vect)
{
    bm::chrono_taker tt(cout, "u32 BLOB read", 1, &timing_map);

    auto ret = bm::read_dump_file(fname, vect);
    
    tt.stop(is_timing);
    
    if (ret != 0)
    {
        std::cerr << "Failed to read file:" << fname << std::endl;
        return 2;
    }
    if (vect.size() == 0)
    {
        std::cerr << "Empty input file:" << fname << std::endl;
        return 3;
    }
    return 0;
}

// convert unsigned vector to sparse format
//
static
int convert_u32(const std::vector<unsigned>& u32, sparse_vector_u32& sv)
{
    BM_DECLARE_TEMP_BLOCK(tb)
    bm::chrono_taker tt(cout, "u32 array to sparse vector transposition conversion", 1, &timing_map);
    
    sv.import(&u32[0], (unsigned)u32.size());
    sv.optimize(tb);
    return 0;
}



int main(int argc, char *argv[])
{
    if (argc < 3)
    {
        show_help();
        return 1;
    }
    
    try
    {
        auto ret = parse_args(argc, argv);
        if (ret != 0)
            return ret;
  
        if (!bv_in_file.empty())
        {
            load_bv(bv_in_file, bv_inp);
        }
        
        if (!sv_in_file.empty())
        {
            auto res = load_sv(sv_in_file, sv_u32_in);
            if (res != 0)
            {
                return res;
            }
            sv_u32_in_flag = true;
        }

        if (!rsc64_in_file.empty())
        {
            auto res = load_rsc64(rsc64_in_file, rsc_u64_in);
            if (res != 0)
            {
                return res;
            }
            rsc_u64_in_flag = true;
        }

        if (!u32_in_file.empty())
        {
            auto res = load_u32(u32_in_file, vect_u32_in);
            if (res != 0)
            {
                return res;
            }
        }
        
        if (is_diag)  // diagnostics required
        {
            if (sv_u32_in_flag)  // input sparse vector loaded
            {
                std::cout << "Input sparse vector statistics:" << std::endl;
                bm::print_svector_stat(std::cout, sv_u32_in, true);
                std::cout << std::endl;
            }
            
            if (rsc_u64_in_flag)
            {
                std::cout << "Input rsc 64-bit sparse vector statistics:" << std::endl;
                bm::print_svector_stat(std::cout, rsc_u64_in, true);
                std::cout << std::endl;
            }

            if (!vect_u32_in.empty())
            {
                std::cout << "Input u32 raw vector size = "
                          << vect_u32_in.size() << " elements."
                          << std::endl;
            }
        }
        
        if (is_timing && sv_u32_in_flag)
        {
            sparse_vector_u32::bvector_type bv_mask;
            sparse_vector_u32::bvector_type bv_out_control;
            sparse_vector_u32::bvector_type bv_out;
            if (bv_inp.any())
                bv_mask = bv_inp;
            else
                bv_mask.set_range(1, 7000000);
            std::cout << "remap bit-vector count = " << bv_mask.count() << std::endl;

            //std::this_thread::sleep_for(std::chrono::seconds(25));
            
            bm::set2set_11_transform<sparse_vector_u32> set2set;
            
            {
                bm::chrono_taker tt(cout, "set2set transform one-by-one(control) remap", 1, &timing_map);
                typename sparse_vector_u32::bvector_type::enumerator en(bv_mask.first());
                for (; en.valid(); ++en)
                {
                    bm::id_t in = *en;
                    bm::id_t out;
                    bool found = set2set.remap(in, sv_u32_in, out);
                    if (found)
                    {
                        bv_out_control.set(out);
                    }
                }

            }


            {
                bm::chrono_taker tt(cout, "set2set transform remap", 1, &timing_map);
                set2set.remap(bv_mask, sv_u32_in, bv_out);
            }
            
            int res = bv_out_control.compare(bv_out);
            if (res != 0)
            {
                std::cerr << "ERROR:" << " control remap mismatch." << std::endl;
                return 1;
            }
            
        }
        
        
        if (!sv_out_file.empty()) // request to make new SV compressed file
        {
            if (!vect_u32_in.empty())
            {
                auto res = convert_u32(vect_u32_in, sv_u32_out);
                if (res != 0)
                    return res;
               
                if (is_diag)  // diagnostics requested
                {
                    std::cout << "Output sparse vector statistics:" << std::endl;
                    bm::print_svector_stat(std::cout, sv_u32_out);
                    std::cout << std::endl;
                }
                
                size_t sv_blob_size = 0;
                
                bm::chrono_taker tt(cout, "sparse vector BLOB save", 1, &timing_map);
                res = bm::file_save_svector(sv_u32_out, sv_out_file, &sv_blob_size);
                tt.stop(is_timing);
                
                if (res != 0)
                {
                    std::cerr << "Failed to save sparse vector file: " << sv_out_file << std::endl;
                    return res;
                }
                if (is_diag)
                    std::cout << "Output sparse vector BLOB size: " << sv_blob_size << std::endl;
            }
            
            if (sv_u32_in_flag) // input data is ready as a sparse vector
            {
                
            }
            
        } // if sv_out_file

        if (!u32_out_file.empty()) // request to de-compressed bmsv file
        {
            if (!sv_u32_in.empty())
            {
                vect_u32_out.resize(sv_u32_in.size());
                {
                    bm::chrono_taker tt(cout, "sparse vector decode", 1, &timing_map);
                    sv_u32_in.decode(&vect_u32_out[0], 0, sv_u32_in.size(), false);
                    tt.stop(is_timing);
                }
                {
                    bm::chrono_taker tt(cout, "u32 vector write", 1, &timing_map);
                    std::ofstream fout(u32_out_file.c_str(), std::ios::binary);
                    if (!fout.good())
                    {
                        std::cerr << "Cannot open file " << u32_out_file << std::endl;
                        return 1;
                    }
                    const char* buf = (const char*)&vect_u32_out[0];
                    fout.write(buf, std::streamsize(vect_u32_out.size() * sizeof(unsigned)));
                    if (!fout.good())
                    {
                        return 2;
                    }
                    fout.close();
                }
            }
        } // if u32_out_file
        
        if (is_diag)
        {
            // diagnostics comparisons
            
            // in/out sparse vectors
            if (!sv_u32_in.empty() && !sv_u32_out.empty())
            {
                bm::chrono_taker tt(cout, "sparse vectors in/out comparison", 1, &timing_map);
                bool eq = sv_u32_in.equal(sv_u32_out);
                if (!eq)
                {
                    std::cerr << "ERROR: input sparse vector is different from output." << std::endl;
                }
            }


            // input sparse compare to input raw
            if (!sv_u32_in.empty() && !vect_u32_in.empty())
            {
                if (sv_u32_in.size() != vect_u32_in.size())
                {
                    std::cerr << "ERROR: input sparse vector size is different from input raw array." << std::endl;
                }
                else
                {
                    bm::chrono_taker tt(cout, "sparse vector in/raw comparison", 1, &timing_map);
                    int res = bm::svector_check(sv_u32_in, vect_u32_in);
                    if (res != 0)
                    {
                        std::cerr << "ERROR: input sparse vector is different from input raw array." << std::endl;
                    }
                }
            }

            // input sparse compare to output raw
            if (!sv_u32_in.empty() && !vect_u32_out.empty())
            {
                bm::chrono_taker tt(cout, "sparse vector in/raw comparison", 1, &timing_map);
                int res = bm::svector_check(sv_u32_in, vect_u32_out);
                if (res != 0)
                {
                    std::cerr << "ERROR: input sparse vector is different from input raw array." << std::endl;
                }
            }
            
            if (!vect_u32_in.empty() && !sv_u32_out.empty())
            {
                bm::chrono_taker tt(cout, "raw in to sparse vector out comparison", 1, &timing_map);
                int res = bm::svector_check(sv_u32_out, vect_u32_in);
                if (res != 0)
                {
                    std::cerr << "ERROR: input raw array is different from output sparse vector." << std::endl;
                }
                
            }
        }
        
        if (is_timing)  // print all collected timings
        {
            std::cout << std::endl << "Timings (ms):" << std::endl;
            bm::chrono_taker<>::print_duration_map(cout, timing_map);
        }
    }
    catch (std::exception& ex)
    {
        std::cerr << "Error:" << ex.what() << std::endl;
        return 1;
    }

    return 0;
}


#ifdef _MSC_VER
#pragma warning( pop )
#endif



