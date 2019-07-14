/*
Copyright(c) 2002-2019 Anatoliy Kuznetsov(anatoliy_kuznetsov at yahoo.com)

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
    std::cout
      << "BitMagic Inverted List Compression Test (c) 2019"            << std::endl
      << "-u32in u32-input-file       -- raw 32-bit unsigned int file" << std::endl
      << "-bvout bvect-output-file    -- bit-vector compressed out file" << std::endl
      << "-bvin bvect-input-file      -- bit-vector compressed in file" << std::endl
      << std::endl
      << "-check                      -- check/validate input list" << std::endl
      << "-verify                     -- verify compressed version " << std::endl
      << "-diag (-d)                  -- print statistics/diagnostics info" << std::endl
      << "-timing (-t)                -- evaluate timing/duration of operations" << std::endl
      ;
}




// Arguments
//
std::string  bv_in_file;
std::string  bv_out_file;

std::string  sv_in_file;
std::string  sv_out_file;
std::string  u32_in_file;
std::string  u32_out_file;

bool         is_diag = false;
bool         is_timing = false;
bool         is_check = false;
bool         is_verify = false;




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
        
        if (arg == "-c" || arg == "-check")
        {
            if (i + 1 < argc)
            {
                is_check = true;
            }
            continue;
        }
        if (arg == "-v" || arg == "-verify")
        {
            if (i + 1 < argc)
            {
                is_verify = true;
            }
            continue;
        }

        if (arg == "-bvout" || arg == "--bvout")
        {
            if (i + 1 < argc)
            {
                bv_out_file = argv[++i];
            }
            else
            {
                std::cerr << "Error: -bvout requires file name" << std::endl;
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
                std::cerr << "Error: -bvin requires file name" << std::endl;
                return 1;
            }
            continue;
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
typedef bm::rsc_sparse_vector<unsigned, sparse_vector_u32> rsc_sparse_vector_u32;


bm::chrono_taker::duration_map_type  timing_map;


/// Read 32-bit vector size-prefix format (length:0, 1, 2, 3, ....)
///
template<class VT>
int io_read_u32_coll(std::ifstream& fin, VT& vec)
{
    typedef typename VT::value_type  value_type;
    vec.resize(0);
    if (!fin.good())
        return -1;
    value_type len;
    fin.read((char*) &len, std::streamsize(sizeof(len)));
    if (!fin.good())
        return -1;
    if (!len)
        return -2; // 0-len detected (broken file)
    
    vec.resize(len);
    fin.read((char*) &vec[0], std::streamsize(len*sizeof(value_type)));
    if (!fin.good())
        return -1;

    return 0;
}

/// Check if input vector is monotonously sorted (true inverted list)
///
template<typename VT>
int validate_inp_vec(const VT& vec)
{
    auto sz = vec.size();
    if (!sz)
        return -1;
    auto i_prev = vec[0];
    for (typename VT::size_type i = 1; i < sz; ++i)
    {
        auto val = vec[i];
        if (val <= i_prev)
            return -2;
        i_prev = val;
    }
    return 0;
}

template<typename VT, typename BV>
int compare_vect(const VT& vec, const BV& bv)
{
    if (vec.size() != bv.count())
    {
        return -1;
    }
    typename BV::enumerator en = bv.first();
    typename VT::const_iterator it = vec.begin();
    for (; en.valid(); ++en, ++it)
    {
        if (*en != *it)
        {
            return -1;
        }
    }
    return 0;
}



///  convert vector into bit-vector and append to file
///
template<typename VT>
void write_as_bvector(std::ofstream& bv_file,
                     const VT& vec,
                     bm::serializer<bm::bvector<> >& bvs,
                     bm::serializer<bm::bvector<> >::buffer& sbuf)
{
    bm::bvector<> bv;
    bv.set(&vec[0], bm::bvector<>::size_type(vec.size()), bm::BM_SORTED);

    bvs.optimize_serialize_destroy(bv, sbuf);

    unsigned bv_size = (unsigned)sbuf.size();
    bv_file.write((char*)&bv_size, sizeof(bv_size));
    bv_file.write((char*)sbuf.data(), (std::streamsize)sbuf.size());
    if (!bv_file.good())
        throw std::runtime_error("Error write to bvect out file");
}

/// read and desrialize bit-bector from the dump file
///
static
int read_bvector(std::ifstream& bv_file,
                  bm::bvector<>& bv,
                  bm::serializer<bm::bvector<> >::buffer& sbuf)
{
    if (!bv_file.good())
        return -1;
    unsigned len;
    bv_file.read((char*) &len, std::streamsize(sizeof(len)));
    if (!bv_file.good())
        return -1;
    if (!len)
        return -2; // 0-len detected (broken file)
    
    sbuf.resize(len);
    bv_file.read((char*) sbuf.data(), std::streamsize(len));
    if (!bv_file.good())
        return -1;
    
    bm::deserialize(bv, sbuf.data());

    return 0;
}

/// read the input collection sequence, write using various compression schemes
///
static
void compress_inv_dump_file(const std::string& fname,
                            const std::string& bv_out_fname,
                            bool  check)
{
    bm::id64_t total_ints = 0;

    cout << "Reading input collection: " << fname << endl;
    if (!bv_out_fname.empty())
        cout << "Writing to BV collection: " << bv_out_fname << endl;
    else
        cout << "NO BV collection specified" << endl;

    vector<unsigned> vec;
    std::ifstream fin(fname.c_str(), std::ios::in | std::ios::binary);
    if (!fin.good())
    {
        throw std::runtime_error("Cannot open input file");
    }
    
    std::ofstream bv_file;
    if (!bv_out_fname.empty())
    {
        bv_file.open(bv_out_fname, std::ios::out | std::ios::binary);
        if (!bv_file.good())
            throw std::runtime_error("Cannot open bvect out file");
    }
    
    
    fin.seekg(0, std::ios::end);
    std::streamsize fsize = fin.tellg();

    fin.seekg(0, std::ios::beg);
    
    // initialize serializer
    //
    bm::serializer<bm::bvector<> > bvs;
    bvs.byte_order_serialization(false);
    bvs.gap_length_serialization(false);
    
    bm::serializer<bm::bvector<> >::buffer sbuf; // resizable memory buffer

    // main loop to read sample vectors
    //
    size_t i;
    for (i = 0; true; ++i)
    {
        int ret = io_read_u32_coll(fin, vec);
        if (ret != 0)
            throw std::runtime_error("Error reading input file");
        if (check)
        {
            ret = validate_inp_vec(vec); // check if we have sorted unique set
            if (ret != 0)
                throw std::runtime_error("Input vector validation failed");
        }
        total_ints += vec.size(); // remember the total size of the collection
        
        // serialize and save as a bit-vector size0:<BLOB0>, size1:<BLOB1>...N
        //
        if (!bv_out_fname.empty())
        {
            write_as_bvector(bv_file, vec, bvs, sbuf);
        }
        
        std::streamsize fpos_curr = fin.tellg();
        if (fpos_curr == fsize)
            break;
        
        cout << "\r" << fpos_curr << "/" << fsize
             << " ( size=" << vec.size() << " )              "
             << flush;
    } // for i
    
    cout << endl;
    cout << "Total vectors=" << i << endl;
    cout << "Total ints=" << total_ints << endl;
    
    if (!bv_out_fname.empty())
    {
        bm::id64_t bv_size = (bm::id64_t)bv_file.tellp();
        cout << "BV size = " << bv_size << endl;
        // calculate bits per int compression ratio corrected not to account
        // for size/length words prefixing the vectors
        double bv_bits_per_int = double(bv_size * 8ull - (i*sizeof(unsigned))) / double(total_ints);
        cout << "BV Bits per/int = " << std::setprecision(3) << bv_bits_per_int << endl;
        bv_file.close();
    }
}

/// read the input collection sequence, write using various compression schemes
///
static
void verify_inv_dump_file(const std::string& fname,
                          const std::string& bv_in_fname,
                          bool  check)
{
    bm::id64_t total_ints = 0;

    cout << "Reading input collection: " << fname << endl;
    if (!bv_in_fname.empty())
        cout << "Reading BV collection: " << bv_in_fname << endl;
    else
        cout << "NO BV collection specified" << endl;


    vector<unsigned> vec;
    std::ifstream fin(fname.c_str(), std::ios::in | std::ios::binary);
    if (!fin.good())
    {
        throw std::runtime_error("Cannot open input file");
    }
    
    std::ifstream bv_file;
    if (!bv_in_fname.empty())
    {
        bv_file.open(bv_in_fname, std::ios::in | std::ios::binary);
        if (!bv_file.good())
            throw std::runtime_error("Cannot open bvect dump file");
    }
    
    
    fin.seekg(0, std::ios::end);
    std::streamsize fsize = fin.tellg();

    fin.seekg(0, std::ios::beg);
    
    // initialize serializer
    //
    
    bm::serializer<bm::bvector<> >::buffer sbuf; // resizable memory buffer

    // main loop to read sample vectors
    //
    size_t i;
    for (i = 0; true; ++i)
    {
        int ret = io_read_u32_coll(fin, vec);
        if (ret != 0)
            throw std::runtime_error("Error reading input file");
        if (check)
        {
            ret = validate_inp_vec(vec); // check if we have sorted unique set
            if (ret != 0)
                throw std::runtime_error("Input vector validation failed");
        }
        total_ints += vec.size(); // remember the total size of the collection
        
        // serialize and save as a bit-vector size0:<BLOB0>, size1:<BLOB1>...N
        //
        if (!bv_in_fname.empty())
        {
            bm::bvector<> bv;
            read_bvector(bv_file, bv, sbuf);
            int cmp = compare_vect(vec, bv);
            if (cmp != 0)
            {
                throw std::runtime_error("Vector comparison failed");
            }
        }
        
        std::streamsize fpos_curr = fin.tellg();
        if (fpos_curr == fsize)
            break;
        
        cout << "\r" << fpos_curr << "/" << fsize
             << " ( size=" << vec.size() << " )              "
             << flush;
    } // for i
    
    cout << endl;
    cout << "Verification complete." << endl;
    cout << "Total vectors=" << i << endl;
    cout << "Total ints=" << total_ints << endl;
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
        
        if (!u32_in_file.empty())
        {
            if (!is_verify)
            {
            cout << "compression" << endl;
                compress_inv_dump_file(u32_in_file, bv_out_file, is_check);
            }
            else
            {
                cout << "Verification." << endl;
                verify_inv_dump_file(u32_in_file, bv_in_file, is_check);
            }
        }
        
        
        if (is_timing)  // print all collected timings
        {
            std::cout << std::endl << "Timings (ms):" << std::endl;
            bm::chrono_taker::print_duration_map(timing_map);
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



