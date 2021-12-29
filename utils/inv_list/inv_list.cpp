/*
Copyright(c) 2002-2020 Anatoliy Kuznetsov(anatoliy_kuznetsov at yahoo.com)

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

/** \example inv_list.cpp
  Utility to compress test sets of inverted lists (Gov2 corpus)
  and study compression characteristics and different comression levels
*/

/*! \file inv_list.cpp
    \brief Utility to compress test sets of inverted lists
*/

#include <iostream>
#include <chrono>
#include <thread>
#include <time.h>
#include <stdio.h>
#include <cstdlib>


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
#include "bmundef.h" /* clear the pre-proc defines from BM */

using namespace std;

static
void show_help()
{
    std::cout
        << "BitMagic Inverted List Compression Test (c) 2019" << std::endl
        << "-u32in u32-input-file       -- raw 32-bit unsigned int file" << std::endl
        << "-bvout bvect-output-file    -- bit-vector compressed out file" << std::endl
        << "-svout svect-output-file    -- bit-transposed sparse vectors out file" << std::endl
        << "-bvin bvect-input-file      -- bit-vector compressed in file" << std::endl
        << std::endl
        << "-level N                    -- compression level to use (up to 5)"     << std::endl
        << "-silent (-s)                -- no progress print or messages"          << std::endl
        << "-verify                     -- verify compressed version "             << std::endl
        << "-decode                     -- run decode test (in-memory)"            << std::endl
        << "-diag (-d)                  -- print statistics/diagnostics info"      << std::endl
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
bool         is_verify = false;
bool         is_silent = false;
bool         is_decode = false;

unsigned     c_level = bm::set_compression_default;


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
        
        if (arg == "-v" || arg == "-verify")
        {
            if (i + 1 < argc)
            {
                is_verify = true;
            }
            continue;
        }
        if (arg == "-decode")
        {
            if (i + 1 < argc)
            {
                is_decode = true;
            }
            continue;
        }
        if (arg == "-l" || arg == "-level")
        {
            if (i + 1 < argc)
            {
                const char* lvl = argv[++i];
                char *end;
                c_level = (unsigned) std::strtoul(lvl, &end, 10);
                if (errno == ERANGE)
                {
                    std::cerr << "Error parsing -level: range error for "
                              << lvl<< std::endl;
                    return 1;
                }
            }
            else
            {
                std::cerr << "Error: -level requires compression level number" << std::endl;
                return 1;
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


        if (arg == "-silent" || arg == "--silent" || arg == "-s" || arg == "--s")
            is_silent = true;
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


bm::chrono_taker<>::duration_map_type  timing_map;


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
/// along the way in computes a minimal delta between values
///
template<typename VT>
int validate_inp_vec(const VT& vec,
                     typename VT::value_type& min_delta,
                     typename VT::value_type& min_delta_cnt
                     )
{
    typename VT::value_type md_cnt = min_delta_cnt = 0;
    auto sz = vec.size();
    if (!sz)
        return -1;
    auto i_prev = vec[0];
    typename VT::value_type md = ~(typename VT::value_type(0)); // max possible uint
    for (typename VT::size_type i = 1; i < sz; ++i)
    {
        auto val = vec[i];
        if (val <= i_prev)
            return -2;
        typename VT::value_type d1 = val - i_prev;
        if (d1 < md)
        {
            md_cnt = 0; md = d1;
        }
        md_cnt += (d1 == md);
        i_prev = val;
    }
    min_delta = md;
    min_delta_cnt = md_cnt;
    return 0;
}

/// Verification check if integer vector is equivalent to a bit-vector
///
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
            return -1;
    }
    return 0;
}

/// Debug utility to detect super sparse bit-vectors
/// which probably get bad compression rate
///
template<typename BV>
bool is_super_sparse(const BV& bv)
{
    typename BV::statistics st;
    bv.calc_stat(&st);
    auto bc = bv.count();
    auto blocks_count = st.gap_blocks + st.bit_blocks;
    if (bc <= blocks_count)
        return true;
    auto bc_parity = blocks_count * 6;
    return (bc <= bc_parity);
}


///  convert vector into bit-vector and append to the file
///
/// @return true if vector was detected as very low cardinality
///
template<typename VT>
bool write_as_bvector(std::ofstream& bv_file,
                     const VT& vec,
                     bm::serializer<bm::bvector<> >& bvs,
                     bm::serializer<bm::bvector<> >::buffer& sbuf)
{
    BM_DECLARE_TEMP_BLOCK(tb)
    bm::bvector<> bv;
    bv.set(&vec[0], bm::bvector<>::size_type(vec.size()), bm::BM_SORTED);

    bv.optimize(tb);
    bvs.serialize(bv, sbuf);

    unsigned bv_size = (unsigned)sbuf.size();
    bv_file.write((char*)&bv_size, sizeof(bv_size));
    bv_file.write((char*)sbuf.data(), (std::streamsize)sbuf.size());
    if (!bv_file.good())
        throw std::runtime_error("Error write to bvect out file");
    return false;
}

///  convert vector into delta coded bit-transposed vector and append to the file
///
template<typename VT>
void write_as_svector(std::ofstream& sv_file,
                     const VT& vec,
                     unsigned min_delta,
                     bm::sparse_vector_serial_layout<sparse_vector_u32>& sv_lay)
{

    sparse_vector_u32 sv;
    bm::id_t prev = vec[0];

    {
        sparse_vector_u32::back_insert_iterator sv_bi = sv.get_back_inserter();
        sv_bi.add(min_delta);
        sv_bi.add(prev);

        for (unsigned k = 1; k < vec.size(); ++k)
        {
            bm::id_t curr = vec[k];
            bm::id_t delta = curr - prev;
            if (delta < min_delta)
                throw std::runtime_error("Input vector validation delta error");
            delta -= min_delta;
            sv_bi.add(delta);
            prev = curr;
        } // for i
        sv_bi.flush();
    }

    BM_DECLARE_TEMP_BLOCK(tb)
    sv.optimize(tb);
    bm::sparse_vector_serialize(sv, sv_lay, tb);

    unsigned sv_size = (unsigned)sv_lay.size();
    sv_file.write((char*)&sv_size, sizeof(sv_size));
    sv_file.write((char*)sv_lay.data(), (std::streamsize)sv_lay.size());
    if (!sv_file.good())
        throw std::runtime_error("Error write to bvect out file");

}


///  convert vector into delta coded bit-transposed vector and append to the file
///
template<typename VT>
void write_as_rsc_svector(std::ofstream& sv_file,
                          const VT& vec,
                          unsigned min_delta,
                          bm::sparse_vector_serial_layout<rsc_sparse_vector_u32>& sv_lay)
{

    sparse_vector_u32 sv(bm::use_null);
    bm::id_t prev = vec[0];

    {
        sparse_vector_u32::back_insert_iterator sv_bi = sv.get_back_inserter();
        sv_bi.add(min_delta);
        sv_bi.add(prev);

        for (unsigned k = 1; k < vec.size(); ++k)
        {
            bm::id_t curr = vec[k];
            bm::id_t delta = curr - prev;
            if (delta < min_delta)
                throw std::runtime_error("Input vector validation delta error");
            delta -= min_delta;
            if (delta)
                sv_bi.add(delta);
            else
                sv_bi.add_null();
            prev = curr;
        } // for i
        sv_bi.flush();
    }

    BM_DECLARE_TEMP_BLOCK(tb)
 
    rsc_sparse_vector_u32 csv; // compressed sparse vector
    csv.load_from(sv); // load rank-select-compacted (rsc) sparse vector
    csv.optimize(tb);

    bm::sparse_vector_serialize(csv, sv_lay, tb);

    unsigned sv_size = (unsigned)sv_lay.size();
    sv_file.write((char*)&sv_size, sizeof(sv_size));
    sv_file.write((char*)sv_lay.data(), (std::streamsize)sv_lay.size());
    if (!sv_file.good())
        throw std::runtime_error("Error write to bvect out file");

}



/// read the input collection sequence, write using various compression schemes
///
static
void compress_inv_dump_file(const std::string& fname,
                            const std::string& bv_out_fname,
                            const std::string& sv_out_fname)
{
    bm::id64_t total_ints = 0;
    bm::id64_t total_low_card = 0;
    bm::id64_t total_low_card_size = 0;
    bm::id64_t min_delta_ints = 0;
    bm::id64_t sv_size = 0;
    bm::id64_t rsc_diff_size = 0;
    bm::id64_t sv_cnt = 0;

    cout << "Reading input collection: " << fname << endl;
    if (!bv_out_fname.empty())
        cout << "Writing to BV collection: " << bv_out_fname << endl;
    else
        cout << "NO BV collection specified" << endl;
    if (!sv_out_fname.empty())
        cout << "Writing to SV collection: " << sv_out_fname << endl;
    else
        cout << "NO SV collection specified" << endl;

    
    cout << "Compression level: " << c_level << endl;

    bm::chrono_taker tt1(std::cout, "1. Convert collection", 1, &timing_map);

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
    std::ofstream sv_file;
    if (!sv_out_fname.empty())
    {
        sv_file.open(sv_out_fname, std::ios::out | std::ios::binary);
        if (!sv_file.good())
            throw std::runtime_error("Cannot open svect out file");
    }

    
    fin.seekg(0, std::ios::end);
    std::streamsize fsize = fin.tellg();

    fin.seekg(0, std::ios::beg);
    
    // initialize serializer
    //
    bm::serializer<bm::bvector<> > bvs;
    bvs.byte_order_serialization(false);
    bvs.gap_length_serialization(false);
    bvs.set_compression_level(c_level);
    
    // serialization target for sparse vector
    bm::sparse_vector_serial_layout<sparse_vector_u32> sv_lay;
    bm::sparse_vector_serial_layout<rsc_sparse_vector_u32> csv_lay;
    
    bm::serializer<bm::bvector<> >::buffer sbuf; // resizable memory buffer

    // main loop to read sample vectors
    //
    bm::id64_t i;
    for (i = 0; true; ++i)
    {
        int ret = io_read_u32_coll(fin, vec);
        if (ret != 0)
            throw std::runtime_error("Error reading input file");
        unsigned min_delta, min_delta_cnt;

        {
            ret = validate_inp_vec(vec, min_delta, min_delta_cnt); // check if we have sorted unique set
            if (ret != 0)
                throw std::runtime_error("Input vector validation failed");
            if (!min_delta || !min_delta_cnt)
                throw std::runtime_error("Input vector validation delta error");
            
            min_delta_ints += min_delta_cnt;
        }
        total_ints += vec.size(); // remember the total size of the collection
        
        // serialize and save as a bit-vector size0:<BLOB0>, size1:<BLOB1>...N
        //
        bool is_low_card = false;
        if (!bv_out_fname.empty())
        {
            is_low_card = write_as_bvector(bv_file, vec, bvs, sbuf);
            if (is_low_card)
            {
                total_low_card_size += sbuf.size();
                ++total_low_card;
            }
        }
        
        // commented out experimental (and very slow code) to evaluate rank-select compression
        #if 0
        if (!sv_out_fname.empty())
        {
            //write_as_svector(sv_file, vec, min_delta, sv_lay);
            write_as_rsc_svector(sv_file, vec, min_delta, csv_lay);
            sv_size += sv_lay.size();

                /*
                rsc_sparse_vector_u32 csv; // compressed sparse vector
                csv.load_from(sv); // load rank-select-compacted (rsc) sparse vector
                
                BM_DECLARE_TEMP_BLOCK(tb)
                csv.optimize(tb);
                bm::sparse_vector_serial_layout<rsc_sparse_vector_u32> sv_lay;
                bm::sparse_vector_serialize(csv, sv_lay, tb);
                if (sv_lay.size() < sbuf.size())
                {
                    rsc_diff = sbuf.size() - sv_lay.size();
                    rsc_diff_size += rsc_diff;
                    sv_size += sv_lay.size();
                    sv_cnt++;
                }
                */
        }
        #endif
        
        
        std::streamsize fpos_curr = fin.tellg();
        if (fpos_curr == fsize)
            break;
        
        if (!is_silent)
        {
            cout << "\r" << i << "   " << fpos_curr << " / " << fsize
                 << " ( size=" << vec.size() << " ) " << (is_low_card ? " * " : "    ")
                 << " sv=" << sv_cnt << " rsc_diff=" << rsc_diff_size
                 << flush;
        }

    } // for i
    
    // print statistics about test set
    
    cout << endl;
    cout << "Total vectors=" << i << endl;
    cout << "  lo-card=" << total_low_card << endl;
    cout << "  lo-card size = " << total_low_card_size << endl;
    cout << "  SV cnt = " << sv_cnt << endl;
    cout << "  SV size = " << sv_size << endl;
    cout << "  RSC diff = " << rsc_diff_size << endl;
    cout << "Total ints=" << total_ints << endl;
    cout << "  min-deltas = " << min_delta_ints << endl;
    {
        double min_delta_ratio = double(min_delta_ints) / double(total_ints);
        cout << "  min delta ratio = " << std::setprecision(3) << min_delta_ratio << endl;
    }
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
    
    if (!sv_out_fname.empty())
    {
        sv_size = (bm::id64_t)sv_file.tellp();
        cout << "SV size = " << sv_size << endl;
        // calculate bits per int compression ratio corrected not to account
        // for size/length words prefixing the vectors
        double sv_bits_per_int = double(sv_size * 8ull - (i*sizeof(unsigned))) / double(total_ints);
        cout << "SV Bits per/int = " << std::setprecision(3) << sv_bits_per_int << endl;
        
        sv_file.close();
    }

}

// ------------------------------------------------------------------------

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
    
    sbuf.resize(len, false); // resize without content preservation
    bv_file.read((char*) sbuf.data(), std::streamsize(len));
    if (!bv_file.good())
        return -1;
    
    bm::deserialize(bv, sbuf.data());

    return 0;
}


/// read the input collection sequence and dump file, verify correctness
///
static
void verify_inv_dump_file(const std::string& fname,
                          const std::string& bv_in_fname)
{
    bm::id64_t total_ints = 0;


    bm::chrono_taker tt1(std::cout, "2. Verify collection", 1, &timing_map);


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
    std::streamsize fsize = 0;
    if (!bv_in_fname.empty())
    {
        bv_file.open(bv_in_fname, std::ios::in | std::ios::binary);
        if (!bv_file.good())
            throw std::runtime_error("Cannot open bvect dump file");
        fin.seekg(0, std::ios::end);
        fsize = fin.tellg();
        fin.seekg(0, std::ios::beg);
    }
    
    
    // initialize serializer
    //
    
    bm::serializer<bm::bvector<> >::buffer sbuf; // resizable memory buffer

    // main loop to read sample vectors
    //
    bm::id64_t i;
    for (i = 0; true; ++i)
    {
        int ret = io_read_u32_coll(fin, vec);
        if (ret != 0)
            throw std::runtime_error("Error reading input file");

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
        
        if (!is_silent)
        {
            cout << "\r" << fpos_curr << "/" << fsize
                 << " ( size=" << vec.size() << " )              "
                 << flush;
        }
    } // for i
    
    cout << endl;
    cout << "Verification complete." << endl;
    cout << "Total vectors=" << i << endl;
    cout << "Total ints=" << total_ints << endl;
}

/// read and decode the compressed dump file
///
static
void decode_test_dump_file(const std::string& bv_in_fname)
{
    bm::chrono_taker tt1(std::cout, "3. Decode collection", 1, &timing_map);

    std::ifstream bv_file;
    std::streamsize fsize;
    if (!bv_in_fname.empty())
    {
        bv_file.open(bv_in_fname, std::ios::in | std::ios::binary);
        if (!bv_file.good())
            throw std::runtime_error("Cannot open bvect dump file");
        bv_file.seekg(0, std::ios::end);
        fsize = bv_file.tellg();
        bv_file.seekg(0, std::ios::beg);
    }
    else
    {
        throw std::runtime_error("Cannot open bvect dump file");
    }
        

    bm::serializer<bm::bvector<> >::buffer sbuf; // resizable memory buffer

    // main loop to read sample vectors
    //
    bm::id64_t i;
    for (i = 0; true; ++i)
    {
        // serialize and save as a bit-vector size0:<BLOB0>, size1:<BLOB1>...N
        //
        if (!bv_in_fname.empty())
        {
            bm::bvector<> bv;
            read_bvector(bv_file, bv, sbuf);
        }

        std::streamsize fpos_curr = bv_file.tellg();
        if (fpos_curr == fsize)
            break;

        if (!is_silent)
        {
            cout << "\r" << fpos_curr << "/" << fsize
                 << flush;
        }
    } // for i

    cout << endl;
    cout << "Decode complete." << endl;
    cout << "Total vectors=" << i << endl;
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
                cout << "Compression" << endl;
                compress_inv_dump_file(u32_in_file, bv_out_file, sv_out_file);
            }
            else
            {
                if (is_verify)
                {
                    cout << "Verification." << endl;
                    verify_inv_dump_file(u32_in_file, bv_in_file);
                }
            }
        }
        if (is_decode)
        {
            cout << "Decode test." << endl;
            decode_test_dump_file(bv_in_file);
        }

                
        if (is_timing)  // print all collected timings
        {
            std::cout << std::endl << "Timings (ms):" << std::endl;
            bm::chrono_taker<>::print_duration_map(std::cout, timing_map);
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



