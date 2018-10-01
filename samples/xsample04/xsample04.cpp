/*
Copyright(c) 2018 Anatoliy Kuznetsov(anatoliy_kuznetsov at yahoo.com)

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

/** \example xsample04.cpp
 
  \sa bm::sparse_vector
  \sa bm::rsc_sparse_vector
*/

/*! \file xsample04.cpp
    \brief Example: DNA compression

*/

#include <iostream>
#include <sstream>
#include <chrono>
#include <regex>
#include <time.h>
#include <stdio.h>

#include <stdexcept>
#include <vector>
#include <chrono>
#include <map>
#include <utility>
#include <algorithm>

#include "bm.h"
#include "bmalgo.h"
#include "bmserial.h"
#include "bmrandom.h"
#include "bmsparsevec.h"
#include "bmsparsevec_compr.h"
#include "bmsparsevec_algo.h"
#include "bmsparsevec_serial.h"
#include "bmalgo_similarity.h"
#include "bmsparsevec_util.h"


#include "bmdbg.h"
#include "bmtimer.h"

static
void show_help()
{
    std::cerr
        << "BitMagic DNA Compression Sample (c) 2018" << std::endl
        << "-fa   file-name            -- input FASTA file" << std::endl
        << "-diag                      -- run diagnostics"  << std::endl
        << "-timing                    -- collect timings"  << std::endl
      ;
}




// Arguments
//
std::string  ifa_name;
bool         is_diag = false;
bool         is_timing = false;
bool         is_bench = false;

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
        
        if (arg == "-fa" || arg == "--fa")
        {
            if (i + 1 < argc)
            {
                ifa_name = argv[++i];
            }
            else
            {
                std::cerr << "Error: -fa requires file name" << std::endl;
                return 1;
            }
            continue;
        }

        if (arg == "-diag" || arg == "--diag" || arg == "-d" || arg == "--d")
            is_diag = true;
        if (arg == "-timing" || arg == "--timing" || arg == "-t" || arg == "--t")
            is_timing = true;
        if (arg == "-bench" || arg == "--bench" || arg == "-b" || arg == "--b")
            is_bench = true;

    } // for i
    return 0;
}


// Global types
//
typedef std::map<std::string, unsigned>                     freq_map;
typedef std::vector<std::pair<unsigned, std::string> >      dict_vect;

typedef bm::sparse_vector<unsigned, bm::bvector<> >         sparse_vector_u32;
typedef bm::rsc_sparse_vector<unsigned, sparse_vector_u32 > rsc_sparse_vector_u32;
typedef std::vector<std::pair<unsigned, unsigned> >         vector_pairs;

// ----------------------------------------------------------------------------

bm::chrono_taker::duration_map_type  timing_map;

// FASTA format parser
static
int load_FASTA(const std::string& fname, std::vector<char>& seq_vect)
{
    bm::chrono_taker tt1("1. parse FASTA", 1, &timing_map);
    
    seq_vect.resize(0);

    std::ifstream fin(fname.c_str(), std::ios::in);
    if (!fin.good())
        return -1;
    
    std::string line;
    
    for (unsigned i = 0; std::getline(fin, line); ++i)
    {
        if (line.empty() ||
            line.front() == '>')
            continue;
        
        for (std::string::iterator it = line.begin(); it != line.end(); ++it)
        {
            seq_vect.push_back(*it);
        }
    } // for
    
    std::cout << std::endl;
    std::cout << "FASTA sequence size=" << seq_vect.size() << std::endl;
    return 0;
}

static
void build_frequency_map_k2(const std::vector<char>& seq_vect, freq_map& fmap)
{
    char c1, c2;
    std::string k2;
    for (unsigned i = 0; i < seq_vect.size(); ++i)
    {
        c1 = seq_vect[i];
        ++i;
        c2 = (i < seq_vect.size()) ? seq_vect[i] : 'A';
        k2 = c1;
        k2.push_back(c2);
        
        fmap[k2]++;
    }
}

static
void build_frequency_map_k3(const std::vector<char>& seq_vect, freq_map& fmap)
{
    char c1, c2, c3;
    std::string k3;
    for (unsigned i = 0; i < seq_vect.size(); ++i)
    {
        c1 = seq_vect[i];
        ++i;
        c2 = (i < seq_vect.size()) ? seq_vect[i] : 'A';
        ++i;
        c3 = (i < seq_vect.size()) ? seq_vect[i] : 'G';
        
        k3 = c1; k3.push_back(c2); k3.push_back(c3);
        
        fmap[k3]++;
    }
}


static
void build_dict(const freq_map& fmap, dict_vect& dict)
{
    dict.resize(0);
    
    freq_map::const_iterator it = fmap.begin();
    freq_map::const_iterator it_end = fmap.end();
    for (; it != it_end; ++it)
        dict.push_back(std::pair<unsigned, std::string>(it->second, it->first));
    std::sort(dict.begin(), dict.end(),
            [](std::pair<unsigned, std::string> a, std::pair<unsigned, std::string> b) {
        return a > b;
    });
}


static
void print_map(const freq_map& fmap, const dict_vect& dict)
{
    dict_vect::const_iterator it = dict.begin();
    dict_vect::const_iterator it_end = dict.end();
    
    std::cout << "Dictionary size = " << dict.size() << std::endl;
    for (unsigned i = 0; it != it_end; ++it, ++i)
    {
        const std::string& key = it->second;
        unsigned cnt = it->first;
        std::cout << i << ": " << key << "=" << cnt << std::endl;
    }
    std::cout << std::endl;
}

inline
unsigned find_id(const dict_vect& dict, const std::string& key)
{
    dict_vect::const_iterator it = dict.begin();
    dict_vect::const_iterator it_end = dict.end();
    for (unsigned i = 0; it != it_end; ++it, ++i)
    {
        const std::string& dkey = it->second;
        if (key == dkey)
            return i;
    }
    throw std::runtime_error("Item not found!");
}


static
void build_sv_k2(const std::vector<char>& seq_vect,
                 const freq_map&          fmap,
                 const dict_vect&         dict,
                 sparse_vector_u32&       sv_k2)
{
    char c1, c2;
    std::string k2;
    unsigned sv_cnt = 0;
    unsigned prev_id = 0;
    for (unsigned i = 0; i < seq_vect.size(); ++i, ++sv_cnt)
    {
        c1 = seq_vect[i];
        ++i;
        c2 = (i < seq_vect.size()) ? seq_vect[i] : 'A';

        k2 = c1;
        k2.push_back(c2);
        
        unsigned dict_id = find_id(dict, k2);
        if (dict_id)
        {
            sv_k2.set(sv_cnt, dict_id);
        }
    }
}

static
void build_sv_k3(const std::vector<char>& seq_vect,
                 const freq_map&          fmap,
                 const dict_vect&         dict,
                 sparse_vector_u32&       sv_k3)
{
    char c1, c2, c3;
    std::string k3;
    unsigned sv_cnt = 0;
    unsigned prev_id = 0;
    for (unsigned i = 0; i < seq_vect.size(); ++i, ++sv_cnt)
    {
        c1 = seq_vect[i];
        ++i;
        c2 = (i < seq_vect.size()) ? seq_vect[i] : 'A';
        ++i;
        c3 = (i < seq_vect.size()) ? seq_vect[i] : 'G';
        k3 = c1;
        k3.push_back(c2);
        k3.push_back(c3);
        
        unsigned dict_id = find_id(dict, k3);
        if (dict_id)
            sv_k3.set(sv_cnt, dict_id);
    }
}



int main(int argc, char *argv[])
{
    if (argc < 3)
    {
        show_help();
        return 1;
    }
    
    std::vector<char> seq_vect;
    freq_map          fmap_k2;
    dict_vect         dict_k2;
    freq_map          fmap_k3;
    dict_vect         dict_k3;
    
    sparse_vector_u32  sv_k2(bm::use_null);
    rsc_sparse_vector_u32 csv_k2;

    sparse_vector_u32  sv_k3(bm::use_null);
    rsc_sparse_vector_u32 csv_k3;

    try
    {
        auto ret = parse_args(argc, argv);
        if (ret != 0)
            return ret;

        if (!ifa_name.empty())
        {
            auto res = load_FASTA(ifa_name, seq_vect);
            if (res != 0)
            {
                return res;
            }
        }
        if (!seq_vect.empty())
        {
            build_frequency_map_k2(seq_vect, fmap_k2);
            build_dict(fmap_k2, dict_k2);

            build_frequency_map_k3(seq_vect, fmap_k3);
            build_dict(fmap_k3, dict_k3);

            
            print_map(fmap_k2, dict_k2);
            print_map(fmap_k3, dict_k3);
            
            build_sv_k2(seq_vect, fmap_k2, dict_k2, sv_k2);
            csv_k2.load_from(sv_k2);
            
            sv_k2.optimize();
            csv_k2.optimize();

            build_sv_k3(seq_vect, fmap_k3, dict_k3, sv_k3);
            csv_k3.load_from(sv_k3);

            sv_k3.optimize();
            csv_k3.optimize();

            file_save_svector(csv_k3, "dna_k3.csv");
        }

        if (is_diag) // run set of benchmarks
        {
            if (!sv_k2.empty())
            {
                std::cout << std::endl
                          << "sparse vector K2 statistics:"
                          << std::endl;
                bm::print_svector_stat(sv_k2, true);
            }
            if (!csv_k2.empty())
            {
                std::cout << std::endl
                          << "RSC sparse vector K2 statistics:"
                          << std::endl;
                bm::print_svector_stat(csv_k2, true);
            }
            if (!sv_k3.empty())
            {
                std::cout << std::endl
                          << "sparse vector K3 statistics:"
                          << std::endl;
                bm::print_svector_stat(sv_k3, true);
            }
            if (!csv_k3.empty())
            {
                std::cout << std::endl
                          << "RSC sparse vector K3 statistics:"
                          << std::endl;
                bm::print_svector_stat(csv_k3, true);
            }


        }


        if (is_bench) // run set of benchmarks
        {
        }

        if (is_timing)  // print all collected timings
        {
            std::cout << std::endl << "Performance:" << std::endl;
            bm::chrono_taker::print_duration_map(timing_map, bm::chrono_taker::ct_ops_per_sec);
        }
    }
    catch (std::exception& ex)
    {
        std::cerr << "Error:" << ex.what() << std::endl;
        return 1;
    }

    return 0;
}

