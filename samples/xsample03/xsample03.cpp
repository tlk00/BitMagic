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

/** \example xsample03.cpp
   Seach for human mutation (SNP) in within chr1.
   Benchmark comaprison of different methods
 
  \sa bm::sparse_vector
  \sa bm::rsc_sparse_vector
  \sa bm::sparse_vector_scanner
*/

/*! \file xsample03.cpp
    \brief Example: SNP search in human genome
 
   Brief description of used method:
   1. Parse SNP chromosome report and extract information about SNP number and
      location in the chromosome
   2. Use this information to build bit-transposed sparse_vector<>
      where vector position matches chromosome position and SNP ids (aka rsid)
      is kept as a bit-transposed matrix
   3. Build rank-select compressed sparse vector, dropping all NULL columns
      (this data format is pretty sparse, since number of SNPs is significantly
       less than number of chromosome bases (1:5 or less)
       Use memory report to understand memory footprint for each form of storage
   4. Run benchmarks searching for 500 randomly selected SNPs using
      - bm::sparse_vector<>
      - bm::rsc_sparse_vector<>
      - std::vector<pair<unsigned, unsigned> >
 
   This example should be useful for construction of compressed columnar
   tables with parallel search capabilities.
 
*/

#include <iostream>
#include <sstream>
#include <chrono>
#include <regex>
#include <time.h>
#include <stdio.h>
#include <vector>
#include <chrono>
#include <map>
#include <utility>

using namespace std;

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
#include "bmundef.h" /* clear the pre-proc defines from BM */

static
void show_help()
{
    std::cerr
        << "BitMagic SNP Search Sample Utility (c) 2018" << std::endl
        << "-isnp   file-name            -- input set file (SNP FASTA) to parse" << std::endl
        << "-svout  spase vector output  -- sparse vector name to save" << std::endl
        << "-rscout rs-compressed spase vector output  -- name to save" << std::endl
        << "-svin   sparse vector input   -- sparse vector file name to load " << std::endl
        << "-rscin  rs-compressed sparse vector input   -- file name to load " << std::endl
        << "-diag                        -- run diagnostics"                  << std::endl
        << "-timing                      -- collect timings"                  << std::endl
      ;
}




// Arguments
//
std::string  sv_out_name;
std::string  rsc_out_name;
std::string  sv_in_name;
std::string  rsc_in_name;
std::string  isnp_name;
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
        
        if (arg == "-svout" || arg == "--svout")
        {
            if (i + 1 < argc)
            {
                sv_out_name = argv[++i];
            }
            else
            {
                std::cerr << "Error: -svout requires file name" << std::endl;
                return 1;
            }
            continue;
        }
        if (arg == "-rscout" || arg == "--rscout")
        {
            if (i + 1 < argc)
            {
                rsc_out_name = argv[++i];
            }
            else
            {
                std::cerr << "Error: -rscout requires file name" << std::endl;
                return 1;
            }
            continue;
        }

        if (arg == "-svin" || arg == "--svin")
        {
            if (i + 1 < argc)
            {
                sv_in_name = argv[++i];
            }
            else
            {
                std::cerr << "Error: -svin requires file name" << std::endl;
                return 1;
            }
            continue;
        }

        if (arg == "-rscin" || arg == "--rscin")
        {
            if (i + 1 < argc)
            {
                rsc_in_name = argv[++i];
            }
            else
            {
                std::cerr << "Error: -rscin requires file name" << std::endl;
                return 1;
            }
            continue;
        }

        if (arg == "-isnp" || arg == "--isnp" || arg == "-snp" || arg == "--snp")
        {
            if (i + 1 < argc)
            {
                isnp_name = argv[++i];
            }
            else
            {
                std::cerr << "Error: -isnp requires file name" << std::endl;
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
typedef bm::sparse_vector<unsigned, bm::bvector<> >         sparse_vector_u32;
typedef bm::rsc_sparse_vector<unsigned, sparse_vector_u32 > rsc_sparse_vector_u32;
typedef std::vector<std::pair<unsigned, unsigned> >         vector_pairs;

// ----------------------------------------------------------------------------

bm::chrono_taker<>::duration_map_type  timing_map;

// SNP report format parser (extraction and transformation)
// Parser extracts SNP id (rsid) and positio on genome to build
// sparse vector where index (position in vector) metches position on the
// chromosome 1.
//
static
int load_snp_report(const std::string& fname, sparse_vector_u32& sv)
{
    bm::chrono_taker tt1(cout, "1. parse input SNP chr report", 1, &timing_map);

    std::ifstream fin(fname.c_str(), std::ios::in);
    if (!fin.good())
        return -1;

    unsigned rs_id, rs_pos;
    size_t idx;

    std::string line;
    std::string delim = " \t";

    std::regex reg("\\s+");
    std::sregex_token_iterator it_end;

    bm::bvector<> bv_rs; 
    bv_rs.init();

    unsigned rs_cnt = 0;
    for (unsigned i = 0; std::getline(fin, line); ++i)
    {
        if (line.empty() ||
            !isdigit(line.front())
            )
            continue;

        // regex based tokenizer
        std::sregex_token_iterator it(line.begin(), line.end(), reg, -1);
        std::vector<std::string> line_vec(it, it_end);
        if (line_vec.empty())
            continue; 
        
        // parse columns of interest
        try
        {
            rs_id = unsigned(std::stoul(line_vec.at(0), &idx));
            
            if (bv_rs.test(rs_id))
            {
                continue;
            }
            rs_pos = unsigned(std::stoul(line_vec.at(11), &idx));

            bv_rs.set_bit_no_check(rs_id);
            sv.set(rs_pos, rs_id);

            ++rs_cnt;
        }
        catch (std::exception& /*ex*/)
        {
            continue; // detailed disgnostics commented out
            // error detected, because some columns are sometimes missing
            // just ignore it
            //
            /*
            std::cerr << ex.what() << "; ";
            std::cerr << "rs=" << line_vec.at(0) << " pos=" << line_vec.at(11) << std::endl;
            continue;
            */
        }
        if (rs_cnt % (4 * 1024) == 0)
            std::cout << "\r" << rs_cnt << " / " << i; // PROGRESS report
    } // for i

    std::cout << std::endl;
    std::cout << "SNP count=" << rs_cnt << std::endl;
    return 0;
}

// Generate random subset of random values from a sparse vector
//
static
void generate_random_subset(const sparse_vector_u32&  sv, std::vector<unsigned>& vect, unsigned count)
{
    const sparse_vector_u32::bvector_type* bv_null = sv.get_null_bvector();

    bm::random_subset<bm::bvector<> > rand_sampler;
    bm::bvector<> bv_sample;
    rand_sampler.sample(bv_sample, *bv_null, count);

    bm::bvector<>::enumerator en = bv_sample.first();
    for (; en.valid(); ++en)
    {
        unsigned idx = *en;
        unsigned v = sv[idx];
        vect.push_back(v);
    }
}

// build std::vector of pairs (position to rs)
//
static
void build_vector_pairs(const sparse_vector_u32& sv, vector_pairs& vp)
{
    sparse_vector_u32::const_iterator it = sv.begin();
    sparse_vector_u32::const_iterator it_end = sv.end();
    
    for (; it != it_end; ++it)
    {
        if (!it.is_null())
        {
            std::pair<unsigned, unsigned> pos2rs = std::make_pair(it.pos(), it.value());
            vp.push_back(pos2rs);
        }
    }
}

// O(N) -- O(N/2) linear search in vector of pairs (position - rsid)
//
static
bool search_vector_pairs(const vector_pairs& vp, unsigned rs_id, unsigned& pos)
{
    for (unsigned i = 0; i < vp.size(); ++i)
    {
        if (vp[i].second == rs_id)
        {
            pos = vp[i].first;
            return true;
        }
    }
    return false;
}

// SNP search benchmark
// Search for SNPs using different data structures (Bitmagic and STL)
//
// An extra step is verification, to make sure all methods are consistent
//
static
void run_benchmark(const sparse_vector_u32& sv, const rsc_sparse_vector_u32& csv)
{
    const unsigned rs_sample_count = 2000;

    std::vector<unsigned> rs_vect;
    generate_random_subset(sv, rs_vect, rs_sample_count);
    if (rs_vect.empty())
    {
        std::cerr << "Benchmark subset empty!" << std::endl;
        return;
    }
    
    // build traditional sparse vector
    vector_pairs vp;
    build_vector_pairs(sv, vp);
    
    // search result bit-vectors
    //
    bm::bvector<> bv_found1;
    bm::bvector<> bv_found2;
    bm::bvector<> bv_found3;

    bv_found1.init(); bv_found2.init(); bv_found3.init();// pre-initialize vectors

    if (!sv.empty())
    {
        bm::chrono_taker tt1(cout, "09. rs search (sv)", unsigned(rs_vect.size()), &timing_map);
        
        // scanner class
        // it's important to keep scanner class outside the loop to avoid
        // unnecessary re-allocs and construction costs
        //

        bm::sparse_vector_scanner<sparse_vector_u32> scanner;

        for (unsigned i = 0; i < rs_vect.size(); ++i)
        {
            unsigned rs_id = rs_vect[i];
            unsigned rs_pos;
            bool found = scanner.find_eq(sv, rs_id, rs_pos);

            if (found)
            {
                bv_found1.set_bit_no_check(rs_pos);
            }
            else
            {
                std::cout << "Error: rs_id = " << rs_id << " not found!" << std::endl;
            }
        } // for
    }

    if (!csv.empty())
    {
        bm::chrono_taker tt1(cout, "10. rs search (rsc_sv)", unsigned(rs_vect.size()), &timing_map);

        bm::sparse_vector_scanner<rsc_sparse_vector_u32> scanner; // scanner class

        for (unsigned i = 0; i < rs_vect.size(); ++i)
        {
            unsigned rs_id = rs_vect[i];
            unsigned rs_pos;
            bool found = scanner.find_eq(csv, rs_id, rs_pos);

            if (found)
            {
                bv_found2.set_bit_no_check(rs_pos);
            }
            else
            {
                std::cout << "rs_id = " << rs_id << " not found!" << std::endl;
            }
        } // for
    }

    if (vp.size())
    {
        bm::chrono_taker tt1(cout, "11. rs search (std::vector<>)", unsigned(rs_vect.size()), &timing_map);

        for (unsigned i = 0; i < rs_vect.size(); ++i)
        {
            unsigned rs_id = rs_vect[i];
            unsigned rs_pos;
            bool found = search_vector_pairs(vp, rs_id, rs_pos);

            if (found)
            {
                bv_found3.set_bit_no_check(rs_pos);
            }
            else
            {
                std::cout << "rs_id = " << rs_id << " not found!" << std::endl;
            }
        } // for
    }

    // compare results from various methods (check integrity)
    int res = bv_found1.compare(bv_found2);
    if (res != 0)
    {
        std::cerr << "Error: search discrepancy (sparse search) detected!" << std::endl;
    }
    res = bv_found1.compare(bv_found3);
    if (res != 0)
    {
        std::cerr << "Error: search discrepancy (std::vector<>) detected!" << std::endl;
    }

}


int main(int argc, char *argv[])
{
    if (argc < 3)
    {
        show_help();
        return 1;
    }

    sparse_vector_u32  sv(bm::use_null);
    rsc_sparse_vector_u32 csv;

    try
    {
        auto ret = parse_args(argc, argv);
        if (ret != 0)
            return ret;

        if (!isnp_name.empty())
        {
            auto res = load_snp_report(isnp_name, sv);
            if (res != 0)
            {
                return res;
            }
        }
        if (!sv_in_name.empty())
        {
            bm::chrono_taker tt1(cout, "02. Load sparse vector", 1, &timing_map);
            file_load_svector(sv, sv_in_name);
        }
        
        // load rs-compressed sparse vector
        if (!rsc_in_name.empty())
        {
            bm::chrono_taker tt1(cout, "03. Load rsc sparse vector", 1, &timing_map);
            file_load_svector(csv, rsc_in_name);
        }
        
        // if rs-compressed vector is not available - build it on the fly
        if (csv.empty() && !sv.empty())
        {
            sparse_vector_u32  sv2(bm::use_null);
            {
                bm::chrono_taker tt1(cout, "04. rs compress sparse vector", 1, &timing_map);
                csv.load_from(sv);
            }
            {
                bm::chrono_taker tt1(cout, "05. rs de-compress sparse vector", 1, &timing_map);
                csv.load_to(sv2);
            }
            
            if (!sv.equal(sv2)) // diagnostics check (just in case)
            {
                std::cerr << "Error:  rs-compressed vector check failed!" << std::endl;
                return 1;
            }
        }
        
        // save SV vector
        if (!sv_out_name.empty() && !sv.empty())
        {
            bm::chrono_taker tt1(cout, "07. Save sparse vector", 1, &timing_map);
            sv.optimize();
            file_save_svector(sv, sv_out_name);
        }

        // save RS sparse vector
        if (!rsc_out_name.empty() && !csv.empty())
        {
            bm::chrono_taker tt1(cout, "08. Save RS sparse vector", 1, &timing_map);
            csv.optimize();
            file_save_svector(csv, rsc_out_name);
        }
        
        if (is_diag) // print memory diagnostics
        {
            if (!sv.empty())
            {
                std::cout << std::endl
                          << "sparse vector statistics:"
                          << std::endl;
                bm::print_svector_stat(std::cout, sv, true);
            }
            if (!csv.empty())
            {
                std::cout << std::endl
                          << "RS compressed sparse vector statistics:"
                          << std::endl;
                bm::print_svector_stat(std::cout, csv, true);
            }
        }

        if (is_bench) // run set of benchmarks
        {
            run_benchmark(sv, csv);
        }

        if (is_timing)  // print all collected timings
        {
            std::cout << std::endl << "Performance:" << std::endl;
            bm::chrono_taker<>::print_duration_map(cout, timing_map, bm::chrono_taker<>::ct_ops_per_sec);
        }
    }
    catch (std::exception& ex)
    {
        std::cerr << "Error:" << ex.what() << std::endl;
        return 1;
    }

    return 0;
}

