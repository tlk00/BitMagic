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

/** \example xsample05.cpp
*/

/*! \file xsample05.cpp
    \brief Example: Example on how to use bit-transposed string sparse vector

    Illustrates how to build a sparse vector, serialize it to disk,
    load back and do search or binary hybrid search.
 
  \sa bm::str_sparse_vector<>
   
*/

#include <iostream>
#include <sstream>
#include <chrono>
#include <regex>
#include <time.h>
#include <stdio.h>


#include <vector>
#include <map>
#include <chrono>
#include <map>
#include <utility>
#include <algorithm>
#include <random>

using namespace std;

//#define BMAVX2OPT

#include "bm.h"
#include "bmalgo.h"
#include "bmserial.h"
#include "bmrandom.h"
#include "bmstrsparsevec.h"
#include "bmsparsevec_algo.h"
#include "bmsparsevec_serial.h"
#include "bmalgo_similarity.h"
#include "bmsparsevec_util.h"


#include "bmdbg.h"
#include "bmtimer.h"
#include "bmundef.h" /* clear the pre-proc defines from BM */

/// Print help
static
void show_help()
{
    std::cerr
        << "BitMagic Dictionary Search Sample (c) 2018" << std::endl
        << "-idict  file-name            -- input set file to parse" << std::endl
        << "-svout  spase vector output  -- sparse vector name to save" << std::endl
        << "-svin   sparse vector input  -- sparse vector file name to load " << std::endl
        << "-remap                       -- re-mapping of string characters " << std::endl
        << "-xor                         -- use XOR compression filtering" << std::endl
        << "-diag                        -- run diagnostics"                  << std::endl
        << "-bench                       -- run benchmarks"                   << std::endl
        << "-timing                      -- collect timings"                  << std::endl
      ;
}


// Arguments
//
std::string  sv_out_name;
std::string  sv_in_name;
std::string  i_name;
bool         is_diag = false;
bool         is_timing = false;
bool         is_bench = false;
bool         is_remap = false;
bool         is_xor = false;

// Parse command line arguments
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

        if (arg == "-idict" || arg == "--idict" )
        {
            if (i + 1 < argc)
            {
                i_name = argv[++i];
            }
            else
            {
                std::cerr << "Error: -idict requires file name" << std::endl;
                return 1;
            }
            continue;
        }

        if (arg == "-remap" || arg == "--remap" || arg == "-r" || arg == "--r")
        {
            is_remap = true;
            continue;
        }
        if (arg == "-xor" || arg == "--xor" || arg == "-x" || arg == "--x")
        {
            is_xor = true;
            continue;
        }
        if (arg == "-diag" || arg == "--diag" || arg == "-d" || arg == "--d")
        {
            is_diag = true;
            continue;
        }
        if (arg == "-timing" || arg == "--timing" || arg == "-t" || arg == "--t")
        {
            is_timing = true;
            continue;
        }
        if (arg == "-bench" || arg == "--bench" || arg == "-b" || arg == "--b")
        {
            is_bench = true;
            continue;
        }
        
        std::cerr << "Error: unknown argument: " << arg << std::endl;
        return 1;


    } // for i
    return 0;
}


// Global types
//
typedef bm::str_sparse_vector<char, bm::bvector<>, 64>  str_sparse_vect;
typedef vector<string>  string_vector;



bm::chrono_taker<>::duration_map_type  timing_map;

/// Parse the input file and extract dictionary values.
/// 
static
int load_dict_report(const std::string& fname, string_vector& str_vec)
{
    bm::chrono_taker tt1(cout, "1. parse input data ", 1, &timing_map);

    std::ifstream fin(fname.c_str(), std::ios::in);
    if (!fin.good())
        return -1;

    std::string line;
    std::regex reg("[|]");
    std::sregex_token_iterator it_end;
    
    string trim_chars("\" ");

    for (unsigned i = 0; std::getline(fin, line); ++i)
    {
        if (line.empty() || !isdigit(line.front()))
            continue;

        // regex based tokenizer
        std::sregex_token_iterator it(line.begin(), line.end(), reg, -1);
        std::vector<std::string> line_vec(it, it_end);
        if (line_vec.empty())
            continue;
        try
        {
            // trip the extra chars
            string& col13 = line_vec.at(13);
            col13.erase(0, col13.find_first_not_of(trim_chars));
            col13.erase(col13.find_last_not_of(trim_chars) + 1);

            if (!col13.empty())
                str_vec.emplace_back(col13);
        }
        catch (std::exception&)
        {
            // just ignore (not ideal, but ok for a sketch ETL)
        }
        if (i % 10000 == 0)
        {
            cout << "\rReading input file: " << i << flush;
        }
    } // for
    cout << endl;
    return 0;
}

/// Compare STL vector with bit-transposed container to check correctness
///
static
void check_sparse(const str_sparse_vect& str_sv, const string_vector& str_vec)
{
    if (str_vec.size() != str_sv.size())
        throw runtime_error("Error. size() comparison failed!");
    string s;
    for (str_sparse_vect::size_type i = 0; i < str_sv.size(); ++i)
    {
        str_sv.get(i, s);
        const string& s_control = str_vec[i];
        if (s != s_control)
        {
            cout << "idx=" << i <<  s << "!=" << s_control << endl;
            throw runtime_error("Error. element comparison failed!");
        }
    } // for
    std::cout << "Check ok. Dictionary size = " << str_sv.size() << std:: endl;
}

const unsigned benchmark_max = 50000;  // benchmark sampling size

/// Sample a few random terms out of collection
static
void pick_benchmark_set(string_vector& bench_vec, string_vector& bench_vec_not_found, const string_vector& str_vec)
{
    std::random_device rand_dev;
    std::mt19937 gen(rand_dev()); // mersenne_twister_engine 
    std::uniform_int_distribution<unsigned> rand_dis(0, unsigned(str_vec.size()-1)); // generate uniform numebrs for [1, vector_max]

    bm::bvector<> bv;
    
    bench_vec.resize(0);
    for (unsigned i = 0; i < benchmark_max; ++i)
    {
        unsigned idx;
        while (true)
        {
            idx = unsigned(rand_dis(gen));
            if (bv.test(idx))  // make sure benchmark example is not repeated
                continue;
            if (idx < str_vec.size())
                bench_vec.push_back(str_vec[idx]);
            break;
        }
        bv.set(idx); // mark as set
        
        // generate not-found case by modifying a letter in an existing sample
        {
            string str_nf = str_vec[idx];
            string::reverse_iterator rit = str_nf.rbegin();
            string::reverse_iterator rit_end = str_nf.rend();
            for (; rit != rit_end; ++rit)
            {
                char ch = *rit;
                int a = rand() % 26 + int('A'); // generate random letter
                ch = char(a);
                *rit = ch;
                auto it = std::lower_bound(str_vec.begin(), str_vec.end(), str_nf);
                if (it == str_vec.end() || *it != str_nf) // check if not found
                {
                    bench_vec_not_found.push_back(str_nf);
                    break;
                }
            } // for rit
        }
         
    } // for
    cout << endl;
}

static
void run_benchmark(const str_sparse_vect& str_sv, const string_vector& str_vec)
{
    string_vector bench_vec;
    string_vector bench_vec_not_found; // vector for impossible dictionary items

    pick_benchmark_set(bench_vec, bench_vec_not_found, str_vec);
    
    bm::bvector<> bv1, bv2, bv3, bv4;

    cout << "Picked " << bench_vec.size() << " / " 
         << bench_vec_not_found.size() << " samples. Running benchmarks." 
         << endl;
    
    unsigned bench_size = unsigned(bench_vec.size());
    {
        {
            bm::chrono_taker tt1(cout, "3.  std::lower_bound() search", bench_size, &timing_map);
            for (const string& term : bench_vec)
            {
                auto it = std::lower_bound(str_vec.begin(), str_vec.end(), term);
                if (it != str_vec.end())
                {
                    string_vector::size_type idx =
                        string_vector::size_type(std::distance(str_vec.begin(), it));
                    bv1.set(unsigned(idx));
                }
            } // for
        }
        {
            bm::chrono_taker tt2(cout, "3a. std::lower_bound() search (empty)", bench_size, &timing_map);
            for (const string& term : bench_vec_not_found)
            {
                auto p = std::lower_bound(str_vec.begin(), str_vec.end(), term);
                (void) p;
            } // for
        }
    }
    
    {
        // construct std::map<>  (RB-tree)
        std::map<string, unsigned> str_map;
        for (string_vector::size_type i = 0; i < str_vec.size(); ++i)
        {
            const string& s = str_vec[i];
            str_map[s] = unsigned(i);
        } // for
        {
            bm::chrono_taker tt1(cout, "4.  std::map<> search", bench_size, &timing_map);
            for (const string& term : bench_vec)
            {
                auto it = str_map.find(term);
                if (it != str_map.end())
                {
                    bv2.set(unsigned(it->second));
                }
            } // for
        }
        {
            bm::chrono_taker tt2(cout, "4a. std::map<> search (empty)", bench_size, &timing_map);
            for (const string& term : bench_vec_not_found)
            {
                auto it = str_map.find(term);
                if (it != str_map.end())
                {
                    cerr << "empty search returned value..." << endl;
                }
            } // for
        }
    }

    {
        bm::sparse_vector_scanner<str_sparse_vect> scanner;
        {
            bm::chrono_taker tt1(cout, "5.  bm::sparse_vector_scanner<> search", bench_size, &timing_map);
            for (const string& term : bench_vec)
            {
                unsigned pos;
                bool found = scanner.find_eq_str(str_sv, term.c_str(), pos);
                if (found)
                {
                    bv3.set(pos);
                }
            } // for
        }
        {
            bm::chrono_taker tt1(cout, "5a. bm::sparse_vector_scanner<> search (empty)", bench_size, &timing_map);
            for (const string& term : bench_vec_not_found)
            {
                unsigned pos;
                bool found = scanner.find_eq_str(str_sv, term.c_str(), pos);
                if (found)
                {
                    cerr << "scanner empty search returned value..." << endl;
                }
            } // for
        }
    }

    {
        bm::sparse_vector_scanner<str_sparse_vect> scanner;
        scanner.bind(str_sv, true); // attach SV as permanent search parameter to share cached values

        {
            bm::chrono_taker tt1(cout, "6.  bm::sparse_vector_scanner<> binary search", bench_size, &timing_map);
            for (const string& term : bench_vec)
            {
                unsigned pos;
                bool found = scanner.bfind_eq_str(term.c_str(), pos);
                if (found)
                {
                    bv4.set(pos);
                }
            } // for
        }
        {
            bm::chrono_taker tt2(cout, "6a. bm::sparse_vector_scanner<> binary search (empty)", bench_size, &timing_map);
            for (const string& term : bench_vec_not_found)
            {
                unsigned pos;
                bool found = scanner.bfind_eq_str(term.c_str(), pos);
                if (found)
                {
                    cerr << "scanner empty search returned value..." << endl;
                }
            } // for
        }
    }

    // various integrity checks
    //
    int cmp = bv1.compare(bv2);
    if (cmp != 0)
        throw runtime_error("Error. RB-search mismatch!");
    cmp = bv1.compare(bv3);
    if (cmp != 0)
        throw runtime_error("Error. scanner mismatch!");

    cmp = bv1.compare(bv4);
    if (cmp != 0)
        throw runtime_error("Error. binary scanner mismatch!");

    if (bv1.count() != bench_size)
        throw runtime_error("Error. Search result missing elements!");


}


int main(int argc, char *argv[])
{
    if (argc < 3)
    {
        show_help();
        return 1;
    }

    string_vector      str_vec;  // dictionary vector (STL)
    str_sparse_vect    str_sv;   // bit-transposed sparse vector

    try
    {
        auto ret = parse_args(argc, argv);
        if (ret != 0)
        {
            return ret;
        }
        if (!i_name.empty())
        {
            auto res = load_dict_report(i_name, str_vec);
            if (res != 0)
            {
                return res;
            }
            cout << "Loaded " << str_vec.size() << " dictionary names." << endl;
            
            std::sort(str_vec.begin(), str_vec.end());
        }
        
        if (str_vec.size()) // load the sparse vector
        {
            bm::chrono_taker tt1(cout, "2. build sparse vector", 1, &timing_map);
            {
                // use insert iterator to load vector (faster than push-back)
                //
                auto bi = str_sv.get_back_inserter();
                for (const string& term : str_vec)
                    bi = term;
                bi.flush();
            }
            
            // build remapped (dense) succinct vector
            // (this should be final), no more edits in this form
            if (is_remap)
            {
                str_sparse_vect    str_sv_remap;
                str_sv_remap.remap_from(str_sv);
                str_sv.swap(str_sv_remap);
            }
            
            BM_DECLARE_TEMP_BLOCK(tb)
            str_sv.optimize(tb); // memory optimization after load
        }

        if (!sv_in_name.empty())
        {
            {
                bm::chrono_taker tt1(cout, "8. Load sparse vector", 1, &timing_map);
                file_load_svector(str_sv, sv_in_name);
            }
            if (str_sv.empty())
            {
                cout << "Input vector empty!" << endl;
                exit(1);
            }
            if (str_vec.empty())
            {
                for (str_sparse_vect::size_type i = 0; i < str_sv.size(); ++i)
                {
                    string s;
                    str_sv.get(i, s);
                    str_vec.emplace_back(std::move(s));
                } // for
            }
        }

        // save SV vector to file
        if (!sv_out_name.empty() && !str_sv.empty())
        {
            bm::chrono_taker tt1(cout, "7. Save sparse vector", 1, &timing_map);
            file_save_svector(str_sv, sv_out_name, 0, is_xor);

            str_sparse_vect    str_sv_control;
            file_load_svector(str_sv_control, sv_out_name);
            bool eq = str_sv.equal(str_sv_control);
            if (!eq)
            {
                cerr << "Serialization control failed" << endl;
                assert(0); exit(1);
            }
        }
        

        
        if (is_diag)
        {
            if (!str_sv.empty())
            {
                print_svector_stat(cout,str_sv, false);
            }
            
            if (str_vec.size())
            {
                size_t total_size = 0;
                for (const string& term : str_vec)
                {
                    total_size += term.size();
                }
                cout << "String dictionary size: "
                     << total_size / 1024 << "KB (" << total_size / (1024*1024) << "MB)"
                     << endl;
            }
            
            if (str_sv.size() && str_vec.size())
            {
                cout << "Run full comparison check..." << endl;
                check_sparse(str_sv, str_vec); // run a paranoiya check
                cout << "Ok" << endl;
            }
        }

        if (is_bench) // run set of benchmarks
        {
            run_benchmark(str_sv, str_vec);
        }

        if (is_timing)  // print all collected timings
        {
            std::cout << std::endl << "Performance:" << std::endl;
            bm::chrono_taker<>::print_duration_map(cout, timing_map, bm::chrono_taker<>::ct_time);
        }
    }
    catch (std::exception& ex)
    {
        std::cerr << "Error:" << ex.what() << std::endl;
        return 1;
    }

    return 0;
}

