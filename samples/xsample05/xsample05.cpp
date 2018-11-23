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
    \brief Example:
 
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

using namespace std;

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

/// Print help
static
void show_help()
{
    std::cerr
        << "BitMagic Dictionary Search Sample (c) 2018" << std::endl
        << "-idict  file-name            -- input set file to parse" << std::endl
        << "-svout  spase vector output  -- sparse vector name to save" << std::endl
        << "-svin   sparse vector input  -- sparse vector file name to load " << std::endl
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
typedef bm::str_sparse_vector<char, bm::bvector<>, 64>  str_sparse_vect;
typedef vector<string>  string_vector;



bm::chrono_taker::duration_map_type  timing_map;

/// Parse the input file and extract dictionary values.
/// 
static
int load_dict_report(const std::string& fname, string_vector& str_vec)
{
    bm::chrono_taker tt1("1. parse input data ", 1, &timing_map);

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

/// Compare STL vector with bit-transposed container to check correcness
///
static
void check_sparse(const str_sparse_vect& str_sv, const string_vector& str_vec)
{
    if (str_vec.size() != str_sv.size())
        throw runtime_error("Error. size() comparison failed!");
    for (str_sparse_vect::size_type i = 0; i < str_sv.size(); ++i)
    {
        string s;
        str_sv.get(i, s);
        const string& s_control = str_vec[i];
        if (s != s_control)
            throw runtime_error("Error. element comparison failed!");
    } // for
}

const unsigned benchmark_max = 1500;  // benchmark sampling size

/// Sample a few random terms out of collection
static
void pick_benchmark_set(string_vector& bench_vec, const string_vector& str_vec)
{
    bm::bvector<> bv;
    
    bench_vec.resize(0);
    for (unsigned i = 0; i < benchmark_max;)
    {
        unsigned idx = unsigned(rand() % str_vec.size());
        if (bv.test(idx))  // make sure benchmark example is not repeated
            continue;
        if (idx < str_vec.size())
            bench_vec.push_back(str_vec[idx]);
        else
            continue;
        bv.set(idx); // mark as set
        ++i;
    } // for
}

static
void run_benchmark(const str_sparse_vect& str_sv, const string_vector& str_vec)
{
    string_vector bench_vec;
    pick_benchmark_set(bench_vec, str_vec);
    
    bm::bvector<> bv1, bv2, bv3;

    cout << "Picked " << bench_vec.size() << " samples. Running benchmarks." << endl;
    
    unsigned bench_size = unsigned(bench_vec.size());
    {
        bm::chrono_taker tt1("3. lower_bound search", bench_size, &timing_map);
        for (const string& term : bench_vec)
        {
            auto it = std::lower_bound(str_vec.begin(), str_vec.end(), term);
            if (it != str_vec.end())
            {
                string_vector::size_type idx =
                  string_vector::size_type(std::distance(str_vec.begin(), it));
                /*
                const string& value = str_vec[idx];
                if (value != term)
                    throw runtime_error("Error. Incorrect search!");
                */
                bv1.set(unsigned(idx));
            }
        } // for
    }
    
    {
        // construct std::map<>  (which is most likely an RB-tree)
        std::map<string, unsigned> str_map;
        for (string_vector::size_type i = 0; i < str_vec.size(); ++i)
        {
            const string& s = str_vec[i];
            str_map[s] = unsigned(i);
        } // for
        
        bm::chrono_taker tt1("4. std::map<> search", bench_size, &timing_map);
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
        bm::sparse_vector_scanner<str_sparse_vect> scanner;
        bm::chrono_taker tt1("5. bm::sparse_vector_scanner<> search", bench_size, &timing_map);
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
    
    // various integrity checks
    //
    int cmp = bv1.compare(bv2);
    if (cmp != 0)
        throw runtime_error("Error. RB-search mismatch!");
    cmp = bv1.compare(bv3);
    if (cmp != 0)
        throw runtime_error("Error. scanner mismatch!");

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
            bm::chrono_taker tt1("2. build sparse vector", 1, &timing_map);

            for (const string& term : str_vec)
            {
                str_sv.push_back(term);
            }
            BM_DECLARE_TEMP_BLOCK(tb)
            str_sv.optimize(tb); // memory optimization after load
        }
        
        // save SV vector to file
        if (!sv_out_name.empty() && !str_sv.empty())
        {
            bm::chrono_taker tt1("7. Save sparse vector", 1, &timing_map);
            file_save_svector(str_sv, sv_out_name);
        }
        
        if (!sv_in_name.empty())
        {
            {
                bm::chrono_taker tt1("8. Load sparse vector", 1, &timing_map);
                file_load_svector(str_sv, sv_in_name);
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

        
        if (is_diag)
        {
            if (str_sv.size())
                print_svector_stat(str_sv, true);
            
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
            bm::chrono_taker::print_duration_map(timing_map, bm::chrono_taker::ct_all);
        }
    }
    catch (std::exception& ex)
    {
        std::cerr << "Error:" << ex.what() << std::endl;
        return 1;
    }

    return 0;
}

