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
        << "BitMagic Range Index Search Sample (c) 2019" << std::endl
        << "-i      file-name            -- input set file to parse" << std::endl
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

        if (arg == "-i" || arg == "--i" )
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
typedef bm::sparse_vector<unsigned, bm::bvector<> >     sparse_vector_u32;
typedef bm::str_sparse_vector<char, bm::bvector<>, 16>  str_sparse_vect;


/// Columnar structure for succinct vectors
///
struct RangeIndex
{
    str_sparse_vect    id_from_v;
    sparse_vector_u32  range_from_v;
    sparse_vector_u32  len_v;
    str_sparse_vect    id_to_v;
};


bm::chrono_taker::duration_map_type  timing_map;

/// Parse the input tab delimited file and load it into vectors
/// 
static
int load_range_idx(const std::string& fname, RangeIndex& range_idx)
{
    bm::chrono_taker tt1("1. parse input data ", 1, &timing_map);

    std::ifstream fin(fname.c_str(), std::ios::in);
    if (!fin.good())
        return -1;

    std::string line;
    std::regex reg("\t|, |\\s+"); // break on TAB or comma+space or spaces
    std::sregex_token_iterator it_end;
    
    string trim_chars(", ");
    
    // use back-insert iterators for sparse vectors
    //
    str_sparse_vect::back_insert_iterator id_from_bi = range_idx.id_from_v.get_back_inserter();
    sparse_vector_u32::back_insert_iterator range_from_bi = range_idx.range_from_v.get_back_inserter();
    sparse_vector_u32::back_insert_iterator len_bi = range_idx.len_v.get_back_inserter();
    str_sparse_vect::back_insert_iterator id_to_bi = range_idx.id_to_v.get_back_inserter();

    // read the input file
    //
    for (unsigned i = 0; std::getline(fin, line); ++i)
    {
        if (line.empty())
            continue;

        // regex based tokenizer
        std::sregex_token_iterator it(line.begin(), line.end(), reg, -1);
        std::vector<std::string> line_vec(it, it_end);
        if (line_vec.empty() || line_vec.size() != 4)
            continue;
        
        try
        {
            const string& id_from = line_vec.at(0);
            const string& r_from = line_vec.at(1);
            const string& r_to = line_vec.at(2);
            const string& id_to = line_vec.at(3);

            size_t p;
            unsigned from = unsigned(std::stoul(r_from, &p));
            unsigned to = unsigned(std::stoul(r_to, &p));
            if (to < from)
            {
                std::cerr << "Incorrect index for:" << std::endl;
                std::cerr << id_from << " " << r_from << " " << r_to << " " << id_to << std::endl;
                continue; // need better handling
            }
            
            id_from_bi.add(id_from.c_str());
            range_from_bi.add(from);
            unsigned len = to - from; // only need to-from diff
            len_bi.add(len);
            id_to_bi.add(id_to.c_str());
            
            //std::cout << id_from << " " << from << " " << len << " " << id_to << std::endl;
            
            
/*
            for (auto s : line_vec)
            {
                std::cout << "'" << s << "' ";
            }
            std::cout << std::endl;
*/
        /*
            // trip the extra chars
            string& col13 = line_vec.at(13);
            col13.erase(0, col13.find_first_not_of(trim_chars));
            col13.erase(col13.find_last_not_of(trim_chars) + 1);

            if (!col13.empty())
                str_vec.emplace_back(col13); */
        }
        catch (std::exception&)
        {
            // just ignore (not ideal, but ok for a sketch ETL)
        }
        
        if ((i & 0xFFFF) == 0)
        {
            cout << "\rReading input file: " << i << flush;
        }
    } // for
    
    // flush all back-insert iterators
    //
    id_from_bi.flush();
    range_from_bi.flush();
    len_bi.flush();
    id_to_bi.flush();
    
    
    cout << endl;
    return 0;
}


/// optimization-compression of range index
///
static
void optimize_range_idx(RangeIndex& range_idx)
{
    bm::chrono_taker tt1("2. optimize index ", 1, &timing_map);
    {
        str_sparse_vect v_tmp;
        v_tmp.remap_from(range_idx.id_from_v);
        range_idx.id_from_v.swap(v_tmp);
    }
    {
        str_sparse_vect v_tmp;
        v_tmp.remap_from(range_idx.id_to_v);
        range_idx.id_to_v.swap(v_tmp);
    }

    BM_DECLARE_TEMP_BLOCK(tb)
    range_idx.id_from_v.optimize(tb);
    range_idx.range_from_v.optimize(tb);
    range_idx.len_v.optimize(tb);
    range_idx.id_to_v.optimize(tb);
}

/// save vector to a file
///
template<typename SV>
int save_index_vector(const string& fn, const SV& sv)
{
    bm::sparse_vector_serial_layout<SV> sv_lay;

    BM_DECLARE_TEMP_BLOCK(tb)
    bm::sparse_vector_serialize(sv, sv_lay, tb);
    {
        std::ofstream fout(fn.c_str(), std::ios::binary);
        if (!fout.good())
            return -1;
        const char* buf = (char*)sv_lay.buf();
        fout.write(buf, (unsigned)sv_lay.size());
        if (!fout.good())
            return -1;
        fout.close();
    }
    return 0;
}

/// save range index
///
static
int save_range_idx(const string& fname, const RangeIndex& range_idx)
{
    bm::chrono_taker tt1("3. save sparse vectors", 1, &timing_map);

    int res;
    {
        std::string fn = fname;
        fn.append(".id_from.sv");
        res = save_index_vector(fn, range_idx.id_from_v);
        if (res)
            return res;
    }
    {
        std::string fn = fname;
        fn.append(".range_from.sv");
        res = save_index_vector(fn, range_idx.range_from_v);
        if (res)
            return res;
    }
    {
        std::string fn = fname;
        fn.append(".len.sv");
        res = save_index_vector(fn, range_idx.len_v);
        if (res)
            return res;
    }
    {
        std::string fn = fname;
        fn.append(".id_to.sv");
        res = save_index_vector(fn, range_idx.id_to_v);
        if (res)
            return res;
    }
    return res;
}


/// load range index
///
static
int file_load_range_idx(const string& fname, RangeIndex& range_idx)
{
    bm::chrono_taker tt1("4. load sparse vector", 1, &timing_map);

    int res;
    {
        std::string fn = fname;
        fn.append(".id_from.sv");
        res = bm::file_load_svector(range_idx.id_from_v, fn);
        if (res)
            return res;
    }
    {
        std::string fn = fname;
        fn.append(".range_from.sv");
        res = bm::file_load_svector(range_idx.range_from_v, fn);
        if (res)
            return res;
    }
    {
        std::string fn = fname;
        fn.append(".len.sv");
        res = bm::file_load_svector(range_idx.len_v, fn);
        if (res)
            return res;
    }
    {
        std::string fn = fname;
        fn.append(".id_to.sv");
        res = bm::file_load_svector(range_idx.id_to_v, fn);
        if (res)
            return res;
    }
    return res;
}

const unsigned benchmark_max = 15000;  // benchmark sampling size

static
void run_search(const RangeIndex& range_idx, const string& search_key)
{
    bm::sparse_vector_scanner<str_sparse_vect> scanner;
    scanner.bind(range_idx.id_from_v, true); // attach SV as permanent search parameter to share cached values


    {
        bm::chrono_taker tt1("5.  bm::sparse_vector_scanner<> binary search", benchmark_max, &timing_map);
        str_sparse_vect::size_type pos;
        bool found = scanner.bfind_eq_str(search_key.c_str(), pos);
        string id_from, id_to;
        unsigned from, to;
        
        for (unsigned i = 0; i < benchmark_max; ++i)
        {
            while (found)
            {
                range_idx.id_from_v.get(pos, id_from);
                if (id_from != search_key)
                    break;
                from = range_idx.range_from_v[pos];
                to = range_idx.len_v[pos];
                to += from;
                range_idx.id_to_v.get(pos, id_to);
                
                //cout << id_from << "\t" << from << "\t" << to << "\t" << id_to << endl;
                ++pos;
                found = pos < range_idx.id_from_v.size();
            } // while
        }
    }

}



int main(int argc, char *argv[])
{
    if (argc < 3)
    {
        show_help();
        return 1;
    }
    
    RangeIndex r_idx;

    try
    {
        auto ret = parse_args(argc, argv);
        if (ret != 0)
        {
            return ret;
        }
        if (!i_name.empty())
        {
            auto res = load_range_idx(i_name, r_idx);
            if (res != 0)
            {
                return res;
            }
            optimize_range_idx(r_idx);
        }
        if (!sv_out_name.empty())
        {
            save_range_idx(sv_out_name, r_idx);
        }


        if (!sv_in_name.empty())
        {
            file_load_range_idx(sv_in_name, r_idx);
        }
        

        if (is_diag)
        {
/*
            std::cout << "ID_FROM:" << std::endl;
            print_svector_stat(r_idx.id_from_v, false);
            std::cout << "RANGE_FROM:" << std::endl;
            print_svector_stat(r_idx.range_from_v, false);
            std::cout << "LEN:" << std::endl;
            print_svector_stat(r_idx.len_v, false);
*/

            //print_svector_xor_stat(r_idx.id_to_v);
/*
            std::cout << "ID_TO:" << std::endl;
            print_svector_stat(r_idx.id_to_v, false);
*/
        }
        
        run_search(r_idx, "CUU0");

        if (is_timing)  // print all collected timings
        {
            std::cout << std::endl << "Performance:" << std::endl;
            bm::chrono_taker::print_duration_map(timing_map, bm::chrono_taker::ct_time);
        }
    }
    catch (std::exception& ex)
    {
        std::cerr << "Error:" << ex.what() << std::endl;
        return 1;
    }

    return 0;
}

