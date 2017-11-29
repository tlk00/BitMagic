/*
Copyright(c) 2002-2017 Anatoliy Kuznetsov(anatoliy_kuznetsov at yahoo.com)

Permission is hereby granted, free of charge, to any person
obtaining a copy of this software and associated documentation
files (the "Software"), to deal in the Software without restriction,
including without limitation the rights to use, copy, modify, merge,
publish, distribute, sublicense, and/or sell copies of the Software,
and to permit persons to whom the Software is furnished to do so,
subject to the following conditions:

The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
OTHER DEALINGS IN THE SOFTWARE.

You have to explicitly mention BitMagic project in any derivative product,
its WEB Site, published materials, articles or any other work derived from this
project or based on our code or know-how.

For more information please visit:  http://bitmagic.io

*/

#include <iostream>
#include <sstream>
#include <chrono>
#include <regex>
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
#include "bmrandom.h"
#include "bmsparsevec.h"
#include "bmsparsevec_algo.h"
#include "bmsparsevec_serial.h"
#include "bmalgo_similarity.h"


#include "bmdbg.h"
#include "bmtimer.h"


void show_help()
{
    std::cerr
      << "BitMagic Link Test Utility (c) 2017"                     << std::endl
      << "-iset   file-name          -- input set file to parse"   << std::endl
      << "-lin    link-input-file    -- test pairs to load"        << std::endl
      << "-lmout  link-matrix-name   -- Link Matrix name to save"  << std::endl
      << "-lmin   link-matrix-name   -- Link Matrix name to load"  << std::endl
      << "-diag                      -- run diagnostics"           << std::endl
      << "-timing                    -- collect timings"           << std::endl
      << "-bench                     -- run benchmark"             << std::endl
      ;
}




// Arguments
//
std::string  ln_in_file;
std::string  lm_out_name;
std::string  lm_in_name;
std::string  iset_name;
bool         is_diag = false;
bool         is_timing = false;
bool         is_bench = false;


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
        
        if (arg == "-lin" || arg == "--lin")
        {
            if (i + 1 < argc)
            {
                ln_in_file = argv[++i];
            }
            else
            {
                std::cerr << "Error: -lin requires file name" << std::endl;
                return 1;
            }
            continue;
        }
        if (arg == "-lmout" || arg == "--lmout")
        {
            if (i + 1 < argc)
            {
                lm_out_name = argv[++i];
            }
            else
            {
                std::cerr << "Error: -lmout requires file name" << std::endl;
                return 1;
            }
            continue;
        }
        if (arg == "-lmin" || arg == "--lmin")
        {
            if (i + 1 < argc)
            {
                lm_in_name = argv[++i];
            }
            else
            {
                std::cerr << "Error: -lmin requires file name" << std::endl;
                return 1;
            }
            continue;
        }
        if (arg == "-iset" || arg == "--iset")
        {
            if (i + 1 < argc)
            {
                iset_name = argv[++i];
            }
            else
            {
                std::cerr << "Error: -iset requires file name" << std::endl;
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

// ------------------------------------------------------------------
// small utility functions
//

void print_vector(const std::vector<unsigned>& vect)
{
    std::cout << vect.size() << ": [";
    for (size_t i = 0; i < vect.size(); ++i)
    {
        std::cout << vect[i] << ", ";
    }
    std::cout << "]" << std::endl;
}
void convert_to_delta(std::vector<unsigned>& vect)
{
    for (size_t k = vect.size()-1; k >= 1; --k)
    {
        vect[k] -= vect[k-1];
        --vect[k];
    } // for k
}
void convert_from_delta(std::vector<unsigned>& vect)
{
    for (unsigned j = 1; j < vect.size(); ++j)
    {
        vect[j] += vect[j-1] + 1;
    } // for j
}



// Globals
//
typedef bm::sparse_vector<unsigned, bm::bvector<> > sparse_vector_u32;




/// compressed sparse vector
///
struct compress_svector
{
    sparse_vector_u32::bvector_type   bv_descr;  ///< bit-vector descriptor
    sparse_vector_u32                 sv_stor;   ///< sparse-vector store
    sparse_vector_u32::bvector_type::blocks_count blocks_bc;
    
    void load_from(const sparse_vector_u32& sv);
    bool equal(const sparse_vector_u32& sv) const;
    
    sparse_vector_u32::value_type get(unsigned i) const;
    void optimize_gap_size();
    
    void count_blocks();
    
};

void compress_svector::load_from(const sparse_vector_u32& sv)
{
    bm::compute_nonzero_bvector(sv, bv_descr);
    sparse_vector_u32::bvector_type::counted_enumerator enc = bv_descr.first();
    for (;enc.valid(); ++enc)
    {
        unsigned idx = *enc;
        unsigned pos = enc.count();
        
        unsigned v = sv.get(idx);
        
        if (v == 0)
        {
            std::cerr << "Compress vector build error" << std::endl;
            exit(1);
        }
        
        sv_stor.set(pos, v);
    } // for en
    
    BM_DECLARE_TEMP_BLOCK(tb)
    bv_descr.optimize(tb);
    sv_stor.optimize(tb);
}

void compress_svector::count_blocks()
{
    bv_descr.running_count_blocks(&blocks_bc); // compute popcount skip list
}

void compress_svector::optimize_gap_size()
{
    bv_descr.optimize_gap_size();
    sv_stor.optimize_gap_size();
}

sparse_vector_u32::value_type compress_svector::get(unsigned i) const
{
    sparse_vector_u32::value_type v = 0;
    bool found = bv_descr.test(i);
    if (!found)
        return v;
    
    unsigned pos = bv_descr.count_to(i, blocks_bc);
    v = sv_stor.get(pos);
    return v;
}

bool compress_svector::equal(const sparse_vector_u32& sv) const
{
    for (unsigned i = 0; i < sv.size(); ++i)
    {
        unsigned v0 = sv.get(i);
        unsigned v1 = get(i);
        if (v0 != v1)
            return false;
        if (i % 10000 == 0)
        {
            std::cout << "\r" << i << " / " << sv.size() << std::flush;
        }
    }
    std::cout << std::endl;
    return true;
}


/// Sparse matrix storage for experiments with precalculated ER joins
///
struct link_matrix
{
    bm::bvector<> bv_from;
    bm::bvector<> bv_to;

    bm::bvector<> bv_11_flag;


    sparse_vector_u32  sv_11;      ///< vector of 1x1 relationships
    sparse_vector_u32  sv_offs;    ///< vector of location offsets
    sparse_vector_u32  sv_sz;      ///< vector of sizes
    sparse_vector_u32  sv_1m;      ///< vector of 1xM relationships
    
    unsigned off;  /// accumulated offset for vectors
    
    compress_svector  sv_11_c;  ///< compressed vector of 1x1 relationships
    compress_svector  sv_offs_c;///< compressed vector of offsets
    
    link_matrix()
     : bv_from(bm::BM_GAP),
       bv_to(bm::BM_GAP),
       bv_11_flag(bm::BM_GAP),
       off(0)
    {}
    
    /// run memory optimization/compression
    void optimize();
    
    void add_vector(unsigned id_from, std::vector<unsigned>& vect);
    bool get_vector(unsigned id_from, std::vector<unsigned>& vect) const;
    
    /// print statistics
    void print_stat() const;
    
    /// Save link matrix as a collection of pre-calculated files
    void save(const std::string& base_name) const;
    
    /// Load link matrix from a collection of files
    void load(const std::string& base_name);
};



void link_matrix::optimize()
{
    BM_DECLARE_TEMP_BLOCK(tb)
    sparse_vector_u32::statistics st;

    bv_from.optimize(tb);
    bv_to.optimize(tb);
    bv_11_flag.optimize(tb);
    
    sv_11.optimize(tb, bm::bvector<>::opt_compress, &st);
    sv_offs.optimize(tb, bm::bvector<>::opt_compress, &st);
    sv_sz.optimize(tb, bm::bvector<>::opt_compress, &st);
    sv_1m.optimize(tb, bm::bvector<>::opt_compress, &st);
}

void link_matrix::print_stat() const
{
    std::cout << "\nsv 11 statistics:" << std::endl;
    std::cout << "-----------------" << std::endl;
    bm::print_svector_stat(sv_11, false);
    std::cout << "\nsv offs statistics:" << std::endl;
    std::cout << "-----------------" << std::endl;
    bm::print_svector_stat(sv_offs, false);
    std::cout << "\nsv size statistics:" << std::endl;
    std::cout << "-----------------" << std::endl;
    bm::print_svector_stat(sv_sz, false);
    std::cout << "\nsv 1M statistics:" << std::endl;
    std::cout << "-----------------" << std::endl;
    bm::print_svector_stat(sv_1m, false);


    std::cout << "\nbvector-from statistics:" << std::endl;
    std::cout << "-----------------" << std::endl;
    bm::print_bvector_stat(bv_from);
    
    std::cout << "\nbvector-to statistics:" << std::endl;
    std::cout << "-----------------" << std::endl;
    bm::print_bvector_stat(bv_to);


    std::cout << "\nsv 11 C statistics:" << std::endl;
    std::cout << "-----------------" << std::endl;
    bm::print_svector_stat(sv_11_c.sv_stor, false);
    bm::print_bvector_stat(sv_11_c.bv_descr);

    std::cout << "\nsv OFFS C statistics:" << std::endl;
    std::cout << "-----------------" << std::endl;
    bm::print_svector_stat(sv_offs_c.sv_stor, false);
    bm::print_bvector_stat(sv_offs_c.bv_descr);

}

void link_matrix::save(const std::string& base_name) const
{
    {
    std::string bv_from_fname = base_name;
    bv_from_fname.append("_bvfrom.bv");
    SaveBVector(bv_from_fname.c_str(), bv_from);
    }

    {
    std::string bv_to_fname = base_name;
    bv_to_fname.append("_bvto.bv");
    SaveBVector(bv_to_fname.c_str(), bv_to);
    }
    
    {
    std::string bv_11f_fname = base_name;
    bv_11f_fname.append("_bv11f.bv");
    SaveBVector(bv_11f_fname.c_str(), bv_11_flag);
    }

    {
    std::string sv11_fname = base_name;
    sv11_fname.append("_sv11.sv");
    int res = file_save_svector(sv_11, sv11_fname);
    if (res != 0)
    {
        std::cerr << "File save error." << std::endl;
    }
    }

    {
    std::string sv11c_fname = base_name;
    sv11c_fname.append("_sv11-c.sv");
    int res = file_save_svector(sv_11_c.sv_stor, sv11c_fname);
    if (res != 0)
    {
        std::cerr << "File save error." << std::endl;
    }
    sv11c_fname = base_name;
    sv11c_fname.append("_sv11-c.bv");
    SaveBVector(sv11c_fname.c_str(), sv_11_c.bv_descr);
    }


    {
    std::string svoffs_fname = base_name;
    svoffs_fname.append("_svoffs.sv");
    int res = file_save_svector(sv_offs, svoffs_fname);
    if (res != 0)
    {
        std::cerr << "File save error." << std::endl;
    }
    }
    
    {
    std::string svoffsc_fname = base_name;
    svoffsc_fname.append("_svoffs-c.sv");
    int res = file_save_svector(sv_offs_c.sv_stor, svoffsc_fname);
    if (res != 0)
    {
        std::cerr << "File save error." << std::endl;
    }
    svoffsc_fname = base_name;
    svoffsc_fname.append("_svoffs-c.bv");
    SaveBVector(svoffsc_fname.c_str(), sv_offs_c.bv_descr);
    }

    
    {
    std::string svsz_fname = base_name;
    svsz_fname.append("_svsize.sv");
    int res = file_save_svector(sv_sz, svsz_fname);
    if (res != 0)
    {
        std::cerr << "File save error." << std::endl;
    }
    }

    {
    std::string sv1m_fname = base_name;
    sv1m_fname.append("_sv1m.sv");
    int res = file_save_svector(sv_1m, sv1m_fname);
    if (res != 0)
    {
        std::cerr << "File save error." << std::endl;
    }
    }
}

void link_matrix::load(const std::string& base_name)
{

    {
    std::string bv_from_fname = base_name;
    bv_from_fname.append("_bvfrom.bv");
    LoadBVector(bv_from_fname.c_str(), bv_from);
    //bv_from.optimize_gap_size();
    }

    {
    std::string bv_to_fname = base_name;
    bv_to_fname.append("_bvto.bv");
    LoadBVector(bv_to_fname.c_str(), bv_to);
    //bv_to.optimize_gap_size();
    }

    {
    std::string bv_11f_fname = base_name;
    bv_11f_fname.append("_bv11f.bv");
    LoadBVector(bv_11f_fname.c_str(), bv_11_flag);
    //bv_11_flag.optimize_gap_size();
    }

    /*
    {
    std::string sv11_fname = base_name;
    sv11_fname.append("_sv11.sv");
    int res = file_load_svector(sv_11, sv11_fname);
    if (res != 0)
    {
        std::cerr << "File load error." << std::endl;
    }
    }
    */
    
    {
    std::string sv11c_fname = base_name;
    sv11c_fname.append("_sv11-c.sv");
    int res = file_load_svector(sv_11_c.sv_stor, sv11c_fname);
    if (res != 0)
    {
        std::cerr << "File load error." << std::endl;
    }
    sv11c_fname = base_name;
    sv11c_fname.append("_sv11-c.bv");
    LoadBVector(sv11c_fname.c_str(), sv_11_c.bv_descr);
    
    sv_11_c.count_blocks();
    
    //sv_11_c.optimize_gap_size();
    }

    /*
    {
    std::string svoffs_fname = base_name;
    svoffs_fname.append("_svoffs.sv");
    int res = file_load_svector(sv_offs, svoffs_fname);
    if (res != 0)
    {
        std::cerr << "File load error." << std::endl;
    }
    }
    */
    
    {
    std::string svoffsc_fname = base_name;
    svoffsc_fname.append("_svoffs-c.sv");
    int res = file_load_svector(sv_offs_c.sv_stor, svoffsc_fname);
    if (res != 0)
    {
        std::cerr << "File load error." << std::endl;
    }
    svoffsc_fname = base_name;
    svoffsc_fname.append("_svoffs-c.bv");
    LoadBVector(svoffsc_fname.c_str(), sv_offs_c.bv_descr);
    
    sv_offs_c.count_blocks();
    
    //sv_offs_c.optimize_gap_size();
    }

    
    {
    std::string svsz_fname = base_name;
    svsz_fname.append("_svsize.sv");
    int res = file_load_svector(sv_sz, svsz_fname);
    if (res != 0)
    {
        std::cerr << "File load error." << std::endl;
    }
    //sv_sz.optimize_gap_size();
    }

    {
    std::string sv1m_fname = base_name;
    sv1m_fname.append("_sv1m.sv");
    int res = file_load_svector(sv_1m, sv1m_fname);
    if (res != 0)
    {
        std::cerr << "File load error." << std::endl;
    }
    //sv_1m.optimize_gap_size();
    }
    
    // ---------------------------------
    // build compress vectors
    /*
    std::cout << "Building compressed vector (sv_11) " << std::endl;
    sv_11_c.load_from(sv_11);
    std::cout << "OK " << std::endl;
    
    std::cout << "Building compressed vector (sv_offs) " << std::endl;
    sv_offs_c.load_from(sv_offs);
    std::cout << "OK " << std::endl;
    */
    
    /*
    std::cout << "Validating compressed vector " << std::endl;
    bool check = sv_11_c.equal(sv_11);
    if (!check)
    {
        std::cerr << "Compressed sparse comparison failed" << std::endl;
        exit(1);
    }
    */
}


bool link_matrix::get_vector(unsigned id_from, std::vector<unsigned>& vect) const
{
/*
    if (!bv_from.test(id_from))
    {
        std::cerr << "id from not found" << std::endl;
        exit(1);
        vect.resize(0);
        return false;
    }
*/
    unsigned vc;
    //std::cout << id_from << ", " << std::flush;
    if (bv_11_flag.test(id_from))
    {
        vc = sv_11_c.get(id_from);
        /*
        unsigned v = sv_11.get(id_from);
        
        if (v == 0)
        {
            std::cerr << "Single value find error" << std::endl;
            exit(1);
            vect.resize(0);
            return false;
        }
     
        if (v != vc)
        {
            std::cerr << "Single value get error" << std::endl;
            exit(1);
        }
        */
     
        vect.push_back(vc);
        return true;
    }
    
    
    unsigned offs_c = sv_offs_c.get(id_from);
/*
    unsigned offs = sv_offs.get(id_from);
    if (offs_c != offs)
    {
        std::cerr << "offset check error" << std::endl;
        exit(1);
    }
*/
    unsigned sz = sv_sz.get(id_from);
    vect.resize(sz);
    if (sz==0)
    {
        return false;
        std::cerr << "vector size error" << std::endl;
        exit(1);
    }
    
    vc = sv_11_c.get(id_from);
/*
    unsigned v = sv_11.get(id_from);
    if (v != vc)
    {
        std::cerr << "Array 0 value get error" << std::endl;
        exit(1);
    }
*/
    vect[0] = vc;
    for (unsigned j = 1; j < sz; ++j)
    {
        unsigned v = sv_1m.get(j + offs_c - 1);
        vect[j] = v;
    }
    
    //convert_from_delta(vect);

    return true;
}

void link_matrix::add_vector(unsigned id_from, std::vector<unsigned>& vect)
{
    size_t sz = vect.size();
    if (sz == 0) // nothing to do
        return;
    
    if (bv_from.test(id_from))
    {
        std::cerr << "Duplicate unsorted id = " << id_from << std::endl;
        exit(1);
    }
    
    bv_from.set(id_from);
    
    if (sz == 1)  // 1x1 relationship
    {
        sv_11.set(id_from, vect[0]);
        bv_to.set(vect[0]);
        bv_11_flag.set(id_from);
        return;
    }
    
    // add all mapped elements to the target vector
    bm::combine_or(bv_to, vect.begin(), vect.end());

    //convert_to_delta(vect);

    sv_offs.set(id_from, off);
    sv_sz.set(id_from, sz);  // save true size
    sz--;
    sv_11.set(id_from, vect[0]);
    sv_1m.import(&vect[1], sz, off);
    
    off += sz;
}

link_matrix            link_m;
bm::chrono_taker::duration_map_type  timing_map;


// parse the input re-mapping vector
//
int load_bv(const std::string& fname, bm::bvector<>& bv)
{
    std::string line;

    std::string regex_str = "<Id>[0-9]+</Id>";
    std::regex reg1(regex_str);

    std::string regex_str2 = "[0-9]+";
    std::regex reg2(regex_str2);


    std::ifstream fin(fname.c_str(), std::ios::in);
    if (!fin.good())
    {
        return -1;
    }

    unsigned id = 0;
    std::vector<unsigned> vbuf;

    for (unsigned i = 0; std::getline(fin, line); i++)
    {
        if (line.empty())
            continue;

        unsigned id_from = 0;
        unsigned id_to = 0;

        std::sregex_iterator it(line.begin(), line.end(), reg1);
        std::sregex_iterator it_end;
        for (unsigned j = 0; it != it_end; ++it, ++j)
        {
            std::smatch match = *it;
            std::string ms = match.str();

            if (ms.empty())
                continue;

            std::sregex_iterator it1(ms.begin(), ms.end(), reg2);
            std::sregex_iterator it1_end;
            for (; it1 != it1_end; ++it1)
            {
                std::smatch match2 = *it1;
                std::string ms2 = match2.str();
                unsigned long id = std::stoul(ms2);
                if (id)
                {
                    bv.set(id);
                }
                
            } // for

        } // for j
    } // for i
    bv.optimize();

    std::cout << "input vector count=" << bv.count() << std::endl;
    return 0;
}


// load sparse_vector from a file
//
int load_ln(const std::string& fname, link_matrix& lm)
{
    std::string line;
    
    std::string regex_str = "[0-9]+";
    std::regex reg1(regex_str);
    
    std::ifstream fin(fname.c_str(), std::ios::in);
    if (!fin.good())
    {
        return -1;
    }
    
    unsigned id = 0;
    std::vector<unsigned> vbuf;
    
    
    for (unsigned i = 0; std::getline(fin, line); i++)
    {
        if (line.empty())
            continue;
        
        unsigned id_from = 0;
        unsigned id_to = 0;
        
        std::sregex_iterator it(line.begin(), line.end(), reg1);
        std::sregex_iterator it_end;
        for (unsigned j = 0; it != it_end; ++it, ++j)
        {
            std::smatch match = *it;
            std::string ms = match.str();
            
            unsigned long ul = std::stoul(ms);
            switch (j)
            {
            case 0: id_from = (unsigned) ul; break;
            case 1: id_to = (unsigned) ul; break;
            default:
                break;
            }
            
        } // for
        
        // mini-FSM for vector accumulation
        if (id == 0)
        {
            id = id_from;
            vbuf.push_back(id_to);
        }
        else
        if (id == id_from)
        {
            vbuf.push_back(id_to);
        }
        else
        {
            if (vbuf.size() > 1)
                std::sort(vbuf.begin(), vbuf.end());

            lm.add_vector(id, vbuf);
            /*
            {
                std::vector<unsigned> vcheck;
                bool found = lm.get_vector(id, vcheck);
                if (!found)
                {
                    std::cerr << "Vector not found! " << id << std::endl;
                    exit(1);
                }
                if (vbuf.size() != vcheck.size())
                {
                    std::cerr << "Vector size mismatch! " << id << std::endl;
                    exit(1);
                }
                //convert_from_delta(vbuf);
                for (unsigned k = 0; k < vbuf.size(); ++k)
                {
                    if (vcheck[k] != vbuf[k])
                    {
                        std::cerr << "Vector  mismatch. " << id
                                  << " " << k
                                  << std::endl;
                        
                        
                        print_vector(vbuf);
                        print_vector(vcheck);
                        exit(1);
                    }
                } // for
            }
            */
            
            vbuf.resize(0);
            
            id = id_from;
            vbuf.push_back(id_to);
        }
        
        if ((i != 0) && (i % 1000000) == 0)
        {
            std::cout << "\r" << i << std::flush;
            lm.optimize();
        }
    } // for getline()
    
    if (vbuf.size())
    {
        std::sort(vbuf.begin(), vbuf.end());
        lm.add_vector(id, vbuf);
    }
    lm.optimize();
    
    std::cout << "Building compressed vector (sv_11) " << std::endl;
    lm.sv_11_c.load_from(lm.sv_11);
    std::cout << "OK " << std::endl;
    
    std::cout << "Building compressed vector (sv_offs) " << std::endl;
    lm.sv_offs_c.load_from(lm.sv_offs);
    std::cout << "OK " << std::endl;


    std::cout << std::endl;
    

    return 0;
}

unsigned benchmark_ops = 1000;

// id remapping benchmark
//
void run_benchmark(link_matrix& lm)
{
    if (!lm.bv_from.any())  // nothing to do
        return;
    BM_DECLARE_TEMP_BLOCK(tb)

    std::vector<bm::bvector<> > sample_vectors;
    bm::bvector<>  bv_sample_to;
    bm::bvector<>  bv_res;

    {
    bm::chrono_taker tt1("3. generation of remapping samples", benchmark_ops, &timing_map);
    
        bm::random_subset<bm::bvector<> > rsampler;

        for (unsigned i = 0; i < benchmark_ops; ++i)
        {
            unsigned msize = rand() % 10000;
            if (msize < 5) msize = 5;

            {
                bm::bvector<> bv_subset;
                rsampler.sample(bv_subset, lm.bv_from, msize);
                bv_subset.optimize(tb);
                sample_vectors.emplace_back(std::move(bv_subset));
            }
        } // for i
        
        rsampler.sample(bv_sample_to, lm.bv_to, 10000); // sample subset from the other end of join
    }

    {
    bm::chrono_taker tt1("4. remapping", benchmark_ops, &timing_map);
    
        std::vector<unsigned> vect;
        bm::bvector<>         bv_remap;
    
        for (unsigned i = 0; i < benchmark_ops; ++i)
        {
            const bm::bvector<>& bv = sample_vectors[i];
            
            bm::bvector<>::enumerator en = bv.first();
            for (unsigned k = 0; en.valid(); ++k, ++en)
            {
                unsigned id = *en;
                
                lm.get_vector(id, vect);
                if (vect.size())
                {
                    bm::combine_or(bv_remap, vect.begin(), vect.end());
                }
            } // for k
            
            bv_res = bv_remap;
            bv_res &= bv_sample_to;
        } // for i
    }


}

// id remapping benchmark
//
void remap(const link_matrix& lm, const bm::bvector<>& bv_in)
{
    if (!lm.bv_from.any())  // nothing to do
        return;
    BM_DECLARE_TEMP_BLOCK(tb)


    {
        bm::chrono_taker tt1("5. input remap", 1, &timing_map);

        std::vector<unsigned> vect;
        bm::bvector<>         bv_remap;

        const bm::bvector<>& bv = bv_in;

        bm::bvector<>::enumerator en = bv.first();
        for (unsigned k = 0; en.valid(); ++k, ++en)
        {
            unsigned id = *en;

            lm.get_vector(id, vect);
            if (vect.size())
            {
                bm::combine_or(bv_remap, vect.begin(), vect.end());
            }
        } // for k

        
        std::cout << "Remap count =" << bv_remap.count() << std::endl;
    }
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

        bm::bvector<> bv_inp(bm::BM_GAP);

        
        if (!ln_in_file.empty())
        {
            bm::chrono_taker tt1("0. link pairs load", 1, &timing_map);
            auto res = load_ln(ln_in_file, link_m);
            if (res != 0)
            {
                return res;
            }
        }

        if (!iset_name.empty())
        {
            auto res = load_bv(iset_name, bv_inp);
            if (res != 0)
            {
                return res;
            }
        }
        
        if (!lm_in_name.empty())
        {
            bm::chrono_taker tt1("2. matrix load", 1, &timing_map);
            link_m.load(lm_in_name);
        }
        
        if (!lm_out_name.empty())
        {
            bm::chrono_taker tt1("1. matrix save", 1, &timing_map);
            link_m.save(lm_out_name);
        }


        if (is_diag)
        {
            link_m.print_stat();
        }

        if (is_bench)
        {
            run_benchmark(link_m);
        }

        if (bv_inp.any())
        {
            remap(link_m, bv_inp);
        }

        
        
        if (is_timing)  // print all collected timings
        {
            std::cout << std::endl << "Timings (ms):" << std::endl;
            bm::chrono_taker::print_duration_map(timing_map);
            
            std::cout << std::endl << "Performance (ops/sec):" << std::endl;
            bm::chrono_taker::print_duration_map(timing_map, bm::chrono_taker::ct_ops_per_sec);

        }
        
        getchar();
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



