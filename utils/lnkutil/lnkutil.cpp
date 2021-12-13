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

using namespace std;

//#define BMAVX2OPT

#include "bm.h"
#include "bmalgo.h"
#include "bmserial.h"
#include "bmrandom.h"
#include "bmsparsevec.h"
#include "bmsparsevec_algo.h"
#include "bmsparsevec_serial.h"
#include "bmalgo_similarity.h"
#include "bmsparsevec_util.h"


#include "bmdbg.h"
#include "bmtimer.h"
#include "bmundef.h" /* clear the pre-proc defines from BM */

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

// ----------------------------------------------------------------------------

/// raw data read accumulator for sparse matrix data
///
struct sm_accum
{
    typedef std::vector< std::vector<unsigned> >     dense_matr_type;
    typedef bm::sv_addr_resolver<sparse_vector_u32>  sparse_addr_resolver_type;
    typedef sparse_addr_resolver_type::bvector_type  bvector_type;
    typedef bm::serializer<bvector_type>             serializer_type;
    typedef serializer_type::buffer                  buffer_type;
    typedef std::vector<buffer_type>                 buffer_collection_type;

    void add_pair(bm::id_t id_from, bm::id_t id_to);
    void optimize(bool all);

    void release();

    const bvector_type& get_addr_bvector() const { return addr_resolver.get_bvector(); }
    const sparse_addr_resolver_type& get_resolver() const { return addr_resolver; }
    buffer_type& get_buffer(bm::id_t addr_idx) { return buffer_storage.at(addr_idx); }

private:
    void optimize(bm::id_t id_from, unsigned cut_off);

private:
    dense_matr_type             link_storage;
    buffer_collection_type      buffer_storage;
    sparse_addr_resolver_type   addr_resolver;
};

void sm_accum::optimize(bool all)
{
    unsigned cut_off = all ? 0 : 512;
    const bvector_type& bv_set = addr_resolver.get_bvector();
    bvector_type::enumerator en = bv_set.first();
    for (; en.valid(); ++en)
    {
        bm::id_t id_from = *en;
        optimize(id_from, cut_off);
    }
    if (all)
    {
        addr_resolver.optimize();
    }
}

void sm_accum::release()
{
    link_storage.resize(0);

    buffer_storage.resize(0);
    buffer_storage.shrink_to_fit();
}

void sm_accum::optimize(bm::id_t id_from, unsigned cut_off)
{
    bm::id_t addr_idx;
    bool found = addr_resolver.resolve(id_from, &addr_idx);
    if (found)
    {
        std::vector<unsigned>& vu = link_storage.at(addr_idx);
        if (vu.size())
        {
            buffer_type& bv_buf = buffer_storage.at(addr_idx);

            bvector_type bv;
            BM_DECLARE_TEMP_BLOCK(tb)

            size_t bv_buf_size = bv_buf.size();
            if (bv_buf_size && (bv_buf_size >= cut_off))
            {
                bm::deserialize(bv, bv_buf.buf());
                bv_buf.resize(0);
            }
            bm::combine_or(bv, vu.begin(), vu.end());
            vu.resize(0);
            vu.shrink_to_fit();

            bm::serializer<bvector_type> bvs(tb);
            bvs.set_compression_level(4);
            bvs.gap_length_serialization(false);
            bvs.byte_order_serialization(false);

            bvector_type::statistics st;
            bv.optimize(tb, bvector_type::opt_compress, &st);

            bvs.serialize(bv, bv_buf, &st);
        }
    }
}

void sm_accum::add_pair(bm::id_t id_from, bm::id_t id_to)
{
    bm::id_t addr_idx;
    bool found = addr_resolver.resolve(id_from, &addr_idx);
    if (found)
    {
        std::vector<unsigned>& vu = link_storage.at(addr_idx);
        vu.push_back(id_to);
    }
    else
    {
        addr_resolver.set(id_from);
        found = addr_resolver.resolve(id_from, &addr_idx);
        assert(found);
        if (addr_idx >= link_storage.size())
        {
            link_storage.resize(addr_idx + 1);
            buffer_storage.resize(addr_idx + 1);
        }

        std::vector<unsigned>& vu = link_storage.at(addr_idx);
        vu.push_back(id_to);
    }
}




/// Accumulator for unsorted pairs
///
struct pairs_accum
{
    typedef std::pair<unsigned, unsigned> pair_u32_type;
    typedef std::vector<pair_u32_type>    pair_vector_type;

    void add_pair(unsigned first, unsigned second);
    void sort();

    pair_vector_type  pair_vec;
};

void pairs_accum::add_pair(unsigned first, unsigned second)
{
    pair_vec.push_back(std::make_pair(first, second));
}

void pairs_accum::sort()
{
    std::sort(pair_vec.begin(), pair_vec.end()); // sort by first AND second
}



// -----------------------------------------------------------------------------

/// compressed sparse vector
///
struct compress_svector
{
    typedef bm::bvps_addr_resolver<bm::bvector<> > address_resolver_type;
    
    void load_from(const sparse_vector_u32& sv);
    bool equal(const sparse_vector_u32& sv) const;

    
    sparse_vector_u32::value_type get(unsigned i) const;
    void optimize_gap_size();
    
    void count_blocks();
    unsigned size() const { return bv_ares.get_bvector().count(); }
    
    address_resolver_type             bv_ares;   /// bit-vector, prefix sum address resolver
    sparse_vector_u32                 sv_stor;   ///< sparse-vector store
};

void compress_svector::load_from(const sparse_vector_u32& sv)
{
    bm::sparse_vector_scanner<sparse_vector_u32> scanner;
    bm::bvector<>& bv_descr = bv_ares.get_bvector();
    scanner.find_nonzero(sv, bv_descr);
    bv_ares.sync();

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
    bv_ares.optimize(tb);
    sv_stor.optimize(tb);
}

void compress_svector::count_blocks()
{
    bv_ares.sync();
}

void compress_svector::optimize_gap_size()
{
    sv_stor.optimize_gap_size();
}

sparse_vector_u32::value_type compress_svector::get(unsigned i) const
{
    sparse_vector_u32::value_type v = 0;
    bm::id_t pos;
    bool found = bv_ares.get(i, &pos);
    if (!found)
        return v;
    
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
    typedef bm::bvector<>    bvector_type;

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

    bm::compressed_buffer_collection<bvector_type>  buf_coll;  ///< compressed buffers 

    link_matrix()
     : bv_from(bm::BM_GAP),
       bv_to(bm::BM_GAP),
       bv_11_flag(bm::BM_GAP),
       off(0)
    {}
    
    /// run memory optimization/compression
    void optimize();
    
    void add_vector(unsigned id_from, std::vector<unsigned>& vect);
    void add_vector(unsigned id_from, const bm::bvector<>& bv);
    bool get_vector(unsigned id_from, std::vector<unsigned>& vect) const;
    bool combine_or(unsigned id_from, bm::bvector<>& bv) const;
    
    /// print statistics
    void print_stat() const;
    
    /// Save link matrix as a collection of pre-calculated files
    void save(const std::string& base_name) const;
    
    /// Load link matrix from a collection of files
    void load(const std::string& base_name);

    /// Load from accumulated matrix
    ///
    void load_from_acc_matrix(sm_accum& sm_acc);
};

void link_matrix::load_from_acc_matrix(sm_accum& sm_acc)
{
    const sm_accum::sparse_addr_resolver_type& addr_resolver = sm_acc.get_resolver();
    bvector_type bv;

    const bvector_type& bv_addr = sm_acc.get_addr_bvector();
    bvector_type::enumerator en = bv_addr.first();
    for (; en.valid(); ++en)
    {
        bm::id_t id_from = *en;
        bm::id_t addr_idx;
        bool found = addr_resolver.resolve(id_from, &addr_idx);
        if (found)
        {
            sm_accum::buffer_type& buf = sm_acc.get_buffer(addr_idx);
            if (buf.size() < 60)  // TODO: need a better criteria
            {
                bm::deserialize(bv, buf.buf());
                add_vector(id_from, bv);
                bv.clear(true);
            }
            else
            {
                buf_coll.move_buffer(id_from, buf);
                if (bv_from.test(id_from))
                {
                    std::cerr << "Duplicate unsorted id = " << id_from << std::endl;
                    exit(1);
                }

                bv_from.set(id_from);
            }
            buf.release();
        }
    } // for en

    buf_coll.sync();
    buf_coll.optimize();

    sm_acc.release();
}


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
    
    buf_coll.optimize(tb);
}

void link_matrix::print_stat() const
{
    std::cout << "\nsv 11 statistics:" << std::endl;
    std::cout << "-----------------" << std::endl;
    bm::print_svector_stat(cout, sv_11, false);
    std::cout << "\nsv offs statistics:" << std::endl;
    std::cout << "-----------------" << std::endl;
    bm::print_svector_stat(cout, sv_offs, false);
    std::cout << "\nsv size statistics:" << std::endl;
    std::cout << "-----------------" << std::endl;
    bm::print_svector_stat(cout, sv_sz, false);
    std::cout << "\nsv 1M statistics:" << std::endl;
    std::cout << "-----------------" << std::endl;
    bm::print_svector_stat(cout, sv_1m, false);


    std::cout << "\nbvector-from statistics:" << std::endl;
    std::cout << "-----------------" << std::endl;
    bm::print_bvector_stat(cout, bv_from);
    
    std::cout << "\nbvector-to statistics:" << std::endl;
    std::cout << "-----------------" << std::endl;
    bm::print_bvector_stat(cout, bv_to);


    std::cout << "\nsv 11 C statistics:" << std::endl;
    std::cout << "-----------------" << std::endl;
    bm::print_svector_stat(cout, sv_11_c.sv_stor, false);
    bm::print_bvector_stat(cout, sv_11_c.bv_ares.get_bvector());

    std::cout << "\nsv OFFS C statistics:" << std::endl;
    std::cout << "-----------------" << std::endl;
    bm::print_svector_stat(cout, sv_offs_c.sv_stor, false);
    bm::print_bvector_stat(cout, sv_offs_c.bv_ares.get_bvector());

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
    SaveBVector(sv11c_fname.c_str(), sv_11_c.bv_ares.get_bvector());
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
        const bm::bvector<> & bv = sv_offs_c.bv_ares.get_bvector();
        SaveBVector(svoffsc_fname.c_str(), bv);
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

    {
        std::string cbc_fname = base_name;
        cbc_fname.append("_bvcoll.cbc");
        int res = file_save_compressed_collection(buf_coll, cbc_fname);
        if (res != 0)
        {
            std::cerr << "compressed collection save error." << std::endl;
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

    {
        sv11c_fname = base_name;
        sv11c_fname.append("_sv11-c.bv");
        bm::bvector<>& bv_ares_11 = sv_11_c.bv_ares.get_bvector();
        LoadBVector(sv11c_fname.c_str(), bv_ares_11);

        sv_11_c.count_blocks();
    }

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

    {
        svoffsc_fname = base_name;
        svoffsc_fname.append("_svoffs-c.bv");
        bm::bvector<>& bv_ares_11 = sv_offs_c.bv_ares.get_bvector();

        LoadBVector(svoffsc_fname.c_str(), bv_ares_11);

        sv_offs_c.count_blocks();
    }
    }

    
    {
    std::string svsz_fname = base_name;
    svsz_fname.append("_svsize.sv");
    int res = file_load_svector(sv_sz, svsz_fname);
    if (res != 0)
    {
        std::cerr << "File load error." << std::endl;
    }
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

    // load buffer collection
    {
        std::string cbc_fname = base_name;
        cbc_fname.append("_bvcoll.cbc");
        int res = file_load_compressed_collection(buf_coll, cbc_fname);
        if (res != 0)
        {
            std::cerr << "compressed collection load error." << std::endl;
        }
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


extern "C" {
    int bit_visitor_func(void* handle_ptr, bm::id_t bit_idx)
    {
        std::vector<bm::id_t>* vp = (std::vector<bm::id_t>*)handle_ptr;
        vp->push_back(bit_idx);
        return 0;
    }
} // extern C



bool link_matrix::get_vector(unsigned id_from, std::vector<unsigned>& vect) const
{
    vect.resize(0);
/*
    if (!bv_from.test(id_from))
    {
        std::cerr << "id from not found" << std::endl;
        exit(1);
        vect.resize(0);
        return false;
    }
*/
    // check if vector is in the buffer storage
    //
    bm::id_t bv_addr;
    bool cbc_found = buf_coll.resolve(id_from, &bv_addr);
    if (cbc_found)
    {
        BM_DECLARE_TEMP_BLOCK(tb)

        const bm::compressed_buffer_collection<bvector_type>::buffer_type& bv_buf = buf_coll.get(bv_addr);
        bvector_type bv_tmp;
        bm::deserialize(bv_tmp, bv_buf.buf(), tb);
        bm::visit_each_bit(bv_tmp, (void*)&vect, bit_visitor_func);
        return true;
    }

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

   // std::cout << "Extract size = " << sz << std::endl;

    vect[0] = vc;
    sv_1m.extract(&vect[1], sz - 1, offs_c);


    //convert_from_delta(vect);

    return true;
}

bool link_matrix::combine_or(unsigned id_from, bm::bvector<>& bv) const
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
    // check if vector is in the buffer storage
    //
    bm::id_t bv_addr;
    bool cbc_found = buf_coll.resolve(id_from, &bv_addr);
    if (cbc_found)
    {
        BM_DECLARE_TEMP_BLOCK(tb)

        const bm::compressed_buffer_collection<bvector_type>::buffer_type& bv_buf = buf_coll.get(bv_addr);
        bm::operation_deserializer<bm::bvector<> > od;
        od.deserialize(bv, bv_buf.buf(), tb, bm::set_OR);
        return true;
    }

    unsigned vc;
    if (bv_11_flag.test(id_from))
    {
        vc = sv_11_c.get(id_from);
        bv.set_bit(vc);
        return true;
    }


    unsigned offs_c = sv_offs_c.get(id_from);
    unsigned sz = sv_sz.get(id_from);

    std::vector<unsigned> vect(sz);
    if (sz == 0)
    {
        std::cerr << "vector size error" << std::endl;
        exit(1);

        return false;
    }

    vc = sv_11_c.get(id_from);

    vect[0] = vc;
    sv_1m.extract(&vect[1], sz - 1, offs_c);

    bm::combine_or(bv, vect.begin(), vect.end());

    //convert_from_delta(vect);

    return true;
}



void link_matrix::add_vector(unsigned id_from, const bm::bvector<>& bv)
{
    std::vector<unsigned> vect;

    vect.resize(bv.count());
    std::copy(bv.first(), bv.end(), vect.begin());
    add_vector(id_from, vect);
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
    sv_sz.set(id_from, unsigned(sz));  // save true size
    sz--;
    sv_11.set(id_from, vect[0]);
    sv_1m.import(&vect[1], unsigned(sz), off);
    
    off += unsigned(sz);
}

link_matrix            link_m;
bm::chrono_taker<>::duration_map_type  timing_map;


// parse the input re-mapping vector
//
int load_bv(const std::string& fname, bm::bvector<>& bv)
{
    bm::chrono_taker tt1(cout, "1. parse input bit-vector", 1, &timing_map);

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

    std::vector<unsigned> vbuf;

    for (unsigned i = 0; std::getline(fin, line); i++)
    {
        if (line.empty())
            continue;

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



// load un-sorted pairs from a file
//
int load_ln_unsorted(const std::string& fname, link_matrix& lm)
{
    sm_accum   macc;

    std::string line;

    std::string regex_str = "[0-9]+";
    std::regex reg1(regex_str);

    std::ifstream fin(fname.c_str(), std::ios::in);
    if (!fin.good())
    {
        return -1;
    }


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
            case 0: id_from = (unsigned)ul; break;
            case 1: id_to = (unsigned)ul; break;
            default:
                break;
            }

        } // for

        macc.add_pair(id_from, id_to);


        if ((i != 0) && (i % 100000) == 0)
        {
            std::cout << "\r" << i << std::flush;
        }
        /*
        if ((i != 0) && (i % 20000000) == 0)
        {
            std::cout << "\r Optimization at " << i << std::endl;
            macc.optimize(false); // "fast" memory compression
            std::cout << "\r" << i << std::flush;
        }
        */
    } // for getline()

    macc.optimize(true); // full garbage collection

    std::cout << "Sort buffer loading ok. " << std::endl;
    
    lm.load_from_acc_matrix(macc);

    std::cout << "Sparse matrix loading ok. " << std::endl;

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
    bm::chrono_taker tt1(cout, "3. generation of remapping samples", benchmark_ops, &timing_map);
    
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
    bm::chrono_taker tt1(cout, "4. remapping", benchmark_ops, &timing_map);
    
        //std::vector<unsigned> vect;
        bm::bvector<>         bv_remap;
    
        for (unsigned i = 0; i < benchmark_ops; ++i)
        {
            const bm::bvector<>& bv = sample_vectors[i];
            
            bm::bvector<>::enumerator en = bv.first();
            for (unsigned k = 0; en.valid(); ++k, ++en)
            {
                unsigned id = *en;

                bool found = lm.combine_or(id, bv_remap);

                if (!found)
                {
                    std::cerr << "row not found! " << id << std::endl;
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

    {
        bm::chrono_taker tt1(cout, "5. input remap", 1, &timing_map);

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
            bm::chrono_taker tt1(cout, "0. link pairs load", 1, &timing_map);
            auto res = load_ln_unsorted(ln_in_file, link_m);
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
            bm::chrono_taker tt1(cout, "2. matrix load", 1, &timing_map);
            link_m.load(lm_in_name);
        }
        
        if (!lm_out_name.empty())
        {
            bm::chrono_taker tt1(cout, "1. matrix save", 1, &timing_map);
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
            std::cout << std::endl << "Performance (ops/sec):" << std::endl;
            bm::chrono_taker<>::print_duration_map(cout, timing_map, bm::chrono_taker<>::ct_ops_per_sec);
        }
        
        //getchar();
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



