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

/** \example xsample07.cpp
    Use of bvector<> for k-mer fingerprint K should be short,
    no minimizers here
*/

/*! \file xsample07.cpp
    \brief Example: Use of bvector<> for k-mer fingerprint
    K should be short, no minimizers here
*/

#include <assert.h>
#include <stdlib.h>

#include <iostream>
#include <vector>
#include <algorithm>
#include <utility>

#include <future>
#include <thread>
#include <mutex>


//#include "bm.h"
#include "bm64.h"  // use 48-bit vectors
#include "bmserial.h"
#include "bmaggregator.h"

// BitMagic utilities for debug and timings
#include "bmdbg.h"
#include "bmtimer.h"

#include "dna_finger.h"

using namespace std;



// Arguments
//
std::string  ifa_name;
std::string  ikd_name;
bool         is_diag = false;
bool         is_timing = false;
bool         is_bench = false;
unsigned     ik_size = 8;
unsigned     parallel_jobs = 4;

#include "cmd_args.h"


// Globals
//
bm::chrono_taker::duration_map_type  timing_map;
DNA_FingerprintScanner               dna_scanner;

// Global types
//
//typedef bm::sparse_vector<bm::id64_t, bm::bvector<> > svector_u64;
typedef std::vector<char>                             vector_char_type;




/// really simple FASTA file parser
///
static
int load_FASTA(const std::string& fname, vector_char_type& seq_vect)
{
    bm::chrono_taker tt1("1. Parse FASTA", 1, &timing_map);

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
            seq_vect.push_back(*it);
    } // for
    return 0;
}

/// Translate DNA letter to integer code
///
inline
unsigned DNA2int(char DNA_bp)
{
    switch (DNA_bp)
    {
    case 'A':
        return 0; // 00
    case 'T':
        return 1; // 01
    case 'G':
        return 2; // 10
    case 'C':
        return 3; // 11
    default:
        return 0;
    }
}


/// Calculate k-mer as an unsigned long integer
///
///
/// @return true - if k-mer is "true" (not 'NNNNNN')
///
inline
bool get_kmer_code(const char* dna,
                  size_t pos, unsigned k_size,
                  bm::id64_t& k_mer)
{
    assert(k_size > 2 && k_size <= 32);

    // generate k-mer
    //
    bm::id64_t k_acc = 0;
    unsigned shift = 0;
    for (size_t i = 0; i < k_size; ++i)
    {
        char bp = dna[pos+i];
        if (!bp) // boundary: string terminated - ignore short k-mer
            return false;
        if (bp == 'N')
            return false; // 'N' containing k-mers are ignored (for simplicity)
        bm::id64_t dna_code = DNA2int(bp);
        k_acc |= (dna_code << shift); // accumulate new code within 64-bit accum
        shift += 2; // each DNA base pair needs 2-bits to store
    } // for i
    k_mer = k_acc;
    return true;
}


/// Translate integer code to DNA letter
///
inline
char int2DNA(unsigned code)
{
    static char lut[] = { 'A', 'T', 'G', 'C', 'N', '$' };
    if (code < 5)
        return lut[code];
    assert(0);
    return 'N';
}

/// Translate k-mer code into ATGC DNA string
///
/// @param dna    - target string
/// @param k_mer  - k-mer code
/// @param k_size -
inline
void translate_kmer(std::string& dna, bm::id64_t k_mer, unsigned k_size)
{
    dna.resize(k_size);
    for (size_t i = 0; i < k_size; ++i)
    {
        unsigned dna_code = unsigned(k_mer & 3);
        char bp = int2DNA(dna_code);
        dna[i] = bp;
        k_mer >>= 2;
    } // for i
    assert(!k_mer);
}




/// QA function to validate if reverse k-mer decode gives the same string
inline
void validate_k_mer(const char* dna,
                    size_t pos, unsigned k_size,
                    bm::id64_t k_mer)
{
    for (size_t i = 0; i < k_size; ++i)
    {
        char bp = dna[pos+i];
        unsigned dna_code = unsigned(k_mer & 3ull);
        char bp_c = int2DNA(dna_code);
        if (bp != bp_c)
        {
            if (bp == 'N' && bp_c != 'A')
            {
                cerr << bp << " " << bp_c << endl;
                cerr << "Error! N code mismatch at pos = " << pos+i
                     << endl;
                exit(1);
            }
        }
        k_mer >>= 2;
    } // for i
    if (k_mer)
    {
        cerr << "Error! non-zero k-mer remainder at  " << pos << endl;
        exit(1);
    }
}

template<class VECT>
void sort_unique(VECT& vect)
{
    std::sort(vect.begin(), vect.end());
    auto last = std::unique(vect.begin(), vect.end());
    vect.erase(last, vect.end());
}


template<typename BV>
void generate_k_mer_bvector(BV& bv,
                            const vector_char_type& seq_vect,
                            unsigned k_size,
                            bool check)
{
    const bm::id64_t chunk_size = 400000000;
//cerr << "buffer memory consumption:" << chunk_size * sizeof(bm::id64_t);
    bm::chrono_taker tt1("2. Generate k-mers", 1, &timing_map);

    bv.clear();
    bv.init(); // need to explicitly init to use bvector<>::set_bit_no_check()
    if (seq_vect.empty())
        return;
    const char* dna_str = &seq_vect[0];

    vector_char_type::size_type k_cnt = 0;
    std::vector<bm::id64_t> k_buf;
    k_buf.reserve(chunk_size);

    {
        //vector_char_type::size_type opt_cnt = chunk_size;
        typename BV::bulk_insert_iterator bit(bv);
        vector_char_type::size_type dna_sz = seq_vect.size();
        for (vector_char_type::size_type pos = 0; pos < dna_sz; ++pos)
        {
            bm::id64_t k_mer_code;
            bool valid = get_kmer_code(dna_str, pos, k_size, k_mer_code);
            if (valid)
            {
                ++k_cnt;
                if (check)
                    validate_k_mer(dna_str, pos, k_size, k_mer_code);
                //bit = k_mer_code;
                //bv.set_bit_no_check(k_mer_code);
                k_buf.push_back(k_mer_code);

                //if (!(--opt_cnt))
                /*
                if (k_buf.size() == chunk_size)
                {
                    sort_unique(k_buf);
                    if (k_buf.size())
                    {
                        //bv.set(&k_buf[0], k_buf.size(), bm::BM_SORTED);
                        k_buf.resize(0);
                        //bv.optimize(); // periodically re-optimize to save memory
                    }

                    //opt_cnt = chunk_size;
                    float pcnt = float(pos) / float(dna_sz);
                    pcnt *= 100;
                    cout << "\r" << unsigned(pcnt) << "% of " << dna_sz
                         << " (" << (pos+1) <<")    "
                         << flush;
                }
                */
            }
        } // for pos

        if (k_buf.size())
        {
            sort_unique(k_buf);
            cout << "Unique k-mers: " << k_buf.size() << endl;

        }

        //bit.flush();
    }
    bv.optimize();
    cout << "\r valid k-mers:" << k_cnt << "                         " << endl;
}



int main(int argc, char *argv[])
{
    vector_char_type  seq_vect; // read FASTA sequence
    bm::bvector<>     bv_kmers; //
    //svector_u64       sv_kmers;

    try
    {
        auto ret = parse_args(argc, argv);
        if (ret != 0)
            return ret;

        if (!ifa_name.empty()) // FASTA file load
        {
            auto res = load_FASTA(ifa_name, seq_vect);
            if (res != 0)
                return res;
            std::cout << "FASTA sequence size=" << seq_vect.size() << std::endl;
        }

        if (seq_vect.size())
        {
            cout << "k-mer generation for k=" << ik_size << endl;

            generate_k_mer_bvector(bv_kmers, seq_vect, ik_size, is_diag);

            cout << "Found " << bv_kmers.count() << " k-mers." << endl;

            if (is_diag)
            {
                bm::print_bvector_stat(bv_kmers);
                size_t blob_size = bm::compute_serialization_size(bv_kmers);
                cout << "DNA 2-bit coded FASTA size=" << seq_vect.size()/4 << endl;
                cout << "           Compressed size=" << blob_size << endl;
            }
        }

        if (!ikd_name.empty())
        {
            bm::chrono_taker tt1("3. k-mer serialization and save", 1, &timing_map);

            bm::SaveBVector(ikd_name.c_str(), bv_kmers);
        }

        if (seq_vect.size())
        {
            bm::chrono_taker tt1("4. Build DNA fingerprints (bulk, parallel)", 1, &timing_map);
            dna_scanner.BuildParallel(seq_vect, parallel_jobs);

            //dna_scanner.Build(seq_vect);
        }

#if 0
        if (seq_vect.size())
        {
            cout << " Searchging for k-mers..." << endl;
            bm::chrono_taker tt1("5. Fingerprint search", 1, &timing_map);

            std::string kmer_str;
            std::vector<bm::bvector<>::size_type> km_search;
            //std::vector<bm::bvector<>::size_type> km_search_control;
            bm::bvector<>::size_type cnt = 0;
            bm::bvector<>::enumerator en = bv_kmers.first();
            for ( ;en.valid(); ++en, ++cnt)
            {
                auto k_mer_code = *en;
                translate_kmer(kmer_str, k_mer_code, ik_size);

                // find list of sequence positions where k-mer is found
                //
                dna_scanner.FindAggFused(kmer_str, km_search);

//                dna_scanner.Find(kmer_str, km_search);
/*
                if (!km_search.size())
                {
                    bool b = bv_kmers.test(k_mer_code);
                    if (!b)
                    {
                        cout << "enumerator test failure!" << endl;
                    }

                    cerr << "Error! incorrect search " << k_mer_code
                         << " " << kmer_str << endl;
                    dna_scanner.Find(kmer_str, km_search_control);
                    cout << "control size=" << km_search_control.size() << endl;

                    const char* s = &seq_vect[0];
                    const char* pch = ::strstr(s, kmer_str.c_str());
                    if (!pch)
                    {
                        cout << "::strstr() not found ..." << endl;
                    }
                    else
                    {
                        auto pos = pch - s;
                        cout << "strstr() found at:" << pos << endl;
                    }
                    exit(1);
                }
*/
                if ((cnt % 1000) == 0)
                {
                    cout << "\r" << cnt << flush;
                }
                if (cnt == 5000)
                    break;
            } // for en
        }
#endif

        if (is_timing)
        {
            std::cout << std::endl << "Performance:" << std::endl;
            bm::chrono_taker::print_duration_map(timing_map, bm::chrono_taker::ct_time);
        }

    }
    catch(std::exception& ex)
    {
        std::cerr << ex.what() << std::endl;
        return 1;
    }
    
    

    return 0;
}

