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
#include <stdatomic.h>

//#include "bm.h"
#include "bm64.h"  // use 48-bit vectors
#include "bmalgo.h"
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



// Global types
//
//typedef bm::sparse_vector<bm::id64_t, bm::bvector<> > svector_u64;
typedef std::vector<char>                             vector_char_type;
typedef DNA_FingerprintScanner<bm::bvector<> >        dna_scanner_type;


// Globals
//
bm::chrono_taker::duration_map_type     timing_map;
dna_scanner_type                        dna_scanner;
std::atomic_ullong                      k_mer_progress_count(0);


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

/**
    This function turns each k-mer into an integer number and encodes it
    in a bit-vector (presense vector)
    The natural limitation here is that integer has to be less tha 48-bits
    (limitations of bm::bvector<>)
    This method build a presense k-mer fingerprint vector which can be
    used for Jaccard distance comparison.

    @param bv - [out] - target bit-vector
    @param seq_vect - [out] DNA sequence vector
    @param k-size   - dimention for k-mer generation
 */
template<typename BV>
void generate_k_mer_bvector(BV& bv,
                            const vector_char_type& seq_vect,
                            unsigned k_size,
                            bool check)
{
    const bm::id64_t chunk_size = 400000000;
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
        typename BV::bulk_insert_iterator bit(bv);
        vector_char_type::size_type dna_sz = seq_vect.size();
        for (vector_char_type::size_type pos = 0; pos < dna_sz; ++pos)
        {
            bm::id64_t k_mer_code;
            bool valid = get_kmer_code(dna_str, pos, k_size, k_mer_code);
            if (!valid)
                continue;
            ++k_cnt;

            if (check)
                validate_k_mer(dna_str, pos, k_size, k_mer_code);

            // generated k-mer codes are accumulated in buffer for sorting
            k_buf.push_back(k_mer_code);
            if (k_buf.size() == chunk_size) // soring point
            {
                sort_unique(k_buf);
                if (k_buf.size())
                {
                    bv.set(&k_buf[0], k_buf.size(), bm::BM_SORTED);
                    k_buf.resize(0);
                    bv.optimize(); // periodically re-optimize to save memory
                }

                float pcnt = float(pos) / float(dna_sz);
                pcnt *= 100;
                cout << "\r" << unsigned(pcnt) << "% of " << dna_sz
                     << " (" << (pos+1) <<")    "
                     << flush;
            }
        } // for pos

        if (k_buf.size()) // add last incomplete chunk here
        {
            sort_unique(k_buf);
            bv.set(&k_buf[0], k_buf.size(), bm::BM_SORTED);

            cout << "Unique k-mers: " << k_buf.size() << endl;
        }
        bit.flush();
    }
    bv.optimize();
    cout << "\rValid k-mers:" << k_cnt << "                         " << endl;
}

template<typename BV>
void count_kmers(const BV& bv_kmers)
{
    std::string kmer_str;
    typename BV::size_type cnt = 0; // progress report counter
    typename BV::enumerator en = bv_kmers.first();
    for ( ;en.valid(); ++en, ++cnt)
    {
        auto k_mer_code = *en;
        translate_kmer(kmer_str, k_mer_code, ik_size); // translate k-mer code to string

        // find list of sequence positions where k-mer is found
        //
        //dna_scanner.Find(kmer_str, km_search);

        typename BV::size_type count = dna_scanner.FindCount(kmer_str);
        (void) count;
/*
        if (count != km_search.size())
        {
            cerr << "Count failure at: " << cnt << endl;
            cerr << count << " " << km_search.size() << endl;
            exit(1);
        }
*/
        if ((cnt % 1000) == 0)
        {
            cout << "\r" << cnt << flush;
        }

    } // for en
}

/**
    k-mer counting job functor class

    Functor received its range of k-mers in the presence-absense
    bit-vector then follows it to run the search-counting algorithm
    using DNA fingerprints common for all job functors.

    bm::aggregator<> cannot be shared across threads,
    so functor creates its own
 */
template<typename DNA_Scan>
class Counting_JobFunctor
{
public:
    typedef typename DNA_Scan::bvector_type  bvector_type;
    typedef typename bvector_type::size_type size_type;

    /// constructor
    ///
    Counting_JobFunctor(const DNA_Scan&     parent_scanner,
                        const bvector_type& bv_kmers,
                        size_type           from,
                        size_type           to)
        : m_parent_scanner(parent_scanner), m_bv_kmers(bv_kmers),
          m_from(from), m_to(to)
    {}

    /// copy-ctor
    Counting_JobFunctor(const Counting_JobFunctor& func)
        : m_parent_scanner(func.m_parent_scanner), m_bv_kmers(func.m_bv_kmers),
          m_from(func.m_from), m_to(func.m_to)
    {}

    // Main functor method
    void operator() ()
    {
        std::string kmer_str;
        size_type cnt = 0; // progress report counter

        bvector_type bv_search_res;

        typename bvector_type::enumerator en = m_bv_kmers.get_enumerator(m_from);
        for ( ;en.valid(); ++en, ++cnt)
        {
            auto k_mer_code = *en;
            if (k_mer_code > m_to)
                break;
            translate_kmer(kmer_str, k_mer_code, ik_size); // translate k-mer code to string

            // setup the aggregator to perform search
            //
            m_Agg.reset();
            m_Agg.set_compute_count(true); // disable full search, only count
            for (size_t i = 0; i < kmer_str.size(); ++i)
            {
                const bvector_type& bv_mask = m_parent_scanner.GetVector(kmer_str[i]);
                m_Agg.add(&bv_mask);
            }

            m_Agg.combine_shift_right_and(bv_search_res);

            // Note we get search count from the Aggregator, not from search
            // result vector, which will be empty,
            // because we set_compute_count(true)
            //
            size_type search_count = m_Agg.count();
            (void) search_count;

            k_mer_progress_count.fetch_add(1);

        } // for en
    }

private:
    const DNA_Scan&                    m_parent_scanner;
    const bvector_type&                m_bv_kmers;
    size_type                          m_from;
    size_type                          m_to;
    typename DNA_Scan::aggregator_type m_Agg;
};

template<typename BV>
void count_kmers_parallel(const BV& bv_kmers, unsigned concurrency)
{
    typedef typename BV::size_type bv_size_type;
    typedef std::vector<std::pair<bv_size_type, bv_size_type> > bv_ranges_vector;

    bv_ranges_vector pair_vect;

    bv_size_type cnt = bv_kmers.count();
    bv_size_type split_rank = cnt / concurrency; // target population count per job

    if (split_rank < concurrency || concurrency < 2)
    {
        count_kmers(bv_kmers); // run single threaded
        return;
    }

    // run split algorithm to determine equal weight ranges for parallel
    // processing
    bm::rank_range_split(bv_kmers, split_rank, pair_vect);

    // Create parallel async tasks running on a range of source sequence
    //
    std::vector<std::future<void> > futures;
    futures.reserve(concurrency);

    for (size_t k = 0; k < pair_vect.size(); ++k)
    {
        futures.emplace_back(std::async(std::launch::async,
            Counting_JobFunctor<dna_scanner_type>(dna_scanner, bv_kmers,
                                pair_vect[k].first,
                                pair_vect[k].second)));
    } // for k

    // wait for all jobs to finish, print progress report
    //
    for (auto& e : futures)
    {
        unsigned long long c_prev = 0;
        while(1)
        {
            std::future_status status = e.wait_for(std::chrono::seconds(60));
            if (status == std::future_status::ready)
                break;

            // progress report (entertainment)
            //
            unsigned long long c = k_mer_progress_count;
            auto delta = c - c_prev;
            c_prev = c;

            auto remain_cnt = cnt - c;
            auto remain_min = remain_cnt / delta;
            cout << "\r" << c << ": progress per minute=" << delta;
            if (remain_min < 120)
            {
                 cout << " wait for " << remain_min << "m     " << flush;
            }
            else
            {
                auto remain_h = remain_min / 60;
                cout << " wait for " << remain_h << "h     " << flush;
            }
        } // while

    } // for
    cout << endl;

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

        cout << "concurrency=" << parallel_jobs << endl;

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
        }


        if (seq_vect.size())
        {
            cout << " Searching for k-mer counts..." << endl;
            bm::chrono_taker tt1("5. k-mer counting", 1, &timing_map);

            count_kmers_parallel(bv_kmers, parallel_jobs);

#if 0
            std::string kmer_str;
            std::vector<bm::bvector<>::size_type> km_search;
            bm::bvector<>::size_type cnt = 0;
            bm::bvector<>::enumerator en = bv_kmers.first();
            for ( ;en.valid(); ++en, ++cnt)
            {
                auto k_mer_code = *en;
                translate_kmer(kmer_str, k_mer_code, ik_size);

                // find list of sequence positions where k-mer is found
                //
                dna_scanner.Find(kmer_str, km_search);

                bm::bvector<>::size_type count = dna_scanner.FindCount(kmer_str);

                if (count != km_search.size())
                {
                    cerr << "Count failure at: " << cnt << endl;
                    cerr << count << " " << km_search.size() << endl;
                    exit(1);
                }

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
            } // for en
#endif
        }


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

