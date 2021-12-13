/*
Copyright(c) 2020 Anatoliy Kuznetsov(anatoliy_kuznetsov at yahoo.com)

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
    no minimizers are used here (the approach can be used to create a
    scheme with minimizers).

    This example:
    - loads a FASTA file (single DNA molecule is expected)
    - generates K-mers for the specified K
      builds a fingerprint bit-vector (presence-absence vector)
    - runs k-mer counting (fastest method uses sort based algorithm)
      builds a k-mer frequency compressed sparse vector (TF-vector)
    - shows how to use TF vector to pick top N% of most frequent k-mers
      and exclude them (as over-represented)

    \sa bm::sparse_vector
    \sa bm::rsc_sparse_vector
    \sa bm::rank_range_split
    \sa bm::rsc_sparse_vector<>::merge_not_null

*/

/*! \file xsample07.cpp
    \brief Example: Use of bvector<> for k-mer fingerprint
    K should be short, no minimizers here
*/

#include <assert.h>
#include <stdlib.h>

#include <iostream>
#include <vector>
#include <map>
#include <algorithm>
#include <utility>

#include <future>
#include <thread>
#include <mutex>
#include <atomic>

#include "bm64.h"  // use 48-bit vectors
#include "bmalgo.h"
#include "bmserial.h"
#include "bmaggregator.h"
#include "bmsparsevec_compr.h"
#include "bmsparsevec_algo.h"
#include "bmundef.h" /* clear the pre-proc defines from BM */

// BitMagic utilities for debug and timings
#include "bmdbg.h"
#include "bmtimer.h"

#include "dna_finger.h"

using namespace std;



// Arguments
//
std::string  ifa_name;
std::string  ikd_name;
std::string  ikd_counts_name;
std::string  kh_name;
std::string  ikd_rep_name;
std::string  ikd_freq_name;
bool         is_diag = false;
bool         is_timing = false;
bool         is_bench = false;
unsigned     ik_size = 8;
unsigned     parallel_jobs = 4;
unsigned     f_percent = 5; // percent of k-mers we try to clear as over-represented

#include "cmd_args.h"



// Global types
//
typedef std::vector<char>                             vector_char_type;
typedef DNA_FingerprintScanner<bm::bvector<> >        dna_scanner_type;
typedef bm::sparse_vector<unsigned, bm::bvector<> >   sparse_vector_u32;
typedef bm::rsc_sparse_vector<unsigned, sparse_vector_u32 > rsc_sparse_vector_u32;
typedef std::map<unsigned, unsigned>                  histogram_map_u32;


// Global vars
//
bm::chrono_taker<>::duration_map_type     timing_map;
dna_scanner_type                          dna_scanner;
std::atomic_ullong                        k_mer_progress_count(0);


/// really simple FASTA parser (one entry per file)
///
static
int load_FASTA(const std::string& fname, vector_char_type& seq_vect)
{
    bm::chrono_taker tt1(cout, "1. Parse FASTA", 1, &timing_map);

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

inline
bool get_DNA_code(char bp, bm::id64_t& dna_code)
{
    switch (bp)
    {
    case 'A':
        dna_code = 0; // 00
        break;
    case 'T':
        dna_code = 1; // 01
        break;
    case 'G':
        dna_code = 2; // 10
        break;
    case 'C':
        dna_code = 3; // 11
        break;
    default: // ambiguity codes are ignored (for simplicity)
        return false;
    }
    return true;
}

/// Calculate k-mer as an unsigned long integer
///
/// @return true - if k-mer is "true" (not 'NNNNNN')
///
inline
bool get_kmer_code(const char* dna,
                  size_t pos, unsigned k_size,
                  bm::id64_t& k_mer)
{
    // generate k-mer
    //
    bm::id64_t k_acc = 0;
    unsigned shift = 0;
    dna += pos;
    for (size_t i = 0; i < k_size; ++i)
    {
        char bp = dna[i];
        bm::id64_t dna_code;
        bool valid = get_DNA_code(bp, dna_code);
        if (!valid)
            return false;
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
void translate_kmer(std::string& dna, bm::id64_t kmer_code, unsigned k_size)
{
    dna.resize(k_size);
    for (size_t i = 0; i < k_size; ++i)
    {
        unsigned dna_code = unsigned(kmer_code & 3);
        char bp = int2DNA(dna_code);
        dna[i] = bp;
        kmer_code >>= 2;
    } // for i
    assert(!kmer_code);
}


/// QA function to validate if reverse k-mer decode gives the same string
///
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

/// Auxiliary function to do sort+unique on a vactor of ints
/// removes duplicate elements
///
template<typename VECT>
void sort_unique(VECT& vect)
{
    std::sort(vect.begin(), vect.end());
    auto last = std::unique(vect.begin(), vect.end());
    vect.erase(last, vect.end());
}


/// Auxiliary function to do sort+unique on a vactor of ints and save results
/// in a counts vector
///
template<typename VECT, typename COUNT_VECT>
void sort_count(VECT& vect, COUNT_VECT& cvect)
{
    if (!vect.size())
        return;
    std::sort(vect.begin(), vect.end());
    typename VECT::value_type prev = vect[0];
    typename COUNT_VECT::value_type cnt = 1;
    auto vsize = vect.size();
    size_t i = 1;
    for (; i < vsize; ++i)
    {
        auto v = vect[i];
        if (v == prev)
        {
            ++cnt;
            continue;
        }
        cvect.inc_not_null(prev, cnt);
        prev = v; cnt = 1;
    } // for i

    cvect.inc_not_null(prev, cnt);
    assert(cvect.in_sync());
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
    bm::chrono_taker tt1(cout, "2. Generate k-mers", 1, &timing_map);

    bv.clear();
    bv.init(); // need to explicitly init to use bvector<>::set_bit_no_check()
    if (seq_vect.empty())
        return;
    const char* dna_str = &seq_vect[0];

    std::vector<bm::id64_t> k_buf;
    k_buf.reserve(chunk_size);

    {
        bm::id64_t k_mer_code;
        vector_char_type::size_type dna_sz = seq_vect.size()-(k_size-1);
        vector_char_type::size_type pos = 0;
        bool valid = false;
        for (; pos < dna_sz; ++pos)
        {
            valid = get_kmer_code(dna_str, pos, k_size, k_mer_code);
            if (valid)
            {
                k_buf.push_back(k_mer_code);
                break;
            }
        } // for pos

        const unsigned k_shift = (k_size-1) * 2;
        if (valid)
        {
            for (++pos; pos < dna_sz; ++pos)
            {
                bm::id64_t bp_code;
                valid = get_DNA_code(dna_str[pos + (k_size - 1)],  bp_code);
                if (!valid)
                {
                    pos += k_size; // wind fwrd to the next BP char
                    for (; pos < dna_sz; ++pos) // search for the next valid k-mer
                    {
                        valid = get_kmer_code(dna_str, pos, k_size, k_mer_code);
                        if (valid)
                        {
                            k_buf.push_back(k_mer_code);
                            break;
                        }
                    }
                    continue;
                }
                // shift out the previous base pair code, OR the new arrival
                k_mer_code = ((k_mer_code >> 2) | (bp_code << k_shift));

                // generated k-mer codes are accumulated in buffer for sorting
                k_buf.push_back(k_mer_code);

                if (check)
                {
                    validate_k_mer(dna_str, pos, k_size, k_mer_code);
                    bm::id64_t k_check;
                    valid = get_kmer_code(dna_str, pos, k_size, k_check);
                    assert(valid);
                    assert(k_check == k_mer_code);
                }

                if (k_buf.size() >= chunk_size) // sorting check.point
                {
                    sort_unique(k_buf);
                    if (k_buf.size())
                    {
                        bv.set(&k_buf[0], k_buf.size(), bm::BM_SORTED); // fast bulk set
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
        }

        if (k_buf.size()) // add last incomplete chunk here
        {
            sort_unique(k_buf);
            bv.set(&k_buf[0], k_buf.size(), bm::BM_SORTED); // fast bulk set

            cout << "Unique k-mers: " << k_buf.size() << endl;
        }
    }
    bv.optimize();
}

/// k-mer counting algorithm using reference sequence,
/// regenerates k-mer codes, sorts them and counts
///
inline
void count_kmers(const vector_char_type& seq_vect,
                 unsigned k_size,
                 rsc_sparse_vector_u32& kmer_counts)
{
    const bm::id64_t chunk_size = 400000000;
    if (seq_vect.empty())
        return;
    const char* dna_str = &seq_vect[0];
    std::vector<bm::id64_t> k_buf;
    k_buf.reserve(chunk_size);

    bm::id64_t k_mer_code;
    vector_char_type::size_type dna_sz = seq_vect.size()-(k_size-1);
    vector_char_type::size_type pos = 0;
    bool valid = false;
    for (; pos < dna_sz; ++pos)
    {
        valid = get_kmer_code(dna_str, pos, k_size, k_mer_code);
        if (valid)
        {
            k_buf.push_back(k_mer_code);
            break;
        }
    } // for pos
    const unsigned k_shift = (k_size-1) * 2;
    if (valid)
    {
        for (++pos; pos < dna_sz; ++pos)
        {
            bm::id64_t bp_code;
            valid = get_DNA_code(dna_str[pos + (k_size - 1)],  bp_code);
            if (!valid)
            {
                pos += k_size; // wind fwrd to the next BP char
                for (; pos < dna_sz; ++pos) // search for the next valid k-mer
                {
                    valid = get_kmer_code(dna_str, pos, k_size, k_mer_code);
                    if (valid)
                    {
                        k_buf.push_back(k_mer_code);
                        break;
                    }
                }
                continue;
            }
            // shift out the previous base pair code, OR the new arrival
            k_mer_code = ((k_mer_code >> 2) | (bp_code << k_shift));
            // generated k-mer codes are accumulated in buffer for sorting
            k_buf.push_back(k_mer_code);

            if (k_buf.size() == chunk_size) // sorting point
            {
                sort_count(k_buf, kmer_counts);
                k_buf.resize(0);

                float pcnt = float(pos) / float(dna_sz);
                pcnt *= 100;
                cout << "\r" << unsigned(pcnt) << "% of " << dna_sz
                     << " (" << (pos+1) <<")    "
                     << flush;
            }

        } // for pos
    }
    sort_count(k_buf, kmer_counts);
}

/// Functor to process job batch (task)
///
template<typename BV>
class SortCounting_JobFunctor
{
public:
    typedef BV                                 bvector_type;
    typedef typename bvector_type::size_type   size_type;

    /// constructor
    ///
    SortCounting_JobFunctor(const vector_char_type& seq_vect,
                            unsigned            k_size,
                            size_type           from,
                            size_type           to,
                            rsc_sparse_vector_u32& kmer_counts)
        : m_seq_vect(seq_vect), m_k_size(k_size), m_from(from), m_to(to),
          m_kmer_counts(kmer_counts)
    {}

    SortCounting_JobFunctor(const SortCounting_JobFunctor& func)
        : m_seq_vect(func.m_seq_vect), m_k_size(func.m_k_size),
          m_from(func.m_from), m_to(func.m_to),
          m_kmer_counts(func.m_kmer_counts)
    {}

    /// Main logic (functor)
    void operator() ()
    {
        const bvector_type* bv_null = m_kmer_counts.get_null_bvector();
        rsc_sparse_vector_u32 kmer_counts_part(*bv_null);
        kmer_counts_part.sync();

        const bm::id64_t chunk_size = 2000000;
        if (m_seq_vect.empty())
            return;

        const char* dna_str = &m_seq_vect[0];
        std::vector<bm::id64_t> k_buf;
        k_buf.reserve(chunk_size);

        const auto k_size = m_k_size;

        bm::id64_t k_mer_code(0);
        vector_char_type::size_type dna_sz = m_seq_vect.size()-(m_k_size-1);
        vector_char_type::size_type pos = 0;
        bool valid = false;
        for (; pos < dna_sz; ++pos)
        {
            valid = get_kmer_code(dna_str, pos, k_size, k_mer_code);
            if (valid)
            {
                if (k_mer_code >= m_from && k_mer_code <= m_to)
                    k_buf.push_back(k_mer_code);
                break;
            }
        } // for pos

        const unsigned k_shift = (k_size-1) * 2;
        if (valid)
        {
            for (++pos; pos < dna_sz; ++pos)
            {
                bm::id64_t bp_code;
                valid = get_DNA_code(dna_str[pos + (k_size - 1)],  bp_code);
                if (!valid)
                {
                    pos += k_size; // wind fwrd to the next BP char
                    for (; pos < dna_sz; ++pos) // search for the next valid k-mer
                    {
                        valid = get_kmer_code(dna_str, pos, k_size, k_mer_code);
                        if (valid)
                        {
                            if (k_mer_code >= m_from && k_mer_code <= m_to)
                                k_buf.push_back(k_mer_code);
                            break;
                        }
                    }
                    continue;
                }
                // shift out the previous base pair code, OR the new arrival
                k_mer_code = ((k_mer_code >> 2) | (bp_code << k_shift));
                // generated k-mer codes are accumulated in buffer for sorting
                if (k_mer_code >= m_from && k_mer_code <= m_to)
                {
                    k_buf.push_back(k_mer_code);
                    if (k_buf.size() == chunk_size) // sorting point
                    {
                        sort_count(k_buf, kmer_counts_part);
                        k_buf.resize(0);
                    }
                }
            } // for pos
        }

        sort_count(k_buf, kmer_counts_part);

        // merge results
        {
            static std::mutex   mtx_counts_lock;
            std::lock_guard<std::mutex> guard(mtx_counts_lock);

            // merge_not_null fast merge for non-overlapping ranges of
            // rsc_sparse_vector
            m_kmer_counts.merge_not_null(kmer_counts_part);
        }
    }

private:
    const vector_char_type&   m_seq_vect;
    unsigned                  m_k_size;
    size_type                 m_from;
    size_type                 m_to;
    rsc_sparse_vector_u32&    m_kmer_counts;
};

/// MT k-mer counting
///
template<typename BV>
void count_kmers_parallel(const BV& bv_kmers,
                          const vector_char_type& seq_vect,
                          rsc_sparse_vector_u32& kmer_counts,
                          unsigned k_size,
                          unsigned concurrency)
{
    typedef typename BV::size_type bv_size_type;
    typedef std::vector<std::pair<bv_size_type, bv_size_type> > bv_ranges_vector;

    bv_ranges_vector pair_vect;

    bv_size_type cnt = bv_kmers.count();
    bv_size_type split_rank = cnt / concurrency; // target population count per job

    if (split_rank < concurrency || concurrency < 2)
    {
        count_kmers(seq_vect, k_size, kmer_counts); // run single threaded
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
            SortCounting_JobFunctor<BV>(seq_vect, ik_size,
                                pair_vect[k].first,
                                pair_vect[k].second,
                                kmer_counts)));
    } // for k

    // wait for all jobs to finish, print progress report
    //
    for (auto& e : futures)
    {
        unsigned m = 0;
        while(1)
        {
            std::future_status status = e.wait_for(std::chrono::seconds(60));
            if (status == std::future_status::ready)
                break;
            cout << "\r" << ++m << " min" << flush;
        } // while
    } // for
    cout << endl;
}


/// k-mer counting method using Bitap algorithm for occurence search
///  this method is significantly slower than direct regeneration of k-mer
///  codes and sorting count
///
template<typename BV>
void count_kmers(const BV& bv_kmers, rsc_sparse_vector_u32& kmer_counts)
{
    std::string kmer_str;
    typename BV::size_type cnt = 0; // progress report counter
    typename BV::enumerator en = bv_kmers.first();
    for ( ;en.valid(); ++en, ++cnt)
    {
        auto k_mer_code = *en;
        translate_kmer(kmer_str, k_mer_code, ik_size); // translate k-mer code to string

        // find list of sequence positions where k-mer is found
        // (uncomment if you need a full search)
        /*
        typename BV::size_type bv_count;
        {
            std::vector<typename BV::size_type> km_search;
            dna_scanner.Find(kmer_str, km_search);
            bv_count = km_search.size();
        }
        */

        typename BV::size_type count = dna_scanner.FindCount(kmer_str);
        assert(count);
        kmer_counts.set(k_mer_code, (unsigned)count);

        if ((cnt % 1000) == 0)
            cout << "\r" << cnt << flush;

    } // for en
}

/**
    k-mer counting job functor class using bm::aggregator<>

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
                        size_type           to,
                        rsc_sparse_vector_u32& kmer_counts)
        : m_parent_scanner(parent_scanner), m_bv_kmers(bv_kmers),
          m_from(from), m_to(to), m_kmer_counts(kmer_counts)
    {}

    /// copy-ctor
    Counting_JobFunctor(const Counting_JobFunctor& func)
        : m_parent_scanner(func.m_parent_scanner), m_bv_kmers(func.m_bv_kmers),
          m_from(func.m_from), m_to(func.m_to), m_kmer_counts(func.m_kmer_counts)
    {}

    /// Main logic (functor)
    void operator() ()
    {
        std::string kmer_str;
        size_type cnt = 0; // progress report counter

        static std::mutex   mtx_counts_lock;
        bvector_type        bv_search_res;

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

            // counts are shared across threads, use locked access
            // to save the results
            // TODO: implement results buffering to avoid mutex overhead
            {
                std::lock_guard<std::mutex> guard(mtx_counts_lock);
                m_kmer_counts.set(k_mer_code, unsigned(search_count));
                assert(m_kmer_counts.in_sync());
            }

            k_mer_progress_count.fetch_add(1);

        } // for en
    }

private:
    const DNA_Scan&                    m_parent_scanner;
    const bvector_type&                m_bv_kmers;
    size_type                          m_from;
    size_type                          m_to;
    typename DNA_Scan::aggregator_type m_Agg;

    rsc_sparse_vector_u32&             m_kmer_counts;
};

/**
    Runs k-mer counting in parallel
*/
template<typename BV>
void count_kmers_parallel(const BV& bv_kmers,
                          rsc_sparse_vector_u32& kmer_counts,
                          unsigned concurrency)
{
    typedef typename BV::size_type bv_size_type;
    typedef std::vector<std::pair<bv_size_type, bv_size_type> > bv_ranges_vector;

    bv_ranges_vector pair_vect;

    bv_size_type cnt = bv_kmers.count();
    bv_size_type split_rank = cnt / concurrency; // target population count per job

    if (split_rank < concurrency || concurrency < 2)
    {
        count_kmers(bv_kmers, kmer_counts); // run single threaded
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
                                pair_vect[k].second,
                                kmer_counts)));
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


/// Compute a map of how often each k-mer frequency is observed in the
/// k-mer counts vector
/// @param hmap - [out] histogram map
/// @param kmer_counts - [in] kmer counts vector
///
static
void compute_kmer_histogram(histogram_map_u32& hmap,
                            const rsc_sparse_vector_u32& kmer_counts)
{
    const rsc_sparse_vector_u32::bvector_type* bv_null =
                                    kmer_counts.get_null_bvector();
    auto en = bv_null->first();
    for (; en.valid(); ++en)
    {
        auto kmer_code = *en;
        auto kmer_count = kmer_counts.get(kmer_code);
        auto mit = hmap.find(kmer_count);
        if (mit == hmap.end())
            hmap[kmer_count] = 1;
        else
            mit->second++;
    } // for
}

/// Save TSV report of k-mer frequences
/// (reverse sorted, most frequent k-mers first)
///
static
void report_hmap(const string& fname, const histogram_map_u32& hmap)
{
    ofstream outf;
    outf.open(fname, ios::out | ios::trunc );

    outf << "kmer count \t number of kmers\n";

    auto it = hmap.rbegin(); auto it_end = hmap.rend();
    for (; it != it_end; ++it)
    {
        outf << it->first << "\t" << it->second << endl;
    }
}

/// Create vector, representing subset of k-mers of high frequency
///
/// @param frequent_bv[out] - bit-vector of frequent k-mers (subset of all k-mers)
/// @param hmap - histogram map of all k-mers
/// @param kmer_counts - kmer frequency(counts) vector
/// @param percent - percent of frequent k-mers to build a subset (5%)
///   percent here is of total number of k-mers (not percent of all occurences)
/// @param k_size - K mer size
///
template<typename BV>
void compute_frequent_kmers(BV& frequent_bv,
                            const histogram_map_u32& hmap,
                            const rsc_sparse_vector_u32& kmer_counts,
                            unsigned percent,
                            unsigned k_size)
{
    (void)k_size;
    frequent_bv.clear();

    if (!percent)
        return;

    // scanner class for fast search for values in sparse vector
    bm::sparse_vector_scanner<rsc_sparse_vector_u32> scanner;
    BV bv_found;  // search results vector

    const rsc_sparse_vector_u32::bvector_type* bv_null =
                                    kmer_counts.get_null_bvector();

    auto total_kmers = bv_null->count();
    bm::id64_t target_f_count = (total_kmers * percent) / 100; // how many frequent k-mers we need to pick

    auto it = hmap.rbegin();
    auto it_end = hmap.rend();
    for (; it != it_end; ++it)
    {
        auto kmer_count = it->first;

        scanner.find_eq(kmer_counts, kmer_count, bv_found); // seach for all values == 25
        auto found_cnt = bv_found.count();
        (void)found_cnt;
        assert(found_cnt);
        {
            bm::bvector<>::enumerator en = bv_found.first();
            for (; en.valid(); ++en)
            {
                auto kmer_code = *en;
                unsigned k_count = kmer_counts.get(kmer_code);

                if (k_count == 1) // unique k-mer, ignore
                    continue;

                if (it->second == 1)
                {
                    assert(k_count ==  kmer_count);
                }
                frequent_bv.set(kmer_code);
                if (kmer_count >= target_f_count)
                {
                    frequent_bv.optimize();
                    return;
                }
                target_f_count -= 1;

            } // for en
        }
    } // for it
    frequent_bv.optimize();
}


int main(int argc, char *argv[])
{
    vector_char_type       seq_vect; // read FASTA sequence
    bm::bvector<>          bv_kmers; // k-mer presense(-absence) vector

    try
    {
        auto ret = parse_args(argc, argv);
        if (ret != 0)
        {
            cerr << "cmd-line parse error. " << endl;
            return ret;
        }

        cout << "concurrency=" << parallel_jobs << endl;

        if (!ifa_name.empty()) // FASTA file load
        {
            // limitation: loads a single molecule only
            //
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
                bm::print_bvector_stat(cout,bv_kmers);
                size_t blob_size = bm::compute_serialization_size(bv_kmers);
                cout << "DNA 2-bit coded FASTA size=" << seq_vect.size()/4 << endl;
                cout << "           Compressed size=" << blob_size << endl;
            }
        }

        if (!ikd_name.empty())
        {
            bm::chrono_taker tt1(cout, "3. k-mer serialization and save", 1, &timing_map);

            bm::SaveBVector(ikd_name.c_str(), bv_kmers);
        }

        if (seq_vect.size())
        {
            bm::chrono_taker tt1(cout, "4. Build DNA fingerprints (bulk, parallel)", 1, &timing_map);
            dna_scanner.BuildParallel(seq_vect, parallel_jobs);
        }

        if (seq_vect.size() &&
           (!ikd_counts_name.empty() || !ikd_rep_name.empty()))
        {
            rsc_sparse_vector_u32  rsc_kmer_counts(bv_kmers); // rank-select sparse vector for counts
            rsc_kmer_counts.sync();
            rsc_sparse_vector_u32  rsc_kmer_counts2(bv_kmers);
            rsc_kmer_counts2.sync();


            cout << " Searching for k-mer counts..." << endl;

            if (is_diag) // compute reference counts using slower algorithm for verification
            {
                cout << " ... using bm::aggregator<>" << endl;
                bm::chrono_taker tt1(cout, "5a. k-mer counting (bm::aggregator<>)", 1, &timing_map);

                count_kmers_parallel(bv_kmers, rsc_kmer_counts, parallel_jobs);
                //count_kmers(bv_kmers, rsc_kmer_counts);

                rsc_kmer_counts.optimize();
            }

            {
                cout << " ... using std::sort() and count" << endl;
                bm::chrono_taker tt1(cout, "5. k-mer counting std::sort()", 1, &timing_map);

                count_kmers_parallel(bv_kmers, seq_vect,
                                     rsc_kmer_counts2, ik_size, parallel_jobs);

                rsc_kmer_counts2.optimize();

                if (!ikd_counts_name.empty())
                {
                    int res = bm::file_save_svector(rsc_kmer_counts2, ikd_counts_name);
                    if (res)
                    {
                        cerr << "Error: Count vector save failed!" << endl;
                        exit(1);
                    }
                }

                if (is_diag) // verification
                {
                    bool eq = rsc_kmer_counts.equal(rsc_kmer_counts2);
                    if(!eq)
                    {
                        rsc_sparse_vector_u32::size_type idx;
                        bool found = bm::sparse_vector_find_first_mismatch(rsc_kmer_counts,
                                                              rsc_kmer_counts2,
                                                              idx);
                        auto v1 = rsc_kmer_counts.get(idx);
                        auto v2 = rsc_kmer_counts2.get(idx);
                        cerr << "Mismatch at idx=" << idx << " v1=" << v1 << " v2=" << v2 << endl;
                        assert(found); (void)found;
                        assert(eq);
                        cerr << "Integrity check failed!" << endl;
                        exit(1);
                    }
                }
            }

            // build histogram of k-mer counts
            // as a map of kmer_count -> number of occurences of k-mer
            //

            histogram_map_u32 hmap;
            {
                bm::chrono_taker tt1(cout, "6. build histogram of k-mer frequencies", 1, &timing_map);
                compute_kmer_histogram(hmap, rsc_kmer_counts2);
            }
            if (!kh_name.empty())
            {
                report_hmap(kh_name, hmap);
            }

            // here we build a build-vector of frequent (say top 5%) k-mers
            // (if needed we can exclude it because they likely represent repeats
            //
            bm::bvector<> bv_freq(bm::BM_GAP);
            {
                bm::chrono_taker tt1(cout, "7. Build vector of frequent k-mers", 1, &timing_map);
                compute_frequent_kmers(bv_freq, hmap,
                                       rsc_kmer_counts2, f_percent, ik_size);
            }
            if (!ikd_freq_name.empty())
            {
                bm::SaveBVector(ikd_freq_name.c_str(), bv_freq);
            }

            cout << "Found frequent k-mers: " << bv_freq.count() << endl;

            if (!ikd_rep_name.empty())
            {
                // exclude frequent k-mers (logical SUBtraction)
                //
                bv_kmers.bit_sub(bv_freq);
                bv_kmers.optimize();

                bm::SaveBVector(ikd_rep_name.c_str(), bv_kmers);
            }
        }

        if (is_timing)
        {
            std::cout << std::endl << "Performance:" << std::endl;
            bm::chrono_taker<>::print_duration_map(cout, timing_map, bm::chrono_taker<>::ct_time);
        }

    }
    catch(std::exception& ex)
    {
        std::cerr << ex.what() << std::endl;
        return 1;
    }
    
    

    return 0;
}

