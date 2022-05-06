/*
Copyright(c) 2002-2020 Anatoliy Kuznetsov(anatoliy_kuznetsov at yahoo.com)

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

/** \example xsample07a.cpp
    Use of bvector<> for k-mer fingerprint K should be short,
    no minimizers here
*/

/*! \file xsample07a.cpp
    \brief Example: Use of bvector<> for k-mer fingerprint
    K should be short, no minimizers here (k < 24)

    Example loads FASTA file (large multi-molecule file is expected,
    builds a collection of k-mers for each molecule and runs
    clusterization algorithm on the input collection using
    set intersect (logical AND) as a similarity measure.

    Example uses std::async for running parallel jobs.
*/

#include <assert.h>
#include <stdlib.h>
#include <math.h>

#include <iostream>
#include <vector>
#include <list>
#include <map>
#include <algorithm>
#include <utility>
#include <memory>

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
#include "bmrandom.h"


// BitMagic utilities for debug and timings
#include "bmdbg.h"
#include "bmtimer.h"
#include "bmundef.h" /* clear the pre-proc defines from BM */

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
typedef bm::bvector<>                                 bvector_type;
typedef std::vector<char>                             vector_char_type;
typedef bm::dynamic_heap_matrix<unsigned, bm::bvector<>::allocator_type> distance_matrix_type;
typedef std::vector<std::unique_ptr<bvector_type> >   bvector_ptr_vector_type;
typedef bvector_type::size_type                       bv_size_type;
typedef std::vector<std::pair<bv_size_type, bv_size_type> > bv_ranges_vector;


// Global vars
//
bm::chrono_taker<>::duration_map_type     timing_map;


/// wait for any opening in a list of futures
///    used to schedule parallel tasks with CPU overbooking control
///
template<typename FV>
void wait_for_slot(FV& futures, unsigned* parallel_cnt, unsigned concurrency)
{
    do
    {
        for (auto e = futures.begin(); e != futures.end(); ++e)
        {
            std::future_status status = e->wait_for(std::chrono::milliseconds(100));
            if (status == std::future_status::ready)
            {
                (*parallel_cnt) -= 1;
                futures.erase(e);
                break;
            }
        } // for e
    } while (*parallel_cnt >= concurrency);
}

/// Collection of sequences and k-mer fingerprint vectors
///
class CSequenceColl
{
public:
    typedef std::vector<unsigned char> buffer_type;
public:

    CSequenceColl()
    {}
    CSequenceColl(const CSequenceColl&) = delete;

    void add_sequence(const string& acc, vector_char_type* seq_ptr)
    {
        m_acc.push_back(acc);
        m_seqs.emplace_back(seq_ptr);
    }

    void set_buffer(size_t i, const buffer_type& buf)
    {
        unique_ptr<buffer_type> buf_ptr(new buffer_type(buf));
        {
            static std::mutex   mtx_counts_lock;
            std::lock_guard<std::mutex> guard(mtx_counts_lock);

            if (m_kmer_bufs.size() <= i)
                m_kmer_bufs.resize(i+1);
            m_kmer_bufs[i].reset(buf_ptr.release());
        }
    }
    void sync_buffers_size()
    {
        m_kmer_bufs.resize(this->size());
    }

    size_t size() const
        { assert(m_seqs.size() == m_acc.size()); return m_seqs.size(); }

    const string& get_acc(size_t i) const { return m_acc[i]; }
    const vector_char_type& get_sequence(size_t i) const { return *(m_seqs[i]); }

    size_t seq_size(size_t i) const { return m_seqs[i]->size(); }

    size_t total_seq_size() const
    {
        size_t sum = 0;
        for (size_t i = 0; i < m_seqs.size(); ++i)
            sum += seq_size(i);
        return sum;
    }

    ///
    size_t buf_size() const { return m_kmer_bufs.size(); }

    /// Get k-mer vector BLOB size
    size_t get_buf_size(size_t i) const { return m_kmer_bufs[i]->size(); }

    /// Get k-mer BLOB pointer
    const unsigned char* get_buf(size_t i) const
    {
        const buffer_type* p =  m_kmer_bufs[i].get();
        if (!p)
            return 0;
        return p->data();
    }

    /// Deserialize group of k-mer fingerprint vectors
    ///
    void deserialize_k_mers(bvector_ptr_vector_type& k_mers_vect,
                            const bm::bvector<>& bv_req,
                            bm::bvector<>::size_type bv_req_count) const;



private:
    vector<unique_ptr<vector_char_type> > m_seqs;
    vector<string>                        m_acc;
    vector<unique_ptr<buffer_type> >      m_kmer_bufs;
};

void CSequenceColl::deserialize_k_mers(bvector_ptr_vector_type& k_mers_vect,
                                       const bm::bvector<>& bv_req,
                                bm::bvector<>::size_type bv_req_count) const
{
    std::list<std::future<void> > futures;

    k_mers_vect.resize(0);
    k_mers_vect.reserve(bv_req_count);

    bm::bvector<>::enumerator en(bv_req.first());
    for (; en.valid(); ++en)
    {
        auto i_idx = *en;
        bm::bvector<>* bv = new bm::bvector<>(); // k-mer fingerprint
        k_mers_vect.emplace_back(bv);
        const unsigned char* buf = this->get_buf(i_idx);
        if (!buf)
            continue;
        futures.emplace_back(                      // async decompress
            std::async(std::launch::async,
            [bv, buf]()
            { BM_DECLARE_TEMP_BLOCK(tb); bm::deserialize(*bv, buf, tb); }
            ));

    } // for en
    for (auto& e : futures)
        e.wait();
}


// -----------------------------------------------------------------------
//

/// Load multi-sequence FASTA
///
static
int load_FASTA(const std::string& fname, CSequenceColl& seq_coll)
{
    unique_ptr<vector_char_type> seq_vect(new vector_char_type());
    std::string line, acc;

    std::ifstream fin(fname.c_str(), std::ios::in);
    if (!fin.good())
        return -1;
    for (size_t i = 0; std::getline(fin, line); ++i)
    {
        if (line.empty())
            continue;

        if (line.front() == '>') // defline
        {
            if (!acc.empty())
            {
                seq_vect->shrink_to_fit();
                seq_coll.add_sequence(acc, seq_vect.release());
                acc.resize(0);
                seq_vect.reset(new vector_char_type());
            }

            std::size_t pos = line.find_first_of(":");
            if (pos == std::string::npos) // not found
            {
                acc = line;
            }
            else
            {
                acc = line.substr(1, pos-1);
            }
            continue;
        }
        for (std::string::iterator it = line.begin(); it != line.end(); ++it)
            seq_vect->push_back(*it);
    } // for

    if (!acc.empty())
    {
        seq_vect->shrink_to_fit();
        seq_coll.add_sequence(acc, seq_vect.release());
    }

    cout << "\r                            \r" << endl;
    return 0;
}

/// save k-mer vectors to a file
static
void save_kmer_buffers(const std::string& fname, const CSequenceColl& seq_coll)
{
    char magic_ch = '\t';
    std::ofstream bfile (fname, std::ios::out | std::ios::binary);
    if (!bfile.good())
    {
        std::cerr << "Cannot open file for write: " << fname << std::endl;
        exit(1);
    }

    // save collection size
    size_t sz = seq_coll.size();
    bfile.write((char*)&sz, std::streamsize(sizeof(sz)));

    // save the collection elements
    //
    for (size_t i = 0; i < sz; ++i)
    {
        size_t buf_size = 0;
        const unsigned char* buf = seq_coll.get_buf(i);
        if (!buf)
        {
            bfile.write((char*)&buf_size, std::streamsize(sizeof(buf_size)));
            continue;
        }
        buf_size = seq_coll.get_buf_size(i);
        bfile.write((char*)&buf_size, std::streamsize(sizeof(buf_size)));
        if (buf_size)
        {
            bfile.write((char*)buf, std::streamsize(buf_size));
            bfile.write((char*)&magic_ch, 1);
        }
    } // for i
}

/// Load k-mer vectors
///
static
void load_kmer_buffers(const std::string& fname, CSequenceColl& seq_coll)
{
    char magic_ch = '\t';
    std::ifstream bfile (fname, std::ios::in | std::ios::binary);
    if (!bfile.good())
    {
        std::cerr << "Cannot open file for read: " << fname << std::endl;
        exit(1);
    }

    // save collection size
    size_t sz;
    bfile.read((char*)&sz, std::streamsize(sizeof(sz)));

    CSequenceColl::buffer_type buf;

    // load the collection elements
    //
    for (size_t i = 0; i < sz; ++i)
    {
        size_t buf_size = 0;
        bfile.read((char*)&buf_size, std::streamsize(sizeof(buf_size)));
        if (buf_size)
        {
            buf.resize(buf_size);
            bfile.read((char*) buf.data(), std::streamsize(buf_size));
            char control_ch = 0;
            bfile.read((char*)&control_ch, 1);
            if (control_ch != magic_ch)
            {
                cerr << "Error: read failure!" << endl;
                exit(1);
            }
            seq_coll.set_buffer(i, buf);
        }

    } // for i
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
    static char lut[] = { 'A', 'T', 'G', 'C', 'N', 'M', '$' };
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
    @param k_buf    - sort buffer for generated k-mers
    @param chunk_size - sort buffer size (number of k-mers per sort)
 */
template<typename BV>
void generate_k_mer_bvector(BV& bv,
                            const vector_char_type& seq_vect,
                            unsigned k_size,
                            std::vector<bm::id64_t>& k_buf,
                            const bm::id64_t chunk_size = 400000000
                            )
{
    bv.clear();
    bv.init(); // need to explicitly init to use bvector<>::set_bit_no_check()
    if (seq_vect.empty())
        return;
    const char* dna_str = &seq_vect[0];

    k_buf.reserve(size_t(chunk_size));
    k_buf.resize(0);

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

                if (k_buf.size() == chunk_size) // soring check.point
                {
                    std::sort(k_buf.begin(), k_buf.end());
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
            std::sort(k_buf.begin(), k_buf.end());
            bv.set(&k_buf[0], k_buf.size(), bm::BM_SORTED); // fast bulk set
        }
    }
}

std::atomic_ullong                      k_mer_progress_count(0);

static
void generate_k_mers(CSequenceColl& seq_coll, unsigned k_size,
                     size_t from, size_t to)
{
    assert(from <= to);
    if (!seq_coll.size() || (from >= seq_coll.size()))
        return;

    std::vector<bm::id64_t> k_buf; // sort buffer
    BM_DECLARE_TEMP_BLOCK(tb)

    CSequenceColl::buffer_type buf;
    typedef bm::bvector<>::allocator_type        allocator_type;
    typedef allocator_type::allocator_pool_type  allocator_pool_type;
    allocator_pool_type  pool; // local pool for blocks

    bm::bvector<> bv;
    bm::bvector<>::mem_pool_guard mp_guard_bv; // memory pool reduces allocation calls to heap
    mp_guard_bv.assign_if_not_set(pool, bv);

    if (!to || to >= seq_coll.size())
        to = seq_coll.size()-1;

    bm::serializer<bm::bvector<> > bvs; // serializer object
    bvs.set_bookmarks(false);

    unsigned cnt = 0;
    for (size_t i = from; i <= to; ++i)
    {
        const vector_char_type& seq_vect = seq_coll.get_sequence(i);
        generate_k_mer_bvector(bv, seq_vect, k_size, k_buf);

        // serialize the vector
        //
        typename bm::bvector<>::statistics st;
        bv.optimize(tb, bm::bvector<>::opt_compress, &st);

        buf.resize(st.max_serialize_mem);

        size_t blob_size = bvs.serialize(bv, &buf[0], buf.size());
        buf.resize(blob_size);

        seq_coll.set_buffer(i, buf);

        // local progress report counter is just to avoid atomic too often
        ++cnt;
        if (cnt >= 100)
        {
            k_mer_progress_count.fetch_add(cnt);
            cnt = 0;
        }

    } // for i
}

static
void generate_k_mers_parallel(CSequenceColl& seq_coll, unsigned k_size,
                              unsigned concurrency)
{
    size_t total_seq_size = seq_coll.total_seq_size();

    if (!concurrency)
        concurrency = 1;

    size_t batch_size = total_seq_size / concurrency;
    if (!batch_size)
        batch_size = total_seq_size;
    std::list<std::future<void> > futures;

    for (size_t from = 0; from <= seq_coll.size(); )
    {
        size_t to = from;
        for (size_t to_pick = 0; to < seq_coll.size(); ++to)
        {
            to_pick += seq_coll.seq_size(to);
            if (to_pick >= batch_size)
                break;
        } // for

        futures.emplace_back(
            std::async(std::launch::async,
            [&seq_coll, k_size, from, to]() { generate_k_mers(seq_coll, k_size, from, to); }
            ));

        from = to+1;
    } // for from

    // wait for all jobs to finish, print progress report
    //
    unsigned long long cnt = seq_coll.size();
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

// -----------------------------------------------------------------------

/// Group (clustrer) of sequences
///
class CSeqGroup
{
public:
    CSeqGroup(bm::id64_t lead_id = ~0ull)
        : m_lead_id(lead_id), m_bv_members(bm::BM_GAP)
    {
        m_bv_members.set(lead_id);
    }

    /// set id for the group representative
    void set_lead(bm::id64_t lead_id)
        { add_member(m_lead_id = lead_id); }
    /// Get lead id
    bm::id64_t get_lead() const { return m_lead_id; }

    /// check is cluster is non-empty
    bool is_assigned() { return m_lead_id != ~0ull; }

    /// add a member to the group
    void add_member(bm::id64_t id) { m_bv_members.set_bit_no_check(id); }
    void add_member(bm::id64_t id, const bm::bvector<>& bv_kmer)
    {
        m_bv_members.set_bit_no_check(id);
        m_bv_kmer_union |= bv_kmer;
    }
    void add_member_sync(bm::id64_t id, const bm::bvector<>& bv_kmer)
    {
        std::lock_guard<std::mutex> guard(mtx_add_member_lock);
        m_bv_members.set_bit_no_check(id); // a bit faster than set()
        m_bv_kmer_union |= bv_kmer;
    }

    void merge_member_sync(bm::bvector<>& bv_seq, bm::bvector<>& bv_kmer)
    {
        std::lock_guard<std::mutex> guard(mtx_add_member_lock);
        m_bv_members.merge(bv_seq);
        m_bv_kmer_union.merge(bv_kmer);
    }

    bm::id64_t count_and_union_sync(const bm::bvector<>& bv)
    {
        std::lock_guard<std::mutex> guard(mtx_add_member_lock);
        return bm::count_and(bv, m_bv_kmer_union);
    }

    void clear_member(bm::id64_t id) { m_bv_members.set(id, false); }

    bm::bvector<>& get_rep() { return m_bv_rep; }
    const bm::bvector<>& get_rep() const { return m_bv_rep; }

    const bm::bvector<>& get_members() const { return m_bv_members; }
    bm::bvector<>& get_members() { return m_bv_members; }

    bm::bvector<>& get_kmer_union() { return m_bv_kmer_union; }
    const bm::bvector<>& get_kmer_union() const { return m_bv_kmer_union; }

    friend class CSeqClusters;

private:
    bm::id64_t    m_lead_id;       ///< groups lead vector ID
    bm::bvector<> m_bv_members;    ///< vector of all group member ids
    std::mutex    mtx_add_member_lock; ///< concurrency protection for add_member_sync()

    bm::bvector<> m_bv_rep;        ///< group's representative vector
    bm::bvector<> m_bv_kmer_union; ///< union of all k-mers of members
};

// -----------------------------------------------------------------------
//
//

class CSeqClusters
{
public:
    typedef std::vector<std::unique_ptr<CSeqGroup> > groups_vector_type;
public:
    CSeqClusters()
    {}
    CSeqClusters(const CSeqClusters&)=delete;

    void add_group(CSeqGroup* sg) { m_seq_groups.emplace_back(sg); }

    /// memebers moved into their own group
    void take_group(bm::bvector<>& bv_members);

    /// Acquire all groups from another cluster collection
    ///
    void merge_from(CSeqClusters& sc);

    /// Remove groups which turned empty after clusterization
    void clear_empty_groups();

    /// Compute union of all cluster group members
    const bm::bvector<>& union_all_groups();

    /// Resolve duplicate membership between groups
    ///
    void resolve_duplicates(const CSequenceColl& seq_coll);

    /// Find the best representatives in all cluster groups
    /// the criteria is maximum absolute similarity to all members
    ///
    void elect_leaders(const CSequenceColl& seq_coll,
                       unsigned concurrency);

    /// calculate avg cluster population count
    bm::id64_t compute_avg_count() const;


    size_t groups_size() const { return m_seq_groups.size(); }
    CSeqGroup* get_group(size_t idx) { return m_seq_groups[idx].get(); }

    /// print clusterization report
    void print_summary(const char* title) const;

private:
    bm::bvector<>         m_all_members; ///< Union of all group members
    groups_vector_type    m_seq_groups; ///< vector of all formed clusters
    bm::aggregator<bvector_type> agg; ///< fast aggregator for set UNION ops
};

void CSeqClusters::clear_empty_groups()
{
    for (groups_vector_type::iterator it = m_seq_groups.begin();
         it != m_seq_groups.end(); )
    {
        CSeqGroup* sg = it->get();
        bm::bvector<>& bv_mem = sg->get_members();
        auto cnt = bv_mem.count();
        if (cnt < 2) // practically empty group
            it = m_seq_groups.erase(it);
        else
            ++it;
    } // for
}

void CSeqClusters::take_group(bm::bvector<>& bv_members)
{
    bm::id64_t lead_id = bv_members.get_first();
    CSeqGroup* sg = new CSeqGroup(lead_id);
    sg->get_members().swap(bv_members); // move members to the cluster
    add_group(sg);
}

const bm::bvector<>& CSeqClusters::union_all_groups()
{
    agg.set_optimization();
    m_all_members.clear();
    for (groups_vector_type::const_iterator it = m_seq_groups.begin();
         it != m_seq_groups.end(); ++it)
    {
        const CSeqGroup* sg = it->get();
        agg.add(&sg->get_members());
    } // for

    agg.combine_or(m_all_members); // run UNION of all member vectors
    agg.reset();

    return m_all_members;
}

bm::id64_t CSeqClusters::compute_avg_count() const
{
    bm::id64_t cnt = 0;
    for (groups_vector_type::const_iterator it = m_seq_groups.begin();
         it != m_seq_groups.end(); ++it)
    {
        const CSeqGroup* sg = it->get();
        cnt += sg->get_members().count();
    }
    cnt = cnt / m_seq_groups.size();
    return cnt;
}


void CSeqClusters::merge_from(CSeqClusters& sc)
{
    for (auto it = sc.m_seq_groups.begin(); it != sc.m_seq_groups.end(); ++it)
        m_seq_groups.emplace_back(it->release());
    sc.m_seq_groups.clear();
}


void resolve_duplicates(CSeqGroup& seq_group1,
                        CSeqGroup& seq_group2,
                        const CSequenceColl& seq_coll);

void CSeqClusters::resolve_duplicates(const CSequenceColl& seq_coll)
{
    for (size_t i = 0; i < m_seq_groups.size(); ++i)
    {
        CSeqGroup* sg1 = m_seq_groups[i].get();
        for (size_t j = 0; j < i; ++j)
        {
            CSeqGroup* sg2 = m_seq_groups[j].get();
            // resolve pairwise conflicts between cluster groups
            ::resolve_duplicates(*sg1, *sg2, seq_coll);
        } // for j
    } // for i
}

/// Compute similarity distances for one row/vector (1:N) of distance matrix
///
static
void compute_and_sim_row(
        unsigned* row,
        const bm::bvector<>* bv_i,
        size_t i,
        const std::vector<std::unique_ptr<bvector_type> >& k_mers_vect)
{
    size_t j;
    for (j = 0; j < i; ++j)
    {
        const bm::bvector<>* bv_j = k_mers_vect[j].get();
        bm::id64_t and_cnt = bm::count_and(*bv_i, *bv_j);
        row[j] = unsigned(and_cnt);
    }
    auto cnt = bv_i->count();
    row[j] = unsigned(cnt);
}

/// Compute similarity distances matrix (COUNT(AND(a, b))
///
static
void compute_and_sim(distance_matrix_type& dm,
                     const CSequenceColl& seq_coll,
                     const bm::bvector<>& bv_mem,
                     bm::bvector<>::size_type bv_mem_cnt,
                     unsigned concurrency)
{
    if (concurrency < 1)
        concurrency = 1;
    auto N = bv_mem_cnt;

    const unsigned k_max_electors = 500;
    bm::bvector<> bv_sub; // subset vector
    if (N > k_max_electors) // TODO: parameterize the
    {
        bm::random_subset<bm::bvector<> > rsub;
        rsub.sample(bv_sub, bv_mem, k_max_electors); // pick random sunset
    }


    bvector_ptr_vector_type k_mers_vect;

    // materialize list of vectors used in distance calculation
    //
    seq_coll.deserialize_k_mers(k_mers_vect, bv_mem, N);

    std::list<std::future<void> > futures;
    unsigned parallel_cnt = 0; // number of jobs in flight at a time

    size_t i;
    for (i = 0; i < N; ++i)
    {
        const bm::bvector<>* bv_i = k_mers_vect[i].get();
        unsigned* row = dm.row(i);
        do
        {
            if (parallel_cnt < concurrency)
            {
                futures.emplace_back(
                    std::async(std::launch::async,
                    [row, bv_i, i, &k_mers_vect]()
                        { compute_and_sim_row(row, bv_i, i, k_mers_vect); }
                    ));
                ++parallel_cnt;
                break;
            }
            else
            {
                // wait for an async() slot to open (overbooking control)
                ::wait_for_slot(futures, &parallel_cnt, concurrency);
            }
        } while(1);
    } // for i

    // sync point
    for (auto& e : futures)
        e.wait();

    dm.replicate_triange(); // copy to full simetrical distances from triangle
}

/// Compute union (Universe) of all k-mers in the cluster group
/// Implemented as a OR of all k-mer fingerprints
static
void compute_seq_group_union(CSeqGroup&           seq_group,
                             const CSequenceColl& seq_coll)
{
    BM_DECLARE_TEMP_BLOCK(tb)

    bm::bvector<>& bv_kmer_union = seq_group.get_kmer_union();
    bv_kmer_union.clear();

    bm::bvector<>& bv_all_members = seq_group.get_members();
    for(bm::bvector<>::enumerator en(bv_all_members) ;en.valid(); ++en)
    {
        auto idx = *en;
        const unsigned char* buf = seq_coll.get_buf(idx);
        if (!buf)
            continue;

        bm::deserialize(bv_kmer_union, buf, tb); // deserialize is an OR operation
    } // for en
    bv_kmer_union.optimize(tb);
}

void CSeqClusters::elect_leaders(const CSequenceColl& seq_coll,
                                 unsigned concurrency)
{
    bm::operation_deserializer<bm::bvector<> > od;
    bm::random_subset<bm::bvector<> > rsub;

    std::list<std::future<void> > futures;

    for (size_t k = 0; k < m_seq_groups.size(); ++k)
    {
        CSeqGroup* sg = m_seq_groups[k].get();
        bm::bvector<>& bv_all_members = sg->get_members();
        bv_all_members.set(sg->get_lead());
        auto N = bv_all_members.count();
        auto all_members_count = N; (void) all_members_count;

        // determine the size of the "electoral colledge" on available concurrency
        unsigned k_max_electors = 200 * unsigned(log2(concurrency));
        if (k_max_electors < 500)
            k_max_electors = 500;

        const bm::bvector<>* bv_mem = &bv_all_members; // vector of electors
        bm::bvector<> bv_sub; // subset vector
        if (N > k_max_electors) // TODO: parameterize the
        {
            // pick a random sub-set as the "electoral colledge"
            rsub.sample(bv_sub, bv_all_members, k_max_electors);
            bv_sub.set(sg->get_lead()); // current leader always takes part
            bv_mem = &bv_sub;
            N = bv_sub.count();
        }

        // NxN distance matrix between members
        //
        distance_matrix_type dm(N, N);
        dm.init();
        dm.set_zero();

        // compute triangular distance matrix
        //
        compute_and_sim(dm, seq_coll, *bv_mem, N, concurrency);

        // leader election based on maximum similarity to other cluster
        // elements
        //
        bm::id64_t best_score = 0;
        bm::id64_t leader_idx = sg->get_lead();
        assert(bv_mem->test(leader_idx));
        auto rank = bv_mem->count_range(0, leader_idx);
        assert(rank);
        dm.sum(best_score, rank-1);  // use sum or all similarities as a score for the leader
        bm::id64_t old_leader_idx = leader_idx;

        bm::bvector<>::enumerator en(bv_mem->first());
        for (size_t i = 0; en.valid(); ++en, ++i)
        {
            bm::id64_t cand_score;
            dm.sum(cand_score, i);  // use sum or all similarities as a score
            if (cand_score > best_score) // better candidate for a leader
            {
                best_score = cand_score;
                leader_idx = *en;
            }
        } // for en

        if (leader_idx != old_leader_idx) // found a new leader
        {
            sg->set_lead(leader_idx);
            const unsigned char* buf = seq_coll.get_buf(leader_idx);
            assert(buf);
            bm::bvector<>* bv = &sg->get_rep();

            // async replace the leader representative k-mer vector
            //
            futures.emplace_back(
                std::async(std::launch::async,
                [bv, buf]() { bv->clear(); bm::deserialize(*bv, buf); }
                ));
            // re-compute k-mer UNION of all vectors
            futures.emplace_back(
                std::async(std::launch::async,
                [sg, &seq_coll]() { compute_seq_group_union(*sg, seq_coll); }
                ));
        }

    } // for i

    for (auto& e : futures) // wait for all forked tasks
        e.wait();

}


void CSeqClusters::print_summary(const char* title) const
{
    cout << title << endl;
    for (size_t i = 0; i < m_seq_groups.size(); ++i)
    {
        const CSeqGroup* sg = m_seq_groups[i].get();
        const bm::bvector<>& bv_mem = sg->get_members();
        cout << sg->get_lead() << ": "
             << bv_mem.count() << endl;
    }
    cout << "-----------\nTotal: " <<  m_seq_groups.size() << endl << endl;
}

///
static
void compute_group(CSeqGroup& seq_group,
                   const CSequenceColl& seq_coll,
                   const bm::bvector<>& bv_exceptions,
                   float similarity_cut_off)
{
    assert(similarity_cut_off < 1);
    assert(seq_group.is_assigned());

    auto sz = seq_coll.buf_size();
    size_t lead_id = seq_group.get_lead();
    if (lead_id >= sz)
        return;

    const unsigned char* buf = seq_coll.get_buf(lead_id);
    if (!buf)
        return;

    bm::bvector<>& bv = seq_group.get_rep();
    bm::deserialize(bv, buf);

    auto i_cnt = bv.count();
    // approximate number of k-mers we consider similar
    float similarity_target = float(i_cnt * float(similarity_cut_off));


    bm::operation_deserializer<bm::bvector<> > od;

    bool found = false;
    for (size_t i = 0; i < sz; ++i)
    {
        bool is_except = bv_exceptions.test(i);
        if (is_except)
            continue;
        buf = seq_coll.get_buf(i);
        if (!buf)
            continue;
        // constant deserializer AND just to count the product
        // without actual deserialization (from the compressed BLOB)
        //
        bm::id64_t and_cnt = od.deserialize(bv, buf, 0, bm::set_COUNT_AND);

        if (and_cnt && (float(and_cnt) > similarity_target)) // similar enough to be a candidate
        {
            seq_group.add_member(i);
            found = true;
        }
    } // for i

    if (!found)
        seq_group.clear_member(lead_id);
}

/// Resolve duplicate members between two groups
void resolve_duplicates(CSeqGroup& seq_group1,
                        CSeqGroup& seq_group2,
                        const CSequenceColl& seq_coll)
{
    if (&seq_group1 == & seq_group2) // self check
        return;

    bm::bvector<>& bv1 = seq_group1.get_members();
    bm::bvector<>& bv2 = seq_group2.get_members();

    // build intersect between group members
    bm::bvector<> bv_and;
    bv_and.bit_and(bv1, bv2, bm::bvector<>::opt_none);

    if (bv_and.any()) // double membership detected
    {
        bm::bvector<>& bv_rep1 = seq_group1.get_rep();
        bm::bvector<>& bv_rep2 = seq_group2.get_rep();

        bm::operation_deserializer<bm::bvector<> > od;

        // evaluate each double member for best membership placement
        //
        for (bm::bvector<>::enumerator en(bv_and); en.valid(); ++en)
        {
            auto idx = *en;
            auto lead_idx1 = seq_group1.get_lead();
            auto lead_idx2 = seq_group2.get_lead();
            assert(lead_idx1 != lead_idx2);

            if (idx == lead_idx1)
            {
                seq_group2.clear_member(idx);
                continue;
            }
            if (idx == lead_idx2)
            {
                seq_group1.clear_member(idx);
                continue;
            }

            const unsigned char* buf = seq_coll.get_buf(idx);
            assert(buf);
            // resolve conflict based on max.absolute similarity
            //
            bm::id64_t and_cnt1 = od.deserialize(bv_rep1, buf, 0, bm::set_COUNT_AND);
            bm::id64_t and_cnt2 = od.deserialize(bv_rep2, buf, 0, bm::set_COUNT_AND);
            if (and_cnt1 >= and_cnt2)
                seq_group2.clear_member(idx);
            else
                seq_group1.clear_member(idx);

        } // for
    }
}

/// Utility class to accumulate cahnges to cluster
/// before commiting it (mutex syncronous operation)
///
struct CKMerAcc
{
    CKMerAcc(size_t sz)
    : bv_members(sz), bv_kmers(sz)
    {}

    void add(size_t cluster_id,
             bm::bvector<>::size_type m_id, bm::bvector<>& bv_kmer)
    {
        bm::bvector<>* bv_m = bv_members[cluster_id].get();
        bm::bvector<>* bv_k = bv_kmers[cluster_id].get();
        if (!bv_m)
        {
            assert(!bv_kmers[cluster_id].get());
            bv_members[cluster_id].reset(new bm::bvector<>(bm::BM_GAP));
            bv_kmers[cluster_id].reset(new bm::bvector<>(bm::BM_GAP));
            bv_m = bv_members[cluster_id].get();
            bv_k = bv_kmers[cluster_id].get();
        }
        bv_m->set(m_id);
        bv_k->merge(bv_kmer);
    }

    bvector_ptr_vector_type bv_members;
    bvector_ptr_vector_type bv_kmers;
};

/// Compute AND similarity to all available clusters assign to the most similar
/// using cluster representative
///
static
void assign_to_best_cluster(CSeqClusters& cluster_groups,
                            const CSequenceColl& seq_coll,
                            const bm::bvector<>& bv_seq_ids,
                            bm::bvector<>::size_type seq_id_from,
                            bm::bvector<>::size_type seq_id_to)
{
    BM_DECLARE_TEMP_BLOCK(tb)
    bm::bvector<>::allocator_pool_type pool;

    CKMerAcc acc(cluster_groups.groups_size());

    bm::bvector<>::enumerator en(bv_seq_ids);
    en.go_to(seq_id_from);

    for ( ;en.valid(); ++en)
    {
        auto seq_id = *en;
        if (seq_id > seq_id_to)
            break;
        const unsigned char* buf = seq_coll.get_buf(seq_id);
        if (!buf)
            continue;

        bm::bvector<> bv_k_mer;
        bv_k_mer.set_allocator_pool(&pool); // for faster memory recycling

        bm::deserialize(bv_k_mer, buf, tb);

        bm::bvector<>::size_type best_score(0);
        bm::bvector<>::size_type cluster_idx(~0ull);

        // analyse candidate's similarity to all clusters via representative
        //
        for (size_t i = 0; i < cluster_groups.groups_size(); ++i)
        {
            CSeqGroup* sg = cluster_groups.get_group(i);
            //  - COUNT(AND) similarity to the representative of each cluster
            //
            bm::bvector<>& bv_rep_k_mer = sg->get_rep();
            auto rep_and_cnt = bm::count_and(bv_k_mer, bv_rep_k_mer);
            if (rep_and_cnt > best_score)
            {
                cluster_idx = i; best_score = rep_and_cnt;
            }
        } // for i

        if (cluster_idx != ~0ull) // cluster association via representative
        {
            acc.add(cluster_idx, seq_id, bv_k_mer);
        }
    } // for all seq-ids in the range

    // merge all accumulated cluster assignmnets all at once
    //
    for (size_t i = 0; i < cluster_groups.groups_size(); ++i)
    {
        bm::bvector<>* bv_m = acc.bv_members[i].get();
        if (!bv_m)
            continue;
        bm::bvector<>* bv_k = acc.bv_kmers[i].get();
        CSeqGroup* sg = cluster_groups.get_group(i);
        sg->merge_member_sync(*bv_m, *bv_k);
    } // for

}

/// Compute AND similarity to all available clusters assign to the most similar
/// using UNION of k-mers in the cluster
/// This is a more relaxed assignmnet, used when representative does not work
static
void assign_to_best_cluster_union(CSeqClusters& cluster_groups,
                            const CSequenceColl& seq_coll,
                            const bm::bvector<>& bv_seq_ids,
                            bm::bvector<>::size_type seq_id_from,
                            bm::bvector<>::size_type seq_id_to)
{
    BM_DECLARE_TEMP_BLOCK(tb)
    bm::bvector<>::allocator_pool_type pool;

    bm::bvector<>::enumerator en(bv_seq_ids);
    en.go_to(seq_id_from);

    for ( ;en.valid(); ++en)
    {
        auto seq_id = *en;
        if (seq_id > seq_id_to)
            break;
        const unsigned char* buf = seq_coll.get_buf(seq_id);
        if (!buf)
            continue;

        bm::bvector<> bv_k_mer;
        bv_k_mer.set_allocator_pool(&pool); // for faster memory recycling

        bm::deserialize(bv_k_mer, buf, tb);

        bm::bvector<>::size_type best_score(0);
        bm::bvector<>::size_type cluster_idx(~0ull);

        {
            // try to extend search using UNION of all cluster k-mers
            //
            best_score = 0;
            for (size_t i = 0; i < cluster_groups.groups_size(); ++i)
            {
                CSeqGroup* sg = cluster_groups.get_group(i);
                //  - COUNT(AND) similarity to the representative of each cluster
                //
                bm::id64_t uni_and_cnt = sg->count_and_union_sync(bv_k_mer);
                if (uni_and_cnt > best_score)
                {
                    cluster_idx = i; best_score = uni_and_cnt;
                }
            } // for i
            if (cluster_idx != ~0ull) // cluster association via representative
            {
                CSeqGroup* sg = cluster_groups.get_group(cluster_idx);
                sg->add_member_sync(seq_id, bv_k_mer);
            }
        }
    } // for all seq-ids in the range

}


/// Pick random sequences as cluster seed elements, try attach
/// initial sequences based on weighted similarity measure
///
static
void compute_random_clusters(CSeqClusters& cluster_groups,
                              const CSequenceColl& seq_coll,
                              const bm::bvector<>& bv_total,
                              bm::random_subset<bvector_type>& rsub,
                              unsigned num_clust,
                              float similarity_cut_off,
                              unsigned concurrency)
{
    bm::bvector<> bv_rsub; // random subset of sequences
    rsub.sample(bv_rsub, bv_total, num_clust); // pick random sequences as seeds

    std::list<std::future<void> > futures;
    unsigned parallel_cnt = 0;

    for (bm::bvector<>::enumerator en=bv_rsub.first(); en.valid(); ++en)
    {
        auto idx = *en;
        CSeqGroup* sg = new CSeqGroup(idx);
        cluster_groups.add_group(sg);

        do
        {
            if (parallel_cnt < concurrency)
            {
                futures.emplace_back(
                    std::async(std::launch::async,
                    [&seq_coll, sg, &bv_rsub, similarity_cut_off]()
                        { compute_group(*sg, seq_coll, bv_rsub, similarity_cut_off); }
                    ));
                ++parallel_cnt;
                break;
            }
            else
            {
                // wait for an async() slot to open
                ::wait_for_slot(futures, &parallel_cnt, concurrency);
            }
        } while(1);
    } // for en

    // wait for completion of initial cluster group formation
    for (auto& e : futures)
        e.wait();
}

static
void compute_jaccard_clusters(CSeqClusters& seq_clusters,
                              const CSequenceColl& seq_coll,
                              unsigned num_clust,
                              float similarity_cut_off,
                              unsigned concurrency)
{
    assert(similarity_cut_off < 1);
    if (!seq_coll.buf_size())
        return; // nothing to do

    bm::random_subset<bm::bvector<> > rsub; // sub-set getter

    bm::bvector<> bv_total;
    bv_total.set_range(0, seq_coll.buf_size());

    bm::id64_t rcount = 0;

    const unsigned max_pass = 3;
    for (unsigned pass = 0; pass < max_pass; ++pass)
    {
        CSeqClusters cluster_groups;

        compute_random_clusters(cluster_groups, seq_coll, bv_total, rsub,
                                num_clust, similarity_cut_off, concurrency);

        // remove possible empty clusters (inital seeds were picked at random)
        //
        cluster_groups.clear_empty_groups();

        cluster_groups.resolve_duplicates(seq_coll);

        cluster_groups.clear_empty_groups();

        // print summary after the initial formation of cluster groups
        //
        cluster_groups.print_summary("Inital cluster formations:");

        // re-elect representatives
        cluster_groups.elect_leaders(seq_coll, concurrency);

        cluster_groups.print_summary("After lead re-election:");

        // sub-set of sequence ids already distributed into clusters
        bm::bvector<>::size_type total_count = bv_total.count();
        {
            cout << " total = " << total_count << endl;
            const bm::bvector<>& bv_clust = cluster_groups.union_all_groups();
            cout << " clustered = " << bv_clust.count() << endl;

            bv_total -= bv_clust; // exlude all already clustered
            total_count = bv_total.count();
            cout << " remain = " << total_count << endl;
        }

        if (!total_count)
            break;

        std::list<std::future<void> > futures;

        // run split algorithm to determine approximately equal ranges
        // for parallel processing
        bv_ranges_vector pair_vect;
        bm::bvector<>::size_type split_rank = total_count / concurrency; // target population count per job
        bm::rank_range_split(bv_total, split_rank, pair_vect);
        assert(pair_vect.size());
        for (size_t k = 0; k < pair_vect.size(); ++k)
        {
            auto seq_id_from = pair_vect[k].first;
            auto seq_id_to = pair_vect[k].second;
            futures.emplace_back(
                std::async(std::launch::async,
                [&cluster_groups, &seq_coll, &bv_total, seq_id_from, seq_id_to]()
                    { assign_to_best_cluster(cluster_groups, seq_coll, bv_total, seq_id_from, seq_id_to); }
                ));
        }
        for (auto& e : futures) // sync point
            e.wait();

        cluster_groups.print_summary("Clusters after phase 2 recruitment");

        // check if there are sequences not yet belonging to any cluster
        {
            const bm::bvector<>& bv_clust = cluster_groups.union_all_groups();
            bv_total -= bv_clust; // exlude all already clustered
            rcount = bv_total.count();
            if (rcount)
            {
                cout << "Undistributed sequences = " << rcount << endl;
            }
            else
            {
                seq_clusters.merge_from(cluster_groups);
                break;
            }
        }
        seq_clusters.merge_from(cluster_groups);
        bm::id64_t avg_group_count = seq_clusters.compute_avg_count();
        if (rcount < avg_group_count) // not worth another pass ?
        {
            seq_clusters.take_group(bv_total);
            break;
        }

        cout << "PASS=" << (pass+1) << endl << endl;
    } // for pass

    // try to assign to the global pool of clusters using UNION
    // which relaxes assignmnet

    if (rcount)
    {
        {
            const bm::bvector<>& bv_clust = seq_clusters.union_all_groups();
            cout << endl << " clustered = " << bv_clust.count() << endl;

            bv_total -= bv_clust; // exlude all already clustered
            rcount = bv_total.count();
            cout << " remain = " << rcount << endl;
        }

        if (rcount)
        {
            bv_ranges_vector pair_vect;
            bm::bvector<>::size_type split_rank = rcount / concurrency;
            bm::rank_range_split(bv_total, split_rank, pair_vect);

            std::list<std::future<void> > futures;

            for (size_t k = 0; k < pair_vect.size(); ++k)
            {
                auto seq_id_from = pair_vect[k].first;
                auto seq_id_to = pair_vect[k].second;
                futures.emplace_back(
                    std::async(std::launch::async,
                    [&seq_clusters, &seq_coll, &bv_total, seq_id_from, seq_id_to]()
                        { assign_to_best_cluster_union(seq_clusters, seq_coll, bv_total, seq_id_from, seq_id_to); }
                    ));
            }
            for (auto& e : futures) // sync point
                e.wait();
            {
                const bm::bvector<>& bv_clust = seq_clusters.union_all_groups();
                cout << endl << " clustered = " << bv_clust.count() << endl;

                bv_total -= bv_clust; // exlude all already clustered
                rcount = bv_total.count();
                cout << " remain = " << rcount << endl;
            }
        }
    }

    if (rcount)
    {
        seq_clusters.take_group(bv_total);
    }

    seq_clusters.clear_empty_groups();
    seq_clusters.resolve_duplicates(seq_coll);
    seq_clusters.clear_empty_groups();

    seq_clusters.print_summary("Final clusters summary:");
}


int main(int argc, char *argv[])
{
    CSequenceColl          seq_coll;

    try
    {
        auto ret = parse_args(argc, argv);
        if (ret != 0)
        {
            cerr << "cmd-line parse error. " << endl;
            return ret;
        }

        if (!ifa_name.empty()) // FASTA file load
        {
            bm::chrono_taker tt1(cout, "1. Load FASTA", 1, &timing_map);

            // limitation: loads a single molecule only
            //
            auto res = load_FASTA(ifa_name, seq_coll);
            if (res != 0)
                return res;
        }

        cout << "Sequences size = " << seq_coll.size() << endl;

        if (ik_size && !ifa_name.empty())
        {
            {
                bm::chrono_taker tt1(cout, "2. Generate k-mers", 1, &timing_map);
                seq_coll.sync_buffers_size();
                generate_k_mers_parallel(seq_coll, ik_size, parallel_jobs);
            }

            if (!ikd_name.empty())
            {
                bm::chrono_taker tt1(cout, "3. Save k-mers", 1, &timing_map);
                save_kmer_buffers(ikd_name, seq_coll);
            }
        }

        if (ik_size && ifa_name.empty() && !ikd_name.empty())
        {
            {
            bm::chrono_taker tt1(cout, "4. Load k-mers", 1, &timing_map);
            load_kmer_buffers(ikd_name, seq_coll);
            }

            if (seq_coll.buf_size())
            {
                CSeqClusters seq_clusters;
                bm::chrono_taker tt1(cout, "5. k-mer similarity clustering", 1, &timing_map);
                compute_jaccard_clusters(seq_clusters, seq_coll, 10, 0.05f, parallel_jobs);
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

