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

/** \example xsample04a.cpp

*/

/*! \file xsample04a.cpp
    \brief Example: DNA index construction

*/

#include <iostream>
#include <sstream>
#include <regex>
#include <time.h>
#include <stdio.h>

#include <stdexcept>
#include <memory>
#include <vector>

#include <future>
#include <thread>
#include <mutex>

#include "bm.h"

#include "bmdbg.h"
#include "bmtimer.h"
#include "bmundef.h" /* clear the pre-proc defines from BM */

using namespace std;

static
void show_help()
{
    std::cerr
        << "BitMagic DNA Index Build Sample (c) 2018" << std::endl
        << "-fa   file-name     -- input FASTA file" << std::endl
        << "-j    number        -- number of parallel jobs to run" << std::endl
        << "-timing             -- collect timings"  << std::endl
      ;
}




// Arguments
//
std::string  ifa_name;
bool         is_timing = false;
unsigned     parallel_jobs = 4;

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
        if (arg == "-fa" || arg == "--fa")
        {
            if (i + 1 < argc)
            {
                ifa_name = argv[++i];
            }
            else
            {
                std::cerr << "Error: -fa requires file name" << std::endl;
                return 1;
            }
            continue;
        }
        if (arg == "-j" || arg == "--j")
        {
            if (i + 1 < argc)
            {
                parallel_jobs = unsigned(::atoi(argv[++i]));
            }
            else
            {
                std::cerr << "Error: -j requires number of jobs" << std::endl;
                return 1;
            }
            continue;
        }

        if (arg == "-timing" || arg == "--timing" || arg == "-t" || arg == "--t")
            is_timing = true;

    } // for i
    return 0;
}



// ----------------------------------------------------------------------------

bm::chrono_taker<>::duration_map_type  timing_map;

// FASTA format parser
static
int load_FASTA(const std::string& fname, std::vector<char>& seq_vect)
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



/**
    Utility for keeping all DNA finger print vectors and search
    using various techniques
*/
class DNA_FingerprintScanner
{
public:
    enum { eA = 0, eC, eG, eT, eN, eEnd };

    DNA_FingerprintScanner() {}

    /// Build fingerprint bit-vectors from the original sequence
    ///
    void Build(const vector<char>& sequence)
    {
        bm::bvector<>::insert_iterator iA = m_FPrintBV[eA].inserter();
        bm::bvector<>::insert_iterator iC = m_FPrintBV[eC].inserter();
        bm::bvector<>::insert_iterator iG = m_FPrintBV[eG].inserter();
        bm::bvector<>::insert_iterator iT = m_FPrintBV[eT].inserter();
        bm::bvector<>::insert_iterator iN = m_FPrintBV[eN].inserter();

        for (size_t i = 0; i < sequence.size(); ++i)
        {
            unsigned pos = unsigned(i);
            switch (sequence[i])
            {
            case 'A':
                iA = pos;
                break;
            case 'C':
                iC = pos;
                break;
            case 'G':
                iG = pos;
                break;
            case 'T':
                iT = pos;
                break;
            case 'N':
                iN = pos;
                break;
            default:
                break;
            }
        }
    }

    /// Build index using bulk insert iterator
    ///
    void BuildBulk(const vector<char>& sequence)
    {
        bm::bvector<>::bulk_insert_iterator iA(m_FPrintBV[eA], bm::BM_SORTED);
        bm::bvector<>::bulk_insert_iterator iC(m_FPrintBV[eC], bm::BM_SORTED);
        bm::bvector<>::bulk_insert_iterator iG(m_FPrintBV[eG], bm::BM_SORTED);
        bm::bvector<>::bulk_insert_iterator iT(m_FPrintBV[eT], bm::BM_SORTED);
        bm::bvector<>::bulk_insert_iterator iN(m_FPrintBV[eN], bm::BM_SORTED);

        for (size_t i = 0; i < sequence.size(); ++i)
        {
            unsigned pos = unsigned(i);
            switch (sequence[i])
            {
            case 'A':
                iA = pos;
                break;
            case 'C':
                iC = pos;
                break;
            case 'G':
                iG = pos;
                break;
            case 'T':
                iT = pos;
                break;
            case 'N':
                iN = pos;
                break;
            default:
                break;
            }
        }
    }


    /// Build fingerprint bit-vectors using bulk insert iterator and parallel
    /// processing
    ///
    void BuildParallel(const vector<char>& sequence, unsigned threads)
    {
        struct Func
        {
            DNA_FingerprintScanner*      target_idx;
            const std::vector<char>*     src_sequence;

            Func(DNA_FingerprintScanner* idx, const vector<char>& src) 
                : target_idx(idx), src_sequence(&src) {}

            void operator() (size_t from, size_t to)
            {
                const vector<char>& sequence = *src_sequence;
                bm::bvector<> bvA, bvT, bvG, bvC, bvN;
                
                {
                    bm::bvector<>::bulk_insert_iterator iA(bvA, bm::BM_SORTED);
                    bm::bvector<>::bulk_insert_iterator iC(bvC, bm::BM_SORTED);
                    bm::bvector<>::bulk_insert_iterator iG(bvG, bm::BM_SORTED);
                    bm::bvector<>::bulk_insert_iterator iT(bvT, bm::BM_SORTED);
                    bm::bvector<>::bulk_insert_iterator iN(bvN, bm::BM_SORTED);
                    for (size_t i = from; i < sequence.size() && (i < to); ++i)
                    {
                        unsigned pos = unsigned(i);
                        switch (sequence[i])
                        {
                        case 'A':
                            iA = pos;
                            break;
                        case 'C':
                            iC = pos;
                            break;
                        case 'G':
                            iG = pos;
                            break;
                        case 'T':
                            iT = pos;
                            break;
                        case 'N':
                            iN = pos;
                            break;
                        default:
                            break;
                        }
                    } // for
                    // Bulk insert iterator keeps an buffer, which has to be
                    // flushed, before all bits appear in the target vector
                    //
                    iA.flush();
                    iC.flush();
                    iT.flush();
                    iG.flush();
                    iN.flush();
                }
                // merge results of parallel processing back to index
                target_idx->MergeVector('A', bvA);
                target_idx->MergeVector('T', bvT);
                target_idx->MergeVector('G', bvG);
                target_idx->MergeVector('C', bvC);
                target_idx->MergeVector('N', bvN);
            }
        };
        
        if (threads <= 1)
        {
            BuildBulk(sequence);
            return;
        }

        // Create parallel async tasks running on a range of source sequence
        //
        std::vector<std::future<void> > futures;
        futures.reserve(8);
        unsigned range = unsigned(sequence.size() / threads);

        for (unsigned k = 0; k < sequence.size(); k += range)
        {
            futures.emplace_back(std::async(std::launch::async,
                                 Func(this, sequence),  k, k + range));
        }

        // wait for all tasks
        for (auto& e : futures)
        {
            e.wait();
        }
    }

    /// Thread sync bit-vector merge
    ///
    void MergeVector(char letter, bm::bvector<>& bv)
    {
        static std::mutex mtx_A;
        static std::mutex mtx_T;
        static std::mutex mtx_G;
        static std::mutex mtx_C;
        static std::mutex mtx_N;

        switch (letter)
        {
        case 'A':
        {
            std::lock_guard<std::mutex> guard(mtx_A);
            m_FPrintBV[eA].merge(bv);
        }
        break;
        case 'C':
        {
            std::lock_guard<std::mutex> guard(mtx_C);
            m_FPrintBV[eC].merge(bv);
        }
        break;
        case 'G':
        {
            std::lock_guard<std::mutex> guard(mtx_G);
            m_FPrintBV[eG].merge(bv);
        }
        break;
        case 'T':
        {
            std::lock_guard<std::mutex> guard(mtx_T);
            m_FPrintBV[eT].merge(bv);
        }
        break;
        case 'N':
        {
            std::lock_guard<std::mutex> guard(mtx_N);
            m_FPrintBV[eN].merge(bv);
        }
        break;
        default:
            break;
        }

    }

    /// Return fingerprint bit-vector
    const bm::bvector<>& GetVector(char letter) const
    {
        switch (letter)
        {
        case 'A':
            return m_FPrintBV[eA];
        case 'C':
            return m_FPrintBV[eC];
        case 'G':
            return m_FPrintBV[eG];
        case 'T':
            return m_FPrintBV[eT];
        case 'N':
            return m_FPrintBV[eN];
        default:
            break;
        }
        throw runtime_error("Error. Invalid letter!");
    }

private:
    bm::bvector<>   m_FPrintBV[eEnd];
};

/// Check correctness of indexes constructed using different methods
///
static
void fingerprint_compare(const DNA_FingerprintScanner& idx1,
                         const DNA_FingerprintScanner& idx2)
{
    std::vector<char> letters {'A', 'T', 'G', 'C'};
    for (char base : letters)
    {
        const bm::bvector<>& bv1 = idx1.GetVector(base);
        const bm::bvector<>& bv2 = idx2.GetVector(base);
        
        int cmp = bv1.compare(bv2);
        if (cmp != 0)
        {
            throw runtime_error(string("Fingerprint mismatch for:") + string(1, base));
        }
    } // for
}



int main(int argc, char *argv[])
{
    if (argc < 3)
    {
        show_help();
        return 1;
    }

    std::vector<char> seq_vect;

    try
    {
        auto ret = parse_args(argc, argv);
        if (ret != 0)
            return ret;

        DNA_FingerprintScanner idx1;
        DNA_FingerprintScanner idx2;

        if (!ifa_name.empty())
        {
            auto res = load_FASTA(ifa_name, seq_vect);
            if (res != 0)
                return res;
            std::cout << "FASTA sequence size=" << seq_vect.size() << std::endl;
            
            {
                bm::chrono_taker tt1(cout, "2. Build DNA index", 1, &timing_map);
                idx1.Build(seq_vect);
            }
            
            if (parallel_jobs > 0)
            {
                std::cout << "jobs = " << parallel_jobs << std::endl;
                bm::chrono_taker tt1(cout, "3. Build DNA index (bulk, parallel)", 1, &timing_map);
                idx2.BuildParallel(seq_vect, parallel_jobs);
            }

            // compare results (correctness check)
            //
            fingerprint_compare(idx1, idx2);
        }

        if (is_timing)  // print all collected timings
        {
            std::cout << std::endl << "Performance:" << std::endl;
            bm::chrono_taker<>::print_duration_map(cout, timing_map, bm::chrono_taker<>::ct_all);
        }
    }
    catch (std::exception& ex)
    {
        std::cerr << "Error:" << ex.what() << std::endl;
        return 1;
    }

    return 0;
}
