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

/** \example xsample09.cpp

    Then it computes variants of density histograms of using
    presense/absense (NOT NULL bit-vector) and bm::bvector<>::count_range()
    function.

    @sa bm::bvector
    @sa bm::bvector::find_reverse
    @sa bm::bvector::find

    @sa bm::rank_range_split

    @sa bm::sparse_vector
    @sa bm::rsc_sparse_vector
    @sa bm::rsc_sparse_vector::get_back_inserter
*/

/*! \file xsample09.cpp
    \brief Example: Use succinct vectors for histogram construction
*/


#include <iostream>
#include <memory>
#include <vector>
#include <random>
#include <algorithm>
#include <stdexcept>

using namespace std;

#include "bm.h"
#include "bmtimer.h"
#include "bmsparsevec.h"
#include "bmsparsevec_compr.h"

#include "bmdbg.h"

#include "bmundef.h" /* clear the pre-proc defines from BM */

// ----------------------------------------------------
// Global parameters and types
// ----------------------------------------------------

const unsigned  test_size = 250000000;  // number of events (ints) to generate
const unsigned  sampling_interval = 2500;   // size of the histogram sampling

typedef bm::bvector<>                                       bvector_type;
typedef bvector_type::size_type                             bv_size_type;
typedef bm::sparse_vector<unsigned, bvector_type>           sparse_vector_u32;
typedef bm::rsc_sparse_vector<unsigned, sparse_vector_u32>  rsc_sparse_vector_u32;

typedef std::vector<std::pair<bv_size_type, bv_size_type> > bv_ranges_vector;



// timing storage for benchmarking
bm::chrono_taker<>::duration_map_type  timing_map;


/// Generate a test RSC vector with a randomly distributed values
/// imitating distribution density of genome variations
/// it adds a huge area of empty in the middle to simulate chr centromere
static
void generate_test_data(rsc_sparse_vector_u32& csv, unsigned size)
{
    BM_DECLARE_TEMP_BLOCK(tb)
    unsigned cnt = 0;

    // part one: dense
    {
        auto bit = csv.get_back_inserter(); // fastest way to add the data
        for (unsigned i = 0; i < size/3; )
        {
            bit = cnt++;
            unsigned nc = (unsigned)(rand() % 16);
            bit.add_null(nc);

            i+=nc+1;
        } // for
        bit.flush();
    }
    // part two: empty + dense again
    {
        auto bit = csv.get_back_inserter(); // fastest way to add the data
        bit.add_null(size/3); // add huge emty space in the middle

        for (unsigned i = 0; i < size; )
        {
            bit = cnt++;
            unsigned nc = (unsigned)(rand() % 16);
            bit.add_null(nc);

            i+=nc+1;
        } // for
        bit.flush();
    }

    csv.optimize(tb); // memory compression
    csv.sync();     // construct the Rank-Select access index
}

///  generate list of random indexes (locations) to read histogram values
///
static
void generate_access_samples(std::vector<bvector_type::size_type> &sample_vec,
                             unsigned size)
{
    for (unsigned i = 0; i < size; )
    {
        unsigned sample = (unsigned)(rand() % 256); // random
        sample_vec.push_back(sample);
        i += sample;
    }
    std::random_device rd;
    std::mt19937 g(rd());

    // even more random (unordered)
    std::shuffle(sample_vec.begin(), sample_vec.end(), g);
}

/// Compute histogram as a SV vector using fixed sampling interval
///
static
void compute_historgam(sparse_vector_u32& hist_sv,
                       const rsc_sparse_vector_u32& csv,
                       sparse_vector_u32::size_type sampling_size)
{
    assert(sampling_size);
    assert(sampling_size < csv.size());

    BM_DECLARE_TEMP_BLOCK(tb)

    sparse_vector_u32::size_type from = 0;
    sparse_vector_u32::size_type to = sampling_size-1;

    // get NOT NULL vector
    const rsc_sparse_vector_u32::bvector_type* bv_null = csv.get_null_bvector();
    assert(bv_null);

    {
        auto bit = hist_sv.get_back_inserter();
        auto sz = csv.size();
        do
        {
            auto cnt = bv_null->count_range(from, to); // closed interval [from..to]
            bit = cnt;
            from += sampling_size;
            to = from + sampling_size - 1;
        } while (from < sz);
        bit.flush();
    }
    hist_sv.optimize(tb); // memory compaction
}

/// Compute histogram as a RSC vector using fixed sampling interval.
/// Histogram values are stored as "true" interval start coordinates and
/// it is a more flexible scheme if we eventually decide to use adaptive
/// sampling (variable step).
///
static
void compute_rsc_historgam(rsc_sparse_vector_u32& hist_rsc,
                           const rsc_sparse_vector_u32& csv,
                           sparse_vector_u32::size_type sampling_size)
{
    assert(sampling_size);
    assert(sampling_size < csv.size());

    {
    // important! bm::use_null if vector is later converted to RSC
    sparse_vector_u32 hist_sv(bm::use_null);

    rsc_sparse_vector_u32::size_type from = 0;
    rsc_sparse_vector_u32::size_type to = sampling_size-1;

    const rsc_sparse_vector_u32::bvector_type* bv_null = csv.get_null_bvector();
    assert(bv_null);

    auto sz = csv.size();
    do
    {
        auto cnt = bv_null->count_range(from, to); // closed interval [from..to]
        hist_sv.set(from, cnt); // assign historgam value at the interval start

        from += sampling_size;
        to = from + sampling_size - 1;
    } while (from < sz);

    hist_rsc.load_from(hist_sv);
    }

    BM_DECLARE_TEMP_BLOCK(tb)
    hist_rsc.optimize(tb);
    hist_rsc.sync();
}

/// Adaptive histogram identifies number of not NULL elements (events)
/// and varies the size of the histogram bin trying to make sure all
/// bins (but last) are the same weight
///
static
void compute_adaptive_rsc_histogram(rsc_sparse_vector_u32& hist_rsc,
                                const rsc_sparse_vector_u32& csv,
                                sparse_vector_u32::size_type sampling_size)
{
    assert(sampling_size);
    BM_DECLARE_TEMP_BLOCK(tb)

    const bvector_type* bv_null = csv.get_null_bvector();
    bv_size_type split_rank = sampling_size;

    bv_ranges_vector pair_vect;
    bm::rank_range_split(*bv_null, split_rank, pair_vect);

    size_t sz = pair_vect.size();
    for (size_t k = 0; k < sz; ++k)
    {
        const auto& p = pair_vect[k]; // [from..to] sampling rank interval
        hist_rsc.push_back(p.first, 0); // split_rank);
    } // for
    hist_rsc.optimize(tb);
    hist_rsc.sync();
}

/// Some test to confirm correctness
///
static
void verify_histograms(const rsc_sparse_vector_u32& hist_rsc,
                       const sparse_vector_u32& hist_sv,
                       sparse_vector_u32::size_type sampling_size)
{
    const rsc_sparse_vector_u32::bvector_type* bv_null = hist_rsc.get_null_bvector();
    auto en = bv_null->get_enumerator(0);
    for (;en.valid(); ++en)
    {
        auto idx = *en;
        rsc_sparse_vector_u32::size_type sv_idx = (idx / sampling_size);
        auto v1 = hist_rsc[idx];
        auto v2 = hist_sv[sv_idx];
        if (v1 != v2)
        {
            cerr << "Discrepancy at:" << idx << endl;
            exit(1);
        }
    } // for
}


/// Access benchmark 1
///
/// uses regular bit-transposed sparse vector to read
/// histogram values in random order. It relies on fixed inetrval sampling.
///
static
unsigned long long access_bench1(const sparse_vector_u32& hist_sv,
                    const std::vector<bvector_type::size_type>& sample_vec,
                    unsigned sampling_size)
{
    unsigned long long sum = 0;
    for (size_t i = 0; i < sample_vec.size(); ++i)
    {
        auto idx = sample_vec[i];
        idx = idx / sampling_size; // interval is calculated by division
        auto v = hist_sv[idx];
        sum += v;
    }
    return sum;
}

/// Access benchmark 2
///
/// Uses Rank-Select bit-transposed vector to read histogram values
/// Sampling interval can be non-fixed (variadic, adaptive sampling).
/// Method finds the interval start and value using RSC container not NULL vector
static
unsigned long long access_bench2(const rsc_sparse_vector_u32&   hist_rsc,
                    const std::vector<bvector_type::size_type>& sample_vec)
{
    const bvector_type* bv_null = hist_rsc.get_null_bvector();
    assert(bv_null);

    unsigned long long sum = 0;
    for (size_t i = 0; i < sample_vec.size(); ++i)
    {
        auto idx = sample_vec[i];

        // search back into container not NULL bit-vector to find sampling
        // interval start position
        bvector_type::size_type pos;
        bool found = bv_null->find_reverse(idx, pos);
        assert(found); (void)found;

        auto v = hist_rsc[pos];
        sum += v;
    }
    return sum;
}


static
void access_bench3(const rsc_sparse_vector_u32&   hist_rsc,
                   const std::vector<bvector_type::size_type>& sample_vec,
                   const rsc_sparse_vector_u32&                rsc_data)
{
    bvector_type::size_type last;
    // find the last element in the data-set
    {
        const bvector_type* bv_null = rsc_data.get_null_bvector();
        bool found = bv_null->find_reverse(last);
        assert(found); (void)found;
    }

    const bvector_type* bv_null = hist_rsc.get_null_bvector();
    assert(bv_null);

    for (size_t i = 0; i < sample_vec.size(); ++i)
    {
        auto idx = sample_vec[i];

        // search back into container not NULL bit-vector to find sampling
        // interval start position
        bvector_type::size_type pos_start, pos_end;
        bool found = bv_null->find_reverse(idx, pos_start);
        assert(found);

        found = bv_null->find(idx+1, pos_end);
        if (!found)
            pos_end = last;

        // please note that we don't need number of events here
        // (as it is always the same) we need the interval [start......end]
        //
        // (end - start) / sampling_rank defines the average density
        //
    }
}


int main(void)
{
    try
    {
        rsc_sparse_vector_u32                rsc_test;
        bv_size_type                          data_set_size = 0;

        sparse_vector_u32                    hist1;
        unsigned                             hist1_avg = 0;

        rsc_sparse_vector_u32                hist2;
        rsc_sparse_vector_u32                hist3;

        std::vector<bvector_type::size_type> svec;

        generate_test_data(rsc_test, test_size);

        {
            const bvector_type* bv_null = rsc_test.get_null_bvector();
            assert(bv_null);
            data_set_size = bv_null->count(); // number of elements in the data set
            cout << "Data set size:                             " << rsc_test.size() << endl;
            cout << "Number of elements/events in the data set: " << data_set_size << endl;
        }

        {
            bm::chrono_taker tt1(cout, "01. Histogram 1 construction (SV)) ", 1, &timing_map);
            compute_historgam(hist1, rsc_test, sampling_interval);
        }
        // explore some statistics on SV histogram 1
        {
            sparse_vector_u32::statistics st;
            hist1.calc_stat(&st);

            cout << "Histogram 1 (SV)" << endl;
            cout << "  size: " << hist1.size() << endl;
            cout << "  RAM size: " << st.memory_used << endl;
            {
                bm::chrono_taker tt1(cout, "05. Serialize and save the histogram 1 (SV)", 1, &timing_map);
                size_t serialized_size = 0;
                int res = bm::file_save_svector(hist1, "hist1.sv", &serialized_size, true);
                if (res!=0)
                    cerr << "Failed to save!" << endl;
                else
                    cout << "  file size: " << serialized_size << endl;
            }

            // calculate sum and average value in historgam 1
            //

            unsigned h1_sum = 0;
            auto it = hist1.begin();
            auto it_end = hist1.end();
            for (; it != it_end; ++it)
                h1_sum += *it;
            assert (h1_sum == data_set_size);
            hist1_avg = h1_sum / hist1.size();
        }

        {
            bm::chrono_taker tt1(cout, "02. Histogram 2 construction (RSC) ", 1, &timing_map);
            compute_rsc_historgam(hist2, rsc_test, sampling_interval);
        }
        {
            rsc_sparse_vector_u32::statistics st;
            hist2.calc_stat(&st);

            cout << "Histogram 2 (RSC)" << endl;
            cout << "  size: " << hist2.size() << endl;
            cout << "  NOT NULL count: " << hist2.get_null_bvector()->count() << endl;
            cout << "  RAM size: " << st.memory_used << endl;

            {
                bm::chrono_taker tt1(cout, "06. Serialize and save histogram 2(RSC)", 1, &timing_map);
                size_t serialized_size = 0;
                int res = bm::file_save_svector(hist2, "hist2.sv", &serialized_size, true);
                if (res!=0)
                    cerr << "Failed to save!" << endl;
                else
                    cout << "  file size: " << serialized_size << endl;
            }
        }

        {
            bm::chrono_taker tt1(cout, "03. Histogram 3 (adaptive) construction (RSC) ", 1, &timing_map);
            compute_adaptive_rsc_histogram(hist3, rsc_test, hist1_avg);
        }

        {
            rsc_sparse_vector_u32::statistics st;
            hist3.calc_stat(&st);

            cout << "Histogram 3 adaptive (RSC)" << endl;
            cout << "  size: " << hist3.size() << endl;
            cout << "  NOT NULL count: " << hist3.get_null_bvector()->count() << endl;
            cout << "  RAM size: " << st.memory_used << endl;
            {
                bm::chrono_taker tt1(cout, "07. Serialize and save the adaptive histogram 3 (RSC)", 1, &timing_map);
                size_t serialized_size = 0;
                int res = bm::file_save_svector(hist3, "hist3.sv", &serialized_size, true);
                if (res!=0)
                    cerr << "Failed to save!" << endl;
                else
                    cout << "  file size: " << serialized_size << endl;
            }
        }

        generate_access_samples(svec, test_size);
        cout << endl;
        cout << "Access sample size = " << svec.size() << endl;

        {
            bm::chrono_taker tt1(cout, "04. Verification", 1, &timing_map);
            verify_histograms(hist2, hist1, sampling_interval);
        }


        unsigned long long sum1(0), sum2(0);
        {
            bm::chrono_taker tt1(cout, "08. Random access test H1 (SV)", 1, &timing_map);
            sum1 = access_bench1(hist1, svec, sampling_interval);
        }

        {
            bm::chrono_taker tt1(cout, "09. Random access test H2 (RSC)", 1, &timing_map);
            sum2 = access_bench2(hist2, svec);
        }

        if (sum1 != sum2) // paranoiya check
        {
            cerr << "Control sum discrepancy!" << endl;
            return 1;
        }

        {
            bm::chrono_taker tt1(cout, "10. Random access test H3 adaptive (RSC)", 1, &timing_map);
            access_bench3(hist3, svec, rsc_test);
        }

        cout << endl;
        bm::chrono_taker<>::print_duration_map(cout, timing_map, bm::chrono_taker<>::ct_ops_per_sec);

    }
    catch(std::exception& ex)
    {
        std::cerr << ex.what() << std::endl;
        return 1;
    }

    return 0;
}

