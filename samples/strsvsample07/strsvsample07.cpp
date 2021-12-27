/*
Copyright(c) 2002-2021 Anatoliy Kuznetsov(anatoliy_kuznetsov at yahoo.com)

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

/** \example strsvsample07.cpp
    Succinct container scanner search using pipeline to run thousands of searches faster
    one by one scans. scanner::pipeline uses variuous cache and algorithmic optimization techniques
    to run bulk searches faster.

  \sa bm::str_sparse_vector
  \sa bm::sparse_vector_scanner
  \sa bm::sparse_vector_scanner::pipeline
*/

/*! \file strsvsample07.cpp
    \brief Example: Succinct container for strings, bulk search using scanner pipeline
*/

#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include <cassert>
#include <thread>

#include "bm.h"
#include "bmstrsparsevec.h"
#include "bmsparsevec_algo.h"

#include "bmtimer.h"
#include "bmdbg.h"

#include "bmundef.h" /* clear the pre-proc defines from BM */

using namespace std;

typedef bm::bvector<> bvector_type;
typedef bm::str_sparse_vector<char, bvector_type, 8> str_sv_type;

/**

    Test data generation.
    max_coll - defines the number of string variants
    repeat_factor - how often strings should be duplicated (to simulate the compressable collections),
            higher repeat_factor produces more compressable vector.
 */
static
void GenerateTestData(std::vector<string>& str_coll,
                      str_sv_type&      str_sv,
                      unsigned max_coll = 8000000,
                      unsigned repeat_factor=10)
{
    // use back inserter to fill in succinct vector (it is faster than push_back)
    auto bi(str_sv.get_back_inserter());

    string str;
    for (unsigned i = 10; i < max_coll; i+= (rand()&0xF))
    {
        switch (i & 0xF)
        {
        case 0: str = "AB"; break;
        case 1: str = "GTx"; break;
        case 2: str = "cnv"; break;
        default: str = "AbY11"; break;
        }
        str.append(to_string(i));

        for (unsigned k = 0; k < repeat_factor; ++k)
        {
            str_coll.emplace_back(str);
            bi = str; // feed into SV back-inserter
        }
    } // for i
    bi.flush();
}

bool is_diag = true; ///< Flag to print the SV diagnostics

/// Rudimentary cmd-args parser
static
void parse_args(int argc, char *argv[])
{
    for (int i = 1; i < argc; ++i)
    {
        std::string arg = argv[i];
        if (arg == "-nodiag")
        {
            is_diag = false;
            continue;
        }
    } // for i
}




int main(int argc, char *argv[])
{
    try
    {
        parse_args(argc, argv); 


        std::vector<string> str_coll;
        str_sv_type str_sv0(bm::use_null); // sparse-succinct vector

        cout << "Generating the test data... " << flush;
        GenerateTestData(str_coll, str_sv0);
        str_sv_type str_sv1(str_sv0); // make a copy of the original vector

        cout << "OK" << endl;


        {
            cout << "Remapping the data to create compressed vector " << flush;
            BM_DECLARE_TEMP_BLOCK(tb)
            // apply char frequency remapping compression
            // (should not modify after that)
            str_sv0.remap();
            str_sv0.optimize(tb); // optimize the succinct vector

            cout << "OK" << endl;
        }

        // we have two succinct vectors str_sv0 (remapped and optimized)
        // and str_sv1 - original after construction
        //  - print statistics to take a look into details
        //

        if (is_diag)
        {
            cout << "\nStatistics on generated SV:" << endl;
            bm::print_svector_stat(cout, str_sv1);
            // diagnostics print to see the details of succinct structures
            cout << "\nStatistics on remapped/optimized SV:" << endl;
            bm::print_svector_stat(cout, str_sv0);
            cout << endl << endl;
        }

        // create a random sampling of strings to search
        //
        unsigned test_runs = 10000;
        std::vector<string> str_test_coll;
        for (bvector_type::size_type i = 0; i < test_runs; ++i)
        {
            bvector_type::size_type idx = (unsigned) rand() % test_runs;
            if (idx >= test_runs)
                idx = test_runs/2;
            str_test_coll.push_back(str_coll[idx]);
        }
        assert(str_test_coll.size() == test_runs);

        // -------------------------------------------------------------

        std::vector<unique_ptr<bvector_type> > res_vec1;
        bm::sparse_vector_scanner<str_sv_type> scanner;


        cout << "Running benchmark tests.." << endl;

        for (int pass = 0; pass < 2; pass++)
        {
            cout << "PASS = " << pass << ((pass==0) ? " -- remap/optimized" : " -- NOT remapped") << endl;

            res_vec1.resize(0);
            const str_sv_type* str_sv = (pass==0) ? &str_sv0 : &str_sv1;

            // Run experiment 1, search sparse vector in a loop one-by-one
            // This is not a slow method, scanner uses various optimizations
            // (SIMD, "search space prunning" to be efficient)
            {
                bm::chrono_taker tt(cout, "scanner<>::find_eq_str()", test_runs);

                for (bvector_type::size_type i = 0; i < test_runs; ++i)
                {
                    const string& str = str_test_coll[i];
                    unique_ptr<bvector_type> bv_res(new bvector_type);
                    scanner.find_eq_str(*str_sv, str.c_str(), *bv_res);
                    res_vec1.emplace_back(unique_ptr<bvector_type>(bv_res.release()));
                } // for
            }

            // There is a faster way to do the same.
            // if we know a bulk of our searches upfront we can use pipeline
            // schedule to run the search group.
            // Scanner pipeline will anayse the set and try to build a more
            // optimal search plan, taking into account CPU cache optimization,
            // resuse of compressed bit-blocks (inetrnal details) and other
            // factors between the search data set items and the sparce vector
            //

            // pipeline object to run the bulk search
            //
            bm::sparse_vector_scanner<str_sv_type>::pipeline<> pipe1(*str_sv);

            // batch_size instructs how many search to run at once
            // batch_size=0 and this parameter will be identified automatically
            //
            // batch size essentially depends on CPU cache size and it is
            // sometimes difficult to determine without trying.
            // batch_size=0 will try to use euristics for CPU L2 = 256KB,
            // it may be "good enough", but for best results it is best
            // to run trial runs with typical values may be (2 to 20)
            //
            pipe1.options().batch_size = test_runs;
            {
                bm::chrono_taker tt(cout, "scanner::pipeline find_eq_str()", test_runs);

                // add all the search items to the pipeline
                for (size_t i = 0; i < test_runs; ++i)
                {
                    const string& str = str_test_coll[i];
                    pipe1.add(str.c_str());
                }
                pipe1.complete(); // finish the pipeline construction with this call

                scanner.find_eq_str(pipe1); // run the search

                // at this point we have the results (in the pipeline object itself)
            }

            // This example sets the search range for only a part of
            // the sv vector using the mask bit-vector.
            // Since mask vector narrows down the search space - it is faster
            //
            bm::bvector<> bv_mask;
            bv_mask.set_range(0, str_sv->size()/3); // only first 1/3 of elements will be searched


            // scanner is configured to support result vectors AND search masking
            //
            typedef bm::agg_run_options<true, false, true> scanner_custom_mask_opt;
            bm::sparse_vector_scanner<str_sv_type>::pipeline<scanner_custom_mask_opt> pipe1_and(*str_sv);

            pipe1_and.set_search_mask(&bv_mask); // associate search mask with the pipeline
            pipe1_and.options().batch_size = test_runs;

            {
                bm::chrono_taker tt(cout, "scanner::pipeline+MASK find_eq_str()", test_runs);

                // add all the search items to the pipeline
                for (size_t i = 0; i < test_runs; ++i)
                {
                    const string& str = str_test_coll[i];
                    pipe1_and.add(str.c_str()); // this will search within defined mask
                }
                pipe1_and.complete();
                scanner.find_eq_str(pipe1_and); // run the search
            }


            // This variant would build the same pipeline but configure
            // it differently using pipe2.options()
            //
            // The idea here is that sometimes we don't actually need the result
            // vectors, but only need population count from this
            // pipeline can be configured to do that without actually
            // materializing result bit-vectors
            //

            // pipeline configuration passed via template parameter
            // instructs to drop results and provide only counts

            bm::sparse_vector_scanner<str_sv_type>::pipeline<bm::agg_opt_only_counts> pipe2(*str_sv);
            pipe1.options().batch_size = test_runs;

            {
                bm::chrono_taker tt(cout, "scanner::pipeline find_eq_str()-count()", test_runs);

                for (size_t i = 0; i < test_runs; ++i)
                {
                    const string& str = str_test_coll[i];
                    pipe2.add(str.c_str());
                }
                pipe2.complete(); // finish the pipeline construction with this call

                scanner.find_eq_str(pipe2); // run the search pipeline

                // at this point we have the population count
                // ... see below how to use it
            }

            // results from two pileline runs are ready at this point
            // we now get access to it
            //

            bvector_type bv_or_total;

            {
                // please note that returned buffer results vectors are NOT STL
                // vector<> type, but a thin wrapper over C style objects
                // (for the reasons of portability)
                //
                auto& res_vect = pipe1.get_bv_res_vector(); // results from pipeline 1
                auto& res_vect_and = pipe1_and.get_bv_res_vector();
                auto& cnt_vect = pipe2.get_bv_count_vector(); // counts from pileine 2

                assert(res_vect.size() == cnt_vect.size());

                // iterate over results, run some checks...
                size_t res_sz = res_vect.size();
                for (size_t i = 0; i < res_sz; ++i)
                {
                    const bvector_type* bv1 = res_vec1[i].get();
                    assert(bv1);
                    const bvector_type* bv = res_vect[i];
                    assert(bv);
                    bool match = bv1->equal(*bv); // quick check
                    assert(match);

                    auto c = cnt_vect[i];
                    auto cnt = bv->count();
                    (void)cnt; (void)c; // to silence unused warnings (relese)
                    assert(cnt == c); // check if counts match

                    bv_or_total |= *bv; // accumulate OR
                    {
                    auto c_and = bm::count_and(*bv, bv_mask);
                    const bvector_type* bv_and = res_vect_and[i];
                    if (bv_and)
                    {
                        auto c1 = bv_and->count();
                        assert(c1 == c_and); (void)c1; (void)c_and;
                        bvector_type bv_m;
                        bv_m.bit_and(*bv1, bv_mask);
                        match = bv_and->equal(bv_m);
                        assert(match);
                    }
                    else
                    {
                        assert(!c_and);
                    }


                    }
                }
            }

            // Here we define a csutom pipeline run policy which disables
            // both intermediate results and population counting for it...
            //
            // instead it aggregates the results into one UNION ALL (OR) vector
            // which simulates a huge
            // field1 IN ('value1', 'value2', 'value3', .... ) SQL expression
            //

            typedef bm::agg_run_options<false, false> scanner_custom_opt;
            bm::sparse_vector_scanner<str_sv_type>::pipeline<scanner_custom_opt> pipe3(*str_sv);
            pipe1.options().batch_size = test_runs;

            bvector_type bv_or;
            pipe3.set_or_target(&bv_or); // Assign OR aggregation target

            {
                bm::chrono_taker tt(cout, "scanner::pipeline find_eq_str()-OR()", test_runs);

                for (size_t i = 0; i < test_runs; ++i)
                {
                    const string& str = str_test_coll[i];
                    pipe3.add(str.c_str());
                }
                pipe3.complete(); // finish the pipeline construction with this call

                scanner.find_eq_str(pipe3); // run the search pipeline
            }
            bool match = bv_or.equal(bv_or_total);
            if (!match)
            {
                cerr << "OR vector mismatch!" << endl;
                exit(1);
            }




            cout << endl;

            // second pass will run the same benchmarks, only using original
            // non-remapped vector to see the effects of additional compression
            // on performance of scanner searches

        } // for pass


    }
    catch(std::exception& ex)
    {
        std::cerr << ex.what() << std::endl;
        return 1;
    }
    

    return 0;
}

