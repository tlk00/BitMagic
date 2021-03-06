 /*
Copyright(c) 2002-2018 Anatoliy Kuznetsov(anatoliy_kuznetsov at yahoo.com)

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


//#define BMSSE2OPT
//#define BMSSE42OPT
//#define BMAVX2OPT
//#define BM_USE_EXPLICIT_TEMP
//#define BM_USE_GCC_BUILD

#define BMXORCOMP
#define BM_NONSTANDARD_EXTENTIONS

#include <stdio.h>
#include <stdlib.h>
#undef NDEBUG
#include <cassert>
#include <time.h>
#include <math.h>
#include <string.h>

#include <iostream>
#include <iomanip>
#include <utility>
#include <memory>
#include <random>
#include <algorithm>
#include <stdarg.h>
#include <mutex>

#include <functional>

#include <bm.h>
#include <bmalgo.h>
#include <bmxor.h>
#include <bmaggregator.h>
#include <bmutil.h>
#include <bmserial.h>
#include <bmrandom.h>
#include <bmvmin.h>
#include <bmbmatrix.h>
#include <bmintervals.h>
#include <bmsparsevec.h>
#include <bmsparsevec_algo.h>
#include <bmsparsevec_serial.h>
#include <bmalgo_similarity.h>
#include <bmsparsevec_util.h>
#include <bmsparsevec_compr.h>
#include <bmstrsparsevec.h>
#include <bmtimer.h>
#include <bmtask.h>
#include <bmsparsevec_parallel.h>
#include <bmthreadpool.h>

using namespace bm;
using namespace std;

#include <limits.h>

#include <vector>

// ----------------------------------------------------
// Globals
//

unsigned                                num_threads = 4;

bm::chrono_taker::duration_map_type     timing_map;

// ----------------------------------------------------


typedef bm::bvector<> bvect;

typedef std::vector<unsigned> ref_vect_type;

typedef bm::sparse_vector<unsigned, bvect > sparse_vector_u32;
typedef bm::sparse_vector<unsigned long long, bvect > sparse_vector_u64;
typedef bm::rsc_sparse_vector<unsigned, sparse_vector_u32> rsc_sparse_vector_u32;
typedef bm::rsc_sparse_vector<unsigned long long, sparse_vector_u64> rsc_sparse_vector_u64;


typedef
bm::thread_pool<bm::task_description*, bm::spin_lock<bm::pad0_struct> > pool_spin_type;
typedef
bm::thread_pool<bm::task_description*, std::mutex > pool_mutex_type;



static
void TestSimpleQueue()
{
    cout << " ----------------------------- TestSimpleQueue()" << endl;

    typedef bm::simple_queue<unsigned, bvect::allocator_type, true> squeue_u32;

    {
        squeue_u32 sq;
        squeue_u32::value_type val;

        assert(sq.empty());
        assert(sq.size()==0);

        bool b = sq.try_push(0);
        assert(!b);
        b = sq.try_pop(val);
        assert(!b);
    }

    {
        squeue_u32 sq;
        squeue_u32::value_type val;

        sq.reserve(1);
        bool b;
        b = sq.try_push(0);
        assert(b);
        auto v = sq.front();
        assert(v == 0);
        assert(!sq.empty());
        b = sq.try_push(1);
        assert(b);
        b = sq.try_push(2);
        assert(!b);
        sq.push(2); // try again with resize
        auto sz = sq.size();
        assert(sz == 3);

        for (unsigned i = 0; i < 3; ++i)
        {
            b = sq.try_pop(val);
            assert(b);
            assert(val == i);
        } // for
        sz = sq.size();
        assert(sz == 0);
        b = sq.try_pop(val);
        assert(!b);
        assert(sq.empty());

        b = sq.try_push(0);
        assert(b);
        b = sq.try_pop(val);
        assert(b);
        assert(val == 0);
    }

    cout << " ----------------------------- TestSimpleQueue() OK" << endl;
}


static
void Check_SimModel(const bm::xor_sim_model<bvect>& sim1,
                    const bm::xor_sim_model<bvect>& sim2)
{
    bool eq = sim1.bv_blocks.equal(sim2.bv_blocks);
    if (!eq)
    {
        cerr << "Sim model Blocks vectors does not match! (BV)" << std::endl;
        assert(0); exit(1);
    }
    eq = sim1.matr.equal(sim2.matr);
    if (!eq)
    {
        cerr << "Sim model Blocks vectors does not match! (Matr)" << std::endl;
        std::this_thread::sleep_for (std::chrono::seconds(1));
        eq = sim1.matr.equal(sim2.matr);
        if (eq)
        {
            cerr << "Race!" << endl;
        }

        auto cols = sim1.matr.cols();
        auto rows = sim1.matr.rows();
        for (unsigned i = 0; i < rows; ++i)
        {
            for (unsigned j = i+1; j < cols; ++j)
            {
                auto v1 = sim1.matr.get(i,j);
                auto v2 = sim2.matr.get(i, j);
                if (!(v1 == v2))
                {
                    std::cerr << " Mismatch at: " << i << ":" << j << endl;
                    if (v1.chain_size != v2.chain_size)
                        cerr << "Chain size mismatch" << endl;
                    std::cerr << "chains=" << v1.chain_size << " " << v2.chain_size << endl;

                    if (v1.nb != v2.nb)
                        std::cerr << "NB mismatch!" << endl;
                    if (v1.match != v2.match)
                        std::cerr << "XOR match type  mismatch!" << endl;
                    cerr << "Match type = " << v1.match << endl;

                    auto chain_size = v1.chain_size;
                    assert(chain_size < 64);
                    for (unsigned k = 0; k < chain_size; ++k)
                        if (v1.ref_idx[k] != v2.ref_idx[k])
                            std::cerr << "refs " << v1.ref_idx[k] << " " << v2.ref_idx[k] << " ";
                    std::cerr << endl;
                    for (unsigned k = 0; k < chain_size; ++k)
                        if (v1.xor_d64[k] != v2.xor_d64[k])
                            std::cerr << "D64 " << v1.xor_d64[k] << " " << v2.xor_d64[k];
                    std::cerr << endl;

                }
            } // j
        } // i

        assert(0); exit(1);
    }
}

static
void generate_df_svu32_set(sparse_vector_u32& sv1i, 
                           sparse_vector_u32& sv2i, 
                           sparse_vector_u32& sv3i, 
                           sparse_vector_u32& sv4i,
                           bvect::size_type size_max)
{
    auto bi1(sv1i.get_back_inserter());
    auto bi2(sv2i.get_back_inserter());
    auto bi3(sv3i.get_back_inserter());
    auto bi4(sv4i.get_back_inserter());

    for (bvect::size_type i = 0; i < size_max; i += 2)
    {
        bi1 = 4 + 8 + 16 + 32 + 64 + 128 + 512 + 1024 + 2048;
        bi1 = 0;
        bi2 = 2 + 8;
        bi2 = 0;
        bi3 = 1 + 0;
        bi3.add_null();
        bi4 = rand() % 256000;
        bi4.add_null();

    } // for

    bi1.flush();
    bi2.flush();
    bi3.flush();
    bi4.flush();

}

template<typename PoolType>
void TestParallelSV_Serial(const char* test_label,
                           bool        check_sim,
                           bool        serialize,
                           bool        deserialize)
{
    cout << " ----------------------------- TestParallelSV_Serial() - "
         << test_label << endl;
    {
        bm::sparse_vector_serializer<sparse_vector_u32> sv_serializer;
        bm::sparse_vector_deserializer<sparse_vector_u32> sv_deserial;

        sparse_vector_u32 sv1i, sv2i, sv3i(bm::use_null), sv4i(bm::use_null);
        sparse_vector_u32 sv1o, sv2o, sv3o(bm::use_null), sv4o(bm::use_null);

        bm::sparse_vector_serial_layout<sparse_vector_u32> sv_lay1, sv_lay2, sv_lay3, sv_lay4;
        bm::sparse_vector_serial_layout<sparse_vector_u32> sv_lay1p, sv_lay2p, sv_lay3p, sv_lay4p;

        typedef PoolType pool_type;

        {
        cout << "Sample set generation..." << flush;

        const bvect::size_type idmax = 50000000;
        generate_df_svu32_set(sv1i, sv2i, sv3i, sv4i, idmax);

        cout << "Ok" << endl;
        }

        bm::xor_sim_params       xor_search_params;

        pool_type tpool;  // our thread pool here (no threads created yet)
        tpool.start(num_threads); // start the threads

        bm::chrono_taker tt1(test_label, 1, &timing_map);

        for (unsigned pass = 0; pass < 2; ++pass)
        {
            bm::sparse_vector_serializer<sparse_vector_u32>::bv_ref_vector_type bv_ref;
            // add references in reverse(!) order
            bv_ref.add_vectors(sv4i.get_bmatrix());
            bv_ref.add_vectors(sv3i.get_bmatrix());
            bv_ref.add_vectors(sv2i.get_bmatrix());
            bv_ref.add_vectors(sv1i.get_bmatrix());

            sv_serializer.set_xor_ref(&bv_ref);
            assert(sv_serializer.is_xor_ref());

            bm::sparse_vector_serializer<sparse_vector_u32>::xor_sim_model_type sim_model;

            {
                bm::compute_sim_matrix_plan_builder<bvect> pbuilder;
                bm::compute_sim_matrix_plan_builder<bvect>::task_batch tbatch;

                pbuilder.build_plan(tbatch, sim_model,
                                    bv_ref, xor_search_params);

                {
                    bm::thread_pool_executor<pool_type> exec;
                    exec.run(tpool, tbatch, true);
                }
                
                if (check_sim)
                {
                    bm::serializer<bvect>::xor_sim_model_type sim_model_c;
                    sv_serializer.compute_sim_model(sim_model_c, bv_ref, xor_search_params);
                    Check_SimModel(sim_model_c, sim_model);
                }
            }

            if (serialize)
            {
                sv_serializer.set_sim_model(&sim_model);

                bm::sv_serialization_plan_builder<sparse_vector_u32> sbuilder;
                sbuilder.set_xor_ref(&bv_ref);
                sbuilder.set_sim_model(&sim_model);

                {
                    bm::sv_serialization_plan_builder<sparse_vector_u32>::task_batch tbatch;
                    sbuilder.build_plan(tbatch, sv_lay1p, sv1i);
                }

                sv_serializer.serialize(sv1i, sv_lay1);
                {
                    const bvect::size_type* cstat = sv_serializer.get_bv_serializer().get_compression_stat();
                    assert(cstat[bm::set_block_ref_eq]>0);
                }

                sv_serializer.serialize(sv2i, sv_lay2);
                {
                    const bvect::size_type* cstat = sv_serializer.get_bv_serializer().get_compression_stat();
                    assert(cstat[bm::set_block_ref_eq]>=1 || cstat[bm::set_block_xor_ref32] >= 1);
                }

                sv_serializer.serialize(sv3i, sv_lay3);
                {
                    const bvect::size_type* cstat = sv_serializer.get_bv_serializer().get_compression_stat();
                    assert(cstat[bm::set_block_ref_eq]>=1 || cstat[bm::set_block_xor_ref32] >= 1);
                }

                sv_serializer.serialize(sv4i, sv_lay4);
                {
                    //const bvect::size_type* cstat = sv_serializer.get_bv_serializer().get_compression_stat();
                }
            }

            // ----------

            if (deserialize)
            {
                bm::sparse_vector_deserializer<sparse_vector_u32>::bv_ref_vector_type bv_ref_d;

                const unsigned char* buf = sv_lay1.buf();
                auto sz2 = sv_lay1.size();
                assert(sz2); (void)sz2;

                sv_deserial.deserialize_structure(sv1o, sv_lay1.buf());
                sv_deserial.deserialize_structure(sv2o, sv_lay2.buf());
                sv_deserial.deserialize_structure(sv3o, sv_lay3.buf());
                sv_deserial.deserialize_structure(sv4o, sv_lay4.buf());

                bv_ref_d.add_vectors(sv4o.get_bmatrix());
                bv_ref_d.add_vectors(sv3o.get_bmatrix());
                bv_ref_d.add_vectors(sv2o.get_bmatrix());
                bv_ref_d.add_vectors(sv1o.get_bmatrix());

                sv_deserial.set_xor_ref(&bv_ref_d);

                sv_deserial.deserialize(sv1o, buf, false);
                bool eq = sv1i.equal(sv1o);
                assert(eq);

                buf = sv_lay2.buf();
                sz2 = sv_lay2.size();

                sv_deserial.deserialize(sv2o, buf, false);
                eq = sv2i.equal(sv2o);
                assert(eq);

                buf = sv_lay3.buf();
                sz2 = sv_lay3.size();

                sv_deserial.deserialize(sv3o, buf, false);
                eq = sv3i.equal(sv3o);
                assert(eq);

                buf = sv_lay4.buf();

                sv_deserial.deserialize(sv4o, buf, false);
                eq = sv4i.equal(sv4o);
                assert(eq);
                sv_deserial.set_xor_ref(0); // unset
            }

            if (pass == 0)
            {
                sv1i.optimize();
                sv2i.optimize();
                sv3i.optimize();
                sv4i.optimize();
            }

        } // for pass

        tpool.set_stop_mode(pool_type::stop_when_done);
        tpool.join();
    }

    cout << " ----------------------------- TestParallelSV_Serial() OK" << endl;
}

static
void show_help()
{
    std::cerr
        << "BitMagic C++ stress test for parallel algorithms " << endl
        << "-h                      - help" << endl
        << "-j <number-of-threads>  - set number of threads to run" << endl
        ;
}


static
int parse_args(int argc, char *argv[])
{
    for (int i = 1; i < argc; ++i)
    {
        std::string arg = argv[i];
        if ((arg == "-h") || (arg == "--h"))
        {
            show_help();
            return 1;
        }
        if (arg == "-j")
        {
            if (i + 1 < argc)
            {
                const char* nt = argv[++i];
                char *end;
                num_threads = (unsigned) std::strtoul(nt, &end, 10);
            }
            else
            {
                std::cerr << "Error: -j requires number of threads" << std::endl;
                return 1;
            }
            continue;
        }
    } // for
    return 0;
}

int main(int argc, char *argv[])
{
    cout << bm::_copyright<true>::_p << endl;

    typedef std::function<int(void)> job_func_type;


    std::function<int(void)> f1([]()
        {
            cout << "func-type1" << endl; return 0;
        });
    int a = 10;
    int b = 20;
    std::function<int(void)> f2([a, b]()
        {
            cout << "func-type2 " << a << " " << b << endl; return 0;
        });
    cout << sizeof(f2) << endl;
/*
    std::vector<job_func_type> job_vect;
    job_vect.push_back(f1);
    job_vect.push_back(f2);

    for (size_t i = 0; i < job_vect.size(); ++i)
    {
        int r = job_vect[i]();
        cout << "returned: " << r << endl;
    }
*/
    std::deque<job_func_type> job_queue;
    job_queue.push_back(f1);
    job_queue.push_back(f2);

    while(!job_queue.empty()) {
        auto func = job_queue.front();
        int r = func();
        cout << "returned: " << r << endl;
        job_queue.pop_front();
    } // while


//    f1();
//    f2();



    return 0;


    {
        auto ret = parse_args(argc, argv);
        if (ret != 0)
            return ret;
    }
    cout << "Threads:" << num_threads << endl;

    bm::chrono_taker tt("TOTAL", 1);
    try
    {
        cout << endl;


        TestSimpleQueue();


        TestParallelSV_Serial<pool_spin_type>(
            "001-s. XOR filter (check)", true, false, false);

        TestParallelSV_Serial<pool_spin_type>(
            "002-s. XOR filter (no-check)", false, false, false);

        TestParallelSV_Serial<pool_spin_type>(
            "003-s. XOR filter/serialization (no-check)", false, true, false);

        TestParallelSV_Serial<pool_spin_type>(
            "004-s. XOR filter/serialization/deserialization (no-check)",
            false, true, true);


        TestParallelSV_Serial<pool_mutex_type>(
            "001-m. XOR filter (check)", true, false, false);

        TestParallelSV_Serial<pool_mutex_type>(
            "002-m. XOR filter (no-check)", false, false, false);

        TestParallelSV_Serial<pool_mutex_type>(
            "003-m. XOR filter/serialization (no-check)", false, true, false);

        TestParallelSV_Serial<pool_mutex_type>(
            "004-m. XOR filter/serialization/deserialization (no-check)",
            false, true, true);

        {
            std::cout << std::endl << "Performance:" << std::endl;
            bm::chrono_taker::print_duration_map(timing_map, bm::chrono_taker::ct_time);
        }

        cout << endl;
    }
    catch (std::exception& ex)
    {
        cerr << ex.what() << endl;
        exit(1);
    }
    return 0;

}

