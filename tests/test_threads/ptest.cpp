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

typedef bm::bvector<> bvect;

typedef std::vector<unsigned> ref_vect_type;

typedef bm::sparse_vector<unsigned, bvect > sparse_vector_u32;
typedef bm::sparse_vector<unsigned long long, bvect > sparse_vector_u64;
typedef bm::rsc_sparse_vector<unsigned, sparse_vector_u32> rsc_sparse_vector_u32;
typedef bm::rsc_sparse_vector<unsigned long long, sparse_vector_u64> rsc_sparse_vector_u64;


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
void TestParallelSV_Serial()
{
    cout << " ----------------------------- TestParallelSV_Serial()" << endl;
    {
        bm::sparse_vector_serializer<sparse_vector_u32> sv_serializer;
        bm::sparse_vector_deserializer<sparse_vector_u32> sv_deserial;

        sparse_vector_u32 sv1i, sv2i, sv3i(bm::use_null);
        sparse_vector_u32 sv1o, sv2o, sv3o(bm::use_null);

        bm::sparse_vector_serial_layout<sparse_vector_u32> sv_lay1, sv_lay2, sv_lay3;

        typedef
        bm::thread_pool<bm::task_description*, bm::spin_lock<bm::pad0_struct> > pool_type;
        pool_type tpool;  // our thread pool here (no threads created yet)
        tpool.start(4); // start the threads

        cout << "Sample set generation..." << flush;
        bvect::size_type idmax = 50000000;
        for (bvect::size_type i = 0; i < idmax; i+=2)
        {
            sv1i[i] = 4+8+16+32+64+128+512+1024+2048;
            sv2i[i] = 2+8;
            sv3i[i] = 1+0;
        }
        cout << "Ok" << endl;

        for (unsigned pass = 0; pass < 2; ++pass)
        {
            bm::sparse_vector_serializer<sparse_vector_u32>::bv_ref_vector_type bv_ref;
            // add references in reverse(!) order
            bv_ref.add_vectors(sv3i.get_bmatrix());
            bv_ref.add_vectors(sv2i.get_bmatrix());
            bv_ref.add_vectors(sv1i.get_bmatrix());

            sv_serializer.set_xor_ref(&bv_ref);
            assert(sv_serializer.is_xor_ref());

            bm::sparse_vector_serializer<sparse_vector_u32>::xor_sim_model_type sim_model;

            {
                bm::compute_sim_matrix_plan_builder<bvect> pbuilder;
                bm::compute_sim_matrix_plan_builder<bvect>::task_batch tbatch;
                bm::xor_sim_params       xor_search_params;

                pbuilder.build_plan(tbatch, sim_model,
                                    bv_ref, xor_search_params);

                {
                    bm::thread_pool_executor<pool_type> exec;
                    exec.run(tpool, tbatch, true);
                }
                {
                    bm::serializer<bvect>::xor_sim_model_type sim_model_c;
                    sv_serializer.compute_sim_model(sim_model_c, bv_ref, xor_search_params);
                    Check_SimModel(sim_model_c, sim_model);
                }
            }

            sv_serializer.set_sim_model(&sim_model);

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

            // ----------


            bm::sparse_vector_deserializer<sparse_vector_u32>::bv_ref_vector_type bv_ref_d;

            const unsigned char* buf = sv_lay1.buf();
            auto sz2 = sv_lay1.size();


            sv_deserial.deserialize_structure(sv1o, sv_lay1.buf());
            sv_deserial.deserialize_structure(sv2o, sv_lay2.buf());
            sv_deserial.deserialize_structure(sv3o, sv_lay3.buf());

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


            sv_deserial.set_xor_ref(0); // unset

            sv1i.optimize();
            sv2i.optimize();
            sv3i.optimize();

        } // for pass
        
        tpool.set_stop_mode(pool_type::stop_when_done);
        tpool.join();

    }

    cout << " ----------------------------- TestParallelSV_Serial() OK" << endl;
}


int main(void)
{
    cout << bm::_copyright<true>::_p << endl;
//    ptest();

    bm::chrono_taker tt("TOTAL", 1);
    try
    {
        cout << endl;

        TestParallelSV_Serial();

        cout << endl;
    }
    catch (std::exception& ex)
    {
        cerr << ex.what() << endl;
        exit(1);
    }
    return 0;

}
