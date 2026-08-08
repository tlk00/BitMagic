/*
Copyright(c) 2019 Anatoliy Kuznetsov(anatoliy_kuznetsov at yahoo.com)

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

#include <bitset>
#include <iostream>
#include <time.h>
#include <stdio.h>
#include <sstream>
#include <cassert>

//#define BMSSE2OPT
//#define BMSSE42OPT
//#define BMAVX2OPT

#ifdef _MSC_VER
#pragma warning( push )
#pragma warning( disable : 4996)
#endif

#include <vector>
#include <random>
#include <memory>

#include "bm64.h"
#include "bmintervals.h"
#include "bmsparsevec.h"
#include "bmstrsparsevec.h"
#include "bmsparsevec_algo.h"
#include "bmsparsevec_serial.h"
#include "bmtimer.h"
#include "bmdbg.h"

#include "bmsparsevec_float.h"
#include "bmsparsevec_compr.h"

using namespace std;
using namespace bm;

#include "../stress64/gena.h"
#include "../stress64/test_util.h"

//#include "bmdbg.h"

#include <math.h>

using namespace std;

const unsigned int BSIZE = 150000000;
const unsigned int REPEATS = 300;

//typedef  bitset<BSIZE>  test_bitset;


std::random_device rand_dev;
std::mt19937 gen(rand_dev()); // mersenne_twister_engine 
std::uniform_int_distribution<> rand_dis(0, BSIZE); // generate uniform numebrs for [1, vector_max]

bm::chrono_taker<>::duration_map_type  timing_map;

typedef bm::bvector<> bvect64;
typedef std::vector<bm::id64_t> ref_vect;

static
void bvector64_VerySparse_SetDestroyCycle()
{
    ref_vect vect;
    generate_vect_simpl0(vect);
    
    cout << "bvector64_VerySparse_SetDestroyCycle..." << endl;
    {
        bm::chrono_taker tt(std::cout, "bvector64_VerySparse_SetDestroyCycle", 1, &timing_map);
        const unsigned max_try = REPEATS;
        for (unsigned i = 0; i < max_try; ++i)
        {
            bvect64 bv0;
            load_BV_set_ref(std::cout, bv0, vect, false);
            if ((i & 0xF) == 0)
                cout << "\r" << i << "/" << max_try << flush;
        }
    }
    cout << "\r                       " << endl;
}

static
void bvector64_VerySparse_RefAccessCycle()
{
    {
        bvect64 bv0;
        ref_vect vect;
        generate_vect_simpl0(vect);
        load_BV_set_ref(std::cout, bv0, vect, false);

        cout << "bvector64_VerySparse_RefAccessCycle..." << endl;
        {
            bm::chrono_taker tt(std::cout, "bvector64_VerySparse_RefAccessCycle", 1, &timing_map);
            const unsigned max_try = REPEATS;
            for (unsigned i = 0; i < max_try; ++i)
            {
                compare_BV_set_ref(bv0, vect, false);
                if ((i & 0xF) == 0)
                    cout << "\r" << i << "/" << max_try << flush;
            }
        }
        cout << "\r                       " << endl;
        cout << "bvector64_VerySparse_EnumAccessCycle..." << endl;
        {
            bm::chrono_taker tt(std::cout, "bvector64_VerySparse_EnumAccessCycle", 1, &timing_map);
            const unsigned max_try = REPEATS;
            for (unsigned i = 0; i < max_try; ++i)
            {
                compare_BV(bv0, vect, false);
                if ((i & 0xF) == 0)
                    cout << "\r" << i << "/" << max_try << flush;
            }
        }
        cout << "\r                       " << flush;
    }
}

static
void bvector64_Serialization()
{
    {
        bvect64 bv0, bv1;
        {
            ref_vect vect0, vect1;
            
            generate_vect_simpl0(vect0);
            generate_vect48(vect1);
            
            load_BV_set_ref(cout,bv0, vect0, false);
            load_BV_set_ref(cout,bv1, vect1, false);
        }

        bm::serializer<bvect64> bv_ser;
        bm::serializer<bvect64>::buffer sbuf0, sbuf1;
        bv_ser.serialize(bv0, sbuf0, 0);
        bv_ser.serialize(bv1, sbuf1, 0);

        bvect64 bv_tc;
        bv_tc.bit_or(bv0, bv1, bvect64::opt_none);

        // interity check
        {
            bvect64 bv_t0;
            bm::deserialize(bv_t0, sbuf0.buf());
            int res = bv0.compare(bv_t0);
            assert(res==0);
            bvect64 bv_t1;
            bm::deserialize(bv_t1, sbuf1.buf());
            res = bv1.compare(bv_t1);
            assert(res==0);
        }
        {
            bvect64 bv_t;
            operation_deserializer<bvect64> od;
            od.deserialize(bv_t,
                                                         sbuf0.buf(),
                                                         0,
                                                         set_OR);
            od.deserialize(bv_t,
                                                         sbuf1.buf(),
                                                         0,
                                                         set_OR);
            int cmp = bv_tc.compare(bv_t);
            if (cmp)
            {
                std::cerr << "Serialization intergity check failed!" << std::endl;
                assert(0); exit(1);
            }
        }
        {
            bvect64 bv_t;
            bm::deserialize(bv_t, sbuf0.buf());
            int res = bv0.compare(bv_t);
            assert(res == 0); (void)res;
            
            bvect64 bv_t_c(bv_t);
            
            bm::deserialize(bv_t, sbuf1.buf());
            int cmp = bv_tc.compare(bv_t);
            if (cmp)
            {
                cmp = bv_tc.compare(bv_t);
                bm::deserialize(bv_t_c, sbuf1.buf());

                std::cerr << "Serialization intergity check failed!" << cmp << std::endl;
                assert(0); exit(1);
            }
        }
        cout << "bvector64_DeserializationCycle..." << endl;
        {
            bm::chrono_taker tt(std::cout, "bvector64_DeserializationCycle", 1, &timing_map);
            const unsigned max_try = REPEATS;
            for (unsigned i = 0; i < max_try; ++i)
            {
                bvect64 bv_t;
                bm::deserialize(bv_t, sbuf0.buf());
                bm::deserialize(bv_t, sbuf1.buf());
                if ((i & 0xF) == 0)
                    cout << "\r" << i << "/" << max_try << flush;
            }
        }
        cout << "\r                       " << endl;
    }
}


typedef bm::sparse_vector_float<bm::sparse_vector<unsigned int, bvect64>> sparseVecFloat;
typedef bm::sparse_vector<unsigned int, bvect64> sparse_vec_u32;
typedef bm::sparse_vector_float<bm::rsc_sparse_vector<unsigned int, sparse_vec_u32>> sparseVecFloatRSC;

//Finds all values in range [from, to] in a given std::vector<float> and flipts the corresponding bits in bv_out
inline
void in_range_vect(const std::vector<float>& fv, float from, float to, sparseVecFloat::bvector_type &bv_out)
{
    if(from > to) std::swap(from, to);
    for (sparseVecFloat::size_type i = 0; i < fv.size(); i++)
    {
        if (fv[i] >= from && fv[i] <= to)
            bv_out.set(i);
    } // for
}

//Finds all values in range [from, to] in a given sparse_vector_float using a const_iterator and flips the corresponding bits in bv_out
inline
void in_range_const(const sparseVecFloat& sv, float from, float to, sparseVecFloat::bvector_type& bv_out)
{
    sparseVecFloat::const_iterator ci = sv.begin();
    if (from > to) std::swap(from, to);
    for (; ci.valid(); ++ci)
    {
        if (auto v = ci.value(); (v >= from && v <= to))
            bv_out.set(ci.pos());
    }
}

struct svf_range_for_each_sparse_func
{
    typedef sparseVecFloat::size_type size_type;

    svf_range_for_each_sparse_func(float from, float to,
                                   sparseVecFloat::bvector_type& bv_out)
        : from_(from),
          to_(to),
          bv_out_(bv_out)
    {}

    int operator()(float v, bool is_null, size_type idx)
    {
        if (!is_null && v >= from_ && v <= to_)
            bv_out_.set(idx);
        return 0;
    }

    float from_;
    float to_;
    sparseVecFloat::bvector_type& bv_out_;
};

inline
void in_range_for_each_sparse(const sparseVecFloat& sv, float from, float to,
                              sparseVecFloat::bvector_type& bv_out)
{
    if (from > to) std::swap(from, to);
    svf_range_for_each_sparse_func func(from, to, bv_out);
    bm::for_each_sparse(sv, func);
}

void TestSVFScanner()
{
    BM_DECLARE_TEMP_BLOCK(tb)

    typedef sparseVecFloat::bvector_type bvect_type;
    typedef bvect_type::allocator_type::allocator_pool_type bvect_pool_type;
    bvect_pool_type bv_pool;

    sparseVecFloat::size_type N = 20000000;
    std::random_device rd;
    //std::mt19937 gen(rd());

    float upper = 1000000.0f;
    float lower = -1000000.0f;
    std::uniform_real_distribution<float> dis(lower, upper);

    std::vector<float> linData(N);

    for(sparseVecFloat::size_type i = 0; i < N/2; i++)
        linData[i] = -1.0f * (float)i * 0.00123f;
    for(sparseVecFloat::size_type i = 0; i < N/2; i++)
        linData[i+N/2] = (float)i * 0.00123f;

    sparseVecFloat testSVF;
    testSVF.import(linData.data(), N);
    testSVF.optimize(tb);

    unsigned int tests = 1000;

    {
        sparseVecFloat::bvector_type xorSV;
        sparseVecFloat::bvector_type xorVect;
        sparseVecFloat::bvector_type xorConst;
        sparseVecFloat::bvector_type xorForEachSparse;

        sparseVecFloat::bvector_type bv_range;
        sparseVecFloat::bvector_type::mem_pool_guard bv_range_guard(bv_pool, bv_range);

        std::vector<float> fromVect(tests);
        std::vector<float> toVect(tests);
        for (unsigned int i = 0; i < tests; i++)
        {
            fromVect[i] = dis(gen);
            toVect[i] = dis(gen);
        }

        {
            bm::chrono_taker<> tt(cout, "SVF with Linear Data find values in range with scanner", tests);
            bm::sparse_vector_scanner<sparseVecFloat> scan;
            for (unsigned int i = 0; i < tests; i++)
            {
                float from = fromVect[i];
                float to   = toVect[i];

                scan.find_range_float(testSVF, from, to, bv_range);

                xorSV ^= bv_range;
                bv_range.clear(true);
            }
        }

        {
            bm::chrono_taker<> tt(cout, "std::vector<float> with Linear Data find values in range", tests);
            for (unsigned int i = 0; i < tests; i++)
            {
                float from = fromVect[i];
                float to   = toVect[i];

                in_range_vect(linData, from, to, bv_range);
                xorVect ^= bv_range;
                bv_range.clear(true);
            }
        }

        {
            bm::chrono_taker<> tt(cout, "SVF with Linear Data find values in range with Const Iterator", tests);
            for (unsigned int i = 0; i < tests; i++)
            {
                float from = fromVect[i];
                float to   = toVect[i];

                in_range_const(testSVF, from, to, bv_range);
                xorConst ^= bv_range;
                bv_range.clear(true);
            }
        }

        {
            bm::chrono_taker<> tt(cout, "SVF with Linear Data find values in range with bm::for_each_sparse", tests);
            for (unsigned int i = 0; i < tests; i++)
            {
                float from = fromVect[i];
                float to   = toVect[i];

                in_range_for_each_sparse(testSVF, from, to, bv_range);
                xorForEachSparse ^= bv_range;
                bv_range.clear(true);
            }
        }

        bool range_eq_vector = (xorSV == xorVect);
        bool range_eq_const  = (xorSV == xorConst);
        bool range_eq_for_each_sparse = (xorSV == xorForEachSparse);

        if (!range_eq_vector || !range_eq_const || !range_eq_for_each_sparse)
        {
            cerr << "Linear: MISMATCH" << endl;
            exit(1);
        }
    }

    testSVF.clear();
    std::vector<float> randData(N);

    for (sparseVecFloat::size_type i = 0; i < N; ++i)
    {
        randData[i] = dis(gen);
    }

    testSVF.import(randData.data(), N);
    testSVF.optimize(tb);

    {
        sparseVecFloat::bvector_type xorSV;
        sparseVecFloat::bvector_type xorVect;
        sparseVecFloat::bvector_type xorConst;
        sparseVecFloat::bvector_type xorForEachSparse;

        sparseVecFloat::bvector_type bv_range;
        sparseVecFloat::bvector_type::mem_pool_guard bv_range_guard(bv_pool, bv_range);

        std::vector<float> fromVect(tests);
        std::vector<float> toVect(tests);
        for (unsigned int i = 0; i < tests; i++)
        {
            fromVect[i] = dis(gen);
            toVect[i] = dis(gen);
        }

        {
            bm::chrono_taker<> tt(cout, "SVF with Random Data find values in range with scanner", tests);
            bm::sparse_vector_scanner<sparseVecFloat> scan;

            for (unsigned int i = 0; i < tests; i++)
            {
                float from = fromVect[i];
                float to   = toVect[i];

                scan.find_range_float(testSVF, from, to, bv_range);
                xorSV ^= bv_range;
                bv_range.clear(true);
            }
        }

        {
            bm::chrono_taker<> tt(cout, "std::vector<float> with Random Data find values in range", tests);
            for (unsigned int i = 0; i < tests; i++)
            {
                float from = fromVect[i];
                float to   = toVect[i];

                in_range_vect(randData, from, to, bv_range);
                xorVect ^= bv_range;
                bv_range.clear(true);
            }
        }

        {
            bm::chrono_taker<> tt(cout, "SVF with Random Data find values in range with Const Iterator", tests);
            for (unsigned int i = 0; i < tests; i++)
            {
                float from = fromVect[i];
                float to   = toVect[i];

                in_range_const(testSVF, from, to, bv_range);
                xorConst ^= bv_range;
                bv_range.clear(true);
            }
        }

        {
            bm::chrono_taker<> tt(cout, "SVF with Random Data find values in range with bm::for_each_sparse", tests);
            for (unsigned int i = 0; i < tests; i++)
            {
                float from = fromVect[i];
                float to   = toVect[i];

                in_range_for_each_sparse(testSVF, from, to, bv_range);
                xorForEachSparse ^= bv_range;
                bv_range.clear(true);
            }
        }

        bool range_eq_vector = (xorSV == xorVect);
        bool range_eq_const  = (xorSV == xorConst);
        bool range_eq_for_each_sparse = (xorSV == xorForEachSparse);

        if (!range_eq_vector || !range_eq_const || !range_eq_for_each_sparse)
        {
            cerr << "Random: MISMATCH" << endl;
            exit(1);
        }
    }

    testSVF.clear();
    std::vector<float> skewData(N);

    for (sparseVecFloat::size_type i = 19000000; i < N; ++i)
    {
        skewData[i] = dis(gen);
    }

    testSVF.import(skewData.data(), N);
    testSVF.optimize(tb);

    {
        sparseVecFloat::bvector_type xorSV;
        sparseVecFloat::bvector_type xorVect;
        sparseVecFloat::bvector_type xorConst;
        sparseVecFloat::bvector_type xorForEachSparse;

        sparseVecFloat::bvector_type bv_range;
        sparseVecFloat::bvector_type::mem_pool_guard bv_range_guard(bv_pool, bv_range);

        std::vector<float> fromVect(tests);
        std::vector<float> toVect(tests);
        for (unsigned int i = 0; i < tests; i++)
        {
            fromVect[i] = dis(gen);
            toVect[i] = dis(gen);
        }

        {
            bm::chrono_taker<> tt(cout, "SVF with Skewed Data find values in range with scanner", tests);
            bm::sparse_vector_scanner<sparseVecFloat> scan;
            for (unsigned int i = 0; i < tests; i++)
            {
                float from = fromVect[i];
                float to   = toVect[i];

                scan.find_range_float(testSVF, from, to, bv_range);
                xorSV ^= bv_range;
                bv_range.clear(true);
            }
        }

        {
            bm::chrono_taker<> tt(cout, "std::vector<float> with Skewed Data find values in range", tests);
            for (unsigned int i = 0; i < tests; i++)
            {
                float from = fromVect[i];
                float to   = toVect[i];

                in_range_vect(skewData, from, to, bv_range);
                xorVect ^= bv_range;
                bv_range.clear(true);
            }
        }

        {
            bm::chrono_taker<> tt(cout, "SVF with Skewed Data find values in range with Const Iterator", tests);
            for (unsigned int i = 0; i < tests; i++)
            {
                float from = fromVect[i];
                float to   = toVect[i];

                in_range_const(testSVF, from, to, bv_range);
                xorConst ^= bv_range;
                bv_range.clear(true);
            }
        }

        {
            bm::chrono_taker<> tt(cout, "SVF with Skewed Data find values in range with bm::for_each_sparse", tests);
            for (unsigned int i = 0; i < tests; i++)
            {
                float from = fromVect[i];
                float to   = toVect[i];

                in_range_for_each_sparse(testSVF, from, to, bv_range);
                xorForEachSparse ^= bv_range;
                bv_range.clear(true);
            }
        }

        bool range_eq_vector = (xorSV == xorVect);
        bool range_eq_const  = (xorSV == xorConst);
        bool range_eq_for_each_sparse = (xorSV == xorForEachSparse);

        if (!range_eq_vector || !range_eq_const || !range_eq_for_each_sparse)
        {
            cerr << "Skewed: MISMATCH" << endl;
            exit(1);
        }
    }
}


//-----------------------------------------------------------------------------------------------


//Finds all values in range [from, to] in a given std::vector<float> and flipts the corresponding bits in bv_out
inline
void in_range_vect_rsc(const std::vector<float>& fv, float from, float to, sparseVecFloatRSC::bvector_type &bv_out)
{
    if (from > to) std::swap(from, to);
    for (sparseVecFloatRSC::size_type i = 0; i < fv.size(); i++)
    {
        if(fv[i] >= from && fv[i] <= to)
            bv_out.set(i);
    }
}

//Finds all values in range [from, to] in a given sparse_vector_float which uses a rsc sparse vector
//using a const_iterator and flips the corresponding bits in bv_out
inline
void in_range_const_rsc(const sparseVecFloatRSC& sv, float from, float to, sparseVecFloatRSC::bvector_type &bv_out)
{
    if (from > to) std::swap(from, to);
    sparseVecFloatRSC::const_iterator ci = sv.begin();
    for (; ci.valid(); ++ci)
    {
        if (auto v = ci.value(); v >= from && v <= to)
            bv_out.set(ci.pos());
    }
}

struct svf_rsc_range_for_each_sparse_func
{
    typedef sparseVecFloatRSC::size_type size_type;

    svf_rsc_range_for_each_sparse_func(float from, float to,
                                       sparseVecFloatRSC::bvector_type& bv_out)
        : from_(from),
          to_(to),
          bv_out_(bv_out)
    {}

    int operator()(float v, bool is_null, size_type idx)
    {
        if (!is_null && v >= from_ && v <= to_)
            bv_out_.set(idx);
        return 0;
    }

    float from_;
    float to_;
    sparseVecFloatRSC::bvector_type& bv_out_;
};

inline
void in_range_for_each_sparse_rsc(const sparseVecFloatRSC& sv, float from, float to,
                                  sparseVecFloatRSC::bvector_type& bv_out)
{
    if (from > to) std::swap(from, to);
    svf_rsc_range_for_each_sparse_func func(from, to, bv_out);
    bm::for_each_sparse(sv, func);
}

// -------------------------------------------------------------------

void TestSVFScannerRSC()
{
    BM_DECLARE_TEMP_BLOCK(tb)

    typedef sparseVecFloatRSC::bvector_type bvect_type;
    typedef bvect_type::allocator_type::allocator_pool_type bvect_pool_type;
    bvect_pool_type bv_rsc_pool;

    sparseVecFloatRSC::size_type N = 20000000;
    std::random_device rd;

    float upper = 1000000.0f;
    float lower = -1000000.0f;
    std::uniform_real_distribution<float> dis(lower, upper);
    std::uniform_real_distribution<float> null_chance(0.0f, 1.0f);

    std::vector<float> linData(N);

    for(sparseVecFloatRSC::size_type i = 0; i < N/2; i++)
    {
        if (null_chance(gen) >= 0.35f)
        {
            linData[i] = -1.0f * (float)i * 0.00123f;
        }
        else
        {
            linData[i] = std::numeric_limits<float>::quiet_NaN();
        }
    }
    for(sparseVecFloatRSC::size_type i = 0; i < N/2; i++)
    {
        if (null_chance(gen) >= 0.35f)
        {
            linData[i+N/2] = (float)i * 0.00123f;
        }
        else
        {
            linData[i] = std::numeric_limits<float>::quiet_NaN();
        }
    }

    sparseVecFloatRSC testSVF;
    testSVF.import(linData.data(), N);
    testSVF.optimize(tb);
    testSVF.sync(true, true);

    unsigned int tests = 1000;

    {
        sparseVecFloatRSC::bvector_type xorRSC;
        sparseVecFloatRSC::bvector_type xorVect;
        sparseVecFloatRSC::bvector_type xorConst;
        sparseVecFloatRSC::bvector_type xorForEachSparse;

        sparseVecFloatRSC::bvector_type bv_range;
        sparseVecFloatRSC::bvector_type::mem_pool_guard bv_range_guard(bv_rsc_pool, bv_range);

        std::vector<float> fromVect(tests);
        std::vector<float> toVect(tests);
        for (unsigned int i = 0; i < tests; i++)
        {
            fromVect[i] = dis(gen);
            toVect[i] = dis(gen);
        }

        {
            bm::chrono_taker<> tt(cout, "SVF RSC with Linear Data find values in range with scanner", tests);
            bm::sparse_vector_scanner<sparseVecFloatRSC> scan;
            for (unsigned int i = 0; i < tests; i++)
            {
                float from = fromVect[i];
                float to   = toVect[i];

                scan.find_range_float(testSVF, from, to, bv_range);
                xorRSC ^= bv_range;
                bv_range.clear(true);
            }
        }

        {
            bm::chrono_taker<> tt(cout, "std::vector<float> with Linear Data find values in range", tests);
            for (unsigned int i = 0; i < tests; i++)
            {
                float from = fromVect[i];
                float to   = toVect[i];

                in_range_vect_rsc(linData, from, to, bv_range);
                xorVect ^= bv_range;
                bv_range.clear(true);
            }
        }

        {
            bm::chrono_taker<> tt(cout, "SVF RSC with Linear Data find values in range with Const Iterator", tests);
            for (unsigned int i = 0; i < tests; i++)
            {
                float from = fromVect[i];
                float to   = toVect[i];

                in_range_const_rsc(testSVF, from, to, bv_range);
                xorConst ^= bv_range;
                bv_range.clear(true);
            }
        }

        {
            bm::chrono_taker<> tt(cout, "SVF RSC with Linear Data find values in range with bm::for_each_sparse", tests);
            for (unsigned int i = 0; i < tests; i++)
            {
                float from = fromVect[i];
                float to   = toVect[i];

                in_range_for_each_sparse_rsc(testSVF, from, to, bv_range);
                xorForEachSparse ^= bv_range;
                bv_range.clear(true);
            }
        }

        bool range_eq_vector = (xorRSC == xorVect);
        bool range_eq_const  = (xorRSC == xorConst);
        bool range_eq_for_each_sparse = (xorRSC == xorForEachSparse);

        if (!range_eq_vector)
        {
            cerr << "LinearRSC: MISMATCH Vect" << endl;
        }
        if (!range_eq_const)
        {
            cerr << "LinearRSC: MISMATCH Const" << endl;
        }
        if (!range_eq_for_each_sparse)
        {
            cerr << "LinearRSC: MISMATCH for_each_sparse" << endl;
        }
        if (!range_eq_vector || !range_eq_const || !range_eq_for_each_sparse)
        {
            exit(1);
        }
    }

    testSVF.clear();
    std::vector<float> randData(N);

    for (sparseVecFloatRSC::size_type i = 0; i < N; ++i)
    {
        if (null_chance(gen) >= 0.35f)
        {
            randData[i] = dis(gen);
        }
        else
        {
            randData[i] = std::numeric_limits<float>::quiet_NaN();
        }
    }

    testSVF.import(randData.data(), N);
    testSVF.optimize(tb);
    testSVF.sync(true, true);

    //Using a completely random dataset

    {
        sparseVecFloatRSC::bvector_type xorRSC;
        sparseVecFloatRSC::bvector_type xorVect;
        sparseVecFloatRSC::bvector_type xorConst;
        sparseVecFloatRSC::bvector_type xorForEachSparse;

        sparseVecFloatRSC::bvector_type bv_range;
        sparseVecFloatRSC::bvector_type::mem_pool_guard bv_range_guard(bv_rsc_pool, bv_range);

        std::vector<float> fromVect(tests);
        std::vector<float> toVect(tests);
        for (unsigned int i = 0; i < tests; i++)
        {
            fromVect[i] = dis(gen);
            toVect[i] = dis(gen);
        }

        {
            bm::chrono_taker<> tt(cout, "SVF RSC with Random Data find values in range with scanner", tests);
            bm::sparse_vector_scanner<sparseVecFloatRSC> scan;
            for (unsigned int i = 0; i < tests; i++)
            {
                float from = fromVect[i];
                float to   = toVect[i];

                scan.find_range_float(testSVF, from, to, bv_range);
                xorRSC ^= bv_range;
                bv_range.clear(true);
            }
        }

        {
            bm::chrono_taker<> tt(cout, "std::vector<float> with Random Data find values in range", tests);
            for (unsigned int i = 0; i < tests; i++)
            {
                float from = fromVect[i];
                float to   = toVect[i];

                in_range_vect_rsc(randData, from, to, bv_range);
                xorVect ^= bv_range;
                bv_range.clear(true);
            }
        }

        {
            bm::chrono_taker<> tt(cout, "SVF RSC with Random Data find values in range with Const Iterator", tests);
            for (unsigned int i = 0; i < tests; i++)
            {
                float from = fromVect[i];
                float to   = toVect[i];

                in_range_const_rsc(testSVF, from, to, bv_range);
                xorConst ^= bv_range;
                bv_range.clear(true);
            }
        }

        {
            bm::chrono_taker<> tt(cout, "SVF RSC with Random Data find values in range with bm::for_each_sparse", tests);
            for (unsigned int i = 0; i < tests; i++)
            {
                float from = fromVect[i];
                float to   = toVect[i];

                in_range_for_each_sparse_rsc(testSVF, from, to, bv_range);
                xorForEachSparse ^= bv_range;
                bv_range.clear(true);
            }
        }

        bool range_eq_vector = (xorRSC == xorVect);
        bool range_eq_const  = (xorRSC == xorConst);
        bool range_eq_for_each_sparse = (xorRSC == xorForEachSparse);

        if (!range_eq_vector)
        {
            cerr << "RandomRSC: MISMATCH Vect" << endl;
        }
        if (!range_eq_const)
        {
            cerr << "RandomRSC: MISMATCH Const" << endl;
        }
        if (!range_eq_for_each_sparse)
        {
            cerr << "RandomRSC: MISMATCH for_each_sparse" << endl;
        }
        if (!range_eq_vector || !range_eq_const || !range_eq_for_each_sparse)
        {
            exit(1);
        }
    }

    testSVF.clear();
    std::vector<float> skewData(N);
    for (sparseVecFloatRSC::size_type i = 0; i < 19000000; ++i)
    {
        skewData[i] = std::numeric_limits<float>::quiet_NaN();
    }
    for (sparseVecFloatRSC::size_type i = 19000000; i < N; ++i)
    {
        if (null_chance(gen) >= 0.35f)
        {
            skewData[i] = dis(gen);
        }
        else
        {
            skewData[i] = std::numeric_limits<float>::quiet_NaN();
        }
    }

    testSVF.import(skewData.data(), N);
    testSVF.optimize(tb);
    testSVF.sync(true, true);

    {
        sparseVecFloatRSC::bvector_type xorRSC;
        sparseVecFloatRSC::bvector_type xorVect;
        sparseVecFloatRSC::bvector_type xorConst;
        sparseVecFloatRSC::bvector_type xorForEachSparse;

        sparseVecFloatRSC::bvector_type bv_range;
        sparseVecFloatRSC::bvector_type::mem_pool_guard bv_range_guard(bv_rsc_pool, bv_range);

        std::vector<float> fromVect(tests);
        std::vector<float> toVect(tests);
        for (unsigned int i = 0; i < tests; i++)
        {
            fromVect[i] = dis(gen);
            toVect[i] = dis(gen);
        }

        {
            bm::chrono_taker<> tt(cout, "SVF RSC with Skewed Data find values in range with scanner", tests);
            bm::sparse_vector_scanner<sparseVecFloatRSC> scan;
            for (unsigned int i = 0; i < tests; i++)
            {
                float from = fromVect[i];
                float to   = toVect[i];

                scan.find_range_float(testSVF, from, to, bv_range);
                xorRSC ^= bv_range;
                bv_range.clear(true);
            }
        }

        {
            bm::chrono_taker<> tt(cout, "std::vector<float> with Skewed Data find values in range", tests);
            for (unsigned int i = 0; i < tests; i++)
            {
                float from = fromVect[i];
                float to   = toVect[i];

                in_range_vect_rsc(skewData, from, to, bv_range);
                xorVect ^= bv_range;
                bv_range.clear(true);
            }
        }

        {
            bm::chrono_taker<> tt(cout, "SVF RSC with Skewed Data find values in range with Const Iterator", tests);
            for (unsigned int i = 0; i < tests; i++)
            {
                float from = fromVect[i];
                float to   = toVect[i];

                in_range_const_rsc(testSVF, from, to, bv_range);
                xorConst ^= bv_range;
                bv_range.clear(true);
            }
        }

        {
            bm::chrono_taker<> tt(cout, "SVF RSC with Skewed Data find values in range with bm::for_each_sparse", tests);
            for (unsigned int i = 0; i < tests; i++)
            {
                float from = fromVect[i];
                float to   = toVect[i];

                in_range_for_each_sparse_rsc(testSVF, from, to, bv_range);
                xorForEachSparse ^= bv_range;
                bv_range.clear(true);
            }
        }

        bool range_eq_vector = (xorRSC == xorVect);
        bool range_eq_const  = (xorRSC == xorConst);
        bool range_eq_for_each_sparse = (xorRSC == xorForEachSparse);

        if (!range_eq_vector)
        {
            cerr << "SkewedRSC: MISMATCH Vect" << endl;
        }
        if (!range_eq_const)
        {
            cerr << "SkewedRSC: MISMATCH Const" << endl;
        }
        if (!range_eq_for_each_sparse)
        {
            cerr << "SkewedRSC: MISMATCH for_each_sparse" << endl;
        }
        if (!range_eq_vector || !range_eq_const || !range_eq_for_each_sparse)
        {
            exit(1);
        }
    }
}

void TestSVFComparison()
{
    BM_DECLARE_TEMP_BLOCK(tb)

    typedef sparseVecFloat::bvector_type bvect_type;
    typedef bvect_type::allocator_type::allocator_pool_type bvect_pool_type;
    bvect_pool_type bv_pool;

    typedef sparseVecFloatRSC::bvector_type bvect_rsc_type;
    typedef bvect_rsc_type::allocator_type::allocator_pool_type bvect_rsc_pool_type;
    bvect_rsc_pool_type bv_rsc_pool;
    
    sparseVecFloat::size_type N = 20000000;
    std::random_device rd;

    float upper = 15000.0f;
    float lower = -15000.0f;
    std::uniform_real_distribution<float> dis(lower, upper);
    std::uniform_real_distribution<float> null_chance(0.0f, 1.0f);
    
    std::vector<float> linData(N);

    for(sparseVecFloatRSC::size_type i = 0; i < N/2; i++)
    {
        if (null_chance(gen) >= 0.35f)
        {
            linData[i] = -1.0f * (float)i * 0.00123f;
        }
        else
        {
            linData[i] = std::numeric_limits<float>::quiet_NaN();
        }
    }
    for(sparseVecFloatRSC::size_type i = 0; i < N/2; i++)
    {
        if (null_chance(gen) >= 0.35f)
        {
            linData[i+N/2] = (float)i * 0.00123f;
        }
        else
        {
            linData[i] = std::numeric_limits<float>::quiet_NaN();
        }
    }
    unsigned int tests = 1000;
    std::vector<float> fromVect(tests);
    std::vector<float> toVect(tests);
    for (unsigned int i = 0; i < tests; i++)
    {
        fromVect[i] = dis(gen);
        toVect[i] = dis(gen);
    }
    
    sparseVecFloat svf(bm::use_null);
    svf.import(linData.data(), N);
    {
        sparseVecFloat::bvector_type xorSVF;
        sparseVecFloat::bvector_type xorConst;
        sparseVecFloat::bvector_type bv_range;
        sparseVecFloat::bvector_type::mem_pool_guard bv_range_guard(bv_pool, bv_range);
        
        {
            bm::chrono_taker<> tt(cout, "Unoptimized SVF with random data in range with Scanner", tests);
            bm::sparse_vector_scanner<sparseVecFloat> scan;
            for (unsigned int i = 0; i < tests; i++)
            {
                float from = fromVect[i];
                float to   = toVect[i];
                
                scan.find_range_float(svf, from, to, bv_range);
                xorSVF ^= bv_range;
                bv_range.clear(true);
            }
        }
        
        {
            bm::chrono_taker<> tt(cout, "Unoptimized SVF with random data in range with Const Iterator", tests);
            for (unsigned int i = 0; i < tests; i++)
            {
                float from = fromVect[i];
                float to   = toVect[i];
                
                in_range_const(svf, from, to, bv_range);
                xorConst ^= bv_range;
                bv_range.clear(true);
            }
        }
        
        if(xorSVF != xorConst){
            cerr << "SVF Non-optimized Scanner and Const Iterator do not match" << endl;
            exit(1);
        }
    }
    
    svf.optimize(tb);
    
    {
        sparseVecFloat::bvector_type xorSVF;
        sparseVecFloat::bvector_type xorConst;
        sparseVecFloat::bvector_type bv_range;
        sparseVecFloat::bvector_type::mem_pool_guard bv_range_guard(bv_pool, bv_range);
        
        {
            bm::chrono_taker<> tt(cout, "Optimized SVF with random data in range with Scanner", tests);
            bm::sparse_vector_scanner<sparseVecFloat> scan;
            for (unsigned int i = 0; i < tests; i++)
            {
                float from = fromVect[i];
                float to   = toVect[i];
                
                scan.find_range_float(svf, from, to, bv_range);
                xorSVF ^= bv_range;
                bv_range.clear(true);
            }
        }
        
        {
            bm::chrono_taker<> tt(cout, "Optimized SVF with random data in range with Const Iterator", tests);
            for (unsigned int i = 0; i < tests; i++)
            {
                float from = fromVect[i];
                float to   = toVect[i];
                
                in_range_const(svf, from, to, bv_range);
                xorConst ^= bv_range;
                bv_range.clear(true);
            }
        }
        
        if(xorSVF != xorConst){
            cerr << "SVF Optimized Scanner and Const Iterator do not match" << endl;
            exit(1);
        }
    }
    
    svf.freeze();
    {
        sparseVecFloat::bvector_type xorSVF;
        sparseVecFloat::bvector_type xorConst;
        sparseVecFloat::bvector_type bv_range;
        sparseVecFloat::bvector_type::mem_pool_guard bv_range_guard(bv_pool, bv_range);
        {
            bm::chrono_taker<> tt(cout, "Optimized and Frozen SVF with random data in range with Scanner", tests);
            bm::sparse_vector_scanner<sparseVecFloat> scan;
            for (unsigned int i = 0; i < tests; i++)
            {
                float from = fromVect[i];
                float to   = toVect[i];
                
                scan.find_range_float(svf, from, to, bv_range);
                xorSVF ^= bv_range;
                bv_range.clear(true);
            }
        }
        
        {
            bm::chrono_taker<> tt(cout, "Optimized and Frozen SVF with random data in range with Const Iterator", tests);
            for (unsigned int i = 0; i < tests; i++)
            {
                float from = fromVect[i];
                float to   = toVect[i];
                
                in_range_const(svf, from, to, bv_range);
                xorConst ^= bv_range;
                bv_range.clear(true);
            }
        }
        
        if(xorSVF != xorConst){
            cerr << "SVF Optimized and Frozen Scanner and Const Iterator do not match" << endl;
            exit(1);
        }
    }
    svf.clear();
    
    sparseVecFloatRSC rscSVF;
    rscSVF.import(linData.data(), N);
    rscSVF.sync(true, true);
    {
        sparseVecFloatRSC::bvector_type xorRSC;
        sparseVecFloatRSC::bvector_type xorConst;
        sparseVecFloatRSC::bvector_type bv_range;
        sparseVecFloatRSC::bvector_type::mem_pool_guard bv_range_guard(bv_rsc_pool, bv_range);
        
        {
            bm::chrono_taker<> tt(cout, "Unoptimized RSC SVF with random data in range with Scanner", tests);
            bm::sparse_vector_scanner<sparseVecFloatRSC> scan;
            for (unsigned int i = 0; i < tests; i++)
            {
                float from = fromVect[i];
                float to   = toVect[i];
                
                scan.find_range_float(rscSVF, from, to, bv_range);
                xorRSC ^= bv_range;
                bv_range.clear(true);
            }
        }
        
        {
            bm::chrono_taker<> tt(cout, "Unoptimized RSC SVF with random data in range with Const Iterator", tests);
            for (unsigned int i = 0; i < tests; i++)
            {
                float from = fromVect[i];
                float to   = toVect[i];
                
                in_range_const_rsc(rscSVF, from, to, bv_range);
                xorConst ^= bv_range;
                bv_range.clear(true);
            }
        }
        
        if(xorRSC != xorConst){
            cerr << "SVF RSC Non-Optimized Scanner and Const Iterator do not match" << endl;
            exit(1);
        }
    }
    
    rscSVF.optimize();
    rscSVF.sync(true, true);
    
    {
        sparseVecFloatRSC::bvector_type xorRSC;
        sparseVecFloatRSC::bvector_type xorConst;
        sparseVecFloatRSC::bvector_type bv_range;
        sparseVecFloatRSC::bvector_type::mem_pool_guard bv_range_guard(bv_rsc_pool, bv_range);
        
        {
            bm::chrono_taker<> tt(cout, "Optimized RSC SVF with random data in range with Scanner", tests);
            bm::sparse_vector_scanner<sparseVecFloatRSC> scan;
            for (unsigned int i = 0; i < tests; i++)
            {
                float from = fromVect[i];
                float to   = toVect[i];
                
                scan.find_range_float(rscSVF, from, to, bv_range);
                xorRSC ^= bv_range;
                bv_range.clear(true);
            }
        }
        
        {
            bm::chrono_taker<> tt(cout, "Optimized RSC SVF with random data in range with Const Iterator", tests);
            for (unsigned int i = 0; i < tests; i++)
            {
                float from = fromVect[i];
                float to   = toVect[i];
                
                in_range_const_rsc(rscSVF, from, to, bv_range);
                xorConst ^= bv_range;
                bv_range.clear(true);
            }
        }
        
        if(xorRSC != xorConst){
            cerr << "SVF RSC Optimized Scanner and Const Iterator do not match" << endl;
            exit(1);
        }
    }
    
    rscSVF.freeze();
    rscSVF.sync(true, true);
    
    {
        sparseVecFloatRSC::bvector_type xorRSC;
        sparseVecFloatRSC::bvector_type xorConst;
        sparseVecFloatRSC::bvector_type bv_range;
        sparseVecFloatRSC::bvector_type::mem_pool_guard bv_range_guard(bv_rsc_pool, bv_range);
        
        {
            bm::chrono_taker<> tt(cout, "Optimized and Frozen RSC SVF with random data in range with Scanner", tests);
            bm::sparse_vector_scanner<sparseVecFloatRSC> scan;
            for (unsigned int i = 0; i < tests; i++)
            {
                float from = fromVect[i];
                float to   = toVect[i];
                
                scan.find_range_float(rscSVF, from, to, bv_range);
                xorRSC ^= bv_range;
                bv_range.clear(true);
            }
        }
        
        {
            bm::chrono_taker<> tt(cout, "Optimized and Frozen RSC SVF with random data in range with Const Iterator", tests);
            for (unsigned int i = 0; i < tests; i++)
            {
                float from = fromVect[i];
                float to   = toVect[i];
                
                in_range_const_rsc(rscSVF, from, to, bv_range);
                xorConst ^= bv_range;
                bv_range.clear(true);
            }
        }
        
        if(xorRSC != xorConst){
            cerr << "SVF RSC Optimized and Frozen Scanner and Const Iterator do not match" << endl;
            exit(1);
        }
    }
}

void TestSVFScannerSpike()
{
    unsigned int N = 200000000;
    
    sparseVecFloat testSVF;
    
    float upper = 15.0f;
    float lower = 5.0f;
    std::uniform_real_distribution<float> flatDis(lower, upper);
    
    unsigned int u = 500000;
    unsigned int l = 100000;
    std::uniform_int_distribution<unsigned int> distDis(l, u);
    unsigned int distance = distDis(gen);
    unsigned int spikeDist = 1000;
    
    upper=1050.0f;
    lower=950.0f;
    std::uniform_real_distribution<float> spikeDis(lower, upper);
    
    sparseVecFloat::bvector_type correct;
    std::vector<float> testVect;
    
    for (sparseVecFloat::size_type i = 0; i < N; ++i)
    {
        if (spikeDist > 0) {
            float toAdd = spikeDis(gen);
            testSVF.push_back(toAdd);
            testVect.push_back(toAdd);
            
            spikeDist--;
            
            correct.set(i);
            
            if (spikeDist == 0) {
                distance = distDis(gen);
                testSVF.optimize();
            }
        }
        else if (distance > 0) {
            float toAdd = flatDis(gen);
            testSVF.push_back(toAdd);
            testVect.push_back(toAdd);
            distance--;
            
            if (distance == 0) {
                spikeDist = 1000;
            }
        }
    }
    testSVF.optimize();
    
    typedef sparseVecFloat::bvector_type bvect_type;
    typedef bvect_type::allocator_type::allocator_pool_type bvect_pool_type;
    bvect_pool_type bv_pool;

    sparseVecFloat::bvector_type scanResult;
    sparseVecFloat::bvector_type::mem_pool_guard scan_result_guard(bv_pool, scanResult);
    bm::sparse_vector_scanner<sparseVecFloat> scan;
    
    unsigned int numTests = 1000;
    
    {
        bm::chrono_taker<> tt(cout, "Scanner time to find spikes", numTests);
        for(unsigned int i = 0; i < numTests; i++)
        {
            scan.find_gt_float(testSVF, 900.0f, scanResult);
        }
    }
    
    if(scanResult != correct)
    {
        std::cerr << "Incorrect vector" << std::endl;
        exit(0);
    }
}

int main(void)
{
    bvector64_VerySparse_SetDestroyCycle();
    bvector64_VerySparse_RefAccessCycle();
    bvector64_Serialization();
    
    TestSVFScanner();
    cout << endl;

    TestSVFScannerRSC();
    cout << endl;
    
    TestSVFComparison();
    cout << endl;
    
    TestSVFScannerSpike();
    cout << endl;
    
    std::cout << std::endl << "Performance:" << std::endl;
    bm::chrono_taker<>::print_duration_map(std::cout, timing_map, bm::chrono_taker<>::ct_all);

    return 0;
}


#ifdef _MSC_VER
#pragma warning( pop )
#endif



