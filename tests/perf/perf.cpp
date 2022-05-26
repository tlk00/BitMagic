/*
Copyright(c) 2002-2017 Anatoliy Kuznetsov(anatoliy_kuznetsov at yahoo.com)

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

//#define BMWASMSIMDOPT
//#define BMSSE2OPT
//#define BMSSE42OPT
//#define BMAVX2OPT
//#define BM_USE_GCC_BUILD

#define BM_NONSTANDARD_EXTENTIONS
//#define BM64ADDR
#ifdef _MSC_VER
#pragma warning( push )
#pragma warning( disable : 4996)
#endif

#include <vector>
#include <random>
#include <memory>
#include <algorithm>

#include "bm.h"
#include "bmalgo.h"
#include "bmintervals.h"
#include "bmaggregator.h"
#include "bmserial.h"
#include "bmsparsevec.h"
#include "bmsparsevec_algo.h"
#include "bmsparsevec_serial.h"
#include "bmstrsparsevec.h"
#include "bmsparsevec_compr.h"
#include "bmrandom.h"

#include "bmtimer.h"

#include "bmdbg.h"

#include <math.h>

using namespace std;

const unsigned int BSIZE = 130000000;
const unsigned int REPEATS = 1200;

typedef  bitset<BSIZE>  test_bitset;

unsigned platform_test = 1;

std::random_device rand_dev;
std::mt19937 gen(rand_dev()); // mersenne_twister_engine 
std::uniform_int_distribution<> rand_dis(0, BSIZE); // generate uniform numebrs for [1, vector_max]

typedef bm::bvector<> bvect;
typedef bm::str_sparse_vector<char, bvect, 8> str_svect_type;


// generate pseudo-random bit-vector, mix of compressed/non-compressed blocks
//
static
void generate_bvector(bvect& bv, unsigned vector_max = 40000000, bool optimize = true)
{
    unsigned i, j;
    for (i = 0; i < vector_max;)
    {
        // generate bit-blocks
        {
            bvect::bulk_insert_iterator iit(bv);
            for (j = 0; j < 65535 * 10; i += 10, j++)
            {
                iit = i;
            }
        }
        if (i > vector_max)
            break;
        // generate GAP (compressed) blocks
        for (j = 0; j < 65535; i += 120, j++)
        {
            unsigned len = (unsigned)rand() % 64;
            bv.set_range(i, i + len);
            i += len;
            if (i > vector_max)
                break;
        }
    }
    if (optimize)
        bv.optimize();
}


static
void SimpleFillSets(test_bitset* bset, 
                       bvect& bv,
                       unsigned min, 
                       unsigned max,
                       unsigned fill_factor,
                       bool set_flag=true)
{
    bv.init();
    for (unsigned i = min; i < max; i+=fill_factor)
    {
        if (bset)
            (*bset)[i] = set_flag;
        if (set_flag)
            bv.set_bit_no_check(i);
        else
            bv.clear_bit_no_check(i);
    } // for i
}


//
// Interval filling.
// 111........111111........111111..........11111111.......1111111...
//
static
void FillSetsIntervals(test_bitset* bset,
                       bvect& bv,
                       unsigned min = 0, 
                       unsigned max = BSIZE,
                       unsigned fill_factor = 10,
                       bool set_flag=true)
{
    while(fill_factor==0)
    {
        fill_factor=(unsigned)rand()%10;
    }

    unsigned i, j;
    unsigned factor = 10 * fill_factor;
    for (i = min; i < max; ++i)
    {
        unsigned len, end; 

        do
        {
            len = unsigned(rand()) % factor;
            end = i+len;
            
        } while (end >= max);
        for (j = i; j < end; ++j)
        {
            if (set_flag)
            {
                if (bset)
                    bset->set(j, true);
                bv[j]= true;
            }
            else
            {
                if (bset)
                    bset->set(j, false);
                bv[j] = false;
            }
                           
        } // j
        i = end;
        len = (unsigned)rand() % 10;

        i+=len;
        {
            for(unsigned k=0; k < 1000 && i < max; k+=3,i+=3)
            {
                if (set_flag)
                {
                    if (bset)
                        bset->set(i, true);
                    bv[i] = true;
                }
                else
                {
                    if (bset)
                        bset->set(j, false);
                    bv[j] = false;
                }
            }
        }

    } // for i

}

static
void generate_sparse_bvector(bvect& bv,
                             unsigned min = 0,
                             unsigned max = BSIZE,
                             unsigned fill_factor = 65536)
{
    bvect::bulk_insert_iterator iit(bv);
    unsigned ff = fill_factor / 10;
    for (unsigned i = min; i < max; i+= ff)
    {
        //bv.set(i);
        iit = i;
        ff += ff / 2;
        if (ff > fill_factor)
            ff = fill_factor / 10;
    }
    iit.flush();
}


static
void GenerateTestCollection(std::vector<bvect>* target,
                            unsigned count = 30,
                            unsigned vector_max = 40000000,
                            bool optimize = true)
{
    assert(target);
    bvect bv_common; // sub-vector common for all collection
    generate_sparse_bvector(bv_common, vector_max/10, vector_max, 250000);
    
    unsigned cnt1 = (count / 2);
    
    unsigned i = 0;
    
    for (i = 0; i < cnt1; ++i)
    {
        std::unique_ptr<bvect> bv (new bvect);
        generate_bvector(*bv, vector_max, optimize);
        *bv |= bv_common;
        if (optimize)
            bv->optimize();
        target->push_back(std::move(*bv));
    } // for
    
    unsigned fill_factor = 10;
    for (; i < count; ++i)
    {
        std::unique_ptr<bvect> bv (new bvect);
        
        FillSetsIntervals(0, *bv, vector_max/ 10, vector_max, fill_factor);
        *bv |= bv_common;

        target->push_back(std::move(*bv));
    } // for
}

static
void GenerateTestStrCollection(std::vector<string>& str_coll, unsigned max_coll)
{
    string prefix = "az";
    string str;
    for (unsigned i = 0; i < max_coll; ++i)
    {
        str = prefix;
        str.append(to_string(i));
        str_coll.emplace_back(str);

        if (i % 1024 == 0) // generate new prefix
        {
            prefix.clear();
            unsigned prefix_len = (unsigned)rand() % 5;
            for (unsigned j = 0; j < prefix_len; ++j)
            {
                char cch = char('a' + (unsigned)rand() % 26);
                prefix.push_back(cch);
            } // for j
        }
    } // for i
}


static
void MemCpyTest()
{
    unsigned* m1 = new unsigned[BSIZE/32];
    unsigned* m2 = new unsigned[BSIZE/32];
    
    unsigned int i,j;

    if (!platform_test)
    {
    bm::chrono_taker<std::ostream> tt(cout, "Memory ADD transfer test", REPEATS * 4);
    for (i = 0; i < REPEATS*4; ++i)
    {
        for (j = 0; j < BSIZE/32; j+=4)
        {
            m1[j+0] += m2[j+0];
            m1[j+1] += m2[j+1];
            m1[j+2] += m2[j+2];
            m1[j+3] += m2[j+3];
        }
    }
    }
    
    if (!platform_test)
    {
    bm::chrono_taker<std::ostream> tt(cout, "memcpy transfer test", REPEATS * 4);
    for (i = 0; i < REPEATS*4; ++i)
    {
        memcpy(m1, m2, BSIZE/32 * sizeof(unsigned));
    }
    }
    
    delete [] m1;
    delete [] m2;
}

static
void BitCountTest()
{
    {
    bvect*  bv = new bvect();
    test_bitset*  bset = new test_bitset();
    unsigned value = 0;

    FillSetsIntervals(bset, *bv, 0, BSIZE, 10);
    {
        auto c0 = bset->count();
        auto c1 = bv->count();
        assert(c0 == c1);
        if (c0 != c1)
        {
            cerr << "Fill sets integrity failed!" << endl;
            exit(1);
        }
    }

    //if (!platform_test)
    {
        bm::chrono_taker<std::ostream> tt(cout, "bvector<>::count() ", REPEATS*10);
        for (unsigned i = 0; i < REPEATS*10; ++i)
        {
            value+=bv->count();
        }
    }


    volatile unsigned* p = &value;
    unsigned c1;
    c1 = value = 0;

    if (!platform_test)
    {
    bm::chrono_taker<std::ostream> tt(cout, "BitCount. Random bitvector (STL)", REPEATS*2);
    for (unsigned i = 0; i < REPEATS*2; ++i)
    {    
        value += (unsigned)bset->count();
    }
    }

    c1 = *p;
    c1 = value = 0;
    stringstream s;
    s << value << c1; // to fool the optimization

    bv->freeze();

    {
        bm::chrono_taker<std::ostream> tt(cout, "bvector<>::count() (RO)", REPEATS*10);
        for (unsigned i = 0; i < REPEATS*10; ++i)
        {
            value+=bv->count();
        }
    }


    delete bset;
    delete bv;

    }
}

static
void CheckBitList(const unsigned* bl1, unsigned bl1_cnt,
                  const unsigned* bl2, unsigned bl2_cnt)
{
    assert(bl1_cnt == bl2_cnt);
    if (bl1_cnt != bl2_cnt)
    {
        cerr << "Check list count failed!" << endl;
        exit(1);
    }
    for (unsigned i = 0; i < bl1_cnt; ++i)
    {
        assert(bl1[i] == bl2[i]);
        if (bl1[i] != bl2[i])
        {
            cerr << "BitList check failed!" << endl;
            exit(1);
        }
    }
}

static
void BitForEachTest()
{
    // setup the test data
    //
    unsigned value_to = 65536 * 64;
    size_t value = 0;
    unsigned* test_arr = new unsigned[value_to];
    {
        unsigned j=0;
        for (; j < value_to/3; ++j) // sparse
        {
            test_arr[j] = (1u << 25) | (1u << 7) | 2u | 4u;
            test_arr[++j] = 0;
        }
        for (; j < 2 * (value_to/3) ; ++j) // medium
            test_arr[j] = j;

        for (; j < value_to ; ++j) // dense
            test_arr[j] = ~0u & ~((7u << 25) | (19u << 7) | 2u | 4u);
    }

    // quick test
    {
        unsigned bit_list0[32];
        unsigned bit_list1[32];
        unsigned bit_list2[32];
        for (unsigned j = 0; j < value_to; ++j)
        {
            auto c0 = bm::bitscan_nibble(test_arr[j], bit_list0);
            auto c1 = bm::bitscan_popcnt(test_arr[j], bit_list1);
            auto c2 = bm::bitscan_bsf(test_arr[j], bit_list2);
            CheckBitList(bit_list0, c0, bit_list1, c1);
            CheckBitList(bit_list0, c0, bit_list2, c2);
        }
    }

    bm::id64_t sum1(0);
    if (platform_test)
    {
        unsigned bit_list[32];
        bm::chrono_taker<std::ostream> tt(cout, "BitScan-nibble (switch based)", REPEATS);

        for (unsigned i = 0; i < REPEATS; ++i)
        {
            for (unsigned j = 0; j < value_to; ++j)
            {
                sum1 += bm::bitscan_nibble(test_arr[j], bit_list);
            }
        }
    }
    #ifdef BM_NONSTANDARD_EXTENTIONS
    #ifdef __GNUC__
    bm::id64_t sum2(0);
    if (platform_test)
    {
        unsigned bit_list[32];
        bm::chrono_taker<std::ostream> tt(cout, "BitScan-nibble (GCC goto)", REPEATS);

        for (unsigned i = 0; i < REPEATS; ++i)
        {
            for (unsigned j = 0; j < value_to; ++j)
            {
                sum2 += bm::bitscan_nibble_gcc(test_arr[j], bit_list);
            }
        }
    }
    assert(sum1 == sum2);
    #endif
    #endif

    bm::id64_t sum3(0);
    {
        unsigned bit_list[32];
        bm::chrono_taker<std::ostream> tt(cout, "BitScan-POPCNT ", REPEATS);
        for (unsigned i = 0; i < REPEATS; ++i)
        {
            for (unsigned j = 0; j < value_to; j+=2)
            {
                sum3 += bm::bitscan_popcnt(test_arr[j], bit_list);
                sum3 += bm::bitscan_popcnt(test_arr[j+1], bit_list);
            }
        }
    }
    assert(sum1 == sum3);
    sum3 = 0;
    {
        unsigned bit_list[64];
        bm::chrono_taker<std::ostream> tt(cout, "BitScan-POPCNT-64 ", REPEATS);
        for (unsigned i = 0; i < REPEATS; ++i)
        {
            for (unsigned j = 0; j < value_to; j+=2)
            {
                bm::id64_t w0 = test_arr[j];
                bm::id64_t w1 = test_arr[j+1];
                bm::id64_t w = w0 | (w1 << 32);
                sum3 += bm::bitscan_popcnt64(w, bit_list);
            }
        }
    }
    assert(sum1 == sum3);


    bm::id64_t sum4(0);
    {
        unsigned bit_list[32];
        bm::chrono_taker<std::ostream> tt(cout, "BitScan-BSF ", REPEATS);
        for (unsigned i = 0; i < REPEATS; ++i)
        {
            for (unsigned j = 0; j < value_to; j+=2)
            {
                sum4 += bm::bitscan_bsf(test_arr[j], bit_list);
                sum4 += bm::bitscan_bsf(test_arr[j+1], bit_list);
            }
        }
    }
    assert(sum1 == sum4);
    sum4 = 0;

    {
        unsigned bit_list[32];
        bm::chrono_taker<std::ostream> tt(cout, "BitScan-BSF-64 ", REPEATS);
        for (unsigned i = 0; i < REPEATS; ++i)
        {
            for (unsigned j = 0; j < value_to; j+=2)
            {
                bm::id64_t w0 = test_arr[j];
                bm::id64_t w1 = test_arr[j+1];
                bm::id64_t w = w0 | (w1 << 32);
                sum4 += bm::bitscan_bsf64(w, bit_list);
            }
        }
    }
    assert(sum1 == sum4);

    char buf[256];
    sprintf(buf, "%i", (int)value); // to fool some smart compilers like ICC


    delete [] test_arr;
}

static
void WordSelectTest()
{
    const unsigned test_size = 50000000;

    std::vector<bm::id64_t> vect_v(test_size);
    std::vector<unsigned> vect_cnt(test_size);
    std::vector<std::pair<unsigned, unsigned> > vect_cnt32(test_size);
    std::vector<unsigned> vect_r1(test_size);
    std::vector<unsigned> vect_r2(test_size);
    std::vector<std::pair<unsigned, unsigned> > vect_r2_32(test_size);

    for (unsigned i = 0; i < test_size; ++i)
    {
        unsigned w;
        bm::id64_t w64;
        do
        {
            w = unsigned(rand());
            if (i % 3 == 0)
                w64 = w << rand()%2048;
            else
            if (i % 4 == 0)
                w64 = bm::id64_t(w) | (bm::id64_t(w) << 32);
            else
                w64 = bm::id64_t(w) | (0xFFAull << 33ull);
        } while (!w64);
        
        if (i % 3 == 0)
            w64 = w << rand()%32;
        else
        if (i % 4 == 0)
            w64 = w | (bm::id64_t(w) << 32);
        else
            w64 = w | (1ull << 60);
        
        unsigned bc = bm::word_bitcount64(w64);
        
        vect_v[i] = w64;
        vect_cnt[i] = bc;
        unsigned w0 = unsigned(w64);
        unsigned w1 = unsigned(w64 >> 32);
        vect_cnt32[i].first = bm::word_bitcount64(w0);
        vect_cnt32[i].second = bm::word_bitcount64(w1);
    }
    
    {
        bm::chrono_taker<std::ostream> tt(cout, "select64 linear", 1);
        for (unsigned i = 0; i < vect_v.size(); ++i)
        {
            bm::id64_t w64 = vect_v[i];
            unsigned bc = vect_cnt[i];
            if (bc)
            {
                for (unsigned j = 1; j <= bc; ++j)
                {
                    unsigned idx = bm::word_select64_linear(w64, j);
                    vect_r1[i] = idx;
                }
            }
            else
            {
                vect_r1[i] = 0;
            }
        }
    }

    {
        bm::chrono_taker<std::ostream> tt(cout, "select64 bitscan_popcnt", 1);
        for (unsigned i = 0; i < vect_v.size(); ++i)
        {
            bm::id64_t w64 = vect_v[i];
            unsigned bc = vect_cnt[i];
            if (bc)
            {
                for (unsigned j = 1; j <= bc; ++j)
                {
                    unsigned idx = bm::word_select64_bitscan_popcnt(w64, j);
                    vect_r2[i] = idx;
                }
            }
            else
            {
                vect_r2[i] = 0;
            }
        }
    }
    {
        bm::chrono_taker<std::ostream> tt(cout, "select64 bitscan_tz", 1);
        for (unsigned i = 0; i < vect_v.size(); ++i)
        {
            bm::id64_t w64 = vect_v[i];
            unsigned bc = vect_cnt[i];
            if (bc)
            {
                for (unsigned j = 1; j <= bc; ++j)
                {
                    unsigned idx = bm::word_select64_bitscan_tz(w64, j);
                    vect_r2[i] = idx;
                }
            }
            else
            {
                vect_r2[i] = 0;
            }
        }
    }

    {
        bm::chrono_taker<std::ostream> tt(cout, "select32 bitscan_popcnt", 1);
        for (unsigned i = 0; i < vect_v.size(); ++i)
        {
            bm::id64_t w64 = vect_v[i];
            unsigned w0 = unsigned(w64);
            unsigned w1 = unsigned(w64 >> 32);
            unsigned bc0 = vect_cnt32[i].first;
            unsigned bc1 = vect_cnt32[i].second;
            if (bc0)
            {
                for (unsigned j = 1; j <= bc0; ++j)
                {
                    unsigned idx = bm::word_select32_bitscan_popcnt(w0, j);
                    vect_r2_32[i].first = idx;
                }
            }
            else
            {
                vect_r2_32[i].first = 0;
            }
            if (bc1)
            {
                for (unsigned j = 1; j <= bc1; ++j)
                {
                    unsigned idx = bm::word_select32_bitscan_popcnt(w1, j);
                    vect_r2_32[i].second = idx;
                }
            }
            else
            {
                vect_r2_32[i].second = 0;
            }
        }
    }

    {
        bm::chrono_taker<std::ostream> tt(cout, "select32 bitscan_tz", 1);
        for (unsigned i = 0; i < vect_v.size(); ++i)
        {
            bm::id64_t w64 = vect_v[i];
            unsigned w0 = unsigned(w64);
            unsigned w1 = unsigned(w64 >> 32);
            unsigned bc0 = vect_cnt32[i].first;
            unsigned bc1 = vect_cnt32[i].second;
            if (bc0)
            {
                for (unsigned j = 1; j <= bc0; ++j)
                {
                    unsigned idx = bm::word_select32_bitscan_tz(w0, j);
                    vect_r2_32[i].first = idx;
                }
            }
            else
            {
                vect_r2_32[i].first = 0;
            }
            if (bc1)
            {
                for (unsigned j = 1; j <= bc1; ++j)
                {
                    unsigned idx = bm::word_select32_bitscan_tz(w1, j);
                    vect_r2_32[i].second = idx;
                }
            }
            else
            {
                vect_r2_32[i].second = 0;
            }
        }
    }


#ifdef BMBMI1OPT
    std::vector<unsigned> vect_r3(test_size);
    std::vector<unsigned> vect_r4(test_size);

    {
        bm::chrono_taker tt(cout, "select64 BMI1 lead-zero", 1);
        for (unsigned i = 0; i < vect_v.size(); ++i)
        {
            bm::id64_t w64 = vect_v[i];
            unsigned bc = vect_cnt[i];
            if (bc)
            {
                for (unsigned j = 1; j <= bc; ++j)
                {
                    unsigned idx = bm::bmi1_select64_lz(w64, j);
                    vect_r3[i] = idx;
                }
            }
            else
            {
                vect_r3[i] = 0;
            }
        }
    }
    
    {
        bm::chrono_taker<std::ostream> tt(cout, "select64 BMI1 bitscan", 1);
        for (unsigned i = 0; i < vect_v.size(); ++i)
        {
            bm::id64_t w64 = vect_v[i];
            unsigned bc = vect_cnt[i];
            if (bc)
            {
                for (unsigned j = 1; j <= bc; ++j)
                {
                    unsigned idx = bm::bmi1_select64_tz(w64, j);
                    vect_r4[i] = idx;
                }
            }
            else
            {
                vect_r4[i] = 0;
            }
        }
    }
#endif

#ifdef BMBMI2OPT
    std::vector<unsigned> vect_r5(test_size);

    {
        bm::chrono_taker<std::ostream> tt(cout, "select64 BMI2 pdep", 1);
        for (unsigned i = 0; i < vect_v.size(); ++i)
        {
            bm::id64_t w64 = vect_v[i];
            unsigned bc = vect_cnt[i];
            if (bc)
            {
                for (unsigned j = 1; j <= bc; ++j)
                {
                    unsigned idx = bm::bmi2_select64_pdep(w64, j);
                    vect_r5[i] = idx;
                }
            }
            else
            {
                vect_r5[i] = 0;
            }
        }
    }
    
#endif


    // validation
    //
    for (unsigned i = 0; i < vect_v.size(); ++i)
    {
        auto r1 = vect_r1[i];
        auto r2 = vect_r2[i];
        
        if (r1 != r2)
        {
            std::cerr << "WordSelect64 error(1) at: " << i << std::endl;
            exit(1);
        }
#ifdef BMBMI1OPT
        auto r3 = vect_r3[i];
        auto r4 = vect_r4[i];

        if (r1 != r3)
        {
            std::cerr << "WordSelect64 BMI1 error(3) at: " << i << std::endl;
            exit(1);
        }
        if (r1 != r4)
        {
            std::cerr << "WordSelect64 BMI1 error(4) at: " << i << std::endl;
            exit(1);
        }
#endif

#ifdef BMBMI2OPT
        auto r5 = vect_r5[i];

        if (r1 != r5)
        {
            std::cerr << "WordSelect64 BMI2 error(5) at: " << i << std::endl;
            exit(1);
        }
#endif

    }

}

static
void BitCountSparseTest()
{
    bvect*  bv = new bvect();
    test_bitset*  bset = new test_bitset();
    size_t value = 0, c1, value2(0);
    size_t value_t = 0, value2_t = 0;
    volatile size_t* p = &value;

    SimpleFillSets(bset, *bv, 0, BSIZE, 130);
    bv->set_range(65536, 65536*3); // create a dense area of FULL blocks
/*
    {
        bm::chrono_taker<std::ostream> tt(cout, "BitCount: Sparse bitset ", REPEATS*10);
        for (unsigned i = 0; i < REPEATS*10; ++i)
        {    
            value += bv->count();
        }
    }
*/
    if (!platform_test)
    {
        bm::chrono_taker<std::ostream> tt(cout, "BitCount: Sparse bitset (STL)", REPEATS*10);
        for (unsigned int i = 0; i < REPEATS*10; ++i)
        {    
            value += bset->count();
        }
    }

    c1 = *p;
    value = c1 = 0;

    bvect*  bv_c = new bvect(*bv);


    std::vector<bvect::size_type> sample_vec;
    sample_vec.reserve(BSIZE*2);
    {
        auto en = bv->get_enumerator(0);
        auto prev = *en;
        for (; en.valid(); ++en)
        {
            auto v = *en;
            sample_vec.push_back(v);
            if (prev+1 != v)
                sample_vec.push_back(prev);
            prev = v;
        }
        std::random_device rd;
        std::mt19937 g(rd());
        std::shuffle(sample_vec.begin(), sample_vec.end(), g);
    }


    BM_DECLARE_TEMP_BLOCK(tb)
    bv->optimize(tb);
    {
    bvect::statistics st;
    bv->calc_stat(&st);
    assert(st.gap_blocks);
    assert(!st.bit_blocks);
    }

    {
    bvect::statistics st_c;
    bv_c->calc_stat(&st_c);
    assert(st_c.gap_blocks == 0);
    assert(st_c.bit_blocks > 0);
    }

/*
    {
        bm::chrono_taker<std::ostream> tt(cout, "BitCount: (GAP) Sparse bitset", REPEATS*100);
        for (unsigned i = 0; i < REPEATS*100; ++i)
        {    
            value3 += bv->count();
        }
    }
*/

    std::unique_ptr<bvect::rs_index_type> bc_arr(new bvect::rs_index_type());
    bv->build_rs_index(bc_arr.get());

    std::unique_ptr<bvect::rs_index_type> bc_arr_c(new bvect::rs_index_type());
    bv_c->build_rs_index(bc_arr_c.get());



    unsigned rs_max = 100;


    {
        bm::chrono_taker<std::ostream> tt(cout, "count_to(): (RS-index) (GAP) bitset", rs_max);
        for (unsigned j = 0; j < rs_max; j++)
        {
            for (size_t i = 0; i < sample_vec.size(); ++i)
            {
                auto idx = sample_vec[i];
                value += bv->count_to_test(idx, *bc_arr);
            }
        }
    }

    {
        bm::chrono_taker<std::ostream> tt(cout, "count_to_test(): (RS-index) (GAP) bitset", rs_max);
        for (unsigned j = 0; j < rs_max; j++)
        {
            for (size_t i = 0; i < sample_vec.size(); ++i)
            {
                auto idx = sample_vec[i];
                value_t += bv->count_to_test(idx, *bc_arr);
            }
        }
    }

    {
        bm::chrono_taker<std::ostream> tt(cout, "count_to(): (RS-index) (BITS) bitset", rs_max);
        for (unsigned j = 0; j < rs_max; j++)
        {
            for (size_t i = 0; i < sample_vec.size(); ++i)
            {
                auto idx = sample_vec[i];
                value2 += bv_c->count_to(idx, *bc_arr_c);
            }
        }
        assert(value == value2);
        assert(value);
    }

    {
        bm::chrono_taker<std::ostream> tt(cout, "count_to_test(): (RS-index) (BITS) bitset", rs_max);
        for (unsigned j = 0; j < rs_max; j++)
        {
            for (size_t i = 0; i < sample_vec.size(); ++i)
            {
                auto idx = sample_vec[i];
                value2_t += bv_c->count_to_test(idx, *bc_arr_c);
            }
        }

        assert(value_t == value2_t);
        assert(value_t);
    }


    // validation
    {
        bm::chrono_taker<std::ostream> tt(cout, "count_to(): (RS-index) (GAP+BITS) bitset", rs_max);
        for (unsigned j = 0; j < rs_max; j++)
        {
            for (size_t i = 0; i < sample_vec.size(); ++i)
            {
                auto idx = sample_vec[i];
                auto cnt   = bv->count_to(idx, *bc_arr);
                auto cnt_c = bv_c->count_to(idx, *bc_arr_c);

                if(cnt != cnt_c)
                {
                    std::cerr << "Error in count_to()!" << std::endl;
                    assert(0); exit(1);
                }
            }
        }
    }

    bv->freeze();
    bv_c->freeze();

    {
        bm::chrono_taker<std::ostream> tt(cout, "count_to(): (RS-index) (GAP+BITS) bitset (RO) ", rs_max);
        for (unsigned j = 0; j < rs_max; j++)
        {
            for (size_t i = 0; i < sample_vec.size(); ++i)
            {
                auto idx = sample_vec[i];
                auto cnt   = bv->count_to(idx, *bc_arr);
                auto cnt_c = bv_c->count_to(idx, *bc_arr_c);

                if(cnt != cnt_c)
                {
                    std::cerr << "Error in count_to() (RO)!" << std::endl;
                    assert(0); exit(1);
                }
            }
        }
    }

    delete bv;
    delete bv_c;
    delete bset;
}

static
void BitTestSparseTest()
{
    unique_ptr<bvect>  bv0(new bvect());
    unique_ptr<bvect>  bv1(new bvect());
    unique_ptr<bvect>  bv2(new bvect());
    unique_ptr<test_bitset>  bset0(new test_bitset());
    unique_ptr<test_bitset>  bset1(new test_bitset());
    unique_ptr<test_bitset>  bset2(new test_bitset());

    const unsigned repeats = REPEATS * 300000;

    size_t value = 0;

    SimpleFillSets(bset0.get(), *bv0, 0, BSIZE, 9530);
    SimpleFillSets(bset1.get(), *bv1, 0, BSIZE, 1000);
    SimpleFillSets(bset2.get(), *bv2, 0, BSIZE, 120);

    std::vector<unsigned> idx;
    idx.resize(repeats);
    for (unsigned i = 0; i < repeats; ++i)
    {
        idx.push_back(unsigned(rand_dis(gen)));
    }

    {
        bm::chrono_taker<std::ostream> tt(cout, "BitTest: bvector<>::test() (BIT) ", repeats);
        for (unsigned i = 0; i < repeats; ++i)
        {
            unsigned id = idx[i];
            value += bv0->test(id);
            value += bv1->test(id);
            value += bv2->test(id);
        }
    }


    size_t value2 = 0;

    BM_DECLARE_TEMP_BLOCK(tb)
    bv0->optimize(tb);
    bv1->optimize(tb);
    bv2->optimize(tb);

    {
        bm::chrono_taker<std::ostream> tt(cout, "BitTest: bvector<>::test() (GAP) ", repeats);
        for (unsigned i = 0; i < repeats; ++i)
        {
            unsigned id = idx[i];
            value2 += bv0->test(id);
            value2 += bv1->test(id);
            value2 += bv2->test(id);
        }
    }
    assert(value == value2);

}

static
void BitSetConditionalTest()
{
    const unsigned repeats = 10;
    std::random_device rd;
    std::mt19937 g(rd());

    std::vector<bvect::size_type> idx_vect;
    for (bvect::size_type i = 0; i < bm::id_max-256; i+=100)
        idx_vect.push_back(i);

    std::shuffle(idx_vect.begin(), idx_vect.end(), g);

    {
        bvect bv2(bm::BM_BIT);
        bv2.invert();
        bvect bv3(bm::BM_GAP);
        bv3.invert();

        {
            bm::chrono_taker<std::ostream> tt(cout, "BitSetConditional: bvector<>::set_bit_conditional() (BIT) ", repeats);
            for (bvect::size_type i = 0; i < idx_vect.size(); ++i)
            {
                auto idx = idx_vect[i];
                bool changed = bv2.set_bit_conditional(idx, false, true);
                if (!changed) {
                    cout << "Conditional bit set failed. " << i << endl;
                    exit(1);
                }
                changed = bv2.set_bit_conditional(idx, false, true);
                if (changed) {
                    cout << "Conditional bit set failed. " << i << endl;
                    exit(1);
                }
            } // for
        }
        unsigned c = 0;
        {
            bv2.invert();
            auto cnt = bv2.count();
            c+=cnt;
            assert(cnt == idx_vect.size());
            bv2.clear(true);
        }
        {
            bm::chrono_taker<std::ostream> tt(cout, "BitSetConditional: bvector<>::set_bit_conditional() (GAP) ", repeats);
            for (bvect::size_type i = 0; i < idx_vect.size(); ++i)
            {
                auto idx = idx_vect[i];
                bool changed = bv3.set_bit_conditional(idx, false, true);
                if (!changed) {
                    cout << "Conditional bit set failed. " << i << endl;
                    exit(1);
                }
                changed = bv3.set_bit_conditional(idx, false, true);
                if (changed) {
                    cout << "Conditional bit set failed. " << i << endl;
                    exit(1);
                }
            } // for
        }
        {
            bv3.invert();
            auto cnt = bv3.count();
            c+=cnt;
            assert(cnt == idx_vect.size());
        }
        if (!c) // this is to fool compiler from de-optimizing
            cout << "\r";
    }
}


static
void EnumeratorGoToTest()
{
    unique_ptr<bvect>  bv0(new bvect());
    unique_ptr<bvect>  bv1(new bvect());
    unique_ptr<bvect>  bv2(new bvect());
    unique_ptr<test_bitset>  bset0(new test_bitset());
    unique_ptr<test_bitset>  bset1(new test_bitset());
    unique_ptr<test_bitset>  bset2(new test_bitset());

    const unsigned repeats = REPEATS * 300000;

    bm::id64_t value = 0;

    SimpleFillSets(bset0.get(), *bv0, 0, BSIZE, 512);
    SimpleFillSets(bset1.get(), *bv1, 0, BSIZE, 256);
    SimpleFillSets(bset2.get(), *bv2, 0, BSIZE, 120);

    std::vector<unsigned> idx;
    idx.resize(repeats);
    for (unsigned i = 0; i < repeats; ++i)
    {
        idx.push_back(unsigned(rand_dis(gen)));
    }


    {
        bm::chrono_taker<std::ostream> tt(cout, "Enumerator at BIT pos:  ", repeats);
        for (unsigned i = 0; i < repeats; ++i)
        {
            unsigned id = idx[i];
            bvect::enumerator en0 = bv0->get_enumerator(id);
            bvect::enumerator en1 = bv1->get_enumerator(id);
            bvect::enumerator en2 = bv2->get_enumerator(id);

            value += *en0;
            value += *en1;
            value += *en2;
        }
    }

    bm::id64_t value2 = 0;

    BM_DECLARE_TEMP_BLOCK(tb)
    bv0->optimize(tb);
    bv1->optimize(tb);
    bv2->optimize(tb);

    {
        bm::chrono_taker<std::ostream> tt(cout, "Enumerator at GAP pos: ", repeats);
        for (unsigned i = 0; i < repeats; ++i)
        {
            unsigned id = idx[i];
            bvect::enumerator en0 = bv0->get_enumerator(id);
            bvect::enumerator en1 = bv1->get_enumerator(id);
            bvect::enumerator en2 = bv2->get_enumerator(id);

            value2 += *en0;
            value2 += *en1;
            value2 += *en2;
        }
    }
    assert(value == value2);
}


static
void BitCompareTest()
{
    {
    bvect*  bv1 = new bvect();
    bvect*  bv2 = new bvect();
    test_bitset*  bset = new test_bitset();
    int value = 0;

    SimpleFillSets(bset, *bv1, 0, BSIZE, 10);
    SimpleFillSets(bset, *bv2, 0, BSIZE, 10);

    {
    bm::chrono_taker<std::ostream> tt(cout, "BitCompare: Random bitvector", REPEATS*10);
    for (unsigned int i = 0; i < REPEATS*10; ++i)
    {    
        value+=bv1->compare(*bv2);
    }
    }

    delete bset;
    delete bv1;
    delete bv2;

    }

    if (platform_test) return;

    unsigned cnt = REPEATS * 100000;
    unsigned* arr1 = new unsigned[cnt];
    unsigned* arr2 = new unsigned[cnt];

    unsigned i;
    for (i = 0; i < cnt; ++i)
    {
        if ((rand() % 10) == 0)
        {
            arr1[i] = 0;
        }
        else 
        {
            arr1[i] = unsigned(rand());
            arr2[i] = unsigned(rand());
        }
    }

    {
    bm::chrono_taker<std::ostream> tt(cout, "wordcmp complex: Random words comparison", cnt);

    for (i = 0; i < cnt; ++i)
    {    
        int res2 = bm::wordcmp(arr1[i], arr2[i]);
        int res = bm::wordcmp0(arr1[i], arr2[i]);

        if (res != res2)
        {
            cerr << "Incorrect result ! " << arr1[i] 
                 << "<=>" << arr2[i] << " res=" << res <<
                 endl;
            exit(1);
        }
    }
    }

    int c = 0;
    volatile void* p = &c;

    {
    bm::chrono_taker<std::ostream> tt(cout, "wordcmp0. Random words comparison", cnt);
    for (i = 0; i < cnt; ++i)
    {    
        c += bm::wordcmp0(arr1[i], arr2[i]);
    }
    }


    {
    bm::chrono_taker<std::ostream> tt(cout, "wordcmp. Random words comparison", cnt);
    for (i = 0; i < cnt; ++i)
    {    
        c += bm::wordcmp(arr1[i], arr2[i]);
    }
    }

    c = 0;;

    delete [] arr1;
    delete [] arr2;

    char buf[256];
    sprintf(buf, "%p", p);
}

unsigned long long g_acc = 0;
extern "C" {
    static
    int bit_visitor_func(void* handle_ptr, bm::id_t bit_idx)
    {
        std::vector<bvect::size_type>* vp = (std::vector<bm::id_t>*)handle_ptr;
        (void)vp;
        //vp->push_back(bit_idx);
        g_acc += bit_idx;
        return 0;
    }
} // extern C


static
void FindTest()
{
    bvect                 bv1, bv2, bv3, bv4;
    bvect                 bv_empty;
    test_bitset*  bset = new test_bitset();
    
    bv_empty[1] = true;
    bv_empty[1] = false;
    bv_empty[100000000] = true;
    bv_empty[100000000] = false;
    
    SimpleFillSets(bset, bv1, 0, BSIZE/2, 3);
    SimpleFillSets(bset, bv1, 0, BSIZE/2, 4);
    SimpleFillSets(bset, bv1, 0, BSIZE/2, 5);
    
    FillSetsIntervals(bset, bv2, 0, BSIZE/2, 8);
    FillSetsIntervals(bset, bv3, 0, BSIZE/2, 12);
    FillSetsIntervals(bset, bv4, 0, BSIZE/2, 120);
    
    unsigned i;
    unsigned pos_sum = 0;
    {
        bm::chrono_taker<std::ostream> tt(cout, "bvector<>::find_reverse()", REPEATS*100);
        for (i = 0; i < REPEATS*100; ++i)
        {
            bm::id_t pos;
            bool found;
            found = bv1.find_reverse(pos);
            if (found)
                pos_sum += pos;
            found = bv2.find_reverse(pos);
            if (found)
                pos_sum += pos;
            found = bv3.find_reverse(pos);
            if (found)
                pos_sum += pos;
            found = bv4.find_reverse(pos);
            if (found)
                pos_sum += pos;
            found = bv_empty.find_reverse(pos);
            if (!found)
                pos_sum += pos;
            else
            {
                cerr << "incorrect find result!" << endl;
                exit(1);
            }
        }
    }
    char cbuf[256];
    sprintf(cbuf, "%i ", pos_sum); // attempt to avoid agressive optmizations

    {
        bm::chrono_taker<std::ostream> tt(cout, "bvector<>::find()", REPEATS*100);
        for (i = 0; i < REPEATS*100; ++i)
        {
            bm::id_t pos;
            bool found;
            found = bv1.find(0, pos);
            if (found)
                pos_sum += pos;
            found = bv2.find(0, pos);
            if (found)
                pos_sum += pos;
            found = bv3.find(0, pos);
            if (found)
                pos_sum += pos;
            found = bv4.find(0, pos);
            if (found)
                pos_sum += pos;
            found = bv_empty.find(0, pos);
            if (!found)
                pos_sum += pos;
            else
            {
                cerr << "incorrect find result!" << endl;
                exit(1);
            }
        }
    }

    delete bset;
}

static
void EnumeratorTest()
{
    bvect                 bv1, bv2, bv3, bv4;
    bvect::statistics     st1, st2, st3, st4;;
    std::vector<bvect::size_type> v1,  v2,  v3,  v4;
    
    {
        test_bitset*  bset = 0;// new test_bitset();

        SimpleFillSets(bset, bv1, 0, BSIZE / 2, 3);
        SimpleFillSets(bset, bv1, 0, BSIZE / 2, 4);
        SimpleFillSets(bset, bv1, 0, BSIZE / 2, 5);

        //delete bset;
    }

    FillSetsIntervals(nullptr, bv2, 0, BSIZE/2, 8);
    FillSetsIntervals(nullptr, bv3, 0, BSIZE/2, 12);
    FillSetsIntervals(nullptr, bv4, 0, BSIZE/2, 120);

    bv1.calc_stat(&st1);
    bv2.calc_stat(&st2);
    bv3.calc_stat(&st3);
    bv4.calc_stat(&st4);


    unsigned i;

    {
        unsigned long long acc = 0;
        bm::chrono_taker<std::ostream> tt(cout, "bvector<>::enumerator", REPEATS/10);
        for (i = 0; i < REPEATS/10; ++i)
        {
            {
                bvect::enumerator en = bv1.first();
                for (;en.valid();++en)
                {
                    auto v = *en;
                    acc += v;
                }
            }
            {
                bvect::enumerator en = bv2.first();
                for (;en.valid();++en)
                {
                    auto v = *en;
                    acc ^= v;
                }
            }
            {
                bvect::enumerator en = bv3.first();
                for (;en.valid();++en)
                {
                    auto v = *en;
                    acc += v;
                }
            }
            {
                bvect::enumerator en = bv4.first();
                for (;en.valid();++en)
                {
                    auto v = *en;
                    acc ^= v;
                }
            }

        } // for REPEATS
    }
    


    // -----------------------------------------------
    {
        unsigned long long acc = 0;

        bm::chrono_taker<std::ostream> tt(cout, "bvector<>::get_next()", REPEATS/10);
        for (i = 0; i < REPEATS/10; ++i)
        {
            if (bv1.any())
            {
                unsigned v = bv1.get_first();
                do
                {
                    v = bv1.get_next(v);
                    acc += v;
                    //v1.push_back(v);
                } while(v);
            }

            if (bv2.any())
            {
                unsigned v = bv2.get_first();
                do
                {
                    v = bv2.get_next(v);
                    acc ^= v;
                } while(v);
            }
            if (bv3.any())
            {
                unsigned v = bv3.get_first();
                do
                {
                    v = bv3.get_next(v);
                    acc += v;
                } while(v);
            }
            if (bv4.any())
            {
                unsigned v = bv4.get_first();
                do
                {
                    v = bv4.get_next(v);
                    acc ^= v;
                } while(v);
            }
        } // for REPEATS
    }

    // -----------------------------------------------
    
    
    {
        bm::chrono_taker<std::ostream> tt(cout, "bm::visit_each_bit()", REPEATS/10);
        for (i = 0; i < REPEATS/10; ++i)
        {
            bm::visit_each_bit(bv1, (void*)&v1, bit_visitor_func);
            bm::visit_each_bit(bv2, (void*)&v2, bit_visitor_func);
            bm::visit_each_bit(bv3, (void*)&v3, bit_visitor_func);
            bm::visit_each_bit(bv4, (void*)&v4, bit_visitor_func);

        } // for REPEATS

    }

    bv1.freeze();
    bv2.freeze();
    bv3.freeze();
    bv4.freeze();

    {
        bm::chrono_taker<std::ostream> tt(cout, "bm::visit_each_bit() (RO)", REPEATS/10);
        for (i = 0; i < REPEATS/10; ++i)
        {
            bm::visit_each_bit(bv1, (void*)&v1, bit_visitor_func);
            bm::visit_each_bit(bv2, (void*)&v2, bit_visitor_func);
            bm::visit_each_bit(bv3, (void*)&v3, bit_visitor_func);
            bm::visit_each_bit(bv4, (void*)&v4, bit_visitor_func);
        } // for REPEATS
    }

    
    // -----------------------------------------------


}

static
void EnumeratorTestGAP()
{
    bvect*  bv = new bvect();
    test_bitset*  bset = new test_bitset();
    unsigned i;

    SimpleFillSets(bset, *bv, 0, BSIZE, 2500);
    bv->optimize();

    unsigned k = 1;
    unsigned cnt = 0;
    {
    bm::chrono_taker<std::ostream> tt(cout, "bvector<>::get_next() (GAP) ", REPEATS*10*(k+1));
    for (i = 0; i < REPEATS*10*(k+1); ++i)
    {
        if (bv->any())
        {
            unsigned v = bv->get_first();
            do
            {
                v = bv->get_next(v);
                cnt += v;
            } while (v);
        }
    }
    }

    unsigned v = 0;
    {
    bm::chrono_taker<std::ostream> tt(cout, "bvector<>::enumerator (GAP) ", REPEATS*10*(k+1));
    for (i = 0; i < REPEATS*10*(k+1); ++i)
    {    
        bvect::enumerator en = bv->first();
        bvect::enumerator bend = bv->end();

        while (en < bend)
        {
            v += *en;
            ++en;
        }
    }
    }
    stringstream s;
    s << v << endl; // attempt to fool optimization

    bv->freeze();
    {
    bm::chrono_taker<std::ostream> tt(cout, "bvector<>::enumerator (GAP) (RO)", REPEATS*10*(k+1));
    for (i = 0; i < REPEATS*10*(k+1); ++i)
    {
        bvect::enumerator en = bv->first();
        bvect::enumerator bend = bv->end();

        while (en < bend)
        {
            v += *en;
            ++en;
        }
    }
    }
    s << v << endl; // attempt to fool optimization

    char buf[256];
    sprintf(buf, "%i", cnt); // to fool some smart compilers like ICC

    delete bv;
    delete bset;
}

static
void SerializationTest()
{
    bvect bv_sparse;
    // stack declaration of temp block instead of re-allocation makes things faster
    BM_DECLARE_TEMP_BLOCK(tb)


    // prepare a test bitset with a small number of bits set somewhere
    // far from the beginning
    for (unsigned i = 0; i < 5; ++i)
    {
        bv_sparse[BSIZE/2 + i * 3] = true;		
    }
    bv_sparse[100] = true;
    bv_sparse[70000] = true;
    bv_sparse[200000] = true;

    bv_sparse.optimize(tb);

    unsigned cnt = bv_sparse.count();
    bvect::statistics st;
    bv_sparse.calc_stat(&st);
    unsigned char*  buf = new unsigned char[st.max_serialize_mem];

    size_t len, id_size;
    len = id_size = 0;
    {
    bm::chrono_taker<std::ostream> tt(cout, "Small bvector serialization", REPEATS*70000);
    for (unsigned i = 0; i < REPEATS*70000; ++i)
    {
        len += bm::serialize(bv_sparse, buf, tb, bm::BM_NO_BYTE_ORDER|bm::BM_NO_GAP_LENGTH);
        id_size += cnt * (unsigned)sizeof(unsigned);
    }
    }
    
    delete [] buf; buf = 0;
        
    bvect*  bv = new bvect();
    //test_bitset*  bset = new test_bitset();
    unsigned value = 0;

    SimpleFillSets(nullptr, *bv, 0, BSIZE, 4);
    
    cnt = bv->count();
    bv->calc_stat(&st);
    buf = new unsigned char[st.max_serialize_mem];
    
    {
    bm::chrono_taker<std::ostream> tt(cout, "Large bvector serialization", REPEATS/3);
    for (unsigned i = 0; i < REPEATS/3; ++i)
    {
        len += bm::serialize(*bv, buf, tb, bm::BM_NO_BYTE_ORDER|bm::BM_NO_GAP_LENGTH);
        id_size += cnt * (unsigned)sizeof(unsigned);
    }
    }


    char cbuf[256] = {0, };
    sprintf(cbuf, "%u", value);
    /*
    cout << cbuf << " " << id_size << " " << len << " " << value << endl;
    */
        
    delete bv;
    delete [] buf;
}

static
void InvertTest()
{
    bvect*  bv = new bvect();
    test_bitset*  bset = new test_bitset();
    unsigned i;
    //unsigned value = 0;

    SimpleFillSets(bset, *bv, 0, BSIZE, 2500);
    {
    bm::chrono_taker<std::ostream> tt(cout, "Invert bvector", REPEATS*4);
    for (i = 0; i < REPEATS*4; ++i)
    {
        bv->flip();    
    }
    }

    if (!platform_test)
    {
    bm::chrono_taker<std::ostream> tt(cout, "Invert bvector (STL)", REPEATS*4);
    for (i = 0; i < REPEATS*4; ++i)
    {
        bset->flip();    
    }
    }

    delete bv;
    delete bset;
}

static
void OrTest()
{
    bvect bv1, bv2;
    bvect bvt1, bvt2, bvt3;
    generate_bvector(bv1, 40000000, false);
    generate_bvector(bv2, 40000000, false);
    
    {
    bm::chrono_taker<std::ostream> tt(cout, "OR-optimize (2 operand) bvector test", REPEATS*4);
    for (unsigned i = 0; i < REPEATS*4; ++i)
    {
        bvt1 = bv1;
        bvt1 |= bv2;
        bvt1.optimize();
    }
    }
    
    {
    bm::chrono_taker<std::ostream> tt(cout, "OR-optimize (3 operand) bvector test", REPEATS*4);
    for (unsigned i = 0; i < REPEATS*4; ++i)
    {
        bvt2.bit_or(bv1, bv2, bvect::opt_compress);
    }
    }

    int cmp = bvt1.compare(bvt2);
    if (cmp)
    {
        cerr << "Error: OR mismatch!" << endl;
        exit(1);
    }

    bv1.freeze();
    bv2.freeze();

    {
    bm::chrono_taker<std::ostream> tt(cout, "OR-optimize (3 operand) bvector test (RO)", REPEATS*4);
    for (unsigned i = 0; i < REPEATS*4; ++i)
    {
        bvt3.bit_or(bv1, bv2, bvect::opt_compress);
    }
    }

    bool eq = bvt1.equal(bvt3);
    if (!eq)
    {
        cerr << "Error: RO OR mismatch!" << endl;
        exit(1);
    }

}

static
void AndTest()
{
    bvect*  bv1 = new bvect();
    test_bitset*  bset1 = new test_bitset();
    test_bitset*  bset2 = new test_bitset();
    bvect*  bv2 = new bvect();
    unsigned i;
    //unsigned value = 0;

    SimpleFillSets(bset1, *bv1, 0, BSIZE, 100);
    SimpleFillSets(bset1, *bv2, 0, BSIZE, 100);

    {
    bm::chrono_taker<std::ostream> tt(cout, "AND bvector test", REPEATS*4);
    for (i = 0; i < REPEATS*4; ++i)
    {
        *bv1 &= *bv2;
    }
    }

    if (!platform_test)
    {
    bm::chrono_taker<std::ostream> tt(cout, "AND bvector test(STL)", REPEATS*4);
    for (i = 0; i < REPEATS*4; ++i)
    {
        *bset1 &= *bset2;
    }
    }


    delete bv1;
    delete bv2;

    delete bset1;
    delete bset2;
}

static
void OrAndTest()
{
    bvect bv_res1, bv_res2, bv_res3;
    bvect*  bv1 = new bvect();
    test_bitset*  bset1 = new test_bitset();
    test_bitset*  bset2 = new test_bitset();
    bvect*  bv2 = new bvect();
    unsigned i;
    //unsigned value = 0;

    SimpleFillSets(bset1, *bv1, 0, BSIZE, 100);
    SimpleFillSets(bset1, *bv2, 0, BSIZE, 100);

    {
    bm::chrono_taker<std::ostream> tt(cout, "AND+OR bvector test", REPEATS*4);
    for (i = 0; i < REPEATS*4; ++i)
    {
        bvect bv;
        bv.bit_and(*bv1, *bv2);
        bv_res1 |= bv;
    }
    }

    {
    bm::chrono_taker<std::ostream> tt(cout, "AND_OR( fused) bvector test", REPEATS*4);
    for (i = 0; i < REPEATS*4; ++i)
    {
        bv_res2.bit_or_and(*bv1, *bv2);
    }
    }

    bool b = bv_res1.equal(bv_res2);
    if (!b)
    {
        cerr << "OR-AND test failed!" << endl;
        exit(1);
    }

    bv1->freeze();
    bv2->freeze();

    {
    bm::chrono_taker<std::ostream> tt(cout, "AND_OR( fused) bvector test (RO)", REPEATS*4);
    for (i = 0; i < REPEATS*4; ++i)
    {
        bv_res3.bit_or_and(*bv1, *bv2);
    }
    }

    b = bv_res1.equal(bv_res3);
    if (!b)
    {
        cerr << "OR-AND (RO) test failed!" << endl;
        exit(1);
    }

    delete bv1;
    delete bv2;

    delete bset1;
    delete bset2;
}


static
void XorTest()
{
    bvect*  bv1 = new bvect();
    test_bitset*  bset1 = new test_bitset();
    test_bitset*  bset2 = new test_bitset();
    bvect*  bv2 = new bvect();
    unsigned i;
    //unsigned value = 0;

    SimpleFillSets(bset1, *bv1, 0, BSIZE, 100);
    SimpleFillSets(bset1, *bv2, 0, BSIZE, 100);

    {
        bm::chrono_taker<std::ostream> tt(cout, "XOR bvector test", REPEATS * 10);
        for (i = 0; i < REPEATS * 4; ++i)
        {
            *bv1 ^= *bv2;
        }
    }

    bv2->freeze();

    {
        bm::chrono_taker<std::ostream> tt(cout, "XOR bvector test (RO)", REPEATS * 10);
        for (i = 0; i < REPEATS * 4; ++i)
        {
            *bv1 ^= *bv2;
        }
    }

    if (!platform_test)
    {
        bm::chrono_taker<std::ostream> tt(cout, "XOR bvector test(STL)", REPEATS * 10);
        for (i = 0; i < REPEATS * 4; ++i)
        {
            *bset1 ^= *bset2;
        }
    }

    delete bv1;
    delete bv2;

    delete bset1;
    delete bset2;
}

static
void SubTest()
{
    bvect bv1, bv2;
    bvect bvt0, bvt1, bvt2, bvt3;
    generate_bvector(bv1, 40000000, false);
    generate_bvector(bv2, 40000000, false);
    
    {
    bm::chrono_taker tt(cout, "AND-NOT bvector test", REPEATS*4);
    for (unsigned i = 0; i < REPEATS*4; ++i)
    {
        bvect bv_tmp(bv2);
        bv_tmp.invert();
        
        bvt0 = bv1;
        bvt0 &= bv_tmp;
        bvt0.optimize();
    }

    }
    
    {
    bm::chrono_taker tt(cout, "SUB-optimize (2 operand) bvector test", REPEATS*4);
    for (unsigned i = 0; i < REPEATS*4; ++i)
    {
        bvt1 = bv1;
        bvt1 -= bv2;
        bvt1.optimize();
    }
    }
    
    {
    bm::chrono_taker tt(cout, "SUB-optimize (3 operand) bvector test", REPEATS*4);
    for (unsigned i = 0; i < REPEATS*4; ++i)
    {
        bvt2.bit_sub(bv1, bv2, bvect::opt_compress);
    }
    }

    bv1.freeze();
    bv2.freeze();
    {
    bm::chrono_taker tt(cout, "SUB-optimize (3 operand) bvector test (RO)", REPEATS*4);
    for (unsigned i = 0; i < REPEATS*4; ++i)
    {
        bvt3.bit_sub(bv1, bv2, bvect::opt_compress);
    }
    }


    int cmp = bvt1.compare(bvt2);
    if (cmp)
    {
        cerr << "Error: (1)SUB mismatch!" << endl;
        exit(1);
    }
    cmp = bvt0.compare(bvt2);
    if (cmp)
    {
        cerr << "Error: (2)SUB mismatch!" << endl;
        exit(1);
    }
    cmp = bvt0.compare(bvt2);
    if (cmp)
    {
        cerr << "Error: (2)SUB mismatch (RO)!" << endl;
        exit(1);
    }

}

static
void XorCountTest()
{
    bvect*  bv1 = new bvect();
    bvect*  bv2 = new bvect();
    test_bitset*  bset1 = new test_bitset();
    test_bitset*  bset2 = new test_bitset();
    unsigned i;

    SimpleFillSets(bset1, *bv1, 0, BSIZE, 400);
    SimpleFillSets(bset2, *bv2, 0, BSIZE, 500);

    unsigned count1 = 0;
    unsigned count2 = 0;
    unsigned count3 = 0;
    unsigned test_count = 0;

    if (!platform_test)
    {
    bvect bv_tmp;
    bm::chrono_taker tt(cout, "XOR COUNT bvector test with TEMP vector", REPEATS*10);
    for (i = 0; i < REPEATS*4; ++i)
    {
        bv_tmp.clear();
        bv_tmp |= *bv1;
        bv_tmp ^= *bv2;
        count1 += bv_tmp.count();
    }
    }

    if (!platform_test)
    {
    test_bitset*  bset_tmp = new test_bitset();
    bm::chrono_taker tt(cout, "XOR COUNT bvector test with TEMP vector (STL)", REPEATS*10);
    for (i = 0; i < REPEATS*4; ++i)
    {
        bset_tmp->reset();
        *bset_tmp |= *bset1;
        *bset_tmp ^= *bset2;
        test_count += (unsigned)bset_tmp->count();
    }
    }


    {
    bm::chrono_taker tt(cout, "XOR COUNT bvector test", REPEATS*10);
    for (i = 0; i < REPEATS*4; ++i)
    {
        count2 += bm::count_xor(*bv1, *bv2);
    }
    }

//    assert(count1 == count2);

    
    if (!platform_test)    
        if (count1 != count2)
        {
            cout << "Check failed !" << endl;
            cout << count1 << " " << count2 << " " << test_count << endl;
            exit(1);
        }
    count1 = count2 = 0;
    
    // -----------------------------------------
    if (!platform_test)
    {
        cout << "One optimized vector" << endl;
    }
    BM_DECLARE_TEMP_BLOCK(tb)
    bv2->optimize(tb);

    if (!platform_test)
    {
    bvect bv_tmp;
    bm::chrono_taker tt(cout, "XOR COUNT bvector test with TEMP vector", REPEATS*10);
    for (i = 0; i < REPEATS*4; ++i)
    {
        bv_tmp.clear(false);
        bv_tmp |= *bv1;
        bv_tmp ^= *bv2;
        count1 += (unsigned)bv_tmp.count();
    }
    }

    {
    bm::chrono_taker tt(cout, "XOR COUNT bvector(opt1) test", REPEATS*10);
    for (i = 0; i < REPEATS*4; ++i)
    {
        count2 += (unsigned)bm::count_xor(*bv1, *bv2);
    }
    }

    if (!platform_test)
        if (count1 != count2)
        {
            cout << "Check failed !" << endl;
            exit(1);
        }
    count1 = count2 = 0;

    // -----------------------------------------
    if (!platform_test)
    {
        cout << "Both vectors optimized" << endl;
    }
    bv1->optimize(tb);
    //bv1->stat();
    if (!platform_test)
    {
    bvect bv_tmp;
    bm::chrono_taker tt(cout, "XOR COUNT bvector test with TEMP vector", REPEATS*10);
    for (i = 0; i < REPEATS*4; ++i)
    {
        bv_tmp.clear(false);
        bv_tmp |= *bv1;
        bv_tmp ^= *bv2;
        count1 += (unsigned)bv_tmp.count();
    }
    }

    {
    bm::chrono_taker tt(cout, "XOR COUNT bvector(opt2) test", REPEATS*10);
    for (i = 0; i < REPEATS*4; ++i)
    {
        count2 += (unsigned)bm::count_xor(*bv1, *bv2);
    }
    }
    if (!platform_test)
        if (count1 != count2)
        {
            cout << "Check failed !" << endl;
            exit(1);
        }

    bv1->freeze();
    bv2->freeze();

    {
    bm::chrono_taker tt(cout, "XOR COUNT bvector(opt2) test (RO)", REPEATS*10);
    for (i = 0; i < REPEATS*4; ++i)
    {
        count3 += (unsigned)bm::count_xor(*bv1, *bv2);
    }
    }
    assert(count2 == count3);



    delete bv1;
    delete bv2;
    
    delete bset1;
    delete bset2;    
}

static
void AndCountTest()
{
    bvect*  bv1 = new bvect();
    bvect*  bv2 = new bvect();
    test_bitset*  bset1 = new test_bitset();
    test_bitset*  bset2 = new test_bitset();
    unsigned i;

    SimpleFillSets(bset1, *bv1, 0, BSIZE, 400);
    SimpleFillSets(bset2, *bv2, 0, BSIZE, 500);

    unsigned count1 = 0;
    unsigned count2 = 0;
    unsigned count3 = 0;
    unsigned test_count = 0;

    if (!platform_test)
    {
        bvect bv_tmp;
        bm::chrono_taker tt(cout, "AND COUNT bvector test with TEMP vector", REPEATS * 10);
        for (i = 0; i < REPEATS * 4; ++i)
        {
            bv_tmp.clear(false);
            bv_tmp |= *bv1;
            bv_tmp &= *bv2;
            count1 += bv_tmp.count();
        }
    }

    if (!platform_test)
    {
        test_bitset*  bset_tmp = new test_bitset();
        bm::chrono_taker tt(cout, "AND COUNT bvector test with TEMP vector (STL)", REPEATS * 10);
        for (i = 0; i < REPEATS * 4; ++i)
        {
            bset_tmp->reset();
            *bset_tmp |= *bset1;
            *bset_tmp &= *bset2;
            test_count += (unsigned)bset_tmp->count();
        }
    }


    {
        bm::chrono_taker tt(cout, "AND COUNT bvector test", REPEATS * 10);
        for (i = 0; i < REPEATS * 4; ++i)
        {
            count2 += bm::count_and(*bv1, *bv2);
        }
    }


    if (!platform_test)
        if (count1 != count2)
        {
            cout << "Check failed !" << endl;
            cout << count1 << " " << count2 << " " << test_count << endl;
            exit(1);
        }
    count1 = count2 = 0;

    // -----------------------------------------
    if (!platform_test)
    {
        cout << "One optimized vector" << endl;
    }
    BM_DECLARE_TEMP_BLOCK(tb)
        bv2->optimize(tb);

    if (!platform_test)
    {
        bvect bv_tmp;
        bm::chrono_taker tt(cout, "AND COUNT bvector test with TEMP vector", REPEATS * 10);
        for (i = 0; i < REPEATS * 4; ++i)
        {
            bv_tmp.clear(false);
            bv_tmp |= *bv1;
            bv_tmp &= *bv2;
            count1 += (unsigned)bv_tmp.count();
        }
    }

    {
        bm::chrono_taker tt(cout, "AND COUNT bvector test", REPEATS * 10);
        for (i = 0; i < REPEATS * 4; ++i)
        {
            count2 += (unsigned)bm::count_and(*bv1, *bv2);
        }
    }

    if (!platform_test)
        if (count1 != count2)
        {
            cout << "Check failed !" << endl;
            exit(1);
        }
    count1 = count2 = 0;

    // -----------------------------------------
    if (!platform_test)
    {
        cout << "Both vectors optimized" << endl;
    }
    bv1->optimize(tb);
    //bv1->stat();
    if (!platform_test)
    {
        bvect bv_tmp;
        bm::chrono_taker tt(cout, "AND COUNT bvector test with TEMP vector", REPEATS * 10);
        for (i = 0; i < REPEATS * 4; ++i)
        {
            bv_tmp.clear(false);
            bv_tmp |= *bv1;
            bv_tmp &= *bv2;
            count1 += (unsigned)bv_tmp.count();
        }
    }

    {
        bm::chrono_taker tt(cout, "AND COUNT bvector(opt) test", REPEATS * 10);
        for (i = 0; i < REPEATS * 4; ++i)
        {
            count2 += (unsigned)bm::count_and(*bv1, *bv2);
        }
    }
    if (!platform_test)
        if (count1 != count2)
        {
            cout << "Check failed !" << endl;
            exit(1);
        }

    bv1->freeze();
    bv2->freeze();

    {
        bm::chrono_taker tt(cout, "AND COUNT bvector(opt) test (RO)", REPEATS * 10);
        for (i = 0; i < REPEATS * 4; ++i)
        {
            count3 += (unsigned)bm::count_and(*bv1, *bv2);
        }
    }
    assert(count3 == count2);


    count1 = count2 = 0;


    delete bv1;
    delete bv2;

    delete bset1;
    delete bset2;
}


static
void TI_MetricTest()
{
    bvect*  bv1 = new bvect();
    bvect*  bv2 = new bvect();
    test_bitset*  bset1 = new test_bitset();
    test_bitset*  bset2 = new test_bitset();
    unsigned i;

    SimpleFillSets(bset1, *bv1, 0, BSIZE, 500);
    SimpleFillSets(bset2, *bv2, 0, BSIZE, 250);

    unsigned count1 = 0;
    unsigned count2 = 0;
    unsigned countA=0, countB=0, test_countA=0, test_countB=0;
    unsigned test_count = 0;
    double ti1=0, ti2=0, ti3=0;
    double di1=0, di2=0;
    {
    bm::chrono_taker tt(cout, "Tversky Index bvector test vector", REPEATS);
    for (i = 0; i < REPEATS; ++i)
    {
        count1 = bm::count_and(*bv1, *bv2);
        
        countA = bm::count_sub(*bv1, *bv2);
        countB = bm::count_sub(*bv2, *bv1);
        
        ti1 = double(count1) / double(0.4*countA + 0.5*countB + count1);
    }
    }


    if (!platform_test)
    {
    test_bitset*  bset_tmp = new test_bitset();
    double test_dice = 0;
    bm::chrono_taker tt(cout, "Dice bvector test with TEMP vector(STL)", REPEATS);
    for (i = 0; i < REPEATS; ++i)
    {
        bset_tmp->reset();
        *bset_tmp |= *bset1;
        *bset_tmp &= *bset2;
        test_count += (unsigned)bset_tmp->count();
        
        test_countA += (unsigned)bset1->count();
        test_countB += (unsigned)bset2->count();
        
        test_countA += (unsigned)bset1->count();
        test_countB += (unsigned)bset2->count();
        
        test_dice += double(2*test_count) / double(test_countA + test_countB);
    }
    }


    {
    bm::distance_metric_descriptor dmd[3];
    dmd[0].metric = bm::COUNT_AND;
    dmd[1].metric = bm::COUNT_SUB_AB;
    dmd[2].metric = bm::COUNT_SUB_BA;    
    
    bm::chrono_taker tt(cout, "Tversky Index bvector test (pipeline)", REPEATS);
    for (i = 0; i < REPEATS; ++i)
    {
        bm::distance_operation(*bv1, *bv2, &dmd[0], (&dmd[0])+3);
                
        ti2 = double(dmd[0].result) / double(0.4*dmd[1].result + 0.5*dmd[2].result + dmd[0].result);
        
        dmd[0].result = dmd[1].result = dmd[2].result = 0;
    }
    }

    
    if (fabs(ti2 - ti1) > 0.1)
    {
        cout << "Check failed ! error=" << fabs(ti2 - ti1) << endl;
        cout << ti1 << " " << ti2 << endl;
        exit(1);
    }
    count1 = count2 = 0;

    // -----------------------------------------
    if (!platform_test)
    {
        cout << "One optimized vector" << endl;
    }
    BM_DECLARE_TEMP_BLOCK(tb)
    bv2->optimize(tb);
    bv1->count(); // trying to fool the CPU cache

    
    {
    bm::chrono_taker tt(cout, "Dice metric bvector test", REPEATS);
    for (i = 0; i < REPEATS; ++i)
    {
        count1 = bm::count_and(*bv1, *bv2);
        
        countA = bm::count_sub(*bv1, *bv2);
        countB = bm::count_sub(*bv2, *bv1);
        
        di1 = double(count1) / double(0.4*countA + 0.5*countB + count1);
    }
    }

    {
        bvect bv1_ro(*bv1, bm::finalization::READONLY);
        bvect bv2_ro(*bv2, bm::finalization::READONLY);

        {
        bm::chrono_taker tt(cout, "Dice metric bvector test (RO)", REPEATS);
        for (i = 0; i < REPEATS; ++i)
        {
            count1 = bm::count_and(bv1_ro, bv2_ro);

            countA = bm::count_sub(bv1_ro, bv2_ro);
            countB = bm::count_sub(bv2_ro, bv1_ro);

            di2 = double(count1) / double(0.4*countA + 0.5*countB + count1);
        }
        }
        if (fabs(di2 - di1) > 0.1)
        {
            cout << "Dice RO Check failed !" << endl;
            cout << di1 << " " << di2 << endl;
            exit(1);
        }

    }



    {
    bm::distance_metric_descriptor dmd[3];
    dmd[0].metric = bm::COUNT_AND;
    dmd[1].metric = bm::COUNT_SUB_AB;
    dmd[2].metric = bm::COUNT_SUB_BA;    
    
    bm::chrono_taker tt(cout, "Tversky Index bvector test(pipeline)", REPEATS);
    for (i = 0; i < REPEATS; ++i)
    {
        bm::distance_operation(*bv1, *bv2, &dmd[0], (&dmd[0])+3);
                
        ti2 = double(dmd[0].result) / double(0.4*dmd[1].result + 0.5*dmd[2].result + dmd[0].result);
        
        dmd[0].result = dmd[1].result = dmd[2].result = 0;
    }
    }


    if (fabs(ti2 - ti1) > 0.1)
    {
        cout << "Check failed !" << endl;
        cout << ti1 << " " << ti2 << endl;
        exit(1);
    }
    count1 = count2 = 0;
    count1 = count2 = 0;

    // -----------------------------------------
    if (!platform_test)
    {
        cout << "Both vectors optimized" << endl;
    }
    bv1->optimize(tb);

    {
    bm::chrono_taker tt(cout, "Tversky index bvector test", REPEATS);
    for (i = 0; i < REPEATS; ++i)
    {
        count1 = bm::count_and(*bv1, *bv2);
        
        countA = bm::count_sub(*bv1, *bv2);
        countB = bm::count_sub(*bv2, *bv1);
        
        ti1 = double(count1) / double(0.4*countA + 0.5*countB + count1);
    }
    }

    {
    bm::distance_metric_descriptor dmd[3];
    dmd[0].metric = bm::COUNT_AND;
    dmd[1].metric = bm::COUNT_SUB_AB;
    dmd[2].metric = bm::COUNT_SUB_BA;    
    
    bm::chrono_taker tt(cout, "Tversky Index bvector test (pipeline)", REPEATS);
    for (i = 0; i < REPEATS; ++i)
    {
        bm::distance_operation(*bv1, *bv2, &dmd[0], (&dmd[0])+3);
                
        ti2 = double(dmd[0].result) / double(0.4*dmd[1].result + 0.5*dmd[2].result + dmd[0].result);
        
        dmd[0].result = dmd[1].result = dmd[2].result = 0;
    }
    }

    if (fabs(ti2 - ti1) > 0.1)
    {
        cout << "Check failed !" << endl;
        cout << ti1 << " " << ti2 << endl;
        exit(1);
    }

    bv1->freeze();
    bv2->freeze();


    {
    bm::distance_metric_descriptor dmd[3];
    dmd[0].metric = bm::COUNT_AND;
    dmd[1].metric = bm::COUNT_SUB_AB;
    dmd[2].metric = bm::COUNT_SUB_BA;

    bm::chrono_taker tt(cout, "Tversky Index bvector test (pipeline) (RO)", REPEATS);
    for (i = 0; i < REPEATS; ++i)
    {
        bm::distance_operation(*bv1, *bv2, &dmd[0], (&dmd[0])+3);

        ti3 = double(dmd[0].result) / double(0.4*dmd[1].result + 0.5*dmd[2].result + dmd[0].result);

        dmd[0].result = dmd[1].result = dmd[2].result = 0;
    }
    }

    if (fabs(ti3 - ti1) > 0.1)
    {
        cout << "Check failed ! (RO)" << endl;
        cout << ti1 << " " << ti3 << endl;
        exit(1);
    }



    delete bv1;
    delete bv2;
    
    delete bset1;
    delete bset2;    
}

#if 0
static
void BitBlockTransposeTest()
{
/*
    bm::word_t BM_VECT_ALIGN block1[bm::set_block_size] BM_VECT_ALIGN_ATTR = { 0, };

    for (unsigned i = 0; i < bm::set_block_size; ++i)
    {
        block1[i] = 1 | (1 << 5) | (7 << 15) | (3 << 22);
    }
*/
    unsigned   BM_VECT_ALIGN tmatrix1[32][bm::set_block_plain_size] BM_VECT_ALIGN_ATTR;

    const unsigned blocks_count = 70000;
    bm::word_t* blocks[blocks_count];
    for (unsigned k = 0; k < blocks_count; ++k)
    {
        blocks[k] = bm::block_allocator::allocate(bm::set_block_size, 0);
        for (unsigned i = 0; i < bm::set_block_size; ++i)
        {
            blocks[k][i] = 1 | (1 << 5) | (7 << 15) | (3 << 22);
        }
    }

    unsigned cnt=0;
/*
    {
    TimeTaker tt("Bit-block transpose.", REPEATS*1000);
    for (unsigned i = 0; i < REPEATS*1000; ++i)
    {
        bm::bit_block_transpose(block1, tmatrix1);
    }
    }

    {
    TimeTaker tt("Bit-block trestore.", REPEATS*1000);
    for (unsigned i = 0; i < REPEATS*1000; ++i)
    {
        bm::bit_block_trestore(tmatrix1, block2);
        cnt += block2[10];

    }
    }

    {
    TimeTaker tt("Bit-block transpose distance.", REPEATS*1000);
    unsigned distance[bm::set_block_plain_cnt][bm::set_block_plain_cnt];
    for (unsigned i = 0; i < REPEATS*1000; ++i)
    {
        bm::bit_block_tmatrix_distance(tmatrix1, distance);
        cnt += distance[1][1];
    }
    }
    printf("", cnt);
*/

    unsigned d2[bm::set_block_plain_cnt][bm::set_block_plain_cnt];
    {
    TimeTaker tt("Bit-block transpose+distance", 100000);
    unsigned distance[bm::set_block_plain_cnt][bm::set_block_plain_cnt];
    unsigned idx = 0;
    for (unsigned i = 0; i < 100000; ++i)
    {
        bm::vect_bit_transpose<unsigned, 
                               bm::set_block_plain_cnt, 
                               bm::set_block_plain_size>
                               (blocks[idx], bm::set_block_size, tmatrix1);
        bm::tmatrix_distance<unsigned, 
                             bm::set_block_plain_cnt, 
                             bm::set_block_plain_size>
                             (tmatrix1, distance);
    
        cnt += distance[1][1];
        ++idx;
        if (idx >= blocks_count) idx = 0;
        memcpy(d2, distance, sizeof(distance));
    }

    }
    
    char cbuf[256];
    sprintf(cbuf, "%i %i", cnt, d2[10][10]);

    for (unsigned i = 0; i < blocks_count; ++i)
    {
        bm::block_allocator::deallocate(blocks[i], 0);
    }

}
#endif

static
void BitBlockRotateTest()
{
    bm::word_t blk0[bm::set_block_size] = { 0 };
    bm::word_t blk1[bm::set_block_size] = { 0 };
    unsigned i;
    unsigned repeats = 20000000;

    for (i = 0; i < bm::set_block_size; ++i)
    {
        blk0[i] = blk1[i] = unsigned(rand());
    }

    {
        bm::chrono_taker tt(cout, "Bit-block left rotate 1", repeats);
        for (i = 0; i < repeats; ++i)
        {
            bm::bit_block_rotate_left_1(blk0);
        }
    }
    {
        bm::chrono_taker tt(cout, "Bit-block left rotate 1 unrolled", repeats);
        for (i = 0; i < repeats; ++i)
        {
            bm::bit_block_rotate_left_1_unr(blk1);
        }
    }

    for (i = 0; i < bm::set_block_size; ++i)
    {
        if (blk0[i] != blk1[i])
        {
            cerr << "Steress Cyclic rotate check failed" << endl;
            exit(1);
        }
    }

}

static
void BitBlockShiftTest()
{
    bm::word_t BM_VECT_ALIGN blk0[bm::set_block_size] BM_VECT_ALIGN_ATTR = { 0 };
    bm::word_t BM_VECT_ALIGN blk1[bm::set_block_size] BM_VECT_ALIGN_ATTR = { 0 };
    unsigned i;
    unsigned repeats = 20000000;
    unsigned acc0, acc1;

    for (i = 0; i < bm::set_block_size; ++i)
    {
        blk0[i] = blk1[i] = unsigned(rand());
    }

    {
        bm::chrono_taker tt(cout, "Bit-block shift-r(1)", repeats);
        {
            for (i = 0; i < repeats; ++i)
            {
                bm::bit_block_shift_r1(blk0, &acc0, 0);
            }
        }
    }

    {
        bm::chrono_taker tt(cout, "Bit-block shift-r(1) unrolled", repeats);
        for (i = 0; i < repeats; ++i)
        {
            bm::bit_block_shift_r1_unr(blk1, &acc1, 0);
        }
    }

    for (i = 0; i < bm::set_block_size; ++i)
    {
        if (blk0[i] != blk1[i])
        {
            cerr << "Stress SHIFT-r(1) check failed" << endl;
            exit(1);
        }
    }
    
    
    for (i = 0; i < bm::set_block_size; ++i)
    {
        blk0[i] = blk1[i] = unsigned(rand());
    }

    {
        bm::chrono_taker tt(cout, "Bit-block shift-l(1)", repeats);
        for (i = 0; i < repeats; ++i)
        {
            bm::bit_block_shift_l1(blk0, &acc0, 0);
        }
    }

    {
        bm::chrono_taker tt(cout, "Bit-block shift-l(1) unrolled", repeats);
        for (i = 0; i < repeats; ++i)
        {
            bm::bit_block_shift_l1_unr(blk1, &acc1, 0);
        }
    }

    for (i = 0; i < bm::set_block_size; ++i)
    {
        if (blk0[i] != blk1[i])
        {
            cerr << "Stress SHIFT-l(1) check failed" << endl;
            exit(1);
        }
    }


}


inline
void ptest()
{
    bvect*  bv_small = new bvect(bm::BM_GAP);
    bvect*  bv_large = new bvect(bm::BM_GAP);

    test_bitset*  bset = new test_bitset();

    FillSetsIntervals(bset, *bv_large, 0, 2000000000, 10);

    for (unsigned i = 0; i < 2000; ++i)
    {
        bv_small->set(i*10000);
    }

    {
    bm::chrono_taker tt(cout, "Operation &= test", REPEATS * 10);
    unsigned count = 0;
    for (unsigned i = 0; i < REPEATS*10; ++i)
    {
        bvect t1(bm::BM_GAP);
        t1 = *bv_small;
        t1 &= *bv_large;
        count += t1.count();
    }
    }



    {
    bm::chrono_taker tt(cout, "Operation &= with enumerator test", REPEATS * 10);
    unsigned count = 0;
    for (unsigned i = 0; i < REPEATS*10; ++i)
    {
        bvect t1(bm::BM_GAP);
        bvect t2(bm::BM_GAP);
        t1 = *bv_small;
        
        for (bvect::enumerator it = t1.first(); it != t1.end(); ++it) {
            if ((*bv_large)[*it]) {
                t2.set_bit(*it);
            }
        }
        count += t2.count();
    }
    }


}


typedef bm::sparse_vector<unsigned, bvect> svect;
typedef bm::sparse_vector<unsigned, bvect> sparse_vector_u32;
typedef bm::sparse_vector<int, bvect> sparse_vector_i32;

typedef bm::sparse_vector<unsigned, bvect > sparse_vector_u32;
typedef bm::sparse_vector<unsigned long long, bvect > sparse_vector_u64;
typedef bm::rsc_sparse_vector<unsigned, sparse_vector_u32> rsc_sparse_vector_u32;


// create a benchmark svector with a few dufferent distribution patterns
//
template<class SV>
void FillSparseIntervals(SV& sv)
{
    sv.resize(250000000);
    unsigned i;
    for (i = 256000; i < 712000 * 2; ++i)
    {
        if constexpr (SV::is_signed())
            sv.set(i, -0xFFE);
        else
            sv.set(i, 0xFFE);
    }
    for (i = 712000 * 3; i < 712000 * 5; ++i)
    {
        sv.set(i, (typename SV::value_type)i);
    }
    for (i = 180000000; i < 190000000; ++i)
    {
        sv.set(i, (typename SV::value_type)rand() % 128000);
    }
    for (i = 200000000; i < 210000000; ++i)
    {
        if constexpr (SV::is_signed())
            sv.set(i, (typename SV::value_type)(0 - (rand() % 128000)));
        else
            sv.set(i, (typename SV::value_type)rand() % 128000);
    }
}


template<class SV>
void FillSparseNullVector(SV& sv, typename SV::size_type size,
                          unsigned data_size, unsigned null_factor)
{
    typename SV::size_type i = 0;
    for (; i < size; i+= null_factor)
    {
        typename SV::size_type k = 0;
        for (;k < data_size; ++k, ++i)
        {
            sv[i] = 0x0DFA;
        } // for k
    } // for i
}

template<class SV>
void FillRandomSparseNullVector(SV& sv, typename SV::size_type size,
                          unsigned data_size, unsigned null_factor)
{
    typename SV::size_type i = 0;
    if (null_factor)
    {
        for (; i < size; i+= null_factor)
        {
            typename SV::size_type k = 0;
            for (;k < data_size; ++k, ++i)
                sv.push_back(i, unsigned(rand()) & 0xFFF);
        } // for i
    }
    else
    {
        auto bi = sv.get_back_inserter();
        for (; i < size; i++)
            bi = unsigned(rand()) & 0xFFF;
    }
}


static
void SparseVectorAccessTest()
{
    std::vector<unsigned> target, target1, target2;
    svect   sv1;

    FillSparseIntervals(sv1);
    BM_DECLARE_TEMP_BLOCK(tb)
    sv1.optimize(tb);

    {
        svect sv2, sv3;
        {
            bm::chrono_taker tt(cout, "sparse_vector random element assignment test", REPEATS / 10);
            for (unsigned i = 0; i < REPEATS / 10; ++i)
            {
                for (unsigned j = 256000; j < 19000000 / 2; ++j)
                {
                    sv2.set(j, 0xF8A);
                }
            }
        }

        {
            bm::chrono_taker tt(cout, "sparse_vector back_inserter test", REPEATS / 10);
            for (unsigned i = 0; i < REPEATS / 10; ++i)
            {
                {
                    sv3.resize(256000);
                    svect::back_insert_iterator bi(sv3.get_back_inserter());
                    for (unsigned j = 256000; j < 19000000 / 2; ++j)
                    {
                        *bi = 0xF8A;
                    }
                }
            }
        }

        // check 
        //
        if (!sv2.equal(sv3))
        {
            std::cerr << "Error! sparse_vector back_insert mismatch."
                << std::endl;
            std::cerr << "sv2.size()=" << sv2.size() << std::endl;
            std::cerr << "sv3.size()=" << sv3.size() << std::endl;

            for (unsigned i = 0; i < sv2.size(); ++i)
            {
                unsigned v2 = sv2[i];
                unsigned v3 = sv3[i];
                if (v2 != v3)
                {
                    std::cerr << "mismatch at: " << i
                        << " v2=" << v2 << " v3=" << v3
                        << std::endl;
                    exit(1);
                }
            }
            exit(1);
        }
    }


    unsigned long long cnt = 0;
    
    unsigned gather_from = 256000;
    unsigned gather_to = 19000000/2;
    std::vector<unsigned> idx;
    for (unsigned j = gather_from; j < gather_to; ++j)
    {
        idx.push_back(j);
    }

    bm::id64_t sum1 = 0;
    {
        bm::chrono_taker tt(cout, "sparse_vector random element access test", REPEATS/10 );
        for (unsigned i = 0; i < REPEATS/10; ++i)
        {
            for (unsigned j = gather_from; j < gather_to; ++j)
                sum1 += sv1[j];
        }
    }

    std::vector<unsigned> target_v;
    target_v.resize(idx.size());
    {
        bm::chrono_taker tt(cout, "sparse_vectot<>::gather() UNSORTED ", REPEATS/5 );
        for (unsigned i = 0; i < REPEATS/5; ++i)
        {
            sv1.gather(target_v.data(), idx.data(), unsigned(idx.size()), bm::BM_UNSORTED);
        }
    }
    bm::id64_t sum2 = 0;
    for (size_t i = 0; i < idx.size(); ++i)
    {
        sum2 += target_v[i];
        auto idx_v = idx[i];
        auto v = sv1[idx_v]; (void) v;
        assert(v == target_v[i]);
    }
    if (sum1 != sum2)
    {
        cout << "\r";
    }

    {
        bm::chrono_taker tt(cout, "sparse_vector<>::decode()", REPEATS / 5);
        auto from = gather_from;
        for (unsigned i = 0; i < REPEATS; ++i)
        {
            auto dsize = sv1.decode(target_v.data(), gather_from, (unsigned)idx.size(), (i == 0));
            from += (dsize % 123);
        }
    }

    target2.resize(idx.size());

    unsigned long long cnt1 = 0, cnt2 = 0;
    {
        bm::chrono_taker tt(cout, "sparse_vector const_iterator test", REPEATS );
        for (unsigned i = 0; i < REPEATS; ++i)
        {
            auto it = sv1.begin();
            auto it_end = sv1.end();
            auto sz = target2.size();
            for (unsigned k = 0; it != it_end && k < sz; ++it, ++k)
            {
                auto v = *it;
                cnt1 += target2[k] = v;
            }
        }
    }

    // check
    //
    size_t sz = min(target1.size(), target2.size());
    for (unsigned j = 0; j < sz; ++j)
    {
        if (target1[j] != target2[j])
        {
            std::cerr << "Error! sparse_vector mismatch at: " << j
            << " t1 = " << target1[j] << " t2 = " << target2[j]
            << " sv[] = " << sv1[j]
            << std::endl;
            exit(1);
        }
    }

    sv1.freeze();
    target_v.resize(idx.size());

    {
        bm::chrono_taker tt(cout, "sparse_vectot<>::gather() UNSORTED (RO) ", REPEATS/5 );
        for (unsigned i = 0; i < REPEATS/5; ++i)
        {
            sv1.gather(target_v.data(), idx.data(), unsigned(idx.size()), bm::BM_UNSORTED);
        }
    }

    {
        bm::chrono_taker tt(cout, "sparse_vector<>::decode() (RO)", REPEATS / 5);
        auto from = gather_from;
        for (unsigned i = 0; i < REPEATS ; ++i)
        {
            auto dsize = sv1.decode(target_v.data(), gather_from, (unsigned)idx.size(), (i == 0));
            from += (dsize % 123);
        }
    }

    {
        bm::chrono_taker tt(cout, "sparse_vector const_iterator test (RO)", REPEATS );
        for (unsigned i = 0; i < REPEATS; ++i)
        {
            auto it = sv1.begin();
            auto it_end = sv1.end();
            auto tsz = target2.size();
            for (unsigned k = 0; it != it_end && k < tsz; ++it, ++k)
            {
                auto v = *it;
                cnt2 += target2[k] = v;
            }
        }
    }

    assert(cnt1 == cnt2);



    
    char buf[256];
    sprintf(buf, "%i", (int)cnt); // to fool some smart compilers like ICC

}

static
void SparseVectorSignedAccessTest()
{
    std::vector<int> target, target1, target2;
    sparse_vector_i32   sv1;

    FillSparseIntervals(sv1);
    BM_DECLARE_TEMP_BLOCK(tb)
    sv1.optimize(tb);

    {
        sparse_vector_i32 sv2, sv3;
        {
            bm::chrono_taker tt(cout, "sparse_vector<int> random element assignment test", REPEATS / 10);
            for (unsigned i = 0; i < REPEATS / 10; ++i)
            {
                for (unsigned j = 256000; j < 19000000 / 2; ++j)
                {
                    sv2.set(j, -0xF8A);
                }
            }
        }

        {
            bm::chrono_taker tt(cout, "sparse_vector<int> back_inserter test", REPEATS / 10);
            for (unsigned i = 0; i < REPEATS / 10; ++i)
            {
                {
                    sv3.resize(256000);
                    sparse_vector_i32::back_insert_iterator bi(sv3.get_back_inserter());
                    for (unsigned j = 256000; j < 19000000 / 2; ++j)
                    {
                        *bi = -0xF8A;
                    }
                }
            }
        }

        // check
        //
        if (!sv2.equal(sv3))
        {
            std::cerr << "Error! sparse_vector<int> back_insert mismatch."
                << std::endl;
            std::cerr << "sv2.size()=" << sv2.size() << std::endl;
            std::cerr << "sv3.size()=" << sv3.size() << std::endl;

            for (unsigned i = 0; i < sv2.size(); ++i)
            {
                const auto v2 = sv2[i];
                const auto v3 = sv3[i];
                if (v2 != v3)
                {
                    std::cerr << "mismatch at: " << i
                        << " v2=" << v2 << " v3=" << v3
                        << std::endl;
                    exit(1);
                }
            }
            exit(1);
        }
    }

    unsigned long long cnt = 0;

    unsigned gather_from = 256000;
    unsigned gather_to = 19000000/2;
    std::vector<unsigned> idx;
    for (unsigned j = gather_from; j < gather_to; ++j)
    {
        idx.push_back(j);
    }

    long long sum1 = 0;
    {
        bm::chrono_taker tt(cout, "sparse_vector<int> random element access test", REPEATS/10 );
        for (unsigned i = 0; i < REPEATS/10; ++i)
        {
            for (unsigned j = gather_from; j < gather_to; ++j)
                sum1 += sv1[j];
        }
    }

    std::vector<int> target_v;
    target_v.resize(idx.size());
    {
        bm::chrono_taker tt(cout, "sparse_vectot<int>::gather() UNSORTED ", REPEATS/5 );
        for (unsigned i = 0; i < REPEATS/10; ++i)
        {
            sv1.gather(target_v.data(), idx.data(), unsigned(idx.size()), bm::BM_UNSORTED);
        }
    }
    long long sum2 = 0;
    for (size_t i = 0; i < idx.size(); ++i)
    {
        sum2 += target_v[i];
        const auto idx_v = idx[i];
        const auto v = sv1[idx_v]; (void)v;
        assert(v == target_v[i]);
    }
    if (sum1 != sum2)
    {
        cout << "\r";
    }

    {
        bm::chrono_taker tt(cout, "sparse_vector<int>::decode()", REPEATS / 5);
        auto from = gather_from;
        for (unsigned i = 0; i < REPEATS / 10; ++i)
        {
            auto dsize = sv1.decode(target_v.data(), gather_from, (unsigned)idx.size(), (i == 0));
            from += (dsize % 123);
        }
    }
    target_v.resize(0);
    target_v.shrink_to_fit();


    {
        bm::chrono_taker tt(cout, "sparse_vector<int> const_iterator test", REPEATS );
        for (unsigned i = 0; i < REPEATS/10; ++i)
        {
            auto it = sv1.begin();
            auto it_end = sv1.end();
            auto sz = target2.size();
            for (unsigned k = 0; it != it_end && k < sz; ++it, ++k)
            {
                auto v = *it;
                target2[k] = v;
            }
        }
    }


    // check
    //
    size_t sz = min(target1.size(), target2.size());
    for (unsigned j = 0; j < sz; ++j)
    {
        if (target1[j] != target2[j])
        {
            std::cerr << "Error! sparse_vector<int> mismatch at: " << j
            << " t1 = " << target1[j] << " t2 = " << target2[j]
            << " sv[] = " << sv1[j]
            << std::endl;
            exit(1);
        }
    }


    char buf[256];
    sprintf(buf, "%i", (int)cnt); // to fool some smart compilers like ICC

}


static
void RSC_SparseVectorFillTest()
{
    bvect bv;// { 10, 20, 100, 200 };

    generate_bvector(bv, 4000000);
    bvect::size_type first, last, mid;
    bv.find_range(first, last);
    mid = first + ((last - first) / 4);


    rsc_sparse_vector_u32 csv1;
    rsc_sparse_vector_u32 csv2(bv);

    {
        bm::chrono_taker tt(cout, "rsc_sparse_vector() set values", REPEATS*1);

        bvect::enumerator en = bv.get_enumerator(mid);
        for (;en.valid(); ++en)
        {
            auto idx = *en;
            csv1.set(idx, 40);
        }
        en.go_to(0);
        for (;en.valid(); ++en)
        {
            auto idx = *en;
            if (idx >= mid)
                break;
            csv1.set(idx, 40);
        }
    }

    {
        bm::chrono_taker tt(cout, "rsc_sparse_vector() set values (rs-index)", REPEATS*1);
        csv2.sync();

        bvect::enumerator en = bv.get_enumerator(mid);
        for (;en.valid(); ++en)
        {
            auto idx = *en;
            csv2.set(idx, 40);
        }

        en.go_to(0);
        for (;en.valid(); ++en)
        {
            auto idx = *en;
            if (idx >= mid)
                break;
            csv2.set(idx, 40);
        }

    }

    bool eq = csv1.equal(csv2);
    if (!eq)
    {
        cerr << "Error: rsc_sparse_vector() set values check failed" << endl;
        assert(0); exit(1);
    }

}

static
void RSC_SparseVectorRandomAccesTest()
{
    const unsigned test_size = 250000000;
    BM_DECLARE_TEMP_BLOCK(tb)

    rsc_sparse_vector_u32   sv1(bm::use_null);
    rsc_sparse_vector_u32   sv2(bm::use_null);


    FillRandomSparseNullVector(sv1, test_size, 100, 153);
    FillRandomSparseNullVector(sv2, test_size, 0, 0);
    sv2.optimize();
    sv2.sync();

    {
        bm::chrono_taker tt(cout, "rsc_sparse_vector<>::sync() (BIT)", REPEATS*10 );
        for (unsigned i = 0; i < REPEATS*10; ++i)
        {
            sv1.sync(true);
        }
    }
    sv1.sync(false);

    //bm::print_svector_stat(cout, sv2);

    std::vector<rsc_sparse_vector_u32::size_type> test_idx;
    std::vector<rsc_sparse_vector_u32::size_type> idx_buf_vec;
    std::vector<rsc_sparse_vector_u32::size_type> test_arr;

    unsigned idx = 0;
    test_idx.push_back(0);
    for (unsigned i = 0; i < 300*65536; ++i)
    {
        idx = i*65535;
        if (idx < sv1.size())
            test_idx.push_back(idx);
        idx = i*65535;
        if (idx < sv1.size())
            test_idx.push_back(idx);
        idx = (unsigned)rand() % test_size;
        if (idx < sv1.size())
            test_idx.push_back(idx);
    }

    idx_buf_vec.resize(test_idx.size());
    test_arr.resize(test_idx.size());

    unsigned long long sum1 = 0, sum2 = 0, sum1d = 0, sum1g = 0, sum1g_ro = 0;

    {
        bm::chrono_taker tt(cout, "rsc_sparse_vector<>::is_null()+get() (DENSE)", REPEATS*10 );
        for (unsigned i = 0; i < test_idx.size(); ++i)
        {
            idx = test_idx[i];
            if (!sv2.is_null(idx))
            {
                unsigned v = sv2.get(idx);
                sum1d += v;
            }
        }
    }

    {
        bm::chrono_taker tt(cout, "rsc_sparse_vector<>::is_null()+get() (BIT)", REPEATS*10 );
        for (unsigned i = 0; i < test_idx.size(); ++i)
        {
            idx = test_idx[i];
            if (!sv1.is_null(idx))
            {
                unsigned v = sv1.get(idx);
                sum1 += v;
            }
        }
    }

    {
        bm::chrono_taker tt(cout, "rsc_sparse_vector<>::try_get_sync() (BIT)", REPEATS*10 );
        for (unsigned i = 0; i < test_idx.size(); ++i)
        {
            idx = test_idx[i];
            unsigned v;
            if (sv1.try_get_sync(idx, v))
                sum2 += v;
        }
    }

    {
        {
        bm::chrono_taker tt(cout, "rsc_sparse_vector<>::gather() (BIT)", REPEATS*10 );
        unsigned sz = (unsigned)test_idx.size();
        sz = sv1.gather(test_arr.data(), test_idx.data(), idx_buf_vec.data(), sz, bm::BM_UNKNOWN);
        }
        // validate
        for (unsigned i = 0; i < test_idx.size(); ++i)
        {
            unsigned v = test_arr[i];
            sum1g += v;
            unsigned v1;
            if (sv1.try_get_sync(test_idx[i], v1))
            {
                assert(v1 == v);
            }
            else
            {
                assert(!v);
            }
        }
    }

    if (sum1 != sum2)
    {
        cerr << "Error! RSC random access check failed!" << endl;
        assert(0);exit(1);
    }
    if (sum1 != sum1g)
    {
        cerr << "Error! RSC random access check failed (gather)!" << endl;
        assert(0);exit(1);
    }

    unsigned long long sum3 = 0, sum4 = 0, sum5 = 0;
    sv1.optimize(tb);

    {
        bm::chrono_taker tt(cout, "rsc_sparse_vector<>::sync() (GAP)", REPEATS*10 );
        for (unsigned i = 0; i < REPEATS*10; ++i)
        {
            sv1.sync(true);
        }
    }

    {
        bm::chrono_taker tt(cout, "rsc_sparse_vector<>::is_null()+get() (GAP)", REPEATS*10 );
        for (unsigned i = 0; i < test_idx.size(); ++i)
        {
            idx = test_idx[i];
            if (!sv1.is_null(idx))
            {
                unsigned v = sv1.get(idx);
                sum3 += v;
            }
        }
    }

    {
        bm::chrono_taker tt(cout, "rsc_sparse_vector<>::try_get_sync() (GAP)", REPEATS*10 );
        for (unsigned i = 0; i < test_idx.size(); ++i)
        {
            idx = test_idx[i];
            unsigned v;
            if (sv1.try_get_sync(idx, v))
                sum4 += v;
        }
    }
    {
        bm::chrono_taker tt(cout, "rsc_sparse_vector<>::try_get() (GAP)", REPEATS * 10);
        for (unsigned i = 0; i < test_idx.size(); ++i)
        {
            idx = test_idx[i];
            unsigned v;
            if (sv1.try_get(idx, v))
                sum5 += v;
        }
    }

    {
        sum1g = 0;
        {
        bm::chrono_taker tt(cout, "rsc_sparse_vector<>::gather() (GAP)", REPEATS*10 );
        unsigned sz = (unsigned)test_idx.size();
        sz = sv1.gather(test_arr.data(), test_idx.data(), idx_buf_vec.data(), sz, bm::BM_UNKNOWN);
        }
        // validate
        for (unsigned i = 0; i < test_idx.size(); ++i)
        {
            unsigned v = test_arr[i];
            sum1g += v;
            unsigned v1;
            if (sv1.try_get_sync(test_idx[i], v1))
            {
                assert(v1 == v);
            }
            else
            {
                assert(!v);
            }
        }
    }


    sv1.freeze();
    {
        sum1g_ro = 0;
        {
        bm::chrono_taker tt(cout, "rsc_sparse_vector<>::gather() (GAP) (RO)", REPEATS*10 );
        unsigned sz = (unsigned)test_idx.size();
        sz = sv1.gather(test_arr.data(), test_idx.data(), idx_buf_vec.data(), sz, bm::BM_UNKNOWN);
        }

        // validate
        for (unsigned i = 0; i < test_idx.size(); ++i)
        {
            unsigned v = test_arr[i];
            sum1g_ro += v;
            unsigned v1;
            if (sv1.try_get_sync(test_idx[i], v1))
            {
                assert(v1 == v);
            }
            else
            {
                assert(!v);
            }
        }
    }

    if (sum3 != sum4 || sum2 != sum3 || sum4 != sum5 || sum5 != sum1g)
    {
        cerr << "Error! RSC random access check failed!" << endl;
        assert(0);exit(1);
    }


}


static
void RSC_SparseVectorAccesTest()
{
    BM_DECLARE_TEMP_BLOCK(tb)

    unsigned test_size = 250000000;
    unsigned dec_size = 65536 * 3;

    std::vector<unsigned> vect, vect1, vect2;
    vect.resize(dec_size);
    vect1.resize(dec_size);
    vect2.resize(dec_size);


    sparse_vector_u32   sv1(bm::use_null);
    rsc_sparse_vector_u32::size_type sz, sz1;
    (void)sz1;

    FillSparseNullVector(sv1, test_size, 2, 150);
    {
        rsc_sparse_vector_u32 csv1;
        csv1.load_from(sv1);
        csv1.optimize(tb);
        sv1.clear();

        {
            bm::chrono_taker tt(cout, "rsc_sparse_vector()::decode() test (sparse)", REPEATS*10 );
            unsigned from = 0;
            for (unsigned i = 0; i < REPEATS*10; ++i)
            {
                sz = csv1.decode(&vect[0], from, dec_size);
                assert(sz); (void)sz;
                from += 3;//rand() % dec_size;
                if (from > csv1.size())
                    from = 0;
            } // for
        }
        {
            bm::chrono_taker tt(cout, "rsc_sparse_vector()::decode_buf() test (sparse)", REPEATS*10 );
            unsigned from = 0;
            for (unsigned i = 0; i < REPEATS*10; ++i)
            {
                sz1 = csv1.decode_buf(&vect1[0], &vect2[0], from, dec_size);
                assert(sz); (void)sz;
                from += 3;//rand() % dec_size;
                if (from > csv1.size())
                    from = 0;
            } // for
        }
        assert (sz == sz1);
        for (unsigned i = 0; i < sz; ++i)
        {
            auto v = vect[i];
            auto v1 = vect1[i];
            if (v != v1)
            {
                std::cerr << "decode integrity check failed" << std::endl;
                assert(0);  exit(1);
            }
        }

    }

    FillSparseNullVector(sv1, test_size, 64, 5);
    {
        rsc_sparse_vector_u32 csv1;
        csv1.load_from(sv1);
        csv1.optimize(tb);
        sv1.clear();

        {
            bm::chrono_taker tt(cout, "rsc_sparse_vector<>::decode() test (dense)", REPEATS*10 );
            unsigned from = 0;
            for (unsigned i = 0; i < REPEATS*10; ++i)
            {
                sz = csv1.decode(&vect[0], from, dec_size);

                assert(sz);(void)sz;
                from += 7; // rand() % dec_size;
                if (from > csv1.size())
                    from = 0;
            } // for
        }
        {
            bm::chrono_taker tt(cout, "rsc_sparse_vector()::decode_buf() test (dense)", REPEATS*10 );
            unsigned from = 0;
            for (unsigned i = 0; i < REPEATS*10; ++i)
            {
                sz1 = csv1.decode_buf(&vect1[0], &vect2[0], from, dec_size);
                assert(sz); (void)sz;
                from += 7; // rand() % dec_size;
                if (from > csv1.size())
                    from = 0;
            } // for
        }
        assert (sz == sz1);
        for (unsigned i = 0; i < sz; ++i)
        {
            assert(vect[i] == vect1[i]);
        }

    }



}




static
void OptimizeTest()
{
    std::vector<bvect> bv_coll;
    std::vector<bvect> bv_coll_control;

    GenerateTestCollection(&bv_coll, 100, 80000000, false);
    bv_coll_control = bv_coll;
    
    BM_DECLARE_TEMP_BLOCK(tb)
    {
        bm::chrono_taker tt(cout, "bvector<>::optimize() ", 1);
        for (unsigned k = 0; k < bv_coll.size(); ++k)
        {
            bv_coll[k].optimize(tb);
        }
    }
    for (unsigned k = 0; k < bv_coll.size(); ++k)
    {
        const bvect& bv1 = bv_coll[k];
        const bvect& bv2 = bv_coll_control[k];
        int cmp = bv1.compare(bv2);
        if (cmp != 0)
        {
            std::cerr << "Optimization error!" << endl;
            exit(1);
        }
    }

    
}

static
void AggregatorTest()
{
    int res;
    bvect* bv_arr[128] = { 0, };
    bvect* bv_arr2[128] = { 0, };
    bm::aggregator<bvect> agg;
    
    std::vector<bvect> bv_coll;
    GenerateTestCollection(&bv_coll, 25, 80000000);
    std::vector<bvect> bv_coll2;
    GenerateTestCollection(&bv_coll2, 7, 50000000);

    if (!bv_coll.size())
        return;
    
    std::vector<bvect>& bvc = bv_coll;
    for (unsigned k = 0; k < bv_coll.size(); ++k)
    {
        bv_arr[k] = &(bvc[k]);
    } // for
    std::vector<bvect>& bvc2 = bv_coll2;
    for (unsigned k = 0; k < bv_coll2.size(); ++k)
    {
        bv_arr2[k] = &(bvc2[k]);
    } // for

    std::unique_ptr<bvect> bv_target1(new bvect);
    std::unique_ptr<bvect> bv_target2(new bvect);

    {
    bm::chrono_taker tt(cout, "Horizontal aggregator OR ", REPEATS);
    for (unsigned i = 0; i < REPEATS; ++i)
    {
        agg.combine_or_horizontal(*bv_target1, bv_arr, unsigned(bv_coll.size()));
    } // for
    }

    {
    bm::chrono_taker tt(cout, "aggregator OR", REPEATS);
    for (unsigned i = 0; i < REPEATS; ++i)
    {
        agg.combine_or(*bv_target2, bv_arr, unsigned(bv_coll.size()));
    } // for
    }

    res = bv_target1->compare(*bv_target2);
    if (res != 0)
    {
        std::cerr << "Error: Aggregator OR integrity failed." << std::endl;
        exit(1);
    }


    // ------------------------------------------------------------------
    {
    bm::chrono_taker tt(cout, "Horizontal aggregator AND", REPEATS);
    for (unsigned i = 0; i < REPEATS; ++i)
    {
        agg.combine_and_horizontal(*bv_target1, bv_arr, unsigned(bv_coll.size()));
    } // for
    }

    {
    bm::chrono_taker tt(cout, "aggregator AND", REPEATS);
    for (unsigned i = 0; i < REPEATS; ++i)
    {
        agg.combine_and(*bv_target2, bv_arr, unsigned(bv_coll.size()));
    } // for
    }
    auto and_cnt = bv_target1->count();

    res = bv_target1->compare(*bv_target2);
    if (res != 0)
    {
        std::cerr << "Error: Aggregator AND integrity failed." << std::endl;
        exit(1);
    }

    {
    bm::chrono_taker tt(cout, "Horizontal aggregator AND-SUB", REPEATS);
    for (unsigned i = 0; i < REPEATS; ++i)
    {
        agg.combine_and_sub_horizontal(*bv_target1,
                                        bv_arr, unsigned(bv_coll.size()),
                                        bv_arr2, unsigned(bv_coll2.size())
                                        );
    } // for
    }
    
    auto and_sub_cnt = bv_target1->count();
    if (and_sub_cnt > and_cnt)
    {
        std::cerr << "Error: Aggregator count AND-SUB integrity failed." << std::endl;
        exit(1);
    }

    {
    bm::chrono_taker tt(cout, "aggregator AND-SUB", REPEATS);
    for (unsigned i = 0; i < REPEATS; ++i)
    {
        agg.combine_and_sub(*bv_target2,
                            bv_arr, unsigned(bv_coll.size()),
                            bv_arr2, unsigned(bv_coll2.size()),
                            false
                            );
    } // for
    }
    
    res = bv_target1->compare(*bv_target2);
    if (res != 0)
    {
        std::cerr << "Error: Aggregator AND-SUB integrity failed." << std::endl;
        exit(1);
    }


    for (size_t i = 0; i < bv_coll.size(); ++i)
        bv_arr[i]->freeze();
    for (size_t i = 0; i < bv_coll2.size(); ++i)
        bv_arr2[i]->freeze();


    {
    bm::chrono_taker tt(cout, "aggregator AND-SUB (RO)", REPEATS);
    for (unsigned i = 0; i < REPEATS; ++i)
    {
        agg.combine_and_sub(*bv_target2,
                            bv_arr, unsigned(bv_coll.size()),
                            bv_arr2, unsigned(bv_coll2.size()),
                            false
                            );
    } // for
    }
    res = bv_target1->compare(*bv_target2);
    if (res != 0)
    {
        std::cerr << "Error: Aggregator AND-SUB integrity failed. (RO)" << std::endl;
        exit(1);
    }

    //std::cout << bv_target1->count() << std::endl;
}


static
void BvectorShiftTest()
{

    {
        std::vector<bvect> bv_coll;
        GenerateTestCollection(&bv_coll, 25, 80000000);

        if (!bv_coll.size())
            return;

        {
            bm::chrono_taker tt(cout, "bvector<>::shift_right() ", REPEATS);
            for (unsigned i = 0; i < REPEATS; ++i)
            {
                for (unsigned k = 0; k < bv_coll.size(); ++k)
                {
                    bv_coll[k].shift_right();
                } // for
            } // for
        }
    }

    {
        std::vector<bvect> bv_coll;
        GenerateTestCollection(&bv_coll, 25, 80000000);

        if (!bv_coll.size())
            return;

        {
            bm::chrono_taker tt(cout, "bvector<>::shift_left() ", REPEATS);
            for (unsigned i = 0; i < REPEATS; ++i)
            {
                for (unsigned k = 0; k < bv_coll.size(); ++k)
                {
                    bv_coll[k].shift_left();
                } // for
            } // for
        }
    }


    bvect mask_bv; // mask vector
    mask_bv.init();
    generate_bvector(mask_bv, 75000000, false); // mask is shorter on both ends

    std::vector<bvect> bv_coll1;
    GenerateTestCollection(&bv_coll1, 25, 80000000);
    
    std::vector<bvect> bv_coll2;
    GenerateTestCollection(&bv_coll2, 25, 80000000);
    

    {
        {
            bm::chrono_taker tt(cout, "bvector<>::shift_right()+AND ", REPEATS);
            for (unsigned i = 0; i < REPEATS; ++i)
            {
                bvect bv(mask_bv);
                for (unsigned k = 0; k < bv_coll1.size(); ++k)
                {
                    bv.shift_right();
                    bv &= bv_coll1[k];
                } // for
            } // for
        }
    }

    {
        bm::aggregator<bvect> agg;
        
        {
            bm::chrono_taker tt(cout, "aggregator::shift_right_and() ", REPEATS);
            agg.add(&mask_bv);
            for (unsigned k = 0; k < bv_coll2.size(); ++k)
            {
                agg.add(&bv_coll2[k]);
            }
            for (unsigned i = 0; i < REPEATS; ++i)
            {
                bvect bv;
                agg.combine_shift_right_and(bv);
            }
        }
    }
/*
    // check correcness
    //
    for (unsigned k = 0; k < bv_coll1.size(); ++k)
    {
        const bvect& bv1 = bv_coll1[k];
        const bvect& bv2 = bv_coll2[k];
        auto cmp = bv1.compare(bv2);
        if (cmp != 0)
        {
            cerr << "Mismatch error!" << endl;
            exit(1);
        }
    } // for
*/

    //std::cout << bv_target1->count() << std::endl;
}



static
void Set2SetTransformTest()
{
    svect   sv;
    
    FillSparseIntervals(sv);
    BM_DECLARE_TEMP_BLOCK(tb)
    sv.optimize(tb);
    
    bvect bv_values;
    bvect bv_non_values;
    bvect bv_sample;
    bvect bv_non_sample;

    {
        bm::sparse_vector_scanner<svect> scanner;
        scanner.find_nonzero(sv, bv_values);
        scanner.find_zero(sv, bv_non_values);
    
        unsigned non_v_count = bv_non_values.count() / 20;
        if (non_v_count > 3000000)
            non_v_count = 3000000;

        if (non_v_count)
        {
            bm::random_subset<bvect> rand_sampler;
            rand_sampler.sample(bv_sample, bv_values, 6000000);
            rand_sampler.sample(bv_non_sample, bv_non_values, non_v_count);
            //cout << "zero mix = " << bv_non_sample.count() << endl;
            bv_sample |= bv_non_sample; // add some missing values
        }
    }
    
    bm::set2set_11_transform<svect> set2set;
    int cnt = 0;
    {
        bm::chrono_taker tt(cout, "set2set_11_transform::run()", REPEATS/10);
        for (unsigned i = 0; i < REPEATS/10; ++i)
        {
            bvect bv_out;
            set2set.run(bv_sample, sv, bv_out);
            cnt += bv_out.any();
        }
    }

    bv_sample.freeze();
    sv.freeze();

    int cnt2 = 0;
    {
        bm::chrono_taker tt(cout, "set2set_11_transform::run() (RO)", REPEATS/10);
        for (unsigned i = 0; i < REPEATS/10; ++i)
        {
            bvect bv_out;
            set2set.run(bv_sample, sv, bv_out);
            cnt2 += bv_out.any();
        }
    }

    assert(cnt == cnt2);


    /*
    {
    TimeTaker tt("set2set_11_transform::one_pass_run", REPEATS/10);

        for (unsigned i = 0; i < REPEATS/10; ++i)
        {
            bvect bv_out;
            set2set.one_pass_run(bv_sample, sv, bv_out);

            cnt += bv_out.any();
        }
    }
    */
    char buf[256];
    sprintf(buf, "%i", (int)cnt); // to fool some smart compilers like ICC
}

static
void RangeCopyTest()
{
    const unsigned vect_max = BSIZE;
    bvect bv;
    generate_bvector(bv, vect_max);

    {
        bm::chrono_taker tt(cout, "bvector<>::copy_range()", REPEATS * 25);
        for (unsigned i = 0; i < REPEATS * 25; ++i)
        {
            unsigned from = vect_max / 4;
            from = from * (unsigned)(rand() % 3);
            unsigned to = from + 1 + (unsigned)(65536 * rand() % 5);
            bvect bv_cp;
            bv_cp.copy_range(bv, from, to);
        } // for
    }
    {
        bm::chrono_taker tt(cout, "bvector<>:: copy range constructor", REPEATS * 25);
        for (unsigned i = 0; i < REPEATS * 25; ++i)
        {
            unsigned from = vect_max / 4;
            from = from * (unsigned)(rand() % 3);
            unsigned to = from + 1 + (unsigned)(65536 * rand() % 5);
            bvect bv_cp(bv, from, to);
        } // for
    }

    {
        bm::chrono_taker tt(cout, "copy range with AND", REPEATS * 25);
        for (unsigned i = 0; i < REPEATS * 25; ++i)
        {
            unsigned from = vect_max / 4;
            from = from * (unsigned)(rand() % 3);
            unsigned to = from + 1 + (unsigned)(65536 * rand() % 5);
            bvect bv_cp;
            bv_cp.set_range(from, to);
            bv_cp &= bv;
        } // for
    }
}

/// Reference (naive) interval detector based on population counting 
/// and boundaries tests
///
template<typename BV>
bool test_interval(const BV& bv,
        typename BV::size_type left, typename BV::size_type right) noexcept
{
    if (left > right)
        bm::xor_swap(left, right); // make sure left <= right
    bool is_left(0), is_right(0);
    if (left) // check left-1 bit (if exists)
        is_left = bv.test(left - 1);
    if ((is_left == false) && (right < bm::id_max - 1))
        is_right = bv.test(right + 1); // check [...right] range condition
    if (is_left == false && is_right == false)
    {
        typename BV::size_type cnt = bv.count_range(left, right);
        if (cnt == (1 + right - left))
            return true;
    }
    return false;
}

static
void IntervalsTest()
{
    const unsigned vect_max = BSIZE * 2;
    bvect bv, bv_inv;

    // generate the test vector
    {
        bool b;
        bvect::size_type istart(0), ilen(0);
        bvect::bulk_insert_iterator iit(bv);
        for (istart = 0; istart < vect_max; )
        {
            for (bvect::size_type i = istart; i <= (istart+ilen); ++i)
            {
                iit = i;
            } // for i
            iit.flush();

            ilen += 1;
            b = bv.test(istart + ilen); (void)b;
            assert(!b);
            istart += (ilen + 2);
            if (ilen > 1024)
                ilen = 0;
            
        } // for istart
    }
    bv_inv = bv;
    bv_inv.invert();
    //auto cnt = bm::count_intervals(bv);
    //cout << cnt << endl;

    const char* msg = "bvector<>::is_all_one_range() (BITS)";
    const char* msg2 = "bvector<>::any_range() (BITS)";
    const char* msg3 = "bvector<>::count_range() (BITS)";
    const char* msg4 = "bvector<>::is_interval() (BITS)";
    const char* msg5 = "reference_is_interval() (BITS)";
    const char* msg6 = "bvector<>::find_interval_start() (BITS)";
    const char* msg7 = "bvector<>::find_interval_end() (BITS)";
    const char* msg8 = "interval_enumerator<> (BITS)";
    const char* msg9 = "bvector<>::enumerator (BITS)";

    for (unsigned pass = 0; pass < 2; ++pass)
    {
        {
            bm::chrono_taker tt(cout, msg3, REPEATS * 1);
            bvect::size_type istart(0), ilen(0);
            for (istart = 0; istart < vect_max; )
            {
                for (bvect::size_type i = istart; i <= (istart + ilen); ++i)
                {
                    auto cnt = bv.count_range(istart, i);
                    if (!cnt)
                    {
                        cerr << "Errro: count_range test failed! (1)" << endl;
                        assert(0); exit(1);
                    }
                    auto icnt = bv_inv.count_range(istart, i);
                    if (icnt)
                    {
                        cerr << "Errro: count_range test failed! (2)" << endl;
                        assert(0); exit(1);
                    }
                } // for i
                ilen += 1;
                istart += (ilen + 2);
                if (ilen > 1024)
                    ilen = 0;
            } // for istart
        }

        {
            bm::chrono_taker tt(cout, msg, REPEATS * 1);
            bvect::size_type istart(0), ilen(0);
            for (istart = 0; istart < vect_max; )
            {
                for (bvect::size_type i = istart; i <= (istart+ilen); ++i)
                {
                    bool all_one = bv.is_all_one_range(istart, i);
                    if (!all_one)
                    {
                        cerr << "Errro: is_all_one_range test failed! (1)" << endl;
                        assert(0); exit(1);
                    }
                    all_one = bv_inv.is_all_one_range(istart, i);
                    if (all_one)
                    {
                        cerr << "Errro: is_all_one_range test failed! (2)" << endl;
                        assert(0); exit(1);
                    }
                } // for i
                ilen += 1;
                istart += (ilen + 2);
                if (ilen > 1024)
                    ilen = 0;
            } // for istart
        }

        {
            bm::chrono_taker tt(cout, msg4, REPEATS * 1);
            bvect::size_type istart(0), ilen(0);
            for (istart = 0; istart < vect_max; )
            {
                for (bvect::size_type i = istart; i <= (istart + ilen); ++i)
                {
                    bool is_int = bm::is_interval(bv, istart, i);
                    if (!is_int)
                    {
                        if (i == istart + ilen)
                        {
                            cerr << "Error: is_interval test failed! (1)" << endl;
                            is_int = bm::is_interval(bv, istart, i);
                            assert(0); exit(1);
                        }
                    }
                    else
                    {
                        assert(i == istart + ilen);
                    }

                    is_int = bm::is_interval(bv_inv, istart, i);
                    if (is_int)
                    {
                        cerr << "Errro: is_interval test failed! (2)" << endl;
                        assert(0); exit(1);
                    }
                } // for i
                ilen += 1;
                istart += (ilen + 2);
                if (ilen > 1024)
                    ilen = 0;
            } // for istart
        }

        {
            bm::chrono_taker tt(cout, msg5, REPEATS * 1);
            bvect::size_type istart(0), ilen(0);
            for (istart = 0; istart < vect_max; )
            {
                for (bvect::size_type i = istart; i <= (istart + ilen); ++i)
                {
                    bool is_int = test_interval(bv, istart, i);
                    if (!is_int)
                    {
                        if (i == istart + ilen)
                        {
                            cerr << "Error: test_interval test failed! (1)" << endl;
                            assert(0); exit(1);
                        }
                    }
                    else
                    {
                        assert(i == istart + ilen);
                    }
                    is_int = test_interval(bv_inv, istart, i);
                    if (is_int)
                    {
                        cerr << "Errro: test_interval test failed! (2)" << endl;
                        assert(0); exit(1);
                    }
                } // for i
                ilen += 1;
                istart += (ilen + 2);
                if (ilen > 1024)
                    ilen = 0;
            } // for istart
        }

        {
            bm::chrono_taker tt(cout, msg2, REPEATS * 1);
            bvect::size_type istart(0), ilen(0);
            for (istart = 0; istart < vect_max; )
            {
                for (bvect::size_type i = istart; i <= (istart+ilen); ++i)
                {
                    bool any_one = bv.any_range(istart, i);
                    if (!any_one)
                    {
                        cerr << "Errro: any_range test failed! (1)" << endl;
                        assert(0); exit(1);
                    }
                    any_one = bv_inv.any_range(istart, i);
                    if (any_one)
                    {
                        cerr << "Errro: any_range test failed! (2)" << endl;
                        assert(0); exit(1);
                    }

                } // for i
                ilen += 1;
                istart += (ilen + 2);
                if (ilen > 1024)
                    ilen = 0;
            } // for istart
        }

        {
            bm::chrono_taker tt(cout, msg6, REPEATS * 1);
            bvect::size_type istart(0), ilen(0);
            for (istart = 0; istart < vect_max; )
            {
                for (bvect::size_type i = istart; i <= (istart + ilen); ++i)
                {
                    bvect::size_type pos;
                    auto diff = i - istart;
                    if (!diff)
                        continue;
                    bool b = bm::find_interval_start(bv, istart+diff, pos);
                    assert(b); (void)b;
                    assert(pos == istart);
                } // for i
                ilen += 1;
                istart += (ilen + 2);
                if (ilen > 1024)
                    ilen = 0;
            } // for istart
        }

        {
            bm::chrono_taker tt(cout, msg7, REPEATS * 1);
            bvect::size_type istart(0), ilen(0);
            for (istart = 0; istart < vect_max; )
            {
                for (bvect::size_type i = istart; i <= (istart + ilen); ++i)
                {
                    bvect::size_type pos;
                    auto diff = i - istart;
                    if (!diff)
                        continue;
                    bool b = bm::find_interval_end(bv, istart+diff, pos);
                    assert(b); (void)b;
                    assert(pos == istart+ilen);
                } // for i
                ilen += 1;
                istart += (ilen + 2);
                if (ilen > 1024)
                    ilen = 0;
            } // for istart
        }

        bvect::size_type cnt_c = bv.count();
        {
            bm::chrono_taker tt(cout, msg8, REPEATS * 10);
   
            for (unsigned i = 0; i < REPEATS * 10; ++i)
            {
                bvect::size_type cnt = 0;
                bm::interval_enumerator<bvect> ien(bv);
                if (ien.valid())
                {
                    do {
                        auto s = ien.start();
                        auto e = ien.end();
                        cnt += e - s + 1;
                    } while (ien.advance());
                }
                assert(cnt == cnt_c);
                if (cnt != cnt_c)
                {
                    cerr << "Count mismatch!" << endl;
                    exit(1);
                }
            }
        }
        {
            bvect::size_type sum = 0;
            bm::chrono_taker tt(cout, msg9, REPEATS/4);

            for (unsigned i = 0; i < REPEATS/4; ++i)
            {
                bvect::size_type cnt = 0;
                bvect::enumerator en = bv.get_enumerator(0);
                while (en.valid())
                {
                    auto v = *en; 
                    sum += v;
                    ++cnt;
                    ++en;
                }
                assert(cnt == cnt_c);
            }
            char buf[256];
            sprintf(buf, "%u", sum); // this is to prevent unwanted optimizations by some compilers

        }

        bv.optimize();
        bv_inv.optimize();

        msg = "bvector<>::is_all_one_range() (BITS+GAPS)";
        msg2 = "bvector<>::any_range() (BITS+GAPS)";
        msg3 = "bvector<>::count_range() (BITS+GAPS)";
        msg4 = "bvector<>::is_interval() (BITS+GAPS)";
        msg5 = "reference_is_interval() (BITS+GAPS)";
        msg6 = "bvector<>::find_interval_start() (BITS+GAPS)";
        msg7 = "bvector<>::find_interval_end() (BITS+GAPS)";
        msg8 = "interval_enumerator<> (BITS+GAPS)";
        msg9 = "bvector<>::enumerator (BITS+GAPS)";
    } // for pass

}

static
void RankCompressionTest()
{
    bvect bv_i1, bv_s1, bv11, bv12, bv21, bv22;
    bvect bv_i2, bv_s2;
    bvect bv11_s, bv12_s, bv21_s, bv22_s;

    generate_bvector(bv_i1);
    generate_bvector(bv_i2);
    bv_i2.optimize();

    bm::random_subset<bvect> rsub;
    rsub.sample(bv_s1, bv_i1, 1000);
    rsub.sample(bv_s2, bv_i2, 1000);
    bv_s2.optimize();

    bm::rank_compressor<bvect> rc;

    std::unique_ptr<bvect::rs_index_type> bc1(new bvect::rs_index_type());
    bv_i1.build_rs_index(bc1.get());
    std::unique_ptr<bvect::rs_index_type> bc2(new bvect::rs_index_type());
    bv_i2.build_rs_index(bc2.get());

    {
        bm::chrono_taker tt(cout, "Rank compression test", REPEATS * 10);
        for (unsigned i = 0; i < REPEATS * 10; ++i)
        {
            rc.compress(bv11, bv_i1, bv_s1);
            rc.compress(bv12, bv_i2, bv_s2);
        } // for
    }
    {
        bm::chrono_taker tt(cout, "Rank compression (by source) test", REPEATS * 10);
        for (unsigned i = 0; i < REPEATS * 10; ++i)
        {
            rc.compress_by_source(bv21, bv_i1, *bc1, bv_s1);
            rc.compress_by_source(bv22, bv_i2, *bc2, bv_s2);
        } // for
    }
    
    {
        bm::chrono_taker tt(cout, "Rank decompression test", REPEATS * 10);
        for (unsigned i = 0; i < REPEATS * 10; ++i)
        {
            rc.decompress(bv11_s, bv_i1, bv11);
            rc.decompress(bv12_s, bv_i2, bv12);
            rc.decompress(bv21_s, bv_i1, bv21);
            rc.decompress(bv22_s, bv_i2, bv22);
        } // for
    }

    int cmp;
    cmp = bv11_s.compare(bv_s1);
    if (cmp != 0)
    {
        std::cerr << "1. Rank check failed!" << std::endl;
        exit(1);
    }
    cmp = bv12_s.compare(bv_s2);
    if (cmp != 0)
    {
        std::cerr << "2. Rank check failed!" << std::endl;
        exit(1);
    }
    cmp = bv21_s.compare(bv_s1);
    if (cmp != 0)
    {
        std::cerr << "3. Rank check failed!" << std::endl;
        exit(1);
    }
    cmp = bv22_s.compare(bv_s2);
    if (cmp != 0)
    {
        std::cerr << "4. Rank check failed!" << std::endl;
        exit(1);
    }
    
}

static
void generate_serialization_test_set(sparse_vector_u32&   sv,
                                     unsigned vector_max = BSIZE)
{
    sparse_vector_u32::back_insert_iterator bi(sv.get_back_inserter());

    unsigned v = 0;
    for (unsigned i = 0; i < vector_max; ++i)
    {
        unsigned plato = (unsigned)rand() % 16;
        for (unsigned j = 0; i < vector_max && j < plato; ++i, ++j)
        {
            *bi = v;
        } // for j
        if (++v > 100000)
            v = 0;
        unsigned nulls = (unsigned)rand() % 16;
        if (nulls)
            bi.add_null(nulls);
        i += nulls;
    } // for i

    sv.optimize();
}


template<typename VECT, typename SVECT>
void generate_scanner_test_set(VECT&          vect,
                               bvect&         bv_null,
                               SVECT&         sv,
                               unsigned vector_max = BSIZE)
{
    auto bi(sv.get_back_inserter());

    vect.resize(vector_max);
    bv_null.reset();

    for (unsigned i = 0; i < vector_max; ++i)
    {
        unsigned vg = unsigned(rand_dis(gen));
        typename VECT::value_type v;
        if constexpr (std::is_signed<typename VECT::value_type>::value)
        {
            v = (typename VECT::value_type)vg;
            if (v & 1)
                v = -v;
        }
        else
            v = vg;

        vect[i] = v;
        bv_null[i] = true; // not NULL(assigned) element
        *bi = v; // push back an element to sparse vector
        if (i % 64 == 0)
        {
            bi.add_null(5);  // add 5 unassigned elements using back inserter
            i += 5;  // insert a small NULL plate (unassigned values)
        }
    } // for
    sv.optimize();
}

template<typename VECT>
void vector_search(const VECT&                vect,
                   const bvect&               bv_null,
                   typename VECT::value_type  value,
                   bvect&                     bv_res)
{
    bv_res.init(); // always use init() if set_bit_no_check()
    for (size_t i = 0; i < vect.size(); ++i)
    {
        if (vect[i] == value)
            bv_res.set_bit_no_check((bm::id_t)i);
    } // for
    bv_res &= bv_null; // correct results to only include non-NULL values
}

template<typename VECT>
void vector_search_GT(const VECT& vect,
                      const bvect&                 bv_null,
                      typename VECT::value_type    value,
                      bvect&                       bv_res)
{
    bv_res.init(); // always use init() if set_bit_no_check()
    auto sz = vect.size();
    for (size_t i = 0; i < sz; ++i)
    {
        if (vect[i] > value)
            bv_res.set_bit_no_check((bm::id_t)i);
    } // for
    bv_res &= bv_null; // correct results to only include non-NULL values
}


template<typename VECT>
void vector_search_GT_sorted(const VECT& vect,
                             const bvect&                 bv_null,
                             typename VECT::value_type    value,
                             bvect&                       bv_res)
{
    bv_res.init(); // always use init() if set_bit_no_check()
    auto it = std::lower_bound(vect.begin(), vect.end(), value);

    if (it == vect.end())
        return;

    for (; it < vect.end(); ++it)
    {
        auto v = *it;
        if (v > value)
        {
            auto idx_from = it - vect.begin();
            bv_res.set_range(bvect::size_type(idx_from), bvect::size_type(vect.size()-1));
            break;
        }
    }
    bv_res &= bv_null; // correct results to only include non-NULL values
}




template<typename VECT>
void generate_search_samples(VECT& search_vect,
                             size_t search_size,
                             const VECT& src_vect)
{
    assert(src_vect.size() >= search_size);
    search_vect.reserve(search_size);
    size_t step = src_vect.size() / search_size;
    for (size_t i = 0; i < search_size; i+= step)
    {
        auto v = src_vect[i];
        search_vect.push_back(v);
    } // for
}


static
void SparseVectorScannerCmpTest()
{
    std::vector<unsigned> vect;
    bvect bv_null;
    sparse_vector_u32 sv(bm::use_null);

    generate_scanner_test_set(vect, bv_null, sv, BSIZE/4);

    // generate a search vector for benchmarking
    std::vector<unsigned> search_vect;
    generate_search_samples(search_vect, REPEATS, vect);

    bm::sparse_vector_scanner<sparse_vector_u32> scanner;

    bvect bv_res1, bv_res2, bv_res3, bv_res4;

    unsigned search_repeats = REPEATS;

    // run GT validation
    {
        bm::chrono_taker tt(cout, "GT validation", search_repeats);
        for (unsigned i = 0; i < search_repeats; ++i)
        {
            unsigned vs = search_vect[i];
            bvect bv1, bv2;
            vector_search_GT(vect, bv_null, vs, bv1);

            scanner.find_gt_horizontal(sv, vs, bv2, true /*NULL correct*/);
            bool eq = bv1.equal(bv2);
            if (!eq)
            {
                cerr << "ERROR: Failed GT integrity checks!" << endl;
                assert(0); exit(1);
            }
        } // for
    }


    {
        bm::chrono_taker tt(cout, "std::vector<> GT scan ", search_repeats);
        for (unsigned i = 0; i < search_repeats; ++i)
        {
            unsigned vs = search_vect[i];
            vector_search_GT(vect, bv_null, vs, bv_res1);
            bv_res1.clear();
        } // for
    }


    std::sort(vect.begin(), vect.end());
    sparse_vector_u32 sv2;
    {
        auto bi(sv2.get_back_inserter());
        for (size_t i=0; i < vect.size(); ++i)
            bi = vect[i];
        bi.flush();
    }
    sv2.optimize();

    {
        bvect bv_all;
        bv_all.set_range(0, sv2.size()-1);
        bm::chrono_taker tt(cout, "GT validation (sorted)", search_repeats);
        for (unsigned i = 0; i < search_repeats; ++i)
        {
            unsigned vs = search_vect[i];
            bvect bv1, bv2, bv3;
            vector_search_GT(vect, bv_all, vs, bv1);
            vector_search_GT_sorted(vect, bv_all, vs, bv2);

            scanner.find_gt_horizontal(sv2, vs, bv3, false /*no NULL correct*/);
            bool eq = bv1.equal(bv2);
            if (!eq)
            {
                cerr << "ERROR: Failed GT (sorted) integrity checks (1)!" << endl;
                assert(0); exit(1);
            }
            eq = bv1.equal(bv3);
            if (!eq)
            {
                cerr << "ERROR: Failed GT (sorted) integrity checks (2)!" << endl;
                assert(0); exit(1);
            }
        } // for
    }

    {
        bm::chrono_taker tt(cout, "std::vector<> GT scan (lower_bound) ", search_repeats);
        for (unsigned i = 0; i < search_repeats; ++i)
        {
            unsigned vs = search_vect[i];
            vector_search_GT_sorted(vect, bv_null, vs, bv_res3);
        } // for
    }


    vect.resize(0);
    vect.shrink_to_fit();

    {
        bm::chrono_taker tt(cout, "horizontal sparse vector scanner find_GT()", search_repeats);
        for (unsigned i = 0; i < search_repeats; ++i)
        {
            unsigned vs = search_vect[i];
            scanner.find_gt_horizontal(sv, vs, bv_res2, true /*NULL correct*/);
        } // for
    }

    {
        bm::chrono_taker tt(cout, "horizontal sparse vector scanner find_GT() (sorted)", search_repeats);
        for (unsigned i = 0; i < search_repeats; ++i)
        {
            unsigned vs = search_vect[i];
            scanner.find_gt_horizontal(sv2, vs, bv_res4, true /*NULL correct*/);
        } // for
    }


}

static
void SparseVectorScannerSignedCmpTest()
{
    std::vector<int> vect;
    bvect bv_null;
    sparse_vector_i32 sv(bm::use_null);

    generate_scanner_test_set(vect, bv_null, sv, BSIZE/4);

    // generate a search vector for benchmarking
    std::vector<int> search_vect;
    generate_search_samples(search_vect, REPEATS, vect);

    bm::sparse_vector_scanner<sparse_vector_i32> scanner;

    bvect bv_res1, bv_res2, bv_res3, bv_res4;

    unsigned search_repeats = REPEATS;

    // run GT validation
    {
        bm::chrono_taker tt(cout, "GT (signed)validation", search_repeats);
        for (unsigned i = 0; i < search_repeats; ++i)
        {
            auto vs = search_vect[i];
            bvect bv1, bv2;
            vector_search_GT(vect, bv_null, vs, bv1);

            scanner.find_gt_horizontal(sv, vs, bv2, true /*NULL correct*/);
            bool eq = bv1.equal(bv2);
            if (!eq)
            {
                cerr << "ERROR: Failed GT integrity checks!" << endl;
                assert(0); exit(1);
            }
        } // for
    }


    {
        bm::chrono_taker tt(cout, "std::vector<> GT (signed) scan ", search_repeats);
        for (unsigned i = 0; i < search_repeats; ++i)
        {
            auto vs = search_vect[i];
            vector_search_GT(vect, bv_null, vs, bv_res1);
            bv_res1.clear();
        } // for
    }

    std::sort(vect.begin(), vect.end());
    sparse_vector_i32 sv2;
    {
        auto bi(sv2.get_back_inserter());
        for (size_t i=0; i < vect.size(); ++i)
            bi = vect[i];
        bi.flush();
    }
    sv2.optimize();

    {
        bvect bv_all;
        bv_all.set_range(0, sv2.size()-1);
        bm::chrono_taker tt(cout, "GT validation (signed)(sorted)", search_repeats);
        for (unsigned i = 0; i < search_repeats; ++i)
        {
            auto vs = search_vect[i];
            bvect bv1, bv2, bv3;
            vector_search_GT(vect, bv_all, vs, bv1);
            vector_search_GT_sorted(vect, bv_all, vs, bv2);

            scanner.find_gt_horizontal(sv2, vs, bv3, false /*no NULL correct*/);
            bool eq = bv1.equal(bv2);
            if (!eq)
            {
                cerr << "ERROR: Failed GT (sorted) integrity checks (1)!" << endl;
                assert(0); exit(1);
            }
            eq = bv1.equal(bv3);
            if (!eq)
            {
                cerr << "ERROR: Failed GT (sorted) integrity checks (2)!" << endl;
                assert(0); exit(1);
            }
        } // for
    }

    {
        bm::chrono_taker tt(cout, "std::vector<> GT (signed) scan (lower_bound) ", search_repeats);
        for (unsigned i = 0; i < search_repeats; ++i)
        {
            auto vs = search_vect[i];
            vector_search_GT_sorted(vect, bv_null, vs, bv_res3);
        } // for
    }


    vect.resize(0);
    vect.shrink_to_fit();

    {
        bm::chrono_taker tt(cout, "horizontal sparse vector scanner (signed) find_GT()", search_repeats);
        for (unsigned i = 0; i < search_repeats; ++i)
        {
            auto vs = search_vect[i];
            scanner.find_gt_horizontal(sv, vs, bv_res2, true /*NULL correct*/);
        } // for
    }

    {
        bm::chrono_taker tt(cout, "horizontal sparse vector scanner find_GT()(signed) (sorted)", search_repeats);
        for (unsigned i = 0; i < search_repeats; ++i)
        {
            auto vs = search_vect[i];
            scanner.find_gt_horizontal(sv2, vs, bv_res4, true /*NULL correct*/);
        } // for
    }

}


static
void SparseVectorScannerTest()
{
    std::vector<unsigned> vect;
    bvect bv_null;
    sparse_vector_u32 sv(bm::use_null);
    
    generate_scanner_test_set(vect, bv_null, sv, BSIZE/4);

    // generate a search vector for benchmarking
    std::vector<unsigned> search_vect;
    {
        bm::bvector<> bv_tmp;
        search_vect.reserve(REPEATS);
        for (unsigned i = 0; i < REPEATS;)
        {
            bm::id_t idx = bm::id_t(rand_dis(gen));
            if (!bv_tmp.test(idx)) // check if number is unique
            {
                search_vect.push_back(idx);
                bv_tmp[idx] = 1;
                ++i;
            }
        }
    }

    bm::sparse_vector_scanner<sparse_vector_u32> scanner;

    bvect bv_res1, bv_res2, bv_res3, bv_res4;

    unsigned search_repeats = REPEATS;
    {
        bm::chrono_taker tt(cout, "std::vector<> scan ", search_repeats);
        for (unsigned i = 0; i < search_repeats; ++i)
        {
            unsigned vs = search_vect[i];
            vector_search(vect, bv_null, vs, bv_res1);
        } // for
    }
    vect.resize(0);
    vect.shrink_to_fit();

    {
    bm::chrono_taker tt(cout, "horizontal sparse vector scanner find_eq()", search_repeats);
    for (unsigned i = 0; i < search_repeats; ++i)
    {
        {
        bvect bv1;
        scanner.find_eq_with_nulls_horizontal(sv, search_vect[i], bv1);
        bv_res2 |= bv1;
        }
    } // for
    scanner.correct_nulls(sv, bv_res2);
    }
    
    {
    bm::chrono_taker tt(cout, "sparse vector scanner find_eq() ", search_repeats);
    {
        scanner.find_eq(sv, search_vect.begin(), search_vect.end(), bv_res3);
    } // for
    }
    
    int res = bv_res3.compare(bv_res1);
    if (res != 0)
    {
        std::cerr << "1. Sparse scanner integrity check failed!" << std::endl;
        exit(1);
    }

    res = bv_res2.compare(bv_res1);
    if (res != 0)
    {
        std::cerr << "2. Sparse scanner integrity check failed!" << std::endl;
        exit(1);
    }

    sv.freeze();
    {
    bm::chrono_taker tt(cout, "sparse vector scanner find_eq() (RO)", search_repeats);
    {
        scanner.find_eq(sv, search_vect.begin(), search_vect.end(), bv_res4);
    } // for
    }

    res = bv_res4.compare(bv_res1);
    if (res != 0)
    {
        std::cerr << "3. Sparse scanner integrity check failed! (RO)" << std::endl;
        exit(1);
    }


}



inline
void GeneratePipelineTestData(std::vector<string>& str_coll,
                              str_svect_type&      str_sv,
                              unsigned max_coll = 8000000,
                              unsigned repeat_factor=10)
{
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
            bi = str;
        }
    } // for i
    bi.flush();
}

static
void SparseVectorPipelineScannerTest()
{
   const unsigned max_coll = 8000000;
   std::vector<string> str_coll;
   str_svect_type      str_sv;


   GeneratePipelineTestData(str_coll, str_sv, max_coll, 10);

    str_sv.remap();
    str_sv.optimize();


    //bm::print_svector_stat(str_sv);

    unsigned test_runs = 10000;
    std::vector<string> str_test_coll;
    for (bvect::size_type i = 0; i < test_runs; ++i)
    {
        bvect::size_type idx = (unsigned) rand() % test_runs;
        if (idx >= test_runs)
            idx = test_runs/2;
        str_test_coll.push_back(str_coll[idx]);
    }
    assert(str_test_coll.size() == test_runs);

    unsigned search_repeats = test_runs;

    std::vector<unique_ptr<bvect> > res_vec1;
    bm::sparse_vector_scanner<str_svect_type> scanner;

    {
        bm::chrono_taker tt(cout, "scanner::find_eq_str()", search_repeats);

        for (bvect::size_type i = 0; i < test_runs; ++i)
        {
            const string& str = str_test_coll[i];
            str_svect_type::bvector_type* bv_res(new bvect);
            scanner.find_eq_str(str_sv, str.c_str(), *bv_res);
            res_vec1.emplace_back(unique_ptr<bvect>(bv_res));
        } // for
    }

    bm::sparse_vector_scanner<str_svect_type>::pipeline<> pipe(str_sv);
    {
        bm::chrono_taker tt(cout, "scanner::pipeline find_eq_str()", search_repeats);

        for (bvect::size_type i = 0; i < test_runs; ++i)
        {
            const string& str = str_test_coll[i];
            pipe.add(str.c_str());
        }
        pipe.complete(); // finish the pipeline construction with this call
        scanner.find_eq_str(pipe); // run the search pipeline
    }

    bm::sparse_vector_scanner<str_svect_type>::pipeline<bm::agg_opt_only_counts> pipe2(str_sv);
    {
        bm::chrono_taker tt(cout, "scanner::pipeline find_eq_str()-count()", search_repeats);

        for (bvect::size_type i = 0; i < test_runs; ++i)
        {
            const string& str = str_test_coll[i];
            pipe2.add(str.c_str());
        }
        pipe2.complete(); // finish the pipeline construction with this call
        scanner.find_eq_str(pipe2); // run the search pipeline
    }

    unsigned c1 = 0;
    bvect bv_or_acc;
    {
        auto& res_vect = pipe.get_bv_res_vector();
        auto& cnt_vect = pipe2.get_bv_count_vector();

        for (size_t i = 0; i < res_vect.size(); ++i)
        {
            const bvect* bv1 = res_vec1[i].get();
            const auto* bv = res_vect[i];
            assert(bv);
            bool match = bv1->equal(*bv);
            assert(match); (void)match;
            auto c = cnt_vect[i];
            auto cnt = bv->count();
            assert(cnt == c);
            c1 += cnt; c1+= c;
            bv_or_acc |= *bv;
        }
    }
    if (!c1)
        cout << "\r";

    bvect bv_or;
    typedef bm::agg_run_options<false, false> scanner_custom_opt;
    bm::sparse_vector_scanner<str_svect_type>::pipeline<scanner_custom_opt> pipe3(str_sv);
    {
        bm::chrono_taker tt(cout, "scanner::pipeline find_eq_str()-OR", search_repeats);
        pipe3.set_or_target(&bv_or); // Assign OR aggregation target

        for (bvect::size_type i = 0; i < test_runs; ++i)
        {
            const string& str = str_test_coll[i];
            pipe3.add(str.c_str());
        }
        pipe3.complete(); // finish the pipeline construction with this call
        scanner.find_eq_str(pipe3); // run the search pipeline
    }
    bool match = bv_or.equal(bv_or_acc);
    if (!match)
    {
        cerr << "scanner::pipeline<>-OR vector mismatch!" << endl;
        assert(0);
        exit(1);
    }

    bv_or.clear(true);

    typedef bm::agg_run_options<false, true> scanner_custom_opt4;
    bm::sparse_vector_scanner<str_svect_type>::pipeline<scanner_custom_opt4> pipe4(str_sv);
    {
        bm::chrono_taker tt(cout, "scanner::pipeline find_eq_str()-count-OR", search_repeats);
        pipe4.set_or_target(&bv_or); // Assign OR aggregation target

        for (bvect::size_type i = 0; i < test_runs; ++i)
        {
            const string& str = str_test_coll[i];
            pipe4.add(str.c_str());
        }
        pipe4.complete(); // finish the pipeline construction with this call
        scanner.find_eq_str(pipe4); // run the search pipeline
    }

    match = bv_or.equal(bv_or_acc);
    if (!match)
    {
        cerr << "scanner::pipeline<>-count-OR vector mismatch!" << endl;
        assert(0);
        exit(1);
    }

    bv_or.clear(true);
    str_sv.freeze();

    bm::sparse_vector_scanner<str_svect_type>::pipeline<scanner_custom_opt4> pipe5(str_sv);
    {
        bm::chrono_taker tt(cout, "scanner::pipeline find_eq_str()-count-OR (RO)", search_repeats);
        pipe5.set_or_target(&bv_or); // Assign OR aggregation target

        for (bvect::size_type i = 0; i < test_runs; ++i)
        {
            const string& str = str_test_coll[i];
            pipe5.add(str.c_str());
        }
        pipe5.complete(); // finish the pipeline construction with this call
        scanner.find_eq_str(pipe5); // run the search pipeline
    }

    match = bv_or.equal(bv_or_acc);
    if (!match)
    {
        cerr << "scanner::pipeline<>-count-OR vector mismatch (RO)!" << endl;
        assert(0);
        exit(1);
    }



    {
        auto& res_vect = pipe.get_bv_res_vector();
        auto& cnt_vect = pipe4.get_bv_count_vector();

        for (size_t i = 0; i < res_vect.size(); ++i)
        {
            const bvect* bv1 = res_vec1[i].get();
            auto c = cnt_vect[i];
            auto cnt = bv1->count();
            assert(cnt == c);
            if (cnt != c)
            {
                cerr << "Failed count check!" << endl;
                exit(1);
            }
        }
    }

}


static
void SparseVectorSerializationTest()
{
    const unsigned char* buf;
    bool eq;
    size_t sz1, sz2;

    sparse_vector_u32 sv1(bm::use_null);
    sparse_vector_u32 sv2(bm::use_null);
    sparse_vector_u32 sv3(bm::use_null);

    unsigned sv_size = BSIZE;
    if constexpr (sizeof(void*) == 4)
    {
        sv_size = sv_size / 2;
    }

    generate_serialization_test_set(sv1, sv_size);

    bm::sparse_vector_serial_layout<sparse_vector_u32> sv_lay;

    bm::sparse_vector_serializer<sparse_vector_u32> sv_serializer;
    bm::sparse_vector_deserializer<sparse_vector_u32> sv_deserial;


    {
        {
            bm::chrono_taker tt(cout, "bm::sparse_vector<> serialization XOR disabled ", 1);

            sv_serializer.set_xor_ref(false); // disable XOR compression
            sv_serializer.serialize(sv1, sv_lay);
        }

        buf = sv_lay.buf();
        sz1 = sv_lay.size();

        sv_deserial.deserialize(sv2, buf);

        eq = sv1.equal(sv2);
        if (!eq)
        {
            cerr << "Error: SparseVectorSerializationTest() integrity failure! (1)" << endl;
            sparse_vector_u32::size_type pos;

            bool f = bm::sparse_vector_find_first_mismatch(sv1, sv2, pos);
            assert(f);
            cerr << "Mismatch at: " << pos << " " << f << endl;

            sv_deserial.deserialize(sv2, buf);

            exit(1);
        }
        sv2.resize(0);
    }

    {
        bm::chrono_taker tt(cout, "bm::sparse_vector<> serialization XOR enabled ", 1);
        sv_serializer.set_xor_ref(true); // enable XOR compression
        sv_serializer.serialize(sv1, sv_lay);
    }

    buf = sv_lay.buf();
    sz2 = sv_lay.size();
    sv_deserial.deserialize(sv3, buf);
    eq = sv1.equal(sv3);
    if (!eq)
    {
        cerr << "Error: SparseVectorSerializationTest() integrity failure! (2)" << endl;
        sparse_vector_u32::size_type pos;
        bool f = bm::sparse_vector_find_first_mismatch(sv1, sv3, pos);
        assert(f); (void)f;
        cerr << "Mismatch at: " << pos << endl;

        sv_deserial.deserialize(sv3, buf);

        exit(1);
    }

    if (sz2 >= sz1)
    {
        cerr << "XOR negative compression!" << endl;
        assert(0);
    }
    else
    {
        //cout << "sz1 = " << sz1 << " gain=" << (sz1 - sz2) << endl;
    }

}

static
void SparseVectorRangeDeserializationTest()
{
    std::vector<unsigned> vect;
    bvect bv_null;
    sparse_vector_u32 sv1(bm::use_null);
    sparse_vector_u32 sv2(bm::use_null);

    unsigned sv_size = BSIZE;
    if (bm::conditional<sizeof(void*) == 4>::test())
    {
        sv_size = sv_size / 2;
    }
    generate_scanner_test_set(vect, bv_null, sv1, sv_size);

    bm::sparse_vector_deserializer<sparse_vector_u32> sv_deserial;
    bm::sparse_vector_serializer<sparse_vector_u32> sv_serializer;
    sv_serializer.set_bookmarks(false);

    bm::sparse_vector_serial_layout<sparse_vector_u32> sv_lay;
    const unsigned char* buf;

    sv_serializer.serialize(sv1, sv_lay);

    buf = sv_lay.buf();
    
    {
        bm::chrono_taker tt(cout, "bm::sparse_vector<> Range Deserialization() - NO bookmarks ", 1);
        for (unsigned i = 0; i < 15; ++i)
        {
            sv_deserial.deserialize(sv2, buf, 0, 65536 * 2);
            sv_deserial.deserialize(sv2, buf, sv_size / 4, sv_size / 2);
            sv_deserial.deserialize(sv2, buf, sv_size / 2, (65536 * 2) + (sv_size / 2));
        }
    }

    assert(sv1.size() == sv2.size());
    // validation
    {
        sparse_vector_u32::size_type to = (65536 * 2) + (sv_size / 2);
        for (sparse_vector_u32::size_type from = sv_size / 2; from <= to; ++from)
        {
            auto v1 = sv1[from];
            auto v2 = sv2[from];
            if (v1 != v2)
            {
                cerr << "Range extraction failed!" << endl;
                exit(1);
            }
        }
    }


    // book-mark enabled serialization
    //
    sv_serializer.set_bookmarks(true, 64);
    sv_serializer.serialize(sv1, sv_lay);

    buf = sv_lay.buf();

    {
        bm::chrono_taker tt(cout, "bm::sparse_vector<> Range Deserialization() - WITH bookmarks ", 1);
        for (unsigned i = 0; i < 15; ++i)
        {
            sv_deserial.deserialize_range(sv2, buf, 0, 65536 * 2);
            sv_deserial.deserialize_range(sv2, buf, sv_size / 4, sv_size / 2);
            sv_deserial.deserialize_range(sv2, buf, sv_size / 2, (65536 * 2) + (sv_size / 2));
        }
    }

    // validation
    //
    {
        sparse_vector_u32::size_type to = (65536 * 2) + (sv_size / 2);
        for (sparse_vector_u32::size_type from = sv_size / 2; from <= to; ++from)
        {
            auto v1 = sv1[from];
            auto v2 = sv2[from];
            if (v1 != v2)
            {
                cerr << "Bookmark Range extraction failed!" << endl;
                exit(1);
            }
        }
    }

}






static
void StrSparseVectorTest()
{
    unsigned max_coll = 30000000;
    if constexpr (sizeof(void*) == 4)
    {
        max_coll = max_coll / 2;
    }

   std::vector<string> str_coll;
   str_svect_type str_sv;

   GenerateTestStrCollection(str_coll, max_coll);

    {
       bm::chrono_taker tt(cout, "bm::str_sparse_vector<>::push_back() ", 1);
       for (auto str : str_coll)
       {
           str_sv.push_back(str);
       }
    }
    str_sv.optimize();

    {
       str_svect_type str_sv0;
       {
           bm::chrono_taker tt(cout, "bm::str_sparse_vector<>::back_insert_iterator ", 1);
           str_svect_type::back_insert_iterator bi = str_sv0.get_back_inserter();
           for (auto str : str_coll)
           {
               bi = str;
           }
           bi.flush();
       }
       bool eq = str_sv.equal(str_sv0);
       if (!eq)
       {
            cerr << "Bit-transposed str_sv comparison failed(1)!" << endl;
            exit(1);
       }

       str_svect_type str_sv1;
       {
           bm::chrono_taker tt(cout, "bm::str_sparse_vector<>::remap ", 1);
           str_sv1.remap_from(str_sv0);
       }
       str_sv1.optimize();

       {
            str_svect_type::const_iterator it = str_sv.begin();
            str_svect_type::const_iterator it_end = str_sv.end();
            str_svect_type::const_iterator it1 = str_sv1.begin();
            for (; it < it_end; ++it, ++it1)
            {
                assert(it1.valid());
                const char* s = *it;
                const char* s1 = *it1;
                int cmp = ::strcmp(s, s1);
                if (cmp != 0)
                {
                    cerr << "Bit-transposed(remap) str_sv comparison failed(2)!" << endl;
                    exit(1);
                }
            }
       }
    }

    {
    string str;

        bm::chrono_taker tt(cout, "bm::str_sparse_vector<> - random access ", 1);
        for (unsigned i = 0; i < str_sv.size(); ++i)
        {
            str_sv.get(i, str);
            const std::string& sc = str_coll[i];
            if (str != sc)
            {
                cerr << "String random access check failure!" << endl;
                exit(1);
            }
        }
    }

    {
        bm::chrono_taker tt(cout, "bm::str_sparse_vector<>::const_iterator ", 1);
        str_svect_type::const_iterator it = str_sv.begin();
        str_svect_type::const_iterator it_end = str_sv.end();

        for (unsigned i=0; it < it_end; ++it, ++i)
        {
            const char* s = *it;
            const std::string& sc = str_coll[i];
            int cmp = ::strcmp(s, sc.c_str());
            if (cmp != 0)
            {
                cerr << "String random access check failure!" << endl;
                exit(1);
            }
        }
    }

    {
        bm::chrono_taker tt(cout, "bm::str_sparse_vector<>::const_iterator (remap)", 1);
        str_svect_type::const_iterator it = str_sv.begin();
        str_svect_type::const_iterator it_end = str_sv.end();

        for (unsigned i=0; it < it_end; ++it, ++i)
        {
            const char* s = *it;
            const std::string& sc = str_coll[i];
            int cmp = ::strcmp(s, sc.c_str());
            if (cmp != 0)
            {
                cerr << "String random access check failure!" << endl;
                exit(1);
            }
        }
    }



    // ---------
    // test on sorted collection

    std::sort(str_coll.begin(), str_coll.end());

    std::vector<string> str_coll_test1;
    std::vector<string> str_coll_test2;

    const unsigned pick_factor = 5;
    for (unsigned i = 0; i < unsigned(str_coll.size()); i+=pick_factor)
    {
        const string& s = str_coll[i];
        str_coll_test1.push_back(s);
        string s_imposs = s;

        if (rand() & 1)
            s_imposs.append("ATGCCCqz");
        else
            s_imposs = "1200ATGCCCqz";

        str_coll_test2.push_back(s_imposs);
    }
    std::random_device rd;
    std::mt19937 g(rd());

    std::shuffle(str_coll_test1.begin(), str_coll_test1.end(), g);


   str_svect_type str_sv_srt;
   {
       str_svect_type::back_insert_iterator bi = str_sv_srt.get_back_inserter();
       for (auto str : str_coll)
           bi = str;
       bi.flush();
   }
   str_sv_srt.remap();
   str_sv_srt.optimize();
   str_sv_srt.freeze();
/*
    {
       unsigned k =0;
       for (auto str : str_coll)
       {
           string s;
           str_sv_srt.get(k, s);
           assert(str == s);
           ++k;
       }
    }
*/


    bm::id64_t  f_sum4_0{0}, f_sum1{0}, f_sum2_4{0}, f_sum2_8{0}, f_sum2_16{0},
                f_sum2_32{0}, f_sum2_64{0};
    unsigned pos2;

    {
        bm::chrono_taker tt(cout, "std::lower_bound()", 1);

        for (unsigned i = 0; i < unsigned(str_coll_test1.size()); ++i)
        {
            const string& s = str_coll_test1[i];
            auto it = std::lower_bound(str_coll.begin(), str_coll.end(), s);
            assert(it != str_coll.end());
            auto pos1 = it - str_coll.begin();
            f_sum1 += (unsigned)pos1;
        }
    }

    {
    bm::sparse_vector_scanner<str_svect_type, 4> scanner_4;
        bm::chrono_taker tt(cout, "bm::sparse_vector_scanner<>::bfind_eq_str() [no-bind] [4] (sort/remap/ro)", 1);
        for (unsigned i = 0; i < unsigned(str_coll_test1.size()); ++i)
        {
            const string& s = str_coll_test1[i];
            bool found2 = scanner_4.bfind_eq_str(str_sv_srt, s.c_str(), pos2);
            if (!found2)
            {
                cerr << "String bfind_eq_str() failure!" << endl;
                cerr << s << endl;
                bm::file_save_svector(str_sv_srt, "perf.str_sv");
                cout << "debug dump created: " << "perf.str_sv" << endl;
                assert(0); exit(1);
            }
            f_sum4_0 += pos2;
        }
        assert(f_sum4_0 == f_sum1);
    }

    {
    bm::sparse_vector_scanner<str_svect_type, 4> scanner_4;
    scanner_4.bind(str_sv_srt, true); // bind sorted vector

        bm::chrono_taker tt(cout, "bm::sparse_vector_scanner<>::bfind_eq_str() [4] (sort/remap/ro)", 1);
        for (unsigned i = 0; i < unsigned(str_coll_test1.size()); ++i)
        {
            const string& s = str_coll_test1[i];
            bool found2 = scanner_4.bfind_eq_str(s.c_str(), pos2);
            if (!found2)
            {
                cerr << "String bfind_eq_str() failure!" << endl;
                assert(0); exit(1);
            }
            f_sum2_4 += pos2;
        }
        assert(f_sum2_4 == f_sum1);
    }

    {
    bm::sparse_vector_scanner<str_svect_type, 8> scanner_8;
    scanner_8.bind(str_sv_srt, true); // bind sorted vector

        bm::chrono_taker tt(cout, "bm::sparse_vector_scanner<>::bfind_eq_str() [8] (sort/remap/ro)", 1);
        for (unsigned i = 0; i < unsigned(str_coll_test1.size()); ++i)
        {
            const string& s = str_coll_test1[i];
            bool found2 = scanner_8.bfind_eq_str(s.c_str(), pos2);
            if (!found2)
            {
                cerr << "String bfind_eq_str() failure!" << endl;
                assert(0); exit(1);
            }
            f_sum2_8 += pos2;
        }
        assert(f_sum2_8 == f_sum1);
    }

    {
    bm::sparse_vector_scanner<str_svect_type, 16> scanner_16;
    scanner_16.bind(str_sv_srt, true); // bind sorted vector

        bm::chrono_taker tt(cout, "bm::sparse_vector_scanner<>::bfind_eq_str() [16] (sort/remap/ro)", 1);
        for (unsigned i = 0; i < unsigned(str_coll_test1.size()); ++i)
        {
            const string& s = str_coll_test1[i];
            bool found2 = scanner_16.bfind_eq_str(s.c_str(), pos2);
            if (!found2)
            {
                cerr << "String bfind_eq_str() failure!" << endl;
                assert(0); exit(1);
            }
            f_sum2_16 += pos2;
        }
        assert(f_sum2_16 == f_sum1);
    }

    {
    bm::sparse_vector_scanner<str_svect_type, 32> scanner_32;
    scanner_32.bind(str_sv_srt, true); // bind sorted vector

        bm::chrono_taker tt(cout, "bm::sparse_vector_scanner<>::bfind_eq_str() [32] (sort/remap/ro)", 1);
        for (unsigned i = 0; i < unsigned(str_coll_test1.size()); ++i)
        {
            const string& s = str_coll_test1[i];
            bool found2 = scanner_32.bfind_eq_str(s.c_str(), pos2);
            if (!found2)
            {
                cerr << "String bfind_eq_str() failure!" << endl;
                assert(0); exit(1);
            }
            f_sum2_32 += pos2;
        }
        assert(f_sum2_32 == f_sum1);
    }

    {
    bm::sparse_vector_scanner<str_svect_type, 64> scanner_64;
    scanner_64.bind(str_sv_srt, true); // bind sorted vector

        bm::chrono_taker tt(cout, "bm::sparse_vector_scanner<>::bfind_eq_str() [64] (sort/remap/ro)", 1);
        for (unsigned i = 0; i < unsigned(str_coll_test1.size()); ++i)
        {
            const string& s = str_coll_test1[i];
            bool found2 = scanner_64.bfind_eq_str(s.c_str(), pos2);
            if (!found2)
            {
                cerr << "String bfind_eq_str() failure!" << endl;
                assert(0); exit(1);
            }
            f_sum2_64 += pos2;
        }
        assert(f_sum2_64 == f_sum1);
    }

    if (f_sum1 != f_sum4_0 || f_sum1 != f_sum2_4 || f_sum1 != f_sum2_8 ||
        f_sum1 != f_sum2_16 || f_sum1 != f_sum2_32 || f_sum1 != f_sum2_64)
    {
        cerr << "Search positions sum mismatch" << endl;
        assert(0); exit(1);
    }

    {
    bm::sparse_vector_scanner<str_svect_type, 16> scanner_16;
    scanner_16.bind(str_sv_srt, true); // bind sorted vector

        bm::chrono_taker tt(cout, "bm::bm::sparse_vector_scanner<>::bfind_eq_str() [16] (sort/remap/ro) (not-found)", 1);

        for (unsigned i = 0; i < unsigned(str_coll_test2.size()); ++i)
        {
            const string& s = str_coll_test2[i];
            bool found2 = scanner_16.bfind_eq_str(s.c_str(), pos2);
            if (found2)
            {
                cerr << "String bfind_eq_str() failure!" << endl;
                assert(0); exit(1);
            }
        }
    }

}


/// Random numbers test
template<typename V>
unsigned generate_inter_test_linear(V* arr, unsigned inc, unsigned target_size)
{
    V maskFF = (V)~V(0u);
    if (inc < 2 || inc > 65535)
        inc = 1;

    unsigned start = 1;
    unsigned sz = 0;
    while (sz < target_size)
    {
        arr[sz++] = V(start);
        if (inc + start >= maskFF)
        {
            arr[sz++] = maskFF;
            break;
        }
        start += inc;
        if (start < arr[sz - 1])
            break;

    } // while
    return sz;
}

/// Random numbers test
template<typename V>
unsigned generate_inter_test(V* arr, unsigned inc_factor, unsigned target_size)
{
    V maskFF = (V)~V(0u);

    if (inc_factor < 2)
        inc_factor = 65535;

    unsigned start = (unsigned)rand() % 256;
    if (!start)
        start = 1;
    unsigned sz = 0;
    while (sz < target_size)
    {
        arr[sz++] = V(start);

        unsigned inc = unsigned(rand()) % inc_factor;
        if (!inc)
            inc = 1;
        start += inc;
        if (start >= maskFF)
            break;
    } // while

    for (unsigned i = 1; i < sz; ++i)
    {
        if (arr[i - 1] >= arr[i])
        {
            return i;
        }
    }
    return sz;
}


static
void InterpolativeCodingTest()
{
    unsigned char buf[1024 * 200] = { 0, };

    const unsigned code_repeats = 500000;
    const unsigned test_size = 12000;

    vector<unsigned> sa; sa.resize(test_size);
    vector<unsigned> da; da.resize(test_size);

    bm::word_t* src_arr = &sa[0];
    bm::word_t* dst_arr = &da[0];
    unsigned sz;
    unsigned inc = (unsigned)rand() % (65536 * 256);
    sz = generate_inter_test_linear(src_arr, inc, test_size);
    assert(sz);
    assert(src_arr[0]);
    {
        bm::encoder enc(buf, sizeof(buf));
        bm::bit_out<bm::encoder> bout(enc);

        bout.bic_encode_u32_cm(src_arr, sz - 1, 0, src_arr[sz - 1]);
        bout.flush();
        auto ssz = enc.size();
        assert(ssz < sizeof(buf));
        (void)ssz;
    }

    {
        bm::chrono_taker tt(cout, "bic_decode_u32_cm() ", 1);

        for (unsigned k = 0; k < code_repeats; ++k)
        {
            bm::decoder dec(buf);
            bm::bit_in<bm::decoder> bin(dec);

            bin.bic_decode_u32_cm(&dst_arr[0], sz - 1, 0, src_arr[sz - 1]);
            dst_arr[sz - 1] = src_arr[sz - 1];
            for (unsigned i = 0; i < sz; ++i)
            {
                assert(src_arr[i] == dst_arr[i]);
                if (i)
                {
                    assert(src_arr[i - 1] < src_arr[i]);
                }
            }
        } // for k
    }
}



int main(void)
{
    cout << bm::_copyright<true>::_p << endl;
    cout << "SIMD code = " << bm::simd_version() << endl;
    #if defined (BM64OPT)
        cout << " 64-bit optimizations: ON" << endl;
    #else
        cout << " 64-bit optimizations: OFF" << endl;
    #endif

//    ptest();

    bm::chrono_taker tt(cout, "TOTAL", 1);
    try
    {
        cout << endl;

        MemCpyTest();
        cout << endl;

        BitCountTest();
        cout << endl;

        BitCountSparseTest();
        cout << endl;

        BitTestSparseTest();
        cout << endl;

        BitForEachTest();
        cout << endl;

        WordSelectTest();
        cout << endl;


        BitSetConditionalTest();
        cout << endl;


        BitCompareTest();
        cout << endl;

        OptimizeTest();
        cout << endl;

        FindTest();
        cout << endl;

        BitBlockRotateTest();
        BitBlockShiftTest();
        cout << endl;

        InterpolativeCodingTest();
        cout << endl;

        EnumeratorTest();
        EnumeratorTestGAP();
        EnumeratorGoToTest();
        cout << endl;

        BvectorShiftTest();
        cout << endl;

        RangeCopyTest();
        cout << endl;

        IntervalsTest();
        cout << endl;

        AggregatorTest();
        cout << endl;

        OrTest();
        AndTest();
        XorTest();
        SubTest();

        cout << endl;
        OrAndTest();
        cout << endl;

        XorCountTest();
        AndCountTest();
        TI_MetricTest();
        cout << endl;

        InvertTest();
        cout << endl;

        Set2SetTransformTest();
        cout << endl;

        SerializationTest();
        cout << endl;

        SparseVectorAccessTest();
        cout << endl;

        SparseVectorSignedAccessTest();
        cout << endl;

        SparseVectorScannerTest();
        cout << endl;

        SparseVectorScannerCmpTest();
        cout << endl;

        SparseVectorScannerSignedCmpTest();
        cout << endl;

        SparseVectorPipelineScannerTest();
        cout << endl;

        SparseVectorSerializationTest();
        SparseVectorRangeDeserializationTest();
        cout << endl;

        RSC_SparseVectorFillTest();

        RSC_SparseVectorAccesTest();

        RSC_SparseVectorRandomAccesTest();
        cout << endl;

        RankCompressionTest();
        cout << endl;

        StrSparseVectorTest();
        cout << endl;

    }
    catch (std::exception& ex)
    {
        cerr << ex.what() << endl;
        exit(1);
    }
    return 0;
}


#ifdef _MSC_VER
#pragma warning( pop )
#endif



