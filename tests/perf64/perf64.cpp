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

#define BM64ADDR
#include "bm.h"
#include <bmsparsevec.h>
#include <bmsparsevec_algo.h>
#include <bmsparsevec_serial.h>
#include "bmtimer.h"
#include "bmdbg.h"

using namespace std;

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

bm::chrono_taker::duration_map_type  timing_map;

typedef bm::bvector<> bvect64;
typedef std::vector<bm::id64_t> ref_vect;

static
void bvector64_VerySparse_SetDestroyCycle()
{
    ref_vect vect;
    generate_vect_simpl0(vect);
    
    cout << "bvector64_VerySparse_SetDestroyCycle..." << endl;
    {
        bm::chrono_taker tt("bvector64_VerySparse_SetDestroyCycle", 1, &timing_map);
        const unsigned max_try = REPEATS;
        for (unsigned i = 0; i < max_try; ++i)
        {
            bvect64 bv0;
            load_BV_set_ref(bv0, vect, false);
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
        load_BV_set_ref(bv0, vect, false);

        cout << "bvector64_VerySparse_RefAccessCycle..." << endl;
        {
            bm::chrono_taker tt("bvector64_VerySparse_RefAccessCycle", 1, &timing_map);
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
            bm::chrono_taker tt("bvector64_VerySparse_EnumAccessCycle", 1, &timing_map);
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
int main(void)
{
    bvector64_VerySparse_SetDestroyCycle();
    bvector64_VerySparse_RefAccessCycle();

    std::cout << std::endl << "Performance:" << std::endl;
    bm::chrono_taker::print_duration_map(timing_map, bm::chrono_taker::ct_all);

    return 0;
}


#ifdef _MSC_VER
#pragma warning( pop )
#endif



