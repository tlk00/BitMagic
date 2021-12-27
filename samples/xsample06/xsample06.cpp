/*
Copyright(c) 2002-2019 Anatoliy Kuznetsov(anatoliy_kuznetsov at yahoo.com)

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

/** \example xsample06.cpp
  Use of sparse vector for compressed DNA strings (2 bits per bp)
  Construction of comparison functions on bit-transposed compressed
  containers, benchmarking.
 
  \sa bm::sparse_vector
  \sa bm::sparse_vector_find_first_mismatch

*/

/*! \file xsample06.cpp
    \brief Example: Use of sparse vector for compressed DNA strings
*/

#include <assert.h>
#include <stdlib.h>

#include <iostream>
#include <vector>
#include <algorithm>
#include <utility>

#include "bm.h"
#include "bmsparsevec.h"
#include "bmsparsevec_algo.h"
#include "bmsparsevec_serial.h"

// BitMagic utilities for debug and timings
#include "bmdbg.h"
#include "bmtimer.h"
#include "bmundef.h" /* clear the pre-proc defines from BM */

using namespace std;


/// Translate DNA letter to integer code
/// Please note this function uses extended alphabet ATGC plus 'N' and '$'
/// Ns are used to indicate ambiguity and $ is part of Burrows-Wheeler transform
/// 
inline unsigned DNA2int(char DNA_bp)
{
    switch (DNA_bp)
    {
    case 'A':
        return 0; // 000
    case 'T':
        return 1; // 001
    case 'G':
        return 2; // 010
    case 'C':
        return 3; // 011
    case 'N':
        return 4; // 100
    case '$':
        return 5; // 101
    default:
        assert(0);
        return 0;
    }
}
/// Translate integer code to DNA letter
///
inline char int2DNA(unsigned code)
{
    static char lut[] = { 'A', 'T', 'G', 'C', 'N', '$' };
    if (code < 5)
        return lut[code];
    assert(0);
    return 'N';
}

/// Print sparse vector 
template<typename SV> void PrintSV(const SV& sv)
{
    cout << "size() = " << sv.size() << " : ";
    auto it = sv.begin(); auto it_end = sv.end();
    for (; it != it_end; ++it)
    {
        auto v = *it;
        char bp = int2DNA(v);
        cout << bp;
    }
    cout << endl;
}

/// Test sparse vector against reference vector 
/// (for QA purposes)
///
template<typename SV, typename VECT> 
void TestSV(const SV& sv, const VECT& vect)
{
    auto it = sv.begin(); auto it_end = sv.end();
    auto itv = vect.begin(); auto itv_end = vect.end();
    for (; it != it_end && itv != itv_end; ++it, ++itv)
    {
        auto v = *it;
        char bp = int2DNA(v);
        char bpv = *itv;
        assert(bp == bpv);
        if (bp != bpv)
        {
            cerr << "Error: reference vector mismatch!" << endl;
            exit(1);
        }
    }
    if (it != it_end || itv != itv_end)
    {
        cerr << "Error: reference size mismatch!" << endl;
        exit(1);
    }
}


typedef bm::sparse_vector<unsigned, bm::bvector<> > svector_u32;
typedef std::vector<char>                           vector_char_type;
typedef std::vector<std::pair<unsigned, unsigned> > vector_pairs_type;

/*
    Lexicographical compare of two sparse vectors using bit-plane
    mismatch search
*/
static
int compare_sv(const svector_u32& sv1, const svector_u32& sv2)
{
    svector_u32::size_type pos;
    bool found = bm::sparse_vector_find_first_mismatch(sv1, sv2, pos);
    if (found)
    {
        auto v1 = sv1[pos];
        auto v2 = sv2[pos];
        if (v1 < v2)
            return -1;
        return 1;
    }
    return 0; // EQ
}

/*
    Lexicographical compare of two sparse vectors using
    iterators. This variant is quite slow, because sparse vector iterators
    are heavy classes doing reverse transformation on the fly
*/
static
int compare_sv_it(const svector_u32& sv1, const svector_u32& sv2)
{
    svector_u32::const_iterator it1 = sv1.begin();
    svector_u32::const_iterator it2 = sv2.begin();
    svector_u32::size_type sz = sv1.size();
    for (svector_u32::size_type i = 0; i < sz; ++i, ++it1, ++it2)
    {
        auto v1 = *it1; auto v2 = *it2;
        if (v1 == v2)
            continue;
        if (v1 < v2)
            return -1;
        return 1;
    }
    return 0;
}

/*
    Utility function to generate random vector of DNA chars (4 letters)
*/
static
void generate_DNA_vector(svector_u32& sv, vector_char_type& vect, unsigned sz)
{
    svector_u32::back_insert_iterator bi = sv.get_back_inserter();
    for (unsigned i = 0; i < sz; ++i)
    {
        unsigned code = rand() % 4; // generate code between 0 and 3
        char bp = int2DNA(code);
        assert(bp == 'A' || bp == 'T' || bp == 'G' || bp == 'C');
        vect.push_back(bp);
        bi = code;
    }
    bi.flush();
}

/*
    Utility function to add symulated centromer (block of Ns in the middle)
*/
static 
void add_centromer_Ns(svector_u32& sv, vector_char_type& vect, unsigned csz)
{
    assert(csz);

    svector_u32::size_type sz = sv.size();
    assert(sz == svector_u32::size_type(vect.size()));
    assert(csz < sz);

    svector_u32::size_type half = sz / 2; // center position
    svector_u32::size_type c_half = csz / 2; 
    svector_u32::size_type c_start = half - c_half; // simulated centromere start
    svector_u32::size_type c_end = half + c_half; // simulated centromere end

    // fill-in simulated centromere 'NNNNNN...NNNN' 
    for (svector_u32::size_type i = c_start; i < c_end; ++i)
    {
        vect[i] = 'N';
        sv[i] = DNA2int('N'); // 4
    }
}

/*
 Generate large vectors to simulate human chr1
*/
static
void generate_big_case(svector_u32& sv, vector_char_type& vect)
{
    const unsigned case_size =  200000000; // 200M bps
    const unsigned centr_size =   7000000; // 7M centromere Ns

    generate_DNA_vector(sv, vect, case_size);
    add_centromer_Ns(sv, vect, centr_size);

    sv.optimize(); // this will compress the sparse (N) plane

    TestSV(sv, vect); // paranoiya check
}

/*
    Generate benchmark set of mismatches
*/
static
void generate_mismatches(vector_pairs_type& vect_m,
                         const vector_char_type& vect,
                         unsigned                m_count)
{
    vect_m.resize(0);
    if (!vect.size())
        return;

    vector_char_type::size_type sz = vect.size();
    vector_char_type::size_type delta = sz / m_count;

    for (vector_char_type::size_type i = 0; i < sz; i += delta)
    {
        vector_char_type::value_type v1 = vect[i];
        unsigned code = rand() % 4;
        vector_char_type::value_type v2 = int2DNA(code);
        if (v2 == v1)
            continue;
        vect_m.push_back(vector_pairs_type::value_type(unsigned(i), code));
    } // for i

    // add some extra with a distrubution skewed to the beginning
    //
    for (vector_char_type::size_type i = 1; i < sz / 4; i += (rand()%(1024 * 10)))
    {
        vector_char_type::value_type v1 = vect[i];
        unsigned code = rand() % 4;
        vector_char_type::value_type v2 = int2DNA(code);
        if (v2 == v1)
            continue;
        vect_m.push_back(vector_pairs_type::value_type(unsigned(i), code));
    } // for i

}

// -------------------------------------------------------------------

const unsigned rep = 20000;
bm::chrono_taker<>::duration_map_type  timing_map;


/*
    Mismatch search benchmark for bm::sparse_vector_find_first_mismatch
*/
static
void test_mismatch_search(const svector_u32& sv, const vector_pairs_type& vect_m)
{
    svector_u32 sv1(sv); // copy sv
    svector_u32::size_type sz = (svector_u32::size_type)vect_m.size();

    unsigned cnt = 0;

    bm::chrono_taker tt1(cout, "1. SV mismatch", sz, &timing_map);

    for (unsigned k = 0; k < sz; ++k)
    {
        unsigned i = vect_m[k].first;
        svector_u32::value_type v1 = sv1[i];
        svector_u32::value_type v2 = vect_m[k].second;
        assert(v1 != v2);

        sv1.set(i, v2); // change the copy vector

        svector_u32::size_type pos;
        bool found = bm::sparse_vector_find_first_mismatch(sv, sv1, pos);
        cnt += found;
        assert(found);
        assert(pos == i);

        sv1.set(pos, v1); // restore old value
    } // for
    assert(cnt == vect_m.size());
}

/*
    Mismatch search benchmark for std::vector::const_iterator
*/
static
void test_vect_mismatch_search(const vector_char_type& vect,
                          const vector_pairs_type& vect_m)
{
    vector_char_type vect1(vect);
    unsigned sz = (unsigned)vect_m.size();

    unsigned cnt = 0;

    bm::chrono_taker tt1(cout, "2. STL iterator", sz, &timing_map);

    for (unsigned k = 0; k < sz; ++k)
    {
        unsigned i = vect_m[k].first;
        vector_char_type::value_type v1 = vect[i];
        vector_char_type::value_type v2 = int2DNA(vect_m[k].second);
        assert(v1 != v2);

        vect1[i] = v2; // change the copy vector

        // possible to use std::mismatch() here
        bool found = false;
        vector_char_type::size_type pos;
        auto it = vect.begin(); auto it_end = vect.end();
        auto it1 = vect1.begin();
        for (; it != it_end; ++it, ++it1)
        {
            if (*it != *it1)
            {
                found = true; 
                pos = vector_char_type::size_type(it - vect.begin());
                break;
            }
        } // for it
        cnt += found;
        assert(found);
        assert(pos == i);
        vect1[i] = v1; // restore old value
    } // for
    assert(cnt == vect_m.size());
}


/*
    Mismatch search benchmark for std::vector using strncmp()
*/
static
void test_strcmp(const vector_char_type& vect, const vector_pairs_type& vect_m)
{
    vector_char_type vect1(vect);
    unsigned sz = (unsigned)vect_m.size();
    unsigned cnt = 0;

    bm::chrono_taker tt1(cout, "3. strcmp() test ", sz, &timing_map);

    unsigned vs = unsigned(vect.size());

    for (unsigned k = 0; k < sz; ++k)
    {
        unsigned i = vect_m[k].first;
        vector_char_type::value_type v1 = vect[i];
        assert(vect_m[k].second < 4);
        vector_char_type::value_type v2 = int2DNA(vect_m[k].second);
        assert(v1 != v2);

        vect1[i] = v2; // change the copy vector

        const char* s1 = vect.data();
        const char* s2 = vect1.data();

        int res = ::strncmp(s1, s2, vs);
        cnt += bool(res);
        assert(res);
        vect1[i] = v1; // restore
    } // for
    assert(cnt == vect_m.size());
}

/*
    Sparse vector compare benchmark
*/
static
void test_sv_cmp(const svector_u32& sv, const vector_pairs_type& vect_m)
{
    svector_u32 sv1(sv);
    unsigned sz = (unsigned)vect_m.size();
    unsigned cnt = 0;

    bm::chrono_taker tt1(cout, "4. sv compare", sz, &timing_map);

    for (unsigned k = 0; k < sz; ++k)
    {
        unsigned i = vect_m[k].first;
        svector_u32::value_type v1 = sv1[i];
        svector_u32::value_type v2 = vect_m[k].second;
        assert(v1 != v2);

        sv1.set(i, v2); // change the copy vector

        int res = compare_sv(sv, sv1);
        assert(res != 0);
        cnt += bool(res);

        sv1.set(i, v1); // restore
    } // for
    assert(cnt == vect_m.size());
}

/*
    Sparse vector compare benchmark using iterators (very slow)
*/
inline
void test_sv_cmp_it(const svector_u32& sv, const vector_pairs_type& vect_m)
{
    svector_u32 sv1(sv);
    unsigned sz = (unsigned)vect_m.size();
    unsigned cnt = 0;

    bm::chrono_taker tt1(cout, "4. sv-cmp-it()", sz, &timing_map);

    for (unsigned k = 0; k < sz; ++k)
    {
        unsigned i = vect_m[k].first;
        svector_u32::value_type v1 = sv1[i];
        svector_u32::value_type v2 = vect_m[k].second;
        assert(v1 != v2);

        sv1.set(i, v2); // change the copy vector

        int res = compare_sv_it(sv, sv1);
        cnt += bool(res);

        assert(res != 0);

        sv1.set(i, v1); // restore
    } // for
    assert(cnt == vect_m.size());
}

// -------------------------------------------------------------------

int main(void)
{
    try
    {
        // generate short DNA vector as STL vector and BitMagic 
        // bit-transposed container, print the results (as a demo)
        // we don't have to re-load it every time,
        // bit transposed container is serializable
        {
            const char* s = "ATGTCNNNNNTATA";
            svector_u32 sv;
            {
                svector_u32::back_insert_iterator bi = sv.get_back_inserter();
                for (unsigned i = 0; s[i]; ++i)
                {
                    bi = DNA2int(s[i]);
                }
                bi.flush();
            }
            sv.optimize(); // this will compress the sparse (N) plane
            PrintSV(sv);
        }

        // run a battery of benchmarks for variants of mismatch searches
        // and compare functions
        {
            svector_u32 sv1;
            vector_char_type dna_vect;
            vector_pairs_type vect_m;

            cout << "Generate benchmark test data." << endl;
            generate_big_case(sv1, dna_vect);

            generate_mismatches(vect_m, dna_vect, rep);
            cout << "generated mismatches=" << vect_m.size() << endl;

            cout << "SV mismatch search test" << endl;
            test_mismatch_search(sv1, vect_m);

            cout << "STL interator mismatch test" << endl;
            test_vect_mismatch_search(dna_vect, vect_m);

            // test compare functions
            //

            cout << "strncmp() test" << endl;
            test_strcmp(dna_vect, vect_m);

            cout << "SV compare test" << endl;
            test_sv_cmp(sv1, vect_m);

            // really sow, not worth testing
            //
            #if 0
            cout << "SV compare iterators test" << endl;
            test_sv_cmp_it(sv1, vect_m);
            #endif

        }

        std::cout << std::endl << "Performance:" << std::endl;
        bm::chrono_taker<>::print_duration_map(cout, timing_map, bm::chrono_taker<>::ct_time);

    }
    catch(std::exception& ex)
    {
        std::cerr << ex.what() << std::endl;
        return 1;
    }
    
    

    return 0;
}

