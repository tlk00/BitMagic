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

/** \example svsample09.cpp
  Use of sparse vector for compressed DNA strings
 
  \sa bm::sparse_vector
  \sa bm::sparse_vector_find_first_mismatch

*/

/*! \file svsample09.cpp
    \brief Example: Use of sparse vector for compressed DNA strings
*/

#include <assert.h>
#include <stdlib.h>

#include <iostream>
#include <vector>
#include <algorithm>

#include "bm.h"
#include "bmsparsevec.h"
#include "bmsparsevec_algo.h"
#include "bmsparsevec_serial.h"

// BitMagic utilities for debug and timings
#include "bmdbg.h"
#include "bmtimer.h"


using namespace std;


/// Translate DNA letter to integer code
/// Please note this function uses extended alphabet ATGC plus N$
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

/// Generate large vectors to approximately simulate human chr1
///
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

bm::chrono_taker::duration_map_type  timing_map;

const unsigned repeats = 20000;

static
void test_mismatch_search(const svector_u32& sv)
{
    svector_u32 sv1(sv); // copy sv

    bm::chrono_taker tt1("1. SV mismatch search ", repeats, &timing_map);
    svector_u32::size_type sz = sv.size();
    svector_u32::size_type delta = sz / repeats;

    for (svector_u32::size_type i = 0; i < sz; i+=delta)
    {
        svector_u32::value_type v1 = sv1[i];
        svector_u32::value_type v2 = rand() % 4; // generate random code
        sv1.set(i, v2); // change the copy vector

        svector_u32::size_type pos;
        bool found = bm::sparse_vector_find_first_mismatch(sv, sv1, pos);
        if (found)
        {
            assert(pos == i);
            if (pos != i)
            {
                cerr << "Search error!" << endl; 
                exit(1);
            }
            sv1.set(i, v1); // restore old value
        }
        else
        {
            assert(v1 == v2);
            if (v1 != v2)
            {
                cerr << "Search mismatch error!" << endl;
                exit(1);
            }
        }
        //cout << "\r" << i << flush;
    } // for i
    //cout << endl;
}

static
void test_mismatch_search(const vector_char_type& vect)
{
    vector_char_type vect1(vect); 

    bm::chrono_taker tt1("2. STL iterator search ", repeats, &timing_map);
    vector_char_type::size_type sz = vect.size();
    vector_char_type::size_type delta = sz / repeats;

    for (vector_char_type::size_type i = 0; i < sz; i += delta)
    {
        vector_char_type::value_type v1 = vect[i];
        unsigned code = rand() % 4; 
        vector_char_type::value_type v2 = int2DNA(code); 

        vect1[i] = v2; // change the copy vector

        // possible to use std::mismatch here
        bool found = false;
        vector_char_type::size_type pos;
        auto it = vect.begin(); auto it_end = vect.end();
        auto it1 = vect1.begin();
        for (; it != it_end; ++it, ++it1)
        {
            if (*it != *it1)
            {
                found = true; 
                pos = it - vect.begin();
                break;
            }
        } // for it

        if (found)
        {
            //cout << "+" << flush;
            assert(pos == i);
            if (pos != i)
            {
                cerr << "Search error!" << endl;
                exit(1);
            }
            vect1[i] = v1; // restore old value

        }
        else
        {
            //cout << "-" << flush;
            assert(v1 == v2);
            if (v1 != v2)
            {
                cerr << "Search mismatch error!" << endl;
                exit(1);
            }
        }
        //cout << "\r" << i << flush;
    } // for i
    //cout << endl;

}


static
void test_strcmp(const vector_char_type& vect)
{
    vector_char_type vect1(vect);

    bm::chrono_taker tt1("3. strcmp() test ", repeats, &timing_map);
    vector_char_type::size_type sz = vect.size();
    vector_char_type::size_type delta = sz / repeats;

    for (vector_char_type::size_type i = 0; i < sz; i += delta)
    {
        vector_char_type::value_type v1 = vect[i];
        unsigned code = rand() % 4;
        vector_char_type::value_type v2 = int2DNA(code);

        vect1[i] = v2; // change the copy vector

        const char* s1 = vect.data();
        const char* s2 = vect1.data();

        int res = ::strncmp(s1, s2, sz);
        if (res == 0)
        {
            assert(v1 == v2);
            if (v1 != v2)
            {
                cerr << "Search mismatch error!" << endl;
                exit(1);
            }
        }
        else
        {
            if (v1 == v2)
            {
                cerr << "Search mismatch error!" << endl;
                exit(1);
            }
        }
        vect1[i] = v1; // restore
        //cout << "\r" << i << flush;
    } // for i
    //cout << endl;
}

static
void test_sv_cmp(const svector_u32& sv)
{
    svector_u32 sv1(sv);

    bm::chrono_taker tt1("4. sv-cmp() test ", repeats, &timing_map);
    svector_u32::size_type sz = sv.size();
    svector_u32::size_type delta = sz / repeats;

    for (svector_u32::size_type i = 0; i < sz; i += delta)
    {
        svector_u32::value_type v1 = sv[i];
        svector_u32::value_type v2 = rand() % 4;

        sv1[i] = v2; // change the copy vector

        int res = compare_sv(sv, sv1);
        if (res == 0)
        {
            assert(v1 == v2);
            if (v1 != v2)
            {
                cerr << "Search mismatch error!" << endl;
                exit(1);
            }
        }
        else
        {
            if (v1 == v2)
            {
                cerr << "Search mismatch error!" << endl;
                exit(1);
            }
            sv1[i] = v1;
            //sv1.set(i, v1); // restore
        }
        //cout << "\r" << i << flush;
    } // for i
    //cout << endl;
}

int main(void)
{
    try
    {
        // generate short DNA vector as STL vector and BitMagic 
        // bit-transposed container
        {
            svector_u32 sv1;
            vector_char_type dna_vect;

            generate_DNA_vector(sv1, dna_vect, 20);
            add_centromer_Ns(sv1, dna_vect, 6);

            sv1.optimize(); // this will compress the sparse (N) plane

            PrintSV(sv1);
            TestSV(sv1, dna_vect);
        }

        {
            svector_u32 sv1;
            vector_char_type dna_vect;

            cout << "Generate test data." << endl;
            generate_big_case(sv1, dna_vect);

            cout << "SV mismatch search test" << endl;
            test_mismatch_search(sv1);

            cout << "STL interator mismatch test" << endl;
            test_mismatch_search(dna_vect);

            cout << "::strncmp test" << endl;
            test_strcmp(dna_vect);

            cout << "SV compare test" << endl;
            test_sv_cmp(sv1);
        }

        std::cout << std::endl << "Performance:" << std::endl;
        bm::chrono_taker::print_duration_map(timing_map, bm::chrono_taker::ct_time);

    }
    catch(std::exception& ex)
    {
        std::cerr << ex.what() << std::endl;
        return 1;
    }
    
    

    return 0;
}

