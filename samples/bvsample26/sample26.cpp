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

/** \example sample26.cpp
    Immutable bit-vectors.
    This example generates sparse bit-vectors and shows how to make them use less memory
    via optimization (after optimization bit-vector remains mutable) and freezing it into immutable,
    where memory all edit reservations are getting dropped.
    Freezing also allocates memory blocks together, reducing the heap fragmentation.

    \sa bm::bvector::optimize
    \sa bm::bvector::freeze
    \sa bm::bvector::is_ro
*/

/*! \file sample26.cpp
    \brief Example: bvector<> with immutability (read-only)
*/
#include <stdlib.h>
#include <iostream>
#include <cassert>

#include "bm.h"
#include "bmundef.h" /* clear the pre-proc defines from BM */

using namespace std;

const unsigned MAX_VALUE = 1000000;

/// Fill bit-vectors with values using dense and sparse distrubutions
///
static
void fill_bvector(bm::bvector<>* bv1, bm::bvector<>* bv2)
{
    const unsigned fill_factor1 = 2500;
    const unsigned fill_factor2 = 150;
    for (unsigned i = 0; i < MAX_VALUE; ++i)
        if (unsigned(rand()) % fill_factor1)
            bv1->set(i);
    for (unsigned i = 0; i < MAX_VALUE; i+=fill_factor2)
        bv2->set(i);
}

static
void print_statistics(const bm::bvector<>& bv)
{
    bm::bvector<>::statistics st;
    bv.calc_stat(&st);

    cout << "Bit blocks        :" << st.bit_blocks << endl;
    cout << "GAP blocks        :" << st.gap_blocks << endl;
    cout << "Memory used       :" << st.memory_used << endl;
    cout << "Memory overhead   :" << st.gap_cap_overhead << endl;
    cout << "Max.serialize mem :" << st.max_serialize_mem << endl << endl;;
}


int main(void)
{
    try
    {
        bm::bvector<>   bv1;
        bm::bvector<>   bv2;

        fill_bvector(&bv1, &bv2);  // Fill bvectors (dense and sparse)

        cout << "Statistics before memory optimization" << endl;
        cout << "-------------------------------------" << endl << endl;
        print_statistics(bv1);
        print_statistics(bv2);


        {
            // optimization scratch buffer on stack to avoid dynamic allocation
            BM_DECLARE_TEMP_BLOCK(tb)

            bv1.optimize(tb);
            bv2.optimize(tb);
        }

        cout << "Statistics after memory optimization" << endl;
        cout << "-------------------------------------" << endl;
        print_statistics(bv1);
        print_statistics(bv2);

        // two ways to turn vector read-only:
        // 1. construct a new vector as a READ_ONLY copy (drop the old one)
        // 2. Use bm::bvector<>::freeze()

        bm::bvector<>   bv1_ro(bv1, bm::finalization::READONLY);
        assert(bv1_ro.is_ro());
        bv2.freeze();
        cout << bv2.is_ro() << endl; // 1 - yes it is now READ_ONLY


        cout << "Statistics after freezing" << endl;
        cout << "-------------------------------------" << endl;
        print_statistics(bv1_ro);
        print_statistics(bv2);

        // interface-wise, immutable bit-vectors remain just the same bvector<>
        // bit you CANNOT call non-cost modification functions on it
        // (it would be undefined behavior) don't shoot yourself in the leg
        //

        // READ-ONLY status is preserved over the copy construction and
        // assignment
        {
        // copy construct of RO bit-vector (result is read-only vector vector)
        bm::bvector<>   bv11(bv1_ro);
        assert(bv11.is_ro());
        bool eq = bv11.equal(bv1_ro);
        assert(eq); (void) eq;
        }


        // if you can rewet it back here is how to do it
        // construct bv11 as READ-WRITE from READONLY
        {
        assert(bv2.is_ro()); // bv2 is immutable

        bm::bvector<>   bv21(bv2, bm::finalization::READWRITE);
        bv2.swap(bv21);
        assert(!bv2.is_ro()); // bv2 is now a mutable vector again

        bv2.set(0); // we can do chnage operations at this point
        }




    }
    catch(std::exception& ex)
    {
        std::cerr << ex.what() << std::endl;
    }

    return 0;
}

