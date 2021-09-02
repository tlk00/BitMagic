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

/** \example bv3vlogic.cpp
     Example demonstrates variety three-valued logic.
     https://en.wikipedia.org/wiki/Three-valued_logic
 
     <a href="https://en.wikipedia.org/wiki/Three-valued_logic">Three-valued logic</a>
*/

/*! \file bv3vlogic.cpp
    \brief Example: Kleene algebra operations

    BitMagic implements 3-value (Kleene) logic using two separate bit-vectors:
    bit-vector of values and bit-vector of knowns.

    bit-vector of vlaues contains 1s in the positions of true
    bit-vector of NULLs contains 1s where values in known (or set)

    The convention is that unknown elements (0s in the NULL vector) MUST NOT
    have 1s in the corresponding positions of the value bit-vector so
    "unknown TRUE" is not a correct situation.
*/


#include <iostream>
#include <vector>
#include <cassert>

#include "bm.h"
#include "bm3vl.h"
#include "bmundef.h" /* clear the pre-proc defines from BM */

using namespace std;

/**
    Print 3-value vector
 */
static
void PrintKleeneVector(const bm::bvector<>& bv_v, const bm::bvector<>& bv_null)
{
    bm::bvector<>::enumerator en_n = bv_null.first();
    auto prev = *en_n;
    if (prev > 0)
        prev = 0;

    for ( ;en_n.valid(); ++en_n)
    {
        auto curr = *en_n;
        for (auto i = prev; i < curr; ++i)
            cout << i << ": NULL" << endl;
        bool v = bv_v.test(curr);
        cout << curr << ": " << (v ? "true" : "false") << endl;
        prev = curr + 1;
    } // for en_n

    cout << endl;

}

/**
    This demo shows how to use bm::set_value_kleene and bm::get_value_kleene
    functions to set values into pair of vectors (value vector and knowns vector)

    Kleene algebra operates on 3 values, which are by convention read as:
    -1 (known false), 0 (unknown), 1 (known true).

    bm::set_value_kleene takes a pair of vectors, position and an int value
    to set value in a pair of bit-vectors representing value and knowns

    Please note that this is the easy but relatively slow method to init the
    vectors since because it uses random access initialization.
 */
static
void Set3VL_ValueDemo()
{
    bm::bvector<> bv_v;     // (true/false) values bit-vector
    bm::bvector<> bv_null;  // (known/unknown (or NULL) values bit-vector

    int v = 0; // start with the unknown
    for (unsigned i = 0; i < 10; ++i)
    {
        bm::set_value_kleene(bv_v, bv_null, i, v);
        auto v1 = bm::get_value_kleene(bv_v, bv_null, i);
        assert(v == v1); (void) v1;
        v += 1;
        if (v > 1)
            v = -1;
    }

    BM_DECLARE_TEMP_BLOCK(tb)
    bv_v.optimize(tb);
    bv_null.optimize(tb);

    PrintKleeneVector(bv_v, bv_null);
}

/**
    Faster way to initialize Kleene bit-vectors via bulk_insert_iterator
 */
static
void Set3VL_ValueDemo2()
{
    bm::bvector<> bv_v;     // (true/false) values bit-vector
    bm::bvector<> bv_null;  // (known/unknown (or NULL) values bit-vector

    // use insert iterators to load the vectors
    {
        bm::bvector<>::bulk_insert_iterator iit_v = bv_v.inserter();
        bm::bvector<>::bulk_insert_iterator iit_n = bv_null.inserter();

        for (unsigned i = 0; i < 13; i+=3)
        {
            if (i & 1) // add only true values as indexes of set bits via insert iterator
            {
                iit_v = i;
            }
            // set the known bit for both true AND false values
            iit_n = i;
        }
        // flush the insert iterators to empty the temp.buffers
        //    it is best to do it explicitly (it can throw exceptions)
        iit_v.flush();
        iit_n.flush();
    }

    // init used to guarantee that value bit-vector would not contain
    // "unknown" true values, which is a requirement of BM library to be able to
    // correct 3-value logical ANDs and ORs
    //
    // If you are absolutely sure that you did not set true values for unknowns
    // then you do not need this call
    //
    bm::init_kleene(bv_v, bv_null);

    // optimize bit-vectors after initalization
    //
    BM_DECLARE_TEMP_BLOCK(tb)
    bv_v.optimize(tb);
    bv_null.optimize(tb);

    PrintKleeneVector(bv_v, bv_null);
}

/**
    Generate Kleene vector (as two bit-vectors)
*/
static
void GenerateDemoVector(bm::bvector<>& bv_v, bm::bvector<>& bv_null)
{
    int v = 0; // start with the unknown
    for (unsigned i = 0; i < 10; ++i)
    {
        bm::set_value_kleene(bv_v, bv_null, i, v);
        v += 1;
        if (v > 1)
            v = -1;
    } // for
    // optimize bit-vectors after initalization
    //
    BM_DECLARE_TEMP_BLOCK(tb)
    bv_v.optimize(tb);
    bv_null.optimize(tb);
}

/**
    Demo for 3-value logic (Kleene) NOT
 */
static
void Set3VL_InvertDemo()
{
    bm::bvector<> bv_v;     // (true/false) values bit-vector
    bm::bvector<> bv_null;  // (known/unknown (or NULL) values bit-vector

    GenerateDemoVector(bv_v, bv_null);

    cout << "Input vector:" << endl;
    PrintKleeneVector(bv_v, bv_null);

    bm::invert_kleene(bv_v, bv_null); // 3-value logic NOT

    cout << "Inverted vector:" << endl;
    PrintKleeneVector(bv_v, bv_null);
}


/**
    Demo for 3-value logic (Kleene) AND

    Kleene algebra AND
    produces known FALSE when known FALSE meets UNKNOWN (false)
 */
static
void Set3VL_AndDemo()
{
    bm::bvector<> bv_v1;     // (true/false) values bit-vector
    bm::bvector<> bv_null1;  // (known/unknown (or NULL) values bit-vector

    GenerateDemoVector(bv_v1, bv_null1);

    bm::bvector<> bv_v2;     // (true/false) values bit-vector
    bm::bvector<> bv_null2;  // (known/unknown (or NULL) values bit-vector

    bm::set_value_kleene(bv_v2, bv_null2, 0, 0); // idx = 0 (unknown)
    bm::set_value_kleene(bv_v2, bv_null2, 1, 1); // idx = 1 (known true)
    bm::set_value_kleene(bv_v2, bv_null2, 2, -1); // idx = 2 (known false)


    bm::bvector<> bv_v_t, bv_null_t;
    // bv_v_t := bv_v2 & bv_v1
    bm::and_kleene(bv_v_t, bv_null_t, bv_v2, bv_null2, bv_v1, bv_null1); // 3-value logic AND

    // bv_v2 and bv_null2 are modified in place:
    // bv_v2 &= bv_v1
    bm::and_kleene(bv_v2, bv_null2, bv_v1, bv_null1); // 3-value logic AND

    bool b = bv_v_t.equal(bv_v2);
    assert(b);
    b = bv_null_t.equal(bv_null2);
    assert(b);

    cout << "AND vector:" << endl;
    PrintKleeneVector(bv_v2, bv_null2);

}

/**
    Demo for 3-value logic (Kleene) OR

    Kleene algebra OR
    produces known TRUE when known TRUE meets UNKNOWN (false)
 */
static
void Set3VL_ORDemo()
{
    bm::bvector<> bv_v1;     // (true/false) values bit-vector
    bm::bvector<> bv_null1;  // (known/unknown (or NULL) values bit-vector

    GenerateDemoVector(bv_v1, bv_null1);

    bm::bvector<> bv_v2;     // (true/false) values bit-vector
    bm::bvector<> bv_null2;  // (known/unknown (or NULL) values bit-vector

    bm::set_value_kleene(bv_v2, bv_null2, 0, 1); // idx = 0 (known true)
    bm::set_value_kleene(bv_v2, bv_null2, 1, 0); // idx = 1 (NULL)
    bm::set_value_kleene(bv_v2, bv_null2, 2, -1); // idx = 2 (known false)

    bm::bvector<> bv_v_t, bv_null_t;
    // bv_v_t := bv_v2 | bv_v1
    bm::or_kleene(bv_v_t, bv_null_t, bv_v2, bv_null2, bv_v1, bv_null1); // 3-value logic AND

    // bv_v2 and bv_null2 are modified in place:
    // bv_v2 |= bv_v1
    bm::or_kleene(bv_v2, bv_null2, bv_v1, bv_null1); // 3-value logic OR

    bool b = bv_v_t.equal(bv_v2);
    assert(b);
    b = bv_null_t.equal(bv_null2);
    assert(b);


    cout << "OR vector:" << endl;
    PrintKleeneVector(bv_v2, bv_null2);

}



int main(void)
{
    try
    {
        cout << endl << "3VL Set values:" << endl << endl;
        Set3VL_ValueDemo();
        Set3VL_ValueDemo2();

        cout << endl << "3VL Invert vector:" << endl << endl;
        Set3VL_InvertDemo();

        cout << endl << "3VL AND:" << endl << endl;
        Set3VL_AndDemo();

        cout << endl << "3VL OR:" << endl << endl;
        Set3VL_ORDemo();
    }
    catch(std::exception& ex)
    {
        std::cerr << ex.what() << std::endl;
    }

    return 0;
}
