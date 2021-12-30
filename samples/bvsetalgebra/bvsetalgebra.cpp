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

/** \example bvsetalgebra.cpp
     Example demonstrates variety of algebra of sets operations.
     http://bitmagic.io/set-algebra.html
 
     <a href="http://bitmagic.io/set-algebra.html">Algebra Of Sets</a>

  \sa bvector

  \sa bm::bvector::bit_or
  \sa bm::bvector::bit_and
  \sa bm::bvector::bit_xor
  \sa bm::bvector::bit_sub

  \sa bm::aggregator
  \sa bm::operation_deserializer

  \sa bm::combine_and
  \sa bm::combine_and_sorted
  \sa bm::combine_sub
  \sa bm::combine_or
  \sa bm::combine_xor

  \sa sample7.cpp

*/

/*! \file bvsetalgebra.cpp
    \brief Example: algebra of sets operations
*/


#include <assert.h>
#include <iostream>
#include <vector>

#include "bm.h"
#include "bmalgo.h"
#include "bmserial.h"
#include "bmaggregator.h"
#include "bmundef.h" /* clear the pre-proc defines from BM */

using namespace std;

// utility function to print a set
static
void print_bvector(const bm::bvector<>& bv)
{
    bm::bvector<>::size_type cnt = 0;
    bm::bvector<>::enumerator en = bv.first();
    for (; en.valid() && cnt < 10; ++en, ++cnt)
        cout << *en << ", ";
    if (cnt == 10)
        cout << " ...";
    cout << "(size = "<< bv.size() << ")" << endl;
}

// utility function to create serialized bit-vector BLOB
static
void make_BLOB(vector<unsigned char>& target_buf, bm::bvector<>& bv)
{
    BM_DECLARE_TEMP_BLOCK(tb)
    bm::serializer<bm::bvector<> > bvs(tb);
    bvs.set_compression_level(4);
    
    bv.optimize(tb, bm::bvector<>::opt_compress); // memory compression

    bm::serializer<bm::bvector<> >::buffer sbuf;
    bvs.serialize(bv, sbuf, 0);
    target_buf.resize(sbuf.size());
    ::memcpy(target_buf.data(), sbuf.buf(), sbuf.size());
}


// -------------------------------------------------------------
// Demo for Set Union (OR) operations
//
static
void DemoOR()
{
    typedef bm::bvector<>::size_type size_type;
    bm::operation_deserializer<bm::bvector<> > od;

    // bit-vector set union operation: bv_A |= bv_B
    {
        bm::bvector<>   bv_A { 1, 2, 3 };
        bm::bvector<>   bv_B { 1, 2, 4 };
        bv_A.bit_or(bv_B);
        
        print_bvector(bv_A); // 1, 2, 3, 4
    }
    // same, but sizes are set, observe size gets extended up
    {
        bm::bvector<>   bv_A { 1, 2, 3 };
        bm::bvector<>   bv_B { 1, 2, 4 };
        bv_A.resize(5);
        bv_B.resize(10);

        bv_A.bit_or(bv_B);
        
        print_bvector(bv_A); // 1, 2, 3, 4 (size = 10)
    }
    // 3-operand OR: bv_T = bv_A | bv_B
    {
        bm::bvector<>   bv_T;
        bm::bvector<>   bv_A { 1, 2, 3 };
        bm::bvector<>   bv_B { 1, 2, 4 };

        bv_T.bit_or(bv_A, bv_B, bm::bvector<>::opt_compress);
        
        print_bvector(bv_T); // 1, 2, 3, 4 (size = 10)
    }
    
    
    // merge operation is a logical equivalent of OR
    // except it can destroy the source vector to borrow memory blocks from it
    // (this is faster, especially in multi-threaded cases)
    {
        bm::bvector<>   bv_A { 1, 2, 3 };
        bm::bvector<>   bv_B { 1, 2, 4 };
        bv_A.merge(bv_B);
        
        print_bvector(bv_A); // 1, 2, 3, 4 (size = 10)
    }

    // bit-vector set union operation (opcode interpeter mode)
    // maybe useful for building query interpetors
    {
        bm::bvector<>   bv_A { 1, 2, 3 };
        bm::bvector<>   bv_B { 1, 2, 4 };
        bv_A.combine_operation(bv_B, bm::BM_OR);
        
        print_bvector(bv_A); // 1, 2, 3, 4
    }

    // Set union between bit-vector and STL container
    {
        bm::bvector<>      bv_A { 1, 2, 3 };
        vector<size_type>   vect_B { 1, 2, 4 };
        
        bm::combine_or(bv_A, vect_B.begin(), vect_B.end());
        print_bvector(bv_A); // 1, 2, 3, 4
    }
    
    // Set union between bit-vector and C-array.
    // This tends to be faster then "combine_or()" especially on sorted vectors
    // and in SIMD enabled configurations
    {
        bm::bvector<>      bv_A { 1, 2, 3 };
        vector<size_type>   vect_B { 1, 2, 4 };
        
        const size_type* arr = &vect_B[0];
        bv_A.set(arr, unsigned(vect_B.size()), bm::BM_SORTED); // sorted - fastest
        print_bvector(bv_A); // 1, 2, 3, 4
    }

    // Set union between bit-vector and a serialized bit-vector BLOB
    // (created on the fly)
    {
        bm::bvector<>   bv_A { 1, 2, 3 };
        vector<unsigned char> blob;
        {
            bm::bvector<>   bv_B { 1, 2, 4 };
            make_BLOB(blob, bv_B);
        }
        od.deserialize(bv_A, blob.data(), bm::set_OR);
        print_bvector(bv_A); // 1, 2, 3, 4
    }
    
    // Union of many sets with bm::aggegator<>
    // target := A OR B OR C
    //
    // This method is best when we have multiple vectors at hands, aggregator
    // is capable of doing it faster, than pair by pair OR
    {
        bm::bvector<>    bv_T; // target vector
        
        bm::bvector<>    bv_A { 1, 2 };
        bm::bvector<>    bv_B { 2, 3 };
        bm::bvector<>    bv_C { 3, 4 };
        
        bm::aggregator<bm::bvector<> > agg;
        agg.set_optimization(); // perform on-the-fly optimization of result
        
        // attach vectors to group 0 for OR operation
        agg.add(&bv_A);
        agg.add(&bv_B);
        agg.add(&bv_C);
        
        agg.combine_or(bv_T);
        
        agg.reset(); // reset the aggregator parameters
        
        print_bvector(bv_T); // 1, 2, 3, 4
    }
    
}


// -------------------------------------------------------------
// Demo for Set Intersect (AND) operations
//
static
void DemoAND()
{
    typedef bm::bvector<>::size_type size_type;
    bm::operation_deserializer<bm::bvector<> > od;

    // bit-vector set intersect operation: bv_A &= bv_B
    {
        bm::bvector<>   bv_A { 1, 2, 3 };
        bm::bvector<>   bv_B { 1, 2, 4 };
        bv_A.bit_and(bv_B);
        
        print_bvector(bv_A); // 1, 2
    }
    // same, but sizes are set, observe size gets extended up
    {
        bm::bvector<>   bv_A { 1, 2, 3 };
        bm::bvector<>   bv_B { 1, 2, 4 };
        bv_A.resize(5);
        bv_B.resize(10);

        bv_A.bit_and(bv_B);
        
        print_bvector(bv_A); // 1, 2 (size = 10)
    }
    // 3-operand AND: bv_T = bv_A & bv_B
    {
        bm::bvector<>   bv_T;
        bm::bvector<>   bv_A { 1, 2, 3 };
        bm::bvector<>   bv_B { 1, 2, 4 };
        bv_T.bit_and(bv_A, bv_B, bm::bvector<>::opt_compress);
        
        print_bvector(bv_T); // 1, 2
    }

    // bit-vector set intesect operation (opcode interpeter mode)
    // maybe useful for building query interpetors
    {
        bm::bvector<>   bv_A { 1, 2, 3 };
        bm::bvector<>   bv_B { 1, 2, 4 };
        bv_A.combine_operation(bv_B, bm::BM_AND);
        
        print_bvector(bv_A); // 1, 2
    }

    // Set Intersect between bit-vector and STL container
    {
        bm::bvector<>      bv_A { 1, 2, 3 };
        vector<unsigned>   vect_B { 1, 2, 4 };
        
        bm::combine_and(bv_A, vect_B.begin(), vect_B.end());
        print_bvector(bv_A); // 1, 2
    }
    
    // Set Intersect between bit-vector and C-array.
    // This may be faster then "combine_and()" especially on sorted vectors
    {
        bm::bvector<>      bv_A { 1, 2, 3 };
        vector<size_type>   vect_B { 1, 2, 4 };
        
        const size_type* arr = &vect_B[0];
        bv_A.keep(arr, size_type(vect_B.size()), bm::BM_SORTED); // sorted - fastest
        print_bvector(bv_A); // 1, 2
    }

    // Set Intersect between bit-vector and a serialized bit-vector BLOB
    //
    {
        bm::bvector<>   bv_A { 1, 2, 3 };
        vector<unsigned char> blob;
        {
            bm::bvector<>   bv_B { 1, 2, 4 };
            make_BLOB(blob, bv_B);
        }
        od.deserialize(bv_A, blob.data(), bm::set_AND);
        print_bvector(bv_A); // 1, 2
    }
    
    // Intersection of many sets with bm::aggegator<> (find common subset)
    // target := A AND B AND C
    //
    // This method is best when we have multiple vectors at hands, aggregator
    // is capable of doing it faster, than pair by pair AND
    {
        bm::bvector<>    bv_T; // target vector
        
        bm::bvector<>    bv_A { 1, 2 };
        bm::bvector<>    bv_B { 1, 2, 3 };
        bm::bvector<>    bv_C { 1, 2, 3, 4 };
        
        bm::aggregator<bm::bvector<> > agg;
        agg.set_optimization(); // perform on-the-fly optimization of result
        
        // attach vectors to group 0 for OR operation
        agg.add(&bv_A);
        agg.add(&bv_B);
        agg.add(&bv_C);
        
        agg.combine_and(bv_T);
        
        agg.reset(); // reset the aggregator parameters
        
        print_bvector(bv_T); // 1, 2
    }
}

// -------------------------------------------------------------
// Demo for XOR operations
//
static
void DemoXOR()
{
    bm::operation_deserializer<bm::bvector<> > od;

    // bit-vector xor operation: bv_A ^= bv_B
    {
        bm::bvector<>   bv_A { 1, 2, 3 };
        bm::bvector<>   bv_B { 1, 2, 4 };
        bv_A.bit_xor(bv_B);
        
        print_bvector(bv_A); // 3, 4
    }
    // same, but sizes are set, observe size gets extended up
    {
        bm::bvector<>   bv_A { 1, 2, 3 };
        bm::bvector<>   bv_B { 1, 2, 4 };
        bv_A.resize(5);
        bv_B.resize(10);

        bv_A.bit_xor(bv_B);
        
        print_bvector(bv_A); // 3, 4 (size = 10)
    }
    // 3-operand XOR: bv_T = bv_A ^ bv_B
    {
        bm::bvector<>   bv_T;
        bm::bvector<>   bv_A { 1, 2, 3 };
        bm::bvector<>   bv_B { 1, 2, 4 };
        bv_T.bit_xor(bv_A, bv_B, bm::bvector<>::opt_compress);
        
        print_bvector(bv_T); // 3, 4
    }

    // bit-vector xor operation (opcode interpeter mode)
    // maybe useful for building query interpetors
    {
        bm::bvector<>   bv_A { 1, 2, 3 };
        bm::bvector<>   bv_B { 1, 2, 4 };
        bv_A.combine_operation(bv_B, bm::BM_XOR);
        
        print_bvector(bv_A); // 3, 4
    }

    // xor between bit-vector and STL container
    {
        bm::bvector<>      bv_A { 1, 2, 3 };
        vector<unsigned>   vect_B { 1, 2, 4 };
        
        bm::combine_xor(bv_A, vect_B.begin(), vect_B.end());
        print_bvector(bv_A); // 3, 4
    }

    // xor between bit-vector and a serialized bit-vector BLOB
    //
    {
        bm::bvector<>   bv_A { 1, 2, 3 };
        vector<unsigned char> blob;
        {
            bm::bvector<>   bv_B { 1, 2, 4 };
            make_BLOB(blob, bv_B);
        }
        od.deserialize(bv_A, blob.data(), bm::set_XOR);
        print_bvector(bv_A); // 3, 4
    }
}


// -------------------------------------------------------------
// Demo for Set Substract (AND NOT) operations
//
static
void DemoSUB()
{
    typedef bm::bvector<>::size_type size_type;
    bm::operation_deserializer<bm::bvector<> > od;

    // bit-vector set union operation: bv_A -= bv_B
    {
        bm::bvector<>   bv_A { 1, 2, 3 };
        bm::bvector<>   bv_B { 1, 2, 4 };
        bv_A.bit_sub(bv_B);
        
        print_bvector(bv_A); // 3
    }
    // same, but sizes are set, observe size gets extended up
    {
        bm::bvector<>   bv_A { 1, 2, 3 };
        bm::bvector<>   bv_B { 1, 2, 4 };
        bv_A.resize(5);
        bv_B.resize(10);

        bv_A.bit_sub(bv_B);
        
        print_bvector(bv_A); // 3 (size = 10)
    }
    
    // 3-operand SUB: bv_T = bv_A - bv_B
    {
        bm::bvector<>   bv_T;
        bm::bvector<>   bv_A { 1, 2, 3 };
        bm::bvector<>   bv_B { 1, 2, 4 };
        bv_T.bit_sub(bv_A, bv_B, bm::bvector<>::opt_compress);
        
        print_bvector(bv_T); // 3
    }

    // bit-vector minus operation (opcode interpeter mode)
    // maybe useful for building query interpetors
    {
        bm::bvector<>   bv_A { 1, 2, 3 };
        bm::bvector<>   bv_B { 1, 2, 4 };
        bv_A.combine_operation(bv_B, bm::BM_SUB);
        
        print_bvector(bv_A); // 3
    }

    // and not between bit-vector and STL container
    {
        bm::bvector<>      bv_A { 1, 2, 3 };
        vector<size_type>   vect_B { 1, 2, 4 };
        
        bm::combine_sub(bv_A, vect_B.begin(), vect_B.end());
        print_bvector(bv_A); // 3
    }

    // Set Intersect between bit-vector and C-array.
    // This may be faster then "combine_and()" especially on sorted vectors
    {
        bm::bvector<>      bv_A { 1, 2, 3 };
        vector<size_type>   vect_B { 1, 2, 4 };
        
        const size_type* arr = &vect_B[0];
        bv_A.clear(arr, unsigned(vect_B.size()), bm::BM_SORTED); // sorted - fastest
        print_bvector(bv_A); // 3
    }

    // Set union between bit-vector and a serialized bit-vector BLOB
    //
    {
        bm::bvector<>   bv_A { 1, 2, 3 };
        vector<unsigned char> blob;
        {
            bm::bvector<>   bv_B { 1, 2, 4 };
            make_BLOB(blob, bv_B);
        }
        od.deserialize(bv_A, blob.data(), bm::set_SUB);
        print_bvector(bv_A); // 3
    }

    // Subtraction of many sets with bm::aggegator<>
    // target := (target SUB A) OR (target SUB B) OR (target SUB C)
    //
    {
        bm::bvector<>    bv_T; // target vector
        
        bm::bvector<>    bv_A { 1, 2, 3, 4 };
        bm::bvector<>    bv_B { 1, 2 };
        bm::bvector<>    bv_C { 1, 2, 4 };
        
        bm::aggregator<bm::bvector<> > agg;
        agg.set_optimization(); // perform on-the-fly optimization of result
        
        // here we are really using AND SUB operation
        // where group 0 is all ANDed and group 1 SUBtracted from the result
        // group 1 is only 1 vector, so AND part will be no-op
        //
        agg.add(&bv_A, 0); // add to group 0 (subtraction source)
        
        agg.add(&bv_B, 1); // add to group 1 (subtraction arguments)
        agg.add(&bv_C, 1);
        
        agg.combine_and_sub(bv_T);
        
        agg.reset(); // reset the aggregator parameters
        
        print_bvector(bv_T); // 3
    }
}

// -------------------------------------------------------------
// Demo for Set Intersect (AND) operations fused with OR
//
static
void DemoAND_OR()
{
    // bit-vector set intersect with fused OR: bv_R |= bv_A & bv_B
    {
        bm::bvector<>   bv_R { 0, 5 };

        bm::bvector<>   bv_A { 1, 2, 3 };
        bm::bvector<>   bv_B { 1, 2, 4 };
        bv_R.bit_or_and(bv_A, bv_B);

        print_bvector(bv_R); // 0, 1, 2, 5
    }

    // same, but sizes are set, observe size gets extended up
    {
        bm::bvector<>   bv_R { 0, 5 };
        bv_R.resize(12);

        bm::bvector<>   bv_A { 1, 2, 3 };
        bm::bvector<>   bv_B { 1, 2, 4 };
        bv_A.resize(5);
        bv_B.resize(10);

        bv_R.bit_or_and(bv_A, bv_B);

        print_bvector(bv_R); // 0, 1, 2, 5 (size = 12)
    }

    // AND OR with aggregator pipeline

    bm::aggregator<bm::bvector<> > agg;

    // Aggregator pipeline can efficiently run a group of AND (or AND-SUB)
    // operations on bit-vectors, assuming that groups of ANDs are comprised
    // of essentially the same super-set of vectors which are going to be
    // reused (cached) by the CPU, thus run more efficintly
    //
    // In this case pipeline is configured to ignore the sub-results of AND
    // operations and just UNION all set-results into a target bit-vector
    // (logical OR)
    //
    // The formula here is:
    // bv_TARGET = bv_TARGET OR (BV1 AND BV2) OR (BV1 & BV3) & (BV2 & BV3)
    //
    // This pipeline configuration avoids allocations of temporary vectors
    // and CPU cache-friendly (thus fast)
    //
    bm::aggregator<bm::bvector<> >::pipeline<bm::agg_opt_disable_bvects_and_counts> agg_pipe;

    bm::bvector<>   bv_T { 0, 65536 };

    bm::bvector<>   bv_A { 1, 2, 3 };
    bm::bvector<>   bv_B { 1, 2  };
    bm::bvector<>   bv_C { 2, 3, 256  };

    // feed input queries into pipeline
    {
        bm::aggregator<bm::bvector<> >::arg_groups* args;

        // Define AND group 1 as: (A & B)
        args = agg_pipe.add();
        args->add(&bv_A, 0); // 0 - means AND arg group
        args->add(&bv_B, 0);

        // Define AND group 2 as: (B & C)
        args = agg_pipe.add();
        args->add(&bv_B, 0);
        args->add(&bv_C, 0);
    }

    // T |= (A & B) | (B & C)
    //
    agg_pipe.set_or_target(&bv_T);

    agg_pipe.complete(); // finish the pipeline configuration

    // We use AND-MINUS operation but since we did not add any subtraction
    // vectors the result will be just AND
    //
    agg.combine_and_sub(agg_pipe); // T |= (A & B) | (B & C)

    print_bvector(bv_T); // 0, 1, 2, 65536
}



// -------------------------------------------------------------
// Demo for Set Invert (NOT)
//
static
void DemoINV()
{
    // bit-vector invert operation
    // by default it inverts the whole 32-bit space
    {
        bm::bvector<>   bv_A { 4, 5, 6  };
        bv_A.invert();
        
        print_bvector(bv_A); // 0, 1, 2, 3, 7 ...
    }
    
    // bit-vector invert operation
    // it is size bound, inverts within set limits
    {
        bm::bvector<>   bv_A { 4, 5, 6  };
        bv_A.resize(7);
        bv_A.invert();
        
        print_bvector(bv_A); // 0, 1, 2, 3, size = 7
    }
}


// -------------------------------------------------------------
// Demo for AND-SUB
//  AND-SUB implements a search pattern "all this but not that"
//
static
void DemoAND_SUB()
{
    bm::aggregator<bm::bvector<> > agg;

    // Operation on two groups of vectors using aggregator
    // 1. Group 0 - find common subset (Set Intersect / AND)
    // 2. Group 1 - find union (OR) of the group and SUBtract it from #1
    //
    // target := (A AND D AND ...) AND NOT (B OR C OR ...)
    //
    bm::bvector<>    bv_A { 1, 2, 3, 4 };
    bm::bvector<>    bv_B { 1, 2 };
    bm::bvector<>    bv_C { 1, 2, 4 };
    bm::bvector<>    bv_D { 0, 2, 3, 4, 5 };

    {
        bm::bvector<>    bv_T; // target vector

        agg.set_optimization(); // perform on-the-fly optimization of result
        
        // here we are really using AND SUB operation
        // where group 0 is all ANDed and group 1 SUBtracted from the result
        //
        agg.add(&bv_A, 0); // add to group 0 for AND
        agg.add(&bv_D, 0); //

        agg.add(&bv_B, 1); // add to group 1 SUB tract from group 0 result
        agg.add(&bv_C, 1);
        
        agg.combine_and_sub(bv_T);
        
        agg.reset(); // reset the aggregator parameters
        
        print_bvector(bv_T); // 3
    }

    cout << endl;

    // Next example shows how to run a group of AND-SUB queries
    // using aggregator::pipeline
    //
    // Pipeline runs multiple search formulas, consiting of essentially the
    // same super-set bit-vectors (vectors are re-used).
    // This reuse opens an oppotinity to run groups of queries faster
    // because of CPU cache reuse
    //
    {
        // default pipeline configuration will produce intermediate results
        // for all pipeline searches
        //
        bm::aggregator<bm::bvector<> >::pipeline<> agg_pipe;
        bm::bvector<> bv_T;

        // Configure 3 different AND-SUB queries based on the same set of vectors
        //
        // Note that if we don't add SUBtraction vectors it will work just as
        // AND formula (AND_SET MINUS empty_set == AND_SET)
        {
            bm::aggregator<bm::bvector<> >::arg_groups* args;

            args = agg_pipe.add();
            // T1 := (A AND B) MINUS D
            args->add(&bv_A, 0); // 0 - means AND arg group
            args->add(&bv_B, 0);
            args->add(&bv_D, 1); // 1 - means SUB arg group

            // T2 := (A AND D) MINUS B
            args = agg_pipe.add();
            args->add(&bv_A, 0);
            args->add(&bv_D, 0);
            args->add(&bv_B, 1); // 1 - means SUB arg group

            // T3 := D MINUS (C OR B)
            // or
            // T3 := (D AND NOT C) OR (D AND NOT B)
            args = agg_pipe.add();
            args->add(&bv_D, 0);
            args->add(&bv_C, 1); // 1 - means SUB arg group
            args->add(&bv_B, 1);
        }
        // this is optional, pipeline can run OR aggregate (UNION ALL)
        // on all searches in the pipeline
        //
        // T = T OR T1 OR T2 OR T3
        //
        agg_pipe.set_or_target(&bv_T);

        agg_pipe.complete(); // finish the configuration

        agg.combine_and_sub(agg_pipe); // run the pipeline search

        // now we iterate the results and print it all
        // we expect 3 output vectors, one for each pipeline AND-SUB group
        //
        auto& res_vect = agg_pipe.get_bv_res_vector();
        assert(res_vect.size()==3); // we expect 3 results
        for (size_t i = 0; i < res_vect.size(); ++i)
        {
            const bm::bvector<>* bv = res_vect[i];
            if (bv)
                print_bvector(*bv); // 1
                                    // 3, 4
                                    // 0, 3, 5
            else
                cout << "Empty result" << endl; // possible outcome as a "not found"
        } // for i

        // OR target contains the UNION of all pipeline AND-SUB searches
        //
        print_bvector(bv_T); // 0, 1, 3, 4, 5,
    }

}



int main(void)
{
    try
    {
        cout << endl << "Set Union (OR) demo" << endl << endl;
        DemoOR();
        
        cout << endl << "Set Intersect (AND) demo" << endl << endl;
        DemoAND();

        cout << endl << "XOR demo" << endl << endl;
        DemoXOR();

        cout << endl << "Set Minus (SUB/AND-NOT) demo" << endl << endl;
        DemoSUB();
        
        cout << endl << "Set Invert (NOT) demo" << endl << endl;
        DemoINV();
        
        cout << endl << "Set AND-SUB demo" << endl << endl;
        DemoAND_SUB();

        cout << endl << "Set AND-OR demo" << endl << endl;
        DemoAND_OR();
    }
    catch(std::exception& ex)
    {
        std::cerr << ex.what() << std::endl;
    }

    return 0;
}
