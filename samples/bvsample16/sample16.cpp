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

/** \example sample16.cpp
Example for finding first and last bits in bit-vector (dynamic range).

Ranges of bit-vectors can be used to find probability of intersection.
For instance, in some corner cases AND product can be predicted empty if vectors
belong to different ranges.

    @sa bm::aggregator
*/

/*! \file sample16.cpp
    \brief Example: how to use bm::aggregator<> for logical operations

    bm::aggregator<> uses cache blocking techniques and bandwidth optimizations
    to do logical operations (OR, AND, AND-SUB) faster, than if we do it by
    combining bit-vectors one by one, sequentially.
*/

#include <stdlib.h>
#include <iostream>
#include <vector>
#include <memory>

#include "bm.h"
#include "bmaggregator.h"
#include "bmundef.h" /* clear the pre-proc defines from BM */

// Utility template function used to print container
template<class T> void PrintContainer(T first, T last)
{
    if (first == last)
        std::cout << "<EMPTY SET>";
    else
        for (; first != last; ++first)
            std::cout << *first << ";";
    std::cout << std::endl;
}


const unsigned max_vectors = 10;

int main(void)
{
    try
    {
        // declare standalone aggregator for logical operations
        bm::aggregator<bm::bvector<> > agg;
        
        std::cout << "AGRUMENT (GROUP 0) SETS:" << std::endl;
        // make vector of bit-vectors, set some bits
        std::vector<std::unique_ptr<bm::bvector<> > > vect;
        for (unsigned i = 0; i < max_vectors; ++i)
        {
            std::unique_ptr<bm::bvector<>> bv(new bm::bvector<>());
            bv->set(i);
            bv->set(i+1);
            bv->set(10000);
            bv->set(20000);
            PrintContainer(bv->first(), bv->end());
            bv->optimize();
            vect.push_back(std::move(bv));
        }
        std::cout << std::endl;
        
        try
        {
            // in a loop we add all agruments to the aggregator
            // (aggregator does not take ownership of pointers it receives)
            //
            for (unsigned i = 0; i < vect.size(); ++i)
            {
                agg.add(vect[i].get());
            }
            
            bm::bvector<> bv_res; // target vector for aggregation

            agg.combine_or(bv_res);  // perform logical OR on a vector group
            
            std::cout << "OR:" << std::endl;
            PrintContainer(bv_res.first(), bv_res.end()); // 0, 1, 2, 3, 4... 10000, 20000
            
            // since we did not call bm::aggregator::reset() the vector group
            // still remains attached
            
            std::cout << "AND:" << std::endl;
            agg.combine_and(bv_res);  // AND on the same vector group
            PrintContainer(bv_res.first(), bv_res.end()); // 10000, 20000

            agg.reset(); // reset the aggregator
            
            // fused logical AND MINUS example
            // AND(vector-group-0) SUBstract (AND NOT) vector-group-1
            for (unsigned i = 0; i < vect.size(); ++i)
            {
                const unsigned group0 = 0; // vector group 0
                agg.add(vect[i].get(), group0);
            }
            // Note that we do not set vector group 1 yet, so operation just
            // runs as a regular AND MINUS an empty-set

            agg.combine_and_sub(bv_res);  // AND-SUB vector group1 and group2
            
            std::cout << "AND-SUB(empty):" << std::endl;
            PrintContainer(bv_res.first(), bv_res.end()); // 10000, 20000
            
            bm::bvector<> bv_not { 10, 10000 };
            agg.add(&bv_not, 1); // add to vector-group-1

            agg.combine_and_sub(bv_res);  // AND-SUB vector group0 minus group1
            
            std::cout << "AND-SUB:" << std::endl;
            PrintContainer(bv_res.first(), bv_res.end()); // 20000
            
            agg.reset();

        }
        catch(std::exception& ex)
        {
            std::cerr << ex.what() << std::endl;
            agg.reset();  // reset
            throw;
        }

    }
    catch(std::exception& ex)
    {
        std::cerr << ex.what() << std::endl;
        return 1;
    }
    

    return 0;
}


