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

/** \example sample9.cpp

  Example demonstrates binary distance metrics. BM library can compute 
  various distances without materializing a product set.
 
  \sa bm::count_and 
  \sa bm::count_xor
  \sa bm::count_sub
  \sa bm::count_or
  \sa bm::distance_operation

   For more information please visit:  http://bmagic.sourceforge.net
*/

/*! \file sample9.cpp
    \brief Example: bvector<> binary similarity / distance algorithms
*/


#include <iostream>

#include "bm.h"
#include "bmalgo.h"
#include "bmundef.h" /* clear the pre-proc defines from BM */

using namespace std;


int main(void)
{
    try
    {
        bm::bvector<>   bv1;
        bm::bvector<>   bv2;

        bv1[10] = true;
        bv1[100] = true;
        bv1[10000] = true;
        
        bv2[10] = true;
        bv2[100] = true;
        bv2[20000] = true;
        bv2[30000] = true;
        
        // Hamming distance:
        
        auto hamming = bm::count_xor(bv1, bv2);
        
        cout << "Hamming distance = " << hamming << endl;

        // Dice distance using basic distance functions
        
        double dice =
            double(2 * bm::count_and(bv1, bv2))/double(bv1.count() + bv2.count());
        
        cout << "Dice distance = " << dice << endl;
        
        
        // Dice distance, can be computed using "metric pipeline"
        
        bm::distance_metric_descriptor dmd[3];
        dmd[0].metric = bm::COUNT_AND;
        dmd[1].metric = bm::COUNT_A;
        dmd[2].metric = bm::COUNT_B;
                
        bm::distance_operation(bv1, bv2, dmd, dmd+3);
        
        double dice_p =
            double(2 *  dmd[0].result) / double(dmd[1].result + dmd[2].result);
        
        cout << "Dice distance (pipeline) = " << dice_p << endl;
    }
    catch(std::exception& ex)
    {
        std::cerr << ex.what() << std::endl;
        return 1;
    }
    return 0;
}

