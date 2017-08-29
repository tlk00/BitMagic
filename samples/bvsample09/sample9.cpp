/*
Copyright(c) 2002-2005 Anatoliy Kuznetsov(anatoliy_kuznetsov at yahoo.com)

Permission is hereby granted, free of charge, to any person 
obtaining a copy of this software and associated documentation 
files (the "Software"), to deal in the Software without restriction, 
including without limitation the rights to use, copy, modify, merge, 
publish, distribute, sublicense, and/or sell copies of the Software, 
and to permit persons to whom the Software is furnished to do so, 
subject to the following conditions:

The above copyright notice and this permission notice shall be included 
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, 
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, 
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR 
OTHER DEALINGS IN THE SOFTWARE.
*/

/** \example sample9.cpp

    Example demonstrates binary distance metrics.
 
  \sa bm::count_and 
  \sa bm::count_xor
  \sa bm::count_sub
  \sa bm::count_or
  \sa bm::distance_operation

   For more information please visit:  http://bmagic.sourceforge.net

*/

#include <iostream>

#include "bm.h"
#include "bmalgo.h"


using namespace std;


int main(void)
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
    
    unsigned hamming = bm::count_xor(bv1, bv2);
    
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
    
    return 0;
}

