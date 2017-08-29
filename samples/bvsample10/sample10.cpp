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

/** \example sample10.cpp
  Example of how to get random subset of a bit-vector
 
  For more information please visit:  http://bmagic.sourceforge.net

  \sa bm::random_subset  
 */


#include <iostream>
#include "bm.h"
#include "bmrandom.h"

using namespace std;

template<class T> void PrintContainer(T first, T last)
{
    if (first == last)
        cout << "<EMPTY SET>";
    else
        for(;first != last; ++first)
            cout << *first << ";";
    cout << endl;
}


int main(void)
{
    bm::bvector<>   bv;    

    // -----------------------------------------------
    // set some bits
    //
    unsigned i;
    for (i = 0; i < 30000; i+=10)
    {
        bv.set(i);
    }

    for (i = 300000; i < 400000; i+=100)
    {
        bv.set(i);
    }

    bm::bvector<>   bvsubset1;
    bm::bvector<>   bvsubset2;

    // random sampler instance can be shared between calls
    //
    bm::random_subset<bm::bvector<> > rand_sampler;
    rand_sampler.sample(bvsubset1, bv, 20);
    rand_sampler.sample(bvsubset2, bv, 20);
 
    PrintContainer(bvsubset1.first(), bvsubset1.end());
    cout << endl;
    PrintContainer(bvsubset2.first(), bvsubset2.end());

    return 0;
}

