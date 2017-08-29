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

/** \example sample5.cpp
 Example demonstrates using enumerators - the fastest way to retrieve 
 indexes of 1 bits from the bitvector. This approach works faster than
 get_first/get_next functions.
 
  \sa bm::bvector<>::enumerator 

For more information please visit:  http://bmagic.sourceforge.net

*/

#include <iostream>
#include <algorithm>
#include "bm.h"

using namespace std;

void Print(unsigned n)
{
    cout << n << endl;;
}

int main(void)
{
    bm::bvector<>   bv;    

    bv[10] = true;
    bv[100] = true;
    bv[10000] = true;

    bm::bvector<>::enumerator en = bv.first();
    bm::bvector<>::enumerator en_end = bv.end();

    while (en < en_end)
    {
        cout << *en << endl;
        ++en;  // Fastest way to increment enumerator
    }

    en = bv.first();

    // This is not the fastest way to do the job, because for_each 
    // often will try to calculate difference between iterators,
    // which is expensive for enumerators.
    // But it can be useful for some STL loyal applications. 

    std::for_each(en, en_end, Print);

    return 0;
}
