/*
Copyright(c) 2002-2017 Anatoliy Kuznetsov(anatoliy_kuznetsov at yahoo.com)

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

You have to explicitly mention BitMagic project in any derivative product,
its WEB Site, published materials, articles or any other work derived from this
project or based on our code or know-how.

For more information please visit:  http://bitmagic.io

*/

/** \example sample5.cpp
 Example demonstrates using enumerators - the fastest way to retrieve 
 indexes of 1 bits from the bitvector. This approach works faster than
 get_first()/get_next() functions.
 
  \sa bm::bvector<>::enumerator 
  \sa bm::bvector<>::first()
  \sa bm::bvector<>::end()
  \sa bm::bvector<>::get_enumerator()
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
    try
    {
        bm::bvector<>   bv;

        bv[10] = true;
        bv[100] = true;
        bv[10000] = true;
        bv[65536] = true;
        bv[65537] = true;
        bv[65538] = true;
        bv[65540] = true;

        bm::bvector<>::enumerator en = bv.first();
        bm::bvector<>::enumerator en_end = bv.end();

        while (en < en_end)
        {
            cout << *en << ", ";
            ++en;  // Fastest way to increment enumerator
        }
        cout << endl;

        en = bv.first();

        // This is not the fastest way to do the job, because for_each
        // often will try to calculate difference between iterators,
        // which is expensive for enumerators.
        // But it can be useful for some STL loyal applications.

        std::for_each(en, en_end, Print);
        cout << endl;

        // example to illustrate random positioning of enumerator 
        // go to a random bit number, enumerator automatically finds the available bit
        //
        en.go_to(65537);
        for (; en.valid(); ++en)
        {
            cout << *en << ", ";
        }
        cout << endl;
    }
    catch(std::exception& ex)
    {
        std::cerr << ex.what() << std::endl;
        return 1;
    }

    return 0;
}
