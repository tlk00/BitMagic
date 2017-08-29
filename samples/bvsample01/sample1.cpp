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

/** \example sample1.cpp
  Example of how to use bvector template class to set
  bits and then retrieve indexes of ON bits.
 
  For more information please visit:  http://bmagic.sourceforge.net

  \sa bm::bvector<>::get_next() 
  \sa bm::bvector<>::get_first() 
  \sa bm::bvector<>::set()
  \sa bm::bvector<>::count() 
  \sa bm::bvector<>::clear()
  
 */

#include <iostream>
#include "bm.h"

using namespace std;

int main(void)
{
    bm::bvector<>   bv;    // Bitvector variable declaration.

    cout << bv.count() << endl;

    // Set some bits.

    bv.set(10);
    bv.set(100);
    bv.set(1000000);

    // New bitvector's count.

    cout << bv.count() << endl;


    // Print the bitvector.

    unsigned value = bv.get_first();
    do
    {
        cout << value;
        value = bv.get_next(value);
        if (value)
        {
            cout << ",";
        }
        else
        {
            break;
        }
    } while(1);

    cout << endl;

    bv.clear();   // Clean up.

    cout << bv.count() << endl;

    // We also can use operators to set-clear bits;

    bv[10] = true;
    bv[100] = true;
    bv[10000] = true;

    cout << bv.count() << endl;

    if (bv[10])
    {
        bv[10] = false;
    }

    cout << bv.count() << endl;

    return 0;
}

