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

/** \example sample2.cpp
  Example demonstrates using set operations AND, OR, XOR, etc.
  For more information please visit:  http://bmagic.sourceforge.net
 
*/


#include <iostream>
#include "bm.h"

using namespace std;

void print_bvector(const bm::bvector<>& bv)
{
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
}

int main(void)
{
    bm::bvector<>   bv1;    
    bm::bvector<>   bv2;
    bm::bvector<>   bv3;

    bv1.set(10);
    bv1.set(100);
    bv1.set(1000000);


    bv2.set(10);
    bv2.set(100);

    // Logical AND operation on bv2 (bv1 is the argument)
    // bv2 = bv2 AND bv1

    bv3 = bv1 & bv2;
    print_bvector(bv3);

    bv2 &= bv1;  // You also can use: bv2.bit_and(bv1);
    print_bvector(bv2);
    
    // bv2 = bv2 OR bv1

    bv3 = bv1 | bv2;
    print_bvector(bv3);

    bv2 |= bv1;  //  You can also use: bv2.bit_or(bv1);
    print_bvector(bv2);

    
    bv1.set(1000000, false);
    
    // bv2 = bv2 SUB bv1

    bv3 = bv2 - bv1;
    print_bvector(bv3);

    bv2 -= bv1;   // You can also use: bv2.bit_sub(bv1);
    print_bvector(bv2);

    // bv2 XOR bv1

    bv3 = bv2 ^ bv1;
    print_bvector(bv3);

    bv2 ^= bv1;  // You can also use: bv2.bit_xor(bv1);
    print_bvector(bv2);

    // For lexicographical comparison there is set of overloaded
    // operators and function compare.

    if (bv2 == bv3)
    {
        cerr << "Equivalent. Comparison result = " 
             << bv2.compare(bv3) << endl;
    }
    else
    {
        cout << "Error." << endl;
        return 1;
    }


    return 0;
}
