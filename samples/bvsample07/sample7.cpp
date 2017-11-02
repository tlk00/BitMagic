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

/** \example sample7.cpp

  Example how to use logical operations between arrays and bit-vectors
 
  \sa bm::combine_and()
  \sa bm::combine_and_sorted()
  \sa bm::combine_sub()
  \sa bm::combine_or()
  \sa bm::combine_xor()

*/

#include <iostream>
#include <algorithm>
#include <vector>
#include <list>

using std::vector;
using std::list;

// This example requires STL compatibility
#ifdef BM_NO_STL
# undef BM_NO_STL
#endif

#include "bm.h"
#include "bmalgo.h"

using namespace std;

void Print(unsigned n)
{
    cout << n << endl;;
}

// Utility template function used to print container
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
    try
    {
        bm::bvector<>   bv;
        bv[10] = true;
        bv[100] = true;
        bv[10000] = true;

        // initialize unsorted, fairly random array for an experiment
        // it even allowes duplicates (see 12)
        //
        unsigned int arr[] = {2, 10000, 5, 12, 255, 12, 300};

        
        cout << "Source set 1:";
        PrintContainer(bv.first(), bv.end());
        
        cout << "Source set 2:";
        PrintContainer(&arr[0], &arr[0] + (sizeof(arr)/sizeof(arr[0])));
        
        // AND operation between bit-vector and a plain array
        // expect one result: 10000
        // please note, that array in this case comes unsorted
        //
        bm::combine_and(bv, &arr[0], &arr[0] + (sizeof(arr)/sizeof(arr[0])));
        cout << "Result 1(AND): ";
        PrintContainer(bv.first(), bv.end());

        
        // re-initalize the bit-vector
        bv.clear();
        bv[10] = true;
        bv[100] = true;
        bv[10000] = true;

        // OR operation to merge bit-vector and array
        // please note that it naturally works as sort-unique for the array
        //
        bm::combine_or(bv, &arr[0], &arr[0] + (sizeof(arr)/sizeof(arr[0])));
        cout << "Result 2(OR): ";
        PrintContainer(bv.first(), bv.end());

        // sort the array, using STL sort method
        // combine operation on sorted arrays tend to be faster
        //
        std::sort(&arr[0], &arr[0] + (sizeof(arr)/sizeof(arr[0])));
        
        // AND on sorted array is faster
        //
        bm::combine_and_sorted(bv, &arr[0], &arr[0] + (sizeof(arr)/sizeof(arr[0])));
        cout << "Result 3(AND): ";
        PrintContainer(bv.first(), bv.end());

        // SUB (AND NOT or MINUS) also works faster on sorted input
        // the result should be an EMPTY set
        bm::combine_sub(bv, &arr[0], &arr[0] + (sizeof(arr)/sizeof(arr[0])));
        cout << "Result 4(MINUS): ";
        PrintContainer(bv.first(), bv.end());

    }
    catch(std::exception& ex)
    {
        std::cerr << ex.what() << std::endl;
        return 1;
    }
    
    
    return 0;
}

