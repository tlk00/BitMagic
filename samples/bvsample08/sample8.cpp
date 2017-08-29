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

/** \example sample8.cpp

    Example demonstrates some STL compatability set operations using
    set iterators.
 
  \sa bm::bvector<>::enumerator 
  \sa bm::bvector<>::insert_iterator

   For more information please visit:  http://bmagic.sourceforge.net

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
    bm::bvector<>   bv;    

    bv[10] = true;
    bv[100] = true;
    bv[10000] = true;
    
    cout << "Source set:";
    PrintContainer(bv.first(), bv.end());
    
    // copy all bitset information into STL vector using copy algorithm
    {
        vector<unsigned> vect;
        vect.resize(bv.count());
        std::copy(bv.first(), bv.end(), vect.begin());
        cout << "Vector:";
        PrintContainer(vect.begin(), vect.end());
    }

    // doing the same with the help of back_inserter

    {
        list<unsigned> lst;
        std::copy(bv.first(), bv.end(), std::back_inserter(lst));
        cout << "List:";
        PrintContainer(lst.begin(), lst.end());
    }

    {
        vector<unsigned>   vect;
        vector<unsigned>   res1, res2, res3;
        
        vect.push_back(100);
        vect.push_back(15);
        vect.push_back(150);
        
        cout << "Argument vector for set operations:";
        PrintContainer(vect.begin(), vect.end());
        
        // set should be ordered by < to make set algorithms possible
        std::sort(vect.begin(), vect.end());
        cout << endl;
        
        std::set_union(bv.first(), bv.end(),
                       vect.begin(), vect.end(),
                       std::back_inserter(res1)); //10;15;100;150;10000
        cout << "Set union:" << endl;
        PrintContainer(res1.begin(), res1.end());
        
        std::set_intersection(bv.first(), bv.end(),
                              vect.begin(), vect.end(),
                              std::back_inserter(res2));  // 100
        cout << "Set intersection:" << endl;
        PrintContainer(res2.begin(), res2.end());

        vector<unsigned>::const_iterator it1 = vect.begin();
        vector<unsigned>::const_iterator it2 = vect.end();
        bm::bvector<>::enumerator en = bv.first();
        bm::bvector<>::enumerator en2= bv.end();
        
        std::set_difference(en, en2,
                            it1, it2,
                            std::back_inserter(res3));  // 10;10000

        cout << "Set diff:" << endl;
        PrintContainer(res3.begin(), res3.end());
        
    }

    // Using bvector<>::insert_iterator to set bits
    {
        bm::bvector<> bv1;
        std::vector<unsigned> vect;
        
        vect.push_back(300);
        vect.push_back(200);
        vect.push_back(275);
        vect.push_back(200);
        
        cout << endl << "Source vector:";
        PrintContainer(vect.begin(), vect.end()); // 300;200;275;200;
        
        // The "side effect" of this operation is that we sorted
        // the input sequence and eliminated duplicates
        
        std::copy(vect.begin(), vect.end(), bv1.inserter());
        cout << "Bitset:";
        
        PrintContainer(bv1.first(), bv1.end());  // 200;275;300
    }
    
    
    return 0;
}

