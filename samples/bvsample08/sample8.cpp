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

/** \example sample8.cpp

  Example for STL interoperability and set operations with iterators.
 
  \sa bm::bvector<>::enumerator 
  \sa bm::bvector<>::insert_iterator
*/

/*! \file sample8.cpp
    \brief Example: bvector<> - STL interoperability
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
#include "bmundef.h" /* clear the pre-proc defines from BM */

using namespace std;

inline
void Print(bm::bvector<>::size_type n)
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
    typedef bm::bvector<>::size_type bm_size_type;
    try
    {
        bm::bvector<>   bv;

        bv[10] = true;
        bv[100] = true;
        bv[10000] = true;
        
        cout << "Source set:";
        PrintContainer(bv.first(), bv.end());
        
        // copy all bitset information into STL vector using copy algorithm
        {
            vector<bm_size_type> vect;
            vect.resize(bv.count());
            std::copy(bv.first(), bv.end(), vect.begin());
            cout << "Vector:";
            PrintContainer(vect.begin(), vect.end());
        }

        // doing the same with the help of back_inserter

        {
            list<bm_size_type> lst;
            std::copy(bv.first(), bv.end(), std::back_inserter(lst));
            cout << "List:";
            PrintContainer(lst.begin(), lst.end());
        }

        {
            vector<bm_size_type>   vect;
            vector<bm_size_type>   res1, res2, res3;
            
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

            vector<bm_size_type>::const_iterator it1 = vect.begin();
            vector<bm_size_type>::const_iterator it2 = vect.end();
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
            std::vector<bm_size_type> vect;
            
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
    }
    catch(std::exception& ex)
    {
        std::cerr << ex.what() << std::endl;
        return 1;
    }
    
    
    return 0;
}

