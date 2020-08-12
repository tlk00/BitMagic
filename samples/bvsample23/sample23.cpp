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

/** \example sample23.cpp

  Example on intervals enumarator.

  BitMagic bvector<> implements high performance operation on RLE
  coded bit-vectors, transparently supporting all logical operations
  like intersections.Serialization uses compressive encoding
  (binary interpolative codes) to efficiently store collections 
  of intervals.

  This example illustrates use of inetrval_enumerator<> to interpret
  bit-vector as a sequence of 011110 ranges/intervals.
 
  \sa bm::bvector::set_range
  \sa bm::interval_enumerator

  \sa sample22.cpp
  \sa bvintervals
*/

/*! \file sample23.cpp
    \brief Example: interval_enumerator<> - interator class for intervals
*/

#include <iostream>
#include <assert.h>

#include "bm.h"
#include "bmintervals.h"
#include "bmundef.h" /* clear the pre-proc defines from BM */

using namespace std;

typedef bm::interval_enumerator<bm::bvector<> > interval_enumerator_type;

int main(void)
{
    try
    {        
        bm::bvector<> bv(bm::BM_GAP); // use RLE compressed vector from the start
        
        bv.set_range(10, 10);
        bv.set_range(100, 110); // sets a range of 1s [100, 110] .....11111....
        bv.set_range(777, 888);
        bv.set_range(65536, 65536);

        bv.optimize(); 


        {
            interval_enumerator_type ien(bv);
            if (ien.valid())
            {
                do
                { 
                    cout << "[" << ien.start() << ".." << ien.end() << "]";
                } while (ien.advance());
                cout << endl;
            }
        }

        // case 2: slightly less efficient, but more STL iterator conformant 
        // 
        {
            interval_enumerator_type ien(bv);
            interval_enumerator_type ien_end;

            // please note prefix increment "++en" is way more efficient 
            // than possible postfix notation of "ien++" (no temp.copy)
            //
            for (; ien != ien_end; ++ien)
            {
                // also uses pair notation to represent interval
                cout << "[" << (*ien).first << ".." << (*ien).second << "]";
            }
            cout << endl;
        }


        {
            // construct enumerator to start from position 102 and 
            // extend start (to 100)
            interval_enumerator_type ien(bv, 102, true);
            interval_enumerator_type ien_end;
            for (; ien != ien_end; ++ien)
            {
                cout << "[" << ien.get().first << ".." << ien.get().second << "]";
            }
            cout << endl;

            ien.go_to(105, false); // now re-position enumerator to position 105 without extend start
            for (; ien != ien_end; ++ien)
            {
                cout << "[" << ien.get().first << ".." << ien.get().second << "]";
            }
            cout << endl;

            // now re-position enumerator to position 105 without extend start
            ien.go_to(105, false); 
            for (; ien != ien_end; ++ien)
            {
                cout << "[" << ien.get().first << ".." << ien.get().second << "]";
            }
            cout << endl;

            // now re-position enumerator to position 115 wit extend start
            // target position is not in the interval, so it should find the
            // next available one automatically
            //
            ien.go_to(115, true); // extend_left=true but it should not matter
            for (; ien != ien_end; ++ien)
            {
                cout << "[" << ien.get().first << ".." << ien.get().second << "]";
            }
            cout << endl;

            // go beyond end
            ien.go_to(1150000, true);
            if (ien.valid())
            {
                assert(0);
            }
            else
            {
                cout << "EMPTY" << endl; 
            }
        }

    }
    catch(std::exception& ex)
    {
        std::cerr << ex.what() << std::endl;
        return 1;
    }
        
    return 0;
}

