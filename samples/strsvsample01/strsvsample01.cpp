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

/** \example strsvsample01.cpp
  Example of how to use bm::str_sparse_vector<> - succinct container for
  bit-transposed string collections
 
  \sa bm::str_sparse_vector
*/

/*! \file strsvsample01.cpp
    \brief Example: str_sparse_vector<> set values, optimize memory
*/

#include <iostream>
#include <string>
#include <vector>
#include <cassert>

#include "bm.h"
#include "bmstrsparsevec.h"
#include "bmundef.h" /* clear the pre-proc defines from BM */

using namespace std;

typedef bm::bvector<> bvector_type;

// define the sparse vector type for 'char' type using bvector as
// a container of bits for bit-transposed planes
// 32 - is maximum string length for this container.
//      Memory allocation is dynamic using sparse techniques, so this number
//      just defines the max capacity.
//
typedef bm::str_sparse_vector<char, bvector_type, 32> str_sv_type;

int main(void)
{
    try
    {
        str_sv_type str_sv;
        
        const char* s0 = "asz1234";
        std::string str1 = "aqw1234";
        std::string str3 = "54z";
        std::string str00 = "00";
        
        str_sv.set(0, s0);       // set a C-string
        str_sv.push_back(str1);  // add from an STL string
        
        str_sv.assign(3, str3);  // assign element 3 from an STL string
        
        // please note that container automatically resizes
        std::cout << "sv size()=" << str_sv.size() << endl; // 4
        
        str_sv.insert(0, str00); // insert element

        str_sv.erase(0);         // erase element 0


        // optimize runs compression in all bit-plains
        {
            BM_DECLARE_TEMP_BLOCK(tb)
            str_sv.optimize(tb);
        }
        
        // print out the container content using [] reference
        //
        for (str_sv_type::size_type i = 0; i < str_sv.size(); ++i)
        {
            const char* s = str_sv[i];
            cout << i << ":" << s << endl;
        } // for i
        cout << "----" << endl;
        
        // print content using const_iterator
        //
        // const iterator has a bulk implementation,
        // not like random access const_iterator is transposing a whole set
        // of elements at once, it is faster for bulk traverse
        //
        {
            str_sv_type::const_iterator it = str_sv.begin();
            str_sv_type::const_iterator it_end = str_sv.end();
            for (; it != it_end; ++it)
            {
                cout << *it << endl;
            } // for it
        }

        // ---------------------------------------------

        str_sv_type str_sv2(bm::use_null); // NULL-able string vector

        str_sv2.set(1, "s1");
        str_sv2.push_back("s2"); 
        str_sv2.push_back_null();
        str_sv2.push_back("s3");

        std::string s;
        bool found = str_sv2.try_get(0, s); // FALSE (NULL value)
        assert(!found);  (void)found;
        found = str_sv2.try_get(2, s); // TRUE, "s2"
        assert(found);
        std::cout << s << std::endl;

        // print out vector content with NULL check
        {
            str_sv_type::const_iterator it = str_sv2.begin();
            it.go_to(1); // position to idx=1, skipping the element 0
            for (; it.valid(); ++it)
            {
                if (it.is_null())
                    cout << "NULL, ";
                else
                    cout << *it << ", ";
            } // for it.valid()
            cout << endl;
        }
    }
    catch(std::exception& ex)
    {
        std::cerr << ex.what() << std::endl;
        return 1;
    }
    

    return 0;
}

