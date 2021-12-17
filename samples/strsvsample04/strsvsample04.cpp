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

/** \example strsvsample04.cpp
  Example of how to use bm::str_sparse_vector<> - succinct container for
  bit-transposed string collections with NULL (unassigned) values support
 
  \sa bm::str_sparse_vector
  \sa bm::str_sparse_vector::set_null
  \sa bm::str_sparse_vector::push_back
  \sa bm::str_sparse_vector::const_iterator
  \sa bm::str_sparse_vector::optimize

*/

/*! \file strsvsample04.cpp
    \brief Example: str_sparse_vector<> how to work with NULL values
*/

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

#include "bm.h"
#include "bmstrsparsevec.h"
#include "bmundef.h" /* clear the pre-proc defines from BM */

using namespace std;

typedef bm::bvector<> bvector_type;
typedef bm::str_sparse_vector<char, bvector_type, 32> str_sv_type;

int main(void)
{
    try
    {
        // construct container with NULL (unassigned support)
        str_sv_type str_sv(bm::use_null);
        
        const char* s0 = "asz1234";
        std::string str1 = "aqw1234";
        std::string str3 = "54z";
        std::string str00 = "00";
        
        str_sv.set(0, s0);       // set a C-string
        str_sv.push_back(str1);  // add from an STL string
        
        str_sv.assign(3, str3);  // assign element 3 from an STL string
        
        // back-insert iterator can be used to add NULL values
        {
            auto bi = str_sv.get_back_inserter();
            bi = "456";    // add regular value
            bi.add_null(); // add NULL value explicitly
            bi = (const char*)0; // add NULL via 0 ptr
            
            bi.flush(); // flush insert iterator
        }
        
        str_sv.set_null(str_sv.size()); // set a NULL element with automatic resize
        
        // please note that container automatically resizes
        std::cout << "sv size()=" << str_sv.size() << endl; // 4
        
        
        // print out the container content using [] reference
        //
        for (str_sv_type::size_type i = 0; i < str_sv.size(); ++i)
        {
            if (str_sv[i].is_null())
                cout << i << ":NULL" << endl;
            else
            {
                const char* s = str_sv[i];
                cout << i << ":" << s << endl;
            }
        } // for i
        cout << endl;
        
        // clear vector elements in [4, 7] closed interval
        // and set elements to NULL
        //
        str_sv.clear_range(4, 7, true);
        
        
        // print content using const_iterator which also supports is_null()
        //
        {
            str_sv_type::const_iterator it = str_sv.begin();
            str_sv_type::const_iterator it_end = str_sv.end();
            for (; it != it_end; ++it)
            {
                if (it.is_null())
                    cout << "NULL" << endl;
                else
                    cout << *it << endl;
            } // for it
        }

        // bulk clear and set to NULL using bit-vector as an index
        // (faster than random access)
        {
            bvector_type bv_idx { 0, 3, 4, 5 }; // index vector of elements
            str_sv.set_null(bv_idx);
            str_sv.optimize(); // recompress blocks to free some memory
        }
        cout << endl;
        {
            str_sv_type::const_iterator it = str_sv.begin();
            str_sv_type::const_iterator it_end = str_sv.end();
            for (; it != it_end; ++it)
            {
                std::string_view s = it.get_string_view();
                if (s.data() == 0)
                    cout << "NULL" << endl;
                else
                    cout << s << endl;
            } // for it
        }

    }
    catch(std::exception& ex)
    {
        std::cerr << ex.what() << std::endl;
        return 1;
    }
    

    return 0;
}

