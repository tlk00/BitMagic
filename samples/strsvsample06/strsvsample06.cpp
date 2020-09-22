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

/** \example strsvsample06.cpp
    Succinct container iterator for (substring), search, mismatch

  \sa bm::str_sparse_vector
  \sa bm::sparse_vector_find_mismatch
  \sa bm::sparse_vector_find_first_mismatch
  \sa bm::sparse_vector_scanner
*/

/*! \file strsvsample06.cpp
    \brief Example: Succinct container iterator for (substring), search, mismatch
*/

#include <iostream>
#include <string>
#include <vector>
#include <cassert>

#include "bm.h"
#include "bmstrsparsevec.h"
#include "bmsparsevec_algo.h"
#include "bmundef.h" /* clear the pre-proc defines from BM */

using namespace std;

typedef bm::bvector<> bvector_type;
typedef bm::str_sparse_vector<char, bvector_type, 32> str_sv_type;

int main(void)
{
    try
    {
        str_sv_type str_sv(bm::use_null); // construct it as a NULL-able vector

        // back inserted is the fastest way to load the container, use it
        // where necessary
        {
            auto iit = str_sv.get_back_inserter();
            iit = "240-423-0567";
            iit = "420-123-5078";
            iit.add_null();
            iit = "301-223-1234";
            iit = "131-423-1235";
            iit.add_null();
            iit = "301-113-6535";

            iit.flush();
        }
        BM_DECLARE_TEMP_BLOCK(tb)
        str_sv.optimize(tb); // optimize the vector

        // print the content using iterator
        {
            auto it = str_sv.begin();
            auto it_end = str_sv.end();
            for (; it != it_end; ++it)
            {
                if (it.is_null())
                    cout << "NULL";
                else
                    cout << *it;
                cout << endl;
            }
        }

        // print all the NULL values in the vector
        //
        {
            const bvector_type* bv_not_null = str_sv.get_null_bvector();
            assert(bv_not_null);
            bvector_type bv_tmp(*bv_not_null);
            bv_tmp.resize(str_sv.size());
            bv_tmp.invert();
            if (bv_tmp.any())
            {
                auto en = bv_tmp.get_enumerator(0);
                cout << "List of NULL values: ";
                for (; en.valid(); ++en)
                {
                    cout << *en << ", ";
                }
                cout << endl;
            }
        }



        // use substring iterators
        // substring performs partial transposition and can be faster
        // than the full iterator with extract ans substr
        //
        {
            cout << endl << "--------------------------\n";
            auto it1 = str_sv.begin();
            it1.set_substr(0, 3);

            auto it2 = str_sv.begin();
            it2.set_substr(4, 3);

            auto it3 = str_sv.begin();
            it3.set_substr(8); // from 8 to end

            auto it_end = str_sv.end();
            for (; it1 != it_end; ++it1, ++it2, ++it3)
            {
                if (it1.is_null())
                    cout << "NULL";
                else
                    cout << "(" << *it1 << ")" << *it2 << "." << *it3;
                cout << endl;
            } // for
        }

        str_sv_type str_sv1(str_sv); // construct a copy



        // see also sparse_vector_find_mismatch
        str_sv_type::size_type idx;
        bool b = bm::sparse_vector_find_first_mismatch(str_sv, str_sv1, idx);
        assert(!b); // should not find at this point

        str_sv1.set(2, "788-125-6789");
        b = bm::sparse_vector_find_first_mismatch(str_sv, str_sv1, idx);
        if (b)
        {
            cout << endl;
            cout << "Mismatch at: " << idx << endl;
            cout << "value=" << str_sv1[idx] << endl;
        }
        cout << endl;

        // use scanner to perform vector search
        //
        // NOTE: re-use scanner for multiple searches for best performance
        //
        bm::sparse_vector_scanner<str_sv_type> str_scan;

        {
            bvector_type bv_res;
            str_scan.find_eq_str(str_sv1, "250-113-6535", bv_res);
            {
                auto en = bv_res.get_enumerator(0); // from the start (pos=0)
                for (;en.valid(); ++en)
                {
                    idx = *en;
                    cout << idx << ": " << str_sv1[idx] << endl;
                } // for
            }

            cout << "Prefix search:" << endl;
            str_scan.find_eq_str_prefix(str_sv1, "301", bv_res);
            cout << "Found: " << bv_res.count() << endl;

            auto en = bv_res.get_enumerator(0); // from the start (pos=0)
            for (;en.valid(); ++en)
            {
                idx = *en;
                cout << idx << ": " << str_sv1[idx] << endl;
            } // for


        }



    }
    catch(std::exception& ex)
    {
        std::cerr << ex.what() << std::endl;
        return 1;
    }
    

    return 0;
}

