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

/** \example strsvsample03.cpp
  Example of how to use bm::str_sparse_vector<> - succinct container for
  bit-transposed string collections
 
  \sa bm::str_sparse_vector
*/

/*! \file strsvsample03.cpp
    \brief Example: str_sparse_vector<> back insert iterator example
 
    This example loads sparse vector from an STL container uses character re-mapping
    to compress, serialize and save container to disk.
    Example also illustrates how to check memory footprint.
*/

#include <assert.h>
#include <iostream>
#include <string>
#include <vector>
#include <random>
#include <algorithm>
#include <fstream>


#include "bm.h"
#include "bmstrsparsevec.h"
#include "bmsparsevec_serial.h"
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


// generate collection of strings from integers and shuffle it
//
static
void generate_string_set(vector<string>& str_vec)
{
    const unsigned max_coll = 50000;
   
    str_vec.resize(0);
    string str;
    for (unsigned i = 10; i < max_coll; i += (unsigned)rand() % 3)
    {
        str = to_string(i);
        str_vec.emplace_back(str);
    } // for i
    
    // shuffle the data set
    //
    std::random_device rd;
    std::mt19937       g(rd());
    std::shuffle(str_vec.begin(), str_vec.end(), g);
}


int main(void)
{
    try
    {
        str_sv_type str_sv;

        vector<string> str_vec;
        generate_string_set(str_vec);
        std::sort(str_vec.begin(), str_vec.end()); // sort the input vector


        // load sparse vector from an STL container
        //
        {
            size_t vect_size = 0;   // approx std::vector<string> memory usage
            str_sv_type str_sv_tmp; // temp vector
            {
                str_sv_type::back_insert_iterator bi =
                                        str_sv_tmp.get_back_inserter();
                for (auto str : str_vec)
                {
                    bi = str;
                    
                    // some approximate estimate of std::string element cost
                    //
                    size_t str_size = str.size() + sizeof(str);
                    vect_size += str_size;
                }
                
                // it is important to use flush, because back inserter is
                // buffering data. Of cause it flashes automatically on
                // destruction but explicit flush is somewhat better
                // because of possible exception is thrown here and not from
                // destructor.
                //
                
                bi.flush();
                
                cout << "STL vector<string> approx.memory consumption:"
                     << vect_size << endl;
            }
            
            // calculate memory footprint
            //
            str_sv_type::statistics st;
            str_sv_tmp.calc_stat(&st);
            
            cout << "Used memory: " << st.memory_used << std::endl;
            
            
            // final step is re-mapping, which increses chances for
            // good memory compression.
            // A side-effect here is that remapping makes container
            // effectively read-only.
            //
            str_sv.remap_from(str_sv_tmp);

            BM_DECLARE_TEMP_BLOCK(tb)
            str_sv.optimize(tb); // optimize the vector
            
            str_sv.calc_stat(&st);
            cout << "Used memory after remap and optimization: "
                 << st.memory_used
                 << std::endl;


            // a slightly different way to do the reampped loading
            // iterator in this case is used to collect character frequency
            // tables and do remap on flush() (which is a bit faster)
            //
            // another option here is to turn vector into read-only mode
            // which runs an embedded memory defragmentation algorithm
            //
            {
                str_sv_type str_sv1;
                {
                str_sv_type::back_insert_iterator bi =
                                        str_sv1.get_back_inserter();
                bi.set_remap(true);
                for (auto str : str_vec)
                    bi = str;
                bi.flush();
                }


                str_sv1.optimize(tb);

                // freeze is used to turn the vector into read-only (immutable)
                //
                str_sv1.freeze();

                assert(str_sv1.is_ro());
                bool eq = str_sv1.equal(str_sv);
                assert(eq); (void)eq;

                str_sv1.calc_stat(&st);
                cout << "Used memory after remap / optimization / freeze: "
                     << st.memory_used
                     << std::endl;

                // construct a vector with remap table derived from another
                // we can just add the same data into the new vector
                // or part of the same data (projection) or reorderd data.
                //
                // the key thing is that we have to guarantee that the
                // remapping tables will be the same
                //
                {
                str_sv_type str_sv2(str_sv1, bm::remap_setup::COPY_RTABLES);
                {
                    unsigned cnt = 0;
                    str_sv_type::back_insert_iterator bi =
                                            str_sv2.get_back_inserter();
                    for (auto str : str_vec)
                    {
                        bi = str;
                        if (++cnt >= 10)
                            break; // pick first 10
                    }
                    bi.flush();
                }
                str_sv2.optimize();
                cout << "size2=" << str_sv2.size() << endl; // 10
                }

            }


        }
        
        // serialize and save
        //
        {
            std::string fname = "test.sv";
            bm::sparse_vector_serial_layout<str_sv_type> sv_lay;
            
            BM_DECLARE_TEMP_BLOCK(tb)
            bm::sparse_vector_serialize(str_sv, sv_lay, tb);

            std::ofstream fout(fname.c_str(), std::ios::binary);
            if (!fout.good())
            {
                return -1;
            }
            const char* buf = (char*)sv_lay.buf();
            fout.write(buf, (unsigned)sv_lay.size());
            if (!fout.good())
            {
                return -1;
            }
            fout.close();
            
            cout << "Saved size: " << sv_lay.size() << endl;
        }

    }
    catch(std::exception& ex)
    {
        std::cerr << ex.what() << std::endl;
        return 1;
    }
    

    return 0;
}

