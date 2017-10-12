/*
Copyright(c) 2017 Anatoliy Kuznetsov(anatoliy_kuznetsov at yahoo.com)

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

/** \example svsample01.cpp
  Example of how to use bm::sparse_vector<> template class to set values
 
  For more information please visit:  http://bmagic.sourceforge.net

  \sa bm::sparse_vector<>::import
  \sa bm::sparse_vector<>::at
  \sa bm::sparse_vector<>::optimize
  \sa bm::sparse_vector<>::size
 
 */

#include <iostream>
#include "bmsparsevec.h"

using namespace std;

int main(void)
{
    try
    {
        bm::sparse_vector<unsigned, bm::bvector<> > sv1;
        
        unsigned arr[3] = {1,2,3};
        sv1.import(arr, 3); // import from a C-style array (fastest way to populate)

        // optimize memory allocation of sparse vector
        {
            BM_DECLARE_TEMP_BLOCK(tb)
            sv1.optimize(tb);
        }
        
        cout << "sv1.size() = " << sv1.size() << endl;
        cout << "sv[]:";
        
        // print the vector elements using direc access operator
        for (unsigned i = 0; i < sv1.size(); ++i)
        {
            cout << sv1.at(i) << ",";
        }
        cout << endl;
    }
    catch(std::exception& ex)
    {
        std::cerr << ex.what() << std::endl;
        return 1;
    }
    

    return 0;
}

