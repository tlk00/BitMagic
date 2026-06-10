/*
Copyright(c) 2026 Anatoliy Kuznetsov(tolikkuznetsov66 at gmail.com)

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

/** \example svfsample01.cpp
  Example of how to use bm::sparse_vector_float template class to set or
  import values.

  bm::sparse_vector_float is a wrapper of bm::sparse_vector<> that takes in floats and converts them into
  a bitvector of signs, and 2 sparse_vectors of the floats exponent and mantissa
 
  \sa bm::sparse_vector_float::push_back
  \sa bm::sparse_vector_float::set
  \sa bm::sparse_vector_float::size
  \sa bm::sparse_vector_float::empty
  \sa bm::sparse_vector_float::get
  \sa bm::sparse_vector_float::import
  \sa bm::sparse_vector_float::swap
  \sa bm::sparse_vector_float::optimize
*/

/*! \file svfsample01.cpp
    \brief Example: sparse_vector_float container set values
*/

#include <iostream>
#include <fstream>
#include <cassert>

#include <bm.h>
#include <bmsparsevec_float.h>
#include <bmsparsevec_float_serial.h>

typedef bm::sparse_vector_float<bm::bvector<>> sparse_vec_float;

void Demo1(){

    sparse_vec_float svf1;

    //initially the sparse_vector_float is empty
    std::cout << "svf1.empty() = " << svf1.empty() << std::endl;

    float toAdd[] = {1.0123, 2.468, 340000.56};
    
    //you can add to the end of the svf with push_back
    svf1.push_back(toAdd[0]);

    //get values with .get()
    std::cout << "svf1.get(0) = " << svf1.get(0) << std::endl;
    std::cout << "toAdd[0] = " << toAdd[0] << std::endl;

    //You can set any element in the svf to a float value you want
    svf1.set(1, toAdd[1]);
    svf1.set(100, toAdd[2]);

    std::cout << "svf1.get(1) = " << svf1.get(1) << std::endl;
    std::cout << "toAdd[1] = " << toAdd[1] << std::endl;
    std::cout << "svf1.get(100) = " << svf1.get(100) << std::endl;
    std::cout << "toAdd[2] = " << toAdd[2] << std::endl;

    //if a value inserted is greater than the size then the values between the last value and the set value
    //are set to 0
    std::cout << "svf1.get(50) = " << svf1.get(50) << std::endl;

    //check the size with .size()
    std::cout << "svf1.size() = " << svf1.size();
    
}

void Demo2(){

    sparse_vec_float svf1;
    sparse_vec_float svf2;
    float toAdd1[] = {1.0123, 2.468, 340000.56};
    float toAdd2[] = {7.000, 89000.01, 324.5006};

    //you can import entire arrays into svfs
    svf1.import(toAdd1, 3);
    svf2.import(toAdd2, 3);

    // optimize memory allocation of sparse vector float
    BM_DECLARE_TEMP_BLOCK(tb)
    svf1.optimize(tb);
    svf2.optimize(tb);

    std::cout << "svf1.size() = " << svf1.size() << std::endl;
    std::cout << "svf2.size() = " << svf2.size() << std::endl;

    for(int i = 0; i < svf1.size(); i++){
        std::cout << "svf1.get(" << i << ") = " << svf1.get(i) << std::endl;
        std::cout << "toAdd1[" << i << "] = " << toAdd1[i] << std::endl;

        std::cout << "svf2.get(" << i << ") = " << svf2.get(i) << std::endl;
        std::cout << "toAdd2[" << i << "] = " << toAdd2[i] << std::endl;
    }

    //swap two svfs with swap
    svf1.swap(svf2);
    std::cout << "Swapped" << std::endl;

    for(int i = 0; i < svf1.size(); i++){
        std::cout << "svf1.get(" << i << ") = " << svf1.get(i) << std::endl;
        std::cout << "toAdd1[" << i << "] = " << toAdd1[i] << std::endl;

        std::cout << "svf2.get(" << i << ") = " << svf2.get(i) << std::endl;
        std::cout << "toAdd2[" << i << "] = " << toAdd2[i] << std::endl;
    }

}

int main(void){
    try
    {
        // Demo1 for general methods
        Demo1();

        std::cout << std::endl << std::endl;

        //Demo2 for importing an array, optimizing it, and swapping 2 sparse_vector_floats
        Demo2();
    }
    catch(std::exception& ex)
    {
        std::cerr << ex.what() << std::endl;
        return 1;
    }
    return 0;
}
