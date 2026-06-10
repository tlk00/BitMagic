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

/** \example svfsample02.cpp
  Example of how to use bm::sparse_vector_float template class to serialize bm::sparse_vector_float
  with bm_sparse_vector_float_serialized and how to use bm::sparse_vector_float's const_iterator

  bm::sparse_vector_float is a wrapper of bm::sparse_vector<> that takes in floats and converts them into
  a bitvector of signs, and 2 sparse_vectors of the floats exponent and mantissa
 
  \sa bm::sparse_vector_float::clear
  \sa bm::sparse_vector_float::clear_range
  \sa bm::sparse_vector_float::equal
  \sa bm::sparse_vector_float::compare

  \sa bm::sparse_vector_float::join
  \sa bm::sparse_vector_float::merge
  \sa bm::sparse_vector_float::copy_range
*/

/*! \file svfsample03.cpp
    \brief Example: sparse_vector_float container comparison and interaction between containers
*/

#include <iostream>
#include <fstream>
#include <cassert>

#include <bm.h>
#include <bmsparsevec_float.h>
#include <bmsparsevec_float_serial.h>

void Demo1(){
    float toAdd[] = {1.0123, -2.468, 340000.56, -7008.0, 0.900102};

    bm::sparse_vector_float<bm::bvector<>> svf1;
    svf1.import(toAdd, 5);
    BM_DECLARE_TEMP_BLOCK(tb)
    svf1.optimize(tb);

    //You can remove all elements in a sparse_vector_float using the clear() method
    std::cout << "svf1.size() = " << svf1.size() << std::endl;
    svf1.clear();
    std::cout << "svf1.size() = " << svf1.size() << std::endl;

    svf1.import(toAdd, 5);

    //clear_range(left, right) sets all values between left and right inclusive to 0
    svf1.clear_range(1, 3);
    std::cout << "svf1.size() = " << svf1.size();
    for(int i = 0; i < 5; i++){
        std::cout << "svf1.get(" << i << ") = " << svf1.get(i) << std::endl;
        std::cout << "toAdd[" << i << "] = " << toAdd[i] << std::endl;
    }

    svf1.clear();
    svf1.import(toAdd, 5);

    //Copy constructor
    bm::sparse_vector_float<bm::bvector<>> svf2(svf1);

    //equal() checks if two svf's have the same values in them
    std::cout << "svf1.equal(svf2) = " << svf1.equal(svf2) << std::endl;
    svf2.set(1, 0.0);
    std::cout << "svf1.equal(svf2) = " << svf1.equal(svf2) << std::endl;

    //compare(index, compareTo, margin) allows you to check how a certain index compares to some float
    //A -1 means the given float is larger than the index float
    //A 0 means the two are equal, or within the margin of error of equality given (by default epsilon)
    //A 1 means the given float is smaller than the index float
    std::cout << "svf1.compare(1, 2.468) = " << svf1.compare(1, 2.468) << std::endl;
    std::cout << "svf1.compare(1, -2.468) = " << svf1.compare(1, -2.468) << std::endl;
    std::cout << "svf1.compare(1, -3) = " << svf1.compare(1, -3) << std::endl;
}

void Demo2(){
    bm::sparse_vector_float<bm::bvector<>> svf1;
    bm::sparse_vector_float<bm::bvector<>> svf2;
    float toAdd1[] = {1.0123, -2.468, 0.0, 0.0, 0.0, 1.5};
    float toAdd2[] = {0.0, 0.0, 0.0, -7008.0, 0.900102, 2.5};
    svf1.import(toAdd1, 6);
    svf2.import(toAdd2, 6);
    BM_DECLARE_TEMP_BLOCK(tb)
    svf1.optimize(tb);
    svf2.optimize(tb);

    //You can combine 2 sparse_vector_floats using the join method
    //Be careful when joining, as when joining 2 floats together, if one of them is not 0 then it is possible
    //to make an NaN float if the exponent becomes all 1's
    svf1.join(svf2);
    for(int i = 0; i < svf1.size(); i++){
        std::cout << "svf1.get(" << i << ") = " << svf1.get(i) << std::endl;
    }

    svf1.clear();
    svf2.clear();
    svf1.import(toAdd1, 6);
    svf2.import(toAdd2, 6);
    svf1.optimize(tb);
    svf2.optimize(tb);

    //You can also use the merge operation which is different from join as 
    //it borrows data from the source vector, so it gets modified.
    svf1.merge(svf2);
    for(int i = 0; i < svf1.size(); i++){
        std::cout << "svf1.get(" << i << ") = " << svf1.get(i) << std::endl;
    }


    svf1.clear();
    svf2.clear();
    svf1.import(toAdd1, 6);
    svf2.import(toAdd2, 6);
    svf1.optimize(tb);
    svf2.optimize(tb);

    //copy_range(left, right) copies the indeces in the given range, inclusive
    //all indeces outside that range are set to 0
    svf1.copy_range(svf2, 2, 4);
    for(int i = 0; i < svf1.size(); i++){
        std::cout << "svf1.get(" << i << ") = " << svf1.get(i) << std::endl;
    }
}

int main(void){
    try
    {
        //Demo1 for clearing and comparison operations
        Demo1();

        std::cout << std::endl << std::endl;

        //Demo2 for combining sparse_vector_float operations
        Demo2();
    }
    catch(std::exception& ex)
    {
        std::cerr << ex.what() << std::endl;
        return 1;
    }
    return 0;
}
