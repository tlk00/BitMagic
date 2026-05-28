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
  \sa bm::sparse_vector_float::begin()
  \sa bm::sparse_vector_float::end()
*/

/*! \file svfsample01.cpp
    \brief Example: sparse_vector<> container set values
*/

#include <iostream>
#include <fstream>
#include <cassert>

#include <bm.h>
#include <bmsparsevec_float.h>
#include <bmsparsevec_float_serial.h>

void Demo1(){

    bm::sparse_vector_float svf1;

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

    bm::sparse_vector_float svf1;
    bm::sparse_vector_float svf2;
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

void Demo3(){
    bm::sparse_vector_float svf1;

    float toAdd[] = {1.0123, 2.468, 340000.56};
    svf1.import(toAdd, 3);

    BM_DECLARE_TEMP_BLOCK(tb)
    svf1.optimize(tb);

    // Serializing a sparse_vector_float is a little different to sparse_vector<> and bitvector<>
    //in order to serialize the sparse_vector_float, create a sparse_vector_float_serialized class
    //in this class call .serialize(sparse_vector_float) to serialize
    //or .deserialize(sparse_vector_float) to deserialize
    bm::sparse_vector_float_serialized svf1Serial;
    svf1Serial.serialize(svf1);

    // .size() returns the size in bytes used instead of number of elements
    std::cout << "svf1Serial.size() = " << svf1Serial.size() << std::endl;

    svf1Serial.deserialize(svf1);
    
    std::cout << "svf1.size() = " << svf1.size();

    for(int i = 0; i < svf1.size(); i++){
        std::cout << "svf1.get(" << i << ") = " << svf1.get(i) << std::endl;
        std::cout << "toAdd[" << i << "] = " << toAdd[i] << std::endl;
    }
}

void Demo4(){
    
    bm::sparse_vector_float svf1;
    bm::sparse_vector_float svf2;
    float toAdd1[] = {1.0123, 2.468, 340000.56};
    float toAdd2[] = {7.000, 89000.01, 324.5006};

    //you can import entire arrays into svfs
    svf1.import(toAdd1, 3);
    svf2.import(toAdd2, 3);

    // optimize memory allocation of sparse vector float
    BM_DECLARE_TEMP_BLOCK(tb)
    svf1.optimize(tb);
    svf2.optimize(tb);
    
    // You can create a const_iterator for a svf in multiple ways
    
    // From the start of a svf: bm::sparse_vector_float::const_iterator itFromSV(&svf1);
    // You can also use .begin()
    bm::sparse_vector_float::const_iterator itBegin = svf1.begin();
    // From a specific element of a svf: bm::sparse_vector_float::const_iterator itFromSVPos(&svf1, 1);
    // From the last element of a svf: bm::sparse_vector_float::const_iterator itEnd   = svf1.end();
    // Copying another const_iterator: bm::sparse_vector_float::const_iterator itCopy(itBegin);
    
    // operator* and value() return the current element
    std::cout << "*itBegin = " << *itBegin << std::endl;
    std::cout << "itBegin.value() = " << itBegin.value() << std::endl;
    std::cout << "toAdd[0] = " << toAdd1[0] << std::endl;

    // pos() returns the current position of the iterator
    std::cout << "itBegin.pos() = " << itBegin.pos() << std::endl;

    // valid() returns true if the const_iterator is in a valid position, not bm::id_max
    // invalidate() invalidates the const_iterator
    std::cout << "itBegin.valid() = " << itBegin.valid() << std::endl;
    

    // operator== and operator!= returns true/false if 2 const_iterators are at the same position 
    // of the same svf
    bm::sparse_vector_float::const_iterator itA = svf1.begin();
    bm::sparse_vector_float::const_iterator itB = svf1.begin();
    std::cout << "itA == itB = " << (itA == itB) << std::endl;
    std::cout << "itA != itB = " << (itA != itB) << std::endl;

    //You can use .advance() or ++ as a prefix or postfix to increase the position of the const_iterator by 1
    itB.advance();
    ++itB;
    std::cout << "itB.pos() = " << itB.pos() << std::endl;

    // the < <= > >= operators check the position difference between 2 const_iterators
    // They do not check that the 2 const_iterators are based on the same svf
    std::cout << "itA < itB = " << (itA < itB) << std::endl;
    std::cout << "itA <= itB = " << (itA <= itB) << std::endl;
    std::cout << "itA <= itA = " << (itA <= itB) << std::endl;
    std::cout << "itA > itB = " << (itA > itB) << std::endl;
    std::cout << "itA >= itB = " << (itA >= itB) << std::endl;
    std::cout << "itA >= itA = " << (itA >= itA) << std::endl;

    // You can use go_to() in order to advance to any position in the const_iterator, forwards or backwards
    itBegin.go_to(2);
    std::cout << "itBegin.pos() = " << itBegin.pos() << std::endl;

    itBegin.go_to(0);
    std::cout << "itBegin.pos() = " << itBegin.pos() << std::endl;

}

int main(void){
    try
    {
        // Demo1 for general methods
        Demo1();

        std::cout << std::endl << std::endl;

        //Demo2 for importing an array, optimizing it, and swapping 2 sparse_vector_floats
        Demo2();

        std::cout << std::endl << std::endl;

        //Demo3 for serializing a sparse_vector_float
        Demo3();

        std::cout << std::endl << std::endl;

        //Demo4 for sparse_vector const_iterator methods
        Demo4();
    }
    catch(std::exception& ex)
    {
        std::cerr << ex.what() << std::endl;
        return 1;
    }
    return 0;
}
