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
  Example of how to use bm::sparse_vector_float template class to serialize and deserialize a 
  bm::sparse_vector_float and how to use bm::sparse_vector_float's const_iterator

  bm::sparse_vector_float is a wrapper of bm::sparse_vector<> that takes in floats and converts them into
  a bitvector of signs, and 2 sparse_vectors of the floats exponent and mantissa
 
  \sa bm::sparse_vector_float_serialize
  \sa bm::sparse_vector_float_deserialize
  \sa bm::sparse_vector_float_serial_layout::size
  \sa bm::sparse_vector_float_serial_layout::capacity
  \sa bm::sparse_vector_float_serial_layout::buf
  \sa bm::sparse_vector_float_serial_layout::reserve
  \sa bm::sparse_vector_float_serial_layout::resize
  \sa bm::sparse_vector_float_serial_layout::shrink
  \sa bm::sparse_vector_float_serial_layout::freemem
  \sa bm::sparse_vector_float_deserializer::deserialize
  \sa bm::sparse_vector_float_deserializer::deserialize(mask)
  \sa bm::sparse_vector_float_deserializer::deserialize_range
  
  \sa bm::sparse_vector_float::begin
  \sa bm::sparse_vector_float::end
  \sa bm::sparse_vector_float::const_iterator::operator*
  \sa bm::sparse_vector_float::const_iterator::operator++
  \sa bm::sparse_vector_float::const_iterator::&operator++
  \sa bm::sparse_vector_float::const_iterator::operator=
  \sa bm::sparse_vector_float::const_iterator::operator==
  \sa bm::sparse_vector_float::const_iterator::operator!=
  \sa bm::sparse_vector_float::const_iterator::operator<
  \sa bm::sparse_vector_float::const_iterator::operator<=
  \sa bm::sparse_vector_float::const_iterator::operator>
  \sa bm::sparse_vector_float::const_iterator::operator>=
  \sa bm::sparse_vector_float::const_iterator::advance
  \sa bm::sparse_vector_float::const_iterator::go_to
  \sa bm::sparse_vector_float::const_iterator::value
  \sa bm::sparse_vector_float::const_iterator::pos
  \sa bm::sparse_vector_float::const_iterator::valid
  \sa bm::sparse_vector_float::const_iterator::invalidate
*/

/*! \file svfsample02.cpp
    \brief Example: sparse_vector_float container serialization and const_iterator
*/

#include <iostream>
#include <fstream>
#include <cassert>

#include <bm.h>
#include <bmsparsevec_float.h>
#include <bmsparsevec_float_serial.h>

typedef bm::sparse_vector_float<bm::sparse_vector<unsigned int, bm::bvector<>>> sparseVecFloat;

void Demo1(){


    //Creating and optimizing a sparse_vector_float
    float toAdd[] = {1.0123, -2.468, 340000.56};
    sparseVecFloat svf1;
    svf1.import(toAdd, 3);
    BM_DECLARE_TEMP_BLOCK(tb)
    svf1.optimize(tb);

    //To serialize a sparse_vector_float you need to create a sparse_vector_float_serial_layout
    bm::sparse_vector_float_serial_layout<sparseVecFloat> layout1;
    
    //You can then serialize by calling sparse_vector_float's serialize method giving it the svf to serialize and the
    //layout to serialize into
    bm::sparse_vector_float_serialize(svf1, layout1);

    //The size() method returns the size of the buf in layout in bytes
    std::cout << layout1.size() << std::endl;
    //The capacity() method returns how many bytes are reservered for the buf in the layout
    std::cout << layout1.capacity() << std::endl;

    //The reserve(size) method allocates more capacity to the layout if more capacity is given to it than is currently used
    //layout1.reserve(layout1.capacity() + 100);
    //The resize(size) methods reserves() the size given to it and then changes the size to be the size given
    //layout1.resize(layout1.size() + 100);
    //The shrink(size) method decreases the size of the layout
    //layout1.shrink(layout1.size()-100);

    /*
        In order to deserialize you must get the buf from the layout, and call the deserialize function
        Giving it the buffer and the svf. This will clear the given svf before deserializing
        buf consists of 4 parts:
            Header - locations of the other 3 parts
            sign buf - the serialized sign bvector
            exponent buf - the serialized exponents sparse_vector
            mantissa buf - the serialized mantissas sparse_vector
    */
    const unsigned char* buf = layout1.buf();
    //const unsigned char* buf = testLayout.data(); - you can also use .data() to get the buf
    bm::sparse_vector_float_deserialize(svf1, buf);
    
    //Deserializing an svf layout returns the values into the svf
    std::cout << svf1.size() << std::endl;
    for(int i = 0; i < 3; i++){
        std::cout << "svf1.get(" << i << ") = " << svf1.get(i) << std::endl;
        std::cout << "toAdd[" << i << "] = " << toAdd[i] << std::endl;
    }

    //The freemem() method frees the buf and resets the layout
    layout1.freemem();



    sparseVecFloat svf2;
    int N = 10000;
    for(int i = 0; i < N; i++){
        float f = i * 0.000123;
        svf2.push_back(f);
    }

    svf2.optimize(tb);
    bm::sparse_vector_float_serial_layout<sparseVecFloat> layout2;
    bm::sparse_vector_float_serialize(svf2, layout2);


    buf = layout2.buf();
    //Instead of just calling deserialize, you can create a deserializer to customize the deserialization
    bm::sparse_vector_float_deserializer<sparseVecFloat> svfD;
    sparseVecFloat svf2_restored;

    //This deserializer is what gets called in the deserialize(testSVF, buf) function
    //Using a customized deserializer you can instead not clear the svf you deserialize into
    //svfD.deserialize(svf2, buf, true);

    //deserialize_range will be given the svf to deserialize into, the buf, the indeces of the original
    //svf that was serialized to deserialize, and then whether to clear the new svf before deserializing into it
    svfD.deserialize_range(svf2_restored, buf, 300, 400, true);
    
    //it will deserialize into the given indeces of the new SVF, not at the start or end of it
    int errorCount = 0;
    for (int i = 300; i < 400; i++) {
        float f = i * 0.000123;
        if (std::fabs(svf2_restored.get(i) - f) > 0.001f){
            errorCount++;
        }
    }
    std::cout << errorCount << std::endl;
    assert(errorCount == 0);



    bm::bvector<> mask_bv;
    int maskIndices[] = {0, 1, 50, 100, 500, 999, 5000, 9999};
    int maskSize = sizeof(maskIndices) / sizeof(maskIndices[0]);
    for (int i = 0; i < maskSize; i++)
        mask_bv.set(maskIndices[i]);
    
    //You can also deserialize a mask of indeces
   sparseVecFloat svf2_masked;
    svfD.deserialize(svf2_masked, buf, mask_bv);

    //The given indeces in the mask will be deserialized into the given svf, and those not in it will be set to 0
    for (int i = 0; i < maskSize; i++) {
        int idx   = maskIndices[i];
        float f   = idx * 0.000123f;
        std::cout << "svf2_masked.get(" << idx << ") = " << svf2_masked.get(idx) << std::endl;
        std::cout << idx << " * " << 0.000123f << " = " << f << std::endl;
    }
    std::cout << "testSVF.get(2) = " << svf2_masked.get(2) << std::endl;
}

void Demo2(){
    
    sparseVecFloat svf1;
    sparseVecFloat svf2;
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
    sparseVecFloat::const_iterator itBegin = svf1.begin();
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
    sparseVecFloat::const_iterator itA = svf1.begin();
    sparseVecFloat::const_iterator itB = svf1.begin();
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
        //Demo1 for serializing a sparse_vector_float
        Demo1();

        std::cout << std::endl << std::endl;

        //Demo2 for sparse_vector const_iterator methods
        Demo2();
    }
    catch(std::exception& ex)
    {
        std::cerr << ex.what() << std::endl;
        return 1;
    }
    return 0;
}
