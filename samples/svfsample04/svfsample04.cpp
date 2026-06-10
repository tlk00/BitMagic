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

/** \example svfsample04.cpp
  Example of how to use bm::sparse_vector_float template class to extract elements into a normal array and how to
  use a back insert iterator

  bm::sparse_vector_float is a wrapper of bm::sparse_vector<> that takes in floats and converts them into
  a bitvector of signs, and 2 sparse_vectors of the floats exponent and mantissa
 
  \sa bm::sparse_vector_float::decode
  \sa bm::sparse_vector_float::extract
  \sa bm::sparse_vector_float::extract_range
  \sa bm::sparse_vector_float::gather

  \sa bm::sparse_vector_float::back_insert_iterator::operator=
  \sa bm::sparse_vector_float::back_insert_iterator::add
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
    auto floatEq = [](float a, float b) {
        return std::fabs(a - b) < 0.001f;
    };

    int N = 128000;
    float m = 0.5f;
    bm::sparse_vector_float<bm::bvector<>> svf1;
    std::vector<float> temp(N*2);

    for(int i = 0; i < N; i++){
        temp[i] = (i * 0.001) * m;
    }
    for(int i = N; i < N*2; i++){
        temp[i] = -1*(i * 0.001) * m;
    }
    BM_DECLARE_TEMP_BLOCK(tb)
    svf1.import(temp.data(), N*2);
    svf1.optimize(tb);


    //You can take out the items in an sparse_vector_float by using decode and giving them an array to put
    //the data into, a vector in this example, where to start, and how many elements to extract, and whether the
    //array is preinitialized, considered true by default
    std::vector<float> svfExtract(N*2);
    //svf1.extract(svfExtract.data() svf1.size()); - You can also use extract, but that just calls decode
    svf1.decode(svfExtract.data(), 0, N*2, true);

    int errorCount = 0;
    for (int i = 0; i < N*2; i++) {
        if (!floatEq(svfExtract[i], temp[i])){
            errorCount++;
        }
    }
    assert(errorCount == 0);

    //You can take out a portion of the svf into an array by using extract_range and giving it the
    //array to extract into, how many elements to extract, and where to start extracting
    std::vector<float> svfExtractRange(48000);
    svf1.extract_range(svfExtractRange.data(), 48000, 16000);

    errorCount = 0;
    for (int i = 16000; i < 64000; i++) {
        if (!floatEq(svfExtractRange[i-16000], temp[i])){
            errorCount++;
        }
    }
    assert(errorCount == 0);


    bm::id_t gatherIndeces[1024];
    for(int i = 0; i < 1024; i++){
        gatherIndeces[i] = rand() % 128000;
    }

    //The last way to extract data into an array is to use gather, which takes in an array of which indeces to
    //extract and extracts them into the given array, and how many of those indeces to extract. You can also
    //let it know if the indeces are sorted or not, which can improve performance
    std::vector<float> svfGather(1024);
    svf1.gather(svfGather.data(), gatherIndeces, 1024, bm::BM_UNKNOWN);

    errorCount = 0;
    for (int i = 0; i < 1024; i++) {
        if (!floatEq(svfGather[i], temp[gatherIndeces[i]])){
            errorCount++;
        }
    }
    assert(errorCount == 0);

}

void Demo2(){
    bm::sparse_vector_float<bm::bvector<>> svf1;

    //You can create a back insert iterator for sparse_vector_float like this
    bm::sparse_vector_float<bm::bvector<>>::back_insert_iterator testBI(&svf1);

    //You can add elements to the sparse_vector using add or the = operator
    //This calles push_back on the svf
    testBI.add(1.0023);
    testBI.add(400005.6);
    testBI=78.9;

    std::cout << "svf1.size() = " << svf1.size() << std::endl;
    std::cout << "svf1.get(0) = " << svf1.get(0) << std::endl;
    std::cout << "svf1.get(1) = " << svf1.get(1) << std::endl;
    std::cout << "svf1.get(2) = " << svf1.get(2) << std::endl;

    //You can copy a back insert iterator or move it like this, the move disconnects the old inserter
    bm::sparse_vector_float<bm::bvector<>>::back_insert_iterator testBI2(testBI);
    bm::sparse_vector_float<bm::bvector<>>::back_insert_iterator testBI3(std::move(testBI));
    testBI2=100;
    testBI3=12345.6789;
    std::cout << "svf1.size() = " << svf1.size() << std::endl;
    std::cout << "svf1.get(3) = " << svf1.get(3) << std::endl;
    std::cout << "svf1.get(4) = " << svf1.get(4) << std::endl;
}

int main(void){
    try
    {
        //Demo1 for extraction operations
        Demo1();

        std::cout << std::endl << std::endl;

        //Demo2 for back inserter operations
        Demo2();
    }
    catch(std::exception& ex)
    {
        std::cerr << ex.what() << std::endl;
        return 1;
    }
    return 0;
}
