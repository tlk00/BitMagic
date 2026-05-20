#include <iostream>
#include <fstream>
#include <cassert>

#include <bm.h>
#include <bmsparsevec_float.h>


int main(int argc, char *argv[]){
    SparseVecFloatTests();
}

void SparseVecFloatTests(){

    SparseVecFloatGeneralTests();

    /*Tests to be added:
        - const_iterator tests
        - operator

        - swap
    */
}

void SparseVecFloatGeneralTests(){

    float toAdd[] = {1.0123, 2.468, 340000.56};

    bm::sparse_vector_float testSVF;

    assert(testSVF.empty());

    testSVF.push_back(toAdd[0]);

    assert(!testSVF.empty());
    assert(testSVF.size() == 1);
    assert(testSVF.get(0) - toAdd[0] == 0);

    testSVF.push_back(toAdd[1]);
    testSVF.push_back(toAdd[2]);

    assert(testSVF.size() == 3);
    assert(testSVF.get(0) - toAdd[0] == 0);
    assert(testSVF.get(1) - toAdd[1] == 0);
    assert(testSVF.get(2) - toAdd[2] == 0);
}
