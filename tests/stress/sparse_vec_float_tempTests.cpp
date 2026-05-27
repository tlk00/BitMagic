#include <iostream>
#include <fstream>
#include <cassert>

#include <bm.h>
#include <bmsparsevec_float_serial.h>

void SparseVecFloatTests();
void SparseVecFloatGeneralTests();
void SparseVecFloatConstIteratorTests();
void SparseVecFloatImportTest();
void SparseVecFloatSerializeTest();

int main(int argc, char *argv[]){
    SparseVecFloatTests();
}

void SparseVecFloatTests(){

    SparseVecFloatGeneralTests();
    SparseVecFloatConstIteratorTests();
    SparseVecFloatImportTest();
    std::cout << "Sparse Vector Float Tests Complete" << std::endl;
}

void SparseVecFloatConstIteratorTests(){
    
    float toAdd[] = {1.0123, 2.468, 340000.56};

    bm::sparse_vector_float testSVF;
    testSVF.import(toAdd, 3);

    auto floatEq = [](float a, float b) {
        return std::fabs(a - b) < 0.001f;
    };
    
    // --- construction ---
    bm::sparse_vector_float::const_iterator defaultit;
    bm::sparse_vector_float::const_iterator itFromSV(&testSVF);
    bm::sparse_vector_float::const_iterator itFromSVPos(&testSVF, 1);
    bm::sparse_vector_float::const_iterator itBegin = testSVF.begin();
    bm::sparse_vector_float::const_iterator itEnd   = testSVF.end();
    bm::sparse_vector_float::const_iterator itCopy(itBegin);
    
    // --- operator* and value() ---
    assert(floatEq(*itBegin, toAdd[0]));
    assert(floatEq(itBegin.value(), toAdd[0]));
    assert(floatEq(*itFromSVPos, toAdd[1]));

    // --- valid() ---
    assert(itBegin.valid());
    assert(!itEnd.valid());
    
    // --- invalidate() ---
    bm::sparse_vector_float::const_iterator itInvalid = testSVF.begin();
    itInvalid.invalidate();
    assert(!itInvalid.valid());

    // --- pos() ---
    assert(itBegin.pos() == 0);
    assert(itFromSVPos.pos() == 1);

    // --- operator== and operator!= ---
    bm::sparse_vector_float::const_iterator itA = testSVF.begin();
    bm::sparse_vector_float::const_iterator itB = testSVF.begin();
    assert(itA == itB);
    assert(!(itA != itB));
    assert(itA != itEnd);
    assert(!(itA == itEnd));

    // --- operator< <= > >= ---
    bm::sparse_vector_float::const_iterator itFirst  = testSVF.begin();
    bm::sparse_vector_float::const_iterator itSecond(&testSVF, 1);
    assert(itFirst  <  itSecond);
    assert(itFirst  <= itSecond);
    assert(itFirst  <= itFirst);
    assert(itSecond >  itFirst);
    assert(itSecond >= itFirst);
    assert(itFirst  >= itFirst);

    // --- prefix operator++ ---
    bm::sparse_vector_float::const_iterator itPre = testSVF.begin();
    ++itPre;
    assert(itPre.pos() == 1);
    assert(floatEq(*itPre, toAdd[1]));
    ++itPre;
    assert(itPre.pos() == 2);
    assert(floatEq(*itPre, toAdd[2]));

    // --- postfix operator++ ---
    bm::sparse_vector_float::const_iterator itPost = testSVF.begin();
    bm::sparse_vector_float::const_iterator itPostOld = itPost++;
    assert(itPostOld.pos() == 0);
    assert(itPost.pos() == 1);
    assert(floatEq(*itPostOld, toAdd[0]));
    assert(floatEq(*itPost, toAdd[1]));

    // --- advance() ---
    bm::sparse_vector_float::const_iterator itAdv = testSVF.begin();
    assert(floatEq(itAdv.value(), toAdd[0]));
    bool stillValid = itAdv.advance();
    assert(stillValid);
    assert(itAdv.pos() == 1);
    assert(floatEq(itAdv.value(), toAdd[1]));
    stillValid = itAdv.advance();
    assert(stillValid);
    assert(itAdv.pos() == 2);
    assert(floatEq(itAdv.value(), toAdd[2]));
    bool pastEnd = itAdv.advance();
    assert(!pastEnd);
    assert(!itAdv.valid());

    // --- go_to() ---
    bm::sparse_vector_float::const_iterator itGoto = testSVF.begin();
    itGoto.go_to(2);
    assert(itGoto.pos() == 2);
    assert(floatEq(itGoto.value(), toAdd[2]));
    itGoto.go_to(0);
    assert(itGoto.pos() == 0);
    assert(floatEq(itGoto.value(), toAdd[0]));
    itGoto.go_to(1);
    assert(itGoto.pos() == 1);
    assert(floatEq(itGoto.value(), toAdd[1]));

    // --- is_null() ---
    bm::sparse_vector_float::const_iterator itNull = testSVF.begin();
    assert(!itNull.is_null());

    // --- full iteration ---
    int idx = 0;
    for (auto it = testSVF.begin(); it != testSVF.end(); ++it, ++idx) {
        assert(floatEq(*it, toAdd[idx]));
    }
    assert(idx == 3);
}

void SparseVecFloatImportTest(){
    int N = 128000;
    float m = 0.5f;
    bm::sparse_vector_float testSVF;
    std::vector<float> temp(N);

    for(int i = 0; i < N; i++){
        temp[i] = (i * 0.001) * m;
    }
    
    testSVF.import(temp.data(), N);

    int errorCount;
    for(int i = 0; i < N; i++){
        float err = std::fabs(temp[i] - testSVF.get(i));
        if (err > 0.0f) {
            errorCount++;
        }
    }

    assert(errorCount == 0);

    BM_DECLARE_TEMP_BLOCK(tb)
    testSVF.optimize();

    errorCount = 0;
    for(int i = 0; i < N; i++){
        float err = std::fabs(temp[i] - testSVF.get(i));
        if (err > 0.0f) {
            errorCount++;
        }
    }

    assert(errorCount == 0);
}

void SparseVecFloatGeneralTests(){

    auto floatEq = [](float a, float b) {
        return std::fabs(a - b) < 0.001f;
    };

    float toAdd[] = {1.0123, 2.468, 340000.56};

    bm::sparse_vector_float testSVF;
    

    assert(testSVF.empty());

    testSVF.push_back(toAdd[0]);

    assert(!testSVF.empty());
    assert(testSVF.size() == 1);
    assert(floatEq(testSVF.get(0), toAdd[0]));

    testSVF.push_back(toAdd[1]);
    testSVF.push_back(toAdd[2]);

    assert(testSVF.size() == 3);
    assert(floatEq(testSVF.get(0), toAdd[0]));
    assert(floatEq(testSVF.get(1), toAdd[1]));
    assert(floatEq(testSVF.get(2), toAdd[2]));

    bm::sparse_vector_float testSVF2;
    testSVF2.import(toAdd, 3);
    assert(testSVF  == testSVF2);
    assert(!(testSVF != testSVF2));

    float toAdd2[] = {9.0f, 8.0f, 7.0f};
    bm::sparse_vector_float testSVF3;
    testSVF3.import(toAdd2, 3);
    assert(testSVF  != testSVF3);
    assert(!(testSVF == testSVF3));

    bm::sparse_vector_float testSVFAssigned;
    testSVFAssigned = testSVF;
    assert(testSVFAssigned == testSVF);
    assert(floatEq(testSVFAssigned.get(0), toAdd[0]));
    assert(floatEq(testSVFAssigned.get(1), toAdd[1]));
    assert(floatEq(testSVFAssigned.get(2), toAdd[2]));
    assert(testSVFAssigned.size() == testSVF.size());

    testSVF.set(1, 8.258f);
    assert(!floatEq(testSVF.get(1), toAdd[1]));
    assert(floatEq(testSVF.get(1), 8.258f));

    testSVF.set(100, 100.001f);
    assert(testSVF.size() == 101);
    assert(floatEq(testSVF.get(100), 100.001f));
    assert(floatEq(testSVF.get(50), 0.0f));

    bm::sparse_vector_float svA;
    bm::sparse_vector_float svB;
    float aVals[] = {1.0f, 2.0f, 3.0f};
    float bVals[] = {4.0f, 5.0f, 6.0f};
    svA.import(aVals, 3);
    svB.import(bVals, 3);

    svA.swap(svB);

    assert(floatEq(svA.get(0), bVals[0]));
    assert(floatEq(svA.get(1), bVals[1]));
    assert(floatEq(svA.get(2), bVals[2]));
    assert(floatEq(svB.get(0), aVals[0]));
    assert(floatEq(svB.get(1), aVals[1]));
    assert(floatEq(svB.get(2), aVals[2]));

    svA.swap(svB);
    assert(floatEq(svA.get(0), aVals[0]));
    assert(floatEq(svA.get(1), aVals[1]));
    assert(floatEq(svA.get(2), aVals[2]));
}

void SparseVecFloatSerializeTest(){
    auto floatEq = [](float a, float b) {
        return std::fabs(a - b) < 0.001f;
    };

    float toAdd[] = {1.0123, 2.468, 340000.56};

    bm::sparse_vector_float testSVF;
    testSVF.import(toAdd, 3);

    bm::sparse_vector_float_serialized testSVFSerial;
    testSVFSerial.serialize(testSVF);

    testSVFSerial.deserialize(testSVF);

    assert(testSVF.size() == 3);
    assert(floatEq(testSVF.get(0), toAdd[0]));
    assert(floatEq(testSVF.get(1), toAdd[1]));
    assert(floatEq(testSVF.get(2), toAdd[2]));
}
