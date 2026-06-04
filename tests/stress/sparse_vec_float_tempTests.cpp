#include <iostream>
#include <fstream>
#include <cassert>

#include <bm.h>
#include <bmsparsevec_float.h>
//#include <bmsparsevec_float_serial.h>
#include <random>


//----------------------------------------------------------------

void SparseVecFloatConstIteratorTests(){
    
    float toAdd[] = {1.0123, 2.468, 340000.56};

    bm::sparse_vector_float<bm::bvector<>> testSVF;
    testSVF.import(toAdd, 3);

    auto floatEq = [](float a, float b) {
        return std::fabs(a - b) < 0.001f;
    };
    
    // --- construction ---
    bm::sparse_vector_float<bm::bvector<>>::const_iterator defaultit;
    bm::sparse_vector_float<bm::bvector<>>::const_iterator itFromSV(&testSVF);
    bm::sparse_vector_float<bm::bvector<>>::const_iterator itFromSVPos(&testSVF, 1);
    bm::sparse_vector_float<bm::bvector<>>::const_iterator itBegin = testSVF.begin();
    bm::sparse_vector_float<bm::bvector<>>::const_iterator itEnd   = testSVF.end();
    bm::sparse_vector_float<bm::bvector<>>::const_iterator itCopy(itBegin);
    
    // --- operator* and value() ---
    assert(floatEq(*itBegin, toAdd[0]));
    assert(floatEq(itBegin.value(), toAdd[0]));
    assert(floatEq(*itFromSVPos, toAdd[1]));

    // --- valid() ---
    assert(itBegin.valid());
    assert(!itEnd.valid());
    
    // --- invalidate() ---
    bm::sparse_vector_float<bm::bvector<>>::const_iterator itInvalid = testSVF.begin();
    itInvalid.invalidate();
    assert(!itInvalid.valid());

    // --- pos() ---
    assert(itBegin.pos() == 0);
    assert(itFromSVPos.pos() == 1);

    // --- operator== and operator!= ---
    bm::sparse_vector_float<bm::bvector<>>::const_iterator itA = testSVF.begin();
    bm::sparse_vector_float<bm::bvector<>>::const_iterator itB = testSVF.begin();
    assert(itA == itB);
    assert(!(itA != itB));
    assert(itA != itEnd);
    assert(!(itA == itEnd));

    // --- operator< <= > >= ---
    bm::sparse_vector_float<bm::bvector<>>::const_iterator itFirst  = testSVF.begin();
    bm::sparse_vector_float<bm::bvector<>>::const_iterator itSecond(&testSVF, 1);
    assert(itFirst  <  itSecond);
    assert(itFirst  <= itSecond);
    assert(itFirst  <= itFirst);
    assert(itSecond >  itFirst);
    assert(itSecond >= itFirst);
    assert(itFirst  >= itFirst);

    // --- prefix operator++ ---
    bm::sparse_vector_float<bm::bvector<>>::const_iterator itPre = testSVF.begin();
    ++itPre;
    assert(itPre.pos() == 1);
    assert(floatEq(*itPre, toAdd[1]));
    ++itPre;
    assert(itPre.pos() == 2);
    assert(floatEq(*itPre, toAdd[2]));

    // --- postfix operator++ ---
    bm::sparse_vector_float<bm::bvector<>>::const_iterator itPost = testSVF.begin();
    bm::sparse_vector_float<bm::bvector<>>::const_iterator itPostOld = itPost++;
    assert(itPostOld.pos() == 0);
    assert(itPost.pos() == 1);
    assert(floatEq(*itPostOld, toAdd[0]));
    assert(floatEq(*itPost, toAdd[1]));

    // --- advance() ---
    bm::sparse_vector_float<bm::bvector<>>::const_iterator itAdv = testSVF.begin();
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
    bm::sparse_vector_float<bm::bvector<>>::const_iterator itGoto = testSVF.begin();
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
    bm::sparse_vector_float<bm::bvector<>>::const_iterator itNull = testSVF.begin();
    assert(!itNull.is_null());

    // --- full iteration ---
    int idx = 0;
    for (auto it = testSVF.begin(); it != testSVF.end(); ++it, ++idx) {
        assert(floatEq(*it, toAdd[idx]));
    }
    assert(idx == 3);
}

//----------------------------------------------------------------

void SparseVecFloatImportTest(){
    int N = 128000;
    float m = 0.5f;
    bm::sparse_vector_float<bm::bvector<>> testSVF;
    std::vector<float> temp(N*2);

    for(int i = 0; i < N; i++){
        temp[i] = (i * 0.001) * m;
    }
    for(int i = N; i < N*2; i++){
        temp[i] = -1*(i * 0.001) * m;
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
    testSVF.optimize(tb);
    
    errorCount = 0;
    for(int i = 0; i < N; i++){
        float err = std::fabs(temp[i] - testSVF.get(i));
        if (err > 0.0f) {
            errorCount++;
        }
    }
    
    assert(errorCount == 0);
}

//----------------------------------------------------------------

void SparseVecFloatGeneralTests(){
    auto floatEq = [](float a, float b) {
        return std::fabs(a - b) < 0.001f;
    };

    float toAdd[] = {1.0123, -2.468, 340000.56};

    bm::sparse_vector_float<bm::bvector<>> testSVF;
    
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

    bm::sparse_vector_float<bm::bvector<>> testSVF2;
    testSVF2.import(toAdd, 3);
    assert(testSVF  == testSVF2);
    assert(!(testSVF != testSVF2));

    float toAdd2[] = {9.0f, -8.0f, -7.0f};
    bm::sparse_vector_float<bm::bvector<>> testSVF3;
    testSVF3.import(toAdd2, 3);
    assert(testSVF  != testSVF3);
    assert(!(testSVF == testSVF3));

    bm::sparse_vector_float<bm::bvector<>> testSVFAssigned;
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

    bm::sparse_vector_float<bm::bvector<>> svA;
    bm::sparse_vector_float<bm::bvector<>> svB;
    float aVals[] = {1.0f, -2.0f, 3.0f};
    float bVals[] = {-4.0f, -5.0f, 6.0f};
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

//----------------------------------------------------------------
/*
void SparseVecFloatSerializeTest(){
    auto floatEq = [](float a, float b) {
        return std::fabs(a - b) < 0.001f;
    };

    float toAdd[] = {1.0123, -2.468, 340000.56};

    bm::sparse_vector_float<bm::bvector<>> testSVF;
    testSVF.import(toAdd, 3);
    BM_DECLARE_TEMP_BLOCK(tb)
    testSVF.optimize(tb);

    bm::sparse_vector_float_serialized<bm::bvector<>> testSVFSerial;
    testSVFSerial.serialize(testSVF);

    testSVFSerial.deserialize(testSVF);

    assert(testSVF.size() == 3);
    assert(floatEq(testSVF.get(0), toAdd[0]));
    assert(floatEq(testSVF.get(1), toAdd[1]));
    assert(floatEq(testSVF.get(2), toAdd[2]));
}*/

//----------------------------------------------------------------

void SparseVecFloatRangeTests(){
    auto floatEq = [](float a, float b) {
        return std::fabs(a - b) < 0.001f;
    };

    float toAdd[] = {1.0123, -2.468, 340000.56, -7008.0, 0.900102};

    bm::sparse_vector_float<bm::bvector<>> testSVF;
    testSVF.import(toAdd, 5);
    BM_DECLARE_TEMP_BLOCK(tb)
    testSVF.optimize(tb);

    assert(testSVF.size() == 5);
    testSVF.clear();
    assert(testSVF.size() == 0);

    testSVF.import(toAdd, 5);
    testSVF.clear_range(1, 3);
    assert(testSVF.size() == 5);
    assert(floatEq(testSVF.get(0), toAdd[0]));
    assert(floatEq(testSVF.get(1), 0.0));
    assert(floatEq(testSVF.get(2), 0.0));
    assert(floatEq(testSVF.get(3), 0.0));
    assert(floatEq(testSVF.get(4), toAdd[4]));

    testSVF.clear();
    testSVF.import(toAdd, 5);

    bm::sparse_vector_float<bm::bvector<>> testSVF2(testSVF);
    assert(testSVF.equal(testSVF2));

    testSVF2.set(1, 0.0);
    assert(!testSVF.equal(testSVF2));

    assert(testSVF.compare(1, 2.468) == -1);
    assert(testSVF.compare(1, -2.468) == 0);
    assert(testSVF.compare(1, -3) == 1);

    bm::sparse_vector_float<bm::bvector<>> svf1;
    bm::sparse_vector_float<bm::bvector<>> svf2;
    float toAdd1[] = {1.0123, -2.468, 0.0, 0.0, 0.0, 1.5};
    float toAdd2[] = {0.0, 0.0, 0.0, -7008.0, 0.900102, 2.5};
    svf1.import(toAdd1, 6);
    svf2.import(toAdd2, 6);
    svf1.optimize(tb);
    svf2.optimize(tb);
    
    svf1.join(svf2);
    assert(svf1.size() == 6);
    assert(floatEq(svf1.get(0), toAdd[0]));
    assert(floatEq(svf1.get(1), toAdd[1]));
    assert(floatEq(svf1.get(2), 0.0));
    assert(floatEq(svf1.get(3), toAdd[3]));
    assert(floatEq(svf1.get(4), toAdd[4]));
    assert(svf1.get(5) != svf1.get(5));

    svf1.clear();
    svf2.clear();
    svf1.import(toAdd1, 6);
    svf2.import(toAdd2, 6);
    svf1.optimize(tb);
    svf2.optimize(tb);

    svf1.merge(svf2);
    assert(svf1.size() == 6);
    assert(floatEq(svf1.get(0), toAdd[0]));
    assert(floatEq(svf1.get(1), toAdd[1]));
    assert(floatEq(svf1.get(2), 0.0));
    assert(floatEq(svf1.get(3), toAdd[3]));
    assert(floatEq(svf1.get(4), toAdd[4]));
    assert(svf1.get(5) != svf1.get(5));


    svf1.clear();
    svf2.clear();
    svf1.import(toAdd1, 6);
    svf2.import(toAdd2, 6);
    svf1.optimize(tb);
    svf2.optimize(tb);

    svf1.copy_range(svf2, 2, 4);
    assert(svf1.size() == 6);
    assert(floatEq(svf1.get(0), 0.0));
    assert(floatEq(svf1.get(1), 0.0));
    assert(floatEq(svf1.get(2), 0.0));
    assert(floatEq(svf1.get(3), toAdd[3]));
    assert(floatEq(svf1.get(4), toAdd[4]));
    assert(floatEq(svf1.get(5), 0.0));
}

//----------------------------------------------------------------

void SparseVecFloatExtractionTests(){

    auto floatEq = [](float a, float b) {
        return std::fabs(a - b) < 0.001f;
    };

    int N = 128000;
    float m = 0.5f;
    bm::sparse_vector_float<bm::bvector<>> testSVF;
    std::vector<float> temp(N*2);

    for(int i = 0; i < N; i++){
        temp[i] = (i * 0.001) * m;
    }
    for(int i = N; i < N*2; i++){
        temp[i] = -1*(i * 0.001) * m;
    }

    testSVF.import(temp.data(), N*2);
    testSVF.optimize();

    std::vector<float> testExtract(N*2);
    testSVF.decode(testExtract.data(), 0, N*2);

    int errorCount = 0;
    for (int i = 0; i < N*2; i++) {
        if (!floatEq(testExtract[i], temp[i])){
            errorCount++;
        }
    }
    assert(errorCount == 0);

    
    std::vector<float> testExtractRange(48000);
    testSVF.extract_range(testExtractRange.data(), 48000, 16000);

    errorCount = 0;
    for (int i = 16000; i < 64000; i++) {
        if (!floatEq(testExtractRange[i-16000], temp[i])){
            errorCount++;
        }
    }
    assert(errorCount == 0);

    bm::id_t gatherIndeces[1024];
    for(int i = 0; i < 1024; i++){
        gatherIndeces[i] = rand() % 128000;
    }

    std::vector<float> testGather(1024);
    testSVF.gather(testGather.data(), gatherIndeces, 1024, bm::BM_UNKNOWN);

    errorCount = 0;
    for (int i = 0; i < 1024; i++) {
        if (!floatEq(testGather[i], temp[gatherIndeces[i]])){
            errorCount++;
            std::cout << "Mismatch at sample " << i 
                  << " (Index " << gatherIndeces[i] << "): "
                  << temp[gatherIndeces[i]] << " vs " << testGather[i] 
                  << std::endl;
        }
    }
    assert(errorCount == 0);
}

//----------------------------------------------------------------
//perf

void SparseVecFloatStressTests(){
    //Random data import and get test
    {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<float> dist(-1000000.0f, 1000000.0f);

        int N = 1000000;
        std::vector<float> data(N);
        for (int i = 0; i < N; i++)
            data[i] = dist(gen);

        bm::sparse_vector_float<bm::bvector<>> sv;
        sv.import(data.data(), N);

        BM_DECLARE_TEMP_BLOCK(tb)
        sv.optimize(tb);

        int errorCount = 0;
        for (int i = 0; i < N; i++) {
            float restored = sv.get(i);
            if (std::fabs(restored - data[i]) > 0.001f){
                errorCount++;
            }
        }
        assert(errorCount == 0);
        assert((int)sv.size() == N);
    }

    //Random data push back and iterator tests
    {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<float> dist(-500.0f, 500.0f);

        int N = 500000;
        std::vector<float> data(N);
        bm::sparse_vector_float<bm::bvector<>> sv;
        for (int i = 0; i < N; i++){
            data[i] = dist(gen);
            sv.push_back(data[i]);
        }

        BM_DECLARE_TEMP_BLOCK(tb)
        sv.optimize(tb);

        // validate via iterator
        int idx = 0;
        int errorCount = 0;
        for (auto it = sv.begin(); it != sv.end(); ++it, ++idx) {
            if (std::fabs(*it - data[idx]) > 0.001f)
                errorCount++;
        }
        assert(errorCount == 0);
        assert(idx == N);
    }

    //Decaying Sinusoid Test
    {
        int N = 2000000;
        std::vector<float> data(N);
        const double e = std::exp(1.0);
        for (int i = 0; i < N; i++) {
            float t = 0.001f * i;
            data[i] = 100.0f * (float)pow(e, -1.0 * t)
                      * (cos(5.0f * t + 1.0f) + sin(5.0f * t + 1.0f));
        }

        bm::sparse_vector_float<bm::bvector<>> sv;
        sv.import(data.data(), N);
        BM_DECLARE_TEMP_BLOCK(tb)
        sv.optimize(tb);

        int errorCount = 0;
        for (int i = 0; i < N; i++) {
            if (std::fabs(sv.get(i) - data[i]) > 0.001f)
                errorCount++;
        }
        assert((int)sv.size() == N);
    }

}
/*
void SparseVecFloatSerialStressTests(){
    //Serializing then deserializing a large random data set
    {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<float> dist(-1000000.0f, 1000000.0f);

        int N = 1000000;
        std::vector<float> data(N);
        for (int i = 0; i < N; i++)
            data[i] = dist(gen);

        bm::sparse_vector_float<bm::bvector<>> sv;
        sv.import(data.data(), N);
        BM_DECLARE_TEMP_BLOCK(tb)
        sv.optimize(tb);

        // serialize
        bm::sparse_vector_float_serialized<bm::bvector<>> serial;
        serial.serialize(sv);

        // deserialize into a fresh vector
        bm::sparse_vector_float<bm::bvector<>> sv_restored;
        serial.deserialize(sv_restored);

        // validate
        assert((int)sv_restored.size() == N);
        int errorCount = 0;
        for (int i = 0; i < N; i++) {
            if (std::fabs(sv_restored.get(i) - data[i]) > 0.001f)
                errorCount++;
        }
        assert(errorCount == 0);
    }

    // Decaying Sinusoid Compression Test
    {
        int N = 2000000;
        std::vector<float> data(N);
        const double e = std::exp(1.0);
        for (int i = 0; i < N; i++) {
            float t = 0.001f * i;
            data[i] = 100.0f * (float)pow(e, -1.0 * t)
                      * (cos(5.0f * t + 1.0f) + sin(5.0f * t + 1.0f));
        }

        bm::sparse_vector_float<bm::bvector<>> sv;
        sv.import(data.data(), N);
        BM_DECLARE_TEMP_BLOCK(tb)
        sv.optimize(tb);

        bm::sparse_vector_float_serialized<bm::bvector<>> serial;
        serial.serialize(sv);

        // deserialize and validate
        bm::sparse_vector_float<bm::bvector<>> sv_restored;
        serial.deserialize(sv_restored);

        assert((int)sv_restored.size() == N);
        int errorCount = 0;
        for (int i = 0; i < N; i++) {
            if (std::fabs(sv_restored.get(i) - data[i]) > 0.001f)
                errorCount++;
        }
        assert(errorCount == 0);

        // check that serialized size is smaller than raw size
        size_t rawSize        = N * sizeof(float);
        size_t serializedSize = serial.size();
        assert(serializedSize < rawSize);
    }
}*/

void SparseVecFloatTests(){

    SparseVecFloatGeneralTests();
    SparseVecFloatConstIteratorTests();
    SparseVecFloatImportTest();

    SparseVecFloatRangeTests();
    SparseVecFloatExtractionTests();
    std::cout << "Sparse Vector Float Tests Complete" << std::endl;
}

int main(int argc, char *argv[]){
    SparseVecFloatTests();
    //SparseVecFloatStressTests();
    //SparseVecFloatSerialStressTests();
}
