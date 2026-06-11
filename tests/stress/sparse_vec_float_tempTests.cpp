#include <iostream>
#include <fstream>
#include <cassert>

#include <bm.h>
#include <bmsparsevec_float.h>
#include <bmsparsevec_float_serial.h>
#include <bmsparsevec_algo.h>
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

void SparseVecFloatSerializeTest(){
    auto floatEq = [](float a, float b) {
        return std::fabs(a - b) < 0.001f;
    };

    float toAdd[] = {1.0123, -2.468, 340000.56};

    bm::sparse_vector_float<bm::bvector<>> testSVF;
    testSVF.import(toAdd, 3);
    BM_DECLARE_TEMP_BLOCK(tb)
    testSVF.optimize(tb);

    bm::sparse_vector_float_serial_layout<bm::sparse_vector_float<bm::bvector<>>> testLayout;
    
    bm::sparse_vector_float_serialize(testSVF, testLayout);

    const unsigned char* buf = testLayout.buf();
    bm::sparse_vector_float_deserialize(testSVF, buf);
    
    assert(testSVF.size() == 3);
    assert(floatEq(testSVF.get(0), toAdd[0]));
    assert(floatEq(testSVF.get(1), toAdd[1]));
    assert(floatEq(testSVF.get(2), toAdd[2]));

    bm::sparse_vector_float<bm::bvector<>> testSVF2;
    int N = 10000;
    for(int i = 0; i < N; i++){
        float f = i * 0.000123;
        testSVF2.push_back(f);
    }

    testSVF2.optimize(tb);
    bm::sparse_vector_float_serial_layout<bm::sparse_vector_float<bm::bvector<>>> testLayout2;
    bm::sparse_vector_float_serialize(testSVF2, testLayout2);

    buf = testLayout2.buf();
    bm::sparse_vector_float_deserializer<bm::sparse_vector_float<bm::bvector<>>> testDeserializer;
    bm::sparse_vector_float<bm::bvector<>> testSVF2_restored;
    testDeserializer.deserialize_range(testSVF2_restored, buf, 300, 400, true);
    
    int errorCount = 0;
    for (int i = 300; i <= 400; i++) {
        float f = i * 0.000123;
        if (!floatEq(testSVF2_restored.get(i), f)){
            errorCount++;
        }
    }
    assert(errorCount == 0);

    bm::bvector<> mask_bv;
    int maskIndices[] = {0, 1, 50, 100, 500, 999, 5000, 9999};
    int maskSize = sizeof(maskIndices) / sizeof(maskIndices[0]);
    for (int i = 0; i < maskSize; i++)
        mask_bv.set(maskIndices[i]);
    
    bm::sparse_vector_float<bm::bvector<>> testSVF2_masked;
    testDeserializer.deserialize(testSVF2_masked, buf, mask_bv);

    errorCount = 0;
    for (int i = 0; i < maskSize; i++) {
        int idx   = maskIndices[i];
        float f   = idx * 0.000123f;
        if (!floatEq(testSVF2_masked.get(idx), f))
            errorCount++;
    }
    assert(errorCount == 0);

    assert(floatEq(testSVF2_masked.get(2),    0.0f));
    assert(floatEq(testSVF2_masked.get(200),  0.0f));
    assert(floatEq(testSVF2_masked.get(1000), 0.0f));
}

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
        }
    }
    assert(errorCount == 0);
}

void SparseVecFloatBackInsertTests(){

    auto floatEq = [](float a, float b) {
        return std::fabs(a - b) < 0.001f;
    };

    bm::sparse_vector_float<bm::bvector<>> testSVF;

    bm::sparse_vector_float<bm::bvector<>>::back_insert_iterator testBI(&testSVF);

    testBI.add(1.0023);
    assert(testSVF.size() == 1);
    assert(floatEq(testSVF.get(0), 1.0023));

    testBI.add(400005.6);
    assert(testSVF.size() == 2);
    assert(floatEq(testSVF.get(1), 400005.6));

    bm::sparse_vector_float<bm::bvector<>>::back_insert_iterator testBI2(testBI);
    testBI2=78.9;
    assert(testSVF.size() == 3);
    assert(floatEq(testSVF.get(2), 78.9));

    bm::sparse_vector_float<bm::bvector<>>::back_insert_iterator testBI3(std::move(testBI));
    testBI3=12345.6789;
    assert(testSVF.size() == 4);
    assert(floatEq(testSVF.get(3), 12345.6789));

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

        bm::sparse_vector_float_serial_layout<bm::sparse_vector_float<bm::bvector<>>> testLayout;
        bm::sparse_vector_float_serialize(sv, testLayout);

        // deserialize into a fresh vector
        bm::sparse_vector_float<bm::bvector<>> sv_restored;
        const unsigned char* buf = testLayout.buf();
        bm::sparse_vector_float_deserialize(sv_restored, buf);

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

        bm::sparse_vector_float_serial_layout<bm::sparse_vector_float<bm::bvector<>>> testLayout;
        bm::sparse_vector_float_serialize(sv, testLayout);

        // deserialize into a fresh vector
        bm::sparse_vector_float<bm::bvector<>> sv_restored;
        const unsigned char* buf = testLayout.buf();
        bm::sparse_vector_float_deserialize(sv_restored, buf);

        assert((int)sv_restored.size() == N);
        int errorCount = 0;
        for (int i = 0; i < N; i++) {
            if (std::fabs(sv_restored.get(i) - data[i]) > 0.001f)
                errorCount++;
        }
        assert(errorCount == 0);

        // check that serialized size is smaller than raw size
        size_t rawSize        = N * sizeof(float);
        size_t serializedSize = testLayout.size();
        assert(serializedSize < rawSize);
    }
}

void SparseVecFloatTests(){

    SparseVecFloatGeneralTests();
    SparseVecFloatConstIteratorTests();
    SparseVecFloatImportTest();
    SparseVecFloatSerializeTest();
    SparseVecFloatRangeTests();
    SparseVecFloatExtractionTests();
    SparseVecFloatBackInsertTests();
    std::cout << "Sparse Vector Float Tests Complete" << std::endl;
}

struct SplitFloat {
    unsigned int sign;
    unsigned int exponent;
    unsigned int mantissa;
};

SplitFloat split_float(float f) {
    unsigned int bits;
    std::memcpy(&bits, &f, sizeof(float));
    return {
        (bits >> 31) & 0x1,
        (bits >> 23) & 0xFF,
        bits         & 0x7FFFFF
    };
}

typedef bm::sparse_vector_float<bm::bvector<>> sparse_vector_float;
typedef bm::sparse_vector_scanner<bm::sparse_vector<unsigned int, bm::bvector<>>> svfScanner;

void in_pos_range(sparse_vector_float sv, SplitFloat fromSplit, SplitFloat toSplit, bm::bvector<> &bv_out){
    typename sparse_vector_float::size_type sz = sv.size();

    svfScanner scan;
    scan.find_range(sv.exponents_, fromSplit.exponent, toSplit.exponent, bv_out);
    bv_out -= sv.signs_;

    if(bv_out.none()) return;

    bm::bvector<> onBounds;
    scan.find_eq(sv.exponents_, fromSplit.exponent, onBounds);
    onBounds &= bv_out;
    

    if (onBounds.any()) {
        bm::bvector<> invalid_min_mantissas;
        scan.find_lt(sv.mantissas_, fromSplit.mantissa, invalid_min_mantissas);
        
        onBounds &= invalid_min_mantissas;
        bv_out -= onBounds;
    }

    onBounds.clear();
    scan.find_eq(sv.exponents_, toSplit.exponent, onBounds);
    onBounds &= bv_out;

    if (onBounds.any()) {
        bm::bvector<> invalid_max_mantissas;
        scan.find_gt(sv.mantissas_, toSplit.mantissa, invalid_max_mantissas);
        
        onBounds &= invalid_max_mantissas;
        bv_out -= onBounds;
    }
}

void in_neg_range(sparse_vector_float sv, SplitFloat fromSplit, SplitFloat toSplit, bm::bvector<> &bv_out){
    typename sparse_vector_float::size_type sz = sv.size();

    svfScanner scan;
    scan.find_range(sv.exponents_, fromSplit.exponent, toSplit.exponent, bv_out);
    bv_out &= sv.signs_;
    
    if(bv_out.none()) return;

    bm::bvector<> onBounds;
    scan.find_eq(sv.exponents_, fromSplit.exponent, onBounds);
    onBounds &= bv_out;

    if (onBounds.any()) {
        bm::bvector<> invalid_min_mantissas;
        scan.find_gt(sv.mantissas_, fromSplit.mantissa, invalid_min_mantissas);
        
        onBounds &= invalid_min_mantissas;
        bv_out -= onBounds;
    }

    onBounds.clear();
    scan.find_eq(sv.exponents_, toSplit.exponent, onBounds);
    onBounds &= bv_out;

    if (onBounds.any()) {
        bm::bvector<> invalid_max_mantissas;
        scan.find_lt(sv.mantissas_, toSplit.mantissa, invalid_max_mantissas);
        
        onBounds &= invalid_max_mantissas;
        bv_out -= onBounds;
    }
}

void in_arb_range(sparse_vector_float sv, SplitFloat fromSplit, SplitFloat toSplit, bm::bvector<> &bv_out){
    typename sparse_vector_float::size_type sz = sv.size();

    svfScanner scan;
    bm::bvector<> neg_exp_range;
    
    scan.find_range(sv.exponents_, 0, fromSplit.exponent, neg_exp_range);
    neg_exp_range &= sv.signs_;

    if (neg_exp_range.any()) {
        bm::bvector<> neg_bounds;
        scan.find_eq(sv.exponents_, fromSplit.exponent, neg_bounds);
        neg_bounds &= neg_exp_range;

        if (neg_bounds.any()) {
            bm::bvector<> invalid_neg_mantissas;
            scan.find_gt(sv.mantissas_, fromSplit.mantissa, invalid_neg_mantissas);
            
            neg_bounds &= invalid_neg_mantissas;
            neg_exp_range -= neg_bounds;
        }
    }

    bm::bvector<> pos_exp_range;
    
    scan.find_range(sv.exponents_, 0, toSplit.exponent, pos_exp_range);
    pos_exp_range -= sv.signs_;

    if (pos_exp_range.any()) {
        bm::bvector<> pos_bounds;
        scan.find_eq(sv.exponents_, toSplit.exponent, pos_bounds);
        pos_bounds &= pos_exp_range;

        if (pos_bounds.any()) {
            bm::bvector<> invalid_pos_mantissas;
            scan.find_gt(sv.mantissas_, toSplit.mantissa, invalid_pos_mantissas);
            
            pos_bounds &= invalid_pos_mantissas;
            pos_exp_range -= pos_bounds; // Drop the invalid ones
        }
    }

    bv_out = std::move(neg_exp_range);
    bv_out |= pos_exp_range;
}

void in_range(sparse_vector_float sv, float from, float to, bm::bvector<> &bv_out){
    
    if(from > to)
        std::swap(to, from);

    SplitFloat fromSplit = split_float(from);
    SplitFloat toSplit = split_float(to);

    bv_out.clear();

    if(fromSplit.sign == 0 && toSplit.sign == 0){
        in_pos_range(sv, fromSplit, toSplit, bv_out);
    }else if(fromSplit.sign == 1 && toSplit.sign == 1){
        in_neg_range(sv, fromSplit, toSplit, bv_out);
    }else{
        in_arb_range(sv, fromSplit, toSplit, bv_out);
    }
}

void TestScanner(){
    typedef bm::sparse_vector_float<bm::bvector<>> sparse_vector_float;
    BM_DECLARE_TEMP_BLOCK(tb)

    auto floatEq = [](float a, float b) {
        return std::fabs(a - b) < 0.0001f;
    };

    auto print_time = [](const char* test_name, std::chrono::high_resolution_clock::time_point start,
                                                std::chrono::high_resolution_clock::time_point end)
    {
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        std::cout << test_name << ": (" << duration / 1000.0 << " ms)\n";
    };

    int N = 20000000;
    std::vector<float> temp(N);

    for(int i = 0; i < N/2; i++)
        temp[i] = -1.0f * i * 0.00123f;
    for(int i = 0; i < N/2; i++)
        temp[i+N/2] = i * 0.00123f;

    sparse_vector_float testSVF;
    testSVF.import(temp.data(), N);
    testSVF.optimize(tb);

    bm::bvector<> bv_out;

    // --- Test 1: large positive range ---
    {
        float from = 100.0f;
        float to   = 500.0f;

        auto start = std::chrono::high_resolution_clock::now();
        in_range(testSVF, from, to, bv_out);
        auto end = std::chrono::high_resolution_clock::now();
        print_time("Test 1 - large positive range [100, 500]", start, end);
        std::cout << "hits: " << bv_out.count() << "\n";

        bm::bvector<>::enumerator en     = bv_out.first();
        bm::bvector<>::enumerator en_end = bv_out.end();
        for (; en != en_end; ++en)
        {
            float v = temp[*en];
            assert(v >= from && v <= to);
        }
        for (int i = 0; i < N; i++)
        {
            if (temp[i] >= from && temp[i] <= to){
                assert(bv_out.test(i));
            }
        }
        bv_out.clear();
    }

    // --- Test 2: small positive range ---
    {
        float from = 1.0f;
        float to   = 2.0f;

        auto start = std::chrono::high_resolution_clock::now();
        in_range(testSVF, from, to, bv_out);
        auto end = std::chrono::high_resolution_clock::now();
        print_time("Test 2 - small positive range [1, 2]", start, end);
        std::cout << "hits: " << bv_out.count() << "\n";

        bm::bvector<>::enumerator en     = bv_out.first();
        bm::bvector<>::enumerator en_end = bv_out.end();
        for (; en != en_end; ++en)
        {
            float v = temp[*en];
            assert(v >= from && v <= to);
        }
        for (int i = 0; i < N; i++)
        {
            if (temp[i] >= from && temp[i] <= to)
                assert(bv_out.test(i));
        }
        bv_out.clear();
    }

    // --- Test 3: large negative range ---
    {
        float from = -500.0f;
        float to   = -100.0f;

        auto start = std::chrono::high_resolution_clock::now();
        in_range(testSVF, from, to, bv_out);
        auto end = std::chrono::high_resolution_clock::now();
        print_time("Test 3 - large negative range [-500, -100]", start, end);
        std::cout << "hits: " << bv_out.count() << "\n";

        bm::bvector<>::enumerator en     = bv_out.first();
        bm::bvector<>::enumerator en_end = bv_out.end();
        for (; en != en_end; ++en)
        {
            float v = temp[*en];
            assert(v >= from && v <= to);
        }
        for (int i = 0; i < N; i++)
        {
            if (temp[i] >= from && temp[i] <= to)
                assert(bv_out.test(i));
        }
        bv_out.clear();
    }

    // --- Test 4: small negative range ---
    {
        float from = -2.0f;
        float to   = -1.0f;

        auto start = std::chrono::high_resolution_clock::now();
        in_range(testSVF, from, to, bv_out);
        auto end = std::chrono::high_resolution_clock::now();
        print_time("Test 4 - small negative range [-2, -1]", start, end);
        std::cout << "hits: " << bv_out.count() << "\n";

        bm::bvector<>::enumerator en     = bv_out.first();
        bm::bvector<>::enumerator en_end = bv_out.end();
        for (; en != en_end; ++en)
        {
            float v = temp[*en];
            assert(v >= from && v <= to);
        }
        for (int i = 0; i < N; i++)
        {
            if (temp[i] >= from && temp[i] <= to)
                assert(bv_out.test(i));
        }
        bv_out.clear();
    }

    // --- Test 5: range crossing zero ---
    {
        float from = -10.0f;
        float to   =  10.0f;

        auto start = std::chrono::high_resolution_clock::now();
        in_range(testSVF, from, to, bv_out);
        auto end = std::chrono::high_resolution_clock::now();
        print_time("Test 5 - range crossing zero [-10, 10]", start, end);
        std::cout << "hits: " << bv_out.count() << "\n";

        bm::bvector<>::enumerator en     = bv_out.first();
        bm::bvector<>::enumerator en_end = bv_out.end();
        for (; en != en_end; ++en)
        {
            float v = temp[*en];
            assert(v >= from && v <= to);
        }
        for (int i = 0; i < N; i++)
        {
            if (temp[i] >= from && temp[i] <= to)
                assert(bv_out.test(i));
        }
        bv_out.clear();
    }

    // --- Test 6: range with no results ---
    {
        float from = 99999.0f;
        float to   = 99999.9f;

        auto start = std::chrono::high_resolution_clock::now();
        in_range(testSVF, from, to, bv_out);
        auto end = std::chrono::high_resolution_clock::now();
        print_time("Test 6 - empty range [99999, 99999.9]", start, end);
        std::cout << "hits: " << bv_out.count() << "\n";

        assert(bv_out.count() == 0);
        bv_out.clear();
    }

    // --- Test 7: tiny range near zero ---
    {
        float from = -0.001f;
        float to   =  0.001f;

        auto start = std::chrono::high_resolution_clock::now();
        in_range(testSVF, from, to, bv_out);
        auto end = std::chrono::high_resolution_clock::now();
        print_time("Test 7 - tiny range near zero [-0.001, 0.001]", start, end);
        std::cout << "hits: " << bv_out.count() << "\n";

        bm::bvector<>::enumerator en     = bv_out.first();
        bm::bvector<>::enumerator en_end = bv_out.end();
        for (; en != en_end; ++en)
        {
            float v = temp[*en];
            assert(v >= from && v <= to);
        }
        for (int i = 0; i < N; i++)
        {
            if (temp[i] >= from && temp[i] <= to)
                assert(bv_out.test(i));
        }
        bv_out.clear();
    }

    // --- Test 8: inverted range ---
    {
        float from = 10.0f;
        float to   =  5.0f;

        auto start = std::chrono::high_resolution_clock::now();
        in_range(testSVF, from, to, bv_out);
        auto end = std::chrono::high_resolution_clock::now();
        print_time("Test 8 - inverted range (10, 5) -> swapped to [5, 10]", start, end);
        std::cout << "hits: " << bv_out.count() << "\n";

        bm::bvector<>::enumerator en     = bv_out.first();
        bm::bvector<>::enumerator en_end = bv_out.end();
        for (; en != en_end; ++en)
        {
            float v = temp[*en];
            assert(v >= 5.0f && v <= 10.0f);
        }
        bv_out.clear();
    }

    // --- also time a naive float comparison for reference ---
    {
        float from = 100.0f;
        float to   = 500.0f;

        auto start = std::chrono::high_resolution_clock::now();
        bm::bvector<> naive_out;
        for (int i = 0; i < N; i++)
        {
            if (temp[i] >= from && temp[i] <= to)
                naive_out.set(i);
        }
        auto end = std::chrono::high_resolution_clock::now();
        print_time("Naive std::vector float comparison [100, 500] (reference)", start, end);
        std::cout << "  hits: " << naive_out.count() << "\n";
    }
}

int main(int argc, char *argv[]){
    //SparseVecFloatTests();
    //SparseVecFloatStressTests();
    //SparseVecFloatSerialStressTests();
    TestScanner();
}
