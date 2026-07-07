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

#include <iostream>
#include <fstream>
#include <cassert>
#include <sstream>
#include <string>

#include <bm.h>
#include "bmstrsparsevec.h"
#include <bmsparsevec_float.h>
#include <bmsparsevec_float_serial.h>
#include <bmsparsevec_algo.h>


typedef bm::sparse_vector<unsigned int, bm::bvector<>> sparseVecUint;
typedef bm::sparse_vector_float<sparseVecUint> sparseVecFloat;
typedef bm::str_sparse_vector<char, bm::bvector<>, 32> sparseVecString;
typedef bm::sparse_vector_scanner<sparseVecFloat> sparseVecScanner;

void fillSparseVecs(sparseVecString &sv_day,
                    sparseVecFloat &sv_open,
                    sparseVecFloat &sv_high,
                    sparseVecFloat &sv_low,
                    sparseVecFloat &sv_pct_change,
                    sparseVecFloat &sv_close,
                    sparseVecUint &sv_volume,
                    std::string filename)
{
    std::ifstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "Could not open file: " << filename << "\n";
        return;
    }
    
    std::vector<float> o;
    std::vector<float> h;
    std::vector<float> l;
    std::vector<float> c;
    std::vector<unsigned int> v;

    std::string line;
    bool isHeader = true;
    unsigned int dataIdx = 0;
    float prev_close = 0.0f;

    while (std::getline(file, line))
    {
        if (line.empty()) continue;
        
        if (isHeader)
        {
            isHeader = false;
            continue;
        }

        std::stringstream ss(line);
        
        // Target scalar variables for space/tab delimited processing
        std::string date_str, time_str;
        float open_val, high_val, low_val, close_val;
        unsigned int volume_val;
        int extra_val;

        // If your file separates date and time with a space,
        // extracting them sequentially keeps everything cleanly aligned.
        if (ss >> date_str >> time_str >> open_val >> high_val >> low_val >> close_val >> volume_val >> extra_val)
        {
            std::string full_timestamp = date_str + " " + time_str;
            
            sv_day.push_back(full_timestamp.c_str());
            o.push_back(open_val);
            h.push_back(high_val);
            l.push_back(low_val);
            c.push_back(close_val);
            
            if (dataIdx == 0)
            {
                sv_pct_change.push_back(0.0f);
            }
            else
            {
                float pct = ((close_val - prev_close) / prev_close) * 100.0f;
                sv_pct_change.push_back(pct);
            }
            
            v.push_back(volume_val);
            
            prev_close = close_val;
            dataIdx++;
        }
    }
    
    sv_open.import(o.data(), (unsigned)o.size());
    sv_high.import(h.data(), (unsigned)h.size());
    sv_low.import(o.data(), (unsigned)o.size());
    sv_close.import(h.data(), (unsigned)h.size());
    sv_volume.import(v.data(), (unsigned)v.size());
    
    BM_DECLARE_TEMP_BLOCK(tb)
    
    sv_day.optimize(tb);
    sv_open.optimize(tb);
    sv_high.optimize(tb);
    sv_low.optimize(tb);
    sv_close.optimize(tb);
    sv_volume.optimize(tb);
    
    sv_day.freeze();
    sv_open.freeze();
    sv_high.freeze();
    sv_low.freeze();
    sv_close.freeze();
    sv_volume.freeze();
}

sparseVecFloat::bvector_type find_inverse_pct_matches(const sparseVecFloat& eur_pct, const sparseVecFloat& jpy_pct,
                                                      float target_pct = 1.0f)
{
    sparseVecFloat::bvector_type eur_scan;
    sparseVecFloat::bvector_type jpy_scan;
    
    sparseVecScanner scanner;
    
    scanner.find_ge_float(eur_pct, target_pct, eur_scan);
    scanner.find_le_float(jpy_pct, -1 * target_pct, jpy_scan);
    
    eur_scan &= jpy_scan;
    
    return eur_scan;
}

int main(void){
    try
    {

        sparseVecString  eur_day;
        sparseVecFloat   eur_open;
        sparseVecFloat   eur_high;
        sparseVecFloat   eur_low;
        sparseVecFloat   eur_pct_change;
        sparseVecFloat   eur_close;
        sparseVecUint    eur_volume;

        fillSparseVecs(eur_day, eur_open, eur_high, eur_low, eur_pct_change, eur_close, eur_volume, "EURUSD_H1.csv");
        
        

        sparseVecString  jpy_day;
        sparseVecFloat   jpy_open;
        sparseVecFloat   jpy_high;
        sparseVecFloat   jpy_low;
        sparseVecFloat   jpy_pct_change;
        sparseVecFloat   jpy_close;
        sparseVecUint    jpy_volume;

        fillSparseVecs(jpy_day, jpy_open, jpy_high, jpy_low, jpy_pct_change, jpy_close, jpy_volume, "USDJPY_H1.csv");

        {
            if (eur_day.size() != jpy_day.size() || eur_pct_change.size() != jpy_pct_change.size())
            {
                std::cerr << "[Validation Error] Vector sizes do not match!\n"
                          << "EUR rows: " << eur_day.size() << " | JPY rows: " << jpy_day.size() << "\n";
            }

            bool timelines_identical = eur_day.equal(jpy_day);
            if (!timelines_identical)
            {
                std::cerr << "[Validation Error] Timestamps are misaligned across files!\n";
            }
        }
        
        {
            std::cout << "Storing All of this data in std::vector<float>'s size: " << sizeof(float) * jpy_close.size() * 10 << " bytes" << std::endl;
            
            sparseVecFloat::statistics stat_eur_open, stat_eur_high, stat_eur_low, stat_eur_pct, stat_eur_close;
            sparseVecFloat::statistics stat_jpy_open, stat_jpy_high, stat_jpy_low, stat_jpy_pct, stat_jpy_close;

            // 2. Invoke calc_stat via pointers (&) to gather inner bit-plane info
            eur_open.calc_stat(&stat_eur_open);
            eur_high.calc_stat(&stat_eur_high);
            eur_low.calc_stat(&stat_eur_low);
            eur_pct_change.calc_stat(&stat_eur_pct);
            eur_close.calc_stat(&stat_eur_close);

            jpy_open.calc_stat(&stat_jpy_open);
            jpy_high.calc_stat(&stat_jpy_high);
            jpy_low.calc_stat(&stat_jpy_low);
            jpy_pct_change.calc_stat(&stat_jpy_pct);
            jpy_close.calc_stat(&stat_jpy_close);

            // 3. Compute structural totals
            size_t total_eur_mem = stat_eur_open.memory_used + stat_eur_high.memory_used +
                                   stat_eur_low.memory_used  + stat_eur_pct.memory_used +
                                   stat_eur_close.memory_used;

            size_t total_jpy_mem = stat_jpy_open.memory_used + stat_jpy_high.memory_used +
                                   stat_jpy_low.memory_used  + stat_jpy_pct.memory_used +
                                   stat_jpy_close.memory_used;

            // 4. Print clean, padded console output
            const int w = 16;
            
            std::cout << std::setw(20) << "Vector Name"
                        << " | " << std::setw(w) << "Memory Used (B)";

            // EUR
            std::cout << std::setw(20) << "eur_open" << " | " << std::setw(w) << stat_eur_open.memory_used << "\n";
            std::cout << std::setw(20) << "eur_high" << " | " << std::setw(w) << stat_eur_high.memory_used << "\n";
            std::cout << std::setw(20) << "eur_low"  << " | " << std::setw(w) << stat_eur_low.memory_used  << "\n";
            std::cout << std::setw(20) << "eur_pct_change" << " | " << std::setw(w) << stat_eur_pct.memory_used << "\n";
            std::cout << std::setw(20) << "eur_close" << " | " << std::setw(w) << stat_eur_close.memory_used << "\n";
            std::cout << "--------------------------------------------------------\n";
            
            // JPY
            std::cout << std::setw(20) << "jpy_open" << " | " << std::setw(w) << stat_jpy_open.memory_used << "\n";
            std::cout << std::setw(20) << "jpy_high" << " | " << std::setw(w) << stat_jpy_high.memory_used << "\n";
            std::cout << std::setw(20) << "jpy_low"  << " | " << std::setw(w) << stat_jpy_low.memory_used  << "\n";
            std::cout << std::setw(20) << "jpy_pct_change" << " | " << std::setw(w) << stat_jpy_pct.memory_used << "\n";
            std::cout << std::setw(20) << "jpy_close" << " | " << std::setw(w) << stat_jpy_close.memory_used << "\n";
            
            std::cout << "\bTotal EUR Float Matrix Footprint: " << total_eur_mem << " bytes\n";
            std::cout << "Total JPY Float Matrix Footprint: " << total_jpy_mem << " bytes\n" << std::endl;
            
            
        }
        
        sparseVecFloat::bvector_type sameChanges = find_inverse_pct_matches(eur_pct_change, jpy_pct_change, 1.0f);
        
        sparseVecFloat::bvector_type::enumerator parser = sameChanges.first();
        
        std::cout << "--- Found " << sameChanges.count() << " Same Percent Changes ---\n";
        
        int w = 18;
        
        while (parser.valid())
        {
            // 1. Extract the current matching row index
            unsigned int idx = *parser;

            // 2. Use the index to pull the parallel data out of your containers
            std::string timestamp = eur_day[idx].get();
            float eur_change     = eur_pct_change.get(idx);
            float jpy_change     = jpy_pct_change.get(idx);

            // 3. Print or process your matched data point
            std::cout << "Row [" << idx << "] | Time: " << timestamp << "\n"
                      << " ├─ EUR Pct: " << std::setw(w) << std::left << (std::to_string(eur_change) + "%")
                      << " | Current Close: " << std::setw(w) << std::left << eur_close.get(idx)
                      << " | Last Close: " << eur_close.get(idx-1) << "\n"
                      << " └─ JPY Pct: " << std::setw(w) << std::left << (std::to_string(jpy_change) + "%")
                      << " | Current Close: " << std::setw(w) << std::left << jpy_close.get(idx)
                      << " | Last Close: " << jpy_close.get(idx-1) << std::endl << std::endl;

            // 4. Advance the parser to the next matching bit index
            ++parser;
        }
        
    }
    catch(std::exception& ex)
    {
        std::cerr << ex.what() << std::endl;
        return 1;
    }
    return 0;
}
