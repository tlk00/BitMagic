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

//Reads the given file, either EUR to USD or USD to JPY in this case and stores the data
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
        
        std::string date_str, time_str;
        float open_val, high_val, low_val, close_val;
        unsigned int volume_val;
        int extra_val;

        //Extracts the data from each row and stores in a vector
        if (ss >> date_str >> time_str >> open_val >> high_val >> low_val >> close_val >> volume_val >> extra_val)
        {
            std::string full_timestamp = date_str + " " + time_str;
            
            sv_day.push_back(full_timestamp.c_str());
            o.push_back(open_val);
            h.push_back(high_val);
            l.push_back(low_val);
            c.push_back(close_val);
            
            //Percent change is the easiest way to see a simultaneous change in value
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
    
    //move information from vectors to sparse_vector_float's
    sv_open.import(o.data(), (unsigned)o.size());
    sv_high.import(h.data(), (unsigned)h.size());
    sv_low.import(o.data(), (unsigned)o.size());
    sv_close.import(h.data(), (unsigned)h.size());
    sv_volume.import(v.data(), (unsigned)v.size());
    
    //Optimizes and freezes the data since we will not be gathering any more data
    BM_DECLARE_TEMP_BLOCK(tb)
    
    sv_day.optimize(tb);
    sv_open.optimize(tb);
    sv_high.optimize(tb);
    sv_low.optimize(tb);
    sv_close.optimize(tb);
    sv_pct_change.optimize(tb);
    sv_volume.optimize(tb);
    
    sv_day.freeze();
    sv_open.freeze();
    sv_high.freeze();
    sv_low.freeze();
    sv_close.freeze();
    sv_pct_change.freeze();
    sv_volume.freeze();
}

//Finds when the USD price fluctuated by comparing JPY and EUR prices
sparseVecFloat::bvector_type find_inverse_pct_matches(const sparseVecFloat& eur_pct, const sparseVecFloat& jpy_pct,
                                                      float target_pct)
{
    sparseVecFloat::bvector_type eur_scan;
    sparseVecFloat::bvector_type jpy_scan;
    
    /*Finds every percent greater than or equal to target_pct
     -target_pct for jpy since the data is stored in USD to JPY instead of EUR to USD, so it is an inverse
     */
    sparseVecScanner scanner;
    scanner.find_ge_float(eur_pct, target_pct, eur_scan);
    scanner.find_le_float(jpy_pct, -1 * target_pct, jpy_scan);
    
    //and operation checks removes any fluctuations in EUR or JPY independely, so any percent change is likely due to USD price fluctuation
    eur_scan &= jpy_scan;
    
    return eur_scan;
}

int main(void){
    try
    {
        //Stores data from EUR dataset
        sparseVecString  eur_day;
        sparseVecFloat   eur_open;
        sparseVecFloat   eur_high;
        sparseVecFloat   eur_low;
        sparseVecFloat   eur_pct_change;
        sparseVecFloat   eur_close;
        sparseVecUint    eur_volume;

        fillSparseVecs(eur_day, eur_open, eur_high, eur_low, eur_pct_change, eur_close, eur_volume, "EURUSD_H1.csv");
        
        
        //Stores data from JPY dataset
        sparseVecString  jpy_day;
        sparseVecFloat   jpy_open;
        sparseVecFloat   jpy_high;
        sparseVecFloat   jpy_low;
        sparseVecFloat   jpy_pct_change;
        sparseVecFloat   jpy_close;
        sparseVecUint    jpy_volume;

        fillSparseVecs(jpy_day, jpy_open, jpy_high, jpy_low, jpy_pct_change, jpy_close, jpy_volume, "USDJPY_H1.csv");

        //makes sure they are the same size
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
        
        //Calculates the statistics of the stored variables
        {
            std::cout << "Storing All of this data in std::vector<float>'s size: " << sizeof(float) * jpy_close.size() * 10 << " bytes" << std::endl;
            
            sparseVecFloat::statistics stat_eur_open, stat_eur_high, stat_eur_low, stat_eur_pct, stat_eur_close;
            sparseVecFloat::statistics stat_jpy_open, stat_jpy_high, stat_jpy_low, stat_jpy_pct, stat_jpy_close;

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

            size_t total_eur_mem = stat_eur_open.memory_used + stat_eur_high.memory_used +
                                   stat_eur_low.memory_used  + stat_eur_pct.memory_used +
                                   stat_eur_close.memory_used;

            size_t total_jpy_mem = stat_jpy_open.memory_used + stat_jpy_high.memory_used +
                                   stat_jpy_low.memory_used  + stat_jpy_pct.memory_used +
                                   stat_jpy_close.memory_used;

            const int w = 16;
            
            std::cout << std::setw(w) << "Vector Name" << " | " << std::setw(w-1) << "Memory Used (B)" << std::endl;

            // EUR
            std::cout << std::setw(w) << "eur_open" << " | " << stat_eur_open.memory_used << "\n";
            std::cout << std::setw(w) << "eur_high" << " | " << stat_eur_high.memory_used << "\n";
            std::cout << std::setw(w) << "eur_low"  << " | " << stat_eur_low.memory_used  << "\n";
            std::cout << std::setw(w) << "eur_pct_change" << " | " << stat_eur_pct.memory_used << "\n";
            std::cout << std::setw(w) << "eur_close" << " | " << stat_eur_close.memory_used << "\n";
            std::cout << "--------------------------------------------------------" << std::endl;
            
            // JPY
            std::cout << std::setw(w) << "jpy_open" << " | " << stat_jpy_open.memory_used << "\n";
            std::cout << std::setw(w) << "jpy_high" << " | " << stat_jpy_high.memory_used << "\n";
            std::cout << std::setw(w) << "jpy_low"  << " | " << stat_jpy_low.memory_used  << "\n";
            std::cout << std::setw(w) << "jpy_pct_change" << " | " << stat_jpy_pct.memory_used << "\n";
            std::cout << std::setw(w) << "jpy_close" << " | " << stat_jpy_close.memory_used << "\n";
            std::cout << "--------------------------------------------------------" << std::endl;
            
            std::cout << "Total EUR sparse_vector_float memory usage: " << total_eur_mem << " bytes\n";
            std::cout << "Total JPY sparse_vector_float memory usage: " << total_jpy_mem << " bytes\n" << std::endl;
            std::cout << "Total sparse_vector_float memory usage: " << total_jpy_mem + total_eur_mem << " bytes\n" << std::endl;
            
        }
        
        //Finds when the USD fluctuated by 1 percent down
        sparseVecFloat::bvector_type sameChanges = find_inverse_pct_matches(eur_pct_change, jpy_pct_change, 1.0f);
        
        sparseVecFloat::bvector_type::enumerator parser = sameChanges.first();
        
        std::cout << "\n\n--- Found " << sameChanges.count() << " Same Percent Changes ---\n";
        
        int w = 18;
        
        //Prints the day and time the USD fluctuated
        while (parser.valid())
        {
            unsigned int idx = *parser;

            std::string timestamp = eur_day[idx].get();
            float eur_change     = eur_pct_change.get(idx);
            float jpy_change     = jpy_pct_change.get(idx);

            std::cout << "Row [" << idx << "] | Time: " << timestamp << "\n"
                      << " ├─ EUR Pct: " << std::setw(w) << std::left << (std::to_string(eur_change) + "%")
                      << " | Current Close: " << std::setw(w) << std::left << eur_close.get(idx)
                      << " | Last Close: " << eur_close.get(idx-1) << "\n"
                      << " └─ JPY Pct: " << std::setw(w) << std::left << (std::to_string(jpy_change) + "%")
                      << " | Current Close: " << std::setw(w) << std::left << jpy_close.get(idx)
                      << " | Last Close: " << jpy_close.get(idx-1) << std::endl << std::endl;

            ++parser;
        }
        
        
        //Shows the serialized size of all the sparse_vector_floats
        {
            //Create layouts
            bm::sparse_vector_float_serial_layout<sparseVecFloat> layout_eur_open, layout_eur_high, layout_eur_low, layout_eur_pct, layout_eur_close;
            bm::sparse_vector_float_serial_layout<sparseVecFloat> layout_jpy_open, layout_jpy_high, layout_jpy_low, layout_jpy_pct, layout_jpy_close;

            //Serializes the EUR data
            bm::sparse_vector_float_serialize(eur_open,       layout_eur_open);
            bm::sparse_vector_float_serialize(eur_high,       layout_eur_high);
            bm::sparse_vector_float_serialize(eur_low,        layout_eur_low);
            bm::sparse_vector_float_serialize(eur_pct_change, layout_eur_pct);
            bm::sparse_vector_float_serialize(eur_close,      layout_eur_close);

            //Serializes the JPY data
            bm::sparse_vector_float_serialize(jpy_open,       layout_jpy_open);
            bm::sparse_vector_float_serialize(jpy_high,       layout_jpy_high);
            bm::sparse_vector_float_serialize(jpy_low,        layout_jpy_low);
            bm::sparse_vector_float_serialize(jpy_pct_change, layout_jpy_pct);
            bm::sparse_vector_float_serialize(jpy_close,      layout_jpy_close);

            size_t size_eur_open  = layout_eur_open.size();
            size_t size_eur_high  = layout_eur_high.size();
            size_t size_eur_low   = layout_eur_low.size();
            size_t size_eur_pct   = layout_eur_pct.size();
            size_t size_eur_close = layout_eur_close.size();

            size_t size_jpy_open  = layout_jpy_open.size();
            size_t size_jpy_high  = layout_jpy_high.size();
            size_t size_jpy_low   = layout_jpy_low.size();
            size_t size_jpy_pct   = layout_jpy_pct.size();
            size_t size_jpy_close = layout_jpy_close.size();

            size_t total_serialized_bytes = size_eur_open + size_eur_high + size_eur_low + size_eur_pct + size_eur_close +
                                            size_jpy_open + size_jpy_high + size_jpy_low + size_jpy_pct + size_jpy_close;

            //Prints the size in bytes of the serialized data
            std::cout << "\n\n         SERIALIZED DATA BLOB SIZES (ON-DISK)          " << std::endl << std::endl;
            std::cout << "EUR Open Serialized Size:       " << size_eur_open  << " bytes\n";
            std::cout << "EUR High Serialized Size:       " << size_eur_high  << " bytes\n";
            std::cout << "EUR Low Serialized Size:        " << size_eur_low   << " bytes\n";
            std::cout << "EUR Pct Change Serialized Size: " << size_eur_pct   << " bytes\n";
            std::cout << "EUR Close Serialized Size:      " << size_eur_close << " bytes\n";
            std::cout << "--------------------------------------------------------" << std::endl;
            std::cout << "JPY Open Serialized Size:       " << size_jpy_open  << " bytes\n";
            std::cout << "JPY High Serialized Size:       " << size_jpy_high  << " bytes\n";
            std::cout << "JPY Low Serialized Size:        " << size_jpy_low   << " bytes\n";
            std::cout << "JPY Pct Change Serialized Size: " << size_jpy_pct   << " bytes\n";
            std::cout << "JPY Close Serialized Size:      " << size_jpy_close << " bytes\n";
            std::cout << "--------------------------------------------------------" << std::endl;
            std::cout << "COMBINED TOTAL SERIALIZED SIZE: " << total_serialized_bytes << " bytes" << std::endl;
        }
    }
    catch(std::exception& ex)
    {
        std::cerr << ex.what() << std::endl;
        return 1;
    }
    return 0;
}
