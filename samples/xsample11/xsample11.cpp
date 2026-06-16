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

    std::string line;
    int rowIdx = 0;
    unsigned int dataIdx = 0;

    while (std::getline(file, line))
    {
        if (line.empty()) continue;
        
        if (rowIdx == 0)
        {
            rowIdx++;
            continue;
        }

        // everything after is data
        std::stringstream ss(line);
        std::string cell;
        int colIdx = 0;

        while (std::getline(ss, cell, ','))
        {
            try
            {
                switch (colIdx)
                {
                    case 0: sv_day.push_back(cell.c_str());                      break;
                    case 1: sv_open.push_back(std::stof(cell));                  break;
                    case 2: sv_high.push_back(std::stof(cell));                  break;
                    case 3: sv_low.push_back(std::stof(cell));                   break;
                    case 4:
                        sv_close.push_back(std::stof(cell));
                        if (dataIdx != 0)
                        {
                            float pct = sv_close.get(dataIdx) - sv_close.get(dataIdx-1);
                            pct /= sv_close.get(dataIdx-1);
                            pct *= 100;
                            sv_pct_change.push_back(pct);
                        }
                        break;
                    case 5: sv_volume.push_back((unsigned int)std::stoul(cell)); break;
                    default: break;
                }
                if(dataIdx == 0){
                    sv_pct_change.push_back(0.0f);
                }
            }
            catch (...)
            {

            }
            colIdx++;
        }

        dataIdx++;
        rowIdx++;
    }
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

        fillSparseVecs(eur_day, eur_open, eur_high, eur_low, eur_pct_change, eur_close, eur_volume, "EURUSD_D1.csv");

        sparseVecString  jpy_day;
        sparseVecFloat   jpy_open;
        sparseVecFloat   jpy_high;
        sparseVecFloat   jpy_low;
        sparseVecFloat   jpy_pct_change;
        sparseVecFloat   jpy_close;
        sparseVecUint    jpy_volume;

        fillSparseVecs(jpy_day, jpy_open, jpy_high, jpy_low, jpy_pct_change, jpy_close, jpy_volume, "USDJPY_D1.csv");

        sparseVecScanner scanner;


        
    }
    catch(std::exception& ex)
    {
        std::cerr << ex.what() << std::endl;
        return 1;
    }
    return 0;
}
