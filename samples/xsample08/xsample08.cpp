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

/** \example xsample08.cpp

  Example on intervals and how to use it for layout calculation


  \sa bm::bvector::set_range
  \sa bm::bvector::any_range
  \sa bm::interval_enumerator

  \sa sample22.cpp
  \sa sample23.cpp
  \sa bvintervals
*/

/*! \file sample23.cpp
    \brief Example: interval_enumerator<> - interator class for intervals
*/

#include <iostream>
#include <utility>
#include <vector>
#include <memory>
#include <cassert>

#include "bm.h"
#include "bmintervals.h"

using namespace std;

typedef bm::interval_enumerator<bm::bvector<> > interval_enumerator_type;

typedef std::vector<std::unique_ptr<bm::bvector<> > > layout_vector_type;

struct data_model
{
    layout_vector_type lv_;  ///< layout vector
    layout_vector_type stv_; ///< strand vector
};

void set_strand(data_model& dm, size_t plane, unsigned start, bool strand)
{
    if (!strand)
        return;
    while (dm.stv_.size() <= plane)
        dm.stv_.emplace_back(std::unique_ptr<bm::bvector<> >(nullptr));

    bm::bvector<>* bv;
    bv = dm.stv_[plane].get();
    if (!bv)
    {
        bv = new bm::bvector<>(bm::BM_GAP);
        dm.stv_[plane] = std::unique_ptr<bm::bvector<> >(bv);
    }
    bv->set_bit(start);
}

void add_object(data_model& dm, unsigned start, unsigned end, bool strand)
{
    assert(start <= end);

    bm::bvector<>* bv; // layout plane vector

    for (size_t i  = 0; i < dm.lv_.size(); ++i)
    {
        bv = dm.lv_[i].get();
        if (!bv)
        {
            bv = new bm::bvector<>(bm::BM_GAP);
            dm.lv_[i] = std::unique_ptr<bm::bvector<> >(bv);
            bv->set_range(start, end);
            set_strand(dm, i, start, strand);
            return;
        }
        if (!bv->any_range(start, end))
        {
            bv->set_range(start, end);
            set_strand(dm, i, start, strand);
            return;
        }
    } // for i

    // not found, make new plane
    //
    bv = new bm::bvector<>(bm::BM_GAP);
    dm.lv_.emplace_back(std::unique_ptr<bm::bvector<> >(bv));
    bv->set_range(start, end);
    set_strand(dm, dm.lv_.size()-1, start, strand);
}

void print_model(const data_model& dm)
{
    const bm::bvector<>* bv; // layout
    const bm::bvector<>* bv_strand; // strand
    for (size_t i  = 0; i < dm.lv_.size(); ++i)
    {
        bv = dm.lv_[i].get();
        if (bv)
        {
            bv_strand = i < dm.stv_.size() ? dm.stv_[i].get() : nullptr;
            interval_enumerator_type ien(*bv);
            if (ien.valid())
            {
                bm::bvector<>::size_type spaces = 0;
                do
                {
                    auto st = ien.start();

                    char ch_strand = '<';
                    if (bv_strand && bv_strand->test(st))
                        ch_strand = '>';

                    auto end = ien.end();
                    bool first = true;
                    for (; spaces < st; ++spaces)
                        cout << " ";
                    for (; st <= end; ++st)
                    {
                        if (st == end)
                            cout << ch_strand;
                        else
                            cout << (first ? ch_strand : '.');
                        first = false;
                    } // for
                    spaces = end+1;
                } while (ien.advance());
                cout << endl;
            }
        }
    } // for
}


int main(void)
{
    try
    {
        data_model dm;

        add_object(dm, 0, 0, false);
        add_object(dm, 5, 10, false);
        add_object(dm, 4, 70, true);
        add_object(dm, 15, 20, true);
        add_object(dm, 20, 30, true);
        add_object(dm, 16, 21, false);

        print_model(dm);
    }
    catch(std::exception& ex)
    {
        std::cerr << ex.what() << std::endl;
        return 1;
    }
        
    return 0;
}

