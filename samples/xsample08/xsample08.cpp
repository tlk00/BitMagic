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

  Example on intervals and how to use it for layout calculation.
  As a use case this example uses genomics visualization for
  features mapped into genomic coordinates.

  It is also illustartes vector model using coordinate ranges
  or feature vectors. Various properties of the initial model acn be
  dropped (sliced) to improve memory efficiency, better storage
  or network transfer.

  This example does NOT do serialization of models (which is possible)
  for the clarity of the sample code.

  \sa bm::bvector::set_range
  \sa bm::bvector::any_range
  \sa bm::bvector::copy_range
  \sa bm::interval_enumerator
  \sa bm::rsc_sparse_vector
  \sa bm::rsc_sparse_vector::copy_range
  \sa bm::find_interval_start
  \sa bm::find_interval_end

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
#include "bmsparsevec_compr.h"
#include "bmundef.h" /* clear the pre-proc defines from BM */

using namespace std;

typedef bm::interval_enumerator<bm::bvector<> > interval_enumerator_type;
typedef std::vector<std::unique_ptr<bm::bvector<> > > layout_vector_type;

typedef bm::sparse_vector<unsigned char, bm::bvector<> >       sparse_vector_u8;
typedef bm::rsc_sparse_vector<unsigned char, sparse_vector_u8> rsc_vector_u8;
typedef std::vector<std::unique_ptr<rsc_vector_u8> >  starnds_vector_type;


// -------------------------------------------------------------------

/// Data frame object, sued to buid succinct data model
///
///
struct data_model
{
    /// Optimize memory layoput, build index for faster read access
    ///
    void optimize();

    void add_layout(size_t plane, bm::bvector<>* bv);
    void add_strand(size_t plane, rsc_vector_u8* strand);

    layout_vector_type     layout_v;   ///< layout vector
    starnds_vector_type    strand_v;  ///< strand planes vector
};

void data_model::optimize()
{
    BM_DECLARE_TEMP_BLOCK(tb); // explicit temp for faster optimization
    for (size_t i  = 0; i < layout_v.size(); ++i)
    {
        auto bv = layout_v[i].get();
        if (bv)
            bv->optimize(tb); // memory optimization
    } // for i
    for (size_t i  = 0; i < strand_v.size(); ++i)
    {
        auto strand_plane = strand_v[i].get();
        if (strand_plane)
        {
            strand_plane->optimize(tb);
            strand_plane->sync(); // build rank-select idx (faster read access)
        }
    } // for i
}

void data_model::add_layout(size_t plane, bm::bvector<>* bv)
{
    unique_ptr<bm::bvector<> > ap(bv);
    if (layout_v.size() == plane) // push back requested
    {
        layout_v.emplace_back(move(ap));
    }
    else
    {
        while (layout_v.size() < plane) // this is crude resize() but it would do
            layout_v.emplace_back(new bm::bvector<>(bm::BM_GAP));
        layout_v[plane] = std::move(ap);
    }
}

void data_model::add_strand(size_t plane, rsc_vector_u8* strand)
{
    unique_ptr<rsc_vector_u8 > ap(strand);
    if (strand_v.size() == plane) // push back requested
    {
        strand_v.emplace_back(move(ap));
    }
    else
    {
        while (strand_v.size() < plane) // this is crude resize() but it would do
            strand_v.emplace_back(new rsc_vector_u8());
        strand_v[plane] = std::move(ap);
    }
}


// -------------------------------------------------------------------

static
void set_feature_strand(data_model& dm, size_t   plane,
                        bm::bvector<>::size_type pos,
                        unsigned char strand)
{
    if (!strand)
        return;
    while (dm.strand_v.size() <= plane) // add planes
    {
        std::unique_ptr<rsc_vector_u8> p2(new rsc_vector_u8());
        dm.strand_v.emplace_back(move(p2));
    }

    rsc_vector_u8* strand_plane = dm.strand_v[plane].get();
    if (!strand_plane)
    {
        strand_plane = new rsc_vector_u8();
        dm.strand_v[plane] = unique_ptr<rsc_vector_u8 >(strand_plane);
    }
    assert(strand_plane->is_null(pos));
    strand_plane->set(pos, strand);
}

/// Register new object in the data model: [start..end] + strand
///
static
void add_object(data_model& dm,
                unsigned start, unsigned end,
                unsigned char strand)
{
    assert(start <= end);

    bm::bvector<>* bv; // layout plane vector

    for (size_t i  = 0; i < dm.layout_v.size(); ++i)
    {
        bv = dm.layout_v[i].get();
        if (!bv)
        {
            bv = new bm::bvector<>(bm::BM_GAP);
            dm.layout_v[i] = unique_ptr<bm::bvector<> >(bv);
            // bv just created (empty) no need to do range check
            bv->set_range(start, end);
            set_feature_strand(dm, i, start, strand);

            return;
        }
        if (!bv->any_range(start, end)) // check if layout space is not used
        {
            bv->set_range(start, end);  // add [start..end] coordinates
            // set strand at the start of feature
            set_feature_strand(dm, i, start, strand);

            return;
        }
    } // for i

    // not found, make new plane
    //
    bv = new bm::bvector<>(bm::BM_GAP);
    dm.layout_v.emplace_back(std::unique_ptr<bm::bvector<> >(bv));
    bv->set_range(start, end);
    set_feature_strand(dm, dm.layout_v.size()-1, start, strand);
}

/// Data model splicer
///
static
void splice_model(data_model& dm_target, const data_model& dm,
                  bm::bvector<>::size_type start,
                  bm::bvector<>::size_type end,
                  bool copy_strands)
{
    const bm::bvector<>* bv; // layout
    const rsc_vector_u8* strand_plane;

    size_t t_plane = 0;
    for (size_t i  = 0; i < dm.layout_v.size(); ++i)
    {
        bv = dm.layout_v[i].get();
        if (bv)
        {
            bm::bvector<>::size_type start_pos;
            bm::bvector<>::size_type end_pos;

            bool found = bm::find_interval_start(*bv, start, start_pos);
            if (!found)
                start_pos = start;
            found = bm::find_interval_end(*bv, end, end_pos);
            if (!found)
                end_pos = end;

            unique_ptr<bm::bvector<>> bv_ptr(new bm::bvector<>(bm::BM_GAP));
            bv_ptr->copy_range(*bv, start_pos, end_pos);
            if (bv_ptr->any()) // copy range may have ended as empty
            {
                dm_target.add_layout(t_plane, bv_ptr.release());

                // slice the strands plane (if requested)
                //
                if (copy_strands)
                {
                    if (i < dm.strand_v.size())
                    {
                        strand_plane = dm.strand_v[i].get();
                        if (strand_plane)
                        {
                            unique_ptr<rsc_vector_u8> strand_ptr(new rsc_vector_u8());
                            strand_ptr->copy_range(*strand_plane, start_pos, end_pos);
                            dm_target.add_strand(t_plane, strand_ptr.release());
                        }
                    }
                }
                ++t_plane;
            } // if any()

        } // if bv
    } // for i

}


/// This is ASCII art "renderer" for the data model.
/// illustrates how to manipulate succinct data model to create graphics
///
static
void print_model(const data_model& dm)
{
    const bm::bvector<>* bv; // layout
    const rsc_vector_u8* strand_plane;

    // Sequence on top is for purely decorative purposes
    cout <<
    "-------------------------------------------------------------------------"
    << endl <<
    "ATGTTAGCCCGCGCATATTATATATGTAGCGTATTAAGCGDGGAGATTACCCTTGCATTAGGTTANNNNNNNN"
    << endl <<
    "-------------------------------------------------------------------------"
    << endl;

    for (size_t i  = 0; i < dm.layout_v.size(); ++i)
    {
        bv = dm.layout_v[i].get();
        if (bv)
        {
            strand_plane = i < dm.strand_v.size() ? dm.strand_v[i].get() : nullptr;
            interval_enumerator_type ien(*bv);
            if (ien.valid())
            {
                bm::bvector<>::size_type spaces = 0;
                do
                {
                    auto st = ien.start(); auto end = ien.end();
                    char ch_strand = '?';
                    if (strand_plane)
                    {
                        auto strand = strand_plane->get(st);
                        switch (strand)
                        {
                        case 0: ch_strand = '>'; break; // positive
                        case 1: ch_strand = '<'; break; // negative
                        default: break; // unknown strand
                        }
                    }
                    for (; spaces < st; ++spaces)
                        cout << " ";
                    for (bool first = true; st <= end; ++st, first = false)
                    {
                        if (st == end)
                            cout << ch_strand;
                        else
                            cout << (first ? ch_strand : '.');
                    } // for
                    spaces = end+1;
                } while (ien.advance());
                cout << endl;
            }
        }
    } // for
}

enum Strand { positive=0, negative=1, unknown=2 };



int main(void)
{
    try
    {
        data_model dm;

        // build the data model using succinct vectors
        //
        add_object(dm, 0, 0,   negative);
        add_object(dm, 5, 10,  positive);
        add_object(dm, 4, 70,  negative);
        add_object(dm, 15, 20, negative);
        add_object(dm, 20, 30, positive);
        add_object(dm, 16, 21, unknown);

        dm.optimize(); // run compression and build access index

        // View the model using toy ASCII art renderer
        //
        print_model(dm);

        // create a model splice for [5..10] range
        // plus drop strand property (renderer will assume unknown)
        //
        {
            data_model dm_splice;

            splice_model(dm_splice, dm, 5, 10, false);
            dm_splice.optimize();

            cout << endl;
            print_model(dm_splice);
        }

        // create a model splice for [5..10] range
        // now WITH strand property
        //
        {
            data_model dm_splice;

            splice_model(dm_splice, dm, 5, 10, true);
            dm_splice.optimize();

            cout << endl;
            print_model(dm_splice);
        }
    }
    catch(std::exception& ex)
    {
        std::cerr << ex.what() << std::endl;
        return 1;
    }
        
    return 0;
}

