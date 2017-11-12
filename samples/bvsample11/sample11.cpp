/*
Copyright(c) 2002-2017 Anatoliy Kuznetsov(anatoliy_kuznetsov at yahoo.com)

Permission is hereby granted, free of charge, to any person
obtaining a copy of this software and associated documentation
files (the "Software"), to deal in the Software without restriction,
including without limitation the rights to use, copy, modify, merge,
publish, distribute, sublicense, and/or sell copies of the Software,
and to permit persons to whom the Software is furnished to do so,
subject to the following conditions:

The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
OTHER DEALINGS IN THE SOFTWARE.

You have to explicitly mention BitMagic project in any derivative product,
its WEB Site, published materials, articles or any other work derived from this
project or based on our code or know-how.

For more information please visit:  http://bitmagic.io

*/

/** \example sample11.cpp
  Demo on memory consumption control techniques
 
 */


#include <iostream>
#include <memory>
#include <map>
#include <vector>

#include "bm.h"
#include "bmrandom.h"

// ----------------------------------------------------
// Global parameters
// ----------------------------------------------------

// Number of vectors generated for the test
const unsigned index_size = 1000000;

// Dynamic range for constructed sets
const unsigned max_size = 2000000;

// Number of bits per one vector
const unsigned bits_per_vect = 5;

// bit-vector type for this example
typedef  bm::bvector<>   TBVector;



/* BitMagic provides two GAP length tables for situations when we have
   standard or embarassingly sparse vectors.
   bm::gap_len_table - default standard
   bm::gap_len_table_min - option for smaller vectors
 
   Here we define an alternative table for very sparse vectors
*/
template<bool T> struct gap_len_table_sparse
{
    static const bm::gap_word_t _len[bm::gap_levels];
};
template<bool T>
const bm::gap_word_t gap_len_table_sparse<T>::_len[bm::gap_levels] =
                                    { 8, 32, 128, 512 };


// simple bit-vector class factory for the project
//
TBVector* construct_bvector()
{
    // in this example we plan to keep lots of vectors in memory, thus
    // use parameters to minimize memory consumption
    //
    TBVector* bv =
        new TBVector(bm::BM_GAP,  // use GAP compressed mode
                     gap_len_table_sparse<true>::_len, // custom lens for super sparse vectors
                     max_size  // limit the maximum size
                   );
    return bv;
}

// Generic utility to destroy map of pointers
template<typename TM>
void destroy_map(TM& id_map)
{
    for (typename TM::iterator it = id_map.begin();
         it != id_map.end();
         ++it)
    {
        typename TM::mapped_type mp = it->second;
        delete mp;
    } // for
    id_map.clear();
}

// ------------------------------------------------------------------
// Sample data structures
// ------------------------------------------------------------------

// Sample index structure to keep a map of in-memory bit-vectors
//
struct bv_index
{
    typedef std::map<unsigned, TBVector*> map_type;
    ~bv_index()
    {
        destroy_map(idx_);
    }
    map_type  idx_;
};

// Sample index structure to keep map of in-memory serialized/compressed bit-vectors
//
struct bvs_index
{
    typedef std::vector<unsigned char>           buffer_type;
    typedef std::map<unsigned, buffer_type>      map_type;
    
    map_type  idx_;
};

// --------------------------------------------------------------------


// set bits in a vector using two methods picked at random
// one method will generate a plato of non-random integers,
// the other adds ints randomly without following any system
//
void generate_random_vector(TBVector* bv)
{
    unsigned method = rand()%2; // pick a generation method
    if (method == 0) // generate a incremental linear sequence at random location
    {
        unsigned seed_id = rand() % max_size;
        for (unsigned i  = seed_id; i < bits_per_vect; ++i)
        {
            if (i >= max_size)
                break;
            bv->set_bit(i);
        } // for i
    }
    else // generate a few random bits
    {
        for (unsigned i  = 0; i < bits_per_vect; ++i)
        {
            unsigned id = rand() % max_size;
            if (i >= max_size) // paranoiya check
                break;
            bv->set_bit(id);
        } // for i
    }
}

// generate map of bit-vectors, each filled with just a few bits
//
void generate_bv_index(bv_index& bvi)
{
    for (unsigned i = 0; i < index_size; ++i)
    {
        std::auto_ptr<TBVector> ap(construct_bvector());
        
        generate_random_vector(ap.get());
        bvi.idx_[i] = ap.release();
    }
}

// calculate memory footprint for in memory index
//
size_t calc_memory_footprint(const bv_index& bvi)
{
    size_t mem_total = 0;
    for (bv_index::map_type::const_iterator it = bvi.idx_.begin();
         it != bvi.idx_.end();
         ++it)
    {
        const TBVector* mp = it->second;
        TBVector::statistics st;
        mp->calc_stat(&st);
        mem_total += st.memory_used;
    } // for
    
    return mem_total;
}


int main(void)
{
    try
    {
        bv_index bvi;
        
        generate_bv_index(bvi);
        
        size_t mem_total = calc_memory_footprint(bvi);
        
        std::cout << "bvector index memory footprint = " << mem_total << std::endl;
        
        
        getchar();

    }
    catch(std::exception& ex)
    {
        std::cerr << ex.what() << std::endl;
        return 1;
    }

    return 0;
}

