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
  Demo and a benchmark on memory consumption control and logical operation
 
*/


#include <iostream>
#include <memory>
#include <map>
#include <vector>
#include <chrono>


#include "bm.h"
#include "bmalgo.h"
#include "bmtimer.h"
#include "bmserial.h"
#include "bmsparsevec.h"

#include "bmsparsevec_algo.h"
#include "bmsparsevec_serial.h"
#include "bmalgo_similarity.h"

#include "bmdbg.h"


// ----------------------------------------------------
// Global parameters and types
// ----------------------------------------------------

// Number of vectors generated for the test
const unsigned index_size = 1000000;

// Dynamic range for constructed sets
const unsigned max_size = 2000000;

// Number of bits per one vector
const unsigned bits_per_vect = 5;

// benchmark operation count
const unsigned benchmark_ops = 1000;

// subset of vectors used as a sample
const unsigned sample_cnt = 250;

// index values to extract
const unsigned result_set_cnt = 200;


// bit-vector type for this example
typedef  bm::bvector<>   TBVector;


// timing storage for benchmarking
bm::chrono_taker::duration_map_type  timing_map;



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

// Sample index structure to keep map of in-memory vector<unsigned int>
//
struct vect_index
{
    typedef std::vector<unsigned int>           buffer_type;
    typedef std::map<unsigned, buffer_type>      map_type;
    
    map_type  idx_;
};


// Sample index structure as in-memory sparse_vector
//
struct sparse_vect_index
{
    struct vect_addr
    {
        unsigned offset;
        unsigned size;
    };
    
    typedef bm::sparse_vector<unsigned, bm::bvector<> > sparse_vector_type;
    typedef std::map<unsigned, vect_addr>               map_type;
    
    sparse_vector_type  sv_storage_;
    map_type            idx_;
    
};


// --------------------------------------------------------------------


// set bits in a vector using two methods picked at random
// one method will generate a plato of non-random integers,
// the other adds ints randomly without following any system
//
void generate_random_vector(TBVector* bv)
{
    unsigned method = rand() % 5; // pick a generation method
    if (method == 0) // generate a incremental linear sequence at random location
    {
        unsigned seed_id = rand() % max_size;
        for (unsigned i = seed_id; i < seed_id+bits_per_vect; ++i)
        {
            if (i >= max_size)
                break;
            bv->set_bit(i);
        } // for i
    }
    else
    if (method == 1)
    {
        unsigned seed_id = rand() % max_size;
        unsigned id = seed_id;
        for (unsigned i = 0; i < bits_per_vect; ++i)
        {
            if (id >= max_size)
                break;
            bv->set_bit(id);
            id += (rand() % 10);
            if (i >= max_size)
                break;
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
        
        if (!ap->any()) // integrity check
        {
            // this should never happen
            std::cerr << "Warning. Empty vector generated!" << std::endl;
        }

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
        
        mem_total += sizeof(void*);
    } // for
    
    return mem_total;
}

// convert bit-vector index to bit-vector serialized index
//
size_t convert_bv2bvs(const bv_index& bvi, bvs_index& bvs)
{
    size_t  mem_total = 0;
    std::vector<unsigned char> buf;      // prepare a temporary buffer
    buf.reserve(1024);

    // bit-vector serializer
    // (keep it out of the serialization loop to minimize buffers re-allocations)
    //
    bm::serializer<TBVector> bvsr;
    bvsr.byte_order_serialization(false);
    bvsr.gap_length_serialization(false);
    bvsr.set_compression_level(4);
    
    for (bv_index::map_type::const_iterator it = bvi.idx_.begin();
         it != bvi.idx_.end();
         ++it)
    {
        unsigned id = it->first;
        const TBVector* bvp = it->second;
        
        TBVector::statistics st;
        bvp->calc_stat(&st);  // calculate max. serialized size
        
        buf.resize(st.max_serialize_mem); // prepare the temp buffer
        
        // run serialization, actual serialization size is expacted to be smaller
        //
        unsigned bvs_size = bvsr.serialize(*bvp, buf.data(), st.max_serialize_mem);
        
        // move from temp serialization buffer to compressed in-memory index
        //
        bvs_index::buffer_type& vbuf = bvs.idx_[id];
        vbuf.resize(bvs_size);
        ::memcpy(vbuf.data(), buf.data(), bvs_size);
        
        mem_total += bvs_size;
        
        mem_total += sizeof(std::vector<unsigned char>::size_type);

        
        // paranoia check compare source and desirialized vectors
        //
        #ifdef DEBUG
        {
            TBVector bv1;
            bm::deserialize(bv1, vbuf.data());
            if (bv1.compare(*bvp) !=0 )
            {
                std::cerr << "Deserialization check failed!" << std::endl;
                exit(1);
            }
        }
        #endif

    } // for
    
    return mem_total;
}

// convert bit-vector index to vector<usingned>
//
size_t convert_bv2vect(const bv_index& bvi, vect_index& vidx)
{
    size_t  mem_total = 0;
    
    for (bv_index::map_type::const_iterator it = bvi.idx_.begin();
         it != bvi.idx_.end();
         ++it)
    {
        unsigned id = it->first;
        const TBVector* bvp = it->second;
        
        unsigned count = bvp->count();  // population count
        
        vect_index::buffer_type& vect = vidx.idx_[id];
        vect.resize(count);

        for (TBVector::enumerator en = bvp->first(); en.valid(); ++en)
        {
            vect.push_back(*en);
        }
        
        mem_total +=
            sizeof(vect_index::buffer_type::value_type) * vect.size() +
            sizeof(vect_index::buffer_type::size_type);

    } // for
    return mem_total;
}

// convert bit-vector index to bm::sparse_vector
//
size_t convert_bv2sv(const bv_index& bvi, sparse_vect_index& sv_idx)
{
    size_t  mem_total = 0;
    
    std::vector<unsigned> vect;
    
    unsigned sv_pos = 0; // current position in sparse vector
    
    for (bv_index::map_type::const_iterator it = bvi.idx_.begin();
         it != bvi.idx_.end();
         ++it)
    {
        unsigned id = it->first;
        const TBVector* bvp = it->second;
        
        // convert into a plain vector first
        //
        unsigned count = bvp->count();
        vect.resize(count);
        for (TBVector::enumerator en = bvp->first(); en.valid(); ++en)
        {
            vect.push_back(*en);
        }
        
        // convert into delta-vector
        //
        for (unsigned k  = 1; k < count; ++k)
        {
            vect[k] -= vect[k-1];
        }
        
        // merge to sparse vector
        sparse_vect_index::vect_addr vaddr;
        vaddr.offset = sv_pos;
        vaddr.size = (unsigned)vect.size();
        
        sv_idx.sv_storage_.import(vect.data(), vaddr.size, vaddr.offset);
        sv_pos += vaddr.size;
        
        sv_idx.idx_[id] = vaddr;
        
        //std::cout << "offs=" << sv_pos << sv_idx.sv_storage_.size() << " " << std::flush;
        
    } // for
    
    // optimize sparse vector storage, compute memory consumption
    {
        sparse_vect_index::sparse_vector_type::statistics st;
        
        BM_DECLARE_TEMP_BLOCK(tb)
        sv_idx.sv_storage_.optimize(tb, TBVector::opt_compress, &st);
        mem_total += st.memory_used;
    }
    //mem_total += sv_idx.idx_.size() * sizeof(sparse_vect_index::vect_addr);
    //bm::print_svector_stat(sv_idx.sv_storage_, false);

    return mem_total;
}


// speed test for in-memory bit vectors
// benchmark performs a mix of logical operations
//
void speed_test_bv_index(const bv_index& bvi)
{
    TBVector bv_join; // OR join vector
    
    bm::chrono_taker tt1("1. BitVector index operations", 1, &timing_map);

    // join all vectors using OR operation
    for (bv_index::map_type::const_iterator it = bvi.idx_.begin();
         it != bvi.idx_.end();
         ++it)
    {
        const TBVector* bvp = it->second;
        bv_join |= *bvp;
    } // for
    tt1.add_repeats(bvi.idx_.size());
    
    // a group of random vectors from the index map, compute OR
    // then compute AND with the join vector
    //
    TBVector bv_res(bm::BM_GAP);
    std::vector<unsigned> result_set;
    result_set.reserve(result_set_cnt); // memory reservation to avoid reallocs
    
    for (unsigned i = 0; i < benchmark_ops; ++i)
    {
        bv_res.clear(true); // free all blocks
        result_set.resize(0);
        
        for (unsigned j = 0; j < sample_cnt; ++j)
        {
            unsigned id = rand() % index_size;
            bv_index::map_type::const_iterator it = bvi.idx_.find(id);
            if (it == bvi.idx_.end())
                continue;
            const TBVector& bv = *(it->second);
            bv_res |= bv;
        }
        
        bv_res &= bv_join;
        
        // enumerate the final result set, extract first N elements
        //
        TBVector::enumerator en = bv_res.first();
        for (unsigned k = 0; en.valid() && k < result_set_cnt; ++k)
        {
            result_set.push_back(*en);
        }

    } // for i
    tt1.add_repeats(benchmark_ops * sample_cnt + 1);

    std::cout << std::endl;
}

// speed test for in-memory serialized bit vectors
// this function uses bm::operation_deserializer
// to perform logical operation between a BLOB and bvector<> in memory
// and avoids extra decompression overhead
//
void speed_test_bvs_index(const bvs_index& bvs)
{
    TBVector bv_join; // OR join vector
    
    BM_DECLARE_TEMP_BLOCK(tb)

    bm::operation_deserializer<TBVector> des;
    
    bm::chrono_taker tt1("2. BitVector BLOB operations", 1, &timing_map);

    // join all vectors using OR operation
    for (bvs_index::map_type::const_iterator it = bvs.idx_.begin();
         it != bvs.idx_.end();
         ++it)
    {
        const bvs_index::buffer_type& svect = it->second;
        if (svect.size() == 0)
        {
            std::cerr << "Error! Empty buffer detected." << std::endl;
            exit(1);
        }
        const unsigned char* buf = it->second.data();
        
        des.deserialize(bv_join, buf, tb, bm::set_OR);
    } // for
    tt1.add_repeats(bvs.idx_.size());

    // a group of random vectors from the index map, compute OR
    // then compute AND with the join vector
    //
    TBVector bv_res(bm::BM_GAP);
    std::vector<unsigned> result_set;
    result_set.reserve(result_set_cnt); // memory reservation to avoid reallocs
    
    for (unsigned i = 0; i < benchmark_ops; ++i)
    {
        bv_res.clear(true); // free all blocks
        result_set.resize(0);
        
        for (unsigned j = 0; j < sample_cnt; ++j)
        {
            unsigned id = rand() % index_size;
            bvs_index::map_type::const_iterator it = bvs.idx_.find(id);
            if (it == bvs.idx_.end())
                continue;
            
            const unsigned char* buf = it->second.data();
            des.deserialize(bv_res, buf, tb, bm::set_OR);
        } // for j
        
        bv_res &= bv_join;
        
        // enumerate the final result set, extract first N elements
        //
        TBVector::enumerator en = bv_res.first();
        for (unsigned k = 0; en.valid() && k < result_set_cnt; ++k)
        {
            result_set.push_back(*en);
        }

    } // for i
    tt1.add_repeats(benchmark_ops * sample_cnt + 1);


    std::cout << std::endl;
}

void speed_test_vect_index(const vect_index& vecti)
{
    TBVector bv_join; // OR join vector
    
    bm::chrono_taker tt1("3. vector<unsigned> operations", 1, &timing_map);

    // join all vectors using OR operation
    for (vect_index::map_type::const_iterator it = vecti.idx_.begin();
         it != vecti.idx_.end();
         ++it)
    {
        const vect_index::buffer_type& vect = it->second;
        if (vect.size() == 0)
        {
            std::cerr << "Error! Empty vector detected." << std::endl;
            exit(1);
        }
        
        bm::combine_or(bv_join, vect.begin(), vect.end());
    } // for
    tt1.add_repeats(vecti.idx_.size());


    // a group of random vectors from the index map, compute OR
    // then compute AND with the join vector
    //
    TBVector bv_res(bm::BM_GAP);
    std::vector<unsigned> result_set;
    result_set.reserve(result_set_cnt); // memory reservation to avoid reallocs
    
    for (unsigned i = 0; i < benchmark_ops; ++i)
    {
        bv_res.clear(true); // free all blocks
        result_set.resize(0);
        
        for (unsigned j = 0; j < sample_cnt; ++j)
        {
            unsigned id = rand() % index_size;
            vect_index::map_type::const_iterator it = vecti.idx_.find(id);
            if (it == vecti.idx_.end())
                continue;
            
            const vect_index::buffer_type& vect = it->second;
            
            bm::combine_or(bv_join, vect.begin(), vect.end());
        } // for j
        
        bv_res &= bv_join;
        
        // enumerate the final result set, extract first N elements
        //
        TBVector::enumerator en = bv_res.first();
        for (unsigned k = 0; en.valid() && k < result_set_cnt; ++k)
        {
            result_set.push_back(*en);
        }

    } // for i
    tt1.add_repeats(benchmark_ops * sample_cnt + 1);

    std::cout << std::endl;
}



int main(void)
{
    try
    {
        bv_index  bvi; // regular in-memory index id to bvector<>
        bvs_index bvs; // compressed in-memory index id to bvector<> BLOB
        vect_index vecti; // index based on plain uncompressed vector<unsigned>
        sparse_vect_index svi; // all ids in a sparse vector
        
        generate_bv_index(bvi);
        
        size_t bv_mem_total = calc_memory_footprint(bvi);
        size_t bv_mem_total_MB = bv_mem_total / (1024*1024);
        
        std::cout << "bvector index memory footprint = "
                  << bv_mem_total << " (" << bv_mem_total_MB << "MB)"
                  << std::endl;
        
        size_t bvs_mem_total = convert_bv2bvs(bvi, bvs);
        size_t bvs_mem_total_MB = bvs_mem_total / (1024*1024);

        std::cout << "bvector BLOB index memory footprint = "
                  << bvs_mem_total << " (" << bvs_mem_total_MB << "MB)"
                  << std::endl;

        size_t vecti_mem_total = convert_bv2vect(bvi, vecti);
        size_t vecti_mem_total_MB = vecti_mem_total / (1024*1024);

        std::cout << "vector<unsigned> index memory footprint = "
                  << vecti_mem_total << " (" << vecti_mem_total_MB << "MB)"
                  << std::endl;

        size_t svi_mem_total = convert_bv2sv(bvi, svi);
        size_t svi_mem_total_MB = svi_mem_total / (1024*1024);

        std::cout << "bm::sparse_vector index memory footprint = "
                  << svi_mem_total << " (" << svi_mem_total_MB << "MB)"
                  << std::endl;


        speed_test_bv_index(bvi);
        speed_test_bvs_index(bvs);
        speed_test_vect_index(vecti);

        std::cout << std::endl << "Performance (ops/sec):" << std::endl;
        bm::chrono_taker::print_duration_map(timing_map, bm::chrono_taker::ct_ops_per_sec);


        //getchar();

    }
    catch(std::exception& ex)
    {
        std::cerr << ex.what() << std::endl;
        return 1;
    }

    return 0;
}

