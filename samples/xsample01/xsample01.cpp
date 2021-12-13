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

/** \example xsample01.cpp
  Demo and a benchmark on memory consumption control and logical operation
*/
/*! \file xsample01.cpp
    \brief Example: Example: memory consumption techniques
*/


#include <iostream>
#include <memory>
#include <map>
#include <vector>
#include <chrono>
#include <algorithm>
#include <stdexcept>

using namespace std;

#include "bm.h"
#include "bmalgo.h"
#include "bmtimer.h"
#include "bmserial.h"
#include "bmsparsevec.h"

#include "bmsparsevec_algo.h"
#include "bmsparsevec_serial.h"
#include "bmalgo_similarity.h"

#include "bmdbg.h"
#include "bmundef.h" /* clear the pre-proc defines from BM */

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
bm::chrono_taker<>::duration_map_type  timing_map;



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
static
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
    typedef std::map<unsigned, buffer_type>     map_type;
    
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
    
    typedef bm::sparse_vector<unsigned, bm::bvector<> >   sparse_vector_type;
    typedef std::map<unsigned, vect_addr>                 map_type;
    typedef std::vector< std::pair<uint64_t, unsigned> >  delta_sum_map_type;

    
    void get_vector(unsigned id, std::vector<unsigned>& vect) const;
    
    
    sparse_vector_type  sv_storage_;
    sparse_vector_type  sv_storage1_;
    map_type            idx_;
};

void sparse_vect_index::get_vector(unsigned id, std::vector<unsigned>& vect) const
{
    map_type::const_iterator it = idx_.find(id);
    if (it != idx_.end())
    {
        const sparse_vect_index::vect_addr& vaddr = it->second;
        vect.resize(vaddr.size+1);
        vect[0] = sv_storage1_.get(id);
        for (unsigned j = 1; j < vect.size(); ++j)
        {
            unsigned a = sv_storage_.get(j + vaddr.offset - 1);
            a += (vect[j-1] + 1);
            vect[j] = a;
        } // for j
    }
    else
    {
        vect.resize(0);
    }
}

// --------------------------------------------------------------------


// set bits in a vector using various methods picked at random
///
// one method will generate a plato of non-random integers,
// another random integers of near neighbors
// the other adds ints randomly without following any system
//
static
void generate_random_vector(TBVector* bv)
{
    unsigned method = rand() % 5; // pick a generation method
    if (method == 0) // generate a incremental linear sequence at random location
    {
        unsigned seed_id = unsigned(rand()) % max_size;
        for (unsigned i = seed_id; i < seed_id+bits_per_vect; ++i)
        {
            if (i >= max_size)
                break;
            bv->set_bit(i);
        } // for i
    }
    else
    if (method == 1) // generate near neighbors
    {
        unsigned seed_id = unsigned(rand()) % max_size;
        unsigned id = seed_id;
        for (unsigned i = 0; i < bits_per_vect; ++i)
        {
            if (id >= max_size)
                break;
            bv->set_bit(id);
            id += (rand() % 10);
            if (id >= max_size)
                id = unsigned(rand()) % max_size;
        } // for i
    }
    else // generate completely random bits
    {
        for (unsigned i  = 0; i < bits_per_vect; ++i)
        {
            unsigned id = unsigned(rand()) % max_size;
            if (i >= max_size) // paranoiya check
                break;
            bv->set_bit(id);
        } // for i
    }
}

// generate map of bit-vectors, each filled with just a few bits
//
static
void generate_bv_index(bv_index& bvi)
{
    for (unsigned i = 0; i < index_size; ++i)
    {
        std::unique_ptr<TBVector> ap(construct_bvector());
        
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
static
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
static
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
                throw std::runtime_error("deserialization check failed");
            }
        }
        #endif

    } // for
    
    return mem_total;
}

// convert bit-vector index to vector<usingned>
//
static
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

static
void bv2delta(const TBVector& bv, std::vector<unsigned>& vect)
{
    // convert into a plain vector first
    //
    vect.resize(0);
    for (TBVector::enumerator en = bv.first(); en.valid(); ++en)
    {
        vect.push_back(*en);
    }

    // convert into delta-vector
    //
    {
        for (size_t k  = vect.size()-1; k >= 1; --k)
        {
            vect[k] -= vect[k-1];
            --vect[k];
        } // for
    }
}

// convert bit-vector index to bm::sparse_vector
//
static
size_t convert_bv2sv(const bv_index& bvi, sparse_vect_index& sv_idx)
{
    size_t  mem_total = 0;
    
    std::vector<unsigned> vect;
    
    sparse_vect_index::delta_sum_map_type  delta_map;
    
    for (bv_index::map_type::const_iterator it = bvi.idx_.begin();
         it != bvi.idx_.end();
         ++it)
    {
        unsigned id = it->first;
        const TBVector* bvp = it->second;
        
        bv2delta(*bvp, vect);
        
        // compute sum of the delta-vector elements add to the sort map
        {
            uint64_t sum = 0;
            for (unsigned k  = 1; k < vect.size(); ++k)
            {
                sum += vect[k];
            } // for
            delta_map.push_back(std::make_pair(sum, id));
        }
        
    } // for
    
    // sort by "enthropy" (sort of)
    //
    std::sort(delta_map.begin(), delta_map.end());
    if (delta_map.size() != bvi.idx_.size()) // paranoia check
    {
        throw std::runtime_error("delta map size is incorrect");
    }
    
    unsigned sv_pos = 0; // current position in sparse vector
    for (unsigned j = 0; j < delta_map.size(); ++j)
    {
        unsigned id = delta_map[j].second;
        
        bv_index::map_type::const_iterator it = bvi.idx_.find(id);
        if (it == bvi.idx_.end())
            continue;
        const TBVector& bv = *(it->second);
        
        // convert into a plain delta vector again
        bv2delta(bv, vect);

        sparse_vect_index::vect_addr vaddr;
        vaddr.offset = sv_pos;
        vaddr.size = (unsigned)(vect.size() - 1);
        
        sv_idx.sv_storage1_.set(id, vect[0]);
        
        if (vaddr.size)
        {
            sv_idx.sv_storage_.import(&vect[1], vaddr.size, vaddr.offset);
            sv_pos += vaddr.size;
        }
        
        sv_idx.idx_[id] = vaddr;
    } // for

    
    // optimize sparse vector storage, compute memory consumption
    {
        sparse_vect_index::sparse_vector_type::statistics st;
        
        BM_DECLARE_TEMP_BLOCK(tb)
        sv_idx.sv_storage_.optimize(tb, TBVector::opt_compress, &st);
        mem_total += st.memory_used;
        sv_idx.sv_storage1_.optimize(tb, TBVector::opt_compress, &st);
        mem_total += st.memory_used;
    }
    
    // check
    for (bv_index::map_type::const_iterator it = bvi.idx_.begin();
         it != bvi.idx_.end();
         ++it)
    {
        unsigned id = it->first;
        const TBVector* bvp = it->second;
        
        // convert into a plain vector first
        //
        vect.resize(0);
        for (TBVector::enumerator en = bvp->first(); en.valid(); ++en)
        {
            vect.push_back(*en);
        }
        
        
        std::vector<unsigned> svect;
        sv_idx.get_vector(id, svect);
        if (svect.size() != vect.size())
        {
            std::cerr << "Size check failed! id = " << id
                      << "size() = " << svect.size()
                      << std::endl;
            throw std::runtime_error("sparse vector content check failed");
        }
        
        for (unsigned k = 0; k < vect.size(); ++k)
        {
            if (vect[k] != svect[k])
            {
                std::cerr << "SV content check failed! id = " << id
                          <<  " i=" << k << std::endl;
                for (unsigned h = 0; h < vect.size(); ++h)
                {
                    std::cout << "[" << vect[h] << "=" << svect[h] << "], ";
                } // for h
                std::cout << std::endl;
                throw std::runtime_error("sparse vector content check failed");
            }
        } // for k

    } // for
    
    #ifdef DEBUG
    bm::print_svector_stat(sv_idx.sv_storage_, true);
    bm::print_svector_stat(sv_idx.sv_storage1_, true);
    #endif

    return mem_total;
}


// speed test for in-memory bit vectors
// benchmark performs a mix of logical operations
//
static
void speed_test_bv_index(const bv_index& bvi)
{
    TBVector bv_join; // OR join vector
    
    bm::chrono_taker tt1(cout, "1. bm::bvector<> index", 1, &timing_map);

    // join all vectors using OR operation
    for (bv_index::map_type::const_iterator it = bvi.idx_.begin();
         it != bvi.idx_.end();
         ++it)
    {
        const TBVector* bvp = it->second;
        bv_join |= *bvp;
    } // for
    bv_join.optimize();
    
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
            unsigned id = unsigned(rand()) % index_size;
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
        for (unsigned k = 0; en.valid() && k < result_set_cnt; ++k, ++en)
        {
            result_set.push_back(*en);
        }

    } // for i
    
    tt1.add_repeats(benchmark_ops + 1);
}

// speed test for in-memory serialized bit vectors
// this function uses bm::operation_deserializer
// to perform logical operation between a BLOB and bvector<> in memory
// and avoids extra decompression overhead
//
static
void speed_test_bvs_index(const bvs_index& bvs)
{
    TBVector bv_join; // OR join vector
    
    BM_DECLARE_TEMP_BLOCK(tb)

    bm::operation_deserializer<TBVector> des;
    
    bm::chrono_taker tt1(cout, "2. serialized bvector", 1, &timing_map);

    // join all vectors using OR operation
    for (bvs_index::map_type::const_iterator it = bvs.idx_.begin();
         it != bvs.idx_.end();
         ++it)
    {
        const bvs_index::buffer_type& svect = it->second;
        if (svect.size() == 0)
        {
           throw std::runtime_error("empty buffer error");
        }
        const unsigned char* buf = it->second.data();
        
        des.deserialize(bv_join, buf, tb, bm::set_OR);
    } // for
    bv_join.optimize();

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
            unsigned id = unsigned(rand()) % index_size;
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
        for (unsigned k = 0; en.valid() && k < result_set_cnt; ++k, ++en)
        {
            result_set.push_back(*en);
        }
    } // for i
    
    tt1.add_repeats(benchmark_ops + 1);
}

static
void speed_test_vect_index(const vect_index& vecti)
{
    TBVector bv_join; // OR join vector
    
    bm::chrono_taker tt1(cout, "3. std::vector<unsigned> ", 1, &timing_map);

    // join all vectors using OR operation
    for (vect_index::map_type::const_iterator it = vecti.idx_.begin();
         it != vecti.idx_.end();
         ++it)
    {
        const vect_index::buffer_type& vect = it->second;
        if (vect.size() == 0)
        {
            throw std::runtime_error("empty buffer error");
        }
        
        bm::combine_or(bv_join, vect.begin(), vect.end());
    } // for
    bv_join.optimize();


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
            unsigned id = unsigned(rand()) % index_size;
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
        for (unsigned k = 0; en.valid() && k < result_set_cnt; ++k, ++en)
        {
            result_set.push_back(*en);
        }

    } // for i
    
    tt1.add_repeats(benchmark_ops + 1);
}

static
void speed_test_sv_index(const sparse_vect_index& svi)
{
    TBVector bv_join; // OR join vector
    
    bm::chrono_taker tt1(cout, "4. bm::sparse_vector<unsigned> ", 1, &timing_map);
    
    std::vector<unsigned> vect;

    // join all vectors using OR operation
    for (sparse_vect_index::map_type::const_iterator it = svi.idx_.begin();
         it != svi.idx_.end();
         ++it)
    {
        unsigned id = it->first;
        svi.get_vector(id, vect);
        
        bm::combine_or(bv_join, vect.begin(), vect.end());
    } // for
    bv_join.optimize();


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
            unsigned id = unsigned(rand()) % index_size;
            svi.get_vector(id, vect);
            if (vect.size() == 0)
                continue;
            
            bm::combine_or(bv_join, vect.begin(), vect.end());
        } // for j
        
        bv_res &= bv_join;
        
        // enumerate the final result set, extract first N elements
        //
        TBVector::enumerator en = bv_res.first();
        for (unsigned k = 0; en.valid() && k < result_set_cnt; ++k, ++en)
        {
            result_set.push_back(*en);
        }

    } // for i

    tt1.add_repeats(benchmark_ops + 1);
}



int main(void)
{
    try
    {
        bv_index  bvi; // regular in-memory index id to bvector<>
        bvs_index bvs; // compressed in-memory index id to bvector<> BLOB
        vect_index vecti; // index based on plain uncompressed vector<unsigned>
        sparse_vect_index svi; // all ids in a sparse vector
        
        // experiments generation, measuring memory footprints
        //
        generate_bv_index(bvi);
        
        size_t bv_mem_total = calc_memory_footprint(bvi);
        size_t bv_mem_total_MB = bv_mem_total / (1024*1024);
        
        std::cout << "bm::bvector<> memory footprint = "
                  << bv_mem_total << " (" << bv_mem_total_MB << "MB)"
                  << std::endl;
        
        size_t bvs_mem_total = convert_bv2bvs(bvi, bvs);
        size_t bvs_mem_total_MB = bvs_mem_total / (1024*1024);

        std::cout << "bm::bvector<> BLOB memory footprint = "
                  << bvs_mem_total << " (" << bvs_mem_total_MB << "MB)"
                  << std::endl;

        size_t vecti_mem_total = convert_bv2vect(bvi, vecti);
        size_t vecti_mem_total_MB = vecti_mem_total / (1024*1024);

        std::cout << "std::vector<unsigned> memory footprint = "
                  << vecti_mem_total << " (" << vecti_mem_total_MB << "MB)"
                  << std::endl;

        size_t svi_mem_total = convert_bv2sv(bvi, svi);
        size_t svi_mem_total_MB = svi_mem_total / (1024*1024);

        std::cout << "bm::sparse_vector<> memory footprint = "
                  << svi_mem_total << " (" << svi_mem_total_MB << "MB)"
                  << std::endl;

        // run performance tests
        //

        speed_test_bv_index(bvi);
        speed_test_bvs_index(bvs);
        speed_test_vect_index(vecti);
        speed_test_sv_index(svi);


        std::cout << std::endl << "Performance (ops/sec):" << std::endl;
        bm::chrono_taker<>::print_duration_map(cout, timing_map, bm::chrono_taker<>::ct_ops_per_sec);

        //getchar();  // uncomment to check memory consumption at the OS level

    }
    catch(std::exception& ex)
    {
        std::cerr << ex.what() << std::endl;
        return 1;
    }

    return 0;
}

