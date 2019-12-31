/*
Copyright(c) 2002-2018 Anatoliy Kuznetsov(anatoliy_kuznetsov at yahoo.com)

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


//#define BMSSE2OPT
//#define BMSSE42OPT
//#define BMAVX2OPT
//#define BM_USE_EXPLICIT_TEMP
//#define BM_USE_GCC_BUILD

#define BMXORCOMP

#include <stdio.h>
#include <stdlib.h>
#undef NDEBUG
#include <cassert>
#include <time.h>
#include <math.h>
#include <string.h>

#include <iostream>
#include <iomanip>
#include <utility>
#include <memory>
#include <random>
#include <algorithm>
#include <stdarg.h>  

#include <bm.h>
#include <bmalgo.h>
#include <bmxor.h>
#include <bmaggregator.h>
#include <bmutil.h>
#include <bmserial.h>
#include <bmrandom.h>
#include <bmvmin.h>
#include <bmbmatrix.h>
#include <bmsparsevec.h>
#include <bmsparsevec_algo.h>
#include <bmsparsevec_serial.h>
#include <bmalgo_similarity.h>
#include <bmsparsevec_util.h>
#include <bmsparsevec_compr.h>
#include <bmstrsparsevec.h>
#include <bmtimer.h>

using namespace bm;
using namespace std;

#include "rlebtv.h"
#include <encoding.h>
#include <limits.h>

#include <bmdbg.h>

#include <vector>


#define POOL_SIZE 5000

//#define MEM_POOL


template<class T> T* pool_allocate(T** pool, int& i, size_t n)
{
    return i ? pool[i--] : (T*) ::malloc(n * sizeof(T));
}

inline void* pool_allocate2(void** pool, int& i, size_t n)
{
    return i ? pool[i--] : malloc(n * sizeof(void*));
}



template<class T> void pool_free(T** pool, int& i, T* p)
{
    i < POOL_SIZE ? (free(p),(void*)0) : pool[++i]=p;
}


class pool_block_allocator
{
public:

    static bm::word_t* allocate(size_t n, const void *)
    {
        int *idx = 0;
        bm::word_t** pool = 0;

        switch (n)
        {
        case bm::set_block_size:
            idx = &bit_blocks_idx_;
            pool = free_bit_blocks_;
            break;

        case 64:
            idx = &gap_blocks_idx0_;
            pool = gap_bit_blocks0_;
            break;

        case 128:
            idx = &gap_blocks_idx1_;
            pool = gap_bit_blocks1_;
            break;
        
        case 256:
            idx = &gap_blocks_idx2_;
            pool = gap_bit_blocks2_;
            break;

        case 512:
            idx = &gap_blocks_idx3_;
            pool = gap_bit_blocks3_;
            break;

        default:
            assert(0);
        }

        return pool_allocate(pool, *idx, n);
    }

    static void deallocate(bm::word_t* p, size_t n)
    {
        int *idx = 0;
        bm::word_t** pool = 0;

        switch (n)
        {
        case bm::set_block_size:
            idx = &bit_blocks_idx_;
            pool = free_bit_blocks_;
            break;

        case 64:
            idx = &gap_blocks_idx0_;
            pool = gap_bit_blocks0_;
            break;

        case 128:
            idx = &gap_blocks_idx1_;
            pool = gap_bit_blocks1_;
            break;
        
        case 256:
            idx = &gap_blocks_idx2_;
            pool = gap_bit_blocks2_;
            break;

        case 512:
            idx = &gap_blocks_idx3_;
            pool = gap_bit_blocks3_;
            break;

        default:
            assert(0);
        }

        pool_free(pool, *idx, p);
    }

private:
    static bm::word_t* free_bit_blocks_[];
    static int         bit_blocks_idx_;

    static bm::word_t* gap_bit_blocks0_[];
    static int         gap_blocks_idx0_;

    static bm::word_t* gap_bit_blocks1_[];
    static int         gap_blocks_idx1_;

    static bm::word_t* gap_bit_blocks2_[];
    static int         gap_blocks_idx2_;

    static bm::word_t* gap_bit_blocks3_[];
    static int         gap_blocks_idx3_;
};

bm::word_t* pool_block_allocator::free_bit_blocks_[POOL_SIZE];
int pool_block_allocator::bit_blocks_idx_ = 0;

bm::word_t* pool_block_allocator::gap_bit_blocks0_[POOL_SIZE];
int pool_block_allocator::gap_blocks_idx0_ = 0;

bm::word_t* pool_block_allocator::gap_bit_blocks1_[POOL_SIZE];
int pool_block_allocator::gap_blocks_idx1_ = 0;

bm::word_t* pool_block_allocator::gap_bit_blocks2_[POOL_SIZE];
int pool_block_allocator::gap_blocks_idx2_ = 0;

bm::word_t* pool_block_allocator::gap_bit_blocks3_[POOL_SIZE];
int pool_block_allocator::gap_blocks_idx3_ = 0;




class pool_ptr_allocator
{
public:
    static void* allocate(size_t n, const void *)
    {
        return pool_allocate2(free_ptr_blocks_, ptr_blocks_idx_, n);
    }
    static void deallocate(void* p, size_t)
    {
        pool_free(free_ptr_blocks_, ptr_blocks_idx_, p);
    }
private:
    static void*  free_ptr_blocks_[];
    static int    ptr_blocks_idx_;
};

void* pool_ptr_allocator::free_ptr_blocks_[POOL_SIZE];
int pool_ptr_allocator::ptr_blocks_idx_ = 0;

#if defined(BMSSE2OPT) || defined(BMSSE42OPT) || defined(BMAVX2OPT) || defined(BMAVX512OPT)
#else
# define MEM_DEBUG
#endif
 
#ifdef MEM_DEBUG


class dbg_block_allocator
{
public:
static size_t na_;
static size_t nf_;

    static bm::word_t* allocate(size_t n, const void *)
    {
        ++na_;
        assert(n);
        bm::word_t* p =
            (bm::word_t*) ::malloc((n+1) * sizeof(bm::word_t));
        if (!p)
        {
            std::cerr << "ERROR Failed allocation!" << endl;
            assert(0);exit(1);
        }
        *p = (bm::word_t)n;
        return ++p;
    }

    static void deallocate(bm::word_t* p, size_t n)
    {
        ++nf_;
        --p;
        if (*p != n)
        {
            printf("Block memory deallocation ERROR! n = %i (expected %i)\n", (int)n, (int)*p);
            assert(0);exit(1);
        }
        ::free(p);
    }

    static size_t balance()
    {
        return nf_ - na_;
    }
};

size_t dbg_block_allocator::na_ = 0;
size_t dbg_block_allocator::nf_ = 0;

class dbg_ptr_allocator
{
public:
static size_t na_;
static size_t nf_;

    static void* allocate(size_t n, const void *)
    {
        ++na_;
        assert(sizeof(size_t) == sizeof(void*));
        void* p = ::malloc((n+1) * sizeof(void*));
        if (!p)
        {
            std::cerr << "ERROR! Failed allocation!" << endl;
            exit(1);
        }
        size_t* s = (size_t*) p;
        *s = n;
        return (void*)++s;
    }

    static void deallocate(void* p, size_t n)
    {
        ++nf_;
        size_t* s = (size_t*) p;
        --s;
        if(*s != n)
        {
            printf("Ptr memory deallocation ERROR!\n");
            assert(0);
            exit(1);
        }
        ::free(s);
    }

    static size_t balance()
    {
        return nf_ - na_;
    }

};

size_t dbg_ptr_allocator::na_ = 0;
size_t dbg_ptr_allocator::nf_ = 0;


typedef mem_alloc<dbg_block_allocator, dbg_ptr_allocator, alloc_pool<dbg_block_allocator, dbg_ptr_allocator> > dbg_alloc;

typedef bm::bvector<dbg_alloc> bvect;
typedef bm::bvector_mini<dbg_block_allocator> bvect_mini;
typedef bm::rs_index<dbg_alloc> rs_ind;

#else

#ifdef MEM_POOL

typedef mem_alloc<pool_block_allocator, pool_ptr_allocator> pool_alloc;
typedef bm::bvector<pool_alloc> bvect;
typedef bm::bvector_mini<bm::block_allocator> bvect_mini;
typedef bm::rs_index<pool_block_allocator> rs_ind;


#else

typedef bm::bvector<> bvect;
typedef bm::bvector_mini<bm::block_allocator> bvect_mini;
typedef bm::rs_index<> rs_ind;

#endif

#endif

typedef bm::sparse_vector<unsigned, bvect > sparse_vector_u32;
typedef bm::sparse_vector<unsigned long long, bvect > sparse_vector_u64;
typedef bm::rsc_sparse_vector<unsigned, sparse_vector_u32> rsc_sparse_vector_u32;

//const unsigned BITVECT_SIZE = 100000000 * 8;

// This this setting program will consume around 150M of RAM
const unsigned BITVECT_SIZE = 100000000 * 2;

const unsigned ITERATIONS = 180000;
//const unsigned PROGRESS_PRINT = 2000000;



template<class BV>
void DetailedCompareBVectors(const BV& bv1, const BV& bv2)
{
    bvect::counted_enumerator en1 = bv1.first();
    bvect::counted_enumerator en2 = bv2.first();
    
    for (; en1.valid(); ++en1)
    {
        assert(en2.valid());
        
        bm::id_t i1 = *en1;
        bm::id_t i2 = *en2;
        
        if (i1 != i2)
        {
            unsigned nb1 = unsigned(i1 >>  bm::set_block_shift);
            unsigned nb2 = unsigned(i2 >>  bm::set_block_shift);
            unsigned ii1 = nb1 >> bm::set_array_shift;
            unsigned jj1 = nb1 &  bm::set_array_mask;
            unsigned ii2 = nb2 >> bm::set_array_shift;
            unsigned jj2 = nb2 &  bm::set_array_mask;

            std::cerr << "Difference detected at: position="
                      << i1 << " nb=" << nb1
                      << "[" << ii1 << ", " << jj1 << "]"
                      " other position = " << i2 <<
                      " nb=" << nb2
                      << "[" << ii2 << ", " << jj2 << "]"
                      << std::endl;
            std::cerr << " en1.count()=" << en1.count() << " en2.count()=" << en2.count()
                      << std::endl;
            
            exit(1);
            return;
        }
        ++en2;
    } // for

    bool eq = bv1.equal(bv2);
    if (!eq)
    {
        cerr << "EQ (1-2) discrepancy! " << endl;
        exit(1);
    }
    int cmp = bv1.compare(bv2);
    if (cmp != 0)
    {
        cerr << "Compare (1-2) discrepancy! " << cmp << endl;
        exit(1);
    }
    cmp = bv2.compare(bv1);
    if (cmp != 0)
    {
        cerr << "Compare (2-1) discrepancy! " << cmp << endl;
        exit(1);
    }

    cout << "Detailed compare OK (no difference)." << endl;
}

void CheckVectors(bvect_mini &bvect_min, 
                  bvect      &bvect_full,
                  unsigned size,
                  bool     detailed = true);

void generate_bvector(bvect& bv, unsigned vector_max = 40000000, bool optimize=true);

static
unsigned random_minmax(unsigned min, unsigned max)
{
    unsigned r = (unsigned(rand()) << 16u) | unsigned(rand());
    return r % (max-min) + min;
}

template<typename BV>
void TestFindDiff(const BV& bv1, BV& bv2)
{
    bool f;
    typename BV::size_type pos, pos_c, pos_l;
    f = bv1.find_first_mismatch(bv2, pos);
    bvect bv_x;
    bv_x.bit_xor(bv1, bv2, bvect::opt_compress);
    if (!f)
    {
        auto a = bv_x.any();
        assert(!a);
        return;
    }
    else // found
    {
        bool f2 = bv1.find_first_mismatch(bv2, pos_l, pos);
        assert(f2 == f);
        assert(pos_l == pos);
        if (pos)
        {
            f2 = bv1.find_first_mismatch(bv2, pos_l, pos-1);
            assert(!f2);
        }
    }
    bool cf = bv_x.find(pos_c);
    assert(f == cf);
    assert(pos == pos_c);

    f = bv2.find_first_mismatch(bv1, pos);
    assert(f == cf);
    assert(pos == pos_c);
}


static
void FillSets(bvect_mini* bvect_min, 
              bvect* bvect_full,
              unsigned min, 
              unsigned max,
              unsigned fill_factor)
{
    unsigned i;
    unsigned id;

    //Random filling
    if(fill_factor == 0)
    {
        unsigned n_id = (max - min) / 100;
        printf("random filling : %i\n", n_id);
        for (i = 0; i < n_id; i++)
        {
            id = random_minmax(min, max);
            bvect_min->set_bit(id);
            bvect_full->set_bit(id);
        }
        cout << endl;
    }
    else
    {
        printf("fill_factor random filling : factor = %i\n", fill_factor);

        for(i = 0; i < fill_factor; i++)
        {
            unsigned k = unsigned(rand()) % 10;
            if (k == 0)
                k+=2;

            //Calculate start
            unsigned start = min + (max - min) / (fill_factor * k);

            //Randomize start
            start += random_minmax(1, (max - min) / (fill_factor * 10));

            if (start > max)
            {
                start = min;
            }
            
            //Calculate end 
            unsigned end = start + (max - start) / (fill_factor *2);

            //Randomize end
            end -= random_minmax(1, (max - start) / (fill_factor * 10));

            if (end > max )
            {
                end = max;
            }

            bvect::bulk_insert_iterator iit = bvect_full->inserter();

            if (fill_factor > 1)
            {
                for(; start < end;)
                {
                    unsigned r = unsigned(rand()) % 8;

                    if (r > 7)
                    {
                        unsigned inc = unsigned(rand()) % 3;
                        ++inc;
                        unsigned end2 = start + rand() % 1000;
                        if (end2 > end)
                            end2 = end;
                        while (start < end2)
                        {
                            bvect_min->set_bit(start);
                            iit = start;
                            start += inc;
                        }
                        continue;
                    }

                    if (r)
                    {
                        bvect_min->set_bit(start);
                        iit = start;
                        ++start;
                    }
                    else
                    {
                        start+=r;
                        bvect_min->set_bit(start);
                        iit = start;
                    }
                }
            }
            else
            {
                unsigned c = unsigned(rand()) % 15;
                if (c == 0)
                    ++c;
                for(; start < end; ++start)
                {
                    bvect_min->set_bit(start);
                    iit = start;
                    if (start % c)
                    {
                        start += c;
                    }
                }
            }
            cout << endl;

        }
    }
}

//
// Interval filling.
// 111........111111........111111..........11111111.......1111111...
//

static
void FillSetsIntervals(bvect_mini* bvect_min, 
              bvect& bvect_full,
              unsigned min, 
              unsigned max,
              unsigned fill_factor,
              bool set_flag=true)
{

    while(fill_factor==0)
    {
        fill_factor=rand()%10;
    }
    bvect_full.init();

    cout << "Intervals filling. Factor=" 
         <<  fill_factor << endl << endl;

    unsigned i, j;
    unsigned factor = 70 * fill_factor;
    for (i = min; i < max; ++i)
    {
        unsigned len, end; 

        do
        {
            len = unsigned(rand()) % factor;
            end = i+len;
            
        } while (end >= max);
        if (i < end)
        {
            bvect_full.set_range(i, end-1, set_flag);
        }
       
        for (j = i; j < end; ++j)
        {
            if (set_flag)
            {
                if (bvect_min)
                    bvect_min->set_bit(j);
                //bvect_full.set_bit(j);
            }
            else
            {
                if (bvect_min)
                    bvect_min->clear_bit(j);
                //bvect_full.clear_bit(j);
            }

                           
        } // j


        i = end;


        len = unsigned(rand()) % (factor* 10 * bm::gap_max_bits);
        if (len % 2)
        {
            len *= unsigned(rand()) % (factor * 10);
        }

        i+=len;

        if ( (len % 6) == 0)  
        {
            for(unsigned k=0; k < 1000 && i < max; k+=3,i+=3)
            {
                if (set_flag)
                {
                    if (bvect_min)
                        bvect_min->set_bit(i);
                    bvect_full.set_bit_no_check(i);
                }
                else
                {
                    if (bvect_min)
                        bvect_min->clear_bit(j);
                    bvect_full.clear_bit(j);
                }
            }
        }
    } // for i

}

static
void FillSetClearIntervals(bvect_mini* bvect_min, 
              bvect* bvect_full,
              unsigned min, 
              unsigned max,
              unsigned fill_factor)
{
    FillSetsIntervals(bvect_min, *bvect_full, min, max, fill_factor, true);
    FillSetsIntervals(bvect_min, *bvect_full, min, max, fill_factor, false);
}

static
void FillSetsRandomOne(bvect_mini* bvect_min, 
                       bvect* bvect_full,
                       unsigned min, 
                       unsigned max)
{
    unsigned range = max - min;
    unsigned bit_idx = unsigned(rand()) % range;
    bvect_min->set_bit(bit_idx);
    bvect_full->set_bit(bit_idx);
    cout << "Bit_idx=" << bit_idx << endl;
}

static
void FillSetsRandom(bvect_mini* bvect_min, 
              bvect* bvect_full,
              unsigned min, 
              unsigned max,
              unsigned fill_factor)
{
    bvect_full->init();
    unsigned diap = max - min;
    unsigned count;

    switch (fill_factor)
    {
    case 0:
        count = diap / 1000;
        break;
    case 1:
        count = diap / 255;
        break;
    default:
        count = diap / 10;
        break;
    }

    for (unsigned i = 0; i < count; ++i)
    {
        unsigned bn = unsigned(rand()) % count;
        bn += min;

        if (bn > max)
        {
            bn = max;
        }
        bvect_min->set_bit(bn);
        bvect_full->set_bit_no_check(bn);
    }
    cout << "Ok" << endl;

}

static
void FillSetsRegular(bvect_mini* bvect_min,
                     bvect* bvect_full,
              unsigned /*min*/,
              unsigned max,
              unsigned /*fill_factor*/)
{
    bvect::bulk_insert_iterator iit = bvect_full->inserter();

    unsigned step = rand() % 4;
    if (step < 2) ++step;
    for (unsigned i = 0; i < max; i+=step)
    {
        bvect_min->set_bit(i);
        iit = i;
        //bvect_full->set_bit_no_check(i);
    }
    cout << "Ok" << endl;
}




//
//  Quasi random filling with choosing randomizing method.
//
//
static
void FillSetsRandomMethod(bvect_mini* bvect_min, 
                          bvect* bvect_full,
                          unsigned min, 
                          unsigned max,
                          int optimize = 0,
                          int method = -1)
{
    if (method == -1)
    {
        method = rand() % 7;
    }
    unsigned factor;
///method = 3;
    switch (method)
    {

    case 0:
        cout << "Random filling: method - FillSets - factor(0)" << endl;
        FillSets(bvect_min, bvect_full, min, max, 0);
        break;

    case 1:
        cout << "Random filling: method - FillSets - factor(random)" << endl;
        factor = rand()%3;
        FillSets(bvect_min, bvect_full, min, max, factor?factor:1);
        break;

    case 2:
        cout << "Random filling: method - Set-Clear Intervals - factor(random)" << endl;
        factor = rand()%10;
        FillSetClearIntervals(bvect_min, bvect_full, min, max, factor);
        break;
    case 3:
        cout << "Random filling: method - FillRandom - factor(random)" << endl;
        factor = rand()%3;
        FillSetsRandom(bvect_min, bvect_full, min, max, factor?factor:1);
        break;
    case 4:
        cout << "Random set one bit" << endl;
        FillSetsRandomOne(bvect_min, bvect_full, min, max);
        break;
    case 5:
        cout << "Regular pattern filling" << endl;
        FillSetsRegular(bvect_min, bvect_full, min, max, 2);
        break;
    default:
        cout << "Random filling: method - Set Intervals - factor(random)" << endl;
        factor = rand()%10;
        FillSetsIntervals(bvect_min, *bvect_full, min, max, factor);
        break;

    } // switch

    if (optimize && (method <= 1))
    {
        cout << "Vector optimization..." << flush;
        BM_DECLARE_TEMP_BLOCK(tb)
        bvect_full->optimize(tb);
        cout << "OK" << endl;
    }
}

static
void print_bv(const bvect& bv)
{
    std::cout << bv.count() << ": ";
    bvect::enumerator en = bv.first();
    for (; en.valid(); ++en)
    {
        std::cout << *en << ", ";
    }
    std::cout << std::endl;
}

// reference SHIFT right
static
void ShiftRight(bvect*  bv, unsigned shift)
{
    bvect bv_tmp;
    {
        bvect::bulk_insert_iterator bi = bv_tmp.inserter();
        bvect::enumerator en = bv->first();
        for (; en.valid(); ++en)
        {
            unsigned v = *en;
            unsigned new_v = v + shift;
            if (new_v < v || new_v == bm::id_max) // check overflow
            {}
            else
            {
                bi = new_v;
            }
        }
    }
    bv->swap(bv_tmp);
}

// Reference bit insert
static
void BVectorInsert(bvect*  bv, unsigned pos, bool value)
{
    bvect bv_tmp;
    if (pos)
        bv_tmp.copy_range(*bv, 0, pos-1);
    
    {
        bvect::bulk_insert_iterator bi = bv_tmp.inserter();
        bvect::enumerator en = bv->first();
        for (; en.valid(); ++en)
        {
            unsigned v = *en;
            if (v < pos)
                continue;
            unsigned new_v = v + 1;
            if (new_v < v || new_v == bm::id_max) // check overflow
            {}
            else
            {
                bi = new_v;
            }
        }
    }
    bv->swap(bv_tmp);
    bv->set(pos, value);
}

// Reference bit erase
static
void BVectorErase(bvect*  bv, unsigned pos)
{
    bvect bv_tmp;
    if (pos)
        bv_tmp.copy_range(*bv, 0, pos-1);
    
    {
        bvect::bulk_insert_iterator bi = bv_tmp.inserter();
        bvect::enumerator en = bv->first();
        en.go_to(pos+1);
        for (; en.valid(); ++en)
        {
            unsigned v = *en;
            assert(v > pos);
            bi = v -1;
        }
    }
    bv->swap(bv_tmp);
}



// do logical operation through serialization
static
unsigned SerializationOperation(bvect*             bv_target,
                                /*const*/ bvect&   bv1,
                                /*const*/ bvect&   bv2,
                                set_operation      op,
                                bool               check_reverse=false)
{
    bvect bv_tmp;
    if (!bv_target)
    {
        bv_target = &bv_tmp;
    }

    if (op == set_COUNT_SUB_AB ||
        op == set_COUNT_SUB_BA)
    {
        check_reverse = false;
    }

    // serialize input vectors
    bvect::statistics *st1_op, *st2_op;
    st1_op = new bvect::statistics;
    st2_op = new bvect::statistics;

    BM_DECLARE_TEMP_BLOCK(tb)
    bv1.optimize(tb, bvect::opt_compress, st1_op);
    bv2.optimize(tb, bvect::opt_compress, st2_op);


   struct bvect::statistics st1, st2;
   bv1.calc_stat(&st1);
   bv2.calc_stat(&st2);


   if (st1.max_serialize_mem > st1_op->max_serialize_mem)
   {
       cout << "Error: Optimize failed to compute max_serialize_mem" << endl;
       cout << "calc_stat=" << st1.max_serialize_mem << endl;
       cout << "optimize=" << st1_op->max_serialize_mem << endl;
       assert(0);
       exit(1);
   }
   if (st2.max_serialize_mem > st2_op->max_serialize_mem)
   {
       cout << "Error:Optimize failed to compute max_serialize_mem" << endl;
       cout << "calc_stat=" << st2.max_serialize_mem << endl;
       cout << "optimize=" << st2_op->max_serialize_mem << endl;
       assert(0); exit(1);
   }

   delete st1_op;
   delete st2_op;

   unsigned char* smem1 = new unsigned char[st1.max_serialize_mem];
   unsigned char* smem2 = new unsigned char[st2.max_serialize_mem];

   size_t slen1 = bm::serialize(bv1, smem1, tb);
   size_t slen2 = bm::serialize(bv2, smem2, tb);

   if (slen1 > st1.max_serialize_mem || slen2 > st2.max_serialize_mem)
   {
       cout << "Serialization override detected!" << endl;
       exit(1);
   }

   operation_deserializer<bvect> od;

   auto count = od.deserialize(*bv_target, smem1, nullptr, set_ASSIGN);
   cout << slen1 << " " << slen2 << endl;
   int res = bv1.compare(*bv_target);
   if (res != 0)
   {
       cout << "---------------------------------- " << endl;
       cout << "bv1.count()=" << bv1.count() << endl;
       print_stat(bv1);
       cout << "---------------------------------- " << endl;
       cout << "bv_target.count()=" << bv_target->count() << endl;
       print_stat(*bv_target);
       
       bv_target->bit_xor(bv1);
       cout << "First diff=" << bv_target->get_first() << endl;
       cout << "set_ASSIGN 1 failed!" << endl;
       exit (1);
   }
   cout << "Deserialization ASSIGN into bv1 OK" << endl;

   {
       bvect* bv_tmp2 = new bvect();
       bm::deserialize(*bv_tmp2, smem1);
       if (*bv_tmp2 != bv1)
       {
           cout << "Deserialize NOT equal to Operation deserialize!" << endl;
           assert(0); exit(1);
       }
       delete bv_tmp2;
   }


   cout << "Operation deserialization... " << op << endl;

    count=
       od.deserialize(*bv_target,
                                                  smem2,
                                                  0,
                                                  op);
    cout << "OK" << endl;

    // check if operation was ok
    {
        bvect bv_agg, bv_agg2;
        
        bm::aggregator<bvect> agg;
        agg.set_optimization();
        
        bvect* agg_list[10];
        bvect* agg_list2[10];
        agg_list[0] = &bv1;
        agg_list[1] = &bv2;
        agg_list2[0] = &bv2;


        bool agg_check = false;

        bvect bvt(bv1);
        switch(op)
        {
        case bm::set_OR:
            {
                bvect bvc(bv1);
                bvc |= bv2;
                bvect bv_merge1(bv1);
                bvect bv_merge2(bv2);
                bv_merge1.merge(bv_merge2);
                
                if (bv_merge1 != bvc)
                {
                    cerr << "Merge(OR) check error!" << endl;
                    exit(1);
                }
                // 2-way
                {
                    bvect bvt1;
                    bvt1.bit_or(bv1, bv2, bvect::opt_none);
                    if (bvt1 != bvc)
                    {
                        cerr << "1. OR 2-way check error!" << endl;
                        exit(1);
                    }
                    bvect bvt2;
                    bvt2.bit_or(bv2, bv1, bvect::opt_compress);
                    if (bvt2 != bvc)
                    {
                        cerr << "2. OR 2-way check error!" << endl;
                        exit(1);
                    }
                }
            }
            bvt |= bv2;
            agg.combine_or(bv_agg, agg_list, 2);
            agg_check = true;
            break;
        case bm::set_XOR:
            bvt ^= bv2;
            // 2-way
            {
                bvect bvc(bv1);
                bvc ^= bv2;
                
                bvect bvt1;
                bvt1.bit_xor(bv1, bv2, bvect::opt_none);
                if (bvt1 != bvc)
                {
                    cerr << "1. XOR 2-way check error!" << endl;
                    cerr << "Run detailed check (1)..." << endl;
                    DetailedCompareBVectors(bvt1, bvc);

                    bvect bvt2;
                    bvt2.bit_xor(bv2, bv1, bvect::opt_compress);
                    if (bvt2 != bvc)
                    {
                        cerr << "Reverse XOR 2-way check failed too!" << endl;
                    }

                    cerr << "Run detailed check (2)..." << endl;
                    DetailedCompareBVectors(bvt2, bvc);

                    exit(1);
                }
                bvect bvt2;
                bvt2.bit_xor(bv2, bv1, bvect::opt_compress);
                if (bvt2 != bvc)
                {
                    cerr << "2. XOR 2-way check error!" << endl;
                    exit(1);
                }
            }
            break;
        case bm::set_AND:
            bvt &= bv2;
            agg.combine_and(bv_agg, agg_list, 2);
            agg.combine_and_sub(bv_agg2, agg_list, 2, 0, 0, false);
            {
                if (bv_agg != bv_agg2)
                {
                    cerr << "Error: Aggregator AND - AND-SUB(0) comparison failed!" << endl;
                    exit(1);
                }
            }
            agg_check = true;
            
            // 2-way
            {
                bvect bvc(bv1);
                bvc &= bv2;
                
                bvect bvt1;
                bvt1.bit_and(bv1, bv2, bvect::opt_none);
                if (bvt1 != bvc)
                {
                    cerr << "1. AND 2-way check error!" << endl;
                    //print_bv(bvt1);
                    //cout << "control:" << endl;
                    //print_bv(bvc);
                    assert(0);
                    exit(1);
                }
                bvect bvt2;
                bvt2.bit_and(bv2, bv1, bvect::opt_compress);
                if (bvt2 != bvc)
                {
                    cerr << "2. AND 2-way check error!" << endl;
                    exit(1);
                }
            }
            break;
        case bm::set_SUB:
            bvt -= bv2;
            agg.combine_and_sub(bv_agg, agg_list, 1, agg_list2, 1, false);
            {
                bvect bv_h;
                agg.combine_and_sub_horizontal(bv_h, agg_list, 1, agg_list2, 1);
                if (bv_agg != bv_h)
                {
                    cerr << "Error: Aggregator Horz-AND-SUB comparison failed!" << endl;
                    exit(1);
                }
            }
            agg_check = true;
            // 2-way
            {
                bvect bvc1(bv1);
                bvect bvc2(bv2);
                bvc1 -= bv2;
                bvc2 -= bv1;
                
                bvect bvt1;
                bvt1.bit_sub(bv1, bv2, bvect::opt_compress);
                if (bvt1 != bvc1)
                {
                    DetailedCompareBVectors(bvt1, bvc1);
                    cerr << "1. SUB 2-way check error!" << endl;
                    exit(1);
                }
                bvect bvt2;
                bvt2.bit_sub(bv2, bv1, bvect::opt_compress);
                if (bvt2 != bvc2)
                {
                    cerr << "2. SUB 2-way check error!" << endl;
                    exit(1);
                }
            }

            break;
        default:
            goto no_compare;
        }
        if (bvt.compare(*bv_target) != 0)
        {
            cout << "Direct Serial operation comparison failed!" << endl;
            exit(1);
        }
        if (agg_check && bvt.compare(bv_agg) != 0)
        {
            cerr << "Error: Aggregator operation comparison failed!" << endl;
            exit(1);
        }

        no_compare:
        ;

    }
/*
    if (op == bm::set_AND || op == bm::set_OR || op == bm::set_XOR || op == bm::set_SUB)
    {
        cout << "3 way operation check... " << op << endl;
        operation_deserializer<bvect>::deserialize(*bv_target,
                                                   bv1,
                                                   smem2,
                                                   0,
                                                   op);
        cout << "OK" << endl;

        bvect bvt(bv1);
        switch(op)
        {
        case bm::set_OR:
            bvt |= bv2;
            break;
        case bm::set_XOR:
            bvt ^= bv2;
            break;
        case bm::set_AND:
            bvt &= bv2;
            break;
        case bm::set_SUB:
            bvt -= bv2;
            break;
        default:
            goto no_compare2;
        }
        if (bvt.compare(*bv_target) != 0)
        {
            cout << "3-way Serial operation comparison failed!" << endl;
            exit(1);
        }
        no_compare2:
        ;
    }
*/

   if (check_reverse)
   {
        cout << "Reverse check... " << endl;
        bvect bv_tmp2(BM_GAP);
        od.deserialize(bv_tmp2, smem2, 0, set_ASSIGN);
        res = bv_tmp2.compare(bv2);
        if (res != 0)
        {
            cout << "set_ASSIGN failed 2! " << endl;
            exit(1);
        }
        if (bv_tmp2 != bv2)
        {
            cout << "set_ASSIGN failed 2-1! " << endl;
            exit(1);
        }
        cout << "Deserialization assign to bv_tmp2 OK" << endl;
        unsigned count_rev =
        od.deserialize(bv_tmp2, smem1, 0, op);
        if (count != count_rev)
        {
//            print_stat(bv1);
/*
            unsigned c = count_or(bv1, bv2);
            cout << "Correct count=" << c << endl;

            c = count_or(bv2, bv1);
            cout << "Correct count=" << c << endl;

            bv1 |= bv2;
            cout << "Count3 = " << bv1.count() << endl;;
*/
            //SaveBVector("err1.bv", bv1);
            //SaveBVector("err2.bv", bv2);

            

            cout << "Operation=" << op << endl;

            cout << "Serialization operation reverse check failed"
                 << " count = " << count 
                 << " count rev= " << count_rev
                 << endl;
            cout << "See bvector dumps: err1.bv, err2.bv" << endl;

            exit(1);
        }

   }

   delete [] smem1;
   delete [] smem2;

   return count;
}

static
void SerializationOperation2Test(bvect*        bv_target,
                                 bvect&        bv1,
                                 bvect&        bv2,
                                 unsigned      predicted_count,
                                 set_operation op_count,
                                 set_operation op_combine)
{
    bv_target->clear(true);
    cout << "Serialization operation count..." << endl;

    unsigned scount1 = SerializationOperation(0, 
                                              bv1,
                                              bv2,
                                              op_count,
                                              true //reverse check
                                            );
    cout << "Serialization operation count OK." << endl;

    cout << "Serialization operation. " << endl;
    unsigned scount2 = SerializationOperation(bv_target,
        bv1,
        bv2,
        op_combine);
    scount2 = bv_target->count();
    if (predicted_count != scount2 || scount1 != scount2)
    {
        cout << "Serialization count != predicted" << endl
            << " predicted=" << predicted_count
            << " scount1=" << scount1
            << " scount2=" << scount2
            << endl;

        cout << endl << "target:" << endl;
        print_stat(*bv_target);
        cout << endl << endl << "Reference" << endl;
        if (op_combine == set_OR)
        {
            bv1 |= bv2;
            if (bv1 != *bv_target)
            {
                cout << "Comparison OR error!" << endl;
            }
            cout << "OR operation count=" << bv1.count() << endl;
            print_stat(bv1);
        }
        else
            if (op_combine == set_AND)
            {
                bv1 &= bv2;
                print_stat(bv1);
            }

        exit(1);
    }
    cout << "OK" << endl;
}

static
void print_mv(const bvect_mini &bvect_min, unsigned size)
{
    unsigned i;
    for (i = 0; i < size; ++i)
    {
        bool bflag = bvect_min.is_bit_true(i) != 0;

        if (bflag)
            printf("1");
        else
            printf("0");
        if ((i % 31) == 0 && (i != 0))
            printf(".");
    }

    printf("\n");
}

static
void print_gap(const gap_vector& gap_vect, unsigned /*size*/)
{
    const gap_word_t *buf = gap_vect.get_buf();
    unsigned len = gap_length(buf);
    printf("[%i:]", *buf++ & 1);

    for (unsigned i = 1; i < len; ++i)
    {
        printf("%i,", *buf++);
    }

    printf("\n");
}


static
void CheckGAPMin(const gap_vector& gapv, const bvect_mini& bvect_min, unsigned len)
{
    int last_bit = -1;
    for (unsigned i = 0; i < len; ++i)
    {
        int bit1 = (gapv.is_bit_true(i) == 1);
        int bit2 = (bvect_min.is_bit_true(i) != 0);
        if (bit1 != bit2)
        {
            cout << "Bit comparison failed. " << "Bit N=" << i << endl;
            assert(0);
            exit(1);
        }
        if (bvect_min.is_bit_true(i))
        {
            last_bit = (int)i;
        }
    }
    unsigned glast;
    bool found = gapv.get_last(&glast);
    if (!found && last_bit != -1)
    {
        cout << "Gap last search failed. " << "Bit=" << last_bit << endl;
        assert(0);
        exit(1);
    }

    if (found && last_bit == -1)
    {
        cout << "Gap last search ok but should failed. " << "Bit=" <<
             glast << endl;
        assert(0);
        exit(1);
    }

    if (last_bit != (int)glast)
    {
        if (!found && last_bit == -1)
        {}
        else
        {
            cout << "Gap last search discrepancy:" << " found Bit=" <<
                 glast << " last=" << last_bit << endl;
            assert(0);
            exit(1);
        }
    }
}

static
void CheckIntervals(const bvect& bv, unsigned /*max_bit*/)
{
    unsigned cnt0 = count_intervals(bv);
    unsigned cnt1 = 1;
    //bool bit_prev = bv.test(0);

    unsigned cnt2 = 0;
    bvect::enumerator en = bv.first();
    if (!en.valid())
    {
        cnt2 = 1;
    }
    else
    {
        if (*en > 0)
            ++cnt2;
        unsigned prev = *en;
        for (++en; en.valid(); ++en)
        {
            if (++prev == *en)
            {
            }
            else
            {
                cnt2 += 2;
                prev = *en;
            }
        }
        cnt2 += 2;
    }
/*
    for (unsigned i = 1; i < max_bit; ++i)
    {
        bool bit = bv.test(i);
        cnt1 += bit_prev ^ bit;
        bit_prev = bit;
    }
*/
    if (cnt0 != cnt2)
    {
        cout << "CheckIntervals error. " << "bm count=" << cnt0
             << " Control = " << cnt1 << endl;
        exit(1);
    }
}

template<class T> void CheckCountGapRange(const T& vect,
                                       unsigned left,
                                       unsigned right,
                                       unsigned* block_count_arr=0)
{
    unsigned cnt1 = vect.count_range(left, right, block_count_arr);
    unsigned cnt2 = 0;
    for (unsigned i = left; i <= right; ++i)
    {
        if (vect.test(i))
        {
            ++cnt2;
        }
    }
    if (cnt1 != cnt2)
    {
        cout << "Bitcount range failed!" << "left=" << left
             << " right=" << right << endl
             << "count_range()=" << cnt1
             << " check=" << cnt2;
        exit(1);
    }
}

template<typename T>
bool FindRank(const T& bv, bm::id_t rank, bm::id_t from, bm::id_t& pos)
{
    assert(rank);
    bool res = false;
    
    typename T::enumerator en = bv.get_enumerator(from);
    
    typename T::enumerator en2 = bv.get_enumerator(from);
    if (rank > 1)
        en2.skip_to_rank(rank);
    if (!en.valid())
        return false;
    for (; en.valid(); ++en)
    {
        rank -= en.valid();
        if (rank == 0)
        {
            pos = *en;
            res = true;
            break;
        }
    } // for en
    
    bm::id_t pos2 = *en2;
    if (pos != pos2)
    {
        cerr << "FindRank enumerator::skip() failed: "
        << "pos=" << pos << " skip()pos=" << pos2
        << " from=" << from
        << endl;
        exit(1);
    }
    
    return res;
}

inline
void CheckRangeCopy(const bvect& bv, unsigned from, unsigned to)
{
    bm::id_t f1, l1, f2, l2;
    
    bvect bv_cp(bv, from, to);
    bvect bv_cp2(bv, to, from); // swapped interval copy is legal

    bvect bv_control;
    bv_control.set_range(from, to);

    bv_control &= bv;
    
    int res = bv_control.compare(bv_cp);
    if (res != 0)
    {
        bool found1 =  bv_cp.find_range(f1, l1);
        bool found2 =  bv_control.find_range(f2, l2);
        
        bvect bv_cp3(bv, from, to);
        
        cerr << "Error: bvector<>::range_copy() failed. from=" << from << " to=" << to << endl;
        if (found1)
            cerr << " range copy from=" << f1 << " to=" << l1 << endl;
        if (found2)
            cerr << " control    from=" << f2 << " to=" << l2 << endl;
        exit(1);
    }

    int res2 = bv_control.compare(bv_cp2);
    if (res2 != 0)
    {
        bool found1 = bv_cp2.find_range(f1, l1);
        bool found2 = bv_control.find_range(f2, l2);

        cerr << "Error: reversed bvector<>::range_copy() failed. from=" << from << " to=" << to << endl;
        if (found1)
            cerr << " range copy from=" << f1 << " to=" << l1 << endl;
        if (found2)
            cerr << " control    from=" << f2 << " to=" << l2 << endl;
        exit(1);
    }

    bool found1 =  bv_cp.find_range(f1, l1);
    bool found2 =  bv_control.find_range(f2, l2);
    if (found1 != found2)
    {
        cerr << "Error: Dynamic range integrity check." << endl;
        exit(1);
    }
    if (found1)
    {
        if (f1 != f2 || l1 != l2)
        {
            cerr << "Error: bvector<>::range_copy() failed (dynamic range check). from=" << from << " to=" << to << endl;
            cerr << " range copy from=" << f1 << " to=" << l1 << endl;
            cerr << " control    from=" << f2 << " to=" << l2 << endl;
            exit(1);
        }
    }

}


template<class T> void CheckCountRange(const T& vect, 
                                       unsigned left, 
                                       unsigned right)
{
    unsigned cnt1 = vect.count_range(left, right);
    unsigned cnt2 = 0;

    typename T::enumerator en = vect.get_enumerator(left);
    for (; en.valid(); ++en)
    {
        if (*en > right)
            break;
        cnt2 += en.valid();
    }
    if (cnt1 != cnt2)
    {
        cout << "2. Bitcount range failed!" << "left=" << left
             << " right=" << right << endl
             << "count_range()=" << cnt1 
             << " check=" << cnt2;
        exit(1);
    }
    
    CheckRangeCopy(vect, left, right);
    
    bvect::rs_index_type bc_arr;
    vect.build_rs_index(&bc_arr);
    
    // run a cycle to check count_to()
    //
    //for (unsigned i = 0; i <= right; ++i)
    {
        if (left > right)
            swap(left, right);

        unsigned cnt_to_r = vect.count_to(right, bc_arr);
        cnt1 = vect.count_range(left, right, bc_arr);
        unsigned cnt_to_l = left ? vect.count_to(left - 1, bc_arr) : 0;
                 cnt2 = cnt_to_r - cnt_to_l;
        if (cnt1 != cnt2)
        {
            cout << "Bitcount range TO failed!" << " left=" << left
                 << " right=" << right << endl
                 << " count_range()=" << cnt1
                 << " check=" << cnt2;
            exit(1);
        }
        
        bm::id_t range = 1 + right - left;
        if (cnt1 > range)
        {
            cerr << "Impossible count_range detected!" << endl;
            exit(1);
        }
        
        
        if (cnt1) // check if we can reverse the search (rank)
        {
            unsigned pos, pos1;
            pos = pos1 = 0;
            bool rf = vect.find_rank(cnt1, left, pos);
            bool rf1 = vect.find_rank(cnt1, left, pos1, bc_arr);
            if (!rf || !rf1)
            {
                cerr << "1. find_rank() failed!" << " left=" << left
                     << " right=" << right
                     << " count_range()=" << cnt1
                     << " pos=" << pos
                     << " pos1=" << pos1
                     << " range=" << range
                     << endl;
                
                unsigned pos2;
                bool rf2 = FindRank(vect, cnt1, left, pos2);
                if (!rf2)
                {
                    cerr << "Debug FindRank failed!" << endl;
                }
                else
                {
                    cerr << " rank=" << pos2 << endl;
                }
                
                cerr << "detailed bug search..." << endl;
                for (unsigned k = 1; k <= cnt1; ++k)
                {
                    rf = vect.find_rank(k, left, pos);
                    rf2 = FindRank(vect, k, left, pos2);
                    if (rf != rf2 || pos != pos2)
                    {
                        rf = vect.find_rank(k, left, pos);
                        
                        cerr << "Failed for rank=" << k << endl;
                        cerr << rf << " " << rf2 << endl;
                        cerr << "pos = " << pos << " pos2 = " << pos2 << endl;
                        exit(1);
                    }
                }
                exit(1);
            }
            assert(pos == pos1);
            
            if (left > 0)
            {
                unsigned pos3;
                T bv1(vect, left, bm::id_max-1);
                std::unique_ptr<bvect::rs_index_type> bc_arr2(new bvect::rs_index_type);
                bv1.build_rs_index(bc_arr2.get());

                bool rf3 = bv1.select(cnt1, pos3, *bc_arr2);
                assert(rf3 == rf);
                assert(pos == pos3);
            }
            
            if (right != pos)
            {
                unsigned pos2;
                bool rf2 = FindRank(vect, cnt1, left, pos2);
                assert(rf2);
                // check if we found zero-tail
                //auto cnt3 = vect.count_range(pos+1, right, block_count_arr);
                if (pos2 != pos)
                {
                    rf = vect.find_rank(cnt1, left, pos);
                    cout << "2. find_rank() check failed! \n" << "left=" << left
                         << " right=" << right
                         << " count_range()=" << cnt1
                         << " pos=" << pos
                         << " rank = " << pos2
                         << endl;
                    exit(1);
                }
            }
        }
    }
}

static
unsigned BitCountChange(unsigned word)
{
    unsigned count = 1;
    unsigned bit_prev = word & 1;
    word >>= 1;
    for (unsigned i = 1; i < 32; ++i)
    {
        unsigned bit = word & 1;
        count += bit ^ bit_prev;
        bit_prev = bit;
        word >>= 1;
    }
    return count;
}

static
void DetailedCheckVectors(const bvect      &bv1,
                          const bvect      &bv2)
{
    bvect::enumerator en1 = bv1.first();
    bvect::enumerator en2 = bv2.first();
    
    while (en1.valid())
    {
        if (!en2.valid())
        {
            cout << "Second vector - invalid enumerator at:" << *en1;
            return;
        }
        if (*en1 != *en2)
        {
            cout << "Discrepancy at bit position: " << *en1;
            cout << " second vector is at:" << *en2;
            return;
        }
        ++en1;
        ++en2;
    }
    cout << "Detailed check OK" << endl; // Hmmm... Why?
}

static
void DetailedCheckVectors(const bvect_mini &bvect_min, 
                          const bvect      &bvect_full,
                          unsigned size)
{
    cout << "Detailed check" << endl;

    //bvect_full.stat();

    // detailed bit by bit comparison. Paranoia check.

    unsigned i;
    for (i = 0; i < size; ++i)
    {
        bool bv_m_flag = bvect_min.is_bit_true(i) != 0; 
        bool bv_f_flag = bvect_full.get_bit(i) != 0;

        if (bv_m_flag != bv_f_flag)
        {
            printf("Bit %u is non conformant. vect_min=%i vect_full=%i\n",
                i, (int)bv_m_flag, (int)bv_f_flag);

            cout << "Non-conformant block number is: " << unsigned(i >>  bm::set_block_shift) << endl;
            exit(1);
        }
    }
    
    printf("\n detailed check ok.\n");

}

static
void CompareEnumerators(const bvect::enumerator& en1, const bvect::enumerator& en2)
{
    if (!en1.valid() && !en2.valid())
        return;
    bool fsm_equal = en1.compare_state(en2);
    if (!fsm_equal)
    {
        cerr << "Enumerators FSM comparison failed" << endl;
        exit(1);
    }
}

// find last set bit by scan (not optimal)
//
static
bool FindLastBit(const bvect& bv, bm::id_t& last_pos)
{
    bvect::enumerator en = bv.first();
    if (!en.valid())
        return false;
    for (; en.valid(); ++en)
    {
        last_pos = *en;
    }
    return true;
}


// vectors comparison check

void CheckVectors(bvect_mini &bvect_min, 
                  bvect      &bvect_full,
                  unsigned size,
                  bool     detailed)
{
    cout << "\nVectors checking...bits to compare = " << size << endl;

    cout << "Bitcount summary : " << endl;
    unsigned min_count = bvect_min.bit_count();
    cout << "minvector count = " << min_count << endl;
    unsigned count = bvect_full.count();
    unsigned full_count = bvect_full.recalc_count();
    cout << "fullvector re-count = " << full_count << endl;
    
    if (min_count != full_count)
    {
        cout << "fullvector count = " << count << endl;
        cout << "Count comparison failed !!!!" << endl;
        print_stat(bvect_full);
        DetailedCheckVectors(bvect_min, bvect_full, size);
        assert(0);
        exit(1);  
    } 

    if (full_count)
    {
        bool any = bvect_full.any();
        if (!any)
        {
            cout << "Anycheck failed!" << endl;
            exit(1);
        }
    }

    // find_last check
    {
        bm::id_t pos1 = 0;
        bm::id_t pos2 = 0;
        bool last_found1 = FindLastBit(bvect_full, pos1);
        bool last_found2 = bvect_full.find_reverse(pos2);
        
        assert(last_found1 == last_found2);
        if (last_found1)
        {
            assert(pos1 == pos2);
        }
    }
    
    if (!detailed)
        return;
    
    // get_next comparison
    cout << "Positive bits comparison..." << flush;
    unsigned nb_min = bvect_min.get_first();
    unsigned nb_ful = bvect_full.get_first();

    bvect::counted_enumerator en = bvect_full.first();
    unsigned nb_en = *en;
    bvect::enumerator en1 = bvect_full.get_enumerator(*en);
    if (nb_min != nb_ful)
    {
         cout << "!!!! First bit comparison failed. Full id = " 
              << nb_ful << " Min id = " << nb_min 
              << endl;

         bool bit_f = bvect_full.get_bit(nb_ful);
         cout << "Full vector'd bit #" << nb_ful << "is:" 
              << bit_f << endl;

         bool bit_m = (bvect_min.is_bit_true(nb_min) == 1);
         cout << "Min vector'd bit #" << nb_min << "is:" 
              << bit_m << endl;


         print_stat(bvect_full);

         DetailedCheckVectors(bvect_min, bvect_full, size);

         exit(1);
    }
    CompareEnumerators(en, en1);

    if (full_count)
    {
       unsigned bit_count = 1;
       unsigned en_prev = nb_en;

       do
       {
           nb_min = bvect_min.get_next(nb_min);
           if (nb_min == 0)
           {
               break;
           }

           en_prev = nb_en;
           ++en;
           ++en1;
           CompareEnumerators(en, en1);

           if ((bit_count % 10 == 0) || (bit_count % 128 == 0))
           {
                bvect::enumerator en2 = bvect_full.get_enumerator(*en);
                CompareEnumerators(en, en2);
           }

           nb_en = *en;
           ++bit_count;

           if (nb_en != nb_min)
           {
               nb_ful = bvect_full.get_next(en_prev);
               cout << "!!!!! next bit comparison failed. Full id = " 
                    << nb_ful << " Min id = " << nb_min 
                    << " Enumerator = " << nb_en
                    << endl;

     //          bvect_full.stat();

     //          DetailedCheckVectors(bvect_min, bvect_full, size);

               exit(1);
           }
       } while (en.valid());
       if (bit_count != min_count)
       {
           cout << " Bit count failed."  
                << " min = " << min_count 
                << " bit = " << bit_count 
                << endl;
           exit(1);
       }
    }

    cout << "OK" << endl;

    return;
}

static
void DynamicMatrixTest()
{
    cout << "---------------------------- DynamicMatrixTest() test" << endl;
    
    {
        bm::dynamic_heap_matrix<unsigned, bvect::allocator_type> matr(0, 0);
        matr.resize(3, bm::set_sub_array_size);
        matr.set_zero();
        
        {
            unsigned* r = matr.row(1);
            for (unsigned i = 0; i < matr.cols(); ++i)
            {
                r[i] = i;
            }
        }
        
        matr.resize(matr.rows()+1, matr.cols());
        {
            unsigned* r = matr.row(3);
            for (unsigned i = 0; i < matr.cols(); ++i)
            {
                r[i] = i;
            }
        }

        {
            const unsigned* r = matr.row(1);
            for (unsigned i = 0; i < matr.cols(); ++i)
            {
                assert(r[i] == i);
            }
        }


    }
    
    cout << "---------------------------- DynamicMatrixTest() test OK" << endl;
}

static
void RSIndexTest()
{
    cout << "---------------------------- RSIndexTest() test" << endl;

    {
        rs_ind rsi;
        
        rsi.resize(bm::set_sub_array_size*5);
        rsi.resize_effective_super_blocks(3);
        
        rsi.set_super_block(0, 1);
        rsi.set_super_block(1, 2);
        rsi.set_super_block(2, 3);
        rsi.set_null_super_block(3);
        rsi.set_full_super_block(4);

        auto sb_size = rsi.super_block_size();
        assert(sb_size == 5);
        
        auto rc = rsi.get_super_block_count(0);
        assert(rc == 1);
        rc = rsi.get_super_block_count(1);
        assert(rc == 2);
        rc = rsi.get_super_block_count(2);
        assert(rc == 3);
        rc = rsi.get_super_block_count(256);
        assert(rc == 0);


        auto bc = rsi.get_super_block_rcount(0);
        assert(bc == 1);
        bc = rsi.get_super_block_rcount(1);
        assert(bc == 3);
        bc = rsi.get_super_block_rcount(2);
        assert(bc == 6);
        
        unsigned i = rsi.find_super_block(1);
        assert(i == 0);
        i = rsi.find_super_block(2);
        assert(i == 1);
        i = rsi.find_super_block(3);
        assert(i == 1);
        i = rsi.find_super_block(4);
        assert(i == 2);
        i = rsi.find_super_block(200);
        assert(i == 4);
/*
        i = rsi.find_super_block(bm::id_max);
        assert(i == 5);
*/
    }
    
    {
        unsigned bcount[bm::set_sub_array_size];
        unsigned sub_count1[bm::set_sub_array_size];
        unsigned sub_count2[bm::set_sub_array_size];
        for (unsigned i = 0; i < bm::set_sub_array_size; ++i)
        {
            bcount[i] = sub_count1[i] = sub_count2[i] = 0;
        } // for
        bcount[0] = 1;
        bcount[255] = 2;
        
        sub_count1[0] = 1;          // sub-range 0
        sub_count1[255] = 0;        // sub-3
        sub_count2[0] = 1 << 16;    // sub-2
        sub_count2[255] = 1 << 16;  // sub 2,3

        
        rs_ind rsi;
        // -----------------------------------------
        rsi.resize(bm::set_sub_array_size*4);
        rsi.resize_effective_super_blocks(2);
        rsi.set_total(bm::set_sub_array_size*4);
        
        
        rsi.set_null_super_block(0);
        auto tcnt = rsi.count();
        assert(tcnt == 0);
        rsi.register_super_block(1, &bcount[0], &sub_count1[0]);
        tcnt = rsi.count();
        assert(tcnt == 3);
        rsi.register_super_block(2, &bcount[0], &sub_count2[0]);
        tcnt = rsi.count();
        assert(tcnt == 6);
        rsi.set_full_super_block(3);
        tcnt = rsi.count();
        assert(tcnt == 6 + bm::set_sub_array_size * 65536);

        unsigned i = rsi.find_super_block(1);
        assert(i == 1);
        i = rsi.find_super_block(3);
        assert(i == 1);
        i = rsi.find_super_block(4);
        assert(i == 2);
        i = rsi.find_super_block(400);
        assert(i == 3);
//        i = rsi.find_super_block(bm::id_max);
//        assert(i == rsi.super_block_size());
        
        unsigned bc;
        rs_ind::size_type rbc;
        for (unsigned nb = 0; nb < bm::set_sub_array_size; ++nb)
        {
            bc = rsi.count(nb);
            assert(bc == 0);
            rbc = rsi.rcount(nb);
            assert(!rbc);
        }
        bc = rsi.count(bm::set_sub_array_size);
        assert(bc == 1);
        rbc = rsi.rcount(bm::set_sub_array_size);
        assert(rbc == 1);
        
        bc = rsi.count(bm::set_sub_array_size+1);
        assert(bc == 0);
        rbc = rsi.rcount(bm::set_sub_array_size+1);
        assert(rbc == 1);

        bc = rsi.count(bm::set_sub_array_size + 255);
        assert(bc == 2);
        rbc = rsi.rcount(bm::set_sub_array_size + 255);
        assert(rbc == 3);

        bc = rsi.count(bm::set_sub_array_size*3);
        assert(bc == 65536);
        rbc = rsi.rcount(bm::set_sub_array_size*3);
        assert(rbc == 65536+6);
        rbc = rsi.rcount(bm::set_sub_array_size*3 + 1);
        assert(rbc == 65536+6 + 65536);
        
        
        // ==========================
        {
            auto nb = rsi.find(1);
            assert(nb == 256);
            
            nb = rsi.find(2);
            assert(nb == bm::set_sub_array_size+255);
            nb = rsi.find(3);
            assert(nb == bm::set_sub_array_size+255);

            nb = rsi.find(4);
            assert(nb == bm::set_sub_array_size+255+1);

            nb = rsi.find(65536);
            assert(nb == 3*bm::set_sub_array_size+0);
            nb = rsi.find(65536*2);
            assert(nb == 3*bm::set_sub_array_size+1);
            nb = rsi.find(65536*3);
            assert(nb == 3*bm::set_sub_array_size+2);
        }
        // ==========================

        {
            bool b;
            rs_ind::size_type rank;
            rs_ind::block_idx_type nb;
            bm::gap_word_t sub_range;

            rank = 1;
            b = rsi.find(&rank, &nb, &sub_range);
            assert(b);
            assert(nb == 256);
            assert(sub_range == 0);
            assert(rank == 1);

            rank = 2;
            b = rsi.find(&rank, &nb, &sub_range);
            assert(b);
            assert(nb == bm::set_sub_array_size+255);
            assert(sub_range == bm::rs3_border1 + 1);
            assert(rank == 1);

            rank = 3;
            b = rsi.find(&rank, &nb, &sub_range);
            assert(b);
            assert(nb == bm::set_sub_array_size+255);
            assert(sub_range == bm::rs3_border1 + 1);
            assert(rank == 2);

            rank = 4;
            b = rsi.find(&rank, &nb, &sub_range);
            assert(b);
            assert(nb == bm::set_sub_array_size+255+1);
            assert(sub_range == bm::rs3_border0 + 1);
            assert(rank == 1);

            rank = 5;
            b = rsi.find(&rank, &nb, &sub_range);
            assert(b);
            assert(nb == bm::set_sub_array_size+256+255);
            assert(sub_range == bm::rs3_border0 + 1);
            assert(rank == 1);

            rank = 6;
            b = rsi.find(&rank, &nb, &sub_range);
            assert(b);
            assert(nb == bm::set_sub_array_size+256+255);
            assert(sub_range == bm::rs3_border1 + 1);
            assert(rank == 1);

            rank = 65536;
            b = rsi.find(&rank, &nb, &sub_range);
            assert(b);
            assert(nb == 3*bm::set_sub_array_size+0);
            assert(sub_range == bm::rs3_border1 + 1);
            assert(rank == 65536 - 6 - bm::rs3_border1);

            rank = 65536 + 7;
            b = rsi.find(&rank, &nb, &sub_range);
            assert(b);
            assert(nb == 3*bm::set_sub_array_size+1);
            assert(sub_range == 0);
            assert(rank == 1);

            rank = bm::id_max;
            b = rsi.find(&rank, &nb, &sub_range);
            assert(!b);


        }

    }
    

    cout << "---------------------------- RSIndexTest() test OK" << endl;
}


static
void ClearAllTest()
{
    bvect     bvect_full;

    for (unsigned i = 0; i < 100000; ++i)
    {
        bvect_full.set_bit(i);
    }
    BM_DECLARE_TEMP_BLOCK(tb)
    bvect_full.optimize(tb);
    bvect_full.clear();

    print_stat(bvect_full);

    unsigned count = bvect_full.count();
    assert(count == 0);
    print_stat(bvect_full);
}

static
void WordCmpTest()
{
    cout << "---------------------------- WordCmp test" << endl;

    for (int i = 0; i < 10000000; ++i)
    {
        unsigned w1 = unsigned(rand());
        unsigned w2 = unsigned(rand());
        int res = wordcmp0(w1, w2);
        int res2 = wordcmp(w1, w2);
        if (res != res2)
        {
            printf("WordCmp failed !\n");
            exit(1);
        }

        res = wordcmp0((unsigned)0U, (unsigned)w2);
        res2 = wordcmp((unsigned)0U, (unsigned)w2);

        if (res != res2)
        {
            printf("WordCmp 0 test failed !\n");
            exit(1);
        }

        res = wordcmp0((unsigned)~0U, (unsigned)w2);
        res2 = wordcmp((unsigned)~0U, (unsigned)w2);

        if (res != res2)
        {
            printf("WordCmp ~0 test failed !\n");
            exit(1);
        }

        res = wordcmp0((unsigned)w2, (unsigned)0);
        res2 = wordcmp((unsigned)w2, (unsigned)0);

        if (res != res2)
        {
            printf("WordCmp 0-2 test failed !\n");
            exit(1);
        }

    }

    cout << "Ok." << endl;
}


static
void TestBlockCountChange()
{
    cout << "---------------------------- CountChange test" << endl;

#ifdef VECT_BLOCK_CHANGE

    unsigned i, c, cc;
    bm::word_t BM_VECT_ALIGN blk[bm::set_block_size] BM_VECT_ALIGN_ATTR = { 0 };

    for (i = 0; i < bm::set_block_size; ++i)
        blk[i] = 0;
    
    c = VECT_BLOCK_CHANGE(blk);
    cc = bm::bit_block_change32(blk);
    assert(c == 1);
    assert(c == cc);

    blk[0] = 1;
    c = VECT_BLOCK_CHANGE(blk);
    cc = bm::bit_block_change32(blk);
    assert(c == 2);
    assert(c == cc);

    blk[0] = 0xFF;
    c = VECT_BLOCK_CHANGE(blk);
    cc = bm::bit_block_change32(blk);
    assert(c == 2);
    assert(c == cc);

    blk[0] = ~0u;
    c = VECT_BLOCK_CHANGE(blk);
    cc = bm::bit_block_change32(blk);
    assert(c == 2);
    assert(c == cc);

    blk[0] = blk[1] = blk[2] = blk[3] = 2;
    c = VECT_BLOCK_CHANGE(blk);
    cc = bm::bit_block_change32(blk);
    assert(c == cc);
    
    blk[4] = blk[5] = blk[6] = blk[7] = 2;
    c = VECT_BLOCK_CHANGE(blk);
    cc = bm::bit_block_change32(blk);
    assert(c == cc);

    {
    for (i = 0; i < bm::set_block_size; ++i)
        blk[i] = 2;

    c = VECT_BLOCK_CHANGE(blk);
    cc = bm::bit_block_change32(blk);
    assert(c == cc);
    }
    
    {
    for (i = 0; i < bm::set_block_size; ++i)
        blk[i] = 1u << 31;

    c = VECT_BLOCK_CHANGE(blk);
    cc = bm::bit_block_change32(blk);
    assert(c == cc);
    }

    {
    for (i = 0; i < bm::set_block_size; ++i)
        blk[i] = ~0u << 30;

    c = VECT_BLOCK_CHANGE(blk);
    cc = bm::bit_block_change32(blk);
    assert(c == cc);
    }

    cout << "Block change stress..." << endl;
    {
    std::chrono::time_point<std::chrono::steady_clock> s;
    std::chrono::time_point<std::chrono::steady_clock> f;
    s = std::chrono::steady_clock::now();

    
        unsigned k_max = (1u << 31) / 4;
        for (unsigned k = 0; k <= k_max; ++k)
        {
            for (i = 0; i < bm::set_block_size; ++i)
                blk[i] = k;
            c = VECT_BLOCK_CHANGE(blk);
            cc = bm::bit_block_change32(blk);
            assert(c == cc);
            
            if (k % 100000 == 0)
            {
                f = std::chrono::steady_clock::now();
                auto diff = f - s;
                auto d = std::chrono::duration <double, std::milli> (diff).count();

                cout << "\r" << k << " / " << k_max << " (" << d << "ms)" << flush;
                
                s = std::chrono::steady_clock::now();
            }
        } // for k
    }
    cout << endl;

#endif

    
    cout << "---------------------------- CountChange test OK" << endl;
}


inline
bm::id64_t bit_block_calc_xor_change_digest(
                        const bm::word_t*  block,
                        const bm::word_t*  xor_block,
                        block_waves_xor_descr&  x_descr)
{
    unsigned bgain;
    bm::compute_complexity_descr(block, x_descr);
    return bm::compute_xor_complexity_descr(block, xor_block, x_descr, bgain);
}

static
void Check_XOR_Product(const bm::word_t*  block,
                       const bm::word_t*  xor_block,
                       bm::id64_t         digest)
{
    assert(digest);

    bm::word_t BM_VECT_ALIGN t_blk1[bm::set_block_size] BM_VECT_ALIGN_ATTR = { 0 };
    bm::word_t BM_VECT_ALIGN t_blk2[bm::set_block_size] BM_VECT_ALIGN_ATTR = { 0 };

    bm::bit_block_xor(t_blk1, block, xor_block, digest);
    bm::bit_block_xor(t_blk2, t_blk1, xor_block, digest);

    unsigned cnt = bm::bit_block_xor_count(block, t_blk2);
    assert(cnt == 0); // identically restored
}


static
void TestBlockCountXORChange()
{
    cout << "---------------------------- TestBlockCountXORChange() test" << endl;
    unsigned i;
    bm::id64_t d64;
    bm::block_waves_xor_descr x_descr;

    {
        bm::word_t BM_VECT_ALIGN blk[bm::set_block_size] BM_VECT_ALIGN_ATTR = { 0 };
        bm::word_t BM_VECT_ALIGN blk_xor[bm::set_block_size] BM_VECT_ALIGN_ATTR = { 0 };

        for (i = 0; i < bm::set_block_size; ++i)
            blk[i] = blk_xor[i] = 0;

        d64 = bit_block_calc_xor_change_digest(blk, blk_xor, x_descr);
        assert(d64 == ~0ull);
        for (unsigned k = 0; k < bm::block_waves; ++k)
        {
            assert(x_descr.sb_change[k] == 1);
            assert(x_descr.sb_xor_change[k] == 1);
        } // for k

        blk[0] = 1;
        d64 = bit_block_calc_xor_change_digest(blk, blk_xor, x_descr);
        assert(d64);
        assert(x_descr.sb_change[0] == 2);
        assert(x_descr.sb_xor_change[0] == 2);
        for (unsigned k = 1; k < bm::block_waves; ++k)
        {
            assert(x_descr.sb_change[k] == 1);
            assert(x_descr.sb_xor_change[k] == 1);
        } // for k


        blk[0] = 1; blk_xor[0] = 1;
        d64 = bit_block_calc_xor_change_digest(blk, blk_xor, x_descr);
        cout << x_descr.sb_xor_change[0] << endl;
        assert(x_descr.sb_change[0] == 2);
        // next assert hides non-critical discrepancy between SIMD versions
        assert(x_descr.sb_xor_change[0] == 1 || x_descr.sb_xor_change[0] == 0);
        assert(d64 == ~0ull);
        for (unsigned k = 1; k < bm::block_waves; ++k)
        {
            assert(x_descr.sb_change[k] == 1);
            assert(x_descr.sb_xor_change[k] == 1);
        } // for k

        Check_XOR_Product(blk, blk_xor, d64);

        blk[0] = (1 << 10) | (1 << 12); blk_xor[0] = (1 << 11);
        unsigned off = (60 * bm::set_block_digest_wave_size);
        blk[off] = (1 << 10) | (1 << 12);

        d64 = bit_block_calc_xor_change_digest(blk, blk_xor, x_descr);
        assert(x_descr.sb_change[0] == 5);
        assert(x_descr.sb_xor_change[0] == 3);
        assert((d64 & 1));

        Check_XOR_Product(blk, blk_xor, d64);

        for (unsigned k = 1; k < bm::block_waves; ++k)
        {
            if (k!= 60)
            {
                assert(x_descr.sb_change[k] == 1);
                assert(x_descr.sb_xor_change[k] == 1);
            }
            else
            {
                assert(x_descr.sb_change[60] == 5);
                assert(x_descr.sb_xor_change[60] == 5);
            }
        } // for k

        blk_xor[off] = (1 << 10) | (1 << 11) | (1 << 12);
        d64 = bit_block_calc_xor_change_digest(blk, blk_xor, x_descr);
        assert(x_descr.sb_change[0] == 5);
        assert(x_descr.sb_xor_change[0] == 3);
        assert((d64 & 1) && (d64 & (1ull << 60)));

        for (unsigned k = 1; k < bm::block_waves; ++k)
        {
            if (k!= 60)
            {
                assert(x_descr.sb_change[k] == 1);
                assert(x_descr.sb_xor_change[k] == 1);
            }
            else
            {
                assert(x_descr.sb_change[60] == 5);
                assert(x_descr.sb_xor_change[60] == 3);
            }
        } // for k

        Check_XOR_Product(blk, blk_xor, d64);

    }

    cout << "---------------------------- TestBlockCountXORChange() test OK" << endl;
}


/*!
   \brief Converts bit block to GAP.
   \param dest - Destinatio GAP buffer.
   \param src - Source bitblock buffer.
   \param bits - Number of bits to convert.
   \param dest_len - length of the dest. buffer.
   \return  New length of GAP block or 0 if conversion failed
   (insufficicent space).

   @ingroup gapfunc
*/
template<typename T>
unsigned bit_convert_to_gap(T*  dest,
                            const unsigned*  src,
                            bm::id_t bits,
                            unsigned dest_len)
{
    T*  pcurr = dest;
    T*  end = dest + dest_len;
    unsigned bitval = (*src) & 1u;
    *pcurr = (T)bitval;

    ++pcurr;
    *pcurr = 0;
    unsigned bit_idx = 0;
    unsigned bitval_next;

    unsigned val = *src;

    do
    {
        // We can fast pace if *src == 0 or *src = 0xffffffff

        while (val == 0 || val == 0xffffffff)
        {
           bitval_next = val ? 1 : 0;
           if (bitval != bitval_next)
           {
               *pcurr++ = (T)(bit_idx-1);
               assert((pcurr-1) == (dest+1) || *(pcurr-1) > *(pcurr-2));
               if (pcurr >= end)
               {
                   return 0; // OUT of memory
               }
               bitval = bitval_next;
           }
           bit_idx += unsigned(sizeof(*src) * 8);
           if (bit_idx >= bits)
           {
               goto complete;
           }
           ++src;
           val = *src;
        }

        unsigned mask = 1;
        while (mask)
        {
            // Now plain bitshifting. TODO: Optimization wanted.

            bitval_next = val & mask ? 1 : 0;
            if (bitval != bitval_next)
            {
                *pcurr++ = (T)(bit_idx-1);
                assert((pcurr-1) == (dest+1) || *(pcurr-1) > *(pcurr-2));
                bitval = bitval_next;
                if (pcurr >= end)
                    return 0; // OUT of memory
            }
            mask <<= 1;
            ++bit_idx;
        } // while mask

        if (bit_idx >= bits)
            goto complete;

        ++src;
        val = *src;

    } while(1);

complete:
    *pcurr = (T)(bit_idx-1);
    unsigned len = (unsigned)(pcurr - dest);
    *dest = (T)((*dest & 7) + (len << 3));
    return len;
}


static
void TestBlockToGAP()
{
    cout << "---------------------------- TestBlockToGAP" << endl;
    
    unsigned i;
    bm::word_t BM_VECT_ALIGN blk[bm::set_block_size] BM_VECT_ALIGN_ATTR = { 0 };

    for (i = 0; i < bm::set_block_size; ++i)
        blk[i] = 0;
    
    {
       gap_vector gapv1(0);
       gap_vector gapv2(0);
       gap_word_t* gap_buf1 = gapv1.get_buf();
       *gap_buf1 = 0;
       unsigned len1 = bit_convert_to_gap(gap_buf1, blk, bm::gap_max_bits, bm::gap_max_buff_len);
       gap_word_t* gap_buf2 = gapv2.get_buf();
       unsigned len2 = bm::bit_to_gap(gap_buf2, blk, bm::gap_max_buff_len);
       print_gap(gapv2, 100);
       assert(len1 == len2);
       int cmp = bm::gapcmp(gap_buf1, gap_buf2);
       assert(cmp == 0);

       unsigned pos;
       bool f = bm::gap_find_first_diff(gap_buf1, gap_buf2, &pos);
       assert(!f);
    }

    unsigned test_arr[] = { 1, 2, 3, (~0u << 4), ~0u, (~0u >> 1), (~0u >> 2) };
    unsigned arr_size = sizeof(test_arr)/sizeof(test_arr[0]);
    for (unsigned k = 0; k < bm::set_block_size; ++k)
    {
        for (i = 0; i < bm::set_block_size; ++i)
            blk[i] = 0;
        for (i = 0; i < arr_size; ++i)
        {
            blk[k] = test_arr[i];
            
           gap_vector gapv1(0);
           gap_vector gapv2(0);
           gap_word_t* gap_buf1 = gapv1.get_buf();
           *gap_buf1 = 0;
           unsigned len1 = bit_convert_to_gap(gap_buf1, blk, bm::gap_max_bits, bm::gap_max_buff_len);
           print_gap(gapv1, 100);
           gap_word_t* gap_buf2 = gapv2.get_buf();
           unsigned len2 = bm::bit_to_gap(gap_buf2, blk, bm::gap_max_buff_len);
           assert(len1 == len2);
           assert(len1);
           print_gap(gapv2, 100);
           int cmp = bm::gapcmp(gap_buf1, gap_buf2);
           assert(cmp == 0);

           unsigned pos;
           bool f = bm::gap_find_first_diff(gap_buf1, gap_buf2, &pos);
           assert(!f);
        }
    } // for k
    
    cout << "Test arr - ok" << endl;
    
    {
       gap_vector gapv1(0);
       gap_vector gapv2(0);
       gap_word_t* gap_buf1 = gapv1.get_buf();
       *gap_buf1 = 0;
       gap_word_t* gap_buf2 = gapv2.get_buf();
       *gap_buf2 = 0;

        unsigned mask = 1u;
        for (unsigned k = 0; k < bm::set_block_size; ++k)
        {
            for (i = 0; i < bm::set_block_size; ++i)
                blk[i] = 0;
            
            blk[k] = mask;
            mask <<= 1;
            if (!mask)
                mask = 1u;
            
            *gap_buf1 = 0;
            *gap_buf2 = 0;
            unsigned len1 = bit_convert_to_gap(gap_buf1, blk, bm::gap_max_bits, bm::gap_max_buff_len);
            unsigned len2 = bm::bit_to_gap(gap_buf2, blk, bm::gap_max_buff_len);
            assert(len1);
            assert(len1 == len2);
            if (len1)
            {
               int cmp = bm::gapcmp(gap_buf1, gap_buf2);
               assert(cmp == 0);

               unsigned pos;
               bool f = bm::gap_find_first_diff(gap_buf1, gap_buf2, &pos);
               assert(!f);
            }
        }
    }
    
    cout << "mask shift - ok" << endl;

    {
       gap_vector gapv1(0);
       gap_vector gapv2(0);
       gap_word_t* gap_buf1 = gapv1.get_buf();
       *gap_buf1 = 0;
       gap_word_t* gap_buf2 = gapv2.get_buf();
       *gap_buf2 = 0;

        unsigned max_try = 1000000;
        for (unsigned k = 0; k < max_try; ++k)
        {
            for (i = 0; i < bm::set_block_size; ++i)
                blk[i] = 0;
            
            unsigned idx = (unsigned)rand() % bm::set_block_size;
            blk[idx] |= (unsigned)rand();
            idx = (unsigned)rand() % bm::set_block_size;
            blk[idx] |= (unsigned)rand();
            idx = (unsigned)rand() % bm::set_block_size;
            blk[idx] |= (unsigned)rand();
            idx = (unsigned)rand() % bm::set_block_size;
            blk[idx] |= (unsigned)rand();

            idx = (unsigned)rand() % bm::set_block_size;
            blk[idx] |= ~0u;
            idx = (unsigned)rand() % bm::set_block_size;
            blk[idx] |= (~0u << (rand()%31));

            *gap_buf1 = 0;
            *gap_buf2 = 0;
            unsigned len1 = bit_convert_to_gap(gap_buf1, blk, bm::gap_max_bits, bm::gap_max_buff_len);
            unsigned len2 = bm::bit_to_gap(gap_buf2, blk, bm::gap_max_buff_len);
            assert(len1);
            assert(len1 == len2);
            if (len1)
            {
               int cmp = bm::gapcmp(gap_buf1, gap_buf2);
               assert(cmp == 0);
               unsigned pos;
               bool f = bm::gap_find_first_diff(gap_buf1, gap_buf2, &pos);
               assert(!f);
            }
        }
    }


    cout << "---------------------------- TestBlockToGAP  OK" << endl;

}


static
void ShiftRotateTest()
{
    cout << "---------------------------- ShiftRotate test" << endl;

    bm::word_t BM_VECT_ALIGN blk0[bm::set_block_size] BM_VECT_ALIGN_ATTR = { 0 };
    bm::word_t BM_VECT_ALIGN blk1[bm::set_block_size] BM_VECT_ALIGN_ATTR = { 0 };
    unsigned i;

    for (i = 0; i < bm::set_block_size; ++i)
    {
        blk0[i] = blk1[i] = 1;
    }

    bm::bit_block_rotate_left_1(blk0);
    bm::bit_block_rotate_left_1_unr(blk1);

    for (i = 0; i < bm::set_block_size; ++i)
    {
        if (blk0[i] != 2 || blk0[i] != blk1[i])
        {
            cerr << "Cyclic rotate check failed" << endl;
            exit(1);
        }
    }

    bm::bit_block_rotate_left_1(blk0);
    bm::bit_block_rotate_left_1_unr(blk1);

    for (i = 0; i < bm::set_block_size; ++i)
    {
        if (blk0[i] != 4 || blk0[i] != blk1[i])
        {
            cerr << "Cyclic rotate check failed" << endl;
            exit(1);
        }
    }

    for (i = 0; i < bm::set_block_size; ++i)
    {
        blk0[i] = blk1[i] = 1u << 31;
    }
    bm::bit_block_rotate_left_1(blk0);
    bm::bit_block_rotate_left_1_unr(blk1);

    for (i = 0; i < bm::set_block_size; ++i)
    {
        if (blk0[i] != 1 || blk0[i] != blk1[i])
        {
            cerr << "Cyclic rotate check failed" << endl;
            exit(1);
        }
    }

    for (i = 0; i < bm::set_block_size; ++i)
    {
        blk0[i] = blk1[i] = unsigned(rand());
    }

    for (unsigned j = 0; j < bm::set_block_size * 32; ++j)
    {
        bm::bit_block_rotate_left_1(blk0);
        bm::bit_block_rotate_left_1_unr(blk1);

        for (i = 0; i < bm::set_block_size; ++i)
        {
            if (blk0[i] != blk1[i])
            {
                cerr << "Stress Cyclic rotate check failed" << endl;
                exit(1);
            }
        }
    }
    
    // SHIFT-R tests
    //

    unsigned acc0, acc1;

    for (i = 0; i < bm::set_block_size; ++i)
    {
        blk0[i] = blk1[i] = 1;
    }

    bm::bit_block_shift_r1(blk0, &acc0, 0);
    bm::bit_block_shift_r1_unr(blk1, &acc1, 0);

    for (i = 0; i < bm::set_block_size; ++i)
    {
        if (blk0[i] != blk1[i])
        {
            cerr << "1. SHIFT-r check failed" << endl;
            exit(1);
        }
        assert(blk0[i] == 2);
    }


    for (i = 0; i < bm::set_block_size; ++i)
    {
        blk0[i] = blk1[i] = (1u << 31);
    }

    bm::bit_block_shift_r1(blk0, &acc0, 0);
    bm::bit_block_shift_r1_unr(blk1, &acc1, 0);

    for (i = 0; i < bm::set_block_size; ++i)
    {
        if (blk0[i] != blk1[i])
        {
            cerr << "2. SHIFT-r check failed" << endl;
            exit(1);
        }
    }

    for (i = 0; i < bm::set_block_size; ++i)
    {
        blk0[i] = blk1[i] = unsigned(rand());
    }

    for (unsigned j = 0; j < bm::set_block_size * 32; ++j)
    {
        bm::bit_block_shift_r1(blk0, &acc0, 0);
        bm::bit_block_shift_r1_unr(blk1, &acc1, 0);
        
        assert(bool(acc0) == bool(acc1));

        for (i = 0; i < bm::set_block_size; ++i)
        {
            if (blk0[i] != blk1[i])
            {
                cerr << "Stress SHIFT-r check failed" << endl;
                exit(1);
            }
        }
    }

    cout << "---------------------------- ShiftRotate test OK" << endl;
}

static
void BlockBitInsertTest()
{
    cout << "---------------------------- BlockBitInsertTest test" << endl;
    
    bm::word_t BM_VECT_ALIGN blk0[bm::set_block_size] BM_VECT_ALIGN_ATTR = { 0 };
    bm::word_t BM_VECT_ALIGN blk1[bm::set_block_size] BM_VECT_ALIGN_ATTR = { 0 };
    
    unsigned i;
    for (i = 0; i < bm::set_block_size; ++i)
        blk0[i] = blk1[i] = 0;
    
    {
        bm::word_t co = bm::bit_block_insert(blk0, 0, 1);
        assert(co == 0);
        assert(blk0[0]==1u);
        co = bm::bit_block_insert(blk0, 0, 1);
        assert(co == 0);
        assert(blk0[0]==3u);
        
        blk0[bm::set_block_size-1] = ~0u;
        co = bm::bit_block_insert(blk0, 0, 0);
        assert(co == 1);
        assert(blk0[0]==(3u << 1));
        assert(blk0[bm::set_block_size-1] == (~0u << 1));
        
        blk0[0] = ~0u;
        co = bm::bit_block_insert(blk0, 1, 0);
        assert(co == 1);
        assert(blk0[0]==(~0u & ~(1u << 1)));
        assert(blk0[bm::set_block_size-1] == (~0u << 2));
        
        blk0[0] = ~0u;
        co = bm::bit_block_insert(blk0, 31, 0);
        assert(co == 1);
        assert(blk0[0]==(~0u >> 1));
        assert(blk0[1]==3);
        assert(blk0[bm::set_block_size-1] == (~0u << 3));

        blk0[0] = 0u;
        co = bm::bit_block_insert(blk0, 31, 1);
        assert(co == 1);
        assert(blk0[0]==(1u << 31));
        assert(blk0[1]==(3u << 1));
        assert(blk0[bm::set_block_size-1] == (~0u << 4));
    }
    
    cout << "bit-insert stress 0..." << endl;
    {
        for (i = 0; i < bm::set_block_size; ++i)
            blk0[i] = blk1[i] = 0;
        
        blk0[0] = 1;
        for (i = 0; i < 65536; ++i)
        {
            bm::word_t co = bm::bit_block_insert(blk0, i, 0);
            assert(co == 0 || i == 65535);
            if (i < 65535)
            {
                unsigned t = bm::test_bit(blk0, i+1);
                assert(t);
            }
            for (unsigned k = 0; k < 65536; ++k)
            {
                if (k != i+1)
                {
                    unsigned t = bm::test_bit(blk0, k);
                    assert(!t);
                }
            } // for k
        } // for i
        for (i = 0; i < bm::set_block_size; ++i)
        {
            if (blk0[i] != blk1[i])
            {
                cerr << "Stress insert(0) failed" << endl;
                exit(1);
            }
        }
    }
    cout << "OK" << endl;

    cout << "bit-insert stress 1..." << endl;
    {
        for (i = 0; i < bm::set_block_size; ++i)
            blk0[i] = ~0u;
        
        for (i = 0; i < 65536; ++i)
        {
            bm::word_t co = bm::bit_block_insert(blk0, i, 0);
            assert(co == 1 || i == 65535);
            for (unsigned k = 0; k < 65536; ++k)
            {
                unsigned t = bm::test_bit(blk0, k);
                if (k <= i)
                {
                    assert(!t);
                }
                else
                {
                    assert(t);
                }
            } // for k
        } // for i
        for (i = 0; i < bm::set_block_size; ++i)
        {
            if (blk0[i] != blk1[i])
            {
                cerr << "Stress insert(1) failed" << endl;
                exit(1);
            }
        }
    }
    cout << "OK" << endl;

    
    cout << "---------------------------- BlockBitInsertTest test OK" << endl;
}


static
void BlockBitEraseTest()
{
    cout << "---------------------------- BlockBitEraseTest test" << endl;
    bm::word_t BM_VECT_ALIGN blk0[bm::set_block_size] BM_VECT_ALIGN_ATTR = { 0 };
    bm::word_t BM_VECT_ALIGN blk1[bm::set_block_size] BM_VECT_ALIGN_ATTR = { 0 };

    bm::bit_block_set(blk0, 0);
    bm::bit_block_set(blk1, 0);
    
    {
        blk0[0] = 1;
        bm::bit_block_erase(blk0, 0, true);
        assert(blk0[0] == 0);
        assert(blk0[bm::set_block_size-1] == (1u << 31u));
        
        blk0[0] = ~0u;
        blk0[1] = 1u;
        blk0[bm::set_block_size-1] = (1u << 31u);
        bm::bit_block_erase(blk0, 0, false);
        assert(blk0[0] == ~0u);
        assert(blk0[1] == 0);
        assert(blk0[bm::set_block_size-1] == (1u << 30u));
        
        bm::bit_block_set(blk0, 0);
        blk0[1] = 15u; // ..01111
        bm::bit_block_erase(blk0, 31, false);
        assert(blk0[1] == 7u); // ..0111
        assert(blk0[0] == 1u << 31u);

        bm::bit_block_erase(blk0, 32, false);
        assert(blk0[1] == 3u); // ..011

        blk0[1] = 15u; // ..01111
        bm::bit_block_erase(blk0, 33, false);
        assert(blk0[1] == 0b111);
    }
    
    {
        bm::bit_block_set(blk0, ~0u);
        bm::word_t acc;
        bm::bit_block_shift_l1_unr(blk0, &acc, true);
        for (unsigned i = 0; i < bm::set_block_size; ++i)
        {
            assert(blk0[i] == ~0u);
        }
        bm::bit_block_erase(blk0, 0, false);
        for (unsigned i = 0; i < bm::set_block_size-1; ++i)
        {
            assert(blk0[i] == ~0u);
        }
        assert(blk0[bm::set_block_size-1] == ((~0u) >> 1) );

        bm::bit_block_shift_l1_unr(blk0, &acc, false);
        for (unsigned i = 0; i < bm::set_block_size-1; ++i)
        {
            assert(blk0[i] == ~0u);
        }
        assert(blk0[bm::set_block_size-1] == ((~0u) >> 2) );
        
        bm::bit_block_erase(blk0, 0, true);
        for (unsigned i = 0; i < bm::set_block_size-1; ++i)
        {
            assert(blk0[i] == ~0u);
        }
        assert(blk0[bm::set_block_size-1] == ((~0u >> 3) |  (1u << 31)));
    }

    cout << "bit-insert-erase stress 0..." << endl;
    unsigned c = 0;
    {
    bm::bit_block_set(blk0, 0);
    bm::bit_block_set(blk1, 0);
    
    blk0[bm::set_block_size-1] = (1u << 31);
    for (unsigned i = 65535; i != 0; --i)
    {
        unsigned t = bm::test_bit(blk0, i);
        assert(t);
        bm::bit_block_erase(blk0, 0, false);
        unsigned cnt = bm::bit_block_count(blk0);
        c += cnt;
        if (cnt != 1)
        {
            unsigned control = 0;
            for (unsigned k = 0; k < 65536; ++k)
            {
                t = bm::test_bit(blk0, k);
                control += t;
            }

            t = bm::test_bit(blk0, i-1);
            assert(t);
            cerr << t << " " << control << endl;
            cerr << "i=" << i << " cnt=" << cnt;
            cerr << " CNT==1 failed!" << endl;
            assert(cnt == 1);
            exit(1);
        }
    } // for i

    std::cout << c << endl;

    bm::bit_block_set(blk0, 0);

    {
    blk0[bm::set_block_size-1] = (1u << 31);
    unsigned j = 0;
    for (unsigned i = 65535; i != 0; --i, ++j)
    {
        unsigned t = bm::test_bit(blk0, i);
        assert(t);
        bm::bit_block_erase(blk0, j, false);
        if (i <= j)
        {
            t = bm::test_bit(blk0, j);
            assert(!t);
            t = bm::test_bit(blk0, i);
            assert(t);

            auto cnt = bm::bit_block_count(blk0);
            assert(cnt);
            
            bm::bit_block_erase(blk0, i, false);
            t = bm::test_bit(blk0, i);
            assert(!t);
            cnt = bm::bit_block_count(blk0);
            assert(!cnt);
            break;
        }
        auto cnt = bm::bit_block_count(blk0);
        assert(cnt == 1);
    } // for i
    }
    
    {
        bm::bit_block_set(blk0, 0);
        for (unsigned i = 0; i < 65535; ++i)
        {
            for(unsigned j = i; j < 65535; ++j)
            {
                bm::bit_block_set(blk0, 0);
                unsigned bitcount = j - i + 1;
                assert(i + bitcount < 65536);
                bm::or_bit_block(blk0, i, bitcount);
                if (bitcount == 1)
                {
                    break;
                }
                unsigned bc = bm::bit_block_count(blk0);
                for (unsigned k = i + bitcount/2; k < i+bitcount; ++k)
                {
                    auto t = bm::test_bit(blk0, k);
                    assert(t);

                    bm::bit_block_erase(blk0, k, false);
                    
                    unsigned bc2 = bm::bit_block_count(blk0);
                    assert(bc2 == bc-1);
                    --bc;
                } // for k
                
            } // for j
        } // for i
    }


    }
    cout << "ok" << endl;

    cout << "---------------------------- BlockBitEraseTest test OK" << endl;
}

static
void EmptyBVTest()
{
    cout << "---------------------------- Empty bvector test" << endl;

    {
        bvect bv1;
        bvect bv2;
        
        bvect bv3(bv1 & bv2);
        bvect bv4 = (bv1 & bv2);
        
        std::vector< bvect > v;
        v.push_back(bvect());
    }

    {
        bvect  bv1;
        bm::id_t cnt = bv1.count_range(0, 10);
        if (cnt)
        {
            cerr << "Failed count_range()" << endl;
            exit(1);
        }
        bool b = bv1.test(0);
        if (b)
        {
            cerr << "Failed test" << endl;
            exit(1);
        }
        
        b = bv1.any();
        if (b)
        {
            cerr << "Failed any" << endl;
            exit(1);
        }
        
        bv1.set_bit(0);
        if (!bv1.any())
        {
            cerr << "Failed set_bit" << endl;
            exit(1);
        }
    }
    {
        bvect  bv1;
        bool b = bv1.set_bit_and(0, false);
        if (bv1.any() || b)
        {
            cerr << "Failed set_bit" << endl;
            exit(1);
        }
    }
    {
        bvect  bv1;
        bv1.set_range(0, 1, false);
        if (bv1.any())
        {
            cerr << "Failed set_range" << endl;
            exit(1);
        }
        bv1.set_range(0, 1, true);
        if (bv1.count()!=2)
        {
            cerr << "Failed set_range(0,1)" << endl;
            exit(1);
        }
    }
    {
        bvect  bv1;
        bv1.clear_bit(0);
        if (bv1.any())
        {
            cerr << "Failed clear_bit" << endl;
            exit(1);
        }
    }
    {
        bvect  bv1;
        bv1.clear();
        if (bv1.any())
        {
            cerr << "Failed clear()" << endl;
            exit(1);
        }
    }
    {
        bvect  bv1;
        bv1.invert();
        if (!bv1.any())
        {
            cerr << "Failed invert()" << endl;
            exit(1);
        }
    }
    {
        bvect  bv1;
        bvect  bv2;
        bv1.swap(bv2);
        if (bv1.any() || bv2.any())
        {
            cerr << "Failed swap()" << endl;
            exit(1);
        }
    }
    {
        bvect  bv1;
        if (bv1.get_first() != 0)
        {
            cerr << "Failed get_first()" << endl;
            exit(1);
        }
    }
    {
        bvect  bv1;
        if (bv1.extract_next(0) != 0)
        {
            cerr << "Failed extract_next()" << endl;
            exit(1);
        }
    }
    {
        bvect  bv1;
        bvect::statistics st;
        bv1.calc_stat(&st);
        if (st.memory_used == 0)
        {
            cerr << "Failed calc_stat()" << endl;
            exit(1);
        }
    }
    {
        bvect  bv1;
        bvect  bv2;
        bvect  bv3;
        bv1.bit_or(bv2);
        if (bv1.any())
        {
            cerr << "Failed bit_or()" << endl;
            exit(1);
        }
        bv2.set_bit(100000000);
        bv1.bit_or(bv2);
        if (!bv1.any())
        {
            cerr << "Failed bit_or()" << endl;
            exit(1);
        }
        bv1.bit_or(bv3);
        if (bv1.count()!=1)
        {
            cerr << "Failed bit_or()" << endl;
            exit(1);
        }
    }
    
    {
        bvect  bv1;
        bvect  bv2;
        bv1.bit_and(bv2);
        if (bv1.any())
        {
            cerr << "Failed bit_and()" << endl;
            exit(1);
        }
        bv2.set_bit(100000000);
        bv1.bit_and(bv2);
        if (bv1.any())
        {
            cerr << "Failed bit_and()" << endl;
            exit(1);
        }
        bv2.bit_and(bv1);
        if (bv2.count()!=0)
        {
            cerr << "Failed bit_and()" << endl;
            exit(1);
        }
    }
    
    {
        bvect  bv1;
        bvect::statistics st1;
        bv1.optimize(0, bvect::opt_compress, &st1);
        if (st1.memory_used == 0)
        {
            cerr << "Failed calc_stat()" << endl;
            exit(1);
        }
        bv1.optimize_gap_size();
    }
    
    {
        bvect  bv1;
        bvect  bv2;
        
        int r = bv1.compare(bv2);
        if (r != 0)
        {
            cerr << "Failed compare()" << endl;
            exit(1);
            
        }
        bv2.set_bit(1000);
        r = bv1.compare(bv2);
        if (r == 0)
        {
            cerr << "Failed compare()" << endl;
            exit(1);
            
        }
        r = bv2.compare(bv1);
        if (r == 0)
        {
            cerr << "Failed compare()" << endl;
            exit(1);
        }
    }
    
    {
        bvect  bv1;
        bvect::enumerator en1 = bv1.first();
        bvect::enumerator en2 = bv1.get_enumerator(65535);
        CompareEnumerators(en1, en2);
        if (en1.valid() || en2.valid())
        {
            cerr << "failed first enumerator" << endl;
            exit(1);
        }
    }

    {
        bvect bv1;
        bv1.resize(0);

        bv1.invert();
        assert(!bv1.any());
        assert(bv1.size()==0);
    }
    
    cout << "---------------------------- Empty bvector test OK" << endl;
    
}




static
void BasicFunctionalityTest()
{
    cout << "---------------------------- Basic functinality test" << endl;

    assert(ITERATIONS < BITVECT_SIZE);

    bvect_mini     bvect_min(BITVECT_SIZE);
    bvect          bvect_full;
    bvect          bvect_full1;
    bvect::rs_index_type bc_arr;
    bvect::rs_index_type bc_arr1;

    printf("\nBasic functionality test.\n");
    
    {
        bvect::rs_index_type rs_idx;
        for (unsigned i = 0; i < ITERATIONS; ++i)
        {
            rs_idx.resize(i);
        }
    }

    // filling vectors with regular values
    
    cout << "test data generation... " << endl;
    unsigned i;
    for (i = 0; i < ITERATIONS; ++i)
    {
        bvect_min.set_bit(i);
        bvect_full.set_bit(i);

        bvect_full.build_rs_index(&bc_arr);
        
        bm::id_t pos1, pos2, pos3, pos4;
        auto rf1 = FindRank(bvect_full, i+1, 0, pos1);
        auto rf2 = bvect_full.find_rank(i+1, 0, pos2);
        auto rf3 = bvect_full.find_rank(i+1, 0, pos3);
        auto rf4 = bvect_full.select(i+1, pos4, bc_arr);

        assert(rf1);
        assert(rf2);
        assert(rf3);
        assert(rf4);
        if (pos1 != pos2 || pos1 != pos3 || pos1 != pos4)
        {
            rf2 = bvect_full.find_rank(i+1, 0, pos2);
            cerr << "1.Rank check error!\n"
                 << " pos1 = " << pos1
                 << " pos2 = " << pos2
                 << " pos3 = " << pos3
                 << " pos4 = " << pos4
                 << endl;
                 ;
            exit(1);
        }

    }

    bvect_full1.set_range(0, ITERATIONS-1);

    bvect_full1.build_rs_index(&bc_arr1);

    cout << "Rank check 2" << endl;

    for (i = 0; i < ITERATIONS; ++i)
    {
        bm::id_t pos1, pos2, pos3, pos4;
        auto rf1 = FindRank(bvect_full1, i+1, 0, pos1);
        auto rf2 = bvect_full1.find_rank(i+1, 0, pos2);
        auto rf3 = bvect_full1.find_rank(i+1, 0, pos3);
        auto rf4 = bvect_full1.select(i+1, pos4, bc_arr1);
        assert(rf1);
        assert(rf2);
        assert(rf3);
        assert(rf4);
        if (pos1 != pos2 || pos1 != pos3 || pos1 != pos4)
        {
            rf2 = bvect_full1.find_rank(i+1, 0, pos2);
            cerr << "2.Rank check error!\n"
                 << " pos1 = " << pos1
                 << " pos2 = " << pos2
                 << " pos3 = " << pos3
                 << " pos4 = " << pos4
                 << endl;
                 ;
            exit(1);
        }
        if (i % 1000 == 0)
        {
            cout << "\r" << i << " / " << ITERATIONS << flush;
        }
    }
    cout << endl;

    CheckCountRange(bvect_full, 0, ITERATIONS);
    CheckCountRange(bvect_full, 10, ITERATIONS+10);
    CheckCountRange(bvect_full1, 0, ITERATIONS);
    CheckCountRange(bvect_full1, ITERATIONS-10, ITERATIONS+10);
    CheckCountRange(bvect_full1, 10, ITERATIONS+10);

    if (bvect_full1 != bvect_full)
    {
        cout << "set_range failed!" << endl;
        print_stat(bvect_full1);
        exit(1);
    }

    print_stat(bvect_full);
    print_stat(bvect_full1);

    // checking the results
    unsigned count_min = 0;

    for (i = 0; i < ITERATIONS; ++i)
    {
        if (bvect_min.is_bit_true(i))
            ++count_min;
    }

    cout << "Rank check 3" << endl;
    for (i = 0; i < ITERATIONS; ++i)
    {
        CheckCountRange(bvect_full, i, ITERATIONS);
        CheckCountRange(bvect_full1, i, ITERATIONS);
    }

    cout << "Rank check 4" << endl;
    for (i = ITERATIONS; i > 0; --i)
    {
        CheckCountRange(bvect_full, 0, i);
        CheckCountRange(bvect_full1, 0, i);
    }
    
    unsigned count_full = bvect_full.count();

    if (count_min == count_full)
    {
        printf("simple count test ok.\n");
    }
    else
    {
        printf("simple count test failed count_min = %i  count_full = %i\n", 
               count_min, count_full);
        exit(1);
    }


    // detailed vectors verification

    CheckVectors(bvect_min, bvect_full, ITERATIONS);

    // now clearning

    for (i = 0; i < ITERATIONS; i+=2)
    {
        bvect_min.clear_bit(i);
        bvect_full.clear_bit(i);
        bvect_full1.set_range(i, i, false);
    }

    CheckVectors(bvect_min, bvect_full, ITERATIONS);
    CheckVectors(bvect_min, bvect_full1, ITERATIONS);

    for (i = 0; i < ITERATIONS; ++i)
    {
        bvect_min.clear_bit(i);
    }
    bvect_full.clear();

    CheckVectors(bvect_min, bvect_full, ITERATIONS);

    cout << "Random step filling" << endl;

    for (i = rand()%10; i < ITERATIONS; i+=rand()%10)
    {
        bvect_min.clear_bit(i);
        bvect_full.clear_bit(i);
    }
    
    CheckVectors(bvect_min, bvect_full, ITERATIONS);

    bvect bv1;
    bvect bv2;

    bv1[10] = true;
    bv1[1000] = true;

    bv2[200] = bv2[700] = bv2[500] = true;

    bv1.swap(bv2);

    if (bv1.count() != 3)
    {
        cout << "Swap test failed!" << endl;
        exit(1);
    }

    if (bv2.count() != 2)
    {
        cout << "Swap test failed!" << endl;
        exit(1);
    }

    {
        //bm::standard_alloc_pool pool;
        bvect::allocator_pool_type pool;
        bvect bv3, bv4;
        bv3.set_allocator_pool(&pool);
        bv3.set(10, true);
        bv4.set(10, true);
        bv4.set(10, false);
        bv3 &= bv4;
    }
    {
        bvect bv(100);
        bv.set_range(150, 151);
        assert(bv.size() == 152);
        bv.set_bit(160);
        assert(bv.size() == 161);
    }
}


static
void generate_test_vectors(std::vector<bm::id_t> &v1,
                           std::vector<bm::id_t> &v2,
                           std::vector<bm::id_t> &v3,
                           unsigned vector_max)
{
    bm::id_t j;
    for (j = 0; j < vector_max; j += 2)
        v1.push_back(j);
    for (j = 0; j < vector_max; j += 5)
        v2.push_back(j);
    for (j = 0; j < vector_max; j += 120)
        v3.push_back(j);
}


static
void BvectorBulkSetTest()
{
    cout << "---------------------------- Bvector BULK set test" << endl;

    
    {
        unsigned ids[] = { 0 };
        
        bvect bv1, bv2, bv3;
        {
        bvect::bulk_insert_iterator iit = bv3.inserter();
        for (unsigned i = 0; i < sizeof(ids)/sizeof(ids[0]); ++i)
        {
            bv1.set(ids[i]);
            iit = ids[i];
        }
        }
        bv2.set(&ids[0], sizeof(ids)/sizeof(ids[0]));
        
        int cmp = bv1.compare(bv2);
        assert(cmp==0);
        cmp = bv1.compare(bv3);
        assert(cmp==0);
        
        bv2.set(0);
        bv2.keep(&ids[0], sizeof(ids)/sizeof(ids[0]));
        assert(bv2.count()==1);
        assert(bv2.test(0));
    }

    {
        unsigned ids[] = {65535, bm::id_max };
        unsigned cnt;
        bvect bv2;
        bv2.set(&ids[0], sizeof(ids)/sizeof(ids[0]));
        
        cnt = bv2.count();
        cout << cnt << endl;

        assert(cnt == 1);
        assert(bv2.test(ids[0]));
    }

    // set bits in FULL vector
    {
        unsigned ids[] = {0, 10, 65535, bm::id_max-1, bm::id_max };
        unsigned cnt;
        bvect bv2;
        bv2.invert();
        struct bvect::statistics st1, st2;
        bv2.calc_stat(&st1);

        bv2.set(&ids[0], sizeof(ids)/sizeof(ids[0]));

        bv2.calc_stat(&st2);
        assert(st1.bit_blocks == st2.bit_blocks);
        assert(st1.gap_blocks == st2.gap_blocks);

        cnt = bv2.count();
        cout << cnt << endl;

        assert(cnt == bm::id_max);
    }

    // test correct sizing
    {
        unsigned ids[] = {65536 };
        bvect bv1;
        bv1.resize(10);

        bv1.set(&ids[0], sizeof(ids)/sizeof(ids[0]));
        assert(bv1.size()==65536+1);
        bv1.keep(&ids[0], sizeof(ids)/sizeof(ids[0]));
        cout << bv1.size() << endl;
        assert(bv1.size()==65536+1);
    }

    {
        unsigned ids[] = {65536, 1280000, 65535 };
        bvect bv1, bv2;

        for (unsigned i = 0; i < sizeof(ids)/sizeof(ids[0]); ++i)
            bv1.set(ids[i]);
        
        bv2.set(&ids[0], sizeof(ids)/sizeof(ids[0]));
        int cmp = bv1.compare(bv2);
        assert(cmp==0);
        
        bv2.keep(&ids[0], sizeof(ids)/sizeof(ids[0]));
        cmp = bv1.compare(bv2);
        assert(cmp==0);
    }


    {
    unsigned ids[] = { 0, 1, 2, 3, 4, 5, 256, 1024, 1028, 256000 };
    
    bvect bv1, bv2;
    for (unsigned i = 0; i < sizeof(ids)/sizeof(ids[0]); ++i)
        bv1.set(ids[i]);
    bv2.set(&ids[0], sizeof(ids)/sizeof(ids[0]));
    
    DetailedCompareBVectors(bv1, bv2);

    int cmp = bv1.compare(bv2);
    assert(cmp==0);
    
    {
        bvect bv_inv;
        bv_inv.flip();
        unsigned keep_cnt = sizeof(ids)/sizeof(ids[0]);
        bv_inv.keep(&ids[0], keep_cnt, bm::BM_SORTED);
        unsigned cnt_inv2 = bv_inv.count();
        assert(cnt_inv2 == keep_cnt);
        for (unsigned i = 0; i < sizeof(ids)/sizeof(ids[0]); ++i)
        {
            assert(bv_inv.test(ids[i]));
        }
    }
    
    bvect bv3, bv4, bv5;
    bv3.invert();
    bv4.invert();
    bv5.invert();

    {
        bvect::bulk_insert_iterator iit = bv5.inserter();
        for (unsigned i = 0; i < sizeof(ids)/sizeof(ids[0]); ++i)
        {
            bv3.set(ids[i]);
            iit = ids[i];
        }
    }
    bv4.set(&ids[0], sizeof(ids)/sizeof(ids[0]));
    cmp = bv3.compare(bv4);
    assert(cmp==0);
    cmp = bv4.compare(bv3);
    assert(cmp==0);

    cmp = bv3.compare(bv5);
    assert(cmp==0);
    }
    
    {
    unsigned vector_max = 4000000;
    std::vector<bm::id_t> v1, v2, v3;
    generate_test_vectors(v1, v2, v3, vector_max);

    for (unsigned k = 0; k < 2; ++ k)
    {
        bvect bvu, bvuc;
        bvect bv1, bv2, bv3, bv11;
        bvect bv1c, bv2c, bv3c;
        
        {
        bvect::bulk_insert_iterator iit(bv11);
        for (unsigned i = 0; i < v1.size(); ++i)
        {
            bv1c.set(v1[i]);
            iit = v1[i];
        }
        iit.flush();
        }
        
        for (unsigned i = 0; i < v2.size(); ++i)
            bv2c.set(v2[i]);
        for (unsigned i = 0; i < v3.size(); ++i)
            bv3c.set(v3[i]);
        
        // union of 3 vectors
        bvuc = bv1c;
        bvuc |= bv2c;
        bvuc |= bv3c;

        bv1.set(&v1[0], unsigned(v1.size()));
        bv2.set(&v2[0], unsigned(v2.size()));
        bv3.set(&v3[0], unsigned(v3.size()));

        // imported union of 3 vectors
        bvu.set(&v1[0], unsigned(v1.size()));
        bvu.set(&v2[0], unsigned(v2.size()));
        bvu.set(&v3[0], unsigned(v3.size()));

        cout << bv1.count() << " " << bv1c.count() << endl;
        
        int cmp;
        cmp = bv1c.compare(bv1);
        if (cmp != 0)
        {
            DetailedCompareBVectors(bv1, bv1c);
        }
        assert(cmp==0);
        cmp = bv2c.compare(bv2);
        assert(cmp==0);
        cmp = bv3c.compare(bv3);
        assert(cmp==0);
        cmp = bv1.compare(bv11);
        assert(cmp==0);

        cmp = bvuc.compare(bvu);
        assert(cmp == 0);
        
        {
            std::random_device rd;
            std::mt19937 g(rd());
            
            std::shuffle(v1.begin(), v1.end(), g);
            std::shuffle(v2.begin(), v2.end(), g);
            std::shuffle(v3.begin(), v3.end(), g);
        }
    }
    
    
    }
    
    cout << "Bulk bvector<>::set() stress.." << endl;
    {
        unsigned vector_max =40000000;
        unsigned delta_max = 65537;

        bvect bv1, bv2;
        bvect bv1c;
        std::vector<bm::id_t> v1;

        for (unsigned delta = 1; delta < delta_max; ++delta)
        {
            v1.resize(0);
            bvect::bulk_insert_iterator iit(bv2);
            for (unsigned i = 0; i < vector_max; i+=delta)
            {
                v1.push_back(i);
                iit = i;
            }
            iit.flush();
            
            bv1.set(&v1[0], unsigned(v1.size()));
            bm::combine_or(bv1c, v1.begin(), v1.end());
            
            int cmp = bv1.compare(bv1c);
            if (cmp!=0)
            {
                cerr << "1.Failed bulk set test at delta=" << delta << endl;
                exit(1);
            }
            cmp = bv1.compare(bv2);
            if (cmp!=0)
            {
                cerr << "2.Failed bulk set test at delta=" << delta << endl;
                exit(1);
            }
            
            // test AND/keep
            {
                bvect bv3(bv1);
                bvect bv4;
                bv3.keep(&v1[0], unsigned(v1.size()));
                bv4.set(&v1[0], unsigned(v1.size()));
                cmp = bv3.compare(bv4);
                if (cmp!=0)
                {
                    cerr << "3.Failed keep() test at delta=" << delta << endl;
                    exit(1);
                }
                if (v1.size())
                {
                    bv3.keep(&v1[0], 1);
                    assert(bv3.count()==1);
                    assert(bv3.test(v1[0]));
                }
            }

            bv1.clear();
            bv1c.clear();
            bv2.clear();

            if (delta % 256 == 0)
            {
                cout << "\r" << delta << "/" << delta_max << flush;
            }
            
        } // for delta
        cout << endl;
    }
    
    
    cout << "---------------------------- Bvector BULK set test OK" << endl;
}



static
void RankFindTest()
{
    cout << "---------------------------- Find Rank test" << endl;
    
    {
    bvect bv1;
    bv1[30] = true;
    bv1[65534] = true;

    bvect::rs_index_type bc_arr1;
    bv1.build_rs_index(&bc_arr1);

    bool rf1, rf2, rf3;
    bm::id_t pos, pos1;
    rf1 = bv1.find_rank(1, 20, pos);
    rf3 = bv1.find_rank(1, 20, pos1, bc_arr1);
    assert(rf1);
    assert(rf3);
    assert(pos == 30);
    assert(pos1 == 30);

    rf2 = bv1.find_rank(2, 30, pos);
    rf3 = bv1.find_rank(2, 30, pos1, bc_arr1);
    assert(rf2);
    assert(rf3);
    assert(pos == 65534);
    assert(pos1 == 65534);
    }
    
    cout << "Find Rank test stress 1\n" << endl;
    
    {
        const unsigned max_size = 2000000;
        bvect bv1;
        for (unsigned i = 0; i < max_size;)
        {
            bv1.set(i);
            i += rand()%5;
        }
        bvect::rs_index_type bc_arr1;
        bv1.build_rs_index(&bc_arr1);

        
        for (unsigned i = 0; i < max_size; ++i)
        {
            bool rf1, rf3;
            bm::id_t pos, pos1;
            
            rf1 = bv1.find_rank(0, i, pos);
            rf3 = bv1.find_rank(0, i, pos1, bc_arr1);
            assert (rf1 == rf3);
            if (rf1)
            {
                if (pos != pos1)
                {
                    cerr << "Rank cmp test failed! i=" << i
                         << " pos="  << pos
                         << " pos1=" << pos1
                         << endl;
                    exit(1);
                }
            }
            
            rf1 = bv1.find_rank(i, max_size-i, pos);
            rf3 = bv1.find_rank(i, max_size-i, pos1, bc_arr1);
            assert (rf1 == rf3);
            if (rf1)
            {
                if (pos != pos1)
                {
                    cerr << "Rank cmp test failed! i=" << i
                         << " pos="  << pos
                         << " pos1=" << pos1
                         << endl;
                    exit(1);
                }
            }
            if (i % 100 == 0)
                cout << "\r" << i << "/" << max_size << flush;
        } // for
        cout << endl;
    }
    
    cout << "---------------------------- Find Rank test OK" << endl;
}

static
void BvectorIncTest()
{
    cout << "---------------------------- Bvector inc test" << endl;
    
    {
    bvect bv1;
    bool b;
    
    b = bv1.inc(0);
    assert(!b);
    b = bv1.inc(0);
    assert(b);
    
    b = bv1.inc(10);
    assert(!b);
    b = bv1.inc(10);
    assert(b);
    }
    
    {
    bvect bv1(BM_GAP);
    bool b;
    
    assert(bv1.count()==0);
    
    b = bv1.inc(0);
    assert(!b);
    cout << bv1.count() << endl;
    assert(bv1.count()==1);
    b = bv1.inc(0);
    assert(b);
    assert(bv1.count()==0);

    b = bv1.inc(10);
    assert(!b);
    b = bv1.inc(10);
    assert(b);
    }

    {
    bvect bv1(BM_GAP);
    bool b;

    bv1.flip();
    
    b = bv1.inc(0);
    assert(b);
    b = bv1.inc(0);
    assert(!b);
    
    b = bv1.inc(10);
    assert(b);
    b = bv1.inc(10);
    assert(!b);
    }

    cout << "---------------------------- Bvector inc test OK" << endl;
}

// -----------------------------------------------------------------------

static
void optimize_fill(bvect& bv, bvect::size_type base, unsigned inc,
                   bvect::size_type max_bits = bm::gap_max_bits,
                   bool value = true)
{
    for (bvect::size_type i = 0; i < bm::set_sub_array_size; ++i)
    {
        bvect::size_type base_idx = i * bm::gap_max_bits;
        for (bvect::size_type j = base; j < max_bits; j += inc)
        {
            bv.set(base_idx + j, value);
        } // for j
    } // for i
}

static
void OptimizeTest()
{
    cout << "---------------------------- Bvector Optimize test" << endl;
    BM_DECLARE_TEMP_BLOCK(tb)

    {
        bvect bv;
        optimize_fill(bv, 0, 1, bm::gap_max_bits, true);
        
        bvect::statistics st1;
        bv.calc_stat(&st1);
        
        assert(st1.bit_blocks == bm::set_sub_array_size);
        assert(st1.gap_blocks == 0);
        assert(st1.ptr_sub_blocks == 1);
        
        bv.optimize(tb, bvect::opt_compress, &st1);

        assert(st1.bit_blocks == 0);
        assert(st1.gap_blocks == 0);
        assert(st1.ptr_sub_blocks == 0);
        
        bv.calc_stat(&st1);
        
        assert(st1.bit_blocks == 0);
        assert(st1.gap_blocks == 0);
        assert(st1.ptr_sub_blocks == 0);
    }
    
    {
        bvect bv;
        optimize_fill(bv, 0, 100, bm::gap_max_bits, true);
        optimize_fill(bv, 0, 100, bm::gap_max_bits, false);

        bvect::statistics st1;
        bv.calc_stat(&st1);
        
        assert(st1.bit_blocks == bm::set_sub_array_size);
        assert(st1.gap_blocks == 0);
        assert(st1.ptr_sub_blocks == 1);
        
        bv.optimize(tb, bvect::opt_compress, &st1);

        assert(st1.bit_blocks == 0);
        assert(st1.gap_blocks == 0);
        assert(st1.ptr_sub_blocks == 0);
    }


    {
        bvect bv(BM_GAP);
        optimize_fill(bv, 0, 1, bm::gap_max_bits, true);
        
        bvect::statistics st1;
        bv.calc_stat(&st1);
        
        assert(st1.bit_blocks == 0);
        assert(st1.gap_blocks == bm::set_sub_array_size);
        assert(st1.ptr_sub_blocks == 1);
        
        bv.optimize(tb, bvect::opt_compress, &st1);

        assert(st1.bit_blocks == 0);
        assert(st1.gap_blocks == 0);
        assert(st1.ptr_sub_blocks == 0);
        
        bv.calc_stat(&st1);
        
        assert(st1.bit_blocks == 0);
        assert(st1.gap_blocks == 0);
        assert(st1.ptr_sub_blocks == 0);
    }
    
    {
        bvect bv(BM_GAP);
        optimize_fill(bv, 0, 1000, bm::gap_max_bits, true);
        optimize_fill(bv, 0, 1000, bm::gap_max_bits, false);

        bvect::statistics st1;
        bv.calc_stat(&st1);
        
        assert(st1.bit_blocks == 0);
        assert(st1.gap_blocks == bm::set_sub_array_size);
        assert(st1.gaps_by_level[0] == 0);
        assert(st1.gaps_by_level[1] == bm::set_sub_array_size);
        assert(st1.ptr_sub_blocks == 1);
        
        bv.optimize(tb, bvect::opt_compress, &st1);

        assert(st1.bit_blocks == 0);
        assert(st1.gap_blocks == 0);
        assert(st1.ptr_sub_blocks == 0);
    }
    
    {
        bvect bv(BM_GAP);
        optimize_fill(bv, 0, 1000, bm::gap_max_bits, true);
        optimize_fill(bv, 1, 1, bm::gap_max_bits, false);

        bvect::statistics st1;
        bv.calc_stat(&st1);
        
        assert(st1.bit_blocks == 0);
        assert(st1.gap_blocks == bm::set_sub_array_size);
        assert(st1.gaps_by_level[0] == 0);
        assert(st1.gaps_by_level[1] == bm::set_sub_array_size);
        assert(st1.ptr_sub_blocks == 1);
        
        bv.optimize(tb, bvect::opt_compress, &st1);

        assert(st1.bit_blocks == 0);
        assert(st1.gap_blocks == bm::set_sub_array_size);
        assert(st1.ptr_sub_blocks == 1);
        assert(st1.gaps_by_level[0] == bm::set_sub_array_size);
        assert(st1.gaps_by_level[1] == 0);
    }


    cout << "---------------------------- Bvector Optimize test OK" << endl;
}


// -----------------------------------------------------------------------

static
void generate_sparse_bvector(bvect& bv,
                             unsigned min = 0,
                             unsigned max = 40000000,
                             unsigned fill_factor = 65536)
{
    bvect::bulk_insert_iterator iit(bv);
    unsigned ff = fill_factor / 10;
    for (unsigned i = min; i < max; i+= ff)
    {
        //bv.set(i);
        iit = i;
        ff += ff / 2;
        if (ff > fill_factor)
            ff = fill_factor / 10;
    }
    iit.flush();
}


static
void GenerateShiftTestCollection(std::vector<bvect>* target,
                            unsigned count = 30,
                            unsigned vector_max = 40000000,
                            bool optimize = true)
{
    assert(target);
    bvect bv_common; // sub-vector common for all collection
    generate_sparse_bvector(bv_common, vector_max/10, vector_max, 250000);
    
    unsigned cnt1 = (count / 2);
    
    unsigned i = 0;
    
    for (i = 0; i < cnt1; ++i)
    {
        std::unique_ptr<bvect> bv (new bvect);
        generate_bvector(*bv, vector_max, optimize);
        *bv |= bv_common;
        if (optimize)
            bv->optimize();
        target->push_back(std::move(*bv));
    } // for
    
    unsigned fill_factor = 10;
    for (; i < count; ++i)
    {
        std::unique_ptr<bvect> bv (new bvect);
        
        FillSetsIntervals(0, *bv, vector_max/ 10, vector_max, fill_factor);
        *bv |= bv_common;

        target->push_back(std::move(*bv));
    } // for
}



static
void BvectorShiftTest()
{
    cout << "---------------------------- Bvector SHIFT test" << endl;


    {
    bvect bv;
    
    bv.set(bm::id_max-1);
    bv.shift_right();
    print_bv(bv);
    assert(bv.count()==0);
    }
    
    {
    bvect bv(BM_GAP);
    
    bv.set(bm::id_max-1);
    bv.optimize();
    bv.shift_right();
    assert(bv.count()==0);
    }
    
    {
        bvect bv;
        
        bv.set();
        auto cnt1 = bv.count();
        bv.shift_right();
        auto cnt2 = bv.count();
        assert(cnt1-1 == cnt2);
        bool b = bv.test(0);
        assert(!b);
    }
    {
        bvect bv;
        bv.set();
        bv.set(bm::gap_max_bits * bm::set_sub_array_size, false);
        
        auto cnt1 = bv.count();
        bv.shift_left();
        auto cnt2 = bv.count();
        assert(cnt1-1 == cnt2);
        bool b = bv.test(bm::id_max-1);
        assert(!b);
        b = bv.test(bm::id_max-2);
        assert(b);
        b = bv.test(bm::gap_max_bits * bm::set_sub_array_size);
        assert(b);
        b = bv.test(bm::gap_max_bits * bm::set_sub_array_size - 1);
        assert(!b);
    }
    
    {
        bvect bv;
        
        bv.set();
        auto cnt1 = bv.count();
        bv.shift_left();
        auto cnt2 = bv.count();
        assert(cnt1-1 == cnt2);
        bool b = bv.test(bm::id_max-1);
        assert(!b);
    }



    {
    bvect bv;

    bv.set(0);
    bv.set(65535);
    bv.set(bm::id_max-1);
    bvect bv1(bv);
    
    bvect bv2(bm::BM_GAP);
    bv2 = bv;
    bv2.optimize();

    ShiftRight(&bv, 1);
    assert(bv.count() == 2);
    assert(bv.test(1));
    assert(bv.test(65536));

    bv1.shift_right();
    print_bv(bv1);
    int cmp = bv.compare(bv1);
    assert(cmp == 0);
    
    bv2.shift_right();
    print_bv(bv2);
    cmp = bv.compare(bv2);
    assert(cmp == 0);
    struct bvect::statistics st;
    bv2.calc_stat(&st);
    assert(st.gap_blocks >= 2);

    }
    
    {
    bvect bv(BM_GAP);
    struct bvect::statistics st;
    
    for (unsigned i = 0; i < 65536; ++i)
    {
        bv.set(i);
    }
    bv.calc_stat(&st);
    assert(st.gap_blocks == 1);
    
    auto cnt = bv.count();
    bv.shift_right();
    assert(bv.test(0)==0);
    assert(bv.count() == cnt);
    
    bv.calc_stat(&st);
    auto bcnt = st.bit_blocks + st.gap_blocks;
    assert(bcnt == 2);
    assert(st.gap_blocks);
    for (unsigned i = 0+1; i < 65536+1; ++i)
    {
        assert(bv.test(i));
    }
    }


    {
    bvect bv;
    bv.invert();
    unsigned cnt = bv.count();
    bool carry_over = bv.shift_right();
    assert(carry_over);
    unsigned cnt1 = bv.count();
    assert(cnt1 == cnt - 1);
    assert(bv.test(0)==0);
    assert(bv.test(1)==1);

    struct bvect::statistics st;
    bv.calc_stat(&st);
    auto bcnt = st.bit_blocks + st.gap_blocks;
    assert(bcnt == 2);

    }
    
    {
    bvect bv;
    
    bv.set(0);
    bv.set(65535);
    bv.set(66000);
    bv.optimize();
    bvect bv1(bv);
    ShiftRight(&bv, 1);
    bv1.shift_right();
    int cmp = bv.compare(bv1);
    
    assert(cmp == 0);
    }


    {
    bvect bv { 1 };
    bool carry_over = bv.shift_left();
    assert(!carry_over);
    unsigned idx = bv.get_first();
    assert(idx == 0);
    carry_over = bv.shift_left();
    assert(carry_over);
    idx = bv.get_first();
    std::cout << idx << endl;
    assert(idx == 0);
    assert(bv.count()==0);
    }
    
    {
    bvect bv { 4278190080 };
    bv.shift_left();
    unsigned idx = bv.get_first();
    assert(idx == 4278190080-1);
    bv.shift_left();
    idx = bv.get_first();
    assert(idx == 4278190080-2);
    }
    
    {
    bvect bv { 4278190080 };
    bv.optimize();
    bv.shift_left();
    unsigned idx = bv.get_first();
    assert(idx == 4278190080-1);
    bv.shift_left();
    idx = bv.get_first();
    assert(idx == 4278190080-2);
    }

    {
    std::cout << "Shift-L stress (1 bit shift)..\n";
    unsigned start = bm::id_max-1;
    bvect bv;
    bv.set(start);

   struct bvect::statistics st;
   bv.calc_stat(&st);
   auto bcnt = st.bit_blocks + st.gap_blocks;
   assert(bcnt == 1);
   
    std::chrono::time_point<std::chrono::steady_clock> s;
    std::chrono::time_point<std::chrono::steady_clock> f;
    
    s = std::chrono::steady_clock::now();

    for( ; start; --start)
    {
        bool carry_over = bv.shift_left();
        if (carry_over)
        {
            cout << "CO at " << start << endl;
            assert(bv.count()==0);
            assert(start == bm::id_max-1);
            break;
        }
        /*
        unsigned idx = bv.get_first();
        if(idx != start-1)
        {
            cerr << bv.count() << endl;
            cerr << "Shift-L Failed at idx=" << idx << " != " << start << endl;
            exit(1);
        }
        */

        if ((start % (1024 * 1024)) == 0)
        {
            f = std::chrono::steady_clock::now();
            auto diff = f - s;
            auto d = std::chrono::duration <double, std::milli> (diff).count();
            cout << "\r" << start << " (" << d << ") " << flush;

            unsigned idx = bv.get_first();
            assert(idx == start-1);

            bv.calc_stat(&st);
            bcnt = st.bit_blocks + st.gap_blocks;
            assert(bcnt == 1);
            
            s = std::chrono::steady_clock::now();
        }
    }
    cout << "ok.\n";
    }

    {
    std::cout << "Shift-R stress (1 bit shift)..\n";
    unsigned start = 0;
    bvect bv;
    bv.set(start);

   struct bvect::statistics st;
   bv.calc_stat(&st);
   auto bcnt = st.bit_blocks + st.gap_blocks;
   assert(bcnt == 1);
   
    std::chrono::time_point<std::chrono::steady_clock> s;
    std::chrono::time_point<std::chrono::steady_clock> f;
    
    s = std::chrono::steady_clock::now();

    while(1)
    {
        bool carry_over = bv.shift_right();
        if (carry_over)
        {
            cout << "CO at " << start << endl;
            assert(bv.count()==0);
            assert(start == bm::id_max-1);
            break;
        }

        if ((start % (1024 * 1024)) == 0)
        {
            f = std::chrono::steady_clock::now();
            auto diff = f - s;
            auto d = std::chrono::duration <double, std::milli> (diff).count();
            cout << "\r" << start << " (" << d << ") " << flush;
            
            unsigned idx = bv.get_first();
            assert(idx-1 == start);
            
            bv.calc_stat(&st);
            bcnt = st.bit_blocks + st.gap_blocks;
            assert(bcnt == 1);
 
            s = std::chrono::steady_clock::now();
        }
        ++start;
    }
    cout << "ok.\n";
    }

    {
        std::cout << "Shift-R stress (large vector shift)..\n";
        bvect bv;
        generate_bvector(bv);
        bvect bv_control(bv);
        
        unsigned max_shifts = 10000;
        for (unsigned i = 0; i < max_shifts; ++i)
        {
            ShiftRight(&bv_control, 1);
            bv.shift_right();
            int cmp = bv.compare(bv_control);
            assert(cmp==0);
            if ((i % 16) == 0)
            {
                cout << "\r" << i << "/" << max_shifts << flush;
            }
        }
    }
    cout << "ok.\n";


    // stress test for shifting aggregator
    //
    cout << "Aggregator based SHIT-R tests..." << endl;
    {
        const unsigned int REPEATS = 300;

        bvect mask_bv; // mask vector
        mask_bv.init();
        generate_bvector(mask_bv, 75000000, false); // mask is shorter on both ends

        std::vector<bvect> bv_coll1;
        GenerateShiftTestCollection(&bv_coll1, 25, 80000000);
        
        {
            bm::aggregator<bvect> agg;
            agg.add(&mask_bv);
            for (unsigned k = 0; k < bv_coll1.size(); ++k)
            {
                agg.add(&bv_coll1[k]);
            }

            for (unsigned i = 0; i < REPEATS; ++i)
            {
                bvect bv1(mask_bv);
                for (unsigned k = 0; k < bv_coll1.size(); ++k)
                {
                    bv1.shift_right();
                    bv1 &= bv_coll1[k];
                } // for
                
                bvect bv2;
                agg.combine_shift_right_and(bv2);
                int cmp = bv1.compare(bv2);
                if (cmp != 0)
                {
                    cerr << "Shift-R compare failure!" << endl;
                    exit(1);
                }
            } // for
        }
    }


    cout << "---------------------------- Bvector SHIFT test OK" << endl;
}

static
void BvectorInsertTest()
{
    cout << "---------------------------- Bvector INSERT test" << endl;
    
    {
        bvect bv { 1, 2, 3 };
        bvect bv_c { 2, 3, 4 };
        bvect bv1(bv);
        BVectorInsert(&bv, 0, false);
        int cmp = bv.compare(bv_c);
        assert(cmp == 0);
        
        bv1.insert(0, false);
        cmp = bv1.compare(bv_c);
        assert(cmp == 0);
    }
    
    {
        bvect bv { 1, 2, 3 };
        bvect bv_c { 0, 2, 3, 4 };
        bvect bv1(bv);
        BVectorInsert(&bv, 0, true);
        int cmp = bv.compare(bv_c);
        assert(cmp == 0);

        bv1.insert(0, true);
        cmp = bv1.compare(bv_c);
        assert(cmp == 0);
    }
    
    {
        bvect bv { 1, 20, 65535 };
        bvect bv_c { 1, 20, 65535, 65536 };
        bvect bv1(bv);
        BVectorInsert(&bv, 65535, true);
        int cmp = bv.compare(bv_c);
        assert(cmp == 0);

        bv1.insert(65535, true);
        print_bv(bv1);
        cmp = bv1.compare(bv_c);
        assert(cmp == 0);
    }
    
    // bit-vector insert checks
    {
        bvect bv;
        bv.resize(10);
        bv.insert(120303030, true);
        assert(bv.test(120303030));
        assert(bv.count()==1);
        assert(bv.size() == 120303030+1);
    }
    
    {
        bvect bv { 120303030u, 120303031u };
        bvect bv1(bv);
        BVectorInsert(&bv, 120303031u, true);
        bv1.insert(120303031u, true);
        int cmp = bv1.compare(bv);
        assert(cmp==0);
        bv.optimize();
        bv1.optimize();
        BVectorInsert(&bv, 120303031u, true);
        bv1.insert(120303031u, true);
        cmp = bv1.compare(bv);
        assert(cmp==0);
    }
    
    {
        bvect bv, bv1;
        bv.set(10);
        bv.set_range(1203030u, 1203030u+65535u*10u);
        bv1 = bv;
        BVectorInsert(&bv, 1203030u+10, false);
        bv1.insert(1203030u+10, false);
        int cmp = bv1.compare(bv);
        assert(cmp==0);
    }
    
    {
        std::cout << "INSERT stress (large vector insert)..\n";
        bvect bv;
        generate_bvector(bv, 40000000);
        bvect bv_control(bv);
        
        unsigned max_shifts = 10000;
        for (unsigned i = 0; i < max_shifts; ++i)
        {
            bvect bv2(bm::BM_GAP);
            bv2 = bv_control;
            bv2.optimize();
            
            unsigned i_pos = rand()%40000000;
            
            BVectorInsert(&bv_control, i_pos, i & 1u);
            bv.insert(i_pos, i & 1u);
            int cmp = bv.compare(bv_control);
            if (cmp != 0)
            {
                DetailedCompareBVectors(bv, bv_control);
            }
            assert(cmp==0);
            
            bv2.insert(i_pos, i & 1u);
            cmp = bv2.compare(bv_control);
            assert(cmp==0);
            
            if ((i % 16) == 0)
            {
                cout << "\r" << i << "/" << max_shifts << flush;
            }
        } // for i
    }
    cout << "ok.\n";


    cout << "---------------------------- Bvector INSERT test OK" << endl;
}

static
void BvectorEraseTest()
{
    cout << "---------------------------- Bvector ERASE test" << endl;
    
    {
        bvect bv;
        bv.erase(100);
        assert(!bv.any());
    }
    
    {
        bvect bv { 1, 2, 3 };
        bvect bv_c { 1, 2 };
        bv.erase(1);
        print_bv(bv);
        int cmp = bv.compare(bv_c);
        assert(cmp == 0);
    }

    {
        bvect bv {100, 65536 };
        bvect bv_c(bv);
        bv.optimize();
        bv.erase(99);
        print_bv(bv);
        BVectorErase(&bv_c, 99);
        
        assert(bv.test(99));
        assert(bv.test(65535));
        assert(bv.count()==2);
        int cmp = bv.compare(bv_c);
        assert(!cmp);
    }
    
    {
        bvect bv;
        bv.set_range(65536, 65536 + 65536);
        bvect bv_c(bv);

        unsigned cnt1 = bv.count();
        bv.optimize();
        bv.erase(0);
        BVectorErase(&bv_c, 0);
        
        unsigned cnt2 = bv.count();
        assert(cnt1 == cnt2);
        unsigned cnt3 = bv.count_range(65535, 65535 + 65536);
        assert(cnt3 == cnt1);
        
        struct bvect::statistics st;
        bv.calc_stat(&st);
        assert(st.bit_blocks == 1);
        int cmp = bv.compare(bv_c);
        assert(!cmp);
    }
    
    {
        bvect bv;
        bv.set_range(65536, 65536 + 65535);
        unsigned cnt1 = bv.count();
        bv.optimize();
        assert(cnt1 == bv.count());
        bv.erase(0);
        unsigned cnt2 = bv.count();
        assert(cnt1 == cnt2);
        unsigned cnt3 = bv.count_range(65535, 65535 + 65535);
        assert(cnt3 == cnt1);
    }
    
    {
        bvect bv;
        bv.invert();
        unsigned cnt1 = bv.count();
        bv.erase(65536);
        unsigned cnt2 = bv.count();
        cout << cnt1 << " " << cnt2 << endl;
        assert(cnt1 == (cnt2 + 1));
        assert(!bv.test(bm::id_max-1));
        
        struct bvect::statistics st;
        bv.calc_stat(&st);
        assert(st.bit_blocks == 2);

    }
    
    // test how emty blocks get deallocated on left shift
    {
        unsigned start = 1000000;
        bvect bv;
        bv.set(start);
        unsigned finish = 10;
        for(;true;)
        {
            bv.erase(finish);
            --start;
            if (start == finish)
            {
                assert(bv.test(start));
                bv.erase(finish);
                assert(!bv.test(start));

                struct bvect::statistics st;
                bv.calc_stat(&st);
                assert(st.bit_blocks == 1);

                break;
            }
            assert(bv.test(start));
            unsigned cnt = bv.count();
            assert(cnt == 1);
            
            struct bvect::statistics st;
            bv.calc_stat(&st);
            assert(st.bit_blocks == 1);
        } // for
    }
    
    cout << "bit erase stress test" << endl;
    {
        std::random_device rd;

        bvect bv;
        generate_bvector(bv, 750000000, false);
        bvect bv2(bv);
        
        bvect bv_c(bv);
        
        unsigned max_erase = 256;
        
        for(unsigned k = 0; k < max_erase; ++k)
        {
            bm::id_t pos;
            unsigned from = rd();
            bool b = bv.find(from, pos);
            if (!b)
            {
                bool any = bv.any();
                if (!any)
                    break;
                pos = 0;
            }
            bv.erase(pos);
            bv2.erase(pos);
            BVectorErase(&bv_c, pos);
            
            int cmp = bv.compare(bv_c);
            if (cmp != 0)
            {
                cerr << "Erase test failed! at pos=" << pos << endl;
                exit(1);
            }
            cmp = bv2.compare(bv_c);
            if (cmp != 0)
            {
                cerr << "2. Erase test failed! at pos=" << pos << endl;
                exit(1);
            }

            
            if ((k % 4) == 0)
            {
                cout << "\r" << k << "/" << max_erase << flush;
                bv.optimize();
            }

        } // for
    }
    cout << "\nOK" << endl;

    cout << "---------------------------- Bvector ERASE test OK" << endl;
}


// -----------------------------------------------------------------------

static
void TestRandomSubset(const bvect& bv, bm::random_subset<bvect>& rsub)
{
    bvect bv_subset;
    bvect::size_type bcnt = bv.count();

    bvect::size_type samples[] =
      { 0, 1, 2, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, bcnt / 5, bcnt / 4, bcnt / 3, bcnt / 2, (bcnt * 2)/3, bcnt };
    bvect::size_type samples_size = sizeof(samples)/sizeof(*samples);

    printf("Taking random sub-sets: ");
    
    for (unsigned i = 0; i < samples_size; ++i)
    {
        unsigned sample_count = samples[i];
        printf(" %u, ", sample_count);
        rsub.sample(bv_subset, bv, sample_count);
        if (sample_count > bcnt)
            sample_count = bcnt;

        if (sample_count != bv_subset.count())
        {
            printf("\nRandom subset failed! sample_count = %u result_count=%u\n", 
                   sample_count,
                   bv_subset.count());
            exit(1);
        }
        {
            bvect bv_subset_copy(bv_subset);
            bvect bv_set_copy(bv);
            bv_set_copy.invert();
            bv_subset_copy -= bv_set_copy;
            int res = bv_subset_copy.compare(bv_subset);
            if (res != 0)
            {
                printf("\nRandom subset failed! inverted set MINUS error! \n");
                exit(1);
            }
        }

        bv_subset -= bv;
        if (bv_subset.count() != 0)
        {
            printf("\nRandom subset failed! Extra bits set! \n");
            exit(1);
        }
        
        
    }
    printf("\n");
}

static
void SimpleRandomFillTest()
{
    assert(ITERATIONS < BITVECT_SIZE);

    bm::random_subset<bvect> rsub;

    cout << "-------------------------- SimpleRandomFillTest" << endl;

    printf("Test for Random inverted subset.");

    {
        bvect bv;
        bv.invert();
        TestRandomSubset(bv, rsub);
    }


    {
    printf("Simple random fill test 1.");
    bvect_mini   bvect_min(BITVECT_SIZE);
    bvect      bvect_full;
    bvect_full.set_new_blocks_strat(bm::BM_BIT);


    unsigned iter = ITERATIONS / 5;

    printf("\nSimple Random fill test ITERATIONS = %i\n", iter);

    bvect_min.set_bit(0);
    bvect_full.set_bit(0);

    unsigned i;
    for (i = 0; i < iter; ++i)
    {
        unsigned num = unsigned(::rand()) % iter;
        bvect_min.set_bit(num);
        bvect_full.set_bit(num);
        if ((i % 1000) == 0) cout << "." << flush;
        CheckCountRange(bvect_full, 0, num);
        CheckCountRange(bvect_full, num, num+iter);
    }

    CheckVectors(bvect_min, bvect_full, iter);
    CheckCountRange(bvect_full, 0, iter);

    TestRandomSubset(bvect_full, rsub);

    printf("Simple random fill test 2.");

    for(i = 0; i < iter; ++i)
    {
        unsigned num = unsigned(::rand()) % iter;
        bvect_min.clear_bit(num);
        bvect_full.clear_bit(num);
    }

    CheckVectors(bvect_min, bvect_full, iter);
    }


    {
    printf("\nSimple random fill test 3.\n");
    bvect_mini   bvect_min(BITVECT_SIZE);
    bvect      bvect_full(bm::BM_GAP);


    unsigned iter = ITERATIONS;

    printf("\nSimple Random fill test ITERATIONS = %i\n", iter);

    unsigned i;
    for(i = 0; i < iter; ++i)
    {
        unsigned num = unsigned(::rand()) % iter;
        bvect_min.set_bit(num);
        bvect_full.set_bit(num);
        CheckCountRange(bvect_full, 0, 65535);
        CheckCountRange(bvect_full, 0, num);
        CheckCountRange(bvect_full, num, num+iter);
    }

    CheckVectors(bvect_min, bvect_full, iter);

    TestRandomSubset(bvect_full, rsub);

    printf("Simple random fill test 4.");

    for(i = 0; i < iter; ++i)
    {
        unsigned num = unsigned(rand()) % iter;
        bvect_min.clear_bit(num);
        bvect_full.clear_bit(num);
        CheckCountRange(bvect_full, 0, num);
        CheckCountRange(bvect_full, num, num+iter);
    }

    CheckVectors(bvect_min, bvect_full, iter);
    CheckCountRange(bvect_full, 0, iter);

    TestRandomSubset(bvect_full, rsub);
    }

}



static
void RangeRandomFillTest()
{
    assert(ITERATIONS < BITVECT_SIZE);

    cout << "----------------------------------- RangeRandomFillTest" << endl;

    {
    bvect_mini   bvect_min(BITVECT_SIZE);
    bvect     bvect_full;

    printf("Range Random fill test\n");

    bvect::size_type min = BITVECT_SIZE / 2;
    bvect::size_type max = BITVECT_SIZE / 2 + ITERATIONS;
    if (max > BITVECT_SIZE) 
        max = BITVECT_SIZE - 1;

    FillSets(&bvect_min, &bvect_full, min, max, 0);

    CheckVectors(bvect_min, bvect_full, BITVECT_SIZE);
    CheckCountRange(bvect_full, min, max);

    }

    
    {
    bvect_mini   bvect_min(BITVECT_SIZE);
    bvect     bvect_full;

    printf("Range Random fill test\n");

    bvect::size_type min = BITVECT_SIZE / 2;
    bvect::size_type max = BITVECT_SIZE / 2 + ITERATIONS;
    if (max > BITVECT_SIZE) 
        max = BITVECT_SIZE - 1;

    FillSetsIntervals(&bvect_min, bvect_full, min, max, 4);

    CheckVectors(bvect_min, bvect_full, BITVECT_SIZE);
    CheckCountRange(bvect_full, min, max);
    }
}

static
void RangeCopyTest()
{
    cout << "----------------------------------- RangeCopyTest" << endl;
    {
        const unsigned to_max = 65536 * bm::set_sub_array_size + 10;
        cout << "Basic range-copy test" << endl;
        bvect     bvect1
        { 10, 20, 21, 100, 65535, 65536, 100000, to_max/2, to_max-1, to_max };
        

        CheckRangeCopy(bvect1, 0, 0);
        CheckRangeCopy(bvect1, 10, 10);
        CheckRangeCopy(bvect1, 15, 15);
        CheckRangeCopy(bvect1, 65535, 65535);
        CheckRangeCopy(bvect1, 65536, 65536);
        CheckRangeCopy(bvect1, 65535, 65536);

        for (unsigned k = 0; k < 2; ++k)
        {
            cout << "Pass " << k << "-0" << endl;
            for (unsigned i = 0; i < to_max; ++i)
            {
                CheckRangeCopy(bvect1, i, to_max);
            }
            cout << "Pass " << k << "-1" << endl;
            for (unsigned i = to_max-1; i > 0; --i)
            {
                CheckRangeCopy(bvect1, 0, i);
            }
            cout << "Pass " << k << "-2" << endl;
            auto to = to_max;
            for (unsigned i = 0; i != to_max; ++i, --to)
            {
                CheckRangeCopy(bvect1, i, to_max);
            }
            bvect1.optimize();
        } // for k
        cout << "OK" << endl;
    }

    {
        cout << "Inverted vector stress test" << endl;
        bvect     bvect1;
        bvect1.invert();

        {
            const bvect::size_type to_max = bm::gap_max_bits * bm::set_sub_array_size;

            cout << "T1" << endl;
            auto to = to_max;

            for (unsigned i = 0; i < to; ++i)
            {
                CheckRangeCopy(bvect1, i, to);
            }
            
            cout << "T2" << endl;
            to = to_max;
            for (unsigned i = to; i > 0; --i)
            {
                CheckRangeCopy(bvect1, 0, i);
            }
            
            cout << "T3" << endl;
            to = to_max;
            for (unsigned i = 0; i != to; ++i, --to)
            {
                CheckRangeCopy(bvect1, i, to);
            }
            cout << "T4" << endl;
            to = bm::id_max-1 - to_max - 100;
            for (unsigned i = to; i < bm::id_max; ++i)
            {
                CheckRangeCopy(bvect1, i, bm::id_max);
                if ((i & 0xFFFF) == 0)
                    cout << "\r" << i << flush;
            }
            cout << endl;
            to = bm::id_max-1 - to_max - 100;
            for (unsigned i = to; i < bm::id_max-(65536 * 3); i+=65536 * 2)
            {
                bvect1.set(i, false);
            }
            for (unsigned k = 0; k < 2; ++k)
            {
                cout << "T5 pass=" << k << endl;
                to = bm::id_max-1 - to_max - (bm::gap_max_bits/2);
                for (unsigned i = to; i < bm::id_max; ++i)
                {
                    CheckRangeCopy(bvect1, i, bm::id_max);
                    if ((i & 0xFFFF) == 0)
                        cout << "\r" << i << flush;
                }
                bvect1.optimize();
            }
        }
    }
    
    

    cout << "----------------------------------- RangeCopyTest OK" << endl;
}



static
void AndOperationsTest(bool detailed)
{
    assert(ITERATIONS < BITVECT_SIZE);

    cout << "----------------------------------- AndOperationTest" << endl;

    {

    bvect_mini   bvect_min1(256);
    bvect_mini   bvect_min2(256);
    bvect        bvect_full1;
    bvect        bvect_full2;

    bvect_full1.set_new_blocks_strat(bm::BM_GAP);
    bvect_full2.set_new_blocks_strat(bm::BM_GAP);



    printf("AND test\n");

    bvect_min1.set_bit(1);
    bvect_min1.set_bit(12);
    bvect_min1.set_bit(13);

    bvect_min2.set_bit(12);
    bvect_min2.set_bit(13);

    bvect_min1.combine_and(bvect_min2);

    bvect_full1.set_bit(1);
    bvect_full1.set_bit(12);
    bvect_full1.set_bit(13);

    bvect_full2.set_bit(12);
    bvect_full2.set_bit(13);

    bm::id_t predicted_count = bm::count_and(bvect_full1, bvect_full2);

    bm::id_t predicted_any = bm::any_and(bvect_full1, bvect_full2);
    if (predicted_any == 0 && predicted_count != 0)
    {
        cout << "Predicted any error!" << endl;
        exit(1);
    }

    bvect    bv_target_s;
    SerializationOperation2Test(&bv_target_s,
                                bvect_full1,
                                bvect_full2,
                                predicted_count,
                                set_COUNT_AND,
                                set_AND);


    bvect_full1.bit_and(bvect_full2);

    bm::id_t count = bvect_full1.count();
    if (count != predicted_count)
    {
        cout << "Predicted count error!" << endl;
        exit(1);
    }

    CheckVectors(bvect_min1, bvect_full1, 256, detailed);
    CheckVectors(bvect_min1, bv_target_s, 256, detailed);
    CheckCountRange(bvect_full1, 0, 256);

    }
    
    {
        bvect        bvect1;
        bvect        bvect2 { 256, 165535 };
        bvect        bvect_control { 256  };
        bvect1.set_range(0, 100000);

        bvect2.optimize();

        bvect1 &= bvect2;
        int res = bvect1.compare(bvect_control);
        assert(res==0);
    }
    
    {
        bvect        bvect1 { 1, 2, 3};
        bvect        bvect2 { 256, 165535 };
        bvect        bvect_control;
        bvect1.optimize();

        bvect1 &= bvect2;
        int res = bvect1.compare(bvect_control);
        assert(res==0);
    }


    {

    bvect_mini   bvect_min1(BITVECT_SIZE);
    bvect_mini   bvect_min2(BITVECT_SIZE);
    bvect        bvect_full1;
    bvect        bvect_full2;


    printf("AND test stage 1.\n");

    for (unsigned i = 0; i < 112; ++i)
    {
        bvect_min1.set_bit(i);
        bvect_full1.set_bit(i);

        bvect_min2.set_bit(i);
        bvect_full2.set_bit(i);
    }

    CheckVectors(bvect_min1, bvect_full1, BITVECT_SIZE/10+10, detailed);
    CheckCountRange(bvect_full1, 0, BITVECT_SIZE/10+10);

//    FillSets(&bvect_min1, &bvect_full1, 1, BITVECT_SIZE/7, 0);
//    FillSets(&bvect_min2, &bvect_full2, 1, BITVECT_SIZE/7, 0);

    bvect_min1.combine_and(bvect_min2);

    bm::id_t predicted_count = bm::count_and(bvect_full1,bvect_full2);
    bm::id_t predicted_any = bm::any_and(bvect_full1, bvect_full2);
    if (predicted_any == 0 && predicted_count != 0)
    {
        cout << "Predicted any error!" << endl;
        exit(1);
    }

    bvect    bv_target_s;
    SerializationOperation2Test(&bv_target_s,
                                bvect_full1,
                                bvect_full2,
                                predicted_count,
                                set_COUNT_AND,
                                set_AND);

    bvect_full1.bit_and(bvect_full2);

    bm::id_t count = bvect_full1.count();
    if (count != predicted_count)
    {
        cout << "Predicted count error!" << endl;
        exit(1);
    }

    CheckVectors(bvect_min1, bvect_full1, BITVECT_SIZE/10+10, detailed);
    CheckVectors(bvect_min1, bv_target_s, BITVECT_SIZE/10+10, detailed);
    CheckCountRange(bvect_full1, 0, BITVECT_SIZE/10+10);

    }


    {

    bvect_mini   bvect_min1(BITVECT_SIZE);
    bvect_mini   bvect_min2(BITVECT_SIZE);
    bvect        bvect_full1;
    bvect        bvect_full2;

    bvect_full1.set_new_blocks_strat(bm::BM_GAP);
    bvect_full2.set_new_blocks_strat(bm::BM_GAP);

    printf("AND test stage 2.\n");


    FillSets(&bvect_min1, &bvect_full1, 1, BITVECT_SIZE/7, 0);
    FillSets(&bvect_min2, &bvect_full2, 1, BITVECT_SIZE/7, 0);

    bm::id_t predicted_count = bm::count_and(bvect_full1,bvect_full2);
    bm::id_t predicted_any = bm::any_and(bvect_full1, bvect_full2);
    if (predicted_any == 0 && predicted_count != 0)
    {
        cout << "Predicted any error!" << endl;
        exit(1);
    }

    bvect    bv_target_s;
    SerializationOperation2Test(&bv_target_s,
                                bvect_full1,
                                bvect_full2,
                                predicted_count,
                                set_COUNT_AND,
                                set_AND);

    bvect_min1.combine_and(bvect_min2);

    bvect_full1.bit_and(bvect_full2);

    bm::id_t count = bvect_full1.count();
    if (count != predicted_count)
    {
        cout << "Predicted count error!" << endl;
        print_stat(bvect_full1);
        exit(1);
    }

    CheckVectors(bvect_min1, bvect_full1, BITVECT_SIZE/10+10, detailed);
    CheckVectors(bvect_min1, bv_target_s, BITVECT_SIZE/10+10, detailed);
    CheckCountRange(bvect_full1, 0, BITVECT_SIZE/10+10);

    }

    {

    bvect_mini   bvect_min1(BITVECT_SIZE);
    bvect_mini   bvect_min2(BITVECT_SIZE);
    bvect        bvect_full1;
    bvect        bvect_full2;

    bvect_full1.set_new_blocks_strat(bm::BM_BIT);
    bvect_full2.set_new_blocks_strat(bm::BM_BIT);

    cout << "------------------------------" << endl;
    printf("AND test stage 3.\n");


    FillSets(&bvect_min1, &bvect_full1, 1, BITVECT_SIZE/5, 2);
    FillSets(&bvect_min2, &bvect_full2, 1, BITVECT_SIZE/5, 2);

    bvect_min1.combine_and(bvect_min2);

    bm::id_t predicted_count = bm::count_and(bvect_full1, bvect_full2);
    bm::id_t predicted_any = bm::any_and(bvect_full1, bvect_full2);
    if (predicted_any == 0 && predicted_count != 0)
    {
        cout << "Predicted any error!" << endl;
        exit(1);
    }
    
    bvect    bv_target_s;
    SerializationOperation2Test(&bv_target_s,
                                bvect_full1,
                                bvect_full2,
                                predicted_count,
                                set_COUNT_AND,
                                set_AND);

    bvect_full1.bit_and(bvect_full2);

    bm::id_t count = bvect_full1.count();
    if (count != predicted_count)
    {
        cout << "Predicted count error!" << endl;
        exit(1);
    }

    {
        bvect::size_type pos;
        bool f = bvect_full1.find_first_mismatch(bv_target_s, pos);
        if (f)
        {
            cerr << "Mismatch at position = " << pos << endl;
            assert(0); exit(1);
        }
    }

    CheckVectors(bvect_min1, bvect_full1, BITVECT_SIZE, detailed);
    CheckVectors(bvect_min1, bv_target_s, BITVECT_SIZE, detailed);
    CheckCountRange(bvect_full1, 0, BITVECT_SIZE);

    BM_DECLARE_TEMP_BLOCK(tb)
    bvect_full1.optimize(tb);
    CheckVectors(bvect_min1, bvect_full1, BITVECT_SIZE, detailed);
    CheckCountRange(bvect_full1, 0, BITVECT_SIZE);
    CheckCountRange(bvect_full1, BITVECT_SIZE/2, BITVECT_SIZE);

    }

    printf("AND test stage 4. combine_and_sorted\n");
    {
    unsigned ids[] = {0, 1, 2, 3, 10, 65535, 65536, 65535*2, 65535*3};
    unsigned to_add = sizeof(ids)/sizeof(unsigned);
    bvect        bvect_full1;
    bvect        bvect_full2;    
    bvect_mini   bvect_min1(BITVECT_SIZE);
    bvect_mini   bvect_min2(BITVECT_SIZE);

    bvect_full1.set_new_blocks_strat(bm::BM_GAP);
    bvect_full2.set_new_blocks_strat(bm::BM_GAP);
    
    for (unsigned i = 2; i < to_add; ++i)
    {
        bvect_full1.set(ids[i]);
        bvect_min1.set_bit(ids[i]);
        bvect_full2.set(ids[i]);
        bvect_min2.set_bit(ids[i]);
    }
    
    unsigned* first = ids;
    unsigned* last = ids + to_add;
    
    bvect_min1.combine_and(bvect_min2);

    bm::combine_and_sorted(bvect_full1, first, last);
    CheckVectors(bvect_min1, bvect_full1, BITVECT_SIZE, detailed);
    }
    
    {
    bvect        bvect1 { 1, 10, 12 };
    bvect        bvect2 { 2, 15, 165535 };
    
    bvect1 &= bvect2;
    
    bvect::statistics st;
    bvect1.calc_stat(&st);
    if (st.bit_blocks != 0 || st.gap_blocks != 0)
    {
        cerr << "Error: AND-optimization reduction failed!" << endl;
        exit(1);
    }
    }
    
    
    // ------------------------------------------
    // 2-way AND
    //
    {
        bvect        bv1 { 0, 1 };
        bvect        bv2 ;
        bv2.bit_and(bv1, bv2, bvect::opt_compress);
        int cmp = bv2.any();
        assert(cmp == 0);
    }

    {
        bvect        bv1 { 0, 1 };
        bvect        bv2 { 1, 3 };
        bvect bv1c(bv1);
        bv1c.bit_and(bv2);

        bvect bv;
        bv.bit_and(bv1, bv2, bvect::opt_compress);
        int cmp = bv.compare(bv1c);
        assert(cmp == 0);
        struct bvect::statistics st1;
        bv.calc_stat(&st1);
        assert(!st1.bit_blocks);
        assert(st1.gap_blocks == 1);
    }
    
    {
        bvect        bv1 { 0, 1 };
        bvect        bv2;
        for (unsigned i = 2; i < 65536; ++i)
            bv2.set(i);
        
        bvect bv1c(bv1);
        bv1c.bit_and(bv2);

        bvect bv;
        bv.bit_and(bv1, bv2, bvect::opt_none); // should detect 0 automatically
        int cmp = bv.compare(bv1c);
        assert(cmp == 0);
        struct bvect::statistics st1;
        bv.calc_stat(&st1);
        assert(!st1.bit_blocks);
        assert(!st1.gap_blocks);
    }
    
    {
        bvect        bv1 { 0, 1 };
        bvect        bv2 { 1 };
        bv1.clear(0); bv1.clear(1);
        bv2.clear(1);
        
        bvect bv;
        bv.bit_and(bv1, bv2, bvect::opt_none); // should detect empty automatically

        struct bvect::statistics st1;
        bv.calc_stat(&st1);
        assert(!st1.bit_blocks);
        assert(!st1.gap_blocks);
        assert(!st1.ptr_sub_blocks);
    }


    {
        bvect        bvect_full1;
        bvect        bvect_full2;
        bvect_full1.invert();
        bvect_full2.set();
        
        {
        bvect    bv_target_s;
        bv_target_s.bit_and(bvect_full1, bvect_full2, bvect::opt_none);
        int cmp = bv_target_s.compare(bvect_full1);
        assert(cmp == 0);
        }
    }
    {
        bvect        bv1;
        bvect        bv2;
        bv1.invert();
        bv2.set();
        bv2.set(bm::id_max/2, false);
        
        {
            bvect    bv_s;
            bv_s.bit_and(bv1, bv2, bvect::opt_none);
            int cmp = bv_s.compare(bv2);
            assert(cmp == 0);
        }
        {
            bvect    bv_s;
            bv_s.bit_and(bv2, bv1, bvect::opt_none);
            int cmp = bv_s.compare(bv2);
            assert(cmp == 0);
        }
        bv1 &= bv2;
        int cmp = bv1.compare(bv2);
        assert(cmp == 0);
    }

    {
        bvect        bv1 { 0, 1 };
        bvect        bv2 { 1, 3 };
        bv2.optimize();
        bvect bv1c(bv1);
        bv1c.bit_and(bv2);

        {
            bvect bv;
            bv.bit_and(bv1, bv2, bvect::opt_compress);
            int cmp = bv.compare(bv1c);
            assert(cmp == 0);
        }
        bv1.optimize();
        {
            bvect bv;
            bv.bit_and(bv1, bv2, bvect::opt_compress);
            int cmp = bv.compare(bv1c);
            assert(cmp == 0);
        }
        bv2.clear();
        bv2.invert();
        {
            bvect bv;
            bv.bit_and(bv1, bv2, bvect::opt_compress);
            int cmp = bv.compare(bv1);
            assert(cmp == 0);
        }
    }
    
    cout << "----------------------------------- AndOperationTest OK" << endl;

}

static
void OrOperationsTest(bool detailed)
{
    assert(ITERATIONS < BITVECT_SIZE);

    cout << "----------------------------------- OrOperationTest" << endl;
    
    {

    bvect_mini   bvect_min1(256);
    bvect_mini   bvect_min2(256);
    bvect        bvect_full1;
    bvect        bvect_full2;

    bvect_full1.set_new_blocks_strat(bm::BM_GAP);
    bvect_full2.set_new_blocks_strat(bm::BM_GAP);



    printf("OR test\n");

    bvect_min1.set_bit(1);
    bvect_min1.set_bit(12);
    bvect_min1.set_bit(13);

    bvect_min2.set_bit(12);
    bvect_min2.set_bit(13);

    bvect_min1.combine_or(bvect_min2);

    bvect_full1.set_bit(1);
    bvect_full1.set_bit(12);
    bvect_full1.set_bit(13);

    bvect_full2.set_bit(12);
    bvect_full2.set_bit(13);
    
    bm::id_t predicted_count = bm::count_or(bvect_full1, bvect_full2);    
    bm::id_t predicted_any = bm::any_or(bvect_full1, bvect_full2);
    if (predicted_any == 0 && predicted_count != 0)
    {
        cout << "Predicted any error!" << endl;
        exit(1);
    }

    bvect    bv_target_s;
    SerializationOperation2Test(&bv_target_s,
                                bvect_full1,
                                bvect_full2,
                                predicted_count,
                                set_COUNT_OR,
                                set_OR);


    bvect_full1.bit_or(bvect_full2);

    bm::id_t count = bvect_full1.count();
    if (count != predicted_count)
    {
        cout << "Predicted count error!" << endl;
        cout << predicted_count << " " << count << endl;
        print_stat(bvect_full1);
        exit(1);
    }

    {
        bvect::size_type pos;
        bool f = bvect_full1.find_first_mismatch(bv_target_s, pos);
        if (f)
        {
            cerr << "Mismatch found pos=" << pos << endl;
            assert(0); exit(1);
        }
    }

    CheckVectors(bvect_min1, bvect_full1, 256, detailed);
    CheckVectors(bvect_min1, bv_target_s, 256, detailed);
    CheckCountRange(bvect_full1, 0, 256);
    CheckCountRange(bvect_full1, 128, 256);
    }

    {

    bvect_mini   bvect_min1(BITVECT_SIZE);
    bvect_mini   bvect_min2(BITVECT_SIZE);
    bvect        bvect_full1;
    bvect        bvect_full2;

    bvect_full1.set_new_blocks_strat(bm::BM_GAP);
    bvect_full2.set_new_blocks_strat(bm::BM_GAP);

    printf("OR test stage 2.\n");


    FillSets(&bvect_min1, &bvect_full1, 1, BITVECT_SIZE/7, 0);
    FillSets(&bvect_min2, &bvect_full2, 1, BITVECT_SIZE/7, 0);

    bvect_min1.combine_or(bvect_min2);

    bm::id_t predicted_count = bm::count_or(bvect_full1, bvect_full2);    
    bm::id_t predicted_any = bm::any_or(bvect_full1, bvect_full2);
    if (predicted_any == 0 && predicted_count != 0)
    {
        cout << "Predicted any error!" << endl;
        exit(1);
    }

    bvect    bv_target_s;
    SerializationOperation2Test(&bv_target_s,
                                bvect_full1,
                                bvect_full2,
                                predicted_count,
                                set_COUNT_OR,
                                set_OR);


    bvect_full1.bit_or(bvect_full2);

    bm::id_t count = bvect_full1.count();
    if (count != predicted_count)
    {
        cout << "Predicted count error!" << endl;
        exit(1);
    }

    {
        bvect::size_type pos;
        bool f = bvect_full1.find_first_mismatch(bv_target_s, pos);
        if (f)
        {
            cerr << "Mismatch found pos=" << pos << endl;
            assert(0); exit(1);
        }
    }

    CheckVectors(bvect_min1, bvect_full1, BITVECT_SIZE/10+10, detailed);
    CheckVectors(bvect_min1, bvect_full1, BITVECT_SIZE/10+10, detailed);
    CheckCountRange(bvect_full1, 0, BITVECT_SIZE/10+10);

    }
    
    {
        bvect bv1;
        bvect bv2;
        bv1.flip(); bv2.flip();
        unsigned cnt1 = bv1.count();
        bv1.bit_or(bv2);
        unsigned cnt2 = bv1.count();
        assert(cnt1 == cnt2);
        struct bvect::statistics st;
        bv1.calc_stat(&st);
        auto bcnt = st.bit_blocks + st.gap_blocks;
        assert(bcnt == 1);
    }

    {

    bvect_mini   bvect_min1(BITVECT_SIZE);
    bvect_mini   bvect_min2(BITVECT_SIZE);
    bvect        bvect_full1;
    bvect        bvect_full2;

    bvect_full1.set_new_blocks_strat(bm::BM_BIT);
    bvect_full2.set_new_blocks_strat(bm::BM_BIT);

    cout << "------------------------------" << endl;
    printf("OR test stage 3.\n");


    FillSets(&bvect_min1, &bvect_full1, 1, BITVECT_SIZE/5, 2);
    FillSets(&bvect_min2, &bvect_full2, 1, BITVECT_SIZE/5, 2);

    bvect_min1.combine_or(bvect_min2);
    unsigned mcnt = bvect_min1.bit_count();

    cout << mcnt << endl;
    
    bm::id_t predicted_count = bm::count_or(bvect_full1, bvect_full2);    
    cout << predicted_count << endl;
    bm::id_t predicted_any = bm::any_or(bvect_full1, bvect_full2);
    if (predicted_any == 0 && predicted_count != 0)
    {
        cout << "Predicted any error!" << endl;
        exit(1);
    }

    bvect    bv_target_s;
    SerializationOperation2Test(&bv_target_s,
                                bvect_full1,
                                bvect_full2,
                                predicted_count,
                                set_COUNT_OR,
                                set_OR);

    bvect_full1.bit_or(bvect_full2);

    bm::id_t count = bvect_full1.count();
    if (count != predicted_count)
    {
        cout << "Predicted count error!" << endl;
        exit(1);
    }

    CheckVectors(bvect_min1, bvect_full1, BITVECT_SIZE);

    BM_DECLARE_TEMP_BLOCK(tb)
    bvect_full1.optimize(tb);

    {
        bvect::size_type pos;
        bool f = bvect_full1.find_first_mismatch(bv_target_s, pos);
        if (f)
        {
            cerr << "Mismatch found pos=" << pos << endl;
            assert(0); exit(1);
        }
    }

    CheckVectors(bvect_min1, bvect_full1, BITVECT_SIZE, detailed);
    CheckVectors(bvect_min1, bv_target_s, BITVECT_SIZE, detailed);
    CheckCountRange(bvect_full1, 0, BITVECT_SIZE);


    }
    
    cout << "Testing combine_or" << endl;
    
    {
    
    bvect        bvect_full1;
    bvect        bvect_full2;
    bvect_mini   bvect_min1(BITVECT_SIZE);
    
    bvect_full1.set_new_blocks_strat(bm::BM_GAP);
    bvect_full2.set_new_blocks_strat(bm::BM_GAP);

    unsigned ids[10000];
    unsigned to_add = 10000;
    
    unsigned bn = 0;
    for (unsigned i = 0; i < to_add; ++i)
    {
        ids[i] = bn;
        bvect_full2.set(bn);
        bvect_min1.set_bit(bn);
        bn += 15;
    }
    
    unsigned* first = ids;
    unsigned* last = ids + to_add;
    
    bm::combine_or(bvect_full1, first, last);

    CheckVectors(bvect_min1, bvect_full1, BITVECT_SIZE);
    
    bm::combine_or(bvect_full1, first, last);
    CheckVectors(bvect_min1, bvect_full1, BITVECT_SIZE);
    
    }
    
    
    {
    unsigned ids[] = {0, 65536, 65535, 65535*3, 65535*2, 10};
    unsigned to_add = sizeof(ids)/sizeof(unsigned);
    bvect        bvect_full1;
    bvect        bvect_full2;    
    bvect_mini   bvect_min1(BITVECT_SIZE);

    bvect_full1.set_new_blocks_strat(bm::BM_GAP);
    bvect_full2.set_new_blocks_strat(bm::BM_GAP);
    
    unsigned bn = 0;
    for (unsigned i = 0; i < to_add; ++i)
    {
        ids[i] = bn;
        bvect_full2.set(bn);
        bvect_min1.set_bit(bn);
        bn += 15;
    }
    
    unsigned* first = ids;
    unsigned* last = ids + to_add;
    
    bm::combine_or(bvect_full1, first, last);
    CheckVectors(bvect_min1, bvect_full1, BITVECT_SIZE);

    bm::combine_or(bvect_full1, first, last);
    CheckVectors(bvect_min1, bvect_full1, BITVECT_SIZE);    
    }
    
    {
    bvect        bv0;
    bvect        bv1 { 0, 36500 };
    bvect        bv2 { 128000, bm::id_max-1 };
    
    bvect        bvc(bv1);
    bvc |= bv2;

    bv1.merge(bv2);
    int cmp = bv1.compare(bvc);
    assert(cmp==0);

    struct bvect::statistics st2;
    bv2.calc_stat(&st2);
    auto bcnt = st2.bit_blocks + st2.gap_blocks;
    assert(bcnt == 0);
    
    
    bv0.merge(bv1);
    struct bvect::statistics st1;
    bv1.calc_stat(&st1);
    bcnt = st1.bit_blocks + st1.gap_blocks;
    assert(bcnt == 0);

    }

    // ------------------------------------------
    // 2-way OR
    //
    {
        bvect        bv1 { 0, 1 };
        bvect        bv2;
        bv2.bit_or(bv1, bv2, bvect::opt_compress);
        int cmp = bv1.compare(bv2);
        assert(cmp == 0);
    }

    {
        bvect        bv1 { 0, 1 };
        bvect        bv2 { 2, 3 };
        bvect bv1c(bv1);
        bv1c.bit_or(bv2);

        bvect bv;
        bv.bit_or(bv1, bv2, bvect::opt_compress);
        int cmp = bv.compare(bv1c);
        assert(cmp == 0);
        struct bvect::statistics st1;
        bv.calc_stat(&st1);
        assert(!st1.bit_blocks);
        assert(st1.gap_blocks == 1);
    }
    
    {
        bvect        bv1 { 0, 1 };
        bvect        bv2;
        for (unsigned i = 2; i < 65536; ++i)
            bv2.set(i);
        
        bvect bv1c(bv1);
        bv1c.bit_or(bv2);

        bvect bv;
        bv.bit_or(bv1, bv2, bvect::opt_none); // should detect FULL automatically
        int cmp = bv.compare(bv1c);
        assert(cmp == 0);
        struct bvect::statistics st1;
        bv.calc_stat(&st1);
        assert(!st1.bit_blocks);
        assert(!st1.gap_blocks);
    }
    
    {
        bvect        bv1 { 0, 1 };
        bvect        bv2 { 1 };
        bv1.clear(0); bv1.clear(1);
        bv2.clear(1);
        
        bvect bv;
        bv.bit_or(bv1, bv2, bvect::opt_none); // should detect FULL automatically

        struct bvect::statistics st1;
        bv.calc_stat(&st1);
        assert(!st1.bit_blocks);
        assert(!st1.gap_blocks);
        assert(!st1.ptr_sub_blocks);
    }

    {
        bvect        bv1;
        bvect        bv2;
        bv1.invert(); bv2.set();
        bv1 |= bv2;
        bvect bv;
        bv.bit_or(bv1, bv2, bvect::opt_compress);
        int cmp = bv.compare(bv2);
        assert(cmp == 0);
        cmp = bv.compare(bv1);
        assert(cmp == 0);
        bvect        bv3 { 10, bm::id_max - 1 };
        bv1 |= bv3;
        cmp = bv.compare(bv1);
        assert(cmp == 0);
    }
    {
        bvect        bv1;
        bvect        bv2;
        bv1.invert();
        bv2.set();
        bv2.set(bm::id_max/2, false);
        
        {
            bvect    bv_s;
            bv_s.bit_or(bv1, bv2, bvect::opt_none);
            int cmp = bv_s.compare(bv1);
            assert(cmp == 0);
        }
        {
            bvect    bv_s;
            bv_s.bit_or(bv2, bv1, bvect::opt_none);
            auto cnt = bv_s.count();
            assert(cnt == bm::id_max);
            bool b = bv_s.test(bm::id_max/2);
            assert(b);
            int cmp = bv_s.compare(bv1);
            assert(cmp == 0);
        }
        bv2 |= bv1;
        int cmp = bv1.compare(bv2);
        assert(cmp == 0);
    }
    
    {
        bvect        bv1 { 0, 1 };
        bvect        bv2 { 2, 3 };
        bv2.optimize();
        bvect bv1c(bv1);
        bv1c.bit_or(bv2);

        {
            bvect bv;
            bv.bit_or(bv1, bv2, bvect::opt_compress);
            int cmp = bv.compare(bv1c);
            assert(cmp == 0);
        }
        bv1.optimize();
        {
            bvect bv;
            bv.bit_or(bv1, bv2, bvect::opt_compress);
            int cmp = bv.compare(bv1c);
            assert(cmp == 0);
        }
        bv2.clear();
        bv2.invert();
        {
            bvect bv;
            bv |= bv2;
            int cmp = bv.compare(bv2);
            assert(cmp == 0);
            bv.bit_or(bv1, bv2, bvect::opt_compress);
            cmp = bv.compare(bv2);
            assert(cmp == 0);
        }
    }
    
    
    
    cout << "----------------------------------- OrOperationTest OK" << endl;

}


static
void SubOperationsTest(bool detailed)
{
    assert(ITERATIONS < BITVECT_SIZE);

    cout << "----------------------------------- SubOperationTest" << endl;

    {

    bvect_mini   bvect_min1(256);
    bvect_mini   bvect_min2(256);
    bvect        bvect_full1;
    bvect        bvect_full2;

    bvect_full1.set_new_blocks_strat(bm::BM_GAP);
    bvect_full2.set_new_blocks_strat(bm::BM_GAP);



    printf("SUB test\n");

    bvect_min1.set_bit(1);
    bvect_min1.set_bit(12);
    bvect_min1.set_bit(13);

    bvect_min2.set_bit(12);
    bvect_min2.set_bit(13);

    bvect_min1.combine_sub(bvect_min2);

    bvect_full1.set_bit(1);
    bvect_full1.set_bit(12);
    bvect_full1.set_bit(13);

    bvect_full2.set_bit(12);
    bvect_full2.set_bit(13);

    bm::id_t predicted_count = bm::count_sub(bvect_full1, bvect_full2);
    bm::id_t predicted_any = bm::any_sub(bvect_full1, bvect_full2);
    if (predicted_any == 0 && predicted_count != 0)
    {
        cout << "Predicted any error!" << endl;
        exit(1);
    }

    bvect    bv_target_s;
    SerializationOperation2Test(&bv_target_s,
                                bvect_full1,
                                bvect_full2,
                                predicted_count,
                                set_COUNT_SUB_AB,
                                set_SUB);


    bvect_full1.bit_sub(bvect_full2);
    
    bm::id_t count = bvect_full1.count();
    if (count != predicted_count)
    {
        cout << "Predicted count error!" << endl;
        exit(1);
    }

    {
        bvect::size_type pos;
        bool f = bvect_full1.find_first_mismatch(bv_target_s, pos);
        if (f)
        {
            cerr << "Mismatch found pos=" << pos << endl;
            assert(0); exit(1);
        }
    }

    CheckVectors(bvect_min1, bvect_full1, 256, detailed);
    CheckVectors(bvect_min1, bv_target_s, 256, detailed);
    CheckCountRange(bvect_full1, 0, 256);

    }

    {

    bvect_mini   bvect_min1(BITVECT_SIZE);
    bvect_mini   bvect_min2(BITVECT_SIZE);
    bvect        bvect_full1;
    bvect        bvect_full2;

    bvect_full1.set_new_blocks_strat(bm::BM_GAP);
    bvect_full2.set_new_blocks_strat(bm::BM_GAP);

    printf("SUB test stage 2.\n");


    FillSets(&bvect_min1, &bvect_full1, 1, BITVECT_SIZE/7, 0);
    FillSets(&bvect_min2, &bvect_full2, 1, BITVECT_SIZE/7, 0);

    bvect_min1.combine_sub(bvect_min2);

    bm::id_t predicted_count = bm::count_sub(bvect_full1, bvect_full2);
    bm::id_t predicted_any = bm::any_sub(bvect_full1, bvect_full2);
    if (predicted_any == 0 && predicted_count != 0)
    {
        cout << "Predicted any error!" << endl;
        exit(1);
    }

    bvect    bv_target_s;
    SerializationOperation2Test(&bv_target_s,
                                bvect_full1,
                                bvect_full2,
                                predicted_count,
                                set_COUNT_SUB_AB,
                                set_SUB);

    bvect_full1.bit_sub(bvect_full2);
    
    bm::id_t count = bvect_full1.count();
    if (count != predicted_count)
    {
        cout << "Predicted count error!" << endl;
        cout << predicted_count << " " << count << endl;
        print_stat(bvect_full1);    
        
        exit(1);
    }
    
    {
        bvect::size_type pos;
        bool f = bvect_full1.find_first_mismatch(bv_target_s, pos);
        if (f)
        {
            cerr << "Mismatch found pos=" << pos << endl;
            assert(0); exit(1);
        }
    }

    CheckVectors(bvect_min1, bvect_full1, BITVECT_SIZE/10+10, detailed);
    CheckVectors(bvect_min1, bv_target_s, BITVECT_SIZE/10+10, detailed);
    CheckCountRange(bvect_full1, 0, BITVECT_SIZE/10+10);

    }

    {

    bvect_mini   bvect_min1(BITVECT_SIZE);
    bvect_mini   bvect_min2(BITVECT_SIZE);
    bvect        bvect_full1;
    bvect        bvect_full2;

    bvect_full1.set_new_blocks_strat(bm::BM_BIT);
    bvect_full2.set_new_blocks_strat(bm::BM_BIT);

    cout << "------------------------------" << endl;
    printf("SUB test stage 3.\n");


    FillSets(&bvect_min1, &bvect_full1, 1, BITVECT_SIZE/5, 2);
    FillSets(&bvect_min2, &bvect_full2, 1, BITVECT_SIZE/5, 2);

    bvect_min1.combine_sub(bvect_min2);
    
    bm::id_t predicted_count = bm::count_sub(bvect_full1, bvect_full2);
    bm::id_t predicted_any = bm::any_sub(bvect_full1, bvect_full2);
    if (predicted_any == 0 && predicted_count != 0)
    {
        cout << "Predicted any error!" << endl;
        exit(1);
    }

    bvect    bv_target_s;
    SerializationOperation2Test(&bv_target_s,
                                bvect_full1,
                                bvect_full2,
                                predicted_count,
                                set_COUNT_SUB_AB,
                                set_SUB);

    bvect_full1.bit_sub(bvect_full2);

    bm::id_t count = bvect_full1.count();
    if (count != predicted_count)
    {
        cout << "Predicted count error!" << endl;
        exit(1);
    }


    CheckVectors(bvect_min1, bvect_full1, BITVECT_SIZE);

    BM_DECLARE_TEMP_BLOCK(tb)
    bvect_full1.optimize(tb);
    {
        bvect::size_type pos;
        bool f = bvect_full1.find_first_mismatch(bv_target_s, pos);
        if (f)
        {
            cerr << "Mismatch found pos=" << pos << endl;
            assert(0); exit(1);
        }
    }

    CheckVectors(bvect_min1, bvect_full1, BITVECT_SIZE, detailed);
    CheckVectors(bvect_min1, bv_target_s, BITVECT_SIZE, detailed);
    CheckCountRange(bvect_full1, 0, BITVECT_SIZE);

    }
    
    
    // ------------------------------------------
    // 2-way SUB
    //

    {
        bvect        bv1 { 0, 1 };
        bvect        bv2 { 1, 3 };
        bvect bv1c(bv1);
        bv1c.bit_sub(bv2);

        bvect bv;
        bv.bit_sub(bv1, bv2, bvect::opt_compress);
        int cmp = bv.compare(bv1c);
        assert(cmp == 0);
        struct bvect::statistics st1;
        bv.calc_stat(&st1);
        assert(!st1.bit_blocks);
        assert(st1.gap_blocks == 1);
    }
    
    {
        bvect        bv1 { 0, 1 };
        bvect        bv2;
        for (unsigned i = 0; i < 65536; ++i)
            bv2.set(i);
        
        bvect bv1c(bv1);
        bv1c.bit_sub(bv2);

        bvect bv;
        bv.bit_sub(bv1, bv2, bvect::opt_none); // should detect 0 automatically
        int cmp = bv.compare(bv1c);
        assert(cmp == 0);
        struct bvect::statistics st1;
        bv.calc_stat(&st1);
        assert(!st1.bit_blocks);
        assert(!st1.gap_blocks);
    }
    
    {
        bvect        bv1 { 0, 1 };
        bvect        bv2 { 1 };
        bv1.clear(0); bv1.clear(1);
        bv2.clear(1);
        
        bvect bv;
        bv.bit_or(bv1, bv2, bvect::opt_none); // should detect 0 automatically

        struct bvect::statistics st1;
        bv.calc_stat(&st1);
        assert(!st1.bit_blocks);
        assert(!st1.gap_blocks);
        assert(!st1.ptr_sub_blocks);
    }


    
    {
        bvect        bv1 { 0, 1 };
        bvect        bv2 { 1, 3 };
        bv1.optimize();
        bvect bv1c(bv1);
        bv1c.bit_sub(bv2);

        {
            bvect bv;
            bv.bit_sub(bv1, bv2, bvect::opt_compress);
            int cmp = bv.compare(bv1c);
            assert(cmp == 0);
        }
        bv2.optimize();
        {
            bvect bv;
            bv.bit_sub(bv1, bv2, bvect::opt_compress);
            int cmp = bv.compare(bv1c);
            assert(cmp == 0);
        }
        bv2.clear();
        bv2.invert();
        {
            bvect bv;
            bv.bit_sub(bv1, bv2, bvect::opt_compress);
            assert(!bv.any());
        }
    }
    
    {
        bvect        bvect_full1;
        bvect        bvect_full2;
        bvect_full1.invert();
        bvect_full2.set();
        
        {
        bvect    bv_target_s;
        bv_target_s.bit_sub(bvect_full1, bvect_full2, bvect::opt_none);
        auto b = bv_target_s.none();
        assert(b);
        }
    }

    {
        bvect        bv1;
        bvect        bv2;
        bv1.invert();
        bv2.set();
        bv2.set(0, false);
        
        {
            bvect    bv_s;
            bv_s.bit_sub(bv1, bv2, bvect::opt_none);
            auto cnt = bv_s.count();
            assert(cnt == 1);
            assert(bv1.test(bm::id_max/2));
        }
        {
            bvect    bv_s;
            bv_s.bit_sub(bv2, bv1, bvect::opt_none);
            auto cnt = bv_s.count();
            assert(cnt == 0);
        }
        bv1 -= bv2;
        auto cnt = bv1.count();
        assert(cnt == 1);
        assert(bv1.test(0));
    }

    cout << "----------------------------------- SubOperationTest OK" << endl;
}


static
void XorOperationsTest(bool detailed)
{
    assert(ITERATIONS < BITVECT_SIZE);

    cout << "----------------------------------- XorOperationTest" << endl;
    {

    bvect_mini   bvect_min1(256);
    bvect_mini   bvect_min2(256);
    bvect        bvect_full1;
    bvect        bvect_full2;

    bvect_full1.set_new_blocks_strat(bm::BM_GAP);
    bvect_full2.set_new_blocks_strat(bm::BM_GAP);



    printf("XOR test\n");

    bvect_min1.set_bit(1);
    bvect_min1.set_bit(12);
    bvect_min1.set_bit(13);

    bvect_min2.set_bit(12);
    bvect_min2.set_bit(13);

    bvect_min1.combine_xor(bvect_min2);

    bvect_full1.set_bit(1);
    bvect_full1.set_bit(12);
    bvect_full1.set_bit(13);

    bvect_full2.set_bit(12);
    bvect_full2.set_bit(13);

    bm::id_t predicted_count = bm::count_xor(bvect_full1, bvect_full2);
    bm::id_t predicted_any = bm::any_xor(bvect_full1, bvect_full2);
    if (predicted_any == 0 && predicted_count != 0)
    {
        cout << "Predicted any error!" << endl;
        exit(1);
    }

    bvect    bv_target_s;
    SerializationOperation2Test(&bv_target_s,
                                bvect_full1,
                                bvect_full2,
                                predicted_count,
                                set_COUNT_XOR,
                                set_XOR);


    bvect_full1.bit_xor(bvect_full2);

    bm::id_t count = bvect_full1.count();
    if (count != predicted_count)
    {
        cout << "1.Predicted count error!" << endl;
        exit(1);
    }

    {
        bvect::size_type pos;
        bool f = bvect_full1.find_first_mismatch(bv_target_s, pos);
        if (f)
        {
            cerr << "Mismatch found pos=" << pos << endl;
            assert(0); exit(1);
        }
    }
    CheckVectors(bvect_min1, bvect_full1, 256, detailed);
    CheckVectors(bvect_min1, bv_target_s, 256, detailed);
    CheckCountRange(bvect_full1, 0, 256);
    CheckCountRange(bvect_full1, 128, 256);

    }
    {
        bvect  bvect1;
        bvect_mini  bvect_min1(BITVECT_SIZE);

        bvect  bvect2;
        bvect_mini  bvect_min2(BITVECT_SIZE);


        for (unsigned i = 0; i < 150000; ++i)
        {
            bvect2.set_bit(i);
            bvect_min2.set_bit(i);
        }

        BM_DECLARE_TEMP_BLOCK(tb)
        bvect2.optimize(tb);

        bm::id_t predicted_count = bm::count_xor(bvect1, bvect2);
        bm::id_t predicted_any = bm::any_xor(bvect1, bvect2);
        if (predicted_any == 0 && predicted_count != 0)
        {
            cout << "Predicted any error!" << endl;
            exit(1);
        }

        bvect    bv_target_s;
        SerializationOperation2Test(&bv_target_s,
                                    bvect1,
                                    bvect2,
                                    predicted_count,
                                    set_COUNT_XOR,
                                    set_XOR);

        bvect1.bit_xor(bvect2);
        
        bm::id_t count = bvect1.count();
        if (count != predicted_count)
        {
            cout << "2.Predicted count error!" << endl;
            exit(1);
        }
        
        bvect_min1.combine_xor(bvect_min2);

        {
            bvect::size_type pos;
            bool f = bvect1.find_first_mismatch(bv_target_s, pos);
            if (f)
            {
                cerr << "Mismatch found pos=" << pos << endl;
                assert(0); exit(1);
            }
        }

        CheckVectors(bvect_min1, bvect1, BITVECT_SIZE, detailed);
        CheckVectors(bvect_min1, bv_target_s, BITVECT_SIZE, detailed);
        CheckCountRange(bvect1, 0, BITVECT_SIZE);
    }
    
    {
        bvect bv1;
        bvect bv2;
        bv1.flip();
        bv2.flip();
        bv1.bit_xor(bv2);
        unsigned cnt2 = bv1.count();
        assert(0 == cnt2);
        struct bvect::statistics st;
        bv1.calc_stat(&st);
        auto bcnt = st.bit_blocks + st.gap_blocks;
        assert(bcnt == 0);
    }


    {
        bvect  bvect1;
        bvect_mini  bvect_min1(BITVECT_SIZE);

        bvect  bvect2;
        bvect_mini  bvect_min2(BITVECT_SIZE);


        for (unsigned i = 0; i < 150000; ++i)
        {
            bvect1.set_bit(i);
            bvect_min1.set_bit(i);
        }

        BM_DECLARE_TEMP_BLOCK(tb)
        bvect1.optimize(tb);
        
        bm::id_t predicted_count = bm::count_xor(bvect1, bvect2);
        bm::id_t predicted_any = bm::any_xor(bvect1, bvect2);
        if (predicted_any == 0 && predicted_count != 0)
        {
            cout << "Predicted any error!" << endl;
            exit(1);
        }

        bvect    bv_target_s;
        SerializationOperation2Test(&bv_target_s,
                                    bvect1,
                                    bvect2,
                                    predicted_count,
                                    set_COUNT_XOR,
                                    set_XOR);

        bvect1.bit_xor(bvect2);

        bm::id_t count = bvect1.count();
        if (count != predicted_count)
        {
            cout << "3.Predicted count error!" << endl;
            exit(1);
        }
        
        bvect_min1.combine_xor(bvect_min2);
        {
            bvect::size_type pos;
            bool f = bvect1.find_first_mismatch(bv_target_s, pos);
            if (f)
            {
                cerr << "Mismatch found pos=" << pos << endl;
                assert(0); exit(1);
            }
        }

        CheckVectors(bvect_min1, bvect1, BITVECT_SIZE, detailed);
        CheckVectors(bvect_min1, bv_target_s, BITVECT_SIZE, detailed);
    }


    {
        bvect  bvect1;
        bvect_mini  bvect_min1(BITVECT_SIZE);

        bvect  bvect2;
        bvect_mini  bvect_min2(BITVECT_SIZE);


        for (unsigned i = 0; i < 150000; ++i)
        {
            bvect1.set_bit(i);
            bvect_min1.set_bit(i);
            bvect2.set_bit(i);
            bvect_min2.set_bit(i);
        }

        BM_DECLARE_TEMP_BLOCK(tb)
        bvect1.optimize(tb);
        
        bm::id_t predicted_count = bm::count_xor(bvect1, bvect2);
        bm::id_t predicted_any = bm::any_xor(bvect1, bvect2);
        if (predicted_any == 0 && predicted_count != 0)
        {
            cout << "Predicted any error!" << endl;
            exit(1);
        }

        bvect    bv_target_s;
        SerializationOperation2Test(&bv_target_s,
                                    bvect1,
                                    bvect2,
                                    predicted_count,
                                    set_COUNT_XOR,
                                    set_XOR);

        bvect1.bit_xor(bvect2);

        bm::id_t count = bvect1.count();
        if (count != predicted_count)
        {
            cout << "4.Predicted count error!" << endl;
            cout << count << " " << predicted_count << endl;
            
            exit(1);
        }
        
        bvect_min1.combine_xor(bvect_min2);
        CheckVectors(bvect_min1, bvect1, BITVECT_SIZE, true);
    }



    {

    bvect_mini   bvect_min1(BITVECT_SIZE);
    bvect_mini   bvect_min2(BITVECT_SIZE);
    bvect        bvect_full1;
    bvect        bvect_full2;

    bvect_full1.set_new_blocks_strat(bm::BM_GAP);
    bvect_full2.set_new_blocks_strat(bm::BM_GAP);

    printf("XOR test stage 2.\n");

    FillSets(&bvect_min1, &bvect_full1, 1, BITVECT_SIZE/7, 0);
    FillSets(&bvect_min2, &bvect_full2, 1, BITVECT_SIZE/7, 0);

    bvect_min1.combine_xor(bvect_min2);
    
    bm::id_t predicted_count = bm::count_xor(bvect_full1, bvect_full2);
    bm::id_t predicted_any = bm::any_xor(bvect_full1, bvect_full2);
    if (predicted_any == 0 && predicted_count != 0)
    {
        cout << "Predicted any error!" << endl;
        exit(1);
    }

    bvect    bv_target_s;
    SerializationOperation2Test(&bv_target_s,
                                bvect_full1,
                                bvect_full2,
                                predicted_count,
                                set_COUNT_XOR,
                                set_XOR);


    bvect_full1.bit_xor(bvect_full2);
    
    bm::id_t count = bvect_full1.count();
    if (count != predicted_count)
    {
        cout << "5.Predicted count error!" << endl;
        cout << count << " " << predicted_count << endl;
        print_stat(bvect_full1);
        exit(1);
    }

    {
        bvect::size_type pos;
        bool f = bvect_full1.find_first_mismatch(bv_target_s, pos);
        if (f)
        {
            cerr << "Mismatch found pos=" << pos << endl;
            assert(0); exit(1);
        }
    }

    CheckVectors(bvect_min1, bvect_full1, BITVECT_SIZE/10+10, detailed);
    CheckVectors(bvect_min1, bv_target_s, BITVECT_SIZE/10+10, detailed);
    CheckCountRange(bvect_full1, 0, BITVECT_SIZE/10+10);

    }

    {

    bvect_mini   bvect_min1(BITVECT_SIZE);
    bvect_mini   bvect_min2(BITVECT_SIZE);
    bvect        bvect_full1;
    bvect        bvect_full2;

    bvect_full1.set_new_blocks_strat(bm::BM_BIT);
    bvect_full2.set_new_blocks_strat(bm::BM_BIT);

    cout << "------------------------------" << endl;
    printf("XOR test stage 3.\n");


    FillSets(&bvect_min1, &bvect_full1, 1, BITVECT_SIZE/5, 2);
    FillSets(&bvect_min2, &bvect_full2, 1, BITVECT_SIZE/5, 2);

    bm::id_t predicted_count = bm::count_xor(bvect_full1, bvect_full2);
    bm::id_t predicted_any = bm::any_xor(bvect_full1, bvect_full2);
    if (predicted_any == 0 && predicted_count != 0)
    {
        cout << "Predicted any error!" << endl;
        exit(1);
    }

    bvect    bv_target_s;
    SerializationOperation2Test(&bv_target_s,
                                bvect_full1,
                                bvect_full2,
                                predicted_count,
                                set_COUNT_XOR,
                                set_XOR);

    bvect_min1.combine_xor(bvect_min2);

    bvect_full1.bit_xor(bvect_full2);

    bm::id_t count = bvect_full1.count();
    if (count != predicted_count)
    {
        cout << "6.Predicted count error!" << endl;
        exit(1);
    }


    CheckVectors(bvect_min1, bvect_full1, BITVECT_SIZE);

    BM_DECLARE_TEMP_BLOCK(tb)
    bvect_full1.optimize(tb);
    {
        bvect::size_type pos;
        bool f = bvect_full1.find_first_mismatch(bv_target_s, pos);
        if (f)
        {
            cerr << "Mismatch found pos=" << pos << endl;
            assert(0); exit(1);
        }
    }

    CheckVectors(bvect_min1, bvect_full1, BITVECT_SIZE, detailed);
    CheckVectors(bvect_min1, bv_target_s, BITVECT_SIZE, detailed);
    CheckCountRange(bvect_full1, 0, BITVECT_SIZE);


    }


    cout << "Testing combine_xor" << endl;
    
    {
    
    bvect        bvect_full1;
    bvect        bvect_full2;
    bvect_mini   bvect_min1(BITVECT_SIZE);
    
    bvect_full1.set_new_blocks_strat(bm::BM_GAP);
    bvect_full2.set_new_blocks_strat(bm::BM_GAP);

    unsigned ids[10000];
    unsigned to_add = 10000;
    
    unsigned bn = 0;
    for (unsigned i = 0; i < to_add; ++i)
    {
        ids[i] = bn;
        bvect_full2.set(bn);
        bvect_min1.set_bit(bn);
        bn += 15;
    }
    
    unsigned* first = ids;
    unsigned* last = ids + to_add;
    
    bm::combine_xor(bvect_full1, first, last);

    CheckVectors(bvect_min1, bvect_full1, BITVECT_SIZE);
    
    bm::combine_xor(bvect_full1, first, last);
    if (bvect_full1.count())
    {
        cout << "combine_xor count failed!" << endl;
        exit(1);
    }
    
    }

    {
    
    bvect        bvect_full1;
    bvect        bvect_full2;
    bvect_mini   bvect_min1(BITVECT_SIZE);
    
    bvect_full1.set_new_blocks_strat(bm::BM_GAP);
    bvect_full2.set_new_blocks_strat(bm::BM_GAP);

    unsigned ids[10000]={0,};
    unsigned to_add = 10000;
    
    for (unsigned i = 0; i < to_add; i+=100)
    {
        ids[i] = i;
        bvect_full2.set(i);
        bvect_min1.set_bit(i);
    }
    unsigned* first = ids;
    unsigned* last = ids + to_add;
    
    bm::combine_xor(bvect_full1, first, last);

    CheckVectors(bvect_min1, bvect_full1, BITVECT_SIZE);
    
    bm::combine_xor(bvect_full1, first, last);
    if (bvect_full1.count())
    {
        cout << "combine_xor count failed!" << endl;
        exit(1);
    }
    
    }

    
    {
    unsigned ids[] = {0, 65536, 65535, 65535*3, 65535*2, 10};
    unsigned to_add = sizeof(ids)/sizeof(unsigned);
    bvect        bvect_full1;
    bvect        bvect_full2;    
    bvect_mini   bvect_min1(BITVECT_SIZE);

    bvect_full1.set_new_blocks_strat(bm::BM_BIT);
    bvect_full2.set_new_blocks_strat(bm::BM_BIT);
    
    unsigned bn = 0;
    for (unsigned i = 0; i < to_add; ++i)
    {
        ids[i] = bn;
        bvect_full2.set(bn);
        bvect_min1.set_bit(bn);
        bn += 15;
    }
    
    unsigned* first = ids;
    unsigned* last = ids + to_add;
    
    bm::combine_xor(bvect_full1, first, last);
    CheckVectors(bvect_min1, bvect_full1, BITVECT_SIZE);

    bm::combine_xor(bvect_full1, first, last);
    if (bvect_full1.count())
    {
        cout << "combine_xor count failed!" << endl;
        exit(1);
    }
    }
    
    
    {
    unsigned ids[] = {0, 65536, 65535, 65535*3, 65535*2, 10};
    unsigned to_add = sizeof(ids)/sizeof(unsigned);
    bvect        bvect_full1;
    bvect        bvect_full2;    
    bvect_mini   bvect_min1(BITVECT_SIZE);

    bvect_full1.set_new_blocks_strat(bm::BM_GAP);
    bvect_full2.set_new_blocks_strat(bm::BM_GAP);
    
    unsigned bn = 0;
    for (unsigned i = 0; i < to_add; ++i)
    {
        ids[i] = bn;
        bvect_full2.set(bn);
        bvect_min1.set_bit(bn);
        bn += 15;
    }
    
    unsigned* first = ids;
    unsigned* last = ids + to_add;
    
    bm::combine_xor(bvect_full1, first, last);
    CheckVectors(bvect_min1, bvect_full1, BITVECT_SIZE);

    bm::combine_xor(bvect_full1, first, last);
    if (bvect_full1.count())
    {
        cout << "combine_xor count failed!" << endl;
        exit(1);
    }
    }
    
    
    // ------------------------------------------
    // 2-way XOR
    //
    {
        bvect        bv1 { 0, 1 };
        bvect        bv2;
        bv2.bit_xor(bv1, bv2, bvect::opt_compress);
        int cmp = bv1.compare(bv2);
        assert(cmp == 0);
    }

    {
        bvect        bv1 { 0, 1 };
        bvect        bv2 { 2, 3 };
        bvect bv1c(bv1);
        bv1c.bit_xor(bv2);

        bvect bv;
        bv.bit_xor(bv1, bv2, bvect::opt_compress);
        int cmp = bv.compare(bv1c);
        assert(cmp == 0);
        struct bvect::statistics st1;
        bv.calc_stat(&st1);
        assert(!st1.bit_blocks);
        assert(st1.gap_blocks == 1);
    }
    
    {
        bvect        bv1;
        bvect        bv2;
        for (unsigned i = 2; i < 65536; ++i)
        {
            bv1.set(i);
            bv2.set(i);
        }
        
        bvect bv1c(bv1);
        bv1c.bit_xor(bv2);

        bvect bv;
        bv.bit_xor(bv1, bv2, bvect::opt_none); // should detect 0 automatically
        int cmp = bv.compare(bv1c);
        assert(cmp == 0);
        struct bvect::statistics st1;
        bv.calc_stat(&st1);
        assert(!st1.bit_blocks);
        assert(!st1.gap_blocks);
    }
    
    {
        bvect        bv1 { 0, 1 };
        bvect        bv2 { 1 };
        bv1.clear(0); bv1.clear(1);
        bv2.clear(1);
        
        bvect bv;
        bv.bit_xor(bv1, bv2, bvect::opt_none); // should detect FULL automatically

        struct bvect::statistics st1;
        bv.calc_stat(&st1);
        assert(!st1.bit_blocks);
        assert(!st1.gap_blocks);
        assert(!st1.ptr_sub_blocks);
    }


    
    {
        bvect        bv1 { 0, 1 };
        bvect        bv2 { 2, 3 };
        bv2.optimize();
        bvect bv1c(bv1);
        bv1c.bit_xor(bv2);

        {
            bvect bv;
            bv.bit_xor(bv1, bv2, bvect::opt_compress);
            int cmp = bv.compare(bv1c);
            if (cmp != 0)
            {
                DetailedCompareBVectors(bv, bv1c);
            }
            assert(cmp == 0);
        }
        bv1.optimize();
        {
            bvect bv;
            bv.bit_xor(bv1, bv2, bvect::opt_compress);
            int cmp = bv.compare(bv1c);
            assert(cmp == 0);
        }
        bv1.clear();
        bv2.clear();
        bv2.invert();
        {
            bvect bv;
            bv.bit_xor(bv1, bv2, bvect::opt_compress);
            int cmp = bv.compare(bv2);
            assert(cmp == 0);
        }
    }

    {
        bvect        bvect_full1;
        bvect        bvect_full2;
        bvect_full1.invert();
        bvect_full2.set();
        
        {
        bvect    bv_target_s;
        bv_target_s.bit_xor(bvect_full1, bvect_full2, bvect::opt_none);
        auto b = bv_target_s.none();
        assert(b);
        }
    }

    {
        bvect        bv1;
        bvect        bv2;
        bv1.invert();
        bv2.set();
        bv2.set(bm::id_max/2, false);
        
        {
            bvect    bv_s;
            bv_s.bit_xor(bv1, bv2, bvect::opt_none);
            auto cnt = bv_s.count();
            assert(cnt == 1);
            assert(bv_s.test(bm::id_max/2));
        }
        {
            bvect    bv_s;
            bv_s.bit_xor(bv2, bv1, bvect::opt_none);
            auto cnt = bv_s.count();
            assert(cnt == 1);
            assert(bv_s.test(bm::id_max/2));
        }
        bv1.bit_xor(bv2);
        auto cnt = bv1.count();
        assert(cnt == 1);
        assert(bv1.test(bm::id_max/2));
    }



}

static
void ComparisonTest()
{
    cout << "-------------------------------------- ComparisonTest" << endl;

    bvect_mini   bvect_min1(BITVECT_SIZE);
    bvect_mini   bvect_min2(BITVECT_SIZE);
    bvect        bvect_full1;
    bvect        bvect_full2;
    int res1, res2;

    bvect_full1.set_bit(31); 
    bvect_full2.set_bit(63); 

    res1 = bvect_full1.compare(bvect_full2);
    if (res1 != 1)
    {
        printf("Comparison test failed 1\n");
        exit(1);
    }

    bvect_full1.clear();
    bvect_full2.clear();

    bvect_min1.set_bit(10);
    bvect_min2.set_bit(10);

    bvect_full1.set_bit(10);
    bvect_full2.set_bit(10);

    res1 = bvect_min1.compare(bvect_min2);
    res2 = bvect_full1.compare(bvect_full2);

    if (res1 != res2)
    {
        printf("Comparison test failed 1\n");
        exit(1);
    }

    printf("Comparison 2.\n");

    bvect_min1.set_bit(11);
    bvect_full1.set_bit(11);

    res1 = bvect_min1.compare(bvect_min2);
    res2 = bvect_full1.compare(bvect_full2);

    if (res1 != res2 && res1 != 1)
    {
        printf("Comparison test failed 2\n");
        exit(1);
    }

    res1 = bvect_min2.compare(bvect_min1);
    res2 = bvect_full2.compare(bvect_full1);

    if (res1 != res2 && res1 != -1)
    {
        printf("Comparison test failed 2.1\n");
        exit(1);
    }

    printf("Comparison 3.\n");

    BM_DECLARE_TEMP_BLOCK(tb)
    bvect_full1.optimize(tb);

    res1 = bvect_min1.compare(bvect_min2);
    res2 = bvect_full1.compare(bvect_full2);

    if (res1 != res2 && res1 != 1)
    {
        printf("Comparison test failed 3\n");
        exit(1);
    }

    res1 = bvect_min2.compare(bvect_min1);
    res2 = bvect_full2.compare(bvect_full1);

    if (res1 != res2 && res1 != -1)
    {
        printf("Comparison test failed 3.1\n");
        exit(1);
    }

    printf("Comparison 4.\n");

    bvect_full2.optimize();

    res1 = bvect_min1.compare(bvect_min2);
    res2 = bvect_full1.compare(bvect_full2);

    if (res1 != res2 && res1 != 1)
    {
        printf("Comparison test failed 4\n");
        exit(1);
    }

    res1 = bvect_min2.compare(bvect_min1);
    res2 = bvect_full2.compare(bvect_full1);

    if (res1 != res2 && res1 != -1)
    {
        printf("Comparison test failed 4.1\n");
        exit(1);
    }

    printf("Comparison 5.\n");

    unsigned i;
    for (i = 0; i < 65536; ++i)
    {
        bvect_full1.set_bit(i);
    }

    res1 = bvect_min1.compare(bvect_min2);
    res2 = bvect_full1.compare(bvect_full2);

    if (res1 != res2 && res1 != 1)
    {
        printf("Comparison test failed 5\n");
        exit(1);
    }

    bvect_full1.optimize();

    res1 = bvect_min2.compare(bvect_min1);
    res2 = bvect_full2.compare(bvect_full1);

    if (res1 != res2 && res1 != -1)
    {
        printf("Comparison test failed 5.1\n");
        exit(1);
    }

}


static
void BvectorFindFirstDiffTest()
{
    cout << "-------------------------------------- BvectorFindFirstDiffTest" << endl;

    // empty test
    {
        bvect bv1, bv2;
        TestFindDiff(bv1, bv2);

        bv1.set(0);
        TestFindDiff(bv1, bv2);
        bv2.set(0);
        TestFindDiff(bv1, bv2);
    }
    // test GAP bits
    {
        bvect bv1(bm::BM_GAP), bv2(bm::BM_GAP);

        bv1.set_range(10, 15);
        bv2.set_range(10, 15);

        TestFindDiff(bv1, bv2);
    }
    {
        bvect bv1(bm::BM_GAP), bv2(bm::BM_GAP);

        bv1.set_range(10, 15);
        bv2.set_range(10, 12);

        TestFindDiff(bv1, bv2);
    }
    {
        bvect bv1(bm::BM_GAP), bv2(bm::BM_GAP);

        bv1.set_range(bm::id_max32/2 - 10, bm::id_max32/2 + 15);
        bv2.set_range(bm::id_max32/2 - 10, bm::id_max32/2 + 12);

        TestFindDiff(bv1, bv2);
    }
    // test GAP-bit mix
    {
        bvect bv1(bm::BM_GAP), bv2;

        bv1.set_range(10, 15);
        bv2.set_range(10, 12);

        TestFindDiff(bv1, bv2);
    }
    {
        bvect bv1(bm::BM_GAP), bv2;

        bv1.set_range(bm::id_max32/2 - 10, bm::id_max32/2 + 15);
        bv2.set_range(bm::id_max32/2 - 10, bm::id_max32/2 + 12);

        TestFindDiff(bv1, bv2);
    }

    // test inverted
    {
        bvect bv1, bv2;

        bv1.invert();
        TestFindDiff(bv1, bv2);
        bv2.invert();
        TestFindDiff(bv1, bv2);
        bv2[123456] = false;
        TestFindDiff(bv1, bv2);
    }


    // test bits far
    {
        bvect bv1, bv2;
        bv1.set(bm::id_max32/2);
        TestFindDiff(bv1, bv2);
        bv2.set(bm::id_max32/2);
        TestFindDiff(bv1, bv2);
        bv2.set(bm::id_max32/2+1);
        TestFindDiff(bv1, bv2);
        bv1.optimize();
        TestFindDiff(bv1, bv2);
        bv2.optimize();
        TestFindDiff(bv1, bv2);
    }
    {
        bvect bv1, bv2;
        bv1.set(bm::id_max-1);
        TestFindDiff(bv1, bv2);
        bv1.optimize();
        TestFindDiff(bv1, bv2);
        bv2.set(bm::id_max-1);
        TestFindDiff(bv1, bv2);
        bv2.optimize();
        TestFindDiff(bv1, bv2);
    }

    // test FULL blocks
    {
        bvect bv1, bv2;
        bv1.set_range(0, bm::id_max32/2);
        TestFindDiff(bv1, bv2);
        bv2.set_range(0, bm::id_max32/2);
        TestFindDiff(bv1, bv2);

        bv1[bm::id_max32/2 - 100] = false;
        TestFindDiff(bv1, bv2);
        bv1.optimize();
        TestFindDiff(bv1, bv2);

        bv2[bm::id_max32/2 - 100] = false;
        TestFindDiff(bv1, bv2);
        bv2.optimize();
        TestFindDiff(bv1, bv2);
    }

    cout << "-------------------------------------- BvectorFindFirstDiffTest OK" << endl;
}

static
void DesrializationTest2()
{
   bvect  bvtotal;
   unsigned size = BITVECT_SIZE - 10;
   BM_DECLARE_TEMP_BLOCK(tb)


   bvect  bv1;
   bvect  bv2;
   unsigned i;
   for (i = 10; i < 165536; i+=2)
   {
      bv1.set_bit(i);
   }

   bv1.optimize(tb);
   print_stat(bv1);

   struct bvect::statistics st1;
   bv1.calc_stat(&st1);

   std::vector<unsigned char> sermemv(st1.max_serialize_mem);
   
   size_t slen2 = bm::serialize(bv1, sermemv.data(), tb);
   assert(slen2);
   slen2 = 0;

   bm::deserialize(bvtotal, sermemv.data());
    bvect  bv_target_s;
    operation_deserializer<bvect> od;
    od.deserialize(bv_target_s, sermemv.data(), 0, set_OR);

   bvtotal.optimize(tb);
   int res = bvtotal.compare(bv_target_s);
   if (res != 0)
   {
       cout << "Operation deserialization error 1" << endl;
       exit(1);
   }

   for (i = 55000; i < 165536; ++i)
   {
      bv2.set_bit(i);
   }
   bv2.optimize();
   print_stat(bv2);

   struct bvect::statistics st2;
   bv2.calc_stat(&st2);

   std::vector<unsigned char> sermemv2(st2.max_serialize_mem);

   size_t slen = bm::serialize(bv2, sermemv2.data());
   assert(slen);
   slen = 0;

   bm::deserialize(bvtotal, sermemv2.data());
   print_stat(bvtotal);
   od.deserialize(bv_target_s, sermemv2.data(), 0, set_OR);
    res = bvtotal.compare(bv_target_s);
    if (res != 0)
    {
        cout << "Operation deserialization error 2" << endl;
        assert(0);
        exit(1);
    }

   bm::deserialize(bvtotal, sermemv2.data());
   bm::deserialize(bvtotal, sermemv.data());

    od.deserialize(bv_target_s,
                   sermemv2.data(),
                   0,
                   set_OR);
    od.deserialize(bv_target_s,
                   sermemv2.data(),
                   0,
                   set_OR);

    res = bvtotal.compare(bv_target_s);
    if (res != 0)
    {
        cout << "Deserialization test failed! 3" << endl;
        exit(1);
    }


   bvtotal.clear();
   bv_target_s.clear(false);

   int clcnt = 0;

   unsigned repetitions = 25;
   for (i = 0; i < repetitions; ++i)
   {
        cout << endl << endl << "Deserialization STEP " << i << endl;

        bvect_mini*   bvect_min1= new bvect_mini(size);
        bvect*        bvect_full1= new bvect();

        FillSetsRandomMethod(bvect_min1, bvect_full1, 1, size, 1);

       struct bvect::statistics st;
       bvect_full1->calc_stat(&st);

       std::vector<unsigned char> sermemv1(st.max_serialize_mem);
       slen = bm::serialize(*bvect_full1, sermemv1.data(), tb);

       std::vector<unsigned char> smemv(slen);
       ::memcpy(smemv.data(), sermemv1.data(), slen);

        bm::deserialize(bvtotal, smemv.data());
       
        {
            bvect bv_c;
            bm::deserialize(bv_c, smemv.data());
            res = bv_c.compare(*bvect_full1);
            assert(res == 0);
            
            bvect bv3;
            od.deserialize(bv3,
                           smemv.data(),
                           0,
                           set_OR);
            res = bv3.compare(*bvect_full1);
            assert(res == 0);
        }
       
        od.deserialize(bv_target_s,
                       smemv.data(),
                       0,
                       set_OR);
        res = bvtotal.compare(bv_target_s);
        if (res != 0)
        {
            res = bvtotal.compare(bv_target_s);
            
            unsigned bit_idx = bv_target_s.get_first();
            cout << bit_idx << " " << bv_target_s.get_next(bit_idx) << endl;;
            print_stat(*bvect_full1);
            print_stat(bv_target_s);
            cout << "Operation deserialization error 2" << endl;
            assert(0); exit(1);
        }

       bvtotal.optimize(tb);
       bv_target_s.optimize(tb);

       if (++clcnt == 5)
       {
          clcnt = 0;
          bvtotal.clear();
          bv_target_s.clear();
       }

       delete bvect_min1;
       delete bvect_full1;

   } // for i

}



// ------------------------------------------------------------------------

static
bool agg_shift_right_and(bm::aggregator<bvect>& agg,
                         bvect& bv_target,
                         const bvect* bv, ...)
{
    va_list args;
    va_start(args, bv);
    agg.add(bv);
    
    for (int i = 0; true; ++i)
    {
        const bvect* bv_arg = (const bvect*)va_arg(args, void*);
        if (!bv_arg)
            break;
        agg.add(bv_arg);
    }
    va_end(args);
    
    agg.combine_shift_right_and(bv_target);
    agg.reset();
    return bv_target.any();
}


static
void AggregatorTest()
{
  cout << "---------------------------- Aggregator Test" << endl;

    bvect* bv_arr[128] = { 0, };
    bvect* bv_arr2[128] = { 0, };

    cout << "OR tests..." << endl;
    {
        bm::aggregator<bvect> agg;
        agg.set_optimization();

        {
            bvect bv_target;
            bvect bv1 { 1, 2, 3};
            bvect bv2 { 0, 4, 5};
            
            agg.add(&bv1);
            agg.add(&bv2);
            
            agg.combine_or(bv_target);
            agg.reset();
            
            struct bvect::statistics st;
            bv_target.calc_stat(&st);
            assert (st.gap_blocks == 1);
            assert (st.bit_blocks == 0);
            
            unsigned bc = bv_target.count();
            assert(bc == 6);
        }
        
        // FULL block optimization test
        {
            bvect bv_target;
            bvect bv1;
            bvect bv2;
            bv1.set_range(0, 256);
            bv2.set_range(256, 65535);

            agg.reset();
            agg.add(&bv1);
            agg.add(&bv2);
            
            agg.combine_or(bv_target);
            
            struct bvect::statistics st;
            bv_target.calc_stat(&st);
            assert (st.gap_blocks == 0);
            assert (st.bit_blocks == 0);
            
            unsigned bc = bv_target.count();
            assert(bc == 65536);

            bv_target.set(1);
            agg.reset();
            agg.add(&bv1);
            agg.add(&bv2);
            
            agg.combine_and(bv_target);
            bv_target.calc_stat(&st);
            assert (st.gap_blocks == 1);
            assert (st.bit_blocks == 0);
            bc = bv_target.count();
            assert(bc == 1);
        }

        // 0-block optimization test
        {
            bvect bv_target;
            bvect bv1 { 1 };
            bvect bv2 { 5 };
            
            bv1.set(1, false); bv2.set(5, false);
            
            agg.reset();
            agg.add(&bv1);
            agg.add(&bv2);
            
            agg.combine_or(bv_target);
            agg.reset();
            
            struct bvect::statistics st;
            bv_target.calc_stat(&st);
            assert (st.gap_blocks == 0);
            assert (st.bit_blocks == 0);
            
            unsigned bc = bv_target.count();
            assert(bc == 0);
        }
    }

    bm::aggregator<bvect> agg;
    agg.set_optimization();

    cout << "AND-SUB tests..." << endl;
    {
        bvect bv1, bv2, bv3;
        bvect bv_empty;
        bvect bv_control;

        bv_arr[0] = &bv1;
        agg.combine_or(bv3, bv_arr, 1);
        assert(bv3.count()==0);

        bv3[100] = true;
        bv1.invert();
        bv_control.invert();
        bv_arr[0] = &bv1;
        agg.combine_or(bv3, bv_arr, 1);

        int res = bv_control.compare(bv3);
        assert(res == 0);

        bv_arr[0] = &bv1;
        bv_arr[1] = &bv2;
        agg.combine_or(bv3, bv_arr, 2);

        res = bv_control.compare(bv3);
        assert(res == 0);

        bv2[1000000] = true;
        bv_arr[0] = &bv1;
        bv_arr[1] = &bv2;
        agg.combine_or(bv3, bv_arr, 2);
        res = bv_control.compare(bv3);
        assert(res == 0);
    }
    
    {
        bvect bv1, bv2, bv3;
        bvect bv_empty;
        bvect bv_control;
        int res;

        bv1.invert();
        bv_control.invert();
        bv_arr[0] = &bv1;
        agg.combine_and(bv3, bv_arr, 1);
        res = bv_control.compare(bv3);
        assert(res == 0);
        
        bv2.invert();
        bv2.set_range(100000, 100100);
        bv_arr[0] = &bv1;
        bv_arr[1] = &bv2;
        bv_control.set_range(100000, 100100);
        agg.combine_and(bv3, bv_arr, 2);
        res = bv_control.compare(bv3);
        agg.combine_and_sub(bv3, bv_arr, 2, 0, 0, false);
        res = bv_control.compare(bv3);
        assert(res == 0);

    }

    {
        bvect bv1, bv2, bv3, bv4;
        bvect bv_empty;
        bvect bv_control;
        int res;

        bv1.invert();
        bv2.set_range(200000, 300000);
        bv3.set_range(200000, 300000);
        
        bv_control.set_range(200000, 300000);
        
        bv_arr[0] = &bv1;
        bv_arr[1] = &bv2;
        bv_arr[2] = &bv3;
        
        agg.combine_and(bv4, bv_arr, 3);
        res = bv_control.compare(bv4);
        assert(res == 0);
        agg.combine_and_sub(bv4, bv_arr, 3, 0, 0, false);
        res = bv_control.compare(bv3);
        assert(res == 0);
    }

    {
        bvect bv1, bv2, bv3;
        bvect bv_empty;
        bvect bv_control;
        int res;
        
        bv_arr[0] = &bv1;
        agg.combine_and(bv3, bv_arr, 1);
        assert(bv3.count()==0);
        
        bv1[100] = true;
        bv_arr[0] = &bv1;
        agg.combine_and(bv3, bv_arr, 1);
        assert(bv3.count()==1);
        assert(bv3[100] == true);
        
        bv1[100] = true;
        bv2.invert();
        bv_arr[0] = &bv1;
        bv_arr[1] = &bv2;
        agg.combine_and(bv3, bv_arr, 2);
        assert(bv3.count()==1);
        assert(bv3[100] == true);
        
        bv1.clear();
        bv1.invert();
        bv_control.invert();
        bv_arr[0] = &bv1;
        bv_arr[1] = &bv2;
        agg.combine_and(bv3, bv_arr, 2);
        res = bv_control.compare(bv3);
        assert(res == 0);
        agg.combine_and_sub(bv3, bv_arr, 2, 0, 0, false);
        res = bv_control.compare(bv3);
        assert(res == 0);
    }


    //  ---------------------------
    {
        bvect bv1, bv2, bv3, bv4;
        bvect bv_empty;
        bvect bv_control;
        
        bv1[100] = true;
        bv1[100000] = true;
        bv2[100] = true;
        bv2[100000] = true;
        bv3[100000] = true;

        bv_arr[0] = &bv1;
        bv_arr[1] = &bv2;
        bv_arr2[0] = &bv3;

        agg.combine_and_sub(bv4, bv_arr, 2, bv_arr2, 1, false);
        assert(bv4.count()==1);
        assert(bv4.test(100));
        
        bv3.optimize();
        agg.combine_and_sub(bv4, bv_arr, 2, bv_arr2, 1, false);
        assert(bv4.count()==1);
        assert(bv4.test(100));
        
        bv1.optimize();
        agg.combine_and_sub(bv4, bv_arr, 2, bv_arr2, 1, false);
        assert(bv4.count()==1);
        assert(bv4.test(100));

        bv2.optimize();
        agg.combine_and_sub(bv4, bv_arr, 2, bv_arr2, 1, false);
        assert(bv4.count()==1);
        assert(bv4.test(100));
    }

    {
        bvect bv1, bv2, bv3, bv4;
        bvect bv_empty;
        bvect bv_control;
        
        bv1[100] = true;
        bv1[100000] = true;
        bv2[100] = true;
        bv2[100000] = true;
        
        bv3.invert();

        bv_arr[0] = &bv1;
        bv_arr[1] = &bv2;
        bv_arr2[0] = &bv3;

        agg.combine_and_sub(bv4, bv_arr, 2, bv_arr2, 1, false);
        assert(bv4.count()==0);
        assert(!bv4.any());
    }
    
    // SHIFT-R_AND
    
    cout << "SHIFT-R-AND tests..." << endl;
    
    {
    bvect bv0, bv1, bv2;
    bv1[0] = true;
    bv1[65535]=true;
    
    bv2[1]=true;
    bv2[65536]=true;
    
    agg.add(&bv1); agg.add(&bv2);
    
    agg.combine_shift_right_and(bv0);
    agg.reset();
    bool any = bv0.any();
    
///    bool any = agg.shift_right_and(bv1, bv2);
    assert(any);
    assert(bv0.count()==2);
    assert(bv0.test(1));
    assert(bv0.test(65536));
    }

    {
    bvect bv0, bv1, bv2;
    bv1[0] = true;
    bv1[65535]=true;

    bv2[0]=true;
    bv2[65535]=true;
    
    bool any = agg_shift_right_and(agg, bv0, &bv1, &bv2, 0);
    assert(!any);
    assert(bv0.count()==0);
    }


    {
    bvect bv0, bv1, bv2;
    bv1[0] = true;
    bv1[65535]=true;
    bv1.optimize();

    bv2[1]=true;
    bv2[65536]=true;
    bv2.optimize();

    agg_shift_right_and(agg, bv0, &bv1, &bv2, 0);
    assert(bv0.count()==2);
    assert(bv0.test(1));
    assert(bv0.test(65536));
    }


    {
    bvect bv0, bv1, bv2;
    bv1[65535]=true;
    
    bv2[65536]=true;
    bv2.optimize();
    
    bool any = agg_shift_right_and(agg, bv0, &bv1, &bv2, 0);
    assert(bv0.count()==1);
    assert(bv0.test(65536));
    assert(any);
    struct bvect::statistics st1;
    bv0.calc_stat(&st1);
    auto bcnt = st1.bit_blocks + st1.gap_blocks;
    assert(bcnt == 1);
    }

    
    {
    bvect bv0, bv1, bv2;
    bv1[0] = true;
    bv1[65535]=true;

    bv2.invert();
    
    agg_shift_right_and(agg, bv0, &bv1, &bv2, 0);
    assert(bv0.count()==2);
    assert(bv0.test(1));
    assert(bv0.test(65536));
    }
    
    {
    bvect bv0, bv1, bv2;
    bvect bv1c, bv2c;
    bv1.invert();
    bv2.invert();
    bv1c.invert();
    bv2c.invert();

    bv1c.shift_right();
    bv1c &= bv2c;

    bool any = agg_shift_right_and(agg, bv0, &bv1, &bv2, 0);

    assert(any);
    assert(!bv0.test(0));
    assert(!bv1c.test(0));

    struct bvect::statistics st1;
    bv0.calc_stat(&st1);
    auto bcnt = st1.bit_blocks + st1.gap_blocks;
    cout << bcnt << endl;
    assert(bcnt == 2);
    
    assert(bv0.count()==bv1c.count());
    auto cmp = bv1c.compare(bv0);
    assert(cmp==0);
    }
    
    {
    bvect bv0, bv1, bv2;
    bv1.set_range(0, 65536*4);
    bv2.set_range(0, 65536*4);
    
    bool any = agg_shift_right_and(agg, bv0, &bv1, &bv2, 0);
    
    assert(any);
    assert(!bv0.test(0));
    assert(bv0.count() == 65536*4);

    struct bvect::statistics st1;
    bv0.calc_stat(&st1);
    auto bcnt = st1.bit_blocks + st1.gap_blocks;
    assert(bcnt == 2); // TODO: not critical (optimization) needs a fix
    }

    {
    bvect bv0, bv1, bv2;
    bv1.set_range(0, 65536*4);
    bv2.set_range(0, 65536*2);
    
    bool any = agg_shift_right_and(agg, bv0, &bv1, &bv2, 0);
    
    assert(any);
    assert(!bv0.test(0));
    assert(bv0.count() == 65536*2);

    struct bvect::statistics st1;
    bv0.calc_stat(&st1);
    auto bcnt = st1.bit_blocks + st1.gap_blocks;
    assert(bcnt == 2);
    }


    {
    bvect bv0, bv1, bv2;
    bv1.set_range(0, 65536*4);
    bv2.set_range(65536, 65536+10);
    
    bool any = agg_shift_right_and(agg, bv0, &bv1, &bv2, 0);
    
    assert(any);
    assert(!bv0.test(0));
    cout << bv0.count() << endl;
    assert(bv0.count() == 11);

    struct bvect::statistics st1;
    bv0.calc_stat(&st1);
    auto bcnt = st1.bit_blocks + st1.gap_blocks;
    assert(bcnt == 1);
    }


  cout << "---------------------------- Aggregator Test OK" << endl;
}


static
void StressTestAggregatorOR(unsigned repetitions)
{
  cout << "---------------------------- Aggregator OR Stress Test" << endl;
   unsigned size = BITVECT_SIZE - 10;
   bvect bv_target1, bv_target2;


    unsigned i;
    for (i = 0; i < repetitions; ++i)
    {
        int opt = rand() % 2;
        cout << endl << " - - - - - - - - - - - - AGG OR STRESS STEP " << i << endl;;
        
        switch (rand() % 3)
        {
        case 0:
            size = BITVECT_SIZE / 10;
            break;
        case 1:
            size = BITVECT_SIZE / 2;
            break;
        default:
            size = BITVECT_SIZE - 10;
            break;
        } // switch
        
        unsigned start1 = 0;
        switch (rand() % 3)
        {
        case 1:
            start1 += size / 5;
            break;
        default:
            break;
        }

        unsigned start2 = 0;
        switch (rand() % 3)
        {
        case 1:
            start2 += size / 5;
            break;
        default:
            break;
        }

        bvect_mini   bvect_min1(size);
        bvect bv0, bv1, bv2, bv3, bv4, bv5, bv6, bv7, bv8, bv9;

        // 0 skipped
        FillSetsRandomMethod(&bvect_min1, &bv1, start1, size, opt);
        FillSetsRandomMethod(&bvect_min1, &bv2, start2, size, opt);
        // 3 skipped
        FillSetsRandomMethod(&bvect_min1, &bv5, start1, size, opt);
        FillSetsRandomMethod(&bvect_min1, &bv6, start2, size, opt);
        FillSetsRandomMethod(&bvect_min1, &bv7, start1, size, opt);
        FillSetsRandomMethod(&bvect_min1, &bv8, start2, size, opt);
        FillSetsRandomMethod(&bvect_min1, &bv9, start2, size, opt);
        
        bm::aggregator<bvect> agg;
        agg.set_optimization();
        
        bvect* agg_list[32] = {0, };
        
        agg_list[0] = &bv0;
        agg_list[1] = &bv1;
        agg_list[2] = &bv2;
        agg_list[3] = &bv3;
        agg_list[4] = &bv4;
        agg_list[5] = &bv5;
        agg_list[6] = &bv6;
        agg_list[7] = &bv7;
        agg_list[8] = &bv8;
        agg_list[9] = &bv9;
        
        unsigned cnt = 10;
        agg.combine_or(bv_target1, agg_list, cnt);
        agg.combine_or_horizontal(bv_target2, agg_list, cnt);

        int res = bv_target1.compare(bv_target2);
        if (res!=0)
        {
            cerr << "Error: Aggregator OR check failed!" << endl;
            DetailedCompareBVectors(bv_target1, bv_target2);
            exit(1);
        }
        for (unsigned j = 1; j < cnt; ++j)
        {
            agg.combine_or(bv_target1, agg_list, j);
            agg.combine_or_horizontal(bv_target2, agg_list, j);
            res = bv_target1.compare(bv_target2);
            if (res!=0)
            {
                cerr << "Error: Aggregator OR check failed! 1.laddder step = "
                     << j << endl;
                exit(1);
            }
        }
        
        for (unsigned j = 0; j < cnt; ++j)
        {
            agg.combine_or(bv_target1, agg_list+j, cnt-j);
            agg.combine_or_horizontal(bv_target2, agg_list+j, cnt-j);
            res = bv_target1.compare(bv_target2);
            if (res!=0)
            {
                cerr << "Error: Aggregator OR check failed! 2.laddder step = "
                     << j << endl;
                exit(1);
            }
        }


    } // for i

  cout << "---------------------------- Aggregator OR Stress Test OK" << endl;
}

static
void StressTestAggregatorAND(unsigned repetitions)
{
  cout << "---------------------------- Aggregator AND Stress Test" << endl;
   unsigned size = BITVECT_SIZE - 10;


    unsigned i;
    for (i = 0; i < repetitions; ++i)
    {
        int opt = rand() % 2;
        cout << endl << " - - - - - - - - - - - - AGG AND STRESS STEP " << i << endl;;
        
        switch (rand() % 3)
        {
        case 0:
            size = BITVECT_SIZE / 10;
            break;
        case 1:
            size = BITVECT_SIZE / 2;
            break;
        default:
            size = BITVECT_SIZE - 10;
            break;
        } // switch
        
        unsigned start1 = 0;
        switch (rand() % 3)
        {
        case 1:
            start1 += size / 5;
            break;
        default:
            break;
        }

        unsigned start2 = 0;
        switch (rand() % 3)
        {
        case 1:
            start2 += size / 5;
            break;
        default:
            break;
        }

        bvect_mini   bvect_min1(size);
        bvect bv0, bv1, bv2, bv3, bv4, bv5, bv6, bv7, bv8, bv9;

        // 0 skipped
        FillSetsRandomMethod(&bvect_min1, &bv1, start1, size, opt);
        FillSetsRandomMethod(&bvect_min1, &bv2, start2, size, opt);
        // 3 skipped
        FillSetsRandomMethod(&bvect_min1, &bv5, start1, size, opt);
        FillSetsRandomMethod(&bvect_min1, &bv6, start2, size, opt);
        FillSetsRandomMethod(&bvect_min1, &bv7, start1, size, opt);
        FillSetsRandomMethod(&bvect_min1, &bv8, start2, size, opt);
        FillSetsRandomMethod(&bvect_min1, &bv9, start2, size, opt);
        
        bm::aggregator<bvect> agg;
        agg.set_optimization();
        
        bvect* agg_list[32] = {0, };
        
        agg_list[0] = &bv0;
        agg_list[1] = &bv1;
        agg_list[2] = &bv2;
        agg_list[3] = &bv3;
        agg_list[4] = &bv4;
        agg_list[5] = &bv5;
        agg_list[6] = &bv6;
        agg_list[7] = &bv7;
        agg_list[8] = &bv8;
        agg_list[9] = &bv9;
        
        bvect bv_target1, bv_target2, bv_target3, bv_target4;
        bvect bv_empty;
        
        unsigned cnt = 10;
        agg.combine_and_sub(bv_target3, agg_list, cnt, 0, 0, false);
        agg.combine_and(bv_target1, agg_list, cnt);
        agg.combine_and_horizontal(bv_target2, agg_list, cnt);
        agg.combine_and_sub(bv_empty, agg_list, cnt, agg_list, cnt, false);

        int res = bv_target1.compare(bv_target2);
        if (res!=0)
        {
            cerr << "Error: Aggregator AND check failed!" << endl;
            assert(0);
            exit(1);
        }
        res = bv_target3.compare(bv_target1);
        if (res!=0)
        {
            cerr << "Error: Aggregator AND-SUB(0) check failed!" << endl;
            assert(0);exit(1);
        }
        assert(!bv_empty.any());
        for (unsigned j = 1; j < cnt; ++j)
        {
            agg.combine_and(bv_target1, agg_list, j);
            agg.combine_and_horizontal(bv_target2, agg_list, j);
            agg.combine_and_sub(bv_target3, agg_list, cnt, 0, 0, false);
            agg.combine_and_sub(bv_empty, agg_list, cnt, agg_list, cnt, false);

            
            res = bv_target1.compare(bv_target2);
            if (res!=0)
            {
                cerr << "Error: Aggregator AND check failed! 1.laddder step = "
                     << j << endl;
                assert(0);
                exit(1);
            }
            res = bv_target1.compare(bv_target3);
            if (res!=0)
            {
                cerr << "Error: Aggregator AND-SUB(0) check failed! 1.laddder step = "
                     << j << endl;
                assert(0);
                exit(1);
            }
            assert(!bv_empty.any());
        }
        
        for (unsigned j = 0; j < cnt; ++j)
        {
            if (j == 9)
                cerr << j << endl;
            agg.combine_and(bv_target1, agg_list+j, cnt-j);
            agg.combine_and_horizontal(bv_target2, agg_list+j, cnt-j);
            agg.combine_and_sub(bv_target3, agg_list+j, cnt-j, 0, 0, false);
            agg.combine_and_sub_horizontal(bv_target4, agg_list+j, cnt-j, 0, 0);
            agg.combine_and_sub(bv_empty, agg_list+j, cnt-j, agg_list+j, cnt-j, false);

            res = bv_target1.compare(bv_target2);
            if (res!=0)
            {
                cerr << "Error: Aggregator AND check failed! 2.laddder step = "
                     << j << endl;
                exit(1);
            }
            res = bv_target1.compare(bv_target4);
            if (res!=0)
            {
                cerr << "Error: Aggregator Horz-AND-SUB(0) check failed! 2.laddder step = "
                     << j << endl;
                res = bv_target3.compare(bv_target4);
                if (res == 0)
                {
                    cerr << "Warning. Aggregator AND-SUB ok... \n";
                }
                exit(1);
            }

            res = bv_target1.compare(bv_target3);
            if (res!=0)
            {
                cerr << "Error: Aggregator AND-SUB(0) check failed! 2.laddder step = "
                     << j << endl;
                exit(1);
            }
            assert(!bv_empty.any());
        }


    } // for i

  cout << "---------------------------- Aggregator AND Stress Test OK" << endl;
}


static
void GenerateTestCollection(std::vector<bvect>* target, unsigned count = 30, unsigned vector_max = 40000000)
{
    assert(target);
    bvect bv_common; // sub-vector common for all collection
    bvect_mini bvect_min(vector_max);
    
    FillSetsRandomMethod(&bvect_min, &bv_common, 0, vector_max, 1);
    
    for (unsigned i = 0; i < count; ++i)
    {
        std::unique_ptr<bvect> bv (new bvect);
        FillSetsRandomMethod(&bvect_min, bv.get(), 0, vector_max, 1);
        *bv |= bv_common;
        target->push_back(std::move(*bv));
    } // for
}


static
void StressTestAggregatorShiftAND(unsigned repeats)
{
   cout << "----------------------------StressTestAggregatorShiftAND " << endl;

    unsigned vector_max = 400000000;
    unsigned coll_size = 20;

    for (unsigned r = 0; r < repeats; ++r)
    {
        bvect mask_bv0;
        bvect_mini bvect_min(vector_max);
        FillSetsRandomMethod(&bvect_min, &mask_bv0, 0, vector_max - (vector_max / 5), 1);

        std::vector<bvect> bv_coll1;
        GenerateTestCollection(&bv_coll1, coll_size, vector_max);


        bm::aggregator<bvect> agg;
        agg.set_optimization();

        unsigned shift_repeats = 65536/3;
        for (unsigned i = 0; i < shift_repeats; ++i)
        {
            bvect bv_target0(mask_bv0);
            for (unsigned k = 0; k < bv_coll1.size(); ++k)
            {
                bv_target0.shift_right();
                bv_target0 &= bv_coll1[k];
            } // for
            
            agg.reset();
            agg.add(&mask_bv0);
            for (unsigned k = 0; k < bv_coll1.size(); ++k)
            {
                agg.add(&bv_coll1[k]);
            } // for
            
            bvect bv_target1;
            agg.combine_shift_right_and(bv_target1);
            auto cmp = bv_target1.compare(bv_target0);
            if (cmp != 0)
            {
                cerr << "Error: Mismatch! " << "STEP=" << i << endl;
                DetailedCheckVectors(bv_target0, bv_target1);
                exit(1);
            }
            if (i % 250 == 0)
                cout << "\r" << i << flush;
        } // for
        cout << "\n\n ---------- SHIFT-AND step: " << r << endl;
    } // for
    
   cout << "\n----------------------------StressTestAggregatorShiftAND OK" << endl;

}



static
void StressTest(unsigned repetitions, int set_operation, bool detailed)
{

   unsigned RatioSum = 0;
   unsigned SRatioSum = 0;
   unsigned DeltaSum = 0;
   unsigned SDeltaSum = 0;

   unsigned clear_count = 0;

   bvect  bvtotal;
   bvtotal.set_new_blocks_strat(bm::BM_GAP);

   bm::random_subset<bvect> rsub;


   cout << "----------------------------StressTest" << endl;

   unsigned size = BITVECT_SIZE - 10;
//size = BITVECT_SIZE / 10;
   unsigned i;
   for (i = 0; i < repetitions; ++i)
   {
        cout << endl << " - - - - - - - - - - - - STRESS STEP " << i; 
        switch (set_operation)
        {
        case 0: cout << " [OR]"; break;
        case 1: cout << " [SUB]";break;
        case 2: cout << " [XOR]";break;
        case 3: cout << " [AND]";break;
        default:
            cout << " [RANDOM]";
        }
        cout << endl;

        switch (rand() % 3)
        {
        case 0:
            size = BITVECT_SIZE / 10;
            break;
        case 1:
            size = BITVECT_SIZE / 2;
            break;
        default:
            size = BITVECT_SIZE - 10;
            break;
        } // switch


        bvect_mini*   bvect_min1= new bvect_mini(size);
        bvect_mini*   bvect_min2= new bvect_mini(size);
        bvect*        bvect_full1= new bvect();
        bvect*        bvect_full2= new bvect();

        bvect_full1->set_new_blocks_strat(i&1 ? bm::BM_GAP : bm::BM_BIT);
        bvect_full2->set_new_blocks_strat(i&1 ? bm::BM_GAP : bm::BM_BIT);

        int opt = rand() % 2;

        unsigned start1 = 0;

        switch (rand() % 3)
        {
        case 1:
            start1 += size / 5;
            break;
        default:
            break;
        }

        unsigned start2 = 0;

        switch (rand() % 3)
        {
        case 1:
            start2 += size / 5;
            break;
        default:
            break;
        }
       
        FillSetsRandomMethod(bvect_min1, bvect_full1, start1, size, opt);
        FillSetsRandomMethod(bvect_min2, bvect_full2, start2, size, opt);

        unsigned arr[bm::set_total_blocks]={0,};
        bm::id_t cnt = bvect_full1->count();
        unsigned last_block = bvect_full1->count_blocks(arr);
        unsigned sum = (unsigned)bm::sum_arr(&arr[0], &arr[last_block+1]);

        if (sum != cnt)
        {
            cout << "Error in function count_blocks." << endl;
            cout << "Array sum = " << sum << endl;
            cout << "BitCount = " << cnt << endl;
            cnt = bvect_full1->count();
            for ( i = 0; i <= last_block; ++i)
            {
                if (arr[i])
                {
                    cout << "[" << i << ":" << arr[i] << "]";
                }
            }
            cout << endl;
            cout << "================" << endl;
            print_stat(*bvect_full1);


            exit(1);
        }

        CheckCountRange(*bvect_full1, start1, BITVECT_SIZE);
        CheckIntervals(*bvect_full1, BITVECT_SIZE);


        CheckCountRange(*bvect_full2, start2, BITVECT_SIZE);

        CheckCountRange(*bvect_full1, 0, start1);
        CheckCountRange(*bvect_full2, 0, start2);


        TestRandomSubset(*bvect_full1, rsub);
        TestRandomSubset(*bvect_full2, rsub);

        // test find first difference
        //
        TestFindDiff(*bvect_full1, *bvect_full1);



/*        
        cout << "!!!!!!!!!!!!!!!" << endl;
        CheckVectors(*bvect_min1, *bvect_full1, size);
        cout << "!!!!!!!!!!!!!!!" << endl;
        CheckVectors(*bvect_min2, *bvect_full2, size);
        cout << "!!!!!!!!!!!!!!!" << endl;
 
        
         bvect_full1->stat();
         cout << " --" << endl;
         bvect_full2->stat();
*/

        int operation = rand()%5;
        if (set_operation != -1)
            operation = set_operation;

        switch(operation)
        {
        case 0:
            cout << "Operation OR" << endl;
            bvect_min1->combine_or(*bvect_min2);
            break;

        case 1:
            cout << "Operation SUB" << endl;
            bvect_min1->combine_sub(*bvect_min2);
            break;

        case 2:
            cout << "Operation XOR" << endl;
            bvect_min1->combine_xor(*bvect_min2);
            break;

        default:
            cout << "Operation AND" << endl;
            bvect_min1->combine_and(*bvect_min2);
            break;
        }

        int cres1 = bvect_min1->compare(*bvect_min2);

        delete bvect_min2;

        switch(operation)
        {
        case 0:
            {
            cout << "Operation OR" << endl;

            bm::id_t predicted_count = bm::count_or(*bvect_full1, *bvect_full2);
            bm::id_t predicted_any = bm::any_or(*bvect_full1, *bvect_full2);
            if (predicted_any == 0 && predicted_count != 0)
            {
                cout << "Predicted any error!" << endl;
                exit(1);
            }
            
            bvect    bv_target_s;
            SerializationOperation2Test(&bv_target_s,
                                        *bvect_full1,
                                        *bvect_full2,
                                        predicted_count,
                                        set_COUNT_OR,
                                        set_OR);

            bvect_full1->bit_or(*bvect_full2);
            
            bm::id_t count = bvect_full1->count();

            if (count != predicted_count)
            {
                cout << "Predicted count error!" << endl;
                cout << "Count = " << count << "Predicted count = " << predicted_count << endl;
                exit(1);
            }
            int res = bvect_full1->compare(bv_target_s);
            if (res != 0)
            {
                cout << "Serialization operation failed!" << endl;
                exit(1);
            }
            
            }
            break;

        case 1:
            {
            cout << "Operation SUB" << endl;
            
            bm::id_t predicted_count = bm::count_sub(*bvect_full1, *bvect_full2);
            bm::id_t predicted_any = bm::any_sub(*bvect_full1, *bvect_full2);
            if (predicted_any == 0 && predicted_count != 0)
            {
                cout << "Predicted any error!" << endl;
                exit(1);
            }
            
            bvect    bv_target_s;
            SerializationOperation2Test(&bv_target_s,
                                        *bvect_full1,
                                        *bvect_full2,
                                        predicted_count,
                                        set_COUNT_SUB_AB,
                                        set_SUB);

            bvect_full1->bit_sub(*bvect_full2);
            
            bm::id_t count = bvect_full1->count();

            if (count != predicted_count)
            {
                cout << "Predicted count error!" << endl;
                cout << "Count = " << count << "Predicted count = " << predicted_count << endl;
                exit(1);
            }
            int res = bvect_full1->compare(bv_target_s);
            if (res != 0)
            {
                cout << "Serialization operation failed!" << endl;
                exit(1);
            }
            
            
            }
            break;

        case 2:
            {
            cout << "Operation XOR <<<" << endl;
           
            bm::id_t predicted_count = bm::count_xor(*bvect_full1, *bvect_full2);
            bm::id_t predicted_any = bm::any_xor(*bvect_full1, *bvect_full2);
            if (predicted_any == 0 && predicted_count != 0)
            {
                cout << "Predicted any error!" << endl;
                exit(1);
            }

            bvect    bv_target_s;
            SerializationOperation2Test(&bv_target_s,
                                        *bvect_full1,
                                        *bvect_full2,
                                        predicted_count,
                                        set_COUNT_XOR,
                                        set_XOR);
            
            bvect_full1->bit_xor(*bvect_full2);
            
            bm::id_t count = bvect_full1->count();

            if (count != predicted_count)
            {
                cout << "Predicted count error!" << endl;
                cout << "Count = " << count << "Predicted count = " << predicted_count << endl;
                exit(1);
            }
            int res = bvect_full1->compare(bv_target_s);
            if (res != 0)
            {
                cout << "Serialization operation failed!" << endl;
                exit(1);
            }
            
            }
            
            break;

        default:
            {
            cout << "Operation AND" << endl;

            bm::id_t predicted_count = bm::count_and(*bvect_full1, *bvect_full2);
            bm::id_t predicted_any = bm::any_and(*bvect_full1, *bvect_full2);
            if (predicted_any == 0 && predicted_count != 0)
            {
                cout << "Predicted any error!" << endl;
                exit(1);
            }

            bvect    bv_target_s;
            SerializationOperation2Test(&bv_target_s,
                                        *bvect_full1,
                                        *bvect_full2,
                                        predicted_count,
                                        set_COUNT_AND,
                                        set_AND);

            TestRandomSubset(bv_target_s, rsub);


            bvect bv1(*bvect_full1);


            bvect_full1->bit_and(*bvect_full2);
            bm::id_t count = bvect_full1->count();

            int res = bvect_full1->compare(bv_target_s);
            if (res != 0)
            {
                //SaveBVector("bv1.bv", bv1);
                //SaveBVector("bv2.bv", *bvect_full2);
                cout << "Serialization operation failed!" << endl;
                exit(1);
            }

            if (count != predicted_count)
            {
                cout << "Predicted count error!" << endl;
                cout << "Count = " << count << "Predicted count = " << predicted_count << endl;
                exit(1);
            }

            }
            break;
        }



        cout << "Operation comparison" << endl;
        CheckVectors(*bvect_min1, *bvect_full1, size, detailed);

        int cres2 = bvect_full1->compare(*bvect_full2);

        CheckIntervals(*bvect_full1, BITVECT_SIZE);

        if (cres1 != cres2)
        {
            cout << cres1 << " " << cres2 << endl;
            cout << bvect_full1->get_first() << " " << bvect_full1->count() << endl;
            cout << bvect_full2->get_first() << " " << bvect_full2->count() << endl;

            cout << endl;
            printf("Bitset comparison operation failed.\n");
            assert(0); exit(1);
        }
        if (cres1 == 0) // re-confirm match
        {
            bvect::size_type pos;
            bool f = bvect_full1->find_first_mismatch(*bvect_full2, pos);
            if (f)
            {
                cerr << "Mismatch found pos=" << pos << endl;
                assert(0); exit(1);
            }
        }


        {
            bvect bv1(*bvect_full1);
            unsigned idx = unsigned(rand()) % size;
            bool b = bv1[idx];
            bool changed;
            if (b) 
            {
                changed = bv1.set_bit_conditional(idx, true, false);
                if (changed)
                {
                    cout << "Set bit conditional failed!" << endl;
                    exit(1);
                }
                b = bv1[idx];
                if (!b)
                {
                    cout << "Set bit conditional failed!" << endl;
                    exit(1);
                }

                changed = bv1.set_bit_conditional(idx, false, false);
                if (changed)
                {
                    cout << "Set bit conditional failed!" << endl;
                    exit(1);
                }
                changed = bv1.set_bit_conditional(idx, true, true);
                if (changed)
                {
                    cout << "Set bit conditional failed!" << endl;
                    exit(1);
                }
                changed = bv1.set_bit_conditional(idx, false, true);
                if (!changed)
                {
                    cout << "Set bit conditional failed!" << endl;
                    exit(1);
                }
                b = bv1[idx];
                if (b)
                {
                    cout << "Set bit conditional failed!" << endl;
                    exit(1);
                }
            } 
            else 
            {
                changed = bv1.set_bit_conditional(idx, false, true);
                if (changed)
                {
                    cout << "Set bit conditional failed!" << endl;
                    exit(1);
                }
                changed = bv1.set_bit_conditional(idx, true, false);
                if (!changed)
                {
                    cout << "Set bit conditional failed!" << endl;
                    exit(1);
                }
                b = bv1[idx];
                if (!b)
                {
                    cout << "Set bit conditional failed!" << endl;
                    exit(1);
                }
            }


        }

        delete bvect_full2;


        struct bvect::statistics st1;
        bvect_full1->calc_stat(&st1);
        bvect_full1->optimize();
        bvect_full1->optimize_gap_size();
        struct bvect::statistics st2;
        bvect_full1->calc_stat(&st2);

        unsigned Ratio = unsigned((st2.memory_used * 100)/st1.memory_used);
        RatioSum+=Ratio;
        DeltaSum+=unsigned(st1.memory_used - st2.memory_used);

        cout << "Optimization statistics: " << endl  
             << "   MemUsedBefore=" << st1.memory_used
             << "   MemUsed=" << st2.memory_used 
             << "   Ratio=" << Ratio << "%"
             << "   Delta=" << st1.memory_used - st2.memory_used
             << endl;
                
        cout << "Optimization comparison" << endl;

        CheckVectors(*bvect_min1, *bvect_full1, size, detailed);
        {
            bvect bv(*bvect_full1);
            bvect_full1->set_gap_levels(gap_len_table_min<true>::_len);
            {
                bvect::size_type pos;
                bool f = bv.find_first_mismatch(*bvect_full1, pos);
                if (f)
                {
                    cerr << "Mismatch found pos=" << pos << endl;
                    assert(0); exit(1);
                }
            }
        }
        CheckVectors(*bvect_min1, *bvect_full1, size, detailed);
        CheckIntervals(*bvect_full1, BITVECT_SIZE);

        //CheckCountRange(*bvect_full1, 0, size);


        // Serialization
       
        bm::serializer<bvect> bv_ser;
        bv_ser.gap_length_serialization(false);
        bv_ser.byte_order_serialization(false);
       
        bm::serializer<bvect>::buffer sermem_buf;
       
        bv_ser.serialize(*bvect_full1, sermem_buf, 0);
        unsigned slen = (unsigned)sermem_buf.size();
       
        delete bvect_full1;

        unsigned SRatio = unsigned((slen*100)/st2.memory_used);
        SRatioSum+=SRatio;
        SDeltaSum+=unsigned(st2.memory_used) - slen;


        cout << "Serialized mem_max = " << st2.max_serialize_mem 
             << " size= " << slen 
             << " Ratio=" << SRatio << "%"
             << " Delta=" << st2.memory_used - slen
             << endl;

        bvect*        bvect_full3= new bvect();
       
        bm::serializer<bvect>::buffer new_sermem_buf;
        new_sermem_buf = sermem_buf;
        cout << "Deserialization...";

        bm::deserialize(*bvect_full3, new_sermem_buf.buf());

        bm::deserialize(bvtotal, new_sermem_buf.buf());
        operation_deserializer<bvect> od;
        {
            int res;
            bvect bv_ac;
            bvect bv_a;
            bv_a.invert(); bv_ac.invert();
            bm::deserialize(bv_a, new_sermem_buf.buf());
            res = bv_a.compare(bv_ac);

            od.deserialize(bv_a,
                            new_sermem_buf.buf(),
                            0,
                            set_OR);
            res = bv_a.compare(bv_ac);
            assert(res == 0);
        }

        bvect* bv_target_s=new bvect();
        od.deserialize(*bv_target_s,
                        new_sermem_buf.buf(),
                        0,
                        set_OR);

        cout << "Ok." << endl;
//        delete [] new_sermem;

        cout << "Optimization...";
        bvtotal.optimize();
        cout << "Ok." << endl;

        ++clear_count;

        if (clear_count == 4)
        {
           bvtotal.clear();
           clear_count = 0;
        }

        cout << "Serialization comparison" << endl;

        {
            bvect::size_type pos;
            bool f = bv_target_s->find_first_mismatch(*bvect_full3, pos);
            if (f)
            {
                cerr << "Mismatch found pos=" << pos << endl;
                assert(0); exit(1);
            }
        }
        CheckVectors(*bvect_min1, *bvect_full3, size, detailed);

        int res = bv_target_s->compare(*bvect_full3);
        if (res != 0)
        {
            CheckVectors(*bvect_min1, *bv_target_s, size, true);
        }

        delete bv_target_s;
        delete bvect_min1;
        delete bvect_full3;

    }

    --i;
    cout << "Repetitions:" << i <<
            " AVG optimization ratio:" << RatioSum/i 
         << " AVG Delta:" << DeltaSum/i
         << endl
         << " AVG serialization Ratio:"<< SRatioSum/i
         << " Delta:" << SDeltaSum/i
         << endl;
}

static
void CheckGap2DGap(gap_vector& gapv)
{
   bm::gap_word_t   dgap_buf[bm::gap_max_buff_len+3]; 
   bm::gap_word_t   gap_buf[bm::gap_max_buff_len+3] = {0, };
   
   bm::gap_2_dgap(gapv.get_buf(), dgap_buf);
   bm::dgap_2_gap(dgap_buf, gap_buf);
   
   int c = bm::gapcmp(gap_buf, gapv.get_buf());
   if (c != 0)
   {
        cout << "Gap1: ";
        PrintGap(gapv.get_buf());
        cout << "D-Gap: ";
        PrintGap(dgap_buf);
        cout << "Gap2:";
        PrintGap(gap_buf);
        
        cout << "DGap conversion failed!" << endl;
        exit(1);
   }
} 

static
void GAPCheck()
{
   cout << "-------------------------------------------GAPCheck" << endl;

    {

    gap_vector   gapv(0);
    bvect_mini  bvect_min(bm::gap_max_bits);

    unsigned i;
    for( i = 0; i < 454; ++i)
    {
        bvect_min.set_bit(i);
        gapv.set_bit(i);
    }

    for(i = 0; i < 254; ++i)
    {
        bvect_min.clear_bit(i);
        gapv.clear_bit(i);
    }

    for(i = 5; i < 10; ++i)
    {
        bvect_min.set_bit(i);
        gapv.set_bit(i);
    }

    for( i = 0; i < bm::gap_max_bits; ++i)
    {
        int bit1 = (gapv.is_bit_true(i) == 1);
        int bit2 = (bvect_min.is_bit_true(i) != 0);
        int bit3 = (gapv.test(i) == 1);
        if (bit1 != bit2)
        {
            cout << "problem with bit comparison " << i << endl;
            exit(1);
        }
        if (bit1 != bit3)
        {
            cout << "problem with bit test comparison " << i << endl;
            exit(1);
        }

    }

    }


   {
   gap_vector gapv(1);
   int bit = gapv.is_bit_true(65535);

   if (bit != 1)
   {
      cout << "Bit != 1" << endl;
      exit(1);
   }
   
   unsigned i;
   for (i = 0; i < 65536; ++i)
   {
        bit = gapv.is_bit_true(i);
        if (bit != 1)
        {
            cout << "2.Bit != 1" << endl;
            exit(1);
        }
   }
   unsigned cnt = gapv.count_range(0, 65535);
   if (cnt != 65536)
   {
       cout << "count_range failed:" << cnt << endl;
       exit(1);
   }
   
   CheckCountGapRange(gapv, 10, 20);
   CheckCountGapRange(gapv, 0, 20);

   CheckGap2DGap(gapv);

   printf("gapv 1 check ok\n");
   }

   {
   gap_vector gapv(0);


   int bit = gapv.is_bit_true(65535);
   int bit1 = gapv.test(65535);
   if(bit != 0)
   {
      cout << "Bit != 0" << endl;
      exit(1);
   }
      
   unsigned i;
   for (i = 0; i < 65536; ++i)
   {
        bit = gapv.is_bit_true(i);
        bit1 = gapv.test(i);
        if (bit != 0)
        {
            cout << "2.Bit != 0 bit =" << i << endl;
            exit(1);
        }
        if (bit1 != 0)
        {
            cout << "2.Bit test != 0 bit =" << i << endl;
            exit(1);
        }
   }
   unsigned cnt = gapv.count_range(0, 65535);
   if (cnt != 0)
   {
       cout << "count_range failed:" << cnt << endl;
       exit(1);
   }
   CheckCountGapRange(gapv, 10, 20);
   CheckCountGapRange(gapv, 0, 20);

   CheckGap2DGap(gapv);



   printf("gapv 2 check ok\n");
   }

   {
   gap_vector gapv(0);

   gapv.set_bit(1);
   gapv.set_bit(0);

   gapv.control();
   CheckCountGapRange(gapv, 0, 20);

   int bit = gapv.is_bit_true(0);

   if (bit != 1)
   {
      cout << "Trouble" << endl;
      exit(1);
   }
   
   bit = gapv.is_bit_true(1);
   if (bit != 1)
   {
      cout << "Trouble 2." << endl;
      exit(1);
   }


   bit = gapv.is_bit_true(2);
   if(bit != 0)
   {
      cout << "Trouble 3." << endl;
      exit(1);
   }

   CheckGap2DGap(gapv);

   }

   {
   gap_vector gapv(0);

   gapv.set_bit(0);
   gapv.control();
   gapv.set_bit(1);
   gapv.control();

   gapv.set_bit(4);
   gapv.control();
   gapv.set_bit(5);
   gapv.control();
   CheckCountGapRange(gapv, 4, 5);
   CheckCountGapRange(gapv, 3, 5);

   gapv.set_bit(3);
   CheckCountGapRange(gapv, 3, 3);
   CheckCountGapRange(gapv, 3, 5);

   gapv.control();
   
   int bit = gapv.is_bit_true(0);
   if(bit!=1)
   {
      cout << "Bug" << endl;
   }
   bit = gapv.is_bit_true(1);
   if(bit!=1)
   {
      cout << "Bug2" << endl;
   }

   gapv.control();
   gapv.set_bit(4);
   gapv.control();

   CheckGap2DGap(gapv);


   printf("gapv 3 check ok\n");
   }

   {
        gap_vector gapv(0);
        bvect_mini   bvect_min(bm::gap_max_bits);
        
        cout << "++++++1" << endl;
        print_gap(gapv, 10);
        
        gapv.set_bit(bm::gap_max_bits-1);
        gapv.control();
        print_gap(gapv, 10);


        bvect_min.set_bit(bm::gap_max_bits-1);
        
        cout << "++++++3" << endl;
        
        gapv.set_bit(5);
        print_gap(gapv,15);
        gapv.control();
        bvect_min.set_bit(5);
        
        cout << "++++++4" << endl;

        CheckCountGapRange(gapv, 13, 150);
        gapv.control();
        
        CheckGap2DGap(gapv);


        unsigned i;
        for (i = 0; i < bm::gap_max_bits; ++i)
        {
            if (i == 65535)
                printf("%i\n", i);
            int bit1 = (gapv.is_bit_true(i) == 1);
            int bit2 = (bvect_min.is_bit_true(i) != 0);
            int bit3 = (gapv.test(i) == 1);
            if (bit1 != bit2)
            {
                cout << "problem with bit comparison " << i << endl;
            }
            if (bit1 != bit3)
            {
                cout << "problem with bit test comparison " << i << endl;
            }

        }

        gapv.clear_bit(5);
        bvect_min.clear_bit(5);
        gapv.control();
        
        CheckGap2DGap(gapv);


        for ( i = 0; i < bm::gap_max_bits; ++i)
        {
            if (i == 65535)
                printf("%i\n", i);
            int bit1 = (gapv.is_bit_true(i) == 1);
            int bit2 = (bvect_min.is_bit_true(i) != 0);
            int bit3 = (gapv.test(i) == 1);
            if (bit1 != bit2)
            {
                cout << "2.problem with bit comparison " << i << endl;
            }
            if (bit1 != bit3)
            {
                cout << "2.problem with bit test comparison " << i << endl;
            }
        }
   printf("gapv check 4 ok.\n");
   }

   {
        gap_vector gapv(0);
        bvect_mini   bvect_min(65536);
        
        unsigned i;
        for (i = 10; i > 0; i-=2)
        {
            bvect_min.set_bit(i);
            gapv.set_bit(i);
            gapv.control();
            CheckCountGapRange(gapv, 0, i);
            

            int bit1 = (gapv.is_bit_true(i) == 1);
            int bit2 = (bvect_min.is_bit_true(i) != 0);
            int bit3 = (gapv.test(i) != 0);
            if (bit1 != bit2)
            {
                cout << "3.problem with bit comparison " << i << endl;
            }
            if (bit1 != bit3)
            {
                cout << "3.problem with bit test comparison " << i << endl;
            }
            
            CheckGap2DGap(gapv);
            

        }
        for (i = 0; i < (int)bm::gap_max_bits; ++i)
        {
            int bit1 = (gapv.is_bit_true(i) == 1);
            int bit2 = (bvect_min.is_bit_true(i) != 0);
            int bit3 = (gapv.test(i) == 1);

            if (bit1 != bit2)
            {
                cout << "3.problem with bit comparison " << i << endl;
            }
            if (bit1 != bit3)
            {
                cout << "3.problem with bit test comparison " << i << endl;
            }
        }
   printf("gapv check 5 ok.\n");
   }

   {
        gap_vector gapv(0);
        bvect_mini   bvect_min(bm::gap_max_bits);
        
        int i;
        for (i = 0; i < 25; ++i)
        {
            unsigned id = random_minmax(0, bm::gap_max_bits);
            bvect_min.set_bit(id);
            gapv.set_bit(id);
            gapv.control();
            CheckCountGapRange(gapv, 0, id);
            CheckCountGapRange(gapv, id, 65535);
            
            CheckGap2DGap(gapv);

        }

        for (i = 0; i < (int)bm::gap_max_bits; ++i)
        {
            unsigned idx = unsigned(i);
            int bit1 = (gapv.is_bit_true(idx) == 1);
            int bit2 = (bvect_min.is_bit_true(idx) != 0);
            if (bit1 != bit2)
            {
                cout << "4.problem with bit comparison " << i << endl;
            }
        }

        for (i = bm::gap_max_bits; i < 0; --i)
        {
            unsigned idx = unsigned(i);
            int bit1 = (gapv.is_bit_true(idx) == 1);
            int bit2 = (bvect_min.is_bit_true(idx) != 0);
            if (bit1 != bit2)
            {
                cout << "5.problem with bit comparison " << i << endl;
            }
        }
   printf("gapv check 6 ok.\n");

   }

   printf("gapv random bit set check ok.\n");


   // conversion functions test
   
   {
   // aligned position test
   bvect        bvect_a;

   bvect_a.set_bit(1);
   bvect_a.set_bit(1, false);
   //bvect_a.clear();


   unsigned* buf = (unsigned*) bvect_a.get_blocks_manager().get_block_ptr(0, 0);

   bm::or_bit_block(buf, 0, 4);
   unsigned cnt = bm::bit_block_calc_count_range(buf, 0, 3);
   assert(cnt == 4);
   
   bool bit = bvect_a.get_bit(0);
   assert(bit);
   bit = bvect_a.get_bit(1);
   assert(bit);
   bit = bvect_a.get_bit(2);
   assert(bit);
   bit = bvect_a.get_bit(3);
   assert(bit);
   bit = bvect_a.get_bit(4);
   assert(bit==0);

   bm::or_bit_block(buf, 0, 36); 
   cnt = bm::bit_block_calc_count_range(buf, 0, 35);
   assert(cnt == 36);

   for (unsigned i = 0; i < 36; ++i)
   {
        bit = (bvect_a.get_bit(i) != 0);
        assert(bit);
   }
   bit = (bvect_a.get_bit(36) != 0);
   assert(bit==0);

   unsigned count = bvect_a.recalc_count();
   assert(count == 36);   
   
   cout << "Aligned position test ok." << endl; 

   }


   {
   // unaligned position test
   bvect   bvect_u;

   bvect_u.set_bit(0);
   bvect_u.set_bit(0, false);
//   bvect_u.clear();

   unsigned* buf = (unsigned*) bvect_u.get_blocks_manager().get_block_ptr(0, 0);

   bm::or_bit_block(buf, 5, 32);
   bool bit = (bvect_u.get_bit(4) != 0);
   assert(bit==0);
   unsigned cnt = bm::bit_block_calc_count_range(buf, 5, 5+32-1);
   assert(cnt == 32);
   cnt = bm::bit_block_calc_count_range(buf, 5, 5+32);
   assert(cnt == 32);

   unsigned i;
   for (i = 5; i < 4 + 32; ++i)
   {
        bit = bvect_u.get_bit(i);
        assert(bit);
   }
   unsigned count = bvect_u.recalc_count();
   assert(count==32);

   cout << "Unaligned position ok." << endl;

   } 

   // random test
   {
   cout << "random test" << endl;

   bvect   bvect_r;

   bvect_r.set_bit(0);
   bvect_r.set_bit(0, false);
   //bvect_r.clear();

   for (int i = 0; i < 5000; ++i)
   {
        unsigned* buf = (unsigned*) bvect_r.get_blocks_manager().get_block_ptr(0, 0);
        assert(buf);
        unsigned start = rand() % 65535;
        unsigned end = rand() % 65535;
        if (start > end)
        {
            unsigned tmp = end;
            end = start;
            start = tmp;
        }
        unsigned len = end - start;
        if (len)
        {
           bm::or_bit_block(buf, start, len);
           unsigned cnt = bm::bit_block_calc_count_range(buf, start, end);
           if (cnt != len)
           {
            cout << "random test: count_range comparison failed. " 
                 << " LEN = " << len << " cnt = " << cnt
                 << endl;
                 exit(1);
           }

           unsigned count = bvect_r.recalc_count();

           if (count != len)
           {
            cout << "random test: count comparison failed. " 
                 << " LEN = " << len << " count = " << count
                 << endl;
            assert(0);  exit(1);
           }            

           for (unsigned j = start; j < end; ++j)
           {
                bool bit = bvect_r.get_bit(j);
                if (!bit)
                {
                    cout << "random test: bit comparison failed. bit#" 
                         << i << endl;
                    exit(1);
                } 
           } // for j

        } 
        bvect_r.clear();
        bvect_r.set_bit(0);
        bvect_r.set_bit(0, false);

        if ((i % 100)==0)
        {
            cout << "*" << flush;
        }
   } // for i

   cout << endl << "Random test Ok." << endl;

   }


   // conversion test
 
   cout << "Conversion test" << endl;
    
   {
   
   gap_vector gapv(0);
   bvect   bvect_a;

   gapv.set_bit(0);
   gapv.set_bit(2);
   gapv.set_bit(10);
   gapv.set_bit(11);
   gapv.set_bit(12);
   
   CheckCountGapRange(gapv, 3, 15);

   print_gap(gapv, 100);
   bvect_a.set_bit(0);
   bvect_a.set_bit(0, false);
   //bvect_a.clear();

   unsigned* buf = (unsigned*) bvect_a.get_blocks_manager().get_block_ptr(0, 0);

   gapv.convert_to_bitset(buf);


   unsigned bitcount = bvect_a.recalc_count();


   if (bitcount != 5)
   {
      cout << "test failed: bitcout = " << bitcount << endl;
      exit(1);
   }


   gap_vector gapv1(0);
   gap_word_t* gap_buf = gapv1.get_buf();
   *gap_buf = 0;
   bit_convert_to_gap(gap_buf, buf, bm::gap_max_bits, bm::gap_max_buff_len);
   print_gap(gapv1, 100);

   bitcount = gapv1.bit_count();
   if(bitcount != 5)
   {
      cout << "2.test_failed: bitcout = " << bitcount << endl;
      exit(1);
   }

   printf("conversion test ok.\n");
    
   }

   // gap AND test

   {
   // special case 1: operand is all 1
   gap_vector gapv1(0);
   gapv1.set_bit(2);
   gap_vector gapv2(1); 

   gapv1.combine_and(gapv2.get_buf());
   gapv1.control();
   print_gap(gapv1, 0);

   unsigned count = gapv1.bit_count();
   assert(count == 1);
   int bit = gapv1.is_bit_true(2);
   if(bit == 0)
   {
      cout << "Wrong bit" << endl;
      exit(1);
   }
   CheckCountGapRange(gapv1, 0, 17);

   }

   {
   // special case 2: src is all 1
   gap_vector gapv1(1);
   gap_vector gapv2(0); 
   gapv2.set_bit(2);

   gapv1.combine_and(gapv2.get_buf());
   gapv1.control();
   print_gap(gapv1, 0);

   unsigned count = gapv1.bit_count();
   assert(count == 1);
   int bit = gapv1.is_bit_true(2);
   assert(bit);

   }

   {
   gap_vector gapv;
   gap_vector gapv1(0);

   gapv1.set_bit(3);
   gapv1.set_bit(4);
   print_gap(gapv1, 0);

   gap_vector gapv2(0); 
   gapv2.set_bit(2);
   gapv2.set_bit(3);
   print_gap(gapv2, 0);

   unsigned dsize=0;
   bm::gap_buff_op((gap_word_t*)gapv.get_buf(), 
                         gapv1.get_buf(), 0,
                         gapv2.get_buf(), 0, bm::and_op, 
                         dsize); 
   print_gap(gapv, 0);
   gapv.control();


    int bit1 = (gapv.is_bit_true(3) == 1);
    if(bit1 == 0)
    {
       cout << "Checking failed." << endl;
       exit(0);
    }

   gapv1.combine_or(gapv2);
   print_gap(gapv1, 0);
   gapv1.control();

   }

   {
        printf("gap AND test 1.\n");
        gap_vector gapv1(0);
        gap_vector gapv2(0);
        bvect_mini   bvect_min1(65536);
        bvect_mini   bvect_min2(65536);

        gapv1.set_bit(65535);
        bvect_min1.set_bit(65535);
        gapv1.set_bit(4);
        bvect_min1.set_bit(4);

        gapv2.set_bit(65535);
        bvect_min2.set_bit(65535);
        gapv2.set_bit(3);
        bvect_min2.set_bit(3);
        CheckCountGapRange(gapv2, 3, 65535);

        gapv2.control();

        printf("vect1:"); print_gap(gapv1, 0);
        printf("vect2:");print_gap(gapv2, 0);

        gapv1.combine_and(gapv2.get_buf());
        printf("vect1:");print_gap(gapv1, 0);

        gapv1.control();
        unsigned bit1 = (unsigned)gapv1.is_bit_true(65535u);
        assert(bit1);

        bvect_min1.combine_and(bvect_min2);
        CheckGAPMin(gapv1, bvect_min1, bm::gap_max_bits);
   }

   {
        printf("gap random AND test.\n");
        gap_vector gapv1(0);
        gap_vector gapv2(0);
        bvect_mini   bvect_min1(65536);
        bvect_mini   bvect_min2(65536);
        
        int i;
        for (i = 0; i < 25; ++i)
        {
            unsigned id = random_minmax(0, 65535);
            bvect_min1.set_bit(id);
            gapv1.set_bit(id);
            gapv1.control();
            CheckCountGapRange(gapv1, 0, id);
            CheckCountGapRange(gapv1, id, 65535);
        }
        for (i = 0; i < 25; ++i)
        {
            unsigned id = random_minmax(0, 65535);
            bvect_min2.set_bit(id);
            gapv2.set_bit(id);
            gapv2.control();
        }

        gapv1.combine_and(gapv2.get_buf());
        gapv1.control();
        gapv2.control();
        bvect_min1.combine_and(bvect_min2);

        CheckGAPMin(gapv1, bvect_min1, bm::gap_max_bits);

        printf("gap random AND test ok.\n");

   }

   {
        printf("gap OR test.\n");

        gap_vector gapv1(0);
        gap_vector gapv2(0);

        gapv1.set_bit(2);
        gapv2.set_bit(3);

        gapv1.combine_or(gapv2);
        gapv1.control();
        print_gap(gapv1, 0);   
        int bit1 = (gapv1.is_bit_true(0) == 1);
        assert(bit1==0);
        bit1=(gapv1.is_bit_true(2) == 1);
        assert(bit1);
        bit1=(gapv1.is_bit_true(3) == 1);
        assert(bit1);
   }

   {
        printf("gap XOR test.\n");

        gap_vector gapv1(0);
        gap_vector gapv2(0);

        gapv1.set_bit(2);
        gapv2.set_bit(3);
        gapv1.set_bit(4);
        gapv2.set_bit(4);
        print_gap(gapv1, 0);   
        print_gap(gapv2, 0);   

        gapv1.combine_xor(gapv2);
        gapv1.control();
        print_gap(gapv1, 0);   
        int bit1 = (gapv1.is_bit_true(0) == 0);
        assert(bit1);
        bit1=(gapv1.is_bit_true(2) == 1);
        assert(bit1);
        bit1=(gapv1.is_bit_true(3) == 1);
        assert(bit1);
        bit1=(gapv1.is_bit_true(4) == 0);
        assert(bit1);

   }


   {
        unsigned i;
        printf("gap random OR test.\n");
        gap_vector gapv1(0);
        gap_vector gapv2(0);
        bvect_mini   bvect_min1(bm::gap_max_bits);
        bvect_mini   bvect_min2(bm::gap_max_bits);
        
        for (i = 0; i < 10; ++i)
        {
            unsigned id = random_minmax(0, 100);
            bvect_min1.set_bit(id);
            gapv1.set_bit(id);
            gapv1.control();
        }
        for (i = 0; i < 10; ++i)
        {
            unsigned id = random_minmax(0, 100);
            bvect_min2.set_bit(id);
            gapv2.set_bit(id);
            gapv2.control();
        }

        print_mv(bvect_min1, 64);
        print_mv(bvect_min2, 64);

        gapv1.combine_or(gapv2);
        gapv1.control();
        gapv2.control();
        bvect_min1.combine_or(bvect_min2);

        print_mv(bvect_min1, 64);

        CheckGAPMin(gapv1, bvect_min1, bm::gap_max_bits);

        printf("gap random OR test ok.\n");

   }


   {
        unsigned i;
        printf("gap random SUB test.\n");
        gap_vector gapv1(0);
        gap_vector gapv2(0);
        bvect_mini   bvect_min1(bm::gap_max_bits);
        bvect_mini   bvect_min2(bm::gap_max_bits);
        
        for (i = 0; i < 25; ++i)
        {
            unsigned id = random_minmax(0, 100);
            bvect_min1.set_bit(id);
            gapv1.set_bit(id);
            gapv1.control();
        }
        for (i = 0; i < 25; ++i)
        {
            unsigned id = random_minmax(0, 100);
            bvect_min2.set_bit(id);
            gapv2.set_bit(id);
            gapv2.control();
        }

        print_mv(bvect_min1, 64);
        print_mv(bvect_min2, 64);

        gapv1.combine_sub(gapv2);
        gapv1.control();
        gapv2.control();
        bvect_min1.combine_sub(bvect_min2);

        print_mv(bvect_min1, 64);

        CheckGAPMin(gapv1, bvect_min1, bm::gap_max_bits);

        printf("gap random SUB test ok.\n");
   }

   {
       printf("GAP comparison test.\n");

       gap_vector gapv1(0);
       gap_vector gapv2(0);

       gapv1.set_bit(3);
       gapv2.set_bit(3);

       int res = gapv1.compare(gapv2);
       if (res != 0)
       {
           printf("GAP comparison failed!");
           exit(1);
       }

       gapv1.set_bit(4);
       gapv2.set_bit(4);

       res = gapv1.compare(gapv2);
       if (res != 0)
       {
           printf("GAP comparison failed!");
           exit(1);
       }

       gapv1.set_bit(0);
       gapv1.set_bit(1);

       res = gapv1.compare(gapv2);
       if (res != 1)
       {
           printf("GAP comparison failed!");
           exit(1);
       }

       gapv2.set_bit(0);
       gapv2.set_bit(1);
       res = gapv1.compare(gapv2);
       if (res != 0)
       {
           printf("GAP comparison failed!");
           exit(1);
       }

       gapv1.clear_bit(1);

       res = gapv1.compare(gapv2);
       if (res != -1)
       {
           printf("GAP comparison failed!");
           exit(1);
       }


   }
}

// -----------------------------------------------------------------------------
static
void SimpleGapFillSets(bvect&   bv0,
                       bvect&   bv1,
                       unsigned min,
                       unsigned max,
                       unsigned fill_factor)
{
    bvect::bulk_insert_iterator bii0(bv0);
    bvect::bulk_insert_iterator bii1(bv1);

    for (unsigned i = min; i < max; i += fill_factor)
    {
        bii0 = i;
        bii1 = i;
    } // for i
}

static
void GAPTestStress()
{
    cout << "----------------------------------- GAP test stress " << endl;
    const unsigned BV_SIZE = 65535 * 3;

    for (unsigned ff = 64; ff < 10000; ff++)
    {
        bvect bv0, bv1;
        SimpleGapFillSets(bv0, bv1, 0, BV_SIZE, ff);
        bv1.optimize();
        for (unsigned i = 0; i < BV_SIZE+1; ++i)
        {
            bool b0 = bv0.test(i);
            bool b1 = bv1.test(i);
            if (b0 != b1)
            {
                cerr << "GAP Test Stress failure!" << " FillFactor=" << ff << " bit=" << i << endl;
                exit(1);
            }
        } // for i
        if (ff % 100 == 0)
        {
            cout << "*" << flush;
        }
    } // for j

    cout << "----------------------------------- GAP test stress " << endl;
}

// -----------------------------------------------------------------------------
static
void MutationTest()
{

    cout << "--------------------------------- MutationTest" << endl;
    {
        bvect_mini     bvect_min(BITVECT_SIZE);
        bvect          bvect_full;

        printf("\nMutation test.\n");

        bvect_full.set_new_blocks_strat(bm::BM_GAP);

        bvect_full.set_bit(5);
        bvect_full.set_bit(5);

        bvect_min.set_bit(5);

        bvect_full.set_bit(65535);
        bvect_full.set_bit(65537);
        bvect_min.set_bit(65535);
        bvect_min.set_bit(65537);

        bvect_min.set_bit(100000);
        bvect_full.set_bit(100000);

        // detailed vectors verification
        ::CheckVectors(bvect_min, bvect_full, ITERATIONS, false);

        unsigned i;
        for (i = 5; i < 20000; i += 3)
        {
            bvect_min.set_bit(i);
            bvect_full.set_bit(i);
        }
        ::CheckVectors(bvect_min, bvect_full, ITERATIONS, false);

        for (i = 100000; i < 200000; i += 3)
        {
            bvect_min.set_bit(i);
            bvect_full.set_bit(i);
        }

        ::CheckVectors(bvect_min, bvect_full, 300000);
    }
    // set-clear functionality

    {
        printf("Set-clear functionality test.");

        bvect_mini     bvect_min(BITVECT_SIZE);
        bvect          bvect_full;
        bvect_full.set_new_blocks_strat(bm::BM_GAP);

        unsigned i;
        for (i = 100000; i < 100010; ++i)
        {
            bvect_min.set_bit(i);
            bvect_full.set_bit(i);            
        }
        ::CheckVectors(bvect_min, bvect_full, 300000);

        for (i = 100000; i < 100010; ++i)
        {
            bvect_min.clear_bit(i);
            bvect_full.clear_bit(i);            
        }
        ::CheckVectors(bvect_min, bvect_full, 300000);
        
        bvect_full.optimize();
        CheckVectors(bvect_min, bvect_full, 65536);//max+10);
    }

}

static
void MutationOperationsTest()
{

   cout << "------------------------------------ MutationOperationsTest" << endl;

   printf("Mutation operations test 1.\n");
   {
    bvect_mini   bvect_min1(BITVECT_SIZE);
    bvect_mini   bvect_min2(BITVECT_SIZE);
    bvect        bvect_full1;
    bvect        bvect_full2;

    bvect_full1.set_new_blocks_strat(bm::BM_GAP);
    bvect_full2.set_new_blocks_strat(bm::BM_BIT);

    bvect_full1.set_bit(100);
    bvect_min1.set_bit(100);

    unsigned i;
    for(i = 0; i < 10000; i+=2)
    {
       bvect_full2.set_bit(i);
       bvect_min2.set_bit(i);
    }
    print_stat(bvect_full2);
    CheckVectors(bvect_min2, bvect_full2, 65536, true);
    
    bvect_min1.combine_and(bvect_min2);
    bvect_full1.bit_and(bvect_full2);

    CheckVectors(bvect_min1, bvect_full1, 65536);//max+10);

   }

   printf("Mutation operations test 2.\n");
   {
    unsigned delta = 65536;
    bvect_mini   bvect_min1(BITVECT_SIZE);
    bvect_mini   bvect_min2(BITVECT_SIZE);
    bvect        bvect_full1;
    bvect        bvect_full2;

    bvect_full1.set_new_blocks_strat(bm::BM_GAP);
    bvect_full2.set_new_blocks_strat(bm::BM_GAP);

    unsigned i;
    for(i = 0; i < 1000; i+=1)
    {
       bvect_full1.set_bit(delta+i);
       bvect_min1.set_bit(delta+i);
    }

    for(i = 0; i < 100; i+=2)
    {
       bvect_full2.set_bit(delta+i);
       bvect_min2.set_bit(delta+i);
    }
//    CheckVectors(bvect_min2, bvect_full2, 65536);
    
    bvect_min1.combine_and(bvect_min2);
    bvect_full1.bit_and(bvect_full2);

    CheckVectors(bvect_min1, bvect_full1, 65536);//max+10);
    bvect_full1.optimize();
    CheckVectors(bvect_min1, bvect_full1, 65536);//max+10);

   }

   {
    bvect_mini   bvect_min1(BITVECT_SIZE);
    bvect        bvect_full1;

    bvect_full1.set_bit(3);
    bvect_min1.set_bit(3);

    struct bvect::statistics st;
    bvect_full1.calc_stat(&st);

    // serialization
    
    BM_DECLARE_TEMP_BLOCK(tb)

    unsigned char* sermem = new unsigned char[st.max_serialize_mem];
    size_t slen = bm::serialize(bvect_full1, sermem, tb);
    cout << "BVECTOR SERMEM=" << slen << endl;


    bvect        bvect_full3;
    bm::deserialize(bvect_full3, sermem);
    print_stat(bvect_full3);
    CheckVectors(bvect_min1, bvect_full3, 100, true);
   }


   printf("Mutation operations test 3.\n");
   {
    bvect_mini   bvect_min1(BITVECT_SIZE);
    bvect_mini   bvect_min2(BITVECT_SIZE);
    bvect        bvect_full1;
    bvect        bvect_full2;

    bvect_full1.set_new_blocks_strat(bm::BM_GAP);
    bvect_full2.set_new_blocks_strat(bm::BM_GAP);

   
    unsigned min = BITVECT_SIZE / 2 - ITERATIONS;
    unsigned max = BITVECT_SIZE / 2 + ITERATIONS;
    if (max > BITVECT_SIZE) 
        max = BITVECT_SIZE - 1;

    unsigned len = max - min;

    FillSets(&bvect_min1, &bvect_full1, min, max, 0);
    FillSets(&bvect_min1, &bvect_full1, 0, len, 5);
    printf("Bvect_FULL 1 STAT\n");
    print_stat(bvect_full1);
//    CheckVectors(bvect_min1, bvect_full1, max+10, false);
    FillSets(&bvect_min2, &bvect_full2, min, max, 0);
    FillSets(&bvect_min2, &bvect_full2, 0, len, 0);
    printf("Bvect_FULL 2 STAT\n");
    print_stat(bvect_full2);
//    CheckVectors(bvect_min2, bvect_full2, max+10);
    

    bvect_min1.combine_and(bvect_min2);
    bvect_full1.bit_and(bvect_full2);
    printf("Bvect_FULL 1 STAT after AND\n");
    print_stat(bvect_full1);

    CheckVectors(bvect_min1, bvect_full1, max+10, false);

    struct bvect::statistics st;
    bvect_full1.calc_stat(&st);
    cout << "BVECTOR: GAP=" << st.gap_blocks << " BIT=" << st.bit_blocks 
         << " MEM=" << st.memory_used << " SERMAX=" << st.max_serialize_mem
         << endl;
    cout << "MINIVECT: " << bvect_min1.mem_used() << endl;

    bvect_full1.optimize();
    print_stat(bvect_full1);

    CheckVectors(bvect_min1, bvect_full1, max+10, false);

    bvect_full1.calc_stat(&st);
    cout << "BVECTOR: GAP=" << st.gap_blocks << " BIT=" << st.bit_blocks 
         << " MEM=" << st.memory_used << " SERMAX=" << st.max_serialize_mem
         << endl;
    cout << "MINIVECT: " << bvect_min1.mem_used() << endl;



    // serialization
    
    BM_DECLARE_TEMP_BLOCK(tb)
    bm::serializer<bvect> bv_ser(tb);
    bm::serializer<bvect>::buffer sermem_buf;
    
    bv_ser.serialize(bvect_full1, sermem_buf, 0);
    unsigned slen = (unsigned)sermem_buf.size();

    cout << "BVECTOR SERMEM=" << slen << endl;


    
    bvect        bvect_full3;
    bm::deserialize(bvect_full3, sermem_buf.buf());
    print_stat(bvect_full3);
    CheckVectors(bvect_min1, bvect_full3, max+10, true);
    
    cout << "Copy constructor check." << endl;

    {
    bvect       bvect_full4(bvect_full3);
    print_stat(bvect_full3);
    CheckVectors(bvect_min1, bvect_full4, max+10, true);
    }
    

   }

}

static
void SerializationBufferTest()
{
   cout << " ----------------------------------- Serialization Buffer test" << endl;

    {
        bm::serializer<bvect>::buffer buf1;

        assert(buf1.size() == 0);
        assert(buf1.capacity() == 0);
        
        bm::serializer<bvect>::buffer buf2(100);
        assert(buf2.size() == 0);
        assert(buf2.capacity() != 0);
        
        buf1.reserve(100);
        assert(buf1.capacity() == buf2.capacity());
        
        const unsigned char str[] = "abc";
        buf1.copy_from(str, 3);
        assert(buf1.size() == 3);
        
        {
            const unsigned char* s = buf1.buf();
            assert(s[0] == 'a' && s[1] == 'b' && s[2] == 'c');
        }
        
        buf2 = buf1;
        assert(buf2.size() == buf1.size());

        {
            const unsigned char* s = buf2.buf();
            assert(s[0] == 'a' && s[1] == 'b' && s[2] == 'c');
        }
        buf2.reserve(100000);
        assert(buf2.size() == buf1.size());
        {
            const unsigned char* s = buf2.buf();
            assert(s[0] == 'a' && s[1] == 'b' && s[2] == 'c');
        }
        size_t cap2_1 = buf2.capacity();

        buf2.optimize();
        size_t cap2_2 = buf2.capacity();
        
        assert(cap2_2 < cap2_1);
        assert(buf2.size() == buf1.size());
        
        {
            const unsigned char* s = buf2.buf();
            assert(s[0] == 'a' && s[1] == 'b' && s[2] == 'c');
        }

        bm::serializer<bvect>::buffer buf3(buf2);
        assert(buf3.size() == buf2.size());
        {
            const unsigned char* s = buf3.buf();
            assert(s[0] == 'a' && s[1] == 'b' && s[2] == 'c');
        }
        
        buf3.reinit(1000000);
        assert(buf3.size() == 0);
    }
   

   cout << " ----------------------------------- Serialization Buffer test OK" << endl;
}


static
void SerializationCompressionLevelsTest()
{
   cout << " ----------------------------------- SerializationCompressionLevelsTest()" << endl;
   operation_deserializer<bvect> od;

   {
        BM_DECLARE_TEMP_BLOCK(tb)

        bvect bv(bm::BM_GAP);
        bv.set_range(0, 2);
        bv.set_range(10, 20);
        bv.set_range(100, 200);
        bv.optimize();
       
        bm::serializer<bvect> bv_ser(tb);
        bv_ser.set_compression_level(4); // use elias gamma
        bv_ser.set_bookmarks(true);
       
        bm::serializer<bvect>::buffer sermem_buf;

        bv_ser.serialize(bv, sermem_buf, 0);
       
        const bvect::size_type* cstat = bv_ser.get_compression_stat();
        assert(cstat[bm::set_block_gap_egamma] == 1);
       
        bvect bv2;
        bm::deserialize(bv2, sermem_buf.buf());

        int cmp = bv.compare(bv2);
        assert(cmp == 0);

        bvect bv3;
        od.deserialize(bv3,
                       sermem_buf.buf(),
                       tb,
                       set_OR);

        cmp = bv.compare(bv2);
        assert(cmp == 0);
   }

   {
        BM_DECLARE_TEMP_BLOCK(tb)

        bvect bv(bm::BM_GAP);
        bv.set_range(0, 2);
        bv.set_range(10, 20);
        bv.set_range(100, 200);
        bv.set_range(250, 300);
        bv.optimize();
       
        bm::serializer<bvect> bv_ser(tb);
        bv_ser.set_compression_level(5); // use interpolative encoder
        bv_ser.set_bookmarks(true);

        bm::serializer<bvect>::buffer sermem_buf;

        bv_ser.serialize(bv, sermem_buf, 0);
       
        const bvect::size_type* cstat = bv_ser.get_compression_stat();
        assert(cstat[bm::set_block_gap_bienc] == 1);
       
        bvect bv2;
        bm::deserialize(bv2, sermem_buf.buf());

        int cmp = bv.compare(bv2);
        assert(cmp == 0);
       
        bvect bv3;
        od.deserialize(bv3,
                       sermem_buf.buf(),
                       tb,
                       set_OR);

        cmp = bv.compare(bv2);
        assert(cmp == 0);
   }
   
   {
        BM_DECLARE_TEMP_BLOCK(tb)
        bvect bv { 0, 1, 2, 10, 100, 200 };
        bv.optimize();
       
        bm::serializer<bvect> bv_ser(tb);
        bv_ser.set_compression_level(4); // use elias gamma
        bv_ser.set_bookmarks(true);

        bm::serializer<bvect>::buffer sermem_buf;

        bv_ser.serialize(bv, sermem_buf, 0);
       
        const bvect::size_type* cstat = bv_ser.get_compression_stat();
        assert(cstat[bm::set_block_arrgap_egamma] == 1);
       
        bvect bv2;
        bm::deserialize(bv2, sermem_buf.buf());

        int cmp = bv.compare(bv2);
        assert(cmp == 0);
        bvect bv3;
        od.deserialize(bv3,
                       sermem_buf.buf(),
                       tb,
                       set_OR);
        cmp = bv.compare(bv2);
        assert(cmp == 0);
   }

   {
        BM_DECLARE_TEMP_BLOCK(tb)
        bvect bv { 0, 1, 2, 10, 100, 200 };
        bv.optimize();
       
        bm::serializer<bvect> bv_ser(tb);
        bv_ser.set_compression_level(5); // binary interpolated coding
        bv_ser.set_bookmarks(true);

        bm::serializer<bvect>::buffer sermem_buf;

        bv_ser.serialize(bv, sermem_buf, 0);
       
        const bvect::size_type* cstat = bv_ser.get_compression_stat();
        assert(cstat[bm::set_block_arrgap_bienc_v2] == 1);
       
        bvect bv2;
        bm::deserialize(bv2, sermem_buf.buf());

        int cmp = bv.compare(bv2);
        assert(cmp == 0);
        bvect bv3;
        od.deserialize(bv3,
                       sermem_buf.buf(),
                       tb,
                       set_OR);
        cmp = bv.compare(bv2);
        assert(cmp == 0);
   }
   
   {
        BM_DECLARE_TEMP_BLOCK(tb)
        bvect bv;
        bv.set_range(0, 65535);
        bv.clear_bit(1);
        bv.clear_bit(5);
        bv.clear_bit(10);
        bv.clear_bit(100);
        bv.clear_bit(200);
        bv.clear_bit(250);

        bv.optimize();
       
        bm::serializer<bvect> bv_ser(tb);
        bv_ser.set_compression_level(5); // binary interplated coding
        bv_ser.set_bookmarks(true);

        bm::serializer<bvect>::buffer sermem_buf;

        bv_ser.serialize(bv, sermem_buf, 0);
       
        const bvect::size_type* cstat = bv_ser.get_compression_stat();
        assert(cstat[bm::set_block_arrgap_bienc_inv_v2] == 1);
       
        bvect bv2;
        bm::deserialize(bv2, sermem_buf.buf());

        int cmp = bv.compare(bv2);
        assert(cmp == 0);
        bvect bv3;
        od.deserialize(bv3,
                       sermem_buf.buf(),
                       tb,
                       set_OR);
        cmp = bv.compare(bv2);
        assert(cmp == 0);
   }


   {
        BM_DECLARE_TEMP_BLOCK(tb)
       
        bvect bv(bm::BM_GAP);
        bv.set_range(0, bm::gap_max_bits-1);
        bv.clear_bit(1);
        bv.clear_bit(10);
        bv.clear_bit(100);
        bv.clear_bit(200);
       
        bv.optimize();
       
        bm::serializer<bvect> bv_ser(tb);
        bv_ser.set_compression_level(4); // use elias gamma
        bv_ser.set_bookmarks(true);

        bm::serializer<bvect>::buffer sermem_buf;

        bv_ser.serialize(bv, sermem_buf, 0);
       
        const bvect::size_type* cstat = bv_ser.get_compression_stat();
        assert(cstat[bm::set_block_arrgap_egamma_inv] == 1);
       
        bvect bv2;
        bm::deserialize(bv2, sermem_buf.buf());

        int cmp = bv.compare(bv2);
        assert(cmp == 0);
        bvect bv3;
        od.deserialize(bv3,
                       sermem_buf.buf(),
                       tb,
                       set_OR);
        cmp = bv.compare(bv2);
        assert(cmp == 0);
   }
   
   {
        BM_DECLARE_TEMP_BLOCK(tb)
       
        bvect bv(bm::BM_GAP);
        bv.set(100);
       
        bm::serializer<bvect> bv_ser(tb);
        bv_ser.set_compression_level(4); // use elias gamma
        bv_ser.set_bookmarks(true);

        bm::serializer<bvect>::buffer sermem_buf;

        bv_ser.serialize(bv, sermem_buf, 0);
       
        const bvect::size_type* cstat = bv_ser.get_compression_stat();
        assert(cstat[bm::set_block_bit_1bit] == 1);
       
        bvect bv2;
        bm::deserialize(bv2, sermem_buf.buf());

        int cmp = bv.compare(bv2);
        assert(cmp == 0);
        bvect bv3;
        od.deserialize(bv3,
                       sermem_buf.buf(),
                       tb,
                       set_OR);
        cmp = bv.compare(bv2);
        assert(cmp == 0);
   }

    // --------------------------------------------------------------
    // XOR serialization
    {{
        bvect bv1, bv2;
        bv1[1] = true;
        bv2[1] = true;

        bm::serializer<bvect>::bv_ref_vector_type bv_ref;
        bv_ref.add(&bv1, 1); // idx = 0
        bv_ref.add(&bv2, 5); // idx = 1

        bm::serializer<bvect> bms;
        bms.set_ref_vectors(&bv_ref);
        bms.set_curr_ref_idx(0);
        bms.set_bookmarks(true);

        bm::serializer<bvect>::buffer buf;
        bms.serialize(bv1, buf);

        const bvect::size_type* cstat = bms.get_compression_stat();
        assert(cstat[bm::set_block_ref_eq] == 1);

        bvect bv3;
        bm::deserialize(bv3, buf.buf(), 0, &bv_ref);
        auto eq = bv3.equal(bv2);
        assert(eq);

        bvect bv4;
        od.set_ref_vectors(&bv_ref);
        od.deserialize(bv4,
                       buf.buf(),
                       set_OR);
        eq = bv4.equal(bv2);
        assert(eq);

        bvect bv5;
        od.deserialize_range(bv5, buf.buf(), 0, 100);
        eq = bv5.equal(bv2);
        assert(eq);

    }}

    // --------------------------------------------------------------
    // XOR serialization (2)
    {{
        bvect bv1, bv2;
        for (unsigned i = 0; i < 100; i+=2)
        {
            bv1[i] = true;
            bv2[i+1] = true;
        }

        bm::serializer<bvect>::bv_ref_vector_type bv_ref;
        bv_ref.add(&bv1, 1); // idx = 0
        bv_ref.add(&bv2, 5); // idx = 1

        bm::serializer<bvect> bms;
        bms.set_ref_vectors(&bv_ref);
        bms.set_curr_ref_idx(0);
        bms.set_bookmarks(true);

        bm::serializer<bvect>::buffer buf;
        bms.serialize(bv1, buf);

        const bvect::size_type* cstat = bms.get_compression_stat();
        assert(cstat[bm::set_block_xor_ref32] == 1);

        bvect bv3;
        bm::deserialize(bv3, buf.buf(), 0, &bv_ref);
        auto eq = bv3.equal(bv1);
        assert(eq);

        bvect bv4;
        od.set_ref_vectors(&bv_ref);
        od.deserialize(bv4,
                       buf.buf(),
                       set_OR);
        eq = bv4.equal(bv1);
        assert(eq);

        bvect bv5;
        od.deserialize_range(bv5, buf.buf(), 0, 100);
        eq = bv5.equal(bv1);
        assert(eq);

        bvect bv6;
        bv6.set(0);
        assert(bv6.test(0) == bv1.test(0));

        bm::deserialize(bv6, buf.buf(), 0, &bv_ref);
        eq = bv6.equal(bv1);
        assert(eq);

        bvect bv7(bm::BM_GAP);
        bv7.set(0);
        bm::deserialize(bv7, buf.buf(), 0, &bv_ref);
        eq = bv7.equal(bv1);
        assert(eq);
        struct bvect::statistics st1;
        bv7.calc_stat(&st1);
        assert(!st1.bit_blocks);
        assert(st1.gap_blocks == 1);
    }}

    // --------------------------------------------------------------
    // XOR serialization (3) - GAPs

    {{
        bvect bv1, bv2;
        for (unsigned i = 0; i < 100; i+=2)
        {
            bv1[i] = true;
            bv2[i+1] = true;
        }
        bv1.optimize();
        bv2.optimize();
        {
            struct bvect::statistics st1;
            bv1.calc_stat(&st1);
            assert(!st1.bit_blocks);
            assert(st1.gap_blocks);
            struct bvect::statistics st2;
            bv2.calc_stat(&st2);
            assert(!st2.bit_blocks);
            assert(st2.gap_blocks);
        }

        bm::serializer<bvect>::bv_ref_vector_type bv_ref;
        bv_ref.add(&bv1, 1000000000); // idx = 0
        bv_ref.add(&bv2, 65500); // idx = 1

        bm::serializer<bvect> bms;
        bms.set_ref_vectors(&bv_ref);
        bms.set_curr_ref_idx(0);
        bms.set_bookmarks(true);

        bm::serializer<bvect>::buffer buf;
        bms.serialize(bv1, buf);

        const bvect::size_type* cstat = bms.get_compression_stat();
        assert(cstat[bm::set_block_xor_gap_ref32] == 1);

        bvect bv3;
        bm::deserialize(bv3, buf.buf(), 0, &bv_ref);
        auto eq = bv3.equal(bv1);
        assert(eq);

        bvect bv4;
        od.set_ref_vectors(&bv_ref);
        od.deserialize(bv4,
                       buf.buf(),
                       set_OR);
        eq = bv4.equal(bv1);
        assert(eq);

        bvect bv5;
        od.deserialize_range(bv5, buf.buf(), 0, 100);
        eq = bv5.equal(bv1);
        assert(eq);

        bvect bv6;
        bv6.set(0);
        assert(bv6.test(0) == bv1.test(0));

        bm::deserialize(bv6, buf.buf(), 0, &bv_ref);
        eq = bv6.equal(bv1);
        assert(eq);

        bvect bv7(bm::BM_GAP);
        bv7.set(0);
        bm::deserialize(bv7, buf.buf(), 0, &bv_ref);
        eq = bv7.equal(bv1);
        assert(eq);
        struct bvect::statistics st1;
        bv7.calc_stat(&st1);
        assert(!st1.bit_blocks);
        assert(st1.gap_blocks == 1);
    }}



   // --------------------------------------------------------------
   // bit-block tests
   //

   {
        bvect bv { 100 };
       
        bm::serializer<bvect> bv_ser;
        bv_ser.set_compression_level(4);
        bv_ser.set_bookmarks(true);

        bm::serializer<bvect>::buffer sermem_buf;

        bv_ser.serialize(bv, sermem_buf, 0);
       
        const bvect::size_type* cstat = bv_ser.get_compression_stat();
        assert(cstat[bm::set_block_bit_1bit] == 1);
       
        bvect bv2;
        bm::deserialize(bv2, sermem_buf.buf());
        int cmp = bv.compare(bv2);
        assert(cmp == 0);
        bvect bv3;
        od.deserialize(bv3,
                       sermem_buf.buf(),
                       0, set_OR);
        cmp = bv.compare(bv2);
        assert(cmp == 0);
   }

   {
        bvect bv;
        for (bvect::size_type i = 0; i < 65536; ++i)
            bv.set(i);
        bv.set(100, false);
       
        bm::serializer<bvect> bv_ser;
        bv_ser.set_compression_level(4);
        bv_ser.set_bookmarks(true);

        bm::serializer<bvect>::buffer sermem_buf;

        bv_ser.serialize(bv, sermem_buf, 0);
       
        const bvect::size_type* cstat = bv_ser.get_compression_stat();
        assert(cstat[set_block_arrbit_inv] == 1);
       
        bvect bv2;
        bm::deserialize(bv2, sermem_buf.buf());
        int cmp = bv.compare(bv2);
        assert(cmp == 0);
        bvect bv3;
        od.deserialize(bv3,
                       sermem_buf.buf(),
                       0, set_OR);
        cmp = bv.compare(bv2);
        assert(cmp == 0);
   }

   {
        bvect bv;
        for (bvect::size_type i = 0; i < 22345; i+=2)
            bv.set(1000+i);
       
        bm::serializer<bvect> bv_ser;
        bv_ser.set_compression_level(4);
        bv_ser.set_bookmarks(true);

        bm::serializer<bvect>::buffer sermem_buf;

        bv_ser.serialize(bv, sermem_buf, 0);
       
        const bvect::size_type* cstat = bv_ser.get_compression_stat();
        assert(cstat[bm::set_block_bit_0runs] == 1);
       
        bvect bv2;
        bm::deserialize(bv2, sermem_buf.buf());
        int cmp = bv.compare(bv2);
        assert(cmp == 0);
        bvect bv3;
        od.deserialize(bv3,
                       sermem_buf.buf(),
                       0, set_OR);
        cmp = bv.compare(bv2);
        assert(cmp == 0);
   }


   {
        bvect bv;
        for (bvect::size_type i = 0; i < 22345; i+=2)
            bv.set(1000+i);
       
        bm::serializer<bvect> bv_ser;
        bv_ser.set_compression_level(4);
        bv_ser.set_bookmarks(true);

        bm::serializer<bvect>::buffer sermem_buf;

        bv_ser.serialize(bv, sermem_buf, 0);
       
        const bvect::size_type* cstat = bv_ser.get_compression_stat();
        assert(cstat[bm::set_block_bit_0runs] == 1);
       
        bvect bv2;
        bm::deserialize(bv2, sermem_buf.buf());
        int cmp = bv.compare(bv2);
        assert(cmp == 0);
        bvect bv3;
        od.deserialize(bv3,
                       sermem_buf.buf(),
                       0, set_OR);
        cmp = bv.compare(bv2);
        assert(cmp == 0);
   }

   {
        bvect bv;
        for (bvect::size_type i = 0; i < 145; i+=64)
            bv.set(1000+i);
       
        bm::serializer<bvect> bv_ser;
        bv_ser.set_compression_level(3);
        bv_ser.set_bookmarks(true);

        bm::serializer<bvect>::buffer sermem_buf;

        bv_ser.serialize(bv, sermem_buf, 0);
       
        const bvect::size_type* cstat = bv_ser.get_compression_stat();
        assert(cstat[bm::set_block_arrbit] == 1);
       
        bvect bv2;
        bm::deserialize(bv2, sermem_buf.buf());
        int cmp = bv.compare(bv2);
        assert(cmp == 0);
        bvect bv3;
        od.deserialize(bv3,
                       sermem_buf.buf(),
                       0, set_OR);
        cmp = bv.compare(bv2);
        assert(cmp == 0);
   }


   {
        bvect bv;
        for (bvect::size_type i = 0; i < 545; i+=64)
            bv.set(1000+i);
       
        bm::serializer<bvect> bv_ser;
        bv_ser.set_compression_level(4);
        bv_ser.set_bookmarks(true);

        bm::serializer<bvect>::buffer sermem_buf;

        bv_ser.serialize(bv, sermem_buf, 0);
       
        const bvect::size_type* cstat = bv_ser.get_compression_stat();
        assert(cstat[bm::set_block_arrgap] == 1);
       
        bvect bv2;
        bm::deserialize(bv2, sermem_buf.buf());
        int cmp = bv.compare(bv2);
        assert(cmp == 0);
        bvect bv3;
        od.deserialize(bv3,
                       sermem_buf.buf(),
                       0, set_OR);
        cmp = bv.compare(bv2);
        assert(cmp == 0);
   }

   {
        bvect bv;
        for (bvect::size_type i = 0; i < 1045; i+=64)
        {
            auto target = i + 25;
            for ( ;i < target; ++i)
                bv.set(i);
        }
       
        bm::serializer<bvect> bv_ser;
        bv_ser.set_compression_level(4);
        bv_ser.set_bookmarks(true);

        bm::serializer<bvect>::buffer sermem_buf;

        bv_ser.serialize(bv, sermem_buf, 0);
       
        const bvect::size_type* cstat = bv_ser.get_compression_stat();
        assert(cstat[set_block_gap_egamma] == 1);
       
        bvect bv2;
        bm::deserialize(bv2, sermem_buf.buf());
        int cmp = bv.compare(bv2);
        assert(cmp == 0);
        bvect bv3;
        od.deserialize(bv3,
                       sermem_buf.buf(),
                       0, set_OR);
        cmp = bv.compare(bv2);
        assert(cmp == 0);
   }


   {
        bvect bv;
        for (bvect::size_type i = 0; i < 65536; ++i)
            bv.set(i);
        for (bvect::size_type i = 0; i < 1045; i+=64)
            bv.set(1000+i, false);

        bm::serializer<bvect> bv_ser;
        bv_ser.set_compression_level(4);
        bv_ser.set_bookmarks(true);

        bm::serializer<bvect>::buffer sermem_buf;

        bv_ser.serialize(bv, sermem_buf, 0);
       
        const bvect::size_type* cstat = bv_ser.get_compression_stat();
        assert(cstat[set_block_arrgap_egamma_inv] == 1);
       
        bvect bv2;
        bm::deserialize(bv2, sermem_buf.buf());
        int cmp = bv.compare(bv2);
        assert(cmp == 0);
        bvect bv3;
        od.deserialize(bv3,
                       sermem_buf.buf(),
                       0, set_OR);
        cmp = bv.compare(bv2);
        assert(cmp == 0);
   }

   {
        for (bvect::size_type k = 0; k < 4; ++k)
        {
            bvect bv;
            for (bvect::size_type i = 5; true; )
            {
                auto idx = (k * 1024) + i;
                if (idx >= 65536)
                    break;
                bv.set(idx);
                i += (rand() % 9);
                if (bv.count() > 13000)
                    break;
            }
            auto bc = bv.count();
           
            size_t l4size = 0;
            bvect bv_l4;
            {
                bm::serializer<bvect> bv_ser;
                bv_ser.set_compression_level(4);
                bm::serializer<bvect>::buffer sermem_buf;
                bv_ser.serialize(bv, sermem_buf, 0);
                //const bvect::size_type* cstat = bv_ser.get_compression_stat();
                //assert(cstat[bm::set_block_bit_0runs] == 1);
                l4size= sermem_buf.size();
            }

            bm::serializer<bvect> bv_ser;
            bv_ser.set_compression_level(5); // interpolated binary
            bv_ser.set_bookmarks(true);


            bm::serializer<bvect>::buffer sermem_buf;
            bv_ser.serialize(bv, sermem_buf, 0);
            size_t l5size = sermem_buf.size();
            size_t raw_int = bc * sizeof(bm::word_t);
            assert(raw_int > l5size);
            assert(l5size < l4size);
            cout << " offset = " << k * 1024 << endl;
            cout << "raw = " << raw_int << " l4 = " << l4size << " l5 = " << l5size << "  Diff(5-4)="
                 << l4size - l5size
                 << endl;
           
            const bvect::size_type* cstat = bv_ser.get_compression_stat();
            assert(cstat[bm::set_block_arr_bienc] == 1);
           
            bvect bv2;
            bm::deserialize(bv2, sermem_buf.buf());
            int cmp = bv.compare(bv2);
            assert(cmp == 0);
            bvect bv3;
            od.deserialize(bv3,
                           sermem_buf.buf(),
                           0, set_OR);
            cmp = bv.compare(bv2);
            assert(cmp == 0);
        }
   }


   {
        cout << "Generate large split in the mid-block" << endl;
        for (bvect::size_type k = 0; k < 4; ++k)
        {
            bvect bv;
            for (bvect::size_type i = 5; true; )
            {
                auto idx = i;
                if (idx >= 65536)
                    break;
                bv.set(idx);
                i += (rand() % 9);
                if (bv.count() > 3000)
                    break;
            }
            bvect bv_shift(bv);
            for (bvect::size_type i = 0; i < 10000 * 3; ++i)
            {
                bv_shift.shift_right();
            }
            bv |= bv_shift;


            auto bc = bv.count();
           
            size_t l4size = 0;
            bvect bv_l4;
            {
                bm::serializer<bvect> bv_ser;
                bv_ser.set_compression_level(4);
                bv_ser.set_bookmarks(true);

                bm::serializer<bvect>::buffer sermem_buf;
                bv_ser.serialize(bv, sermem_buf, 0);
                const bvect::size_type* cstat = bv_ser.get_compression_stat();
                assert(cstat[bm::set_block_bit_0runs] == 1 || cstat[bm::set_block_bit_digest0]);
                l4size= sermem_buf.size();
            }

            bm::serializer<bvect> bv_ser;
            bv_ser.set_compression_level(5); // interpolated binary
            bv_ser.set_bookmarks(true);


            bm::serializer<bvect>::buffer sermem_buf;
            bv_ser.serialize(bv, sermem_buf, 0);
            size_t l5size = sermem_buf.size();
            size_t raw_int = bc * sizeof(bm::word_t);
            assert(raw_int >= l5size);
            assert(l5size <= l4size);
            cout << " offset = " << k * 1024 << endl;
            cout << "raw = " << raw_int << " l4 = " << l4size << " l5 = " << l5size << "  Diff(5-4)="
                 << l4size - l5size
                 << endl;
           
            const bvect::size_type* cstat = bv_ser.get_compression_stat();
            assert(cstat[bm::set_block_arr_bienc]==1 ||
                   cstat[bm::set_block_bit_digest0]==1 ||
                   cstat[bm::set_block_bit_0runs]== 1);
           
            bvect bv2;
            bm::deserialize(bv2, sermem_buf.buf());
            int cmp = bv.compare(bv2);
            assert(cmp == 0);
            bvect bv3;
            od.deserialize(bv3,
                           sermem_buf.buf(),
                           0, set_OR);
            cmp = bv.compare(bv2);
            assert(cmp == 0);
        }
   }



   {
        bvect bv;
        bv.set(1); bv.set(1, false);
        bv.set_range(10, 20);
        bv.set_range(100, 200);
        bv.set_range(1000, 2000);
        bv.set_range(2010, 2020);
        bv.set_range(3000, 4020);
        bv.set_range(5000, 6000);
        bv.set_range(6000, 7000);

        bm::serializer<bvect> bv_ser;
        bv_ser.set_compression_level(5);
        bv_ser.set_bookmarks(true);

        bm::serializer<bvect>::buffer sermem_buf;

        bv_ser.serialize(bv, sermem_buf, 0);
       
        const bvect::size_type* cstat = bv_ser.get_compression_stat();
        assert(cstat[bm::set_block_gap_bienc] == 1);
       
        bvect bv2;
        bm::deserialize(bv2, sermem_buf.buf());
        int cmp = bv.compare(bv2);
        assert(cmp == 0);
        bvect bv3;
        od.deserialize(bv3,
                       sermem_buf.buf(),
                       0, set_OR);
        cmp = bv.compare(bv2);
        assert(cmp == 0);
   }

   {
        bvect bv;
        for (bvect::size_type i = 0; i < 65536; ++i)
            bv.set(i);
        for (bvect::size_type i = 0; i < 12045; ++i)
            bv.set(rand()%65535, false);

        bm::serializer<bvect> bv_ser;
        bv_ser.set_compression_level(5);
        bv_ser.set_bookmarks(true);

        bm::serializer<bvect>::buffer sermem_buf;

        bv_ser.serialize(bv, sermem_buf, 0);
       
        const bvect::size_type* cstat = bv_ser.get_compression_stat();
        assert(cstat[bm::set_block_arr_bienc_inv] == 1);
       
        bvect bv2;
        bm::deserialize(bv2, sermem_buf.buf());
        int cmp = bv.compare(bv2);
        assert(cmp == 0);
        bvect bv3;
        od.deserialize(bv3,
                       sermem_buf.buf(),
                       0, set_OR);
        cmp = bv.compare(bv2);
        assert(cmp == 0);
   }

   {
        bvect bv;
        bv.set(100);
        bvect::size_type from = 0;
        for (bvect::size_type i = 0; i < 3000; ++i)
        {
            auto to = from + rand()%15;
            bv.set_range(from, to);
            from = to + rand()%15;
        }

        bm::serializer<bvect> bv_ser;
        bv_ser.set_compression_level(5);
        bv_ser.set_bookmarks(true);

        bm::serializer<bvect>::buffer sermem_buf;
        bv_ser.serialize(bv, sermem_buf, 0);
 
        const bvect::size_type* cstat = bv_ser.get_compression_stat();
        assert(cstat[set_block_bitgap_bienc] == 1);
 
        bvect bv2;
        bm::deserialize(bv2, sermem_buf.buf());
        int cmp = bv.compare(bv2);
        assert(cmp == 0);
        bvect bv3;
        od.deserialize(bv3,
                       sermem_buf.buf(),
                       0, set_OR);
        cmp = bv.compare(bv2);
        assert(cmp == 0);
   }


   cout << " ----------------------------------- SerializationCompressionLevelsTest() OK" << endl;
}

static
void SerializationTest()
{

   cout << " ----------------------------------- SerializationTest" << endl;

   cout << "Compression level test (GAP blocks)" << endl;



   // ------------------------------------------------------------

   cout << "Serialization STEP 0" << endl;

   {
    unsigned size = BITVECT_SIZE/6000;


    bvect_mini*   bvect_min1= new bvect_mini(BITVECT_SIZE);
    bvect*        bvect_full1= new bvect();
    bvect*        bvect_full2= new bvect();
    bvect*        bv_target_s= new bvect();

    bvect_full1->set_new_blocks_strat(bm::BM_BIT);
    bvect_full2->set_new_blocks_strat(bm::BM_BIT);

    for(unsigned i = 0; i < size; ++i)
    {
        bvect_full1->set_bit(i);
        bvect_min1->set_bit(i);
    }

    bvect_full1->optimize();
    CheckVectors(*bvect_min1, *bvect_full1, size, true);



    bvect::statistics st;
    bvect_full1->calc_stat(&st);
    BM_DECLARE_TEMP_BLOCK(tb)
    
    bm::serializer<bvect> bv_ser(tb);
    
    bm::serializer<bvect>::buffer sermem_buf;
    bv_ser.optimize_serialize_destroy(*bvect_full1, sermem_buf);
    unsigned slen = (unsigned)sermem_buf.size();

    cout << "Serialized mem_max = " << st.max_serialize_mem
         << " size= " << slen 
         << " Ratio=" << (slen*100)/st.max_serialize_mem << "%"
         << endl;

    bm::deserialize(*bvect_full2, sermem_buf.buf());

    operation_deserializer<bvect> od;
    od.deserialize(*bv_target_s,
                   sermem_buf.buf(),
                   tb,
                   set_OR);


    CheckVectors(*bvect_min1, *bvect_full2, size, true);
    CheckVectors(*bvect_min1, *bv_target_s, size, true);


    delete bvect_full2;
    delete bvect_min1;
    delete bvect_full1;
    delete bv_target_s;

    }


   {
    unsigned size = BITVECT_SIZE/6000;


    bvect_mini*   bvect_min1= new bvect_mini(BITVECT_SIZE);
    bvect*        bvect_full1= new bvect();
    bvect*        bvect_full2= new bvect();
    bvect*        bv_target_s= new bvect();

    bvect_full1->set_new_blocks_strat(bm::BM_BIT);
    bvect_full2->set_new_blocks_strat(bm::BM_BIT);

        bvect_full1->set_bit(131072);
        bvect_min1->set_bit(131072);
    

    bvect_full1->optimize();

    bvect::statistics st;
    bvect_full1->calc_stat(&st);
    unsigned char* sermem = new unsigned char[st.max_serialize_mem];
    size_t slen = bm::serialize(*bvect_full1, sermem);
    cout << "Serialized mem_max = " << st.max_serialize_mem 
         << " size= " << slen 
         << " Ratio=" << (slen*100)/st.max_serialize_mem << "%"
         << endl;

    bm::deserialize(*bvect_full2, sermem);

    operation_deserializer<bvect> od;
    od.deserialize(*bv_target_s,
                   sermem,
                   0,
                   set_OR);

    delete [] sermem;

    CheckVectors(*bvect_min1, *bvect_full2, size, true);
    CheckVectors(*bvect_min1, *bv_target_s, size, true);

    delete bvect_full2;
    delete bvect_min1;
    delete bvect_full1;
    delete bv_target_s;

    }


    cout << "Serialization STEP 1." << endl;

    {
    bvect_mini   bvect_min1(BITVECT_SIZE);
    bvect        bvect_full1;

    bvect_full1.set_new_blocks_strat(bm::BM_GAP);
   
    unsigned min = BITVECT_SIZE / 2 - ITERATIONS;
    unsigned max = BITVECT_SIZE / 2 + ITERATIONS;
    if (max > BITVECT_SIZE) 
        max = BITVECT_SIZE - 1;

    unsigned len = max - min;

    FillSets(&bvect_min1, &bvect_full1, min, max, 0);
    FillSets(&bvect_min1, &bvect_full1, 0, len, 5);

    // shot some random bits

    unsigned i;
    for (i = 0; i < 10000; ++i)
    {
        unsigned bit = unsigned(rand()) % BITVECT_SIZE;
        bvect_full1.set_bit(bit);
        bvect_min1.set_bit(bit);
    }

    bvect::statistics st;
    bvect_full1.calc_stat(&st);

    unsigned char* sermem = new unsigned char[st.max_serialize_mem];
    print_stat(bvect_full1);
    
    size_t slen = bm::serialize(bvect_full1, sermem);

    cout << "Serialized len = " << slen << endl;

    bvect        bvect_full3;
    bm::deserialize(bvect_full3, sermem);
    bvect*  bv_target_s = new bvect();

    operation_deserializer<bvect> od;
    od.deserialize(*bv_target_s,
                   sermem,
                   0,
                   set_OR);

    CheckVectors(bvect_min1, bvect_full3, max+10, true);
    CheckVectors(bvect_min1, *bv_target_s, max+10, true);

    delete [] sermem;
    delete bv_target_s;

    }


   cout << "Stage 2" << endl;

   {

    bvect_mini*   bvect_min1= new bvect_mini(BITVECT_SIZE);
//    bm::bvect_mini*   bvect_min2= new bm::bvect_mini(BITVECT_SIZE);
    bvect*        bvect_full1= new bvect();
    bvect*        bvect_full2= new bvect();

    bvect_full1->set_new_blocks_strat(bm::BM_GAP);
    bvect_full2->set_new_blocks_strat(bm::BM_GAP);

    FillSetsRandomMethod(bvect_min1, bvect_full1, 1, BITVECT_SIZE-10, 1);
//    FillSetsRandomMethod(bvect_min2, bvect_full2, 1, BITVECT_SIZE-10, 1);

    bvect::statistics st;
    bvect_full1->calc_stat(&st);

    bm::serializer<bvect> bv_ser;
    bm::serializer<bvect>::buffer sermem_buf;
       
    bv_ser.serialize(*bvect_full1, sermem_buf, &st);
    unsigned slen = (unsigned)sermem_buf.size();

    cout << "Serialized mem_max = " << st.max_serialize_mem
         << " size= " << slen 
         << " Ratio=" << (slen*100)/st.max_serialize_mem << "%"
         << endl;
    bm::deserialize(*bvect_full2, sermem_buf.buf());
    bvect*  bv_target_s=new bvect();


    operation_deserializer<bvect> od;
    od.deserialize(*bv_target_s,
                   sermem_buf.buf(),
                   0,
                   set_OR);

    CheckVectors(*bvect_min1, *bvect_full2, BITVECT_SIZE, true);
    CheckVectors(*bvect_min1, *bv_target_s, BITVECT_SIZE, true);

    delete bv_target_s;
    delete bvect_full2;
    delete bvect_min1;
    delete bvect_full1;

    }



   cout << "Stage 3" << endl;

   {

    bvect_mini*   bvect_min1= new bvect_mini(BITVECT_SIZE);
    bvect_mini*   bvect_min2= new bvect_mini(BITVECT_SIZE);
    bvect*        bvect_full1= new bvect();
    bvect*        bvect_full2= new bvect();

    bvect_full1->set_new_blocks_strat(bm::BM_GAP);
    bvect_full2->set_new_blocks_strat(bm::BM_GAP);


    FillSetsRandomMethod(bvect_min1, bvect_full1, 1, BITVECT_SIZE, 1);
    FillSetsRandomMethod(bvect_min2, bvect_full2, 1, BITVECT_SIZE, 1);
    CheckVectors(*bvect_min1, *bvect_full1, BITVECT_SIZE, true);
    CheckVectors(*bvect_min2, *bvect_full2, BITVECT_SIZE, true);


    bvect::statistics st;
    bvect_full1->calc_stat(&st);
    unsigned char* sermem = new unsigned char[st.max_serialize_mem];
    size_t slen = bm::serialize(*bvect_full1, sermem);

    bvect bvt;
    bm::deserialize(bvt, sermem);
    if (bvt != *bvect_full1)
    {
        print_stat(bvt);
        print_stat(*bvect_full1);
        cout << "Error!" << endl;
        exit(1);
    }

    CheckVectors(*bvect_min1, *bvect_full1, BITVECT_SIZE, true);
    CheckVectors(*bvect_min2, *bvect_full2, BITVECT_SIZE, true);

    cout << "Serialized mem_max = " << st.max_serialize_mem 
         << " size= " << slen 
         << " Ratio=" << (slen*100)/st.max_serialize_mem << "%"
         << endl;

    bvect*  bv_target_s=new bvect(*bvect_full2);
    print_stat(*bv_target_s);

    print_stat(*bvect_full2);

    bvect*  bvect_full3= new bvect();
    *bvect_full3 = *bvect_full1;
    *bvect_full3 |= *bvect_full2;
//    CheckVectors(*bvect_min2, *bvect_full3, BITVECT_SIZE, true);


    bm::deserialize(*bvect_full2, sermem);

    operation_deserializer<bvect> od;
    od.deserialize(*bv_target_s,
                   sermem,
                   0,
                   set_OR);
    delete [] sermem;
    
    CheckVectors(*bvect_min1, *bvect_full1, BITVECT_SIZE, true);
//    CheckVectors(*bvect_min1, *bvect_full3, BITVECT_SIZE, true);

    bvect_min2->combine_or(*bvect_min1);
    delete bvect_min1;
    
    if (*bvect_full2 != *bvect_full3)
    {
        print_stat(*bvect_full2);
        print_stat(*bvect_full3);

        cout << "Error!" << endl;
        exit(1);
    }


    CheckVectors(*bvect_min2, *bvect_full2, BITVECT_SIZE, true);
    CheckVectors(*bvect_min2, *bv_target_s, BITVECT_SIZE, true);

    delete bv_target_s;
    delete bvect_full1;
    delete bvect_full2;
    delete bvect_full3;
    delete bvect_min2;    


    }

   cout << "Stage 4. " << endl;

   {
    unsigned size = BITVECT_SIZE/3;


    bvect_mini*   bvect_min1= new bvect_mini(BITVECT_SIZE);
    bvect*        bvect_full1= new bvect();
    bvect*        bvect_full2= new bvect();

    bvect_full1->set_new_blocks_strat(bm::BM_BIT);
    bvect_full2->set_new_blocks_strat(bm::BM_BIT);

    unsigned i;
    for(i = 0; i < 65000; ++i)
    {
        bvect_full1->set_bit(i);
        bvect_min1->set_bit(i);
    }

    for(i = 65536; i < 65536+65000; ++i)
    {
        bvect_full1->set_bit(i);
        bvect_min1->set_bit(i);
    }

    for (i = 65536*2; i < size/6; ++i)
    {
        bvect_full1->set_bit(i);
        bvect_min1->set_bit(i);
    }


    bvect_full1->optimize();

    print_stat(*bvect_full1);

    bvect::statistics st;
    bvect_full1->calc_stat(&st);
    unsigned char* sermem = new unsigned char[st.max_serialize_mem];
    size_t slen = bm::serialize(*bvect_full1, sermem);
    cout << "Serialized mem_max = " << st.max_serialize_mem 
         << " size= " << slen 
         << " Ratio=" << (slen*100)/st.max_serialize_mem << "%"
         << endl;
    
    unsigned char* new_sermem = new unsigned char[st.max_serialize_mem];
    ::memcpy(new_sermem, sermem, slen);

    bvect  bv_target_s(*bvect_full2);

    bm::deserialize(*bvect_full2, new_sermem);
    operation_deserializer<bvect> od;

    od.deserialize(bv_target_s,
                   new_sermem,
                   0,
                   set_OR);

    delete [] sermem;
    delete [] new_sermem;

    CheckVectors(*bvect_min1, *bvect_full2, size, true);
    CheckVectors(*bvect_min1, bv_target_s, size, true);


    delete bvect_full2;
    delete bvect_min1;
    delete bvect_full1;

    }


}

template<typename SV>
void generate_serialization_test_set(SV&   sv,
                                     typename SV::size_type vector_max)
{
    typename SV::back_insert_iterator bi(sv.get_back_inserter());

    unsigned v = 0;
    for (typename SV::size_type i = 0; i < vector_max; ++i)
    {
        unsigned plato = rand() % 32;
        for (unsigned j = 0; i < vector_max && j < plato; ++i, ++j)
        {
            *bi = v;
        } // for j
        if (++v > 100000)
            v = 0;
        unsigned nulls = rand() % 24;
        if (nulls)
            bi.add_null(nulls);
        i += nulls;
    } // for i
}

static
void TestSparseVectorSerialization2()
{
    cout << " ------------------------------ TestSparseVectorSerialization2()" << endl;
    const unsigned int BSIZE = 150000000;

    const unsigned char* buf;
    bool eq;
    size_t sz1, sz2;

    sparse_vector_u32 sv1(bm::use_null);
    sparse_vector_u32 sv2(bm::use_null);
    sparse_vector_u32 sv3(bm::use_null);

    generate_serialization_test_set(sv1, BSIZE);

    bm::sparse_vector_serial_layout<sparse_vector_u32> sv_lay;

    bm::sparse_vector_serializer<sparse_vector_u32> sv_serializer;
    bm::sparse_vector_deserializer<sparse_vector_u32> sv_deserial;

    for (unsigned k = 0; k < 2; ++k)
    {
        {
            {
                sv_serializer.set_xor_ref(false); // disable XOR compression
                sv_serializer.serialize(sv1, sv_lay);
            }

            buf = sv_lay.buf();
            sz1 = sv_lay.size();

            sv_deserial.deserialize(sv2, buf);

            eq = sv1.equal(sv2);
            if (!eq)
            {
                cerr << "Error: SparseVectorSerializationTest() integrity failure! (1)" << endl;
                sparse_vector_u32::size_type pos;

                bool f = bm::sparse_vector_find_first_mismatch(sv1, sv2, pos);
                assert(f);
                cerr << "Mismatch at: " << pos << endl;

                sv_deserial.deserialize(sv2, buf);

                exit(1);
            }
            sv2.resize(0);
        }

        {
            sv_serializer.set_xor_ref(true); // enable XOR compression
            sv_serializer.serialize(sv1, sv_lay);
        }

        buf = sv_lay.buf();
        sz2 = sv_lay.size();
        sv_deserial.deserialize(sv3, buf);
        eq = sv1.equal(sv3);
        if (!eq)
        {
            cerr << "Error: SparseVectorSerializationTest() integrity failure! (2)" << endl;
            sparse_vector_u32::size_type pos;
            bool f = bm::sparse_vector_find_first_mismatch(sv1, sv3, pos);
            assert(f);
            cerr << "Mismatch at: " << pos << endl;

            sv_deserial.deserialize(sv3, buf);

            exit(1);
        }

        if (sz2 > sz1)
        {
            cerr << "XOR negative compression!" << endl;
            assert(0);
        }
        else
        {
            cout << "sz1 = " << sz1 << " gain=" << (sz1 - sz2) << endl;
        }
        sv1.optimize();
    } // for k

    cout << " ------------------------------ TestSparseVectorSerialization2() OK" << endl;

}





static
void GetNextTest()
{
   cout << "-------------------------------------------- GetNextTest" << endl;
   
   cout << "testing bvector<>::find() in bit-mode" << endl;
   
   {
        bvect bv;
        bv.set();
        bvect::size_type pos=1;
        bool b;
        b = bv.find(pos);
        assert(pos == 0);
        assert(b);
       
        b = bv.find(bm::id_max/2, pos);
        assert(b);
        assert(pos == bm::id_max/2);
       
        bv.set(bm::id_max/2, false);
        b = bv.find(bm::id_max/2, pos);
        assert(b);
        assert(pos == bm::id_max/2 + 1);

        b = bv.find_reverse(pos);
        assert(b);
        assert(pos == bm::id_max-1);
       
        bvect::size_type f, l;
        b = bv.find_range(f, l);
        assert(b);
        assert(f == 0);
        assert(l == bm::id_max-1);
   }

   {
       bvect  bv;
       bool found;
       bm::id_t pos;
       found = bv.find(0, pos);
       
       if (found)
       {
           cout << "1. find() failed" << endl;
           exit(1);
       }
       found = bv.find_reverse(pos);
       assert(!found);
       
       bv[0] = true;
       found = bv.find(0, pos);
       if (!found || pos != 0)
       {
           cout << "2. find() failed " << pos << endl;
           exit(1);
       }
       found = bv.find_reverse(pos);
       assert(found && pos == 0);
       
       bv[0] = false;
       found = bv.find_reverse(pos);
       assert(!found);
       bv[0] = true;

       found = bv.find(1, pos);
       if (found)
       {
           cout << "3. find() failed" << endl;
           exit(1);
       }
       
       bv[100000] = true;
       bv[100001] = true;
       found = bv.find(1, pos);
       if (!found || pos != 100000)
       {
           cout << "4. find() failed " << pos << " " << found << endl;
           exit(1);
       }
       found = bv.find_reverse(pos);
       assert(found && pos == 100001);

       found = bv.find(100000, pos);
       if (!found || pos != 100000)
       {
           cout << "5. find() failed " << pos << " " << found << endl;
           exit(1);
       }
       found = bv.find(100001, pos);
       if (!found || pos != 100001)
       {
           cout << "6. find() failed " << pos << " " << found << endl;
           exit(1);
       }
       found = bv.find(100002, pos);
       if (found)
       {
           cout << "7. find() failed "<< endl;
           exit(1);
       }
       bv[100001] = false;
       found = bv.find_reverse(pos);
       assert(found && pos == 100000);
   }

   cout << "testing bvector<>::find() in GAP-mode" << endl;

   {
       bvect  bv(BM_GAP);
       bool found;
       bm::id_t pos;
       found = bv.find(0, pos);
       if (found)
       {
           cout << "1. find() failed" << endl;
           exit(1);
       }
       found = bv.find_reverse(pos);
       assert(!found);

       bv[0] = true;
       found = bv.find(0, pos);
       if (!found || pos != 0)
       {
           cout << "2. find() failed " << pos << endl;
           exit(1);
       }
       found = bv.find_reverse(pos);
       assert(found && pos == 0);

       found = bv.find(1, pos);
       if (found)
       {
           cout << "3. find() failed" << endl;
           exit(1);
       }
       bv[100000] = true;
       bv[100001] = true;
       found = bv.find(1, pos);
       if (!found || pos != 100000)
       {
           cout << "4. find() failed " << pos << " " << found << endl;
           exit(1);
       }
       found = bv.find(100000, pos);
       if (!found || pos != 100000)
       {
           cout << "5. find() failed " << pos << " " << found << endl;
           exit(1);
       }
       found = bv.find(100001, pos);
       if (!found || pos != 100001)
       {
           cout << "6. find() failed " << pos << " " << found << endl;
           exit(1);
       }
       found = bv.find_reverse(pos);
       assert(found && pos == 100001);

       found = bv.find(100002, pos);
       if (found)
       {
           cout << "7. find() failed "<< endl;
           exit(1);
       }
       bv[100001] = false;
       found = bv.find_reverse(pos);
       assert(found && pos == 100000);


   }

   {
       bvect  bv;
       bool found;
       
       bv.set_range(100000, 20000000);
       bm::id_t pos;
       found = bv.find_reverse(pos);
       assert(found && pos == 20000000);

       bv.optimize();
       found = bv.find_reverse(pos);
       assert(found && pos == 20000000);
       
       bv[bm::id_max-1] = true;
       found = bv.find_reverse(pos);
       assert(found && pos == bm::id_max-1);
       
       bv[bm::id_max-1] = false;
       found = bv.find_reverse(pos);
       assert(found && pos == 20000000);

       bv.set_range(100000, 20000000, false);
       found = bv.find_reverse(pos);
       assert(!found);

       found = bv.find(0, pos);
       assert(!found);
   }
   
   {
       bvect  bv;
       bool found;
       bm::id_t pos;
       bv.invert();
       
       found = bv.find_reverse(pos);
       assert(found && pos == bm::id_max-1);

       bv.optimize();
       found = bv.find_reverse(pos);
       assert(found && pos == bm::id_max-1);
       
       bv.invert();
       found = bv.find_reverse(pos);
       assert(!found);
       found = bv.find(0, pos);
       assert(!found);
   }


   int i;
   for(i = 0; i < 2; ++i)
   {
      cout << "Strategy " << i << endl;

   {
      bvect       bvect_full1;
      bvect_mini  bvect_min1(BITVECT_SIZE);

      bvect_full1.set_new_blocks_strat(i ? bm::BM_GAP : bm::BM_BIT);

      bvect_full1.set_bit(0);
      bvect_min1.set_bit(0);


      bvect_full1.set_bit(65536);
      bvect_min1.set_bit(65536);

      unsigned nbit1 = bvect_full1.get_first();
      unsigned nbit2 = bvect_min1.get_first();

      if (nbit1 != nbit2)
      {
         cout << "1. get_first failed() " <<  nbit1 << " " << nbit2 << endl;
         exit(1);
      }
      nbit1 = bvect_full1.get_next(nbit1);
      nbit2 = bvect_min1.get_next(nbit2);
      if ((nbit1 != nbit2) || (nbit1 != 65536))
      {
         cout << "1. get_next failed() " <<  nbit1 << " " << nbit2 << endl;
         exit(1);
      }
   }



   {
      bvect       bvect_full1;
      bvect_mini  bvect_min1(BITVECT_SIZE);
      bvect_full1.set_new_blocks_strat(i ? bm::BM_GAP : bm::BM_BIT);

      bvect_full1.set_bit(65535);
      bvect_min1.set_bit(65535);

      unsigned nbit1 = bvect_full1.get_first();
      unsigned nbit2 = bvect_min1.get_first();

      if ((nbit1 != nbit2) || (nbit1 != 65535))
      {
         cout << "1. get_first failed() " <<  nbit1 << " " << nbit2 << endl;
         exit(1);
      }
      nbit1 = bvect_full1.get_next(nbit1);
      nbit2 = bvect_min1.get_next(nbit2);
      if (nbit1 != nbit2 )
      {
         cout << "1. get_next failed() " <<  nbit1 << " " << nbit2 << endl;
         exit(1);
      }
   }

   {
      cout << "--------------" << endl;
      bvect       bvect_full1;
      bvect_mini  bvect_min1(BITVECT_SIZE);
      bvect_full1.set_new_blocks_strat(i ? bm::BM_GAP : bm::BM_BIT);

      bvect_full1.set_bit(655350);
      bvect_min1.set_bit(655350);

      unsigned nbit1 = bvect_full1.get_first();
      unsigned nbit2 = bvect_min1.get_first();

      if (nbit1 != nbit2 || nbit1 != 655350)
      {
         cout << "1. get_first failed() " <<  nbit1 << " " << nbit2 << endl;
         exit(1);
      }

      nbit1 = bvect_full1.get_next(nbit1);
      nbit2 = bvect_min1.get_next(nbit2);
      if (nbit1 != nbit2)
      {
         cout << "1. get_next failed() " <<  nbit1 << " " << nbit2 << endl;
         exit(1);
      }
   }


   {
   bvect       bvect_full1;
   bvect_mini  bvect_min1(BITVECT_SIZE);

   bvect_full1.set_new_blocks_strat(i ? bm::BM_GAP : bm::BM_BIT);

   bvect_full1.set_bit(256);
   bvect_min1.set_bit(256);

//   bvect_full1.clear_bit(256);
   bvect_full1.set_bit(65536);
   bvect_min1.set_bit(65536);

   unsigned nbit1 = bvect_full1.get_first();
   unsigned nbit2 = bvect_min1.get_first();

   if (nbit1 != nbit2)
   {
      cout << "get_first failed " <<  nbit1 << " " << nbit2 << endl;
      exit(1);
   }

   unsigned last_found = 0;
   while (nbit1)
   {
      cout << nbit1 << endl;
      nbit1 = bvect_full1.get_next(nbit1);
      nbit2 = bvect_min1.get_next(nbit2);
      if (nbit1 != nbit2)
      {
         cout << "get_next failed " <<  nbit1 << " " << nbit2 << endl;
         exit(1);
      }
     if (nbit1)
        last_found = nbit1;
   } // while
   
   unsigned pos = 0;
   bool found = bvect_full1.find_reverse(pos);
   assert(found && pos == last_found);

   }

   
   }// for

}

// Test contributed by Maxim Shemanarev.
static
void MaxSTest()
{
   bvect vec;

   int i, j;
   unsigned id;
   for(i = 0; i < 100; i++)
   {
      int n = rand() % 2000 + 1;
      id = 1;
      for(j = 0; j < n; j++)
      {
         id += rand() % 10 + 1;
         vec.set_bit(id);

      }
      vec.optimize();
      vec.clear();
      fprintf(stderr, ".");
   }
}

static
void CalcBeginMask()
{
    printf("BeginMask:\n");

    unsigned i;
    for (i = 0; i < 32; ++i)
    {
    unsigned mask = 0;

        for(unsigned j = i; j < 32; ++j)
        {
            unsigned nbit  = j; 
            nbit &= bm::set_word_mask;
            bm::word_t  mask1 = (((bm::word_t)1) << j);

            mask |= mask1;
        }

        printf("0x%x, ", mask);
        
    } 
    printf("\n");
}

static
void CalcEndMask()
{
    printf("EndMask:\n");

    unsigned i;
    for (i = 0; i < 32; ++i)
    {
    unsigned mask = 1;

        for(unsigned j = i; j > 0; --j)
        {
            unsigned nbit  = j; 
            nbit &= bm::set_word_mask;
            bm::word_t  mask1 = (((bm::word_t)1) << j);

            mask |= mask1;
        }

        printf("0x%x,", mask);
        
    } 
    printf("\n");
}

static
void EnumeratorTest()
{
    cout << "-------------------------------------------- EnumeratorTest" << endl;

    {
    bvect bvect1;

    bvect1.set_bit(100);
    
    {
        unsigned n = bvect1.get_next(101);
        assert(!n);
    }

    bvect::enumerator en = bvect1.first();
    unsigned n = bvect1.get_next(0);
    
    bvect::enumerator en1 = bvect1.get_enumerator(n);
    if (*en != 100 || n != 100 || *en1 != 100)
    {
        cout << "1.Enumerator error !" << endl;
        exit(1);
    }
    CompareEnumerators(en, en1);

    bvect1.clear_bit(100);

    bvect1.set_bit(2000000000);
    en.go_first();
    n = bvect1.get_next(0);
    en1.go_to(n);
    if (*en != 2000000000 || n != *en || *en1 != *en)
    {
        cout << "2. Enumerator error !" << endl;
        assert(0);
        exit(1);
    }
    CompareEnumerators(en, en1);

    bvect1.optimize();
    en = bvect1.first();
    n = bvect1.get_next(0);
    en1 = bvect1.first();
    en1.go_to(n);
    if (*en != 2000000000 || n != *en || *en1 != *en)
    {
        cout << "2. Enumerator error !" << endl;
        assert(0);
        exit(1);
    }
    CompareEnumerators(en, en1);

    }

    {
        bvect bvect1;
        bvect1.set_bit(0);
        bvect1.set_bit(10);
        bvect1.set_bit(35);
        bvect1.set_bit(1000);
        bvect1.set_bit(2016519);
        bvect1.set_bit(2034779);
        bvect1.set_bit(bm::id_max-1);

        bvect::enumerator en = bvect1.first();

        auto num = bvect1.get_first();

        bvect::enumerator end = bvect1.end();
        while (en < end)
        {
            cout << num << endl;
            bvect::enumerator en1 = bvect1.get_enumerator(num ? num-1 : num);
            if (*en != num || *en != *en1)
            {
                cout << "Enumeration comparison failed !" << 
                        " enumerator = " << *en <<
                        " get_next() = " << num <<
                        " goto enumerator = " << *en1 <<
                        endl;
                exit(1);
            }
            CompareEnumerators(en, en1);
            
            ++en;
            num = bvect1.get_next(num);
            ++en1;
            CompareEnumerators(en, en1);
            {
                auto num2 = num / 2;
                if (num2 < num)
                {
                    if (num2 == 2147483647)
                        cout << "!" << endl;
                    auto idx0 = bvect1.get_next(num2);
                    bvect::enumerator en3 = bvect1.get_enumerator(num2);
                    assert(idx0 == *en3);
                }
            }
        }
        if (num != 0)
        {
            cout << "Enumeration error!" << endl;
            exit(1);
        }
    }

    cout << "FULL bvector enumerator stress test ..." << endl;
    {
        bvect bvect1;
        bvect1.set();

        bvect::enumerator en = bvect1.first();
        unsigned num = bvect1.get_first();
        while (en.valid())
        {
            if (*en != num)
            {
                cout << "Enumeration comparison failed !" << 
                        " enumerator = " << *en <<
                        " get_next() = " << num << endl;
                assert(0);
                exit(1);
            }

            ++en;
            num = bvect1.get_next(num);
            {
                bvect::enumerator en2(&bvect1, num);
                if (*en2 != num)
                {
                    cout << "Enumeration comparison failed !" <<
                            " enumerator = " << *en <<
                            " get_next() = " << num << endl;
                    assert(0);
                    exit(1);
                }
                CompareEnumerators(en, en2);
            }
            if (num > (bm::set_sub_array_size * bm::gap_max_bits * 2))
                break;
        } // while
    }
    cout << "FULL bvector enumerator stress test ... OK" << endl;

    {
        bvect bvect1;

        unsigned i;
        for(i = 0; i < 65536; ++i)
        {
            bvect1.set_bit(i);
        }
        for(i = 65536*10; i < 65536*20; i+=3)
        {
            bvect1.set_bit(i);
        }


        bvect::enumerator en = bvect1.first();
        unsigned num = bvect1.get_first();

        while (en < bvect1.end())
        {
            bvect::enumerator en1 = bvect1.get_enumerator(num);
            if (*en != num || *en != *en1)
            {
                cout << "Enumeration comparison failed !" << 
                        " enumerator = " << *en <<
                        " get_next() = " << num <<
                        " goto enumerator = " << *en1
                        << endl; 
                exit(1);
            }
            ++en;
            num = bvect1.get_next(num);
            if (num == 31)
            {
                num = num + 0;
            }
            ++en1;
            CompareEnumerators(en, en1);
        }
        if (num != 0)
        {
            cout << "Enumeration error!" << endl;
            exit(1);
        }
    }


    {
    bvect bvect1;
    bvect1.set_new_blocks_strat(bm::BM_GAP);
    bvect1.set_bit(100);

    bvect::enumerator en = bvect1.first();
    bvect::enumerator en1 = bvect1.get_enumerator(99);
    if (*en != 100 || *en != *en1)
    {
        cout << "Enumerator error !" << endl;
        exit(1);
    }
    CompareEnumerators(en, en1);

    bvect1.clear_bit(100);

    bvect1.set_bit(2000000);
    en.go_first();
    en1.go_to(10);

    if (*en != 2000000 || *en != *en1)
    {
        cout << "Enumerator error !" << endl;
        exit(1);
    }
    CompareEnumerators(en, en1);
    print_stat(bvect1);
    }

    {
        bvect bvect1;
        bvect1.set_new_blocks_strat(bm::BM_GAP);
        bvect1.set_bit(0);
        bvect1.set_bit(1);
        bvect1.set_bit(10);
        bvect1.set_bit(100);
        bvect1.set_bit(1000);

        bvect::enumerator en = bvect1.first();

        unsigned num = bvect1.get_first();

        while (en < bvect1.end())
        {
            bvect::enumerator en1 = bvect1.get_enumerator(num);
            if (*en != num || *en != *en1)
            {
                cout << "Enumeration comparison failed !" << 
                        " enumerator = " << *en <<
                        " get_next() = " << num << 
                        " goto enumerator = " << *en1 << endl; 
                exit(1);
            }
            CompareEnumerators(en, en1);
            ++en;
            num = bvect1.get_next(num);
            ++en1;
            CompareEnumerators(en, en1);
        }
        if (num != 0)
        {
            cout << "Enumeration error!" << endl;
            exit(1);
        }
    }
}


static
void BlockLevelTest()
{
    bvect  bv;
    bvect  bv2;

    bv.set_new_blocks_strat(bm::BM_GAP);
    bv2.set_new_blocks_strat(bm::BM_GAP);

    unsigned i;
    for (i = 0; i < 500; i+=1)
    {
        bv.set_bit(i);
    }
    print_stat(bv);

    for (i = 0; i < 1000; i+=2)
    {
        bv2.set_bit(i);
    }
    print_stat(bv2);

    struct bvect::statistics st;
    bv2.calc_stat(&st);

    unsigned char* sermem = new unsigned char[st.max_serialize_mem];

    size_t slen = bm::serialize(bv2, sermem);
    assert(slen);
    slen = 0;
    
    bm::deserialize(bv, sermem);

    print_stat(bv);

}

/*
__int64 CalcBitCount64(__int64 b)
{
    b = (b & 0x5555555555555555) + (b >> 1 & 0x5555555555555555);
    b = (b & 0x3333333333333333) + (b >> 2 & 0x3333333333333333);
    b = b + (b >> 4) & 0x0F0F0F0F0F0F0F0F;
    b = b + (b >> 8);
    b = b + (b >> 16);
    b = b + (b >> 32) & 0x0000007F;
    return b;
}


*/

static
void FindNotNullPtrTest()
{
    cout << "----------------------------- FindNotNullPtrTest()" << endl;
    bm::word_t*** arr = 0;
    unsigned arr_size = 895;
    arr = (bm::word_t***)::malloc(sizeof(void*) * arr_size);

    for (unsigned i = 0; i < arr_size; ++i)
    {
        arr[i] = 0;
    }
    bool found;
    unsigned pos;
    found = bm::find_not_null_ptr(arr, 0u, arr_size, &pos);
    assert(!found);
    
    for (unsigned i = arr_size-1; i > 0; --i)
    {
        arr[i] = (bm::word_t**)~0;
        for (unsigned j = 0; j < i; ++j)
        {
            found = bm::find_not_null_ptr(arr, j, arr_size, &pos);
            assert(found);
            assert(pos == i);
        }
    } // for i

    ::free(arr);
    cout << "----------------------------- FindNotNullPtrTest() OK" << endl;
}

// function to return bvector by value to test move semantics
//
static
bvect bvect_test_return()
{
    bvect bv1;
    bvect bv2;

    bv1[100] = true;
    bv1[1000] = true;
    bv2[100] = true;
    bv2[10001] = true;

    if (rand()%2)
    {
        return bv1;
    }
    return (bv1 & bv2);
}



static
void SyntaxTest()
{
    cout << "----------------------------- Syntax test." << endl;
    
    {
        bvect bv1;
        bv1.set_bit(100);
        bv1.set_bit(100 + 10 *65535 * 256);
        {
        bvect bv2(bv1);
        int res = bv2.compare(bv1);
        assert(res == 0);
        }
    }
    {
        bvect bv1;
        
        //bvect::allocator_type a = bv1.get_allocator();

        bvect bv2(bv1);
        bvect bv3;
        bv3.swap(bv1);
        
        bv1[100] = true;
        bool v = bv1[100];
        assert(v);
        v = false;

        bv1[100] = false;

        bv2 |= bv1;
        bv2 &= bv1;
        bv2 ^= bv1;
        bv2 -= bv1;

        bv3 = bv1 | bv2;

        if (bv1 < bv2)
        {
        }
        
        bv3 = bv1 & bv2;
        bv3 = bv1 ^ bv2;

        bvect::reference ref = bv1[10];
        bool bn = !ref;
        bool bn2 = ~ref;
        bv1[10] = bn2;
        bv1[10] = bn;

        bn = bn2 = false;

        ref.flip();

        bvect bvn = ~bv1;
        
        // this should trigger move
        bvect bv4 = bvect_test_return();
        bvect bv41 = bvect_test_return() | bv2;
        bvect bv5(bvect_test_return());
        
        cout << bv4.count() << " " << bv41.count() << " " << bv5.count() << endl;
    }
    
    
    {
        bvect bv0;
        bvect bv1, bv2, bv3;
        bv0 = bv1 | bv2 | bv3;
        assert(bv0.count() == 0);
    }
    {
        bvect bv0;
        bvect bv1, bv2, bv3;
        bv0 = bv1 & bv2 & bv3;
        assert(bv0.count() == 0);
    }
    {
        bvect bv0;
        bvect bv1, bv2, bv3;
        bv0 = bv1 ^ bv2 ^ bv3;
        assert(bv0.count() == 0);
    }
    {
        bvect bv0;
        bvect bv1, bv2, bv3;
        bv0 = bv1 - bv2 - bv3;
        assert(bv0.count() == 0);
    }


    cout << "----------------------------- Syntax test ok." << endl;
}

static
void SetTest()
{
    {
        bvect bv{ 0, 10, 65536, 10000, bm::id_max-1 };
        unsigned cnt = bv.count();
        if (cnt != 5)
        {
            cout << "Brace initialization test failed!." << endl;
            exit(1);
        }
        bvect bv2;
        bv2.set(0).set(10).set(65536).set(10000).set(bm::id_max-1);

        if (bv != bv2)
        {
            cout << "Brace initialization comparison test failed!." << endl;
            exit(1);
        }
    }
    {
        unsigned cnt;
        bvect bv;
        bv.set();

        cnt = bv.count();
        if (cnt != bm::id_max)
        {
            cout << "Set test failed!." << endl;
            assert(0);
            exit(1);
        }

        bv.invert();
        cnt = bv.count();
        if (cnt != 0)
        {
            cout << "Set invert test failed!." << endl;
            exit(1);
        }

        bv.set(0);
        bv.set(bm::id_max - 1);
        cnt = bv.count();

        assert(cnt == 2);

        bv.invert();
        print_stat(bv);
        cnt = bv.count();

        if (cnt != bm::id_max - 2)
        {
            cout << "Set invert test failed!." << endl;
            exit(1);
        }

        bv.clear();
        bv[1] &= true;
        bool v = bv[1];
        if (v)
        {
            cout << "Set &= test failed!" << endl;
            exit(1);
        }


        bv[1] = true;
        bv[1] &= true;
        v = bv[1];
        if (!v)
        {
            cout << "Set &= test failed (2)!" << endl;
            exit(1);
        }
        bv.clear(true);
        bv.invert();
        bv[1] &= true;
        v = bv[1];
        if (!v)
        {
            cout << "Set &= test failed (2)!" << endl;
            exit(1);
        }
    }
    
    {
        bvect bv_full;
        bv_full.invert();
        assert(bv_full.test(bm::id_max/2));
    }
    
    {
        bvect bv1;
        bv1.set(0);
        bv1.set();
        auto cnt1 = bv1.count();
        assert (cnt1 == bm::id_max);
    }
    
    {
        bvect bv1, bv2(BM_GAP);
        bv1.set(0); bv2.set(0);
        bv1.set(bm::id_max-1);bv2.set(bm::id_max-1);
        bv1.set((bm::id_max-1)/2);bv2.set((bm::id_max-1)/2);
        for (unsigned i = 0; i < 2; ++i)
        {
            bv1.set();
            bv2.set();
            auto cnt1 = bv1.count();
            auto cnt2 = bv2.count();
            assert (cnt1 == bm::id_max);
            assert (cnt2 == bm::id_max);
        }
    }
    

    bvect bv2;
    bv2[1] = true;
    bv2[1] = false;
    bvect::statistics stat1;
    bv2.calc_stat(&stat1);
    
    bv2.optimize();

    bvect::statistics stat2;
    bv2.calc_stat(&stat2);

    if (stat2.bit_blocks != 0 || 
        stat2.gap_blocks != 0 ||
        stat1.memory_used <= stat2.memory_used)
    {
        cout << "Optimization memory free test failed (2)!" << endl;
        exit(1);
    }


    {
        bvect bv3;
        bool changed;
        changed = bv3.set_bit_conditional(10, true, true);
        bool v = bv3[10];
        if (v || changed) {
            cout << "Conditional bit set failed." << endl;
            exit(1);
        }
        changed = bv3.set_bit_conditional(10, true, false);
        v = bv3[10];
        if (!v || !changed) {
            cout << "Conditional bit set failed." << endl;
            exit(1);
        }
        changed = bv3.set_bit_conditional(10, false, false);
        v = bv3[10];
        if (!v || changed) {
            cout << "Conditional bit set failed." << endl;
            exit(1);
        }
        changed = bv3.set_bit_conditional(10, false, true);
        v = bv3[10];
        if (v || !changed) {
            cout << "Conditional bit set failed." << endl;
            exit(1);
        }
    }

    {
        bvect bv3(bm::BM_GAP);
        bool changed;
        changed = bv3.set_bit_conditional(10, true, true);
        bool v = bv3[10];
        if (v || changed) {
            cout << "Conditional bit set failed." << endl;
            exit(1);
        }
        changed = bv3.set_bit_conditional(10, true, false);
        v = bv3[10];
        if (!v || !changed) {
            cout << "Conditional bit set failed." << endl;
            exit(1);
        }
        changed = bv3.set_bit_conditional(10, false, false);
        v = bv3[10];
        if (!v || changed) {
            cout << "Conditional bit set failed." << endl;
            exit(1);
        }
        changed = bv3.set_bit_conditional(10, false, true);
        v = bv3[10];
        if (v || !changed) {
            cout << "Conditional bit set failed." << endl;
            exit(1);
        }
    }

    {
        bvect bv3(bm::BM_GAP);
        bv3.invert();
        bv3.optimize();
        bool changed;
        changed = bv3.set_bit_conditional(10, true, true);
        bool v = bv3[10];
        if (!v || changed) {
            cout << "Conditional bit set failed." << endl;
            exit(1);
        }
        changed = bv3.set_bit_conditional(10, true, false);
        v = bv3[10];
        if (!v || changed) {
            cout << "Conditional bit set failed." << endl;
            exit(1);
        }
        changed = bv3.set_bit_conditional(10, false, false);
        v = bv3[10];
        if (!v || changed) {
            cout << "Conditional bit set failed." << endl;
            exit(1);
        }
        changed = bv3.set_bit_conditional(10, false, true);
        v = bv3[10];
        if (v || !changed) {
            cout << "Conditional bit set failed." << endl;
            exit(1);
        }
        changed = bv3.set_bit_conditional(10, true, false);
        v = bv3[10];
        if (!v || !changed) {
            cout << "Conditional bit set failed." << endl;
            exit(1);
        }
    }


    {
        bvect bv(0);
        bv.resize(100);
        bv[10] = true;
        bv.resize(1000001);
        bv[100000] = 1;

        if (bv.size() != 1000001)
        {
            cout << "Resize failed" << endl;
            exit(1);
        }
        if (bv.count() != 2)
        {
            cout << "Resize count failed" << endl;
            exit(1);
        }

        bv.resize(100);
        if (bv.size() != 100)
        {
            cout << "Resize failed" << endl;
            exit(1);
        }
        if (bv.count() != 1)
        {
            cout << "Resize count failed" << endl;
            exit(1);
        }
        
        bv.resize(60000100);
        bv.invert();
        bv.clear(true);


        if (bv.size() != 60000100)
        {
            cout << "Resize failed" << endl;
            exit(1);
        }
        if (bv.count() != 0)
        {
            cout << "Resize count failed" << endl;
            exit(1);
        }
    }
    
    {
        bvect bv(100);
        assert(bv.size()==100);
        bv[10000000] = true;
        assert(bv.size() == 10000001);
        bv.set_bit(10000001);
        assert(bv.size() == 10000002);
    }

}


template<class A, class B> void CompareMiniSet(const A& ms,
                                          const B& bvm)
{
    for (unsigned i = 0; i < bm::set_total_blocks; ++i)
    {
        bool ms_val = ms.test(i)!=0;
        bool bvm_val = bvm.is_bit_true(i)!=0;
        if (ms_val != bvm_val)
        {
            printf("MiniSet comparison error: %u\n",i);
            exit(1);
        }
    }
}

static
void MiniSetTest()
{
    cout << "----------------------- MiniSetTest" << endl;
    {
    bm::miniset<bm::block_allocator, bm::set_total_blocks> ms;
    bvect_mini bvm(bm::set_total_blocks);


    CompareMiniSet(ms, bvm);


    ms.set(1);
    bvm.set_bit(1);

    CompareMiniSet(ms, bvm);

    unsigned i;

    for (i = 1; i < 10; i++)
    {
        ms.set(i);
        bvm.set_bit(i);
    }
    CompareMiniSet(ms, bvm);

    for (i = 1; i < 10; i++)
    {
        ms.set(i, false);
        bvm.clear_bit(i);
    }
    CompareMiniSet(ms, bvm);


    for (i = 1; i < 10; i+=3)
    {
        ms.set(i);
        bvm.set_bit(i);
    }
    CompareMiniSet(ms, bvm);

    for (i = 1; i < 5; i+=3)
    {
        ms.set(i, false);
        bvm.clear_bit(i);
    }
    CompareMiniSet(ms, bvm);
    }


    {
    bm::miniset<bm::block_allocator, bm::set_total_blocks> ms;
    bvect_mini bvm(bm::set_total_blocks);


    ms.set(1);
    bvm.set_bit(1);

    CompareMiniSet(ms, bvm);

    unsigned i;
    for (i = 1; i < bm::set_total_blocks; i+=3)
    {
        ms.set(i);
        bvm.set_bit(i);
    }
    CompareMiniSet(ms, bvm);

    for (i = 1; i < bm::set_total_blocks/2; i+=3)
    {
        ms.set(i, false);
        bvm.clear_bit(i);
    }
    CompareMiniSet(ms, bvm);
    }


    {
    bm::bvmini<bm::set_total_blocks> ms(0);
    bvect_mini bvm(bm::set_total_blocks);


    CompareMiniSet(ms, bvm);


    ms.set(1);
    bvm.set_bit(1);

    CompareMiniSet(ms, bvm);

    unsigned i;

    for (i = 1; i < 10; i++)
    {
        ms.set(i);
        bvm.set_bit(i);
    }
    CompareMiniSet(ms, bvm);

    for (i = 1; i < 10; i++)
    {
        ms.set(i, false);
        bvm.clear_bit(i);
    }
    CompareMiniSet(ms, bvm);


    for (i = 1; i < bm::set_total_blocks; i+=3)
    {
        ms.set(i);
        bvm.set_bit(i);
    }
    CompareMiniSet(ms, bvm);

    for (i = 1; i < bm::set_total_blocks/2; i+=3)
    {
        ms.set(i, false);
        bvm.clear_bit(i);
    }
    CompareMiniSet(ms, bvm);
    }


    {
    bm::miniset<bm::block_allocator, bm::set_total_blocks> ms;
    bvect_mini bvm(bm::set_total_blocks);


    ms.set(1);
    bvm.set_bit(1);

    CompareMiniSet(ms, bvm);

    unsigned i;
    for (i = 1; i < 15; i+=3)
    {
        ms.set(i);
        bvm.set_bit(i);
    }
    CompareMiniSet(ms, bvm);

    for (i = 1; i < 7; i+=3)
    {
        ms.set(i, false);
        bvm.clear_bit(i);
    }
    CompareMiniSet(ms, bvm);
    }


    cout << "----------------------- MiniSetTest ok" << endl;
}

inline
unsigned CalcBitCount32(unsigned b)
{
    b = (b & 0x55555555) + (b >> 1 & 0x55555555);
    b = (b & 0x33333333) + (b >> 2 & 0x33333333);
    b = b + ((b >> 4) & 0x0F0F0F0F);
    b = b + (b >> 8);
    b = b + ((b >> 16) & 0x0000003F);
    return b;
}

static
void PrintGapLevels(const gap_word_t* glevel)
{
    cout << "Gap levels:" << endl;
    unsigned i;
    for (i = 0; i < bm::gap_levels; ++i)
    {
        cout << glevel[i] << ",";
    }
    cout << endl;
}

static
void OptimGAPTest()
{
    gap_word_t    glevel[bm::gap_levels];
    ::memcpy(glevel, gap_len_table<true>::_len, bm::gap_levels * sizeof(gap_word_t));

    {
    gap_word_t  length[] = { 2, 2, 5, 5, 10, 11, 12 };
    unsigned lsize = sizeof(length) / sizeof(gap_word_t);

    bm::improve_gap_levels(length, length + lsize, glevel);

    PrintGapLevels(glevel);
    }

    {
    gap_word_t  length[] = { 3, 5, 15, 15, 100, 110, 120 };
    unsigned lsize = sizeof(length) / sizeof(gap_word_t);

    bm::improve_gap_levels(length, length + lsize, glevel);
    PrintGapLevels(glevel);
    }

    {
    gap_word_t  length[] = { 15, 80, 5, 3, 100, 110, 95 };
    unsigned lsize = sizeof(length) / sizeof(gap_word_t);

    bm::improve_gap_levels(length, length + lsize, glevel);
    PrintGapLevels(glevel);
    }

    {
    gap_word_t  length[] = 
    { 16,30,14,24,14,30,18,14,12,16,8,38,28,4,20,18,28,22,32,14,12,16,10,8,14,18,14,8,
      16,30,8,8,58,28,18,4,26,14,52,12,18,10,14,18,22,18,20,70,12,6,26,6,8,22,12,4,8,8,
      8,54,18,6,8,4,4,10,4,4,4,4,4,6,22,14,38,40,56,50,6,10,8,18,82,16,6,18,20,12,12,
      16,8,14,14,10,16,12,10,16,14,12,18,14,18,34,14,12,18,18,10,20,10,18,8,14,14,22,16,
      10,10,18,8,20,14,10,14,12,12,14,16,16,6,10,14,6,10,10,10,10,12,4,8,8,8,10,10,8,
      8,12,10,10,14,14,14,8,4,4,10,10,4,10,4,8,6,52,104,584,218
    };
    unsigned lsize = sizeof(length) / sizeof(gap_word_t);

    bm::improve_gap_levels(length, length + lsize, glevel);
    PrintGapLevels(glevel);
    }

    {
    gap_word_t  length[] = {
     30,46,26,4,4,68,72,6,10,4,6,14,6,42,198,22,12,4,6,24,12,8,18,4,6,10,6,4,6,6,12,6
    ,6,4,4,78,38,8,52,4,8,10,6,8,8,6,10,4,6,6,4,10,6,8,16,22,28,14,10,10,16,10,20,10
    ,14,12,8,18,4,8,10,6,10,4,6,12,16,12,6,4,8,4,14,14,6,8,4,10,10,8,8,6,8,6,8,4,8,4
    ,8,10,6,4,6 
    };
    unsigned lsize = sizeof(length) / sizeof(gap_word_t);

    bm::improve_gap_levels(length, length + lsize, glevel);
    PrintGapLevels(glevel);

    }

}


static
void BitCountChangeTest()
{
    cout << "---------------------------- BitCountChangeTest " << endl;

    unsigned i;
    for (i = 0xFFFFFFFF; i; i <<= 1)
    {
        unsigned a0 = bm::bit_count_change(i);
        unsigned a1 = BitCountChange(i);

        if (a0 != a1)
        {
            cout << hex
                << "Bit count change test failed!"
                << " i = " << i << " return = "
                << a0 << " check = " << a1
                << endl;
            exit(1);
        }
    }

    cout << "---------------------------- STEP 2 " << endl;

    for (i = 0xFFFFFFFF; i; i >>= 1)
    {
        unsigned a0 = bm::bit_count_change(i);
        unsigned a1 = BitCountChange(i);

        if (a0 != a1)
        {
            cout << "Bit count change test failed!"
                << " i = " << i << " return = "
                << a0 << " check = " << a1
                << endl;
            exit(1);
        }
    }

    cout << "---------------------------- STEP 3 " << endl;

    for (i = 0; i < 0xFFFFFFF; ++i)
    {
        unsigned a0 = bm::bit_count_change(i);
        unsigned a1 = BitCountChange(i);

        if (a0 != a1)
        {
            cout << "Bit count change test failed!"
                << " i = " << i << " return = "
                << a0 << " check = " << a1
                << endl;
            exit(1);
        }
    }
    cout << "!" << endl;

    bm::word_t  BM_VECT_ALIGN arr0[32] BM_VECT_ALIGN_ATTR = { 0, };
    arr0[0] = (bm::word_t)(1 << 31);
    arr0[1] = 1; //(bm::word_t)(1 << 31);

    bm::id_t cnt;

    cnt = bm::bit_count_change(arr0[1]);
    cout << cnt << endl;
    if (cnt != 2)
    {
        cout << "0.count_change() failed " << cnt << endl;
        exit(1);
    }
    
    // check solid block
    {
        BM_DECLARE_TEMP_BLOCK(tb1);
        for (i = 0; i < bm::set_block_size; ++i)
        {
            tb1.b_.w32[i] = 0;
        }
        unsigned gap_count = bm::bit_block_calc_change(tb1);
        assert(gap_count == 1);
        for (i = 0; i < bm::set_block_size; ++i)
        {
            tb1.b_.w32[i] = ~0u;
        }
        gap_count = bm::bit_block_calc_change(tb1);
        assert(gap_count == 1);
    }

    cout << "---------------------------- STEP 4 " << endl;
    
    bvect   bv1;
    cnt = bm::count_intervals(bv1);

    if (cnt != 1)
    {
        cout << "1.count_intervals() failed " << cnt << endl;
        exit(1);
    }
    CheckIntervals(bv1, 65536);

    bv1.invert();

    cnt = count_intervals(bv1);
    cout << "Inverted cnt=" << cnt << endl;

    if (cnt != 2)
    {
        cout << "2.inverted count_intervals() failed " << cnt << endl;
        exit(1);
    }

    bv1.invert();

    for (i = 10; i < 100000; ++i)
    {
        bv1.set(i);
    }

    cnt = count_intervals(bv1);

    if (cnt != 3)
    {
        cout << "3.count_intervals() failed " << cnt << endl;
        exit(1);
    }
    cout << "-----" << endl;
    CheckIntervals(bv1, 65536 * 2);
    cout << "Optmization..." << endl;
    bv1.optimize();
    cnt = count_intervals(bv1);

    if (cnt != 3)
    {
        cout << "4.count_intervals() failed " << cnt << endl;
        exit(1);
    }

    CheckIntervals(bv1, 65536 * 2);
    
    cout << "---------------------------- array GAP test " << endl;

    {
        bm::gap_word_t arr[] = { 0 };

        unsigned gap_count;
        
        gap_count = bit_array_compute_gaps(arr, sizeof(arr)/sizeof(arr[0]));
        if (gap_count != 1)
        {
            cout << "Array gap test failed. 1. " << endl;
            exit(1);
        }

        bm::gap_word_t gap[20] = {0};
        bm::gap_word_t gap_cntrl[20] = {0};

        gap_set_all(gap_cntrl, bm::gap_max_bits, 0);
        for (i = 0; i < sizeof(arr)/sizeof(arr[0]); ++i)
        {
            unsigned is_set;
            gap_set_value(1, gap_cntrl, arr[i], &is_set);
        }
        unsigned gap_l_cntrl = gap_length(gap_cntrl);
        unsigned gap_len = gap_set_array(&gap[0], arr, sizeof(arr)/sizeof(arr[0]));
        unsigned gap_len1 = gap_length(gap);

        if (gap_len != gap_l_cntrl || gap_len1 != gap_l_cntrl)
        {
            cout << "Array gap test failed. 1. " << endl;
            exit(1);
        }
        int cmpres = gapcmp(gap, gap_cntrl);
        if (cmpres != 0)
        {
            cout << "Array gap cmp test failed. 1. " << endl;
            exit(1);
        }
    }

    {
        bm::gap_word_t arr[] = { 65535 };

        unsigned gap_count;

        gap_count = bit_array_compute_gaps(arr, sizeof(arr)/sizeof(arr[0]));
        if (gap_count != 2)
        {
            cout << "Array gap test failed. 1.1 " << endl;
            exit(1);
        }

        bm::gap_word_t gap[20] = {0};
        bm::gap_word_t gap_cntrl[20] = {0};

        gap_set_all(gap_cntrl, bm::gap_max_bits, 0);
        for (i = 0; i < sizeof(arr)/sizeof(arr[0]); ++i)
        {
            unsigned is_set;
            gap_set_value(1, gap_cntrl, arr[i], &is_set);
        }
        unsigned gap_l_cntrl = gap_length(gap_cntrl);

        unsigned gap_len = gap_set_array(&gap[0], arr, sizeof(arr)/sizeof(arr[0]));
        unsigned gap_len1 = gap_length(gap);

        if (gap_len != gap_l_cntrl || gap_len1 != gap_l_cntrl)
        {
            cout << "Array gap test failed. 1.1 " << endl;
            exit(1);
        }
        int cmpres = gapcmp(gap, gap_cntrl);
        if (cmpres != 0)
        {
            cout << "Array gap cmp test failed. 1. " << endl;
            exit(1);
        }
    }

    {
        bm::gap_word_t arr[] = { 0, 65535 };

        unsigned gap_count;

        gap_count = bit_array_compute_gaps(arr, sizeof(arr)/sizeof(arr[0]));
        if (gap_count != 3)
        {
            cout << "Array gap test failed. 1.2 " << endl;
            exit(1);
        }

        bm::gap_word_t gap[20] = {0};
        bm::gap_word_t gap_cntrl[20] = {0};

        gap_set_all(gap_cntrl, bm::gap_max_bits, 0);
        for (i = 0; i < sizeof(arr)/sizeof(arr[0]); ++i)
        {
            unsigned is_set;
            gap_set_value(1, gap_cntrl, arr[i], &is_set);
        }
        unsigned gap_l_cntrl = gap_length(gap_cntrl);

        unsigned gap_len = gap_set_array(&gap[0], arr, sizeof(arr)/sizeof(arr[0]));
        unsigned gap_len1 = gap_length(gap);

        if (gap_len != gap_l_cntrl || gap_len1 != gap_l_cntrl)
        {
            cout << "Array gap test failed. 1.2 " << endl;
            exit(1);
        }
        int cmpres = gapcmp(gap, gap_cntrl);
        if (cmpres != 0)
        {
            cout << "Array gap cmp test failed. 1.2 " << endl;
            exit(1);
        }
    }

    {
        bm::gap_word_t arr[] = { 0, 1, 2, 65534, 65535 };

        unsigned gap_count;

        gap_count = bit_array_compute_gaps(arr, sizeof(arr)/sizeof(arr[0]));
        if (gap_count != 3)
        {
            cout << "Array gap test failed. 1.3 " << endl;
            exit(1);
        }

        bm::gap_word_t gap[20] = {0};
        bm::gap_word_t gap_cntrl[20] = {0};

        gap_set_all(gap_cntrl, bm::gap_max_bits, 0);
        for (i = 0; i < sizeof(arr)/sizeof(arr[0]); ++i)
        {
            unsigned is_set;
            gap_set_value(1, gap_cntrl, arr[i], &is_set);
        }
        unsigned gap_l_cntrl = gap_length(gap_cntrl);

        unsigned gap_len = gap_set_array(&gap[0], arr, sizeof(arr)/sizeof(arr[0]));
        unsigned gap_len1 = gap_length(gap);

        if (gap_len != gap_l_cntrl || gap_len1 != gap_l_cntrl)
        {
            cout << "Array gap test failed. 1.3 " << endl;
            exit(1);
        }
        int cmpres = gapcmp(gap, gap_cntrl);
        if (cmpres != 0)
        {
            cout << "Array gap cmp test failed. 1.3 " << endl;
            exit(1);
        }
    }

    {
        bm::gap_word_t arr[] = { 0, 1, 2 };
        unsigned gap_count;

        gap_count = bit_array_compute_gaps(arr, sizeof(arr)/sizeof(arr[0]));
        if (gap_count != 1)
        {
            cout << "Array gap test failed. 2. " << endl;
            exit(1);
        }
        bm::gap_word_t gap[20] = {0};
        bm::gap_word_t gap_cntrl[20] = {0};

        gap_set_all(gap_cntrl, bm::gap_max_bits, 0);
        for (i = 0; i < sizeof(arr)/sizeof(arr[0]); ++i)
        {
            unsigned is_set;
            gap_set_value(1, gap_cntrl, arr[i], &is_set);
        }
        unsigned gap_l_cntrl = gap_length(gap_cntrl);

        unsigned gap_len = gap_set_array(&gap[0], arr, sizeof(arr)/sizeof(arr[0]));
        unsigned gap_len1 = gap_length(gap);

        if (gap_len != gap_l_cntrl || gap_len1 != gap_l_cntrl)
        {
            cout << "Array gap test failed. 2 " << endl;
            exit(1);
        }
        int cmpres = gapcmp(gap, gap_cntrl);
        if (cmpres != 0)
        {
            cout << "Array gap cmp test failed. 2 " << endl;
            exit(1);
        }

    }

    {
        bm::gap_word_t arr[] = { 1, 2 };
        unsigned gap_count;

        gap_count = bit_array_compute_gaps(arr, sizeof(arr)/sizeof(arr[0]));
        if (gap_count != 2)
        {
            cout << "Array gap test failed. 3. " << endl;
            exit(1);
        }
        bm::gap_word_t gap[20] = {0};
        bm::gap_word_t gap_cntrl[20] = {0};

        gap_set_all(gap_cntrl, bm::gap_max_bits, 0);
        for (i = 0; i < sizeof(arr)/sizeof(arr[0]); ++i)
        {
            unsigned is_set;
            gap_set_value(1, gap_cntrl, arr[i], &is_set);
        }
        unsigned gap_l_cntrl = gap_length(gap_cntrl);

        unsigned gap_len = gap_set_array(&gap[0], arr, sizeof(arr)/sizeof(arr[0]));
        unsigned gap_len1 = gap_length(gap);

        if (gap_len != gap_l_cntrl || gap_len1 != gap_l_cntrl)
        {
            cout << "Array gap test failed. 3 " << endl;
            exit(1);
        }
        int cmpres = gapcmp(gap, gap_cntrl);
        if (cmpres != 0)
        {
            cout << "Array gap cmp test failed. 3 " << endl;
            exit(1);
        }
    }

    {
        bm::gap_word_t arr[] = { 1, 2, 10 };
        unsigned gap_count;

        gap_count = bit_array_compute_gaps(arr, sizeof(arr)/sizeof(arr[0]));
        if (gap_count != 4)
        {
            cout << "Array gap test failed. 4. " << endl;
            exit(1);
        }
        bm::gap_word_t gap[20] = {0};
        bm::gap_word_t gap_cntrl[20] = {0};

        gap_set_all(gap_cntrl, bm::gap_max_bits, 0);
        for ( i = 0; i < sizeof(arr)/sizeof(arr[0]); ++i)
        {
            unsigned is_set;
            gap_set_value(1, gap_cntrl, arr[i], &is_set);
        }
        unsigned gap_l_cntrl = gap_length(gap_cntrl);

        unsigned gap_len = gap_set_array(&gap[0], arr, sizeof(arr)/sizeof(arr[0]));
        unsigned gap_len1 = gap_length(gap);

        if (gap_len != gap_l_cntrl || gap_len1 != gap_l_cntrl)
        {
            cout << "Array gap test failed. 4 " << endl;
            exit(1);
        }
        int cmpres = gapcmp(gap, gap_cntrl);
        if (cmpres != 0)
        {
            cout << "Array gap cmp test failed. 4 " << endl;
            exit(1);
        }
    }

    {
        bm::gap_word_t arr[] = { 1, 2, 10, 11 };
        unsigned gap_count;

        gap_count = bit_array_compute_gaps(arr, sizeof(arr)/sizeof(arr[0]));
        if (gap_count != 4)
        {
            cout << "Array gap test failed. 5. " << endl;
            exit(1);
        }
        bm::gap_word_t gap[20] = {0};
        bm::gap_word_t gap_cntrl[20] = {0};

        gap_set_all(gap_cntrl, bm::gap_max_bits, 0);
        for ( i = 0; i < sizeof(arr)/sizeof(arr[0]); ++i)
        {
            unsigned is_set;
            gap_set_value(1, gap_cntrl, arr[i], &is_set);
        }
        unsigned gap_l_cntrl = gap_length(gap_cntrl);

        unsigned gap_len = gap_set_array(&gap[0], arr, sizeof(arr)/sizeof(arr[0]));
        unsigned gap_len1 = gap_length(gap);

        if (gap_len != gap_l_cntrl || gap_len1 != gap_l_cntrl)
        {
            cout << "Array gap test failed. 5 " << endl;
            exit(1);
        }
        int cmpres = gapcmp(gap, gap_cntrl);
        if (cmpres != 0)
        {
            cout << "Array gap cmp test failed. 5 " << endl;
            exit(1);
        }

    }

    {
        bm::gap_word_t arr[] = { 1, 2, 10, 11, 256 };
        unsigned gap_count;

        gap_count = bit_array_compute_gaps(arr, sizeof(arr)/sizeof(arr[0]));
        if (gap_count != 6)
        {
            cout << "Array gap test failed. 6. " << endl;
            exit(1);
        }
        bm::gap_word_t gap[20] = {0};
        bm::gap_word_t gap_cntrl[20] = {0};

        gap_set_all(gap_cntrl, bm::gap_max_bits, 0);
        for ( i = 0; i < sizeof(arr)/sizeof(arr[0]); ++i)
        {
            unsigned is_set;
            gap_set_value(1, gap_cntrl, arr[i], &is_set);
        }
        unsigned gap_l_cntrl = gap_length(gap_cntrl);

        unsigned gap_len = gap_set_array(&gap[0], arr, sizeof(arr)/sizeof(arr[0]));
        unsigned gap_len1 = gap_length(gap);

        if (gap_len != gap_l_cntrl || gap_len1 != gap_l_cntrl)
        {
            cout << "Array gap test failed. 6 " << endl;
            exit(1);
        }
        int cmpres = gapcmp(gap, gap_cntrl);
        if (cmpres != 0)
        {
            cout << "Array gap cmp test failed. 6 " << endl;
            exit(1);
        }

    }


    cout << "---------------------------- BitCountChangeTest Ok." << endl;
}



static
void DNACompressionTest()
{
    const char seeds[] = 
        { 'A', 'C', 'G', 'T', 'A', 'C', 'G', 'A', 'N', 'A', 'C', 'G' };
    
    const unsigned arr_size = bm::set_block_size*4;
    const unsigned arr_plain_size = arr_size / 8;    
    
    unsigned char BM_VECT_ALIGN block1[arr_size] BM_VECT_ALIGN_ATTR = {0,};

    unsigned char BM_VECT_ALIGN tmatrix1[8][arr_plain_size] BM_VECT_ALIGN_ATTR;
    unsigned BM_VECT_ALIGN distance1[8][8] BM_VECT_ALIGN_ATTR;
    unsigned char pc_vector1[8] = {0,};
    unsigned pc_vector_stat1[bm::ibpc_end];
/*
    unsigned   BM_ALIGN16 tmatrix2[32][bm::set_block_plain_size] BM_ALIGN16ATTR;
    unsigned  
    BM_ALIGN16 distance2[bm::set_block_plain_cnt][bm::set_block_plain_cnt] BM_ALIGN16ATTR;
    unsigned char pc_vector2[32] = {0,};


    unsigned   BM_ALIGN16 tmatrix3[32][bm::set_block_plain_size] BM_ALIGN16ATTR;
    unsigned  
    BM_ALIGN16 distance3[bm::set_block_plain_cnt][bm::set_block_plain_cnt] BM_ALIGN16ATTR;
    unsigned char pc_vector3[32] = {0,};
*/
    
    // generate pseudo-random DNA sequence
    for (unsigned i = 0; i < arr_size; ++i)
    {
        unsigned letter_idx = unsigned(rand()) % unsigned(sizeof(seeds));
        unsigned char l = (unsigned char)seeds[letter_idx];
        unsigned char c = 0;
        switch (l)
        {
        case 'A':
            c = 0; break;
        case 'C':
            c = 1; break;
        case 'G':
            c = 2; break;
        case 'T':
            c = 3; break;
        case 'N':
            c = 4; break;
        default:
            cout << "Alphabet error!" << endl;
            exit(1);
        };
        block1[i] = c;
        //cout << block1[i];
    }
    cout << endl;
        
    bm::vect_bit_transpose<unsigned char, 
                           8, 
                           arr_plain_size>
                           (block1, arr_size, tmatrix1);
    
    bm::tmatrix_distance<unsigned char, 
                         8, 
                         arr_plain_size>
                         (tmatrix1, distance1);
    
    cout << "ALL count=" << sizeof(char)*8*arr_plain_size << endl;
    bm::bit_iblock_make_pcv<unsigned char, 8, arr_plain_size>(distance1, pc_vector1);
    
    bm::bit_iblock_pcv_stat(pc_vector1, pc_vector1 + 8, pc_vector_stat1);
    
    for (unsigned s = 0; s < bm::ibpc_end; ++s)
    {
        switch(s)
        {
        case bm::ibpc_uncompr:
            cout << "Uncompressed: "; 
            break;
        case bm::ibpc_all_zero:
            cout << "    All ZERO: "; 
            break;
        case bm::ibpc_all_one:
            cout << "     All ONE: "; 
            break;
        case bm::ibpc_equiv:
            cout << "       Equiv: "; 
            break;
        case bm::ibpc_close:
            cout << "     Similar: "; 
            break;
        default:
            //cout << "Oops!" << s << " "; 
            break;
        }
        cout << pc_vector_stat1[s] << endl;
    } // for
    

    // print out the pc_vector    
    for (unsigned j = 0; j < 8; ++j)
    {
        unsigned ibpc = pc_vector1[j] & 7;
        unsigned n_row = (pc_vector1[j] >> 3);
        cout << j << ":" << "->" << n_row << " ";
        
        switch(ibpc)
        {
        case bm::ibpc_uncompr:
            cout << "Uncompressed: "; 
            cout << " popcnt=" << distance1[j][j];
            break;
        case bm::ibpc_all_zero:
            cout << "ZERO";
            break;            
        case bm::ibpc_all_one:
            cout << "ONE: "; 
            break;
        case bm::ibpc_equiv:
            cout << "Equiv: "; 
            break;            
        case bm::ibpc_close:
            cout << " Similar: "; 
            cout << " popcnt="  << distance1[j][j]
                 << " Humming=" << distance1[j][n_row];            
            break;
        default:
            assert(0);
        }
        cout << endl;
    }
/*
    cout << endl << "Second round." << endl << endl;

    bm::bit_iblock_reduce(tmatrix1, pc_vector1, pc_vector1+32, tmatrix2);
    bm::tmatrix_distance<unsigned, 
                         bm::set_block_plain_cnt, 
                         bm::set_block_plain_size>
                         (tmatrix2, distance2);    
    
    bm::bit_iblock_make_pcv(distance2, pc_vector2);
    
    // print out the pc_vector    
    for (unsigned j = 0; j < 32; ++j)
    {
        unsigned ibpc = pc_vector2[j] & 7;
        unsigned n_row = (pc_vector2[j] >> 3);
        cout << j << ":" << "->" << n_row << " ";
        
        switch(ibpc)
        {
        case bm::ibpc_uncompr:
            cout << "Uncompressed: "; 
            cout << " popcnt=" << distance2[j][j];
            break;
        case bm::ibpc_all_zero:
            cout << "ZERO";
            break;            
        case bm::ibpc_all_one:
            cout << "ONE: "; 
            break;
        case bm::ibpc_equiv:
            cout << "Equiv: "; 
            break;            
        case bm::ibpc_close:
            cout << " Similar: "; 
            cout << " popcnt="  << distance2[j][j]
                 << " Humming=" << distance2[j][n_row] << endl; 
             {
                const unsigned* r1 = tmatrix2[j];
                for (unsigned i = 0; i < bm::set_block_plain_size; ++i)
                {
                    cout << hex << r1[i] << " ";
                }
                cout << dec << endl << endl;                         
             }           
            break;
        }
        cout << endl;
    }


    cout << endl << "3rd round." << endl << endl;

    bm::bit_iblock_reduce(tmatrix2, pc_vector2, pc_vector2+32, tmatrix3);

    bm::tmatrix_distance<unsigned, 
                         bm::set_block_plain_cnt, 
                         bm::set_block_plain_size>
                         (tmatrix3, distance3);    
    
    bm::bit_iblock_make_pcv(distance3, pc_vector3);
    
    // print out the pc_vector    
    for (unsigned j = 0; j < 32; ++j)
    {
        unsigned ibpc = pc_vector3[j] & 7;
        unsigned n_row = (pc_vector3[j] >> 3);
        cout << j << ":" << "->" << n_row << " ";
        
        switch(ibpc)
        {
        case bm::ibpc_uncompr:
            cout << "Uncompressed: "; 
            cout << " popcnt=" << distance3[j][j];
            break;
        case bm::ibpc_all_zero:
            cout << "ZERO";
            break;            
        case bm::ibpc_all_one:
            cout << "ONE: "; 
            break;
        case bm::ibpc_equiv:
            cout << "Equiv: "; 
            break;            
        case bm::ibpc_close:
            cout << " Similar: "; 
            cout << " popcnt="  << distance3[j][j]
                 << " Humming=" << distance3[j][n_row] << endl; 
             {
                const unsigned* r1 = tmatrix3[j];
                for (unsigned i = 0; i < bm::set_block_plain_size; ++i)
                {
                    cout << hex << r1[i] << " ";
                }
                cout << dec << endl << endl;                         
             }           
            break;
        }
        cout << endl;
    }
*/    
    
}

void BitBlockTransposeTest();

void BitBlockTransposeTest()
{
    DNACompressionTest();
   

    bm::word_t BM_ALIGN16 block1[bm::set_block_size] BM_ALIGN16ATTR = {0,};
    bm::word_t BM_ALIGN16 block2[bm::set_block_size] BM_ALIGN16ATTR = {0xFF,};
    unsigned   BM_ALIGN16 tmatrix1[32][bm::set_block_plain_size] BM_ALIGN16ATTR;


    cout << "---------------------------- BitTransposeTest" << endl;

    cout << "Transpose 1" << endl;

    for (unsigned i = 0; i < bm::set_block_size; ++i)
    {
        block1[i] = 1;
    }

    bm::vect_bit_transpose<unsigned, 
                           bm::set_block_plain_cnt, 
                           bm::set_block_plain_size>
                           (block1, bm::set_block_size, tmatrix1);

    bm::vect_bit_trestore<unsigned, 
                           bm::set_block_plain_cnt, 
                           bm::set_block_plain_size>
                           (tmatrix1, block2);

    for (unsigned i = 0; i < bm::set_block_size; ++i)
    {
        if (block1[i] != block2[i])
        {
            cout << "Bit transpose error! " << i << endl; exit(1);
        }
    }

    {
    unsigned BM_ALIGN16 distance[bm::set_block_plain_cnt][bm::set_block_plain_cnt];
    bm::tmatrix_distance<unsigned, 
                         bm::set_block_plain_cnt, 
                         bm::set_block_plain_size>
                         (tmatrix1, distance);
    
    PrintDistanceMatrix(distance);

    // distance matrix verification:
    {
    for (unsigned i = 0; i < bm::set_block_plain_cnt; ++i)
    {
        const unsigned* row = distance[i];
        for (unsigned j = i; j < bm::set_block_plain_cnt; ++j)
        {
            if (i == j)
            {
                if (distance[0][0] != 2048)
                {
                    cout << "Self distance(bitcount) is incorrect!" << endl;
                    exit(1);
                }
            }
            else
            {
                if (i == 0)
                {
                    if (row[j] != 2048) // max distance
                    {
                        cout << "Incorrect max distance!" << endl; exit(1);
                    }
                }
                else
                {
                    if (row[j] != 0) // max distance
                    {
                        cout << "Incorrect min distance!" << endl; exit(1);
                    }
                }
            }
        }
    }
    }

    }

    cout << "Transpose 2" << endl;

    for (unsigned i = 0; i < bm::set_block_size; ++i)
    {
        block1[i] = 1 | (1 << 17);
    }

    bm::vect_bit_transpose<unsigned, 
                           bm::set_block_plain_cnt, 
                           bm::set_block_plain_size>
                           (block1, bm::set_block_size, tmatrix1);
    bm::vect_bit_trestore<unsigned, 
                           bm::set_block_plain_cnt, 
                           bm::set_block_plain_size>
                           (tmatrix1, block2);


    for (unsigned i = 0; i < bm::set_block_size; ++i)
    {
        if (block1[i] != block2[i])
        {
            cout << "Bit transpose error! " << i << endl; exit(1);
        }
    }

    cout << "Transpose 3" << endl;

    for (unsigned i = 0; i < bm::set_block_size; ++i)
    {
        block1[i] = ~1u;
    }

    bm::vect_bit_transpose<unsigned, 
                           bm::set_block_plain_cnt, 
                           bm::set_block_plain_size>
                           (block1, bm::set_block_size, tmatrix1);
    bm::vect_bit_trestore<unsigned, 
                           bm::set_block_plain_cnt, 
                           bm::set_block_plain_size>
                           (tmatrix1, block2);

    for (unsigned i = 0; i < bm::set_block_size; ++i)
    {
        if (block1[i] != block2[i])
        {
            cout << "Bit transpose error! " << i << endl; exit(1);
        }
    }

    cout << "Transpose 4" << endl;

    for (unsigned i = 0; i < bm::set_block_size; ++i)
    {
        block1[i] = i;
    }

    bm::vect_bit_transpose<unsigned, 
                           bm::set_block_plain_cnt, 
                           bm::set_block_plain_size>
                           (block1, bm::set_block_size, tmatrix1);
    bm::vect_bit_trestore<unsigned, 
                           bm::set_block_plain_cnt, 
                           bm::set_block_plain_size>
                           (tmatrix1, block2);

    for (unsigned i = 0; i < bm::set_block_size; ++i)
    {
        if (block1[i] != block2[i])
        {
            cout << "Bit transpose error! " << i << endl; exit(1);
        }
    }
/*    
    cout << "Transpose 5 - random" << endl;

    for (unsigned c = 0; c < 10000; ++c)
    {
        if ((c % 100) == 0) cout << ".";

        for (unsigned i = 0; i < bm::set_block_size; ++i)
        {
            block1[i] = rand();
        }

        bm::vect_bit_transpose<unsigned, 
                               bm::set_block_plain_cnt, 
                               bm::set_block_plain_size>
                               (block1, bm::set_block_size, tmatrix1);

        bm::vect_bit_trestore<unsigned, 
                               bm::set_block_plain_cnt, 
                               bm::set_block_plain_size>
                               (tmatrix1, block2);


        for (unsigned i = 0; i < bm::set_block_size; ++i)
        {
            if (block1[i] != block2[i])
            {
                cout << "Bit transpose error! " << i << endl; exit(1);
            }
        }
    }
 */   

    cout << "Transpose GAP block 1" << endl;
    
    {
    gap_vector   gapv(0);
    gap_vector   gapv1(9);
    gapv.set_bit(1);
    gapv.set_bit(2);
    gapv.set_bit(10);
    gapv.set_bit(65000);
    

    gap_transpose_engine<bm::gap_word_t, bm::word_t, bm::set_block_size> gte;
    
    if ( bm::conditional<sizeof(gte.tmatrix_) != (2048 * sizeof(bm::word_t))>::test())
    {
        cout << "TMatrix recalculation error!" << sizeof(gte.tmatrix_) << endl;
        exit(1);
    }
//    gte.transpose(gapv.get_buf());//, block1);


    gte.compute_distance_matrix();
    gte.reduce();
    gte.restore();
    
    unsigned glen = *(gapv.get_buf()) >> 3;
    PrintGap(gapv.get_buf());
    PrintDGap((gap_word_t*) block1, glen-1);
    PrintDGapGamma((gap_word_t*) block1, glen-1);
    
    PrintTMatrix(gte.tmatrix_, gte.eff_cols_, true);
    
    //bm::gap_word_t gap_head = *gapv.get_buf();
//    gte.trestore(gap_head, gapv1.get_buf());//, block2);
/*    
    if (gapv.compare(gapv1))
    {
        cout << "GAP block transpose error!" << endl;
        PrintGap(gapv.get_buf());
        PrintGap(gapv1.get_buf());
        exit(1);
    }
*/    
    }

    cout << "Transpose GAP block 2" << endl;

    {
    gap_vector   gapv(0);
    gap_vector   gapv1(0);

    unsigned gcnt = 5;
    for (unsigned i = 0; i < 65500; i+= 50)
    {
        for (unsigned j = 0; j < gcnt ; ++j)
        {
            gapv.set_bit(i);

            if (++i > 65500) 
                break;
        }
        gcnt += 2;
    }

    gap_transpose_engine<bm::gap_word_t, bm::word_t, bm::set_block_size> gte;    
//    gte.transpose(gapv.get_buf());


    gte.compute_distance_matrix();
    gte.reduce();
    gte.restore();
    
    unsigned glen = *(gapv.get_buf()) >> 3;
    cout << glen << endl;

    // bm::gap_word_t gap_head = *gapv.get_buf();
//    gte.trestore(gap_head, gapv1.get_buf());
    
/*
    if (gapv.compare(gapv1))
    {
        cout << "GAP block transpose error!" << endl;
        PrintGap(gapv.get_buf());
        PrintGap(gapv1.get_buf());
        exit(1);
    }
*/

    }
    

    cout << endl << "---------------------------- BitTransposeTest ok" << endl;
}

/*
#define POWER_CHECK(w, mask) \
    (bm::bit_count_table<true>::_count[(w&mask) ^ ((w&mask)-1)])

void BitListTest()
{
    unsigned bits[64] = {0,};

    unsigned w = 0;

    w = (1 << 3) | 1;


    int bn3 = POWER_CHECK(w, 1 << 3) - 1;
    int bn2 = POWER_CHECK(w, 1 << 2) - 1;
    int bn0 = POWER_CHECK(w, 1 << 0) - 1;

    bit_list(w, bits+1);
  
}
*/

static
void ResizeTest()
{
    {{
        bvect bv(0);
        assert(bv.any() == false);
        assert(bv.count() == 0);
    }}

    {{
        bvect bv1(10);
        bvect bv2(bv1);
        assert(bv1.size() == bv2.size());
    }}

    {{
        bvect bv;
        bv.invert();
        bvect::size_type cnt = bv.count();
        assert(cnt == bm::id_max);
        assert(bv.test(bm::id_max-1));
    }}

    {{
        bvect bv(10);
        assert(bv.any() == false);
        assert(bv.count() == 0);
        bv.invert();
        unsigned cnt = bv.count();
        assert(cnt == 10);
    }}

    {{
        bvect bv1(10);
        bv1.set(1);
        bvect bv2(0);

        assert(bv1.size() == 10);
        assert(bv2.size() == 0);
        assert(bv1.count() == 1);
        assert(bv2.count() == 0);
        
        bv1.swap(bv2);

        assert(bv2.size() == 10);
        assert(bv2.count() == 1);
        assert(bv1.size() == 0);
        assert(bv1.count() == 0);
    }}

    {{
        bvect bv1;
        bv1.set(65536);
        bv1.set(100);
        assert(bv1.size() == bm::id_max);
        assert(bv1.count() == 2);
        bv1.resize(101);
        assert(bv1.size() == 101);
        assert(bv1.count() == 1);
        {{
            auto f = bv1.get_first();
            assert(f == 100);
            f = bv1.get_next(f);
            assert(f == 0);
        }}

        bv1.resize(10);
        assert(bv1.size() == 10);
        assert(bv1.count() == 0);
        auto f = bv1.get_first();
        assert(f == 0);
    }}

    {{
        bvect bv;
        print_stat(bv);
        bv.set(100);
        bv.set(65536 + 10);
        print_stat(bv);
        bv.set_range(0, 65536*10, false);
        print_stat(bv);
    }}

    // test logical operations

    {{
        bvect bv1(65536 * 10);
        bvect bv2(65536 * 100);
        bv1.set(5);
        bv2.set(5);
        bv2.set(65536 * 2);
        bv2 &= bv1;
        assert(bv2.size() == 65536 * 100);
        assert(bv2.count() == 1);
        assert(bv2.get_first() == 5);
    }}

    {{
        bvect bv1(10);
        bvect bv2;
        bv1.set(5);
        bv2.set(5);
        bv2.set(65536 * 2);
        bv1 &= bv2;
        assert(bv1.size() == bv2.size());
        assert(bv1.count() == 1);
        assert(bv1.get_first() == 5);
    }}

    {{
        bvect bv1(10);
        bvect bv2;
        bv1.set(5);
        bv2.set(6);
        bv2.set(65536 * 2);
        bv1 |= bv2;
        assert(bv1.size() == bv2.size());
        assert(bv1.count() == 3);
    }}

    // comparison test

    {{
        int cmp;
        bvect bv1(10);
        bvect bv2;
        bv2.set(65536 * 2);

        cmp = bv1.compare(bv2);
        assert(cmp < 0);

        bv1.set(5);
        assert(cmp < 0);
        cmp = bv1.compare(bv2);
        assert(cmp > 0);
        cmp = bv2.compare(bv1);
        assert(cmp < 0);
    }}

    // inserter

    {{
        bvect bv1(10);
        {
            bvect::insert_iterator it(bv1);
            *it = 100 * 65536;
        }
        assert(bv1.size() ==  100 * 65536 + 1);
    }}

    // serialization

    {{
        bvect bv1(10);
        bv1.set(5);
        struct bvect::statistics st1;
        bv1.calc_stat(&st1);

        unsigned char* sermem = new unsigned char[st1.max_serialize_mem];
        size_t slen2 = bm::serialize(bv1, sermem);
        cout << slen2 << endl;

        bvect bv2(0);
        bm::deserialize(bv2, sermem);
        delete [] sermem;

        assert(bv2.size() == 10);
        assert(bv2.count() == 1);
        assert(bv2.get_first() == 5);
    }}

    {{
        bvect bv1(10);
        bv1.set(5);
        unsigned int arg[] = { 10, 65536, 65537, 65538 * 10000 };
        unsigned* it1 = arg;
        unsigned* it2 = arg + 4;
        combine_or(bv1, it1, it2);
        assert(bv1.size() == 65538 * 10000 + 1);
        bvect::enumerator en = bv1.first();
        while (en.valid())
        {
            cout << *en << " ";
            ++en;
        }
    }}
}

static
void VerifyCountRange(const bvect& bv,
                      const bvect::rs_index_type& bc_arr,
                      bvect::size_type from,
                      bvect::size_type to)
{
    for (bvect::size_type i = from; i < to; ++i)
    {
        bvect::size_type cnt1 = bv.count_range(0, i);
        bvect::size_type cnt2 = bv.count_to(i, bc_arr);
        auto cnt3 = bv.count_to_test(i, bc_arr);
        
        assert(cnt1 == cnt2);
        if (cnt1 != cnt2)
        {
            cerr << "VerifyCountRange failed!" << " count_range()=" << cnt1
                << " count_to()=" << cnt2 << endl;
        }
        if (cnt3 != cnt1)
        {
            bool b = bv.test(i);
            if (b)
            {
                cerr << "VerifyCountRange failed! count_to_test()" << cnt3 << " count_range()=" << cnt1
                     << endl;
            }
        }
        
        bvect::size_type cnt4 = bv.count_range(i, to);
        bvect::size_type cnt5 = bv.count_range(i, to, bc_arr);
        if (cnt4 != cnt5)
        {
            cnt5 = bv.count_range(i, to, bc_arr);
            assert(cnt4 == cnt5);
            exit(1);
        }
    } // for
}

static
void CountRangeTest()
{
    cout << "---------------------------- CountRangeTest..." << endl;

    {{
        bvect bv1 { 0, 1 };
        bv1.set(0);
        bv1.set(1);
 
        bvect::rs_index_type bc_arr;
        bv1.build_rs_index(&bc_arr);
        assert(bc_arr.count() == 2);

        assert(bc_arr.count(0) == 2);
        for (bvect::size_type i = 1; i < bm::set_total_blocks; ++i)
        {
            assert(bc_arr.count(i) == 0);
        } // for
 
        VerifyCountRange(bv1, bc_arr, 0, 200000);
 
        bv1.optimize();
        bvect::rs_index_type bc_arr1;
        bv1.build_rs_index(&bc_arr1);
 
        assert(bc_arr.count(0) == 2);
        for (bvect::size_type i = 1; i < bm::set_total_blocks; ++i)
        {
            assert(bc_arr.count(i) == 0);
        } // for

        VerifyCountRange(bv1, bc_arr1, 0, 200000);
    }}

    {{
        bvect bv1 { bm::id_max - 100, bm::id_max - 1 };
        
        bvect::rs_index_type bc_arr;
        bv1.build_rs_index(&bc_arr);
        assert(bc_arr.count() == 2);
        
        assert(bc_arr.rcount(bm::set_total_blocks-1) == 2);
        for (bvect::size_type i = 0; i < bm::set_total_blocks-1; ++i)
        {
            assert(bc_arr.rcount(i) == 0);
        } // for
        
        VerifyCountRange(bv1, bc_arr, 0, 200000);
        
        bv1.optimize();
        bvect::rs_index_type bc_arr1;
        bv1.build_rs_index(&bc_arr1);
        
        assert(bc_arr.rcount(bm::set_total_blocks-1) == 2);
        for (bvect::size_type i = 0; i < bm::set_total_blocks-1; ++i)
        {
            assert(bc_arr1.rcount(i) == 0);
        } // for
        
        VerifyCountRange(bv1, bc_arr1, 0, 200000);
        VerifyCountRange(bv1, bc_arr, bm::id_max-200000, bm::id_max-1);
    }}

    {{
        bvect bv1 { 0, 1, 65535+10, 65535+20, 65535+21, bm::id_max-100};

 
        bvect::rs_index_type bc_arr;
        bv1.build_rs_index(&bc_arr);
        {
            auto cnt1 = bv1.count();
            auto cnt2 = bc_arr.count();
            assert(cnt1 == cnt2);
        }

        assert(bc_arr.rcount(0) == 2);
        assert(bc_arr.rcount(1) == 5);
        auto cnt = bc_arr.rcount(768);
        assert(cnt == 5);
        for (bvect::size_type i = 2; i < bm::set_total_blocks; ++i)
        {
            assert(bc_arr.rcount(i) == 5 || (bc_arr.rcount(i) == 6 && i == bm::set_total_blocks-1));
        } // for
 
        VerifyCountRange(bv1, bc_arr, bm::id_max-1, bm::id_max-1);
        for (unsigned i = 0; i < 2; ++i)
        {
            VerifyCountRange(bv1, bc_arr, 0, 200000);
            VerifyCountRange(bv1, bc_arr, bm::id_max-200000, bm::id_max-1);

            // check within empty region
            VerifyCountRange(bv1, bc_arr, bm::id_max/2-200000, bm::id_max/2+200000);

            bv1.optimize();
        }
    }}

    cout << "check 11111-filled bvector" << endl;
    {{
        bvect bv1;
        for (unsigned i = 0; i <= 200000; ++i)
            bv1.set(i, true);
        
        for (unsigned i = 0; i < 2; ++i)
        {
            bvect::rs_index_type rs_idx;
            bv1.build_rs_index(&rs_idx);

            VerifyCountRange(bv1, rs_idx, 0, 200000);

            bv1.optimize();
        }
    }}

    cout << "check inverted bvector" << endl;
    {{
            bvect bv1;
        
            bv1.invert();

            bvect::rs_index_type bc_arr;
            bv1.build_rs_index(&bc_arr);
            auto cnt1 = bv1.count();
            auto cnt2 = bc_arr.count();
            assert(cnt1 == cnt2);

            VerifyCountRange(bv1, bc_arr, bm::id_max-1, bm::id_max-1);

            VerifyCountRange(bv1, bc_arr, 0, 200000);
            VerifyCountRange(bv1, bc_arr, bm::id_max-200000, bm::id_max-1);
            VerifyCountRange(bv1, bc_arr, bm::id_max/2-200000, bm::id_max/2+200000);
    }}
    
    cout << "---------------------------- CountRangeTest OK" << endl;
}

static
void KeepRangeTest()
{
    std::cout << "--------------------- KeepRangeTest()" << endl;

    {
        bvect bv;
        bv.keep_range(10, 100);
        assert(!bv.any());
        bv.keep_range(0, 0);
        assert(!bv.any());
    }

    {
        bvect bv;
        bv.invert();
        bv.keep_range(10, 20);

        assert(bv.count() == 11);
        assert(bv.count_range(10, 20) == 11);
        bv.keep_range(20, 10);
        assert(bv.count() == 11);
        assert(bv.count_range(10, 20) == 11);

        bv.keep_range(10, 10);
        assert(bv.test(10));
        assert(bv.count() == 1);
        assert(bv.count_range(10, 10) == 1);
    }

    {
        bvect bv{ 10, 256, bm::id_max / 2, bm::id_max - 1 };
        bv.optimize();
        bv.keep_range(bm::id_max / 2 - 100, bm::id_max / 2);
        assert(bv.count() == 1);
        assert(bv.test(bm::id_max / 2));
    }


    std::cout << "--------------------- KeepRangeTest() OK" << endl;
}

static
void ExportTest()
{
    cout << "---------------------------- ExportTest..." << endl;

    {
        char buf[20] = {0,};

        buf[0] = 1;
        buf[1] = 1;
        buf[2]= (char)(1 << 1);

        bvect bv1;
        export_array(bv1, buf + 0, buf + 20);

        auto cnt = bv1.count();
        assert(cnt == 3);
        assert(bv1.test(0));
        assert(bv1.test(8));
        assert(bv1.test(17));
    }

    {
        char buf[65536*10] = {0,};

        buf[0] = 1;
        buf[1] = 1;
        buf[2]= (char)(1 << 1);

        bvect bv1;
        export_array(bv1, buf + 0, buf + 65536*10);

        assert(bv1.count() == 3);
        assert(bv1.test(0));
        assert(bv1.test(8));
        assert(bv1.test(17));
    }

    {
        short buf[20] = {0,};

        buf[0] = 1;
        buf[1] = 1;
        buf[2]= (char)(1 << 1);

        bvect bv1;
        export_array(bv1, buf + 0, buf + 20);

        assert(bv1.count() == 3);
        assert(bv1.test(0));
        assert(bv1.test(16));
        assert(bv1.test(33));
    }

    {
        int buf[20] = {0,};

        buf[0] = 1;
        buf[1] = 1;
        buf[2]= (char)(1 << 1);

        bvect bv1;
        export_array(bv1, buf + 0, buf + 20);

        assert(bv1.count() == 3);
        assert(bv1.test(0));
        assert(bv1.test(32));
        assert(bv1.test(65));
    }


    cout << "---------------------------- ExportTest Ok." << endl;
}

static
void TestRecomb()
{
    bm::word_t b1[bm::set_block_size]= {0,};
    bm::word_t b2[bm::set_block_size]= {0,};
    bm::word_t br[bm::set_block_size]= {0,};
 
    b1[0] = 1;
    b1[1] = 1;
    b2[0] = 1;

    bitblock_get_adapter bga1(b1);
    bitblock_get_adapter bga2(b2);
    bitblock_store_adapter bbsa(br);
    bm::bit_AND<bm::word_t> and_func;
    bit_recomb<bitblock_get_adapter,
               bitblock_get_adapter,
               bm::bit_AND<bm::word_t>,
               bitblock_store_adapter>
           (bga1, bga2,and_func, bbsa);
/*
    bit_recomb(bitblock_get_adapter(b1),
               bitblock_get_adapter(b2),
               bit_AND<bm::word_t>(),
               bitblock_store_adapter(br)
               );

    assert(br[0] == 1);
    for (int i = 1; i < bm::set_block_size; ++i)
    {
        assert(br[i] == 0);
    }

    bitblock_sum_adapter sa;
    bit_recomb(bitblock_get_adapter(b1),
               bitblock_get_adapter(b2),
               bit_COUNT_AND<bm::word_t>(),
               sa
               );
    assert(sa.sum() == 1);
*/
}

static
void BitForEachTest()
{
    cout << "---------------------------- BitForEachTest..." << endl;

    cout << "testing bit_list_4(), bitscan_popcnt().." << endl;
    {
        unsigned bit_list1[32];
        unsigned bit_list2[32];
        unsigned bit_list3[32];


        for (unsigned i = 0; i < 65536*50; ++i)
        {
            unsigned bits1 = bm::bit_list(i, bit_list1);
            unsigned bits2 = bm::bit_list_4(i, bit_list2);
            unsigned bits3 = bm::bitscan_popcnt(i, bit_list3);
            if (bits1 != bits2 || bits1 != bits3)
            {
                cout << "Bit for each test failed bit_cnt criteria!" << endl;
                exit(1);
            }
            for (unsigned j = 0; j < bits1; ++j)
            {
                if (bit_list1[j] != bit_list2[j] || bit_list1[j] != bit_list3[j])
                {
                    cout << "Bit for each check failed for " << j << endl;
                    exit(1);
                }
            }

        } // for
    }
    
    {
        cout << "testing bitscan_popcnt64()..." << endl;

        unsigned char bit_list[64];
        bm::id64_t w = 0;
        unsigned cnt;
        
        cnt = bm::bitscan_popcnt64(w, bit_list);
        if (cnt)
        {
            cout << "bitscan_popcnt64 cnt for 0x00 failed " << cnt << endl;
            exit(1);
        }
        
        w = ~w; // 0xFFFFF...
        
        cnt = bm::bitscan_popcnt64(w, bit_list);
        if (cnt != 64)
        {
            cout << "bitscan_popcnt64 cnt for 0xFFF failed " << cnt << endl;
            exit(1);
        }
        for (unsigned i = 0; i < cnt; ++i)
        {
            if (bit_list[i] != i)
            {
                cout << "bitscan_popcnt64 cnt at " << i << " != " << bit_list[i] << endl;
                exit(1);
            }
        } // for
        
        for (unsigned k = 63; k != 0; --k)
        {
            w <<= 1;
            cnt = bm::bitscan_popcnt64(w, bit_list);
            if (cnt != k)
            {
                cout << "bitscan_popcnt64 cnt for " << w << " cnt=" << cnt << endl;
                exit(1);
            }
            cout << "[" << cnt << "]:";
            for (unsigned i = 0; i < cnt; ++i)
            {
                cout << (unsigned)bit_list[i] << ", ";
            } // for
            cout << endl;
        } // for
    }
    

    cout << "---------------------------- BitForEachTest Ok." << endl;
}

static
void Log2Test()
{
    cout << "---------------------------- Log2 Test..." << endl;

    {
    unsigned l = bm::bit_scan_reverse32(~0u);
    cout << l << endl;
    assert(l == 31);
    l = bm::bit_scan_reverse(~0u);
    cout << l << endl;
    assert(l == 31);
    l = bm::bit_scan_reverse(~0ull);
    cout << l << endl;
    assert(l == 63);
    }

    {
    bm::id64_t v8 = 0x8000000000000000U;
    unsigned l = bm::bit_scan_reverse(v8);
    assert(l == 63);
    v8 = 0x4000000000000000U;
    l = bm::bit_scan_reverse(v8);
    assert(l == 62);
    }

    cout << "Stage 1" << endl;
    for (unsigned  i = 1; i <= 65535; ++i)
    {
        unsigned l1 = bm::ilog2<unsigned short>((unsigned short)i);
        unsigned l2 = iLog2(i);
        unsigned l3 = ilog2_LUT<unsigned short>((unsigned short)i);
        unsigned l4 = bm::bit_scan_reverse(i);
        if (l1 != l2 || l2 != l3 || l2 != l4)
        {
            cout << "Log2 error for " << i << endl;
            cout << l2 << " " << l3 << endl;;
            exit(1);
        }
    }
    cout << "Stage 2" << endl;
    for (unsigned  i = 1; i <= 10000*65535; ++i)
    {
        unsigned l1 = bm::ilog2<unsigned>(i);
        unsigned l2 = iLog2(i);
        unsigned l3 = ilog2_LUT<unsigned>(i);
        unsigned l4 = bm::bit_scan_reverse(i);
        if (l1 != l2 || l2 != l3 || l2 != l4)
        {
            cout << "Log2 error for " << i << endl;
            cout << l2 << " " << l3 << endl;;
            exit(1);
        }
    }
    cout << "Stage 3" << endl;
    unsigned v = 1;
    for (unsigned  i = 1; i <= 31; ++i)
    {
        v |= 1;
        unsigned l1 = bm::ilog2<unsigned>(v);
        unsigned l2 = iLog2(v);
        unsigned l3 = ilog2_LUT<unsigned>(v);
        unsigned l4 = bm::bit_scan_reverse(v);
        if (l1 != l2 || l2 != l3 || l2 != l4)
        {
            cout << "Log2 error for " << i << endl;
            cout << l2 << " " << l3 << endl;;
            exit(1);
        }
        
        bm::id64_t v8 = v;
        v8 <<= 32;
        unsigned l5 = bm::bit_scan_reverse(v8);
        if ((l4 + 32) != l5)
        {
            cout << "Log2 error for " << v8 << " " << v << endl;
            cout << i << " " <<" " << l4 << " " << (l4+32) << " " << l5 << endl;;
            exit(1);
        }
        
        v <<= 1;
    }
    cout << "---------------------------- Log2 Test Ok." << endl;
}


static
void LZCNTTest()
{
    cout << "---------------------------- LZCNT Test..." << endl;

    unsigned bsr;
    unsigned l = bm::count_leading_zeros(0);
    assert(l == 32);

    unsigned t = bm::count_trailing_zeros(0);
    assert(t == 32);

    l = bm::count_leading_zeros(2);
    unsigned bsf = bm::bit_scan_forward32(2);
    assert(bsf == 1);
    t = bm::count_trailing_zeros(2);
    assert(t == 1);

    l = bm::count_leading_zeros(~0u);
    assert(l == 0);
    l = bm::count_leading_zeros(~0u >> 1u);
    assert(l == 1);


    unsigned mask = ~0u;
    for (unsigned i = 1; i; i <<= 1)
    {
        l = bm::count_leading_zeros(i);
        t = bm::count_trailing_zeros(i);
        bsr = bm::bit_scan_reverse32(i);
        bsf = bm::bit_scan_forward32(i);
        assert(bsf == bsr);
        assert(l == 31 - bsf);
        assert(t == bsf);

        l = bm::count_leading_zeros(mask);
        bsf = bm::bit_scan_forward32(mask);
        bsr = bm::bit_scan_reverse32(mask);
        assert(l == 31 - bsr);
        mask >>= 1;
    }

    {
        bm::id64_t w = ~0ull;
        for (unsigned i = 0; i < 63; ++i)
        {
            unsigned lz = bm::count_leading_zeros_u64(w);
            assert(lz == i);
            w >>= 1;
        }
    }

    cout << "---------------------------- LZCNT Test..." << endl;
}

inline
unsigned proxy_bmi1_select64_lz(bm::id64_t val, unsigned rank)
{
#ifdef BMBMI1OPT
    return bmi1_select64_lz(val, rank);
#else
    return bm::word_select64_linear(val, rank);
#endif
}

inline
unsigned proxy_bmi1_select64_tz(bm::id64_t val, unsigned rank)
{
#ifdef BMBMI1OPT
    return bmi1_select64_tz(val, rank);
#else
    return bm::word_select64_linear(val, rank);
#endif
}


inline
unsigned proxy_bmi2_select64_pdep(bm::id64_t val, unsigned rank)
{
#ifdef BMBMI2OPT
    return bmi2_select64_pdep(val, rank);
#else
    return bm::word_select64_linear(val, rank);
#endif
}


// Returns the position of the rank'th 1.  (rank = 0 returns the 1st 1)
// Returns 64 if there are fewer than rank+1 1s.
/*
inline
unsigned select64_pdep_tzcnt(bm::id64_t val, unsigned rank) {
    uint64_t i = 1ull << rank;
    asm("pdep %[val], %[mask], %[val]"
            : [val] "+r" (val)
            : [mask] "r" (i));
    asm("tzcnt %[bit], %[index]"
            : [index] "=r" (i)
            : [bit] "g" (val)
            : "cc");
    return unsigned(i);
}
*/



static
void SelectTest()
{
    cout << "---------------------------- SELECT Test" << endl;
    
    {
        bm::id64_t w64 = 1;
        unsigned idx = bm::word_select64_linear(w64, 1);
        unsigned idx0 = word_select64_bitscan(w64, 1);
        unsigned idx1, idx4, idx3, idx5;
        assert(idx == 0);
        assert(idx0 == idx);
        idx4 = proxy_bmi1_select64_lz(w64, 1);
        assert(idx4 == idx);
        idx5 = proxy_bmi1_select64_tz(w64, 1);
        assert(idx5 == idx);
        
        
        idx3 = proxy_bmi2_select64_pdep(w64, 1);
        std::cerr << idx3 << " " << idx << endl;
        assert(idx3 == idx);

//        idx4 = word_select64_part(w64, 1);
//        assert(idx4 == idx);

/*
        w64 = 1 | (1 << 2) | (1 << 3) | (1 << 4) | (1 << 7);
        w64 = (1ull << 63) | 1;
        idx3 = avx2_select64(w64, 2);
        w64 = (1ull << 63) | 8 | 2;
        idx3 = avx2_select64(w64, 2);

        exit(1);

        w64 = ~0u;
        idx3 = avx2_select64(w64, 1);
        idx3 = avx2_select64(w64, 2);
        w64 = 16 | 2 | 1;
        idx3 = avx2_select64(w64, 1);
        idx3 = avx2_select64(w64, 2);
        idx3 = avx2_select64(w64, 3);
        w64 = w64 << 1;
        idx3 = avx2_select64(w64, 1);
        idx3 = avx2_select64(w64, 2);
        idx3 = avx2_select64(w64, 3);
        idx3 = avx2_select64(w64, 4);
        idx3 = avx2_select64(w64, 5);
        exit(1);
*/

        for (unsigned sel = 1; sel <= 64; ++sel)
        {
            idx = bm::word_select64_linear(~0ull, sel);
            assert(idx == sel-1);
            idx0 = word_select64_bitscan(~0ull, sel);
            assert(idx0 == idx);
            idx4 = proxy_bmi1_select64_lz(~0ull, sel);
            assert(idx4 == idx);
            idx5 = proxy_bmi1_select64_tz(~0ull, sel);
            assert(idx5 == idx);
            idx3 = proxy_bmi2_select64_pdep(~0ull, sel);
            assert(idx3 == idx);
        }
        
        for (idx = 0; w64; w64 <<= 1)
        {
            idx0 = bm::word_select64_linear(w64, 1);
            assert(idx0 == idx);
            idx1 = word_select64_bitscan(w64, 1);
            assert(idx1 == idx0);
            idx4 = proxy_bmi1_select64_lz(w64, 1);
            assert(idx4 == idx);
            idx5 = proxy_bmi1_select64_tz(w64, 1);
            assert(idx5 == idx);
            idx3 = proxy_bmi2_select64_pdep(w64, 1);
            assert(idx3 == idx);

            ++idx;
        }
    }

    {
        cout << "SELECT stress test." << std::endl;
        const unsigned test_size = 1000000 * 10;
        for (unsigned i = 1; i < test_size; ++i)
        {
            bm::id64_t w64 = i;
            bm::id64_t w64_1 = (w64 << 32) | w64;

            unsigned count = bm::word_bitcount64(w64);
            for (unsigned j = 1; j <= count; ++j)
            {
                unsigned idx0 = bm::word_select64_linear(w64, j);
                unsigned idx1 = word_select64_bitscan(w64, j);
                assert(idx0 == idx1);
                unsigned idx4 = proxy_bmi1_select64_lz(w64, j);
                assert(idx4 == idx1);
                unsigned idx5 = proxy_bmi1_select64_tz(w64, j);
                assert(idx5 == idx1);
                unsigned idx3 = proxy_bmi2_select64_pdep(w64, j);
                assert(idx3 == idx1);

            }
            
            count = bm::word_bitcount64(w64_1);
            for (unsigned j = 1; j <= count; ++j)
            {
                unsigned idx0 = bm::word_select64_linear(w64_1, j);
                unsigned idx1 = word_select64_bitscan(w64_1, j);
                assert(idx0 == idx1);
                unsigned idx4 = proxy_bmi1_select64_lz(w64_1, j);
                assert(idx4 == idx1);
                unsigned idx5 = proxy_bmi1_select64_tz(w64_1, j);
                assert(idx5 == idx1);
                unsigned idx3 = proxy_bmi2_select64_pdep(w64_1, j);
                assert(idx3 == idx1);

            }
            
            if (i % 1000000 == 0)
                cout << "\r" << i << std::flush;
        }
    }

    {
        cout << "SELECT bit-block test." << std::endl;
        unsigned cnt;
        
        BM_DECLARE_TEMP_BLOCK(tb1);
        for (unsigned i = 0; i < bm::set_block_size; ++i)
        {
            tb1.b_.w32[i] = ~0u;
        }
        for (unsigned i = 1; i <= 65535; ++i)
        {
            unsigned idx;
            unsigned rank = bm::bit_find_rank(tb1, i, 0, idx);
            assert(rank == 0);
            cnt = bm::bit_block_calc_count_to(tb1, i-1);
            assert(idx == cnt-1);

            for (unsigned j = 65535; j > i; --j)
            {
                cnt = bm::bit_block_calc_count_range(tb1, i, j);
                if (cnt)
                {
                    rank = bm::bit_find_rank(tb1, cnt, i, idx);
                    assert(rank == 0);
                    assert(idx == j);
                }
            }

            if (i % 128 == 0)
                cout << "\r" << i << std::flush;
        }
    }
    
    cout << "\n---------------------------- SELECT Test OK" << endl;
}


static
void BitEncoderTest()
{
    cout << "---------------------------- BitEncoderTest" << endl;
    
    unsigned char buf[1024 * 200] = {0, };
    
    {
        bm::encoder enc(buf, sizeof(buf));
        bm::bit_out<bm::encoder> bout(enc);
        
        unsigned value = 1024 + 3;
        bout.put_bits(value, 32);
        value = 1024 + 5;
        bout.put_bits(value, 32);
        
        bout.flush();
        
        bm::decoder dec(buf);
        bm::bit_in<bm::decoder> bin(dec);
        value = bin.get_bits(32);
        assert(value == 1024 + 3);
        value = bin.get_bits(32);
        assert(value == 1024 + 5);
    }
    
    {
        unsigned bits = 1;
        for (unsigned i = 1; i < (1u << 31u); i <<= 1, bits++)
        {
            bm::encoder enc(buf, sizeof(buf));
            bm::bit_out<bm::encoder> bout(enc);
            
            for (unsigned j = 0; j < 160; ++j)
            {
                bout.put_bits(i, bits);
            }
            bout.flush();
            
            bm::decoder dec(buf);
            bm::bit_in<bm::decoder> bin(dec);
            for (unsigned j = 0; j < 160; ++j)
            {
                unsigned value = bin.get_bits(bits);
                if (value != i)
                {
                    cerr << "Invalid encoding for i=" << i
                         << " value=" << value << " bits=" << bits << endl;
                    exit(1);
                }
            }
            
        } // for
    }
    
    {
        bm::encoder enc(buf, sizeof(buf));
        bm::bit_out<bm::encoder> bout(enc);
        
        for (unsigned i = 1; i < 65536; ++i)
        {
            unsigned bits = bm::bit_scan_reverse(i)+1;
            bout.put_bits(i, bits);
        } // for
        bout.flush();
        
        bm::decoder dec(buf);
        bm::bit_in<bm::decoder> bin(dec);
        for (unsigned i = 1; i < 65536; ++i)
        {
            unsigned bits = bm::bit_scan_reverse(i)+1;
            unsigned value = bin.get_bits(bits);
            if (value != i)
            {
                cerr << "2. Invalid encoding for i=" << i
                     << " value=" << value << " bits=" << bits << endl;
                exit(1);
            }

        } // for
        
    }

    // test for 24-48-bit encode
    {
        bm::encoder enc(buf, sizeof(buf));
        enc.put_24(0xFFFFFF);
        enc.put_24(0xFEAAFE);
        enc.put_48(0xF0FAFFBEEFFFUL);
        enc.put_8(0);

        assert(enc.size() == 6+6+1);

        bm::decoder dec(buf);
        auto v1 = dec.get_24();
        assert(v1 == 0xFFFFFF);
        auto v2 = dec.get_24();
        assert(v2 == 0xFEAAFE);
        auto v3 = dec.get_48();
        assert(v3 == 0xF0FAFFBEEFFFUL);

        unsigned char e = dec.get_8();
        assert(!e);
    }
    
    
    cout << "---------------------------- BitEncoderTest" << endl;
}

/// Random numbers test
template<typename V>
unsigned generate_inter_test_linear(V* arr, unsigned inc, unsigned target_size)
{
    V maskFF = (V)~V(0u);
    
    
    if (inc < 2 || inc > 65535)
        inc = 1;
    
    unsigned start = 1;
    unsigned sz = 0;
    while (sz < target_size)
    {
        arr[sz++] = V(start);
        if (inc + start >= maskFF)
        {
            arr[sz++] = maskFF;
            break;
        }
        start += inc;
        if (start < arr[sz-1])
            break;
        
    } // while
    return sz;
}


/// Random numbers test
template<typename V>
unsigned generate_inter_test(V* arr, unsigned inc_factor, unsigned target_size)
{
    V maskFF = (V)~V(0u);

    if (inc_factor < 2)
        inc_factor = 65535;
    
    unsigned start = rand() % 256;
    if (!start)
        start = 1;
    unsigned sz = 0;
    while (sz < target_size)
    {
        arr[sz++] = V(start);
        
        unsigned inc = unsigned(rand()) % inc_factor;
        if (!inc)
            inc = 1;
        start += inc;
        if (start >= maskFF)
            break;
    } // while

    for(unsigned i = 1; i < sz; ++i)
    {
        if (arr[i-1] >= arr[i])
        {
            return i;
        }
    }
    return sz;
}


static
void InterpolativeCodingTest()
{
    cout << "---------------------------- InterpolativeCodingTest() " << endl;
    
    unsigned char buf[1024 * 200] = {0, };
    unsigned char buf2[1024 * 200] = {0, };
    const bm::gap_word_t arr1[] = { 3, 4, 7, 13, 14, 15, 21, 25, 36, 38, 54, 62 };
    const bm::word_t arr2[] = { 30, 44, 78, 130, 140, 150, 210, 250, 3600, 3800, 540001, 620000258 };

    unsigned sz, sz2;
    {
        bm::encoder enc(buf, sizeof(buf));
        bm::bit_out<bm::encoder> bout(enc);
        
        sz = sizeof(arr1)/sizeof(arr1[0])-1;
        bout.bic_encode_u16_rg(arr1, sz, 0, 62);
        
        bout.flush();
    }
    {
        bm::encoder enc(buf2, sizeof(buf2));
        bm::bit_out<bm::encoder> bout(enc);
        
        sz2 = sizeof(arr2)/sizeof(arr2[0])-1;
        bout.bic_encode_u32_cm(arr2, sz2, 0, 620000258);
        bout.flush();
    }

    {
        decoder dec(buf);
        bm::bit_in<decoder> bin(dec);
        
        bm::gap_word_t arr2c[256] = {0, };
        bin.bic_decode_u16_rg(&arr2c[0], sz, 0, 62);
        for (unsigned i = 0; i < sz; ++i)
        {
            assert(arr1[i] == arr2c[i]);
        }
    }
    {
        decoder dec(buf2);
        bm::bit_in<decoder> bin(dec);
        
        bm::word_t arr2c[256] = {0, };
        bin.bic_decode_u32_cm(&arr2c[0], sz2, 0, 620000258);
        for (unsigned i = 0; i < sz2; ++i)
        {
            assert(arr2[i] == arr2c[i]);
        }
    }

    // -------------------------------------------------------
    cout << "BIC encode-decode U16 unit test" << endl;
    {
        bm::encoder enc(buf, sizeof(buf));
        bm::bit_out<bm::encoder> bout(enc);
        
        sz = sizeof(arr1)/sizeof(arr1[0])-1;
        bout.bic_encode_u16_cm(arr1, sz, 0, 62);
        
        bout.flush();
    }
    {
        decoder dec(buf);
        bm::bit_in<decoder> bin(dec);
        
        bm::gap_word_t arr2c[256] = {0, };
        bin.bic_decode_u16_cm(&arr2c[0], sz, 0, 62);
        for (unsigned i = 0; i < sz; ++i)
        {
            assert(arr1[i] == arr2c[i]);
        }
    }

    // -------------------------------------------------------

    cout << "\nu16 interpolated cm encoding Stress..." << endl;
    {
        const unsigned code_repeats = 1000000;
        const unsigned test_size = 12000;
        vector<gap_word_t> sa; sa.resize(test_size);
        vector<gap_word_t> da; da.resize(test_size);

        bm::gap_word_t* src_arr=&sa[0];
        bm::gap_word_t* dst_arr = &da[0];

        std::chrono::time_point<std::chrono::steady_clock> s;
        std::chrono::time_point<std::chrono::steady_clock> f;
        s = std::chrono::steady_clock::now();

        cout << "  linear pattern" << endl;
        for (unsigned k = 0; k < code_repeats; ++k)
        {
            unsigned inc = rand()%(65536*256);
            if (k == 0)
                inc = 1;
            sz = generate_inter_test_linear(src_arr, inc, test_size);
            assert(sz);
            assert(src_arr[0]);
            {
                bm::encoder enc(buf, sizeof(buf));
                bm::bit_out<bm::encoder> bout(enc);
                
                bout.bic_encode_u16_cm(src_arr, sz-1, 0, src_arr[sz-1]);
                bout.flush();
                auto ssz = enc.size();
                assert(ssz < sizeof(buf));
            }
            {
                decoder dec(buf);
                bm::bit_in<decoder> bin(dec);
                
                bin.bic_decode_u16_cm(&dst_arr[0], sz-1, 0, src_arr[sz-1]);
                dst_arr[sz-1]=src_arr[sz-1];
                for (unsigned i = 0; i < sz; ++i)
                {
                    assert(src_arr[i] == dst_arr[i]);
                    if (i)
                    {
                        assert(src_arr[i-1] < src_arr[i]);
                    }
                }
            }
            if ((k & 0xFFFF) == 0)
            {
                f = std::chrono::steady_clock::now();
                auto diff = f - s;
                auto d = std::chrono::duration <double, std::milli> (diff).count();

                cout << "\r" << k << "-" << code_repeats << " (" << d << "ms)" << flush;
                s = std::chrono::steady_clock::now();
            }
        }
        
        cout << "\n  random pattern" << endl;
        for (unsigned k = 0; k < code_repeats; ++k)
        {
            sz = generate_inter_test(src_arr, k, test_size);
            if (sz < 3)
                continue;

            assert(sz);
            assert(src_arr[0]);
            {
                bm::encoder enc(buf, sizeof(buf));
                bm::bit_out<bm::encoder> bout(enc);
                
                bout.bic_encode_u16_cm(src_arr, sz, 0, src_arr[sz-1]);
                bout.flush();
            }
            {
                decoder dec(buf);
                bm::bit_in<decoder> bin(dec);
                
                bin.bic_decode_u16_cm(&dst_arr[0], sz, 0, src_arr[sz-1]);
                //dst_arr[sz-1]=src_arr[sz-1];
                for (unsigned i = 0; i < sz; ++i)
                {
                    if (i)
                    {
                        assert(src_arr[i-1] < src_arr[i]);
                    }
                    assert(src_arr[i] == dst_arr[i]);
                }
            }
            if ((k & 0xFFF) == 0)
            {
                f = std::chrono::steady_clock::now();
                auto diff = f - s;
                auto d = std::chrono::duration <double, std::milli> (diff).count();

                cout << "\r" << k << "-" << code_repeats << " (" << d << "ms)" << flush;
                s = std::chrono::steady_clock::now();
            }
        } // for k

    }

    // -------------------------------------------------------

    cout << "\nu32 interpolated cm encoding Stress..." << endl;
    {
        const unsigned code_repeats = 1000000;
        const unsigned test_size = 12000;
        vector<unsigned> sa; sa.resize(test_size);
        vector<unsigned> da; da.resize(test_size);

        bm::word_t* src_arr=&sa[0];
        bm::word_t* dst_arr = &da[0];

        std::chrono::time_point<std::chrono::steady_clock> s;
        std::chrono::time_point<std::chrono::steady_clock> f;
        s = std::chrono::steady_clock::now();

        cout << "  linear pattern" << endl;
        for (unsigned k = 0; k < code_repeats; ++k)
        {
            unsigned inc = rand()%(65536*256);
            if (k == 0)
                inc = 1;
            sz = generate_inter_test_linear(src_arr, inc, test_size);
            assert(sz);
            assert(src_arr[0]);
            {
                bm::encoder enc(buf, sizeof(buf));
                bm::bit_out<bm::encoder> bout(enc);
                
                bout.bic_encode_u32_cm(src_arr, sz-1, 0, src_arr[sz-1]);
                bout.flush();
                auto ssz = enc.size();
                assert(ssz < sizeof(buf));
            }
            {
                decoder dec(buf);
                bm::bit_in<decoder> bin(dec);
                
                bin.bic_decode_u32_cm(&dst_arr[0], sz-1, 0, src_arr[sz-1]);
                dst_arr[sz-1]=src_arr[sz-1];
                for (unsigned i = 0; i < sz; ++i)
                {
                    assert(src_arr[i] == dst_arr[i]);
                    if (i)
                    {
                        assert(src_arr[i-1] < src_arr[i]);
                    }
                }
            }
            if ((k & 0xFFFF) == 0)
            {
                f = std::chrono::steady_clock::now();
                auto diff = f - s;
                auto d = std::chrono::duration <double, std::milli> (diff).count();

                cout << "\r" << k << "-" << code_repeats << " (" << d << "ms)" << flush;
                s = std::chrono::steady_clock::now();
            }
        }

        cout << "\n  random pattern" << endl;
        for (unsigned k = 0; k < code_repeats; ++k)
        {
            sz = generate_inter_test(src_arr, k, test_size);
            if (sz < 3)
                continue;

            assert(sz);
            assert(src_arr[0]);
            {
                bm::encoder enc(buf, sizeof(buf));
                bm::bit_out<bm::encoder> bout(enc);
                
                bout.bic_encode_u32_cm(src_arr, sz, 0, src_arr[sz-1]);
                bout.flush();
            }
            {
                decoder dec(buf);
                bm::bit_in<decoder> bin(dec);
                
                bin.bic_decode_u32_cm(&dst_arr[0], sz, 0, src_arr[sz-1]);
                //dst_arr[sz-1]=src_arr[sz-1];
                for (unsigned i = 0; i < sz; ++i)
                {
                    if (i)
                    {
                        assert(src_arr[i-1] < src_arr[i]);
                    }
                    assert(src_arr[i] == dst_arr[i]);
                }
            }
            if ((k & 0xFFF) == 0)
            {
                f = std::chrono::steady_clock::now();
                auto diff = f - s;
                auto d = std::chrono::duration <double, std::milli> (diff).count();

                cout << "\r" << k << "-" << code_repeats << " (" << d << "ms)" << flush;
                s = std::chrono::steady_clock::now();
            }
        } // for k

    }

    cout << "\nu16 interpolated encoding Stress..." << endl;
    {
        const unsigned code_repeats = 1000000;
        const unsigned test_size = 65536;
        vector<bm::gap_word_t> sa; sa.resize(test_size);
        vector<bm::gap_word_t> da; da.resize(test_size);

        bm::gap_word_t* src_arr= &sa[0];
        bm::gap_word_t* dst_arr= &da[0];

        std::chrono::time_point<std::chrono::steady_clock> s;
        std::chrono::time_point<std::chrono::steady_clock> f;
        s = std::chrono::steady_clock::now();

        cout << "  linear pattern" << endl;
//        unsigned inc = rand()%128;
//        sz = generate_inter_test_linear(src_arr, inc);
        for (unsigned k = 0; k < code_repeats; ++k)
        {
            unsigned inc = rand()%128;
            if (k == 0)
                inc = 1;

            sz = generate_inter_test_linear(src_arr, inc, 65536);
            assert(sz);
            assert(src_arr[0]);
            if(src_arr[sz-1]<65535)
               src_arr[sz-1]=65535;
            {
                bm::encoder enc(buf, sizeof(buf));
                bm::bit_out<bm::encoder> bout(enc);
                
                bout.bic_encode_u16_rg(src_arr, sz, 0, 65535);
                bout.flush();
            }
            {
                decoder dec(buf);
                bm::bit_in<decoder> bin(dec);
                
                bin.bic_decode_u16_rg(&dst_arr[0], sz, 0, 65535);
                //dst_arr[sz-1]=65535;
                for (unsigned i = 0; i < sz; ++i)
                {
                    assert(src_arr[i] == dst_arr[i]);
                }
            }
            if ((k & 0xFFFF) == 0)
            {
                f = std::chrono::steady_clock::now();
                auto diff = f - s;
                auto d = std::chrono::duration <double, std::milli> (diff).count();

                cout << "\r" << k << "-" << code_repeats << " (" << d << "ms)" << flush;
                s = std::chrono::steady_clock::now();
            }
        }
        cout << "  random pattern" << endl;

        for (unsigned k = 0; k < code_repeats; ++k)
        {
            sz = generate_inter_test(src_arr, k, 65536);
            if (sz < 3)
                continue;
            assert(sz);
            assert(src_arr[0]);
            assert(src_arr[sz-1]<=65535);
            {
                bm::encoder enc(buf, sizeof(buf));
                bm::bit_out<bm::encoder> bout(enc);
                
                bout.bic_encode_u16_rg(src_arr, sz, 0, 65535);
                bout.flush();
            }
            {
                decoder dec(buf);
                bm::bit_in<decoder> bin(dec);
                
                bin.bic_decode_u16_rg(&dst_arr[0], sz, 0, 65535);
                //dst_arr[sz-1]=65535;
                for (unsigned i = 0; i < sz; ++i)
                {
                    assert(src_arr[i] == dst_arr[i]);
                }
            }
            if ((k & 0xFFFF) == 0)
                cout << "\r" << k << "-" << code_repeats << flush;
        } // for k
        
    }


    cout << "---------------------------- InterpolativeCodingTest() OK " << endl;
}

static
void GammaEncoderTest()
{
    cout << "---------------------------- GammaEncoderTest" << endl;
    
    
    unsigned char buf1[2048 * 4] = {0, };
    
    cout << "Stage 1" << endl;

    {
    encoder enc(buf1, sizeof(buf1));
    typedef bit_out<encoder>  TBitIO;
    bit_out<encoder> bout(enc);
    gamma_encoder<bm::gap_word_t, TBitIO> gamma(bout);     
    gamma(65534);
    }

    {
    decoder dec(buf1);
    typedef bit_in<decoder> TBitIO;
    bit_in<decoder> bin(dec);
    gamma_decoder<bm::gap_word_t, TBitIO> gamma(bin);
    
    gap_word_t value = gamma();
    if (value != 65534)
        {
            cout << "Gamma decoder error for value=" << value << endl;
            exit(1);
        }             
    }


    {
    encoder enc(buf1, sizeof(buf1));
    typedef bit_out<encoder>  TBitIO;
    bit_out<encoder> bout(enc);
    gamma_encoder<bm::gap_word_t, TBitIO> gamma(bout);
     
    for (gap_word_t i = 1; i < 15; ++i)
    {
        gamma(i);
    } 
    }    
    
    {
    decoder dec(buf1);
    typedef bit_in<decoder> TBitIO;
    bit_in<decoder> bin(dec);
    gamma_decoder<bm::gap_word_t, TBitIO> gamma(bin);
    
    for (gap_word_t i = 1; i < 15; ++i)
    {
        gap_word_t value = gamma();
        if (value != i)
        {
            cout << "Gamma decoder error for " << i << " value=" << value << endl;
            exit(1);
        }
    }     
    
    }

    cout << "Stage 2" << endl;

    for (unsigned i = 0; i < 256; ++i)
    {
        gap_word_t short_block[64] = {0,};
        
        {
        encoder enc(buf1, sizeof(buf1));
        typedef bit_out<encoder>  TBitIO;
        bit_out<encoder> bout(enc);
        gamma_encoder<bm::gap_word_t, TBitIO> gamma(bout);
         

        for (unsigned j = 0; j < 64; ++j)
        {
            gap_word_t a = gap_word_t(rand() % 65535);
            if (!a) a = 65535; // 0 is illegal
            gap_word_t value = short_block[j] = a;
            gamma(value);
        } // for
        }

        {
        decoder dec(buf1);
        typedef bit_in<decoder> TBitIO;
        bit_in<decoder> bin(dec);
        gamma_decoder<bm::gap_word_t, TBitIO> gamma(bin);
        
        for (unsigned j = 0; j < 64; ++j)
        {
            gap_word_t value = short_block[j];
            gap_word_t r = gamma();
            if (r != value)
            {
                cout << "Gamma encoding failure for value=" << value << " gamma=" << r << endl;
                exit(1);
            }
        } // for
        }
    }


    cout << "Stage 3" << endl;

    unsigned code_value = 65535;
    for (unsigned i = 0; i < 10000; ++i)
    {
        gap_word_t short_block[1000] = {0,};
        
        {
        encoder enc(buf1, sizeof(buf1));
        typedef bit_out<encoder>  TBitIO;
        bit_out<encoder> bout(enc);
        gamma_encoder<bm::gap_word_t, TBitIO> gamma(bout);
         
        for (unsigned j = 0; j < 1000; ++j)
        {
            gap_word_t a = (unsigned short)code_value;
            if (!a) 
            {
                code_value = a = 65535;
            }

            gap_word_t value = short_block[j] = a;
            gamma(value);
            --code_value;
        } // for
        }

        {
        decoder dec(buf1);
        typedef bit_in<decoder> TBitIO;
        bit_in<decoder> bin(dec);
        gamma_decoder<bm::gap_word_t, TBitIO> gamma(bin);
        
        for (unsigned j = 0; j < 1000; ++j)
        {
            gap_word_t value = short_block[j];
            gap_word_t r = gamma();
            if (r != value)
            {
                cout << "Gamma encoding failure for value=" << value << " gamma=" << r << endl;
                exit(1);
            }
        } // for
        }
    }


    cout << "---------------------------- GammaEncoderTest Ok." << endl;

}

template<class SV, class Vect>
bool CompareSparseVector(const SV& sv, const Vect& vect,
                         bool interval_filled = false,
                         bool detailed = true)
{

    if (vect.size() != sv.size())
    {
        cerr << "Sparse vector size test failed!" << vect.size() << "!=" << sv.size() << endl;
        return false;
    }
    
    if (sv.is_nullable())
    {
        const typename SV::bvector_type* bv_null = sv.get_null_bvector();
        assert(bv_null);
        unsigned non_null_cnt = bv_null->count();
        if (vect.size() != non_null_cnt)
        {
            if (!interval_filled)
            {
                cerr << "NULL vector count failed." << non_null_cnt << " size=" << vect.size() << endl;
                exit(1);
            }
        }
    }

    if (detailed)
    {
        typename SV::const_iterator it = sv.begin();
        typename SV::const_iterator it_end = sv.end();

        for (unsigned i = 0; i < vect.size(); ++i)
        {
            typename Vect::value_type v1 = vect[i];
            typename SV::value_type v2 = sv[i];
            typename SV::value_type v3 = *it;

            int cmp = sv.compare(i, v1);
            assert(cmp == 0);
            if (v1 > 0)
            {
                cmp = sv.compare(i, v1-1);
                assert(cmp > 0);
            }

            if (v1 != v2)
            {
                cerr << "SV discrepancy:" << "sv[" << i << "]=" << v2
                     <<  " vect[" << i << "]=" << v1
                     << endl;
                assert(0);return false;
            }
            if (v1 != v3)
            {
                cerr << "SV discrepancy:" << "sv[" << i << "]=" << v2
                     <<  " *it" << v3
                     << endl;
                return false;
            }
            assert(it < it_end);
            ++it;
        } // for
        if (it != it_end)
        {
            cerr << "sv const_iterator discrepancy!" << endl;
            return false;
        }
    }
    
    // extraction comparison
    if (detailed)
    {
        std::vector<typename SV::value_type> v1(sv.size());
        std::vector<typename SV::value_type> v1r(sv.size());
        sv.extract(&v1[0], sv.size(), 0);
        sv.extract_range(&v1r[0], sv.size(), 0);
        for (unsigned i = 0; i < sv.size(); ++i)
        {
            if (v1r[i] != v1[i] || v1[i] != vect[i])
            {
                cerr << "TestEqualSparseVectors Extract 1 failed at:" << i
                     << " v1[i]=" << v1[i] << " v1r[i]=" << v1r[i]
                     << endl;
                exit(1);
            }
        } // for
    }

    // serialization comparison
    BM_DECLARE_TEMP_BLOCK(tb)
    sparse_vector_serial_layout<SV> sv_lay;
    bm::sparse_vector_serialize<SV>(sv, sv_lay, tb);
    SV sv2;
    const unsigned char* buf = sv_lay.buf();
    int res = bm::sparse_vector_deserialize(sv2, buf, tb);
    if (res != 0)
    {
        cerr << "De-Serialization error" << endl;
        exit(1);
    }
    if (sv.is_nullable() != sv2.is_nullable())
    {
        cerr << "Serialization comparison of two svectors failed (NULL vector)" << endl;
        exit(1);
    }
    const typename SV::bvector_type* bv_null = sv.get_null_bvector();
    const typename SV::bvector_type* bv_null2 = sv.get_null_bvector();
    
    if (bv_null != bv_null2 && (bv_null == 0 || bv_null2 == 0))
    {
        cerr << "Serialization comparison (NUUL vector missing)!" << endl;
        exit(1);
    }
    if (bv_null)
    {
        if (bv_null->compare(*bv_null2) != 0)
        {
            cerr << "Serialization comparison of two svectors (NUUL vectors unmatch)!" << endl;
            exit(1);
        }
    }

    if (!sv.equal(sv2) )
    {
        cerr << "Error: Serialization comparison of two svectors failed!" << endl;
        typename SV::size_type pos;
        bool b = bm::sparse_vector_find_first_mismatch(sv, sv2, pos);
        assert(b);
        cerr << "Mismatch at: " << pos << endl;

        sparse_vector_serial_layout<SV> sv_lay1;
        bm::sparse_vector_serialize<SV>(sv, sv_lay1, tb);
        SV sv3;
        bm::sparse_vector_deserialize(sv3, buf, tb);


        exit(1);
    }
    
    return true;
}

template<class SV>
bool TestEqualSparseVectors(const SV& sv1, const SV& sv2, bool detailed = true)
{
    if (sv1.size() != sv2.size())
    {
        cerr << "TestEqualSparseVectors failed incorrect size" << endl;
        exit(1);
    }
    
    if (sv1.is_nullable() == sv2.is_nullable())
    {
        bool b = sv1.equal(sv2);
        if (!b)
        {
            cerr << "sv1.equal(sv2) failed" << endl;
            return b;
        }
        const typename SV::bvector_type* bv_null1 = sv1.get_null_bvector();
        const typename SV::bvector_type* bv_null2 = sv2.get_null_bvector();
        
        if (bv_null1 != bv_null2)
        {
            int r = bv_null1->compare(*bv_null2);
            if (r != 0)
            {
                cerr << "sparse NULL-vectors comparison failed" << endl;
                exit(1);
            }
        }
    }
    else  // NULLable does not match
    {
        detailed = true; // simple check not possible, use slow, detailed
    }
    

    // test non-offset extraction
    //
    {
        std::vector<unsigned> v1(sv1.size());
        std::vector<unsigned> v1r(sv1.size());
        std::vector<unsigned> v1p(sv1.size());
        
        sv1.extract(&v1[0], sv1.size(), 0);
        sv1.extract_range(&v1r[0], sv1.size(), 0);
        sv1.extract_plains(&v1p[0], sv1.size(), 0);
        
        for (unsigned i = 0; i < sv1.size(); ++i)
        {
            if (v1r[i] != v1[i] || v1p[i] != v1[i])
            {
                cerr << "TestEqualSparseVectors Extract 1 failed at:" << i
                     << " v1[i]=" << v1[i] << " v1r[i]=" << v1r[i] << " v1p[i]=" << v1p[i]
                     << endl;
                exit(1);
            }
        } // for
    }

    // test offset extraction
    //
    {
        std::vector<unsigned> v1(sv1.size());
        std::vector<unsigned> v1r(sv1.size());
        std::vector<unsigned> v1p(sv1.size());
        
        unsigned pos = sv1.size() / 2;
        
        sv1.extract(&v1[0], sv1.size(), pos);
        sv1.extract_range(&v1r[0], sv1.size(), pos);
        sv1.extract_plains(&v1p[0], sv1.size(), pos);
        
        for (unsigned i = 0; i < sv1.size(); ++i)
        {
            if (v1r[i] != v1[i] || v1p[i] != v1[i])
            {
                cerr << "TestEqualSparseVectors Extract 1 failed at:" << i
                     << " v1[i]=" << v1[i] << " v1r[i]=" << v1r[i] << " v1p[i]=" << v1p[i]
                     << endl;
                exit(1);
            }
        } // for
    }

    {
        SV svv1(sv1);
        SV svv2(sv2);
        
        bm::null_support is_null = (sv1.is_nullable() == sv2.is_nullable()) ? bm::use_null : bm::no_null;
        
        bool b = svv1.equal(svv2, is_null);
        if (!b)
        {
            cerr << "Equal, copyctor comparison failed" << endl;
            return b;
        }

        svv1.swap(svv2);
        b = svv1.equal(svv2, is_null);
        if (!b)
        {
            cerr << "Equal, copyctor-swap comparison failed" << endl;
            return b;
        }
    }

    // comparison using elements assignment via reference
    if (detailed)
    {
        SV sv3;
        sv3.resize(sv1.size());
        for (unsigned i = 0; i < sv1.size(); ++i)
        {
            sv3[i] = sv1[i];
            unsigned v1 = sv1[i];
            unsigned v2 = sv3[i];
            if (v1 != v2)
            {
                cerr << "1. sparse_vector reference assignment validation failed" << endl;
                return false;
            }
        }
        bm::null_support is_null = (sv1.is_nullable() == sv3.is_nullable()) ? bm::use_null : bm::no_null;
        bool b = sv1.equal(sv3, is_null);
        if (!b)
        {
            cerr << "2. sparse_vector reference assignment validation failed" << endl;
            return b;
        }
    }
    
    // comparison via const_iterators
    //
    {{
        typename SV::const_iterator it1 = sv1.begin();
        typename SV::const_iterator it2 = sv2.begin();
        typename SV::const_iterator it1_end = sv1.end();
        
        for (; it1 < it1_end; ++it1, ++it2)
        {
            if (*it1 != *it2)
            {
                cerr << "1. sparse_vector::const_iterator validation failed" << endl;
                return false;
            }
        }
    }}

    // comparison through serialization
    //
    {{
        int res;
        bm::sparse_vector_serial_layout<SV> sv_lay;
        bm::sparse_vector_serialize(sv1, sv_lay);
        
        // copy buffer to check if serialization size is actually correct
        const unsigned char* buf = sv_lay.buf();
        size_t buf_size = sv_lay.size();
        
        vector<unsigned char> tmp_buf(buf_size);
        ::memcpy(&tmp_buf[0], buf, buf_size);
        
        SV sv3;
        res = bm::sparse_vector_deserialize(sv3, &tmp_buf[0]);
        if (res != 0)
        {
            cerr << "De-Serialization error in TestEqualSparseVectors()" << endl;
            exit(1);
        }
        
        const typename SV::bvector_type* bv_null1 = sv1.get_null_bvector();
        const typename SV::bvector_type* bv_null2 = sv2.get_null_bvector();
        const typename SV::bvector_type* bv_null3 = sv3.get_null_bvector();
        
        if (bv_null1 && bv_null3)
        {
            int r = bv_null1->compare(*bv_null3);
            if (r != 0)
            {
                cerr << "2. NULL bvectors comparison failed" << endl;
                exit(1);
            }
        }
        if (bv_null1 && bv_null2)
        {
            int r = bv_null1->compare(*bv_null2);
            if (r != 0)
            {
                cerr << "3. NULL bvectors comparison failed" << endl;
                exit(1);
            }
        }

        bm::null_support is_null = (sv1.is_nullable() == sv3.is_nullable()) ? bm::use_null : bm::no_null;
        if (!sv1.equal(sv3, is_null) )
        {
            cerr << "Serialization comparison of two svectors failed (1)" << endl;
            exit(1);
        }
        is_null = (sv2.is_nullable() == sv3.is_nullable()) ? bm::use_null : bm::no_null;
        if (!sv2.equal(sv3, is_null) )
        {
            cerr << "Serialization comparison of two svectors failed (2)" << endl;
            exit(1);
        }
        
    
    }}
    return true;
}

static
void TestBasicMatrix()
{
    cout << "---------------------------- Basic bit-matrix test" << endl;
    
    // construction-destruction
    {
        bm::basic_bmatrix<bvect> bmtr(10);
        bvect* bv = bmtr.construct_row(0);
        assert(bv);
        
        bv->set(10);
        
        // copy content
        //
        bm::basic_bmatrix<bvect> bmtr2(bmtr);
        bvect* bv2 = bmtr2.construct_row(0);
        assert(bv2);
        
        bool b = bv2->test(10);
        assert(b);


        bm::basic_bmatrix<bvect> bmtr4(10);
        {
            bm::basic_bmatrix<bvect> bmtr3(10);
            bmtr3.construct_row(0)->set(110);
            bmtr4 = bmtr3;
            bv = bmtr4.construct_row(0);
            b = bv->test(10);
            assert(!b);
            b = bv->test(110);
            assert(b);
        }
        bmtr4 = bmtr2;
        bv = bmtr4.construct_row(0);
        b = bv->test(10);
        assert(b);
        b = bv->test(110);
        assert(!b);
        
        bm::basic_bmatrix<bvect> bmtr5(2);
        bmtr5 = bmtr2;
        assert(bmtr5.rows()==10);
        
        bm::basic_bmatrix<bvect> bmtr6(11);
        bmtr6.construct_row(0)->set(210);
        bmtr6.swap(bmtr5);
        assert(bmtr6.rows()==10);
        assert(bmtr5.rows()==11);

        bv = bmtr6.construct_row(0);
        b = bv->test(210);
        assert(!b);
        bv = bmtr5.construct_row(0);
        b = bv->test(210);
        assert(b);
    }
    
    
    // octet assignment logic
    {
        bm::basic_bmatrix<bvect> bmtr(32);
        bmtr.set_octet(0, 0, '3');
        bmtr.set_octet(1, 0, 1);
        unsigned char ch;
        ch = bmtr.get_octet(0, 0);
        assert(ch == '3');
        ch = bmtr.get_octet(1, 0);
        assert(ch == 1);
        
        bmtr.optimize();
        
        ch = bmtr.get_octet(0, 0);
        assert(ch == '3');
        ch = bmtr.get_octet(1, 0);
        assert(ch == 1);

    }
    {
        bm::basic_bmatrix<bvect> bmtr(32);
        bmtr.set_octet(0, 0, 1);
        bmtr.set_octet(0, 1, 2);
        bmtr.set_octet(0, 2, 'G');
        bmtr.set_octet(0, 3, 'C');
        unsigned char ch;
        ch = bmtr.get_octet(0, 0);
        assert(ch == 1);
        ch = bmtr.get_octet(0, 1);
        assert(ch == 2);
        ch = bmtr.get_octet(0, 2);
        assert(ch == 'G');
        ch = bmtr.get_octet(0, 3);
        assert(ch == 'C');
        ch = bmtr.get_octet(0, 0);
        assert(ch == 1);
        
        bmtr.optimize();
        ch = bmtr.get_octet(0, 0);
        assert(ch == 1);
        ch = bmtr.get_octet(0, 1);
        assert(ch == 2);
        ch = bmtr.get_octet(0, 2);
        assert(ch == 'G');
        ch = bmtr.get_octet(0, 3);
        assert(ch == 'C');
    }
    
    
    cout << "---------------------------- Basic bit-matrix test OK" << endl;
}



static
void TestSparseVector()
{
    cout << "---------------------------- Bit-plain sparse vector test" << endl;
    BM_DECLARE_TEMP_BLOCK(tb)


    typedef bm::sparse_vector<unsigned,bvect> svector;
    typedef bm::sparse_vector<unsigned long long, bvect> svector64;

    // basic construction (NULL-able vector)
    {{
        bm::sparse_vector<unsigned, bvect> sv1;
        bool n = sv1.is_nullable();
        assert(!n);
        const bvect* bvp = sv1.get_null_bvector();
        assert(bvp==0);
        
        bm::sparse_vector<unsigned, bvect> sv2(bm::use_null);
        n = sv2.is_nullable();
        assert(n);
        
        sv1 = sv2;
        assert(sv1.is_nullable());
        
        bm::sparse_vector<unsigned, bvect> sv3(sv1);
        assert(sv3.is_nullable());
        
        bm::sparse_vector<unsigned, bvect> sv4;
        sv3.swap(sv4);
        assert(sv4.is_nullable());
        assert(!sv3.is_nullable());
        bvp = sv4.get_null_bvector();
        assert(bvp);
    }}

    // reference vector construction for XOR serialization
    {{
        bm::sparse_vector<unsigned, bvect> sv;
        sv.push_back(1);
        sv.push_back(8);

        bm::bv_ref_vector<bvect> ref_vect;
        ref_vect.build(sv.get_bmatrix());

        assert(ref_vect.size() == 2);

        auto idx = ref_vect.find(0);
        assert(idx == 0);
        idx =  ref_vect.find(1);
        assert(idx == ref_vect.not_found());
        idx =  ref_vect.find(2);
        assert(idx == ref_vect.not_found());

        // test for 8 which is 1 << 3
        idx =  ref_vect.find(3);
        assert(idx == 1);

        // test XOR scanner
        //

        bm::xor_scanner<bvect> xscan;
        bm::xor_scanner<bvect>::bv_ref_vector_type r_vect;
        xscan.set_ref_vector(&r_vect);
        r_vect.build(sv.get_bmatrix());

        const bvect* bv_x = sv.get_plain(0);
        const bvect::blocks_manager_type& bman_x = bv_x->get_blocks_manager();
        const bm::word_t* block_x = bman_x.get_block_ptr(0, 0);

        xscan.compute_x_block_stats(block_x);
        assert(xscan.get_x_bc() == 1);
        assert(xscan.get_x_gc() == 2);
        assert(xscan.get_x_block_best() == 1);

        idx = xscan.get_ref_vector().find(0);
        assert(idx == 0);

        bool f = xscan.search_best_xor_mask(block_x,
                                            1, xscan.get_ref_vector().size(),
                                            0, 0, tb);
        assert(!f);
    }}

    // XOR scanner EQ test
    {{
        bm::sparse_vector<unsigned, bvect> sv;
        sv.push_back(9);
        sv.push_back(9);

        bm::xor_scanner<bvect> xscan;
        bm::xor_scanner<bvect>::bv_ref_vector_type r_vect;
        r_vect.build(sv.get_bmatrix());
        xscan.set_ref_vector(&r_vect);

        const bvect* bv_x = sv.get_plain(0);
        const bvect::blocks_manager_type& bman_x = bv_x->get_blocks_manager();
        const bm::word_t* block_x = bman_x.get_block_ptr(0, 0);

        xscan.compute_x_block_stats(block_x);
        assert(xscan.get_x_bc() == 2);
        assert(xscan.get_x_gc() == 2);
        assert(xscan.get_x_block_best() == 2);

        auto idx = xscan.get_ref_vector().find(0);
        assert(idx == 0);

        bool f = xscan.search_best_xor_mask(block_x,
                                            1, xscan.get_ref_vector().size(),
                                            0, 0, tb);
        assert(f);
        idx = xscan.found_ridx();
        assert(idx == 1);
        assert(xscan.get_x_best_metric() == 0); // EQ
        assert(xscan.is_eq_found());
        idx = xscan.get_ref_vector().get_row_idx(idx);
        assert(idx == 3); // matrix row 3
    }}
    
    // basic const_iterator construction
    {{
        bm::sparse_vector<unsigned, bvect> sv1;
        svector::const_iterator it_end;
        svector::const_iterator it = sv1.begin();
        
        assert(!it.valid());
        assert(!it_end.valid());
        assert(it != it_end);
        it.invalidate();
        assert(!it.valid());
        
        it.go_to(1);
        assert(!it.valid());
        
        svector::const_iterator it_end2 = sv1.end();
        assert(!it_end2.valid());
    }}
    
    // test empty vector serialization
    {{
        int res;

        bm::sparse_vector<unsigned, bvect> sv1;
        bm::sparse_vector<unsigned, bvect> sv2;
        bm::sparse_vector_serial_layout<svector> sv_layout;
        bm::sparse_vector_serialize(sv1, sv_layout);

        const unsigned char* buf = sv_layout.buf();
        res = bm::sparse_vector_deserialize(sv2, buf);
        if (res != 0)
        {
            cerr << "De-Serialization error" << endl;
            exit(1);
        }
        if (!sv1.equal(sv2) )
        {
            cerr << "Serialization comparison of empty vectors failed" << endl;
            exit(1);
        }
    }}
    
    // test move construction
    {{
        vector<svector> v;
        v.push_back(svector());
        v.push_back(svector());
        v[0] = svector();
    }}
    
    // test NULL operations
    {{
        bm::sparse_vector<unsigned, bvect> sv1;
        bm::sparse_vector<unsigned, bvect> sv2(bm::use_null);
        sv1.resize(10);
        sv2.resize(10);
        for (unsigned i = 0; i < sv1.size(); ++i)
        {
            assert(!sv1.is_null(i));
            assert(sv2.is_null(i));
        }
        unsigned arr[3] = {1,2,3};
        sv1.import(arr, 3);
        sv2.import(arr, 3);
        assert(!sv1.is_null(0));
        assert(!sv2.is_null(0));
        assert(!sv2.is_null(1));
        assert(!sv2.is_null(2));
        assert(sv2.is_null(3));
        
        assert(sv2.is_null(5));
        sv2.set(5, 123);
        assert(!sv2.is_null(5));
        
        sv2.set_null(5);
        assert(sv2.is_null(5));
        assert(sv2[5].is_null());

        
        bm::sparse_vector<unsigned, bvect> sv3(sv2);
        assert(sv3.is_nullable());
        
        assert(!sv3.is_null(0));
        assert(!sv3.is_null(1));
        assert(!sv3.is_null(2));
        
        sv3.clear_range(0, 1, true);
        assert(sv3.is_null(0));
        assert(sv3.is_null(1));
        assert(!sv3.is_null(2));
        
        
        sv3 = sv1;
        assert(sv3.is_nullable());
        
        sv1.clear();
        assert(!sv1.is_nullable());
        sv2.clear();
        assert(sv2.is_nullable());
    }}
    
    {{
    bm::sparse_vector<unsigned, bvect> sv;
    unsigned arr[3] = {1,2,3};
    sv.import(arr, 3);
    cout << "sv.size() = " << sv.size() << endl;
    cout << "sv[]:";
    for (unsigned i = 0; i < sv.size(); ++i)
    {
        cout << sv.at(i) << ",";
    }
    cout << endl;

    bm::sparse_vector_scanner<bm::sparse_vector<unsigned, bvect> > scanner;
    bvect bv;
    scanner.find_nonzero(sv, bv);
    if (bv.count() != sv.size())
    {
        cerr << "compute_nonzero_bvector test failed" << endl;
        exit(1);
    }
    }}
    
    {{
        bm::sparse_vector<unsigned, bvect> sv;
        sv.push_back(1);
        sv.push_back(1);
        sv.push_back(1);
        unsigned arr[1024];
        
        unsigned esize =  sv.extract(&arr[0], 1024, 0);
        assert(esize == 3);
        assert(arr[0] == 1);
        assert(arr[1] == 1);
        assert(arr[2] == 1);
    }}
    cout << "sv::push_back_null()" << endl;
    {{
        bm::sparse_vector<unsigned, bvect> sv(bm::use_null);
        sv.push_back_null(10);
        auto sz = sv.size();
        assert(sz==10);
        sv.push_back(1);
        sz = sv.size();
        assert(sz==11);
        assert(sv[10] == 1);
        for (bvect::size_type i = 0; i < 10; ++i)
        {
            auto v = sv[i];
            assert(v == 0);
        }
        sv.optimize();
        assert(sv[10] == 1);
        for (bvect::size_type i = 0; i < 10; ++i)
        {
            auto v = sv[i];
            assert(v == 0);
        }
    }}

    cout << "svector Import test..." << endl;
    
    {{
        std::vector<unsigned> vect;
        for (unsigned i = 0; i < 128000; ++i)
        {
            vect.push_back(i);
        }
        
        svector sv;
        sv.import(&vect[0], (unsigned)vect.size());
        bool res = CompareSparseVector(sv, vect);
        if (!res)
        {
            cerr << "0.Bit Plain import test failed" << endl;
            exit(1);
        }
        sv.optimize();
        print_svector_stat(sv);
        res = CompareSparseVector(sv, vect);
        if (!res)
        {
            cerr << "optimized Bit Plain import test failed" << endl;
            exit(1);
        }

        bm::sparse_vector<unsigned, bm::bvector<> > sv_1;
        std::copy(vect.begin(), vect.end(), std::back_inserter(sv_1));
        res = CompareSparseVector(sv_1, vect);
        if (!res)
        {
            cerr << "Bit Plain push_back test failed" << endl;
            exit(1);
        }

        
        bm::sparse_vector<unsigned, bvect>::statistics st;
        sv.calc_stat(&st);
        
        bm::sparse_vector<unsigned, bvect> sv2(sv);
        res = CompareSparseVector(sv2, vect);
        if (!res)
        {
            cerr << "Bit Plain copy-ctor test failed" << endl;
            exit(1);
        }
        
        sv2.clear();
        sv2.import(&vect[0], (unsigned)vect.size());
        res = CompareSparseVector(sv2, vect);
        if (!res)
        {
            cerr << "Bit Plain copy-ctor test failed" << endl;
            exit(1);
        }

        bm::sparse_vector<unsigned, bvect> sv3;
        sv3.set(65536, 10); // set some bit to initiate it
        sv3 = sv;
        res = CompareSparseVector(sv3, vect);
        if (!res)
        {
            cerr << "Bit Plain assignmnet test failed" << endl;
            exit(1);
        }
        
        sv3.clear();
        sv3.import(&vect[0], (unsigned)vect.size());
        res = CompareSparseVector(sv3, vect);
        if (!res)
        {
            cerr << "Bit Plain assignment test failed" << endl;
            exit(1);
        }
    }}

    cout << "Import test 64..." << endl;
    
    {{
        std::vector<unsigned long long> vect;
        for (unsigned long long i = 0; i < 128000; ++i)
        {
            vect.push_back(i << 33);
        }
        
        {{
            svector64 sv;
            unsigned long long v = 3;
            v = v << 33;
            sv.set(10, v);
            unsigned long long v1;
            v1 = sv.get(10);
            if (v != v1)
            {
                cerr << "SV 64-bit set failed" << endl;
                exit(1);
            }
        }}
        
        svector64 sv;
        sv.import(&vect[0], (unsigned)vect.size());
        bool res = CompareSparseVector(sv, vect);
        if (!res)
        {
            cerr << "0.Bit Plain import test failed" << endl;
            exit(1);
        }
        sv.optimize();
        print_svector_stat(sv);
        res = CompareSparseVector(sv, vect);
        if (!res)
        {
            cerr << "optimized Bit Plain import test failed" << endl;
            exit(1);
        }

        svector64 sv_1;
        std::copy(vect.begin(), vect.end(), std::back_inserter(sv_1));
        res = CompareSparseVector(sv_1, vect);
        if (!res)
        {
            cerr << "Bit Plain push_back test failed" << endl;
            exit(1);
        }

        
        svector64::statistics st;
        sv.calc_stat(&st);
        
        svector64 sv2(sv);
        res = CompareSparseVector(sv2, vect);
        if (!res)
        {
            cerr << "Bit Plain copy-ctor test failed" << endl;
            exit(1);
        }
        
        sv2.clear();
        sv2.import(&vect[0], (unsigned)vect.size());
        res = CompareSparseVector(sv2, vect);
        if (!res)
        {
            cerr << "Bit Plain copy-ctor test failed" << endl;
            exit(1);
        }

        svector64 sv3;
        sv3.set(65536, 10); // set some bit to initiate it
        sv3 = sv;
        res = CompareSparseVector(sv3, vect);
        if (!res)
        {
            cerr << "Bit Plain assignmnet test failed" << endl;
            exit(1);
        }
        
        sv3.clear();
        sv3.import(&vect[0], (unsigned)vect.size());
        res = CompareSparseVector(sv3, vect);
        if (!res)
        {
            cerr << "Bit Plain assignment test failed" << endl;
            exit(1);
        }
    }}
    
    cout << "Same value assignment test.." << endl;

    {
        bm::sparse_vector<unsigned, bvect > sv;
        const unsigned max_assign =
                            100 + bm::gap_max_bits * bm::set_sub_array_size;
        {
            bm::sparse_vector<unsigned, bvect>::back_insert_iterator
                                                    bi(sv.get_back_inserter());
            for (unsigned i = 0; i < max_assign; ++i)
            {
                *bi = 1;
            } // for
            bi.flush();
        }
        bm::sparse_vector<unsigned, bvect>::statistics st;
        sv.optimize(0, bvect::opt_compress, &st);
        assert(st.gap_blocks == 1);
        assert(st.bit_blocks == 0);
        for (unsigned i = 0; i < max_assign; ++i)
        {
            auto v = sv[i];
            assert(v == 1);
        } // for
        sv[65536] = 0;
    }
    
    
    cout << "Linear assignment test" << endl;
    {{
    std::vector<unsigned> vect(128000);
    bm::sparse_vector<unsigned, bvect > sv;
    bm::sparse_vector<unsigned, bvect > sv1(bm::use_null);
    bm::sparse_vector<unsigned, bvect > sv2;
    
    {
    bm::sparse_vector<unsigned, bvect>::back_insert_iterator bi(sv2.get_back_inserter());
        for (unsigned i = 0; i < 128000; ++i)
        {
            vect[i] = i;
            sv.set(i, i);
            sv1.set(i, i);
            *bi = i;
        }
        bi.flush();
    }
    

    bool res = CompareSparseVector(sv, vect);
    if (!res)
    {
        cerr << "linear assignment test failed" << endl;
        exit(1);
    }
    res = CompareSparseVector(sv1, vect);
    if (!res)
    {
        cerr << "linear assignment test failed (2)" << endl;
        exit(1);
    }
    res = CompareSparseVector(sv2, vect);
    if (!res)
    {
        cerr << "linear assignment test failed (3 - back_inserter)" << endl;
        exit(1);
    }
    }}
    
    cout << "Extract test" << endl;
    
    {{
    bm::sparse_vector<unsigned, bm::bvector<> > sv;
    for (unsigned i = 0; i < 16; ++i)
    {
        sv.set(i, 8);
    }
    for (unsigned i =32; i < 48; ++i)
    {
        sv.set(i, 255);
    }
        
    std::vector<unsigned> v1(16);
    std::vector<unsigned> v1r(16);
    std::vector<unsigned> v1p(16);
    
    sv.extract(&v1[0], 16, 0);
    sv.extract_range(&v1r[0], 16, 0);
    sv.extract_plains(&v1p[0], 16, 0);
    for (unsigned i = 0; i < 16; ++i)
    {
        if (v1[i] != 8 || v1r[i] != v1[i] || v1p[i] != v1[i])
        {
            cerr << "Extract 1 failed at:" << i << endl;
            exit(1);
        }
    } // for
    
    std::vector<unsigned> v2(10);
    std::vector<unsigned> v2r(10);
    std::vector<unsigned> v2p(10);
    
    sv.extract(&v2[0], 10, 32);
    sv.extract_range(&v2r[0], 10, 32);
    sv.extract_plains(&v2p[0], 10, 32);
        
    for (unsigned i = 0; i < 10; ++i)
    {
        if (v2[i] != 255 || v2r[i] != v2[i] || v2p[i] != v2[i])
        {
            cerr << "Extract 2 failed at:" << i << "=" << v2[i] << " r=" << v2r[i] << " sv=" << sv[i+32] << endl;
            exit(1);
        }
    } // for

    }}
    
    
    // test automatic optimization with back_insert iterator
    {
       bm::sparse_vector<unsigned, bm::bvector<> > sv;
       {
           bm::sparse_vector<unsigned, bm::bvector<>>::back_insert_iterator
                                                    bi = sv.get_back_inserter();
           for (unsigned i = 0; i < 65536; ++i)
           {
                bi = 15;
           } // for
       }
       
       bm::sparse_vector<unsigned, bm::bvector<>>::statistics st;
       
       sv.calc_stat(&st);
       assert(st.bit_blocks == 0);
       assert(st.gap_blocks == 0);
       
       {
           bm::sparse_vector<unsigned, bm::bvector<>>::back_insert_iterator
                                            bi = sv.get_back_inserter();
           for (unsigned i = 0; i < 65536; ++i)
           {
                bi = 15;
           }
       }
       
       sv.calc_stat(&st);
       assert(st.bit_blocks == 0);
       assert(st.gap_blocks == 0);
    }

    
    // test insert() / erase()
    {
        bm::sparse_vector<unsigned, bm::bvector<> > sv1;
        sv1.push_back(1);
        sv1.push_back(2);
        sv1.insert(0, 17);
        sv1.insert(4, 18);

        assert(sv1.size() == 5);
        assert(sv1[0] == 17);
        assert(sv1[1] == 1);
        assert(sv1[2] == 2);
        assert(sv1[4] == 18);
        
        sv1.erase(1);
        
        assert(sv1.size() == 4);
        assert(sv1[0] == 17);
        assert(sv1[1] == 2);
        assert(sv1[3] == 18);
    }
    
    // test lower bound search
    cout << "sparse vector (unsigned) lower_bound()" << endl;
    
    {
        bm::sparse_vector<unsigned, bm::bvector<> > sv1;
        sv1.push_back(1);
        sv1.push_back(2);
        sv1.push_back(2);
        sv1.push_back(2);
        sv1.push_back(20);
        sv1.push_back(2000);

        bvect::size_type pos;
        bool found;
        
        bm::sparse_vector_scanner<bm::sparse_vector<unsigned, bm::bvector<> > > scanner;
        found = scanner.lower_bound(sv1, 0u, pos);
        assert(!found);

        found = scanner.lower_bound(sv1, 1u, pos);
        assert(found);
        assert(pos == 0);

        found = scanner.lower_bound(sv1, 2u, pos);
        assert(found);
        assert(pos == 1);

        found = scanner.lower_bound(sv1, 3u, pos);
        assert(!found);

        found = scanner.lower_bound(sv1, 20u, pos);
        assert(found);

        found = scanner.lower_bound(sv1, 2000u, pos);
        assert(found);
    }
    
    
    
    {{
    cout << "sparse vector inc test" << endl;
    bm::sparse_vector<unsigned, bvect > sv;
    
    for (unsigned i = 1; i < 65536; ++i)
    {
        for (unsigned j = 0; j < 200000; ++j)
        {
            sv.inc(j);
            unsigned v = sv.get(j);
            assert(v == i);
        } // for j
        if ((i % 200) == 0)
            cout << "\r" << i << " / " << 65536 << flush;
    } // for i

    cout <<  "\nOk" << endl;
    }}
    


    {{
    cout << "Dynamic range clipping test 1" << endl;
    bm::sparse_vector<unsigned, bvect > sv;

    unsigned i;
    for (i = 0; i < 16; ++i)
    {
        sv.set(i, i);
    }
    bm::dynamic_range_clip_high(sv, 3);
    
    for (i = 0; i < 16; ++i)
    {
        unsigned v = sv[i];
        if (v > 15)
        {
            cerr << "Clipped Value cmpr failed at:" << i << "=" << v << endl;
            exit(1);
        }
        
    } // for i
    cout << "Ok" << endl;
    
    }}
    
    
    {{
        cout << "Resize test" << endl;
        bm::sparse_vector<unsigned, bvect > sv;
        bm::sparse_vector<unsigned, bvect > sv1(bm::use_null);
        unsigned i;
        for (i = 0; i < 16; ++i)
        {
            sv.set(i, i);
            sv1.set(i, i);
        }
        if (sv.size()!= 16)
        {
            cerr << "1.Incorrect sparse vector size:" << sv.size() << endl;
            exit(1);
        }
        
        
        const bvect* bv_null1 = sv1.get_null_bvector();
        assert(bv_null1);
        if (bv_null1->count() != sv1.size())
        {
            cerr << "1.1. Incorrect sparse vector size() - NOT NULL comparison" << sv1.size() << " " << bv_null1->count() << endl;
        }
        
        sv.resize(10);
        sv1.resize(10);
        if (sv.size()!= 10 || sv1.size() != 10)
        {
            cerr << "2.Incorrect sparse vector size:" << sv.size() << endl;
            exit(1);
        }
        if (bv_null1->count() != sv1.size())
        {
            cerr << "2.1. Incorrect sparse vector size() - NOT NULL comparison" << sv1.size() << " " << bv_null1->count() << endl;
        }

        
        cout << "check values for size()=" << sv.size() << endl;
        for (i = 0; i < sv.size(); ++i)
        {
            unsigned v = sv[i];
            if (v != i)
            {
                cerr << "Wrong sparse vector value: at[" << i << "]=" << v << endl;
                exit(1);
            }
            assert(!sv1[i].is_null());
            v = sv1[i];
            if (v != i)
            {
                cerr << "Wrong null sparse vector value: at[" << i << "]=" << v << endl;
                exit(1);
            }
        }
        
        sv.resize(20);
        sv1.resize(20);
        if (sv.size() != 20 || sv1.size() != 20)
        {
            cerr << "3.Incorrect sparse vector size:" << sv.size() << endl;
            exit(1);
        }
        cout << "check values for size()=" << sv.size() << endl;
        for (i = 0; i < sv.size(); ++i)
        {
            unsigned v = sv[i];
            unsigned v1 = sv1[i];
            
            bool b_null = sv[i].is_null();
            bool b1_null = sv1[i].is_null();
            
            if (i < 10)
            {
                if (v != i || v1 != i)
                {
                    cerr << "Wrong sparse vector value: at[" << i << "]=" << v << endl;
                    exit(1);
                }
                assert(!b_null);
                assert(!b1_null);
            }
            else
            {
                if (v != 0 || v1 != 0)
                {
                    cerr << "Wrong sparse (non-zero) vector value " << v << endl;
                    exit(1);
                }
                assert(!b_null);
                assert(b1_null);
            }
        } // for i
        
        sv.resize(0);
        sv1.resize(0);
        if (sv.size()!= 0 || sv1.size() != 0)
        {
            cerr << "2.Incorrect sparse vector size:" << sv.size() << endl;
            exit(1);
        }
        if (bv_null1->count() != 0)
        {
            cerr << "3. Incorrect sparse vector size() - NOT NULL comparison" << sv1.size() << " " << bv_null1->count() << endl;
        }

        
        sv.resize(65536);
        sv1.resize(65536);
        if (bv_null1->count() != 0)
        {
            cerr << "4. Incorrect sparse vector size() - NOT NULL comparison" << sv1.size() << " " << bv_null1->count() << endl;
        }

        for (i = 0; i < sv.size(); ++i)
        {
            unsigned v = sv[i];
            unsigned v1 = sv1[i];
            if (v || v1)
            {
                if (v != 0)
                {
                    cerr << "Wrong sparse (non-zero) vector value: " << v << endl;
                    exit(1);
                }
            }
            assert(sv1[i].is_null());
        }
    
    }}

    // back insert test
    {{
        bm::sparse_vector<unsigned, bvect > sv1(bm::use_null);
        {
            auto bit = sv1.get_back_inserter();
            bit = 10;
            bit = 11;
            bit.add_null();
            bit = 13;
        }
        assert(sv1.size() == 4);

        assert(sv1.is_null(2));
        assert(!sv1.is_null(0));
        assert(!sv1.is_null(1));
        assert(!sv1.is_null(3));

    }}
    
    {{
    cout << "Dynamic range clipping test 2" << endl;
    bm::sparse_vector<unsigned, bvect > sv;

    unsigned i;
    for (i = 128000; i < 128000 * 3; ++i)
    {
        sv.set(i, 7 + rand() % 256);
    }
    bm::dynamic_range_clip_high(sv, 3);
    
    for (i = 0; i < 128000 * 3; ++i)
    {
        unsigned v = sv[i];
        if (i < 128000)
        {
            if (v != 0)
            {
                cerr << "Value cmpr failed at:" << i << "=" << v << endl;
                exit(1);
            }
        }
        else
        {
            if (v > 15)
            {
                cerr << "Clipped Value cmpr failed at:" << i << "=" << v << endl;
                exit(1);
            }
        }
        
    } // for i
    cout << "Ok" << endl;
    
    }}
    
    
    {{
    cout << "Dynamic range clipping test 3" << endl;
    bm::sparse_vector<unsigned, bvect > sv;

    unsigned i;
    for (i = 0; i <= 16; ++i)
    {
        sv.set(i, 1);
    }
    for (i = 17; i < 32; ++i)
    {
        sv.set(i, 3);
    }


    bm::dynamic_range_clip_low(sv, 3);
    for (i = 0; i < 32; ++i)
    {
        unsigned v = sv[i];
        if (v != 8)
        {
            cerr << "Low Clipped Value cmpr failed at:" << i << "=" << v << endl;
            exit(1);
        }
        
    } // for i
    cout << "Ok" << endl;
    
    }}


    cout << "Test Sparse vector join" << endl;
    {
        bm::sparse_vector<unsigned, bvect> sv1;
        bm::sparse_vector<unsigned, bvect> sv2;
        
        sv1.set(0, 0);
        sv1.set(1, 1);
        sv1.set(2, 2);

        sv2.set(3, 3);
        sv2.set(4, 4);
        sv2.set(5, 5);

        sv1.join(sv2);
        
        if (sv1.size()!=6)
        {
            cerr << "Sparse join size failed:" << sv1.size() << endl;
            exit(1);
        
        }
        for (unsigned i = 0; i < sv1.size(); ++i)
        {
            unsigned v1 = sv1[i];
            if (v1 != i)
            {
                cerr << "Sparse join cmp failed:" << sv1.size() << endl;
                exit(1);
            }
        }
    }

    cout << "Test Sparse vector merge" << endl;
    {
        bm::sparse_vector<unsigned, bvect> sv1;
        bm::sparse_vector<unsigned, bvect> sv2;
        
        sv1.set(0, 0);
        sv1.set(1, 1);
        sv1.set(2, 2);

        sv2.set(3, 3);
        sv2.set(4, 4);
        sv2.set(5, 5);

        sv1.merge(sv2);
        
        if (sv1.size()!=6)
        {
            cerr << "Sparse merge size failed:" << sv1.size() << endl;
            exit(1);
        
        }
        for (unsigned i = 0; i < sv1.size(); ++i)
        {
            unsigned v1 = sv1[i];
            if (v1 != i)
            {
                cerr << "Sparse join cmp failed:" << sv1.size() << endl;
                exit(1);
            }
        }
    }


    cout << "Test Sparse vector join with NULL-able" << endl;
    {
        bm::sparse_vector<unsigned, bvect> sv1;
        bm::sparse_vector<unsigned, bvect> sv2(bm::use_null);

        assert(!sv1.is_nullable());
        
        sv1.set(0, 0);
        sv1.set(1, 1);
        sv1.set(2, 2);

        sv2.set(3, 3);
        sv2.set(4, 4);
        sv2.set(5, 5);

        sv1.join(sv2);
        assert(!sv1.is_nullable());
        
        if (sv1.size()!=6)
        {
            cerr << "Sparse join size failed:" << sv1.size() << endl;
            exit(1);
        
        }
        for (unsigned i = 0; i < sv1.size(); ++i)
        {
            unsigned v1 = sv1[i];
            if (v1 != i)
            {
                cerr << "Sparse join cmp failed:" << sv1.size() << endl;
                exit(1);
            }
            assert(!sv1[i].is_null());
        }
    }

    cout << "Test Sparse vector join NULL-able with not NULL-able" << endl;
    {
        bm::sparse_vector<unsigned, bvect> sv1(bm::use_null);
        bm::sparse_vector<unsigned, bvect> sv2;

        assert(sv1.is_nullable());
        
        //sv1.set(0, 0);
        sv1.set(1, 1);
        sv1.set(2, 2);

        sv2.set(3, 3);
        sv2.set(4, 4);
        sv2.set(5, 5);

        sv1.join(sv2);
        assert(sv1.is_nullable());
        
        if (sv1.size()!=6)
        {
            cerr << "Sparse join size failed:" << sv1.size() << endl;
            exit(1);
        }
        for (unsigned i = 0; i < sv1.size(); ++i)
        {
            unsigned v1 = sv1[i];
            if (v1 != i)
            {
                cerr << "Sparse join cmp failed:" << i << endl;
                exit(1);
            }
            assert(!sv1[i].is_null());
        }
    }

    cout << "Test Sparse vector join NULL-able with NULL-able" << endl;
    {
        bm::sparse_vector<unsigned, bvect> sv1(bm::use_null);
        bm::sparse_vector<unsigned, bvect> sv2(bm::use_null);

        assert(sv1.is_nullable());
        assert(sv2.is_nullable());
        
        //sv1.set(0, 0);
        sv1.set(1, 1);
        sv1.set(2, 2);

        //sv2.set(3, 3);
        sv2.set(4, 4);
        sv2.set(5, 5);

        sv1.join(sv2);
        assert(sv1.is_nullable());
        
        if (sv1.size()!=6)
        {
            cerr << "Sparse join size failed:" << sv1.size() << endl;
            exit(1);
        }
        for (unsigned i = 0; i < sv1.size(); ++i)
        {
            unsigned v1 = sv1[i];
            if (v1 != i)
            {
                if (i == 0 || i == 3) // legitimate test-case exceptions
                {
                }
                else
                {
                    cerr << "Sparse join cmp failed:" << i << endl;
                    exit(1);
                }
            }
            if (sv1[i].is_null())
            {
                assert(i == 0 || i == 3);
            }
        }
    }
    
    cout << "check if optimize keeps the NULL vector" << std::endl;
    {
        bm::sparse_vector<unsigned, bvect> sv(bm::use_null);
        assert(sv.is_nullable());
        sv.optimize();
        assert(sv.is_nullable());
    }
    

    {
        bm::sparse_vector<unsigned, bvect> sv1;
        bm::sparse_vector<unsigned, bvect> sv2;
        bm::sparse_vector<unsigned, bvect> sv3;
        
        unsigned i;
        for (i = 65536; i < 256000; ++i)
        {
            sv1.set(i, 256);
            sv3.set(i, 256);
        }
        for (i = 312000; i < 365636; ++i)
        {
            sv2.set(i, 65536);
            sv3.set(i, 65536);
        }
        
        sv1.optimize();
        sv1.join(sv2);
        
        if (sv1.size() != sv3.size())
        {
            cerr << "Sparse join size failed (2):" << sv1.size() << endl;
            exit(1);
        }
        
        for (i = 0; i < sv1.size(); ++i)
        {
            unsigned v1 = sv1[i];
            unsigned v3 = sv3[i];
            if (v1 != v3)
            {
                cerr << "Sparse join cmp failed (2):" << v1 << "!=" << v3 << endl;
                exit(1);
            }
        } // for i
    }
    cout << "Sparse vector join ok" << endl;
    
    cout << "---------------------------- Bit-plain sparse vector test OK" << endl;
}

static
void TestSparseVectorAlgo()
{
    cout << " -------------------------- TestSparseVectorAlgo()" << endl;

    {
        bm::sparse_vector<unsigned, bvect> sv1;
        bm::sparse_vector<unsigned, bvect> sv2;
        sv1.push_back(1);
        sv1.push_back(1);
        sv1.push_back(1);

        sv2 = sv1;

        bm::sparse_vector<unsigned, bvect>::size_type pos;
        bool f;
        f = bm::sparse_vector_find_first_mismatch(sv1, sv2, pos);
        assert(!f);
        f = bm::sparse_vector_find_first_mismatch(sv2, sv1, pos);
        assert(!f);
        sv2.push_back(4);
        f = bm::sparse_vector_find_first_mismatch(sv1, sv2, pos);
        assert(f);
        assert(pos == 3);
        f = bm::sparse_vector_find_first_mismatch(sv2, sv1, pos);
        assert(f);
        assert(pos == 3);

        sv1.optimize();
        f = bm::sparse_vector_find_first_mismatch(sv2, sv1, pos);
        assert(f);
        assert(pos == 3);

        sv2.optimize();
        f = bm::sparse_vector_find_first_mismatch(sv2, sv1, pos);
        assert(f);
        assert(pos == 3);
    }

    // sparse with NULLs test
    {
        bm::sparse_vector<unsigned, bvect> sv1(bm::use_null);
        bm::sparse_vector<unsigned, bvect> sv2(bm::use_null);
        sv1[100] = 1;
        sv1[1000] = 1;
        sv1[bm::id_max32/2 + 1000] = 1;

        sv2 = sv1;

        bm::sparse_vector<unsigned, bvect>::size_type pos;
        bool f;
        f = bm::sparse_vector_find_first_mismatch(sv1, sv2, pos);
        assert(!f);
        f = bm::sparse_vector_find_first_mismatch(sv2, sv1, pos);
        assert(!f);
        sv2.push_back(4);
        f = bm::sparse_vector_find_first_mismatch(sv1, sv2, pos);
        assert(f);
        assert(pos == bm::id_max32/2 + 1000+1);
        f = bm::sparse_vector_find_first_mismatch(sv2, sv1, pos);
        assert(f);
        assert(pos == bm::id_max32/2 + 1000+1);

        sv1.optimize();
        f = bm::sparse_vector_find_first_mismatch(sv2, sv1, pos);
        assert(f);
        assert(pos == bm::id_max32/2 + 1000+1);

        sv2.optimize();
        f = bm::sparse_vector_find_first_mismatch(sv2, sv1, pos);
        assert(f);
        assert(pos == bm::id_max32/2 + 1000+1);
    }

    {
        bm::sparse_vector<unsigned, bvect>::size_type pos;
        bm::sparse_vector<unsigned, bvect> sv1(bm::use_null);
        bm::sparse_vector<unsigned, bvect> sv2(bm::use_null);

        sv1[1] = 1;
        sv1[2] = 2;
        sv1.set_null(3); // set element 3 to NULL
        sv1[4] = 0;

        sv2 = sv1;

        bool found = bm::sparse_vector_find_first_mismatch(sv1, sv2, pos);
        assert(!found);

        sv2[4] = 10;
        found = bm::sparse_vector_find_first_mismatch(sv1, sv2, pos);
        assert(found);
        assert(pos == 4);

        sv2[3] = 0;
        found = bm::sparse_vector_find_first_mismatch(sv1, sv2, pos);
        assert(found);
        assert(pos == 3);
    }

    {
        bm::sparse_vector<unsigned, bvect>::size_type pos;

        bm::sparse_vector<unsigned, bvect> sv1(bm::use_null);
        bm::sparse_vector<unsigned, bvect> sv2;

        bool found = bm::sparse_vector_find_first_mismatch(sv1, sv2, pos);
        assert(!found);

        sv1[0] = 0;
        sv1[10] = 1;
        sv1[20] = 2;
        sv1.set_null(30); // set element 3 to NULL
        sv1[40] = 0;

        sv2[0] = 0;
        sv2[10] = 1;
        sv2[20] = 2;
        sv2[40] = 0;

        found = bm::sparse_vector_find_first_mismatch(sv1, sv2, pos);
        assert(found);
        assert(pos == 1);

        found = bm::sparse_vector_find_first_mismatch(sv2, sv1, pos);
        assert(found);
        assert(pos == 1);

    }

    cout << " ----- Find mismatches " << endl;

    {
        bvect  bv_m; // mismatch vector
        bm::sparse_vector<unsigned, bvect> sv1;
        bm::sparse_vector<unsigned, bvect> sv2;

        bm::sparse_vector_find_mismatch(bv_m, sv1, sv2, bm::no_null);
        assert(!bv_m.any());

        sv1[0] = 0;
        sv1[10] = 15;
        sv1[20] = 23;

        sv2[0] = 0;
        sv2[10] = 15;
        sv2[20] = 23;

        bm::sparse_vector_find_mismatch(bv_m, sv1, sv2, bm::no_null);
        assert(!bv_m.any());

        sv2[0] = 32;
        bm::sparse_vector_find_mismatch(bv_m, sv1, sv2, bm::no_null);
        assert(bv_m.count()==1);
        {
            bvect bv_c { 0 };
            bool f = bv_m.equal(bv_c);
            assert(f);
        }

        sv1[22] = 255;
        bm::sparse_vector_find_mismatch(bv_m, sv1, sv2, bm::no_null);
        cout << bv_m.count() << endl;
        assert(bv_m.count()==3);
        {
            bvect bv_c { 0,  21, 22 };
            bool f = bv_m.equal(bv_c);
            assert(f);
        }
    }

    {
        bvect  bv_m; // mismatch vector
        bm::sparse_vector<unsigned, bvect> sv1(bm::use_null);
        bm::sparse_vector<unsigned, bvect> sv2(bm::use_null);

        sv1[0] = 0;
        sv1[10] = 15;
        sv1[20] = 23;

        sv2[0] = 0;
        sv2[10] = 15;
        sv2[20] = 23;

        bm::sparse_vector_find_mismatch(bv_m, sv1, sv2, bm::no_null);
        assert(!bv_m.any());
        bm::sparse_vector_find_mismatch(bv_m, sv2, sv1, bm::no_null);
        assert(!bv_m.any());

        sv2[id_max/2] = 0;
        sv2[id_max/2+1] = 256;

        bm::sparse_vector_find_mismatch(bv_m, sv1, sv2, bm::no_null);
        cout << bv_m.count() << endl;
        assert(bv_m.count()==2);
        {
            bvect bv_c { id_max/2,  id_max/2+1 };
            DetailedCompareBVectors(bv_c, bv_m);
            bool f = bv_m.equal(bv_c);
            assert(f);
        }
        bm::sparse_vector_find_mismatch(bv_m, sv2, sv1, bm::no_null);
        cout << bv_m.count() << endl;
        assert(bv_m.count()==2);
        {
            bvect bv_c { id_max/2,  id_max/2+1 };
            DetailedCompareBVectors(bv_c, bv_m);
            bool f = bv_m.equal(bv_c);
            assert(f);
        }
    }


    {
        bvect  bv_m; // mismatch vector
        bm::sparse_vector<unsigned, bvect> sv1;
        bm::sparse_vector<unsigned, bvect> sv2(bm::use_null);

        sv1[0] = 0;
        sv1[1] = 15;
        sv1[2] = 23;

        sv2[0] = 0;
        sv2[1] = 15;
        sv2[2] = 23;

        bm::sparse_vector_find_mismatch(bv_m, sv1, sv2, bm::no_null);
        assert(!bv_m.any());
        bm::sparse_vector_find_mismatch(bv_m, sv1, sv2, bm::use_null);
        cout << bv_m.count() << endl;
        assert(!bv_m.any());
        bm::sparse_vector_find_mismatch(bv_m, sv2, sv1, bm::use_null);
        cout << bv_m.count() << endl;
        assert(!bv_m.any());

        sv2[4] = 10;
        bm::sparse_vector_find_mismatch(bv_m, sv1, sv2, bm::use_null);
        cout << bv_m.count() << endl;
        {
            bvect bv_c { 3, 4 };
            DetailedCompareBVectors(bv_c, bv_m);
            bool f = bv_m.equal(bv_c);
            assert(f);
        }
        bm::sparse_vector_find_mismatch(bv_m, sv2, sv1, bm::use_null);
        cout << bv_m.count() << endl;
        {
            bvect bv_c { 3, 4 };
            DetailedCompareBVectors(bv_c, bv_m);
            bool f = bv_m.equal(bv_c);
            assert(f);
        }

    }

    {
        bvect  bv_m; // mismatch vector
        bm::sparse_vector<unsigned, bvect> sv1;
        bm::sparse_vector<unsigned, bvect> sv2(bm::use_null);

        sv1[0] = 0;
        sv1[1] = 1;
        sv1[2] = 2;

        sv2[0] = 0;
        sv2[1] = 1;
        sv2[2] = 2;

        sv2[3] = 0;
        bm::sparse_vector_find_mismatch(bv_m, sv1, sv2, bm::use_null);
        cout << bv_m.count() << endl;
        {
            bvect bv_c { 3 };
            DetailedCompareBVectors(bv_c, bv_m);
            bool f = bv_m.equal(bv_c);
            assert(f);
        }
        bm::sparse_vector_find_mismatch(bv_m, sv2, sv1, bm::use_null);
        cout << bv_m.count() << endl;
        {
            bvect bv_c { 3 };
            DetailedCompareBVectors(bv_c, bv_m);
            bool f = bv_m.equal(bv_c);
            assert(f);
        }


        sv1[5] = 0;
        bm::sparse_vector_find_mismatch(bv_m, sv2, sv1, bm::use_null);
        cout << bv_m.count() << endl;
        {
            bvect bv_c { 4, 5 };
            DetailedCompareBVectors(bv_c, bv_m);
            bool f = bv_m.equal(bv_c);
            assert(f);
        }
        bm::sparse_vector_find_mismatch(bv_m, sv2, sv1, bm::no_null);
        cout << bv_m.count() << endl;
        assert(!bv_m.any());

        bm::sparse_vector_find_mismatch(bv_m, sv1, sv2, bm::no_null);
        cout << bv_m.count() << endl;
        assert(!bv_m.any());
    }



    cout << " -------------------------- TestSparseVectorAlgo() OK" << endl;
}


static
void TestSparseVector_XOR_Scanner()
{
    cout << " -------------------------- TestSparseVector_XOR_Scanner()" << endl;
    BM_DECLARE_TEMP_BLOCK(tb)

    // XOR scanner EQ test
    {{
        bm::sparse_vector<unsigned, bvect> sv;
        sv.push_back(9);
        sv.push_back(9);

        bm::xor_scanner<bvect> xscan;
        bm::xor_scanner<bvect>::bv_ref_vector_type r_vect;
        r_vect.build(sv.get_bmatrix());
        xscan.set_ref_vector(&r_vect);

        const bvect* bv_x = sv.get_plain(0);
        const bvect::blocks_manager_type& bman_x = bv_x->get_blocks_manager();
        const bm::word_t* block_x = bman_x.get_block_ptr(0, 0);

        xscan.compute_x_block_stats(block_x);

        auto idx = xscan.get_ref_vector().find(0);
        assert(idx == 0);

        bool f = xscan.search_best_xor_mask(block_x,
                                            1, xscan.get_ref_vector().size(),
                                            0, 0, tb);
        assert(f);
        idx = xscan.found_ridx();
        assert(idx == 1);
        assert(xscan.get_x_best_metric() == 0); // EQ
        assert(xscan.is_eq_found());
        idx = xscan.get_ref_vector().get_row_idx(idx);
        assert(idx == 3); // matrix row 3
    }}

    {{
        bm::sparse_vector<unsigned, bvect> sv;
        for (unsigned i = 0; i < 65536; ++i)
            sv.push_back(9);

        bm::xor_scanner<bvect> xscan;
        bm::xor_scanner<bvect>::bv_ref_vector_type r_vect;
        r_vect.build(sv.get_bmatrix());
        xscan.set_ref_vector(&r_vect);

        const bvect* bv_x = sv.get_plain(0);
        const bvect::blocks_manager_type& bman_x = bv_x->get_blocks_manager();
        const bm::word_t* block_x = bman_x.get_block_ptr(0, 0);

        xscan.compute_x_block_stats(block_x);

        auto idx = xscan.get_ref_vector().find(0);
        assert(idx == 0);
        auto sz = xscan.get_ref_vector().size();
        bool f = xscan.search_best_xor_mask(block_x,
                                            1, sz,
                                            0, 0, tb);
        assert(f);
        idx = xscan.found_ridx();
        assert(idx == 1);
        assert(xscan.get_x_best_metric() == 0); // EQ
        assert(xscan.is_eq_found());
        idx = xscan.get_ref_vector().get_row_idx(idx);
        assert(idx == 3); // matrix row 3

    }}

    {{
        bm::sparse_vector<unsigned, bvect> sv;
        for (unsigned i = 0; i < 65536; i+=2)
        {
            sv.push_back(1);
            sv.push_back(8);
        }
        bm::xor_scanner<bvect> xscan;
        bm::xor_scanner<bvect>::bv_ref_vector_type r_vect;
        r_vect.build(sv.get_bmatrix());
        xscan.set_ref_vector(&r_vect);

        const bvect* bv_x = sv.get_plain(0);
        const bvect::blocks_manager_type& bman_x = bv_x->get_blocks_manager();
        const bm::word_t* block_x = bman_x.get_block_ptr(0, 0);

        xscan.compute_x_block_stats(block_x);

        auto idx = xscan.get_ref_vector().find(0);
        assert(idx == 0);
        auto sz = xscan.get_ref_vector().size();
        bool f = xscan.search_best_xor_mask(block_x,
                                            1, sz,
                                            0, 0, tb);
        assert(f);
        idx = xscan.found_ridx();
        assert(idx == 1);
        assert(xscan.get_x_best_metric() == 1);
        assert(!xscan.is_eq_found());
        idx = xscan.get_ref_vector().get_row_idx(idx);
        bm::id64_t d64 = xscan.get_xor_digest();
        assert(d64 == ~bm::id64_t(0));
        assert(idx == 3); // matrix row 3
    }}

    cout << " -------------------------- TestSparseVector_XOR_Scanner() OK" << endl;
}



static
void TestSparseVectorSerial()
{
    cout << "---------------------------- Test sparse vector serializer" << endl;

    bm::sparse_vector_serializer<sparse_vector_u32> sv_ser;
    sv_ser.set_xor_ref(false);

    for (unsigned pass = 0; pass < 2; ++pass)
    {
        // simple test gather for non-NULL vector
        {
            sparse_vector_u32 sv1;
            sparse_vector_u32 sv2;

            for (sparse_vector_u32::size_type i = 0; i < 10; ++i)
                sv1.push_back(i + 1);
            sparse_vector_serial_layout<sparse_vector_u32> sv_lay;
            sv_ser.serialize(sv1, sv_lay);
            const unsigned char* buf = sv_lay.buf();

            bm::sparse_vector_deserializer<sparse_vector_u32> sv_deserial;

            sparse_vector_u32::bvector_type bv_mask;
            bv_mask.set(0);
            bv_mask.set(2);
            sv_deserial.deserialize(sv2, buf, bv_mask);

            assert(sv2.size() == sv1.size());
            assert(sv2.get(0) == 1);
            cout << sv2.get(1) << endl;
            assert(sv2.get(1) == 0);
            assert(sv2.get(2) == 3);

            sparse_vector_u32::statistics st;
            sv2.calc_stat(&st);
            assert(!st.bit_blocks);
            assert(st.gap_blocks);
        }

        // simple test gather for NULL-able vector
        {
            sparse_vector_u32 sv1(bm::use_null);
            sparse_vector_u32 sv2(sv1);

            for (sparse_vector_u32::size_type i = 0; i < 100; i += 2)
            {
                sv1[i] = i + 1;
            }
            sparse_vector_serial_layout<sparse_vector_u32> sv_lay;
            sv_ser.serialize(sv1, sv_lay);
            const unsigned char* buf = sv_lay.buf();

            bm::sparse_vector_deserializer<sparse_vector_u32> sv_deserial;

            {
                sparse_vector_u32::bvector_type bv_mask;
                bv_mask.set(0);
                bv_mask.set(2);
                bv_mask.set(1024); // out of range mask
                sv_deserial.deserialize(sv2, buf, bv_mask);

                assert(sv2.size() == sv1.size());
                assert(sv2.get(0) == 1);
                assert(sv2.get(1) == 0);
                assert(sv2.get(2) == 3);

                const sparse_vector_u32::bvector_type* bv_null = sv2.get_null_bvector();
                auto cnt = bv_null->count();
                assert(cnt == sv1.get_null_bvector()->count());

                sparse_vector_u32::statistics st;
                sv2.calc_stat(&st);
                //assert(!st.bit_blocks);
                assert(st.gap_blocks);
            }
            {
                sparse_vector_u32::bvector_type bv_mask;
                sv_deserial.deserialize(sv2, buf, bv_mask);
                assert(sv2.size() == sv1.size());
                const sparse_vector_u32::bvector_type* bv_null = sv2.get_null_bvector();
                auto cnt = bv_null->count();
                assert(cnt == sv2.get_null_bvector()->count());
                assert(sv2.get(0) == 0);
                assert(sv2.get(1) == 0);
                assert(sv2.get(2) == 0);
            }
        }

        // stress test gather deserialization
        cout << "Gather deserialization stress test..." << endl;
        {
            sparse_vector_u32::size_type from, to;
            sparse_vector_u32 sv1(bm::use_null);
            sparse_vector_u32 sv2(sv1);
            sparse_vector_u32 sv3(sv1);

            from = bm::id_max / 2;
            to = from + 75538;

            unsigned cnt = 0;
            for (sparse_vector_u32::size_type i = from; i < to; ++i, ++cnt)
            {
                if (cnt % 10 == 0)
                    sv1.set_null(i);
                else
                    sv1.set(i, cnt);
            } // for i
            sparse_vector_serial_layout<sparse_vector_u32> sv_lay;
            sv_ser.serialize(sv1, sv_lay);
            const unsigned char* buf = sv_lay.buf();

            {
                bm::sparse_vector_deserializer<sparse_vector_u32> sv_deserial;
                sparse_vector_u32 sv4(bm::use_null);
                sv_deserial.deserialize(sv4, buf);
                {
                    bool is_eq = sv1.equal(sv4);
                    assert(is_eq);
                }


                auto i = from;
                auto j = to;
                for (i = from; i < j; ++i, --j)
                {
                    sparse_vector_u32::bvector_type bv_mask;
                    bv_mask.set_range(i, j);

                    sv_deserial.deserialize(sv2, buf, bv_mask);
                    assert(sv2.size() == sv1.size());

                    sv_deserial.deserialize(sv3, buf, i, j);
                    bool is_eq = sv2.equal(sv3);
                    if (!is_eq)
                    {
                        cerr << "Error: Range deserialization equality failed!" << endl;
                        assert(0); exit(1);
                    }
                    sparse_vector_u32::size_type pos;
                    bool found;

                    sparse_vector_u32 sv_filt(sv1);
                    sv_filt.filter(bv_mask);
                    sparse_vector_u32 sv_range(bm::use_null);
                    sv_range.copy_range(sv1, i, j);

                    is_eq = sv_filt.equal(sv_range);
                    assert(is_eq);

                    //sv3.filter(bv_mask);

                    found = bm::sparse_vector_find_first_mismatch(sv_filt, sv3, pos, bm::no_null);
                    if (found)
                    {
                        found = bm::sparse_vector_find_first_mismatch(sv_filt, sv3, pos, bm::no_null);

                        auto vf = sv_filt.get(pos);
                        auto v1 = sv1.get(pos);
                        auto v3 = sv3.get(pos);

                        cerr << vf << "!=" << v3 << "!=" << v1 << endl;
                        cerr << "Filter Range deserialization mismatch found! at pos=" << pos << endl;
                        cerr << "[" << i << ".." << j << "]" << endl;
                        assert(0); exit(1);
                    }

                    found = bm::sparse_vector_find_first_mismatch(sv_range, sv2, pos, bm::no_null);
                    if (found)
                    {
                        cerr << "Range deserialization mismatch found! at pos=" << pos << endl;
                        cerr << "[" << i << ".." << j << "]" << endl;
                        assert(0); exit(1);
                    }
                    /*
                                    for (auto k = i; k < j; ++k)
                                    {
                                        auto v1 = sv1.get(k);
                                        auto v2 = sv2.get(k);
                                        if (v1 != v2)
                                        {
                                            cerr << "Error:Range deserialization discrepancy!" << endl;
                                            assert(0); exit(1);
                                        }
                                        auto n1 = sv1.is_null(k);
                                        auto n2 = sv2.is_null(k);
                                        if (n1 != n2)
                                        {
                                            cerr << "Error:Range NULL deserialization discrepancy!" << endl;
                                            assert(0); exit(1);
                                        }
                                    } // for k
                    */
                    if (i % 0xFF == 0)
                    {
                        std::cout << "\r" << j - i << flush;
                    }

                } // for i
            }
            cout << "\nOK\n" << endl;
        }

        sv_ser.set_xor_ref(true);
    } // for pass

    cout << "---------------------------- Test sparse vector serializer OK" << endl;
}


static
void TestSparseVectorInserter()
{
    cout << "---------------------------- Test sparse vector inserter" << endl;
    
    {
        sparse_vector_u32 sv1(bm::use_null);
        sparse_vector_u32 sv2(bm::use_null);
        sparse_vector_u32::back_insert_iterator bi2(sv2.get_back_inserter());
        sparse_vector_u32::back_insert_iterator bi3(bi2);
        sparse_vector_u32::back_insert_iterator bi4;

        assert(bi2.empty());
        assert(bi3.empty());
        
        bi4 = bi2;
        assert(bi4.empty());

        for (unsigned i = 0; i < 1280000; ++i)
        {
            if (i % 100 == 0)
            {
                sv1.set_null(i);
                bi2.add_null();
            }
            else
            {
                sv1.set(i, i);
                *bi2 = i;
            }
            assert(!bi2.empty());
        }
        bi2.flush();
        
        if (!sv1.equal(sv2))
        {
            cout << "ERROR! sparse_vector back_insert_iterator mismatch." << endl;
            exit(1);
        }
    }

    {
        sparse_vector_u32 sv1(bm::use_null);
        sparse_vector_u32 sv2(bm::use_null);
        sparse_vector_u32::back_insert_iterator bi2(sv2.get_back_inserter());
        
        for (unsigned i = 0; i < 1280000; ++i)
        {
            if (i % 100 == 0)
            {
                sv1.set_null(i);
                ++i;
                sv1.set_null(i);
                bi2.add_null(2);
            }
            else
            {
                sv1.set(i, i);
                *bi2 = i;
            }
            if (i % 10000 == 0)
            {
                bi2.flush();
            }

        }
        bi2.flush();
        
        if (!sv1.equal(sv2))
        {
            cout << "ERROR! (2)sparse_vector back_insert_iterator mismatch." << endl;
            exit(1);
        }
    }


    cout << "---------------------------- Bit-plain sparse vector inserter OK" << endl;
}

static
void CheckSparseVectorGather(const sparse_vector_u32& sv,
                             unsigned from, unsigned to, unsigned control_value = 0)
{
    assert(sv.size());
    assert (to >= from);
    
    unsigned gather_size = to - from + 1;
    std::vector<unsigned> target_v;
    std::vector<unsigned> target_v_control;
    std::vector<unsigned> idx_v;
    target_v.resize(gather_size);
    target_v_control.resize(gather_size);
    idx_v.reserve(gather_size);
    
    for (unsigned i = from; i <= to; ++i)
    {
        idx_v.push_back(i);
    }
    sv.decode(target_v_control.data(), from, gather_size);


    sv.gather(target_v.data(), idx_v.data(), gather_size, BM_SORTED);
    for (unsigned i = 0; i < gather_size; ++i)
    {
        unsigned vg = target_v[i];
        unsigned vc = target_v_control[i];
        if (vg != vc)
        {
            cerr << "Error! gather/decode control mismatch " << vc << " " << vg
                 << " at=" << i << endl;
            cerr << control_value << endl;
            exit(1);
        }
    }
    
    sv.gather(target_v.data(), idx_v.data(), gather_size, BM_UNSORTED);
    for (unsigned i = 0; i < gather_size; ++i)
    {
        unsigned vg = target_v[i];
        unsigned vc = target_v_control[i];
        if (vg != vc)
        {
            cerr << "Error! gather/decode control mismatch " << vc << " " << vg
                 << " at=" << i << endl;
            cerr << control_value << endl;
            exit(1);
        }
    }

    sv.gather(target_v.data(), idx_v.data(), gather_size, BM_UNKNOWN);
    for (unsigned i = 0; i < gather_size; ++i)
    {
        unsigned vg = target_v[i];
        unsigned vc = target_v_control[i];
        if (vg != vc)
        {
            cerr << "Error! gather/decode control mismatch " << vc << " " << vg
                 << " at=" << i << endl;
            cerr << control_value << endl;
            exit(1);
        }
    }


#if 0
    // detailed check (very slow)
    unsigned k = 0;
    for (unsigned i = from; i <= to; ++i, ++k)
    {
        unsigned v1 = i;
        unsigned v2 = target_v[k];
        if (v1 != v2)
        {
            if (control_value)
            {
                if (v2 != control_value)
                {
                    cerr << "Error! gather control mismatch " << control_value << " " << v2
                         << " at=" << i << endl;
                    exit(1);
                }
            }
            else
            {
                v1 = sv.get(i);
                if (v1 != v2)
                {
                    cerr << "Error! gather mismatch " << v1 << " " << v2
                         << " at=" << i << endl;
                    exit(1);
                }
            }
        }
    } // for
#endif
}

static
void CheckSparseVectorGatherRandom(const sparse_vector_u32& sv,
                                   unsigned gather_size)
{
    assert(sv.size());
    
    if (gather_size == 0)
        gather_size = 1;
    
    std::vector<unsigned> target_v;
    std::vector<unsigned> idx_v;
    target_v.resize(gather_size);
    idx_v.reserve(gather_size);
    
    for (unsigned i = 0; i < gather_size; ++i)
    {
        unsigned r_idx = unsigned(rand()) % (sv.size()-1);
        idx_v.push_back(r_idx);
    }
    
    sv.gather(target_v.data(), idx_v.data(), gather_size, BM_UNSORTED);

    unsigned k = 0;
    for (unsigned i = 0; i < gather_size; ++i, ++k)
    {
        unsigned v1 = sv.get(idx_v[k]);
        unsigned v2 = target_v[k];
        if (v1 != v2)
        {
            {
                cerr << "Error! random gather mismatch " << v1 << " " << v2
                     << " at=" << i << endl;
                exit(1);
            }
        }
    } // for
}


static
void TestSparseVectorGatherDecode()
{
    cout << "---------------------------- Test sparse vector gather decode" << endl;
    
    sparse_vector_u32 sv;
    sparse_vector_u32 sv2;
    sparse_vector_u32 sv3;
    const unsigned control_value = 9;
    {
        cout << " Filling..." << flush;
        sparse_vector_u32::back_insert_iterator bi(sv.get_back_inserter());
        sparse_vector_u32::back_insert_iterator bi2(sv2.get_back_inserter());

        unsigned max_size = 1280000;
        for (unsigned i = 0; i < max_size; ++i)
        {
            *bi = i;
            *bi2 = control_value;
        }
        bi.flush(); bi2.flush();
        sv2.optimize();
        
        for (unsigned i = 0; i < max_size; i+=200)
        {
            sv3[i] = i;
        }
        sv3.optimize();

        cout << "ok" << endl;
    }
    
    {
        cout << "Test 1 (regular pattern)" << endl;
        unsigned probe_to = 100000;
        time_t      start_time = time(0);
        time_t      finish_time;

        for (unsigned i = 0; i < probe_to; ++i)
        {
            CheckSparseVectorGather(sv, i, i);
            CheckSparseVectorGather(sv2, i, i, control_value);
            CheckSparseVectorGather(sv3, i, i);
            unsigned depth = rand() % 30000;
            CheckSparseVectorGather(sv, i, i+depth);
            CheckSparseVectorGather(sv2, i, i+depth, control_value);
            CheckSparseVectorGather(sv3, i, i+depth);
            if (i % 500 == 0)
            {
                finish_time = time(0);
                cout << "\r" << i << "/" << probe_to
                     << " [" << (finish_time - start_time) << "]" << flush;
                start_time = time(0);
            }
        }
        cout << endl;
    }

    {
        cout << "Test 2 (random pattern)" << endl;
        unsigned probe_to = 100000;
        time_t      start_time = time(0);
        time_t      finish_time;

        for (unsigned i = 0; i < probe_to; ++i)
        {
            unsigned gsize = rand()%2024;
            CheckSparseVectorGatherRandom(sv, gsize);
            CheckSparseVectorGatherRandom(sv2, gsize);
            CheckSparseVectorGatherRandom(sv3, gsize);
            if (i % 500 == 0)
            {
                finish_time = time(0);
                cout << "\r" << i << "/" << probe_to
                     << " [" << (finish_time - start_time) << "]" << flush;
                start_time = time(0);
            }
        }
        cout << endl;
    }


    cout << "---------------------------- Test sparse vector gather decode OK" << endl;
}



template<class SV>
void bvector_transform_11(typename SV::bvector_type& bvect_in,
                          const    SV&               sv_brel,
                          typename SV::bvector_type& bvect_out)
{
    bm::set2set_11_transform<SV> bin_trans;
    bin_trans.run(bvect_in, sv_brel, bvect_out);
}

static
void CheckSparseVectorRange(const sparse_vector_u32& sv,
                            unsigned left, unsigned right)
{
    sparse_vector_u32 sv1(bm::use_null);
    sparse_vector_u32 sv2(sv);
    sv1.copy_range(sv, left, right);
    
    if (right >= sv.size())
    {
        right = sv.size()-1;
    }
    
    if (left == right)
    {
        unsigned v1 = sv.get(left);
        unsigned v2 = sv1[right];
        assert(v1 == v2);
        return;
    }
    
    if (left)
    {
        sv2.clear_range(0, left-1, true);
    }
    if (right < sv2.size()-1)
    {
        sv2.clear_range(right+1, sv2.size()-1, true);
    }
    
    bool same = sv2.equal(sv1);
    if (!same)
    {
        cerr << "Hmmm... Range comaprison failed, detailed check..." << endl;
        cerr << "[" << left << ".." << right << "]" << endl;
        for (unsigned i = left; i <= right; ++i)
        {
            unsigned v1 = sv.get(i);
            unsigned v2 = sv1[i];
            if (v1 != v2)
            {
                cerr << "Error! Copy range check failed at:" << i << endl;
                exit(1);
            }
        } // for
        cerr << "detailed check did not find issues. error in test?" << endl;
        exit(1);
    }
}

static
void TestSparseVectorRange()
{
    cout << " ---------------- Sparse vector Range partitioning test" << endl;

    cout << "Basic check" << endl;
    {
        sparse_vector_u32 sv(bm::use_null);
        sv.set(2, 25);
        sv.set(3, 35);
        sv.set(7, 75);
        sv.set(10, 2);
        sv.set(21, 201);
        
        CheckSparseVectorRange(sv, 0, 0);
        CheckSparseVectorRange(sv, 2, 2);
        CheckSparseVectorRange(sv, 7, 10);
    }

    cout << "Stress check 1 (constant)" << endl;
    {
        sparse_vector_u32 sv(bm::use_null);
        const unsigned sv_max = 120000;
        cout << "Filling the vector" << endl;
        for (unsigned i = 0; i < sv_max; ++i)
        {
            sv.push_back(9);
        }
        
        cout << "Phase 1.." << endl;
        for (unsigned i = 0; i < sv_max; ++i)
        {
            CheckSparseVectorRange(sv, 0, i);
            CheckSparseVectorRange(sv, i, sv_max+10);
            cout << "\r" << i << "/" << sv_max << flush;
        }
        cout << endl;
        
        cout << "\nPhase 2.." << endl;
        unsigned k = sv_max;
        for (unsigned i = 0; i < k; ++i, --k)
        {
            CheckSparseVectorRange(sv, i, k);
        }
        
        sv.optimize();
        
        cout << "Phase 3.." << endl;
        for (unsigned i = 0; i < sv_max; ++i)
        {
            CheckSparseVectorRange(sv, 0, i);
            CheckSparseVectorRange(sv, i, sv_max+10);
        }
        
        cout << "Phase 4.." << endl;
        k = sv_max;
        for (unsigned i = 0; i < k; ++i, --k)
        {
            CheckSparseVectorRange(sv, i, k);
        }

    }

    cout << "\nStress check 2 (liner function)" << endl;
    {
        sparse_vector_u32 sv(bm::use_null);
        const unsigned sv_max = 250000;
        cout << "Filling the vector" << endl;
        for (unsigned i = 0; i < sv_max; ++i)
        {
            sv.push_back(i);
        }
        
        cout << "Phase 2-1.." << endl;
        for (unsigned i = 0; i < sv_max; ++i)
        {
            CheckSparseVectorRange(sv, i, i);
            CheckSparseVectorRange(sv, 0, i);
            CheckSparseVectorRange(sv, i, sv_max+10);
            if (i % 256 == 0)
                cout << "\r" << i << "/" << sv_max << flush;
        }
        
        cout << "\nPhase 2-2.." << endl;
        unsigned k = sv_max;
        for (unsigned i = 0; i < k; ++i, --k)
        {
            CheckSparseVectorRange(sv, i, k);
        }
    }
    
    cout << " ---------------- Sparse vector Range partitioning test  OK\n" << endl;
}

static
void CheckSparseVectorFilter(const sparse_vector_u32& sv, unsigned factor)
{
    sparse_vector_u32 sv1(sv);
    
    sparse_vector_u32::bvector_type bv_mask;
    for (unsigned i = 0; i < sv.size(); ++i)
    {
        if (i % factor == 0)
            bv_mask.set(i);
    }
    
    sv1.filter(bv_mask);
    
    for (unsigned i = 0; i < sv.size(); ++i)
    {
        unsigned v = sv.get(i);
        bool is_null = sv.is_null(i);
        unsigned v1 = sv1.get(i);
        bool is_null1 = sv1.is_null(i);

        if (i % factor == 0)
        {
            if (v != v1 || is_null != is_null1)
            {
                cerr << "Error! (1)sparse_vector<>::filter() failed at:" << i << endl;
                exit(1);
            }
        }
        else
        {
            if (v == v1 || is_null == is_null1)
            {
                cerr << "Error! (2)sparse_vector<>::filter() failed at:" << i << endl;
                exit(1);
            }
        }
    }
}

static
void TestSparseVectorFilter()
{
    cout << " ---------------- Sparse vector Filter test" << endl;
    cout << "Basic check" << endl;
    {
        sparse_vector_u32 sv(bm::use_null);
        sv.set(2, 25);
        sv.set(3, 35);
        sv.set(7, 75);
        sv.set(10, 2);
        sv.set(21, 201);
        
        sparse_vector_u32::bvector_type bv_mask { 2, 7 };
        
        sv.filter(bv_mask);
        
        for (unsigned i = 0; i < sv.size(); ++i)
        {
            unsigned v = sv.get(i);
            bool is_null = sv.is_null(i);
            if (i == 2 || i == 7)
            {
                assert(v != 0);
                assert(!is_null);
            }
            else
            {
                assert(v == 0);
                assert(is_null);
            }
        }
        
    }
    
    cout << "Stress check 1" << endl;

    {
        sparse_vector_u32 sv(bm::use_null);
        const unsigned sv_max = 250000;
        cout << "Filling the vector ... " << flush;
        for (unsigned i = 0; i < sv_max; ++i)
        {
            sv.push_back(i);
        }
        sv[0] = 113213;

        
        sparse_vector_u32::bvector_type bv_mask;
        for (unsigned i = 0; i < sv_max; ++i)
        {
            if (i % 2 == 0)
                bv_mask.set(i);
        }
        cout << "done." << endl;
        
        sv.filter(bv_mask);
        for (unsigned i = 0; i < sv_max; ++i)
        {
            unsigned v = sv.get(i);
            bool is_null = sv.is_null(i);
            if (i % 2 == 0)
            {
                assert(v == i || (i == 0 && v == 113213));
                assert(!is_null);
            }
            else
            {
                assert(v == 0);
                assert(is_null);
            }
        }
    }

    cout << "Stress check 2" << endl;

    {
        sparse_vector_u32 sv(bm::use_null);
        const unsigned sv_max = 250000;
        cout << "Filling the vector ... " << flush;
        for (unsigned i = 0; i < sv_max; ++i)
        {
            sv.push_back(i);
        }
        cout << "done" << endl;
        
        const unsigned max_factor = 10000;
        for (unsigned i = 2; i < max_factor; ++i)
        {
            CheckSparseVectorFilter(sv, i);
            if (i % 256 == 0)
                cout << "\r" << i << "/" << max_factor << flush;
        }
        cout << endl;
    }
    
    cout << " ---------------- Sparse vector Filter test OK" << endl;
}



static
void TestSparseVectorTransform()
{
    cout << " ---------------- Test set transformation with sparse vector" << endl;

    {
        sparse_vector_u32 sv(bm::use_null);
        bvect bv_in { 1, 2, 3, 10, 20 };
        bvect bv_out;
        
        bvector_transform_11(bv_in, sv, bv_out);
        assert(bv_out.count() == 0);
        cout << "Transform11 with empty sv - ok" << endl;

        bm::set2set_11_transform<sparse_vector_u32> set2set;
        unsigned to;
        bool found = set2set.remap(0, sv, to);
        assert(!found);
        found = set2set.remap(3, sv, to);
        assert(!found);

    }

    {
        sparse_vector_u32 sv(bm::use_null);

        sv.set(2, 25);
        sv.set(3, 35);
        sv.set(7, 75);
        sv.set(10, 2);
        sv.set(21, 201);

        bm::set2set_11_transform<sparse_vector_u32> set2set;
        unsigned to;
        bool found = set2set.remap(0, sv, to);
        assert(!found);
        found = set2set.remap(3, sv, to);
        assert(found);
        assert(to == 35);

        bvect bv_in { 1, 2, 3, 10, 20 };
        bvect bv_control {25, 35, 2 };

        {
            bvect bv_out;
            bvector_transform_11(bv_in, sv, bv_out);
            int cmp = bv_control.compare(bv_out);
            if (cmp != 0)
            {
                cerr << "Transform11 (1) control comparison failed" << endl;
                exit(1);
            }
            
            sv.optimize();
            bv_out.clear();
            
            bvector_transform_11(bv_in, sv, bv_out);
            cmp = bv_control.compare(bv_out);
            if (cmp != 0)
            {
                cerr << "Transform11 (1, 1) control comparison failed" << endl;
                exit(1);
            }
        }
        
        cout << "Transform11 (1) - ok" << endl;
    }
    {
        sparse_vector_u32 sv;

        sv.set(2, 25);
        sv.set(3, 35);
        sv.set(7, 75);
        sv.set(10, 2);
        sv.set(21, 201);

        bm::set2set_11_transform<sparse_vector_u32> set2set;
        unsigned to;
        bool found = set2set.remap(0, sv, to);
        assert(found);
        assert(to == 0);
        found = set2set.remap(8, sv, to);
        assert(found);
        assert(to == 0);
        found = set2set.remap(3, sv, to);
        assert(found);
        assert(to == 35);


        bvect bv_in { 0, 2, 3, 10};
        bvect bv_control {0, 25, 35, 2 };

        {
            bvect bv_out;
            bvector_transform_11(bv_in, sv, bv_out);
            int cmp = bv_control.compare(bv_out);
            if (cmp != 0)
            {
                cerr << "Transform11 (1-1) control comparison failed" << endl;
                exit(1);
            }
            
            sv.optimize();
            bv_out.clear();
            
            bvector_transform_11(bv_in, sv, bv_out);
            cmp = bv_control.compare(bv_out);
            if (cmp != 0)
            {
                cerr << "Transform11 (1-1, 1) control comparison failed" << endl;
                exit(1);
            }
        }
        
        cout << "Transform11 (1-1) - ok" << endl;
    }

    {
        bvect bv_in, bv_out;
        sparse_vector_u32 sv(bm::use_null);
        
        generate_bvector(bv_in);
        
        {
            bvect::enumerator en = bv_in.first();
            for (;en.valid(); ++en)
            {
                bm::id_t idx = *en;
                sv.set(idx, idx); // 1 to 1 direct
            }
        }
        bvector_transform_11(bv_in, sv, bv_out);
        int cmp = bv_in.compare(bv_out);
        if (cmp != 0)
        {
            cerr << "Transform11 (2) control comparison failed" << endl;
            exit(1);
        }
        
        sv.optimize();
        
        bvector_transform_11(bv_in, sv, bv_out);
        cmp = bv_in.compare(bv_out);
        if (cmp != 0)
        {
            cerr << "Transform11 (2, 2) control comparison failed" << endl;
            exit(1);
        }
        
        cout << "Transform11 (2) - ok" << endl;
    }

    {
        bvect bv_in, bv_out;
        sparse_vector_u32 sv(bm::use_null);
        
        generate_bvector(bv_in);
        
        bvect bv_control;
        {
            bvect::enumerator en = bv_in.first();
            for (;en.valid(); ++en)
            {
                bm::id_t idx = *en;
                bv_control.set(idx + 50000000);
            }
        }
        
        {
            bvect::enumerator en = bv_in.first();
            for (;en.valid(); ++en)
            {
                bm::id_t idx = *en;
                sv.set(idx, idx + 50000000); // 1 to 1 direct with a base shift
            }
        }
        bvector_transform_11(bv_in, sv, bv_out);
        
        int cmp = bv_control.compare(bv_out);
        if (cmp != 0)
        {
            cerr << "Transform11 (3) control comparison failed" << endl;
            exit(1);
        }
        
        sv.optimize();
        
        cmp = bv_control.compare(bv_out);
        if (cmp != 0)
        {
            cerr << "Transform11 (3, 2) control comparison failed" << endl;
            exit(1);
        }
        cout << "Transform11 (3) - ok" << endl;
    }


    {
        bvect bv_in, bv_out;
        sparse_vector_u32 sv(bm::use_null);
        
        generate_bvector(bv_in);
        
        bvect bv_control;
        bv_control.set(50000000);
        
        {
            bvect::enumerator en = bv_in.first();
            for (;en.valid(); ++en)
            {
                bm::id_t idx = *en;
                sv.set(idx, 50000000); // M:1
            }
        }
        bvector_transform_11(bv_in, sv, bv_out);
        
        int cmp = bv_control.compare(bv_out);
        if (cmp != 0)
        {
            cerr << "Transform11 (4) control comparison failed" << endl;
            exit(1);
        }
        
        sv.optimize();
        bvector_transform_11(bv_in, sv, bv_out);

        cmp = bv_control.compare(bv_out);
        if (cmp != 0)
        {
            cerr << "Transform11 (4, 2) control comparison failed" << endl;
            exit(1);
        }
        cout << "Transform11 (4) - ok" << endl;
    }
    
    
    {
        bvect bv_in, bv_out;
        sparse_vector_u32 sv(bm::use_null);
        
        generate_bvector(bv_in);
        bvect bv_control(bv_in);
        {
            bvect::enumerator en = bv_in.first();
            for (;en.valid(); ++en)
            {
                bm::id_t idx = *en;
                sv.set(idx, idx); // same:same
            }
        }
        bvector_transform_11(bv_in, sv, bv_out);
        
        int cmp = bv_control.compare(bv_out);
        if (cmp != 0)
        {
            cerr << "Transform11 (5) control comparison failed" << endl;
            exit(1);
        }
        
        sv.optimize();
        bvector_transform_11(bv_in, sv, bv_out);

        cmp = bv_control.compare(bv_out);
        if (cmp != 0)
        {
            cerr << "Transform11 (5, 2) control comparison failed" << endl;
            exit(1);
        }
        cout << "Transform11 (5) - ok" << endl;
    }


    cout << " --------------- Test set transformation with sparse vector OK" << endl;
}

static
void TestSparseVectorScan()
{
    cout << " --------------- Test sparse_vector<> scan algo" << endl;
    
    bm::sparse_vector_scanner<sparse_vector_u32> scanner;
    bm::sparse_vector_scanner<sparse_vector_u64> scanner_64;
    bm::sparse_vector_scanner<rsc_sparse_vector_u32> rsc_scanner;

    {
        sparse_vector_u32 sv(bm::use_null);
        bvect bv_control;
        scanner.find_eq(sv, 25, bv_control);
        assert(!bv_control.any());
        scanner.invert(sv, bv_control);
        assert(!bv_control.any());
    }

    {
        sparse_vector_u32 sv;
        bvect bv_control;
        for (unsigned i = 0; i < 20; ++i)
        {
            sv.set(i, 0);
        }
        scanner.find_eq(sv, 0, bv_control);
        unsigned found = bv_control.count();
        assert(found == 20);
        scanner.invert(sv, bv_control);
        found = bv_control.count();
        assert(!found);
    }

    {
        cout << endl << "Unique search check" << endl;
        sparse_vector_u32 sv;
        rsc_sparse_vector_u32 csv(bm::use_null);

        bvect bv_control, bv_control2;
        bvect::allocator_pool_type pool;
        bvect::mem_pool_guard(pool, bv_control);
        bvect::mem_pool_guard(pool, bv_control2);

        unsigned sv_size = 1256000;
        {
            sparse_vector_u32::back_insert_iterator bi(sv.get_back_inserter());
            for (unsigned j = 0; j < sv_size; ++j)
            {
                *bi = j;
            }
        }
        csv.load_from(sv);
        
        {
        chrono_taker ct("sparse_vector<> search");

            for (unsigned j = 0; j < sv_size; ++j)
            {
                scanner.find_eq(sv, j, bv_control);
                
                if (bv_control.count()!= 1)
                {
                    cerr << "1. Unique search discrepancy at value=" << j
                         << " count = " << bv_control.count() << endl;
                    exit(1);
                }
                unsigned v1, v2;
                bool b = bv_control.find_range(v1, v2);
                assert(b);
                if (v1 != v2)
                {
                    cerr << "2. Unique search discrepancy at value=" << j
                         << " count = " << bv_control.count() << endl;
                    exit(1);
                }
                
                bm::id_t pos;
                bool found = scanner.find_eq(sv, j, pos);
                if (!found)
                {
                    cerr << "3. Unique search failure at value=" << j
                         << endl;
                    exit(1);
                }
                if (v1 != pos)
                {
                    cerr << "4. Unique search discrepancy at value=" << j
                         << " found = " << pos << endl;
                    exit(1);
                }
                
                rsc_scanner.find_eq(csv, j, bv_control2);
                int res = bv_control.compare(bv_control2);
                if (res != 0)
                {
                    cerr << "RSC scan comparison failed at value =" << j
                    << endl;
                    exit(1);
                }


                if (j % 1000 == 0)
                    cout << "\r" << j << "/" << sv_size << "    " << flush;
            } // for
            cout << endl;
        }
 
        cout << "Unique search OK" << endl;
    }

    {
        cout << "Find EQ test on flat data" << endl;
        bvect::allocator_pool_type pool;
        unsigned max_value = 128000;
        for (unsigned value = 0; value < max_value; ++value)
        {
            sparse_vector_u32 sv;
            sparse_vector_u64 sv_64;
            rsc_sparse_vector_u32 csv;
            
            bvect bv_control, bv_control2;
            bvect::mem_pool_guard(pool, bv_control);

            unsigned sv_size = 67000;
            
            {
                sparse_vector_u32::back_insert_iterator bi(sv.get_back_inserter());
                sparse_vector_u64::back_insert_iterator bi_64(sv_64.get_back_inserter());
                for (unsigned j = 0; j < 67000; ++j)
                {
                    *bi = value;
                    bm::id64_t v64 = value;
                    v64 <<= 32;
                    *bi_64 = v64;
                }
            }
            csv.load_from(sv);
            
            scanner.find_eq(sv, value, bv_control);
            unsigned found = bv_control.count();
            
            if (found != sv_size)
            {
                cerr << "1. sparse_vector<>::find_eq() discrepancy for value=" << value
                     << " count = " << found << endl;
                exit(1);
            }
            
            rsc_scanner.find_eq(csv, value, bv_control2);
            int res = bv_control.compare(bv_control2);
            if (res != 0)
            {
                cerr << "RSC scan comparison failed at value =" << value
                << endl;
                exit(1);
            }

            
            {
                bm::id64_t v64 = value;
                v64 <<= 32;

                scanner_64.find_eq(sv_64, v64, bv_control);
                found = bv_control.count();
                
                if (found != sv_size)
                {
                    cerr << "1. (64) sparse_vector<>::find_eq() discrepancy for value=" << value
                         << " count = " << found << endl;
                    exit(1);
                }
            }
            
            // not found check
            scanner.find_eq(sv, value+1, bv_control);
            if (bv_control.any())
            {
                cerr << "1. sparse_vector<>::find_eq() (any) discrepancy for value=" << value+1
                     << " count = " << bv_control.count() << endl;
                exit(1);
            }
            rsc_scanner.find_eq(csv, value+1, bv_control2);
            res = bv_control.compare(bv_control2);
            if (res != 0)
            {
                cerr << "1. RSC scan comparison failed at value =" << value+1
                << endl;
                exit(1);
            }

            {
            BM_DECLARE_TEMP_BLOCK(tb)
            sv.optimize(tb);
            }
            
            bv_control.clear();
            
            scanner.find_eq(sv, value, bv_control);
            found = bv_control.count();
            
            if (found != sv_size)
            {
                cerr << "2. sparse_vector<>::find_eq() discrepancy for value=" << value
                     << " count = " << found << endl;
                exit(1);
            }
            
            // not found check
            scanner.find_eq(sv, value+1, bv_control);
            if (bv_control.any())
            {
                cerr << "2. sparse_vector<>::find_eq() (any) discrepancy for value=" << value+1
                     << " count = " << bv_control.count() << endl;
                exit(1);
            }

            
            if (value % 256 == 0)
                cout << "\r" << value << "/" << max_value << "    " << flush;
        }
        
        cout << endl << "Flat EQ ok" << endl;
    }
    

    
    cout << " \n--------------- Test sparse_vector<> scan algo OK" << endl;
}

static
void TestCompressedSparseVectorAlgo()
{
    cout << " --------------- TestCompressedSparseVectorAlgo()" << endl;

    {
        rsc_sparse_vector_u32 csv1(bm::use_null);
        rsc_sparse_vector_u32 csv2(bm::use_null);

        bm::sparse_vector<unsigned, bvect>::size_type pos;
        bool f;
        f = bm::sparse_vector_find_first_mismatch(csv1, csv2, pos);
        assert(!f);

        csv1.push_back(10, 10);
        csv1.push_back(11, 10);
        csv1.push_back(200, 0);
        csv1.push_back(300, 0);


        f = bm::sparse_vector_find_first_mismatch(csv1, csv2, pos);
        assert(f);
        assert(pos == 10);

        csv1.sync();
        csv2.sync();

        f = bm::sparse_vector_find_first_mismatch(csv1, csv2, pos);
        assert(f);
        assert(pos == 10);
    }

    {
        rsc_sparse_vector_u32 csv1(bm::use_null);
        rsc_sparse_vector_u32 csv2(bm::use_null);

        csv1.push_back(200, 0);
        csv1.push_back(300, 0);

        csv2.push_back(200, 0);
        csv2.push_back(300, 0);

        bm::sparse_vector<unsigned, bvect>::size_type pos;
        bool f;
        f = bm::sparse_vector_find_first_mismatch(csv1, csv2, pos);
        assert(!f);

        csv2.push_back(400, 0);
        f = bm::sparse_vector_find_first_mismatch(csv1, csv2, pos);
        assert(f);
        assert(pos == 400);
    }

    {
        rsc_sparse_vector_u32 csv1(bm::use_null);
        rsc_sparse_vector_u32 csv2(bm::use_null);

        bm::sparse_vector<unsigned, bvect>::size_type pos;
        bool f;

        csv1.push_back(10, 10);
        csv1.push_back(11, 10);
        csv1.push_back(200, 0);
        csv1.push_back(300, 0);

        csv2 = csv1;
        f = bm::sparse_vector_find_first_mismatch(csv1, csv2, pos);
        assert(!f);
        csv2.push_back(400, 256);

        f = bm::sparse_vector_find_first_mismatch(csv1, csv2, pos);
        assert(f);
        assert(pos == 400);

        csv1.optimize();
        csv2.optimize();

        csv1.sync();
        csv2.sync();

        f = bm::sparse_vector_find_first_mismatch(csv1, csv2, pos);
        assert(f);
        assert(pos == 400);

    }

    // ----------------------------------------------------------


    {
        cout << endl << "Unique mismatch check" << endl;
        sparse_vector_u32 sv1, sv2;
        rsc_sparse_vector_u32 csv1(bm::use_null);
        rsc_sparse_vector_u32 csv2(bm::use_null);


        unsigned sv_size = 525600;
        {
            sparse_vector_u32::back_insert_iterator bi(sv1.get_back_inserter());
            unsigned v = 0; unsigned cnt = 0;
            for (unsigned j = 0; j < sv_size; ++j)
            {
                *bi = v;
                if (++cnt > 256)
                {
                    cnt = 0; ++v;
                }
            }
        }
        csv1.load_from(sv1);

        sv2 = sv1;
        csv2 = csv1;

        bm::sparse_vector<unsigned, bvect>::size_type pos;
        bool f;

        f = bm::sparse_vector_find_first_mismatch(csv1, csv2, pos);
        assert(!f);
        f = bm::sparse_vector_find_first_mismatch(sv1, sv2, pos);
        assert(!f);

        for (unsigned k = 0; k < 4; ++k)
        {
            cout << "PASS = " << k << endl;
            chrono_taker ct("sparse_vector<> unique value mismatch search");

            for (sparse_vector_u32::size_type j = 0; j < sv_size; ++j)
            {
                std::chrono::time_point<std::chrono::steady_clock> st;
                st = std::chrono::steady_clock::now();

                sparse_vector_u32::value_type v2 = sv2[j];
                v2 = ~v2;
                sv2[j] = v2;
                f = bm::sparse_vector_find_first_mismatch(sv1, sv2, pos);
                assert(f);
                assert(pos == j);

                sv2[j] = ~v2; // restore
                f = bm::sparse_vector_find_first_mismatch(sv1, sv2, pos);
                assert(!f);

                v2 = csv2[j];
                v2 = ~v2;
                csv2.set(j, v2);
                f = bm::sparse_vector_find_first_mismatch(csv1, csv2, pos);
                assert(f);
                assert(pos == j);

                csv2.set(j, ~v2); // restore
                f = bm::sparse_vector_find_first_mismatch(csv1, csv2, pos);
                assert(!f);


                if (j % 10000 == 0)
                {
                    sv2.optimize();
                    csv2.optimize();
                    csv2.sync();

                    std::chrono::time_point<std::chrono::steady_clock> f1 = std::chrono::steady_clock::now();
                    auto diff = f1 - st;
                    //auto d = std::chrono::duration <double, std::milli> (diff).count();

                    cout << "\r" << j << "/" << sv_size << " " <<
                            " (" << diff.count() << ")" << flush;
                }
            } // for
            cout << endl;

            switch(k)
            {
            case 0:
                sv1.optimize();
                break;
            case 1:
                sv1.optimize();
                csv1.optimize();
                break;
            case 2:
                sv1.optimize();
                csv1.optimize();
                sv2.optimize();
                break;
            case 3:
                sv1.optimize();
                csv1.optimize();
                sv2.optimize();
                csv2.optimize();
                break;
            default:
                assert(0);
            }
        } // for k

        cout << "Unique search OK" << endl;
    }



    cout << " --------------- TestCompressedSparseVectorAlgo() OK" << endl;
}


static
void TestCompressedSparseVectorScan()
{
    cout << " --------------- Test rsc_sparse_vector<> scan algo" << endl;
    
    bm::sparse_vector_scanner<rsc_sparse_vector_u32> scanner;
    
    {
        rsc_sparse_vector_u32 csv(bm::use_null);
        bvect bv_control;
        scanner.find_eq(csv, 25, bv_control);
        assert(!bv_control.any());
        scanner.invert(csv, bv_control);
        assert(!bv_control.any());
    }
    
    {
        rsc_sparse_vector_u32 csv(bm::use_null);
        bvect bv_res;

        csv.push_back(10, 10);
        csv.push_back(11, 10);
        csv.push_back(200, 0);
        csv.push_back(300, 0);
        
        csv.sync();

        bm::id_t idx;
        for (unsigned i = 0; i < 2; ++i)
        {
            bool found = scanner.find_eq(csv, 10, idx);
            assert(found);
            assert(idx == 10);
            
            scanner.find_eq(csv, 10, bv_res);
            assert(bv_res.count()==2);
            assert(bv_res.test(10));
            assert(bv_res.test(11));
            
            scanner.find_eq(csv, 0, bv_res);
            assert(bv_res.count()==2);
            assert(bv_res.test(200));
            assert(bv_res.test(300));

            found = scanner.find_eq(csv, 0, idx);
            assert(found);
            assert(idx == 200);
            
            csv.optimize();
        } // for
    }

    cout << " --------------- Test rsc_sparse_vector<> scan algo OK" << endl;

}


// fill pseudo-random plato pattern into two vectors
//
template<class SV>
void FillSparseIntervals(std::vector<unsigned>&   vect,
                         SV& svect,
                         unsigned min,
                         unsigned max,
                         unsigned fill_factor)
{
    unsigned diap = max - min;

    unsigned count;


    switch (fill_factor)
    {
    case 0:
        count = diap / 1000;
        break;
    case 1:
        count = diap / 100;
        break;
    default:
        count = diap / 10;
        break;

    }
    
    if (vect.size() < max)
    {
        vect.resize(max + 1);
    }
    if (svect.size() < max)
    {
        svect.resize(max + 1);
    }
    
    unsigned val = 0;
    
    for ( ;min < max; )
    {
        // hi-band interval
        val = rand() % (65535 * 2);
        unsigned i;
        for (i = 0; i < count; ++i)
        {
            vect[min] = val;
            svect.set(min, val);
            ++min;
            if (min > max)
                break;
        } // for i
        
        // gap with all zeroes
        unsigned inc = rand() % 2048;
        min += inc;
        if (min > max)
            break;

        // low band plato
        val = rand() % 8;
        for (i = 0; i < count; ++i)
        {
            vect[min] = val;
            svect.set(min, val);
            ++min;
            if (min > max)
                break;
        } // for i
        
    } // for min
}

static
void TestSparseVector_Stress(unsigned count)
{

    cout << "---------------------------- Bit-plain sparse vector stress" << endl;

    cout << "Interval shift check.\n";
    // interval shift check
    bool detailed_check=true;
    for (unsigned i = 0; i < count; ++i)
    {
        unsigned fill_factor = 0;
        for (unsigned min = 0; min < 10000000; min+= rand()%100000)
        {
            unsigned max = min + (65535 * 10);
            {{
                bm::null_support null_able =
                                (min % 2 == 0) ? bm::no_null : bm::use_null;
                
                std::vector<unsigned> vect;
                bm::sparse_vector<unsigned, bvect > sv(null_able);
            
                FillSparseIntervals(vect, sv, min, max, fill_factor);
                
                bool res = CompareSparseVector(sv, vect, true, detailed_check);
                if (!res)
                {
                    cerr << "sparse-dense vector comparison failed" << endl;
                    exit(1);
                }
                
                sv.optimize();
                res = CompareSparseVector(sv, vect, true);
                if (!res)
                {
                    cerr << "sparse-dense vector comparison failed" << endl;
                    exit(1);
                }
                cout << "\r min=" << min << " max=" << max << " ff= "
                     << fill_factor << flush;
            }}
            if (++fill_factor > 2) fill_factor = 0;
        } // for min
        cout << "\n." << flush;
        detailed_check = false; // all other passes go faster
    } // for i
    cout << endl;

    cout << "--------------------------- Interval shift check Ok" << endl;

    cout << "Join check" << endl;
    for (unsigned i = 0; i < 1; ++i)
    {
        unsigned fill_factor = 0;
        for (unsigned min = 0; min < 10000000; min+= rand()%100000)
        {
            unsigned max = min + (65535 * 10);
            unsigned min2 = max + rand() % 65536;
            unsigned max2 = min2 + (65535 * 10);

            bm::null_support null_able1 =
                            (min % 2 == 0) ? bm::no_null : bm::use_null;

            {{
                std::vector<unsigned> vect1;
                std::vector<unsigned> vect2;
                bm::sparse_vector<unsigned, bvect > sv1(null_able1);
                bm::sparse_vector<unsigned, bvect > sv2(null_able1);
            
                FillSparseIntervals(vect1, sv1, min, max, fill_factor);
                FillSparseIntervals(vect2, sv2, min2, max2, fill_factor);
                
                if (rand()%2)
                {
                    sv1.optimize();
                }
                if (rand()%3 == 0)
                {
                    sv2.optimize();
                }
                bm::sparse_vector<unsigned, bvect > sv3;
                bm::sparse_vector<unsigned, bvect > sv4(sv2);
                
                sv1.join(sv2);
                sv3.join(sv1);
                sv3.join(sv2);
                sv4.join(sv1);
                
                if (sv1.size() != sv3.size() || sv2.size() != sv3.size() || sv3.size() != sv4.size())
                {
                    cerr << "Sparse join size error:" << sv1.size() << endl;
                    exit(1);
                }
                
                
                for ( i = 0; i < sv1.size(); ++i)
                {
                    unsigned v1 = sv1[i];
                    unsigned v3 = sv3[i];
                    unsigned v4 = sv4[i];
                    if (v1 != v3 || v1 != v4)
                    {
                        cerr << "Sparse join cmp failed:" << v1 << "!=" << v3 << "!=" << v4 << endl;
                        exit(1);
                    }
                } // for i
                
                bool b1 = TestEqualSparseVectors(sv1, sv3, false);
                if (!b1)
                {
                    cerr << "Equal 1 comparison failed" << endl;
                    exit(1);
                }
                bool b2 = TestEqualSparseVectors(sv1, sv4, false);
                if (!b2)
                {
                    cerr << "Equal 2 comparison failed" << endl;
                    exit(1);
                }
                
                bool b3 = TestEqualSparseVectors(sv3, sv4, false);
                if (!b3)
                {
                    cerr << "Equal 3 comparison failed" << endl;
                    exit(1);
                }
                
                
                
                cout << "+" << flush;
                
            }}
            ++fill_factor;
            if (fill_factor > 2) fill_factor = 0;
        } // for min
        
        cout << "." << flush;
        
    } // for i
    
    cout << "--------------------------- Join check Ok" << endl;
    
    
    
    cout << "---------------------------- Bit-plain sparse vector stress OK" << endl;
}

static
void TestStrSparseVector()
{
   cout << "---------------------------- Bit-plain STR sparse vector test" << endl;

   typedef str_sparse_vector<char, bvect, 32> str_svect_type;

   {
   str_sparse_vector<char, bvect, 32> str_sv0;
   str_sparse_vector<char, bvect, 32> str_sv1(str_sv0);
   str_sparse_vector<char, bvect, 32> str_sv2;
   str_sv2 = str_sv1;
   
   assert(str_sv1.size() == 0);
   str_sparse_vector<char, bvect, 32> str_sv3(std::move(str_sv0));
   assert(!str_sv3.is_remap());
   }
   
   {
       const char* s0 = "AbC";
       const char* s1 = "jKl";
       char str[256];
       str_sparse_vector<char, bvect, 32> str_sv0;
       int cmp;

        assert(str_sv0.size()==0);
        str_sv0.set(0, s0);
        str_sv0.get(0, str, sizeof(str));
        cmp = ::strcmp(str, s0);
        assert(cmp == 0);
        assert(str_sv0.size()==1);

       str_sv0.set(1, s1);
       str_sv0.get(0, str, sizeof(str));
       cmp = ::strcmp(str, s0);
       assert(cmp == 0);
       
       str_sv0.get(1, str, sizeof(str));
       cmp = ::strcmp(str, s1);
       assert(cmp == 0);
       
       str_sv0.optimize();

       str_sv0.get(0, str, sizeof(str));
       cmp = ::strcmp(str, s0);
       assert(cmp == 0);
       
       str_sv0.get(1, str, sizeof(str));
       cmp = ::strcmp(str, s1);
       assert(cmp == 0);
       
       string str0 = "AtGc";
       str_sv0.assign(3, str0);
       str_sv0.get(3, str, sizeof(str));
       cmp = ::strcmp(str, str0.c_str());
       assert(cmp == 0);
       //auto sz=str_sv0.size();
       assert(str_sv0.size()==4);

       string str1;
       str_sv0.get(3, str1);
       assert(str0 == str1);
       {
           str0 = "TTF";
           str_sv0.assign(3, str0);
           str_sv0.get(3, str, sizeof(str));
           cmp = ::strcmp(str, str0.c_str());
           assert(cmp==0);

           str0.clear();
           str_sv0.assign(3, str0);
           str_sv0.get(3, str, sizeof(str));
           cmp = ::strcmp(str, str0.c_str());
           assert(cmp==0);
       }
       
       // test truncation of input string
       {
          str_sparse_vector<char, bvect, 3> str_sv10;
          const char* s10 = "12345";
          const char* s10c = "123";
           str_sv10.set(1, s10);
           str_sv10.get(1, str, sizeof(str));
           cmp = ::strcmp(str, s10c);
           assert(cmp == 0);
       }
       
       // test string insert
       {
          str_sparse_vector<char, bvect, 3> str_sv10;
          const char* cs0 = "10";
          const char* cs2 = "30";
          const char* cs1 = "200";
          
          str_sv10.push_back(cs0);
          str_sv10.push_back(cs1);
          str_sv10.insert(1, cs2);

           str_sv10.get(0, str, sizeof(str));
           cmp = ::strcmp(str, cs0);
           assert(cmp == 0);
           
           str_sv10.get(1, str, sizeof(str));
           cmp = ::strcmp(str, cs2);
           assert(cmp == 0);

           str_sv10.get(2, str, sizeof(str));
           cmp = ::strcmp(str, cs1);
           assert(cmp == 0);
           
           str_sv10.clear();
           assert(str_sv10.size() == 0);
       }
       
       // test erase
       {
          str_sparse_vector<char, bvect, 3> str_sv10;
          const char* cs0 = "10";
          const char* cs1 = "200";
          const char* cs2 = "30";
          
          str_sv10.push_back(cs0);
          str_sv10.push_back(cs1);
          str_sv10.push_back(cs2);
          
          str_sv10.erase(1);
          assert(str_sv10.size() == 2);

           str_sv10.get(0, str, sizeof(str));
           cmp = ::strcmp(str, cs0);
           assert(cmp == 0);
           str_sv10.get(1, str, sizeof(str));
           cmp = ::strcmp(str, cs2);
           assert(cmp == 0);

          str_sv10.erase(0);
          assert(str_sv10.size() == 1);
       }
       
       // test decode
       {
          str_sparse_vector<char, bvect, 3> str_sv10;
          const char* cs0 = "10";
          const char* cs1 = "200";
          const char* cs2 = "30";
          str_sv10.push_back(cs0);
          str_sv10.push_back(cs1);
          str_sv10.push_back(cs2);
          
          bm::heap_matrix<char, 1024, 64, bvect::allocator_type> hmatr(true);
          
          unsigned d = 0;
          char *s;
          
          d = str_sv10.decode(hmatr, 0, 1);
          s = hmatr.row(0);
          cmp = ::strcmp(s, cs0);
          assert(cmp == 0);
          assert(d == 1);

          d = str_sv10.decode(hmatr, 1, 1);
          s = hmatr.row(0);
          cmp = ::strcmp(s, cs1);
          assert(cmp == 0);
          assert(d == 1);

          d = str_sv10.decode(hmatr, 2, 1);
          s = hmatr.row(0);
          cmp = ::strcmp(s, cs2);
          assert(cmp == 0);
          assert(d == 1);
          
          // decode beyond limit
          d = str_sv10.decode(hmatr, 3, 1);
          s = hmatr.row(0);
          assert(*s == 0);
          assert(d == 0);
          
          d = str_sv10.decode(hmatr, 0, 100000);
          assert(d == str_sv10.size());
          s = hmatr.row(0);
          cmp = ::strcmp(s, cs0);
          assert(cmp == 0);
          s = hmatr.row(1);
          cmp = ::strcmp(s, cs1);
          assert(cmp == 0);
          s = hmatr.row(2);
          cmp = ::strcmp(s, cs2);
          assert(cmp == 0);

          d = str_sv10.decode(hmatr, 1, 100000);
          assert(d == str_sv10.size()-1);
          s = hmatr.row(0);
          cmp = ::strcmp(s, cs1);
          assert(cmp == 0);
          s = hmatr.row(1);
          cmp = ::strcmp(s, cs2);
          assert(cmp == 0);

       }
       
       // test import
       {
          str_sparse_vector<char, bvect, 3> str_sv10;

          const char* cs0 = "@";
          const char* cs1 = " 2";
          const char* cs2 = "034";
          
          bm::heap_matrix<char, 1024, 64, bvect::allocator_type> hmatr(true);
          ::strncpy(hmatr.row(0), cs0, hmatr.cols());
          ::strncpy(hmatr.row(1), cs1, hmatr.cols());
          ::strncpy(hmatr.row(2), cs2, hmatr.cols());
          
          for (unsigned i = 0; i < 3; ++i)
          {
            const char* s = hmatr.row(i);
            cout << s << endl;
          }
          
          str_sv10.import(hmatr, 0, 3);
          assert(str_sv10.size() == 3);
          
           str_sv10.get(0, str, sizeof(str));
           cmp = ::strcmp(str, cs0);
           assert(cmp == 0);
           
           str_sv10.get(1, str, sizeof(str));
           cmp = ::strcmp(str, cs1);
           assert(cmp == 0);

           str_sv10.get(2, str, sizeof(str));
           cmp = ::strcmp(str, cs2);
           assert(cmp == 0);
       }


       // reference test / serialization test
       {
       const char* s = str_sv0[3];
       cmp = ::strcmp(s, str0.c_str());
       assert(cmp == 0);
       str_sv0[3] = "333";
       str_sv0.get(3, str, sizeof(str));

       s = str_sv0[3];
       cmp = ::strcmp(s, "333");
       assert(cmp == 0);
       
       {
           const str_sparse_vector<char, bvect, 32>& ssv = str_sv0;
           str_sparse_vector<char, bvect, 32>::const_reference ref3 = ssv[3];
           s = ref3;
           cmp = ::strcmp(s, "333");
           assert(cmp == 0);
       }
       
        BM_DECLARE_TEMP_BLOCK(tb)
        sparse_vector_serial_layout<str_svect_type> sv_lay;
        bm::sparse_vector_serialize<str_svect_type>(str_sv0, sv_lay, tb);

        str_sparse_vector<char, bvect, 32> str_sv2;
        
        const unsigned char* buf = sv_lay.buf();
        int res = bm::sparse_vector_deserialize(str_sv2, buf, tb);
        if (res != 0)
        {
            cerr << "De-Serialization error" << endl;
            exit(1);
        }
        
        bool eq = str_sv0.equal(str_sv2);
        assert(eq);
        
        str_sparse_vector<char, bvect, 64> str_sv3;   // size increase test
        buf = sv_lay.buf();
        res = bm::sparse_vector_deserialize(str_sv3, buf, tb);
        if (res != 0)
        {
            cerr << "De-Serialization error" << endl;
            exit(1);
        }


       }

   }
    
   {
       str_sparse_vector<char, bvect, 32> str_sv0;
       unsigned str_max = str_sv0.effective_max_str();
       assert(str_max == 0);
       str_sv0[0] = "1";
       str_max = str_sv0.effective_max_str();
       assert(str_max == 1);
       str_sv0[1] = "11";
       str_max = str_sv0.effective_max_str();
       assert(str_max == 2);
       str_sv0[2] = "123";
       str_max = str_sv0.effective_max_str();
       assert(str_max == 3);
       
       str_sv0.clear_range(1, 234567);

       char str[256];
       str_sv0.get(1, str, sizeof(str));
       assert(str[0]==0);
       str_sv0.get(2, str, sizeof(str));
       assert(str[0]==0);
   }
   
   {
       str_sparse_vector<char, bvect, 32> str_sv0;
       str_sv0[0] = "1";
       str_sv0[1] = "11";
       str_sv0[2] = "123";
       
       unsigned pos;
       bm::sparse_vector_scanner<bm::str_sparse_vector<char, bvect, 32> > scanner;
       
       bool found = scanner.find_eq_str(str_sv0, "1", pos);
       assert(found);
       assert(pos == 0);

       found = scanner.find_eq_str(str_sv0, "11", pos);
       assert(found);
       assert(pos == 1);

       found = scanner.find_eq_str(str_sv0, "123", pos);
       assert(found);
       assert(pos == 2);

       found = scanner.find_eq_str(str_sv0, "1234", pos);
       assert(!found);

       found = scanner.find_eq_str(str_sv0, "", pos);
       assert(!found);
   }

   // test basic remappings functions
   {
       str_sparse_vector<char, bvect, 32> str_sv0;
       str_sv0[0] = "1";
       str_sv0[1] = "11";
       str_sv0[2] = "123";
       
       str_sparse_vector<char, bvect, 32>::plain_octet_matrix_type occ_matrix;
       str_sparse_vector<char, bvect, 32>::plain_octet_matrix_type remap_matrix1;
       str_sparse_vector<char, bvect, 32>::plain_octet_matrix_type remap_matrix2;

       str_sv0.calc_octet_stat(occ_matrix);
       str_sv0.build_octet_remap(remap_matrix1, remap_matrix2, occ_matrix);
       
       bool res;
       int cmp;
       char str0[64];
       char str1[64];
       
       res = str_sv0.remap_tosv(&str0[0], 64, "1", remap_matrix2);
       assert(res);
       
       res = str_sv0.remap_fromsv(&str1[0], 64, &str0[0], remap_matrix1);
       assert(res);
       cmp = str_sv0.compare(0, &str1[0]);
       assert(cmp == 0);

       res = str_sv0.remap_tosv(&str0[0], 64, "2", remap_matrix2); // impossible case
       assert(!res);
       
       
       res = str_sv0.remap_tosv(&str0[0], 64, "11", remap_matrix2);
       assert(res);
       
       res = str_sv0.remap_fromsv(&str1[0], 64, &str0[0], remap_matrix1);
       assert(res);
       cmp = str_sv0.compare(1, &str1[0]);
       assert(cmp == 0);

       res = str_sv0.remap_tosv(&str0[0], 64, "123", remap_matrix2);
       assert(res);
       
       res = str_sv0.remap_fromsv(&str1[0], 64, &str0[0], remap_matrix1);
       assert(res);
       cout << str1 << endl;
       cmp = str_sv0.compare(2, &str1[0]);
       assert(cmp == 0);


       res = str_sv0.remap_tosv(&str0[0], 64, "133", remap_matrix2); // impossible case
       assert(!res);
       res = str_sv0.remap_tosv(&str0[0], 64, "1231", remap_matrix2); // impossible case
       assert(!res);
   }
   
   // build-vefify remapped sparse-vector
   {
       str_sparse_vector<char, bvect, 32> str_sv0;
       str_sparse_vector<char, bvect, 32> str_sv1;
       
       str_sv1.remap_from(str_sv0);
       assert(!str_sv1.is_remap());
       
       str_sv0[0] = "1";
       str_sv0[1] = "11";
       str_sv0[2] = "123";

       assert(!str_sv1.is_remap());
       
       str_sv1.remap_from(str_sv0);
       str_sv1.recalc_remap_matrix2();
       
       assert(str_sv1.is_remap());
       assert(str_sv1.size() == str_sv0.size());

       char str[256];
       int cmp;

        str_sv1.get(0, str, sizeof(str));
        cmp = ::strcmp(str, "1");
        assert(cmp==0);
        cmp = str_sv1.compare(0, "1");
        assert(cmp==0);
       
        bm::heap_matrix<char, 1024, 64, bvect::allocator_type> hmatr(true);

        // test remap decoder
        {
          unsigned d = 0;
          char *s;
          
          d = str_sv1.decode(hmatr, 0, 1);
          s = hmatr.row(0);
          cmp = ::strcmp(s, "1");
          assert(cmp == 0);
          assert(d == 1);
        }
       
        str_sv1.get(1, str, sizeof(str));
        cmp = ::strcmp(str, "11");
        assert(cmp==0);
        cmp = str_sv1.compare(1, "11");
        assert(cmp==0);

        str_sv1.get(2, str, sizeof(str));
        cmp = ::strcmp(str, "123");
        assert(cmp==0);

        // test remap decoder
        {
          unsigned d = 0;
          char *s;
          
          d = str_sv1.decode(hmatr, 2, 1);
          s = hmatr.row(0);
          cmp = ::strcmp(s, "123");
          assert(cmp == 0);
          assert(d == 1);
        }

        string s;
        str_sv1.get(2, s);
        cmp = ::strcmp(s.c_str(), "123");
        assert(cmp==0);
       
        {
           string s0 = "113";
           str_sv1.assign(4, s0);
           str_sv1.get(4, s);
           assert(s == s0);
           cout << s << endl;
           cmp = str_sv1.compare(4, "113");
           assert(cmp==0);
        }
       
       {
            bool equal = str_sv1.equal(str_sv0);
            assert(!equal);
           
            str_sparse_vector<char, bvect, 32> str_sv2(str_sv1);
            equal = str_sv1.equal(str_sv2);
            assert(equal);
       }
       {
            str_sparse_vector<char, bvect, 32> str_sv2(str_sv0);
            str_sv2 = str_sv1;
            bool equal = str_sv1.equal(str_sv2);
            assert(equal);
       }
   }

   // scanner search on remapped str vector
   {
       str_sparse_vector<char, bvect, 32> str_sv0;
       str_sparse_vector<char, bvect, 32> str_sv1;
       str_sv0[0] = "1";
       str_sv0[1] = "11";
       str_sv0[2] = "123";

       str_sv1.remap_from(str_sv0);
       
       unsigned pos, pos1;
       bm::sparse_vector_scanner<bm::str_sparse_vector<char, bvect, 32> > scanner;
       bm::sparse_vector_scanner<bm::str_sparse_vector<char, bvect, 32> > scanner1;
       scanner.bind(str_sv1, true);

       
       bool found = scanner.find_eq_str("1", pos);
       assert(found);
       assert(pos == 0);
       found = scanner1.lower_bound_str(str_sv1,"1", pos1);
       assert(found);
       assert(pos == pos1);

       found = scanner.find_eq_str("11", pos);
       assert(found);
       assert(pos == 1);
       found = scanner.bfind_eq_str("11", pos);
       assert(found);
       assert(pos == 1);
       found = scanner1.lower_bound_str(str_sv1, "11", pos1);
       assert(found);
       assert(pos == pos1);

       found = scanner.find_eq_str("123", pos);
       assert(found);
       assert(pos == 2);
       found = scanner.bfind_eq_str("123", pos);
       assert(found);
       assert(pos == 2);
       found = scanner1.lower_bound_str(str_sv1, "123", pos1);
       assert(found);
       assert(pos == pos1);


       found = scanner.find_eq_str("1234", pos);
       assert(!found);
       found = scanner.bfind_eq_str("1234", pos);
       assert(!found);
       // commented out because lower_bound throws exceptions on impossible strings
       /*
       found = scanner1.lower_bound_str(str_sv1,"1234", pos1);
       assert(!found);
       assert(pos1 == str_sv1.size());
       */

       found = scanner.find_eq_str("", pos);
       assert(!found);
   }

    // serialization of remap string vector
    {
       str_sparse_vector<char, bvect, 32> str_sv0;
       str_sparse_vector<char, bvect, 32> str_sv1;
       str_sparse_vector<char, bvect, 32> str_sv2;
       str_sv0[0] = "1";
       str_sv0[1] = "11";
       str_sv0[2] = "123";

       str_sv1.remap_from(str_sv0);

        BM_DECLARE_TEMP_BLOCK(tb)
        sparse_vector_serial_layout<str_svect_type> sv_lay;
        bm::sparse_vector_serialize<str_svect_type>(str_sv1, sv_lay, tb);

        const unsigned char* buf = sv_lay.buf();
        int res = bm::sparse_vector_deserialize(str_sv2, buf, tb);
        if (res != 0)
        {
            cerr << "De-Serialization error!" << endl;
            exit(1);
        }
        
        bool equal = str_sv1.equal(str_sv2);
        assert(equal);
   }
   
   
    // back insert iterator
    //
    {
       str_sparse_vector<char, bvect, 32> str_sv0;
       str_sparse_vector<char, bvect, 32>::back_insert_iterator bit;
       str_sparse_vector<char, bvect, 32>::back_insert_iterator bit2(&str_sv0);
       str_sparse_vector<char, bvect, 32>::back_insert_iterator bit3(bit);
       bit = bit3;
    }
    
    {
       str_sparse_vector<char, bvect, 32> str_sv0(bm::use_null);
       auto bi = str_sv0.get_back_inserter();
       bi = (const char*)0;
       bool b = str_sv0.is_null(0);
       assert(b);
       bi = (const char*)nullptr;
       bi = "123";
       bi.add_null();
       bi.flush();
       
       
       assert(str_sv0.size() == 4);
       b = str_sv0.is_null(0);
       assert(b);
       assert(str_sv0[0].is_null());
       assert(str_sv0.is_null(1));
       assert(!str_sv0.is_null(2));
       
       auto sz = str_sv0.size();
       str_sv0.set_null(sz);
       assert(sz + 1 == str_sv0.size() );
       assert(str_sv0.is_null(sz));
       
    }
    
    // test automatic optimization with back_insert iterator
    {
       str_sparse_vector<char, bvect, 32> str_sv0;
       {
           str_sparse_vector<char, bvect, 32>::back_insert_iterator
                                        bi = str_sv0.get_back_inserter();
           for (unsigned i = 0; i < 65536; ++i)
           {
                bi = "123";
           } // for
       }
       
       str_sparse_vector<char, bvect, 32>::statistics st;
       
       str_sv0.calc_stat(&st);
       assert(st.bit_blocks == 0);
       assert(st.gap_blocks == 0);
       
       {
           str_sparse_vector<char, bvect, 32>::back_insert_iterator
                                        bi = str_sv0.get_back_inserter();
           for (unsigned i = 0; i < 65536; ++i)
           {
                bi = "123";
           }
       }
       
       str_sv0.calc_stat(&st);
       assert(st.bit_blocks == 0);
       assert(st.gap_blocks == 0);
    }
    
    

    // const_iterator / back_inserter
    //
    {
       str_sparse_vector<char, bvect, 32> str_sv0;
       {
        str_sparse_vector<char, bvect, 32>::const_iterator it2(&str_sv0);
        assert(!it2.valid());
        str_sparse_vector<char, bvect, 32>::const_iterator it2c(it2);
        assert(!it2c.valid());
       }
       {
           str_sparse_vector<char, bvect, 32>::back_insert_iterator bi = str_sv0.get_back_inserter();
           bi = "1";
           bi = "11";
           bi = "123";
           
           bi.flush();
       }
    
        {
        str_sparse_vector<char, bvect, 32>::const_iterator it_end;
        assert(!it_end.valid());
        str_sparse_vector<char, bvect, 32>::const_iterator it2(&str_sv0);
        assert(it2.valid());
        const char* s = it2.value();
        int cmp = ::strcmp(s, "1");
        assert(cmp == 0);
        assert(!it2.is_null());

        str_sparse_vector<char, bvect, 32>::const_iterator it3(&str_sv0, 1);
        assert(it3.valid());
        s = it3.value();
        cmp = ::strcmp(s, "11");
        assert(cmp == 0);
        
        str_sparse_vector<char, bvect, 32>::const_iterator it4(&str_sv0, 3);
        assert(!it4.valid());
        it4.go_to(2);
        assert(it4.valid());
        s = it4.value();
        cmp = ::strcmp(s, "123");
        assert(cmp == 0);
        }
    }
    
    
   cout << "---------------------------- Bit-plain STR sparse vector test OK" << endl;
}

static
void TestStrSparseVectorAlgo()
{
    cout << "------------------------------ TestStrSparseVectorAlgo()" << endl;

    {
       str_sparse_vector<char, bvect, 32> str_sv1;
       str_sparse_vector<char, bvect, 32> str_sv2;

       {
           str_sparse_vector<char, bvect, 32>::back_insert_iterator bi = str_sv1.get_back_inserter();
           bi = "123";
           bi = "123";
           bi = "123";

           bi.flush();
       }
       str_sv2 = str_sv1;

       bm::sparse_vector<unsigned, bvect>::size_type pos;
       bool f;
       f = bm::sparse_vector_find_first_mismatch(str_sv1, str_sv2, pos);
       assert(!f);

       str_sv2.push_back("8");
       f = bm::sparse_vector_find_first_mismatch(str_sv1, str_sv2, pos);
       assert(f);
       assert(pos == 3);
       f = bm::sparse_vector_find_first_mismatch(str_sv2, str_sv1, pos);
       assert(f);
       assert(pos == 3);
       str_sv1.optimize();
       f = bm::sparse_vector_find_first_mismatch(str_sv2, str_sv1, pos);
       assert(f);
       assert(pos == 3);
       str_sv2.optimize();
       f = bm::sparse_vector_find_first_mismatch(str_sv2, str_sv1, pos);
       assert(f);
       assert(pos == 3);
    }

    {
       str_sparse_vector<char, bvect, 32> str_sv1(bm::use_null);
       str_sparse_vector<char, bvect, 32> str_sv2(bm::use_null);

       bm::sparse_vector<unsigned, bvect>::size_type pos;
       bool f;
       f = bm::sparse_vector_find_first_mismatch(str_sv1, str_sv2, pos);
       assert(!f);

       str_sv1[1000] = "123";
       str_sv1[10000] = "123";
       str_sv1[100000] = "123";
       str_sv1[1000000] = "123";

       str_sv2 = str_sv1;
       f = bm::sparse_vector_find_first_mismatch(str_sv1, str_sv2, pos);
       assert(!f);
       str_sv1.optimize();
       str_sv2.optimize();

       f = bm::sparse_vector_find_first_mismatch(str_sv1, str_sv2, pos);
       assert(!f);

       str_sv1[10000000] = "9";
       f = bm::sparse_vector_find_first_mismatch(str_sv1, str_sv2, pos);
       assert(f);
       assert(pos == 10000000);

       str_sv1[0] = "A";
       f = bm::sparse_vector_find_first_mismatch(str_sv1, str_sv2, pos);
       assert(f);
       assert(pos == 0);

    }

    cout << "------------------------------ TestStrSparseVectorAlgo() OK" << endl;
}

static
void TestStrSparseVectorSerial()
{
   cout << "---------------------------- TestStrSparseVectorSerial()" << endl;

    {
       str_sparse_vector<char, bvect, 32> str_sv1;
       str_sparse_vector<char, bvect, 32> str_sv2;

       {
           str_sparse_vector<char, bvect, 32>::back_insert_iterator bi = str_sv1.get_back_inserter();
           bi = "1";
           bi = "11";
           bi = "123";

           bi.flush();
       }
        BM_DECLARE_TEMP_BLOCK(tb)
        sparse_vector_serial_layout<str_sparse_vector<char, bvect, 32> > sv_lay;
        bm::sparse_vector_serialize<str_sparse_vector<char, bvect, 32> >(str_sv1, sv_lay, tb);
        const unsigned char* buf = sv_lay.buf();

        sparse_vector_u32::bvector_type bv_mask;
        bv_mask.set(1);
        bv_mask.set(2);
        bv_mask.set(100);
        bm::sparse_vector_deserializer<str_sparse_vector<char, bvect, 32> > sv_deserial;
        sv_deserial.deserialize(str_sv2, buf, bv_mask);

        assert(str_sv1.size() == str_sv2.size());
        char str[256];
        int cmp;

        str_sv2.get(1, str, sizeof(str));
        cmp = ::strcmp(str, "11");
        assert(cmp==0);
        cmp = str_sv2.compare(1, "11");
        assert(cmp==0);
        cmp = str_sv2.compare(2, "123");
        assert(cmp==0);
        cmp = str_sv2.compare(0, "");
        assert(cmp==0);
    }

    cout << "Stress deserialization (AND mask)" << endl;
    {
       str_sparse_vector<char, bvect, 32> str_sv0;
       str_sparse_vector<char, bvect, 32> str_sv1;
       str_sparse_vector<char, bvect, 32> str_sv2;

       {
           str_sparse_vector<char, bvect, 32>::back_insert_iterator bi = str_sv0.get_back_inserter();
           for (unsigned i = 0; i < 100000; ++i)
           {
                bi = "ATGC";
                bi = "GCTA";
                bi = "GCAA";
                bi = "TATA";
           }
       }

       str_sv1.remap_from(str_sv0);
       assert(str_sv1.is_remap());
       str_sv1.optimize();

        BM_DECLARE_TEMP_BLOCK(tb)
        sparse_vector_serial_layout<str_sparse_vector<char, bvect, 32> > sv_lay;
        bm::sparse_vector_serialize<str_sparse_vector<char, bvect, 32> >(str_sv1, sv_lay, tb);
        const unsigned char* buf = sv_lay.buf();

        bm::sparse_vector_deserializer<str_sparse_vector<char, bvect, 32> > sv_deserial;
        char s1[256];
        char s2[256];
        int cmp;

        bvect::size_type from = 100000-1;
        bvect::size_type to = from + 65536;
        for (unsigned i = from; i < to; ++i)
        {
            bvect bv_mask;
            bv_mask.set_range(i, to);
            sv_deserial.deserialize(str_sv2, buf, bv_mask);

            str_sv2.get(1, s2, sizeof(s2));
            assert(s2[0] == 0);

            for (unsigned j = i; j < to; ++j)
            {
                str_sv1.get(j, s1, sizeof(s1));
                str_sv2.get(j, s2, sizeof(s2));
                cmp = ::strcmp(s1, s2);
                assert(cmp==0);

            } // for j

            if ((i & 0xF) == 0)
                cout << "\r" << i << "/" << to << flush;

        } // for

    }
    cout << " ok" << endl;


    cout << "Stress deserialization (AND mask) (use NULL)" << endl;
    {
       str_sparse_vector<char, bvect, 32> str_sv0(bm::use_null);
       str_sparse_vector<char, bvect, 32> str_sv1(bm::use_null);
       str_sparse_vector<char, bvect, 32> str_sv2(bm::use_null);

       {
           str_sparse_vector<char, bvect, 32>::back_insert_iterator bi = str_sv0.get_back_inserter();
           for (unsigned i = 0; i < 100000; ++i)
           {
                bi = "ATGC";
                bi = "GCTA";
                bi = "GCAA";
                bi = "TATA";
                bi.add_null();
           }
       }

       str_sv1.remap_from(str_sv0);
       assert(str_sv1.is_remap());
       str_sv1.optimize();

        BM_DECLARE_TEMP_BLOCK(tb)
        sparse_vector_serial_layout<str_sparse_vector<char, bvect, 32> > sv_lay;
        bm::sparse_vector_serialize<str_sparse_vector<char, bvect, 32> >(str_sv1, sv_lay, tb);
        const unsigned char* buf = sv_lay.buf();

        bm::sparse_vector_deserializer<str_sparse_vector<char, bvect, 32> > sv_deserial;
        char s1[256];
        char s2[256];
        int cmp;

        bvect::size_type from = 100000-1;
        bvect::size_type to = from + 65536;
        for (auto i = from; i < to; ++i)
        {
            bvect bv_mask;
            bv_mask.set_range(i, to);
            sv_deserial.deserialize(str_sv2, buf, bv_mask);

            str_sv2.get(1, s2, sizeof(s2));
            assert(s2[0] == 0);

            for (auto j = i; j < to; ++j)
            {
                str_sv1.get(j, s1, sizeof(s1));
                str_sv2.get(j, s2, sizeof(s2));
                cmp = ::strcmp(s1, s2);
                assert(cmp==0);
                assert(str_sv1.is_null(j) == str_sv2.is_null(j));

            } // for j

            if ((i & 0xF) == 0)
                cout << "\r" << i << "/" << to << flush;

        } // for

    }
    cout << " ok" << endl;


   cout << "---------------------------- TestStrSparseVectorSerial() OK" << endl;
}


typedef str_sparse_vector<char, bvect, 32> str_svect_type;

static
void CompareStrSparseVector(const str_svect_type& str_sv,
                            const vector<string>& str_coll)
{
    assert(str_sv.size() == str_coll.size());
    
    
    string str_h = "z";
    string str_l = "A";

    bm::sparse_vector_scanner<bm::str_sparse_vector<char, bvect, 32> > scanner;

    str_svect_type::const_iterator it = str_sv.begin();
    string str;
    for (unsigned i = 0; i < str_sv.size(); ++i, ++it)
    {
        assert (it.valid());
        assert (it != str_sv.end());
        
        str_sv.get(i, str);
        const string& str_control = str_coll[i];
        if (str != str_control)
        {
            std::cerr << "String mis-match at:" << i << std::endl;
            exit(1);
        }
        {
            const char* s = *it;
            int cmp = ::strcmp(s, str_control.c_str());
            if (cmp != 0)
            {
                cerr << "Iterator comparison failed! " << s << " != " << str_control
                     << endl;
                exit(1);
            }
            str_svect_type::const_iterator it2 = str_sv.get_const_iterator(i);
            assert(it == it2);
            s = *it2;
            cmp = ::strcmp(s, str_control.c_str());
            if (cmp != 0)
            {
                cerr << "2. Iterator comparison failed! " << s << " != " << str_control
                     << endl;
                exit(1);
            }
        }
        int cmp = str_sv.compare(i, str_control.c_str());
        if (cmp != 0)
        {
            std::cerr << "String comparison failure at:" << i << std::endl;
            exit(1);
        }
        if (!str_sv.is_remap()) // re-mapped vectors can give incorrect compare
        {
            cmp = str_sv.compare(i, str_h.c_str());
            if (cmp < 0)
            {
                assert(str < str_h);
            }
            if (cmp > 0)
            {
                assert(str > str_h);
            }

            cmp = str_sv.compare(i, str_l.c_str());
            if (cmp < 0)
            {
                assert(str < str_l);
            }
            if (cmp > 0)
            {
                assert(str > str_l);
            }
        }
        
       unsigned pos;
       bool found = scanner.find_eq_str(str_sv, str_control.c_str(), pos);
       if (!found)
       {
            cerr << "Scanner search failed! " << str_control << endl;
            exit(1);
       }
       assert(pos == i);

        if (i % 100000 == 0)
        {
            cout << "\r" << i << " / " << str_sv.size() << flush;
        }

    } // for
    cout << endl;
}

static 
void GenerateTestStrCollection(std::vector<string>& str_coll, unsigned max_coll)
{
    string prefix = "az";
    string str;
    for (unsigned i = 0; i < max_coll; ++i)
    {
        str = prefix;
        str.append(to_string(i));
        str_coll.emplace_back(str);
        
//        if (i % 1024 == 0) // generate new prefix
        {
            prefix.clear();
            unsigned prefix_len = rand() % 5;
            for (unsigned j = 0; j < prefix_len; ++j)
            {
                char cch = char('a' + rand() % 26);
                prefix.push_back(cch);
            } // for j
        }
    } // for i
}

static
void EraseStrCollection(str_svect_type& str_sv)
{
    std::string s_next, s_curr;
    while (str_sv.size())
    {
        unsigned idx = str_sv.size() / 2;
        unsigned sz = str_sv.size();
        
        if (idx+1 < sz)
        {
            str_sv.get(idx+1, s_next);
        }
        str_sv.erase(idx);
        assert(str_sv.size() == sz-1);
        if (idx+1 < sz)
        {
            str_sv.get(idx, s_curr);
            assert(s_next == s_curr);
        }

    }
}

template<typename SV>
void EraseSVCollection(SV& sv)
{
    typename SV::value_type v_next, v_curr;
    v_next = 0;
    while (sv.size())
    {
        auto idx = sv.size() / 2;
        auto sz = sv.size();
        
        if (idx+1 < sz)
        {
            v_next = sv.get(idx+1);
        }
        sv.erase(idx);
        assert(sv.size() == sz-1);
        if (idx+1 < sz)
        {
            v_curr = sv.get(idx);
            assert(v_next == v_curr);
        }
    }
}


static
void StressTestStrSparseVector()
{
   cout << "---------------------------- Bit-plain STR sparse vector stress test" << endl;
   
   const unsigned max_coll = 2000000;
   std::vector<string> str_coll;
   str_svect_type str_sv;

   GenerateTestStrCollection(str_coll, max_coll);

   cout << "Loading test sparse vector..." << endl;
   {
       str_svect_type::back_insert_iterator bi = str_sv.get_back_inserter();
       for (auto str : str_coll)
       {
           bi = str;
       }
   }

    // -----------------------------------------------------------
    // create sorted collections
    cout << "Sorting str sparse vectors..." << endl;
    vector<string>   str_coll_sorted(str_coll);
    str_svect_type   str_sv_sorted;
    
    std::sort(str_coll_sorted.begin(), str_coll_sorted.end());
    string str_prev;
    for (const string& s : str_coll_sorted)
    {
        if (s != str_prev)
            str_sv_sorted.push_back(s);
        str_prev = s;
    }

    cout << "Build re-mapped vector..." << endl;
    str_svect_type str_sv_remap;
    str_sv_remap.remap_from(str_sv_sorted);
    cout << "Build re-mapped vector... OK" << endl;

    
    // -----------------------------------------------------------

   //print_svector_stat(str_sv);

    cout << "ok. \n Verification..." << endl;

    CompareStrSparseVector(str_sv, str_coll);
    
    cout << "Memory optimization" << endl;
    
    str_sv.optimize();

   print_svector_stat(str_sv, true);

    cout << "ok. \n Verification..." << endl;

    CompareStrSparseVector(str_sv, str_coll);

    cout << "ok. \n Verification of remap vector..." << endl;
    
    CompareStrSparseVector(str_sv_remap, str_coll_sorted);
    
    cout << "Memory optimization" << endl;
    
    str_sv_remap.optimize();

    cout << "ok. \n Verification of remap vector..." << endl;
    
    CompareStrSparseVector(str_sv_remap, str_coll_sorted);


    // serialization check
    //
    cout << "Validate serialization of str-sparse vector..." << endl;
    {
        BM_DECLARE_TEMP_BLOCK(tb)
        sparse_vector_serial_layout<str_svect_type> sv_lay;
        bm::sparse_vector_serialize<str_svect_type>(str_sv, sv_lay, tb);

        str_sparse_vector<char, bvect, 32> str_sv2;
        const unsigned char* buf = sv_lay.buf();
        int res = bm::sparse_vector_deserialize(str_sv2, buf, tb);
        if (res != 0)
        {
            cerr << "De-Serialization error" << endl;
            exit(1);
        }
        CompareStrSparseVector(str_sv2, str_coll);
        bool equal = str_sv.equal(str_sv2);
        assert(equal);
   }
   cout << "Validate serialization of str-sparse vector... OK" << endl;

   cout << "Validate serialization of REMAP str-sparse vector..." << endl;
   {
        BM_DECLARE_TEMP_BLOCK(tb)
        sparse_vector_serial_layout<str_svect_type> sv_lay;
        bm::sparse_vector_serialize<str_svect_type>(str_sv_remap, sv_lay, tb);

        str_sparse_vector<char, bvect, 32> str_sv2;
        const unsigned char* buf = sv_lay.buf();
        int res = bm::sparse_vector_deserialize(str_sv2, buf, tb);
        if (res != 0)
        {
            cerr << "De-Serialization error" << endl;
            exit(1);
        }
        CompareStrSparseVector(str_sv2, str_coll_sorted);
        bool equal = str_sv_remap.equal(str_sv2);
        assert(equal);
   }
   cout << "Validate serialization of REMAP str-sparse vector...OK" << endl;

   // ----------------------------------------------

   cout << "Test common prefix..." << endl;
    {
    const unsigned str_size = 64;
    char str1[str_size];
    char str2[str_size];
    
    unsigned test_size = unsigned(str_coll_sorted.size());
    if (test_size > 20000)
        test_size = 20000;

    for (unsigned i = 0; i < test_size; ++i)
    {
        str_sv_sorted.get(i, &str1[0], str_size);
        for (unsigned j = 0; j < test_size; ++j)
        {
            str_sv_sorted.get(j, &str2[0], str_size);
            unsigned octet_idx = 0;
            for (;true; ++octet_idx)
            {
                if (!str1[octet_idx])
                {
                    if (octet_idx)
                        --octet_idx;
                    break;
                }
                if (!str2[octet_idx])
                {
                    if (octet_idx)
                        --octet_idx;
                    break;
                }
                if (str1[octet_idx] != str2[octet_idx])
                    break;
            }
            unsigned common_prefix = str_sv_sorted.common_prefix_length(i, j);
            if (common_prefix != octet_idx)
            {
                cerr << "Common prefix length mismatch!" <<
                     common_prefix << " != " << octet_idx <<
                     " [" << str1 << "]-[" << str2 << "]" << endl;
                exit(1);
            }
        } // for j
        
        if (i % 512 == 0)
        {
            cout << "\r" << i << " / " << test_size << flush;
        }
    } // for i
    
    }
   cout << "\nTest common prefix...ok." << endl;

   // ----------------------------------------------
   
   cout << "Test sorted search..." << endl;
   
   for (unsigned k = 0; k < 2; ++k)
   {
        bm::sparse_vector_scanner<str_svect_type> scanner;
        bm::sparse_vector_scanner<str_svect_type> scanner2;
        scanner2.bind(str_sv_remap, true); // bind sorted vector

        for (unsigned i = 0; i < unsigned(str_coll_sorted.size()); ++i)
        {
            const string& s = str_coll_sorted[i];
            unsigned pos1, pos2, pos3, pos4;
            bool found1 = scanner.find_eq_str(str_sv_sorted, s.c_str(), pos1);
            if (!found1)
            {
                cerr << "Sorted scan failed at: " << i << " " << s << endl;
                exit(1);
            }
            if (pos1 != i)
            {
                cerr << "Sorted scan position failed at: " << i << "!=" << pos1
                     << " " << s << endl;
                exit(1);
            }
            bool found2 = scanner.bfind_eq_str(str_sv_sorted, s.c_str(), pos2);
            if (!found2)
            {
                cerr << "Error! Sorted binary search failed at: " << i << " value='" << s << "'" << endl;
                //cerr << "Dump file test.sv created." << endl;
                //file_save_svector(str_sv_sorted, "test.sv");
                exit(1);
            }
            if (pos2 != i)
            {
                cerr << "Error! Sorted binary search position mismatch at: " << i << "!=" << pos2
                     << " " << s << endl;
                exit(1);
            }
            bool found4 = scanner.lower_bound_str(str_sv_sorted, s.c_str(), pos4);
            assert(found4);
            assert(pos2 == pos4);
            
            bool found3 = scanner2.bfind_eq_str(s.c_str(), pos3);
            if (!found3)
            {
                cerr << "Error! Sorted-remap binary search failed at: " << i << " value='" << s << "'" << endl;
                exit(1);
            }
            if (pos3 != i)
            {
                cerr << "Error! Sorted-remap binary search position mismatch at: " << i << "!=" << pos2
                     << " value='" << s << "'" << endl;
                exit(1);
            }
            if (i % 65535 == 0)
            {
                cout << "\r" << i << " / " << str_sv_sorted.size() << flush;
            }

        } // for

       str_sv_sorted.optimize();
       cout << "\nPass 2." << endl;
   } // for k
   
   cout << "\nTest sorted search...OK" << endl;

   EraseStrCollection(str_sv_sorted);
   EraseStrCollection(str_sv_remap);
   
   cout << "---------------------------- Bit-plain STR sparse vector stress test OK" << endl;
   cout << endl;
}


static
void TestStrSparseSort()
{
   cout << "---------------------------- Bit-plain STR sparse vector SORT test" << endl;
   const unsigned max_coll = 560000;

   {
       std::vector<string> str_coll;
       str_svect_type      str_sv_sorted;

        // generate sorted vector
        string str;
        for (unsigned i = 10; i < max_coll; i+=10)
        {
            str = to_string(i);
            str_coll.emplace_back(str);
        } // for i
        std::sort(str_coll.begin(), str_coll.end());
        for (const string& s : str_coll)
        {
            str_sv_sorted.push_back(s);
        } // for s
        str_sv_sorted.optimize();
       
        // run lower bound tests
        bm::sparse_vector_scanner<str_svect_type> scanner;

        for (unsigned i = 0; i < max_coll; ++i)
        {
            str = to_string(i);
            
            unsigned pos;
            bool found = scanner.lower_bound_str(str_sv_sorted, str.c_str(), pos);
            string s1;
            if (found)
            {
                str_sv_sorted.get(pos, s1);
                assert(s1 == str);
            }
            
            auto it = std::lower_bound(str_coll.begin(), str_coll.end(), str);
            if (it != str_coll.end())
            {
                unsigned idx = unsigned(it - str_coll.begin());
                const string& s0 = str_coll[idx];
                
                if (s0 == str)
                {
                    assert(found);
                    assert(pos == idx);
                }
                else
                {
                    assert(!found);
                    str_sv_sorted.get(pos, s1);
                    
                    assert(pos == idx);
                }
            }
            if (i % 4096 == 0)
                cout << "\r" << i << "/" << max_coll << flush;

        } // for
        cout << "\n";
    }

    cout << "insertion sort test data generation.." << endl;
    // insertion sort stress test
    {
       std::vector<string> str_coll;
        // generate test values vector
        string str;
        for (unsigned i = 0; i < max_coll; )
        {
            str = to_string(i);
            str_coll.emplace_back(str);
            i += rand() % 3;
        } // for i
        
        // shuffle the data set
        {
            std::random_device rd;
            std::mt19937       g(rd());
            std::shuffle(str_coll.begin(), str_coll.end(), g);
        }

        // insertion sort
        str_svect_type      str_sv_sorted;
        
        cout << "\ninsertion sort..." << endl;
        {
        std::chrono::time_point<std::chrono::steady_clock> st;
        st = std::chrono::steady_clock::now();

            unsigned i = 0;
            bm::sparse_vector_scanner<str_svect_type> scanner;
            for (const string& s : str_coll)
            {
                unsigned pos;
                bool found = scanner.lower_bound_str(str_sv_sorted, s.c_str(), pos);

                auto sz1 = str_sv_sorted.size();
                
                str_sv_sorted.insert(pos, s.c_str());
                
                auto sz2 = str_sv_sorted.size();
                assert(sz1 + 1 == sz2);

                {
                    string str_sv;
                    str_sv_sorted.get(pos, str_sv);
                    assert(s == str_sv);
                }
                
                if (pos)
                {
                    string str_prev;
                    str_sv_sorted.get(pos-1, str_prev);
                    if (str_prev >= s)
                    {
                        cerr << "insertion sort sort order check failed! "
                             << " i = " << i
                             << "s=" << s << " prev=" << str_prev
                             << endl;
                        
                        exit(1);
                    }
                }
                
                {
                    unsigned pos2;
                    found = scanner.lower_bound_str(str_sv_sorted, s.c_str(), pos2);
                    if (!found)
                    {
                        cerr << "control loss at " << i << " " << s << endl;
                        exit(1);
                    }
                    assert(pos == pos2);
                }

                
                if (i % 8096 == 0)
                {
                    std::chrono::time_point<std::chrono::steady_clock> f = std::chrono::steady_clock::now();
                    auto diff = f - st;
                    auto d = std::chrono::duration <double, std::milli> (diff).count();

                    cout << "\r" << i << "/" << max_coll << " (" << d << "ms)" << flush;
                    
                    str_sv_sorted.optimize();
                    
                    st = std::chrono::steady_clock::now();
                }
                ++i;
            } // for s
        }
        cout << endl;
        
        cout << "sort validation.." << endl;
        std::sort(str_coll.begin(), str_coll.end());
        unsigned i = 0;
        string str_prev;
        for (const string& s : str_coll)
        {
            string sv_str;
            str_sv_sorted.get(i, sv_str);
            if (i)
            {
                if (str_prev > sv_str)
                {
                    cerr << "Sort order violation!" << endl;
                    exit(1);
                }
            }
            //cout << s << " = " << sv_str << endl;
            if (s != sv_str)
            {
                cerr << "Sort comparison failed at i=" << i << " s=" << s
                     << " sv_str = " << sv_str << endl;
                
                bm::sparse_vector_scanner<str_svect_type> scanner;
                unsigned pos;
                bool found = scanner.lower_bound_str(str_sv_sorted, s.c_str(), pos);
                
                if (!found)
                {
                    cerr << s << " not found in target." << endl;
                }
                else
                {
                    cerr << s << " is at idx=" << pos << endl;
                }

                exit(1);
            }
            str_prev = sv_str;
            ++i;
        } // for s

        EraseStrCollection(str_sv_sorted);

    }
    
    
    
   cout << "---------------------------- Bit-plain STR sparse vector SORT test OK" << endl;

}

static
void TestSparseSort()
{
   std::cout << "---------------------------- sparse vector SORT test" << endl;
   const unsigned max_coll = 560000;
   typedef bm::sparse_vector<unsigned, bvect > u_svect_type;

   {
       std::vector<unsigned> u_coll;
       u_svect_type          u_sv_sorted;

        // generate sorted vector
        string str;
        for (unsigned i = 10; i < max_coll; i+=10)
        {
            u_coll.emplace_back(i);
        } // for i
        std::sort(u_coll.begin(), u_coll.end());
        for (const unsigned u : u_coll)
        {
            u_sv_sorted.push_back(u);
        } // for s
        u_sv_sorted.optimize();
       
        // run lower bound tests
        bm::sparse_vector_scanner<u_svect_type> scanner;

        for (unsigned i = 0; i < max_coll; ++i)
        {
            bvect::size_type pos;
            bool found = scanner.lower_bound(u_sv_sorted, i, pos);
            unsigned u1;
            if (found)
            {
                u1 = u_sv_sorted[pos];
                assert(u1 == i);
            }
            
            auto it = std::lower_bound(u_coll.begin(), u_coll.end(), i);
            if (it != u_coll.end())
            {
                unsigned idx = unsigned(it - u_coll.begin());
                unsigned u0 = u_coll[idx];
                
                if (u0 == i)
                {
                    assert(found);
                    assert(pos == idx);
                }
                else
                {
                    assert(!found);
                    u1 = u_sv_sorted[pos];
                    
                    assert(pos == idx);
                }
            }
            if (i % 4096 == 0)
                cout << "\r" << i << "/" << max_coll << flush;

        } // for
        cout << "\n";
       
    }
    

    cout << "insertion sort test data generation.." << endl;
    // insertion sort stress test
    {
       std::vector<unsigned> u_coll;
        // generate test values vector
        for (unsigned i = 0; i < max_coll; )
        {
            u_coll.emplace_back(i);
            i += rand() % 3;
        } // for i
        
        // shuffle the data set
        {
            std::random_device rd;
            std::mt19937       g(rd());
            std::shuffle(u_coll.begin(), u_coll.end(), g);
        }

        // insertion sort
        u_svect_type      u_sv_sorted;
        
        cout << "\ninsertion sort..." << endl;
        {
        std::chrono::time_point<std::chrono::steady_clock> st;
        st = std::chrono::steady_clock::now();

            unsigned i = 0;
            bm::sparse_vector_scanner<u_svect_type> scanner;
            for (const unsigned u : u_coll)
            {
                bvect::size_type pos;
                bool found = scanner.lower_bound(u_sv_sorted, u, pos);

                auto sz1 = u_sv_sorted.size();
                
                u_sv_sorted.insert(pos, u);
                
                auto sz2 = u_sv_sorted.size();
                assert(sz1 + 1 == sz2);

                {
                    unsigned u_sv = u_sv_sorted.get(pos);
                    assert(u == u_sv);
                }
                
                if (pos)
                {
                    unsigned u_prev;
                    u_prev = u_sv_sorted.get(pos-1);
                    if (u_prev >= u)
                    {
                        cerr << "insertion sort sort order check failed! "
                             << " i = " << i
                             << "s=" << u << " prev=" << u_prev
                             << endl;
                        assert(0); exit(1);
                    }
                }
                
                {
                    bvect::size_type pos2;
                    found = scanner.lower_bound(u_sv_sorted, u, pos2);
                    if (!found)
                    {
                        cerr << "control loss at " << i << " " << u << endl;
                        assert(0); exit(1);
                    }
                    assert(pos == pos2);
                }

                
                if (i % 8096 == 0)
                {
                    std::chrono::time_point<std::chrono::steady_clock> f = std::chrono::steady_clock::now();
                    auto diff = f - st;
                    auto d = std::chrono::duration <double, std::milli> (diff).count();
                    cout << "\r" << i << "/" << max_coll << " (" << d << "ms)" << flush;
                    
                    u_sv_sorted.optimize();
                    
                    st = std::chrono::steady_clock::now();
                }
                ++i;
            } // for s
        }
        cout << endl;
        
        cout << "sort validation.." << endl;
        std::sort(u_coll.begin(), u_coll.end());
        unsigned i = 0;
        unsigned u_prev = 0;
        for (unsigned u : u_coll)
        {
            unsigned sv_u;
            sv_u = u_sv_sorted.get(i);
            if (i)
            {
                if (u_prev > sv_u)
                {
                    cerr << "Sort order violation!" << endl;
                    assert(0);exit(1);
                }
            }
            //cout << s << " = " << sv_str << endl;
            if (u != sv_u)
            {
                cerr << "Sort comparison failed at i=" << i << " u=" << u
                     << " sv_u = " << sv_u << endl;
                
                bm::sparse_vector_scanner<u_svect_type> scanner;
                bvect::size_type pos;
                bool found = scanner.lower_bound(u_sv_sorted, u, pos);
                
                if (!found)
                {
                    cerr << u << " not found in target." << endl;
                }
                else
                {
                    cerr << u << " is at idx=" << pos << endl;
                }

                exit(1);
            }
            u_prev = sv_u;
            ++i;
        } // for u

        EraseSVCollection(u_sv_sorted);
    }

    
    
   cout << "---------------------------- sparse vector SORT test OK" << endl;

}



inline
void LoadBVDump(const char* filename, const char* filename_out=0, bool validate=false)
{
    ifstream bv_file (filename, ios::in | ios::binary);
    if (!bv_file.good())
    {
        cout << "Cannot open file: " << filename << endl;
        exit(1);
    }

    ofstream* bv_file_out = 0;

    if (filename_out)
    {
        bv_file_out = new ofstream(filename_out, ios::out | ios::binary);
        if (!bv_file_out->good())
        {
            cout << "Cannot create file: " << filename_out << endl;
            exit(1);
        }    
    }


    unsigned buffer_size = 1024*1024;
    unsigned char* buffer = new unsigned char[buffer_size];

    unsigned count = 0;
    clock_t start = clock();
    size_t total_out_size = 0;

    for (;1; ++count)
    {
        unsigned bv_size;
        bv_file.read((char*)&bv_size, sizeof(bv_size));
        if (!bv_file.good())
            break;
        if (bv_size == 0)
        {
            cout << "Warning:Zero vector in dump..." << endl;
            continue;
        }
        if (buffer_size < bv_size)
        {
            delete [] buffer;
            buffer_size = bv_size;
            buffer = new unsigned char[buffer_size];
        }
        bv_file.read((char*)buffer, bv_size);
        {
            bvect bv;
            bm::deserialize(bv, (unsigned char*)buffer);

            bvect::statistics st1;
            bv.calc_stat(&st1);

            if (st1.max_serialize_mem > buffer_size)
            {
                delete [] buffer;
                buffer_size = (unsigned)st1.max_serialize_mem;
                buffer = new unsigned char[buffer_size];
            }

            size_t blob_size = bm::serialize(bv, buffer, BM_NO_GAP_LENGTH|BM_NO_BYTE_ORDER);
            total_out_size += blob_size;

            if (blob_size > bv_size)
            {
                //print_stat(bv);
                //cout << count << ". -" << blob_size-bv_size << endl;
                //exit(1);
            }
            
            if (validate)
            {
                bvect bv_control;
                bm::deserialize(bv_control, (unsigned char*)buffer);
                if (bv_control != bv)
                {
                    cout << "Serialization error!" << endl;
                    exit(1);
                }
            }
            
            if (bv_file_out)
            {
                bv_file_out->write((char*)&blob_size, sizeof(blob_size));
                bv_file_out->write((char*)buffer, (unsigned)blob_size);
            }

        }
        if (count % 1000 == 0)
        {
            cout << count << " out=" << total_out_size << endl;
        }
        //cout << count << ": size=" << bv_size << endl;
    }

    delete [] buffer;
    cout << "Total vectors:" << count << endl;
    cout << "Total out size:" << total_out_size << endl;

    clock_t finish = clock();
    clock_t elapsed_clocks = finish - start;
    double duration = (double)(elapsed_clocks) / CLOCKS_PER_SEC;

    cout << endl
         << "Serialization duration = " << duration 
         << endl;

    bv_file_out->close();
    delete bv_file_out;

}

inline
void GroupByTest(const char* filename, const char* query_filename)
{
    bvect bv_query;

    unsigned count = 0;
    unsigned group_by_count = 0;

    clock_t start = clock();

    // load the query vector
    {
        ifstream bv_file (query_filename, ios::in | ios::binary);
        if (!bv_file.good())
        {
            cout << "Cannot open file: " << query_filename << endl;
            exit(1);
        }
        unsigned buffer_size = 400*1024*1024;
        unsigned char* buffer = new unsigned char[buffer_size];

        unsigned bv_size=0;
        bv_file.read((char*)&bv_size, sizeof(bv_size));
        if (bv_size == 0)
        {
            cout << "Warning:Zero vector in query dump..." << endl;
            return;
        }
        bv_file.read((char*)buffer, bv_size);
        bm::deserialize(bv_query, (unsigned char*)buffer);
        
        delete [] buffer;

    }


    ifstream bv_file (filename, ios::in | ios::binary);
    if (!bv_file.good())
    {
        cout << "Cannot open file: " << filename << endl;
        exit(1);
    }


    unsigned buffer_size = 100*1024*1024;
    unsigned char* buffer = new unsigned char[buffer_size];

    for (;1; ++count)
    {
        unsigned bv_size;
        bv_file.read((char*)&bv_size, sizeof(bv_size));
        if (!bv_file.good())
            break;
        if (bv_size == 0)
        {
            cout << "Warning:Zero vector in dump..." << endl;
            continue;
        }
        if (buffer_size < bv_size)
        {
            delete [] buffer;
            buffer_size = bv_size;
            buffer = new unsigned char[buffer_size];
        }
        bv_file.read((char*)buffer, bv_size);
        bvect bv;
        if (1)
        {
            bv.clear(true);
            bm::deserialize(bv, (unsigned char*)buffer);

            unsigned bc = bm::count_and(bv, bv_query);
            if (bc)
            {
                ++group_by_count;
            }

/*            
            bv &= bv_query;
            if (bv.any())
            {
                ++group_by_count;
            }
*/            

        }
    
    
#if 0
//print_stat(bv_query);
//exit(1);
        {
        bvect bv(BM_GAP);
        operation_deserializer<bvect>::deserialize(bv,
                                                   bv_query,
                                                   (unsigned char*)buffer,
                                                   0,
                                                   bm::set_AND);
        // control			
        if (0)
        {
            bvect bv_control(BM_GAP);
            bm::deserialize(bv_control, (unsigned char*)buffer);
            bv_control &= bv_query;
            if (bv_control != bv)
            {
                cerr << "Group by control failure" << endl;
                cerr << bv_control.count() << endl;
                cerr << bv.count() << endl;
                exit(1);
            }
        }				


        if (bv.any())
        {
            ++group_by_count;
        }
        }
#endif

        if (count % 1000 == 0)
        {
            cout << count << endl;
        }
    }

    delete [] buffer;
    cout << "Total vectors:" << count << endl;
    cout << "Group by vectors:" << group_by_count << endl;

    clock_t finish = clock();
    clock_t elapsed_clocks = finish - start;
    double duration = (double)(elapsed_clocks) / CLOCKS_PER_SEC;

    cout << endl
         << "Test duration = " << duration 
         << endl;
}


inline
void LoadVectors(const char* dir_name, unsigned from, unsigned to)
{
    vector<bvect*>   bv_list;
    vector<unsigned> sz_list;

    size_t total_size = 0;
    size_t total_new_size = 0;

    for(; from <= to; ++from)
    {
        std::stringstream fname_str;
        fname_str << dir_name << "/" << from;
		std::string s = fname_str.str();
        const char* fname = s.c_str();
        
        bvect* bv = new bvect;

        unsigned fsize = 0;
        LoadBVector(fname, *bv, &fsize);
        //bv->optimize();
        //print_stat(*bv);


        // get new size
        size_t blob_size = 0;
        {
        bvect::statistics st1;
        bv->calc_stat(&st1);

        unsigned char* blob = new unsigned char[st1.max_serialize_mem];
        blob_size = bm::serialize(*bv, blob, BM_NO_GAP_LENGTH|BM_NO_BYTE_ORDER);

        if (st1.max_serialize_mem < blob_size)
        {
            printf("BLOB size prediction error!\n");
            exit(1);
        }

        //if (from >= 27)
        {
            bvect bv_control;
            bm::deserialize(bv_control, (unsigned char*)blob);
            if (bv_control != *bv)
            {
                cout << "Serialization error!" << endl;
                exit(1);
            }
        }
                
        delete [] blob;

        }

        cout << fname << "    " 
             << " old=" << fsize << " new=" << blob_size 
             << " diff=" << (int)fsize - (int) blob_size
             << endl;

        bv_list.push_back(bv);
        sz_list.push_back(fsize);

        total_size += fsize;
        total_new_size += blob_size;
    } // for

    cout << "Total size = " << total_size / (1024*1024) << "Mb" << endl;
    cout << "  New size = " << total_new_size / (1024*1024) << "Mb" << endl;
    cout << "Total diff = " << (total_size - total_new_size) / (1024*1024) << "Mb" << endl;

    vector<unsigned char*> bv_blobs;

    cout << "Running serialization benchmark..." << endl;
    {
    clock_t start = clock();

        for (size_t i = 0; i < bv_list.size(); ++i)
        {
            const bvect* bv = bv_list[i];
            bvect::statistics st1;
            bv->calc_stat(&st1);
            unsigned char* blob = new unsigned char[st1.max_serialize_mem*2];
            bv_blobs.push_back(blob);

            for (int j = 0; j < (int)(400/(i?i:1)); ++j)
            {
                //unsigned blob_size = 
                    bm::serialize(*bv, blob);
            }
            // delete [] blob;
        }

    clock_t finish = clock();
    clock_t elapsed_clocks = finish - start;
    double duration = (double)(elapsed_clocks) / CLOCKS_PER_SEC;

    cout << endl
         << "Serialization duration = " << duration 
         << endl;
    }

    cout << "Running de-serialization benchmark..." << endl;
    {
    clock_t start = clock();

        for (size_t i = 0; i < bv_blobs.size(); ++i)
        {
            const unsigned char* blob = bv_blobs[i];
            for (int j = 0; j < (int)(400/(i?i:1)); ++j)
            {
                 bvect bv;
                 bm::deserialize(bv, (unsigned char*)blob);
            }
            // delete [] blob;
        }

    clock_t finish = clock();
    clock_t elapsed_clocks = finish - start;
    double duration = (double)(elapsed_clocks) / CLOCKS_PER_SEC;

    cout << endl
         << "DeSerialization duration = " << duration 
         << endl;
    }




    for (size_t i = 0; i < bv_list.size(); ++i)
    {
        delete bv_list[i];
    }
    for (size_t i = 0; i < bv_blobs.size(); ++i)
    {
        delete [] bv_blobs[i];
    }

}


static
void TestSIMDUtils()
{
    cout << "------------------------ Test SIMD Utils" << endl;
#if defined(BMSSE2OPT)
    unsigned idx;
    cout << "----------------------------> [ SSE2 ]" << endl;
    {
        unsigned short buf[127] = { 65535, 127, 255, 256, 1000, 2000, 2001, 2005, 0xFF, 0, };
        idx = bm::sse2_gap_find(buf, 65535, 1);
        assert(idx == 0);
        idx = bm::sse2_gap_find(buf, 0, 1);
        assert(idx == 0);
        idx = bm::sse2_gap_find(buf, 10, 1);
        assert(idx == 0);
    }

    {
        unsigned short buf[16] = { 60000, 127, 255, 256, 1000, 2000, 2001, 2005, 0xFF, 0, };
        idx = bm::sse2_gap_find(buf, 65535, 1);
        assert(idx == 1);
        idx = bm::sse2_gap_find(buf, 0, 1);
        assert(idx == 0);
        idx = bm::sse2_gap_find(buf, 10, 1);
        assert(idx == 0);
    }

    {
        unsigned short buf[16] = { 10, 65530, 127, 255, 256, 1000, 2000, 2001, 2005, 0xFF, };
        const unsigned vsize = 2;
        idx = bm::sse2_gap_find(buf, 0, vsize);
        assert(idx == 0);
        idx = bm::sse2_gap_find(buf, 65535, vsize);
        assert(idx == vsize);
        for (unsigned i = 0; i < vsize; ++i)
        {
            unsigned short v = buf[i];
            idx = bm::sse2_gap_find(buf, v, vsize);
            assert(idx == i);
            idx = bm::sse2_gap_find(buf, gap_word_t(v - 1), vsize);
            assert(idx == i);
        }
    }

    {
        unsigned short buf[16] = { 10, 256, 65530, 127, 255, 256, 1000, 2000, 2001, 2005, 0xFF, };
        const unsigned vsize = 3;
        idx = bm::sse2_gap_find(buf, 0, vsize);
        assert(idx == 0);
        idx = bm::sse2_gap_find(buf, 65535, vsize);
        assert(idx == vsize);
        for (unsigned i = 0; i < vsize; ++i)
        {
            unsigned short v = buf[i];
            idx = bm::sse2_gap_find(buf, v, vsize);
            assert(idx == i);
            idx = bm::sse2_gap_find(buf, gap_word_t(v - 1), vsize);
            assert(idx == i);
        }
    }
    {
        unsigned short buf[16] = { 10, 256, 258, 65500, 127, 255, 256, 1000, 2000, 2001, 2005, 0xFF, };
        const unsigned vsize = 4;
        idx = bm::sse2_gap_find(buf, 0, vsize);
        assert(idx == 0);
        idx = bm::sse2_gap_find(buf, 65535, vsize);
        assert(idx == vsize);
        for (unsigned i = 0; i < vsize; ++i)
        {
            unsigned short v = buf[i];
            idx = bm::sse2_gap_find(buf, v, vsize);
            assert(idx == i);
            idx = bm::sse2_gap_find(buf, gap_word_t(v - 1), vsize);
            assert(idx == i);
        }
    }

    {
        //        unsigned short buf[16] = { 10, 256, 258, 15525, 64500, 127, 255, 256, 1000, 2000, 2001, 2005, 0xFF, };
        unsigned short buf[16] = { 10, 20, 30, 40, 50, 127, 255, 256, 1000, 2000, 2001, 2005, 0xFF, };
        const unsigned vsize = 5;
        idx = bm::sse2_gap_find(buf, 0, vsize);
        assert(idx == 0);
        idx = bm::sse2_gap_find(buf, 65535, vsize);
        assert(idx == vsize);

        for (unsigned i = 0; i < vsize; ++i)
        {
            unsigned short v = buf[i];
            idx = bm::sse2_gap_find(buf, v, vsize);
            assert(idx == i);
            idx = bm::sse2_gap_find(buf, gap_word_t(v - 1), vsize);
            assert(idx == i);
        }
    }

    {
        unsigned short buf[127] = { 2, 10, 256, 258, 15525, 65530, 127, 255, 256, 1000, 2000, 2001, 2005, 0xFF, };
        const unsigned vsize = 6;
        idx = bm::sse2_gap_find(buf, 0, vsize);
        assert(idx == 0);
        idx = bm::sse2_gap_find(buf, 65535, vsize);
        assert(idx == vsize);
        for (unsigned i = 0; i < vsize; ++i)
        {
            unsigned short v = buf[i];
            idx = bm::sse2_gap_find(buf, v, vsize);
            assert(idx == i);
            idx = bm::sse2_gap_find(buf, gap_word_t(v - 1), vsize);
            assert(idx == i);
        }
    }

    {
        unsigned short buf[16] = { 1, 2, 10, 256, 258, 15525, 65530, 127, 255, 256, 1000, 2000, 2001, 2005, 0xFF, };
        const unsigned vsize = 7;
        idx = bm::sse2_gap_find(buf, 0, vsize);
        assert(idx == 0);
        idx = bm::sse2_gap_find(buf, 65535, vsize);
        assert(idx == vsize);

        for (unsigned i = 0; i < vsize; ++i)
        {
            unsigned short v = buf[i];
            idx = bm::sse2_gap_find(buf, v, vsize);
            assert(idx == i);
            idx = bm::sse2_gap_find(buf, gap_word_t(v - 1), vsize);
            assert(idx == i || (buf[i - 1] == v - 1));
        }
    }
    {
        unsigned short buf[16] = { 1, 2, 10, 256,  258, 1024, 15525, 65530,  127, 255, 256, 0xFF, };
        const unsigned vsize = 8;
        idx = bm::sse2_gap_find(buf, 0, vsize);
        assert(idx == 0);
        idx = bm::sse2_gap_find(buf, 65535, vsize);
        assert(idx == vsize);
        for (unsigned i = 0; i < vsize; ++i)
        {
            unsigned short v = buf[i];
            if (i == vsize - 1)
            {
                assert(v != 65535);
            }
            idx = bm::sse2_gap_find(buf, v, vsize);
            assert(idx == i);
            idx = bm::sse2_gap_find(buf, gap_word_t(v - 1), vsize);
            assert(idx == i || (buf[i - 1] == v - 1));
        }
    }

    {
        unsigned short buf[16] = { 6217, 6300, 6400, 6500,
            6600, 6700, 30584, 40255,
            50256, 60000, 61001, 65255, 65530, 12, 23, 0, };
        const unsigned vsize = 13;
        idx = bm::sse2_gap_find(buf, 0, vsize);
        assert(idx == 0);
        idx = bm::sse2_gap_find(buf, 65535, vsize);
        assert(idx == vsize);

        for (unsigned i = 0; i < vsize; ++i)
        {
            unsigned short v = buf[i];
            if (i == vsize - 1)
            {
                assert(v != 65535);
            }
            idx = bm::sse2_gap_find(buf, v, vsize);
            assert(idx == i);
            idx = bm::sse2_gap_find(buf, gap_word_t(v - 1), vsize);
            assert(idx == i || (buf[i - 1] == v - 1));
        }
    }

    {
        unsigned short buf[16] = { 6217, 6300, 6400, 6500,
            6600, 6700, 30584, 40255,
            50256, 60000, 61001, 65255,
            65256, 65257, 65300, 65530 };
        const unsigned vsize = 16;
        idx = bm::sse2_gap_find(buf, 0, vsize);
        assert(idx == 0);
        idx = bm::sse2_gap_find(buf, 65535, vsize);
        assert(idx == 16);
        for (unsigned i = 0; i < vsize; ++i)
        {
            unsigned short v = buf[i];
            if (i == vsize - 1)
            {
                assert(v != 65535);
            }
            idx = bm::sse2_gap_find(buf, v, vsize);
            assert(idx == i);
            idx = bm::sse2_gap_find(buf, gap_word_t(v - 1), vsize);
            assert(idx == i || (buf[i - 1] == v - 1));
        }
    }
    
    // SSE2 AND block check
    

#endif

#if defined(BMSSE42OPT)
    unsigned idx;
    cout << "----------------------------> [ SSE 4.2 ]" << endl;
    
    {
        BM_DECLARE_TEMP_BLOCK(tb)
        for (unsigned i = 0; i < bm::set_block_size; ++i)
        {
            tb[i] = ~0u;
        }
        bool all_one = sse4_is_all_one((__m128i*)tb);
        assert(all_one);
        
        tb[256] = 1;
        all_one = sse4_is_all_one((__m128i*)tb);
        assert(!all_one);
    }

    {
        BM_DECLARE_TEMP_BLOCK(tb)
        for (unsigned i = 0; i < bm::set_block_size; ++i)
        {
            tb[i] = 0u;
        }
        bool all_z = sse4_is_all_zero((__m128i*)tb);
        assert(all_z);
        
        tb[256] = 1;
        all_z = sse4_is_all_zero((__m128i*)tb);
        assert(!all_z);
    }


    {
        unsigned short buf[127] = { 65535, 127, 255, 256, 1000, 2000, 2001, 2005, 0xFF, 0,  };
        idx = bm::sse4_gap_find(buf, (unsigned short)65535, 1);
        assert(idx == 0);
        idx = bm::sse4_gap_find(buf, 0, 1);
        assert(idx == 0);
        idx = bm::sse4_gap_find(buf, 10, 1);
        assert(idx == 0);
    }

    {
        unsigned short buf[16] = { 60000, 127, 255, 256, 1000, 2000, 2001, 2005, 0xFF, 0, };
        idx = bm::sse4_gap_find(buf, 65535, 1);
        assert(idx == 1);
        idx = bm::sse4_gap_find(buf, 0, 1);
        assert(idx == 0);
        idx = bm::sse4_gap_find(buf, 10, 1);
        assert(idx == 0);
    }

    {
        unsigned short buf[16] = { 10, 65530, 127, 255, 256, 1000, 2000, 2001, 2005, 0xFF, };
        const unsigned vsize = 2;
        idx = bm::sse4_gap_find(buf, 0, vsize);
        assert(idx == 0);
        idx = bm::sse4_gap_find(buf, 65535, vsize);
        assert(idx == vsize);
        for (unsigned i = 0; i < vsize; ++i)
        {
            unsigned short v = buf[i];
            idx = bm::sse4_gap_find(buf, v, vsize);
            assert(idx == i);
            idx = bm::sse4_gap_find(buf, (unsigned short)(v - 1), vsize);
            assert(idx == i);
        }
    }

    {
        unsigned short buf[16] = { 10, 256, 65530, 127, 255, 256, 1000, 2000, 2001, 2005, 0xFF, };
        const unsigned vsize = 3;
        idx = bm::sse4_gap_find(buf, 0, vsize);
        assert(idx == 0);
        idx = bm::sse4_gap_find(buf, 65535, vsize);
        assert(idx == vsize);
        for (unsigned i = 0; i < vsize; ++i)
        {
            unsigned short v = buf[i];
            idx = bm::sse4_gap_find(buf, v, vsize);
            assert(idx == i);
            idx = bm::sse4_gap_find(buf, (unsigned short)(v - 1), vsize);
            assert(idx == i);
        }
    }
    {
        unsigned short buf[16] = { 10, 256, 258, 65500, 127, 255, 256, 1000, 2000, 2001, 2005, 0xFF, };
        const unsigned vsize = 4;
        idx = bm::sse4_gap_find(buf, 0, vsize);
        assert(idx == 0);
        idx = bm::sse4_gap_find(buf, (unsigned short)(65535), vsize);
        assert(idx == vsize);
        for (unsigned i = 0; i < vsize; ++i)
        {
            unsigned short v = buf[i];
            idx = bm::sse4_gap_find(buf, v, vsize);
            assert(idx == i);
            idx = bm::sse4_gap_find(buf, (unsigned short)(v - 1), vsize);
            assert(idx == i);
        }
    }
    
    {
//        unsigned short buf[16] = { 10, 256, 258, 15525, 64500, 127, 255, 256, 1000, 2000, 2001, 2005, 0xFF, };
        unsigned short buf[16] = { 10, 20, 30, 40, 50, 127, 255, 256, 1000, 2000, 2001, 2005, 0xFF, };
        const unsigned vsize = 5;
        idx = bm::sse4_gap_find(buf, 0, vsize);
        assert(idx == 0);
        idx = bm::sse4_gap_find(buf, 65535, vsize);
        assert(idx == vsize);

        for (unsigned i = 0; i < vsize; ++i)
        {
            unsigned short v = buf[i];
            idx = bm::sse4_gap_find(buf, v, vsize);
            assert(idx == i);
            idx = bm::sse4_gap_find(buf, (unsigned short)(v - 1), vsize);
            assert(idx == i);
        }
    }
    
    {
        unsigned short buf[127] = { 2, 10, 256, 258, 15525, 65530, 127, 255, 256, 1000, 2000, 2001, 2005, 0xFF, };
        const unsigned vsize = 6;
        idx = bm::sse4_gap_find(buf, 0, vsize);
        assert(idx == 0);
        idx = bm::sse4_gap_find(buf, 65535, vsize);
        assert(idx == vsize);
        for (unsigned i = 0; i < vsize; ++i)
        {
            unsigned short v = buf[i];
            idx = bm::sse4_gap_find(buf, v, vsize);
            assert(idx == i);
            idx = bm::sse4_gap_find(buf, (unsigned short)(v - 1), vsize);
            assert(idx == i);
        }
    }

    {
        unsigned short buf[16] = { 1, 2, 10, 256, 258, 15525, 65530, 127, 255, 256, 1000, 2000, 2001, 2005, 0xFF, };
        const unsigned vsize = 7;
        idx = bm::sse4_gap_find(buf, 0, vsize);
        assert(idx == 0);
        idx = bm::sse4_gap_find(buf, 65535, vsize);
        assert(idx == vsize);

        for (unsigned i = 0; i < vsize; ++i)
        {
            unsigned short v = buf[i];
            idx = bm::sse4_gap_find(buf, v, vsize);
            assert(idx == i);
            idx = bm::sse4_gap_find(buf, (unsigned short)(v - 1), vsize);
            assert(idx == i || (buf[i-1] == v-1));
        }
    }
    {
    unsigned short buf[16] = { 1, 2, 10, 256,  258, 1024, 15525, 65530,  127, 255, 256, 0xFF, };
    const unsigned vsize = 8;
    idx = bm::sse4_gap_find(buf, 0, vsize);
    assert(idx == 0);
    idx = bm::sse4_gap_find(buf, (unsigned short)(65535), vsize);
    assert(idx == vsize);
    for (unsigned i = 0; i < vsize; ++i)
    {
        unsigned short v = buf[i];
        if (i == vsize - 1)
        {
            assert(v != 65535);
        }
        idx = bm::sse4_gap_find(buf, v, vsize);
        assert(idx == i);
        idx = bm::sse4_gap_find(buf, (unsigned short)(v - 1), vsize);
        assert(idx == i || (buf[i - 1] == v - 1));
    }
    }

    {
        unsigned short buf[16] = { 6217, 6300, 6400, 6500,  
                                   6600, 6700, 30584, 40255, 
                                   50256, 60000, 61001, 65255, 65530, 12, 23, 0,  };
        const unsigned vsize = 13;
        idx = bm::sse4_gap_find(buf, 0, vsize);
        assert(idx == 0);
        idx = bm::sse4_gap_find(buf, (unsigned short)(65535), vsize);
        assert(idx == vsize);

        for (unsigned i = 0; i < vsize; ++i)
        {
            unsigned short v = buf[i];
            if (i == vsize - 1)
            {
                assert(v != 65535);
            }
            idx = bm::sse4_gap_find(buf, v, vsize);
            assert(idx == i);
            idx = bm::sse4_gap_find(buf, (unsigned short)(v - 1), vsize);
            assert(idx == i || (buf[i - 1] == v - 1));
        }
    }

    {
        unsigned short buf[16] = { 6217, 6300, 6400, 6500,
            6600, 6700, 30584, 40255,
            50256, 60000, 61001, 65255, 
            65256, 65257, 65300, 65530  };
        const unsigned vsize = 16;
        idx = bm::sse4_gap_find(buf, 0, vsize);
        assert(idx == 0);
        idx = bm::sse4_gap_find(buf, 65535, vsize);
        assert(idx == 16);
        for (unsigned i = 0; i < vsize; ++i)
        {
            unsigned short v = buf[i];
            if (i == vsize - 1)
            {
                assert(v != 65535);
            }
            idx = bm::sse4_gap_find(buf, v, vsize);
            assert(idx == i);
            idx = bm::sse4_gap_find(buf, (unsigned short)(v - 1), vsize);
            assert(idx == i || (buf[i - 1] == v - 1));
        }
    }
    
    
    // unsigned vector GE search
    
    {
        __m128i vect4 = _mm_set_epi32(-1, int(0x80000000u), 8, 0);
        
        int pos = bm::sse42_cmpge_u32(vect4, 0);
        assert(pos == 0);
        
        pos = bm::sse42_cmpge_u32(vect4, 7);
        assert(pos == 1);

        pos = bm::sse42_cmpge_u32(vect4, 8);
        assert(pos == 1);

        pos = bm::sse42_cmpge_u32(vect4, 0x80000000u);
        assert(pos == 2);

        pos = bm::sse42_cmpge_u32(vect4, ~0u);
        assert(pos == 3);

        pos = bm::sse42_cmpge_u32(vect4, ~0u - 1u);
        assert(pos == 3);
    }
    
    {
        __m128i vect4 = _mm_set_epi32(-1, -1, 8, 8);
        
        int pos = bm::sse42_cmpge_u32(vect4, 0);
        assert(pos == 0);
        
        pos = bm::sse42_cmpge_u32(vect4, 5);
        assert(pos == 0);
        
        pos = bm::sse42_cmpge_u32(vect4, 9);
        assert(pos == 2);

        pos = bm::sse42_cmpge_u32(vect4, ~0u);
        assert(pos == 2);

        pos = bm::sse42_cmpge_u32(vect4, ~0u - 1u);
        assert(pos == 2);

    }
    
    // lower bound SSE42 scan
    
    const unsigned arr_size = 50000;
    unsigned arr[arr_size + 10] = {0, };
    
    for (unsigned i = 0; i < arr_size; ++i)
    {
        arr[i] = 10 + i;
    }
    
    {
        unsigned target, s_idx;
        for (unsigned i = 0; i < arr_size-1; ++i)
        {
            target = arr[i];
            s_idx = bm::sse4_lower_bound_scan_u32(&arr[0], target, 0, arr_size-1);
            assert(s_idx == i);
            s_idx = bm::sse4_lower_bound_scan_u32(&arr[0], target, i, arr_size-1);
            assert(s_idx == i);
            
            target = 1; // not found but lower
            s_idx = bm::sse4_lower_bound_scan_u32(&arr[0], target, 0, arr_size-1);
            assert(s_idx == 0);
            
            target = arr_size * 2; // not found but higher
            s_idx = bm::sse4_lower_bound_scan_u32(&arr[0], target, 0, arr_size-1);
            assert(s_idx == arr_size);
        }

    }

#endif

#if defined(BMAVX2OPT)
    cout << "----------------------------> [ AVX2 ]" << endl;
    
    {
        BM_DECLARE_TEMP_BLOCK(tb)
        for (unsigned i = 0; i < bm::set_block_size; ++i)
        {
            tb[i] = 0;
        }
        bool all_z = avx2_is_all_zero((__m256i*)tb);
        assert(all_z);
        
        tb[256] = 1;
        all_z = avx2_is_all_zero((__m256i*)tb);
        assert(!all_z);
    }

    {
        BM_DECLARE_TEMP_BLOCK(tb)
        for (unsigned i = 0; i < bm::set_block_size; ++i)
        {
            tb[i] = ~0u;
        }
        bool all_one = avx2_is_all_one((__m256i*)tb);
        assert(all_one);
        
        tb[256] = 1;
        all_one = avx2_is_all_one((__m256i*)tb);
        assert(!all_one);
    }

    // AVX2 unsigned vector GE search
    {
        __m256i vect8 = _mm256_set_epi32(-1, -1, int(0x80000000u), 8, 5, 5, 4, 0);

        int pos = bm::avx2_cmpge_u32(vect8, 0);
        assert(pos == 0);

        pos = bm::avx2_cmpge_u32(vect8, 7);
        assert(pos == 4);

        pos = bm::avx2_cmpge_u32(vect8, 8);
        assert(pos == 4);

        pos = bm::avx2_cmpge_u32(vect8, 0x80000000u);
        assert(pos == 5);

        pos = bm::avx2_cmpge_u32(vect8, ~0u);
        assert(pos == 6);

        pos = bm::avx2_cmpge_u32(vect8, ~0u - 1u);
        assert(pos == 6);
    }

    cout << "avx2_cmpge_u32 stress" << endl;
    {
        for (unsigned i = 1; i < ~0u; ++i)
        {
            __m256i vect8 = _mm256_set_epi32(-1, -1, i, 0, 0, 0, 0, 0);
            int pos = bm::avx2_cmpge_u32(vect8, i);
            assert(pos == 5);
        }
    }
    cout << " - ok " << endl;


    {
        // lower bound avx2 scan

        const unsigned arr_size = 50000;
        unsigned arr[arr_size + 10] = { 0, };

        for (unsigned i = 0; i < arr_size; ++i)
        {
            arr[i] = 10 + i;
        }

        {
            unsigned target, s_idx;
            for (unsigned i = 0; i < arr_size - 1; ++i)
            {
                target = arr[i];
                s_idx = bm::avx2_lower_bound_scan_u32(&arr[0], target, 0, arr_size - 1);
                assert(s_idx == i);
                s_idx = bm::avx2_lower_bound_scan_u32(&arr[0], target, i, arr_size - 1);
                assert(s_idx == i);

                target = 1; // not found but lower
                s_idx = bm::avx2_lower_bound_scan_u32(&arr[0], target, 0, arr_size - 1);
                assert(s_idx == 0);

                target = arr_size * 2; // not found but higher
                s_idx = bm::avx2_lower_bound_scan_u32(&arr[0], target, 0, arr_size - 1);
                assert(s_idx == arr_size);
            }
        }
    }

#endif


    cout << "------------------------ Test SIMD Utils OK" << endl;
}

static
void AddressResolverTest()
{
    bm::id_t id_to;
    bool found;

    {
        bvps_addr_resolver<bvect>  ares;
        
        found = ares.resolve(10, &id_to);
        assert(!found);
        assert(id_to == 0);
        
        {
        bvps_addr_resolver<bvect>  ares2(ares);
        
        found = ares2.resolve(10, &id_to);
        assert(!found);
        assert(id_to == 0);
        }
        
        found = ares.resolve(10, &id_to);
        assert(!found);
        assert(id_to == 0);

    }

    {
        bvps_addr_resolver<bvect>  ares;
        
        ares.set(1000);
        ares.set(10000);
        ares.set(100000);

        found = ares.resolve(10, &id_to);
        assert(!found);
        assert(id_to == 0);

        found = ares.resolve(100000, &id_to);
        assert(found);
        assert(id_to == 3);
        
        assert(ares.in_sync() == false);
        
        ares.optimize();
        assert(ares.in_sync() == false);

        ares.sync();
        assert(ares.in_sync());

        found = ares.resolve(100000, &id_to);
        assert(found);
        assert(id_to == 3);

        bvps_addr_resolver<bvect>  ares2(ares);
        bool same = ares.equal(ares2);
        assert(same);
        
        bvps_addr_resolver<bvect>  ares3;
        ares3.move_from(ares2);
        same = ares.equal(ares3);
        assert(same);

    }

    {
        sv_addr_resolver<sparse_vector<bm::id_t, bvect> > ares;
        
        found = ares.resolve(10, &id_to);
        assert(!found);
        assert(id_to == 0);
    }
    
    {
        sv_addr_resolver<sparse_vector<bm::id_t, bvect> > ares;
        
        ares.set(1000);  // 1
        ares.set(10000); // 2
        ares.set(bm::id_max32-1); // 3
        ares.set(5);      // 4
        
        found = ares.resolve(10, &id_to);
        assert(!found);
        assert(id_to == 0);
        
        found = ares.resolve(1000, &id_to);
        assert(found);
        assert(id_to == 1);

        found = ares.resolve(bm::id_max32-1, &id_to);
        assert(found);
        assert(id_to == 3);
        
        ares.optimize();
        
        found = ares.resolve(5, &id_to);
        assert(found);
        assert(id_to == 4);
    }

    
}

// generate pseudo-random bit-vector, mix of blocks
//
void generate_bvector(bvect& bv, unsigned vector_max, bool optimize)
{
    unsigned i, j;
    for (i = 0; i < vector_max;)
    {
        // generate bit-blocks
        for (j = 0; j < 65535*8; i += 10, j++)
        {
            bv.set(i);
        }
        if (i > vector_max)
            break;
        // generate GAP (compressed) blocks
        for (j = 0; j < 65535; i += 120, j++)
        {
            unsigned len = rand() % 64;
            bv.set_range(i, i + len);
            i += len;
            if (i > vector_max)
                break;
        }
    }
    if (optimize)
        bv.optimize();
}

extern "C" {
    static
    int bit_decode_func(void* handle_ptr, bm::id_t bit_idx)
    {
        std::vector<bm::id_t>* vp = (std::vector<bm::id_t>*)handle_ptr;
        vp->push_back(bit_idx);
        return 0;
    }
    
    static
    int bit_decode_func2(void* handle_ptr, bm::id_t bit_idx)
    {
        if (bit_idx > (65536 * 256))
        {
            throw 1;
        }
        std::vector<bm::id_t>* vp = (std::vector<bm::id_t>*)handle_ptr;
        vp->push_back(bit_idx);
        return 0;
    }
} // extern C


static
void RangeForEachTest(bvect::size_type from, bvect::size_type to)
{
    bvect bv;
    bv.set_range(from, to);
    std::vector<bvect::size_type> v;
    bm::visit_each_bit(bv, (void*)&v, bit_decode_func);
    
    assert(v.size() == bv.count());
    assert(v[0] == from);
    
    for (size_t i = 1; i < v.size(); ++i)
    {
        bvect::size_type prev = v[i-1];
        bvect::size_type curr = v[i];
        assert(prev+1 == curr);
    }
    bvect bv_control;
    bm::combine_or(bv_control, v.begin(), v.end());
    int res = bv.compare(bv_control);
    assert(res == 0);
}

static
void BvectorBitForEachTest()
{
    cout << "------------------------ bvector BitForEach Test" << endl;
    int res;
    
    {
        cout << "test empty vector" << endl;
        bvect bv1;
        std::vector<unsigned> v1;
        
        bm::visit_each_bit(bv1, (void*)&v1, bit_decode_func);
        if (v1.size())
        {
            cerr << "1. Failed empty vector decode " << v1.size() << endl;
            exit(1);
        }
        bv1.init();
        bm::visit_each_bit(bv1, (void*)&v1, bit_decode_func);
        if (v1.size())
        {
            cerr << "2. Failed empty vector decode " << v1.size() << endl;
            exit(1);
        }
        
        bv1.set(100000, true);
        bv1.set(100000, false);
        
        bm::visit_each_bit(bv1, (void*)&v1, bit_decode_func);
        if (v1.size())
        {
            cerr << "3. Failed empty vector decode " << v1.size() << endl;
            exit(1);
        }

        bv1.optimize();
        bm::visit_each_bit(bv1, (void*)&v1, bit_decode_func);
        if (v1.size())
        {
            cerr << "3. Failed empty vector decode " << v1.size() << endl;
            exit(1);
        }
    }
    
    {
        bvect bv1 { 0,1,2, 10, 32, 100, 65535,
                            65535+1, 65535+2, 65535+10, 65535+11, 65535+31,
                            20000000 };
        bvect bv2;
        std::vector<unsigned> v1;
        
        bm::visit_each_bit(bv1, (void*)&v1, bit_decode_func);

        {
            for (size_t i = 0; i < v1.size(); ++i)
                cout << v1[i] << ", ";
            cout << endl;
        }

        if (v1.size() != bv1.count())
        {
            std::cerr << "0. Bit for each failed size test. " << v1.size()
                      << " should be " << bv1.count() << std::endl;
            exit(1);
        }
        bm::combine_or(bv2, v1.begin(), v1.end());
        
        res = bv2.compare(bv1);
        if (res != 0)
        {
            std::cerr << "0. Bit for each failed comparison test. " << endl;
            exit(1);
        }
        
        bv1.optimize();
        bv2.reset();
        v1.resize(0);
        bm::visit_each_bit(bv1, (void*)&v1, bit_decode_func);
        
        {
            for (size_t i = 0; i < v1.size(); ++i)
                cout << v1[i] << ", ";
            cout << endl;
        }

        if (v1.size() != bv1.count())
        {
            std::cerr << "1. Bit for each failed size test. " << v1.size()
                      << " should be " << bv1.count() << std::endl;
            exit(1);
        }
        bm::combine_or(bv2, v1.begin(), v1.end());
        
        res = bv2.compare(bv1);
        if (res != 0)
        {
            std::cerr << "1. Bit for each failed comparison test. " << endl;
            exit(1);
        }

    }
    
    {
        bvect bv1, bv2;
        std::vector<unsigned> v1;
        
        generate_bvector(bv1);
        v1.reserve(bv1.count());
        
        bm::visit_each_bit(bv1, (void*)&v1, bit_decode_func);

        if (v1.size() != bv1.count())
        {
            std::cerr << "Bit for each failed size test. " << v1.size()
                      << " should be " << bv1.count() << std::endl;
            exit(1);
        }
        bm::combine_or(bv2, v1.begin(), v1.end());
        
        res = bv2.compare(bv1);
        if (res != 0)
        {
            std::cerr << "Bit for each failed comparison test. " << endl;
            exit(1);
        }

        cout << "for each bit in optimized bit-vector..." << endl;
        
        v1.resize(0);
        bv2.clear(true);
        bv1.optimize();

        bm::visit_each_bit(bv1, (void*)&v1, bit_decode_func);

        if (v1.size() != bv1.count())
        {
            std::cerr << "Bit for each failed size test. " << v1.size()
                      << " should be " << bv1.count() << std::endl;
            exit(1);
        }
        bm::combine_or(bv2, v1.begin(), v1.end());
        
        res = bv2.compare(bv1);
        if (res != 0)
        {
            std::cerr << "Bit for each failed comparison test. " << endl;
            exit(1);
        }
    }
    
    {
        RangeForEachTest(0, 65536);
        RangeForEachTest(65536, 65536+65536);
        RangeForEachTest(bm::id_max/2, bm::id_max/2 + (65536*256));
    }

    {
        bvect bv;
        bv.set();
        
        std::vector<bvect::size_type> v1;
        try
        {
            bm::visit_each_bit(bv, (void*)&v1, bit_decode_func2);
        } catch (...)
        {
        }

    }
    
    
    cout << "------------------------ bvector BitForEach Test OK" << endl;
}

static
void FillTestBuffer(bm::compressed_buffer_collection<bvect>::buffer_type& buf)
{
    unsigned sz_factor = rand() % 10;
    if (!sz_factor)
        sz_factor = 1;
    unsigned size = 65000 + rand() % 65;// (128000 / sz_factor);
    
    buf.resize(size);
    unsigned char* data = buf.data();
    
    for (unsigned i = 0; i < size; ++i)
    {
        data[i] = (unsigned char)i;
    }
}

static
void GenerateCompressedBufferCollection(bm::compressed_buffer_collection<bvect>& cbc)
{
    unsigned sz = rand() % 10000;
    unsigned key = 0;
    unsigned key_factor = rand() % 128;
    if (!key_factor)
        key_factor = 1;
    for (unsigned i = 0; i < sz; ++i)
    {
        {
            bm::compressed_buffer_collection<bvect>::buffer_type buf;
            FillTestBuffer(buf);
            cbc.move_buffer(key, buf);
        }
        key += key_factor;
    } // for
    cbc.sync();
}

static
void TestCompressedCollection()
{
    cout << "------------------------ Compressed collection Test" << endl;
    
    {
        bm::compressed_collection<unsigned, bvect> coll;
        bool added;
        unsigned v;
        
        added = coll.push_back(0, 100);
        assert(added);
        added = coll.push_back(10, 5);
        assert(added);
        added = coll.push_back(150000, 500);
        assert(added);
        
        coll.sync();
        
        bool found;
        bm::id_t idx;
        
        found = coll.resolve(0, &idx);
        assert(found);
        v = coll.get(idx);
        assert(v == 100);

        found = coll.resolve(10, &idx);
        assert(found);
        v = coll.get(idx);
        assert(v == 5);

        found = coll.resolve(150000, &idx);
        assert(found);
        v = coll.get(idx);
        assert(v == 500);

        found = coll.resolve(256, &idx);
        assert(!found);
        
        coll.optimize();
        
        found = coll.resolve(0, &idx);
        assert(found);
        v = coll.get(idx);
        assert(v == 100);

        found = coll.resolve(10, &idx);
        assert(found);
        v = coll.get(idx);
        assert(v == 5);

        found = coll.resolve(150000, &idx);
        assert(found);
        v = coll.get(idx);
        assert(v == 500);

        found = coll.resolve(256, &idx);
        assert(!found);

    }
    
    {
    bm::compressed_buffer_collection<bvect> cbc;
    {
        bm::compressed_buffer_collection<bvect>::buffer_type buf;
        buf.copy_from((unsigned char*)"ABC", 3);
        cbc.move_buffer(10, buf);
    }
    {
        bm::compressed_buffer_collection<bvect>::buffer_type buf;
        buf.copy_from((unsigned char*)"1234", 4);
        cbc.move_buffer(15, buf);
    }
    cbc.sync();
    
    assert(cbc.size() == 2);
    bm::compressed_buffer_collection<bvect>::statistics st;
    cbc.calc_stat(&st);
    
    bm::compressed_collection_serializer<compressed_buffer_collection<bvect> > cbcs;
    bm::compressed_buffer_collection<bvect>::buffer_type sbuf;
    
    cbcs.serialize(cbc, sbuf);
    assert(sbuf.size() > 0);
    
    bm::compressed_buffer_collection<bvect>::buffer_type sbuf2(sbuf);
    bm::compressed_buffer_collection<bvect> cbc2;
    compressed_collection_deserializer<compressed_buffer_collection<bvect> > cbcd;
    cbcd.deserialize(cbc2, sbuf2.buf());
    
    if (!cbc2.equal(cbc))
    {
        std::cerr << "Compressed collection serialization error" << endl;
        exit(1);
    }

    }
    
    {
        cout << "Compressed buffer collection stress." << endl;
        const unsigned test_count = 60;
        for (unsigned i = 0; i < test_count; ++i)
        {
            bm::compressed_buffer_collection<bvect> cbc1;
            bm::compressed_buffer_collection<bvect> cbc2;
            
            GenerateCompressedBufferCollection(cbc1);

            bm::compressed_buffer_collection<bvect>::statistics st;
            cbc1.calc_stat(&st);
            
            bm::compressed_collection_serializer<compressed_buffer_collection<bvect> > cbcs;
            bm::compressed_buffer_collection<bvect>::buffer_type sbuf;
    
            cbcs.serialize(cbc1, sbuf);
            
            bm::compressed_buffer_collection<bvect>::buffer_type sbuf2(sbuf);
            sbuf.release();

            compressed_collection_deserializer<compressed_buffer_collection<bvect> > cbcd;
            cbcd.deserialize(cbc2, sbuf2.buf());
            
            if (!cbc2.equal(cbc1))
            {
                std::cerr << "Compressed collection serialization error at step " << i << endl;
                exit(1);
            }


            cout << "\r" << i << " of " << test_count << flush;
        } // for
        
        cout << endl;
    }
    
    
    cout << "------------------------ Compressed collection Test OK" << endl;
}

static
void TestBlockLast()
{
    cout << " ------------------------------ Test bit-block LAST find" << endl;
    
    {
        bool found;
        unsigned last;
        
        BM_DECLARE_TEMP_BLOCK(tb);
        for (unsigned i = 0; i < bm::set_block_size; ++i)
        {
            tb.b_.w32[i] = 0u;
        }
        found = bm::bit_find_last(tb, &last);
        assert(!found);
        
        tb.b_.w32[0] = 1u;
        found = bm::bit_find_last(tb, &last);
        assert(found);
        assert(last == 0);
        
        for (unsigned j = 0; j < 31; ++j)
        {
            tb.b_.w32[0] = 1u << j;
            found = bm::bit_find_last(tb, &last);
            assert(found);
            assert(last == j);
        }
        tb.b_.w32[0] = 0;
        for (unsigned j = 0; j < 31; ++j)
        {
            tb.b_.w32[0] |= 1u << j;
            found = bm::bit_find_last(tb, &last);
            //cout << "last = " << last << " j = " << j << endl;
            assert(found);
            assert(last == j);
        }
        
        tb.b_.w32[1] = 1u;
        found = bm::bit_find_last(tb, &last);
        cout << "last = " << last << endl;
        assert(found);
        assert(last == 32);

        tb.b_.w32[1] = 1u << 1;
        found = bm::bit_find_last(tb, &last);
        cout << "last = " << last << endl;
        assert(found);
        assert(last == 33);


        tb.b_.w32[bm::set_block_size-1] = 1u << 31;
        found = bm::bit_find_last(tb, &last);
        cout << "last = " << last << endl;
        assert(found);
        assert(last == 65535);

        tb.b_.w32[bm::set_block_size-1] = 1u << 30;
        found = bm::bit_find_last(tb, &last);
        //cout << "last = " << last << " j = " << j << endl;
        assert(found);
        assert(last == 65534);
    }
    cout << "Unit 1 ok." << endl;
    
    {
        bool found;
        unsigned last;
        
        BM_DECLARE_TEMP_BLOCK(tb);
        for (unsigned i = 0; i < bm::set_block_size; ++i)
        {
            tb.b_.w32[i] = 0u;
        }
        for (unsigned i = 0; i < bm::set_block_size; ++i)
        {
            tb.b_.w32[i] = 1u;
            found = bm::bit_find_last(tb, &last);
            assert(found);
            assert(last == (i * 32));
        }
        for (unsigned i = 0; i < bm::set_block_size; ++i)
        {
            tb.b_.w32[i] = 0u;
        }
        for (unsigned i = 0; i < bm::set_block_size; ++i)
        {
            tb.b_.w32[i] = 2;
            found = bm::bit_find_last(tb, &last);
            assert(found);
            assert(last == (i * 32)+1);
        }
    }
    cout << "Unit 2 ok." << endl;

    cout << " ------------------------------ Test bit-block LAST find OK" << endl;
}


static
void TestBlockZero()
{
    cout << " ------------------------------ Test bit-block ZERO" << endl;
    {
        BM_DECLARE_TEMP_BLOCK(tb1);

        unsigned pad = 0xDEAD;
        for (unsigned i = 0; i < bm::set_block_size; ++i)
        {
            tb1.b_.w32[i] = 0;
        }

        auto zero = bm::bit_is_all_zero(tb1);
        assert(zero);
        cout << zero << endl;

        for (unsigned i = 0; i < bm::set_block_size; ++i)
        {
            ::memset(tb1, 0, sizeof(tb1));
            tb1.b_.w32[i] = 1;
            zero = bm::bit_is_all_zero(tb1);
            assert(!zero);
            cout << zero;
        }
        cout << pad << endl;
    }
    cout << "\n ------------------------------ Test bit-block ZERO OK" << endl;

}

static
void TestFindBlockDiff()
{
    cout << " ------------------------------ Test bit_find_first_diff()" << endl;
    {
        unsigned pos;
        bool f;
        BM_DECLARE_TEMP_BLOCK(tb1);
        BM_DECLARE_TEMP_BLOCK(tb2);

        for (unsigned i = 0; i < bm::set_block_size; ++i)
        {
            tb1.b_.w32[i] = 0; tb2.b_.w32[i] = 0;
        }

        f = bm::bit_find_first_diff(tb1, tb2, &pos);
        assert(f ==  false);

        f = bm::bit_find_first(tb1, &pos);
        assert(f ==  false);

        tb2.b_.w32[0] = 1;
        f = bm::bit_find_first_diff(tb1, tb2, &pos);
        assert(f);
        assert(pos == 0);
        f = bm::bit_find_first(tb2, &pos);
        assert(f);
        assert(pos == 0);


        tb2.b_.w32[0] = (1 << 1);
        f = bm::bit_find_first_diff(tb1, tb2, &pos);
        assert(f);
        assert(pos == 1);
        f = bm::bit_find_first(tb2, &pos);
        assert(f);
        assert(pos == 1);

        tb1.b_.w32[0] = (1 << 1);
        f = bm::bit_find_first_diff(tb1, tb2, &pos);
        assert(!f);

        for (unsigned i = 0; i < bm::set_block_size; ++i)
        {
            tb1.b_.w32[i] = 0; tb2.b_.w32[i] = 0;
        }
        tb1.b_.w32[0] = (1<<12); tb2.b_.w32[1] = 18;
        f = bm::bit_find_first_diff(tb1, tb2, &pos);
        assert(f);
        assert(pos == 12);
        f = bm::bit_find_first(tb1, &pos);
        assert(f);
        assert(pos == 12);

        for (unsigned i = 0; i < bm::set_block_size; ++i)
        {
            tb1.b_.w32[i] = 0; tb2.b_.w32[i] = 0;
        }

        bm::set_bit(tb1, 12345);
        f = bm::bit_find_first_diff(tb1, tb2, &pos);
        assert(f);
        assert(pos == 12345);
        f = bm::bit_find_first(tb1, &pos);
        assert(f);
        assert(pos == 12345);

        bm::set_bit(tb2, 12345);
        f = bm::bit_find_first_diff(tb1, tb2, &pos);
        assert(!f);

        for (unsigned i = 0; i < bm::set_block_size; ++i)
        {
            tb1.b_.w32[i] = 0; tb2.b_.w32[i] = 0;
        }

        for (unsigned k = 0; k < 65535; ++k)
        {
            bm::set_bit(tb1, k);
            f = bm::bit_find_first_diff(tb1, tb2, &pos);
            assert(f);
            assert(pos == k);


            for (unsigned j = k+1; j < 65535; ++j)
            {
                bm::set_bit(tb1, j);
                f = bm::bit_find_first_diff(tb1, tb2, &pos);
                assert(f);
                assert(pos == k);
                bm::set_bit(tb2, j);
                f = bm::bit_find_first_diff(tb1, tb2, &pos);
                assert(f);
                assert(pos == k);
                bm::clear_bit(tb1, j);
                bm::clear_bit(tb2, j);
            }

            bm::set_bit(tb2, k);
            f = bm::bit_find_first_diff(tb1, tb2, &pos);
            assert(!f);

            cout << "\r" << k << flush;
        } // for

    }
    cout << "\n ------------------------------ Test bit_find_first_diff() OK" << endl;

}



static
void TestBlockDigest()
{
    cout << " ------------------------------ Test bit-block DIGEST" << endl;

    BM_DECLARE_TEMP_BLOCK(tb1);

    {
        bm::id64_t digest0 = 1;
        bool all_zero;
        all_zero = bm::check_zero_digest(digest0, 0, 0);
        assert(!all_zero);
        for (unsigned i = 0; i < bm::set_block_digest_wave_size * 32; ++i)
        {
            all_zero = bm::check_zero_digest(digest0, 0, i);
            assert(!all_zero);
        }
        for (unsigned i = bm::set_block_digest_wave_size * 32+1; i < 65535; ++i)
        {
            all_zero = bm::check_zero_digest(digest0,
                                        bm::set_block_digest_wave_size * 32+1,
                                        i);
            assert(all_zero);
            
        } // for
    }

    for (unsigned k = 0; k < bm::set_block_size; ++k)
    {
        bm::bit_block_set(tb1, 0);

        tb1.b_.w32[k] = 1;
        bm::id64_t mask1 = bm::widx_to_digest_mask(k);
        bm::id64_t mask2 = bm::calc_block_digest0(tb1);
        assert(mask1 == mask2);
        bm::id64_t mask3 = bm::update_block_digest0(tb1, mask1);
        assert(mask1 == mask3);
        
        assert(mask1);
        unsigned bc = bm::word_bitcount64(mask1);
        assert(bc == 1);

        unsigned first_bit, ffbc;
        bool single_bit_found = bm::bit_find_first_if_1(tb1, &first_bit, mask1);
        assert(single_bit_found);
        unsigned found = bm::bit_find_first(tb1, &ffbc);
        assert(found);
        assert(first_bit == ffbc);

        bm::bit_block_set(tb1, 0);
        mask3 = bm::update_block_digest0(tb1, mask1);
        assert(!mask3);

        single_bit_found = bm::bit_find_first_if_1(tb1, &first_bit, mask1);
        assert(!single_bit_found);

        single_bit_found = bm::bit_find_first_if_1(tb1, &first_bit, mask3);
        assert(!single_bit_found);

        tb1.b_.w32[k] = 3;
        mask3 = bm::calc_block_digest0(tb1); 

        single_bit_found = bm::bit_find_first_if_1(tb1, &first_bit, mask3);
        assert(!single_bit_found);

        tb1.b_.w32[0] = 1;
        tb1.b_.w32[bm::set_block_size/2] = 2;

        mask3 = bm::calc_block_digest0(tb1);

        single_bit_found = bm::bit_find_first_if_1(tb1, &first_bit, mask3);
        assert(!single_bit_found);
    }

    unsigned start = 0;
    unsigned end = bm::set_block_size-1;
    
    while(start <= end)
    {
        bm::bit_block_set(tb1, 0);

        tb1.b_.w32[start] = 1;
        tb1.b_.w32[end] = 1;
        
        bm::id64_t mask_s1 = bm::widx_to_digest_mask(start);
        bm::id64_t mask_e1 = bm::widx_to_digest_mask(end);
        bm::id64_t mask1 = mask_s1 | mask_e1;
        bm::id64_t mask2 = bm::calc_block_digest0(tb1);
        assert(mask1 == mask2);
        bm::id64_t mask3 = bm::update_block_digest0(tb1, mask1);
        assert(mask1 == mask3);
        
        if (mask_s1 != mask_e1)
        {
            unsigned bc = bm::word_bitcount64(mask1);
            assert(bc == 2);
        }
        else
        {
            unsigned bc = bm::word_bitcount64(mask1);
            assert(bc == 1);
        }
        bm::bit_block_set(tb1, 0);
        mask3 = bm::update_block_digest0(tb1, mask1);
        assert(!mask3);
        

        ++start; --end;
    } // while

    cout << " ------------------------------ Test bit-block DIGEST OK" << endl;
}


static
void TestBlockAND()
{
    cout << " ------------------------------ Test bit-block AND" << endl;
    {
        BM_DECLARE_TEMP_BLOCK(tb2);
        BM_DECLARE_TEMP_BLOCK(tb1);
        BM_DECLARE_TEMP_BLOCK(tb0);

        unsigned pad = 0xDEAD;
        for (unsigned i = 0; i < bm::set_block_size; ++i)
        {
            tb1.b_.w32[i] = tb2.b_.w32[i] = 0;
        }

        auto any = bm::bit_block_and(tb1, tb2);
        assert(any == 0);
        tb1.b_.w32[1] = 1;
        any = bm::bit_block_and(tb1, tb2);
        assert(any == 0);
        assert(tb1.b_.w32[1] == 0);
        
        
        tb1.b_.w32[1] = tb2.b_.w32[1] = 1;
        any = bm::bit_block_and(tb1, tb2);
        
        cout << tb1.b_.w32[1] << endl;
        assert(tb1.b_.w32[1] == 1);
        assert(any);
        for (unsigned j = 0; j < 32; ++j)
        {
            tb1.b_.w32[10] = tb2.b_.w32[10] = (1 << j);
            if (tb1.b_.w32[10])
            {
                any = bm::bit_block_and(tb1, tb2);
                assert(any);
            }
        }
        
        for (unsigned i = 0; i < bm::set_block_size; ++i)
        {
            tb1.b_.w32[i] = tb2.b_.w32[i] = 8;
        }
        any = bm::bit_block_and(tb1, tb2);
        assert(any);
        {
//            bm::bit_decode_cache dcache;
            bm::id64_t d1 = ~0ull;
            d1 = bm::bit_block_and(tb1,
                                   tb2,
                                   d1);
            bm::id64_t dc = bm::calc_block_digest0(tb1);
            assert(d1 == dc);
            unsigned bc = bm::word_bitcount64(d1);
            assert(bc == 64);
        }
        for (unsigned i = 0; i < bm::set_block_size; ++i)
        {
            assert(tb1.b_.w32[i] == tb2.b_.w32[i]);
        }


        for (unsigned i = 0; i < bm::set_block_size; ++i)
        {
            tb1.b_.w32[i] = tb2.b_.w32[i] = 0;
        }

        unsigned i, j;
        for (i = 0; i < bm::set_block_size; ++i)
        {
            for (j = 0; j < 32; ++j)
            {
                unsigned v = (1u << j);
                ::memset(tb1, 0, sizeof(tb1));
                ::memset(tb2, 0, sizeof(tb1));
                tb1.b_.w32[i] = tb2.b_.w32[i] = v;
                if (tb1[i])
                {
                    auto any1 = bm::bit_block_and(tb1, tb2);
                    auto all_zero = bm::bit_is_all_zero(tb1.begin());
                    
                    //cout << any1 <<" j=" << j << " i=" << i << " " << tb1[i] << " " << tb2[i] << endl;
                    assert(pad == 0xDEAD);
                    assert(tb1.b_.w32[i] == v);
                    assert((unsigned)(all_zero) != any1);
                    assert(any1);
                    {
                        bm::id64_t d1 = ~0ull;
                        d1 = bm::bit_block_and(tb1,
                                               tb2,
                                               d1);
                        bm::id64_t dc = bm::calc_block_digest0(tb1);
                        assert(d1 == dc);
                        unsigned bc = bm::word_bitcount64(d1);
                        assert(bc == 1);
                    }
                }
            }
            tb1.b_.w32[i] = tb2.b_.w32[i] = 0;
        }
        cout << tb1.b_.w32[0] << pad << endl;


        for (i = 0; i < bm::set_block_size; ++i)
        {
            ::memset(tb1, 0, sizeof(tb1));
            ::memset(tb2, 0, sizeof(tb1));
            
            tb1.b_.w32[i] = tb2.b_.w32[i] = 8u;

            auto any1 = bm::bit_block_and(tb1, tb2);
            assert(tb1.b_.w32[i] == 8u);
            assert(any1);

            ::memset(tb1, 0, sizeof(tb1));
            ::memset(tb2, 0, sizeof(tb1));
            
            tb1.b_.w32[i] = tb2.b_.w32[i] = 8u;

            bm::id64_t d1 = ~0ull;
            d1 = bm::bit_block_and(tb1,
                                   tb2,
                                   d1);
            bm::id64_t dc = bm::calc_block_digest0(tb1);
            assert(d1 == dc);
            unsigned bc = bm::word_bitcount64(d1);
            assert(bc == 1);
            
            d1 = ~0ull;
            d1 = bm::bit_block_and_2way(tb0, tb1, tb2, d1);
            assert(d1 == dc);
        }
        
        

    }
    cout << " ------------------------------ Test bit-block AND  OK" << endl;

}

static
void TestBlockOR()
{
    BM_DECLARE_TEMP_BLOCK(tb3);
    BM_DECLARE_TEMP_BLOCK(tb2);
    BM_DECLARE_TEMP_BLOCK(tb1);
    
    bool all_one;

    cout << " ------------------------------ Test bit-block OR" << endl;

    {
        for (unsigned i = 0; i < bm::set_block_size; ++i)
        {
            tb1.b_.w32[i] = tb3.b_.w32[i] = 0;
            tb2.b_.w32[i] = 8;
        }

        all_one = bm::bit_block_or(tb1, tb2);
        assert(!all_one);
        all_one = bm::bit_block_or_3way(tb3, tb2, tb1);
        assert(!all_one);


        for (unsigned i = 0; i < bm::set_block_size; ++i)
        {
            assert(tb1.b_.w32[i] == 8);
            if (tb1.b_.w32[i] != 8 || tb3.b_.w32[i] != 8)
            {
                cerr << "TestOR failed!" << endl;
                exit(1);
            }
        }
    }
    
    {
        for (unsigned i = 0; i < bm::set_block_size; ++i)
        {
            tb1.b_.w32[i] = ~3u;
            tb2.b_.w32[i] = 3u;
            tb3.b_.w32[i] = 0;
        }

        all_one = bm::bit_block_or(tb1, tb2);
        assert(all_one);
        all_one = bm::bit_block_or_3way(tb3, tb2, tb1);
        assert(all_one);


        for (unsigned i = 0; i < bm::set_block_size; ++i)
        {
            assert(tb1.b_.w32[i] == ~0u);
            assert(tb3.b_.w32[i] == ~0u);
        }
    }

    {
        for (unsigned i = 0; i < bm::set_block_size; ++i)
        {
            tb1.b_.w32[i] = 0;
            tb2.b_.w32[i] = 0;
            tb3.b_.w32[i] = 0;
        }
        for (unsigned i = 0; i < 100; ++i)
        {
            tb1.b_.w32[i] = ~0u;
            tb2.b_.w32[i] = 0;
            tb3.b_.w32[i] = 0;
        }
        for (unsigned i = 100; i < bm::set_block_size; ++i)
        {
            tb1.b_.w32[i] = 0;
            tb2.b_.w32[i] = ~0u;
            tb3.b_.w32[i] = 0;
        }


        all_one = bm::bit_block_or(tb1, tb2);
        assert(all_one);
        all_one = bm::bit_block_or_3way(tb3, tb2, tb1);
        assert(all_one);


        for (unsigned i = 0; i < bm::set_block_size; ++i)
        {
            assert(tb1.b_.w32[i] == ~0u);
            assert(tb3.b_.w32[i] == ~0u);
        }
    }

    
    cout << " ------------------------------ Test bit-block OR  OK" << endl;

}

static
void TestBlockSUB()
{
    cout << " ------------------------------ Test bit-block SUB" << endl;
    {
        BM_DECLARE_TEMP_BLOCK(tb2);
        BM_DECLARE_TEMP_BLOCK(tb1);

        unsigned pad = 0xDEAD;
        for (unsigned i = 0; i < bm::set_block_size; ++i)
        {
            tb1.b_.w32[i] = tb2.b_.w32[i] = 0;
        }

        auto any = bm::bit_block_sub(tb1, tb2);
        assert(any == 0);
        tb1.b_.w32[1] = 1;
        any = bm::bit_block_sub(tb1, tb2);
        assert(any);
        assert(tb1.b_.w32[1] == 1);
        
        
        tb1.b_.w32[1] = tb2.b_.w32[1] = 1;
        any = bm::bit_block_sub(tb1, tb2);
        
        cout << tb1.b_.w32[1] << endl;
        assert(tb1.b_.w32[1] == 0);
        assert(!any);
        for (unsigned j = 0; j < 32; ++j)
        {
            tb1.b_.w32[10] = tb2.b_.w32[10] = (1 << j);
            if (tb1.b_.w32[10])
            {
                any = bm::bit_block_sub(tb1, tb2);
                assert(!any);
            }
        }
        
        for (unsigned i = 0; i < bm::set_block_size; ++i)
        {
            tb1.b_.w32[i] = tb2.b_.w32[i] = 8;
        }
        any = bm::bit_block_sub(tb1, tb2);
        assert(!any);
        for (unsigned i = 0; i < bm::set_block_size; ++i)
        {
            assert(tb1.b_.w32[i] != tb2.b_.w32[i]);
        }

        for (unsigned i = 0; i < bm::set_block_size; ++i)
        {
            tb1.b_.w32[i] = tb2.b_.w32[i] = 0;
        }

        unsigned i, j;
        for (i = 0; i < bm::set_block_size; ++i)
        {
            for (j = 0; j < 32; ++j)
            {
                unsigned v = (1u << j);
                ::memset(tb1, 0, sizeof(tb1));
                ::memset(tb2, 0, sizeof(tb1));
                tb1.b_.w32[i] = tb2.b_.w32[i] = v;
                if (tb1[i])
                {
                    auto any1 = bm::bit_block_sub(tb1, tb2);
                    auto all_zero = bm::bit_is_all_zero(tb1.begin());
                    assert(all_zero);
                    
                    //cout << any1 <<" j=" << j << " i=" << i << " " << tb1[i] << " " << tb2[i] << endl;
                    assert(pad == 0xDEAD);
                    assert(tb1.b_.w32[i] == 0);
                    assert(!any1);
                }
            }
            tb1.b_.w32[i] = tb2.b_.w32[i] = 0;
        }
        cout << tb1.b_.w32[0] << pad << endl;
        
        
        for (i = 0; i < bm::set_block_size; ++i)
        {
            ::memset(tb1, 0, sizeof(tb1));
            ::memset(tb2, 0, sizeof(tb1));
            
            tb1.b_.w32[i] = tb2.b_.w32[i] = 8u;

            auto any1 = bm::bit_block_sub(tb1, tb2);
            assert(tb1.b_.w32[i] == 0);
            assert(!any1);
        }

    }
    cout << " ------------------------------ Test bit-block SUB  OK" << endl;

}


static
void TestRankCompress()
{
    cout << " ------------------------------ Test Rank Compressor " << endl;
    
    int cmp;

    cout << "Step 1" << endl;
    {
        bvect bv1, bv2;
        bvect bv_s { 0, 1,        16 };
        bvect bv_i { 0, 1, 2, 10, 16 };
        bvect bv_sr; // restored vector

        bvect bv_ref { 0, 1, 4 };
        bm::rank_compressor<bvect> rc;

        bvect::rs_index_type bc;
        bv_i.build_rs_index(&bc);

        for (unsigned i = 0; i < 2; ++i)
        {
            rc.compress(bv1, bv_i, bv_s);
            assert(bv1.count() == bv_s.count());
            
            cmp = bv1.compare(bv_ref);
            assert(cmp == 0);
            
            rc.decompress(bv_sr, bv_i, bv1);
            cmp = bv_sr.compare(bv_s);
            assert(cmp == 0);
            
            
            rc.compress_by_source(bv2, bv_i, bc, bv_s);
            assert(bv2.count() == bv_s.count());

            cmp = bv2.compare(bv_ref);
            assert(cmp == 0);
            
            
            bv_i.optimize();
            bv_s.optimize();
        }
    }
    cout << "Step 1 - OK" << endl;

    {
        bvect bv1, bv2;
        bvect bv_s { 0, 100000, 100001,                  1600000, 1600001  };
        bvect bv_i { 0, 100000, 100001, 200000, 1000000, 1600000, 1600001 };
        bvect bv_sr;

        bm::rank_compressor<bvect> rc;

        bvect::rs_index_type bc;
        bv_i.build_rs_index(&bc);

        for (unsigned i = 0; i < 2; ++i)
        {
            rc.compress(bv1, bv_i, bv_s);
            assert(bv1.count() == bv_s.count());

            rc.decompress(bv_sr, bv_i, bv1);
            cmp = bv_sr.compare(bv_s);
            if (cmp != 0)
            {
                DetailedCompareBVectors(bv_sr, bv_s);
            }
            assert(cmp == 0);
 
            rc.compress_by_source(bv2, bv_i, bc, bv_s);
            assert(bv2.count() == bv_s.count());

            cmp = bv2.compare(bv1);
            assert(cmp == 0);

            rc.decompress(bv_sr, bv_i, bv2);
            cmp = bv_sr.compare(bv_s);
            assert(cmp == 0);

 
            bv_i.optimize();
            bv_s.optimize();
        }
    }
    std::cout << "basic test OK..." << std::endl;


    {
        std::cout << "Stress rank compression..." << std::endl;
        bm::rank_compressor<bvect> rc;
        unsigned test_count = 10;
        unsigned bv_size = 1000000;
        for (unsigned i  = 0; i < test_count; ++i)
        {
            if (bv_size > 40000000)
                break;
            cout << "target size = " << bv_size << " " << endl;
            bvect bv_i, bv_s, bv_sr;
            generate_bvector(bv_i, bv_size);
            generate_bvector(bv_s, bv_size);
            bv_i |= bv_s;
            
            assert(bv_i.count() >= bv_s.count());

            bvect::rs_index_type bc;
            bv_i.build_rs_index(&bc);
            
            // quick rank test
            //
            /*
            bm::id_t pos = 5308470;
            bm::id_t r1 = bv_i.count_range(0, pos)-1;
            bm::id_t r2 = bv_i.count_to(pos, bc)-1;
            assert(r1 == r2);
            cout << "i=" << pos << " rank()=" << r1 << endl;
            */

            bvect bv1, bv2;
            
            for (unsigned j = 0; j < 2; ++ j)
            {
                {
                chrono_taker ct("c1");
                rc.compress(bv1, bv_i, bv_s);
                rc.decompress(bv_sr, bv_i, bv1);
                }
                assert(bv1.count() == bv_s.count());
                cmp = bv_sr.compare(bv_s);
                assert(cmp == 0);

                {
                chrono_taker ct("c2");
                rc.compress_by_source(bv2, bv_i, bc, bv_s);
                rc.decompress(bv_sr, bv_i, bv2);
                }
                assert(bv2.count() == bv_s.count());
                cmp = bv_sr.compare(bv_s);
                assert(cmp == 0);

                cmp = bv2.compare(bv1);
                if (cmp!=0)
                {
                    DetailedCompareBVectors(bv1, bv2);
                    exit(1);
                }
                assert(cmp == 0);

                {
                    bm::random_subset<bvect> rsub;
                    bvect bv_subset;
                    rsub.sample(bv_subset, bv_s, 100);
                    {
                    chrono_taker ct("c1-1");
                    rc.compress(bv1, bv_i, bv_subset);
                    rc.decompress(bv_sr, bv_i, bv1);
                    }
                    assert(bv1.count() == bv_subset.count());
                    cmp = bv_sr.compare(bv_subset);
                    assert(cmp == 0);

                    {
                    chrono_taker ct("c2-2");
                    rc.compress_by_source(bv2, bv_i, bc, bv_subset);
                    rc.decompress(bv_sr, bv_i, bv2);
                    }
                    assert(bv2.count() == bv_subset.count());
                    
                    cmp = bv2.compare(bv1);
                    assert(cmp == 0);

                    cmp = bv_sr.compare(bv_subset);
                    assert(cmp == 0);
                }


                bv_i.optimize();
                bv_s.optimize();
            } // for j
            cout << "\n" << i << " of " << test_count << "  " << endl;
            
            bv_size += bv_size;
        } // for i
        std::cout << endl << "Stress rank compression... OK" << std::endl;
    }
    
    
    cout << " ------------------------------ Test Rank Compressor OK " << endl;
}


static
void GenerateSV(sparse_vector_u32&   sv, unsigned strategy = 0)
{
    unsigned max_idx_value = 1000000;
    switch (strategy)
    {
    case 0:
    {
        cout << "SV Ultra sparse generation" << endl;
        for (unsigned i = 0; i < max_idx_value;)
        {
            unsigned v = (rand() * rand()) % 650000;
            sv[i] = v;
            i += 10000 + rand() % 65535;
        }
        break;
    }
    case 1:
    {
        cout << "SV Dense intervals generation 1" << endl;
        for (unsigned i = 0; i < max_idx_value;)
        {
            unsigned v = (rand() * rand()) % 650000;
            for (unsigned j = 0; i < max_idx_value; ++i, ++j)
            {
                sv[i] = v + j;
                if (j > 256)
                    break;
            }
            i += 20000 + rand() % 65535;
        }
        break;
    }
    case 2:
    {
        cout << "SV Dense intervals generation 2" << endl;
        unsigned v = (rand() * rand()) % 650000;
        for (unsigned i = 0; i < max_idx_value/4; ++i)
        {
            sv[i] = v;
        }

        for (unsigned i = 0; i < max_idx_value;)
        {
            v = unsigned(rand() * rand()) % 650000;
            for (unsigned j = 0; i < max_idx_value; ++i, ++j)
            {
                sv[i] = v + i;
                if (j > 256)
                    break;
            }
            i += 30000 + unsigned(rand()) % 65535;
        }
        break;
    }
    case 3:
    {
        cout << "SV random generation" << endl;
        unsigned rand_max = rand() % 300000;
        for (unsigned i = 0; i < rand_max; ++i)
        {
            unsigned v = unsigned(rand() * rand());
            unsigned idx = unsigned(rand()) % max_idx_value;
            sv[idx] = v;
            if (i % 2 == 0)
            {
                sv.clear(idx, true);
            }
        }
        break;
    }
    case 4:
        {
        cout << "SV empty generation" << endl;
        unsigned idx = unsigned(rand()) % max_idx_value;
        sv[idx] = 25557890;
        sv.clear(idx, true);
        }
        break;
    case 5:
        {
        cout << "SV uniform power 2 value generation" << endl;
        unsigned v = 8;//unsigned(rand()) % 64;
        for (unsigned i = 0; i < max_idx_value; ++i)
        {
            sv[i] = v;
        }
        }
        break;
    case 6:
        {
        cout << "SV uniform power 2+1 value generation" << endl;
        unsigned v = 16+1;
        for (unsigned i = 0; i < max_idx_value; ++i)
        {
            sv[i] = v;
        }
        }
        break;
    case 7:
        {
        cout << "SV liner growth value generation" << endl;
        for (unsigned i = 0; i < max_idx_value; ++i)
        {
            sv[i] = i;
        }
        }
        break;
    default:
        break;
    } // switch
    sv.optimize();
}

static
void DetailedCompareSparseVectors(const rsc_sparse_vector_u32& csv,
                                  const sparse_vector_u32&     sv)
{
    sparse_vector_u32   sv_s(bm::use_null);  // sparse vector decompressed
    
    // de-compression test
    csv.load_to(sv_s);
    /*
    if (!sv.equal(sv_s))
    {
        cerr << "compressed vector load_to (decompression) failed!" << endl;
        exit(1);
    }
    */
    

    size_t csv_size = csv.size();
    size_t sv_size = sv.size();
    size_t sv_s_size = sv_s.size();

    const sparse_vector_u32::bvector_type* bv_null_sv = sv.get_null_bvector();
    const sparse_vector_u32::bvector_type* bv_null_sv_s = sv_s.get_null_bvector();
    const sparse_vector_u32::bvector_type* bv_null_csv = csv.get_null_bvector();

    if (csv_size != sv_size || sv_s_size != sv_size)
    {
        assert(bv_null_sv != bv_null_csv);
        
        unsigned cnt_sv = bv_null_sv->count();
        unsigned cnt_sv_s = bv_null_sv_s->count();
        unsigned cnt_csv = bv_null_csv->count();
        
        if (cnt_sv != cnt_csv)
        {
            cerr << "Sparse compressed vector comparison failed (size check):"
                 << "csv.size()=" << csv_size
                 << "sv.size()=" << sv_size
                 << "cnt sv = " << cnt_sv
                 << "cnt csv = " << cnt_csv
                 << endl;
            exit(1);
        }
        if (cnt_sv_s != cnt_csv)
        {
            cerr << "Restored Sparse vector comparison failed (size check):"
                 << "csv.size()=" << csv_size
                 << "sv_s.size()=" << sv_s_size
                 << "cnt sv = " << cnt_sv
                 << "cnt csv = " << cnt_csv
                 << endl;
            exit(1);
        }
    }
    
    for (unsigned i = 0; i < sv_size; ++i)
    {
        bool is_null_sv = sv.is_null(i);
        bool is_null_sv_s = sv_s.is_null(i);
        bool is_null_csv = csv.is_null(i);
        if (is_null_sv != is_null_csv || is_null_sv != is_null_sv_s)
        {
            cerr << "Detailed csv check failed (null mismatch) at i=" << i
                << " sv=" << is_null_sv
                << " sv_s=" << is_null_sv_s
                << " csv=" << is_null_csv
                << endl;
            int cmp = bv_null_sv->compare(*bv_null_csv);
            if (cmp != 0)
            {
                cerr << "1. cmp=" << cmp << endl;
                exit(1);
            }
            cmp = bv_null_sv->compare(*bv_null_sv_s);
            if (cmp != 0)
            {
                cerr << "2. cmp=" << cmp << endl;
                exit(1);
            }

            exit(1);
        }
        
        
        if (!is_null_sv)
        {
            unsigned v1 = sv[i];
            unsigned v1_s = sv_s[i];
            unsigned v2 = csv[i];
            
            if (v1 != v2 || v1_s != v1)
            {
                cerr << "Detailed csv check failed (value mismatch) at i=" << i
                     << " v1=" << v1
                     << " v1_s=" << v1_s
                     << " v2=" << v2
                     << endl;
                exit(1);
            }
        }
    }
    
    {
    BM_DECLARE_TEMP_BLOCK(tb)
    sparse_vector_serial_layout<rsc_sparse_vector_u32> sv_lay;
    bm::sparse_vector_serialize<rsc_sparse_vector_u32>(csv, sv_lay, tb);

    rsc_sparse_vector_u32 csv1;
    const unsigned char* buf = sv_lay.buf();
    bm::sparse_vector_deserialize(csv1, buf, tb);

    if (!csv.equal(csv1))
    {
        cerr << "Compressed sparse vector serialization comparison failed!" << endl;

        rsc_sparse_vector_u32::size_type pos;
        bool b = bm::sparse_vector_find_first_mismatch(csv, csv1, pos);
        assert(b);
        cerr << "Mismatch at: " << pos << endl;

        sparse_vector_serial_layout<rsc_sparse_vector_u32> sv_lay1;
        bm::sparse_vector_serialize<rsc_sparse_vector_u32>(csv, sv_lay1);

        bm::sparse_vector_deserialize(csv1, buf, tb);

        exit(1);
    }
    }
    
}

static
void CheckCompressedDecode(const rsc_sparse_vector_u32& csv,
                           unsigned from, unsigned size)
{
    std::vector<unsigned> vect;
    vect.resize(size);
    
    unsigned sz = csv.decode(&vect[0], from, size);
    unsigned ex_idx = 0;
    for (unsigned i = from; i < from + sz; ++i)
    {
        unsigned v = csv.get(i);
        unsigned vx = vect[ex_idx];
        if (v != vx)
        {
            cerr << "compressed vector decode mismatch from="
                 << from << " idx=" << i
                 << " v=" << v << " vx=" << vx
                 << endl;
            exit(1);
        }
        ++ex_idx;
    }
}

static
void DetailedCheckCompressedDecode(const rsc_sparse_vector_u32& csv)
{
    auto size = csv.size();
    cout << endl;

    {
    unsigned size1 = 100;
    for (unsigned i = 0; i < size1; )
    {
        CheckCompressedDecode(csv, i, size);
        if (i % 128 ==0)
            cout << "\r" << i << "/" << size1 << flush;
        i++;
    }
    }
    cout << endl;

    {
    unsigned size1 = 100000;
    for (unsigned i = 0; i < size1; )
    {
        CheckCompressedDecode(csv, i, size1);
        cout << "\r" << i << "/" << size1 << flush;
        i+=rand()%3;
        size1 -= rand()%5;
    }
    }
    cout << endl;

    {
    unsigned size1 = size;
    for (unsigned i = size-size/2; i < size1; )
    {
        CheckCompressedDecode(csv, i, size1);
        cout << "\r" << i << "/" << size1 << flush;
        i+=(1+i);
    }
    }
    cout << endl;

    for (unsigned i = size-size/2; i < size; )
    {
        CheckCompressedDecode(csv, i, size);
        cout << "\r" << i << "/" << size << flush;
        i += rand() % 25000;
    }
    cout << endl;
    
    for (unsigned i = size-size/2; i < size; )
    {
        if (size <= i)
            break;
        CheckCompressedDecode(csv, i, size);
        cout << "\r" << i << "/" << size << flush;
        i += rand() % 25000;
        size -= rand() % 25000;;
    }
    cout << endl;

}

static
void TestArraysAndBuffers()
{
    cout << " ------------------------------ Test buffers " << endl;
    
    typedef bm::heap_matrix<unsigned, 10, 20, bvect::allocator_type> hmatrix;
    typedef bm::heap_matrix<char, 10, 256, bvect::allocator_type> hmatrix_str;

    {
        hmatrix hm;
        hm.init();
        hm.set_zero();

        for (unsigned i = 0; i < hm.rows(); ++i)
        {
            const unsigned* r = hm.row(i);
            for (unsigned j = 0; j < hm.cols(); ++j)
            {
                assert(r[j] == 0);
            }
        }
    }

    {
        hmatrix hm(true);

        for (unsigned i = 0; i < hm.rows(); ++i)
        {
            unsigned* r = hm.row(i);
            for (unsigned j = 0; j < hm.cols(); ++j)
            {
                r[j] = i;
            }
        }

        for (unsigned i = 0; i < hm.rows(); ++i)
        {
            const unsigned* r = hm.row(i);
            for (unsigned j = 0; j < hm.cols(); ++j)
            {
                assert(r[j] == i);
            }
        }
    }

    {
        hmatrix_str hm(true);
        hm.set_zero();

        for (unsigned i = 0; i < hm.rows(); ++i)
        {
            char* r = hm.row(i);
            ::strncpy(r, "abcd", hm.cols());
        }
        for (unsigned i = 0; i < hm.rows(); ++i)
        {
            const char* r = hm.row(i);
            int c = ::strcmp(r, "abcd");
            assert(c == 0);
        }

    }


    cout << " ------------------------------ Test buffers OK" << endl;
}

static
void TestCompressSparseVector()
{
    cout << " ------------------------------ Test Compressed Sparse Vector " << endl;
    
    {
        rsc_sparse_vector_u32 csv1;
        assert(csv1.size() == 0);
        assert(csv1.equal(csv1));
        rsc_sparse_vector_u32 csv2;
        assert(csv1.equal(csv2));
        rsc_sparse_vector_u32 csv3(csv1);
        assert(csv3.equal(csv2));
    }
    
    {
    cout << "push_back() test" << endl;
    unsigned v, v1;
    
        rsc_sparse_vector_u32 csv1;
        sparse_vector_u32 sv1(bm::use_null);
        
        csv1.push_back(10, 100);
        assert(csv1.size() == 11);
        csv1.push_back(20, 200);
        csv1.push_back(21, 201);
        
        csv1.load_to(sv1);

        v = csv1.at(10);
        assert(v == 100);
        v1 = sv1.at(10);
        assert(v1 == 100);
        
        v = csv1.at(20);
        assert(v == 200);
        v1 = sv1.at(20);
        assert(v1 == 200);

        v = csv1.at(21);
        assert(v == 201);
        v1 = sv1.at(21);
        assert(v1 == 201);

        csv1.sync();

        DetailedCompareSparseVectors(csv1, sv1);
        
        v = csv1.at(10);
        assert(v == 100);
        v = csv1.at(20);
        assert(v == 200);
        v = csv1.at(21);
        assert(v == 201);
        
        csv1.optimize();
        v = csv1.at(10);
        assert(v == 100);
        v = csv1.at(20);
        assert(v == 200);
        v = csv1.at(21);
        assert(v == 201);
        
        rsc_sparse_vector_u32 csv2(csv1);
        bool same = csv2.equal(csv1);
        assert(same);
        
        rsc_sparse_vector_u32 csv3;
        csv3 = ::move(csv2);
        same = csv3.equal(csv1);
        assert(same);
        
        bm::sparse_vector_scanner<rsc_sparse_vector_u32> scanner;
        bm::id_t pos;
        bool found = scanner.find_eq(csv1, 201, pos);
        assert(found);
        assert(pos == 21);

    }
    
    // back inserter tests
    {
        rsc_sparse_vector_u32 csv1;
        {
        rsc_sparse_vector_u32::back_insert_iterator rs_bi = csv1.get_back_inserter();
            rs_bi.add_null();
            rs_bi.add(1);
            rs_bi.add(2);
            rs_bi.flush();
        }
        assert(csv1.size() == 3);
        auto v = csv1.get(0);
        assert(v == 0);
        assert(csv1.is_null(0));
        v = csv1.at(1);
        assert(v == 1);
        v = csv1.at(2);
        assert(v == 2);
    }
    
    {
        rsc_sparse_vector_u32 csv1;
        {
        rsc_sparse_vector_u32::back_insert_iterator rs_bi = csv1.get_back_inserter();
            rs_bi.add(1);
            rs_bi.add(2);
            rs_bi.add_null();
            rs_bi.add(3);
            rs_bi.flush();
        }
        assert(csv1.size() == 4);
        auto v = csv1.at(0);
        assert(v == 1);
        v = csv1.at(1);
        assert(v == 2);
        v = csv1.get(2);
        assert(v == 0);
        assert(csv1.is_null(2));
        v = csv1.at(3);
        assert(v == 3);

        // test copy-range
        {
            rsc_sparse_vector_u32 csv2;
            csv2.copy_range(csv1, 4, 5);
            assert(csv2.size() == 0);

            csv2.copy_range(csv1, 0, 0);
            assert(csv2.size() == 4);
            v = csv2.at(0);
            assert(v == 1);

            csv2.copy_range(csv1, 1, 2);
            assert(csv2.size() == 4);
            v = csv2[0];
            assert(v == 0);

            v = csv2.at(1);
            assert(v == 2);
            v = csv2.get(2);
            assert(v == 0);
            assert(csv2.is_null(2));
        }

    }

    
    {
        rsc_sparse_vector_u32 csv1;
        {
        rsc_sparse_vector_u32::back_insert_iterator rs_bi = csv1.get_back_inserter();
        for (unsigned i = 0; i < 100000; i++)
        {
            if (i&1)
            {
                rs_bi.add_null();
            }
            else
            {
                rs_bi.add(i);
            }
        }
        rs_bi.flush();
        }
        csv1.optimize();
        
        // validation
        for (unsigned i = 0; i < 100000; i++)
        {
            if (i&1)
            {
                assert(csv1.is_null(i));
            }
            else
            {
                assert(!csv1.is_null(i));
                auto v = csv1[i];
                assert(v == i);
            }
        }
    }
    
    // set test
    {
        rsc_sparse_vector_u32 csv;
        csv.set(1, 1);
        assert(csv.is_null(0));
        assert(!csv.is_null(1));
        assert(csv.get(1) == 1);
        
        csv.push_back(10, 11);
        csv.set(11, 12);
        assert(csv.get(11) == 12);
        csv.set(5, 55);
        csv.set(5, 56);

        assert(csv.size() == 12);
        assert(csv.get(1) == 1);
        assert(csv.get(10) == 11);
        assert(csv.get(11) == 12);
        assert(csv.get(5) == 56);
        
        csv.set_null(5);
        assert(csv.is_null(5));
        assert(csv.get(1) == 1);
        assert(csv.get(10) == 11);
        assert(csv.get(11) == 12);
        assert(csv.get(5) == 0);
    }

    // set stress test
    {
        cout << "RSC set stress..." << endl;
        std::vector<std::pair<unsigned, unsigned> > vect;
        rsc_sparse_vector_u32 csv;
        
        const unsigned max_size = 2000000;
        
        cout << "Test set generation." << endl;
        for (unsigned i = 0; i < max_size; i+=2)
        {
            std::pair<unsigned, unsigned> pr(i, i+10);
            vect.push_back(pr);
        } // for
        
        {
            std::random_device rd;
            std::mt19937 g(rd());
            std::shuffle(vect.begin(), vect.end(), g);
        }
        
        cout << "RSC set() " << endl;
        unsigned i = 0;
        for (auto rit = vect.rbegin(); rit != vect.rend(); ++rit)
        {
            std::pair<unsigned, unsigned> pr = *rit;
            csv.set(pr.first, pr.second);
            unsigned v = csv[pr.first];
            assert(v == pr.second);

            if (i % 4096 == 0)
            {
                cout << "\r" << pr.first << "/" << max_size << flush;
                csv.optimize();
            }

            ++i;
        } // for
        
        cout << "\nRSC verification..." << endl;
        
        csv.optimize();
        csv.sync();
        i = 0;
        for (i = 0; i < vect.size(); ++i)
        {
            const std::pair<unsigned, unsigned>& pr = vect[i];
            unsigned v = csv[pr.first];
            assert(v == pr.second);
            if (i % 4096 == 0)
                cout << "\r" << pr.first << "/" << max_size << flush;
        } // for
        
        cout << "\nRSC set null..." << endl;

        i = 0;
        for (auto rit = vect.rbegin(); rit != vect.rend(); ++rit)
        {
            std::pair<unsigned, unsigned> pr = *rit;
            csv.set_null(pr.first);
            assert(csv.is_null(pr.first));
            if (i % 4096 == 0)
            {
                cout << "\r" << i << "/" << max_size << flush;
                csv.optimize();
            }
            ++i;
        } // for


        
        cout << "\nOK" << endl;
    }
    
    
    {
    cout << "decode() test" << endl;
    
        {
        rsc_sparse_vector_u32 csv1;
        
        csv1.push_back(5, 1);
        csv1.push_back(6, 1);
        csv1.push_back(8, 2);
        csv1.push_back(255, 4);

        csv1.sync();

        for (unsigned k = 0; k < 2; ++k)
        {
            CheckCompressedDecode(csv1, 0, 1);
            CheckCompressedDecode(csv1, 0, 2);
            CheckCompressedDecode(csv1, 1, 1);

            CheckCompressedDecode(csv1, 0, 5);
            CheckCompressedDecode(csv1, 0, 6);

            CheckCompressedDecode(csv1, 256, 1);

            for (unsigned i = 0; i < csv1.size(); ++i)
            {
                CheckCompressedDecode(csv1, i, 1);
                CheckCompressedDecode(csv1, i, csv1.size());
            }
            unsigned j = csv1.size();
            for (unsigned i = 0; i < csv1.size(); ++i, --j)
            {
                unsigned size = j - i;
                if (!size)
                    break;
                CheckCompressedDecode(csv1, i, size);
            }

            csv1.optimize();
        }
        }

        {
        rsc_sparse_vector_u32 csv1;
        
        csv1.push_back(bm::id_max-1, 10);
        csv1.sync();

        for (unsigned k = 0; k < 2; ++k)
        {
            for (unsigned i = bm::id_max-20; i < csv1.size(); ++i)
            {
                CheckCompressedDecode(csv1, i, 1);
                CheckCompressedDecode(csv1, i, csv1.size()-i+10);
            }
            csv1.optimize();
        }

        }
    }

    {
    cout << "load() test" << endl;
    unsigned v;
        sparse_vector_u32 sv1(bm::use_null);
        rsc_sparse_vector_u32 csv1;
        rsc_sparse_vector_u32 csv2;

        sv1.set(10, 9);
        sv1.set(20, 200);
        sv1.set(21, 201);
        sv1.set(100, 65535);
        sv1.clear(100, true);
        
        csv1.load_from(sv1);
        assert(csv1.size() == 22);
        csv1.sync();
        
        csv2.push_back(10, 9);
        csv2.push_back(20, 200);
        csv2.push_back(21, 201);
        csv2.sync();


        v = csv1.at(10);
        assert(v == 9);
        v = csv1.at(20);
        assert(v == 200);
        v = csv1.at(21);
        assert(v == 201);
        
        bool same = csv1.equal(csv2);
        assert(same);
        
        DetailedCompareSparseVectors(csv1, sv1);
        
        rsc_sparse_vector_u32 csv4;
        csv4 = std::move(csv1);
        v = csv4.at(10);
        assert(v == 9);
        DetailedCompareSparseVectors(csv4, sv1);
        
        rsc_sparse_vector_u32 csv5(std::move(csv4));
        v = csv5.at(10);
        assert(v == 9);
        DetailedCompareSparseVectors(csv5, sv1);

    }

    {
    cout << "------ Compressed load() stress test" << endl;
    BM_DECLARE_TEMP_BLOCK(tb)
    for (unsigned i = 0; i < 10; ++i)
    {
        cout << "\nPass " << i << endl;
        
        sparse_vector_u32 sv(bm::use_null);
        rsc_sparse_vector_u32 csv1;
        
        GenerateSV(sv, i);
        
        
        csv1.load_from(sv);
        csv1.sync();
        
        cout << "cmp 1...";
        DetailedCompareSparseVectors(csv1, sv);
        DetailedCheckCompressedDecode(csv1);
        cout << "ok" << endl;
        
        cout << "cmp 2...";
        csv1.optimize(tb);
        DetailedCompareSparseVectors(csv1, sv);
        DetailedCheckCompressedDecode(csv1);
        cout << "ok" << endl;

        cout << "cmp 3...";
        csv1.clear();

        sv.optimize(tb);
        rsc_sparse_vector_u32 csv2;
        csv2.load_from(sv);
        DetailedCompareSparseVectors(csv2, sv);
        csv1.sync();
        DetailedCheckCompressedDecode(csv1);

        csv2.optimize(tb);
        csv2.sync();

        DetailedCompareSparseVectors(csv2, sv);
        DetailedCheckCompressedDecode(csv1);
        cout << "ok" << endl;

        cout << "cmp 4...";
        {
        rsc_sparse_vector_u32 csv3(csv2);
        DetailedCompareSparseVectors(csv3, sv);
        }
        cout << "ok" << endl;

        cout << "cmp 5...";
        {
        rsc_sparse_vector_u32 csv4;
        csv4 = std::move(csv2);
        DetailedCompareSparseVectors(csv4, sv);

        rsc_sparse_vector_u32 csv5(std::move(csv4));
        DetailedCompareSparseVectors(csv5, sv);
        }
        cout << "ok" << endl;
    } // for
    cout << "Compressed load() stress test OK" << endl;

    
    }
    
    
    cout << " ------------------------------ Test Compressed Sparse Vector OK" << endl;
}

static
void TestCompressSparseVectorSerial()
{
    cout << " ------------------------------ TestCompressSparseVectorSerial()" << endl;

    {
        rsc_sparse_vector_u32 csv1;
        rsc_sparse_vector_u32 csv2;
        {
        rsc_sparse_vector_u32::back_insert_iterator rs_bi = csv1.get_back_inserter();
            rs_bi.add_null();
            rs_bi.add(1);
            rs_bi.add(2);
            rs_bi.add_null();
            rs_bi.add(4);

            rs_bi.flush();
        }

        BM_DECLARE_TEMP_BLOCK(tb)
        sparse_vector_serial_layout<rsc_sparse_vector_u32> sv_lay;
        bm::sparse_vector_serialize<rsc_sparse_vector_u32>(csv1, sv_lay, tb);
        const unsigned char* buf = sv_lay.buf();

        sparse_vector_u32::bvector_type bv_mask;
        bv_mask.set(1);
        bv_mask.set(2);
        bv_mask.set(100);
        bm::sparse_vector_deserializer<rsc_sparse_vector_u32> sv_deserial;
        sv_deserial.deserialize(csv2, buf, bv_mask);

        assert(csv2.size() == csv1.size());
        assert(csv2.get(1) == 1);
        assert(csv2.get(2) == 2);

        sv_deserial.deserialize(csv2, buf, 1, 2);
        assert(csv2.size() == csv1.size());
        assert(csv2.get(1) == 1);
        assert(csv2.get(2) == 2);

    }

    cout << "\nRSC gather stress test ..." << endl;
    {
        rsc_sparse_vector_u32 csv1;
        rsc_sparse_vector_u32 csv2;
        rsc_sparse_vector_u32 csv3;

        rsc_sparse_vector_u32::size_type from = bm::id_max32 / 2;
        rsc_sparse_vector_u32::size_type to = from + 257 * 65536;

        cout << "   generation... " << endl;
        {
            rsc_sparse_vector_u32::back_insert_iterator rs_bi = csv1.get_back_inserter();
            rs_bi.add_null();
            rs_bi.add(1);
            rs_bi.add(2);
            rs_bi.add_null();
            rs_bi.add(4);
            rs_bi.add_null(from); // add many NULLs

            for (auto i = from; i < to; ++i)
            {
                rs_bi.add(i);
                rs_bi.add_null();
            }
            rs_bi.flush();
        }
        cout << "   generation OK" << endl;
        {
            rsc_sparse_vector_u32 csv_range;
            csv_range.set(1, 10);
            csv_range.copy_range(csv1, from - 10, from-1);
            assert(csv_range.get(1) == 0);

        }
        sparse_vector_serial_layout<rsc_sparse_vector_u32> sv_lay;

        bm::sparse_vector_serializer<rsc_sparse_vector_u32> sv_serializer;
        sv_serializer.set_bookmarks(false);

        for (unsigned pass = 0; pass < 2; ++pass)
        {
            cout << "\nPASS=" << pass << endl;
            sv_serializer.serialize(csv1, sv_lay);

            const unsigned char* buf = sv_lay.buf();

            size_t cnt = 0;
            auto j = to;
            for (auto i = from; i < j; ++cnt, ++i, --j)
            {
                rsc_sparse_vector_u32::bvector_type bv_mask;
                bv_mask.set_range(i, j);
                bm::sparse_vector_deserializer<rsc_sparse_vector_u32> sv_deserial;
                sv_deserial.deserialize(csv2, buf, bv_mask);
                assert(csv2.get(j + 1) == 0);
                csv2.sync();


                {
                    rsc_sparse_vector_u32 csv_range;
                    csv_range.copy_range(csv1, i, j);

                    rsc_sparse_vector_u32::size_type pos;
                    bool f = bm::sparse_vector_find_first_mismatch(csv_range, csv2, pos, bm::no_null);
                    if (f)
                    {
                        auto v2 = csv2.get(pos);
                        auto rv1 = csv_range.get(pos);
                        std::cerr << "Discrepancy at idx=" << pos << endl;
                        std::cerr << v2 << "!=" << rv1 << endl;
                        csv_range.copy_range(csv1, i, j);
                        assert(0);
                    }

                    if (!cnt || (j - i) < 65536)
                    {
                        for (auto i0 = i; i0 < j; ++i0)
                        {
                            auto v1 = csv1[i0];
                            auto v2 = csv2[i0];
                            assert(v1 == v2);
                            assert(csv1.is_null(i0) == csv2.is_null(i0));
                            auto rv1 = csv_range[i0];
                            assert(rv1 == v1);
                            assert(csv_range.is_null(i0) == csv2.is_null(i0));
                        } // for i0
                    }

                }
                sv_deserial.deserialize(csv3, buf, i, j);
                bool eq = csv2.equal(csv3);
                assert(eq);

                cout << "\r  " << (j-i) << "           " << std::flush;

                csv1.sync();

                if (cnt < 512 || (j - i) < 65536/2)
                {
                    continue;
                }

                // gallop to the end
                i += rand() % 65536;
                j -= rand() % 65536;
            } // for i

            cout << "\n bookmarks ON" << endl;
            sv_serializer.set_bookmarks(true);
        } // for pass


    }
    cout << "\nOK" << endl;


    cout << " ------------------------------ TestCompressSparseVectorSerial() OK" << endl;
}

static
void TestHeapVector()
{
    {
        bm::heap_vector<bm::id64_t, bvect::allocator_type> hv;
        hv.add() = ~10ull;
        hv.push_back(~0ull);
        assert(hv[0] == ~10ull);
        assert(hv[1] == ~0ull);
    }
    {
        bm::heap_vector<bm::id64_t, bvect::allocator_type> hv;
        for (unsigned i = 0; i < 65535; ++i)
            hv.push_back(i);
        for (unsigned i = 0; i < 65535; ++i)
        {
            assert(hv[i] == i);
        }
    }

    {
        bm::heap_vector<bvect, bvect::allocator_type> hv;

        bvect& bv0 = hv.add();
        bv0.set(1);
        bvect& bv1 = hv.add();
        bv1.set(2);
        assert(hv.size() == 2);
        assert(hv[1].test(2));
        
        hv.resize(1);
        assert(hv.size() == 1);
        assert(hv[0].test(1));

        hv.resize(2);
        assert(hv.size() == 2);
        assert(!hv[1].any());

        bm::heap_vector<bvect, bvect::allocator_type> hv2(hv);
        assert(hv2.size() == 2);
        assert(hv2[0].test(1));
        assert(!hv2[1].any());

        bm::heap_vector<bvect, bvect::allocator_type> hv3;
        hv3.reserve(10);
        hv3 = hv;
        hv3[1].set(0);
        assert(hv3.size() == 2);
        assert(hv3.at(0).test(1));
        assert(hv3.at(1).any());

        bm::heap_vector<bvect, bvect::allocator_type> hv4;
        hv4.swap(hv3);
        assert(hv3.size() == 0);

    }

}

static
void TestXOR_RefVector()
{
    cout << " ------------------------------ TestXOR_RefVector()" << endl;

    {
        bv_ref_vector<bvect> ref_vect;
        assert(ref_vect.size() == 0);

        bvect bv1, bv2;
        ref_vect.add(&bv1, 10);
        ref_vect.add(&bv2, 15);

        assert(ref_vect.size() == 2);
        assert(ref_vect.get_bv(0) == &bv1);
        assert(ref_vect.get_bv(1) == &bv2);
        assert(ref_vect.get_row_idx(0) == 10);
        assert(ref_vect.get_row_idx(1) == 15);

        size_t idx = ref_vect.find(15);
        assert(idx == 1);
        idx = ref_vect.find(10);
        assert(idx == 0);

        idx = ref_vect.find(100);
        assert(idx == ref_vect.not_found());

        ref_vect.reset();
        assert(ref_vect.size() == 0);
    }

    cout << " ------------------------------ TestXOR_RefVector() OK" << endl;
}

static
void show_help()
{
    std::cerr
        << "BitMagic C++ stress test." << endl
        << "-h                - help" << endl
        << "-llevel (or -ll)      - low level tests" << endl
        << "-support (or -s)      - support containers " << endl
        << "-bvbasic (or -bvb     - bit-vector basic " << endl
        << "-bvser                - bit-vector serialization " << endl
        << "-bvops (-bvo, -bvl)   - bit-vector logical operations" << endl
        << "-bvshift (or -bvs)    - bit-vector shifts " << endl
        << "-rankc (or -rc)       - rank-compress " << endl
        << "-agg (or -aggregator) - bm::aggregator " << endl
        << "-sv                   - test sparse vectors" << endl
        << "-csv                  - test compressed sparse vectors"
        << "-strsv                - test sparse vectors" << endl
        << "-cc                   - test compresses collections" << endl
      ;
}

bool         is_all = true;
bool         is_low_level = false;
bool         is_support = false;
bool         is_bvbasic = false;
bool         is_bvser = false;
bool         is_bvops = false;
bool         is_bvshift = false;
bool         is_rankc = false;
bool         is_agg = false;
bool         is_sv = false;
bool         is_csv = false;
bool         is_str_sv = false;
bool         is_c_coll = false;

static
int parse_args(int argc, char *argv[])
{
    for (int i = 1; i < argc; ++i)
    {
        std::string arg = argv[i];
        if ((arg == "-h"))
        {
            show_help();
            return 1;
        }
        if (arg == "-ll" || arg == "-llevel")
        {
            is_all = false;
            is_low_level = true;
            continue;
        }
        if (arg == "-s" || arg == "-support")
        {
            is_all = false;
            is_support = true;
            continue;
        }
        if (arg == "-bvb" || arg == "-bvbasic")
        {
            is_all = false;
            is_bvbasic = true;
            continue;
        }
        if (arg == "-bvser" || arg == "-ser")
        {
            is_all = false;
            is_bvser = true;
            continue;
        }
        if (arg == "-bvo" || arg == "-bvops" || arg == "-bvl")
        {
            is_all = false;
            is_bvops = true;
            continue;
        }
        if (arg == "-bvs" || arg == "-bvshift")
        {
            is_all = false;
            is_bvshift = true;
            continue;
        }
        if (arg == "-rc" || arg == "-rankc")
        {
            is_all = false;
            is_rankc = true;
            continue;
        }
        if (arg == "-agg" || arg == "-aggregator")
        {
            is_all = false;
            is_agg = true;
            continue;
        }
        if (arg == "-sv")
        {
            is_all = false;
            is_sv = true;
            continue;
        }
        if (arg == "-csv")
        {
            is_all = false;
            is_csv = true;
            continue;
        }

        if (arg == "-strsv" || arg == "-svstr")
        {
            is_all = false;
            is_str_sv = true;
            continue;
        }
        if (arg == "-cc")
        {
            is_all = false;
            is_c_coll = true;
            continue;
        }

    } // for i
    return 0;
}

int main(int argc, char *argv[])
{
    time_t      start_time = time(0);
    time_t      finish_time;

    {
    auto ret = parse_args(argc, argv);
    if (ret != 0)
        return ret;
    }
/*
    BrennerTest("/Users/anatoliykuznetsov/Desktop/dev/git/BitMagic/tests/stress/1.i",
                "/Users/anatoliykuznetsov/Desktop/dev/git/BitMagic/tests/stress/1.bv",
                "/Users/anatoliykuznetsov/Desktop/dev/git/BitMagic/tests/stress/new.bv");
    return 0;
*/

/*
    LoadBVDump("C:/dev/group-by-sets/sets/geo_organization.bvdump", 
               "C:/dev/group-by-sets/sets/geo_organization.bvdump2", 
               false); // not validate!
    exit(1);
*/

/*
    LoadBVDump("C:/dev/group-by-sets/sets/geo_organization.dat", 
               "C:/dev/group-by-sets/sets/geo_organization.bvdump", 
               true); //  validate!

    LoadBVDump("C:/dev/group-by-sets/sets/exec_time_msec_bit.dat", 
               "C:/dev/group-by-sets/sets/exec_time_msec_bit.bvdump", 
               true);

    LoadBVDump("C:/dev/group-by-sets/sets/geo_country.dat", 
               "C:/dev/group-by-sets/sets/geo_country.bvdump", 
               true);
    LoadBVDump("C:/dev/group-by-sets/sets/log_displayeduids.dat", 
               "C:/dev/group-by-sets/sets/log_displayeduids.bvdump", 
               true);

    //LoadVectors("c:/dev/bv_perf", 3, 27);
    exit(1);
*/

//avx2_i32_shift();
//return 0;

//unsigned long long a = 9223372036854775807ULL;
//unsigned long long a = 281474976710655ULL;
//a = a / (65536 * 256);

    if (is_all || is_low_level)
    {
        TestRecomb();

        OptimGAPTest();

        CalcBeginMask();
        CalcEndMask();

        TestSIMDUtils();

        TestArraysAndBuffers();

        TestFindBlockDiff();

        Log2Test();
        FindNotNullPtrTest();

        LZCNTTest();

        SelectTest();

         TestBlockZero();

         TestBlockDigest();

         TestBlockAND();
         TestBlockSUB();
         TestBlockOR();

         TestBlockCountChange();

         TestBlockCountXORChange();

         TestBlockToGAP();

         ShiftRotateTest();

         BlockBitInsertTest();
        
         BlockBitEraseTest();
         TestBlockLast();

         BitForEachTest();
        
         BitCountChangeTest();
         WordCmpTest();
        
        //BitBlockTransposeTest();
    }
    
    if (is_all || is_support)
    {
        TestHeapVector();
        TestXOR_RefVector();

        MiniSetTest();
        BitEncoderTest();
      
        InterpolativeCodingTest();
        GammaEncoderTest();
        GAPCheck();
        SerializationBufferTest();
        TestBasicMatrix();

        DynamicMatrixTest();
        RSIndexTest();
    }

    if (is_all || is_bvbasic)
    {

         ExportTest();
         ResizeTest();

         SyntaxTest();

         SetTest();

         EmptyBVTest();
         ClearAllTest();

         EnumeratorTest();

         CountRangeTest();

         KeepRangeTest();

         BasicFunctionalityTest();

         OptimizeTest();
  
         RankFindTest();

         BvectorBitForEachTest();

         GetNextTest();

         BvectorIncTest();

         BvectorBulkSetTest();

        GAPTestStress();
        
        MaxSTest();

        SimpleRandomFillTest();

        RangeRandomFillTest();

        RangeCopyTest();

        ComparisonTest();

        BvectorFindFirstDiffTest();
        
        MutationTest();
        MutationOperationsTest();
        
        BlockLevelTest();
     }
    
    if (is_all || is_bvser || is_bvbasic)
    {
        SerializationCompressionLevelsTest();
        SerializationTest();
        DesrializationTest2();
    }
    
    if (is_all || is_bvshift)
    {
         BvectorShiftTest();
         BvectorInsertTest();
         BvectorEraseTest();
    }


    if (is_all || is_rankc)
    {
         AddressResolverTest();
         TestRankCompress();
    }

    if (is_all || is_bvops)
    {

        AndOperationsTest(true); // enable detailed check
        OrOperationsTest(true);
        XorOperationsTest(true);
        SubOperationsTest(true);

        StressTest(150, 0, false); // OR - detailed check disabled
        StressTest(150, 3, false); // AND
        StressTest(150, 1, false); // SUB
        StressTest(150, 2, false); // XOR
    }



    if (is_all || is_agg)
    {
         AggregatorTest();
         StressTestAggregatorOR(100);
         StressTestAggregatorAND(100);
         StressTestAggregatorShiftAND(5);

    //     StressTestAggregatorSUB(100);
    }

    if (is_all || is_sv)
    {
        TestSparseVector();

        TestSparseVectorAlgo();

        TestSparseVector_XOR_Scanner();

        TestSparseVectorInserter();

        TestSparseVectorGatherDecode();

        TestSparseVectorSerial();

        TestSparseVectorSerialization2();

        TestSparseVectorTransform();

        TestSparseVectorRange();

        TestSparseVectorFilter();

        TestSparseVectorScan();

        TestSparseSort();
    }

    if (is_all || is_csv)
    {
        TestCompressSparseVector();

        TestCompressedSparseVectorAlgo();
        
        TestCompressSparseVectorSerial();

        TestCompressedSparseVectorScan();

        TestSparseVector_Stress(3);
    }

    if (is_all || is_c_coll)
    {
         TestCompressedCollection();
    }
    
    if (is_all || is_str_sv)
    {
         TestStrSparseVector();

         TestStrSparseVectorAlgo();

         TestStrSparseVectorSerial();

         TestStrSparseSort();

         StressTestStrSparseVector();
    }

    if (is_all || is_bvops)
    {
        StressTest(300, -1, true);
    }

    finish_time = time(0);


    cout << "Test execution time = " << finish_time - start_time << endl;

#ifdef MEM_DEBUG
    cout << "[--------------  Allocation digest -------------------]" << endl;
    cout << "Number of BLOCK allocations = " <<  dbg_block_allocator::na_ << endl;
    cout << "Number of PTR allocations = " <<  dbg_ptr_allocator::na_ << endl << endl;

    if(dbg_block_allocator::balance() != 0)
    {
        cout << "ERROR! Block memory leak! " << endl;
        cout << dbg_block_allocator::balance() << endl;
        exit(1);
    }

    if(dbg_ptr_allocator::balance() != 0)
    {
        cout << "ERROR! Ptr memory leak! " << endl;
        cout << dbg_ptr_allocator::balance() << endl;
        exit(1);
    }
    cout << "[------------  Debug Allocation balance OK ----------]" << endl;
#endif

    return 0;
}


