/*
Copyright(c) 2019 Anatoliy Kuznetsov(anatoliy_kuznetsov at yahoo.com)

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
#include <iterator>
#include <stdarg.h>
#include <vector>
#include <chrono>

#include <bm64.h>
#include <bmrandom.h>
#include <bmaggregator.h>
#include <bmvmin.h>
#include <bmdbg.h>
#include <bmalgo.h>
#include <bmintervals.h>
#include <bmbvimport.h>
#include <bmsparsevec_util.h>
#include <bmtimer.h>

#include <bmsparsevec.h>
#include <bmsparsevec_algo.h>
#include <bmsparsevec_serial.h>
#include <bmsparsevec_compr.h>
#include <bmstrsparsevec.h>


using namespace bm;
using namespace std;


#include "gena.h"
#include "test_util.h"


#define POOL_SIZE 5000

//#define MEM_POOL

bool is_silent = false;

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

#if defined(BMSSE2OPT) || defined(BMSSE42OPT) || defined(BMAVX2OPT) || defined(BMAVX512OPT) || defined(__ARM_NEON__)
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
            (bm::word_t*) ::malloc((n+2) * sizeof(bm::word_t));
        if (!p)
        {
            std::cerr << "ERROR Failed allocation!" << endl;
            assert(0); exit(1);
        }
        size_t* ptr = (size_t*)p;
        *ptr = n;
        p+=2;
        return p;
    }

    static void deallocate(bm::word_t* p, size_t n)
    {
        ++nf_;
        p-=2;
        size_t* ptr = (size_t*)p;
        if (*ptr != n)
        {
            cerr << "Block memory deallocation ERROR! n = " << n << ", block-control = " << *ptr << endl;
            assert(0); exit(1);
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

typedef bm::bvector<dbg_alloc> bvect64;
typedef bm::bvector<dbg_alloc> bvect;
typedef bm::bvector_mini<dbg_block_allocator> bvect_mini;
typedef bm::rs_index<dbg_alloc> rs_ind;

#else

#ifdef MEM_POOL

typedef mem_alloc<pool_block_allocator, pool_ptr_allocator> pool_alloc;
typedef bm::bvector<pool_alloc> bvect64;
typedef bm::bvector<pool_alloc> bvect;
typedef bm::bvector_mini<bm::block_allocator> bvect_mini;
typedef bm::rs_index<pool_block_allocator> rs_ind;


#else

typedef bm::bvector<> bvect64;
typedef bm::bvector<> bvect;
typedef bm::bvector_mini<bm::block_allocator> bvect_mini;
typedef bm::rs_index<> rs_ind;

#endif

#endif

typedef std::vector<bm::id64_t> ref_vect;

typedef bm::sparse_vector<unsigned, bvect > sparse_vector_u32;
typedef bm::sparse_vector<int, bvect > sparse_vector_i32;
typedef bm::sparse_vector<unsigned long long, bvect > sparse_vector_u64;
typedef bm::sparse_vector<signed long long, bvect > sparse_vector_i64;
typedef bm::rsc_sparse_vector<unsigned, sparse_vector_u32> rsc_sparse_vector_u32;
typedef bm::rsc_sparse_vector<int, sparse_vector_i32> rsc_sparse_vector_i32;


static
void SyntaxTest()
{
    cout << "------------------------------------ SyntaxTest()" << endl;
    ref_vect vect;
    generate_vect_simpl0(vect);
    {
        bvect64 bv0;

        load_BV_set_ref(cout, bv0, vect);
        compare_BV_set_ref(bv0, vect);

        
        auto idx = vect.size()-1;
        auto v = vect[idx];
        assert(bv0.test(v));
        
        bvect64 bv1(bv0);
        cout << bv0.count() << endl;
        cout << bv1.count() << endl;        
        compare_BV_set_ref(bv1, vect);
        
        bvect64 bv2;
        bv2 = bv0;
        compare_BV_set_ref(bv1, vect);
        
        bvect64 bv3(std::move(bv2));
        assert(bv2.none());
        compare_BV_set_ref(bv3, vect);
        
        bvect64 bv4;
        bv4.move_from(bv3);
        assert(bv3.none());
        compare_BV_set_ref(bv4, vect);
        
        bv0 &= bv4;
        compare_BV_set_ref(bv0, vect);
        bv0 |= bv4;
        compare_BV_set_ref(bv0, vect);
        bv0 -= bv4;
        assert(bv0.none());
        bv1 ^= bv4;
        assert(bv1.none());
    }
    
    {
        bvect64 bv0, bv1, bv2;

        load_BV_set_ref(cout, bv0, vect);
        bv1 = bv2 = bv0;
        bool b = (bv1 == bv0);
        assert(b);
        b = (bv1 == bv2);
        assert(b);
     
        {
            bvect64 bv3;
            bv3 = bv1 | bv2;
            b = (bv1 == bv3);
            assert(b);
            compare_BV_set_ref(bv3, vect);
        }
        {
            bvect64 bv3;
            bv3 = bv1 & bv2;
            b = (bv1 == bv3);
            assert(b);
            compare_BV_set_ref(bv3, vect);
        }
        {
            bvect64 bv3;
            bv3 = bv1 - bv2;
            assert(bv3.none());
        }
        {
            bvect64 bv3;
            bv3 = bv1 ^ bv2;
            assert(bv3.count() == 0ULL);
        }
    }
    
    {
        bvect64 bv1;
        bvect64::reference ref = bv1[10];
        bool bn = !ref;
        assert(bn);
        bool bn2 = ~ref;
        assert(bn2);
        bv1[10] = bn2;
        bv1[10] = bn;
        bn = bn2 = false;
        assert(!bn);
        ref.flip();
        assert(!bv1.test(10));
        bv1[bm::id_max-1] = 1;
        assert(bv1[bm::id_max-1]);
        auto ref1 = bv1[bm::id_max-1];
        ref1.flip();
        assert(!ref1);
        assert(!bv1.test(bm::id_max-1));
    }
    {
        bvect64 bv1;
        auto ii = bv1.inserter();
        ii = bm::id_max / 2;
        assert(bv1.test(bm::id_max / 2));
    }
    
    cout << "------------------------------------ SyntaxTest() OK" << endl;
}

static
void GenericBVectorTest()
{
    cout << "------------------------------------ GenericBVectorTest()" << endl;
    
    {
        bvect64 bv0;
        ref_vect vect;
        generate_vect_simpl0(vect);
        
        load_BV_set_ref(cout, bv0, vect);
        compare_BV_set_ref(bv0, vect);

        bvect64 bv1(bm::BM_GAP);
        load_BV_set_ref(cout, bv1, vect);
        compare_BV_set_ref(bv1, vect);
        
        int cmp = bv0.compare(bv1);
        assert(cmp == 0);
        
        bvect64::statistics st1, st2, st3;
        bv0.calc_stat(&st1);
        cout << "st1.max_serialize_mem = " << st1.max_serialize_mem << endl;
        assert(st1.max_serialize_mem < 2*(sizeof(bm::word_t)*st1.bit_blocks * 2048));
        assert(st1.ptr_sub_blocks == 6);
        bv0.optimize(0, bvect64::opt_compress, &st2);
        assert(st1.ptr_sub_blocks >= st2.ptr_sub_blocks);
        assert(st2.max_serialize_mem <= st1.max_serialize_mem);
        bv1.calc_stat(&st3);
        assert(st1.ptr_sub_blocks >= st3.ptr_sub_blocks);
        assert(st2.gap_blocks == st3.gap_blocks);
        assert(st3.max_serialize_mem <= st1.max_serialize_mem);
    }

    cout << "------------------------------------ GenericBVectorTest() OK" << endl;
}

static
void ArenaTest()
{
   cout << "----------------------------- ArenaTest() " << endl;

   {
        bm::bv_arena_statistics st;
        bvect bv;

        bv.get_blocks_manager().calc_arena_stat(&st);
        assert(st.gap_blocks_sz == 0);
        assert(st.ptr_sub_blocks_sz == 0);
        assert(st.bit_blocks_sz == 0);

        bv.set(0);

        bv.get_blocks_manager().calc_arena_stat(&st);
        assert(st.gap_blocks_sz == 0);
        assert(st.ptr_sub_blocks_sz == bm::set_sub_array_size);
        assert(st.bit_blocks_sz == bm::set_block_size);

        bv.set(bm::id_max/2);

        bv.get_blocks_manager().calc_arena_stat(&st);
        assert(st.gap_blocks_sz == 0);
        assert(st.ptr_sub_blocks_sz == 2* bm::set_sub_array_size);
        assert(st.bit_blocks_sz == 2 * bm::set_block_size);
   }

   {
        bm::bv_arena_statistics st;
        bvect bv(bm::BM_GAP);

        bv.set(1);
        bv.set(2);

        bv.get_blocks_manager().calc_arena_stat(&st);
        assert(st.gap_blocks_sz == 4);
        assert(st.ptr_sub_blocks_sz == bm::set_sub_array_size);
        assert(st.bit_blocks_sz == 0);

        bv.set(1+bm::id_max/2);

        bv.get_blocks_manager().calc_arena_stat(&st);
        assert(st.gap_blocks_sz == 7);
        assert(st.ptr_sub_blocks_sz == 2*bm::set_sub_array_size);
        assert(st.bit_blocks_sz == 0);
   }

    {
        bm::bv_arena_statistics st;
        bvect bv;

        bv.set(0);
        bv.set(1);
        bv.optimize();

        bv.set(bm::id_max/2);

        bv.get_blocks_manager().calc_arena_stat(&st);
        assert(st.gap_blocks_sz == 3);
        assert(st.ptr_sub_blocks_sz == 2* bm::set_sub_array_size);
        assert(st.bit_blocks_sz == 1 * bm::set_block_size);


        bvect::blocks_manager_type& bman = bv.get_blocks_manager();
        bvect::blocks_manager_type::arena ar;

        bman.alloc_arena(&ar, st, bman.get_allocator());
        bman.free_arena(&ar, bman.get_allocator());

    }


   cout << "----------------------------- ArenaTest() OK" << endl;
}


static
void FreezeTest()
{
    cout << "----------------------------- FreezeTest()" << endl;

    bool eq;
    {
        bvect bv;
        bvect bv_ro(bv, bm::finalization::READONLY);
        assert(!bv.is_ro());
        assert(!bv_ro.is_ro());
    }

    // swap test
    {
        bvect bv {0};
        bvect bv_ro(bv, bm::finalization::READONLY);

        assert(!bv.is_ro());
        assert(bv_ro.is_ro());

        eq = bv.equal(bv_ro);
        assert(eq);

        bv.swap(bv_ro);

        assert(bv.is_ro());
        assert(!bv_ro.is_ro());

        eq = bv.equal(bv_ro);
        assert(eq);
    }

    {
        bvect bv;
        bv.invert();

        {
            bvect bv_ro(bv, bm::finalization::READONLY);

            assert(!bv.is_ro());
            assert(bv_ro.is_ro());

            eq = bv.equal(bv_ro);
            assert(eq);
        }
        bv.optimize();
        {
            bvect bv_ro(bv, bm::finalization::READONLY);

            assert(!bv.is_ro());
            assert(bv_ro.is_ro());

            eq = bv.equal(bv_ro);
            assert(eq);
        }
    }

    {
        bvect bv;
        bv.set_range(256*65536, bm::id_max/2 + 10);

        {
            bvect bv_ro(bv, bm::finalization::READONLY);

            assert(!bv.is_ro());
            assert(bv_ro.is_ro());

            eq = bv.equal(bv_ro);
            assert(eq); // copy-ctor
            {
                bvect bv_ro2(bv_ro, bm::finalization::READONLY);
                assert(bv_ro2.is_ro());
                eq = bv.equal(bv_ro2);
                assert(eq);
            }

            { // assignment
                bvect bv_ro2 { 126 };
                bv_ro2 = bv_ro;
                assert(bv_ro2.is_ro());
                eq = bv.equal(bv_ro2);
                assert(eq);
            }

            { // move ctor/assignment
            bvect bv_ro2 = std::move(bv_ro);
            assert(bv_ro2.is_ro());
            eq = bv.equal(bv_ro2);
            assert(eq);
            bvect bv_ro3;
            bv_ro3 = std::move(bv_ro2);
            assert(bv_ro3.is_ro());
            eq = bv.equal(bv_ro3);
            assert(eq);

            }
        }
    }



    {
        std::vector<bvect::size_type> v1 { 0, 65536, 65536 * 256, bm::id_max/2, bm::id_max-1 };
        for (size_t i = 0; i < v1.size(); ++i)
        {
            auto idx = v1[i];
            bvect bv; bv.set(idx);
            {
                for (int pass = 0; pass < 2; ++pass)
                {
                    bvect bv_ro(bv, bm::finalization::READONLY);
                    assert(!bv.is_ro());
                    assert(bv_ro.is_ro());
                    eq = bv.equal(bv_ro);
                    assert(eq);
                    bvect::size_type pos(bm::id_max);
                    bool found = bv_ro.find(pos);
                    assert(found);
                    assert(pos == idx);

                    { // freezing copyctor
                    bvect bv_ro2(bv_ro, bm::finalization::READONLY);
                    assert(bv_ro2.is_ro());
                    eq = bv.equal(bv_ro2);
                    assert(eq);
                    }

                    { // copyctor
                    bvect bv_ro2(bv_ro);
                    assert(bv_ro2.is_ro());
                    eq = bv.equal(bv_ro2);
                    assert(eq);
                    }
                    { // assignment tests
                    bvect bv_ro2 { 126 };
                    bv_ro2 = bv_ro;
                    assert(bv_ro2.is_ro());
                    eq = bv.equal(bv_ro2);
                    assert(eq);
                    }

                    { // move ctor/assignment
                    bvect bv_ro2 = std::move(bv_ro);
                    assert(bv_ro2.is_ro());
                    eq = bv.equal(bv_ro2);
                    assert(eq);
                    bvect bv_ro3;
                    bv_ro3 = std::move(bv_ro2);
                    assert(bv_ro3.is_ro());
                    eq = bv.equal(bv_ro3);
                    assert(eq);
                    }

                    bv.optimize();
                } // for pass
            }
        } // for i
    }

    {
        bvect bv { 0 };
        {
            for (bvect::size_type i = 0; i < 65535; i+=3)
                bv.set(i);
            bv.set(bm::id_max/2);
            bv.optimize();

            bvect bv_ro(bv, bm::finalization::READONLY);
            assert(!bv.is_ro());
            assert(bv_ro.is_ro());
            eq = bv.equal(bv_ro);
            assert(eq);

            bvect bv_ro2(bv_ro);
            assert(bv_ro2.is_ro());
            eq = bv_ro.equal(bv_ro2);
            assert(eq);
            eq = bv.equal(bv_ro2);
            assert(eq);
        }
    }

    cout << "----------------------------- FreezeTest() ON\n" << endl;
}


static
void SetTest()
{
    cout << "------------------------------------ SetTest()" << endl;
    {
        bvect bv;
        bv.set();
        auto cnt = bv.count();
        assert (cnt == bm::id_max);
    }
    
    {
        bvect bv;
        bv.set();
        bv.set_range(10, bm::id_max, false);
        auto cnt = bv.count();
        assert(cnt == 10);
    }
    {
        bvect bv;
        bv.set_range(bm::id_max-65535, bm::id_max, true);
        auto cnt = bv.count();
        assert(cnt == 65536);
    }

    {
        bvect bv{ 0, 10, 65536, 10000, bm::id_max-1 };
        auto cnt = bv.count();
        assert (cnt == 5);

        bvect bv2;
        bv2.set(0).set(10).set(65536).set(10000).set(bm::id_max-1);

        if (bv != bv2)
        {
            cout << "Brace initialization comparison test failed!." << endl;
            assert(0);exit(1);
        }
    }
    {
        bvect bv;
        bv.set();

        auto cnt = bv.count();
        assert (cnt == bm::id_max);

        bv.invert();
        cnt = bv.count();
        assert (cnt == 0);

        bv.set(0);
        bv.set(bm::id_max - 1);
        cnt = bv.count();
        assert(cnt == 2);

        bv.invert();
        //print_stat(bv);
        cnt = bv.count();
        assert (cnt == bm::id_max - 2);

        bv.clear();
        bv[1] &= true;
        bool v = bv[1];
        assert (!v);

        bv[1] = true;
        bv[1] &= true;
        v = bv[1];
        assert(v);

        bv.clear(true);
        bv.invert();
        bv[1] &= true;
        v = bv[1];
        assert (v);
    }
    {
        bvect bv_full;
        bv_full.invert();
        assert(bv_full.test(bm::id_max/2));
    }
    {
        bvect bv1, bv2(BM_GAP);
        bvect::size_type cnt;
        bv1.set(0); bv2.set(0);
        bv1.set(bm::id_max-1);bv2.set(bm::id_max-1);
        bv1.set((bm::id_max-1)/2);bv2.set((bm::id_max-1)/2);
        for (unsigned i = 0; i < 2; ++i)
        {
            bv1.set();
            bv2.set();
            cnt = bv1.count();
            assert (cnt == bm::id_max);
            cnt = bv2.count();
            assert (cnt == bm::id_max);
        }
    }

    {
        bvect bv2;
        bv2[bm::id_max-1] = true;
        bv2[bm::id_max-1] = false;
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
            assert(0);exit(1);
        }
    }
    
    {
        bvect bv3;
        bool changed;
        changed = bv3.set_bit_conditional(bm::id_max-10, true, true);
        bool v = bv3[10];
        if (v || changed) {
            cout << "Conditional bit set failed." << endl;
            assert(0);exit(1);
        }
        changed = bv3.set_bit_conditional(bm::id_max-10, true, false);
        v = bv3[bm::id_max-10];
        if (!v || !changed) {
            cout << "Conditional bit set failed." << endl;
            assert(0);exit(1);
        }
        changed = bv3.set_bit_conditional(bm::id_max-10, false, false);
        v = bv3[bm::id_max-10];
        if (!v || changed) {
            cout << "Conditional bit set failed." << endl;
            assert(0);exit(1);
        }
        changed = bv3.set_bit_conditional(bm::id_max-10, false, true);
        v = bv3[bm::id_max-10];
        if (v || !changed) {
            cout << "Conditional bit set failed." << endl;
            assert(0);exit(1);
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
        bvect::size_type new_size(1000001);
        bvect bv(0);
        bv.resize(new_size);
        bv[10] = true;
        bv.resize(new_size);
        bv[new_size-1] = 1;

        if (bv.size() != new_size)
        {
            cout << "Resize failed" << endl;
            exit(1);
        }
        if (bv.count() != 2ull)
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
        bv[bm::id_max-1] = true;
        assert(bv.size() == bm::id_max);
        bv.set_bit(bm::id_max-1);
        assert(bv.size() == bm::id_max);
    }

    cout << "------------------------------------ SetTest() OK" << endl;
}

static
void BVImportTest()
{
    cout << "----------------------------- BVImportTest()" << endl;

    {
        unsigned int arr[1] = {0, };
        arr[0] = 0;
        bvect bv;
        bm::bit_import_u32(bv, arr, sizeof(arr)/sizeof(arr[0]), false);
        assert(bv.count() == 0);
        arr[0] = 1;
        bm::bit_import_u32(bv, arr, sizeof(arr)/sizeof(arr[0]), false);
        assert(bv.count() == 1);
        bvect::statistics st;
        bv.calc_stat(&st);
        assert(st.bit_blocks==1);
        assert(st.gap_blocks==0);


        {
        bvect::enumerator en = bv.first();
        assert(en.valid() && *en == 0);
        }
        arr[0] = 1 | (1u << 2) | (1u << 31);
        bm::bit_import_u32(bv, arr, sizeof(arr)/sizeof(arr[0]), true);
        assert(bv.count() == 3);
        {
        bvect::enumerator en = bv.first();
        assert(en.valid() && *en == 0);
        ++en;
        assert(*en == 2);
        ++en;
        assert(*en == 31);
        bv.calc_stat(&st);
        assert(st.bit_blocks==0);
        assert(st.gap_blocks==1);
        }
    }

    {
        unsigned int arr[2048] = {0, };
        arr[0] = 0;
        bvect bv;
        bm::bit_import_u32(bv, arr, sizeof(arr)/sizeof(arr[0]), true);
        assert(bv.count() == 0);
        arr[0] = 1;
        bm::bit_import_u32(bv, arr, sizeof(arr)/sizeof(arr[0]), false);
        assert(bv.count() == 1);
        assert(bv.test(0));

        arr[2047] = 1u << 31;
        bm::bit_import_u32(bv, arr, sizeof(arr)/sizeof(arr[0]), true);
        assert(bv.count() == 2);
        assert(bv.test(0));
        assert(bv.test(65535));
    }

    {
        unsigned int arr[2048 + 10] = {0, };
        arr[0] = 0;
        bvect bv;
        bm::bit_import_u32(bv, arr, sizeof(arr)/sizeof(arr[0]), false);
        assert(bv.count() == 0);
        arr[0] = 1 << 16;
        bm::bit_import_u32(bv, arr, sizeof(arr)/sizeof(arr[0]), false);
        assert(bv.count() == 1);
        assert(bv.test(16));

        arr[2047] = 1u << 31;
        arr[2048] = 1u << 7;
        bm::bit_import_u32(bv, arr, sizeof(arr)/sizeof(arr[0]), false);
        assert(bv.count() == 3);
        assert(bv.test(16));
        assert(bv.test(65535));
        assert(bv.test(65536 + 7));
    }


    cout << "... import stress test" << endl;
    {
        unsigned max = 250000;
        for (unsigned size = 1; size < max; ++size)
        {
            std::vector<unsigned> vect;
            vect.resize(size);
            bvect bv_control;
            {
                bvect::bulk_insert_iterator iit(bv_control);
                for (size_t i = 0; i < vect.size(); ++i)
                {
                    unsigned w = (i & 1) ? (unsigned)rand() : 0;
                    vect[i] = w;
                    bvect::size_type base_idx = bvect::size_type(i * 32);
                    for (unsigned k = 0; (k < 32) && w; ++k)
                    {
                        unsigned mask = (1u << k);
                        if (w & mask) {
                            iit = (base_idx + k);
                        }
                        w &= ~mask;
                    } // for k
                }
                iit.flush();
            }

            bvect bv;
            bm::bit_import_u32(bv, vect.data(), (unsigned)vect.size(), true);
            bool eq = bv.equal(bv_control);
            assert(eq);

            if (size % 256 == 0)
                cout << "\r" << size << " of " << max << "      " << flush;
        } // for size
    }


    cout << "\n----------------------------- BVImportTest() OK" << endl;
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
        bm::export_array(bv1, buf + 0, buf + 20);

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
        unsigned buf[20] = {0,};

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
void ResizeTest()
{
    cout << "---------------------------- ResizeTest()" << endl;
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
        auto cnt = bv.count();
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
        //print_stat(bv);
        bv.set(100);
        bv.set(65536 + 10);
        print_stat(cout,bv);
        bv.set_range(0, 65536*10, false);
        print_stat(cout,bv);
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
    //
    {{
        bvect bv1(10);
        bm::id64_t maxs= bm::id_max - 100;
        {
            bvect::insert_iterator it(bv1);
            *it = 100 * 65536;
            *it = maxs;
        }
        auto sz = bv1.size();
        assert(sz == maxs+1);
    }}

    // serialization
    //
    {{
        const bvect::size_type test_size = bm::id_max - 100;
        bvect bv1(test_size);
        bv1.set(test_size - 1005);
        struct bvect::statistics st1;
        bv1.calc_stat(&st1);

        unsigned char* sermem = new unsigned char[st1.max_serialize_mem];
        size_t slen2 = bm::serialize(bv1, sermem);
        cout << slen2 << endl;

        bvect bv2(0);
        bm::deserialize(bv2, sermem);
        delete [] sermem;

        assert(bv2.size() == test_size);
        assert(bv2.count() == 1);
        auto first = bv2.get_first();
        assert(first == test_size - 1005);
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
    cout << "---------------------------- ResizeTest() OK" << endl;
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
        assert(0);
        exit(1);
    }
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
        auto cnt = bv1.count_range(0, 10);
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
            assert(0);exit(1);
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
        bv2.set_bit(bm::id_max-100);
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
            assert(0);
            exit(1);
            
        }
        bv2.set_bit(bm::id_max-1000);
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
        bvect::enumerator en2 = bv1.get_enumerator(0ULL);
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
void EnumeratorTest()
{
    cout << "-------------------------------------------- EnumeratorTest" << endl;

    {
    bvect bvect1;

    bvect1.set_bit(100);
    
    {
        auto n = bvect1.get_next(101);
        assert(!n);
    }

    bvect::enumerator en = bvect1.first();
    auto n = bvect1.get_next(0);
    
    bvect::enumerator en1 = bvect1.get_enumerator(n);
    if (*en != 100 || n != 100 || *en1 != 100)
    {
        cout << "1.Enumerator error !" << endl;
        exit(1);
    }
    CompareEnumerators(en, en1);

    bvect1.clear_bit(100);

    bvect1.set_bit(bm::id_max - 100);
    en.go_first();
    n = bvect1.get_next(0);
    en1.go_to(n);
    if (*en != bm::id_max - 100 || n != *en || *en1 != *en)
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
    if (*en != bm::id_max - 100 || n != *en || *en1 != *en)
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

    cout << "FULL bvector enumerator stress test (0)..." << endl;
    {
        bvect bvect1;
        bvect1.set();
        
        {
            bvect::enumerator en2(&bvect1, bm::id_max-1);
            ++en2;
            bool b = en2.valid();
            assert(!b);
        }

        bvect::enumerator en = bvect1.first();
        auto num = bvect1.get_first();
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
    cout << "FULL bvector enumerator stress test (0) ... OK" << endl;

    
    cout << "FULL bvector enumerator stress test (1)..." << endl;
    {
        bvect bvect1;
        bvect1.set();
        
        bvect::size_type start_idx = bm::id_max - (bm::set_sub_array_size * bm::gap_max_bits * 2);

        bvect::enumerator en(&bvect1, start_idx);
        while (en.valid())
        {
            bvect::size_type pos;
            bool b = bvect1.find(start_idx, pos);
            if (*en != pos || !b)
            {
                cout << "2. Enumeration comparison failed !" <<
                        " enumerator = " << *en <<
                        " find() = " << pos << endl;
                assert(0);
                exit(1);
            }

            ++en;
            pos = bvect1.get_next(pos);
            if (pos)
            {
                bvect::enumerator en2(&bvect1, pos);
                if (*en2 != pos)
                {
                    cout << "2. Enumeration comparison failed !" <<
                            " enumerator = " << *en <<
                            " get_next() = " << pos << endl;
                    assert(0);
                    exit(1);
                }
                CompareEnumerators(en, en2);
            }
            else
            {
                assert(start_idx == bm::id_max-1);
            }
            start_idx = pos;
        } // while
    }
    cout << "FULL bvector enumerator stress test (1) ... OK" << endl;
    
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
        bvect::size_type num = bvect1.get_first();
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

        bvect1.set_bit(bm::id_max - 100);
        en.go_first();
        en1.go_to(10);

        if ((*en != bm::id_max - 100) || *en != *en1)
        {
            cout << "Enumerator error !" << endl;
            exit(1);
        }
        CompareEnumerators(en, en1);
        print_stat(cout,bvect1);
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

        auto num = bvect1.get_first();

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
void IntervalEnumeratorTest()
{
    cout << "----------------------------- IntervalEnumeratorTest()" << endl;

    bool valid;
    cout << "empty bvector tests" << endl;
    {
        bm::interval_enumerator<bvect> ien;
        valid = ien.valid();
        assert(!valid);
    }

    {
        bvect bv;
        bm::interval_enumerator<bvect> ien(bv);

        valid = ien.valid();
        assert(!valid);
    }


    cout << "inverted bvector tests" << endl;
    {
        bvect bv;
        bv.invert();
        bm::interval_enumerator<bvect> ien(bv);

        valid = ien.valid();
        assert(valid);
        assert(ien.start() == 0);
        assert(ien.end() == bm::id_max-1);
    }

    cout << "GAP bvector tests" << endl;
    {
        bvect bv;
        bv.set_range(0, 33);

        bm::interval_enumerator<bvect> ien(bv);
        valid = ien.valid();
        assert(valid);
        assert(ien.start() == 0);
        assert(ien.end() == 33);

        bv.set_range(bm::id_max/2, bm::id_max/2 + 2);
        ien.go_to(100);
        valid = ien.valid();
        assert(valid);
        assert(ien.start() == bm::id_max/2);
        assert(ien.end() == bm::id_max/2 + 2);

        bv.set_range(bm::id_max-1, bm::id_max-1);
        ien.go_to(bm::id_max-2);
        valid = ien.valid();
        assert(valid);
        assert(ien.start() == bm::id_max-1);
        assert(ien.end() == bm::id_max-1);

        ien.go_to(0);
        valid = ien.valid();
        assert(valid);
        assert(ien.start() == 0);
        assert(ien.end() == 33);

        valid = ien.advance();
        assert(valid);
        assert(ien.start() == bm::id_max/2);
        assert(ien.end() == bm::id_max/2 + 2);

        valid = ien.advance();
        assert(valid);
        assert(ien.start() == bm::id_max-1);
        assert(ien.end() == bm::id_max-1);

        valid = ien.advance();
        assert(!valid);
    }

    {
        bvect bv { 0 };
        bm::interval_enumerator<bvect> ien(bv);

        valid = ien.valid();
        assert(valid);
        assert(ien.start() == 0);
        assert(ien.end() == 0);
    }


    {
        bvect bv { bm::id_max-1};
        bm::interval_enumerator<bvect> ien(bv);

        valid = ien.valid();
        assert(valid);
        assert(ien.start() == bm::id_max-1);
        assert(ien.end() == bm::id_max-1);

    }

    {
        bvect bv { 0, 100, bm::id_max-1 };
        for (unsigned pass = 0; pass < 2; ++pass)
        {
            bm::interval_enumerator<bvect> ien(bv);

            valid = ien.valid();
            assert(valid);
            assert(ien.start() == 0);
            assert(ien.end() == 0);

            valid = ien.advance();
            assert(valid);
            assert(ien.start() == 100);
            assert(ien.end() == 100);

            valid = ien.advance();
            assert(valid);
            assert(ien.start() == bm::id_max-1);
            assert(ien.end() == bm::id_max-1);

            valid = ien.advance();
            assert(!valid);
            bv.optimize();

        } // for pass
    }

    cout << "interval_enumerator +N stress test" << endl;
    {
        unsigned delta_max = 65536;
        for (unsigned inc = 1; inc < delta_max; ++inc)
        {
            if (inc % 256 == 0)
                cout << "\rinc = " << inc << " of " << delta_max << flush;
            bvect bv;
            bvect bv_c;
            bvect::size_type test_max = 65535 * 256;

            for (bvect::size_type i = 0; i < test_max; i+=inc)
                bv.set(i);

            for (unsigned pass = 0; pass < 2; ++pass)
            {
                bm::interval_enumerator<bvect> ien(bv);
                while (ien.valid())
                {
                    auto from = ien.start();
                    auto to = ien.end();
                    bv_c.set_range(from, to);
                    if (!ien.advance())
                        break;
                }
                bool eq = bv.equal(bv_c);
                assert(eq);

                bv.optimize();
                bv_c.clear();
            } // for pass
        } // for inc

    }




    cout << "\n----------------------------- IntervalEnumeratorTest() OK" << endl;
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
        bm::id64_t sub_count1[bm::set_sub_array_size];
        bm::id64_t sub_count2[bm::set_sub_array_size];
        for (unsigned i = 0; i < bm::set_sub_array_size; ++i)
        {
            bcount[i] = 0, sub_count1[i] = sub_count2[i] = 0;
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
            assert(rank == 65536 - 6 - bm::rs3_border1 - 1);

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
void CountRangeTest()
{
    cout << "---------------------------- CountRangeTest..." << endl;

    cout << "Stage 0" << endl;
    {{
        bvect bv1 { 0, 1 };
        
        bvect::rs_index_type bc_arr;
        bv1.build_rs_index(&bc_arr);
        assert(bc_arr.count() == 2);
        
        assert(bc_arr.rcount(512) == 2);
        for (bvect::size_type i = 0; i < bm::set_total_blocks; ++i)
        {
            assert(bc_arr.rcount(i) == 2);
        } // for
        
        VerifyCountRange(bv1, bc_arr, 0, 200000);
        
        bv1.optimize();
        bvect::rs_index_type bc_arr1;
        bv1.build_rs_index(&bc_arr1);
        
        for (bvect::size_type i = 0; i < bm::set_total_blocks; ++i)
        {
            assert(bc_arr1.rcount(i) == 2);
        } // for
        
        VerifyCountRange(bv1, bc_arr1, 0, 200000);
    }}

    cout << "Stage 1" << endl;
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
        VerifyCountRange(bv1, bc_arr1, bm::id_max-2000, bm::id_max-1);
    }}

    cout << "Stage 2" << endl;
    {{
        bvect bv1 { 0, 1, 65535+10, 65535+20, 65535+21, bm::id_max-100};
        
        bvect::rs_index_type bc_arr;
        bv1.build_rs_index(&bc_arr);

        assert(bc_arr.rcount(0) == 2);
        assert(bc_arr.rcount(1) == 5);
        assert(bc_arr.rcount(bm::set_total_blocks-1) == 6);

        for (bvect::size_type i = 2; i < bm::set_total_blocks-1; ++i)
        {
            assert(bc_arr.rcount(i) == 5);
        } // for
        
        VerifyCountRange(bv1, bc_arr, bm::id_max-1, bm::id_max-1);
        for (unsigned i = 0; i < 2; ++i)
        {
            cout << "Pass " << i << endl;
            cout << "1 ..." << endl;
            VerifyCountRange(bv1, bc_arr, 0, 200000);
            cout << "2 ..." << endl;
            VerifyCountRange(bv1, bc_arr, bm::id_max-200, bm::id_max-1);
            // check within empty region
            cout << "3 ..." << endl;
            VerifyCountRange(bv1, bc_arr, bm::id_max/2-200, bm::id_max/2+2000);

            bv1.optimize();
            bv1.build_rs_index(&bc_arr);
        }
    }}

    cout << "Stage 3. check inverted bvector" << endl;
    {{
            bvect bv1;
        
            bv1.invert();

            bvect::rs_index_type bc_arr;
            bv1.build_rs_index(&bc_arr);
            auto cnt1 = bv1.count();
            auto cnt2 = bc_arr.count();
            assert(cnt1 == cnt2);

            cout << "1 ..." << endl;
            VerifyCountRange(bv1, bc_arr, bm::id_max-1, bm::id_max-1);
            cout << "2 ..." << endl;
            VerifyCountRange(bv1, bc_arr, 0, 200000);
            cout << "3 ..." << endl;
            VerifyCountRange(bv1, bc_arr, bm::id_max-2000, bm::id_max-1);
            cout << "4 ..." << endl;
            VerifyCountRange(bv1, bc_arr, bm::id_max/2-2000, bm::id_max/2+2000);
            cout << "Done!" << endl;

    }}
    
    cout << "---------------------------- CountRangeTest OK" << endl;
}

// -----------------------------------------------------------------------

bvect::size_type
from_arr[] = { 0, 0, 0,  0,  7,  1,   0,     65535, 65535,   0,     0,                 bm::id_max/2-2001, bm::id_max-2000};
bvect::size_type
to_arr[]   = { 0, 1, 16, 32, 31, 255, 65535, 65536, 65536*2, 65537, bm::id_max/2-2000, bm::id_max/2+2000, bm::id_max-1};


static
void verify_all_one_ranges(const bvect& bv, bool all_one)
{
    size_t fr_size = sizeof(from_arr) / sizeof(from_arr[0]);
    size_t to_size = sizeof(to_arr) / sizeof(to_arr[0]);

    assert(fr_size == to_size);

    for (unsigned t = 0; t < fr_size; ++t)
    {
        bool one_test, one_test_cnt, any_one_test, is_int;
        bvect::size_type from(from_arr[t]), to(to_arr[t]);

        one_test = bv.is_all_one_range(from, to);
        any_one_test = bv.any_range(from, to);
        is_int = bm::is_interval(bv, from, to);
        if (all_one)
        {
            assert(one_test);
            assert(any_one_test);
            if (is_int)
            {
                assert(from == 0 || bv.test(from-1) == false);
                assert(bv.test(to+1) == false);
            }
        }
        else
        {
            auto cnt = bv.count_range(from, to);
            one_test_cnt = (cnt == to - from + 1);
            assert(one_test_cnt == one_test);
            if (cnt)
            {
                assert(any_one_test);
            }
            else
            {
                assert(!any_one_test);
            }
            if (one_test_cnt)
            {
                if (!is_int)
                {
                    bool l, r;
                    if (to < bm::id_max-1)
                        r = !bv.test(to+1);
                    else
                        r = false;
                    if (from)
                        l = !bv.test(from-1);
                    else
                        l = false;
                    assert(l==false || r==false);
                }
            }
        }

        // [from-1, to] range check
        //
        if (from)
        {
            --from;
            one_test = bv.is_all_one_range(from, to);
            any_one_test = bv.any_range(from, to);
            is_int = bm::is_interval(bv, from, to);
            if (all_one)
            {
                assert(one_test);
                assert(any_one_test);
                if (is_int)
                {
                    assert(from == 0 || bv.test(from-1) == false);
                    assert(bv.test(to+1) == false);
                }
            }
            else
            {
                auto cnt = bv.count_range(from, to);
                one_test_cnt = (cnt == to - from + 1);
                assert(one_test_cnt == one_test);
                if (cnt)
                {
                    assert(any_one_test);
                }
                else
                {
                    assert(!any_one_test);
                }
                if (!is_int && one_test_cnt)
                {
                    bool l, r;
                    if (to < bm::id_max-1)
                        r = bv.test(to+1);
                    else
                        r = false;
                    if (from)
                        l = bv.test(from-1);
                    else
                        l = false;
                    assert(l || r);
                }

            }
            ++from;
        }
        // [from, to+1] range check
        //
        if (to < bm::id_max-1)
        {
            ++to;
            one_test = bv.is_all_one_range(from, to);
            any_one_test = bv.any_range(from, to);
            is_int = bm::is_interval(bv, from, to);
            if (all_one)
            {
                assert(one_test);
                assert(any_one_test);
                if (is_int)
                {
                    assert(from == 0 || bv.test(from-1) == false);
                    assert(bv.test(to+1) == false);
                }
            }
            else
            {
                auto cnt = bv.count_range(from, to);
                one_test_cnt = (cnt == to - from + 1);
                assert(one_test_cnt == one_test);
                if (cnt)
                {
                    assert(any_one_test);
                }
                else
                {
                    assert(!any_one_test);
                }
                if (!is_int && one_test_cnt)
                {
                    bool l, r;
                    if (to < bm::id_max-1)
                        r = bv.test(to+1);
                    else
                        r = false;
                    if (from)
                        l = bv.test(from-1);
                    else
                        l = false;

                    assert(l || r);
                }
            }
            --to;
        }
    } // for t
}

static
void IsAllOneRangeTest()
{
    cout << "---------------------------- IsAllOneRangeTest()" << endl;

    cout << "Check empty bvector" << endl;
    {{
        bvect bv1;
        verify_all_one_ranges(bv1, false);
        bv1.set(0);
        bv1.clear(0);
        verify_all_one_ranges(bv1, false);
        IntervalsCheck(bv1);
    }}

    cout << "Check inverted bvector" << endl;
    {{
        bvect bv1;
        bv1.invert();

        verify_all_one_ranges(bv1, true);
        IntervalsCheck(bv1);
    }}

    cout << "Check set ranges" << endl;
    {{
        size_t fr_size = sizeof(from_arr) / sizeof(from_arr[0]);
        size_t to_size = sizeof(to_arr) / sizeof(to_arr[0]);
        assert(fr_size == to_size);

        for (unsigned t = 0; t < fr_size; ++t)
        {
            bvect::size_type from(from_arr[t]), to(to_arr[t]);

            {
                bvect bv1;
                bvect bv2(bm::BM_GAP);

                if (to - from < 65536*10)
                {
                    for (bvect::size_type i = from; i <= to; ++i)
                        bv1.set(i);
                }
                else
                {
                    bv1.set_range(from, to);
                }
                bv2.set_range(from, to);

                bool one_test1 = bv1.is_all_one_range(from, to);
                bool one_test2 = bv2.is_all_one_range(from, to);
                bool any_one1 = bv1.any_range(from, to);
                bool any_one2 = bv2.any_range(from, to);

                assert(one_test1 == one_test2);
                assert(any_one1 == any_one2);
                auto cnt = bv1.count_range(from, to);
                if (any_one1)
                {
                    assert(cnt);
                }
                else
                {
                    assert(!cnt);
                }

                verify_all_one_ranges(bv1, false);
                verify_all_one_ranges(bv2, false);

                IntervalsCheck(bv1);
                IntervalsCheck(bv2);
                IntervalsEnumeratorCheck(bv1);
                IntervalsEnumeratorCheck(bv2);
            }
        } // for t

    }}

    cout << "---------------------------- IsAllOneRangeTest() OK" << endl;
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
        for (bvect::size_type j = base; j < base+max_bits; j += inc)
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

    bvect::size_type base_idx = bvect::size_type(bm::id_max32)+1;
    {
        bvect bv;
        optimize_fill(bv, base_idx, 1, bm::gap_max_bits, true);
        
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
        optimize_fill(bv, base_idx, 100, bm::gap_max_bits, true);
        optimize_fill(bv, base_idx, 100, bm::gap_max_bits, false);

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
        optimize_fill(bv, base_idx, 1, bm::gap_max_bits, true);
        
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
        optimize_fill(bv, base_idx, 1000, bm::gap_max_bits, true);
        optimize_fill(bv, base_idx, 1000, bm::gap_max_bits, false);

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
        optimize_fill(bv, base_idx+0, 1000, bm::gap_max_bits, true);
        optimize_fill(bv, base_idx+1, 1, bm::gap_max_bits, false);

        bvect::statistics st1;
        bv.calc_stat(&st1);
        
        assert(st1.bit_blocks == 0);
        assert(st1.gap_blocks == bm::set_sub_array_size);
        assert(st1.gaps_by_level[0] == 0);
        assert(st1.gaps_by_level[1] == bm::set_sub_array_size);
        assert(st1.ptr_sub_blocks == 1);
        
        bv.optimize(tb, bvect::opt_compress, &st1);

        assert(st1.bit_blocks == 0);
        assert(st1.gap_blocks == 1);
        assert(st1.ptr_sub_blocks == 1);
        assert(st1.gaps_by_level[0] == 1);
        assert(st1.gaps_by_level[1] == 0);
    }


    cout << "---------------------------- Bvector Optimize test OK" << endl;
}

static
void RankFindTest()
{
    cout << "---------------------------- Find Rank test" << endl;

    bvect::size_type base_idx = bvect::size_type(bm::id_max32)+1;
    assert(base_idx > bm::id_max32);

    {
        bvect bv1;
        bv1[base_idx+30] = true;
        bv1[base_idx+65534] = true;

        for (unsigned i = 0; i < 2; ++i)
        {
            bvect::rs_index_type bc_arr1;
            bv1.build_rs_index(&bc_arr1);

            bool rf1, rf2, rf3;
            bvect::size_type pos, pos1;
            rf1 = bv1.find_rank(1, 20, pos);
            rf3 = bv1.find_rank(1, 20, pos1, bc_arr1);
            assert(rf1);
            assert(rf3);
            assert(pos == base_idx+30);
            assert(pos1 == base_idx+30);

            rf2 = bv1.find_rank(2, base_idx+30, pos);
            rf3 = bv1.find_rank(2, base_idx+30, pos1, bc_arr1);
            assert(rf2);
            assert(rf3);
            assert(pos == base_idx+65534);
            assert(pos1 == base_idx+65534);
            bv1.optimize();
        }
    }
    
    {
        bvect bv1;
        bv1.set_range(base_idx, bm::id_max-1);
        bvect::rs_index_type rsi;
        bv1.build_rs_index(&rsi);
        
        bool rf1, rf3;
        bvect::size_type pos, pos1;
        rf1 = bv1.find_rank(1, 20, pos);
        rf3 = bv1.find_rank(1, 20, pos1, rsi);
        assert(rf1);
        assert(rf3);
        assert(pos == base_idx);
        assert(pos1 == base_idx);
    }

    cout << "Test bvector<>::select()" << endl;
    {
        bvect::size_type r(3*65536*256-1), pos;

        bvect bv;
        bv.invert();
        bv.optimize();

        bvect::rs_index_type rs_idx;
        bv.build_rs_index(&rs_idx);

        bvect::size_type i;
        for (i = 0; i <= r; ++i)
        {
            bvect::size_type ri = bv.rank(i, rs_idx);
            assert(ri == (i+1));
            bool f = bv.select(ri, pos, rs_idx);
            assert(f);
            assert(pos == i);
        }
    }

    
    cout << "Find Rank test stress 1\n" << endl;
    
    {
        const bvect::size_type max_size = base_idx+200000;
        const bvect::size_type min_size = base_idx-20000;
        bvect bv1;
        for (bvect::size_type i = base_idx; i < max_size; i += (unsigned)rand()%5)
        {
            bv1.set(i);
        }
        bvect::rs_index_type bc_arr1;
        bv1.build_rs_index(&bc_arr1);

        
        for (bvect::size_type i = max_size; i > min_size; --i)
        {
            bool rf1, rf3;
            bvect::size_type pos, pos1;
            
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
            if ((i & 0xFFFFULL) == 0)
                cout << "\r" << i << "/" << max_size << flush;
        } // for
        cout << endl;
    }
    
    cout << "---------------------------- Find Rank test OK" << endl;
}

extern "C" {
    static
    int bit_decode_func(void* handle_ptr, bm::id64_t bit_idx)
    {
        std::vector<bvect::size_type>* vp = (std::vector<bvect::size_type>*)handle_ptr;
        vp->push_back(bit_idx);
        return 0;
    }
    
    static
    int bit_decode_func2(void* /*handle_ptr*/, bm::id64_t bit_idx)
    {
        static bm::id64_t prev = 0;
        if (bit_idx <= prev)
        {
            if (prev)
            {
                cerr << "Incorrect loading sequence prev=" << prev <<
                         " next=" << bit_idx << endl;
                assert(0); exit(1);
            }
        }
        else
        {
            assert(prev+1 == bit_idx);
        }
        
        const id64_t limit = 65536ULL * 256 * 2ULL;
        if (bit_idx > limit)
        {
            return -1;
        }
        /*
        std::vector<bvect::size_type>* vp = (std::vector<bvect::size_type>*)handle_ptr;
        vp->push_back(bit_idx);
        */
        prev = bit_idx;
        return 0;
    }
} // extern C

template<class BV>
void VisitorAllRangeTest(const BV& bv, typename BV::size_type step)
{
    cout << ".... VisitorAllRangeTest()";

    typename BV::size_type left, right, next, i_end;
    bool non_empty = bv.find_range(left, right);
    if (!non_empty)
        return;

    auto drange = right - left;
    if (!drange)
        drange = 256;
    if (!step)
    {
        unsigned factor = (unsigned)rand()%32;
        if (!factor) factor = 10;
        step = drange / factor;
    }
    if (!step)
        step = 1;
    cout << "   step=" << step << endl;


    auto pcnt = 4;
    for (auto i = left; i <= right; i+=step, right-=step)
    {
        {
            bvect bv_control;
            bm::bit_vistor_copy_functor<bvect> func(bv_control);
            bm::for_each_bit_range(bv, i, right, func);

            bvect bv2;
            bv2.copy_range(bv, i, right);

            bool eq = bv2.equal(bv_control);
            assert(eq);
        }
        next = bv.get_next(i);
        if (next)
        {
            auto delta = next - i;
            if (delta > 32)
            {
                i += delta / 2;
            }
            else
            if (delta == 1)
            {
                bool f = bm::find_interval_end(bv, next, i_end);
                if (f)
                {
                    delta = i_end - i;
                    if (delta > 4)
                        i += delta / 2;
                }
                else
                {
                    assert(!bv.test(i));
                }
            }
        }
        if (!pcnt)
        {
            cout << "\r" << i << " / " << right << flush;
            pcnt = 4;
        }
        --pcnt;

    } // for i
    cout << endl;
/*
    pcnt = 4;
    for (; left <= right; left+=step, --right)
    {
        {
            bvect bv_control;
            bm::bit_vistor_copy_functor<bvect> func(bv_control);
            bm::for_each_bit_range(bv, left, right, func);

            bvect bv2;
            bv2.copy_range(bv, left, right);

            bool eq = bv2.equal(bv_control);
            assert(eq);
        }
        next = bv.get_next(left);
        if (next)
        {
            auto delta = next - left;
            if (delta > 128)
            {
                left += delta / 2;
            }
            else
            if (delta == 1)
            {
                bool f = bv.find_interval_end(left, i_end);
                if (f)
                {
                    delta = i_end - left;
                    if (delta > 4)
                        left += delta / 2;
                }
            }
        }
        if (!pcnt)
        {
            cout << "\r" << left << " / " << right << flush;
            pcnt = 4;
        }
        --pcnt;
    } // for i
    */
    cout << endl;
}



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
        std::vector<bvect::size_type> v1;
        
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
        
        bv1.set(bm::id_max-1, true);
        bv1.set(bm::id_max-1, false);
        
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
                            20000000, bm::id_max-1 };
        bvect bv2;
        std::vector<bvect::size_type> v1;
        
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
        std::vector<bvect::size_type> v1;
        
        generate_bvector(bv1, bm::id_max32/4, false);
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
void GetNextTest()
{
   cout << "-------------------------------------------- GetNextTest()" << endl;
   
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
       bvect::size_type pos;
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
       bvect::size_type pos;
       bvect::size_type from;
       from = bm::id_max / 2;
       
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
       bv[from] = true;
       bv[from+1] = true;
       found = bv.find(1, pos);
       if (!found || pos != from)
       {
           cout << "4. find() failed " << pos << " " << found << endl;
           exit(1);
       }
       found = bv.find(from, pos);
       if (!found || pos != from)
       {
           cout << "5. find() failed " << pos << " " << found << endl;
           exit(1);
       }
       found = bv.find(from+1, pos);
       if (!found || pos != from+1)
       {
           cout << "6. find() failed " << pos << " " << found << endl;
           exit(1);
       }
       found = bv.find_reverse(pos);
       assert(found && pos == from+1);

       found = bv.find(from+2, pos);
       if (found)
       {
           cout << "7. find() failed "<< endl;
           exit(1);
       }
       bv[from+1] = false;
       found = bv.find_reverse(pos);
       assert(found && pos == from);


   }

   {
       bvect  bv;
       bool found;
       bvect::size_type from, to;
       from = bm::id_max / 2;
       to = from + 20000000;
       
       bv.set_range(from, to);
       bvect::size_type pos;
       found = bv.find_reverse(pos);
       assert(found && pos == to);

       bv.optimize();
       found = bv.find_reverse(pos);
       assert(found && pos == to);
       
       bv[bm::id_max-1] = true;
       found = bv.find_reverse(pos);
       assert(found && pos == bm::id_max-1);
       
       bv[bm::id_max-1] = false;
       found = bv.find_reverse(pos);
       assert(found && pos == to);

       bv.set_range(from, to, false);
       found = bv.find_reverse(pos);
       assert(!found);

       found = bv.find(0, pos);
       assert(!found);
   }
   
   {
       bvect  bv;
       bool found;
       bvect::size_type pos;
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

   const bvect::size_type BITVECT_SIZE = bvect::size_type(bm::id_max32)+100; //100000000 * 2;

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


          bvect_full1.set_bit(BITVECT_SIZE-1);
          bvect_min1.set_bit(BITVECT_SIZE-1);

          auto nbit1 = bvect_full1.get_first();
          auto nbit2 = bvect_min1.get_first();

          if (nbit1 != nbit2)
          {
             cout << "1. get_first failed() " <<  nbit1 << " " << nbit2 << endl;
             exit(1);
          }
          nbit1 = bvect_full1.get_next(nbit1);
          nbit2 = bvect_min1.get_next(nbit2);
          if ((nbit1 != nbit2) || (nbit1 != BITVECT_SIZE-1))
          {
             cout << "2. get_next failed() " <<  nbit1 << " " << nbit2 << endl;
             exit(1);
          }
       }


       {
           bvect       bvect_full1;
           bvect_mini  bvect_min1(BITVECT_SIZE);

           bvect_full1.set_new_blocks_strat(i ? bm::BM_GAP : bm::BM_BIT);

           bvect_full1.set_bit(256);
           bvect_min1.set_bit(256);

           bvect_full1.set_bit(BITVECT_SIZE-65536);
           bvect_min1.set_bit(BITVECT_SIZE-65536);

           bvect::size_type nbit1 = bvect_full1.get_first();
           bvect::size_type nbit2 = bvect_min1.get_first();

           if (nbit1 != nbit2)
           {
              cout << "get_first failed " <<  nbit1 << " " << nbit2 << endl;
              exit(1);
           }

           bvect::size_type last_found = 0;
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
         
           bvect::size_type pos = 0;
           bool found = bvect_full1.find_reverse(pos);
           assert(found && pos == last_found);
       }

   }// for

   cout << "-------------------------------------------- GetNextTest()" << endl;

}

static
void BvectorIncTest()
{
    cout << "---------------------------- Bvector inc test" << endl;
    bvect::size_type pos1 = bvect::size_type(bm::id_max32)+1;
    bvect::size_type pos2 = bvect::size_type(bm::id_max)-1;

    
    
    {
        bvect bv1;
        bool b;
        
        b = bv1.inc(pos1);
        assert(!b);
        assert(bv1.count()==1);
        b = bv1.inc(pos1);
        assert(b);
        assert(bv1.count()==0);

        b = bv1.inc(pos2);
        assert(!b);
        assert(bv1.count()==1);
        b = bv1.inc(pos2);
        assert(b);
        assert(bv1.count()==0);
    }
    
    {
        bvect bv1(BM_GAP);
        bool b;
        
        assert(bv1.count()==0);
        
        b = bv1.inc(pos1);
        assert(!b);
        cout << bv1.count() << endl;
        assert(bv1.count()==1);
        b = bv1.inc(pos1);
        assert(b);
        assert(bv1.count()==0);

        b = bv1.inc(pos2);
        assert(!b);
        b = bv1.inc(pos2);
        assert(b);
    }

    {
        bvect bv1(BM_GAP);
        bool b;

        bv1.flip();
        
        b = bv1.inc(pos1);
        assert(b);
        b = bv1.inc(pos1);
        assert(!b);
        
        b = bv1.inc(pos2);
        assert(b);
        b = bv1.inc(pos2);
        assert(!b);
    }

    cout << "---------------------------- Bvector inc test OK" << endl;
}

template<class BV>
void DetailedCompareBVectors(const BV& bv1, const BV& bv2)
{
    bvect::counted_enumerator en1 = bv1.first();
    bvect::counted_enumerator en2 = bv2.first();
    
    for (; en1.valid(); ++en1)
    {
        assert(en2.valid());
        
        auto i1 = *en1;
        auto i2 = *en2;
        
        if (i1 != i2)
        {
            auto nb1 = (i1 >>  bm::set_block_shift);
            auto nb2 = (i2 >>  bm::set_block_shift);
            auto ii1 = nb1 >> bm::set_array_shift;
            auto jj1 = nb1 &  bm::set_array_mask;
            auto ii2 = nb2 >> bm::set_array_shift;
            auto jj2 = nb2 &  bm::set_array_mask;

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

    int cmp = bv1.compare(bv2);
    if (cmp != 0)
    {
        cerr << "Compare (1-2) discrepancy! " << cmp << endl;
    }
    cmp = bv2.compare(bv1);
    if (cmp != 0)
    {
        cerr << "Compare (2-1) discrepancy! " << cmp << endl;
    }

    cout << "Detailed compare OK (no difference)." << endl;
}



static
void BvectorBulkSetTest()
{
    cout << "---------------------------- Bvector BULK set test" << endl;
    {
        bvect::size_type ids[] = { 0 };
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
        bvect::size_type ids[] = {65535, bm::id_max-1 };
        bvect::size_type cnt;
        bvect bv2;
        bv2.set(&ids[0], sizeof(ids)/sizeof(ids[0]));
 
        cnt = bv2.count();
        cout << cnt << endl;

        assert(cnt == sizeof(ids)/sizeof(ids[0]));
        assert(bv2.test(ids[0]));
        assert(bv2.test(ids[1]));
    }

    // set bits in FULL vector
    {
        bvect::size_type ids[] = {0, 10, 65535, bm::id_max-1, bm::id_max-2 };
        bvect::size_type cnt;
        bvect bv2;
        bv2.invert();
        struct bvect::statistics st1, st2;
        bv2.calc_stat(&st1);

        bv2.set(&ids[0], sizeof(ids)/sizeof(ids[0]));

        bv2.calc_stat(&st2);
        assert(st1.bit_blocks == st2.bit_blocks);
        assert(st1.gap_blocks == st2.gap_blocks);

        cnt = bv2.count();

        assert(cnt == bm::id_max);
    }

    // test correct sizing
    {
        bvect::size_type ids[] = {bm::id_max32+100 };
        bvect bv1;
        bv1.resize(10);

        bv1.set(&ids[0], sizeof(ids)/sizeof(ids[0]));
        assert(bv1.size()==bm::id_max32+100+1);
        bv1.keep(&ids[0], sizeof(ids)/sizeof(ids[0]));
        cout << bv1.size() << endl;
        assert(bv1.size()==bm::id_max32+100+1);
        assert(bv1.count()==sizeof(ids)/sizeof(ids[0]));
    }

    {
        bvect::size_type ids[] = {65536, bm::id_max-1, 65535 };
        bvect bv1, bv2;

        for (size_t i = 0; i < sizeof(ids)/sizeof(ids[0]); ++i)
            bv1.set(ids[i]);
 
        bv2.set(&ids[0], sizeof(ids)/sizeof(ids[0]));
        int cmp = bv1.compare(bv2);
        assert(cmp==0);
 
        bv2.keep(&ids[0], sizeof(ids)/sizeof(ids[0]));
        cmp = bv1.compare(bv2);
        assert(cmp==0);
    }

    {
        bvect::size_type ids[] =
            { 0, 1, 2, 3, 4, 5, 256, 1024, 1028, 256000, bm::id_max-100 };
     
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
            bvect::size_type keep_cnt = sizeof(ids)/sizeof(ids[0]);
            bv_inv.keep(&ids[0], keep_cnt, bm::BM_SORTED);
            auto cnt_inv2 = bv_inv.count();
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

    bvect::size_type vector_from = bvect::size_type(bm::id_max32) - 18*65536;
    bvect::size_type vector_to = bvect::size_type(bm::id_max32) + 128*65536;

    cout << "bvector-set volume test 1" << endl;
    {
        std::vector<bvect::size_type> v1, v2, v3;
        generate_test_vectors(v1, v2, v3, vector_from, vector_to);

        for (unsigned k = 0; k < 2; ++ k)
        {
            bvect bvu, bvuc;
            bvect bv1, bv2, bv3, bv11;
            bvect bv1c, bv2c, bv3c;
     
            {
                bvect::bulk_insert_iterator iit(bv11);
                for (bvect::size_type i = 0; i < v1.size(); ++i)
                {
                    bv1c.set(v1[i]);
                    iit = v1[i];
                }
                iit.flush();
            }
     
            for (bvect::size_type i = 0; i < v2.size(); ++i)
                bv2c.set(v2[i]);
            for (bvect::size_type i = 0; i < v3.size(); ++i)
                bv3c.set(v3[i]);
     
            // union of 3 vectors
            bvuc = bv1c;
            bvuc |= bv2c;
            bvuc |= bv3c;

            bv1.set(&v1[0], bvect::size_type(v1.size()));
            bv2.set(&v2[0], bvect::size_type(v2.size()));
            bv3.set(&v3[0], bvect::size_type(v3.size()));

            // imported union of 3 vectors
            bvu.set(&v1[0], bvect::size_type(v1.size()));
            bvu.set(&v2[0], bvect::size_type(v2.size()));
            bvu.set(&v3[0], bvect::size_type(v3.size()));

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
        assert(vector_from < vector_to);
        unsigned delta_max = 65537;

        bvect bv1, bv2;
        bvect bv1c;
        std::vector<bvect::size_type> v1;

        for (unsigned delta = 1; delta < delta_max; ++delta)
        {
            v1.resize(0);
            bvect::bulk_insert_iterator iit(bv2);
            for (bvect::size_type i = vector_from; i < vector_to; i+=delta)
            {
                v1.push_back(i);
                iit = i;
            }
            iit.flush();
 
            bv1.set(&v1[0], bvect::size_type(v1.size()));
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
                bv3.keep(&v1[0], bvect::size_type(v1.size()));
                bv4.set(&v1[0], bvect::size_type(v1.size()));
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

            //if ((delta & 0xFF) == 0)
                cout << "\r" << delta << "/" << delta_max << flush;
            if (delta < 128)
                delta++;
            else
                delta += delta / 2;

        } // for delta
        cout << endl;
    }

    
    cout << "---------------------------- Bvector BULK set test OK" << endl;
}

static
void GAPTestStress()
{
    cout << "----------------------------------- GAP test stress " << endl;
    const bvect::size_type from = bvect::size_type(bm::id_max32)-65536;
    const bvect::size_type to = bvect::size_type(bm::id_max32)+65536*256;
    assert(from < to);

    unsigned ff_to = 65536*256;
    for (unsigned ff = 64; ff < ff_to; )
    {
        bvect bv0, bv1;
        SimpleGapFillSets(bv0, bv1, from, to, ff);
        bv1.optimize();
        for (bvect::size_type i = from; i < to+1; ++i)
        {
            bool b0 = bv0.test(i);
            bool b1 = bv1.test(i);
            if (b0 != b1)
            {
                cerr << "GAP Test Stress failure!" << " FillFactor=" << ff << " bit=" << i << endl;
                assert(0);
                exit(1);
            }
        } // for i
        //if ((ff & 0xFF) == 0)
        {
            cout << "\r" << ff << "/" << ff_to << flush;
        }
        if (ff <= 64)
            ++ff;
        else
            ff += (ff / 2); 
    } // for j

    cout << "\n----------------------------------- GAP test stress " << endl;
}


// -----------------------------------------------------------------------

static
void TestRandomSubset(const bvect& bv, bm::random_subset<bvect>& rsub, bvect::size_type limit = 0)
{
    bvect bv_subset;
    bvect::size_type bcnt = bv.count();
    if (limit)
        bcnt = limit;

    bvect::size_type samples[] =
      { 0, 1, 2, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, bcnt / 5, bcnt / 4, bcnt / 3, bcnt / 2, (bcnt * 2)/3, bcnt };
    bvect::size_type samples_size = sizeof(samples)/sizeof(*samples);

    cout << "Taking random sub-sets: " << samples_size << endl;
    
    for (bvect::size_type i = 0; i < samples_size; ++i)
    {
        bvect::size_type sample_count = samples[i];
        
        cout << "\r" << i << " / " << samples_size << " take = " << sample_count << flush;
        
        rsub.sample(bv_subset, bv, sample_count);
        bv_subset.optimize();
        if (sample_count > bcnt)
            sample_count = bcnt;

        if (sample_count != bv_subset.count())
        {
            cout << "\nRandom subset failed! sample_count = %u result_count=%u\n"
                 << sample_count
                 << ", " << bv_subset.count()
                 << endl;
            
            exit(1);
        }
        {
            bvect bv_subset_copy(bv_subset);
            {
                bvect bv_set_copy(bv);
                bv_set_copy.invert();
                bv_subset_copy -= bv_set_copy;
            }
            int res = bv_subset_copy.compare(bv_subset);
            if (res != 0)
            {
                printf("\nRandom subset failed! inverted set MINUS error! \n");
                assert(0); exit(1);
            }
        }

        bv_subset -= bv;
        if (bv_subset.count() != 0)
        {
            printf("\nRandom subset failed! Extra bits set! \n");
            assert(0); exit(1);
        }
        
    } // for
    cout << endl;
}

// find last set bit by scan (not optimal)
//
static
bool FindLastBit(const bvect& bv, bvect::size_type& last_pos)
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
inline
void CheckVectors(bvect_mini &bvect_min,
                  bvect      &bvect_full,
                  bvect::size_type size,
                  bool     /*detailed*/)
{
    cout << "\nVectors checking...bits to compare = " << size << endl;

    cout << "Bitcount summary : " << endl;
    bvect::size_type min_count = bvect_min.bit_count();
    cout << "minvector count = " << min_count;
    bvect::size_type count = bvect_full.count();
    bvect::size_type full_count = count;
    cout << "fullvector re-count = " << full_count << endl;
    
    if (min_count != full_count)
    {
        cout << "fullvector count = " << count << endl;
        cout << "Count comparison failed !!!!" << endl;
//        print_stat(bvect_full);
//        DetailedCheckVectors(bvect_min, bvect_full, size);
        assert(0);
        exit(1);
    }

    if (full_count)
    {
        bool any = bvect_full.any();
        if (!any)
        {
            cout << "Anycheck failed!" << endl;
            assert(0);
            exit(1);
        }
    }

    // find_last check
    {
        bvect::size_type pos1 = 0;
        bvect::size_type pos2 = 0;
        bool last_found1 = FindLastBit(bvect_full, pos1);
        bool last_found2 = bvect_full.find_reverse(pos2);
        
        assert(last_found1 == last_found2);
        if (last_found1)
        {
            assert(pos1 == pos2);
        }
    }

    IntervalsCheck(bvect_full);
    IntervalsEnumeratorCheck(bvect_full);

    
    // get_next comparison
    cout << "Positive bits comparison..." << flush;
    bvect::size_type nb_min = bvect_min.get_first();
    bvect::size_type nb_ful = bvect_full.get_first();

    bvect::counted_enumerator en = bvect_full.first();
    bvect::size_type nb_en = *en;
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

//         print_stat(bvect_full);
//         DetailedCheckVectors(bvect_min, bvect_full, size);
         assert(0);
         exit(1);
    }
    CompareEnumerators(en, en1);

    if (full_count)
    {
       bvect::size_type bit_count = 1;
       bvect::size_type en_prev = nb_en;

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

           if (!(bit_count & 0xFF))
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
                // DetailedCheckVectors(bvect_min, bvect_full, size);
               assert(0);
               exit(1);
           }
       } while (en.valid());
       if (bit_count != min_count)
       {
           cout << " Bit count failed."
                << " min = " << min_count
                << " bit = " << bit_count
                << endl;
           assert(0);
           exit(1);
       }
    }

    cout << "OK" << endl;
    return;
}




const unsigned ITERATIONS = 180000;


static
void SimpleRandomFillTest()
{
    std::random_device rd;
    std::mt19937_64 mt_rand(rd());

    bvect::size_type BITVECT_SIZE = bvect::size_type(bm::id_max32) * 2;
    bvect::size_type BITVECT_FROM = bvect::size_type(bm::id_max32) - 65535;
    
    assert(ITERATIONS < BITVECT_SIZE);

    bm::random_subset<bvect> rsub;

    cout << "-------------------------- SimpleRandomFillTest" << endl;

    printf("Test for Random inverted subset.");

    {
        bvect bv;
        bv.invert();
        TestRandomSubset(bv, rsub, 256);
    }


    {
        printf("Simple random fill test 1.");
        bvect_mini   bvect_min(BITVECT_SIZE);
        bvect        bvect_full;
        bvect_full.set_new_blocks_strat(bm::BM_BIT);


        bvect::size_type iter = ITERATIONS / 5;

        cout << "\nSimple Random fill test ITERATIONS = " <<  iter << endl;;

        bvect_min.set_bit(0);
        bvect_full.set_bit(0);

        bvect::size_type i;
        for (i = 0; i < iter; ++i)
        {
            bvect::size_type num = mt_rand() % iter;
            bvect_min.set_bit(BITVECT_FROM + num);
            bvect_full.set_bit(BITVECT_FROM + num);
            if ((i % 1000) == 0) cout << "." << flush;
//            CheckCountRange(bvect_full, 0, num);
//            CheckCountRange(bvect_full, num, num+iter);
        }

        CheckVectors(bvect_min, bvect_full, iter, false);
//        CheckCountRange(bvect_full, 0, iter);

        TestRandomSubset(bvect_full, rsub);

        printf("Simple random fill test 2.");

        for(i = 0; i < iter; ++i)
        {
            bvect::size_type num = mt_rand() % iter;
            bvect_min.clear_bit(BITVECT_FROM + num);
            bvect_full.clear_bit(BITVECT_FROM + num);
        }

        CheckVectors(bvect_min, bvect_full, iter, false);
    }


    {
    printf("\nSimple random fill test 3.\n");
    bvect_mini   bvect_min(BITVECT_SIZE);
    bvect      bvect_full(bm::BM_GAP);


    bvect::size_type iter = ITERATIONS;

    cout << "Simple Random fill test ITERATIONS = " << iter << endl;

    bvect::size_type i;
    for(i = 0; i < iter; ++i)
    {
        bvect::size_type num = mt_rand() % iter;
        bvect_min.set_bit(BITVECT_FROM + num);
        bvect_full.set_bit(BITVECT_FROM + num);
//        CheckCountRange(bvect_full, 0, 65535);
//        CheckCountRange(bvect_full, 0, num);
//        CheckCountRange(bvect_full, num, num+iter);
    }

    CheckVectors(bvect_min, bvect_full, iter, false);

    TestRandomSubset(bvect_full, rsub);

    printf("Simple random fill test 4.");

    for(i = 0; i < iter; ++i)
    {
        bvect::size_type num = mt_rand() % iter;
        bvect_min.clear_bit(BITVECT_FROM + num);
        bvect_full.clear_bit(BITVECT_FROM + num);
//        CheckCountRange(bvect_full, 0, num);
//        CheckCountRange(bvect_full, num, num+iter);
    }

    CheckVectors(bvect_min, bvect_full, iter, false);
//    CheckCountRange(bvect_full, 0, iter);

    TestRandomSubset(bvect_full, rsub);
    }

}


static
void RangeRandomFillTest()
{
    bvect::size_type BITVECT_SIZE = bvect::size_type(bm::id_max32) * 2;
//    bvect::size_type BITVECT_FROM = bvect::size_type(bm::id_max32) - 65535;

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

        FillSets(&bvect_min, &bvect_full, min, max, 0ULL);

        CheckVectors(bvect_min, bvect_full, BITVECT_SIZE, false);
    }
    {
        bvect_mini   bvect_min(BITVECT_SIZE);
        bvect     bvect_full;

        printf("Range Random fill test\n");

        bvect::size_type min = BITVECT_SIZE / 2;
        bvect::size_type max = BITVECT_SIZE / 2 + ITERATIONS;
        if (max > BITVECT_SIZE)
            max = BITVECT_SIZE - 1;

        FillSetsIntervals(&bvect_min, bvect_full, min, max, 4ULL);

        CheckVectors(bvect_min, bvect_full, BITVECT_SIZE, false);
    }
}

inline
void CheckRangeCopy(const bvect& bv, bvect::size_type from, bvect::size_type to)
{
    bvect::size_type f1, l1, f2, l2;
    
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


static
void RangeCopyTest()
{
    cout << "----------------------------------- RangeCopyTest" << endl;

    {
        const bvect::size_type to_max = bvect::size_type(bm::id_max32) + 65536 * bm::set_sub_array_size + 10;
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
            for (bvect::size_type i = 0; i < to_max;)
            {
                CheckRangeCopy(bvect1, i, to_max);
                if (i < 128)
                    ++i;
                else
                {
                    if ((to_max - 256) > i)
                        i += i / 5;
                    else
                        ++i;
                }
                cout << "\r" << i << "/" << to_max << flush;
            }
            cout << "Pass " << k << "-1" << endl;
            for (bvect::size_type i = to_max-1; i > 0; )
            {
                CheckRangeCopy(bvect1, 0, i);
                if (i < 128)
                    --i;
                else
                {
                    i -= i / 3;
                }
                cout << "\r" << i << "/" << to_max << flush;
            }
            
            cout << "Pass " << k << "-2" << endl;
            auto to = to_max;
            for (bvect::size_type i = 0; i < to; )
            {
                CheckRangeCopy(bvect1, i, to_max);
                if (i > 128)
                {
                    i += i / 3;
                    to -= to / 5;
                }
                else
                {
                    i++;
                    --to;
                }
                cout << "\r" << i << "/" << to_max << flush;
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
            const bvect::size_type to_max = bvect::size_type(bm::id_max-1);
            const bvect::size_type from = bvect::size_type(bm::id_max-1)-(200);

            cout << "T1" << endl;
            auto to = to_max;

            for (bvect::size_type i = from; i < to; ++i)
            {
                CheckRangeCopy(bvect1, i, to);
            }
            
            cout << "T2" << endl;
            to = to_max;
            for (bvect::size_type i = to; i > from; --i)
            {
                CheckRangeCopy(bvect1, 0, i);
            }
            
            cout << "T3" << endl;
            to = to_max;
            for (bvect::size_type i = from; i != to; ++i, --to)
            {
                CheckRangeCopy(bvect1, i, to);
            }
            cout << "T4" << endl;
            to = bm::id_max-1 - to_max - 100;
            for (bvect::size_type i = to; i < bm::id_max; ++i)
            {
                CheckRangeCopy(bvect1, i, bm::id_max);
                if ((i & 0xFFFF) == 0)
                    cout << "\r" << i << flush;
            }
            cout << endl;
            to = bm::id_max-1 - to_max - 100;
            for (bvect::size_type i = to; i < bm::id_max-(65536 * 3); i+=65536 * 2)
            {
                bvect1.set(i, false);
            }
            for (bvect::size_type k = 0; k < 2; ++k)
            {
                cout << "T5 pass=" << k << endl;
                to = bm::id_max-1 - to_max - (bm::gap_max_bits/2);
                for (bvect::size_type i = to; i < bm::id_max; ++i)
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
void ComparisonTest()
{
    bvect::size_type BITVECT_SIZE = bvect::size_type(bm::id_max32) * 2;

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

    bvect_min1.set_bit(BITVECT_SIZE/2);
    bvect_min2.set_bit(BITVECT_SIZE/2);

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

    bvect_min1.set_bit(BITVECT_SIZE/2+1);
    bvect_full1.set_bit(BITVECT_SIZE/2+1);

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
    bvect::size_type i;
    for (i = 0; i < 65536; ++i)
    {
        bvect_full1.set_bit(i+BITVECT_SIZE/2);
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
void SerializationTest()
{
    std::random_device rd;
    std::mt19937_64 mt_rand(rd());
    bvect::size_type BITVECT_SIZE = bvect::size_type(bm::id_max32) + (bvect::size_type(bm::id_max32) / 2);

   cout << " ----------------------------------- SerializationTest" << endl;

   cout << "\nSerialization STEP 0" << endl;

   {
    bvect::size_type size = BITVECT_SIZE/6000;


    bvect_mini*   bvect_min1= new bvect_mini(BITVECT_SIZE);
    bvect*        bvect_full1= new bvect();
    bvect*        bvect_full2= new bvect();
    bvect*        bv_target_s= new bvect();

    bvect_full1->set_new_blocks_strat(bm::BM_BIT);
    bvect_full2->set_new_blocks_strat(bm::BM_BIT);

    for(bvect::size_type  i = 0; i < size; ++i)
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
    bv_ser.serialize(*bvect_full1, sermem_buf, &st);
    bvect::size_type slen = sermem_buf.size();

    cout << "Serialized mem_max = " << st.max_serialize_mem
         << " size= " << slen
         << " Ratio=" << (slen*100)/st.max_serialize_mem << "%"
         << endl;

    bm::deserialize(*bvect_full2, sermem_buf.buf());
    operation_deserializer<bvect> od;
    od.deserialize(*bv_target_s, sermem_buf.buf(), tb, set_OR);


    CheckVectors(*bvect_min1, *bvect_full2, size, true);
    CheckVectors(*bvect_min1, *bv_target_s, size, true);


    delete bvect_full2;
    delete bvect_min1;
    delete bvect_full1;
    delete bv_target_s;

    }


   {
    bvect::size_type  size = BITVECT_SIZE/6000;


    bvect_mini*   bvect_min1= new bvect_mini(BITVECT_SIZE);
    bvect*        bvect_full1= new bvect();
    bvect*        bvect_full2= new bvect();
    bvect*        bv_target_s= new bvect();

    bvect_full1->set_new_blocks_strat(bm::BM_BIT);
    bvect_full2->set_new_blocks_strat(bm::BM_BIT);

    bvect_full1->set_bit(BITVECT_SIZE-131072);
    bvect_min1->set_bit(BITVECT_SIZE-131072);
    

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

    CheckVectors(*bvect_min1, *bvect_full2, size, false);
    CheckVectors(*bvect_min1, *bv_target_s, size, false);

    delete bvect_full2;
    delete bvect_min1;
    delete bvect_full1;
    delete bv_target_s;

    }


    cout << "\nSerialization STEP 1." << endl;

    {
    bvect_mini   bvect_min1(BITVECT_SIZE);
    bvect        bvect_full1;

    bvect_full1.set_new_blocks_strat(bm::BM_GAP);
   
    bvect::size_type min = BITVECT_SIZE / 2 - ITERATIONS;
    bvect::size_type max = BITVECT_SIZE / 2 + ITERATIONS;
    if (max > BITVECT_SIZE)
        max = BITVECT_SIZE - 1;

    bvect::size_type len = max - min;

    FillSets(&bvect_min1, &bvect_full1, min, max, 0ull);
    FillSets(&bvect_min1, &bvect_full1, 0ull, len, 5ull);

    // shot some random bits

    bvect::size_type i;
    for (i = 0; i < 10000; ++i)
    {
        bvect::size_type bit = mt_rand() % BITVECT_SIZE;
        bvect_full1.set_bit(bit);
        bvect_min1.set_bit(bit);
    }

    bvect::statistics st;
    bvect_full1.calc_stat(&st);

    unsigned char* sermem = new unsigned char[st.max_serialize_mem];
    //print_stat(bvect_full1);
    
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

    CheckVectors(bvect_min1, bvect_full3, max+10, false);
    CheckVectors(bvect_min1, *bv_target_s, max+10, false);

    delete [] sermem;
    delete bv_target_s;

    }


   cout << "\nStage 2" << endl;

   {

    bvect_mini*   bvect_min1= new bvect_mini(BITVECT_SIZE);
//    bm::bvect_mini*   bvect_min2= new bm::bvect_mini(BITVECT_SIZE);
    bvect*        bvect_full1= new bvect();
    bvect*        bvect_full2= new bvect();

    bvect_full1->set_new_blocks_strat(bm::BM_GAP);
    bvect_full2->set_new_blocks_strat(bm::BM_GAP);

    FillSetsRandomMethod(bvect_min1, bvect_full1, 1ull, BITVECT_SIZE-10, 1);
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
    operation_deserializer<bvect64> od;
    od.deserialize(*bv_target_s,
                   sermem_buf.buf(),
                   0,
                   set_OR);

    CheckVectors(*bvect_min1, *bvect_full2, BITVECT_SIZE, false);
    CheckVectors(*bvect_min1, *bv_target_s, BITVECT_SIZE, false);

    delete bv_target_s;
    delete bvect_full2;
    delete bvect_min1;
    delete bvect_full1;

    }



   cout << "\nStage 3" << endl;

   {

    bvect_mini*   bvect_min1= new bvect_mini(BITVECT_SIZE);
    bvect_mini*   bvect_min2= new bvect_mini(BITVECT_SIZE);
    bvect*        bvect_full1= new bvect();
    bvect*        bvect_full2= new bvect();

    bvect_full1->set_new_blocks_strat(bm::BM_GAP);
    bvect_full2->set_new_blocks_strat(bm::BM_GAP);


    FillSetsRandomMethod(bvect_min1, bvect_full1, 1ull, BITVECT_SIZE, 1);
    FillSetsRandomMethod(bvect_min2, bvect_full2, 1ull, BITVECT_SIZE, 1);
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
        //print_stat(bvt);
        //print_stat(*bvect_full1);
        cout << "Error!" << endl;
        assert(0); exit(1);
    }

    CheckVectors(*bvect_min1, *bvect_full1, BITVECT_SIZE, true);
    CheckVectors(*bvect_min2, *bvect_full2, BITVECT_SIZE, true);

    cout << "Serialized mem_max = " << st.max_serialize_mem
         << " size= " << slen
         << " Ratio=" << (slen*100)/st.max_serialize_mem << "%"
         << endl;

    bvect*  bv_target_s=new bvect(*bvect_full2);
    //print_stat(*bv_target_s);

    //print_stat(*bvect_full2);

    bvect*  bvect_full3= new bvect();
    *bvect_full3 = *bvect_full1;
    *bvect_full3 |= *bvect_full2;
    // bug in test here
    //CheckVectors(*bvect_min2, *bvect_full3, BITVECT_SIZE, true);


    bm::deserialize(*bvect_full2, sermem);

    operation_deserializer<bvect64> od;
    od.deserialize(*bv_target_s,
                   sermem,
                   0,
                   set_OR);
    delete [] sermem;
    
    CheckVectors(*bvect_min1, *bvect_full1, BITVECT_SIZE, true);
    // bug in test, commented out
//    CheckVectors(*bvect_min1, *bvect_full3, BITVECT_SIZE, true);

    bvect_min2->combine_or(*bvect_min1);
    delete bvect_min1;
    
    if (*bvect_full2 != *bvect_full3)
    {
        //print_stat(*bvect_full2);
        //print_stat(*bvect_full3);

        cout << "Error!" << endl;
        assert(0);
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

   cout << "\nStage 4. " << endl;

   {
    bvect::size_type size = BITVECT_SIZE/3;


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

    for (i = 65536; i < 65536+65000; ++i)
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

    print_stat(cout,*bvect_full1);


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

    operation_deserializer<bvect64> od;
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


static
void DesrializationTest2()
{
   cout << "\n-------------------------------------- DesrializationTest2()" << endl;
   bvect::size_type BITVECT_SIZE = bvect::size_type(bm::id_max32) * 3;


   bvect  bv1, bv2;
   bvect  bvtotal;
   bvect::size_type size = BITVECT_SIZE - 10;
   BM_DECLARE_TEMP_BLOCK(tb)
   unsigned repetitions = 10;
   bvect::size_type i;


   for (i = bm::id_max32; i < (bm::id_max32+165536); i+=2)
   {
      bv1.set_bit(i);
   }

   bv1.optimize(tb);

   struct bvect::statistics st1;
   bv1.calc_stat(&st1);

   std::vector<unsigned char> sermemv(st1.max_serialize_mem);
   
   size_t slen2 = bm::serialize(bv1, sermemv.data(), tb);
   assert(slen2);
   slen2 = 0;

   bm::deserialize(bvtotal, sermemv.data());
    bvect  bv_target_s;

    operation_deserializer<bvect64> od;
    od.deserialize(bv_target_s,
                    sermemv.data(),
                    0,
                    set_OR);

   bvtotal.optimize(tb);
   int res = bvtotal.compare(bv_target_s);
   if (res != 0)
   {
       cout << "Operation deserialization error 1" << endl;
       assert(0);
       exit(1);
   }

   for (i = bm::id_max32+55000; i < bm::id_max32+165536; ++i)
   {
      bv2.set_bit(i);
   }
   bv2.optimize();
   print_stat(cout,bv2);

   struct bvect::statistics st2;
   bv2.calc_stat(&st2);

   std::vector<unsigned char> sermemv2(st2.max_serialize_mem);

   size_t slen = bm::serialize(bv2, sermemv2.data());
   assert(slen);
   slen = 0;

   bm::deserialize(bvtotal, sermemv2.data());
   print_stat(cout,bvtotal);
   //operation_deserializer<bvect64> od;

   od.deserialize(bv_target_s,
                  sermemv2.data(),
                  0,
                  set_OR);
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
        assert(0);
        exit(1);
    }


   bvtotal.clear();
   bv_target_s.clear(false);

   int clcnt = 0;

   for (i = 0; i < repetitions; ++i)
   {
        cout << endl << endl << "Deserialization STEP " << i << endl;

        bvect_mini*   bvect_min1= new bvect_mini(size);
        bvect*        bvect_full1= new bvect();

        FillSetsRandomMethod(bvect_min1, bvect_full1, 1ull, size, 1ull);

       struct bvect::statistics st;
       bvect_full1->calc_stat(&st);

       std::vector<unsigned char> sermemv1(st.max_serialize_mem);

       slen = bm::serialize(*bvect_full1, sermemv1.data(), tb);

       std::vector<unsigned char> smemv(slen);
       ::memcpy(smemv.data(), sermemv1.data(), slen);


        bm::deserialize(bvtotal, smemv.data());
        od.deserialize(bv_target_s,
                                                   smemv.data(),
                                                   0,
                                                   set_OR);
        res = bvtotal.compare(bv_target_s);
        if (res != 0)
        {
            bvect::size_type bit_idx = bv_target_s.get_first();
            cout << bit_idx << " " << bv_target_s.get_next(bit_idx) << endl;;
            print_stat(cout,bv_target_s);
            cout << "Operation deserialization error 2" << endl;
            assert(0);
            exit(1);
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

    cout << "de-serialization test 64-bit wide vectors mix.." << endl;
    {
        for (i = 0; i < repetitions; ++i)
        {
            bvect64 bv0;
            bv1.clear();
            {
                ref_vect vect0, vect1;
                
                generate_vect_simpl0(vect0);
                generate_vect48(vect1);
                
                load_BV_set_ref(cout,bv0, vect0, false);
                load_BV_set_ref(cout,bv1, vect1, false);
            }

            bm::serializer<bvect64> bv_ser;
            bm::serializer<bvect64>::buffer sbuf0, sbuf1;
            bv_ser.serialize(bv0, sbuf0, 0);
            bv_ser.serialize(bv1, sbuf1, 0);

            bvect64 bv_tc;
            bv_tc.bit_or(bv0, bv1, bvect64::opt_none);
            {
                bvect64 bv_tc1(bv0);
                bv_tc1 |= bv1;
                int cmp = bv_tc.compare(bv_tc1);
                assert(cmp==0);
            }
            
            bvect64 bv_tc_and;
            bv_tc_and.bit_and(bv0, bv1, bvect64::opt_none);

            // integrity check
            {
                bvect64 bv_t0;
                bm::deserialize(bv_t0, sbuf0.buf());
                res = bv0.compare(bv_t0);
                assert(res==0);
                // pass 2 to make sure deserialization into the same is ok
                bm::deserialize(bv_t0, sbuf0.buf());
                res = bv0.compare(bv_t0);
                assert(res==0);

                
                bvect64 bv_t1;
                bm::deserialize(bv_t1, sbuf1.buf());
                res = bv1.compare(bv_t1);
                assert(res==0);
                bm::deserialize(bv_t1, sbuf1.buf());
                res = bv1.compare(bv_t1);
                assert(res==0);
            }

            {
                bvect64 bv_t0, bv_t1;
                od.deserialize(bv_t0, sbuf0.buf(), 0, set_OR);
                int cmp = bv0.compare(bv_t0);
                assert(cmp == 0);
                od.deserialize(bv_t1, sbuf1.buf(), 0, set_OR);
                cmp = bv1.compare(bv_t1);
                assert(!cmp);
            }

            {
                bvect64 bv_t;
                od.deserialize(bv_t, sbuf1.buf(), 0, set_OR);
                od.deserialize(bv_t, sbuf0.buf(), 0, set_OR);
                int cmp = bv_tc.compare(bv_t);
                if (cmp)
                {
                    std::cerr << "Serialization intergity check failed!" << std::endl;
                    assert(0); exit(1);
                }
            }
            {
                bvect64 bv_t;
                od.deserialize(bv_t, sbuf0.buf(), 0, set_OR);
                od.deserialize(bv_t, sbuf0.buf(), 0, set_AND);
                int cmp = bv0.compare(bv_t);
                if (cmp)
                {
                    std::cerr << "AND Serialization intergity check failed!" << std::endl;
                    assert(0); exit(1);
                }
            }
            {
                bvect64 bv_t;
                od.deserialize(bv_t, sbuf1.buf(), 0, set_OR);
                od.deserialize(bv_t, sbuf1.buf(), 0, set_AND);
                int cmp = bv1.compare(bv_t);
                if (cmp)
                {
                    std::cerr << "AND Serialization intergity check failed!" << std::endl;
                    assert(0); exit(1);
                }
            }

            {
                bvect64 bv_t;
                od.deserialize(bv_t, sbuf0.buf(), 0, set_OR);
                od.deserialize(bv_t, sbuf1.buf(), 0, set_AND);
                int cmp = bv_tc_and.compare(bv_t);
                if (cmp)
                {
                    std::cerr << "AND Serialization intergity check failed!" << std::endl;
                    assert(0); exit(1);
                }
            }

            
            // TODO: add detection that blocks are actually getting freed as a result
            //
            {
                bvect64 bv_t;
                od.deserialize(bv_t, sbuf0.buf(), 0, set_OR);
                od.deserialize(bv_t, sbuf0.buf(), 0, set_XOR);
                auto cnt = bv_t.count();
                if (cnt)
                {
                    std::cerr << "XOR Serialization intergity check failed!" << std::endl;
                    assert(0); exit(1);
                }
            }
            {
                bvect64 bv_t;
                od.deserialize(bv_t, sbuf1.buf(), 0, set_OR);
                od.deserialize(bv_t, sbuf1.buf(), 0, set_XOR);
                auto cnt = bv_t.count();
                if (cnt)
                {
                    std::cerr << "XOR Serialization intergity check failed!" << std::endl;
                    assert(0); exit(1);
                }
            }
            {
                bvect64 bv_t;
                od.deserialize(bv_t, sbuf1.buf(), 0, set_OR);
                od.deserialize(bv_t, sbuf1.buf(), 0, set_SUB);
                auto cnt = bv_t.count();
                if (cnt)
                {
                    std::cerr << "SUB Serialization intergity check failed!" << std::endl;
                    assert(0); exit(1);
                }
            }
            {
                bvect64 bv_t;
                bm::deserialize(bv_t, sbuf0.buf());
                res = bv0.compare(bv_t);
                assert(res==0);
                
                bvect64 bv_t_c(bv_t);
                
                bm::deserialize(bv_t, sbuf1.buf());
                int cmp = bv_tc.compare(bv_t);
                if (cmp)
                {
                    std::cerr << "Serialization intergity check failed!" << std::endl;
                    assert(0); exit(1);
                }
            }
            cout << "\r" << i << "/" << repetitions << flush;
        } // for i
        cout << endl << "OK" << endl;;
    }

   cout << "\n-------------------------------------- DesrializationTest2() OK" << endl;
}

static
void generate_sparse_bv(bvect& bv, bvect::size_type from,
    bvect::size_type to,
    bvect::size_type step = 65536 / 10)
{
    for (bvect::size_type i = from; true; i += step)
    {
        bv.set(i);
        if (to - step < i)
            break;
    }
}

static
void SparseSerializationTest()
{
    cout << " ----------------------------------- SparseSerializationTest()" << endl;
    BM_DECLARE_TEMP_BLOCK(tb)
        operation_deserializer<bvect> od;
    bm::serializer<bvect> bv_ser(tb);
    bm::serializer<bvect>::buffer sermem_buf;
    bm::serializer<bvect>::buffer sermem_buf_shifted;

    std::vector<std::pair<bvect::size_type, bvect::size_type> > ranges;

    ranges.push_back(std::make_pair(0, 65535 * 255));
    ranges.push_back(std::make_pair(0, 65535 * 255 * 2));
    ranges.push_back(std::make_pair(65535 * 5, 65535 * 255 * 2));
    ranges.push_back(std::make_pair(65535 * 255 / 2, 65535 * 255));
    ranges.push_back(std::make_pair(65535 * 255 / 2, 65535 * 255 * 2)); 
    ranges.push_back(std::make_pair(bm::id_max / 2 - 65535 * 255 / 2, bm::id_max / 2 + 65535 * 255 * 2));
    ranges.push_back(std::make_pair(bm::id_max / 2, bm::id_max / 2 + 65535 * 255 * 2));
    ranges.push_back(std::make_pair(bm::id_max - 65535 * 255 * 2, bm::id_max - 1));

    for (size_t k = 0; k < ranges.size(); ++k)
    {
        bvect::size_type from = ranges[k].first;
        bvect::size_type to = ranges[k].second;

        std::cout << "  range [" << from << ", " << to << "]" << std::endl;

        for (unsigned i = 0; i < 2; ++i)
        {
            bvect bv, bv_shifted;

            generate_sparse_bv(bv, from, to, 65536 / 10);
            generate_sparse_bv(bv_shifted, from + 1, to + 1, 65536 / 10);

            //CheckRangeDeserial(bv, from, to);

            auto bv_cnt = bv.count();

            bv_ser.serialize(bv, sermem_buf, 0);
            bv_ser.serialize(bv_shifted, sermem_buf_shifted, 0);

            const bvect::size_type* cstat = bv_ser.get_compression_stat();

            {
                bvect::size_type sb_from = from / (65536 * 256);
                bvect::size_type sb_to = to / (65536 * 256);
                bvect::size_type sb_cnt = sb_to - sb_from + 1;
                assert(cstat[bm::set_sblock_bienc] == sb_cnt || cstat[bm::set_sblock_bienc] == sb_cnt - 1);
            }

            bvect bv2;
            bm::deserialize(bv2, sermem_buf.buf());
            bool eq = bv.equal(bv2);
            assert(eq);

            {
                bvect bv3;
                od.deserialize(bv3, sermem_buf.buf(), tb, set_OR);
                eq = bv3.equal(bv2);
                assert(eq);
            }
            {
                bvect bv3;
                od.deserialize(bv3, sermem_buf.buf(), tb, set_XOR);
                eq = bv3.equal(bv);
                assert(eq);
                od.deserialize(bv3, sermem_buf.buf(), tb, set_XOR);
                assert(bv3.count() == 0);

                bv3 = bv;
                od.deserialize(bv3, sermem_buf_shifted.buf(), tb, set_XOR);
                assert(bv3.count() == 2 * bv_cnt);
            }
            {
                bvect bv3;
                od.deserialize(bv3, sermem_buf.buf(), tb, set_XOR);
                eq = bv3.equal(bv);
                assert(eq);
                od.deserialize(bv3, sermem_buf.buf(), tb, set_XOR);
                assert(bv3.count() == 0);

                bv3 = bv;
                od.deserialize(bv3, sermem_buf_shifted.buf(), tb, set_XOR);
                assert(bv3.count() == 2 * bv_cnt);
            }
            {
                bvect bv3(bv);
                od.deserialize(bv3, sermem_buf_shifted.buf(), tb, set_SUB);
                eq = bv3.equal(bv);
                assert(eq);
                od.deserialize(bv3, sermem_buf.buf(), tb, set_SUB);
                assert(bv3.count() == 0);
            }
            {
                bvect bv3(bv);
                od.deserialize(bv3, sermem_buf.buf(), tb, set_AND);
                eq = bv3.equal(bv);
                assert(eq);
                od.deserialize(bv3, sermem_buf_shifted.buf(), tb, set_AND);
                assert(bv3.count() == 0);
            }


            {
                bvect bv3(bv);
                auto cnt = od.deserialize(bv3, sermem_buf.buf(), tb, set_COUNT_OR);
                assert(cnt == bv_cnt);

                cnt = od.deserialize(bv3, sermem_buf_shifted.buf(), tb, set_COUNT_OR);
                assert(cnt == 2 * bv_cnt);
            }

            {
                bvect bv3;
                auto cnt = od.deserialize(bv3, sermem_buf.buf(), tb, set_COUNT_XOR);
                assert(cnt == bv_cnt);

                cnt = od.deserialize(bv, sermem_buf.buf(), tb, set_COUNT_XOR);
                assert(cnt == 0);

                bv3 = bv_shifted;
                cnt = od.deserialize(bv3, sermem_buf_shifted.buf(), tb, set_COUNT_XOR);
                assert(cnt == 0);
            }
            {
                bvect bv3(bv);
                auto cnt = od.deserialize(bv3, sermem_buf.buf(), tb, set_COUNT_AND);
                assert(cnt == bv_cnt);
                cnt = od.deserialize(bv3, sermem_buf_shifted.buf(), tb, set_COUNT_AND);
                assert(cnt == 0);
            }


            bv.optimize(tb);
            bv_shifted.optimize(tb);
        } // for
    }

    cout << " ----------------------------------- SparseSerializationTest() OK" << endl;
}



// ---------------------------------------------------------------------------

static
void CheckRangeDeserial(const bvect&     bv,
                        bvect::size_type from,
                        bvect::size_type to)
{
    static unsigned bm_distance = 4;
    assert(from < to);

    cout << " Check Range [" << from << ", " << to << "] = " << (to-from) << endl;

    int max_inc = 256;
    if (to - from > 65536)
        max_inc = 1024;


    bool eq;
    bm::operation_deserializer<bvect> od;

    bm::serializer<bvect> bvs;
    bvs.set_bookmarks(false);
    //bvs.set_bookmarks(true, bm_distance++);

    cout << " bookmarks OFF" << endl;

    for (unsigned pass = 0; pass < 2; ++pass)
    {
        cout << "   pass = " << pass << endl;
        bm::serializer<bvect>::buffer buf;
        cout << "   serialize ..." << flush;
        bvs.serialize(bv, buf);
        cout << "OK" << endl;

        bvect bv_r;
        bv_r.copy_range(bv, from, to);
        auto count_r = bv.count_range(from, to);
        auto count = bv_r.count();
        assert(count == count_r);

        {
            bvect bv_c;
            bm::deserialize(bv_c, buf.data());
            eq = bv.equal(bv_c);
            assert(eq);
        }

        {
            bvect bv_x;
            bv_x.bit_xor(bv, bv_r, bvect::opt_compress);
            count_r = bv_x.count_range(from, to);
            assert(!count_r);
        }

        const unsigned char* sdata = buf.data();

        {

            bvect bv_rd_m;
            bv_rd_m.set_range(from, to);
            od.deserialize(bv_rd_m, sdata, 0, bm::set_AND);
            eq = bv_r.equal(bv_rd_m);
            assert(eq);

            bvect bv_rd;
            od.deserialize_range(bv_rd, sdata, from, to);
            eq = bv_r.equal(bv_rd);
            assert(eq);
        }

        bvect::size_type cnt = 0;

        cout << "      start range" << endl;
        {
            bvect::size_type target = from + 65536 * 2;
            if (target > to)
                target = to;

            bvect bv_rd_m;
            bv_rd_m.set_range(from, target);
            for (bvect::size_type i = from; i <= target; ++cnt)
            {
                bv_r.copy_range(bv, i, target);
                bvect bv_rd;
                od.deserialize_range(bv_rd, sdata, i, target);
                eq = bv_r.equal(bv_rd);
                assert(eq);

                od.deserialize(bv_rd_m, sdata, 0, bm::set_AND);
                eq = bv_rd.equal(bv_rd_m);
                assert(eq);

                {
                    bvect bv_rd2;
                    bm::deserialize_range(bv_rd2, sdata, i, target);
                    eq = bv_rd.equal(bv_rd2);
                    assert(eq);
                }

                auto r = target - i;
                cout << "\r      " << r << "      " << flush;
                if (cnt > 128 && (target - i) > 128)
                {
                    {
                        i += bvect::size_type(rand() % max_inc);
                        target -= bvect::size_type(rand() % max_inc);
                    }
                    bv_rd_m.keep_range(i, target);
                    continue;
                }
                else
                {
                    bv_rd_m.set(i, false);
                    bv_rd_m.set(target, false);
                }
                ++i; --target;
            } // for i
        }

        cout << "\r       " << endl;
        cout << "      whole range (randomized)" << endl;

        cnt = 0;
        bvect::size_type j = to;

        for (bvect::size_type i = from; i <= j; ++i, --j, ++cnt)
        {
            bv_r.copy_range(bv, i, j);
            bvect bv_rd;
            od.deserialize_range(bv_rd, buf.data(), i, j);
            eq = bv_r.equal(bv_rd);
            assert(eq);

            bvect bv_rd_m;
            bv_rd_m.set_range(i, j);
            od.deserialize(bv_rd_m, buf.data(), 0, bm::set_AND);
            eq = bv_rd.equal(bv_rd_m);
            assert(eq);

            auto r = j - i;
            cout << "\r      " << r << "      " << flush;
            // turn on random gallop mode
            if (cnt > 100)
            {
                {
                    i = bvect::size_type(i + unsigned(rand() % max_inc));
                    j = bvect::size_type(j - unsigned(rand() % max_inc));
                }
            }
        } // for i-j

        bvs.set_bookmarks(true, bm_distance++);
        cout << "\n bookmarks ON distance=" << (bm_distance-1) << endl;

    } // for pass (bookmarks)

    cout << "\r       " << endl;

}

static
void RangeDeserializationTest()
{
    cout << "\n------------------------------- RangeDeserializationTest()" << endl;

    cout << "======= BV 48-bit sparse " << endl;
    {
        bvect bv3;  // 48-bit super sparse
        ref_vect vect;
        generate_vect_simpl0(vect);
        load_BV_set_ref(cout, bv3, vect);

        CheckRangeDeserial(bv3, 0, 77777);
        CheckRangeDeserial(bv3, (bm::id_max32 - 10), bm::id_max32 + 10);
        CheckRangeDeserial(bv3, bm::id_max48 / 2 + 1, bm::id_max48 / 2 + 32);
        CheckRangeDeserial(bv3, bm::id_max48 - 65536, bm::id_max48 - 1);
    }

    cout << "============ BV sparse vector" << endl;
    {
        std::vector<std::pair<bvect::size_type, bvect::size_type> > ranges;

        ranges.push_back(std::make_pair(0, 65535 * 255));
        ranges.push_back(std::make_pair(0, 65535 * 255 * 3));
        ranges.push_back(std::make_pair(65535 * 5, 65535 * 255 * 2));
        ranges.push_back(std::make_pair(65535 * 255 / 2, 65535 * 255));
        ranges.push_back(std::make_pair(65535 * 255 / 2, 65535 * 255 * 2));
        ranges.push_back(std::make_pair(bm::id_max / 2 - 65535 * 255 / 2, bm::id_max / 2 + 65535 * 255 * 2));
        ranges.push_back(std::make_pair(bm::id_max / 2, bm::id_max / 2 + 65535 * 255 * 2));
        ranges.push_back(std::make_pair(bm::id_max - 65535 * 255 * 2, bm::id_max - 1));

        for (size_t k = 0; k < ranges.size(); ++k)
        {
            bvect bv;  // generated random

            bvect::size_type from = ranges[k].first;
            bvect::size_type to = ranges[k].second;
            std::cout << "- Vector range [" << from << ", " << to << "]" << std::endl;

            generate_sparse_bv(bv, from, to, 65536 / 10);

            CheckRangeDeserial(bv, to - 65536, to);
            CheckRangeDeserial(bv, from, from + 65536);
            auto mid = (to - from) / 2;
            CheckRangeDeserial(bv, mid - 100, mid + 100);
        }
    }


    // generated random
    cout << "======= BV random generated " << endl;
    {
        bvect bv1;  // generated random

        generate_bvector(bv1, bm::id_max32/4, false);
        CheckRangeDeserial(bv1, 0, 5*65536);
        CheckRangeDeserial(bv1, bm::id_max32/4-(8*65536), bm::id_max32/4);

        bv1.optimize();
        CheckRangeDeserial(bv1, 128*65536, 130*65536);
        CheckRangeDeserial(bv1, bm::id_max32/4-(2*65536), bm::id_max32/4);
    }

    cout << "======= BV 48-bit generated " << endl;
    {
        bvect bv4;  // 48-bit
        {
        ref_vect vect0;
        generate_vect48(vect0);
        load_BV_set_ref(cout, bv4, vect0);
        }
        bv4.optimize();
        CheckRangeDeserial(bv4, 0, 64*65536);
        CheckRangeDeserial(bv4, (bm::id_max32-65535), bm::id_max32+10);
        CheckRangeDeserial(bv4, bm::id_max48-75536, bm::id_max48-1);
    }


        cout << "======= BV Empty " << endl;
        {
            bvect bv_e;
            CheckRangeDeserial(bv_e, 0, 256*65536);
            CheckRangeDeserial(bv_e, bm::id_max32/4-(256*65536), bm::id_max32/4);
        }

        // inverted
        cout << "======= BV inverted " << endl;
        {
            bvect bv_i; // inverted
            bv_i.invert();
            CheckRangeDeserial(bv_i, 0, 256*65536);
            CheckRangeDeserial(bv_i, bm::id_max32/4-(256*65536), bm::id_max32/4);
        }



    cout << "\n------------------------------- RangeDeserializationTest() OK" << endl;
}

// ---------------------------------------------------------------------------


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
void ShiftRight(bvect*  bv, bvect::size_type shift)
{
    bvect bv_tmp;
    {
        bvect::bulk_insert_iterator bi = bv_tmp.inserter();
        bvect::enumerator en = bv->first();
        for (; en.valid(); ++en)
        {
            bvect::size_type v = *en;
            bvect::size_type new_v = v + shift;
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



static
void BvectorShiftTest()
{
    cout << "---------------------------- Bvector SHIFT test" << endl;
    bvect::size_type BITVECT_SIZE = bvect::size_type(bm::id_max32) * 3;


    {
        bvect bv;
        
        bv.set(bm::id_max-1);
        bv.shift_right();
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
        auto cnt0 = bv.count();
        bv.set(bm::gap_max_bits * bm::set_sub_array_size, false);
        auto cnt01 = bv.count();
        assert(cnt01 == cnt0 - 1);

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
        auto cnt = bv.count();
        bool carry_over = bv.shift_right();
        assert(carry_over);
        auto cnt1 = bv.count();
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
        auto idx = bv.get_first();
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
        auto idx = bv.get_first();
        assert(idx == 4278190080-1);
        bv.shift_left();
        idx = bv.get_first();
        assert(idx == 4278190080-2);
    }
    
    {
        bvect bv { 4278190080 };
        bv.optimize();
        bv.shift_left();
        auto idx = bv.get_first();
        assert(idx == 4278190080-1);
        bv.shift_left();
        idx = bv.get_first();
        assert(idx == 4278190080-2);
    }

    {
    std::cout << "Shift-L stress (1 bit shift)..\n";
    bvect::size_type start = bm::id_max-1;
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
 #ifdef DEBUG
        bvect::size_type idx = bv.get_first();
        if(idx != start-1)
        {
            cerr << bv.count() << endl;
            cerr << "Shift-L Failed at idx=" << idx << " != " << start << endl;
            exit(1);
        }
 #endif

        if ((start & 0xFF) == 0)
        {
            f = std::chrono::steady_clock::now();
            auto diff = f - s;
            auto d = std::chrono::duration <double, std::milli> (diff).count();
            cout << "\r" << start << " / " << (bm::id_max - start) << " (" << d << ") " << flush;

            auto idx = bv.get_first();
            assert(idx == start-1);

            bv.calc_stat(&st);
            bcnt = st.bit_blocks + st.gap_blocks;
            assert(bcnt == 1);
            
            s = std::chrono::steady_clock::now();
        }
        if (bm::id_max - start > 100000)
            break;
    }
    cout << "ok.\n";
    }

    {
        std::cout << "Shift-R stress (1 bit shift)..\n";
        bvect::size_type start = 0;
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

            if ((start & 0xFF) == 0)
            {
                f = std::chrono::steady_clock::now();
                auto diff = f - s;
                auto d = std::chrono::duration <double, std::milli> (diff).count();
                cout << "\r" << start << " (" << d << ") " << flush;
                
                auto idx = bv.get_first();
                assert(idx-1 == start);
                
                bv.calc_stat(&st);
                bcnt = st.bit_blocks + st.gap_blocks;
                assert(bcnt == 1);
     
                s = std::chrono::steady_clock::now();
            }
            ++start;
            if (start > 100000)
                break;
        }
        cout << "ok.\n";
    }

    {
        std::cout << "Shift-R stress (large vector shift)..\n";
        bvect bv;
        generate_bvector(bv, BITVECT_SIZE, true);
        bvect bv_control(bv);
        
        unsigned max_shifts = 33;
        for (unsigned i = 0; i < max_shifts; ++i)
        {
            ShiftRight(&bv_control, 1);
            bv.shift_right();
            int cmp = bv.compare(bv_control);
            assert(cmp==0);
            
            cout << "\r" << i << "/" << max_shifts << flush;
        }
    }
    cout << "ok.\n";


    // stress test for shifting aggregator
    //
    cout << "Aggregator based SHIFT-R tests..." << endl;
    {
        const unsigned int REPEATS = 300;

        bvect mask_bv; // mask vector
        mask_bv.init();
        generate_bvector(mask_bv, 75000000, false); // mask is shorter on both ends

        std::vector<bvect> bv_coll1;
        GenerateShiftTestCollection(&bv_coll1, 25ULL, 80000000ULL, true);
        
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
                agg.set_compute_count(false);
                agg.combine_shift_right_and(bv2);
                int cmp = bv1.compare(bv2);
                if (cmp != 0)
                {
                    cerr << "Shift-R compare failure!" << endl;
                    exit(1);
                }
                bvect bv3;
                agg.set_compute_count(true);
                agg.combine_shift_right_and(bv3);
                assert(!bv3.any());
                auto cnt = agg.count();
                auto cnt_c = bv1.count();
                assert(cnt == cnt_c);

            } // for
        }
    }

    cout << "---------------------------- Bvector SHIFT test OK" << endl;
}

// Reference bit insert
static
void BVectorInsert(bvect*  bv, bvect::size_type pos, bool value)
{
    bvect bv_tmp;
    if (pos)
        bv_tmp.copy_range(*bv, 0, pos-1);
    
    {
        bvect::bulk_insert_iterator bi = bv_tmp.inserter();
        bvect::enumerator en = bv->first();
        for (; en.valid(); ++en)
        {
            bvect::size_type v = *en;
            if (v < pos)
                continue;
            bvect::size_type new_v = v + 1;
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
        bv.insert(bm::id_max/2, true);
        assert(bv.test(bm::id_max/2));
        assert(bv.count()==1);
        assert(bv.size() == bm::id_max/2+1);
    }
    
    {
        bvect bv { bm::id_max/2, bm::id_max/2+1 };
        bvect bv1(bv);
        BVectorInsert(&bv, bm::id_max/2+1, true);
        bv1.insert(bm::id_max/2+1, true);
        int cmp = bv1.compare(bv);
        assert(cmp==0);
        bv.optimize();
        bv1.optimize();
        BVectorInsert(&bv, bm::id_max/2+1, true);
        bv1.insert(bm::id_max/2+1, true);
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
        generate_bvector(bv, bm::id_max32+65536, true);
        bvect bv_control(bv);
        
        unsigned max_shifts = 100;
        for (unsigned i = 0; i < max_shifts; ++i)
        {
            bvect bv2(bm::BM_GAP);
            bv2 = bv_control;
            bv2.optimize();
            
            bvect::size_type i_pos = (unsigned)rand()%40000000;
            
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
            
            {
                cout << "\r" << i << "/" << max_shifts << flush;
            }
        } // for i
    }
    cout << "ok.\n";


    cout << "---------------------------- Bvector INSERT test OK" << endl;
}

// Reference bit erase
static
void BVectorErase(bvect*  bv, bvect::size_type pos)
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
            bvect::size_type v = *en;
            assert(v > pos);
            bi = v -1;
        }
    }
    bv->swap(bv_tmp);
}


static
void BvectorEraseTest()
{
    cout << "---------------------------- Bvector ERASE test" << endl;
    
    {
        bvect bv;
        bv.erase(bm::id_max/2);
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
        bv.set_range(bm::id_max32+65536, bm::id_max32+65536 + 65536);
        bvect bv_c(bv);

        auto cnt1 = bv.count();
        bv.optimize();
        bv.erase(0);
        BVectorErase(&bv_c, 0);
        
        auto cnt2 = bv.count();
        assert(cnt1 == cnt2);
        auto cnt3 = bv.count_range(bm::id_max32+65535, bm::id_max32+65535 + 65536);
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
        auto cnt1 = bv.count();
        bv.optimize();
        assert(cnt1 == bv.count());
        bv.erase(0);
        auto cnt2 = bv.count();
        assert(cnt1 == cnt2);
        auto cnt3 = bv.count_range(65535, 65535 + 65535);
        assert(cnt3 == cnt1);
    }
    
    {
        bvect bv;
        bv.invert();
        auto cnt1 = bv.count();
        bv.erase(bm::id_max32);
        auto cnt2 = bv.count();
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
            auto cnt = bv.count();
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
            bvect::size_type pos;
            bvect::size_type from = rd();
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


static
void AddressResolverTest()
{
    cout << "---------------------------- AddressResolverTest()" << endl;
    bvect::size_type id_to;
    bool found;

    {
        bm::bvps_addr_resolver<bvect>  ares;
        
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
        bm::bvps_addr_resolver<bvect>  ares;
        
        ares.set(1000);
        ares.set(10000);
        ares.set(bm::id_max-1);

        found = ares.resolve(10, &id_to);
        assert(!found);
        assert(id_to == 0);

        found = ares.resolve(bm::id_max-1, &id_to);
        assert(found);
        assert(id_to == 3);
        
        assert(ares.in_sync() == false);
        
        ares.optimize();
        assert(ares.in_sync() == false);

        ares.sync();
        assert(ares.in_sync());

        found = ares.resolve(bm::id_max-1, &id_to);
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
        bm::sv_addr_resolver<sparse_vector<bm::id_t, bvect> > ares;
        
        found = ares.resolve(10, &id_to);
        assert(!found);
        assert(id_to == 0);
    }

    {
        sv_addr_resolver<sparse_vector<bm::id64_t, bvect> > ares;
        
        ares.set(1000);   // 1
        ares.set(10000);  // 2
        ares.set(bm::id_max-1); // 3
        ares.set(5);      // 4
        
        found = ares.resolve(10, &id_to);
        assert(!found);
        assert(id_to == 0);
        
        found = ares.resolve(1000, &id_to);
        assert(found);
        assert(id_to == 1);

        found = ares.resolve(bm::id_max-1, &id_to);
        assert(found);
        assert(id_to == 3);
        
        ares.optimize();
        
        found = ares.resolve(5, &id_to);
        assert(found);
        assert(id_to == 4);
    }

    cout << "---------------------------- AddressResolverTest() OK" << endl;

}

static
void TestRankCompress()
{
    cout << " ------------------------------ Test Rank Compressor " << endl;
    
    int cmp;

    cout << "Step 1" << endl;
    {
        bvect bv1, bv2;
        bvect bv_s { 0, 1,        bm::id_max/2 };
        bvect bv_i { 0, 1, 2, 10, bm::id_max/2 };
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
            bv_i.build_rs_index(&bc);
        }
    }
    cout << "Step 1 - OK" << endl;

    {
        bvect bv1, bv2;
        bvect bv_s { 0, 100000, 100001,                  1600000, bm::id_max/2  };
        bvect bv_i { 0, 100000, 100001, 200000, 1000000, 1600000, bm::id_max/2 };
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
            bv_i.build_rs_index(&bc);

        }
    }
    std::cout << "basic test OK..." << std::endl;

    {
        std::cout << "\nStress rank compression...\n" << std::endl;
        bm::rank_compressor<bvect> rc;
        unsigned test_count = 2;
        bvect::size_type bv_size = bm::id_max32 * 2;
        for (unsigned i  = 0; i < test_count; ++i)
        {
//            if (bv_size > 40000000)
//                break;
            cout << "target size = " << bv_size << " " << endl;
            bvect bv_i, bv_s, bv_sr;
            generate_bvector(bv_i, bv_size, true);
            generate_bvector(bv_s, bv_size, true);
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
                bm::chrono_taker<std::ostream> ct(cout, "c1");
                rc.compress(bv1, bv_i, bv_s);
                rc.decompress(bv_sr, bv_i, bv1);
                }
                assert(bv1.count() == bv_s.count());
                cmp = bv_sr.compare(bv_s);
                assert(cmp == 0);

                {
                chrono_taker<std::ostream> ct(cout, "c2");
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
                    chrono_taker<std::ostream> ct(cout, "c1-1");
                    rc.compress(bv1, bv_i, bv_subset);
                    rc.decompress(bv_sr, bv_i, bv1);
                    }
                    assert(bv1.count() == bv_subset.count());
                    cmp = bv_sr.compare(bv_subset);
                    assert(cmp == 0);

                    {
                    chrono_taker<std::ostream> ct(cout, "c2-2");
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
                bv_i.build_rs_index(&bc);

            } // for j
            cout << "\n" << i << " of " << test_count << "  " << endl;
            
            bv_size += bv_size;
        } // for i
        std::cout << endl << "Stress rank compression... OK" << std::endl;
    }
    
    
    cout << " ------------------------------ Test Rank Compressor OK " << endl;
}


// do logical operation through serialization
static
bvect::size_type SerializationOperation(bvect*     bv_target,
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
    if (!bv1.is_ro())
        bv1.optimize(tb, bvect::opt_compress, st1_op);
    else
        bv1.calc_stat(st1_op);

    if (!bv2.is_ro())
        bv2.optimize(tb, bvect::opt_compress, st2_op);
    else
        bv2.calc_stat(st2_op);


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

   bm::serializer<bvect64> bv_ser;

    if (rand() & 1) // setup random bookmark set
    {
        unsigned bm_range = (unsigned)rand()%256;
        bv_ser.set_bookmarks(true, bm_range);
        cout << "Bookmark ON at every:" << bm_range << endl;
    }

   size_t slen1 = bv_ser.serialize(bv1, smem1, st1.max_serialize_mem);
   size_t slen2 = bm::serialize(bv2, smem2, tb);

   if (slen1 > st1.max_serialize_mem || slen2 > st2.max_serialize_mem)
   {
       cout << "Serialization override detected!" << endl;
       assert(0);
       exit(1);
   }

   operation_deserializer<bvect64> od;

    // Range extraction test
    {
        bvect::size_type first, last;
        bool found = bv1.find_range(first, last);
        if (found)
        {
            unsigned frac = (unsigned)rand() %4;
            if (!frac)
                frac = 4;
            auto r_part = (last - first) / frac;
            bvect::size_type start, end;
            if (r_part > last)
                start = r_part;
            else
                start = last - r_part;
            end = last;

            bvect bv_r;
            bv_r.copy_range(bv1, start, end);
            bvect bv_od_r;
            od.deserialize_range(bv_od_r, smem1, start, end);

            bool eq;
            eq = bv_r.equal(bv_od_r);
            assert(eq);
        }
    }

   bvect::size_type count =
       od.deserialize(*bv_target,
                      smem1,
                      0,
                      set_ASSIGN);
   //cout << slen1 << " " << slen2 << endl;
   int res = bv1.compare(*bv_target);
   if (res != 0)
   {
       cout << "---------------------------------- " << endl;
       cout << "bv1.count()=" << bv1.count() << endl;
       print_stat(cout,bv1);
       cout << "---------------------------------- " << endl;
       cout << "bv_target.count()=" << bv_target->count() << endl;
       print_stat(cout,*bv_target);
       
       bv_target->bit_xor(bv1);
       cout << "First diff=" << bv_target->get_first() << endl;
       cout << "set_ASSIGN 1 failed!" << endl;
       assert(0); exit (1);
   }

   {
       bvect* bv_tmp2 = new bvect();
       bm::deserialize(*bv_tmp2, smem1);
       if (*bv_tmp2 != bv1)
       {
           cout << "Deserialize NOT equal to Operation deserialize!" << endl;
           assert(0);exit(1);
       }
       delete bv_tmp2;
   }

   cout << "Operation deserialization... " << op << endl;

    count=
       od.deserialize(*bv_target,
                      smem2,
                      0,
                      op);

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

        bvect bvt(bv1, bm::finalization::READWRITE);
        switch(op)
        {
        case bm::set_OR:
            {
                bvect bvc(bv1, bm::finalization::READWRITE);
                bvc |= bv2;
                bvect bv_merge1(bv1, bm::finalization::READWRITE);
                bvect bv_merge2(bv2, bm::finalization::READWRITE);
                bv_merge1.merge(bv_merge2);
                
                if (bv_merge1 != bvc)
                {
                    cerr << "Merge(OR) check error!" << endl;
                    assert(0);
                    exit(1);
                }
                // 2-way
                {
                    bvect bvt1;
                    bvt1.bit_or(bv1, bv2, bvect::opt_none);
                    if (bvt1 != bvc)
                    {
                        cerr << "1. OR 2-way check error!" << endl;
                        assert(0);
                        exit(1);
                    }
                    bvect bvt2;
                    bvt2.bit_or(bv2, bv1, bvect::opt_compress);
                    if (bvt2 != bvc)
                    {
                        cerr << "2. OR 2-way check error!" << endl;
                        assert(0);exit(1);
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
                bvect bvc(bv1, bm::finalization::READWRITE);
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
                    assert(0);
                    exit(1);
                }
                bvect bvt2;
                bvt2.bit_xor(bv2, bv1, bvect::opt_compress);
                if (bvt2 != bvc)
                {
                    cerr << "2. XOR 2-way check error!" << endl;
                    assert(0);
                    exit(1);
                }
            }
            break;
        case bm::set_AND:
            bvt &= bv2;
            agg.combine_and(bv_agg, agg_list, 2);
            agg.combine_and_sub(bv_agg2, agg_list, 2, 0, 0, false);
            {
                if (bv_agg.compare(bv_agg2) != 0)
                {
                    cerr << "Error: Aggregator AND - AND-SUB(0) comparison failed!" << endl;
                    assert(0); exit(1);
                }
            }
            agg_check = true;
            
            // 2-way
            {
                bvect bvc(bv1, bm::finalization::READWRITE);
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
                    assert(0); exit(1);
                }
            }
            break;
        case bm::set_SUB:
            bvt -= bv2;
            agg.combine_and_sub(bv_agg, agg_list, 1, agg_list2, 1, false);
            {
                bvect bv_h;
                agg.combine_and_sub_horizontal(bv_h, agg_list, 1, agg_list2, 1);
                if (bv_agg.compare(bv_h) != 0)
                {
                    cerr << "Error: Aggregator Horz-AND-SUB comparison failed!" << endl;
                    exit(1);
                }
            }
            agg_check = true;
            // 2-way
            {
                bvect bvc1(bv1, bm::finalization::READWRITE);
                bvect bvc2(bv2, bm::finalization::READWRITE);
                bvc1 -= bv2;
                bvc2 -= bv1;
                
                bvect bvt1;
                bvt1.bit_sub(bv1, bv2, bvect::opt_compress);
                if (bvt1 != bvc1)
                {
                    DetailedCompareBVectors(bvt1, bvc1);
                    cerr << "1. SUB 2-way check error!" << endl;
                    assert(0);exit(1);
                }
                bvect bvt2;
                bvt2.bit_sub(bv2, bv1, bvect::opt_compress);
                if (bvt2 != bvc2)
                {
                    cerr << "2. SUB 2-way check error!" << endl;
                    assert(0);exit(1);
                }
            }

            break;
        default:
            goto no_compare;
        }
        if (bvt.compare(*bv_target) != 0)
        {
            cout << "Direct Serial operation comparison failed!" << endl;
            assert(0);exit(1);
        }
        if (agg_check && bvt.compare(bv_agg) != 0)
        {
            cerr << "Error: Aggregator operation comparison failed!" << endl;
            assert(0);exit(1);
        }

        no_compare:
        ;

    }

   if (check_reverse)
   {
        cout << "Reverse check... " << endl;
        bvect bv_tmp2(BM_GAP);
        od.deserialize(bv_tmp2,
                       smem2,
                       0,
                       set_ASSIGN);
        res = bv_tmp2.compare(bv2);
        if (res != 0)
        {
            cout << "set_ASSIGN failed 2! " << endl;
            assert(0);exit(1);
        }
        cout << "Deserialization assign to bv_tmp2 OK" << endl;
        auto count_rev =
        od.deserialize(bv_tmp2,
                       smem1,
                       0,
                       op);
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
            assert(0);
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
                                 bvect::size_type  predicted_count,
                                 set_operation op_count,
                                 set_operation op_combine)
{
    bv_target->clear(true);
    cout << "Serialization operation count..." << endl;

    bvect::size_type scount1;

    bvect bv_ro1(bv1, bm::finalization::READONLY);
    bvect bv_ro2(bv2, bm::finalization::READONLY);
    bvect bv_target2(*bv_target);


    scount1 = SerializationOperation(0,
                                      bv1,
                                      bv2,
                                      op_count,
                                      true //reverse check
                                    );
    cout << "Serialization operation count OK." << endl;


    cout << "Serialization operation. " << endl;
    auto scount2 = SerializationOperation(bv_target,
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
        cout << endl << endl << "Reference" << endl;
        if (op_combine == set_OR)
        {
            bv1 |= bv2;
            if (bv1 != *bv_target)
            {
                cout << "Comparison OR error!" << endl;
            }
            cout << "OR operation count=" << bv1.count() << endl;
        }
        else
            if (op_combine == set_AND)
            {
                bv1 &= bv2;
            }
        assert(0);
        exit(1);
    }
    auto scount3 = SerializationOperation(&bv_target2,
                                          bv_ro1,
                                          bv_ro2,
                                          op_combine);
    scount3 = bv_target2.count();
    assert(scount3 == scount2);
    bool eq = bv_target->equal(bv_target2);
    assert(eq);

    cout << "OK" << endl;
}




static
void AndOperationsTest()
{
    const bvect::size_type BITVECT_SIZE = bm::id_max32 * 2;
    assert(ITERATIONS < BITVECT_SIZE);

    cout << "----------------------------------- AndOperationTest" << endl;

    {
        ref_vect vect;
        generate_vect_simpl0(vect);

        bvect bv0;
        load_BV_set_ref(cout, bv0, vect);
        bvect bv1(bv0);
        bvect bv11(bv0);

        {
            bvect bvt;
            SerializationOperation2Test(&bvt,
                                        bv0,
                                        bv1,
                                        vect.size(),
                                        set_COUNT_AND,
                                        set_AND);
        }
        
        bvect bv_i;
        bv_i.invert();
        bvect bv_ro_i(bv_i, bm::finalization::READONLY);

        for (unsigned i = 0; i < 2; ++i)
        {
            bvect::size_type predicted_count = bm::count_and(bv0, bv1);
            assert(predicted_count == vect.size());
            auto predicted_any = bm::any_and(bv1, bv0);
            if (predicted_any == 0 && predicted_count != 0)
            {
                cout << "Predicted any error!" << endl;
                assert(0);
                exit(1);
            }
            predicted_count = bm::count_and(bv0, bv_i);
            assert(predicted_count == vect.size());

            bv1.bit_and(bv0);
            auto count = bv1.count();
            if (count != predicted_count)
            {
                cout << "Predicted count error!" << endl;
                assert(0);
                exit(1);
            }
            compare_BV_set_ref(bv1, vect);

            
            bv1.bit_and(bv_i);
            compare_BV_set_ref(bv1, vect);

            {
                bvect bv_ro0(bv0, bm::finalization::READONLY);
                auto pcount1 = bm::count_and(bv_ro0, bv_i);
                assert(pcount1 == predicted_count);

                bv11.bit_and(bv_ro0);
                bv11.bit_and(bv_ro_i);
                compare_BV_set_ref(bv11, vect);
            }

            {
                bvect bv_i_c;
                bv_i_c.invert();
                predicted_count = bm::count_and(bv0, bv_i_c);
                assert(predicted_count == vect.size());
                bv_i_c.bit_and(bv0);
                int cmp = bv_i_c.compare(bv0);
                assert(cmp == 0);
                compare_BV_set_ref(bv_i_c, vect);
            }
            
            bv0.optimize();
            bv1.optimize();
        } // for
    }

    {
        cout << "\n48-bit large set intersect test" << endl;
        cout << "generation..." << endl;
        ref_vect vect0;
        generate_vect48(vect0);
        bvect bv0;
        load_BV_set_ref(cout, bv0, vect0);
        compare_BV_set_ref(bv0, vect0);

        ref_vect vect1;
        generate_vect48(vect1);
        bvect bv1;
        load_BV_set_ref(cout, bv1, vect1);
        compare_BV_set_ref(bv1, vect1);
        cout << "ok\n" << endl;

        ref_vect vect_i;
        std::set_intersection(vect0.begin(), vect0.end(),
                              vect1.begin(), vect1.end(),
                              std::back_inserter(vect_i));
        {
            bvect bv_i(bv0);
            bv_i.bit_and(bv1);
            compare_BV_set_ref(bv_i, vect_i);
        }
        bv0.optimize();
        print_bvector_stat(cout,bv0);
        {
            bvect bv_i(bv0);
            bv_i.bit_and(bv1);
            compare_BV_set_ref(bv_i, vect_i);
        }
        bv1.optimize();
        print_bvector_stat(cout,bv1);
        {
            bvect bv_i(bv0);
            bv_i.bit_and(bv1);
            compare_BV_set_ref(bv_i, vect_i);
        }
        bvect bvt;
        SerializationOperation2Test(&bvt,
                                    bv0,
                                    bv1,
                                    vect_i.size(),
                                    set_COUNT_AND,
                                    set_AND);
    }

    {
        bvect_mini   bvect_min1(BITVECT_SIZE);
        bvect_mini   bvect_min2(BITVECT_SIZE);
        bvect        bvect_full1;
        bvect        bvect_full2;

        bvect_full1.set_new_blocks_strat(bm::BM_GAP);
        bvect_full2.set_new_blocks_strat(bm::BM_GAP);



        printf("AND test\n");

        bvect_min1.set_bit(1);
        bvect_min1.set_bit(12);
        bvect_min1.set_bit(bm::id_max32 + 256);

        bvect_min2.set_bit(12);
        bvect_min2.set_bit(bm::id_max32 + 256);

        bvect_min1.combine_and(bvect_min2);

        bvect_full1.set_bit(1);
        bvect_full1.set_bit(12);
        bvect_full1.set_bit(bm::id_max32 + 256);

        bvect_full2.set_bit(12);
        bvect_full2.set_bit(bm::id_max32 + 256);

        auto predicted_count = bm::count_and(bvect_full1, bvect_full2);

        auto predicted_any = bm::any_and(bvect_full1, bvect_full2);
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

        bvect::size_type count = bvect_full1.count();
        if (count != predicted_count)
        {
            cout << "Predicted count error!" << endl;
            exit(1);
        }

        CheckVectors(bvect_min1, bvect_full1, 256, true);
        CheckVectors(bvect_min1, bv_target_s, 256, true);
//        CheckCountRange(bvect_full1, 0, 256);

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

    CheckVectors(bvect_min1, bvect_full1, BITVECT_SIZE/10+10, true);
//    CheckCountRange(bvect_full1, 0, BITVECT_SIZE/10+10);

//    FillSets(&bvect_min1, &bvect_full1, 1, BITVECT_SIZE/7, 0);
//    FillSets(&bvect_min2, &bvect_full2, 1, BITVECT_SIZE/7, 0);

    bvect_min1.combine_and(bvect_min2);

    auto predicted_count = bm::count_and(bvect_full1,bvect_full2);
    auto predicted_any = bm::any_and(bvect_full1, bvect_full2);
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

    auto count = bvect_full1.count();
    if (count != predicted_count)
    {
        cout << "Predicted count error!" << endl;
        exit(1);
    }

    CheckVectors(bvect_min1, bvect_full1, BITVECT_SIZE/10+10, true);
    CheckVectors(bvect_min1, bv_target_s, BITVECT_SIZE/10+10, true);
//    CheckCountRange(bvect_full1, 0, BITVECT_SIZE/10+10);

    }


    {

    bvect_mini   bvect_min1(BITVECT_SIZE);
    bvect_mini   bvect_min2(BITVECT_SIZE);
    bvect        bvect_full1;
    bvect        bvect_full2;

    bvect_full1.set_new_blocks_strat(bm::BM_GAP);
    bvect_full2.set_new_blocks_strat(bm::BM_GAP);

    printf("AND test stage 2.\n");


    FillSets(&bvect_min1, &bvect_full1, BITVECT_SIZE/2, BITVECT_SIZE, 0ull);
    FillSets(&bvect_min2, &bvect_full2, BITVECT_SIZE/2, BITVECT_SIZE, 0ull);

    auto predicted_count = bm::count_and(bvect_full1,bvect_full2);
    auto predicted_any = bm::any_and(bvect_full1, bvect_full2);
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

    auto count = bvect_full1.count();
    if (count != predicted_count)
    {
        cout << "Predicted count error!" << endl;
        print_stat(cout,bvect_full1);
        exit(1);
    }

    CheckVectors(bvect_min1, bvect_full1, BITVECT_SIZE/10+10, true);
    CheckVectors(bvect_min1, bv_target_s, BITVECT_SIZE/10+10, true);
//    CheckCountRange(bvect_full1, 0, BITVECT_SIZE/10+10);

    VisitorAllRangeTest(bv_target_s, 0);

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


    FillSets(&bvect_min1, &bvect_full1, 1ull, BITVECT_SIZE/5, 2ull);
    FillSets(&bvect_min2, &bvect_full2, 1ull, BITVECT_SIZE/5, 2ull);

    bvect_min1.combine_and(bvect_min2);

    auto predicted_count = bm::count_and(bvect_full1, bvect_full2);
    auto predicted_any = bm::any_and(bvect_full1, bvect_full2);
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

    auto count = bvect_full1.count();
    if (count != predicted_count)
    {
        cout << "Predicted count error!" << endl;
        exit(1);
    }

    CheckVectors(bvect_min1, bvect_full1, BITVECT_SIZE, true);
    CheckVectors(bvect_min1, bv_target_s, BITVECT_SIZE, true);

    BM_DECLARE_TEMP_BLOCK(tb)
    bvect_full1.optimize(tb);
    CheckVectors(bvect_min1, bvect_full1, BITVECT_SIZE, true);
    }

    printf("AND test stage 4. combine_and_sorted\n");
    {
        bvect::size_type ids[] = {0, 1, 2, 3, 10, 65535, 65536, 65535*2, 65535*3};
        size_t to_add = sizeof(ids)/sizeof(bvect::size_type);
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
     
        bvect::size_type* first = ids;
        bvect::size_type* last = ids + to_add;
     
        bvect_min1.combine_and(bvect_min2);

        bm::combine_and_sorted(bvect_full1, first, last);
        CheckVectors(bvect_min1, bvect_full1, BITVECT_SIZE, true);
    }
 
    {
        bvect        bvect1 { 1, 10, 12, bm::id_max32 + 256};
        bvect        bvect2 { 2, 15, 165535, bm::id_max32 + 250 };
     
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
void OrOperationsTest()
{
    const bvect::size_type BITVECT_SIZE = bm::id_max32 * 2;
    assert(ITERATIONS < BITVECT_SIZE);

    cout << "----------------------------------- OrOperationTest" << endl;
    
    {
        ref_vect vect;
        generate_vect_simpl0(vect);

        bvect bv0;
        load_BV_set_ref(cout, bv0, vect);
        bvect bv1(bv0);

        {
            bvect bvt;
            SerializationOperation2Test(&bvt,
                                        bv0,
                                        bv1,
                                        vect.size(),
                                        set_COUNT_OR,
                                        set_OR);
        }

        bvect bv_i;
        bv_i.invert();
        auto full_cnt = bm::id_max;
        
        for (unsigned i = 0; i < 2; ++i)
        {
            cout << "Pass " << i << endl;
            bvect::size_type predicted_count = bm::count_or(bv0, bv1);
            assert(predicted_count == vect.size());
            auto predicted_any = bm::any_or(bv1, bv0);
            if (predicted_any == 0 && predicted_count != 0)
            {
                cout << "Predicted any error!" << endl;
                assert(0);
                exit(1);
            }

            predicted_count = bm::count_or(bv0, bv_i);
            assert(predicted_count == full_cnt);

            bv1.bit_or(bv_i);
            {
                int cmp = bv1.compare(bv_i);
                assert(cmp == 0);
            }
            bv1 = bv0;

            {
                bvect bv_i_c;
                bv_i_c.invert();
                predicted_count = bm::count_or(bv0, bv_i_c);
                assert(predicted_count == full_cnt);
                bv_i_c.bit_or(bv0);
                int cmp = bv_i_c.compare(bv_i);
                assert(cmp == 0);
            }

            bv0.optimize();
            bv1.optimize();
        } // for
    }
    
    {
        cout << "\n48-bit large set union test" << endl;
        cout << "generation..." << endl;
        ref_vect vect0;
        generate_vect48(vect0);
        bvect bv0;
        load_BV_set_ref(cout, bv0, vect0);
        compare_BV_set_ref(bv0, vect0);

        ref_vect vect1;
        generate_vect48(vect1);
        bvect bv1;
        load_BV_set_ref(cout, bv1, vect1);
        compare_BV_set_ref(bv1, vect1);
        cout << "ok\n" << endl;

        ref_vect vect_i;
        std::set_union(vect0.begin(), vect0.end(),
                       vect1.begin(), vect1.end(),
                       std::back_inserter(vect_i));
        {
            bvect bv_i(bv0);
            bv_i.bit_or(bv1);
            compare_BV_set_ref(bv_i, vect_i);
        }
        bv0.optimize();
        print_bvector_stat(cout,bv0);
        {
            bvect bv_i(bv0);
            bv_i.bit_or(bv1);
            compare_BV_set_ref(bv_i, vect_i);
        }
        bv1.optimize();
        print_bvector_stat(cout,bv1);
        {
            bvect bv_i(bv0);
            bv_i.bit_or(bv1);
            compare_BV_set_ref(bv_i, vect_i);
        }
        bvect bvt;
        SerializationOperation2Test(&bvt,
                                    bv0,
                                    bv1,
                                    vect_i.size(),
                                    set_COUNT_OR,
                                    set_OR);
    }

    
    {
        bvect_mini   bvect_min1(BITVECT_SIZE);
        bvect_mini   bvect_min2(BITVECT_SIZE);
        bvect        bvect_full1;
        bvect        bvect_full2;

        bvect_full1.set_new_blocks_strat(bm::BM_GAP);
        bvect_full2.set_new_blocks_strat(bm::BM_GAP);

        printf("OR test\n");

        bvect_min1.set_bit(1);
        bvect_min1.set_bit(12);
        bvect_min1.set_bit(bm::id_max32 + 256);

        bvect_min2.set_bit(12);
        bvect_min2.set_bit(bm::id_max32 + 256);

        bvect_min1.combine_or(bvect_min2);

        bvect_full1.set_bit(1);
        bvect_full1.set_bit(12);
        bvect_full1.set_bit(bm::id_max32 + 256);

        bvect_full2.set_bit(12);
        bvect_full2.set_bit(bm::id_max32 + 256);
        
        bvect::size_type predicted_count = bm::count_or(bvect_full1, bvect_full2);
        auto predicted_any = bm::any_or(bvect_full1, bvect_full2);
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

        auto count = bvect_full1.count();
        if (count != predicted_count)
        {
            cout << "Predicted count error!" << endl;
            cout << predicted_count << " " << count << endl;
            print_stat(cout,bvect_full1);
            exit(1);
        }


        CheckVectors(bvect_min1, bvect_full1, 256, true);
        CheckVectors(bvect_min1, bv_target_s, 256, true);
//        CheckCountRange(bvect_full1, 0, 256);
//        CheckCountRange(bvect_full1, 128, 256);
    }

    {
        bvect_mini   bvect_min1(BITVECT_SIZE);
        bvect_mini   bvect_min2(BITVECT_SIZE);
        bvect        bvect_full1;
        bvect        bvect_full2;

        bvect_full1.set_new_blocks_strat(bm::BM_GAP);
        bvect_full2.set_new_blocks_strat(bm::BM_GAP);

        printf("OR test stage 2.\n");


        FillSets(&bvect_min1, &bvect_full1, 1ull, BITVECT_SIZE/7, 0ull);
        FillSets(&bvect_min2, &bvect_full2, 1ull, BITVECT_SIZE/7, 0ull);

        bvect_min1.combine_or(bvect_min2);

        auto predicted_count = bm::count_or(bvect_full1, bvect_full2);
        bvect::size_type predicted_any = bm::any_or(bvect_full1, bvect_full2);
        if (predicted_any == 0 && predicted_count != 0)
        {
            cout << "Predicted any error!" << endl;
            assert(0); exit(1);
        }

        bvect    bv_target_s;
        SerializationOperation2Test(&bv_target_s,
                                    bvect_full1,
                                    bvect_full2,
                                    predicted_count,
                                    set_COUNT_OR,
                                    set_OR);


        bvect_full1.bit_or(bvect_full2);

        auto count = bvect_full1.count();
        if (count != predicted_count)
        {
            cout << "Predicted count error!" << endl;
            exit(1);
        }

        CheckVectors(bvect_min1, bvect_full1, BITVECT_SIZE/10+10, true);
        CheckVectors(bvect_min1, bvect_full1, BITVECT_SIZE/10+10, true);
//        CheckCountRange(bvect_full1, 0, BITVECT_SIZE/10+10);
    }
    
    {
        bvect bv1;
        bvect bv2;
        bv1.flip(); bv2.flip();
        auto cnt1 = bv1.count();
        bv1.bit_or(bv2);
        bvect::size_type cnt2 = bv1.count();
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


        FillSets(&bvect_min1, &bvect_full1, 1ull, BITVECT_SIZE/5, 2ull);
        FillSets(&bvect_min2, &bvect_full2, 1ull, BITVECT_SIZE/5, 2ull);

        bvect_min1.combine_or(bvect_min2);
        auto mcnt = bvect_min1.bit_count();

        cout << mcnt << endl;
        
        auto predicted_count = bm::count_or(bvect_full1, bvect_full2);
        cout << predicted_count << endl;
        auto predicted_any = bm::any_or(bvect_full1, bvect_full2);
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

        auto count = bvect_full1.count();
        if (count != predicted_count)
        {
            cout << "Predicted count error!" << endl;
            exit(1);
        }

        CheckVectors(bvect_min1, bvect_full1, BITVECT_SIZE, true);

        BM_DECLARE_TEMP_BLOCK(tb)
        bvect_full1.optimize(tb);

        CheckVectors(bvect_min1, bvect_full1, BITVECT_SIZE, true);
        CheckVectors(bvect_min1, bv_target_s, BITVECT_SIZE, true);
//        CheckCountRange(bvect_full1, 0, BITVECT_SIZE);
    }
    
    cout << "Testing combine_or" << endl;
    
    {
        bvect        bvect_full1;
        bvect        bvect_full2;
        bvect_mini   bvect_min1(BITVECT_SIZE);
        
        bvect_full1.set_new_blocks_strat(bm::BM_GAP);
        bvect_full2.set_new_blocks_strat(bm::BM_GAP);

        bvect::size_type ids[10000];
        unsigned to_add = 10000;
        
        unsigned bn = bm::id_max32;
        for (unsigned i = 0; i < to_add; ++i)
        {
            ids[i] = bn;
            bvect_full2.set(bn);
            bvect_min1.set_bit(bn);
            bn += 15;
        }
        
        bvect::size_type* first = ids;
        bvect::size_type* last = ids + to_add;
        
        bm::combine_or(bvect_full1, first, last);

        CheckVectors(bvect_min1, bvect_full1, BITVECT_SIZE, true);
        
        bm::combine_or(bvect_full1, first, last);
        CheckVectors(bvect_min1, bvect_full1, BITVECT_SIZE, true);
    }
    
    
    {
        bvect::size_type ids[] = {0, 65536, 65535, 65535*3, 65535*2, 10};
        size_t to_add = sizeof(ids)/sizeof(bvect::size_type);
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
        
        bvect::size_type* first = ids;
        bvect::size_type* last = ids + to_add;
        
        bm::combine_or(bvect_full1, first, last);
        CheckVectors(bvect_min1, bvect_full1, BITVECT_SIZE, true);

        bm::combine_or(bvect_full1, first, last);
        CheckVectors(bvect_min1, bvect_full1, BITVECT_SIZE, true);
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
void XorOperationsTest()
{
    const bvect::size_type BITVECT_SIZE = bm::id_max32 * 2;
    assert(ITERATIONS < BITVECT_SIZE);

    cout << "----------------------------------- XorOperationTest" << endl;
    

    {
        bvect bv_i0;
        bv_i0.invert();
        bvect bv_i1;
        bv_i1.invert();
        auto cnt = bm::count_xor(bv_i0, bv_i1);
        assert(!cnt);
        auto predicted_any = bm::any_xor(bv_i0, bv_i1);
        assert(!predicted_any);
    }
    
    {
        ref_vect vect;
        generate_vect_simpl0(vect);

        bvect bv0;
        load_BV_set_ref(cout,bv0, vect);
        bvect bv1(bv0);
        
        {
            bvect bv3;
            bv3.bit_xor(bv0, bv1, bvect::opt_compress);
            assert(!bv3.count());
            assert(!bv3.any());
            
            bvect bv4;
            int cmp = bv4.compare(bv3);
            assert(cmp == 0);
            
            bv3 = bv0;
            bv3 ^= bv1;
            assert(!bv3.count());
            assert(!bv3.any());
            cmp = bv4.compare(bv3);
            assert(cmp == 0);
        }
        
        bvect bv_i;
        bv_i.invert();
        auto full_cnt = bm::id_max - vect.size();
        
        for (unsigned i = 0; i < 2; ++i)
        {
            bvect::size_type predicted_count = bm::count_xor(bv0, bv1);
            assert(predicted_count == 0);
            auto predicted_any = bm::any_xor(bv1, bv0);
            assert(!predicted_any);

            predicted_count = bm::count_xor(bv0, bv_i);
            assert(predicted_count == full_cnt);
            
            {
                bvect bv3(bv_i);
                clear_BV_set_ref(cout,bv3, vect);
                predicted_count = bm::count_xor(bv0, bv3);
                assert(predicted_count == bm::id_max);
                auto cnt = bv3.count();
                assert(cnt == bm::id_max - vect.size());
                bv3 ^= bv0;
                cnt = bv3.count();
                assert(cnt == bm::id_max);
            }
            
            
            bv0.optimize();
            bv1.optimize();
        } // for
        
    }

    {
        cout << "\n48-bit large set difference test" << endl;
        cout << "generation..." << endl;
        ref_vect vect0;
        generate_vect48(vect0);
        bvect bv0;
        load_BV_set_ref(cout,bv0, vect0);
        compare_BV_set_ref(bv0, vect0);

        ref_vect vect1;
        generate_vect48(vect1);
        bvect bv1;
        load_BV_set_ref(cout,bv1, vect1);
        compare_BV_set_ref(bv1, vect1);
        cout << "ok\n" << endl;

        ref_vect vect_i;
        std::set_symmetric_difference(vect0.begin(), vect0.end(),
                            vect1.begin(), vect1.end(),
                            std::back_inserter(vect_i));
        {
            bvect bv_i(bv0);
            bv_i.bit_xor(bv1);
            compare_BV_set_ref(bv_i, vect_i);
        }
        bv0.optimize();
        print_bvector_stat(cout,bv0);
        {
            bvect bv_i(bv0);
            bv_i.bit_xor(bv1);
            compare_BV_set_ref(bv_i, vect_i);
        }
        bv1.optimize();
        print_bvector_stat(cout,bv1);
        {
            bvect bv_i(bv0);
            bv_i.bit_xor(bv1);
            compare_BV_set_ref(bv_i, vect_i);
        }
        bvect bvt;
        SerializationOperation2Test(&bvt,
                                    bv0,
                                    bv1,
                                    vect_i.size(),
                                    set_COUNT_XOR,
                                    set_XOR);
    }

    

    {
        bvect_mini   bvect_min1(BITVECT_SIZE);
        bvect_mini   bvect_min2(BITVECT_SIZE);
        bvect        bvect_full1;
        bvect        bvect_full2;

        bvect_full1.set_new_blocks_strat(bm::BM_GAP);
        bvect_full2.set_new_blocks_strat(bm::BM_GAP);



        printf("XOR test\n");

        bvect_min1.set_bit(1);
        bvect_min1.set_bit(12);
        bvect_min1.set_bit(bm::id_max32 + 256);

        bvect_min2.set_bit(12);
        bvect_min2.set_bit(bm::id_max32 + 256);

        bvect_min1.combine_xor(bvect_min2);

        bvect_full1.set_bit(1);
        bvect_full1.set_bit(12);
        bvect_full1.set_bit(bm::id_max32 + 256);

        bvect_full2.set_bit(12);
        bvect_full2.set_bit(bm::id_max32 + 256);

        bvect::size_type predicted_count = bm::count_xor(bvect_full1, bvect_full2);
        auto predicted_any = bm::any_xor(bvect_full1, bvect_full2);
        if (predicted_any == 0 && predicted_count != 0)
        {
            cout << "Predicted any error!" << endl;
            assert(0); exit(1);
        }

        bvect    bv_target_s;
        SerializationOperation2Test(&bv_target_s,
                                    bvect_full1,
                                    bvect_full2,
                                    predicted_count,
                                    set_COUNT_XOR,
                                    set_XOR);


        bvect_full1.bit_xor(bvect_full2);

        bvect::size_type count = bvect_full1.count();
        if (count != predicted_count)
        {
            cout << "1.Predicted count error!" << endl;
            exit(1);
        }

        CheckVectors(bvect_min1, bvect_full1, 256, true);
        CheckVectors(bvect_min1, bv_target_s, 256, true);
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

        bvect::size_type predicted_count = bm::count_xor(bvect1, bvect2);
        auto predicted_any = bm::any_xor(bvect1, bvect2);
        if (predicted_any == 0 && predicted_count != 0)
        {
            cout << "Predicted any error!" << endl;
            assert(0); exit(1);
        }

        bvect    bv_target_s;
        SerializationOperation2Test(&bv_target_s,
                                    bvect1,
                                    bvect2,
                                    predicted_count,
                                    set_COUNT_XOR,
                                    set_XOR);

        bvect1.bit_xor(bvect2);
        
        bvect::size_type count = bvect1.count();
        if (count != predicted_count)
        {
            cout << "2.Predicted count error!" << endl;
            exit(1);
        }
        
        bvect_min1.combine_xor(bvect_min2);
        CheckVectors(bvect_min1, bvect1, BITVECT_SIZE, true);
        CheckVectors(bvect_min1, bv_target_s, BITVECT_SIZE, true);
    }
    
    {
        bvect bv1;
        bvect bv2;
        bv1.flip();
        bv2.flip();
        bv1.bit_xor(bv2);
        auto cnt2 = bv1.count();
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
        
        bvect::size_type predicted_count = bm::count_xor(bvect1, bvect2);
        auto predicted_any = bm::any_xor(bvect1, bvect2);
        if (predicted_any == 0 && predicted_count != 0)
        {
            cout << "Predicted any error!" << endl;
            assert(0); exit(1);
        }

        bvect    bv_target_s;
        SerializationOperation2Test(&bv_target_s,
                                    bvect1,
                                    bvect2,
                                    predicted_count,
                                    set_COUNT_XOR,
                                    set_XOR);

        bvect1.bit_xor(bvect2);

        auto count = bvect1.count();
        if (count != predicted_count)
        {
            cout << "3.Predicted count error!" << endl;
            exit(1);
        }
        
        bvect_min1.combine_xor(bvect_min2);
        CheckVectors(bvect_min1, bvect1, BITVECT_SIZE, true);
        CheckVectors(bvect_min1, bv_target_s, BITVECT_SIZE, true);
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
        
        bvect::size_type predicted_count = bm::count_xor(bvect1, bvect2);
        auto predicted_any = bm::any_xor(bvect1, bvect2);
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

        bvect::size_type count = bvect1.count();
        if (count != predicted_count)
        {
            cout << "4.Predicted count error!" << endl;
            cout << count << " " << predicted_count << endl;
            assert(0);
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

        FillSets(&bvect_min1, &bvect_full1, 1ull, BITVECT_SIZE/7, 0ull);
        FillSets(&bvect_min2, &bvect_full2, 1ull, BITVECT_SIZE/7, 0ull);

        bvect_min1.combine_xor(bvect_min2);
        
        bvect::size_type predicted_count = bm::count_xor(bvect_full1, bvect_full2);
        auto predicted_any = bm::any_xor(bvect_full1, bvect_full2);
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
        
        bvect::size_type count = bvect_full1.count();
        if (count != predicted_count)
        {
            cout << "5.Predicted count error!" << endl;
            cout << count << " " << predicted_count << endl;
            print_stat(cout,bvect_full1);
            exit(1);
        }

        CheckVectors(bvect_min1, bvect_full1, BITVECT_SIZE/10+10, true);
        CheckVectors(bvect_min1, bv_target_s, BITVECT_SIZE/10+10, true);

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


        FillSets(&bvect_min1, &bvect_full1, 1ull, BITVECT_SIZE, 2ull);
        FillSets(&bvect_min2, &bvect_full2, 1ull, BITVECT_SIZE, 2ull);

        bvect::size_type predicted_count = bm::count_xor(bvect_full1, bvect_full2);
        auto predicted_any = bm::any_xor(bvect_full1, bvect_full2);
        if (predicted_any == 0 && predicted_count != 0)
        {
            cout << "Predicted any error!" << endl;
            assert(0);
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

        bvect::size_type count = bvect_full1.count();
        if (count != predicted_count)
        {
            cout << "6.Predicted count error!" << endl;
            exit(1);
        }

        CheckVectors(bvect_min1, bvect_full1, BITVECT_SIZE, true);

        BM_DECLARE_TEMP_BLOCK(tb)
        bvect_full1.optimize(tb);
        CheckVectors(bvect_min1, bvect_full1, BITVECT_SIZE, true);
        CheckVectors(bvect_min1, bv_target_s, BITVECT_SIZE, true);
    }


    cout << "Testing combine_xor" << endl;
    
    {
        bvect        bvect_full1;
        bvect        bvect_full2;
        bvect_mini   bvect_min1(BITVECT_SIZE);
        
        bvect_full1.set_new_blocks_strat(bm::BM_GAP);
        bvect_full2.set_new_blocks_strat(bm::BM_GAP);

        bvect::size_type ids[10000];
        unsigned to_add = 10000;
        
        bvect::size_type bn = bm::id_max32;
        for (unsigned i = 0; i < to_add; ++i)
        {
            ids[i] = bn;
            bvect_full2.set(bn);
            bvect_min1.set_bit(bn);
            bn += 15;
        }
        
        bvect::size_type* first = ids;
        bvect::size_type* last = ids + to_add;
        
        bm::combine_xor(bvect_full1, first, last);

        CheckVectors(bvect_min1, bvect_full1, BITVECT_SIZE, true);
        
        bm::combine_xor(bvect_full1, first, last);
        if (bvect_full1.count())
        {
            cout << "combine_xor count failed!" << endl;
            assert(0);
            exit(1);
        }
    }

    {
        bvect        bvect_full1;
        bvect        bvect_full2;
        bvect_mini   bvect_min1(BITVECT_SIZE);
        
        bvect_full1.set_new_blocks_strat(bm::BM_GAP);
        bvect_full2.set_new_blocks_strat(bm::BM_GAP);

        bvect::size_type ids[10000]={0,};
        bvect::size_type to_add = 10000;
        for (bvect::size_type i = 0; i < to_add; i+=100)
        {
            ids[i] = bm::id_max32 + i;
            bvect_full2.set(bm::id_max32 + i);
            bvect_min1.set_bit(bm::id_max32 + i);
        }
        bvect::size_type* first = ids;
        bvect::size_type* last = ids + to_add;
        
        bm::combine_xor(bvect_full1, first, last);

        CheckVectors(bvect_min1, bvect_full1, BITVECT_SIZE, true);
        
        bm::combine_xor(bvect_full1, first, last);
        if (bvect_full1.count())
        {
            cout << "combine_xor count failed!" << endl;
            exit(1);
        }
    }

    
    {
        bvect::size_type ids[] = {0, 65536, 65535, 65535*3, 65535*2, 10};
        unsigned to_add = sizeof(ids)/sizeof(bvect::size_type);
        bvect        bvect_full1;
        bvect        bvect_full2;
        bvect_mini   bvect_min1(BITVECT_SIZE);

        bvect_full1.set_new_blocks_strat(bm::BM_BIT);
        bvect_full2.set_new_blocks_strat(bm::BM_BIT);
        
        bvect::size_type bn = 0;
        for (unsigned i = 0; i < to_add; ++i)
        {
            ids[i] = bn;
            bvect_full2.set(bn);
            bvect_min1.set_bit(bn);
            bn += 15;
        }
        
        bvect::size_type* first = ids;
        bvect::size_type* last = ids + to_add;
        
        bm::combine_xor(bvect_full1, first, last);
        CheckVectors(bvect_min1, bvect_full1, BITVECT_SIZE, true);

        bm::combine_xor(bvect_full1, first, last);
        if (bvect_full1.count())
        {
            cout << "combine_xor count failed!" << endl;
            assert(0); exit(1);
        }
    }
    
    
    {
        bvect::size_type ids[] = {0, 65536, 65535, 65535*3, 65535*2, 10};
        unsigned to_add = sizeof(ids)/sizeof(bvect::size_type);
        bvect        bvect_full1;
        bvect        bvect_full2;
        bvect_mini   bvect_min1(BITVECT_SIZE);

        bvect_full1.set_new_blocks_strat(bm::BM_GAP);
        bvect_full2.set_new_blocks_strat(bm::BM_GAP);
        
        bvect::size_type bn = bm::id_max32;
        for (unsigned i = 0; i < to_add; ++i)
        {
            ids[i] = bn;
            bvect_full2.set(bn);
            bvect_min1.set_bit(bn);
            bn += 15;
        }
        
        bvect::size_type* first = ids;
        bvect::size_type* last = ids + to_add;
        
        bm::combine_xor(bvect_full1, first, last);
        CheckVectors(bvect_min1, bvect_full1, BITVECT_SIZE, true);

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
void SubOperationsTest()
{
    const bvect::size_type BITVECT_SIZE = bm::id_max32 * 2;
    assert(ITERATIONS < BITVECT_SIZE);

    cout << "----------------------------------- SubOperationTest" << endl;

    {
        bvect bv_i0;
        bv_i0.invert();
        bvect bv_i1;
        bv_i1.invert();
        auto cnt = bm::count_sub(bv_i0, bv_i1);
        assert(!cnt);
        auto predicted_any = bm::any_sub(bv_i0, bv_i1);
        assert(!predicted_any);
    }

    {
        ref_vect vect;
        generate_vect_simpl0(vect);

        bvect bv0;
        load_BV_set_ref(cout,bv0, vect);
        bvect bv1(bv0);
        
        {
            bvect bv3;
            bv3.bit_sub(bv0, bv1, bvect::opt_compress);
            assert(!bv3.count());
            assert(!bv3.any());
            
            bvect bv4;
            int cmp = bv4.compare(bv3);
            assert(cmp == 0);
            
            bv3 = bv0;
            bv3 -= bv1;
            assert(!bv3.count());
            assert(!bv3.any());
            cmp = bv4.compare(bv3);
            assert(cmp == 0);
        }
        
        bvect bv_i;
        bv_i.invert();
 
        for (unsigned i = 0; i < 2; ++i)
        {
            bvect::size_type predicted_count = bm::count_sub(bv0, bv1);
            assert(predicted_count == 0);
            auto predicted_any = bm::any_sub(bv1, bv0);
            assert(!predicted_any);

            predicted_count = bm::count_sub(bv0, bv_i);
            assert(predicted_count == 0);
            
            {
                bvect bv3(bv_i);
                clear_BV_set_ref(cout, bv3, vect);
                bvect bv4(bv_i);
                bv4 -= bv0;
                
                int cmp = bv3.compare(bv4);
                assert(cmp == 0);
                
                predicted_count = bm::count_sub(bv4, bv3);
                assert(!predicted_count);
                
                auto pany = bm::any_sub(bv4, bv3);
                assert(!pany);
            }
            
            bv0.optimize();
            bv1.optimize();
        } // for
        
    }

    {
        cout << "\n48-bit large set difference test" << endl;
        cout << "generation..." << endl;
        ref_vect vect0;
        generate_vect48(vect0);
        bvect bv0;
        load_BV_set_ref(cout,bv0, vect0);
        compare_BV_set_ref(bv0, vect0);

        ref_vect vect1;
        generate_vect48(vect1);
        bvect bv1;
        load_BV_set_ref(cout, bv1, vect1);
        compare_BV_set_ref(bv1, vect1);
        cout << "ok\n" << endl;

        ref_vect vect_i;
        std::set_difference(vect0.begin(), vect0.end(),
                            vect1.begin(), vect1.end(),
                            std::back_inserter(vect_i));
        {
            bvect bv_i(bv0);
            bv_i.bit_sub(bv1);
            compare_BV_set_ref(bv_i, vect_i);
        }
        bv0.optimize();
        {
            bvect bv_i(bv0);
            bv_i.bit_sub(bv1);
            compare_BV_set_ref(bv_i, vect_i);
        }
        bv1.optimize();
        {
            bvect bv_i(bv0);
            bv_i.bit_sub(bv1);
            compare_BV_set_ref(bv_i, vect_i);
        }
        bvect bvt;
        SerializationOperation2Test(&bvt,
                                    bv0,
                                    bv1,
                                    vect_i.size(),
                                    set_COUNT_SUB_AB,
                                    set_SUB);
    }



    {
        bvect_mini   bvect_min1(BITVECT_SIZE);
        bvect_mini   bvect_min2(BITVECT_SIZE);
        bvect        bvect_full1;
        bvect        bvect_full2;

        bvect_full1.set_new_blocks_strat(bm::BM_GAP);
        bvect_full2.set_new_blocks_strat(bm::BM_GAP);


        printf("SUB test\n");

        bvect_min1.set_bit(1);
        bvect_min1.set_bit(12);
        bvect_min1.set_bit(bm::id_max32+13);

        bvect_min2.set_bit(12);
        bvect_min2.set_bit(bm::id_max32+13);

        bvect_min1.combine_sub(bvect_min2);

        bvect_full1.set_bit(1);
        bvect_full1.set_bit(12);
        bvect_full1.set_bit(bm::id_max32+13);

        bvect_full2.set_bit(12);
        bvect_full2.set_bit(bm::id_max32+13);

        auto predicted_count = bm::count_sub(bvect_full1, bvect_full2);
        auto predicted_any = bm::any_sub(bvect_full1, bvect_full2);
        if (predicted_any == 0 && predicted_count != 0)
        {
            cout << "Predicted any error!" << endl;
            assert(0);
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
        
        auto count = bvect_full1.count();
        if (count != predicted_count)
        {
            cout << "Predicted count error!" << endl;
            assert(0);exit(1);
        }

        CheckVectors(bvect_min1, bvect_full1, 256, true);
        CheckVectors(bvect_min1, bv_target_s, 256, true);
    }

    {
        bvect_mini   bvect_min1(BITVECT_SIZE);
        bvect_mini   bvect_min2(BITVECT_SIZE);
        bvect        bvect_full1;
        bvect        bvect_full2;

        bvect_full1.set_new_blocks_strat(bm::BM_GAP);
        bvect_full2.set_new_blocks_strat(bm::BM_GAP);

        printf("SUB test stage 2.\n");

        FillSets(&bvect_min1, &bvect_full1, 1ull, BITVECT_SIZE/7, 0ull);
        FillSets(&bvect_min2, &bvect_full2, 1ull, BITVECT_SIZE/7, 0ull);

        bvect_min1.combine_sub(bvect_min2);

        auto predicted_count = bm::count_sub(bvect_full1, bvect_full2);
        auto predicted_any = bm::any_sub(bvect_full1, bvect_full2);
        if (predicted_any == 0 && predicted_count != 0)
        {
            cout << "Predicted any error!" << endl;
            assert(0);exit(1);
        }

        bvect    bv_target_s;
        SerializationOperation2Test(&bv_target_s,
                                    bvect_full1,
                                    bvect_full2,
                                    predicted_count,
                                    set_COUNT_SUB_AB,
                                    set_SUB);

        bvect_full1.bit_sub(bvect_full2);
        
        auto count = bvect_full1.count();
        if (count != predicted_count)
        {
            cout << "Predicted count error!" << endl;
            cout << predicted_count << " " << count << endl;
            print_stat(cout,bvect_full1);
            assert(0);
            exit(1);
        }
        CheckVectors(bvect_min1, bvect_full1, BITVECT_SIZE/10+10, true);
        CheckVectors(bvect_min1, bv_target_s, BITVECT_SIZE/10+10, true);

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


        FillSets(&bvect_min1, &bvect_full1, 1ull, BITVECT_SIZE/5, 2ull);
        FillSets(&bvect_min2, &bvect_full2, 1ull, BITVECT_SIZE/5, 2ull);

        bvect_min1.combine_sub(bvect_min2);
        
        auto predicted_count = bm::count_sub(bvect_full1, bvect_full2);
        auto predicted_any = bm::any_sub(bvect_full1, bvect_full2);
        if (predicted_any == 0 && predicted_count != 0)
        {
            cout << "Predicted any error!" << endl;
            assert(0);exit(1);
        }

        bvect    bv_target_s;
        SerializationOperation2Test(&bv_target_s,
                                    bvect_full1,
                                    bvect_full2,
                                    predicted_count,
                                    set_COUNT_SUB_AB,
                                    set_SUB);

        bvect_full1.bit_sub(bvect_full2);

        auto count = bvect_full1.count();
        if (count != predicted_count)
        {
            cout << "Predicted count error!" << endl;
            exit(1);
        }


        CheckVectors(bvect_min1, bvect_full1, BITVECT_SIZE, true);

        BM_DECLARE_TEMP_BLOCK(tb)
        bvect_full1.optimize(tb);
        CheckVectors(bvect_min1, bvect_full1, BITVECT_SIZE, true);
        CheckVectors(bvect_min1, bv_target_s, BITVECT_SIZE, true);
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

template<typename BV>
void CheckBV_AND_OR(BV& bv_target, const BV& bv1, const BV& bv2)
{
    BV bv_control(bv_target);
    BV bv_t_copy(bv_target);
    {
        BV bv_and;
        bv_and.bit_and(bv1, bv2, bvect::opt_compress);
        bv_t_copy |= bv_and;
    }

    bv_target.bit_or_and(bv1, bv2, bvect::opt_compress);
    bool f;
    typename BV::size_type pos;
    f = bv_target.find_first_mismatch(bv_t_copy, pos);
    if (f)
    {
        cerr << "AND-OR Mismatch at:" << pos << endl;
        unsigned nb = unsigned(pos >> bm::set_block_shift);
        unsigned i,j;
        bm::get_block_coord(nb, i, j);
        cout << "nb=" << nb << " i=" << i << " j=" << j << endl;

        bool v1 = bv_target.test(pos);
        bool vC = bv_t_copy.test(pos);
        cout << "v1=" << v1 << " control=" << vC << endl;

        bv_control.bit_or_and(bv1, bv2, bvect::opt_compress);

        assert(0);
    }

}

static
void AndOrOperationsTest(bool detailed)
{
    (void)detailed;
    cout << "----------------------------------- AndOrOperationTest()" << endl;

    {
        bvect  bvtarget;
        bvect  bv1 { 0, 1 }, bv2 { 1, 3 };
        CheckBV_AND_OR(bvtarget, bv1, bv2);
        assert(bvtarget.count() == 1);
    }

    {
        bvect  bvtarget;
        bvect  bv1, bv2;
        bv1.invert();
        bv2.invert();
        CheckBV_AND_OR(bvtarget, bv1, bv2);
    }
    {
        bvect  bvtarget { 1, 10, 65536 };
        bvect  bv1, bv2;
        bv1.invert();
        bv2.invert();
        CheckBV_AND_OR(bvtarget, bv1, bv2);
    }
    {
        bvect  bvtarget {1, 256, 65536 } ;
        bvect  bv1 { 0, 1 }, bv2 { 1, 3 };
        CheckBV_AND_OR(bvtarget, bv1, bv2);
        auto cnt = bvtarget.count();
        assert(cnt == 3);
    }
    {
        bvect  bvtarget {1, 256, 65536 } ;
        bvect  bv1 { 0, 1 }, bv2 { 1, 3 };
        bvtarget.optimize();
        CheckBV_AND_OR(bvtarget, bv1, bv2);
        auto cnt = bvtarget.count();
        assert(cnt == 3);
    }
    {
        bvect  bvtarget {1, 256, 65536 } ;
        bvect  bv1 { 0, 1 }, bv2 { 1, 3 };
        bv1.optimize();
        CheckBV_AND_OR(bvtarget, bv1, bv2);
        auto cnt = bvtarget.count();
        assert(cnt == 3);
    }
    {
        bvect  bvtarget {1, 256, 65536 } ;
        bvect  bv1 { 0, 1 }, bv2 { 1, 3 };
        bv2.optimize();
        CheckBV_AND_OR(bvtarget, bv1, bv2);
        auto cnt = bvtarget.count();
        assert(cnt == 3);
    }
    {
        bvect  bvtarget {1, 256, 65536 } ;
        bvect  bv1 { 0, 1 }, bv2 { 1, 3 };
        bv1.optimize();
        bv2.optimize();
        CheckBV_AND_OR(bvtarget, bv1, bv2);
        auto cnt = bvtarget.count();
        assert(cnt == 3);
    }
    {
        bvect  bvtarget;
        bvect  bv1 { 0, 1 }, bv2 { 2, 3 };
        CheckBV_AND_OR(bvtarget, bv1, bv2);
        auto cnt = bvtarget.count();
        assert(cnt == 0);
    }
    {
        bvect  bvtarget;
        bvect  bv1, bv2 { 2, 3, bm::id_max/2 };
        bv1.invert();
        CheckBV_AND_OR(bvtarget, bv1, bv2);
        auto cnt = bvtarget.count();
        assert(cnt == 3);
    }
    {
        bvect  bvtarget;
        bvect  bv2, bv1 { 2, 3, bm::id_max/2 };
        bv2.invert();
        CheckBV_AND_OR(bvtarget, bv1, bv2);
        auto cnt = bvtarget.count();
        assert(cnt == 3);
    }


    {
        bvect  bvtarget {1, 256, 65536 } ;
        bvect  bv1 { 0, 1 }, bv2 { 1, 3 };
        bvtarget.optimize();
        bv1.optimize();
        bv2.optimize();
        CheckBV_AND_OR(bvtarget, bv1, bv2);
        auto cnt = bvtarget.count();
        assert(cnt == 3);
    }

    {
        bvect  bvtarget;
        bvect  bv1, bv2;
        bvtarget.set_range(2, 65535);
        bv1.set_range(0,1);
        bv2.set_range(0,1);
        CheckBV_AND_OR(bvtarget, bv1, bv2);
        auto cnt = bvtarget.count();
        assert(cnt == 65536);
        bvect::statistics st;
        bvtarget.calc_stat(&st);
        assert(st.gap_blocks==0 && st.bit_blocks==0);
    }


    cout << "----------------------------------- AndOrOperationTest OK" << endl;
}

static
void TestAND_OR(bm::random_subset<bvect>& rsub,
                bvect::size_type count,
                const bvect& bvect_full1, const bvect& bvect_full2)
{
    cout << "AND-OR tests..." << flush;
    bvect bv_sub1;
    auto sample_count = count / 2;
    if (sample_count)
        rsub.sample(bv_sub1, bvect_full1, sample_count);
    CheckBV_AND_OR(bv_sub1, bvect_full1, bvect_full2);
    if (sample_count)
        rsub.sample(bv_sub1, bvect_full2, sample_count);
    CheckBV_AND_OR(bv_sub1, bvect_full2, bvect_full1);
    bv_sub1 = bvect_full1;
    CheckBV_AND_OR(bv_sub1, bvect_full2, bvect_full1);
    bv_sub1 = bvect_full2;
    CheckBV_AND_OR(bv_sub1, bvect_full1, bvect_full2);
    cout << " OK" << endl;
}


static
void StressTest(unsigned repetitions, int set_operation = -1)
{
    const bvect::size_type BITVECT_SIZE = bm::id_max32 + (bm::id_max32 / 2);

   bvect::size_type RatioSum = 0;
   bvect::size_type SRatioSum = 0;
   bvect::size_type DeltaSum = 0;
   bvect::size_type SDeltaSum = 0;

   bvect::size_type clear_count = 0;

   bvect  bvtotal;
   bvtotal.set_new_blocks_strat(bm::BM_GAP);

   bm::random_subset<bvect> rsub;

    bvect bv_p0, bv_p1;
    {
        ref_vect vect;
        generate_vect_simpl0(vect);
        load_BV_set_ref(cout,bv_p0, vect);
    }
    {
        ref_vect vect;
        generate_vect48(vect);
        load_BV_set_ref(cout,bv_p1, vect);
    }


   cout << "----------------------------StressTest" << endl;

   bvect::size_type size = BITVECT_SIZE - 10;
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

        bvect::size_type start1 = 0;
        if (rand()%3)
            start1 += size / 3;
        else
            start1 += size / 5;
        bvect::size_type start2 = 0;
        if (rand()%3)
            start2 += size / 3;
        else
            start2 += size / 5;

        FillSetsRandomMethod(bvect_min1, bvect_full1, start1, size, opt);
        FillSetsRandomMethod(bvect_min2, bvect_full2, start2, size, opt);

       
        // commented out because it crashes Apple CLang compiler ...
/*
        unsigned arr[bm::set_total_blocks]={0,};
        bvect::size_type cnt = bvect_full1->count();
        bvect::size_type last_block = bvect_full1->count_blocks(arr);

        bvect::size_type sum = bm::sum_arr(&arr[0], &arr[last_block+1]);
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

            assert(0);
            exit(1);
        }
*/
        //CheckCountRange(*bvect_full1, start1, BITVECT_SIZE);
        //CheckIntervals(*bvect_full1, BITVECT_SIZE);


        //CheckCountRange(*bvect_full2, start2, BITVECT_SIZE);

        //CheckCountRange(*bvect_full1, 0, start1);
        //CheckCountRange(*bvect_full2, 0, start2);


        TestRandomSubset(*bvect_full1, rsub);
        TestRandomSubset(*bvect_full2, rsub);

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
                auto predicted_count = bm::count_or(*bvect_full1, *bvect_full2);
                auto predicted_any = bm::any_or(*bvect_full1, *bvect_full2);
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
                auto count = bvect_full1->count();
                if (count != predicted_count)
                {
                    cout << "Predicted count error!" << endl;
                    cout << "Count = " << count << "Predicted count = " << predicted_count << endl;
                    assert(0);
                    exit(1);
                }
                int res = bvect_full1->compare(bv_target_s);
                if (res != 0)
                {
                    cout << "Serialization operation failed!" << endl;
                    assert(0);exit(1);
                }

                TestAND_OR(rsub, predicted_count, *bvect_full1, *bvect_full2);

                // run cross checks with wide-band vectors
                //
                cout << "48-wide vectors checks.." << endl;
                bv_target_s.clear();
                SerializationOperation(&bv_target_s, *bvect_full1, bv_p0, set_OR, true);
                bv_target_s.clear();
                SerializationOperation(&bv_target_s, *bvect_full2, bv_p0, set_OR, true);
                bv_target_s.clear();
                SerializationOperation(&bv_target_s, *bvect_full1, bv_p0, set_OR, true);
                bv_target_s.clear();
                SerializationOperation(&bv_target_s, *bvect_full2, bv_p0, set_OR, true);
                
                bv_target_s.clear();
                SerializationOperation(&bv_target_s, *bvect_full1, bv_p1, set_OR, true);
                bv_target_s.clear();
                SerializationOperation(&bv_target_s, *bvect_full2, bv_p1, set_OR, true);
                bv_target_s.clear();
                SerializationOperation(&bv_target_s, *bvect_full1, bv_p1, set_OR, true);
                bv_target_s.clear();
                SerializationOperation(&bv_target_s, *bvect_full2, bv_p1, set_OR, true);
            }
            break;

        case 1:
            {
            cout << "Operation SUB" << endl;
            
            auto predicted_count = bm::count_sub(*bvect_full1, *bvect_full2);
            auto predicted_any = bm::any_sub(*bvect_full1, *bvect_full2);
            if (predicted_any == 0 && predicted_count != 0)
            {
                cout << "Predicted any error!" << endl;
                assert(0); exit(1);
            }
            
            bvect    bv_target_s;
            SerializationOperation2Test(&bv_target_s,
                                        *bvect_full1,
                                        *bvect_full2,
                                        predicted_count,
                                        set_COUNT_SUB_AB,
                                        set_SUB);

            bvect_full1->bit_sub(*bvect_full2);
            auto count = bvect_full1->count();
            if (count != predicted_count)
            {
                cout << "Predicted count error!" << endl;
                cout << "Count = " << count << "Predicted count = " << predicted_count << endl;
                assert(0); exit(1);
            }
            int res = bvect_full1->compare(bv_target_s);
            if (res != 0)
            {
                cout << "Serialization operation failed!" << endl;
                assert(0); exit(1);
            }
                // run cross checks with wide-band vectors
                //
                cout << "48-wide vectors checks.." << endl;
                bv_target_s.clear();
                SerializationOperation(&bv_target_s, *bvect_full1, bv_p0, set_SUB, true);
                bv_target_s.clear();
                SerializationOperation(&bv_target_s, *bvect_full2, bv_p0, set_SUB, true);
                bv_target_s.clear();
                SerializationOperation(&bv_target_s, *bvect_full1, bv_p0, set_SUB, true);
                bv_target_s.clear();
                SerializationOperation(&bv_target_s, *bvect_full2, bv_p0, set_SUB, true);
                
                bv_target_s.clear();
                SerializationOperation(&bv_target_s, *bvect_full1, bv_p1, set_SUB, true);
                bv_target_s.clear();
                SerializationOperation(&bv_target_s, *bvect_full2, bv_p1, set_SUB, true);
                bv_target_s.clear();
                SerializationOperation(&bv_target_s, *bvect_full1, bv_p1, set_SUB, true);
                bv_target_s.clear();
                SerializationOperation(&bv_target_s, *bvect_full2, bv_p1, set_SUB, true);
            }
            break;

        case 2:
            {
            cout << "Operation XOR <<<" << endl;
           
            auto predicted_count = bm::count_xor(*bvect_full1, *bvect_full2);
            auto predicted_any = bm::any_xor(*bvect_full1, *bvect_full2);
            if (predicted_any == 0 && predicted_count != 0)
            {
                cout << "Predicted any error!" << endl;
                assert(0); exit(1);
            }

            bvect    bv_target_s;
            SerializationOperation2Test(&bv_target_s,
                                        *bvect_full1,
                                        *bvect_full2,
                                        predicted_count,
                                        set_COUNT_XOR,
                                        set_XOR);
            
            bvect_full1->bit_xor(*bvect_full2);
            
            auto count = bvect_full1->count();
            if (count != predicted_count)
            {
                cout << "Predicted count error!" << endl;
                cout << "Count = " << count << "Predicted count = " << predicted_count << endl;
                assert(0); exit(1);
            }
            int res = bvect_full1->compare(bv_target_s);
            if (res != 0)
            {
                cout << "Serialization operation failed!" << endl;
                exit(1);
            }
                // run cross checks with wide-band vectors
                //
                cout << "48-wide vectors checks.." << endl;
                bv_target_s.clear();
                SerializationOperation(&bv_target_s, *bvect_full1, bv_p0, set_XOR, true);
                bv_target_s.clear();
                SerializationOperation(&bv_target_s, *bvect_full2, bv_p0, set_XOR, true);
                bv_target_s.clear();
                SerializationOperation(&bv_target_s, *bvect_full1, bv_p0, set_XOR, true);
                bv_target_s.clear();
                SerializationOperation(&bv_target_s, *bvect_full2, bv_p0, set_XOR, true);
                
                bv_target_s.clear();
                SerializationOperation(&bv_target_s, *bvect_full1, bv_p1, set_XOR, true);
                bv_target_s.clear();
                SerializationOperation(&bv_target_s, *bvect_full2, bv_p1, set_XOR, true);
                bv_target_s.clear();
                SerializationOperation(&bv_target_s, *bvect_full1, bv_p1, set_XOR, true);
                bv_target_s.clear();
                SerializationOperation(&bv_target_s, *bvect_full2, bv_p1, set_XOR, true);
            }
            
            break;

        default:
            {
            cout << "Operation AND" << endl;

            auto predicted_count = bm::count_and(*bvect_full1, *bvect_full2);
            auto predicted_any = bm::any_and(*bvect_full1, *bvect_full2);
            if (predicted_any == 0 && predicted_count != 0)
            {
                cout << "Predicted any error!" << endl;
                assert(0); exit(1);
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
            auto count = bvect_full1->count();

            int res = bvect_full1->compare(bv_target_s);
            if (res != 0)
            {
                //SaveBVector("bv1.bv", bv1);
                //SaveBVector("bv2.bv", *bvect_full2);
                cout << "Serialization operation failed!" << endl;
                assert(0); exit(1);
            }

            if (count != predicted_count)
            {
                cout << "Predicted count error!" << endl;
                cout << "Count = " << count << "Predicted count = " << predicted_count << endl;
                exit(1);
            }
            TestAND_OR(rsub, predicted_count, *bvect_full1, *bvect_full2);

                // run cross checks with wide-band vectors
                //
                cout << "48-wide vectors checks.." << endl;
                bv_target_s.clear();
                SerializationOperation(&bv_target_s, *bvect_full1, bv_p0, set_AND, true);
                bv_target_s.clear();
                SerializationOperation(&bv_target_s, *bvect_full2, bv_p0, set_AND, true);
                bv_target_s.clear();
                SerializationOperation(&bv_target_s, *bvect_full1, bv_p0, set_AND, true);
                bv_target_s.clear();
                SerializationOperation(&bv_target_s, *bvect_full2, bv_p0, set_AND, true);
                
                bv_target_s.clear();
                SerializationOperation(&bv_target_s, *bvect_full1, bv_p1, set_AND, true);
                bv_target_s.clear();
                SerializationOperation(&bv_target_s, *bvect_full2, bv_p1, set_AND, true);
                bv_target_s.clear();
                SerializationOperation(&bv_target_s, *bvect_full1, bv_p1, set_AND, true);
                bv_target_s.clear();
                SerializationOperation(&bv_target_s, *bvect_full2, bv_p1, set_AND, true);
            }
            break;
        }



        cout << "Operation comparison" << endl;
        CheckVectors(*bvect_min1, *bvect_full1, size, true);

        int cres2 = bvect_full1->compare(*bvect_full2);

        //CheckIntervals(*bvect_full1, BITVECT_SIZE);

        if (cres1 != cres2)
        {
            cout << cres1 << " " << cres2 << endl;
            cout << bvect_full1->get_first() << " " << bvect_full1->count() << endl;
            cout << bvect_full2->get_first() << " " << bvect_full2->count() << endl;

           // bvect_full1->stat(1000);
            cout << endl;
           // bvect_full2->stat(1000);
            printf("Bitset comparison operation failed.\n");
            exit(1);
        }

        {
            bvect bv1(*bvect_full1);
            bvect::size_type idx = bvect::size_type(rand()) % size;
            bool b = bv1[idx];
            bool changed;
            if (b)
            {
                changed = bv1.set_bit_conditional(idx, true, false);
                if (changed)
                {
                    cout << "Set bit conditional failed!" << endl;
                    assert(0); exit(1);
                }
                b = bv1[idx];
                if (!b)
                {
                    cout << "Set bit conditional failed!" << endl;
                    assert(0); exit(1);
                }

                changed = bv1.set_bit_conditional(idx, false, false);
                if (changed)
                {
                    cout << "Set bit conditional failed!" << endl;
                    assert(0); exit(1);
                }
                changed = bv1.set_bit_conditional(idx, true, true);
                if (changed)
                {
                    cout << "Set bit conditional failed!" << endl;
                    assert(0); exit(1);
                }
                changed = bv1.set_bit_conditional(idx, false, true);
                if (!changed)
                {
                    cout << "Set bit conditional failed!" << endl;
                    assert(0);exit(1);
                }
                b = bv1[idx];
                if (b)
                {
                    cout << "Set bit conditional failed!" << endl;
                    assert(0); exit(1);
                }
            }
            else
            {
                changed = bv1.set_bit_conditional(idx, false, true);
                if (changed)
                {
                    cout << "Set bit conditional failed!" << endl;
                    assert(0);exit(1);
                }
                changed = bv1.set_bit_conditional(idx, true, false);
                if (!changed)
                {
                    cout << "Set bit conditional failed!" << endl;
                    assert(0);exit(1);
                }
                b = bv1[idx];
                if (!b)
                {
                    cout << "Set bit conditional failed!" << endl;
                    assert(0);exit(1);
                }
            }


        }

        {
            VisitorAllRangeTest(*bvect_full2, 0);
        }

        delete bvect_full2;


        struct bvect::statistics st1;
        bvect_full1->calc_stat(&st1);
        bvect_full1->optimize();
        struct bvect::statistics st2;
        bvect_full1->calc_stat(&st2);

        bvect::size_type Ratio = bvect::size_type((st2.memory_used * 100)/st1.memory_used);
        RatioSum+=Ratio;
        DeltaSum+=unsigned(st1.memory_used - st2.memory_used);

        cout << "Optimization statistics: " << endl
             << "   MemUsedBefore=" << st1.memory_used
             << "   MemUsed=" << st2.memory_used
             << "   Ratio=" << Ratio << "%"
             << "   Delta=" << st1.memory_used - st2.memory_used
             << endl;
       
        cout << "Optimization comparison" << endl;

        CheckVectors(*bvect_min1, *bvect_full1, size, true);

        bvect_full1->set_gap_levels(gap_len_table_min<true>::_len);
        CheckVectors(*bvect_min1, *bvect_full1, size, true);
        //CheckIntervals(*bvect_full1, BITVECT_SIZE);

        //CheckCountRange(*bvect_full1, 0, size);


        // Serialization
       
        bm::serializer<bvect> bv_ser;
        bv_ser.gap_length_serialization(false);
        bv_ser.byte_order_serialization(false);
        bv_ser.set_bookmarks(true, 128);
       
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

        bvect* bv_target_s=new bvect();
        operation_deserializer<bvect64> od;
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

        CheckVectors(*bvect_min1, *bvect_full3, size, true);
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

    cout << " AGG arg pipeline tests" << endl;

    {
        bm::aggregator<bvect>::pipeline agg_pipe;
        {
            bm::aggregator<bvect>::arg_groups* args = agg_pipe.add();
            assert(args);
            args->arg_bv0.push_back(nullptr);
            args->arg_bv1.push_back(nullptr);
        }
        {
            bm::aggregator<bvect>::arg_groups* agrs = agg_pipe.add();
            agrs->arg_bv0.push_back(nullptr);
        }
        auto& arg_vect = agg_pipe.get_args_vector();
        assert(arg_vect.size() == 2);
        assert(arg_vect[0]->arg_bv0.size() == 1);
        assert(arg_vect[0]->arg_bv1.size() == 1);
        assert(arg_vect[1]->arg_bv0.size() == 1);
        assert(arg_vect[1]->arg_bv1.size() == 0);
    }

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
            
            auto bc = bv_target.count();
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
            
            auto bc = bv_target.count();
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
            
            auto bc = bv_target.count();
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

    cout << "AND tests..." << endl;
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
    cout << "AND with range tests..." << endl;
    {
        bvect bv1, bv2, bv3, bv4;
        bvect bv_empty;
        bvect bv_control;
        int res;

        bv1.invert();
        bvect::size_type from = bm::id_max32-200000;
        bvect::size_type to = bm::id_max32+300000;
        bv2.set_range(from, to);
        bv3.set_range(from, to);
        
        bv_control.set_range(from, to);
        
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
    cout << "AND-SUB tests..." << endl;
    {
        bvect bv1, bv2, bv3, bv4;
        bvect bv_empty;
        bvect bv_control;
        
        bv1[100] = true;
        bv1[bm::id_max32+100000] = true;
        bv2[100] = true;
        bv2[bm::id_max32+100000] = true;
        bv3[bm::id_max32+100000] = true;

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
        bv1[bm::id_max32+100000] = true;
        bv2[100] = true;
        bv2[bm::id_max32+100000] = true;
        
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
        bv1[bm::id_max32+65535]=true;

        bv2[0]=true;
        bv2[bm::id_max32+65535]=true;
        
        bool any = agg_shift_right_and(agg, bv0, &bv1, &bv2, 0);
        assert(!any);
        assert(bv0.count()==0);
    }


    {
        bvect bv0, bv1, bv2;
        bv1[0] = true;
        bv1[bm::id_max32+65535]=true;
        bv1.optimize();

        bv2[1]=true;
        bv2[bm::id_max32+65536]=true;
        bv2.optimize();

        agg_shift_right_and(agg, bv0, &bv1, &bv2, 0);
        assert(bv0.count()==2);
        assert(bv0.test(1));
        assert(bv0.test(bm::id_max32+65536));
    }


    {
        bvect bv0, bv1, bv2;
        bv1[bm::id_max32+65535]=true;
        
        bv2[bm::id_max32+65536]=true;
        bv2.optimize();
        
        bool any = agg_shift_right_and(agg, bv0, &bv1, &bv2, 0);
        assert(bv0.count()==1);
        assert(bv0.test(bm::id_max32+65536));
        assert(any);
        struct bvect::statistics st1;
        bv0.calc_stat(&st1);
        auto bcnt = st1.bit_blocks + st1.gap_blocks;
        assert(bcnt == 1);
    }

    
    {
        bvect bv0, bv1, bv2;
        bv1[0] = true;
        bv1[bm::id_max32+65535]=true;

        bv2.invert();
        
        agg_shift_right_and(agg, bv0, &bv1, &bv2, 0);
        assert(bv0.count()==2);
        assert(bv0.test(1));
        assert(bv0.test(bm::id_max32+65536));
    }
    
    /*
    // TODO: optimize this case
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
    */
    
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
  bvect::size_type BITVECT_SIZE =
    bvect::size_type(bm::id_max32) + (bvect::size_type(bm::id_max32) / 2);

  cout << "---------------------------- Aggregator OR Stress Test" << endl;
   bvect::size_type size = BITVECT_SIZE - 10;
   bvect bv_target1, bv_target2;


    unsigned i;
    for (i = 0; i < repetitions; ++i)
    {
        int opt = 1;//rand() % 2;
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
        
        bvect::size_type start1 = 0;
        bvect::size_type start2 = 0;
        if (rand()%3 )
        {
            start1 += size / 2;
            start2 += size / 2;
        }
        bvect bv0, bv1, bv2, bv3, bv4, bv5, bv6, bv7, bv8, bv9;

        {
            bvect_mini   bvect_min1(BITVECT_SIZE);

            // 0 skipped
            FillSetsRandomMethod(&bvect_min1, &bv1, start1, size, opt);
            FillSetsRandomMethod(&bvect_min1, &bv2, start2, size, opt);
            // 3 skipped
            FillSetsRandomMethod(&bvect_min1, &bv5, start1, size, opt);
            FillSetsRandomMethod(&bvect_min1, &bv6, start2, size, opt);
            FillSetsRandomMethod(&bvect_min1, &bv7, start1, size, opt);
        }
        // bv8 and bv9 loaded as wide range vectors
        {
            ref_vect vect0;
            generate_vect48(vect0);
            load_BV_set_ref(cout,bv8, vect0);
            bv8.optimize();
        }
        {
            ref_vect vect0;
            generate_vect48(vect0);
            load_BV_set_ref(cout,bv9, vect0);
            bv9.optimize();
        }
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
            assert(0);exit(1);
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
                assert(0); exit(1);
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
                assert(0); exit(1);
            }
        }


    } // for i

  cout << "---------------------------- Aggregator OR Stress Test OK" << endl;
}


static
void StressTestAggregatorAND(unsigned repetitions)
{
   bvect::size_type BITVECT_SIZE =
      bvect::size_type(bm::id_max32) + (bvect::size_type(bm::id_max32) / 2);

   cout << "---------------------------- Aggregator AND Stress Test" << endl;
   bvect::size_type size = BITVECT_SIZE - 10;

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

        bvect::size_type start1 = 0;
        bvect::size_type start2 = 0;
        if (rand()%3 )
        {
            start1 += size / 2;
            start2 += size / 2;
        }

        bvect bv0, bv1, bv2, bv3, bv4, bv5, bv6, bv7, bv8, bv9;
        {
        bvect_mini   bvect_min1(size);

        // 0 skipped
        FillSetsRandomMethod(&bvect_min1, &bv1, start1, size, opt);
        FillSetsRandomMethod(&bvect_min1, &bv2, start2, size, opt);
        // 3 skipped
        FillSetsRandomMethod(&bvect_min1, &bv5, start1, size, opt);
        FillSetsRandomMethod(&bvect_min1, &bv6, start2, size, opt);
        FillSetsRandomMethod(&bvect_min1, &bv7, start1, size, opt);
        FillSetsRandomMethod(&bvect_min1, &bv8, start2, size, opt);
        FillSetsRandomMethod(&bvect_min1, &bv9, start2, size, opt);
        }
        
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
            assert(0);exit(1);
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
                assert(0);exit(1);
            }

            res = bv_target1.compare(bv_target3);
            if (res!=0)
            {
                cerr << "Error: Aggregator AND-SUB(0) check failed! 2.laddder step = "
                     << j << endl;
                assert(0);exit(1);
            }
            assert(!bv_empty.any());
        }


    } // for i

  cout << "---------------------------- Aggregator AND Stress Test OK" << endl;
}

static
void GenerateTestCollection(std::vector<bvect>* target,
                          unsigned count, bvect::size_type vector_max)
{
    assert(target);
    bvect bv_common; // sub-vector common for all collection
    bvect_mini bvect_min(vector_max);
    
    FillSetsRandomMethod(&bvect_min, &bv_common, 0ull, vector_max, 1);
    
    for (unsigned i = 0; i < count; ++i)
    {
        std::unique_ptr<bvect> bv (new bvect);
        FillSetsRandomMethod(&bvect_min, bv.get(), 0ull, vector_max, 1);
        *bv |= bv_common;
        target->push_back(std::move(*bv));
    } // for
}


static
void StressTestAggregatorShiftAND(unsigned repeats)
{
  bvect::size_type BITVECT_SIZE =
    bvect::size_type(bm::id_max32) + ( bvect::size_type(bm::id_max32) / 2);

   cout << "----------------------------StressTestAggregatorShiftAND " << endl;

    bvect::size_type vector_max = BITVECT_SIZE;
    unsigned coll_size = 20;

    for (unsigned r = 0; r < repeats; ++r)
    {
        bvect mask_bv0;
        {
            bvect_mini bvect_min(vector_max);
            FillSetsRandomMethod(&bvect_min, &mask_bv0, 0ull,
                                 vector_max - (vector_max / 5), 1);
        }
        
        std::vector<bvect> bv_coll1;
        GenerateTestCollection(&bv_coll1, coll_size, vector_max);

        bm::aggregator<bvect> agg;
        agg.set_optimization();

        unsigned shift_repeats = 65536/15;
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
            agg.set_compute_count(false);
            agg.combine_shift_right_and(bv_target1);
            auto cmp = bv_target1.compare(bv_target0);
            if (cmp != 0)
            {
                cerr << "Error: Mismatch! " << "STEP=" << i << endl;
                //DetailedCheckVectors(bv_target0, bv_target1);
                assert(0); exit(1);
            }
            bvect bv_target2;
            agg.set_compute_count(true);
            agg.combine_shift_right_and(bv_target2);
            assert(!bv_target2.any());
            auto cnt = agg.count();
            auto cnt_c = bv_target1.count();
            assert(cnt == cnt_c);

            if (i % 250 == 0)
                cout << "\r" << i << flush;

        } // for
        cout << "\n\n ---------- SHIFT-AND step: " << r << endl;
    } // for
    
   cout << "\n----------------------------StressTestAggregatorShiftAND OK" << endl;

}


static
void TestSparseVector()
{
    cout << "---------------------------- Bit-plain sparse vector test" << endl;

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
        unsigned arr[3] = {1, 2, 3};
        sv1.import(arr, 3);
        sv2.import(arr, 3);
        assert(!sv1.is_null(0));
        assert(!sv2.is_null(0));
        assert(!sv2.is_null(1));
        assert(!sv2.is_null(2));
        assert(sv2.is_null(3));
        
        assert(sv2.is_null(bm::id_max-1));
        sv2.set(bm::id_max-1, 123);
        assert(!sv2.is_null(bm::id_max-1));
        
        sv2.set_null(bm::id_max-1);
        assert(sv2.is_null(bm::id_max-1));
        assert(sv2[bm::id_max-1].is_null());

        
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
        sv2.clear_all(true);
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
        sv.push_back(2);
        sv.push_back(9);
        unsigned arr[1024];
        
        auto esize = sv.extract(&arr[0], 1024, 0);
        assert(esize == 3);
        assert(arr[0] == 1);
        assert(arr[1] == 2);
        assert(arr[2] == 9);
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
        found = scanner.bfind(sv1, 0u, pos);
        assert(!found);

        found = scanner.bfind(sv1, 1u, pos);
        assert(found);
        assert(pos == 0);

        found = scanner.bfind(sv1, 2u, pos);
        assert(found);
        assert(pos == 1);

        found = scanner.bfind(sv1, 3u, pos);
        assert(!found);

        found = scanner.bfind(sv1, 20u, pos);
        assert(found);

        found = scanner.bfind(sv1, 2000u, pos);
        assert(found);
    }


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
            assert(0);exit(1);
        }
        sv.optimize();
        print_svector_stat(cout,sv);
        res = CompareSparseVector(sv, vect);
        if (!res)
        {
            cerr << "optimized Bit Plain import test failed" << endl;
            assert(0);exit(1);
        }

        bm::sparse_vector<unsigned, bm::bvector<> > sv_1;
        std::copy(vect.begin(), vect.end(), std::back_inserter(sv_1));
        res = CompareSparseVector(sv_1, vect);
        if (!res)
        {
            cerr << "Bit Plain push_back test failed" << endl;
            assert(0);exit(1);
        }

        
        bm::sparse_vector<unsigned, bvect>::statistics st;
        sv.calc_stat(&st);
        
        bm::sparse_vector<unsigned, bvect> sv2(sv);
        res = CompareSparseVector(sv2, vect);
        if (!res)
        {
            cerr << "Bit Plain copy-ctor test failed" << endl;
            assert(0);exit(1);
        }
        
        sv2.clear();
        sv2.import(&vect[0], (unsigned)vect.size());
        res = CompareSparseVector(sv2, vect);
        if (!res)
        {
            cerr << "Bit Plain copy-ctor test failed" << endl;
            assert(0);exit(1);
        }

        bm::sparse_vector<unsigned, bvect> sv3;
        sv3.set(bm::id_max/2, 10); // set some bit to initiate it
        sv3 = sv;
        res = CompareSparseVector(sv3, vect);
        if (!res)
        {
            cerr << "Bit Plain assignmnet test failed" << endl;
            assert(0);exit(1);
        }
        
        sv3.clear_all(true);
        sv3.import(&vect[0], (unsigned)vect.size());
        res = CompareSparseVector(sv3, vect);
        if (!res)
        {
            cerr << "Bit Plain assignment test failed" << endl;
            assert(0);exit(1);
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
                assert(0);exit(1);
            }
        }}
        
        svector64 sv;
        sv.import(&vect[0], (unsigned)vect.size());
        bool res = CompareSparseVector(sv, vect);
        if (!res)
        {
            cerr << "0.Bit Plain import test failed" << endl;
            assert(0);exit(1);
        }
        sv.optimize();
        print_svector_stat(cout,sv);
        res = CompareSparseVector(sv, vect);
        if (!res)
        {
            cerr << "optimized Bit Plain import test failed" << endl;
            assert(0);exit(1);
        }

        svector64 sv_1;
        std::copy(vect.begin(), vect.end(), std::back_inserter(sv_1));
        res = CompareSparseVector(sv_1, vect);
        if (!res)
        {
            cerr << "Bit Plain push_back test failed" << endl;
            assert(0);exit(1);
        }

        
        svector64::statistics st;
        sv.calc_stat(&st);
        
        svector64 sv2(sv);
        res = CompareSparseVector(sv2, vect);
        if (!res)
        {
            cerr << "Bit Plain copy-ctor test failed" << endl;
            assert(0);exit(1);
        }
        
        sv2.clear();
        sv2.import(&vect[0], (unsigned)vect.size());
        res = CompareSparseVector(sv2, vect);
        if (!res)
        {
            cerr << "Bit Plain copy-ctor test failed" << endl;
            assert(0);exit(1);
        }

        svector64 sv3;
        sv3.set(65536, 10); // set some bit to initiate it
        sv3 = sv;
        res = CompareSparseVector(sv3, vect);
        if (!res)
        {
            cerr << "Bit Plain assignmnet test failed" << endl;
            assert(0);exit(1);
        }
        
        sv3.clear();
        sv3.import(&vect[0], (unsigned)vect.size());
        res = CompareSparseVector(sv3, vect);
        if (!res)
        {
            cerr << "Bit Plain assignment test failed" << endl;
            assert(0);exit(1);
        }
    }}

    cout << "64-bit sparse assignment test (1)" << endl;
    {{
        bm::sparse_vector<unsigned long long, bvect > sv(bm::use_null);
        ref_vect vect;
        generate_vect_simpl0(vect);
        load_SV_set_ref(&sv, vect);
        compare_SV_set_ref(sv, vect);
        cout << "optimization...";
        sv.optimize();
        cout << "ok" << endl;
        compare_SV_set_ref(sv, vect);
    }}

    cout << "64-bit sparse assignment test (2)" << endl;
    {{
        bm::sparse_vector<unsigned long long, bvect > sv(bm::use_null);
        ref_vect vect;
        generate_vect48(vect);
        load_SV_set_ref(&sv, vect);
        compare_SV_set_ref(sv, vect);
        cout << "optimization...";
        sv.optimize();
        cout << "ok" << endl;
        compare_SV_set_ref(sv, vect);
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
        assert(0);exit(1);
    }
    res = CompareSparseVector(sv1, vect);
    if (!res)
    {
        cerr << "linear assignment test failed (2)" << endl;
        assert(0);exit(1);
    }
    res = CompareSparseVector(sv2, vect);
    if (!res)
    {
        cerr << "linear assignment test failed (3 - back_inserter)" << endl;
        assert(0);exit(1);
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
    sv.extract_planes(&v1p[0], 16, 0);
    for (unsigned i = 0; i < 16; ++i)
    {
        if (v1[i] != 8 || v1r[i] != v1[i] || v1p[i] != v1[i])
        {
            cerr << "Extract 1 failed at:" << i << endl;
            assert(0);exit(1);
        }
    } // for
    
    std::vector<unsigned> v2(10);
    std::vector<unsigned> v2r(10);
    std::vector<unsigned> v2p(10);
    
    sv.extract(&v2[0], 10, 32);
    sv.extract_range(&v2r[0], 10, 32);
    sv.extract_planes(&v2p[0], 10, 32);
        
    for (unsigned i = 0; i < 10; ++i)
    {
        if (v2[i] != 255 || v2r[i] != v2[i] || v2p[i] != v2[i])
        {
            cerr << "Extract 2 failed at:" << i << "=" << v2[i] << " r=" << v2r[i] << " sv=" << sv[i+32] << endl;
            assert(0);exit(1);
        }
    } // for

    }}
    
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


    
    {{
        cout << "sparse vector inc test" << endl;
        bm::sparse_vector<unsigned long long, bvect > sv;
        
        bvect::size_type from = bm::id_max32;
        bvect::size_type to = from + 65536;
        bvect::size_type max_iter = 65536/5;
        for (bvect::size_type i = 1; i < max_iter; ++i)
        {
            for (auto j = from; j < to; ++j)
            {
                sv.inc(j);
                auto v = sv.get(j);
                assert(v == i);
            } // for j
            if ((i & 0xFF) == 0)
            {
                cout << "\r" << i << " / " << max_iter << flush;
                sv.optimize();
            }
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
            assert(0);exit(1);
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
            assert(0);exit(1);
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
                assert(0);exit(1);
            }
            assert(!sv1[i].is_null());
            v = sv1[i];
            if (v != i)
            {
                cerr << "Wrong null sparse vector value: at[" << i << "]=" << v << endl;
                assert(0);exit(1);
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
            assert(0);exit(1);
        }
        if (bv_null1->count() != 0)
        {
            cerr << "3. Incorrect sparse vector size() - NOT NULL comparison" << sv1.size() << " " << bv_null1->count() << endl;
            assert(0);exit(1);
        }

        
        sv.resize(65536);
        sv1.resize(65536);
        if (bv_null1->count() != 0)
        {
            cerr << "4. Incorrect sparse vector size() - NOT NULL comparison" << sv1.size() << " " << bv_null1->count() << endl;
            assert(0);exit(1);
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
                    assert(0);exit(1);
                }
            }
            assert(sv1[i].is_null());
        }
    }}
    
    {{
        cout << "Dynamic range clipping test 2" << endl;
        bm::sparse_vector<unsigned, bvect > sv;

        unsigned i;
        for (i = 128000; i < 128000 * 3; ++i)
        {
            sv.set(i, 7 + (unsigned)rand() % 256);
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
                    assert(0);exit(1);
                }
            }
            else
            {
                if (v > 15)
                {
                    cerr << "Clipped Value cmpr failed at:" << i << "=" << v << endl;
                    assert(0);exit(1);
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
                assert(0); exit(1);
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
        bm::sparse_vector<unsigned, bvect> sv1(bm::use_null);
        bm::sparse_vector<unsigned, bvect> sv2(bm::use_null);

        assert(sv1.is_nullable());
        
        sv1.set(0, 0);
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
                cerr << "Sparse join cmp failed:" << sv1.size() << endl;
                exit(1);
            }
            assert(!sv1[i].is_null());
        }
    }

    cout << "Test Sparse vector join NULL-able with not NULL-able" << endl;
    {
        bm::sparse_vector<unsigned, bvect> sv1(bm::use_null);
        bm::sparse_vector<unsigned, bvect> sv2(bm::use_null);

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
            //assert(sv1[i].is_null());
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
        sv1[bm::id_max32 + 1000] = 1;

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
        assert(pos == bm::id_max32 + 1000+1);
        f = bm::sparse_vector_find_first_mismatch(sv2, sv1, pos);
        assert(f);
        assert(pos == bm::id_max32 + 1000+1);

        sv1.optimize();
        f = bm::sparse_vector_find_first_mismatch(sv2, sv1, pos);
        assert(f);
        assert(pos == bm::id_max32 + 1000+1);

        sv2.optimize();
        f = bm::sparse_vector_find_first_mismatch(sv2, sv1, pos);
        assert(f);
        assert(pos == bm::id_max32 + 1000+1);
    }

    cout << " -------------------------- TestSparseVectorAlgo() OK" << endl;
}

// ---------------------------------------------------------------------


static void TestSignedSparseVector()
{
    cout << " -------------------------- TestSignedSparseVector()" << endl;

    {{
        bm::sparse_vector<int64_t, bvect> sv;
        sv.push_back(-1);
        sv.push_back(1);
        sv.push_back(INT64_MAX);
        sv.push_back(INT64_MIN);

        int64_t v0 = sv.get(0);
        assert(v0 == -1);
        int64_t v1 = sv[1];
        assert(v1 == 1);
        int64_t v2 = sv[2];
        assert(v2 == INT64_MAX);
        int64_t v3 = sv.get(3);
        assert(v3 == INT64_MIN);

        int64_t arr[1024];

        {
        auto esize =  sv.extract(&arr[0], 1024, 0);
        assert(esize == 4);
        assert(arr[0] == -1);
        assert(arr[1] == 1);
        assert(arr[2] == INT64_MAX);
        assert(arr[3] == INT64_MIN);
        }

        for (int pass = 0; pass < 3; ++pass)
        {
            {
            sv[2] = INT_MAX;
            sv[3] = INT_MIN;
            auto esize =  sv.extract(&arr[0], 1024, 0);
            assert(esize == 4);
            assert(arr[0] == -1+pass);
            assert(arr[1] == 1+pass);
            assert(arr[2] == INT_MAX);
            assert(arr[3] == INT_MIN);
            }

            sv[2] = INT64_MAX;
            sv[3] = INT64_MIN;

            std::vector<int64_t> target_v;
            std::vector<bm::sparse_vector<int64_t, bvect>::size_type> idx_v;
            target_v.resize(4);
            for (bm::sparse_vector<int64_t, bvect>::size_type i = 0; i < 5; ++i)
                idx_v.push_back(i);
            sv.gather(target_v.data(), idx_v.data(), 4, pass ? BM_SORTED : BM_UNSORTED);
            assert(target_v[0] == -1+pass);
            assert(target_v[1] == 1+pass);
            assert(target_v[2] == INT64_MAX);
            assert(target_v[3] == INT64_MIN);

            sv.inc(0);
            sv.inc(1);
            if (!pass)
                sv.inc(2);
            sv.inc(3);

            assert(sv.get(0) == 0+pass);
            assert(sv.get(1) == 2+pass);
            v2 = sv.get(2);
            if (!pass)
            {
                assert(sv.get(2) == 0);
            }
            else
            {
                assert(sv.get(2) == INT64_MAX);
            }
            v3 = sv.get(3);
            assert(sv.get(3) == INT64_MIN+1);

            sv.optimize();
        } // for pass
    }}

    // import back
    {{
        std::vector<int64_t> vect {0, 1, -1, INT_MIN, INT_MAX, 0, INT64_MIN+1, INT64_MIN, INT64_MAX };
        bm::sparse_vector<int64_t, bvect> sv;
        sv.import_back(vect.data(), vect.size(), false);

        assert(sv.size() == vect.size());
        for (size_t i = 0; i < vect.size(); ++i)
        {
            auto vc = vect.at(i);
            auto v = sv.at(i);
            assert(v == vc);
        } // for

        sv.resize(0);
        vect.resize(0);

        for (size_t i = 0; i < 65536*256; ++i)
            vect.push_back(int(-i));

        sv.import_back(vect.data(), vect.size(), false);

        assert(sv.size() == vect.size());
        for (size_t i = 0; i < vect.size(); ++i)
        {
            auto vc = vect.at(i);
            auto v = sv.at(i);
            assert(v == vc);
        } // for

    }}

    {{
        bm::sparse_vector<int64_t, bvect > sv;
        {
            auto bi(sv.get_back_inserter());
            *bi = (-1);
            *bi = (1);
            *bi = (INT64_MAX);
            *bi = (INT64_MIN);
            *bi = 0;
            bi.flush();
        }

        int64_t arr[1024];
        auto esize =  sv.extract(&arr[0], 1024, 0);
        assert(esize == 5);
        assert(arr[0] == -1);
        assert(arr[1] == 1);
        assert(arr[2] == INT64_MAX);
        assert(arr[3] == INT64_MIN);
        assert(arr[4] == 0);

    }}

    cout << "svector Import test..." << endl;

    {{
        sparse_vector_i32 sv;
        int arr[3] = {1,-8,3};
        sv.import(arr, 3); // import from a C-style array (fastest way to populate)
        sv.optimize();
        print_svector_stat(cout,sv);

        sparse_vector_i32::statistics st;
        sv.calc_stat(&st);
        assert(st.gap_blocks == 4);
    }}

    {{
        std::vector<int> vect;
        for (int i = 0; i < 128000; ++i)
            vect.push_back(-i);

        sparse_vector_i32 sv;
        sv.import(&vect[0], (bvect::size_type)vect.size());
        bool res = CompareSparseVector(sv, vect);
        if (!res)
        {
            cerr << "0.Bit plane import test failed" << endl;
            assert(0);exit(1);
        }
        sv.optimize();
        print_svector_stat(cout,sv);
        res = CompareSparseVector(sv, vect);
        if (!res)
        {
            cerr << "optimized Bit plane import test failed" << endl;
            assert(0);exit(1);
        }

        sparse_vector_i32 sv_1;
        std::copy(vect.begin(), vect.end(), std::back_inserter(sv_1));
        res = CompareSparseVector(sv_1, vect);
        if (!res)
        {
            cerr << "Bit plane push_back test failed" << endl;
            assert(0);exit(1);
        }

        sparse_vector_i32::statistics st;
        sv.calc_stat(&st);

        sparse_vector_i32 sv2(sv);
        res = CompareSparseVector(sv2, vect);
        if (!res)
        {
            cerr << "Bit plane copy-ctor test failed" << endl;
            assert(0);exit(1);
        }

        sv2.clear();
        sv2.import(&vect[0], (bvect::size_type)vect.size());
        res = CompareSparseVector(sv2, vect);
        if (!res)
        {
            cerr << "Bit plane copy-ctor test failed" << endl;
            assert(0);exit(1);
        }

        sparse_vector_i32 sv3;
        sv3.set(65536, 10); // set some bit to initiate it
        sv3 = sv;
        res = CompareSparseVector(sv3, vect);
        if (!res)
        {
            cerr << "Bit plane assignmnet test failed" << endl;
            exit(1);
        }

        sv3.clear();
        sv3.import(&vect[0], (bvect::size_type)vect.size());
        res = CompareSparseVector(sv3, vect);
        if (!res)
        {
            cerr << "Bit plane assignment test failed" << endl;
            exit(1);
        }
    }}


    cout << "Same value assignment test.." << endl;
    {
        bm::sparse_vector<int64_t, bvect > sv;
        const unsigned max_assign =
                            100 + bm::gap_max_bits * bm::set_sub_array_size;
        {
            auto bi(sv.get_back_inserter());
            for (unsigned i = 0; i < max_assign; ++i)
            {
                *bi = -1;
            } // for
            bi.flush();
        }
        bm::sparse_vector<int64_t, bvect>::statistics st;
        sv.optimize(0, bvect::opt_compress, &st);
        assert(st.gap_blocks == 1);
        assert(st.bit_blocks == 0);

        for (unsigned i = 0; i < max_assign; ++i)
        {
            auto v = sv[i];
            assert(v == -1);
        } // for

        for (auto it = sv.begin(); it < sv.end(); ++it)
        {
            auto v = *it;
            assert(v == -1);
        }
        sv[65536] = 0;
    }


    cout << "Linear assignment test" << endl;
    {{
    std::vector<int64_t> vect(128000);
    typedef std::vector<int64_t>::size_type vect_sz_type;
    typedef bm::sparse_vector<int64_t, bvect >::size_type sv_sz_type;
    bm::sparse_vector<int64_t, bvect > sv;
    bm::sparse_vector<int64_t, bvect > sv1(bm::use_null);
    bm::sparse_vector<int64_t, bvect > sv2;

    {
        auto bi(sv2.get_back_inserter());
        for (int i = 0; i < 128000; ++i)
        {
            vect[vect_sz_type(i)] = -i;
            sv.set(sv_sz_type(i), -i);
            sv1.set(sv_sz_type(i), -i);
            *bi = -i;
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

    cout << " -------------------------- TestSignedSparseVector() OK" << endl;
}



// ---------------------------------------------------------------------

static
void TestSparseVectorSerial()
{
    cout << "---------------------------- Test sparse vector serializer" << endl;

    // simple test gather for non-NULL vector
    {
        sparse_vector_u32 sv1;
        sparse_vector_u32 sv2;


        {
            unsigned k=0;
            for (sparse_vector_u32::size_type i = 0; i < 10; ++i, ++k)
                sv1.push_back(k+1);
        }
        BM_DECLARE_TEMP_BLOCK(tb)
        sparse_vector_serial_layout<sparse_vector_u32> sv_lay;
        bm::sparse_vector_serialize<sparse_vector_u32>(sv1, sv_lay, tb);
        const unsigned char* buf = sv_lay.buf();

        bm::sparse_vector_deserializer<sparse_vector_u32> sv_deserial;

        sparse_vector_u32::bvector_type bv_mask;
        bv_mask.set(0);
        bv_mask.set(2);
        sv_deserial.deserialize(sv2, buf, bv_mask);

        assert(sv2.size() == sv1.size());
        assert(sv2.get(0) == 1);
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

        unsigned k = 0;
        for (sparse_vector_u32::size_type i = 0; i < 100; i+=2, k+=2)
        {
            sv1[i] = k+1;
        }
        BM_DECLARE_TEMP_BLOCK(tb)
        sparse_vector_serial_layout<sparse_vector_u32> sv_lay;
        bm::sparse_vector_serialize<sparse_vector_u32>(sv1, sv_lay, tb);
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
            assert(bv_null);
            //auto cnt = bv_null->count();
            //assert(cnt == 2);

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
            assert(bv_null);
            //auto cnt = bv_null->count();
            //assert(cnt == 0);
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

        from = bm::id_max32 * 2;
        to = from + 75538;

        unsigned cnt = 0;
        for (sparse_vector_u32::size_type i = from; i < to; ++i, ++cnt)
        {
            if (cnt % 10 == 0)
                sv1.set_null(i);
            else
                sv1.set(i, cnt);
        } // for i
        BM_DECLARE_TEMP_BLOCK(tb)
        sparse_vector_serial_layout<sparse_vector_u32> sv_lay;
        bm::sparse_vector_serialize<sparse_vector_u32>(sv1, sv_lay, tb);
        const unsigned char* buf = sv_lay.buf();

        {
            bm::sparse_vector_deserializer<sparse_vector_u32> sv_deserial;

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
                    cerr << "Error: Range deserialiuzation equality failed!" << endl;
                    assert(0); exit(1);
                }

                for (auto k = i; k < j; ++k)
                {
                    auto v1 = sv1[k];
                    auto v2 = sv2[k];
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

                if (i % 0xFF == 0)
                {
                    std::cout << "\r" << j-i << flush;
                }

            } // for i
        }
        cout << "\nOK\n" << endl;
    }


    cout << "---------------------------- Test sparse vector serializer OK" << endl;
}


static
void TestSignedSparseVectorSerial()
{
    cout << "---------------------------- TestSignedSparseVectorSerial()" << endl;

    cout << "  test chain XOR serialization.. " << endl;
    bm::sparse_vector_serializer<sparse_vector_i32> sv_ser;
    bm::sparse_vector_deserializer<sparse_vector_i32> sv_deserial;

    {
        sparse_vector_serial_layout<sparse_vector_i32> sv_lay;

        sv_ser.set_xor_ref(true);

        const unsigned stride_len = 1024;

        for (unsigned pass = 2; pass < 128; ++pass)
        {
            sparse_vector_i32 sv;
            bvect::size_type from = 0;

            unsigned mask = (1u << 1);
            for (unsigned k = 0; k < pass; ++k)
            {
                sparse_vector_i32::value_type v = -int(1u | mask);
                for (unsigned i = 0; i < stride_len; i+=2)
                {
                    bvect::size_type idx = from + i;
                    sv.set(idx, v);
                } // for i
                from += stride_len;
                mask <<= 1;
                if (!mask)
                    mask = 1u << 1;
            } // for k

            sv_ser.serialize(sv, sv_lay);
            const bvect::size_type* cstat = sv_ser.get_bv_serializer().get_compression_stat();
            assert(cstat[bm::set_block_xor_chain]>=1);
            {
                const unsigned char* buf = sv_lay.buf();
                sparse_vector_i32 sv2;
                sv_deserial.deserialize(sv2, buf);

                bool eq = sv.equal(sv2);
                assert(eq);
            }

        } // for pass
    }

    sv_ser.set_xor_ref(false);

    for (unsigned pass = 0; pass < 2; ++pass)
    {
        if (!pass)
        {
            cout << "   XOR ref compression is ON" << endl;
            sv_ser.set_xor_ref(true);
        }
        else
        {
            cout << "   XOR ref compression is OFF" << endl;
            sv_ser.set_xor_ref(false);
        }

        // simple test gather for non-NULL vector
        {
            sparse_vector_i32 sv1;
            sparse_vector_i32 sv2;

            for (sparse_vector_i32::value_type i = 0; i < 10; ++i)
                sv1.push_back(0-(i + 1));
            sparse_vector_serial_layout<sparse_vector_i32> sv_lay;
            sv_ser.serialize(sv1, sv_lay);
            const unsigned char* buf = sv_lay.buf();

            sparse_vector_u32::bvector_type bv_mask;
            bv_mask.set(0);
            bv_mask.set(2);
            sv_deserial.deserialize(sv2, buf, bv_mask);

            assert(sv2.size() == sv1.size());
            assert(sv2.get(0) == -1);
            cout << sv2.get(1) << endl;
            assert(sv2.get(1) == 0);
            assert(sv2.get(2) == -3);

            sparse_vector_i32::statistics st;
            sv2.calc_stat(&st);
            assert(!st.bit_blocks);
            assert(st.gap_blocks);
        }


        // simple test gather for NULL-able vector
        {
            sparse_vector_i32 sv1(bm::use_null);
            sparse_vector_i32 sv2(sv1);

            for (sparse_vector_u32::value_type i = 0; i < 100; i += 2)
            {
                sv1[i] = -int(i + 1);
            }
            sparse_vector_serial_layout<sparse_vector_i32> sv_lay;
            sv_ser.serialize(sv1, sv_lay);
            const unsigned char* buf = sv_lay.buf();

            //bm::sparse_vector_deserializer<sparse_vector_u32> sv_deserial;
            {
                sparse_vector_i32::bvector_type bv_mask;
                bv_mask.set(0);
                bv_mask.set(2);
                bv_mask.set(1024); // out of range mask
                sv_deserial.deserialize(sv2, buf, bv_mask);


                assert(sv2.get(0) == -1);
                assert(sv2.get(1) == 0);
                assert(sv2.get(2) == -3);

                const sparse_vector_u32::bvector_type* bv_null = sv2.get_null_bvector();
                auto cnt = bv_null->count();
                auto cnt1 = sv1.get_null_bvector()->count();
                assert(cnt == 2);
                assert(cnt != cnt1);

                sparse_vector_i32::statistics st;
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
            sparse_vector_i32::size_type from, to;
            sparse_vector_i32 sv1(bm::use_null);
            sparse_vector_i32 sv2(sv1);
            sparse_vector_i32 sv3(sv1);

            from = bm::id_max32 * 3;
            to = from + 75538;

            sparse_vector_i32::size_type cnt = 0;
            for (sparse_vector_u32::size_type i = from; i < to; ++i, ++cnt)
            {
                if (cnt % 10 == 0)
                    sv1.set_null(i);
                else
                    sv1.set(i, sparse_vector_i32::value_type(cnt));
            } // for i
            sv1.sync_size();

            sparse_vector_serial_layout<sparse_vector_i32> sv_lay;

            sv_ser.serialize(sv1, sv_lay);
            const unsigned char* buf = sv_lay.buf();

            {

                //bm::sparse_vector_deserializer<sparse_vector_u32> sv_deserial;
                sparse_vector_i32 sv4(bm::use_null);
                sv_deserial.deserialize(sv4, buf);
                {
                    sparse_vector_u32::size_type midx;
                    bool is_eq = sv1.equal(sv4);
                    if (!is_eq)
                    {
                        bool b = bm::sparse_vector_find_first_mismatch(sv1, sv4, midx);
                        assert(b);
                        cout << "Mismatch at pos = " << midx << endl;
                        auto v1 = sv1[midx];
                        auto v4 = sv4[midx];
                        cout << "v1=" << v1 << " v4=" << v4 << " xor(v1, v4)=" << (v1 ^ v4) << endl;
                    }
                    assert(is_eq);
                }


                auto i = from;
                auto j = to;
                bool is_eq;
                sparse_vector_i32::size_type pos;
                bool found;
                //i = 59982;
                //i = 2147491602; j = 2147551230;

                for (i = from; i < j; ++i, --j)
                {
                    sparse_vector_i32::bvector_type bv_mask;
                    bv_mask.set_range(i, j);

                    sparse_vector_i32 sv_filt(sv1);
                    sv_filt.filter(bv_mask);
                    sparse_vector_i32 sv_range(bm::use_null);
                    sv_range.copy_range(sv1, i, j);

                    is_eq = sv_filt.equal(sv_range);
                    assert(is_eq);

                    sv_deserial.deserialize(sv2, buf, bv_mask);

                    assert(sv2.size() == sv1.size());
                    is_eq = sv2.equal(sv_range);
                    if (!is_eq)
                    {
                        found = bm::sparse_vector_find_first_mismatch(sv2, sv_range, pos, bm::no_null);
                        if (found)
                        {
                            auto vf = sv_filt.get(pos);
                            auto v2 = sv2.get(pos);
                            auto v3 = sv_range.get(pos);

                            cerr << "Mismatch at:" << pos << endl;
                            cerr << vf << "!=" << v2 << "!=" << v3 << endl;
                            cerr << "[i, j] = " << i << ":" << j << endl;

                        }
                        assert(is_eq);
                    }

                    sv_deserial.deserialize(sv3, buf, i, j);
                    {
                        bool b = bm::sparse_vector_find_first_mismatch(sv3, sv_range, pos);
                        if (b)
                        {
                            auto vf = sv_filt.get(pos);
                            auto v2 = sv3.get(pos);
                            auto v3 = sv_range.get(pos);
                            auto xd = v2 ^ v3;

                            cerr << "Mismatch at:" << pos << " xor diff=" << xd << endl;
                            cerr << vf << "!=" << v2 << "!=" << v3 << endl;
                            cerr << "[i, j] = " << i << ":" << j << endl;
                        }
                        assert(!b);
                    }
                    is_eq = sv2.equal(sv3);
                    if (!is_eq)
                    {
                        cerr << "Error: Range deserialization equality failed!" << endl;
                        assert(0); exit(1);
                    }

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

    } // for pass

    cout << "---------------------------- TestSignedSparseVectorSerial()" << endl;
}


// ---------------------------------------------------------------------

typedef bm::sparse_vector<unsigned, bvect > sparse_vector_u32;
typedef bm::sparse_vector<unsigned long long, bvect > sparse_vector_u64;
typedef bm::rsc_sparse_vector<unsigned, sparse_vector_u32> rsc_sparse_vector_u32;
typedef bm::rsc_sparse_vector<unsigned, sparse_vector_u64> rsc_sparse_vector_u64;

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

    cout << "64-bit sparse bulk inserter test (1)" << endl;
    {{
        bm::sparse_vector<unsigned long long, bvect > sv(bm::use_null);
        ref_vect vect;
        generate_vect_simpl0(vect);
        bulk_load_SV_set_ref(&sv, vect);
        compare_SV_set_ref(sv, vect);
        sv.optimize();
        cout << "ok" << endl;
        compare_SV_set_ref(sv, vect);
    }}
    
    cout << "64-bit sparse bulk inserter test (2)" << endl;
    {{
        bm::sparse_vector<unsigned long long, bvect > sv(bm::use_null);
        ref_vect vect;
        generate_vect48(vect);
        bulk_load_SV_set_ref(&sv, vect);
        compare_SV_set_ref(sv, vect);
        sv.optimize();
        cout << "ok" << endl;
        compare_SV_set_ref(sv, vect);
    }}

    cout << "---------------------------- Bit-plain sparse vector inserter OK" << endl;
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

    }

    cout << "\nRSC gather stress test ..." << endl;
    {
        rsc_sparse_vector_u32 csv1;
        rsc_sparse_vector_u32 csv2;

        rsc_sparse_vector_u32::size_type from = bm::id_max32 * 2;
        rsc_sparse_vector_u32::size_type to = from + 75538;

        {
            rsc_sparse_vector_u32::back_insert_iterator rs_bi = csv1.get_back_inserter();
            rs_bi.add_null();
            rs_bi.add(1);
            rs_bi.add(2);
            rs_bi.add_null();
            rs_bi.add(4);
            rs_bi.add_null(from); // add many NULLs

            unsigned k = 0;
            for (auto i = from; i < to; ++i, ++k)
            {
                rs_bi.add(k);
                rs_bi.add_null();
            }
            rs_bi.flush();
        }

        BM_DECLARE_TEMP_BLOCK(tb)
        sparse_vector_serial_layout<rsc_sparse_vector_u32> sv_lay;
        bm::sparse_vector_serialize<rsc_sparse_vector_u32>(csv1, sv_lay, tb);
        const unsigned char* buf = sv_lay.buf();

        auto j = to;
        for (auto i = from; i < j; ++i, --j)
        {
            sparse_vector_u32::bvector_type bv_mask;
            bv_mask.set_range(i, j);
            bm::sparse_vector_deserializer<rsc_sparse_vector_u32> sv_deserial;
            sv_deserial.deserialize(csv2, buf, bv_mask);
            csv2.sync();

            for (auto i0 = i; i0 < j; ++i0)
            {
                auto v1 = csv1[i0];
                auto v2 = csv2[i0];
                assert(v1 == v2);
                assert(csv1.is_null(i0) == csv2.is_null(i0));
            } // for

            cout << "\r" << (j-i) << std::flush;

        } // for i


    }
    cout << "\nOK" << endl;


    cout << " ------------------------------ TestCompressSparseVectorSerial() OK" << endl;
}


static
void TestSparseVectorSerialization2()
{
    cout << " ------------------------------ TestSparseVectorSerialization2()" << endl;
    const unsigned int BSIZE = 250000000;

    const unsigned char* buf;
    bool eq;
    size_t sz1, sz2;


    bm::sparse_vector_serializer<sparse_vector_u32> sv_serializer;
    bm::sparse_vector_deserializer<sparse_vector_u32> sv_deserial;

    cout << "Test data-frame XOR compression" << endl;
    {
        sparse_vector_u32 sv1i, sv2i, sv3i(bm::use_null);
        sparse_vector_u32 sv1o, sv2o, sv3o(bm::use_null);

        bm::sparse_vector_serial_layout<sparse_vector_u32> sv_lay1, sv_lay2, sv_lay3;


        for (unsigned i = 0; i < 65536; i+=2)
        {
            sv1i[i] = 4;
            sv2i[i] = 8;
            sv3i[i] = 0;
        }

        bm::sparse_vector_serializer<sparse_vector_u32>::bv_ref_vector_type bv_ref;
        // add references in reverse(!) order
        bv_ref.add_vectors(sv3i.get_bmatrix());
        bv_ref.add_vectors(sv2i.get_bmatrix());
        bv_ref.add_vectors(sv1i.get_bmatrix());

        assert(bv_ref.size() == 3);

        sv_serializer.set_xor_ref(&bv_ref);
        assert(sv_serializer.is_xor_ref());

        bm::sparse_vector_serializer<sparse_vector_u32>::xor_sim_model_type sim_model;
        xor_sim_params xs_params;
        sv_serializer.compute_sim_model(sim_model, bv_ref, xs_params);
        sv_serializer.set_sim_model(&sim_model);


        sv_serializer.serialize(sv1i, sv_lay1);
        {
            const bvect::size_type* cstat = sv_serializer.get_bv_serializer().get_compression_stat();
            assert(cstat[bm::set_block_ref_eq]==0);
        }
        sv_serializer.serialize(sv2i, sv_lay2);
        {
            const bvect::size_type* cstat = sv_serializer.get_bv_serializer().get_compression_stat();
            assert(cstat[bm::set_block_ref_eq]==1);
        }
        sv_serializer.serialize(sv3i, sv_lay3);
        {
            const bvect::size_type* cstat = sv_serializer.get_bv_serializer().get_compression_stat();
            assert(cstat[bm::set_block_ref_eq]==1);
        }

        // ----------


        bm::sparse_vector_deserializer<sparse_vector_u32>::bv_ref_vector_type bv_ref_d;

        buf = sv_lay1.buf();
        sz2 = sv_lay1.size();


        sv_deserial.deserialize_structure(sv1o, sv_lay1.buf());
        sv_deserial.deserialize_structure(sv2o, sv_lay2.buf());
        sv_deserial.deserialize_structure(sv3o, sv_lay3.buf());

        bv_ref_d.add_vectors(sv3o.get_bmatrix());
        bv_ref_d.add_vectors(sv2o.get_bmatrix());
        bv_ref_d.add_vectors(sv1o.get_bmatrix());

        sv_deserial.set_xor_ref(&bv_ref_d);

        sv_deserial.deserialize(sv1o, buf, false);
        eq = sv1i.equal(sv1o);
        assert(eq);

        buf = sv_lay2.buf();
        sz2 = sv_lay2.size();

        sv_deserial.deserialize(sv2o, buf, false);
        eq = sv2i.equal(sv2o);
        assert(eq);

        buf = sv_lay3.buf();
        sz2 = sv_lay3.size();

        sv_deserial.deserialize(sv3o, buf, false);
        eq = sv3i.equal(sv3o);
        assert(eq);


        sv_deserial.set_xor_ref(0); // unset
    }
    cout << "Test data-frame XOR compression - OK" << endl;

    // -------------------------------------------------


    sparse_vector_u32 sv1(bm::use_null);
    sparse_vector_u32 sv2(bm::use_null);
    sparse_vector_u32 sv3(bm::use_null);

    generate_serialization_test_set(sv1, BSIZE);

    bm::sparse_vector_serial_layout<sparse_vector_u32> sv_lay;


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


// ------------------------------------------------------------------------

static
void CheckSparseVectorGather(const sparse_vector_u32& sv,
                             sparse_vector_u32::size_type from,
                             sparse_vector_u32::size_type to,
                             unsigned control_value = 0)
{
    assert(sv.size());
    assert (to >= from);
    
    sparse_vector_u32::size_type gather_size = to - from + 1;
    std::vector<unsigned> target_v;
    std::vector<unsigned> target_v_control;
    std::vector<sparse_vector_u32::size_type> idx_v;
    target_v.resize(gather_size);
    target_v_control.resize(gather_size);
    idx_v.reserve(gather_size);
    
    for (sparse_vector_u32::size_type i = from; i <= to; ++i)
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
                    assert(0);exit(1);
                }
            }
            else
            {
                v1 = sv.get(i);
                if (v1 != v2)
                {
                    cerr << "Error! gather mismatch " << v1 << " " << v2
                         << " at=" << i << endl;
                    assert(0);exit(1);
                }
            }
        }
    } // for
#endif
}

static
void CheckSparseVectorGatherRandom(const sparse_vector_u32& sv,
                                   sparse_vector_u32::size_type gather_size)
{
    assert(sv.size());
 
    std::random_device rd;
    std::mt19937_64 mt_rand(rd());

    if (gather_size == 0)
        gather_size = 1;
    
    std::vector<unsigned> target_v;
    std::vector<sparse_vector_u32::size_type> idx_v;
    target_v.resize(gather_size);
    idx_v.reserve(gather_size);
    
    for (sparse_vector_u32::size_type i = 0; i < gather_size; ++i)
    {
        sparse_vector_u32::size_type r_idx = mt_rand() % (sv.size()-1);
        idx_v.push_back(r_idx);
    }
    
    sv.gather(target_v.data(), idx_v.data(), gather_size, BM_UNSORTED);

    unsigned k = 0;
    for (sparse_vector_u32::size_type i = 0; i < gather_size; ++i, ++k)
    {
        unsigned v1 = sv.get(idx_v[k]);
        unsigned v2 = target_v[k];
        if (v1 != v2)
        {
            {
                cerr << "Error! random gather mismatch " << v1 << " " << v2
                     << " at=" << i << endl;
                assert(0);
                exit(1);
            }
        }
    } // for
}




static
void TestSparseVectorGatherDecode()
{
    cout << "---------------------------- Test sparse vector gather decode" << endl;
    
    
    {
        unsigned base = 0;
        sparse_vector_u32 sv0;
        sv0[base + 0] = 0;
        sv0[base + 1] = 1;
        sv0[base + 2] = 1;
        sv0[base + 3] = 1;

        for (unsigned pass = 0; pass < 2; ++pass)
        {
            unsigned d[32] = { 25, };

            auto sz = sv0.decode(&d[0], base + 1, 2);
            assert(sz == 2);
            assert(d[0] == 1);
            assert(d[1] == 1);
            sz = sv0.decode(&d[0], base + 2, 2);
            assert(sz == 2);
            assert(d[0] == 1);
            assert(d[1] == 1);

            sv0.optimize();
        }
    }


    cout << "64-bit sparse gather test (1)" << endl;
    {{
        bm::sparse_vector<unsigned long long, bvect > sv(bm::use_null);
        ref_vect vect;
        generate_vect_simpl0(vect);
        load_SV_set_ref(&sv, vect);
        compare_SV_set_ref(sv, vect);

        for (unsigned p = 0; p < 3; ++p)
        {
            cout << "Pass " << p << endl;
            {
                ref_vect vect_d; vect_d.resize(vect.size());
                for (size_t k = 0; k < vect.size(); ++k)
                {
                    auto sz = sv.gather(&vect_d[0], &vect[k], 1, BM_UNSORTED);
                    assert(sz == 1);
                    bvect::size_type idx = vect[k];
                    unsigned long long v = sv[idx];
                    assert(vect_d[0] == v);
                    vect_d[0] = 0;
                    sz = sv.gather(&vect_d[0], &vect[k], 1, BM_SORTED);
                    assert(sz == 1);
                    idx = vect[k];
                    v = sv[idx];
                    assert(vect_d[0] == v);
                    vect_d[0] = 0;
                }
            }
            {
                ref_vect vect_d; vect_d.resize(vect.size());
                auto sz = sv.gather(&vect_d[0], &vect[0], vect.size(), BM_UNSORTED);
                assert(sz == vect.size());
                compare_SV_set_ref(sv, vect_d);
            }
            {
                ref_vect vect_d; vect_d.resize(vect.size());
                auto sz = sv.gather(&vect_d[0], &vect[0], vect.size(), BM_UNKNOWN);
                assert(sz == vect.size());
                compare_SV_set_ref(sv, vect_d);
            }
            {
                ref_vect vect_d; vect_d.resize(vect.size());
                auto sz = sv.gather(&vect_d[0], &vect[0], vect.size(), BM_SORTED);
                assert(sz == vect.size());
                compare_SV_set_ref(sv, vect_d);
            }
            sv.optimize();
        } // for k
        cout << "ok" << endl;
    }}
    
    cout << "64-bit sparse decode test" << endl;
    {{
            bm::sparse_vector<unsigned long long, bvect > sv(bm::use_null);
            bvect::size_type from = bm::id_max - 2048;
            bvect::size_type to = bm::id_max - 1;
            for (auto i = from; i < to; ++i)
            {
                sv[i] = i;
            } // for

            std::vector<unsigned long long> vect;
            vect.resize(to - from + 1);
            bvect::size_type target_sz = 1280;
            for (unsigned pass = 0; pass < 2; ++pass)
            {
                cout << "Pass " << pass << endl;
                for (bvect::size_type k = 1; k < target_sz; ++k)
                {
                    cout << "\r" << k << "/" << target_sz << flush;
                    sv.decode(vect.data(), from, k, true);
                    bvect::size_type fr = from;
                    for (bvect::size_type i = 0; i < k; ++i, ++fr)
                    {
                        auto v = vect[i];
                        assert(v == fr);
                    }
                } // for
                sv.optimize();
                cout << endl;
            } // for pass
    }}



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
            unsigned depth = (unsigned)rand() % 30000;
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
            unsigned gsize = (unsigned)rand()%2024;
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

// ---------------------------------------------------------------------------

template<class SV>
void bvector_transform_11(typename SV::bvector_type& bvect_in,
                          const    SV&               sv_brel,
                          typename SV::bvector_type& bvect_out)
{
    bm::set2set_11_transform<SV> bin_trans;
    bin_trans.run(bvect_in, sv_brel, bvect_out);
}


static
void TestSparseVectorTransform()
{
    cout << " ---------------- Test set transformation with sparse vector" << endl;

    {
        sparse_vector_u64 sv(bm::use_null);
        bvect bv_in { 1, 2, 3, 10, bm::id_max-1 };
        bvect bv_out;
        
        bvector_transform_11(bv_in, sv, bv_out);
        assert(bv_out.count() == 0);
        cout << "Transform11 with empty sv - ok" << endl;

        bm::set2set_11_transform<sparse_vector_u64> set2set;
        bvect::size_type to;
        bool found = set2set.remap(0, sv, to);
        assert(!found);
        found = set2set.remap(3, sv, to);
        assert(!found);
    }

    {
        sparse_vector_u64 sv(bm::use_null);

        sv.set(2, 25);
        sv.set(3, 35);
        sv.set(7, 75);
        sv.set(10, 2);
        sv.set(21, 201);
        sv.set(bm::id_max-1, bm::id_max-2);

        bm::set2set_11_transform<sparse_vector_u64> set2set;
        bvect::size_type to;
        bool found = set2set.remap(0, sv, to);
        assert(!found);
        found = set2set.remap(3, sv, to);
        assert(found);
        assert(to == 35);
        found = set2set.remap(bm::id_max-1, sv, to);
        assert(found);
        assert(to == bm::id_max-2);

        bvect bv_in { 1, 2, 3, 10, 20, bm::id_max-1};
        bvect bv_control {25, 35, 2,  bm::id_max-2};

        {
            bvect bv_out;
            bvector_transform_11(bv_in, sv, bv_out);
            int cmp = bv_control.compare(bv_out);
            if (cmp != 0)
            {
                cerr << "Transform11 (1) control comparison failed" << endl;
                assert(0); exit(1);
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
        sparse_vector_u64 sv;

        sv.set(2, 25);
        sv.set(3, 35);
        sv.set(7, 75);
        sv.set(10, 2);
        sv.set(21, 201);

        bm::set2set_11_transform<sparse_vector_u64> set2set;
        bvect::size_type to;
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
        sparse_vector_u64 sv(bm::use_null);
        
        generate_bvector(bv_in, bm::id_max32/4, false);
        
        {
            bvect::enumerator en = bv_in.first();
            for (;en.valid(); ++en)
            {
                auto idx = *en;
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
        sparse_vector_u64 sv(bm::use_null);
        
        generate_bvector(bv_in, bm::id_max32/4, false);
        
        bvect bv_control;
        {
            bvect::enumerator en = bv_in.first();
            for (;en.valid(); ++en)
            {
                auto idx = *en;
                bv_control.set(idx + bm::id_max32);
            }
        }
        
        {
            bvect::enumerator en = bv_in.first();
            for (;en.valid(); ++en)
            {
                auto idx = *en;
                sv.set(idx, idx + bm::id_max32); // 1 to 1 direct with a base shift
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
        sparse_vector_u64 sv(bm::use_null);
        
        generate_bvector(bv_in, bm::id_max32/4, false);

        bvect bv_control;
        bv_control.set(bm::id_max32-2);
        
        {
            bvect::enumerator en = bv_in.first();
            for (;en.valid(); ++en)
            {
                auto idx = *en;
                sv.set(idx, bm::id_max32-2); // M:1
            }
        }
        bvector_transform_11(bv_in, sv, bv_out);
        
        int cmp = bv_control.compare(bv_out);
        if (cmp != 0)
        {
            cerr << "Transform11 (4) control comparison failed" << endl;
            assert(0); exit(1);
        }
        
        sv.optimize();
        bvector_transform_11(bv_in, sv, bv_out);

        cmp = bv_control.compare(bv_out);
        if (cmp != 0)
        {
            cerr << "Transform11 (4, 2) control comparison failed" << endl;
            assert(0); exit(1);
        }
        cout << "Transform11 (4) - ok" << endl;
    }
    
    
    {
        bvect bv_in, bv_out;
        sparse_vector_u64 sv(bm::use_null);
        
        generate_bvector(bv_in, bm::id_max32/4, false);
        bvect bv_control(bv_in);
        {
            bvect::enumerator en = bv_in.first();
            for (;en.valid(); ++en)
            {
                auto idx = *en;
                sv.set(idx, idx); // same:same
            }
        }
        bvector_transform_11(bv_in, sv, bv_out);
        
        int cmp = bv_control.compare(bv_out);
        if (cmp != 0)
        {
            cerr << "Transform11 (5) control comparison failed" << endl;
            assert(0); exit(1);
        }
        
        sv.optimize();
        bvector_transform_11(bv_in, sv, bv_out);

        cmp = bv_control.compare(bv_out);
        if (cmp != 0)
        {
            cerr << "Transform11 (5, 2) control comparison failed" << endl;
            assert(0); exit(1);
        }
        cout << "Transform11 (5) - ok" << endl;
    }


    cout << " --------------- Test set transformation with sparse vector OK" << endl;
}

// -------------------------------------------------------------------------


template<typename SV>
void CheckSparseVectorRange(const SV& sv,
                            typename SV::size_type left,
                            typename SV::size_type right)
{
    SV sv1(bm::use_null);
    SV sv2(sv);
    sv1.copy_range(sv, left, right);
    
    if (right >= sv.size())
    {
        right = sv.size()-1;
    }
    
    if (left == right)
    {
        auto v1 = sv.get(left);
        auto v2 = sv1[right];
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
        for (bvect::size_type i = left; i <= right; ++i)
        {
            auto v1 = sv.get(i);
            auto v2 = sv1[i];
            if (v1 != v2)
            {
                cerr << "Error! Copy range check failed at:" << i << endl;
                exit(1);
            }
        } // for
        cerr << "detailed check did not find issues. error in test?" << endl;
        assert(0); exit(1);
    }
}

static
void TestSparseVectorRange()
{
    cout << " ---------------- Sparse vector Range partitioning test" << endl;

    cout << "Basic check" << endl;
    {
        sparse_vector_u64 sv(bm::use_null);
        sv.set(2, 25);
        sv.set(3, 35);
        sv.set(7, 75);
        sv.set(10, 2);
        sv.set(21, 201);
        sv.set(bm::id_max32*2, bm::id_max32*10);
        sv.set(bm::id_max32*2+1, bm::id_max32*5);
        sv.set(bm::id_max32*2+3, bm::id_max32);
        sv.set(bm::id_max-2, 25);
        sv.set(bm::id_max-3, 35);
        sv.set(bm::id_max-7, 75);

        CheckSparseVectorRange(sv, 0, 0);
        CheckSparseVectorRange(sv, 2, 2);
        CheckSparseVectorRange(sv, 7, 10);
        CheckSparseVectorRange(sv, bm::id_max32*2, bm::id_max32*2+10);
        CheckSparseVectorRange(sv, bm::id_max-7, bm::id_max-2);
    }

    cout << "Stress check 1 (constant)" << endl;
    {
        sparse_vector_u64 sv(bm::use_null);
        const bvect::size_type sv_from = bm::id_max - 1200001;
        const bvect::size_type sv_max = 120000;
        cout << "Filling the vector" << endl;
        sv.push_back_null(sv_from);
        for (bvect::size_type i = 0; i < sv_max; ++i)
        {
            sv.push_back(9);
        }
        
        cout << "Phase 1.." << endl;
        for (bvect::size_type i = sv_from; i < sv_max; ++i)
        {
            CheckSparseVectorRange(sv, sv_from+0, i);
            CheckSparseVectorRange(sv, sv_from+i, sv_max+10);
            cout << "\r" << i << "/" << sv_max << flush;
        }
        cout << endl;
        
        cout << "\nPhase 2.." << endl;
        bvect::size_type k = sv_max;
        for (bvect::size_type i = sv_from; i < k; ++i, --k)
        {
            CheckSparseVectorRange(sv, i, k);
        }
        
        sv.optimize();
        
        cout << "Phase 3.." << endl;
        for (bvect::size_type i = sv_from; i < sv_max; ++i)
        {
            CheckSparseVectorRange(sv, sv_from+0, i);
            CheckSparseVectorRange(sv, i, sv_max+10);
        }
        
        cout << "Phase 4.." << endl;
        k = sv_max;
        for (bvect::size_type i = sv_from; i < k; ++i, --k)
        {
            CheckSparseVectorRange(sv, i, k);
        }
    }

    cout << "\nStress check 2 (liner function)" << endl;
    {
        sparse_vector_u64 sv(bm::use_null);
        const bvect::size_type sv_max = 250000;
        cout << "Filling the vector" << endl;
        for (bvect::size_type i = 0; i < sv_max; ++i)
        {
            sv.push_back(i);
        }
        
        cout << "Phase 2-1.." << endl;
        for (bvect::size_type i = 0; i < sv_max; ++i)
        {
            CheckSparseVectorRange(sv, i, i);
            CheckSparseVectorRange(sv, 0, i);
            CheckSparseVectorRange(sv, i, sv_max+10);
            if (i % 256 == 0)
                cout << "\r" << i << "/" << sv_max << flush;
        }
        
        cout << "\nPhase 2-2.." << endl;
        bvect::size_type k = sv_max;
        for (bvect::size_type i = 0; i < k; ++i, --k)
        {
            CheckSparseVectorRange(sv, i, k);
        }
    }
    
    cout << " ---------------- Sparse vector Range partitioning test  OK\n" << endl;
}

// ----------------------------------------------------------------------


template <typename SV>
void CheckSparseVectorFilter(const SV& sv, typename SV::size_type factor)
{
    SV sv1(sv);
    
    typename SV::bvector_type bv_mask;
    for (typename SV::size_type i = 0; i < sv.size(); ++i)
    {
        if (i % factor == 0)
            bv_mask.set(i);
    }
    
    sv1.filter(bv_mask);
    
    for (typename SV::size_type i = 0; i < sv.size(); ++i)
    {
        auto v = sv.get(i);
        bool is_null = sv.is_null(i);
        auto v1 = sv1.get(i);
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
        sv.set(bm::id_max32+5, 301);

        sparse_vector_u32::bvector_type bv_mask { 2, 7, 256, bm::id_max32+5 };
        
        sv.filter(bv_mask);
        
        for (bvect::size_type i = 0; i < sv.size(); ++i)
        {
            auto v = sv.get(i);
            bool is_null = sv.is_null(i);
            if (i == 2 || i == 7 || i == bm::id_max32+5)
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
        sparse_vector_u64 sv(bm::use_null);
        const bvect::size_type sv_max = 250000;
        cout << "Filling the vector ... " << flush;
        for (bvect::size_type i = 0; i < sv_max; ++i)
        {
            sv.push_back(i);
        }
        sv[0] = 113213;

        
        sparse_vector_u64::bvector_type bv_mask;
        for (bvect::size_type i = 0; i < sv_max; ++i)
        {
            if (i % 2 == 0)
                bv_mask.set(i);
        }
        cout << "done." << endl;
        
        sv.filter(bv_mask);
        for (bvect::size_type i = 0; i < sv_max; ++i)
        {
            auto v = sv.get(i);
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
        sparse_vector_u64 sv(bm::use_null);
        const bvect::size_type sv_max = 250000;
        cout << "Filling the vector ... " << flush;
        for (bvect::size_type i = 0; i < sv_max; ++i)
        {
            sv.push_back(i);
        }
        cout << "done" << endl;
        const bvect::size_type max_factor = 10000;
        for (bvect::size_type i = 2; i < max_factor; ++i)
        {
            CheckSparseVectorFilter(sv, i);
            if ((i & 0xFF) == 0)
                cout << "\r" << i << "/" << max_factor << flush;
        }
        cout << endl;
    }
    
    cout << " ---------------- Sparse vector Filter test OK" << endl;
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
        auto found = bv_control.count();
        assert(found == 20);
        scanner.invert(sv, bv_control);
        found = bv_control.count();
        assert(!found);
    }
    
    {
        sparse_vector_u32 sv(bm::use_null);
        sv[bm::id_max-1] = 1;
        for(unsigned k = 0; k < 2; ++k)
        {
            bvect::size_type pos;
            bool found = scanner.find_eq(sv, 1, pos);
            assert(found);
            assert(pos == bm::id_max-1);
            sv.optimize();
        }
    }

    {
        cout << endl << "Unique search check" << endl;
        sparse_vector_u32 sv(bm::use_null);
        rsc_sparse_vector_u32 csv(bm::use_null);

        bvect bv_control, bv_control2;
        bvect::allocator_pool_type pool;
        bvect::mem_pool_guard(pool, bv_control);
        bvect::mem_pool_guard(pool, bv_control2);

        cout << "Loading sparse vectors..." << flush;
        unsigned sv_size = 1256;
        bvect::size_type bv_from = bm::id_max - sv_size - 1;
        {
            sparse_vector_u32::back_insert_iterator bi(sv.get_back_inserter());
            bi.add_null(bv_from);
            for (unsigned j = 0; j < sv_size; ++j)
            {
                *bi = j;
            }
        }
        csv.load_from(sv);
        cout << "ok." << endl;

        {
            chrono_taker<std::ostream> ct(cout, "sparse_vector<> search");
            for (unsigned j = 0; j < sv_size; ++j)
            {
                scanner.find_eq(sv, j, bv_control);
                if (bv_control.count() != 1)
                {
                    cerr << "1. Unique search discrepancy at value=" << j
                        << " count = " << bv_control.count() << endl;
                    assert(0);  exit(1);
                }
                bvect::size_type v1, v2;
                bool b = bv_control.find_range(v1, v2);
                assert(b);
                if (v1 != v2)
                {
                    cerr << "2. Unique search discrepancy at value=" << j
                        << " count = " << bv_control.count() << endl;
                    assert(0); exit(1);
                }

                bvect::size_type pos;
                bool found = scanner.find_eq(sv, j, pos);
                if (!found)
                {
                    cerr << "3. Unique search failure at value=" << j
                        << endl;
                    assert(0);  exit(1);
                }
                if (v1 != pos)
                {
                    cerr << "4. Unique search discrepancy at value=" << j
                        << " found = " << pos << endl;
                    assert(0); exit(1);
                }
                assert(pos == bv_from + j);
 

                rsc_scanner.find_eq(csv, j, bv_control2);
                int res = bv_control.compare(bv_control2);
                if (res != 0)
                {
                    cerr << "RSC scan comparison failed at value =" << j
                        << endl;
                    assert(0);  exit(1);
                }

                cout << "\r" << j << "/" << sv_size << "    " << flush;
            } // for
            cout << endl;
        }
        cout << "Unique search OK" << endl;
    }

    {
        cout << "Find EQ test on flat data" << endl;
        bvect::allocator_pool_type pool;
        unsigned max_value = 1280;
        for (unsigned value = 0; value < max_value; ++value)
        {
            sparse_vector_u32 sv(bm::use_null);
            sparse_vector_u64 sv_64(bm::use_null);
            rsc_sparse_vector_u32 csv;

            bvect bv_control, bv_control2;
            bvect::mem_pool_guard(pool, bv_control);

            unsigned sv_size = 67000;
            bvect::size_type sv_from = bm::id_max - sv_size - 10;

            {
                sparse_vector_u32::back_insert_iterator bi(sv.get_back_inserter());
                sparse_vector_u64::back_insert_iterator bi_64(sv_64.get_back_inserter());
                bi.add_null(sv_from);
                bi_64.add_null(sv_from);
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
            auto found = bv_control.count();
            if (found != sv_size)
            {
                cerr << "1. sparse_vector<>::find_eq() discrepancy for value=" << value
                    << " count = " << found << endl;
                assert(0); exit(1);
            }
            auto first = bv_control.get_first();
            assert(first == sv_from);

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
                    assert(0); exit(1);
                }
                first = bv_control.get_first();
                assert(first == sv_from);
            }

            // not found check
            scanner.find_eq(sv, value + 1, bv_control);
            if (bv_control.any())
            {
                cerr << "1. sparse_vector<>::find_eq() (any) discrepancy for value=" << value + 1
                    << " count = " << bv_control.count() << endl;
                assert(0); exit(1);
            }
            rsc_scanner.find_eq(csv, value + 1, bv_control2);
            res = bv_control.compare(bv_control2);
            if (res != 0)
            {
                cerr << "1. RSC scan comparison failed at value =" << value + 1
                    << endl;
                assert(0); exit(1);
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
                assert(0); exit(1);
            }
            first = bv_control.get_first();
            assert(first == sv_from);


            // not found check
            scanner.find_eq(sv, value + 1, bv_control);
            if (bv_control.any())
            {
                cerr << "2. sparse_vector<>::find_eq() (any) discrepancy for value=" << value + 1
                    << " count = " << bv_control.count() << endl;
                assert(0); exit(1);
            }

            cout << "\r" << value << "/" << max_value << "    " << flush;
        }

        cout << endl << "Flat EQ ok" << endl;
    }



    cout << " \n--------------- Test sparse_vector<> scan algo OK" << endl;
}

static
void TestSignedSparseVectorScan()
{
    cout << " --------------- TestSignedSparseVectorScan()" << endl;

    bm::sparse_vector_scanner<sparse_vector_i32> scanner;
    bm::sparse_vector_scanner<sparse_vector_i64> scanner_64;
    bm::sparse_vector_scanner<rsc_sparse_vector_i32> rsc_scanner;

    {
        sparse_vector_i32 sv(bm::use_null);
        bvect bv_control;
        scanner.find_eq(sv, 25, bv_control);
        assert(!bv_control.any());
        scanner.invert(sv, bv_control);
        assert(!bv_control.any());
    }

    {
        sparse_vector_i32 sv;
        bvect bv_control;
        for (unsigned i = 0; i < 20; ++i)
        {
            sv.set(i, 0);
        }
        scanner.find_eq(sv, 0, bv_control);
        bvect::size_type found = bv_control.count();
        assert(found == 20);
        scanner.invert(sv, bv_control);
        found = bv_control.count();
        assert(!found);
    }

    {
        cout << endl << "Unique search check" << endl;
        sparse_vector_i32 sv(bm::use_null);
        rsc_sparse_vector_i32 csv(bm::use_null);

        bvect bv_control, bv_control2;
        bvect::allocator_pool_type pool;
        bvect::mem_pool_guard g1(pool, bv_control);
        bvect::mem_pool_guard g2(pool, bv_control2);

        unsigned sv_size = 1256000;
        {
            sparse_vector_i32::back_insert_iterator bi(sv.get_back_inserter());
            for (unsigned j = 0; j < sv_size; ++j)
            {
                if (j & 1)
                    *bi = int(j);
                else
                    *bi = -int(j);
            }
        }
        csv.load_from(sv);

        {
        chrono_taker<std::ostream> ct(cout, "sparse_vector<> search");

            for (unsigned j = 0; j < sv_size; ++j)
            {
                int search_value;
                if (j & 1)
                    search_value = int(j);
                else
                    search_value = -int(j);
                scanner.find_eq(sv, search_value, bv_control);

                if (bv_control.count()!= 1)
                {
                    cerr << "1. Unique search discrepancy at value=" << j
                         << " count = " << bv_control.count() << endl;
                    assert(0);exit(1);
                }
                bvect::size_type v1, v2;
                bool b = bv_control.find_range(v1, v2);
                assert(b);
                if (v1 != v2)
                {
                    cerr << "2. Unique search discrepancy at value=" << j
                         << " count = " << bv_control.count() << endl;
                    exit(1);
                }

                bvect::size_type pos;
                bool found = scanner.find_eq(sv, search_value, pos);
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

                rsc_scanner.find_eq(csv, search_value, bv_control2);
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
        int max_value = 128000;
        for (int value = 0; value < max_value; ++value)
        {
            sparse_vector_i32 sv(bm::use_null);
            sparse_vector_i64 sv_64(bm::use_null);
            rsc_sparse_vector_i32 csv;

            bvect bv_control, bv_control2;
            bvect::mem_pool_guard g0(pool, bv_control);

            unsigned sv_size = 67000;

            {
                sparse_vector_i32::back_insert_iterator bi(sv.get_back_inserter());
                sparse_vector_i64::back_insert_iterator bi_64(sv_64.get_back_inserter());
                for (unsigned j = 0; j < 67000; ++j)
                {
                    *bi = -value;
                    bm::id64_t v64 = (unsigned)value;
                    v64 <<= 32;
                    *bi_64 = -(signed long long)v64;
                }
            }
            csv.load_from(sv);

            scanner.find_eq(sv, -value, bv_control);
            bvect::size_type found = bv_control.count();

            if (found != sv_size)
            {
                cerr << "1. sparse_vector<>::find_eq() discrepancy for value=" << value
                     << " count = " << found << endl;
                exit(1);
            }

            rsc_scanner.find_eq(csv, -value, bv_control2);
            int res = bv_control.compare(bv_control2);
            if (res != 0)
            {
                cerr << "RSC scan comparison failed at value =" << value
                << endl;
                exit(1);
            }


            {
                bm::id64_t v64 = unsigned(value);
                v64 <<= 32;
                signed long long v64s = -(signed long long)v64;
                scanner_64.find_eq(sv_64, v64s, bv_control);
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

            scanner.find_eq(sv, -value, bv_control);
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



    cout << " \n--------------- TestSignedSparseVectorScan() OK" << endl;
}




//-----------------------------------------------------------------------------------


static
void TestCompressSparseVector()
{
    cout << " ------------------------------ Test Compressed Sparse Vector " << endl;

    {
        rsc_sparse_vector_u32 csv1;
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
        csv1.push_back(20, 200);
        csv1.push_back(bm::id_max-1, 201);

        csv1.load_to(sv1);

        v = csv1.at(10);
        assert(v == 100);
        v1 = sv1.at(10);
        assert(v1 == 100);

        v = csv1.at(20);
        assert(v == 200);
        v1 = sv1.at(20);
        assert(v1 == 200);

        v = csv1.at(bm::id_max - 1);
        assert(v == 201);
        v1 = sv1.at(bm::id_max - 1);
        assert(v1 == 201);

        csv1.sync();

        v = csv1.at(10);
        assert(v == 100);
        v = csv1.at(20);
        assert(v == 200);
        v = csv1.at(bm::id_max - 1);
        assert(v == 201);

        csv1.optimize();
        v = csv1.at(10);
        assert(v == 100);
        v = csv1.at(20);
        assert(v == 200);
        v = csv1.at(bm::id_max - 1);
        assert(v == 201);

        rsc_sparse_vector_u32 csv2(csv1);
        bool same = csv2.equal(csv1);
        assert(same);

        rsc_sparse_vector_u32 csv3;
        csv3 = ::move(csv2);
        same = csv3.equal(csv1);
        assert(same);

        bm::sparse_vector_scanner<rsc_sparse_vector_u32> scanner;
        bvect::size_type pos;
        bool found = scanner.find_eq(csv1, 201, pos);
        assert(found);
        assert(pos == bm::id_max - 1);
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
        auto v0 = csv.get(11);
        assert(v0 == 12);
        csv.set(11, 12); // set it twice
        v0 = csv.get(11);
        assert(csv.get(11) == 12);
        assert(csv.get(0) == 0);
        assert(csv.get(1) == 1);


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


        // gather

        bvect::size_type idx[5] = {11, 10, 5, }, buf[5];
        unsigned v[5];
        csv.gather(&v[0], &idx[0], &buf[0], 3, bm::BM_UNKNOWN);
        assert(v[0] == 12);
        assert(v[1] == 11);
        assert(v[2] == 0);
        assert(buf[0] != bm::id_max);
        assert(buf[1] != bm::id_max);
        assert(buf[2] == bm::id_max);
    }

    cout << "inc() and merge_not_null() tests" << endl;
    {
        rsc_sparse_vector_u32::bvector_type bv { 1, 2, 10, 200, bm::id_max/2, bm::id_max-1 };
        rsc_sparse_vector_u32 csv1(bv);
        rsc_sparse_vector_u32 csv2(bv);

        csv1.sync(); csv2.sync();

        csv1.inc(1);
        csv1.inc(2, 10);

        csv2.set(200, 7);
        csv2.inc(bm::id_max/2);
        csv2.inc(bm::id_max/2, 1);
        csv2.inc(bm::id_max-1, 255);

        csv1.merge_not_null(csv2);

        assert(csv1.in_sync());

        assert(csv1.get(1) == 1);
        assert(csv1.get(2) == 10);
        assert(csv1.get(200) == 7);
        assert(csv1.get(bm::id_max/2) == 2);
        assert(csv1.get(bm::id_max-1) == 255);
    }

    {
        rsc_sparse_vector_u32::bvector_type bv;
        bv.set_range(1, 65536*2);

        for (unsigned i = 1; i < 65536*2; ++i)
        {
            rsc_sparse_vector_u32 csv1(bv);
            rsc_sparse_vector_u32 csv2(bv);
            csv1.sync();

            for (unsigned i0 = 1; i0 < i; ++i0)
            {
                csv1.set(i0, i0);
                csv1.inc(i0, i0);
            }
            for (unsigned i1 = i+1; i1 < 65536*2; ++i1)
            {
                csv2.set(i1, i1);
                csv2.inc(i1, i1);
            }
            csv1.merge_not_null(csv2);
            assert(csv1.in_sync());

            for (unsigned i0 = 1; i0 < i; ++i0)
            {
                assert(csv1.get(i0) == i0*2);
            }
            for (unsigned i1 = i+1; i1 < 65536*2; ++i1)
            {
                assert(csv1.get(i1) == i1*2);
            }
            assert(csv1.get(i) == 0);

        } // for
    }


    // set stress test
    {
        cout << "RSC set stress..." << endl;
        std::vector<std::pair<unsigned, unsigned> > vect;
        rsc_sparse_vector_u32 csv;

        const unsigned max_size = 2000000;

        cout << "Test set generation." << endl;
        for (unsigned i = 0; i < max_size; i += 2)
        {
            std::pair<unsigned, unsigned> pr(i, i + 10);
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

    cout << "random assignmnet in sync() mode...." << endl;
    {
        bvect bv { 10, 20, 100, 200, bm::id_max/2, bm::id_max-1 };

        bvect::size_type first, last, mid;
        bv.find_range(first, last);
        mid = first + ((last - first) / 4);

        rsc_sparse_vector_u32 csv1;
        rsc_sparse_vector_u32 csv2(bv);
        {
            bvect::enumerator en = bv.get_enumerator(mid);
            for (;en.valid(); ++en)
            {
                auto idx = *en;
                csv1.set(idx, 40);
            }
            en.go_to(0);
            for (;en.valid(); ++en)
            {
                auto idx = *en;
                if (idx >= mid)
                    break;
                csv1.set(idx, 40);
            }
        }
        {
            csv2.sync();
            bvect::enumerator en = bv.get_enumerator(mid);
            for (;en.valid(); ++en)
            {
                auto idx = *en;
                csv2.set(idx, 40);
            }

            en.go_to(0);
            for (;en.valid(); ++en)
            {
                auto idx = *en;
                if (idx >= mid)
                    break;
                csv2.set(idx, 40);
            }
        }
        bool eq = csv1.equal(csv2);
        if (!eq)
        {
            cerr << "Error: rsc_sparse_vector() add values check failed" << endl;
            assert(0); exit(1);
        }
    }

    {
        cout << "decode() test" << endl;

        {
            rsc_sparse_vector_u32 csv1;

            csv1.push_back(5, 1);
            csv1.push_back(6, 1);
            csv1.push_back(8, 2);
            csv1.push_back(255, 4);

            csv1.sync(true);

            for (unsigned k = 0; k < 2; ++k)
            {
                CheckCompressedDecode(csv1, 0, 1);
                CheckCompressedDecode(csv1, 0, 2);
                CheckCompressedDecode(csv1, 1, 1);

                CheckCompressedDecode(csv1, 0, 5);
                CheckCompressedDecode(csv1, 0, 6);

                CheckCompressedDecode(csv1, 256, 1);

                for (bvect::size_type i = 0; i < csv1.size(); ++i)
                {
                    CheckCompressedDecode(csv1, i, 1);
                    CheckCompressedDecode(csv1, i, csv1.size());
                }
                bvect::size_type j = csv1.size();
                for (bvect::size_type i = 0; i < csv1.size(); ++i, --j)
                {
                    bvect::size_type size = j - i;
                    if (!size)
                        break;
                    CheckCompressedDecode(csv1, i, size);
                }

                csv1.optimize();
            }
        }

        {
            rsc_sparse_vector_u32 csv1;

            csv1.push_back(bm::id_max - 1, 10);
            csv1.sync();

            for (unsigned k = 0; k < 2; ++k)
            {
                for (bvect::size_type i = bm::id_max - 20; i < csv1.size(); ++i)
                {
                    CheckCompressedDecode(csv1, i, 1);
                    CheckCompressedDecode(csv1, i, csv1.size() - i + 10);
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
            cout << "\nPASS=" << i << endl;

            sparse_vector_u32 sv(bm::use_null);
            rsc_sparse_vector_u32 csv1;

            GenerateSV(sv, i);
            auto sz = sv.size();

            csv1.load_from(sv);
            csv1.sync();

            cout << "cmp 1...";
            DetailedCompareSparseVectors(csv1, sv);
            if (sz < bm::id_max32*2)
                DetailedCheckCompressedDecode(csv1);
            cout << "ok" << endl;

            cout << "cmp 2...";
            csv1.optimize(tb);
            DetailedCompareSparseVectors(csv1, sv);
            if (sz < bm::id_max32*2)
                DetailedCheckCompressedDecode(csv1);
            cout << "ok" << endl;

            cout << "cmp 3...";
            csv1.clear_all(true);

            sv.optimize(tb);
            rsc_sparse_vector_u32 csv2;
            csv2.load_from(sv);
            DetailedCompareSparseVectors(csv2, sv);
            csv1.sync();
            if (sz < bm::id_max32*2)
                DetailedCheckCompressedDecode(csv1);

            csv2.optimize(tb);
            csv2.sync();

            DetailedCompareSparseVectors(csv2, sv);
            if (sz < bm::id_max32*2)
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

// ---------------------------------------------------------------------------

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
        csv.push_back(bm::id_max-1, 0);
        
        csv.sync();

        bvect::size_type idx;
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
            assert(bv_res.test(bm::id_max-1));

            found = scanner.find_eq(csv, 0, idx);
            assert(found);
            assert(idx == 200);
            
            csv.optimize();
        } // for
    }

    cout << " --------------- Test rsc_sparse_vector<> scan algo OK" << endl;

}


static
void TestCompressSparseSignedVector()
{
    cout << " ------------------------------ Test Compressed Sparse Vector " << endl;

    {
        rsc_sparse_vector_i32 csv1;
        assert(csv1.size() == 0);
        assert(csv1.equal(csv1));
        rsc_sparse_vector_i32 csv2;
        assert(csv1.equal(csv2));
        rsc_sparse_vector_i32 csv3(csv1);
        assert(csv3.equal(csv2));
    }


    cout << " set test " << endl;
    {
        rsc_sparse_vector_i32 csv;
        csv.set(1, -1);
        assert(csv.is_null(0));
        assert(!csv.is_null(1));
        auto v = csv.get(1);
        assert(v == -1);

        csv.push_back(10, -11);
        csv.set(11, -12);
        v = csv.get(10);
        assert(v == -11);
        v = csv.get(11);
        assert(v == -12);

        csv.set(5, 55);
        csv.set(5, -56);

        assert(csv.size() == 12);
        assert(csv.get(1) == -1);
        assert(csv.get(10) == -11);
        assert(csv.get(11) == -12);
        assert(csv.get(5) == -56);

        csv.set_null(5);
        assert(csv.is_null(5));
        assert(csv.get(1) == -1);
        assert(csv.get(10) == -11);
        assert(csv.get(11) == -12);
        assert(csv.get(5) == 0);
    }

    {
    cout << "push_back() test" << endl;
    int v, v1;

        rsc_sparse_vector_i32 csv1;
        sparse_vector_i32 sv1(bm::use_null);

        csv1.push_back(10, -100);
        assert(csv1.size() == 11);
        csv1.push_back(20, -200);
        csv1.push_back(21, 201);

        csv1.load_to(sv1);

        v = csv1.at(10);
        assert(v == -100);
        v1 = sv1.at(10);
        assert(v1 == -100);

        v = csv1.at(20);
        assert(v == -200);
        v1 = sv1.at(20);
        assert(v1 == -200);

        v = csv1.at(21);
        assert(v == 201);
        v1 = sv1.at(21);
        assert(v1 == 201);

        csv1.sync();

        DetailedCompareSparseVectors(csv1, sv1);

        v = csv1.at(10);
        assert(v == -100);
        v = csv1.at(20);
        assert(v == -200);
        v = csv1.at(21);
        assert(v == 201);

        csv1.optimize();
        v = csv1.at(10);
        assert(v == -100);
        v = csv1.at(20);
        assert(v == -200);
        v = csv1.at(21);
        assert(v == 201);

        rsc_sparse_vector_i32 csv2(csv1);
        bool same = csv2.equal(csv1);
        assert(same);

        rsc_sparse_vector_i32 csv3;
        csv3 = ::move(csv2);
        same = csv3.equal(csv1);
        assert(same);

        bm::sparse_vector_scanner<rsc_sparse_vector_i32> scanner;
        bm::id64_t pos;
        bool found = scanner.find_eq(csv1, 201, pos);
        assert(found);
        assert(pos == 21);
    }


    cout << "rsc_sparse_vector<>::const_iterator tests" << endl;
    {
        {
        rsc_sparse_vector_i32 csv1;
            {
                rsc_sparse_vector_i32::const_iterator it;
                assert(!it.valid());
            }
        csv1.push_back(0, 100);
        csv1.push_back(2, -200);

        csv1.sync();
            {
                rsc_sparse_vector_i32::const_iterator it(&csv1);
                rsc_sparse_vector_i32::const_iterator it2(it);
                assert(it2.valid());

                assert(it.valid());
                assert(*it == 100);
                bool b = it.advance();
                assert(b);
                assert(*it == 0);
                assert(it.is_null());
                ++it;
                assert(it.valid());
                assert(it.value() == -200);
            }
        }
    }


    cout << " back inserter tests" << endl;
    {
        rsc_sparse_vector_i32 csv1;
        {
        rsc_sparse_vector_i32::back_insert_iterator rs_bi = csv1.get_back_inserter();
            rs_bi.add_null();
            rs_bi.add(1);
            rs_bi.add(-2);
            rs_bi.flush();
        }
        assert(csv1.size() == 3);
        auto v = csv1.get(0);
        assert(v == 0);
        assert(csv1.is_null(0));
        v = csv1.at(1);
        assert(v == 1);
        v = csv1.at(2);
        assert(v == -2);
    }

    {
        rsc_sparse_vector_i32 csv1;
        {
        rsc_sparse_vector_i32::back_insert_iterator rs_bi = csv1.get_back_inserter();
            rs_bi.add(-1);
            rs_bi.add(2);
            rs_bi.add_null();
            rs_bi.add(-3);
            rs_bi.flush();
        }
        assert(csv1.size() == 4);
        auto v = csv1.at(0);
        assert(v == -1);
        v = csv1.at(1);
        assert(v == 2);
        v = csv1.get(2);
        assert(v == 0);
        assert(csv1.is_null(2));
        v = csv1.at(3);
        assert(v == -3);

        // test copy-range
        {
            rsc_sparse_vector_i32 csv2;
            csv2.copy_range(csv1, 4, 5);
            assert(csv2.size() == 0);

            csv2.copy_range(csv1, 0, 0);
            assert(csv2.size() == 4);
            v = csv2.at(0);
            assert(v == -1);

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
        rsc_sparse_vector_i32 csv1;
        {
        rsc_sparse_vector_i32::back_insert_iterator rs_bi = csv1.get_back_inserter();
        for (int i = 0; i < 100000; i++)
        {
            if (i&1)
            {
                rs_bi.add_null();
            }
            else
            {
                if (i & 1)
                    rs_bi.add(-i);
                else
                    rs_bi.add(i);
            }
        }
        rs_bi.flush();
        }
        csv1.optimize();

        // validation
        for (int i = 0; i < 100000; i++)
        {
            if (i&1)
            {
                assert(csv1.is_null(unsigned(i)));
            }
            else
            {
                assert(!csv1.is_null(unsigned(i)));
                auto v = csv1[unsigned(i)];
                if (i & 1)
                {
                    assert(v == -i);
                }
                else
                {
                    assert(v == i);
                }
            }
        }
    }


    {
    cout << "decode() tests" << endl;

        {
            int arr[10];
            int arr1[10];
            int arr2[10];
            rsc_sparse_vector_i32 csv1;

            csv1.push_back(5, 1);
            csv1.push_back(6, -1);
            csv1.push_back(8, 2);

            csv1.push_back(100, -4);
            csv1.sync();

            auto sz = csv1.decode(&arr[0], 100, 1);
            assert(sz==1);
            assert(arr[0] == -4);

            auto sz2 = csv1.decode_buf(&arr1[0], &arr2[0], 100, 1);
            assert(sz2==1);
            assert(arr1[0] == -4);


            csv1.set_null(100);
            csv1.sync(true);

            sz = csv1.decode(&arr[0], 100, 1);
            assert(sz == 0);

            sz2 = 2;
            sz2 = csv1.decode_buf(&arr1[0], &arr2[0], 100, 1);
            if (sz2)
            {
                cout << sz2 << endl;
            }
            assert(sz2==0);
        }

    cout << "inc() and merge_not_null() tests" << endl;
    {
        rsc_sparse_vector_i32::bvector_type bv { 1, 2, 10, 200, bm::id_max/2, bm::id_max-1 };
        rsc_sparse_vector_i32 csv1(bv);
        rsc_sparse_vector_i32 csv2(bv);

        csv1.sync(); csv2.sync();

        csv1.inc(1);
        csv1.inc(2, -10);


        csv2.set(200, -7);
        csv2.inc(200);
        csv2.inc(bm::id_max/2);
        csv2.inc(bm::id_max/2, 1);
        csv2.inc(bm::id_max-1, 255);

        csv1.merge_not_null(csv2);

        assert(csv1.in_sync());

        assert(csv1.get(1) == 1);
        assert(csv1.get(2) == -10);
        int v = csv1.get(200);
        assert(v == -6);
        assert(csv1.get(bm::id_max/2) == 2);
        assert(csv1.get(bm::id_max-1) == 255);

        csv1.inc(200, INT_MAX);
        v = csv1.get(200);
        assert(v == INT_MAX-6);
    }


    cout << "random assignment in sync() mode...." << endl;
    {
        bvect bv { 10, 20, 100, 200, bm::id_max/4 };

        bvect::size_type first, last, mid;
        bv.find_range(first, last);
        mid = first + ((last - first) / 4);

        rsc_sparse_vector_i32 csv1;
        rsc_sparse_vector_i32 csv2(bv);
        {
            bvect::enumerator en = bv.get_enumerator(mid);
            for (;en.valid(); ++en)
            {
                auto idx = *en;
                csv1.set(idx, -(int)idx);
                csv1.inc(idx);
            }
            en.go_to(0);
            for (;en.valid(); ++en)
            {
                auto idx = *en;
                if (idx >= mid)
                    break;
                csv1.set(idx, -int(idx));
                csv1.inc(idx);
            }
            assert(!csv1.in_sync());
        }
        {
            csv2.sync();
            bvect::enumerator en = bv.get_enumerator(mid);
            for (;en.valid(); ++en)
            {
                auto idx = *en;
                csv2.set(idx, -int(idx));
                csv2.inc(idx);
            }

            en.go_to(0);
            for (;en.valid(); ++en)
            {
                auto idx = *en;
                if (idx >= mid)
                    break;
                csv2.set(idx, -int(idx));
                csv2.inc(idx);
            }
            assert(csv2.in_sync());

        }
        bool eq = csv1.equal(csv2);
        if (!eq)
        {
            cerr << "Error: rsc_sparse_vector() add values check failed" << endl;
            assert(0); exit(1);
        }
    }

    cout << "random assignment in sync() mode.... [stress]" << endl;
    {
        bvect bv;
        generate_bvector(bv, 1000000, true);
        bv.optimize();

        bvect::size_type first, last, mid;
        bv.find_range(first, last);
        mid = first + ((last - first) / 4);

        rsc_sparse_vector_i32 csv1;
        rsc_sparse_vector_i32 csv2(bv);
        {
            bvect::enumerator en = bv.get_enumerator(mid);
            for (;en.valid(); ++en)
            {
                auto idx = *en;
                csv1.set(idx, -int(idx & 0xFF));
                csv1.inc(idx);

            }
            csv1.optimize();
            en.go_to(0);
            for (;en.valid(); ++en)
            {
                auto idx = *en;
                if (idx >= mid)
                    break;
                csv1.set(idx, -int(idx & 0xFF));
                csv1.inc(idx);
            }
            csv1.optimize();
        }
        // sync mode
        {
            csv2.sync();
            bvect::enumerator en = bv.get_enumerator(mid);
            for (;en.valid(); ++en)
            {
                auto idx = *en;
                csv2.set(idx, -int(idx & 0xFF));
                csv2.inc(idx);
            }
            assert(csv2.in_sync());
            csv2.optimize();

            en.go_to(0);
            for (;en.valid(); ++en)
            {
                auto idx = *en;
                if (idx >= mid)
                    break;
                csv2.set(idx, -int(idx & 0xFF));
                csv2.inc(idx);
            }
            assert(csv2.in_sync());
            csv2.optimize();

        }
        bool eq = csv1.equal(csv2);
        if (!eq)
        {
            cerr << "Error: rsc_sparse_vector() add values check failed" << endl;
            assert(0); exit(1);
        }
    }

        {
        rsc_sparse_vector_i32 csv1;

        csv1.push_back(5, 1);
        csv1.push_back(6, -1);
        csv1.push_back(8, -2);
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
            unsigned j = (unsigned)csv1.size();
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

    }

/*
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
*/
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
        cout << "\n====================================  PASS: " << i << endl;

        sparse_vector_i32 sv(bm::use_null);
        rsc_sparse_vector_i32 csv1;

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
        rsc_sparse_vector_i32 csv2;
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
        rsc_sparse_vector_i32 csv3(csv2);
        DetailedCompareSparseVectors(csv3, sv);
        }
        cout << "ok" << endl;

        cout << "cmp 5...";
        {
        rsc_sparse_vector_i32 csv4;
        csv4 = std::move(csv2);
        DetailedCompareSparseVectors(csv4, sv);

        rsc_sparse_vector_i32 csv5(std::move(csv4));
        DetailedCompareSparseVectors(csv5, sv);
        }
        cout << "ok" << endl;
    } // for
    cout << "Compressed load() stress test OK" << endl;


    }


    cout << " ------------------------------ Test Compressed Sparse Vector OK" << endl;
}





static
void TestSparseVector_Stress(unsigned count)
{

    cout << "---------------------------- Bit-plain sparse vector stress" << endl;

    cout << "Interval shift check.\n";
    // interval shift check
    for (unsigned i = 0; i < count; ++i)
    {
        unsigned fill_factor = 0;
        for (bvect::size_type min = 0; min < 10000000; min+= (unsigned)rand()%100000)
        {
            bvect::size_type max = min + (65535 * 10);
            {{
                bm::null_support null_able =
                                (min % 2 == 0) ? bm::no_null : bm::use_null;
                
                std::vector<unsigned> vect;
                bm::sparse_vector<unsigned, bvect > sv(null_able);
            
                FillSparseIntervals(vect, sv, min, max, fill_factor);
                
                bool res = CompareSparseVector(sv, vect, true);
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
            }}
            ++fill_factor;
            if (fill_factor > 2) fill_factor = 0;

            std::cout << "\r" << min << std::flush;
        } // for min
        
        cout << "." << flush;
        
    } // for i
    cout << endl;
    cout << "--------------------------- Interval shift check Ok" << endl;

    cout << "Join check" << endl;
    for (unsigned i = 0; i < 1; ++i)
    {
        unsigned fill_factor = 0;
        for (bvect::size_type min = 0; min < 10000000; min+= (unsigned)rand()%100000)
        {
            bvect::size_type max = min + (65535 * 10);
            bvect::size_type min2 = max + (unsigned)rand() % 65536;
            bvect::size_type max2 = min2 + (65535 * 10);

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
                    auto v1 = sv1[i];
                    auto v3 = sv3[i];
                    auto v4 = sv4[i];
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


template<class STR_SV>
void CheckStrSVCompare(const STR_SV& str_sv,
                        typename STR_SV::size_type limit = 0)
{
    if (!limit)
        limit = str_sv.size();

    typename STR_SV::size_type i, j;
    i = j = 0;
    auto it1 = str_sv.begin();
    auto it_end = str_sv.end();
    for (; it1 != it_end; ++it1, ++i)
    {
        const char* s1 = *it1;
        auto it2 = str_sv.begin();
        it2.go_to(i);
        for (j = i; j < limit; ++it2, ++j)
        {
            assert(it2 != it_end);
            const char* s2 = *it2;
            int r2 = ::strcmp(s1, s2);
            int r1 = str_sv.compare(i, j);
            assert (r1 == r2 || (r1 < 0 && r2 < 0) || (r1 > 0 && r2 > 0));
        } // for j
        if ((!is_silent) && ((i & 0xFF) == 0))
            cout << "\r   " << i << " / " << (limit ? limit : str_sv.size()) << flush;
    } // for i
    cout << endl << endl;
}

// -------------------------------------------------------------------------------------------

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
        const char* s_max = "XYz";
        char str[256];
        str_sparse_vector<char, bvect, 32> str_sv0;
        int cmp;

        assert(str_sv0.size() == 0);
        str_sv0.set(0, s0);
        str_sv0.get(0, str, sizeof(str));
        cmp = ::strcmp(str, s0);
        assert(cmp == 0);
        assert(str_sv0.size() == 1);

        str_sv0.set(1, s1);
        str_sv0.get(0, str, sizeof(str));
        cmp = ::strcmp(str, s0);
        assert(cmp == 0);

        str_sv0.set(bm::id_max-1, s_max);
        str_sv0.get(bm::id_max-1, str, sizeof(str));
        cmp = ::strcmp(str, s_max);
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
        assert(str_sv0.size() == bm::id_max);

        string str1;
        str_sv0.get(3, str1);
        assert(str0 == str1);
        {
            str0 = "TTF";
            str_sv0.assign(3, str0);
            str_sv0.get(3, str, sizeof(str));
            cmp = ::strcmp(str, str0.c_str());
            assert(cmp == 0);

            str0.clear();
            str_sv0.assign(3, str0);
            str_sv0.get(3, str, sizeof(str));
            cmp = ::strcmp(str, str0.c_str());
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

            bvect::size_type d = 0;
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
            assert(d == str_sv10.size() - 1);
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

            #if defined(_MSC_VER)
            ::strncpy_s(hmatr.row(0), hmatr.cols(), cs0, hmatr.cols());
            ::strncpy_s(hmatr.row(1), hmatr.cols(), cs1, hmatr.cols());
            ::strncpy_s(hmatr.row(2), hmatr.cols(), cs2, hmatr.cols());
            #else
            ::strncpy(hmatr.row(0),  cs0, hmatr.cols());
            ::strncpy(hmatr.row(1),  cs1, hmatr.cols());
            ::strncpy(hmatr.row(2),  cs2, hmatr.cols());
            #endif
            
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
            auto ref = str_sv0[3];
            const char* s = ref;
            cmp = ::strcmp(s, str0.c_str());
            assert(cmp == 0);
            str_sv0[3] = "333";
            str_sv0.get(3, str, sizeof(str));

            ref = str_sv0[3];
            s = ref;
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
        bvect::size_type str_max = str_sv0.effective_max_str();
        assert(str_max == 33);
        str_sv0[0] = "1";
        str_max = str_sv0.effective_max_str();
        assert(str_max == 33);
        str_sv0[1] = "11";
        str_max = str_sv0.effective_max_str();
        assert(str_max == 33);
        str_sv0[2] = "123";
        str_max = str_sv0.effective_max_str();
        assert(str_max == 33);

        str_sv0.clear_range(1, 234567);

        char str[256];
        str_sv0.get(1, str, sizeof(str));
        assert(str[0] == 0);
        str_sv0.get(2, str, sizeof(str));
        assert(str[0] == 0);
    }

    {
        str_sparse_vector<char, bvect, 32> str_sv0;
        str_sv0[0] = "1";
        str_sv0[1] = "11";
        str_sv0[2] = "123";

        bvect::size_type pos;
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

        str_sparse_vector<char, bvect, 32>::octet_freq_matrix_type occ_matrix;
        str_sparse_vector<char, bvect, 32>::slice_octet_matrix_type remap_matrix1;
        str_sparse_vector<char, bvect, 32>::slice_octet_matrix_type remap_matrix2;

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
        assert(cmp == 0);
        cmp = str_sv1.compare(0, "1");
        assert(cmp == 0);

        bm::heap_matrix<char, 1024, 64, bvect::allocator_type> hmatr(true);

        // test remap decoder
        {
            bvect::size_type d = 0;
            char *s;

            d = str_sv1.decode(hmatr, 0, 1);
            s = hmatr.row(0);
            cmp = ::strcmp(s, "1");
            assert(cmp == 0);
            assert(d == 1);
        }

        str_sv1.get(1, str, sizeof(str));
        cmp = ::strcmp(str, "11");
        assert(cmp == 0);
        cmp = str_sv1.compare(1, "11");
        assert(cmp == 0);

        str_sv1.get(2, str, sizeof(str));
        cmp = ::strcmp(str, "123");
        assert(cmp == 0);

        // test remap decoder
        {
            bvect::size_type d = 0;
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
        assert(cmp == 0);

        {
            string s0 = "113";
            str_sv1.assign(4, s0);
            str_sv1.get(4, s);
            assert(s == s0);
            cout << s << endl;
            cmp = str_sv1.compare(4, "113");
            assert(cmp == 0);
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

        bvect::size_type pos, pos1;
        bm::sparse_vector_scanner<bm::str_sparse_vector<char, bvect, 32> > scanner;
        bm::sparse_vector_scanner<bm::str_sparse_vector<char, bvect, 32> > scanner1;
        scanner.bind(str_sv1, true);


        bool found = scanner.find_eq_str("1", pos);
        assert(found);
        assert(pos == 0);
        found = scanner1.lower_bound_str(str_sv1, "1", pos1);
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
        assert(sz + 1 == str_sv0.size());
        assert(str_sv0.is_null(sz));

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


// -------------------------------------------------------------------------


typedef str_sparse_vector<char, bvect, 32> str_svect_type;


static
void EraseStrCollection(str_svect_type& str_sv)
{
    std::string s_next, s_curr;
    while (str_sv.size())
    {
        bvect::size_type idx = str_sv.size() / 2;
        bvect::size_type sz = str_sv.size();
        
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
            
            bvect::size_type pos;
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
                bvect::size_type idx = bvect::size_type(it - str_coll.begin());
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
            i += (unsigned)rand() % 3;
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

            bvect::size_type i = 0;
            bm::sparse_vector_scanner<str_svect_type> scanner;
            for (const string& s : str_coll)
            {
                bvect::size_type pos;
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
                    bvect::size_type pos2;
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
        bvect::size_type i = 0;
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
                bvect::size_type pos;
                bool found = scanner.lower_bound_str(str_sv_sorted, s.c_str(), pos);
                if (!found)
                {
                    cerr << s << " not found in target." << endl;
                }
                else
                {
                    cerr << s << " is at idx=" << pos << endl;
                }
                assert(0); exit(1);
            }
            str_prev = sv_str;
            ++i;
        } // for s
        EraseStrCollection(str_sv_sorted);
    }
    
    
    
   cout << "---------------------------- Bit-plain STR sparse vector SORT test OK" << endl;

}

template<typename SV>
void EraseSVCollection(SV& sv)
{
    typename SV::value_type v_next, v_curr;
    v_next = v_curr = 0;
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
            bool found = scanner.bfind(u_sv_sorted, i, pos);
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
            i += (unsigned)rand() % 3;
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
                bool found = scanner.bfind(u_sv_sorted, u, pos);

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
                    found = scanner.bfind(u_sv_sorted, u, pos2);
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
                bool found = scanner.bfind(u_sv_sorted, u, pos);
                
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

// ---------------------------------------------------------------------------

static
void TestSignedSparseSort()
{
   std::cout << "---------------------------- TestSignedSparseSort()" << endl;
   const int max_coll = 560000;
   typedef bm::sparse_vector<int, bvect > i_svect_type;

   {
       std::vector<int> i_coll;
       i_svect_type     i_sv_sorted;

        // generate sorted vector
        for (int i = 10; i < max_coll; i+=10)
        {
            i_coll.emplace_back(-i);
        }

        std::sort(i_coll.begin(), i_coll.end());

        for (const auto u : i_coll)
            i_sv_sorted.push_back(u);
        i_sv_sorted.optimize();

        // run lower bound tests
        bm::sparse_vector_scanner<i_svect_type> scanner;

        for (int i = 0; i < max_coll; ++i)
        {
            bvect::size_type pos;
            bool found = scanner.bfind(i_sv_sorted, -i, pos);
            int u1;
            if (found)
            {
                u1 = i_sv_sorted[pos];
                assert(u1 == -i);
            }

            auto it = std::lower_bound(i_coll.begin(), i_coll.end(), -i);
            if (it != i_coll.end())
            {
                auto v = *it;
                if (v == -i)
                {
                    unsigned idx = unsigned(it - i_coll.begin());
                    int u0 = i_coll[idx];

                    if (u0 == -i)
                    {
                        assert(found);
                        assert(pos == idx);
                    }
                    else
                    {
                        assert(!found);
                        u1 = i_sv_sorted[pos];
                        assert(pos == idx);
                    }
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
       std::vector<int> u_coll;
        // generate test values vector
        for (int i = 0; i < max_coll; )
        {
            if (i & 1)
                u_coll.emplace_back(-i);
            else
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
        i_svect_type      i_sv_sorted;

        cout << "\ninsertion sort..." << endl;
        {
        std::chrono::time_point<std::chrono::steady_clock> st;
        st = std::chrono::steady_clock::now();

            unsigned i = 0;
            bm::sparse_vector_scanner<i_svect_type> scanner;
            for (const int u : u_coll)
            {
                bvect::size_type pos;
                bool found = scanner.bfind(i_sv_sorted, u, pos);

                auto sz1 = i_sv_sorted.size();

                i_sv_sorted.insert(pos, u);

                auto sz2 = i_sv_sorted.size();
                assert(sz1 + 1 == sz2);

                {
                    int u_sv = i_sv_sorted.get(pos);
                    assert(u == u_sv);
                }

                if (pos)
                {
                    int u_prev;
                    u_prev = i_sv_sorted.get(pos-1);
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
                    found = scanner.bfind(i_sv_sorted, u, pos2);
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

                    i_sv_sorted.optimize();

                    st = std::chrono::steady_clock::now();
                }
                ++i;
            } // for s
        }
        cout << endl;

        cout << "sort validation.." << endl;
        std::sort(u_coll.begin(), u_coll.end());
        int i = 0;
        int u_prev = 0;
        for (int u : u_coll)
        {
            int sv_u;
            sv_u = i_sv_sorted.get(unsigned(i));
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

                bm::sparse_vector_scanner<i_svect_type> scanner;
                bvect::size_type pos;
                bool found = scanner.bfind(i_sv_sorted, u, pos);

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

        EraseSVCollection(i_sv_sorted);
    }


   cout << "---------------------------- TestSignedSparseSort() OK" << endl;

}


// -----------------------------------------------------------------------

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

   CheckStrSVCompare(str_sv, max_coll / 1000);

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

   //print_svector_stat(cout,str_sv, true);

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
                if (str1[octet_idx] != str2[octet_idx])
                    break;
                if (!str1[octet_idx] || !str2[octet_idx])
                    break;
            }
            if (octet_idx)
            {
                for (unsigned i0 = 0; i0 < octet_idx; ++i0)
                {
                    assert(str1[i0] == str2[i0]);
                }
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
   
   cout << "\n\nTest sorted search..." << endl;
   
   for (unsigned k = 0; k < 2; ++k)
   {
        bm::sparse_vector_scanner<str_svect_type> scanner;
        bm::sparse_vector_scanner<str_svect_type> scanner2;
        scanner2.bind(str_sv_remap, true); // bind sorted vector

        for (unsigned i = 0; i < unsigned(str_coll_sorted.size()); ++i)
        {
            const string& s = str_coll_sorted[i];
            bvect::size_type pos1, pos2, pos3, pos4;

            // validate the compare function
            if (i)
            {
                int res0 = str_sv_remap.compare(0, s.c_str());
                int res1 = str_sv_sorted.compare(0, s.c_str());
                assert(res0 == res1 && res1 < 0);
                res0 = str_sv_remap.compare(i-1, s.c_str());
                res1 = str_sv_sorted.compare(i-1, s.c_str());
                assert(res0 == res1 && res1 < 0);

                if ( i+1 < unsigned(str_coll_sorted.size()))
                {
                    res0 = str_sv_remap.compare(i+1, s.c_str());
                    res1 = str_sv_sorted.compare(i+1, s.c_str());
                    assert(res0 == res1 && res1 > 0);
                }
            }

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
                found2 = scanner.bfind_eq_str(str_sv_sorted, s.c_str(), pos2);
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

inline
void GeneratePipelineTestData(std::vector<string>& str_coll,
                              str_svect_type&      str_sv,
                              unsigned max_coll = 8000000,
                              unsigned repeat_factor=10)
{
    auto bi(str_sv.get_back_inserter());
    string str;
    for (unsigned i = 10; i < max_coll; i+= (rand()&0xF))
    {
        switch (i & 0xF)
        {
        case 0: str = "AB"; break;
        case 1: str = "GTx"; break;
        case 2: str = "cnv"; break;
        default: str = "AbY11"; break;
        }
        str.append(to_string(i));

        for (unsigned k = 0; k < repeat_factor; ++k)
        {
            str_coll.emplace_back(str);
            bi = str;
        }
    } // for i
    bi.flush();
}


static
void TestSparseFindEqStrPipeline()
{
   cout << "\n---------------------------- TestSparseFindEqStrPipeline()" << endl;
   const unsigned max_coll = 8000000;
   std::vector<string> str_coll;
   str_svect_type      str_sv;

   cout << "   generate test set..." << flush;


   GeneratePipelineTestData(str_coll, str_sv, max_coll, 10);

    cout << "remap..." << flush;

    str_sv.remap();
    str_sv.optimize();

   cout << "OK" << endl;

    //bm::print_svector_stat(cout,str_sv);

    unsigned test_runs = 10000;
    std::vector<string> str_test_coll;
    for (bvect::size_type i = 0; i < test_runs; ++i)
    {
        bvect::size_type idx = (unsigned) rand() % test_runs;
        if (idx >= test_runs)
            idx = test_runs/2;
        str_test_coll.push_back(str_coll[idx]);
    }
    assert(str_test_coll.size() == test_runs);


    std::vector<unique_ptr<bvect> > res_vec1;
    bm::sparse_vector_scanner<str_svect_type> scanner;


    {
    std::chrono::time_point<std::chrono::steady_clock> s;
    std::chrono::time_point<std::chrono::steady_clock> f;
    s = std::chrono::steady_clock::now();

        for (bvect::size_type i = 0; i < test_runs; ++i)
        {
            const string& str = str_test_coll[i];

            str_svect_type::bvector_type* bv_res(new bvect);
            scanner.find_eq_str(str_sv, str.c_str(), *bv_res);
            res_vec1.emplace_back(unique_ptr<bvect>(bv_res));
        } // for
    f = std::chrono::steady_clock::now();
    auto diff = f - s;
    auto d = std::chrono::duration <double, std::milli> (diff).count();

    cout << "scanner::find_eq_str()  " << d << "ms" << endl;
    }


    bm::sparse_vector_scanner<str_svect_type>::pipeline<> pipe(str_sv);
    {
    std::chrono::time_point<std::chrono::steady_clock> s;
    std::chrono::time_point<std::chrono::steady_clock> f;
    s = std::chrono::steady_clock::now();

        for (bvect::size_type i = 0; i < test_runs; ++i)
        {
            const string& str = str_test_coll[i];
            pipe.add(str.c_str());
        }
        pipe.complete(); // finish the pipeline construction with this call

        scanner.find_eq_str(pipe); // run the search pipeline

    f = std::chrono::steady_clock::now();
    auto diff = f - s;
    auto d = std::chrono::duration <double, std::milli> (diff).count();
    cout << "scanner::pipeline:  " << d << "ms" << endl;
    }


    bm::sparse_vector_scanner<str_svect_type>::pipeline<bm::agg_opt_only_counts> pipe2(str_sv);
    {
    std::chrono::time_point<std::chrono::steady_clock> s;
    std::chrono::time_point<std::chrono::steady_clock> f;
    s = std::chrono::steady_clock::now();

        for (bvect::size_type i = 0; i < test_runs; ++i)
        {
            const string& str = str_test_coll[i];
            pipe2.add(str.c_str());
        }
        pipe2.complete(); // finish the pipeline construction with this call

        scanner.find_eq_str(pipe2); // run the search pipeline

    f = std::chrono::steady_clock::now();
    auto diff = f - s;
    auto d = std::chrono::duration <double, std::milli> (diff).count();
    cout << "scanner::pipeline::count():  " << d << "ms" << endl;
    }


    cout << "  validation..." << flush;
    {
        auto& res_vect = pipe.get_bv_res_vector();
        auto& cnt_vect = pipe2.get_bv_count_vector();

        for (size_t i = 0; i < res_vect.size(); ++i)
        {
            const bvect* bv1 = res_vec1[i].get();
            const auto* bv = res_vect[i];
            assert(bv);
            bool match = bv1->equal(*bv);
            assert(match);
            auto c = cnt_vect[i];
            auto cnt = bv->count();
            assert(cnt == c);
        }
    }
    cout << "OK" << endl;


   cout << "---------------------------- TestSparseFindEqStrPipeline() OK" << endl;
}


static
void TestCompressSparseGather()
{
    cout << "\n--------------- TestCompressSparseGather()" << endl;

    // DEBUGGING code, reading from a saved memory dump
#if 0
    {
        rsc_sparse_vector_i32 csv1(bm::use_null);
        std::vector<rsc_sparse_vector_i32::size_type>  idx;

        int res = bm::load_vector(idx, "/Users/anatoliykuznetsov/dev/git/BitMagic/tests/stress/gath_idx.vect");
        assert(res == 0);
        res = bm::file_load_svector(csv1, "/Users/anatoliykuznetsov/dev/git/BitMagic/tests/stress/gath_data.csv");
        assert(res==0);

//        unsigned i=33595; unsigned get_idx = 139324;
        unsigned i=469; unsigned get_idx = 24969598;


        {
//            assert(idx[i] == get_idx);
            rsc_sparse_vector_i32::size_type rank;
            bool b = csv1.resolve(get_idx, &rank);
            bool b1 = !csv1.is_null(get_idx);
            assert(b == b1);
        }


        assert(idx[i]==get_idx);
        auto vc = csv1.get(get_idx);


        std::vector<rsc_sparse_vector_i32::size_type>  idx_buf;
        std::vector<rsc_sparse_vector_i32::value_type> vbuf;

        assert(idx.size());
        idx_buf.resize(idx.size());
        vbuf.resize(idx.size());

//        for (unsigned k = 0; k < i; ++k)
        unsigned k = 33594;
        {
            int* arr = vbuf.data() + k;
            const rsc_sparse_vector_i32::size_type* iarr = idx.data()+k;
            rsc_sparse_vector_i32::size_type* iarr_tmp = idx_buf.data() + k;
            size_t gath_sz = idx.size()-k;
/*
            auto sz = csv1.gather(arr, iarr, iarr_tmp,
                                  (rsc_sparse_vector_i32::size_type)gath_sz,
                                  bm::BM_UNSORTED);
*/
            //for (size_t gs = 2; gs <= gath_sz; ++gs)
            {
                auto sz = csv1.gather(arr, iarr, iarr_tmp,
                                      102, //gath_sz,
                                      bm::BM_UNSORTED);
                auto v0 = vbuf.at(k);
                auto get_idx0 = idx[k];
                auto v0c = csv1.at(get_idx0);
                auto t = idx_buf[k];
                if (t != bm::id_max)
                {
                    assert(v0 == v0c);
                }
                auto v = vbuf.at(i);
                cout << sz << " " << flush;
                assert(v == vc);
                //cout << k << " " << flush;
    //            if (v == vc)
    //                break;
            }
        }
        auto v = vbuf.at(i);
        cout << vc << "=" << v << endl;
        assert(v == vc);
        return;
    }
#endif

    unsigned max_passes = 1024;
    unsigned sampling_delta = 3;

    for (unsigned pass = 0; pass < max_passes; ++pass)
    {
        rsc_sparse_vector_i32 csv1(bm::use_null);
        std::vector<int> vect_c;

        rsc_sparse_vector_i32::size_type test_size = 65536 * 256 * 2;
        rsc_sparse_vector_i32::size_type test_start = pass * 65536;
        if (test_start > bm::id_max/4)
            test_start = (unsigned)((pass % 256) * (unsigned(rand()) % 65536));

        {
            auto iit = csv1.get_back_inserter();
            iit.add_null(test_start);

            for (unsigned i = 0; i < test_size; i+=3)
            {
                if (bm::id64_t(csv1.size()) + pass + 1ull > bm::id64_t(bm::id_max))
                    break;
                int v = rand() % 65539;
                if ((v > 0) && ((v % 5) == 0)) v = -v;
                iit = v;
                iit.add_null(pass);

                vect_c.push_back(v);
            } // for i
            iit.flush();
        }

        csv1.optimize();
        csv1.sync();

        if (pass < 10)
        {
            rsc_sparse_vector_i32::size_type rank = 0;
            for (rsc_sparse_vector_i32::size_type i = 0; i < csv1.size(); ++i)
            {
                int v;
                bool b = csv1.try_get(i, v);
                if (b)
                {
                    auto vc = vect_c[rank++];
                    assert(v == vc);
                }
            }
        }

        std::vector<rsc_sparse_vector_i32::size_type>  idx;
        std::vector<rsc_sparse_vector_i32::size_type>  idx_buf;
        std::vector<rsc_sparse_vector_i32::value_type> vbuf;

        rsc_sparse_vector_i32::size_type gather_size = 65536 * 3;

        for (rsc_sparse_vector_i32::size_type i = 0; i < gather_size; i += sampling_delta)
        {
            rsc_sparse_vector_i32::size_type get_idx = test_start - 1025 + i;
            if (get_idx > bm::id_max)
                get_idx = bm::id_max - 1;
            idx.push_back(get_idx);
        }

        assert(idx.size());
        idx_buf.resize(idx.size());
        vbuf.resize(idx.size());

        auto sort_opt = bm::BM_SORTED;
        for (unsigned gpass = 0; gpass < 5; ++gpass)
        {
            auto sz = csv1.gather(vbuf.data(), idx.data(), idx_buf.data(),
                                  (rsc_sparse_vector_i32::size_type)idx.size(),
                                  sort_opt);
            assert(sz == idx.size());
            for (rsc_sparse_vector_i32::size_type i = 0; i < sz; ++i)
            {
                auto get_idx = idx[i];
                auto v = vbuf[i];
                auto t = idx_buf[i];
                if (t == bm::id_max)
                {
                    if (!csv1.is_null(get_idx))
                    {
                        cout << "  i=" << i << " get_idx = " << get_idx
                             << " sz = " << sz
                             << endl;
                        rsc_sparse_vector_i32::size_type rank;
                        bool b = csv1.resolve(get_idx, &rank);
                        cout << "  is_not_NULL=" << b << " rank="
                             << rank-1 << endl;

                        int res = bm::save_vector(idx, "gath_idx.vect");
                        assert(res==0);
                        res = bm::file_save_svector(csv1, "gath_data.csv");
                        assert(res==0);
                        cout << "DBG/Dump creates: " << "gath_data.csv" << "gath_idx.vect" << endl;

                        bool bc = csv1.is_null(get_idx);
                        assert(!b);
                        assert(bc);
                        assert(0);
                    }
                }
                else
                {
                    auto vc = csv1.at(get_idx);
                    auto vc0 = csv1.get(get_idx);
                    assert(vc == vc0);

                    if (v != vc)
                    {
                        rsc_sparse_vector_i32::size_type rank;
                        bool b = csv1.resolve(get_idx, &rank);
                        assert(b);
                        auto vvc = vect_c.at(rank-1);
                        assert(vc == vvc);

                        cout << "  i=" << i << " get_idx = " << get_idx
                             << " sz = " << sz << " vc = " << vc << " v=" << v
                             << endl;

                        int res = bm::save_vector(idx, "gath_idx.vect");
                        assert(res==0);
                        res = bm::file_save_svector(csv1, "gath_data.csv");
                        assert(res==0);
                        cout << "DBG/Dump creates: " << "gath_data.csv" << "gath_idx.vect" << endl;
                    }
                    assert(v == vc);
                }
            } // for i

            // shuffle the indexes
            {
                std::random_device rd;
                std::mt19937       g(rd());
                std::shuffle(idx.begin(), idx.end(), g);
            }
            sort_opt = bm::BM_UNSORTED;
        } // for gpass


        if (!is_silent)
            cout << "\r" << pass << "/" << max_passes << flush;
    } // pass
    cout << endl;

    cout << "--------------- TestCompressSparseGather() OK\n" << endl;
}




static
void TestCompressedSparseVectorAlgo()
{
    cout << " --------------- TestCompressedSparseVectorAlgo() 64-bit" << endl;

    {
        rsc_sparse_vector_u32 csv1(bm::use_null);
        rsc_sparse_vector_u32 csv2(bm::use_null);

        bm::sparse_vector<unsigned, bvect>::size_type pos;
        bool f;
        f = bm::sparse_vector_find_first_mismatch(csv1, csv2, pos);
        assert(!f);

        csv1.push_back(bm::id_max32 + 10, 10);
        csv1.push_back(bm::id_max32 + 11, 10);
        csv1.push_back(bm::id_max32 + 200, 0);
        csv1.push_back(bm::id_max32 + 300, 0);


        f = bm::sparse_vector_find_first_mismatch(csv1, csv2, pos);
        assert(f);
        assert(pos == bm::id_max32 + 10);

        csv1.sync();
        csv2.sync();

        f = bm::sparse_vector_find_first_mismatch(csv1, csv2, pos);
        assert(f);
        assert(pos == bm::id_max32 + 10);
    }

    {
        rsc_sparse_vector_u32 csv1(bm::use_null);
        rsc_sparse_vector_u32 csv2(bm::use_null);

        csv1.push_back(bm::id_max32 + 200, 0);
        csv1.push_back(bm::id_max32 + 300, 0);

        csv2.push_back(bm::id_max32 + 200, 0);
        csv2.push_back(bm::id_max32 + 300, 0);

        bm::sparse_vector<unsigned, bvect>::size_type pos;
        bool f;
        f = bm::sparse_vector_find_first_mismatch(csv1, csv2, pos);
        assert(!f);

        csv2.push_back(bm::id_max32 + 400, 0);
        f = bm::sparse_vector_find_first_mismatch(csv1, csv2, pos);
        assert(f);
        assert(pos == bm::id_max32 + 400);
    }

    {
        rsc_sparse_vector_u32 csv1(bm::use_null);
        rsc_sparse_vector_u32 csv2(bm::use_null);

        bm::sparse_vector<unsigned, bvect>::size_type pos;
        bool f;

        csv1.push_back(bm::id_max32 + 10, 10);
        csv1.push_back(bm::id_max32 + 11, 10);
        csv1.push_back(bm::id_max32 + 200, 0);
        csv1.push_back(bm::id_max32 + 300, 0);

        csv2 = csv1;
        f = bm::sparse_vector_find_first_mismatch(csv1, csv2, pos);
        assert(!f);
        csv2.push_back(bm::id_max32 + 400, 256);

        f = bm::sparse_vector_find_first_mismatch(csv1, csv2, pos);
        assert(f);
        assert(pos == bm::id_max32 + 400);

        csv1.optimize();
        csv2.optimize();

        csv1.sync();
        csv2.sync();

        f = bm::sparse_vector_find_first_mismatch(csv1, csv2, pos);
        assert(f);
        assert(pos == bm::id_max32 + 400);

    }

    // ----------------------------------------------------------


    {
        cout << endl << "Unique mismatch check" << endl;
        sparse_vector_u32 sv1(bm::use_null), sv2(bm::use_null);
        rsc_sparse_vector_u32 csv1(bm::use_null);
        rsc_sparse_vector_u32 csv2(bm::use_null);


        sparse_vector_u32::size_type sv_size = 225600;
        {
            unsigned v = 0; unsigned cnt = 0;
            for (unsigned j = 0; j < sv_size; j+=2)
            {
                sv1[bm::id_max32 + j] = v;
                if (++cnt > 600)
                {
                    cnt = 0; v+=2;
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

        for (unsigned p = 0; p < 4; ++p)
        {
            cout << "PASS = " << p << endl;
            chrono_taker<std::ostream> ct(cout, "sparse_vector<> unique value mismatch search");

            for (sparse_vector_u32::size_type k = 0; k < sv_size; ++k)
            {
                sparse_vector_u32::size_type j = k + bm::id_max32;

                std::chrono::time_point<std::chrono::steady_clock> st;
                st = std::chrono::steady_clock::now();

                bool is_null = sv2.is_null(j);
                if (is_null)
                {
                    continue;
                }

                sparse_vector_u32::value_type v2 = sv2[j];
                v2 = ~v2;
                sv2[j] = v2;
                f = bm::sparse_vector_find_first_mismatch(sv1, sv2, pos);
                assert(f);
                assert(pos == j);

                sv2[j] = ~v2; // restore
                assert(~v2 == sv2.get(j));
                f = bm::sparse_vector_find_first_mismatch(sv1, sv2, pos);
                assert(!f);

                v2 = csv2[j];
                v2 = ~v2;
                csv2.set(j, v2);
                f = bm::sparse_vector_find_first_mismatch(csv1, csv2, pos);
                assert(f);
                assert(pos == j);

                csv2.set(j, ~v2); // restore
                auto vv2 = csv2.get(j);
                assert(~v2 == vv2);
                f = bm::sparse_vector_find_first_mismatch(csv1, csv2, pos);
                assert(!f);


                if (k % 10000 == 0)
                {
                    sv2.optimize();
                    csv2.optimize();
                    csv2.sync();

                    std::chrono::time_point<std::chrono::steady_clock> f1 = std::chrono::steady_clock::now();
                    auto diff = f1 - st;
                    //auto d = std::chrono::duration <double, std::milli> (diff).count();

                    cout << "\r" << k << "/" << sv_size << " " <<
                            " (" << diff.count() << ")" << flush;
                }
            } // for
            cout << endl;

            switch(p)
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
        } // for p

        cout << "Unique search OK" << endl;
    }



    cout << " --------------- TestCompressedSparseVectorAlgo() OK" << endl;
}

template<class SV>
void CheckGTSearch(const SV& sv, typename SV::value_type v,
                   bm::sparse_vector_scanner<SV>& scanner)
{
    bvect bv_res, bv_gt, bv_control;
    bvect bv_ge, bv_ge_control;
    bvect bv_lt, bv_lt_control;
    bvect bv_le, bv_le_control;
    bvect bv_r_0v, bv_r_0v_control;

    scanner.find_gt_horizontal(sv, v, bv_res);
    scanner.find_gt(sv, v, bv_gt);
    scanner.find_ge(sv, v, bv_ge);
    scanner.find_lt(sv, v, bv_lt);
    scanner.find_le(sv, v, bv_le);
    scanner.find_range(sv, 0, v, bv_r_0v);

    {
    bvect bv_r_vv, bv_r_vv_control;
    scanner.find_range(sv, v, v, bv_r_vv);
    scanner.find_eq(sv, v, bv_r_vv_control);
    bool eq = bv_r_vv.equal(bv_r_vv_control);
    if (!eq)
    {
        print_bv(bv_r_vv);
        print_bv(bv_r_vv_control);
        bv_r_vv ^= bv_r_vv_control;
        cout << "diff=" << endl;
        print_bv(bv_r_vv);
        assert(eq);exit(1);
    }
    }

    auto it = sv.begin();
    auto it_end = sv.end();
    for (typename SV::size_type i(0); it != it_end; ++it, ++i)
    {
        if (!it.is_null())
        {
            auto v1 = *it;
            if (v1 > v)
                bv_control.set(i);
            if (v1 >= v)
                bv_ge_control.set(i);
            if (v1 < v)
                bv_lt_control.set(i);
            if (v1 <= v)
                bv_le_control.set(i);
            if (v < 0)
            {
                if (v1 <= 0 && v1 >= v)
                    bv_r_0v_control.set(i);
            }
            else
            {
                if (v1 >= 0 && v1 <= v)
                    bv_r_0v_control.set(i);
            }
        }
    } // for
    bool eq = bv_res.equal(bv_control);
    if (!eq)
    {
        cout << "1. result for v >=" << v << " :" << endl;
        print_bv(bv_res);
        bv_res ^= bv_control;
        cout << "diff=" << endl;
        print_bv(bv_res);
        assert(eq);exit(1);
    }
    eq = bv_res.equal(bv_gt);
    if (!eq)
    {
        cout << "2. result for v >=" << v << " :" << endl;
        print_bv(bv_gt);
        bv_gt ^= bv_control;
        cout << "diff=" << endl;
        print_bv(bv_gt);
        assert(eq);exit(1);
    }
    eq = bv_ge.equal(bv_ge_control);
    if (!eq)
    {
        cout << "3. result for v >=" << v << " :" << endl;
        print_bv(bv_ge);
        bv_ge ^= bv_ge_control;
        cout << "diff=" << endl;
        print_bv(bv_ge);
        assert(eq);exit(1);
    }
    eq = bv_lt.equal(bv_lt_control);
    if (!eq)
    {
        cout << "4. result for v <" << v << " :" << endl;
        print_bv(bv_lt);
        bv_lt ^= bv_lt_control;
        cout << "diff=" << endl;
        print_bv(bv_lt);
        assert(eq);exit(1);
    }
    eq = bv_le.equal(bv_le_control);
    if (!eq)
    {
        cout << "5. result for v <=" << v << " :" << endl;
        print_bv(bv_le);
        bv_le ^= bv_le_control;
        cout << "diff=" << endl;
        print_bv(bv_le);
        assert(eq);exit(1);
    }
    eq = bv_r_0v.equal(bv_r_0v_control);
    if (!eq)
    {
        cout << "6. result for [0, v] " << v << " :" << endl;
        print_bv(bv_r_0v);
        bv_r_0v ^= bv_r_0v_control;
        cout << "diff=" << endl;
        print_bv(bv_r_0v);
        assert(eq);exit(1);
    }
}



static
void TestCompressedSparseVectorScanGT()
{
    cout << " --------------- Test rsc_sparse_vector<> TestSparseVectorScanGT()" << endl;

    bm::sparse_vector_scanner<rsc_sparse_vector_u32> scanner_csv;
    bm::sparse_vector_scanner<sparse_vector_u32> scanner_sv;
    {
    rsc_sparse_vector_u32 csv(bm::use_null);
    sparse_vector_u32 sv(bm::use_null);

        sv.push_back_null(); // 0
        sv.push_back(1);
        sv.push_back(8);     // 2
        sv.push_back(8+7);
        sv.push_back(8+1);   // 4
        sv.push_back_null(); // 5
        sv.push_back(16);
        sv.push_back(0);     // 7
        sv.push_back_null(); // 8
        sv.push_back(1023);   // 9

        csv.load_from(sv);

        for (int pass = 0; pass < 2; ++pass)
        {
            auto it = sv.begin();
            auto it_end = sv.end();
            for (; it != it_end; ++it)
            {
                if (!it.is_null())
                {
                    auto v1 = *it;
                    CheckGTSearch(sv, v1, scanner_sv);
                    CheckGTSearch(csv, v1, scanner_csv);
                }
            } // for
            sv.optimize();
        } // pass
    }

    cout << "  stress test GT on random data " << endl;

    {
        const rsc_sparse_vector_u32::size_type total_size = 55 * 1024 * 1024;
        const unsigned sample_size = 2048;


        rsc_sparse_vector_u32 csv(bm::use_null);
        {
        auto bit = csv.get_back_inserter();
        for (rsc_sparse_vector_u32::size_type i = 0; i < total_size; )
        {
            unsigned v = unsigned(rand()) & 0xFFF;
            unsigned len = unsigned(rand()) % 256;
            for (rsc_sparse_vector_u32::size_type j = 0; j < len; ++j)
                bit = v;
            i += len;
            len = len & 0xF;
            if (len)
                bit.add_null(len);
            i += len;
        } // for i
        bit.add_null(65536 * 4);
        bit = ~0u;

        bit.flush();
        csv.sync();
        cout << " Data generation Ok" << endl;
        }

        bm::random_subset<bvect> rsub;
        bvect bv_subset;
        const bvect* bv_null = csv.get_null_bvector();
        rsub.sample(bv_subset, *bv_null, sample_size);
        bv_subset.set(csv.size()-1);

        for (unsigned pass = 0; pass < 2; ++pass)
        {
            cout << "  PASS = " << pass << " sample count = " << bv_subset.count() << endl;
            auto en = bv_subset.get_enumerator(0);
            for (unsigned cnt = 0;en.valid(); ++en, ++cnt)
            {
                auto i = *en;
                assert (!csv.is_null(i));
                auto v1 = csv.get(i);
                CheckGTSearch(csv, v1, scanner_csv);
                if (cnt & 0xF)
                    cout << "\r" << cnt << "/" << sample_size << flush;
            }

            cout << endl;
            csv.optimize();
        } // for

    }

    cout << " --------------- Test rsc_sparse_vector<> TestSparseVectorScanGT() OK " << endl;
}



// -----------------------------------------------------------------------

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
        bvect::size_type idx;
        
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
        const unsigned test_count = 160;
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
void BvectorFindReverseTest()
{
    cout << "---------------------------- BvectorFindReverseTest()" << endl;

    bool b;
    bvect::size_type pos;


    cout << "Check inverted bvector..." << endl;
    {
        bvect bv;
        bv.flip();
        b = bv.find_reverse(0, pos);
        assert(b);
        assert(pos == 0);
        b = bv.find_reverse(65535, pos);
        assert(b);
        assert(pos == 65535);

        b = bv.find_reverse(bm::id_max-1, pos);
        assert(b);
        assert(pos == bm::id_max-1);
        b = bv.find_reverse(bm::id_max, pos);
        assert(b);
        assert(pos == bm::id_max-1);
    }

    cout << "Check bit bvector..." << endl;
    {
        bvect bv;
        bv[100] = true;
        b = bv.find_reverse(100, pos);
        assert(b);
        assert(pos == 100);

        b = bv.find_reverse(256, pos);
        assert(b);
        assert(pos == 100);

        bv[101] = true;
        b = bv.find_reverse(256, pos);
        assert(b);
        assert(pos == 101);

        bv[65355] = true;

        b = bv.find_reverse(256, pos);
        assert(b);
        assert(pos == 101);

        bv[100] = false;
        bv[101] = false;
        b = bv.find_reverse(256, pos);
        assert(!b);

        b = bv.find_reverse(bm::id_max/2, pos);
        assert(b);
        assert(pos == 65355);

        bv[65355*4] = true;
        b = bv.find_reverse(65355*2, pos);
        assert(b);
        assert(pos == 65355);


    }
    cout << "Check GAP bvector..." << endl;
    {
        bvect bv(bm::BM_GAP);
        bv[100] = true;
        b = bv.find_reverse(100, pos);
        assert(b);
        assert(pos == 100);

        b = bv.find_reverse(256, pos);
        assert(b);
        assert(pos == 100);

        bv[101] = true;
        b = bv.find_reverse(256, pos);
        assert(b);
        assert(pos == 101);

        bv[65355] = true;

        bv[100] = false;
        bv[101] = false;
        b = bv.find_reverse(256, pos);
        assert(!b);

        b = bv.find_reverse(bm::id_max/2, pos);
        assert(b);
        assert(pos == 65355);

        bv[65355*4] = true;
        b = bv.find_reverse(65355*2, pos);
        assert(b);
        assert(pos == 65355);
    }



    cout << "---------------------------- BvectorFindReverseTest() OK" << endl;
}

inline
void PrintStacks(unsigned max_cnt = 10)
{
(void)max_cnt;
#ifdef MEM_DEBUG
    #ifdef BM_STACK_COLL
    unsigned cnt(0);
    for (auto it = g_alloc_trace_map.begin();
        it != g_alloc_trace_map.end() && cnt < max_cnt; ++it)
    {
        cout << "\n--------------------STACK_TRACE: " << cnt++ << endl;
        cout << it->second << endl;
    }
    #endif
#else
    cout << "Stack tracing not enabled (use #define BM_STACK_COLL)" << endl;
#endif
}

static
bool CheckAllocLeaks(bool details = false, bool abort = true)
{
(void)details; (void)abort;
#ifdef MEM_DEBUG
    if (details)
    {
        cout << "[--------------  Allocation digest -------------------]" << endl;
        cout << "Number of BLOCK allocations = " <<  dbg_block_allocator::na_ << endl;
        cout << "Number of PTR allocations = " <<  dbg_ptr_allocator::na_ << endl << endl;
    }

    if (dbg_block_allocator::balance() != 0)
    {
        cout << "ERROR! Block memory leak! " << endl;
        cout << "leaked blocks: " << dbg_block_allocator::balance() << endl;
        PrintStacks();
        if (!abort)
            return true;

        assert(0);exit(1);
    }

    if (dbg_ptr_allocator::balance() != 0)
    {
        cout << "ERROR! Ptr memory leak! " << endl;
        cout << "leaked blocks: " << dbg_ptr_allocator::balance() << endl;
        PrintStacks();
        if (!abort)
            return true;
        assert(0);exit(1);
    }
    cout << "[------------  Debug Allocation balance OK ----------]" << endl;
#endif
    return false;
}



static
void show_help()
{
    std::cout
        << "BitMagic C++ stress test (64-bit vectors)." << endl
        << "-h                - help" << endl
        << "-llevel (or -ll)  - low level tests" << endl
        << "-support (or -s)  - support containers " << endl
        << "-bvbasic (or -bvb - bit-vector basic " << endl
        << "-bvser                - bit-vector serialization " << endl
        << "-bvops (-bvo, -bvl)  - bit-vector logical operations" << endl
        << "-bvshift (or -bvs)- bit-vector shifts " << endl
        << "-rankc (or -rc)   - rank-compress " << endl
        << "-agg (or -aggregator) - bm::aggregator " << endl
        << "-sv                   - test sparse vectors" << endl
        << "-csv                  - test compressed sparse vectors" << endl
        << "-strsv                - test string sparse vectors" << endl
        << "-cc                   - test compresses collections" << endl
        << endl
        << "-onlystress    - run ONLY stress tests " << endl
        << "-nostress      - do NOT run stress tests " << endl
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
bool         is_only_stress = false;
bool         is_nostress = false;

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
        if (arg == "-bvser")
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
        if (arg == "-onlystress")
        {
            is_only_stress = true;
            continue;
        }
        if (arg == "-nostress")
        {
            is_nostress = true;
            continue;
        }

    } // for i
    return 0;
}

#define BM_EXPAND(x)  x ## 1
#define EXPAND(x)     BM_EXPAND(x)


int main(int argc, char *argv[])
{
    time_t      start_time = time(0);
    time_t      finish_time;

#if !defined(BM_ASSERT) || (EXPAND(BM_ASSERT) == 1)
    cerr << "Build error: Test build with undefined BM_ASSERT" << endl;
    exit(1);
#endif

    {
    auto ret = parse_args(argc, argv);
    if (ret != 0)
        return ret;
    }
/*
    {
        typedef str_sparse_vector<char, bvect, 3> str_svect_type;
        str_svect_type sv0(bm::use_null);

        file_load_svector(sv0, "/Volumes/WD-MacOS/svstr-shuf/100M.out.shuf.vec");

        cout << sv0.size() << endl;


        //print_str_svector_stat(sv0);
        //print_svector_stat(sv0);

        size_t samples = 65536;
//        size_t samples = 10000;

        typedef bm::agg_run_options<true, true, true> scanner_custom_mask_opt;
        bm::sparse_vector_scanner<str_svect_type> scanner;

        typedef typename str_svect_type::bvector_type::allocator_type::allocator_pool_type allocator_pool_type;
        allocator_pool_type  pool; // local pool for blocks

        size_t wcnt = 0;
        auto it = sv0.begin();
        auto it_end = sv0.end();
        str_svect_type::bvector_type bv_mask;
        typename str_svect_type::bvector_type::mem_pool_guard mp_guard_bv;
        mp_guard_bv.assign_if_not_set(pool, bv_mask); // pool bv_mask for faster construction
        bv_mask.resize(sv0.size());
        bv_mask.set_range(0, sv0.size());

        for (size_t idx=0; idx < sv0.size(); ++wcnt)
        {
            cout << "\nwindow=" << wcnt << " window start= " << idx << endl;
            bm::chrono_taker tt("scanner::pipeline find_eq_str()", 0);

            bm::sparse_vector_scanner<str_svect_type>::pipeline<scanner_custom_mask_opt> pipe1_and(sv0);
            pipe1_and.options().batch_size = 5;
            pipe1_and.set_search_mask(&bv_mask); // associate search mask with the pipeline
            pipe1_and.set_search_count_limit(2); // only need pairs

            for (auto k = idx; k < (idx + samples); ++k, ++it)
            {
                assert(it != it_end);
                if (sv0.is_null(k))
                    continue;
                const char* s = it.value();
                pipe1_and.add(s); // this will search within defined mask
            } // for k

            if (pipe1_and.size())
            {
                pipe1_and.complete();
                scanner.find_eq_str(pipe1_and); // run the search

                {
                    auto& res_vect = pipe1_and.get_bv_res_vector();

                    // iterate over results, run some checks...
                    size_t res_sz = res_vect.size();
                    for (size_t i = 0; i < res_sz; ++i)
                    {
                        if (const str_svect_type::bvector_type* bv_res = res_vect[i])
                            bv_mask.bit_sub(*bv_res); // subtract all results from search mask
                    }
                }
            } // if size()

            idx += samples;

            bv_mask.keep_range(idx, sv0.size());
            auto cnt = bv_mask.count();
            cout << "bv_mask.count()=" << cnt << endl;

            // correct the NULL vector, which is against the rules since
            // we have an active iterator on the same object
            //
            str_svect_type::bvector_type* bv_null = sv0.get_null_bvect();
            bv_null->bit_and(bv_mask);
        } // for idx
*/

/*
        for (unsigned bs = 0; bs < 1000; ++bs)
        {
            cout << " batch=" << bs << endl;

            typedef bm::agg_run_options<true, false, true> scanner_custom_mask_opt;
            bm::sparse_vector_scanner<str_svect_type>::pipeline<scanner_custom_mask_opt> pipe1_and(sv0);
            pipe1_and.options().batch_size = bs;

            {
                bm::chrono_taker tt("scanner::pipeline find_eq_str()", svect.size());

                // add all the search items to the pipeline
                for (size_t i = 0; i < svect.size(); ++i)
                {
                    const string& str = svect[i];
                    pipe1_and.add(str.c_str()); // this will search within defined mask
                }
                pipe1_and.complete();
                scanner.find_eq_str(pipe1_and); // run the search
            }
        
        }
    }
*/

    // -----------------------------------------------------------------
    
    if (is_all || is_bvbasic)
    {

        SyntaxTest();
        CheckAllocLeaks(false);

        GenericBVectorTest();
        CheckAllocLeaks(false);

        SetTest();
        CheckAllocLeaks(false);

         ArenaTest();
         CheckAllocLeaks(false);

         FreezeTest();
         CheckAllocLeaks(false);

        ExportTest();
        CheckAllocLeaks(false);

        ResizeTest();
        CheckAllocLeaks(false);

        EmptyBVTest();
        CheckAllocLeaks(false);

        EnumeratorTest();
        CheckAllocLeaks(false);

        RSIndexTest();
        CheckAllocLeaks(false);

        CountRangeTest();
        CheckAllocLeaks(false);

        IsAllOneRangeTest();
        CheckAllocLeaks(false);

        KeepRangeTest();
        CheckAllocLeaks(false);

        OptimizeTest();
        CheckAllocLeaks(false);

        BvectorFindReverseTest();
        CheckAllocLeaks(false);

        RankFindTest();
        CheckAllocLeaks(false);

        BvectorBitForEachTest();
        CheckAllocLeaks(false);

        GetNextTest();
        CheckAllocLeaks(false);

        BvectorIncTest();
        CheckAllocLeaks(false);

        BvectorBulkSetTest();
        CheckAllocLeaks(false);

        GAPTestStress();
        CheckAllocLeaks(false);

        SimpleRandomFillTest();
        CheckAllocLeaks(false);

        RangeRandomFillTest();
        CheckAllocLeaks(false);

        RangeCopyTest();
        CheckAllocLeaks(false);

        BvectorFindFirstDiffTest();
        CheckAllocLeaks(false);

        ComparisonTest();
        CheckAllocLeaks(false);

        IntervalEnumeratorTest();
        CheckAllocLeaks(false);

        BVImportTest();
        CheckAllocLeaks(false);
    }
    
    if (is_all || is_bvser || is_bvbasic)
    {
        //SerializationCompressionLevelsTest();
        SerializationTest();
        CheckAllocLeaks(false);

        DesrializationTest2();
        CheckAllocLeaks(false);

        SparseSerializationTest();
        CheckAllocLeaks(false);

        RangeDeserializationTest();
        CheckAllocLeaks(false);

    }

    if (is_all || is_bvshift)
    {
         BvectorShiftTest();
         CheckAllocLeaks(false);

         BvectorInsertTest();
         CheckAllocLeaks(false);

         BvectorEraseTest();
         CheckAllocLeaks(false);
    }

    if (is_all || is_rankc)
    {
         AddressResolverTest();
         TestRankCompress();
    }
    if (is_all || is_bvops)
    {
        if (!is_only_stress)
        {
            AndOperationsTest();
            CheckAllocLeaks(false);

            AndOrOperationsTest(true); // enable detailed check
            CheckAllocLeaks(false);

            OrOperationsTest();
            CheckAllocLeaks(false);

            XorOperationsTest();
            CheckAllocLeaks(false);

            SubOperationsTest();
            CheckAllocLeaks(false);
        }

        if (!is_nostress)
        {
            const unsigned repeats = 20;

            StressTest(repeats, 0); // OR
            CheckAllocLeaks(false);

            StressTest(repeats, 3); // AND
            CheckAllocLeaks(false);

            StressTest(repeats, 1); // SUB
            CheckAllocLeaks(false);

            StressTest(repeats, 2); // XOR
            CheckAllocLeaks(false);
        }
    }


    if (is_all || is_agg)
    {
         AggregatorTest();
         CheckAllocLeaks(false);

         StressTestAggregatorOR(5);
         CheckAllocLeaks(false);

         StressTestAggregatorAND(10);
         CheckAllocLeaks(false);

         StressTestAggregatorShiftAND(5);
         CheckAllocLeaks(false);

    }

    if (is_all || is_sv)
    {
         TestSparseVector();
         CheckAllocLeaks(false);

         TestSignedSparseVector();
         CheckAllocLeaks(false);

         TestSparseVectorAlgo();
         CheckAllocLeaks(false);

         TestSparseVectorInserter();
         CheckAllocLeaks(false);

         TestSparseVectorGatherDecode();
         CheckAllocLeaks(false);

         TestSparseVectorSerial();
         CheckAllocLeaks(false);

         TestSignedSparseVectorSerial();
         CheckAllocLeaks(false);

         TestSparseVectorSerialization2();
         CheckAllocLeaks(false);

         TestSparseVectorTransform();
         CheckAllocLeaks(false);

         TestSparseVectorRange();
         CheckAllocLeaks(false);

         TestSparseVectorFilter();
         CheckAllocLeaks(false);

         TestSparseVectorScan();
         CheckAllocLeaks(false);

         TestSignedSparseVectorScan();
         CheckAllocLeaks(false);

         TestSparseSort();
         CheckAllocLeaks(false);

         TestSignedSparseSort();
         CheckAllocLeaks(false);

         TestSignedSparseSort();
         CheckAllocLeaks(false);
    }

    if (is_all || is_csv)
    {

         TestCompressSparseVector();
         CheckAllocLeaks(false);

         TestCompressSparseGather();
         CheckAllocLeaks(false);

         TestCompressSparseSignedVector();
         CheckAllocLeaks(false);

         TestCompressedSparseVectorAlgo();
         CheckAllocLeaks(false);

         TestCompressedSparseVectorScanGT();
         CheckAllocLeaks(false);

         TestCompressSparseVectorSerial();
         CheckAllocLeaks(false);

         TestCompressedSparseVectorScan();
         CheckAllocLeaks(false);

        if (!is_nostress)
        {
            TestSparseVector_Stress(2);
            CheckAllocLeaks(false);
        }

    }

    if (is_all || is_c_coll)
    {
        TestCompressedCollection();
        CheckAllocLeaks(false);
    }

    if (is_all || is_str_sv)
    {
         TestStrSparseVector();
         CheckAllocLeaks(false);

         TestSparseFindEqStrPipeline();
         CheckAllocLeaks(false);
         TestStrSparseSort();
         CheckAllocLeaks(false);

        if (!is_nostress)
        {
            StressTestStrSparseVector();
            CheckAllocLeaks(false);
        }
    }


    // -----------------------------------------------------------------

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
