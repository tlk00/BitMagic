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

#define BM64ADDR

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

#include <bm.h>
#include <bmrandom.h>
#include <bmaggregator.h>
#include <bmvmin.h>
#include <bmdbg.h>
#include <bmalgo.h>
#include <bmsparsevec_util.h>
#include <bmtimer.h>

#include <bmsparsevec.h>
#include <bmsparsevec_algo.h>
#include <bmsparsevec_serial.h>
//#include <bmsparsevec_compr.h>
//#include <bmstrsparsevec.h>

using namespace bm;
using namespace std;


#include "gena.h"
#include "test_util.h"


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
            exit(1);
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
            assert(0);
            exit(1);
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

static
void SyntaxTest()
{
    cout << "------------------------------------ SyntaxTest()" << endl;
    ref_vect vect;
    generate_vect_simpl0(vect);
    {
        bvect64 bv0;

        load_BV_set_ref(bv0, vect);
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

        load_BV_set_ref(bv0, vect);
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
        
        load_BV_set_ref(bv0, vect);
        compare_BV_set_ref(bv0, vect);

        bvect64 bv1(bm::BM_GAP);
        load_BV_set_ref(bv1, vect);
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
        VerifyCountRange(bv1, bc_arr, bm::id_max-2000, bm::id_max-1);
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
    
    cout << "Find Rank test stress 1\n" << endl;
    
    {
        const bvect::size_type max_size = base_idx+200000;
        const bvect::size_type min_size = base_idx-20000;
        bvect bv1;
        for (bvect::size_type i = base_idx; i < max_size; i += rand()%5)
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
            throw 1;
        }
        /*
        std::vector<bvect::size_type>* vp = (std::vector<bvect::size_type>*)handle_ptr;
        vp->push_back(bit_idx);
        */
        prev = bit_idx;
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
    cout << "minvector count = " << min_count << endl;
    bvect::size_type count = bvect_full.count();
    bvect::size_type full_count = bvect_full.recalc_count();
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
            unsigned num = unsigned(::rand()) % iter;
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
            unsigned num = unsigned(::rand()) % iter;
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
        unsigned num = unsigned(::rand()) % iter;
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
        unsigned num = unsigned(rand()) % iter;
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
//        const bvect::size_type from = bvect::size_type(bm::id_max-1)-(65536*256);
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
   bvect::size_type BITVECT_SIZE = bvect::size_type(bm::id_max32) * 3;

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
    operation_deserializer<bvect>::deserialize(*bv_target_s,
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
    operation_deserializer<bvect>::deserialize(*bv_target_s,
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
    //print_stat(bvect_full1);
    
    size_t slen = bm::serialize(bvect_full1, sermem);

    cout << "Serialized len = " << slen << endl;

    bvect        bvect_full3;
    bm::deserialize(bvect_full3, sermem);
    bvect*  bv_target_s = new bvect();
    operation_deserializer<bvect>::deserialize(*bv_target_s,
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
    operation_deserializer<bvect>::deserialize(*bv_target_s,
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
        exit(1);
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

    operation_deserializer<bvect>::deserialize(*bv_target_s,
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
    operation_deserializer<bvect>::deserialize(bv_target_s,
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
   bvect::size_type BITVECT_SIZE = bvect::size_type(bm::id_max32) * 3;

   bvect  bvtotal;
   bvect::size_type size = BITVECT_SIZE - 10;
   BM_DECLARE_TEMP_BLOCK(tb)


   bvect  bv1;
   bvect  bv2;
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
    operation_deserializer<bvect>::deserialize(bv_target_s,
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
   print_stat(bv2);

   struct bvect::statistics st2;
   bv2.calc_stat(&st2);

   std::vector<unsigned char> sermemv2(st2.max_serialize_mem);

   size_t slen = bm::serialize(bv2, sermemv2.data());
   assert(slen);
   slen = 0;

   bm::deserialize(bvtotal, sermemv2.data());
   print_stat(bvtotal);
   operation_deserializer<bvect>::deserialize(bv_target_s,
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

    operation_deserializer<bvect>::deserialize(bv_target_s,
                                               sermemv2.data(),
                                               0,
                                               set_OR);
    operation_deserializer<bvect>::deserialize(bv_target_s,
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

   unsigned repetitions = 5;
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
        operation_deserializer<bvect>::deserialize(bv_target_s,
                                                   smemv.data(),
                                                   0,
                                                   set_OR);
        res = bvtotal.compare(bv_target_s);
        if (res != 0)
        {
            bvect::size_type bit_idx = bv_target_s.get_first();
            cout << bit_idx << " " << bv_target_s.get_next(bit_idx) << endl;;
            print_stat(bv_target_s);
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
            
            bvect::size_type i_pos = rand()%40000000;
            
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
                bm::chrono_taker ct("c1");
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
       assert(0);
       exit(1);
   }


   bvect::size_type count =
       operation_deserializer<bvect>::deserialize(*bv_target,
                                                  smem1,
                                                  0,
                                                  set_ASSIGN);
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
       assert(0); exit (1);
   }
   cout << "Deserialization ASSIGN into bv1 OK" << endl;

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
       operation_deserializer<bvect>::deserialize(*bv_target,
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
        operation_deserializer<bvect>::deserialize(bv_tmp2,
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
        operation_deserializer<bvect>::deserialize(bv_tmp2,
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
        load_BV_set_ref(bv0, vect);
        bvect bv1(bv0);
        
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
        load_BV_set_ref(bv0, vect0);
        compare_BV_set_ref(bv0, vect0);

        ref_vect vect1;
        generate_vect48(vect1);
        bvect bv1;
        load_BV_set_ref(bv1, vect1);
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
        print_bvector_stat(bv0);
        {
            bvect bv_i(bv0);
            bv_i.bit_and(bv1);
            compare_BV_set_ref(bv_i, vect_i);
        }
        bv1.optimize();
        print_bvector_stat(bv1);
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
        print_stat(bvect_full1);
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
        load_BV_set_ref(bv0, vect);
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
        load_BV_set_ref(bv0, vect0);
        compare_BV_set_ref(bv0, vect0);

        ref_vect vect1;
        generate_vect48(vect1);
        bvect bv1;
        load_BV_set_ref(bv1, vect1);
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
        print_bvector_stat(bv0);
        {
            bvect bv_i(bv0);
            bv_i.bit_or(bv1);
            compare_BV_set_ref(bv_i, vect_i);
        }
        bv1.optimize();
        print_bvector_stat(bv1);
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
            print_stat(bvect_full1);
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
        load_BV_set_ref(bv0, vect);
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
                clear_BV_set_ref(bv3, vect);
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
        load_BV_set_ref(bv0, vect0);
        compare_BV_set_ref(bv0, vect0);

        ref_vect vect1;
        generate_vect48(vect1);
        bvect bv1;
        load_BV_set_ref(bv1, vect1);
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
        print_bvector_stat(bv0);
        {
            bvect bv_i(bv0);
            bv_i.bit_xor(bv1);
            compare_BV_set_ref(bv_i, vect_i);
        }
        bv1.optimize();
        print_bvector_stat(bv1);
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
            print_stat(bvect_full1);
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
        load_BV_set_ref(bv0, vect);
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
                clear_BV_set_ref(bv3, vect);
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
        load_BV_set_ref(bv0, vect0);
        compare_BV_set_ref(bv0, vect0);

        ref_vect vect1;
        generate_vect48(vect1);
        bvect bv1;
        load_BV_set_ref(bv1, vect1);
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
            print_stat(bvect_full1);
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


static
void StressTest(unsigned repetitions, int set_operation = -1)
{
    const bvect::size_type BITVECT_SIZE = bm::id_max32 * 2;

   bvect::size_type RatioSum = 0;
   bvect::size_type SRatioSum = 0;
   bvect::size_type DeltaSum = 0;
   bvect::size_type SDeltaSum = 0;

   bvect::size_type clear_count = 0;

   bvect  bvtotal;
   bvtotal.set_new_blocks_strat(bm::BM_GAP);

   bm::random_subset<bvect> rsub;

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
        switch (rand() % 3)
        {
        case 1:
            start1 += size / 5;
            break;
        default:
            break;
        }

        bvect::size_type start2 = 0;
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

        delete bvect_full2;


        struct bvect::statistics st1;
        bvect_full1->calc_stat(&st1);
        bvect_full1->optimize();
        bvect_full1->optimize_gap_size();
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
        operation_deserializer<bvect>::deserialize(*bv_target_s,
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
  bvect::size_type BITVECT_SIZE = bvect::size_type(bm::id_max32) * 2;

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
        switch (rand() % 3)
        {
        case 1:
            start1 += size / 5;
            break;
        default:
            break;
        }

        bvect::size_type start2 = 0;
        switch (rand() % 3)
        {
        case 1:
            start2 += size / 5;
            break;
        default:
            break;
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
            load_BV_set_ref(bv8, vect0);
            bv8.optimize();
        }
        {
            ref_vect vect0;
            generate_vect48(vect0);
            load_BV_set_ref(bv9, vect0);
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
  bvect::size_type BITVECT_SIZE = bvect::size_type(bm::id_max32) * 2;

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
        switch (rand() % 3)
        {
        case 1:
            start1 += size / 5;
            break;
        default:
            break;
        }

        bvect::size_type start2 = 0;
        switch (rand() % 3)
        {
        case 1:
            start2 += size / 5;
            break;
        default:
            break;
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
  bvect::size_type BITVECT_SIZE = bvect::size_type(bm::id_max32) * 2;

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
            agg.combine_shift_right_and(bv_target1);
            auto cmp = bv_target1.compare(bv_target0);
            if (cmp != 0)
            {
                cerr << "Error: Mismatch! " << "STEP=" << i << endl;
                //DetailedCheckVectors(bv_target0, bv_target1);
                assert(0); exit(1);
            }
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
        print_svector_stat(sv);
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
        
        sv3.clear();
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
        print_svector_stat(sv);
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
    sv.extract_plains(&v1p[0], 16, 0);
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
    sv.extract_plains(&v2p[0], 10, 32);
        
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



// ---------------------------------------------------------------------

typedef bm::sparse_vector<unsigned, bvect > sparse_vector_u32;
typedef bm::sparse_vector<unsigned long long, bvect > sparse_vector_u64;
typedef bm::rsc_sparse_vector<unsigned, sparse_vector_u32> rsc_sparse_vector_u32;

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


static
void show_help()
{
    std::cerr
        << "BitMagic C++ stress test." << endl
        << "-h                - help" << endl
        << "-llevel (or -ll)  - low level tests" << endl
        << "-support (or -s)  - support containers " << endl
        << "-bvbasic (or -bvb - bit-vector basic " << endl
        << "-bvops (-bvo, -bvl)  - bit-vector logical operations" << endl
        << "-bvshift (or -bvs)- bit-vector shifts " << endl
        << "-rankc (or -rc)   - rank-compress " << endl
        << "-agg (or -aggregator) - bm::aggregator " << endl
        << "-sv                   - test sparse vectors" << endl
      ;
}

bool         is_all = true;
bool         is_low_level = false;
bool         is_support = false;
bool         is_bvbasic = false;
bool         is_bvops = false;
bool         is_bvshift = false;
bool         is_rankc = false;
bool         is_agg = false;
bool         is_sv = false;

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

    // -----------------------------------------------------------------
    
    if (is_all || is_bvbasic)
    {
        SyntaxTest();
        GenericBVectorTest();
        SetTest();
        ExportTest();

        ResizeTest();
        EmptyBVTest();
        EnumeratorTest();
        RSIndexTest();

        CountRangeTest();

        OptimizeTest();

        RankFindTest();

        BvectorBitForEachTest();

        GetNextTest();

        BvectorIncTest();

        BvectorBulkSetTest();

        GAPTestStress();
     
        SimpleRandomFillTest();
     
        RangeRandomFillTest();

        RangeCopyTest();

        ComparisonTest();

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
        AndOperationsTest();
        OrOperationsTest();
        XorOperationsTest();
        SubOperationsTest();

        const unsigned repeats = 20;
        StressTest(repeats, 0); // OR
        StressTest(repeats, 3); // AND
        StressTest(repeats, 1); // SUB
        StressTest(repeats, 2); // XOR

    }
    
    if (is_all || is_agg)
    {
         AggregatorTest();

         StressTestAggregatorOR(5);
         StressTestAggregatorAND(10);
         StressTestAggregatorShiftAND(5);

    }

    if (is_all || is_sv)
    {

         TestSparseVector();

         TestSparseVectorInserter();

         TestSparseVectorGatherDecode();
/*
         TestSparseVectorTransform();

         TestSparseVectorRange();

         TestSparseVectorFilter();

         TestSparseVectorScan();

         TestCompressSparseVector();

         TestCompressedSparseVectorScan();

         TestSparseVector_Stress(2);
     
         TestCompressedCollection();

         TestStrSparseVector();

         TestStrSparseSort();

         StressTestStrSparseVector();
*/
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
