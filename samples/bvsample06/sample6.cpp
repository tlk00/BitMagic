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

/** \example sample6.cpp
 This example demonstrates using of custom memory allocators.
 In this case allocator works as a memory checker, counts number of 
 allocations and deallocations to make sure that there is no memory leaks.
*/
/*! \file sample6.cpp
    \brief Example: bvector<> custom memory allocator
*/
#include <iostream>
#include <cassert>

#include "bm.h"
#include "bmundef.h" /* clear the pre-proc defines from BM */

using namespace std;


// Custom allocator keeps number of all alloc/free calls.
// It also reservs the front word of the allocated block and saves
// number of elements allocated. On deallocation it makes sure
// it deallocates the same size as allocated
//
// Please note, that this sample allocator is NOT compatible with SIMD
// optimizations, requiring special address alignment
//

class dbg_block_allocator
{
public:
static unsigned na_;
static unsigned nf_;

    static bm::word_t* allocate(size_t n, const void *)
    {
        ++na_;
        assert(n);
        bm::word_t* p =
            (bm::word_t*) ::malloc((n+1) * sizeof(bm::word_t));
        *p = (unsigned)n;
        return ++p;
    }

    static void deallocate(bm::word_t* p, size_t /* n */)
    {
        ++nf_;
        --p;
        ::free(p);
    }

    static int balance()
    {
        return int(nf_ - na_);
    }
};

unsigned dbg_block_allocator::na_ = 0;
unsigned dbg_block_allocator::nf_ = 0;

class dbg_ptr_allocator
{
public:
static unsigned na_;
static unsigned nf_;

    static void* allocate(size_t n, const void *)
    {
        ++na_;
        assert(sizeof(size_t) == sizeof(void*));
        void* p = ::malloc((n+1) * sizeof(void*));
        size_t* s = (size_t*) p;
        *s = n;
        return (void*)++s;
    }

    static void deallocate(void* p, size_t /* n */)
    {
        ++nf_;
        size_t* s = (size_t*) p;
        --s;
        ::free(s);
    }

    static int balance()
    {
        return int(nf_ - na_);
    }

};

unsigned dbg_ptr_allocator::na_ = 0;
unsigned dbg_ptr_allocator::nf_ = 0;


typedef
bm::mem_alloc<dbg_block_allocator, dbg_ptr_allocator,
            bm::alloc_pool<dbg_block_allocator, dbg_ptr_allocator>> dbg_alloc;

typedef bm::bvector<dbg_alloc> bvect;

int main(void)
{
    try
    {
        {
            bvect bv;

            bv[10] = true;
            bv[100000] = true;
            bv[10000000] = false;
        }

        cout << "Number of BLOCK allocations = " <<  dbg_block_allocator::na_ << endl;
        cout << "Number of PTR allocations = " <<  dbg_ptr_allocator::na_ << endl;
    }
    catch(std::exception& ex)
    {
        std::cerr << ex.what() << std::endl;
        return 1;
    }

    assert(dbg_block_allocator::balance() == 0);
    assert(dbg_ptr_allocator::balance() == 0);

    return 0;
}

  
