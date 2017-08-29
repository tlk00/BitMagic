/*
Copyright(c) 2002-2005 Anatoliy Kuznetsov(anatoliy_kuznetsov at yahoo.com)

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
*/

/** \example sample6.cpp
 This example demonstrates using of custom memory allocators.
 In this case allocator works as a memory checker, counts number of 
 allocations and deallocations to make sure that there is no memory leaks. 

For more information please visit:  http://bmagic.sourceforge.net
*/

#include <iostream>
#include <cassert>
#include "bm.h"

using namespace std;


// Custom allocator keeps number of all alloc/free calls.
// It also reservs the front word of the allocated block and saves
// number of elements allocated. On deallocation it makes sure
// it deallocates the same size as allocated

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
        *p = n;
        return ++p;
    }

    static void deallocate(bm::word_t* p, size_t n)
    {
        ++nf_;
        --p;
        assert(*p == n);
        ::free(p);
    }

    static int balance()
    {
        return nf_ - na_;
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

    static void deallocate(void* p, size_t n)
    {
        ++nf_;
        size_t* s = (size_t*) p;
        --s;
        assert(*s == n);
        ::free(s);
    }

    static int balance()
    {
        return nf_ - na_;
    }

};

unsigned dbg_ptr_allocator::na_ = 0;
unsigned dbg_ptr_allocator::nf_ = 0;


typedef bm::mem_alloc<dbg_block_allocator, dbg_ptr_allocator> dbg_alloc;

typedef bm::bvector<dbg_alloc> bvect;



int main(void)
{
    {
        bvect bv;

        bv[10] = true;
        bv[100000] = true;
        bv[10000000] = false;
    }

    cout << "Number of BLOCK allocations = " <<  dbg_block_allocator::na_ << endl;
    cout << "Number of PTR allocations = " <<  dbg_ptr_allocator::na_ << endl;

    assert(dbg_block_allocator::balance() == 0);
    assert(dbg_ptr_allocator::balance() == 0);

    return 0;
}

  
