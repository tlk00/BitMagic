#ifndef BMCALLOC__H__INCLUDED__
#define BMCALLOC__H__INCLUDED__
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

#include <stdlib.h>
#include "try_throw_catch.h"

#include "bmfunc.h"

namespace libbm
{

#if defined(BMSSE2OPT) || defined(BMSSE42OPT)
#define BM_ALLOC_ALIGN 16
#endif
#if defined(BMAVX2OPT)
#define BM_ALLOC_ALIGN 32
#endif



class block_allocator
{
public:
    static bm::word_t* allocate(size_t n, const void *)
    {
        bm::word_t* ptr;
#if defined(BM_ALLOC_ALIGN)
    #ifdef _MSC_VER
        ptr = (bm::word_t*) ::_aligned_malloc(n * sizeof(bm::word_t), BM_ALLOC_ALIGN);
    #else
        ptr = (bm::word_t*) ::_mm_malloc(n * sizeof(bm::word_t), BM_ALLOC_ALIGN);
    #endif
#else
            ptr = (bm::word_t*) ::malloc(n * sizeof(bm::word_t));
#endif
        if (!ptr)
        {
            BM_THROW( BM_ERR_BADALLOC );
        }
        return ptr;
    }

    static void deallocate(bm::word_t* p, size_t)
    {
#ifdef BM_ALLOC_ALIGN
    # ifdef _MSC_VER
            ::_aligned_free(p);
    #else
            ::_mm_free(p);
    # endif
#else
        ::free(p);
#endif
    }

};

class ptr_allocator
{
public:
    static void* allocate(size_t n, const void *)
    {
        void* ptr = ::malloc(n * sizeof(void*));
        if (!ptr)
        {
            BM_THROW( BM_ERR_BADALLOC );
        }
        return ptr;
    }

    static void deallocate(void* p, size_t)
    {
        ::free(p);
    }
};

template<class BA, class PA> class mem_alloc
{
public:
    typedef BA  block_allocator_type;
    typedef PA  ptr_allocator_type;

public:
    mem_alloc(const BA& block_alloc = BA(), const PA& ptr_alloc = PA())
    : block_alloc_(block_alloc),
      ptr_alloc_(ptr_alloc)
    {}
    
    block_allocator_type get_block_allocator() const 
    { 
        return BA(block_alloc_); 
    }

    ptr_allocator_type get_ptr_allocator() const 
    { 
       return PA(block_alloc_); 
    }

    bm::word_t* alloc_bit_block(unsigned alloc_factor = 1)
    {
        return block_alloc_.allocate(bm::set_block_size * alloc_factor, 0);
    }

    void free_bit_block(bm::word_t* block, unsigned alloc_factor = 1)
    {
		block_alloc_.deallocate(block, bm::set_block_size * alloc_factor);
    }

    bm::gap_word_t* alloc_gap_block(unsigned level, 
                                    const bm::gap_word_t* glevel_len)
    {
        unsigned len = 
            (unsigned)(glevel_len[level] / (sizeof(bm::word_t) / sizeof(bm::gap_word_t)));

        return (bm::gap_word_t*)block_alloc_.allocate(len, 0);
    }

    void free_gap_block(bm::gap_word_t*   block,
                        const bm::gap_word_t* glevel_len)
    {
        unsigned len = bm::gap_capacity(block, glevel_len);
        len /= (unsigned)(sizeof(bm::word_t) / sizeof(bm::gap_word_t));
        block_alloc_.deallocate((bm::word_t*)block, len);        
    }

    void* alloc_ptr(unsigned size = bm::set_array_size)
    {
        return ptr_alloc_.allocate(size, 0);
    }

    void free_ptr(void* p, unsigned size = bm::set_array_size)
    {
        if (p)
            ptr_alloc_.deallocate(p, size);
    }
private:
    BA            block_alloc_;
    PA            ptr_alloc_;
};

typedef mem_alloc<block_allocator, ptr_allocator> standard_allocator;


} // namespace libbm


#endif
