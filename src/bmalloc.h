#ifndef BMALLOC__H__INCLUDED__
#define BMALLOC__H__INCLUDED__
/*
Copyright(c) 2002-2009 Anatoliy Kuznetsov(anatoliy_kuznetsov at yahoo.com)

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

For more information please visit:  http://bmagic.sourceforge.net

*/

#include <stdlib.h>
#include <new>

namespace bm
{


/*! @defgroup alloc Memory Allocation
 *  Memory allocation related units
 *  @ingroup bmagic
 *
 *  @{
 */

/*! 
  @brief Default malloc based bitblock allocator class.

  Functions allocate and deallocate conform to STL allocator specs.
  @ingroup alloc
*/
class block_allocator
{
public:
    /**
    The member function allocates storage for an array of n bm::word_t 
    elements, by calling malloc. 
    @return pointer to the allocated memory. 
    */
    static bm::word_t* allocate(size_t n, const void *)
    {
        bm::word_t* ptr;
#if defined(BMSSE2OPT) || defined(BMSSE42OPT)
# ifdef _MSC_VER
        ptr = (bm::word_t*) ::_aligned_malloc(n * sizeof(bm::word_t), 16);
#else
        ptr = (bm::word_t*) ::_mm_malloc(n * sizeof(bm::word_t), 16);
# endif

#else  
        ptr = (bm::word_t*) ::malloc(n * sizeof(bm::word_t));
#endif
        if (!ptr)
        {
            throw std::bad_alloc();
        }
        return ptr;
    }

    /**
    The member function frees storage for an array of n bm::word_t 
    elements, by calling free. 
    */
    static void deallocate(bm::word_t* p, size_t)
    {
#if defined(BMSSE2OPT) || defined(BMSSE42OPT)
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


// -------------------------------------------------------------------------

/*! @brief Default malloc based bitblock allocator class.

  Functions allocate and deallocate conform to STL allocator specs.
*/
class ptr_allocator
{
public:
    /**
    The member function allocates storage for an array of n void* 
    elements, by calling malloc. 
    @return pointer to the allocated memory. 
    */
    static void* allocate(size_t n, const void *)
    {
        void* ptr = ::malloc(n * sizeof(void*));
        if (!ptr)
        {
            throw std::bad_alloc();
        }
        return ptr;
    }

    /**
    The member function frees storage for an array of n bm::word_t 
    elements, by calling free. 
    */
    static void deallocate(void* p, size_t)
    {
        ::free(p);
    }
};

// -------------------------------------------------------------------------

/*! @brief BM style allocator adapter. 

  Template takes two parameters BA - allocator object for bit blocks and
  PA - allocator object for pointer blocks.
*/
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
    
    /*! @brief Returns copy of the block allocator object
    */
    block_allocator_type get_block_allocator() const 
    { 
        return BA(block_alloc_); 
    }

    /*! @brief Returns copy of the ptr allocator object
    */
    ptr_allocator_type get_ptr_allocator() const 
    { 
       return PA(block_alloc_); 
    }


    /*! @brief Allocates and returns bit block.
        @param alloc_factor 
            indicated how many blocks we want to allocate in chunk
            total allocation is going to be bm::set_block_size * alloc_factor
            Default allocation factor is 1
    */
    bm::word_t* alloc_bit_block(unsigned alloc_factor = 1)
    {
        return block_alloc_.allocate(bm::set_block_size * alloc_factor, 0);
    }

    /*! @brief Frees bit block allocated by alloc_bit_block.
    */
    void free_bit_block(bm::word_t* block, unsigned alloc_factor = 1)
    {
        if (IS_VALID_ADDR(block)) 
            block_alloc_.deallocate(block, bm::set_block_size * alloc_factor);
    }

    /*! @brief Allocates GAP block using bit block allocator (BA).

        GAP blocks in BM library belong to levels. Each level has a 
        correspondent length described in bm::gap_len_table<>.
        
        @param level GAP block level.
    */
    bm::gap_word_t* alloc_gap_block(unsigned level, 
                                    const gap_word_t* glevel_len)
    {
        BM_ASSERT(level < bm::gap_levels);
        unsigned len = 
            (unsigned)(glevel_len[level] / (sizeof(bm::word_t) / sizeof(gap_word_t)));

        return (bm::gap_word_t*)block_alloc_.allocate(len, 0);
    }

    /*! @brief Frees GAP block using bot block allocator (BA)
    */
    void free_gap_block(bm::gap_word_t*   block,
                        const gap_word_t* glevel_len)
    {
        BM_ASSERT(IS_VALID_ADDR((bm::word_t*)block));
         
        unsigned len = gap_capacity(block, glevel_len);
        len /= (unsigned)(sizeof(bm::word_t) / sizeof(bm::gap_word_t));
        block_alloc_.deallocate((bm::word_t*)block, len);        
    }

    /*! @brief Allocates block of pointers.
    */
    void* alloc_ptr(unsigned size = bm::set_array_size)
    {
        return ptr_alloc_.allocate(size, 0);
    }

    /*! @brief Frees block of pointers.
    */
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

/** @} */


} // namespace bm


#endif
