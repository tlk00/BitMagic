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

For more information please visit:  http://bmagic.sourceforge.net

*/

#include <stdio.h>
#include <stdlib.h>
#include "jni.h"
#include "try_throw_catch.h"

#include "bmfunc.h"

#define LOG2(x) (((x) & 8) ? 3 : (((x) & 4) ? 2 : (((x) & 2) ? 1 : 0)))
#define SIZE_OF_JAVA_LONG 8
#define DEADBEEF 0xdeadbeefdeadbeef

struct Jalloc {
//  jlongArray ja;
  jbyteArray ja;
  jobject ref;
  jsize sz;
};

static JavaVM *cached_jvm = 0;

JNIEXPORT jint JNICALL
JNI_OnLoad(JavaVM *jvm, void *reserved)
{
  cached_jvm = 0; // disable java allocations
  return JNI_VERSION_1_8;
}

static JNIEnv* JNU_GetEnv() {
  JNIEnv *env;
  jint rc = cached_jvm->GetEnv((void **)&env, JNI_VERSION_1_8);
  if (rc == JNI_EDETACHED)
    BM_THROW(BM_ERR_DETACHED);
  if (rc == JNI_EVERSION)
    BM_THROW(BM_ERR_JVM_NOT_SUPPORTED);

  return env;
}

namespace libbm
{

//static void * array_alloc(int shift, size_t n) {
//  //printf("Shift: %d\n", shift);
//  if (shift < 0)
//    BM_THROW(BM_ERR_BADALLOC);
//
//  JNIEnv *env = JNU_GetEnv();
//  jsize sz = (n >> shift) + 1 + (sizeof(Jalloc) >> 3) + 1 + 1; // The last one for overrun check
//                                                               //printf("Size to allocate: %d\n", sz);
//  jlongArray ja = env->NewLongArray(sz);
//  if (ja == NULL) {
//    printf("Could not allocate %d units.", n);
//    BM_THROW(BM_ERR_BADALLOC);
//  }
//  
//  jlong *jlongbuf = env->GetLongArrayElements(ja, 0);
//  if (jlongbuf == NULL)
//    BM_THROW(BM_ERR_BADALLOC);
//  
//  // Put the marker
//  jlongbuf[sz - 1] = 1;
//
//  Jalloc *pJalloc = reinterpret_cast<Jalloc*>(jlongbuf);
//  pJalloc->ref = env->NewGlobalRef(ja);
//  if (pJalloc->ref == NULL)
//    BM_THROW(BM_ERR_BADALLOC);
//  pJalloc->ja = ja;
//  pJalloc->sz = sz;
//
//  //printf("Requested %d, allocated %p of size %d, with marker 0x%llx, ref %p, ja %p\n", n, jlongbuf, sz, jlongbuf[sz - 1], pJalloc->ref, ja);
//  return static_cast<void*>(jlongbuf);
//}

//static void jni_free(void *p) {
//  if (p != 0) {
//    jlong *buffer = reinterpret_cast<jlong*>(static_cast<char*>(p) - sizeof(Jalloc));
//    //printf("Freeing %p\n", buffer);
//    Jalloc *pJalloc = reinterpret_cast<Jalloc*>(buffer);
//    if (pJalloc->ref) {
//      JNIEnv *env = JNU_GetEnv();
//      // Check marker
//      //     jsize sz = env->GetArrayLength(pJalloc->ja);
//      jsize sz = pJalloc->sz;
//      switch (buffer[sz - 1]) {
//      case 1L:
//        --buffer[sz - 1];
//        env->DeleteGlobalRef(pJalloc->ref);
//        env->ReleaseLongArrayElements(pJalloc->ja, buffer, 0);
//        printf("Released at %p\n", buffer);
//        break;
//      case 0:
//        printf("Double deletion at %p\n", buffer);
//        break;
//      default:
//        printf("Memory corruption 0x%llx at %p\n", buffer[sz - 1], buffer);
//      }
//    }
//  }
//}
  
static void * array_alloc(int shift, size_t n) {
  //printf("Shift: %d\n", shift);
  if (shift < 0)
    BM_THROW(BM_ERR_BADALLOC);

  JNIEnv *env = JNU_GetEnv();
  jsize sz = n + sizeof(Jalloc) + 1; // The last one for overrun check
  printf("Sizes to allocate: n: %d, Jalloc: %d, total: %d\n", n, sizeof(Jalloc), sz);
  jbyteArray ja = env->NewByteArray(sz);
  if (ja == NULL) {
    printf("Could not allocate %d units.", n);
    BM_THROW(BM_ERR_BADALLOC);
  }

  jbyte *jbuf = env->GetByteArrayElements(ja, 0);
  if (jbuf == NULL)
    BM_THROW(BM_ERR_BADALLOC);

  // Put the marker
  jbuf[sz - 1] = 1;

  Jalloc *pJalloc = reinterpret_cast<Jalloc*>(jbuf);
  pJalloc->ref = env->NewGlobalRef(ja);
  if (pJalloc->ref == NULL)
    BM_THROW(BM_ERR_BADALLOC);
  pJalloc->ja = ja;
  pJalloc->sz = sz;

  //printf("Requested %d, allocated %p of size %d, with marker 0x%llx, ref %p, ja %p\n", n, jlongbuf, sz, jlongbuf[sz - 1], pJalloc->ref, ja);
  return static_cast<void*>(jbuf);
}

static void jni_free(void *p) {
  if (p != 0) {
    jbyte *buffer = reinterpret_cast<jbyte*>(static_cast<char*>(p) - sizeof(Jalloc));
    //printf("Freeing %p\n", buffer);
    Jalloc *pJalloc = reinterpret_cast<Jalloc*>(buffer);
    if (pJalloc->ref) {
      JNIEnv *env = JNU_GetEnv();
      // Check marker
      //     jsize sz = env->GetArrayLength(pJalloc->ja);
      jsize sz = pJalloc->sz;
      switch (buffer[sz - 1]) {
      case 1:
        --buffer[sz - 1];
        env->DeleteGlobalRef(pJalloc->ref);
        env->ReleaseByteArrayElements(pJalloc->ja, buffer, 0);
        printf("Released size %d at %p\n", sz, buffer);
        break;
      case 0:
        printf("Double deletion, size %d at %p\n", sz, buffer);
        break;
      default:
        printf("Memory corruption, size %d,  0x%x at %p\n", sz, buffer[sz - 1], buffer);
      }
    }
  }
}
  
//template<typename T>
//T* jni_alloc(size_t n) {
////  printf("Allocating %d units in %s\n", n, __FUNCSIG__);
//  //printf("Size of T: %d\n", LOG2(sizeof(T)));
//  // Making array of longs to accommodate max value of unsigned int. 
//  // Java long size is 8 bytes
//  int shift = LOG2(SIZE_OF_JAVA_LONG) - LOG2(sizeof(T));
//  void* jbuffer = array_alloc(shift, n);
//  void *mem = static_cast<char*>(jbuffer) + sizeof(Jalloc);
//  fprintf(stderr, "Allocated %d at %p\n", n, mem);
//  return reinterpret_cast<T*>(mem);
//}
//
//template <> 
//void* jni_alloc(size_t n) {
//  //printf("Allocating %d units in %s\n", n, __FUNCSIG__);
//  int shift = LOG2(SIZE_OF_JAVA_LONG) - LOG2(sizeof(void*));
//  void* jbuffer = array_alloc(shift, n);
//  void *mem = static_cast<char*>(jbuffer) + sizeof(Jalloc);
//  fprintf(stderr, "Allocated %d at %p\n", n, mem);
//  return mem;
//}
//



class block_allocator
{
public:
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
      //if (cached_jvm != 0) {
      //  ptr = jni_alloc<bm::word_t>(n);
      //}
      //else
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
#if defined(BMSSE2OPT) || defined(BMSSE42OPT)
# ifdef _MSC_VER
        ::_aligned_free(p);
#else
        ::_mm_free(p);
# endif

#else  
      //if (cached_jvm != 0) {
      //  jni_free(p);
      //}
      //else
        ::free(p);
#endif
    }

};

class ptr_allocator
{
public:
    static void* allocate(size_t n, const void *)
    {
      void* ptr = 0;
      //if (cached_jvm != 0) {
      //  ptr = jni_alloc<void>(n);
      //}
      //else 
        ptr = ::malloc(n * sizeof(void*));

      if (!ptr)
      {
          BM_THROW( BM_ERR_BADALLOC );
      }
      return ptr;
    }

    static void deallocate(void* p, size_t)
    {
      //if (cached_jvm != 0)
      //  jni_free(p);
      //else
        ::free(p);
    }
};

/*
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
*/

/*
template<class BA, class PA> class alloc_pool;
typedef bm::alloc_pool<block_allocator, ptr_allocator> standard_alloc_pool;

template<class BA = block_allocator, class PA = ptr_allocator, class APool = standard_alloc_pool > class mem_alloc;
*/

class pointer_pool_array
{
public:
    enum params
    {
        n_pool_max_size = BM_DEFAULT_POOL_SIZE
    };

    pointer_pool_array() : size_(0)
    {
        allocate_pool(n_pool_max_size);
    }

    pointer_pool_array(const pointer_pool_array&) = delete;
    pointer_pool_array& operator=(const pointer_pool_array&) = delete;

    ~pointer_pool_array()
    {
        free_pool();
    }

    /// Push pointer to the pool (if it is not full)
    ///
    /// @return 0 if pointer is not accepted (pool is full)
    unsigned push(void* ptr)
    {
        if (size_ == n_pool_max_size - 1)
            return 0;
        pool_ptr_[size_++] = ptr;
        return size_;
    }

    /// Get a pointer if there are any vacant
    ///
    void* pop()
    {
        if (size_ == 0)
            return 0;
        return pool_ptr_[--size_];
    }
private:
    void allocate_pool(size_t pool_size)
    {
        pool_ptr_ = (void**)::malloc(sizeof(void*) * pool_size);
        if (!pool_ptr_)
        {
            BM_THROW(BM_ERR_BADALLOC);
        }
    }

    void free_pool()
    {
        ::free(pool_ptr_);
    }
private:
    void**     pool_ptr_;  ///< array of pointers in the pool
    unsigned  size_;                  ///< current size
};

template<class BA, class PA>
class alloc_pool
{
public:
    typedef BA  block_allocator_type;
    typedef PA  ptr_allocator_type;

public:

    alloc_pool() {}
    ~alloc_pool()
    {
        free_pools();
    }

    bm::word_t* alloc_bit_block()
    {
        bm::word_t* ptr = (bm::word_t*)block_pool_.pop();
        if (ptr == 0)
            ptr = block_alloc_.allocate(bm::set_block_size, 0);
        return ptr;
    }
    
    void free_bit_block(bm::word_t* block)
    {
        if (!block_pool_.push(block))
        {
            block_alloc_.deallocate(block, bm::set_block_size);
        }
    }

    void free_pools()
    {
        bm::word_t* block;
        do
        {
            block = (bm::word_t*)block_pool_.pop();
            if (block)
                block_alloc_.deallocate(block, bm::set_block_size);
        } while (block);
    }

protected:
    pointer_pool_array  block_pool_;
    BA                  block_alloc_;
};


template<class BA, class PA, class APool=alloc_pool<BA, PA> >
class mem_alloc
{
public:

    typedef BA      block_allocator_type;
    typedef PA      ptr_allocator_type;
    typedef APool   allocator_pool_type;

public:

    mem_alloc(const BA& block_alloc = BA(), const PA& ptr_alloc = PA())
    : block_alloc_(block_alloc),
      ptr_alloc_(ptr_alloc),
      alloc_pool_p_(0)
    {}

    mem_alloc(const mem_alloc& ma)
        : block_alloc_(ma.block_alloc_),
          ptr_alloc_(ma.ptr_alloc_),
          alloc_pool_p_(0) // do not inherit pool (has to be explicitly defined)
    {}

    mem_alloc& operator=(const mem_alloc& ma)
    {
        block_alloc_ = ma.block_alloc_;
        ptr_alloc_ = ma.ptr_alloc_;
        // alloc_pool_p_ - do not inherit pool (has to be explicitly defined)
        return *this;
    }
    
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

    /*! @brief set pointer to external pool */
    void set_pool(allocator_pool_type* pool)
    {
        alloc_pool_p_ = pool;
    }

    /*! @brief get pointer to allocation pool (if set) */
    allocator_pool_type* get_pool()
    {
        return alloc_pool_p_;
    }

    /*! @brief Allocates and returns bit block.
        @param alloc_factor
            indicated how many blocks we want to allocate in chunk
            total allocation is going to be bm::set_block_size * alloc_factor
            Default allocation factor is 1
    */
    bm::word_t* alloc_bit_block(unsigned alloc_factor = 1)
    {
        if (alloc_pool_p_ && alloc_factor == 1)
            return alloc_pool_p_->alloc_bit_block();
        return block_alloc_.allocate(bm::set_block_size * alloc_factor, 0);
    }

    /*! @brief Frees bit block allocated by alloc_bit_block.
    */
    void free_bit_block(bm::word_t* block, unsigned alloc_factor = 1)
    {
        if (alloc_pool_p_ && alloc_factor == 1)
            alloc_pool_p_->free_bit_block(block);
        else
            block_alloc_.deallocate(block, bm::set_block_size * alloc_factor);
    }

    /*! @brief Allocates GAP block using bit block allocator (BA).

        GAP blocks in BM library belong to levels. Each level has a
        correspondent length described in bm::gap_len_table<>
     
        @param level GAP block level.
        @param glevel_len table of level lengths
    */
    bm::gap_word_t* alloc_gap_block(unsigned level,
                                    const bm::gap_word_t* glevel_len)
    {
        unsigned len =
            (unsigned)(glevel_len[level] / (sizeof(bm::word_t) / sizeof(bm::gap_word_t)));

        return (bm::gap_word_t*)block_alloc_.allocate(len, 0);
    }

    /*! @brief Frees GAP block using bot block allocator (BA)
    */
    void free_gap_block(bm::gap_word_t*   block,
                        const bm::gap_word_t* glevel_len)
    {
        unsigned len = bm::gap_capacity(block, glevel_len);
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
    BA                     block_alloc_;
    PA                     ptr_alloc_;
    allocator_pool_type*   alloc_pool_p_;
};


typedef mem_alloc<block_allocator, ptr_allocator, alloc_pool<block_allocator, ptr_allocator> > standard_allocator;


/// Aligned malloc (unlike classic malloc it throws bad_alloc exception)
///
/// @internal
inline
void* aligned_new_malloc(size_t size)
{
    void* ptr;

#ifdef BM_ALLOC_ALIGN
#ifdef _MSC_VER
    ptr = ::_aligned_malloc(size, BM_ALLOC_ALIGN);
#else
    ptr = ::_mm_malloc(size, BM_ALLOC_ALIGN);
#endif
#else
    ptr = ::malloc(size);
#endif
    if (!ptr)
    {
#ifndef BM_NO_STL
    throw std::bad_alloc();
#else
    BM_THROW(BM_ERR_BADALLOC);
#endif
    }
    return ptr;
}

/// Aligned free (safe to pass NULL ptr)
///
/// @internal
inline
void aligned_free(void* ptr)
{
    if (!ptr)
        return;
#ifdef BM_ALLOC_ALIGN
# ifdef _MSC_VER
    ::_aligned_free(ptr);
#else
    ::_mm_free(ptr);
# endif
#else
    ::free(ptr);
#endif
}


} // namespace libbm


#endif
