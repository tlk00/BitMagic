// Copyright(c) 2002-2009 Anatoliy Kuznetsov(anatoliy_kuznetsov at yahoo.com)
//
// BM library internal header
//
// Set all required preprocessor defines

#include <climits>

// Incorporate appropriate tuneups when the NCBI C++ Toolkit's core
// headers have been included.
//
#ifdef NCBI_ASSERT
#  define BM_ASSERT _ASSERT

#  ifdef HAVE_RESTRICT_CXX
#    define BM_HASRESTRICT
#    define BMRESTRICT NCBI_RESTRICT
#  endif

#  if defined(NCBI_FORCEINLINE)  &&  \
    ( !defined(NCBI_COMPILER_GCC)  ||  NCBI_COMPILER_VERSION >= 400  || \
      defined(__OPTIMIZE__))
#    define BM_HASFORCEINLINE
#    define BMFORCEINLINE NCBI_FORCEINLINE
#  endif
#endif


// macro to define/undefine unaligned memory access (x86, PowerPC)
//
#if defined(__i386) || defined(__x86_64) || defined(__ppc__) || \
	defined(__ppc64__) || defined(_M_IX86) || defined(_M_AMD64) || \
    defined(_M_IX86) || defined(_M_AMD64) || defined(_M_X64) || \
    (defined(_M_MPPC) && !defined(BM_FORBID_UNALIGNED_ACCESS))
#define BM_UNALIGNED_ACCESS_OK 1
#endif

#if defined(_M_IX86) || defined(_M_AMD64) || defined(_M_X64) || \
    defined(__i386) || defined(__x86_64) || defined(_M_AMD64) || \
    defined(BMSSE2OPT) || defined(BMSSE42OPT)
#define BM_x86
#endif

// disable 'register' keyword, which is obsolete in C++11
//
#ifndef BMREGISTER
# define BMREGISTER
#endif


// Enable MSVC 8.0 (2005) specific optimization options
//
#if(_MSC_VER >= 1400)

#  define BM_HASFORCEINLINE
#  ifndef BMRESTRICT
#    define BMRESTRICT __restrict
#  endif
#endif 

#ifdef __GNUG__

#  ifndef BMRESTRICT
#    define BMRESTRICT __restrict__
#  endif

#  ifdef __OPTIMIZE__
#    define BM_NOASSERT
#  endif
#endif

#ifndef BM_ASSERT

# ifndef BM_NOASSERT
#  include <cassert>
#  define BM_ASSERT assert
# else
#  ifndef BM_ASSERT
#    define BM_ASSERT(x)
#  endif
# endif

#endif


#define FULL_BLOCK_REAL_ADDR bm::all_set<true>::_block._p
#define FULL_BLOCK_FAKE_ADDR bm::all_set<true>::_block._p_fullp
#define BLOCK_ADDR_SAN(addr) (addr == FULL_BLOCK_FAKE_ADDR) ? FULL_BLOCK_REAL_ADDR : addr
#define IS_VALID_ADDR(addr) bm::all_set<true>::is_valid_block_addr(addr)
#define IS_FULL_BLOCK(addr) bm::all_set<true>::is_full_block(addr)
#define IS_EMPTY_BLOCK(addr) bool(addr == 0)

// Macro definitions to manipulate bits in pointers
// This trick is based on the fact that pointers allocated by malloc are
// aligned and bit 0 is never set. It means we are safe to use it.
// BM library keeps GAP flag in pointer.



# if ULONG_MAX != 0xffffffff || defined(_WIN64)  // 64-bit

#  define BMPTR_SETBIT0(ptr)   ( ((bm::id64_t)ptr) | 1 )
#  define BMPTR_CLEARBIT0(ptr) ( ((bm::id64_t)ptr) & ~(bm::id64_t)1 )
#  define BMPTR_TESTBIT0(ptr)  ( ((bm::id64_t)ptr) & 1 )

# else // 32-bit

#  define BMPTR_SETBIT0(ptr)   ( ((bm::id_t)ptr) | 1 )
#  define BMPTR_CLEARBIT0(ptr) ( ((bm::id_t)ptr) & ~(bm::id_t)1 )
#  define BMPTR_TESTBIT0(ptr)  ( ((bm::id_t)ptr) & 1 )

# endif

# define BMGAP_PTR(ptr) ((bm::gap_word_t*)BMPTR_CLEARBIT0(ptr))
# define BMSET_PTRGAP(ptr) ptr = (bm::word_t*)BMPTR_SETBIT0(ptr)
# define BM_IS_GAP(ptr) bool(BMPTR_TESTBIT0(ptr)!=0)





#ifdef BM_HASRESTRICT
# ifndef BMRESTRICT
#  define BMRESTRICT restrict
# endif
#else
# ifndef BMRESTRICT
#   define BMRESTRICT
# endif
#endif


#ifdef BM_HASFORCEINLINE
# ifndef BMFORCEINLINE
#  define BMFORCEINLINE __forceinline
# endif
#else
# define BMFORCEINLINE inline
#endif


// --------------------------------
// SSE optmization
//

#if !(defined(BMSSE2OPT) || defined(BMSSE42OPT)) 

# ifndef BM_SET_MMX_GUARD
#  define BM_SET_MMX_GUARD
# endif

#define BM_ALIGN16 
#define BM_ALIGN16ATTR

#else  

# ifndef BM_SET_MMX_GUARD
#  define BM_SET_MMX_GUARD  sse_empty_guard  bm_mmx_guard_;
# endif

#ifdef _MSC_VER

#ifndef BM_ALIGN16
#  define BM_ALIGN16 __declspec(align(16))
#  define BM_ALIGN16ATTR
#endif

# else // GCC

#ifndef BM_ALIGN16
#  define BM_ALIGN16
#  define BM_ALIGN16ATTR __attribute__((aligned(16)))
#endif

#endif

#endif

/*! 
	Define calculates number of 1 bits in 32-bit word.
    @ingroup bitfunc 
*/
#ifndef BM_INCWORD_BITCOUNT

#ifdef BMSSE42OPT

# define BM_INCWORD_BITCOUNT(cnt, w) cnt += _mm_popcnt_u32(w);

#else

# define BM_INCWORD_BITCOUNT(cnt, w) cnt += \
     bm::bit_count_table<true>::_count[(unsigned char)(w)] + \
     bm::bit_count_table<true>::_count[(unsigned char)((w) >> 8)] + \
     bm::bit_count_table<true>::_count[(unsigned char)((w) >> 16)] + \
     bm::bit_count_table<true>::_count[(unsigned char)((w) >> 24)];

#endif

#endif


