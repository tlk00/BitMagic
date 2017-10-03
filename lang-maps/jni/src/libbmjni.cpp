#include "libbm.h"
#include "try_throw_catch.h"
static jmp_buf ex_buf__;

#define BM_NO_STL
#define BMALLOC__H__INCLUDED__

#include "bmconst.h"
#include "jnialloc.h"
#include "bm.h"

#include <new>

typedef libbm::standard_allocator              TBM_Alloc;
typedef bm::bvector<libbm::standard_allocator> TBM_bvector;

// include standard C-library implementation with all allocation should now
// be reimplemented to use Java specific allocator

#include "libbm_impl.cpp"



