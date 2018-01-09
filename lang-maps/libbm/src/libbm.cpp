
#include "libbm.h"
#include "try_throw_catch.h"
static jmp_buf ex_buf__;

#define BM_NO_STL
#define BM_NO_CXX11
#define BMALLOC__H__INCLUDED__

#define BM_ASSERT_THROW(x, xerrcode) if (!(x)) BM_THROW( xerrcode )


#include "bmdef.h"
#include "bmconst.h"
#include "bmcalloc.h"
#include "bm.h"

#include <new>

typedef libbm::standard_allocator              TBM_Alloc;
typedef bm::bvector<libbm::standard_allocator> TBM_bvector;


#include "libbm_impl.cpp"



