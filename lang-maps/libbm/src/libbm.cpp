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

#include "libbm.h"
#include "try_throw_catch.h"

#if defined(_MSC_VER)
__declspec(thread) jmp_buf ex_buf__;
#else
__thread jmp_buf ex_buf__;
#endif

#define BM_NO_STL
#define BM_NO_CXX11
#define BMALLOC__H__INCLUDED__

#define BM_ASSERT_THROW(x, xerrcode) if (!(x)) BM_THROW( xerrcode )


#include "bmdef.h"
#include "bmconst.h"
#include "bmsimd.h"
#include "bmcalloc.h"
#include "bm.h"

#include <new>

typedef libbm::standard_allocator              TBM_Alloc;
typedef bm::bvector<libbm::standard_allocator> TBM_bvector;


#include "libbm_impl.cpp"



