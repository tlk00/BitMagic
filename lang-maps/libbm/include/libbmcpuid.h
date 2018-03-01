#ifndef LIBBMCPUID_INCLUDED_H__
#define LIBBMCPUID_INCLUDED_H__
/*
BitMagic Library License

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

#include <stddef.h>

#define SIMD_SSE     0x1
#define SIMD_SSE2    0x2
#define SIMD_SSE3    0x4
#define SIMD_SSE4_1  0x8
#define SIMD_SSE4_2  0x10
#define SIMD_AVX     0x20
#define SIMD_AVX2    0x40
#define SIMD_AVX512F 0x80

#ifdef _WIN32
#ifdef BMDLLEXPORTS
#    define BM_API_EXPORT __declspec(dllexport)
#else
#    define BM_API_EXPORT __declspec(dllimport)
#endif
#else
#   define BM_API_EXPORT
#endif

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

BM_API_EXPORT unsigned BM_x86_simd(void);

#ifdef __cplusplus
}
#endif /* __cplusplus */


#endif
