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

#include "libbmcpuid.h"
#ifdef _MSC_VER
#include <intrin.h>
#endif

// -------------------------------------------------------------------
//
// code credit for CPU caps identification:
// https://attractivechaos.wordpress.com/2017/09/04/on-cpu-dispatch/
//
// -------------------------------------------------------------------

unsigned BM_x86_simd(void)
{
  unsigned eax, ebx, ecx, edx, flag = 0;
#ifdef _MSC_VER
  int cpuid[4];
  __cpuid(cpuid, 1);
  eax = cpuid[0], ebx = cpuid[1], ecx = cpuid[2], edx = cpuid[3];
#else
  asm volatile("cpuid" : "=a" (eax), "=b" (ebx), "=c" (ecx), "=d" (edx) : "a" (1));
#endif
  if (edx>>25&1) flag |= SIMD_SSE;
  if (edx>>26&1) flag |= SIMD_SSE2;
  if (ecx>>0 &1) flag |= SIMD_SSE3;
  if (ecx>>19&1) flag |= SIMD_SSE4_1;
  if (ecx>>20&1) flag |= SIMD_SSE4_2;
  if (ecx>>28&1) flag |= SIMD_AVX;
  if (ebx>>5 &1) flag |= SIMD_AVX2;
  if (ebx>>16&1) flag |= SIMD_AVX512F;
  return flag;
}


