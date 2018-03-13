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

#include "io_bitmagic_core_SimdUtil0.h"
#include "libbmcpuid.h"
#include <stdio.h>

#define BMJNI_SSE 		  "bmjni_sse"
#define BMJNI_SSE2 		  "bmjni_sse2"
#define BMJNI_SSE3		  "bmjni_sse3"
#define BMJNI_SSE4_1  	"bmjni_sse4_1"
#define BMJNI_SSE4_2  	"bmjni_sse4_2"
#define BMJNI_AVX     	"bmjni_avx"
#define BMJNI_AVX2    	"bmjni_avx2"
#define BMJNI_AVX512F 	"bmjni_avx512f"

const char* BM_simd_libname(void) {
  unsigned simd = BM_x86_simd();
  printf("Simd code: %x\n", simd);
  if (simd & SIMD_AVX512F) 
    return BMJNI_AVX512F;
  else if (simd & SIMD_AVX2)
    return BMJNI_AVX2;
  else if (simd & SIMD_AVX)
    return BMJNI_AVX;
  else if (simd & SIMD_SSE4_2)
    return BMJNI_SSE4_2;
  else if (simd & SIMD_SSE4_1)
    return BMJNI_SSE4_1;
  else if (simd & SIMD_SSE3)
    return BMJNI_SSE3;
  else if (simd & SIMD_SSE2)
    return BMJNI_SSE2;
  else if (simd & SIMD_SSE)
    return BMJNI_SSE;
  else
    return "Unknown simd code, no library name found";  
}

/*
 * Class:     io_bitmagic_core_SimdUtil0
 * Method:    getLibName0
 * Signature: ()Ljava/lang/String;
 */
JNIEXPORT jstring JNICALL Java_io_bitmagic_core_SimdUtil0_getLibName0
  (JNIEnv *env, jclass clazz) {
  const char *libName = BM_simd_libname();
  return (*env)->NewStringUTF(env, libName);
}



