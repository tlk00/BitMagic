#include "org_tlk00_bitmagic_BVector0.h"
#include "libbm.h"

#define exec(x) \
  switch(x) { \
  case BM_ERR_BADALLOC: \
  { \
    jclass ex = env->FindClass("java/lang/RuntimeException"); \
    env->ThrowNew(ex,"Bad alloc in native call: " #x); \
    break; \
  } \
  case BM_ERR_BADARG: \
  { \
    jclass ex = env->FindClass("java/lang/IllegalArgumentException"); \
    env->ThrowNew(ex,"One or more arguments are invalid in: " #x); \
    break; \
  } \
  default: \
    break; \
  }
  
    
JNIEXPORT void JNICALL Java_org_tlk00_bitmagic_BVector0_init0
  (JNIEnv *env, jobject, jlong ptr) {
  exec(BM_init((void*)ptr));
}

/*
 * Class:     org_tlk00_bitmagic_BVector0
 * Method:    create0
 * Signature: (IJ)J
 */
JNIEXPORT jlong JNICALL Java_org_tlk00_bitmagic_BVector0_create0
  (JNIEnv *env, jobject, jint, jlong size) {
  BM_BVHANDLE ptr;
  exec(BM_bvector_construct(&ptr, size));
  return (jlong)ptr;
}

/*
 * Class:     org_tlk00_bitmagic_BVector0
 * Method:    clone0
 * Signature: (J)J
 */
JNIEXPORT jlong JNICALL Java_org_tlk00_bitmagic_BVector0_clone0
  (JNIEnv *env, jobject, jlong ptr) {
  return ptr;
}

/*
 * Class:     org_tlk00_bitmagic_BVector0
 * Method:    dispose0
 * Signature: (J)V
 */
JNIEXPORT void JNICALL Java_org_tlk00_bitmagic_BVector0_dispose0
  (JNIEnv *env, jobject, jlong ptr) {
  exec(BM_bvector_free((BM_BVHANDLE)ptr));
}

/*
 * Class:     org_tlk00_bitmagic_BVector0
 * Method:    set0
 * Signature: (JJZ)V
 */
JNIEXPORT void JNICALL Java_org_tlk00_bitmagic_BVector0_set0
  (JNIEnv *env, jobject, jlong ptr, jlong idx, jboolean bit) {
  exec(BM_bvector_set_bit((BM_BVHANDLE)ptr, (unsigned int)idx, bit ? 1 : 0));
}

/*
 * Class:     org_tlk00_bitmagic_BVector0
 * Method:    get0
 * Signature: (JJ)Z
 */
JNIEXPORT jboolean JNICALL Java_org_tlk00_bitmagic_BVector0_get0
  (JNIEnv *env, jobject, jlong ptr, jlong idx) {
  unsigned int ret;
  exec(BM_bvector_get_bit((BM_BVHANDLE)ptr, (unsigned int)idx, &ret));
  return (jboolean)ret;
}
