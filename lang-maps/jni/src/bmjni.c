#include "io_bitmagic_BVector0.h"
#include "libbm.h"

#include <string.h>
#include <stdio.h>

#define exec(x) \
  switch(x) { \
  case BM_ERR_BADALLOC: \
  { \
    jclass ex = (*env)->FindClass(env, "java/lang/RuntimeException"); \
    (*env)->ThrowNew(env, ex,BM_ERR_BADALLOC_MSG " in native call: " #x); \
    break; \
  } \
  case BM_ERR_BADARG: \
  { \
    jclass ex = (*env)->FindClass(env, "java/lang/IllegalArgumentException"); \
    (*env)->ThrowNew(env, ex, BM_ERR_BADARG_MSG " in: " #x); \
    break; \
  } \
  case BM_ERR_RANGE: \
  { \
    jclass ex = (*env)->FindClass(env, "java/lang/IndexOutOfBoundsException"); \
    (*env)->ThrowNew(env, ex, BM_ERR_RANGE_MSG " in: " #x); \
    break; \
  } \
  default: \
    break; \
  }

JNIEXPORT void JNICALL Java_io_bitmagic_BVector0_init0
  (JNIEnv *env, jobject obj, jlong ptr) {
  exec(BM_init((void*)ptr));
}

/*
 * Class:     io_bitmagic_BVector0
 * Method:    create0
 * Signature: (IJ)J
 */
JNIEXPORT jlong JNICALL Java_io_bitmagic_BVector0_create0
  (JNIEnv *env, jobject obj, jint x, jlong size) {
  BM_BVHANDLE ptr;
  exec(BM_bvector_construct(&ptr, size));
  return (jlong)ptr;
}

/*
 * Class:     io_bitmagic_BVector0
 * Method:    clone0
 * Signature: (J)J
 */
JNIEXPORT jlong JNICALL Java_io_bitmagic_BVector0_copy0
  (JNIEnv *env, jobject obj, jlong ptr) {
  BM_BVHANDLE cptr;
  exec(BM_bvector_construct_copy(&cptr, (BM_BVHANDLE)ptr));
  return (jlong)cptr;
}

/*
 * Class:     io_bitmagic_BVector0
 * Method:    dispose0
 * Signature: (J)V
 */
JNIEXPORT void JNICALL Java_io_bitmagic_BVector0_dispose0
  (JNIEnv *env, jobject obj, jlong ptr) {
  exec(BM_bvector_free((BM_BVHANDLE)ptr));
}

/*
 * Class:     io_bitmagic_BVector0
 * Method:    version0
 * Signature: ()Ljava/lang/String;
 */
JNIEXPORT jstring JNICALL Java_io_bitmagic_BVector0_version0
(JNIEnv *env, jobject obj) {
  int major;
  int minor;
  int patch;
  BM_version(&major, &minor, &patch);
  
  char version[64];
  snprintf(version, 64, "%d.%d.%d", major, minor, patch);
  return (*env)->NewStringUTF(env, version);
}

/*
 * Class:     io_bitmagic_BVector0
 * Method:    copyright0
 * Signature: ()Ljava/lang/String;
 */
JNIEXPORT jstring JNICALL Java_io_bitmagic_BVector0_copyright0
(JNIEnv *env, jobject obj) {
  return (*env)->NewStringUTF(env, BM_version(0, 0, 0));
}

/*
 * Class:     io_bitmagic_BVector0
 * Method:    set0
 * Signature: (JJZ)V
 */
JNIEXPORT void JNICALL Java_io_bitmagic_BVector0_set0
  (JNIEnv *env, jobject obj, jlong ptr, jlong idx, jboolean bit) {
  exec(BM_bvector_set_bit((BM_BVHANDLE)ptr, (unsigned int)idx, bit ? 1 : 0));
}

/*
 * Class:     io_bitmagic_BVector0
 * Method:    get0
 * Signature: (JJ)Z
 */
JNIEXPORT jboolean JNICALL Java_io_bitmagic_BVector0_get0
  (JNIEnv *env, jobject obj, jlong ptr, jlong idx) {
  int ret;
  exec(BM_bvector_get_bit((BM_BVHANDLE)ptr, (unsigned int)idx, &ret));
  return (jboolean)ret;
}

/*
* Class:     io_bitmagic_BVector0
* Method:    getSize0
* Signature: (J)J
*/
JNIEXPORT jlong JNICALL Java_io_bitmagic_BVector0_getSize0
(JNIEnv *env, jclass obj, jlong ptr) {
  unsigned int size;
  exec(BM_bvector_get_size((BM_BVHANDLE)ptr, &size));
  return (jlong)size;
}

/*
* Class:     io_bitmagic_BVector0
* Method:    setSize0
* Signature: (JJ)V
*/
JNIEXPORT void JNICALL Java_io_bitmagic_BVector0_setSize0
(JNIEnv *env, jclass obj, jlong ptr, jlong size) {
  exec(BM_bvector_set_size((BM_BVHANDLE)ptr, size));
}
