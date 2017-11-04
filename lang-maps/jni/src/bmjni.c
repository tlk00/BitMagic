#include "io_bitmagic_core_BVector0.h"
#include "io_bitmagic_core_BVIterator0.h"
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

JNIEXPORT void JNICALL Java_io_bitmagic_core_BVector0_init0
  (JNIEnv *env, jobject obj, jlong ptr) {
  exec(BM_init((void*)ptr));
}

/*
 * Class:     io_bitmagic_core_BVector0
 * Method:    create0
 * Signature: (IJ)J
 */
JNIEXPORT jlong JNICALL Java_io_bitmagic_core_BVector0_create0
  (JNIEnv *env, jobject obj, jint x, jlong size) {
  BM_BVHANDLE ptr;
  exec(BM_bvector_construct(&ptr, size));
  return (jlong)ptr;
}

/*
 * Class:     io_bitmagic_core_BVector0
 * Method:    clone0
 * Signature: (J)J
 */
JNIEXPORT jlong JNICALL Java_io_bitmagic_core_BVector0_copy0
  (JNIEnv *env, jobject obj, jlong ptr) {
  BM_BVHANDLE cptr;
  exec(BM_bvector_construct_copy(&cptr, (BM_BVHANDLE)ptr));
  return (jlong)cptr;
}

/*
 * Class:     io_bitmagic_core_BVector0
 * Method:    dispose0
 * Signature: (J)V
 */
JNIEXPORT void JNICALL Java_io_bitmagic_core_BVector0_dispose0
  (JNIEnv *env, jobject obj, jlong ptr) {
  exec(BM_bvector_free((BM_BVHANDLE)ptr));
}

/*
 * Class:     io_bitmagic_core_BVector0
 * Method:    version0
 * Signature: ()Ljava/lang/String;
 */
JNIEXPORT jstring JNICALL Java_io_bitmagic_core_BVector0_version0
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
 * Class:     io_bitmagic_core_BVector0
 * Method:    copyright0
 * Signature: ()Ljava/lang/String;
 */
JNIEXPORT jstring JNICALL Java_io_bitmagic_core_BVector0_copyright0
(JNIEnv *env, jobject obj) {
  return (*env)->NewStringUTF(env, BM_version(0, 0, 0));
}

/*
* Class:     io_bitmagic_core_BVector0
* Method:    getSize0
* Signature: (J)J
*/
JNIEXPORT jlong JNICALL Java_io_bitmagic_core_BVector0_getSize0
(JNIEnv *env, jclass obj, jlong ptr) {
  unsigned int size;
  exec(BM_bvector_get_size((BM_BVHANDLE)ptr, &size));
  return (jlong)size;
}

/*
* Class:     io_bitmagic_core_BVector0
* Method:    setSize0
* Signature: (JJ)V
*/
JNIEXPORT void JNICALL Java_io_bitmagic_core_BVector0_setSize0
(JNIEnv *env, jclass obj, jlong ptr, jlong size) {
  exec(BM_bvector_set_size((BM_BVHANDLE)ptr, size));
}

/*
 * Class:     io_bitmagic_core_BVector0
 * Method:    set0
 * Signature: (JJZ)V
 */
JNIEXPORT void JNICALL Java_io_bitmagic_core_BVector0_set0
(JNIEnv *env, jobject obj, jlong ptr, jlong idx, jboolean bit) {
  exec(BM_bvector_set_bit((BM_BVHANDLE)ptr, (unsigned int)idx, bit ? 1 : 0));
}

/*
 * Class:     io_bitmagic_core_BVector0
 * Method:    setConditional0
 * Signature: (JJZZ)Z
 */
JNIEXPORT jboolean JNICALL Java_io_bitmagic_core_BVector0_setConditional0
(JNIEnv *env, jclass obj, jlong ptr, jlong idx, jboolean bit, jboolean condition) {
  int changed;
  exec(BM_bvector_set_bit_conditional((BM_BVHANDLE)ptr, (unsigned int)idx, bit ? 1 : 0, condition ? 1 : 0, &changed));
  return (jboolean)changed;
}

/*
 * Class:     io_bitmagic_core_BVector0
 * Method:    flip0
 * Signature: (JJ)V
 */
JNIEXPORT void JNICALL Java_io_bitmagic_core_BVector0_flip0
(JNIEnv *env, jclass obj, jlong ptr, jlong idx) {
  exec(BM_bvector_flip_bit((BM_BVHANDLE)ptr, (unsigned int)idx))
}

/*
 * Class:     io_bitmagic_core_BVector0
 * Method:    setAll0
 * Signature: (J)V
 */
JNIEXPORT void JNICALL Java_io_bitmagic_core_BVector0_setAll0
(JNIEnv *env, jclass obj, jlong ptr) {
  exec(BM_bvector_set((BM_BVHANDLE)ptr));
}

/*
 * Class:     io_bitmagic_core_BVector0
 * Method:    setRange0
 * Signature: (JJJZ)V
 */
JNIEXPORT void JNICALL Java_io_bitmagic_core_BVector0_setRange0
(JNIEnv *env, jclass obj, jlong ptr, jlong left, jlong right, jboolean bit) {
  exec(BM_bvector_set_range((BM_BVHANDLE)ptr, (unsigned int)left, (unsigned int)right, bit ? 1 : 0));
}

/*
 * Class:     io_bitmagic_core_BVector0
 * Method:    invert0
 * Signature: (J)V
 */
JNIEXPORT void JNICALL Java_io_bitmagic_core_BVector0_invert0
(JNIEnv *env, jclass obj, jlong ptr) {
  exec(BM_bvector_invert((BM_BVHANDLE)ptr));
}

/*
 * Class:     io_bitmagic_core_BVector0
 * Method:    clear0
 * Signature: (JI)V
 */
JNIEXPORT void JNICALL Java_io_bitmagic_core_BVector0_clear0
(JNIEnv *env, jclass obj, jlong ptr, jint free_mem) {
  exec(BM_bvector_clear((BM_BVHANDLE)ptr, free_mem));
}

/*
 * Class:     io_bitmagic_core_BVector0
 * Method:    extract0
 * Signature: (JJ)J
 */
JNIEXPORT jlong JNICALL Java_io_bitmagic_core_BVector0_extract0
(JNIEnv *env, jclass obj, jlong ptr, jlong start) {
  unsigned int next;
  exec(BM_bvector_extract_next((BM_BVHANDLE)ptr, (unsigned int)start, &next));
  return next == 0 ? -1 : (jlong)next;
}

/*
 * Class:     io_bitmagic_core_BVector0
 * Method:    get0
 * Signature: (JJ)Z
 */
JNIEXPORT jboolean JNICALL Java_io_bitmagic_core_BVector0_get0
(JNIEnv *env, jobject obj, jlong ptr, jlong idx) {
  int ret;
  exec(BM_bvector_get_bit((BM_BVHANDLE)ptr, (unsigned int)idx, &ret));
  return (jboolean)ret;
}

/*
 * Class:     io_bitmagic_core_BVector0
 * Method:    count0
 * Signature: (J)J
 */
JNIEXPORT jlong JNICALL Java_io_bitmagic_core_BVector0_count0
(JNIEnv *env, jclass obj, jlong ptr) {
  unsigned int count;
  exec(BM_bvector_count((BM_BVHANDLE)ptr, &count));
  return (jlong)count;
}

/*
 * Class:     io_bitmagic_core_BVector0
 * Method:    countInRange0
 * Signature: (JJJ)J
 */
JNIEXPORT jlong JNICALL Java_io_bitmagic_core_BVector0_countInRange0
(JNIEnv *env, jclass obj, jlong ptr, jlong left, jlong right) {
  unsigned int count;
  exec(BM_bvector_count_range((BM_BVHANDLE)ptr, (unsigned int)left, (unsigned int)right, &count));
  return (jlong)count;
}

/*
 * Class:     io_bitmagic_core_BVector0
 * Method:    nonEmpty0
 * Signature: (J)Z
 */
JNIEXPORT jboolean JNICALL Java_io_bitmagic_core_BVector0_nonEmpty0
(JNIEnv *env, jclass obj, jlong ptr) {
  int val;
  exec(BM_bvector_any((BM_BVHANDLE)ptr, &val));
  return (jboolean)val;
}

/*
 * Class:     io_bitmagic_core_BVector0
 * Method:    indexOf0
 * Signature: (JJ)J
 */
JNIEXPORT jlong JNICALL Java_io_bitmagic_core_BVector0_indexOf0
(JNIEnv *env, jclass obj, jlong ptr, jlong start) {
  int found;
  unsigned int pos;
  exec(BM_bvector_find((BM_BVHANDLE)ptr, (unsigned int)start, &pos, &found));
  if (found)
    return (jlong)pos;
  else
    return (jlong)-1;
}

/*
 * Class:     io_bitmagic_core_BVector0
 * Method:    compare0
 * Signature: (JJ)I
 */
JNIEXPORT jint JNICALL Java_io_bitmagic_core_BVector0_compare0
(JNIEnv *env, jclass obj, jlong ptr1, jlong ptr2) {
  int comp;
  exec(BM_bvector_compare((BM_BVHANDLE)ptr1, (BM_BVHANDLE)ptr2, &comp));
  return (jint)comp;
}

/*
 * Class:     io_bitmagic_core_BVector0
 * Method:    optimize0
 * Signature: (JI)Lio/bitmagic/BitVectorStat;
 */
JNIEXPORT jobject JNICALL Java_io_bitmagic_core_BVector0_optimize0
(JNIEnv *env, jclass obj, jlong ptr, jint opt_mode) {
  struct BM_bvector_statistics stat;
  jobject object;
  jmethodID constructor;
  jobject cls;

  exec(BM_bvector_optimize((BM_BVHANDLE)ptr, opt_mode, &stat));
  cls = (*env)->FindClass(env, "io/bitmagic/core/BitVectorStat");
  constructor = (*env)->GetMethodID(env, cls, "<init>", "(JJJJ)V");
  object = (*env)->NewObject(env, cls, constructor, (jlong)stat.bit_blocks, (jlong)stat.gap_blocks, (jlong)stat.max_serialize_mem, (jlong)stat.memory_used);
  return object;
}

/*
 * Class:     io_bitmagic_core_BVector0
 * Method:    calcStat0
 * Signature: (J)Lio/bitmagic/BitVectorStat;
 */
JNIEXPORT jobject JNICALL Java_io_bitmagic_core_BVector0_calcStat0
(JNIEnv *env, jclass obj, jlong ptr) {
  struct BM_bvector_statistics stat;
  jobject object;
  jmethodID constructor;
  jobject cls;
  
  exec(BM_bvector_calc_stat((BM_BVHANDLE)ptr, &stat));
  cls = (*env)->FindClass(env, "io/bitmagic/core/BitVectorStat");
  constructor = (*env)->GetMethodID(env, cls, "<init>", "(JJJJ)V");
  object = (*env)->NewObject(env, cls, constructor, (jlong)stat.bit_blocks, (jlong)stat.gap_blocks, (jlong)stat.max_serialize_mem, (jlong)stat.memory_used);
  return object;
}

/*
 * Class:     io_bitmagic_core_BVector0
 * Method:    operation0
 * Signature: (JJI)V
 */
JNIEXPORT void JNICALL Java_io_bitmagic_core_BVector0_operation0
(JNIEnv *env, jclass obj, jlong dst, jlong src, jint op) {
  exec(BM_bvector_combine_operation((BM_BVHANDLE)dst, (BM_BVHANDLE)src, op));
}

/*
 * Class:     io_bitmagic_core_BVector0
 * Method:    and0
 * Signature: (JJ)V
 */
JNIEXPORT void JNICALL Java_io_bitmagic_core_BVector0_and0
(JNIEnv *env, jclass obj, jlong dst, jlong src) {
  exec(BM_bvector_combine_AND((BM_BVHANDLE)dst, (BM_BVHANDLE)src));
}

/*
 * Class:     io_bitmagic_core_BVector0
 * Method:    or0
 * Signature: (JJ)V
 */
JNIEXPORT void JNICALL Java_io_bitmagic_core_BVector0_or0
(JNIEnv *env, jclass obj, jlong dst, jlong src) {
  exec(BM_bvector_combine_OR((BM_BVHANDLE)dst, (BM_BVHANDLE)src));
}

/*
 * Class:     io_bitmagic_core_BVector0
 * Method:    sub0
 * Signature: (JJ)V
 */
JNIEXPORT void JNICALL Java_io_bitmagic_core_BVector0_sub0
(JNIEnv *env, jclass obj, jlong dst, jlong src) {
  exec(BM_bvector_combine_SUB((BM_BVHANDLE)dst, (BM_BVHANDLE)src));
}

/*
 * Class:     io_bitmagic_core_BVector0
 * Method:    xor0
 * Signature: (JJ)V
 */
JNIEXPORT void JNICALL Java_io_bitmagic_core_BVector0_xor0
(JNIEnv *env, jclass obj, jlong dst, jlong src) {
  exec(BM_bvector_combine_XOR((BM_BVHANDLE)dst, (BM_BVHANDLE)src));
}

/***************************************** BitVector Iterator ****************************************/
/*
* Class:     io_bitmagic_core_BVIterator0
* Method:    create0
* Signature: (J)J
*/
JNIEXPORT jlong JNICALL Java_io_bitmagic_core_BVIterator0_create0
(JNIEnv *env, jclass obj, jlong bvPtr) {
  BM_BVEHANDLE eh;
  exec(BM_bvector_enumerator_construct((BM_BVHANDLE)bvPtr, &eh));
  return (jlong)eh;
}

/*
* Class:     io_bitmagic_core_BVIterator0
* Method:    dispose0
* Signature: (J)V
*/
JNIEXPORT void JNICALL Java_io_bitmagic_core_BVIterator0_dispose0
(JNIEnv *env, jclass obj, jlong ptr) {
  exec(BM_bvector_enumerator_free((BM_BVEHANDLE)ptr));
}

/*
* Class:     io_bitmagic_core_BVIterator0
* Method:    isValid0
* Signature: (J)Z
*/
JNIEXPORT jboolean JNICALL Java_io_bitmagic_core_BVIterator0_isValid0
(JNIEnv *env, jobject obj, jlong ptr) {
  int valid;
  exec(BM_bvector_enumerator_is_valid((BM_BVEHANDLE)ptr, &valid));
  return (jboolean)(valid != 0);
}

/*
* Class:     io_bitmagic_core_BVIterator0
* Method:    get0
* Signature: (J)J
*/
JNIEXPORT jlong JNICALL Java_io_bitmagic_core_BVIterator0_get0
(JNIEnv *env, jobject obj, jlong ptr) {
  unsigned int value;
  exec(BM_bvector_enumerator_get_value((BM_BVEHANDLE)ptr, &value));
  return (jlong)value;
}

/*
* Class:     io_bitmagic_core_BVIterator0
* Method:    next0
* Signature: (J)Z
*/
JNIEXPORT jboolean JNICALL Java_io_bitmagic_core_BVIterator0_next0
(JNIEnv *env, jobject obj, jlong ptr) {
  unsigned int value;
  int valid;
  exec(BM_bvector_enumerator_next((BM_BVEHANDLE)ptr, &valid, &value));
  return (jboolean)(valid != 0);
}
