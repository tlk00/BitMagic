#include "io_bitmagic_core_BVector0.h"
#include "io_bitmagic_core_BVIterator0.h"
#include "libbm.h"

#include <string.h>
#include <stdio.h>

#define exec(x) \
  switch(x) { \
  case BM_OK: \
    break; \
  case BM_ERR_BADALLOC: \
  { \
    jclass ex = (*env)->FindClass(env, "java/lang/RuntimeException"); \
    (*env)->ThrowNew(env, ex, BM_ERR_BADALLOC_MSG " in native call: " #x); \
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
  case BM_ERR_CPU: \
  { \
    jclass ex = (*env)->FindClass(env, "java/lang/RuntimeException"); \
    (*env)->ThrowNew(env, ex, BM_ERR_CPU_MSG " in: " #x); \
    break; \
  } \
  default: \
  { \
    jclass ex = (*env)->FindClass(env, "java/lang/RuntimeException"); \
    (*env)->ThrowNew(env, ex, "Unknown error in native call: " #x); \
    break; \
  }}

#define null_out_of_memory(x) \
  if ((x) == NULL) { \
    jclass ex = (*env)->FindClass(env, "java/lang/OutOfMemoryError"); \
    (*env)->ThrowNew(env, ex, "Out of memory error in native call andArr0()"); \
    return; \
  }

JNIEXPORT void JNICALL Java_io_bitmagic_core_BVector0_init0
  (JNIEnv *env, jclass clazz, jlong ptr) {
  exec(BM_init((void*)ptr));
}

/*
 * Class:     io_bitmagic_core_BVector0
 * Method:    create0
 * Signature: (IJ)J
 */
JNIEXPORT jlong JNICALL Java_io_bitmagic_core_BVector0_create0
  (JNIEnv *env, jclass clazz, jint strategy, jlong size) {
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
  (JNIEnv *env, jclass clazz, jlong ptr) {
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
  (JNIEnv *env, jclass clazz, jlong ptr) {
  exec(BM_bvector_free((BM_BVHANDLE)ptr));
}

/*
 * Class:     io_bitmagic_core_BVector0
 * Method:    version0
 * Signature: ()Ljava/lang/String;
 */
JNIEXPORT jstring JNICALL Java_io_bitmagic_core_BVector0_version0
(JNIEnv *env, jclass clazz) {
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
(JNIEnv *env, jclass clazz) {
  return (*env)->NewStringUTF(env, BM_version(0, 0, 0));
}

/*
* Class:     io_bitmagic_core_BVector0
* Method:    getSize0
* Signature: (J)J
*/
JNIEXPORT jlong JNICALL Java_io_bitmagic_core_BVector0_getSize0
(JNIEnv *env, jclass clazz, jlong ptr) {
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
(JNIEnv *env, jclass clazz, jlong ptr, jlong size) {
  exec(BM_bvector_set_size((BM_BVHANDLE)ptr, size));
}

/*
 * Class:     io_bitmagic_core_BVector0
 * Method:    set0
 * Signature: (JJZ)V
 */
JNIEXPORT void JNICALL Java_io_bitmagic_core_BVector0_set0
(JNIEnv *env, jclass clazz, jlong ptr, jlong idx, jboolean bit) {
  exec(BM_bvector_set_bit((BM_BVHANDLE)ptr, (unsigned int)idx, bit & 1));
}

/*
* Class:     io_bitmagic_core_BVector0
* Method:    inc0
* Signature: (JJ)I
*/
JNIEXPORT jint JNICALL Java_io_bitmagic_core_BVector0_inc0
(JNIEnv *env, jclass clazz, jlong ptr, jlong idx) {
  int co;
  exec(BM_bvector_inc_bit((BM_BVHANDLE)ptr, (unsigned int)idx, &co));
  return (jint)co;
}

/*
 * Class:     io_bitmagic_core_BVector0
 * Method:    setConditional0
 * Signature: (JJZZ)Z
 */
JNIEXPORT jboolean JNICALL Java_io_bitmagic_core_BVector0_setConditional0
(JNIEnv *env, jclass clazz, jlong ptr, jlong idx, jboolean bit, jboolean condition) {
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
(JNIEnv *env, jclass clazz, jlong ptr, jlong idx) {
  exec(BM_bvector_flip_bit((BM_BVHANDLE)ptr, (unsigned int)idx))
}

/*
 * Class:     io_bitmagic_core_BVector0
 * Method:    setAll0
 * Signature: (J)V
 */
JNIEXPORT void JNICALL Java_io_bitmagic_core_BVector0_setAll0
(JNIEnv *env, jclass clazz, jlong ptr) {
  exec(BM_bvector_set((BM_BVHANDLE)ptr));
}

/*
 * Class:     io_bitmagic_core_BVector0
 * Method:    setRange0
 * Signature: (JJJZ)V
 */
JNIEXPORT void JNICALL Java_io_bitmagic_core_BVector0_setRange0
(JNIEnv *env, jclass clazz, jlong ptr, jlong left, jlong right, jboolean bit) {
  exec(BM_bvector_set_range((BM_BVHANDLE)ptr, (unsigned int)left, (unsigned int)right, bit ? 1 : 0));
}

/*
 * Class:     io_bitmagic_core_BVector0
 * Method:    invert0
 * Signature: (J)V
 */
JNIEXPORT void JNICALL Java_io_bitmagic_core_BVector0_invert0
(JNIEnv *env, jclass clazz, jlong ptr) {
  exec(BM_bvector_invert((BM_BVHANDLE)ptr));
}

/*
 * Class:     io_bitmagic_core_BVector0
 * Method:    clear0
 * Signature: (JI)V
 */
JNIEXPORT void JNICALL Java_io_bitmagic_core_BVector0_clear0
(JNIEnv *env, jclass clazz, jlong ptr, jint free_mem) {
  exec(BM_bvector_clear((BM_BVHANDLE)ptr, free_mem));
}

/*
 * Class:     io_bitmagic_core_BVector0
 * Method:    extract0
 * Signature: (JJ)J
 */
JNIEXPORT jlong JNICALL Java_io_bitmagic_core_BVector0_extract0
(JNIEnv *env, jclass clazz, jlong ptr, jlong start) {
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
(JNIEnv *env, jclass clazz, jlong ptr, jlong idx) {
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
(JNIEnv *env, jclass clazz, jlong ptr) {
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
(JNIEnv *env, jclass clazz, jlong ptr, jlong left, jlong right) {
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
(JNIEnv *env, jclass clazz, jlong ptr) {
  int val;
  exec(BM_bvector_any((BM_BVHANDLE)ptr, &val));
  return (jboolean)val;
}

/*
* Class:     io_bitmagic_core_BVector0
* Method:    findFirst0
* Signature: (JJ)J
*/
JNIEXPORT jlong JNICALL Java_io_bitmagic_core_BVector0_findFirst0
(JNIEnv *env, jclass clazz, jlong ptr, jlong start) {
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
* Method:    findLast0
* Signature: (J)J
*/
JNIEXPORT jlong JNICALL Java_io_bitmagic_core_BVector0_findReverse0
(JNIEnv *env, jclass clazz, jlong ptr) {
  int found;
  unsigned int pos;
  exec(BM_bvector_find_reverse((BM_BVHANDLE)ptr, &pos, &found));
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
(JNIEnv *env, jclass clazz, jlong ptr1, jlong ptr2) {
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
(JNIEnv *env, jclass clazz, jlong ptr, jint opt_mode) {
  struct BM_bvector_statistics stat;
  jclass clazzect;
  jmethodID constructor;
  jobject cls;

  exec(BM_bvector_optimize((BM_BVHANDLE)ptr, opt_mode, &stat));
  cls = (*env)->FindClass(env, "io/bitmagic/core/BitVectorStat");
  constructor = (*env)->GetMethodID(env, cls, "<init>", "(JJJJ)V");
  jobject object = (*env)->NewObject(env, cls, constructor, (jlong)stat.bit_blocks, (jlong)stat.gap_blocks, (jlong)stat.max_serialize_mem, (jlong)stat.memory_used);
  return object;
}

/*
 * Class:     io_bitmagic_core_BVector0
 * Method:    calcStat0
 * Signature: (J)Lio/bitmagic/BitVectorStat;
 */
JNIEXPORT jobject JNICALL Java_io_bitmagic_core_BVector0_calcStat0
(JNIEnv *env, jclass clazz, jlong ptr) {
  struct BM_bvector_statistics stat;
  jclass clazzect;
  jmethodID constructor;
  jobject cls;
  
  exec(BM_bvector_calc_stat((BM_BVHANDLE)ptr, &stat));
  cls = (*env)->FindClass(env, "io/bitmagic/core/BitVectorStat");
  constructor = (*env)->GetMethodID(env, cls, "<init>", "(JJJJ)V");
  jobject object = (*env)->NewObject(env, cls, constructor, (jlong)stat.bit_blocks, (jlong)stat.gap_blocks, (jlong)stat.max_serialize_mem, (jlong)stat.memory_used);
  return object;
}

/***************************************** BitVector operations ****************************************/

/*
 * Class:     io_bitmagic_core_BVector0
 * Method:    operation0
 * Signature: (JJI)V
 */
JNIEXPORT void JNICALL Java_io_bitmagic_core_BVector0_operation0
(JNIEnv *env, jclass clazz, jlong dst, jlong src, jint op) {
  exec(BM_bvector_combine_operation((BM_BVHANDLE)dst, (BM_BVHANDLE)src, op));
}

/*
 * Class:     io_bitmagic_core_BVector0
 * Method:    and0
 * Signature: (JJ)V
 */
JNIEXPORT void JNICALL Java_io_bitmagic_core_BVector0_and0
(JNIEnv *env, jclass clazz, jlong dst, jlong src) {
  exec(BM_bvector_combine_AND((BM_BVHANDLE)dst, (BM_BVHANDLE)src));
}

/*
 * Class:     io_bitmagic_core_BVector0
 * Method:    or0
 * Signature: (JJ)V
 */
JNIEXPORT void JNICALL Java_io_bitmagic_core_BVector0_or0
(JNIEnv *env, jclass clazz, jlong dst, jlong src) {
  exec(BM_bvector_combine_OR((BM_BVHANDLE)dst, (BM_BVHANDLE)src));
}

/*
 * Class:     io_bitmagic_core_BVector0
 * Method:    sub0
 * Signature: (JJ)V
 */
JNIEXPORT void JNICALL Java_io_bitmagic_core_BVector0_sub0
(JNIEnv *env, jclass clazz, jlong dst, jlong src) {
  exec(BM_bvector_combine_SUB((BM_BVHANDLE)dst, (BM_BVHANDLE)src));
}

/*
 * Class:     io_bitmagic_core_BVector0
 * Method:    xor0
 * Signature: (JJ)V
 */
JNIEXPORT void JNICALL Java_io_bitmagic_core_BVector0_xor0
(JNIEnv *env, jclass clazz, jlong dst, jlong src) {
  exec(BM_bvector_combine_XOR((BM_BVHANDLE)dst, (BM_BVHANDLE)src));
}

/***************************************** BitVector operations with arrays ****************************************/

/*
* Class:     io_bitmagic_core_BVector0
* Method:    andArr0
* Signature: (J[I)V
*/
JNIEXPORT void JNICALL Java_io_bitmagic_core_BVector0_andArr0
(JNIEnv *env, jclass class, jlong dst, jintArray arr) {
  unsigned int *start = (unsigned int*)(*env)->GetPrimitiveArrayCritical(env, arr, 0);
  null_out_of_memory(start);
  jsize size = (*env)->GetArrayLength(env, arr);
  exec(BM_bvector_combine_AND_arr((BM_BVHANDLE)dst, start, start + size));
  (*env)->ReleasePrimitiveArrayCritical(env, arr, start, JNI_ABORT); // The array is read-only
}

/*
* Class:     io_bitmagic_core_BVector0
* Method:    andArrSorted0
* Signature: (J[I)V
*/
JNIEXPORT void JNICALL Java_io_bitmagic_core_BVector0_andArrSorted0
(JNIEnv *env, jclass class, jlong dst, jintArray arr) {
  unsigned int *start = (unsigned int*)(*env)->GetPrimitiveArrayCritical(env, arr, 0);
  null_out_of_memory(start);
  jsize size = (*env)->GetArrayLength(env, arr);
  exec(BM_bvector_combine_AND_arr_sorted((BM_BVHANDLE)dst, start, start + size));
  (*env)->ReleasePrimitiveArrayCritical(env, arr, start, JNI_ABORT); // The array is read-only
}

/*
* Class:     io_bitmagic_core_BVector0
* Method:    orArr0
* Signature: (J[I)V
*/
JNIEXPORT void JNICALL Java_io_bitmagic_core_BVector0_orArr0
(JNIEnv *env, jclass class, jlong dst, jintArray arr) {
  unsigned int *start = (unsigned int*)(*env)->GetPrimitiveArrayCritical(env, arr, 0);
  null_out_of_memory(start);
  jsize size = (*env)->GetArrayLength(env, arr);
  exec(BM_bvector_combine_OR_arr((BM_BVHANDLE)dst, start, start + size));
  (*env)->ReleasePrimitiveArrayCritical(env, arr, start, JNI_ABORT); // The array is read-only
}

/*
* Class:     io_bitmagic_core_BVector0
* Method:    xorArr0
* Signature: (J[I)V
*/
JNIEXPORT void JNICALL Java_io_bitmagic_core_BVector0_xorArr0
(JNIEnv *env, jclass class, jlong dst, jintArray arr) {
  unsigned int *start = (unsigned int*)(*env)->GetPrimitiveArrayCritical(env, arr, 0);
  null_out_of_memory(start);
  jsize size = (*env)->GetArrayLength(env, arr);
  exec(BM_bvector_combine_XOR_arr((BM_BVHANDLE)dst, start, start + size));
  (*env)->ReleasePrimitiveArrayCritical(env, arr, start, JNI_ABORT); // The array is read-only
}

/*
* Class:     io_bitmagic_core_BVector0
* Method:    subArr0
* Signature: (J[I)V
*/
JNIEXPORT void JNICALL Java_io_bitmagic_core_BVector0_subArr0
(JNIEnv *env, jclass class, jlong dst, jintArray arr) {
  unsigned int *start = (unsigned int*)(*env)->GetPrimitiveArrayCritical(env, arr, 0);
  null_out_of_memory(start);
  jsize size = (*env)->GetArrayLength(env, arr);
  exec(BM_bvector_combine_SUB_arr((BM_BVHANDLE)dst, start, start + size));
  (*env)->ReleasePrimitiveArrayCritical(env, arr, start, JNI_ABORT); // The array is read-only
}

/***************************************** BitVector serialization ****************************************/
/*
* Class:     io_bitmagic_core_BVector0
* Method:    deserialize0
* Signature: (J[B)V
*/
JNIEXPORT void JNICALL Java_io_bitmagic_core_BVector0_deserialize0
(JNIEnv *env, jclass clazz, jlong ptr, jbyteArray ba) {
  void *start = (*env)->GetPrimitiveArrayCritical(env, ba, 0);
  if (start == NULL) {
    jclass ex = (*env)->FindClass(env, "java/lang/OutOfMemoryError");
    (*env)->ThrowNew(env, ex, "Out of memory error in native call deserailize0()");
    return;
  }
  jsize size = (*env)->GetArrayLength(env, ba);
  exec(BM_bvector_deserialize((BM_BVHANDLE)ptr, (const char*)start, size));
  (*env)->ReleasePrimitiveArrayCritical(env, ba, start, JNI_ABORT); // The array is read-only
}

/*
* Class:     io_bitmagic_core_BVector0
* Method:    serialize0
* Signature: (J[B)J
*/
JNIEXPORT jlong JNICALL Java_io_bitmagic_core_BVector0_serialize0
(JNIEnv *env, jclass clazz, jlong ptr, jbyteArray ba) {
  void *start = (*env)->GetPrimitiveArrayCritical(env, ba, 0);
  if (start == NULL) {
    jclass ex = (*env)->FindClass(env, "java/lang/OutOfMemoryError"); 
    (*env)->ThrowNew(env, ex, "Out of memory error in native call serailize0()");
    return 0;
  }
  jsize size = (*env)->GetArrayLength(env, ba);
  size_t actualSize;

  exec(BM_bvector_serialize((BM_BVHANDLE)ptr, start, size, &actualSize));
  (*env)->ReleasePrimitiveArrayCritical(env, ba, start, 0); 
  return actualSize;
}

/***************************************** BitVector Iterator ****************************************/
/*
* Class:     io_bitmagic_core_BVIterator0
* Method:    create0
* Signature: (J)J
*/
JNIEXPORT jlong JNICALL Java_io_bitmagic_core_BVIterator0_create0
(JNIEnv *env, jclass clazz, jlong bvPtr) {
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
(JNIEnv *env, jclass clazz, jlong ptr) {
  exec(BM_bvector_enumerator_free((BM_BVEHANDLE)ptr));
}

/*
* Class:     io_bitmagic_core_BVIterator0
* Method:    isValid0
* Signature: (J)Z
*/
JNIEXPORT jboolean JNICALL Java_io_bitmagic_core_BVIterator0_isValid0
(JNIEnv *env, jclass clazz, jlong ptr) {
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
(JNIEnv *env, jclass clazz, jlong ptr) {
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
(JNIEnv *env, jclass clazz, jlong ptr) {
  unsigned int value;
  int valid;
  exec(BM_bvector_enumerator_next((BM_BVEHANDLE)ptr, &valid, &value));
  return (jboolean)(valid != 0);
}
