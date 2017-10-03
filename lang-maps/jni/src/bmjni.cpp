#include "org_tlk00_bitmagic_BVector0.h"
#include "libbm.h"

extern "C" {

JNIEXPORT void JNICALL Java_org_tlk00_bitmagic_BVector0_init0
  (JNIEnv *, jobject, jlong ptr) {
  BM_init((void*)ptr);
}

/*
 * Class:     org_tlk00_bitmagic_BVector0
 * Method:    create0
 * Signature: (IJ)J
 */
JNIEXPORT jlong JNICALL Java_org_tlk00_bitmagic_BVector0_create0
  (JNIEnv *, jobject, jint, jlong size) {
  BM_BVHANDLE ptr;
  BM_bvector_construct(&ptr, size);
  return (jlong)ptr;
}

/*
 * Class:     org_tlk00_bitmagic_BVector0
 * Method:    clone0
 * Signature: (J)J
 */
JNIEXPORT jlong JNICALL Java_org_tlk00_bitmagic_BVector0_clone0
  (JNIEnv *, jobject, jlong ptr) {
  return ptr;
}

/*
 * Class:     org_tlk00_bitmagic_BVector0
 * Method:    dispose0
 * Signature: (J)V
 */
JNIEXPORT void JNICALL Java_org_tlk00_bitmagic_BVector0_dispose0
  (JNIEnv *, jobject, jlong ptr) {
  BM_bvector_free((BM_BVHANDLE)ptr);
}

/*
 * Class:     org_tlk00_bitmagic_BVector0
 * Method:    set0
 * Signature: (JJ)V
 */
JNIEXPORT void JNICALL Java_org_tlk00_bitmagic_BVector0_set0
  (JNIEnv *, jobject, jlong ptr, jlong bit) {
  //BM_bvector_set_bit((BM_BVHANDLE)ptr, (unsigned int)bit);
}

/*
 * Class:     org_tlk00_bitmagic_BVector0
 * Method:    get0
 * Signature: (JJ)Z
 */
JNIEXPORT jboolean JNICALL Java_org_tlk00_bitmagic_BVector0_get0
  (JNIEnv *, jobject, jlong ptr, jlong bit) {
  unsigned int ret;
  //BM_bvector_get_bit((BM_BVHANDLE)ptr, (unsigned int)bit, &ret);
  return false; //(jboolean)ret;
}
  
} // extern