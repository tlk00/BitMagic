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

/*! \file bmundef.h
    \brief pre-processor un-defines to avoid global space pollution (internal)
*/

#undef BMRESTRICT
#undef BMFORCEINLINE
#undef BMGAP_PTR
#undef BMSET_PTRGAP
#undef BM_IS_GAP
#undef BMPTR_SETBIT0
#undef BMPTR_CLEARBIT0
#undef BMPTR_TESTBIT0
#undef BM_SET_MMX_GUARD
#undef SER_NEXT_GRP
#undef BM_SET_ONE_BLOCKS
#undef DECLARE_TEMP_BLOCK
#undef BM_MM_EMPTY
#undef BM_ASSERT
#undef FULL_BLOCK_ADDR
#undef IS_VALID_ADDR
#undef IS_FULL_BLOCK
#undef IS_EMPTY_BLOCK
#undef BM_INCWORD_BITCOUNT
#undef BM_MINISET_GAPLEN
#undef BM_MINISET_ARRSIZE

#undef BMVECTOPT
#undef VECT_XOR_ARR_2_MASK
#undef VECT_ANDNOT_ARR_2_MASK
#undef VECT_BITCOUNT
#undef VECT_BITCOUNT_AND
#undef VECT_BITCOUNT_OR
#undef VECT_BITCOUNT_XOR
#undef VECT_BITCOUNT_SUB
#undef VECT_INVERT_ARR
#undef VECT_AND_ARR
#undef VECT_OR_ARR
#undef VECT_OR_ARR_3WAY
#undef VECT_OR_ARR_5WAY
#undef VECT_SUB_ARR
#undef VECT_XOR_ARR

#undef VECT_COPY_BLOCK
#undef VECT_SET_BLOCK
#undef VECT_IS_ZERO_BLOCK
#undef VECT_IS_ONE_BLOCK



#undef BM_UNALIGNED_ACCESS_OK
#undef BM_x86


