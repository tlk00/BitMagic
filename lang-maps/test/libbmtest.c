/*
     BitMagic Library C - unit test.
*/


/*
Copyright(c) 2002-2018 Anatoliy Kuznetsov(anatoliy_kuznetsov at yahoo.com)

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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "libbm.h"


static
int check_report_error(int res, const char* msg)
{
    if (res != BM_OK)
    {
        printf("Error: %s. err_code=%i \'%s\'\n", msg, res, BM_error_msg(res));
    }
    return res;
}

#define BMERR_CHECK(x, y) if (check_report_error(x, y) != 0) return x;
#define BMERR_CHECK_GOTO(x, y, z) if (check_report_error(x, y) != 0) goto z;

static
int InitTest()
{
    int res = 0;
    const char* c;
    int major, minor, patch;
    const char* msg;
    int simd_version;
    
    
    res = BM_init(0);
    BMERR_CHECK(res, "BM_init()");
    
    c = BM_version(&major, &minor, &patch);
    BMERR_CHECK(res, "BM_version()");
    
    msg = BM_error_msg(BM_ERR_BADARG);
    if (!msg)
    {
        printf("BM_error_msg() failed \n");
        return 1;
    }
    else
    {
        printf("Test message, ignore. %s\n", msg);
    }
    
    printf ("%s\n", c);
    
    simd_version = BM_simd_version();
    switch(simd_version)
    {
    case BM_SIMD_NO:
        printf("BitMagic vanilla.\n");
        break;
    case BM_SIMD_SSE2:
        printf("BitMagic for SSE2 \n");
        break;
    case BM_SIMD_SSE42:
        printf("BitMagic for SSE4.2 \n");
        break;
    case BM_SIMD_AVX2:
        printf("BitMagic for AVX2 \n");
        break;
    default:
        printf("Unknown SIMD code \n");
        break;
    }

    return 0;
}

static
int ConstructionTest()
{
    int res = 0;
    BM_BVHANDLE bmh = 0;
    
    res = BM_bvector_construct(&bmh, 0);
    BMERR_CHECK(res, "BM_bvector_construct()");

    res = BM_bvector_free(bmh);
    BMERR_CHECK(res, "BM_bvector_free()");

    return res;
}

static
int FreezeTest()
{
    int res = 0;
    int val;
    BM_BVHANDLE bmh = 0;
    BM_BVHANDLE bmh2 = 0;
    BM_BVHANDLE bmh3 = 0;

    res = BM_bvector_construct(&bmh, 0);
    BMERR_CHECK(res, "BM_bvector_construct()");

    res = BM_bvector_is_ro(bmh, &val);
    BMERR_CHECK_GOTO(res, "BM_bvector_is_ro()", free_mem);
    if (val)
    {
        printf("BM_bvector_is_ro incorrect value \n");
        res = 1; goto free_mem;
    }

    res = BM_bvector_set_bit(bmh, 10, BM_TRUE); // vector has to be not empty to freeze
    BMERR_CHECK_GOTO(res, "BM_bvector_set_bit()", free_mem);

    res = BM_bvector_construct_copy_ro(&bmh2, bmh);
    BMERR_CHECK(res, "BM_bvector_construct()");

    res = BM_bvector_is_ro(bmh2, &val);
    BMERR_CHECK_GOTO(res, "BM_bvector_is_ro()", free_mem);
    if (!val)
    {
        printf("BM_bvector_is_ro (bmh2) incorrect value \n");
        res = 1; goto free_mem;
    }


    res = BM_bvector_freeze(bmh);
    BMERR_CHECK_GOTO(res, "BM_bvector_freeze()", free_mem);

    res = BM_bvector_is_ro(bmh, &val);
    BMERR_CHECK_GOTO(res, "BM_bvector_is_ro()", free_mem);
    if (!val)
    {
        printf("BM_bvector_is_ro incorrect value \n");
        res = 1; goto free_mem;
    }

    res = BM_bvector_get_bit(bmh, 10, &val);
    BMERR_CHECK_GOTO(res, "BM_bvector_get_bit()", free_mem);
    if (!val)
    {
        printf("bvector get_bit incorrect value \n");
        res = 1; goto free_mem;
    }

    res = BM_bvector_get_bit(bmh2, 10, &val);
    BMERR_CHECK_GOTO(res, "BM_bvector_get_bit()", free_mem);
    if (!val)
    {
        printf("bvector get_bit incorrect value \n");
        res = 1; goto free_mem;
    }



    res = BM_bvector_construct_copy_rw(&bmh3, bmh2);
    BMERR_CHECK(res, "BM_bvector_construct()");
    if (!bmh3)
    {
        printf("BM_bvector_construct_copy_rw incorrect value \n");
        res = 1; goto free_mem;
    }

    res = BM_bvector_is_ro(bmh3, &val);
    BMERR_CHECK_GOTO(res, "BM_bvector_is_ro()", free_mem);
    printf("%i\n", val);
    if (val)
    {
        printf("BM_bvector_is_ro incorrect value \n");
        res = 1; goto free_mem;
    }
    res = BM_bvector_get_bit(bmh3, 10, &val);
    BMERR_CHECK_GOTO(res, "BM_bvector_get_bit()", free_mem);
    if (!val)
    {
        printf("bvector get_bit incorrect value \n");
        res = 1; goto free_mem;
    }

free_mem:
    res = BM_bvector_free(bmh);
    BMERR_CHECK(res, "BM_bvector_free()");
    res = BM_bvector_free(bmh2);
    BMERR_CHECK(res, "BM_bvector_free()");
    res = BM_bvector_free(bmh3);
    BMERR_CHECK(res, "BM_bvector_free()");

    return res;
}



static
int ResizeTest()
{
    int res = 0;
    BM_BVHANDLE bmh = 0;
    unsigned int size1 = 100000;
    unsigned int size2 = 100000;
    unsigned int size;
    
    
    res = BM_bvector_construct(&bmh, size1);
    BMERR_CHECK(res, "BM_bvector_construct()");
    
    
    res = BM_bvector_get_size(bmh, &size);
    BMERR_CHECK_GOTO(res, "BM_bvector_get_size()", free_mem);
    if (size != size1)
    {
        printf("get size test failed %i\n", size);
        res = 1; goto free_mem;
    }
    /*
    {
    unsigned int capacity;
    res = BM_bvector_get_capacity(bmh, &capacity);
    BMERR_CHECK_GOTO(res, "BM_bvector_get_capacity()", free_mem);
    printf("cap=%i", capacity);
    }
    */
    
    res = BM_bvector_set_size(bmh, size2);
    BMERR_CHECK_GOTO(res, "BM_bvector_get_size", free_mem);
    
    res = BM_bvector_get_size(bmh, &size);
    BMERR_CHECK_GOTO(res, "BM_bvector_get_size()", free_mem);
    if (size != size2)
    {
        printf("bvector get size failed %i\n", size);
        res = 1; goto free_mem;
    }
    
    
    free_mem:
        res = BM_bvector_free(bmh);
        BMERR_CHECK(res, "BM_bvector_free()");
    return res;
}

static
int ConstructionCopyMoveTest()
{
    int res = 0;
    BM_BVHANDLE bmh1 = 0;
    BM_BVHANDLE bmh2 = 0;
    BM_BVHANDLE bmh3 = 0;
    unsigned int count1, count2, count3;
    
    res = BM_bvector_construct(&bmh1, 0);
    BMERR_CHECK(res, "BM_bvector_construct()");
    res = BM_bvector_construct(&bmh2, 0);
    BMERR_CHECK_GOTO(res, "BM_bvector_construct()", free_mem);
    res = BM_bvector_construct(&bmh3, 0);
    BMERR_CHECK_GOTO(res, "BM_bvector_construct()", free_mem);
    
    
    
    res = BM_bvector_set_bit(bmh1, 10, BM_TRUE);
    BMERR_CHECK_GOTO(res, "BM_bvector_set_bit()", free_mem);

    res = BM_bvector_set_bit(bmh2, 100, BM_TRUE);
    BMERR_CHECK_GOTO(res, "BM_bvector_set_bit()", free_mem);
    res = BM_bvector_set_bit(bmh2, 101, BM_TRUE);
    BMERR_CHECK_GOTO(res, "BM_bvector_set_bit()", free_mem);
    
    
    res = BM_bvector_swap(bmh1, bmh2);
    BMERR_CHECK_GOTO(res, "BM_bvector_swap()", free_mem);
    
    res = BM_bvector_count(bmh1, &count1);
    BMERR_CHECK_GOTO(res, "BM_bvector_count()", free_mem);
    if (count1 != 2)
    {
        printf("1. incorrrect count %i \n", count1);
        res = 1; goto free_mem;
    }
    res = BM_bvector_count(bmh2, &count2);
    BMERR_CHECK_GOTO(res, "BM_bvector_count()", free_mem);
    if (count2 != 1)
    {
        printf("2. incorrrect count %i \n", count2);
        res = 1; goto free_mem;
    }
    
    res = BM_bvector_swap(bmh1, bmh3);
    BMERR_CHECK_GOTO(res, "BM_bvector_swap()", free_mem);

    res = BM_bvector_count(bmh3, &count3);
    BMERR_CHECK_GOTO(res, "BM_bvector_count()", free_mem);
    if (count3 != 2)
    {
        printf("3. incorrrect count %i \n", count3);
        res = 1; goto free_mem;
    }
    res = BM_bvector_count(bmh1, &count1);
    BMERR_CHECK_GOTO(res, "BM_bvector_count()", free_mem);
    if (count1 != 0)
    {
        printf("4. incorrrect count %i \n", count1);
        res = 1; goto free_mem;
    }
    
    {
        BM_BVHANDLE bmh4 = 0;
        unsigned count4;
        unsigned pos;
        int cmp, found;
        
        res = BM_bvector_construct_copy(&bmh4, bmh3);
        BMERR_CHECK_GOTO(res, "BM_bvector_construct_copy()", free_mem);
        res = BM_bvector_count(bmh4, &count4);
        BMERR_CHECK_GOTO(res, "BM_bvector_count()", free_mem1);
        if (count4 != 2)
        {
            printf("5. incorrrect count %i \n", count4);
            res = BM_bvector_free(bmh4);
            res = 1; goto free_mem1;
        }
        
        res = BM_bvector_compare(bmh4, bmh3, &cmp);
        BMERR_CHECK_GOTO(res, "BM_bvector_compare()", free_mem1);
        if (cmp != 0)
        {
            printf("5. incorrrect compare result %i \n", cmp);
            res = 1; goto free_mem1;
        }
        res = BM_bvector_find_first_mismatch(bmh4, bmh3, &pos, &found);
        BMERR_CHECK_GOTO(res, "BM_bvector_find_first_mismatch()", free_mem1);
        if (found)
        {
            printf("6. incorrrect find mismatch result %i \n", found);
            res = 1; goto free_mem1;
        }

        
    free_mem1:
        res = BM_bvector_free(bmh4);
    }
    

    free_mem:
        res = BM_bvector_free(bmh1);
        res = BM_bvector_free(bmh2);
        res = BM_bvector_free(bmh3);

    return res;
}


static
int SetGetTest()
{
    int res = 0;
    BM_BVHANDLE bmh1 = 0;
    BM_BVHANDLE bmh2 = 0;
    int val;
    unsigned int count;
    int carry_over;
    
    res = BM_bvector_construct(&bmh1, 200);
    BMERR_CHECK(res, "BM_bvector_construct()");

    res = BM_bvector_construct(&bmh2, 0);
    BMERR_CHECK_GOTO(res, "BM_bvector_construct()", free_mem);

    res = BM_bvector_init(bmh2);
    BMERR_CHECK_GOTO(res, "BM_bvector_init()", free_mem);

    res = BM_bvector_any(bmh1, &val);
    BMERR_CHECK_GOTO(res, "BM_bvector_any()", free_mem);
    if (val)
    {
        printf("bvector any() incorrect value \n");
        res = 1; goto free_mem;
    }
    
    res = BM_bvector_set_bit(bmh1, 10, BM_TRUE);
    BMERR_CHECK_GOTO(res, "BM_bvector_set_bit()", free_mem);

    res = BM_bvector_set_bit_no_check(bmh2, 10);
    BMERR_CHECK_GOTO(res, "BM_bvector_set_bit_no_check()", free_mem);


    res = BM_bvector_any(bmh1, &val);
    BMERR_CHECK_GOTO(res, "BM_bvector_any()", free_mem);
    if (val==0)
    {
        printf("bvector any() incorrect value \n");
        res = 1; goto free_mem;
    }
    
    
    res = BM_bvector_get_bit(bmh1, 10, &val);
    BMERR_CHECK_GOTO(res, "BM_bvector_get_bit()", free_mem);
    if (!val)
    {
        printf("bvector get_bit incorrect value \n");
        res = 1; goto free_mem;
    }

    res = BM_bvector_get_bit(bmh2, 10, &val);
    BMERR_CHECK_GOTO(res, "BM_bvector_get_bit()", free_mem);
    if (!val)
    {
        printf("bvector get_bit incorrect value \n");
        res = 1; goto free_mem;
    }

    res = BM_bvector_inc_bit(bmh2, 10, &carry_over);
    BMERR_CHECK_GOTO(res, "BM_bvector_inc_bit()", free_mem);
    if (!carry_over)
    {
        printf("bvector inc_bit incorrect value \n");
        res = 1; goto free_mem;
    }

    res = BM_bvector_inc_bit(bmh2, 10, &carry_over);
    BMERR_CHECK_GOTO(res, "BM_bvector_inc_bit()", free_mem);
    if (carry_over)
    {
        printf("bvector inc_bit incorrect value \n");
        res = 1; goto free_mem;
    }


    res = BM_bvector_get_bit(bmh1, 0, &val);
    BMERR_CHECK_GOTO(res, "BM_bvector_get_bit()", free_mem);
    if (val)
    {
        printf("bvector get_bit incorrect value \n");
        res = 1; goto free_mem;
    }
    
    res = BM_bvector_flip_bit(bmh1, 0);
    BMERR_CHECK_GOTO(res, "BM_bvector_flip_bit()", free_mem);

    res = BM_bvector_get_bit(bmh1, 0, &val);
    BMERR_CHECK_GOTO(res, "BM_bvector_get_bit()", free_mem);
    if (val == 0)
    {
        printf("bvector get_bit incorrect value \n");
        res = 1; goto free_mem;
    }
    res = BM_bvector_flip_bit(bmh1, 0);
    BMERR_CHECK_GOTO(res, "BM_bvector_flip_bit()", free_mem);
    
    
    res = BM_bvector_count(bmh1, &count);
    BMERR_CHECK_GOTO(res, "BM_bvector_count()", free_mem);
    if (count != 1)
    {
        printf("incorrrect count %i \n", count);
        res = 1; goto free_mem;
    }
    
    res = BM_bvector_set_bit(bmh1, 10, BM_FALSE);
    BMERR_CHECK_GOTO(res, "BM_bvector_set_bit()", free_mem);
    
    res = BM_bvector_get_bit(bmh1, 0, &val);
    BMERR_CHECK_GOTO(res, "BM_bvector_get_bit()", free_mem);
    if (val != BM_FALSE)
    {
        printf("bvector get_bit incorrect value %i\n", val);
        res = 1; goto free_mem;
    }
    res = BM_bvector_count(bmh1, &count);
    BMERR_CHECK_GOTO(res, "BM_bvector_count()", free_mem);
    if (count != 0)
    {
        printf("incorrrect count %i \n", count);
        res = 1; goto free_mem;
    }
    
    {
        int change;
        res = BM_bvector_set_bit_conditional(bmh1, 0, BM_TRUE, BM_FALSE, &change);
        BMERR_CHECK_GOTO(res, "BM_bvector_set_bit_conditional()", free_mem);
        if (!change)
        {
            printf("bvector set_bit_conditional error \n");
            res = 1; goto free_mem;
        }
        res = BM_bvector_count(bmh1, &count);
        BMERR_CHECK_GOTO(res, "BM_bvector_count()", free_mem);
        if (count != 1)
        {
            printf("incorrrect count %i \n", count);
            res = 1; goto free_mem;
        }
        res = BM_bvector_set_bit_conditional(bmh1, 0, BM_TRUE, BM_FALSE, &change);
        BMERR_CHECK_GOTO(res, "BM_bvector_set_bit_conditional()", free_mem);
        if (change)
        {
            printf("bvector set_bit_conditional error \n");
            res = 1; goto free_mem;
        }
    }
    
    res = BM_bvector_set(bmh1);
    BMERR_CHECK_GOTO(res, "BM_bvector_set()", free_mem);

    res = BM_bvector_count(bmh1, &count);
    BMERR_CHECK_GOTO(res, "BM_bvector_count()", free_mem);
    if (count == 0)
    {
        printf("incorrrect count %i \n", count);
        res = 1; goto free_mem;
    }
    res = BM_bvector_clear(bmh1, BM_TRUE);
    BMERR_CHECK_GOTO(res, "BM_bvector_clear()", free_mem);

    res = BM_bvector_count(bmh1, &count);
    BMERR_CHECK_GOTO(res, "BM_bvector_count()", free_mem);
    if (count != 0)
    {
        printf("incorrrect count %i \n", count);
        res = 1; goto free_mem;
    }
    
    
    
    free_mem:
        res = BM_bvector_free(bmh1);
        BMERR_CHECK(res, "bvector free failed");
        res = BM_bvector_free(bmh2);
        BMERR_CHECK(res, "bvector free failed");

    return res;
}

int RangeTest()
{
    int res = 0;
    BM_BVHANDLE bmh = 0;
    unsigned count;
    
    res = BM_bvector_construct(&bmh, 200);
    BMERR_CHECK(res, "BM_bvector_construct()");
    
    res = BM_bvector_set_range(bmh, 10, 20, BM_TRUE);
    BMERR_CHECK_GOTO(res, "BM_bvector_set_range()", free_mem);

    res = BM_bvector_count(bmh, &count);
    BMERR_CHECK_GOTO(res, "BM_bvector_count()", free_mem);
    if (count != 11)
    {
        printf("incorrrect count %i \n", count);
        res = 1; goto free_mem;
    }
    
    res = BM_bvector_count_range(bmh, 10, 20, &count);
    BMERR_CHECK_GOTO(res, "BM_bvector_count_range()", free_mem);
    if (count != 11)
    {
        printf("incorrrect count_range %i \n", count);
        res = 1; goto free_mem;
    }
    
    {
    unsigned int idx;
    int found;
    res = BM_bvector_find_rank(bmh, 10, 11, &idx, &found);
    BMERR_CHECK_GOTO(res, "BM_bvector_find_rank()", free_mem);
    if (idx != 20)
    {
        printf("incorrrect find_rank %i \n", idx);
        res = 1; goto free_mem;
    }
    }
    
    res = BM_bvector_count_range(bmh, 0, 9, &count);
    BMERR_CHECK_GOTO(res, "BM_bvector_count_range()", free_mem);
    if (count != 0)
    {
        printf("incorrrect count_range %i \n", count);
        res = 1; goto free_mem;
    }

    res = BM_bvector_set_range(bmh, 10, 20, BM_FALSE);
    BMERR_CHECK_GOTO(res, "BM_bvector_set_range()", free_mem);

    res = BM_bvector_count(bmh, &count);
    BMERR_CHECK_GOTO(res, "BM_bvector_count()", free_mem);
    if (count != 0)
    {
        printf("incorrrect count %i \n", count);
        res = 1; goto free_mem;
    }

    res = BM_bvector_invert(bmh);
    BMERR_CHECK_GOTO(res, "BM_bvector_invert()", free_mem);
    
    res = BM_bvector_count(bmh, &count);
    BMERR_CHECK_GOTO(res, "BM_bvector_count()", free_mem);
    if (count == 0)
    {
        printf("incorrrect count %i \n", count);
        res = 1; goto free_mem;
    }
    res = BM_bvector_invert(bmh);
    BMERR_CHECK_GOTO(res, "BM_bvector_invert()", free_mem);
    
    res = BM_bvector_count(bmh, &count);
    BMERR_CHECK_GOTO(res, "BM_bvector_count()", free_mem);
    if (count != 0)
    {
        printf("incorrrect count %i \n", count);
        res = 1; goto free_mem;
    }
    
    

    free_mem:
        res = BM_bvector_free(bmh);
        BMERR_CHECK(res, "bvector free failed");

    return res;
}


int GetNextTest()
{
    int res = 0;
    BM_BVHANDLE bmh = 0;
    unsigned int idx;
    int found;

    res = BM_bvector_construct(&bmh, 0);
    BMERR_CHECK(res, "BM_bvector_construct()");
        
                
    res = BM_bvector_get_first(bmh, &idx, &found);
    BMERR_CHECK_GOTO(res, "BM_bvector_get_first()", free_mem);
    if (found || idx != 0)
    {
        printf("1. incorrrect get first found on an empty vector \n");
        res = 1; goto free_mem;
    }
    
    res = BM_bvector_find(bmh, 0, &idx, &found);
    BMERR_CHECK_GOTO(res, "BM_bvector_find()", free_mem);
    if (found)
    {
        printf("1.1 incorrrect find found on an empty vector \n");
        res = 1; goto free_mem;
    }
    res = BM_bvector_find_reverse(bmh, &idx, &found);
    BMERR_CHECK_GOTO(res, "BM_bvector_find_reverse()", free_mem);
    if (found)
    {
        printf("1.2 incorrrect find_reverse found on an empty vector \n");
        res = 1; goto free_mem;
    }

    res = BM_bvector_set_bit(bmh, 0, BM_TRUE);
    BMERR_CHECK_GOTO(res, "BM_bvector_set_bit()", free_mem);
    
    res = BM_bvector_get_first(bmh, &idx, &found);
    BMERR_CHECK_GOTO(res, "BM_bvector_get_first()", free_mem);
    if (!found || idx != 0)
    {
        printf("2. incorrrect get first %i\n", idx);
        res = 1; goto free_mem;
    }
    res = BM_bvector_find(bmh, 0, &idx, &found);
    BMERR_CHECK_GOTO(res, "BM_bvector_find()", free_mem);
    if (!found || idx != 0)
    {
        printf("2.1 incorrrect find found in 0 position \n");
        res = 1; goto free_mem;
    }
    res = BM_bvector_find_reverse(bmh, &idx, &found);
    BMERR_CHECK_GOTO(res, "BM_bvector_find_reverse()", free_mem);
    if (!found || idx != 0)
    {
        printf("1.2 incorrrect find_reverse found in 0 position \n");
        res = 1; goto free_mem;
    }

    res = BM_bvector_get_next(bmh, idx, &idx);
    BMERR_CHECK_GOTO(res, "BM_bvector_get_next()", free_mem);
    if (idx != 0)
    {
        printf("4. incorrrect get next \n");
        res = 1; goto free_mem;
    }

    
    res = BM_bvector_set_range(bmh, 100000, 100002, BM_TRUE);
    BMERR_CHECK_GOTO(res, "BM_bvector_set_range()", free_mem);

    res = BM_bvector_get_next(bmh, idx, &idx);
    BMERR_CHECK_GOTO(res, "BM_bvector_get_next()", free_mem);
    if (idx != 100000)
    {
        printf("4. incorrrect get next \n");
        res = 1; goto free_mem;
    }
    
    res = BM_bvector_extract_next(bmh, 0, &idx);
    BMERR_CHECK_GOTO(res, "BM_bvector_extract_next()", free_mem);
    if (idx != 100000)
    {
        printf("5. incorrrect extract next \n");
        res = 1; goto free_mem;
    }
    
    res = BM_bvector_get_next(bmh, idx, &idx);
    BMERR_CHECK_GOTO(res, "BM_bvector_get_next()", free_mem);
    if (idx != 100001)
    {
        printf("6. incorrrect get next \n");
        res = 1; goto free_mem;
    }
    
    

    free_mem:
        res = BM_bvector_free(bmh);
        BMERR_CHECK(res, "bvector free failed");

    return res;
}

static
int PrintVector(BM_BVHANDLE bmh, unsigned int size)
{
    unsigned int count;
    int res = 0;
    int val;
    
    res = BM_bvector_count(bmh, &count);
    BMERR_CHECK_GOTO(res, "BM_bvector_count()", ret);
    
    printf("%u : ", count);

    for (unsigned int i = 0; i < size; ++i)
    {
        res = BM_bvector_get_bit(bmh, i, &val);
        BMERR_CHECK_GOTO(res, "BM_bvector_get_bit()", ret);
        printf("%u", val);

    } // for
ret:
    printf("\n");
    return res;
}

static
int CompareVectors(BM_BVHANDLE bmh1, BM_BVHANDLE bmh2, int* is_equal)
{
    int res = 0;
    BM_BVEHANDLE bmeh1 = 0;
    BM_BVEHANDLE bmeh2 = 0;
    int valid1, valid2;
    unsigned pos1, pos2;

    res = BM_bvector_enumerator_construct(bmh1, &bmeh1);
    BMERR_CHECK_GOTO(res, "BM_bvector_enumerator_construct()", free_mem);
    res = BM_bvector_enumerator_construct(bmh2, &bmeh2);
    BMERR_CHECK_GOTO(res, "BM_bvector_enumerator_construct()", free_mem);

    res = BM_bvector_enumerator_is_valid(bmeh1, &valid1);
    BMERR_CHECK_GOTO(res, "BM_bvector_enumerator_is_valid()", free_mem);
    res = BM_bvector_enumerator_is_valid(bmeh2, &valid2);
    BMERR_CHECK_GOTO(res, "BM_bvector_enumerator_is_valid()", free_mem);
    
    if (!valid1)
    {
        if (!valid2)
            *is_equal = 1;
        else
            *is_equal = 0;
        return 0;
    }
    res = BM_bvector_enumerator_get_value(bmeh1, &pos1);
    BMERR_CHECK_GOTO(res, "BM_bvector_enumerator_get_value()", free_mem);
    res = BM_bvector_enumerator_get_value(bmeh2, &pos2);
    BMERR_CHECK_GOTO(res, "BM_bvector_enumerator_get_value()", free_mem);
    
    while(valid1 && valid2)
    {
        if (pos1 != pos2)
        {
            *is_equal = 0;
            goto free_mem;
        }
        res = BM_bvector_enumerator_next(bmeh1, &valid1, &pos1);
        BMERR_CHECK_GOTO(res, "BM_bvector_enumerator_next()", free_mem);
        res = BM_bvector_enumerator_next(bmeh2, &valid2, &pos2);
        BMERR_CHECK_GOTO(res, "BM_bvector_enumerator_next()", free_mem);
    }
    
    if (valid1 == valid2)
    {
        *is_equal = 1;
    }

    
free_mem:
    res = BM_bvector_enumerator_free(bmeh1);
    res = BM_bvector_enumerator_free(bmeh2);

    return res;
}


int OperationsTest_AND()
{
    int res = 0;
    BM_BVHANDLE bmh0 = 0;
    BM_BVHANDLE bmh1 = 0;
    BM_BVHANDLE bmh2 = 0;
    unsigned int i;
    unsigned int count1, count2;
    int cmp;

    res = BM_bvector_construct(&bmh0, 0);
    BMERR_CHECK(res, "BM_bvector_construct()");
    res = BM_bvector_construct(&bmh1, 0);
    BMERR_CHECK(res, "BM_bvector_construct()");
    res = BM_bvector_construct(&bmh2, 0);
    BMERR_CHECK_GOTO(res, "BM_bvector_construct()", free_mem);


    for (i = 1; i < 4; ++i)
    {
        res = BM_bvector_set_bit(bmh1, i, BM_TRUE);
        BMERR_CHECK_GOTO(res, "BM_bvector_set_bit()", free_mem);
    } // for

    for (i = 0; i < 3; ++i)
    {
        res = BM_bvector_set_bit(bmh2, i, BM_TRUE);
        BMERR_CHECK_GOTO(res, "BM_bvector_set_bit()", free_mem);
    } // for

    PrintVector(bmh1, 10);
    PrintVector(bmh2, 10);
    printf("AND\n");

    res = BM_bvector_combine_AND_2sc(bmh0, bmh1, bmh2, 1);
    BMERR_CHECK_GOTO(res, "BM_bvector_combine_AND_2sc()", free_mem);

    res = BM_bvector_combine_AND(bmh1, bmh2);
    BMERR_CHECK_GOTO(res, "BM_bvector_combine_AND()", free_mem);

    res = BM_bvector_compare(bmh1, bmh0, &cmp);
    BMERR_CHECK_GOTO(res, "BM_bvector_compare()", free_mem);
    if (cmp != 0)
    {
        printf("0. incorrrect compare result %i \n", cmp);
        res = 1; goto free_mem;
    }

    res = BM_bvector_count(bmh1, &count1);
    BMERR_CHECK_GOTO(res, "BM_bvector_count()", free_mem);
    if (count1 != 2)
    {
        printf("1. incorrrect count %i \n", count1);
        res = 1; goto free_mem;
    }
    res = BM_bvector_count(bmh2, &count2);
    BMERR_CHECK_GOTO(res, "BM_bvector_count()", free_mem);
    if (count2 != 3)
    {
        printf("2. incorrrect count %i \n", count2);
        res = 1; goto free_mem;
    }

    PrintVector(bmh1, 10);
    printf("\n");
    free_mem:
        res = BM_bvector_free(bmh0);
        BMERR_CHECK(res, "bvector free failed");
        res = BM_bvector_free(bmh1);
        BMERR_CHECK(res, "bvector free failed");
        res = BM_bvector_free(bmh2);
        BMERR_CHECK(res, "bvector free failed");

    return res;
}

int OperationsTest_OR()
{
    int res = 0;
    BM_BVHANDLE bmh0 = 0;
    BM_BVHANDLE bmh1 = 0;
    BM_BVHANDLE bmh2 = 0;
    unsigned int i;
    int cmp;

    res = BM_bvector_construct(&bmh0, 0);
    BMERR_CHECK(res, "BM_bvector_construct()");
    res = BM_bvector_construct(&bmh1, 0);
    BMERR_CHECK(res, "BM_bvector_construct()");
    res = BM_bvector_construct(&bmh2, 0);
    BMERR_CHECK_GOTO(res, "BM_bvector_construct()", free_mem);


    for (i = 1; i < 4; ++i)
    {
        res = BM_bvector_set_bit(bmh1, i, BM_TRUE);
        BMERR_CHECK_GOTO(res, "BM_bvector_set_bit()", free_mem);
    } // for

    for (i = 0; i < 3; ++i)
    {
        res = BM_bvector_set_bit(bmh2, i, BM_TRUE);
        BMERR_CHECK_GOTO(res, "BM_bvector_set_bit()", free_mem);
    } // for

    PrintVector(bmh1, 10);
    PrintVector(bmh2, 10);
    printf("OR\n");

    res = BM_bvector_combine_OR_2sc(bmh0, bmh1, bmh2, 1);
    BMERR_CHECK_GOTO(res, "BM_bvector_combine_OR_2sc()", free_mem);

    res = BM_bvector_combine_OR(bmh1, bmh2);
    BMERR_CHECK_GOTO(res, "BM_bvector_combine_OR()", free_mem);

    res = BM_bvector_compare(bmh1, bmh0, &cmp);
    BMERR_CHECK_GOTO(res, "BM_bvector_compare()", free_mem);
    if (cmp != 0)
    {
        printf("0. incorrrect compare result %i \n", cmp);
        res = 1; goto free_mem;
    }

    PrintVector(bmh1, 10);
    printf("\n");

free_mem:
    res = BM_bvector_free(bmh0);
    BMERR_CHECK(res, "bvector free failed");
    res = BM_bvector_free(bmh1);
    BMERR_CHECK(res, "bvector free failed");
    res = BM_bvector_free(bmh2);
    BMERR_CHECK(res, "bvector free failed");

    return res;
}

int OperationsTest_XOR()
{
    int res = 0;
    BM_BVHANDLE bmh0 = 0;
    BM_BVHANDLE bmh1 = 0;
    BM_BVHANDLE bmh2 = 0;
    unsigned int i;
    int cmp;

    res = BM_bvector_construct(&bmh0, 0);
    BMERR_CHECK(res, "BM_bvector_construct()");
    res = BM_bvector_construct(&bmh1, 0);
    BMERR_CHECK(res, "BM_bvector_construct()");
    res = BM_bvector_construct(&bmh2, 0);
    BMERR_CHECK_GOTO(res, "BM_bvector_construct()", free_mem);


    for (i = 1; i < 4; ++i)
    {
        res = BM_bvector_set_bit(bmh1, i, BM_TRUE);
        BMERR_CHECK_GOTO(res, "BM_bvector_set_bit()", free_mem);
    } // for

    for (i = 0; i < 3; ++i)
    {
        res = BM_bvector_set_bit(bmh2, i, BM_TRUE);
        BMERR_CHECK_GOTO(res, "BM_bvector_set_bit()", free_mem);
    } // for

    PrintVector(bmh1, 10);
    PrintVector(bmh2, 10);
    printf("XOR\n");

    res = BM_bvector_combine_XOR_2sc(bmh0, bmh1, bmh2, 1);
    BMERR_CHECK_GOTO(res, "BM_bvector_combine_OR_2sc()", free_mem);

    res = BM_bvector_combine_XOR(bmh1, bmh2);
    BMERR_CHECK_GOTO(res, "BM_bvector_combine_OR()", free_mem);

    res = BM_bvector_compare(bmh1, bmh0, &cmp);
    BMERR_CHECK_GOTO(res, "BM_bvector_compare()", free_mem);
    if (cmp != 0)
    {
        printf("0. incorrrect compare result %i \n", cmp);
        res = 1; goto free_mem;
    }

    PrintVector(bmh1, 10);
    printf("\n");

free_mem:
    res = BM_bvector_free(bmh0);
    BMERR_CHECK(res, "bvector free failed");
    res = BM_bvector_free(bmh1);
    BMERR_CHECK(res, "bvector free failed");
    res = BM_bvector_free(bmh2);
    BMERR_CHECK(res, "bvector free failed");

    return res;
}

int OperationsTest_SUB()
{
    int res = 0;
    BM_BVHANDLE bmh0 = 0;
    BM_BVHANDLE bmh1 = 0;
    BM_BVHANDLE bmh2 = 0;
    unsigned int i;
    int cmp;

    res = BM_bvector_construct(&bmh0, 0);
    BMERR_CHECK(res, "BM_bvector_construct()");
    res = BM_bvector_construct(&bmh1, 0);
    BMERR_CHECK(res, "BM_bvector_construct()");
    res = BM_bvector_construct(&bmh2, 0);
    BMERR_CHECK_GOTO(res, "BM_bvector_construct()", free_mem);


    for (i = 1; i < 4; ++i)
    {
        res = BM_bvector_set_bit(bmh1, i, BM_TRUE);
        BMERR_CHECK_GOTO(res, "BM_bvector_set_bit()", free_mem);
    } // for

    for (i = 0; i < 3; ++i)
    {
        res = BM_bvector_set_bit(bmh2, i, BM_TRUE);
        BMERR_CHECK_GOTO(res, "BM_bvector_set_bit()", free_mem);
    } // for

    PrintVector(bmh1, 10);
    PrintVector(bmh2, 10);
    printf("SUB\n");

    res = BM_bvector_combine_SUB_2sc(bmh0, bmh1, bmh2, 1);
    BMERR_CHECK_GOTO(res, "BM_bvector_combine_SUB_2sc()", free_mem);

    res = BM_bvector_combine_SUB(bmh1, bmh2);
    BMERR_CHECK_GOTO(res, "BM_bvector_combine_SUB()", free_mem);

    res = BM_bvector_compare(bmh1, bmh0, &cmp);
    BMERR_CHECK_GOTO(res, "BM_bvector_compare()", free_mem);
    if (cmp != 0)
    {
        printf("0. incorrrect compare result %i \n", cmp);
        res = 1; goto free_mem;
    }

    PrintVector(bmh1, 10);
    printf("\n");

free_mem:
    res = BM_bvector_free(bmh0);
    BMERR_CHECK(res, "bvector free failed");
    res = BM_bvector_free(bmh1);
    BMERR_CHECK(res, "bvector free failed");
    res = BM_bvector_free(bmh2);
    BMERR_CHECK(res, "bvector free failed");

    return res;
}


int OperationsArrTest()
{
    unsigned int arr1[] = {100, 10, 10000};
    unsigned int arr2[] = {1,   10, 5};
    unsigned int arr3_sorted[] = {1, 5, 10, 10000};
    int res = 0;
    BM_BVHANDLE bmh1 = 0;
    unsigned int count1;
    
    res = BM_bvector_construct(&bmh1, 0);
    BMERR_CHECK(res, "BM_bvector_construct()");
    

    res = BM_bvector_combine_OR_arr(bmh1, &arr1[0], &arr1[0] + (sizeof(arr1)/sizeof(arr1[0])));
    BMERR_CHECK_GOTO(res, "BM_bvector_combine_OR_arr()", free_mem);
    res = BM_bvector_count(bmh1, &count1);
    BMERR_CHECK_GOTO(res, "BM_bvector_count()", free_mem);
    if (count1 != 3)
    {
        printf("1. incorrrect count %i \n", count1);
        res = 1; goto free_mem;
    }

    res = BM_bvector_combine_AND_arr(bmh1, &arr2[0], &arr2[0] + (sizeof(arr2)/sizeof(arr2[0])));
    BMERR_CHECK_GOTO(res, "BM_bvector_combine_AND_arr()", free_mem);
    res = BM_bvector_count(bmh1, &count1);
    BMERR_CHECK_GOTO(res, "BM_bvector_count()", free_mem);
    if (count1 != 1)
    {
        printf("2. incorrrect count %i \n", count1);
        res = 1; goto free_mem;
    }
    
    res = BM_bvector_combine_OR_arr(bmh1, &arr1[0], &arr1[0] + (sizeof(arr1)/sizeof(arr1[0])));
    BMERR_CHECK_GOTO(res, "BM_bvector_combine_OR_arr()", free_mem);
    res = BM_bvector_combine_AND_arr_sorted(bmh1, &arr3_sorted[0],
                                                  &arr3_sorted[0] + (sizeof(arr3_sorted)/sizeof(arr3_sorted[0])));
    BMERR_CHECK_GOTO(res, "BM_bvector_combine_AND_arr_sorted()", free_mem);
    res = BM_bvector_count(bmh1, &count1);
    BMERR_CHECK_GOTO(res, "BM_bvector_count()", free_mem);
    if (count1 != 2)
    {
        printf("3. incorrrect count %i \n", count1);
        res = 1; goto free_mem;
    }

    res = BM_bvector_combine_SUB_arr(bmh1, &arr1[0], &arr1[0] + (sizeof(arr1)/sizeof(arr1[0])));
    BMERR_CHECK_GOTO(res, "BM_bvector_combine_SUB_arr()", free_mem);
    res = BM_bvector_count(bmh1, &count1);
    BMERR_CHECK_GOTO(res, "BM_bvector_count()", free_mem);
    if (count1 != 0)
    {
        printf("4. incorrrect count %i \n", count1);
        res = 1; goto free_mem;
    }
    
    res = BM_bvector_combine_XOR_arr(bmh1, &arr3_sorted[0], &arr3_sorted[0] + (sizeof(arr3_sorted)/sizeof(arr3_sorted[0])));
    BMERR_CHECK_GOTO(res, "BM_bvector_combine_XOR_arr()", free_mem);
    res = BM_bvector_count(bmh1, &count1);
    BMERR_CHECK_GOTO(res, "BM_bvector_count()", free_mem);
    if (count1 != 4)
    {
        printf("5. incorrrect count %i \n", count1);
        res = 1; goto free_mem;
    }

    
    free_mem:
        res = BM_bvector_free(bmh1);
        BMERR_CHECK(res, "bvector free failed");

    return res;
}

int EnumeratorTest()
{
    int res = 0;
    BM_BVHANDLE bmh1 = 0;
    BM_BVEHANDLE bmeh1 = 0;
    BM_BVEHANDLE bmeh2 = 0;
    int valid;
    unsigned int pos;
    
    unsigned int i;

    res = BM_bvector_construct(&bmh1, 0);
    BMERR_CHECK(res, "BM_bvector_construct()");

    res = BM_bvector_enumerator_construct(bmh1, &bmeh1);
    BMERR_CHECK_GOTO(res, "BM_bvector_enumerator_construct()", free_mem);
    res = BM_bvector_enumerator_is_valid(bmeh1, &valid);
    BMERR_CHECK_GOTO(res, "BM_bvector_enumerator_is_valid()", free_mem);
    if (valid)
    {
        printf("1. incorrrect enumerator valid %i \n", valid);
        res = 1; goto free_mem;
    }
    res = BM_bvector_enumerator_free(bmeh1);


    for (i = 1; i < 4; ++i)
    {
        res = BM_bvector_set_bit(bmh1, i, BM_TRUE);
        BMERR_CHECK_GOTO(res, "BM_bvector_set_bit()", free_mem);
    } // for

    res = BM_bvector_enumerator_construct(bmh1, &bmeh1);
    BMERR_CHECK_GOTO(res, "BM_bvector_enumerator_construct()", free_mem);
    
    res = BM_bvector_enumerator_is_valid(bmeh1, &valid);
    BMERR_CHECK_GOTO(res, "BM_bvector_enumerator_is_valid()", free_mem);
    if (!valid)
    {
        printf("2. incorrrect enumerator valid %i \n", valid);
        res = 1; goto free_mem;
    }
    
    res = BM_bvector_enumerator_get_value(bmeh1, &pos);
    BMERR_CHECK_GOTO(res, "BM_bvector_enumerator_get_value()", free_mem);
    if (pos != 1)
    {
        printf("3. incorrrect enumerator traversal position %u \n", pos);
        res = 1; goto free_mem;
    }
    
    res = BM_bvector_enumerator_next(bmeh1, &valid, &pos);
    BMERR_CHECK_GOTO(res, "BM_bvector_enumerator_next()", free_mem);
    if (!valid)
    {
        printf("3. incorrrect enumerator valid %i \n", valid);
        res = 1; goto free_mem;
    }
    if (pos != 2)
    {
        printf("4. incorrrect enumerator traversal position %u \n", pos);
        res = 1; goto free_mem;
    }

    res = BM_bvector_enumerator_construct_from(bmh1, &bmeh2, 3);
    BMERR_CHECK_GOTO(res, "BM_bvector_enumerator_construct_from()", free_mem);

    res = BM_bvector_enumerator_get_value(bmeh2, &pos);
    BMERR_CHECK_GOTO(res, "BM_bvector_enumerator_get_value()", free_mem);
    if (pos != 3)
    {
        printf("4. incorrrect enumerator traversal position %u \n", pos);
        res = 1; goto free_mem;
    }
    
    res = BM_bvector_enumerator_goto(bmeh1, 3, &valid, &pos);
    BMERR_CHECK_GOTO(res, "BM_bvector_enumerator_next()", free_mem);
    if (!valid)
    {
        printf("5. incorrrect enumerator valid %i \n", valid);
        res = 1; goto free_mem;
    }
    if (pos != 3)
    {
        printf("6. incorrrect enumerator goto position %u \n", pos);
        res = 1; goto free_mem;
    }


    free_mem:
        res = BM_bvector_enumerator_free(bmeh2);
        res = BM_bvector_enumerator_free(bmeh1);
        res = BM_bvector_free(bmh1);

    return res;
}


int SerializationTest()
{
    int res = 0;
    BM_BVHANDLE bmh1 = 0;
    BM_BVHANDLE bmh2 = 0;
    BM_BVHANDLE bmh3 = 0;
    char* sbuf1 = 0;
    char* sbuf2 = 0;
    unsigned int i;
    struct BM_bvector_statistics bv_stat;
    size_t blob_size;
    int is_equal;

    res = BM_bvector_construct(&bmh1, 0);
    BMERR_CHECK(res, "BM_bvector_construct()");
    res = BM_bvector_construct(&bmh2, 0);
    BMERR_CHECK_GOTO(res, "BM_bvector_construct()", free_mem);
    res = BM_bvector_construct(&bmh3, 0);
    BMERR_CHECK_GOTO(res, "BM_bvector_construct()", free_mem);


    for (i = 1; i < 4; ++i)
    {
        res = BM_bvector_set_bit(bmh1, i, BM_TRUE);
        BMERR_CHECK_GOTO(res, "BM_bvector_set_bit()", free_mem);
        res = BM_bvector_set_bit(bmh2, i, BM_TRUE);
        BMERR_CHECK_GOTO(res, "BM_bvector_set_bit()", free_mem);
    } // for
    
    res = BM_bvector_optimize(bmh1, 3, &bv_stat);
    BMERR_CHECK_GOTO(res, "BM_bvector_calc_stat()", free_mem);

    sbuf1 = (char*) malloc(bv_stat.max_serialize_mem);
    if (sbuf1 == 0)
    {
        printf("Failed to allocate serialization buffer.\n");
        res = 1; goto free_mem;
    }
    
    res = BM_bvector_serialize(bmh1, sbuf1, bv_stat.max_serialize_mem, &blob_size);
    BMERR_CHECK_GOTO(res, "BM_bvector_serialize()", free_mem);
    
    if (blob_size == 0 || blob_size > bv_stat.max_serialize_mem)
    {
        printf("Failed to serialize correctly.\n");
        res = 1; goto free_mem;
    }
    
    sbuf2 = (char*) malloc(blob_size);
    if (sbuf1 == 0)
    {
        printf("Failed to allocate buffer.\n");
        res = 1; goto free_mem;
    }
    
    memcpy(sbuf2, sbuf1, blob_size); // imitation of I/O
    
    res = BM_bvector_deserialize(bmh3, sbuf2, blob_size);
    BMERR_CHECK_GOTO(res, "BM_bvector_deserialize()", free_mem);
    
    res = CompareVectors(bmh1, bmh2, &is_equal);
    BMERR_CHECK_GOTO(res, "CompareVectors()", free_mem);
    if (is_equal)
    {
        res = CompareVectors(bmh1, bmh3, &is_equal);
        BMERR_CHECK_GOTO(res, "CompareVectors()", free_mem);

        if (!is_equal)
        {
            printf("1.vectors comparison failed!\n");
            
            PrintVector(bmh1, 10);
            PrintVector(bmh2, 10);
            PrintVector(bmh3, 10);
            
            res = 1;
            goto free_mem;
        }
    }
    else
    {
        printf("2.vectors comparison failed!\n");

        PrintVector(bmh1, 10);
        PrintVector(bmh2, 10);
        PrintVector(bmh3, 10);
        
        res = 1;
        goto free_mem;
    }


    free_mem:
        if(sbuf1) free(sbuf1);
        if(sbuf2) free(sbuf2);

        res = BM_bvector_free(bmh1);
        res = BM_bvector_free(bmh2);
        res = BM_bvector_free(bmh3);

    return res;
}




int main(void)
{
    int res = 0;
    
    
    res = InitTest();
    if (res != 0)
    {
        printf("\nInitTest failed!\n");
        return res;
    }
    printf("\n---------------------------------- InitTest OK\n");
    
    res = ConstructionTest();
    if (res != 0)
    {
        printf("\nConstructionTest failed!\n");
        return res;
    }
    
    printf("\n---------------------------------- ConstructionTest OK\n");
    
    res = ResizeTest();
    if (res != 0)
    {
        printf("\nResizeTest failed!\n");
        return res;
    }
    printf("\n---------------------------------- ResizeTest OK\n");
    
    res = SetGetTest();
    if (res != 0)
    {
        printf("\nSetGetTest failed!\n");
        return res;
    }
    printf("\n---------------------------------- SetGetTest OK\n");
    
    res = ConstructionCopyMoveTest();
    if (res != 0)
    {
        printf("\nConstructionCopyMoveTest failed!\n");
        return res;
    }
    
    printf("\n---------------------------------- ConstructionCopyMoveTest OK\n");
    
    
    res = RangeTest();
    if (res != 0)
    {
        printf("\nRangeTest failed!\n");
        return res;
    }
    printf("\n---------------------------------- RangeTest OK\n");


    res = GetNextTest();
    if (res != 0)
    {
        printf("\nGetNextTest failed!\n");
        return res;
    }
    printf("\n---------------------------------- GetNextTest OK\n");

    res = FreezeTest();
    if (res != 0)
    {
        printf("\nFreezeTest failed!\n");
        return res;
    }
    printf("\n---------------------------------- FreezeTest OK\n");

    res = OperationsTest_AND();
    if (res != 0)
    {
        printf("\nOperationsTest AND failed!\n");
        return res;
    }

    res = OperationsTest_OR();
    if (res != 0)
    {
        printf("\nOperationsTest OR failed!\n");
        return res;
    }

    res = OperationsTest_XOR();
    if (res != 0)
    {
        printf("\nOperationsTest XOR failed!\n");
        return res;
    }

    res = OperationsTest_SUB();
    if (res != 0)
    {
        printf("\nOperationsTest SUB failed!\n");
        return res;
    }

    printf("\n---------------------------------- OperationsTest OK\n");

    res = OperationsArrTest();
    if (res != 0)
    {
        printf("\nOperationsArrTest failed!\n");
        return res;
    }
    printf("\n---------------------------------- OperationsArrTest OK\n");


    res = EnumeratorTest();
    if (res != 0)
    {
        printf("\nEnumeratorTest failed!\n");
        return res;
    }
    printf("\n---------------------------------- EnumeratorTest OK\n");


    res = SerializationTest();
    if (res != 0)
    {
        printf("\nSerializationTest failed!\n");
        return res;
    }
    printf("\n---------------------------------- SerializationTest OK\n");


    
    printf("\nlibbm unit test OK\n");
    
    return 0;
}

