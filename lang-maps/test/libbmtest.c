#include <stdio.h>
#include "libbm.h"



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
        int cmp;
        
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
    BM_BVHANDLE bmh = 0;
    int val;
    unsigned int count;
    
    res = BM_bvector_construct(&bmh, 200);
    BMERR_CHECK(res, "BM_bvector_construct()");
    
    res = BM_bvector_any(bmh, &val);
    BMERR_CHECK_GOTO(res, "BM_bvector_any()", free_mem);
    if (val)
    {
        printf("bvector any() incorrect value \n");
        res = 1; goto free_mem;
    }
    
    res = BM_bvector_set_bit(bmh, 10, BM_TRUE);
    BMERR_CHECK_GOTO(res, "BM_bvector_set_bit()", free_mem);
    
    res = BM_bvector_any(bmh, &val);
    BMERR_CHECK_GOTO(res, "BM_bvector_any()", free_mem);
    if (val==0)
    {
        printf("bvector any() incorrect value \n");
        res = 1; goto free_mem;
    }
    
    
    res = BM_bvector_get_bit(bmh, 10, &val);
    BMERR_CHECK_GOTO(res, "BM_bvector_get_bit()", free_mem);
    if (!val)
    {
        printf("bvector get_bit incorrect value \n");
        res = 1; goto free_mem;
    }

    res = BM_bvector_get_bit(bmh, 0, &val);
    BMERR_CHECK_GOTO(res, "BM_bvector_get_bit()", free_mem);
    if (val)
    {
        printf("bvector get_bit incorrect value \n");
        res = 1; goto free_mem;
    }
    
    res = BM_bvector_flip_bit(bmh, 0);
    BMERR_CHECK_GOTO(res, "BM_bvector_flip_bit()", free_mem);

    res = BM_bvector_get_bit(bmh, 0, &val);
    BMERR_CHECK_GOTO(res, "BM_bvector_get_bit()", free_mem);
    if (val == 0)
    {
        printf("bvector get_bit incorrect value \n");
        res = 1; goto free_mem;
    }
    res = BM_bvector_flip_bit(bmh, 0);
    BMERR_CHECK_GOTO(res, "BM_bvector_flip_bit()", free_mem);
    
    
    res = BM_bvector_count(bmh, &count);
    BMERR_CHECK_GOTO(res, "BM_bvector_count()", free_mem);
    if (count != 1)
    {
        printf("incorrrect count %i \n", count);
        res = 1; goto free_mem;
    }
    
    res = BM_bvector_set_bit(bmh, 10, BM_FALSE);
    BMERR_CHECK_GOTO(res, "BM_bvector_set_bit()", free_mem);
    
    res = BM_bvector_get_bit(bmh, 0, &val);
    BMERR_CHECK_GOTO(res, "BM_bvector_get_bit()", free_mem);
    if (val != BM_FALSE)
    {
        printf("bvector get_bit incorrect value %i\n", val);
        res = 1; goto free_mem;
    }
    res = BM_bvector_count(bmh, &count);
    BMERR_CHECK_GOTO(res, "BM_bvector_count()", free_mem);
    if (count != 0)
    {
        printf("incorrrect count %i \n", count);
        res = 1; goto free_mem;
    }
    
    {
        int change;
        res = BM_bvector_set_bit_conditional(bmh, 0, BM_TRUE, BM_FALSE, &change);
        BMERR_CHECK_GOTO(res, "BM_bvector_set_bit_conditional()", free_mem);
        if (!change)
        {
            printf("bvector set_bit_conditional error \n");
            res = 1; goto free_mem;
        }
        res = BM_bvector_count(bmh, &count);
        BMERR_CHECK_GOTO(res, "BM_bvector_count()", free_mem);
        if (count != 1)
        {
            printf("incorrrect count %i \n", count);
            res = 1; goto free_mem;
        }
        res = BM_bvector_set_bit_conditional(bmh, 0, BM_TRUE, BM_FALSE, &change);
        BMERR_CHECK_GOTO(res, "BM_bvector_set_bit_conditional()", free_mem);
        if (change)
        {
            printf("bvector set_bit_conditional error \n");
            res = 1; goto free_mem;
        }
    }
    
    res = BM_bvector_set(bmh);
    BMERR_CHECK_GOTO(res, "BM_bvector_set()", free_mem);

    res = BM_bvector_count(bmh, &count);
    BMERR_CHECK_GOTO(res, "BM_bvector_count()", free_mem);
    if (count == 0)
    {
        printf("incorrrect count %i \n", count);
        res = 1; goto free_mem;
    }
    res = BM_bvector_clear(bmh, BM_TRUE);
    BMERR_CHECK_GOTO(res, "BM_bvector_clear()", free_mem);

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
    
    res = BM_bvector_set_bit(bmh, 0, BM_TRUE);
    BMERR_CHECK_GOTO(res, "BM_bvector_set_bit()", free_mem);
    
    res = BM_bvector_get_first(bmh, &idx, &found);
    BMERR_CHECK_GOTO(res, "BM_bvector_get_first()", free_mem);
    if (!found || idx != 0)
    {
        printf("2. incorrrect get first %i\n", idx);
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
    
    
    printf("\nlibbm unit test OK\n");
    
    return 0;
}

