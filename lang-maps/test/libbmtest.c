#include <stdio.h>
#include "libbm.h"





static
int InitTest()
{
    int res = 0;
    const char* c;
    int major, minor, patch;
    
    const char* msg;
    
    res = BM_init(0);
    if (res != BM_OK)
    {
        printf("BitMagic initialization failed \n");
        return 1;
    }
    
    c = BM_version(&major, &minor, &patch);
    if (!c)
    {
        printf("BM_version() failed \n");
        return 1;
    }
    
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
    
    res = BM_bvector_construct(&bmh, 200);
    if (res != BM_OK)
    {
        printf("bvector construction error \n");
        return res;
    }
    
    res = BM_bvector_free(bmh);
    if (res != BM_OK)
    {
        printf("bvector free error \n");
        return res;
    }
    return 0;
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
    if (res != BM_OK)
    {
        printf("bvector construction error \n");
        return res;
    }
    
    res = BM_bvector_get_size(bmh, &size);
    if (res != BM_OK)
    {
        printf("bvector get size error %s\n", BM_error_msg(res));
        return res;
    }
    if (size != size1)
    {
        printf("bvector get size failed %i\n", size);
    }
    
    res = BM_bvector_set_size(bmh, size2);
    if (res != BM_OK)
    {
        printf("bvector set size error %s\n", BM_error_msg(res));
        return res;
    }
    
    res = BM_bvector_get_size(bmh, &size);
    if (res != BM_OK)
    {
        printf("bvector get size error %s\n", BM_error_msg(res));
        return res;
    }
    if (size != size2)
    {
        printf("bvector get size failed %i\n", size);
    }
    
    
    
    
    res = BM_bvector_free(bmh);
    if (res != BM_OK)
    {
        printf("bvector free error \n");
        return res;
    }
    return 0;
}


static
int SetGetTest()
{
    int res = 0;
    BM_BVHANDLE bmh = 0;
    int val;
    unsigned int count;
    
    res = BM_bvector_construct(&bmh, 200);
    if (res != BM_OK)
    {
        printf("bvector construction error \n");
        return res;
    }
    
    res = BM_bvector_set_bit(bmh, 10, BM_TRUE);
    if (res != BM_OK)
    {
        printf("bvector set_bit error \n");
        return res;
    }
    
    res = BM_bvector_get_bit(bmh, 10, &val);
    if (res != BM_OK)
    {
        printf("bvector get_bit error \n");
        return res;
    }
    if (!val)
    {
        printf("bvector get_bit incorrect value \n");
        return 1;
    }

    res = BM_bvector_get_bit(bmh, 0, &val);
    if (res != BM_OK)
    {
        printf("bvector get_bit error \n");
        return res;
    }
    if (val)
    {
        printf("bvector get_bit incorrect value \n");
        return 1;
    }
    
    res = BM_bvector_count(bmh, &count);
    if (res != BM_OK)
    {
        printf("bvector count error \n");
        return res;
    }
    if (count != 1)
    {
        printf("incorrrect count %i \n", count);
        return res;
    }
    
    res = BM_bvector_set_bit(bmh, 10, BM_FALSE);
    if (res != BM_OK)
    {
        printf("bvector set_bit error \n");
        return res;
    }
    res = BM_bvector_get_bit(bmh, 0, &val);
    if (res != BM_OK)
    {
        printf("bvector get_bit error \n");
        return res;
    }
    if (val != BM_FALSE)
    {
        printf("bvector get_bit incorrect value %i\n", val);
        return 1;
    }
    res = BM_bvector_count(bmh, &count);
    if (res != BM_OK)
    {
        printf("bvector count error \n");
        return res;
    }
    if (count != 0)
    {
        printf("incorrrect count %i \n", count);
        return res;
    }
    
    {
        int change;
        res = BM_bvector_set_bit_conditional(bmh, 0, BM_TRUE, BM_FALSE, &change);
        if (res != BM_OK)
        {
            printf("bvector set_bit error \n");
            return res;
        }
        
        if (!change)
        {
            printf("bvector set_bit_conditional error \n");
            return 10;
        }
        res = BM_bvector_count(bmh, &count);
        if (res != BM_OK)
        {
            printf("bvector count error \n");
            return res;
        }
        if (count != 1)
        {
            printf("incorrrect count %i \n", count);
            return res;
        }
        res = BM_bvector_set_bit_conditional(bmh, 0, BM_TRUE, BM_FALSE, &change);
        if (res != BM_OK)
        {
            printf("bvector set_bit error \n");
            return res;
        }
        if (change)
        {
            printf("bvector set_bit_conditional error \n");
            return 10;
        }
    }
    
    res = BM_bvector_set(bmh);
    if (res != BM_OK)
    {
        printf("bvector set error \n");
        return res;
    }
    res = BM_bvector_count(bmh, &count);
    if (res != BM_OK)
    {
        printf("bvector count error \n");
        return res;
    }
    if (count == 0)
    {
        printf("incorrrect count %i \n", count);
        return res;
    }
    res = BM_bvector_clear(bmh, BM_TRUE);
    if (res != BM_OK)
    {
        printf("bvector clear error \n");
        return res;
    }
    res = BM_bvector_count(bmh, &count);
    if (res != BM_OK)
    {
        printf("bvector count error \n");
        return res;
    }
    if (count != 0)
    {
        printf("incorrrect count %i \n", count);
        return res;
    }
    
    
    
    
    
    res = BM_bvector_free(bmh);
    if (res != BM_OK)
    {
        printf("bvector free error \n");
        return res;
    }
    return 0;

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
    
    
    
    printf("\nlibbm unit test OK\n");
    
    return 0;
}

