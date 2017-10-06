#include <stdio.h>
#include "libbm.h"


static
int InitTest()
{
    int res = 0;
    const char* c;
    unsigned int major, minor, patch;
    
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
    
    return 0;
}

