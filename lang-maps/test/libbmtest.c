#include <stdio.h>
#include "libbm.h"

static
int ConstructionTest()
{
    int res = 0;
    BM_BVHANDLE bmh = 0;
    
    res = BM_bvector_construct(&bmh, 200);
    if (res != BM_OK)
    {
        printf("bvector construction error \n");
        return 1;
    }
    
    res = BM_bvector_free(bmh);
    if (res != BM_OK)
    {
        printf("bvector free error \n");
        return 1;
    }
    return 0;
}

int main(void)
{
    BM_init(0);
    
    ConstructionTest();
    
    return 0;
}

