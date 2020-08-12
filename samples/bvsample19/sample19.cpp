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

/** \example sample19.cpp
  Example of bit-vector merge
  \sa bm::bvector::merge

*/

/*! \file sample19.cpp
    \brief Example: bit-vector merge

    Merge operation is an equivalent of logical OR, except 
    it destroys the argument. This is done to more efficiently 
    transfer the information to the target vector by re-assigning
    memory blocks, instead of doing allocation, copy and logical OR.

    After merge the argument vector will have undefined content, so
    it is better to clear or destroy it.

    Typical use case is to merge vectors after some multi-threaded
    partitioned processing.
*/

#include <stdlib.h>
#include <iostream>

#include "bm.h"
#include "bmundef.h" /* clear the pre-proc defines from BM */

using namespace std;

int main(void)
{
    try
    {
        bm::bvector<>   bv1 { 1, 2, 10};
        bm::bvector<>   bv2 { 10, 100000000 };

        std::cout << "Before merge bv2 count = " << bv2.count() << std::endl; // 2

        std::cout << "Merge." << std::endl;

        bv1.merge(bv2); // merge may change bv2 content (undefined content)

        std::cout << "bv1 count = " << bv1.count() << std::endl; // 4
        std::cout << "bv2 count = " << bv2.count() << std::endl; // undefined result
    }
    catch(std::exception& ex)
    {
        std::cerr << ex.what() << std::endl;
        return 1;
    }
    

    return 0;
}


