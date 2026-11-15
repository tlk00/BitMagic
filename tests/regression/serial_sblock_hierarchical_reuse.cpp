/** @file serial_sblock_hierarchical_reuse.cpp
    @brief Verify allocator propagation and retained hierarchical codec storage.
*/
#undef NDEBUG
#include <cassert>
#ifndef BM_ASSERT
# define BM_ASSERT assert
#endif
#include "../common/serialization_allocator_tests.h"

/** Run hierarchical workspace reuse checks and verify allocation balance.
    @return Zero after all allocation and round-trip assertions pass.
*/
int main()
{
    serialization_allocator_test::CheckHierarchicalStateReuse();
    assert(serialization_allocator_test::tracking_block_allocator::live == 0);
}
