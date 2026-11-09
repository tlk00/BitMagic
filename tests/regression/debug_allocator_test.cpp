/** \file debug_allocator_test.cpp
    Focused scalar/SIMD stress-allocator checks. Build with the same backend
    defines as the stress executable. Optional arguments block-size and ptr-size
    deliberately abort to verify size-mismatch diagnostics in subprocesses.
*/
#undef NDEBUG
#include <cassert>
#include <cstdint>
#include <cstring>
#include <initializer_list>
#include "../common/debug_allocator.h"

/** Check payload alignment, writable bounds, overflow rejection and balance.
    @tparam Alloc Checked block or pointer allocator.
    @tparam ElementSize Bytes per requested allocation unit.
*/
template<class Alloc, size_t ElementSize>
void CheckAllocator()
{
    const size_t before = Alloc::balance();
    for (size_t n : {size_t(0), size_t(1), size_t(3), size_t(2048), size_t(4096)})
    {
        auto p = Alloc::allocate(n, nullptr);
        assert(reinterpret_cast<std::uintptr_t>(p) % Alloc::alignment == 0);
        assert(Alloc::balance() == before + 1);
        std::memset(p, 0xA5, n * ElementSize);
        Alloc::deallocate(p, n);
        assert(Alloc::balance() == before);
    } // for n
    bool caught = false;
    try
    {
        Alloc::allocate(std::numeric_limits<size_t>::max(), nullptr);
    }
    catch (const std::bad_alloc&)
    {
        caught = true;
    }
    assert(caught && Alloc::balance() == before);
}

/** Check aligned block/pointer allocation through containers and pool reuse.
    All owned storage must be released before checking allocator balance.
*/
void CheckContainers()
{
    using Alloc = bm::mem_alloc<dbg_block_allocator, dbg_ptr_allocator,
        bm::alloc_pool<dbg_block_allocator, dbg_ptr_allocator>>;
    {
        Alloc alloc;
        typename Alloc::allocator_pool_type pool;
        alloc.set_pool(&pool);
        bm::word_t* p = alloc.alloc_bit_block();
        assert(reinterpret_cast<std::uintptr_t>(p) % dbg_block_allocator::alignment == 0);
        alloc.free_bit_block(p);
        bm::word_t* reused = alloc.alloc_bit_block();
        assert(reused == p);
        alloc.free_bit_block(reused);
        alloc.set_pool(nullptr);

        bm::bvector<Alloc> bits;
        bits.set(1);
        bits.set(65536);
        bits.optimize();
        bm::bvector<Alloc> copy(bits);
        assert(bits == copy);
    }
    assert(dbg_block_allocator::balance() == 0);
    assert(dbg_ptr_allocator::balance() == 0);
}

/** Run positive checks, or a requested deliberate size mismatch.
    @param argc Argument count; zero or one test selector is accepted.
    @param argv Optional block-size or ptr-size selector.
    @return Zero on success, two for an unknown selector.
*/
int main(int argc, char** argv)
{
    if (argc == 2 && std::strcmp(argv[1], "block-size") == 0)
    {
        auto p = dbg_block_allocator::allocate(2048, nullptr);
        dbg_block_allocator::deallocate(p, 4096);
        return 1;
    }
    if (argc == 2 && std::strcmp(argv[1], "ptr-size") == 0)
    {
        auto p = dbg_ptr_allocator::allocate(3, nullptr);
        dbg_ptr_allocator::deallocate(p, 4);
        return 1;
    }
    if (argc != 1)
    {
        return 2;
    }
    CheckAllocator<dbg_block_allocator, sizeof(bm::word_t)>();
    CheckAllocator<dbg_ptr_allocator, sizeof(void*)>();
    CheckContainers();
#ifdef BM_STACK_COLL
    assert(g_alloc_trace_map.empty());
#endif
    std::puts("Debug allocator checks PASS");
}
