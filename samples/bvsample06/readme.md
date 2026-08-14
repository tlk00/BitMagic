# bvsample06: Custom memory allocators

This example supplies `bm::bvector<>` with custom block and pointer allocators
through `bm::mem_alloc<>`.

The debug allocators count allocation and deallocation calls and reserve a
small header in each allocation. After a vector is created, modified and
destroyed, assertions verify that both allocator balances return to zero.

The implementation is deliberately simple and is meant to explain the
allocator interface, not to serve as a production allocator. In particular,
it does not provide the address alignment required by BitMagic SIMD builds.
