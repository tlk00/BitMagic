#ifndef BM_TEST_DEBUG_ALLOCATOR_H_INCLUDED
#define BM_TEST_DEBUG_ALLOCATOR_H_INCLUDED

/** \file debug_allocator.h
    Shared stress-test allocators with checked deallocation sizes.
    Internal test tooling only. Backend allocation and an alignment-sized
    metadata prefix preserve SIMD alignment, including Mac NEON's alignment.
    Define BM_STACK_COLL to retain allocation backtraces (expensive).
*/
#include "bm.h"
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <mutex>
#include <new>
#include <string>
#include <unordered_map>
#include <utility>

#ifdef BM_STACK_COLL
#include "stacktrace_dbg.h"
inline std::unordered_map<void*, std::string> g_alloc_trace_map;
#endif
inline std::mutex g_trace_lock;

/** @internal Test allocator preserving the backend's alignment and size units.
    @tparam Backend BM block or pointer allocator providing matching allocate/free.
    @tparam ElementSize Bytes in one backend allocation unit.
*/
template<class Backend, size_t ElementSize>
class size_checked_allocator
{
public:
    using pointer = decltype(Backend::allocate(size_t(0), nullptr));
#ifdef BM_ALLOC_ALIGN
    static constexpr size_t alignment = BM_ALLOC_ALIGN; ///< Backend SIMD alignment.
#else
    static constexpr size_t alignment = alignof(std::max_align_t); ///< malloc alignment.
#endif
    static constexpr size_t prefix_bytes =
        ((sizeof(size_t) + alignment - 1) / alignment) * alignment; ///< Padded metadata.
    static inline size_t na_ = 0; ///< Successful allocations.
    static inline size_t nf_ = 0; ///< Successful deallocations.

    /** Allocate aligned payload and record its requested size.
        @param n Payload length in backend units.
        @param hint Passed to the backend without taking ownership.
        @return Owned payload; release with deallocate and the same n.
        @throws std::bad_alloc On size overflow or backend allocation failure.
    */
    static pointer allocate(size_t n, const void* hint);
    /** Verify size and free the original backend allocation.
        @param p Live payload returned by allocate; must not be null.
        @param n Original requested payload length in backend units.
        Terminates on mismatch. Invalid pointers and double frees are not supported.
    */
    static void deallocate(pointer p, size_t n);
    /** @return Outstanding allocations; callers must serialize concurrent access. */
    static size_t balance() BMNOEXCEPT { return na_ - nf_; }
};

template<class Backend, size_t ElementSize>
typename size_checked_allocator<Backend, ElementSize>::pointer
size_checked_allocator<Backend, ElementSize>::allocate(size_t n, const void* hint)
{
    static_assert(prefix_bytes % ElementSize == 0, "Metadata must use whole units");
    static_assert(alignment >= alignof(size_t), "Metadata must be aligned");
    constexpr size_t prefix_units = prefix_bytes / ElementSize;
    if (n > (std::numeric_limits<size_t>::max() - prefix_bytes) / ElementSize)
    {
        throw std::bad_alloc();
    }
#ifdef BM_STACK_COLL
    std::string stack;
    get_stacktrace(stack);
#endif
    std::lock_guard<std::mutex> lock(g_trace_lock);
    auto* base = reinterpret_cast<unsigned char*>(Backend::allocate(n + prefix_units, hint));
    auto* payload = base + prefix_bytes;
    *reinterpret_cast<size_t*>(payload - sizeof(size_t)) = n;
#ifdef BM_STACK_COLL
    try
    {
        g_alloc_trace_map.emplace(base, std::move(stack));
    }
    catch (...)
    {
        Backend::deallocate(reinterpret_cast<pointer>(base), n + prefix_units);
        throw;
    }
#endif
    ++na_;
    return reinterpret_cast<pointer>(payload);
}

template<class Backend, size_t ElementSize>
void size_checked_allocator<Backend, ElementSize>::deallocate(pointer p, size_t n)
{
    std::lock_guard<std::mutex> lock(g_trace_lock);
    auto* payload = reinterpret_cast<unsigned char*>(p);
    const size_t recorded = *reinterpret_cast<const size_t*>(payload - sizeof(size_t));
    if (recorded != n)
    {
        std::fprintf(stderr,
            "Debug allocator deallocation ERROR! pointer=%p unit=%zu n=%zu recorded=%zu\n",
            static_cast<void*>(p), ElementSize, n, recorded);
#ifdef BM_STACK_COLL
        std::string stack;
        get_stacktrace(stack);
        std::fputs(stack.c_str(), stderr);
#endif
        std::abort();
    }
    auto* base = payload - prefix_bytes;
#ifdef BM_STACK_COLL
    g_alloc_trace_map.erase(base);
#endif
    Backend::deallocate(reinterpret_cast<pointer>(base), n + prefix_bytes / ElementSize);
    ++nf_;
}

using dbg_block_allocator = size_checked_allocator<bm::block_allocator, sizeof(bm::word_t)>;
using dbg_ptr_allocator = size_checked_allocator<bm::ptr_allocator, sizeof(void*)>;

#endif // BM_TEST_DEBUG_ALLOCATOR_H_INCLUDED
