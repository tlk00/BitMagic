/** Standalone generic hash traversal checks, including reset and storage shapes.
    Build with -I../../src from utils/inv_list. Also supports -DBM64ADDR.
*/
#include <cassert>
#include <cstring>
#include <vector>
#include <stdexcept>
#include <utility>
#include "bm.h"
#include "bmbvhash.h"

/** @internal Test functor retaining canonical bytes for an independent oracle. */
struct recording_hash
{
    static const unsigned digest_size = 1; ///< Dummy digest output length.
    std::vector<unsigned char> bytes; ///< Captured stream for inspection.
    unsigned resets = 0; ///< Number of explicit resets.
    /** Clear captured bytes and record a new run. */
    void reset()
    {
        bytes.clear();
        ++resets;
    }
    /** Append size borrowed bytes from data to the test capture. */
    void update(const void* data, size_t size)
    {
        const auto* source = static_cast<const unsigned char*>(data);
        bytes.insert(bytes.end(), source, source + size);
    }
    /** Write one dummy byte to caller-owned output. */
    void final(unsigned char* output) { *output = 0; }
};

/** @internal Probe each functor method's independent exception specification. */
template<bool ResetNoexcept, bool UpdateNoexcept, bool FinalNoexcept>
struct exception_hash
{
    static const unsigned digest_size = 1; ///< Dummy digest size.
    unsigned fail_phase = 0; ///< Inject failure in reset=1, update=2 or final=3.
    /** Reset this dummy stream, optionally throwing the injected reset failure. */
    void reset() noexcept(ResetNoexcept)
    {
        if constexpr (!ResetNoexcept)
        {
            if (fail_phase == 1)
            {
                throw std::runtime_error("reset failure");
            }
        }
    }
    /** Ignore size borrowed bytes at data, optionally failing this update. */
    void update(const void* data, size_t size) noexcept(UpdateNoexcept)
    {
        (void)data;
        (void)size;
        if constexpr (!UpdateNoexcept)
        {
            if (fail_phase == 2)
            {
                throw std::runtime_error("update failure");
            }
        }
    }
    /** Write a dummy digest to output, optionally failing finalization. */
    void final(unsigned char* output) noexcept(FinalNoexcept)
    {
        if constexpr (!FinalNoexcept)
        {
            if (fail_phase == 3)
            {
                throw std::runtime_error("final failure");
            }
        }
        *output = 0;
    }
};

using NoThrowHash = bm::bvector_hash<bm::bvector<>, exception_hash<true, true, true>>;
using ThrowResetHash = bm::bvector_hash<bm::bvector<>, exception_hash<false, true, true>>;
using ThrowUpdateHash = bm::bvector_hash<bm::bvector<>, exception_hash<true, false, true>>;
using ThrowFinalHash = bm::bvector_hash<bm::bvector<>, exception_hash<true, true, false>>;
static_assert(noexcept(std::declval<NoThrowHash&>().compute(std::declval<const bm::bvector<>&>(), nullptr)),
              "All non-throwing hash operations make compute noexcept");
static_assert(!noexcept(std::declval<ThrowResetHash&>().compute(std::declval<const bm::bvector<>&>(), nullptr)),
              "Reset failures must propagate");
static_assert(!noexcept(std::declval<ThrowUpdateHash&>().compute(std::declval<const bm::bvector<>&>(), nullptr)),
              "Update failures must propagate");
static_assert(!noexcept(std::declval<ThrowFinalHash&>().compute(std::declval<const bm::bvector<>&>(), nullptr)),
              "Final failures must propagate");
static_assert(noexcept(std::declval<ThrowUpdateHash&>().init()), "Init only depends on reset");
static_assert(!noexcept(std::declval<ThrowResetHash&>().init()), "Init propagates reset failures");

/** Inject phase into Hash, verify propagation and successful reuse after failure. */
template<class Hash>
void test_hash_exception(unsigned phase)
{
    Hash hash;
    bm::bvector<> bv;
    bv.set(3);
    bm::bvector_hash<bm::bvector<>, Hash> traversal(hash, bv.get_allocator());
    unsigned char digest = 0;
    hash.fail_phase = phase;
    bool caught = false;
    try
    {
        traversal.compute(bv, &digest);
    }
    catch (const std::runtime_error&)
    {
        caught = true;
    }
    assert(caught);
    hash.fail_phase = 0;
    traversal.compute(bv, &digest);
}

/** @internal Allocator probe verifying two-block workspace ownership. */
struct recording_allocator : bm::bvector<>::allocator_type
{
    using base_type = bm::bvector<>::allocator_type; ///< Normal allocator implementation.
    static unsigned allocated; ///< Charged bit blocks.
    static unsigned freed; ///< Released bit blocks.
    /** Allocate factor aligned blocks and record the allocation count. */
    bm::word_t* alloc_bit_block(unsigned factor = 1);
    /** Free block with its original factor and record the release count. */
    void free_bit_block(bm::word_t* block, size_t factor = 1);
};
unsigned recording_allocator::allocated = 0;
unsigned recording_allocator::freed = 0;

bm::word_t* recording_allocator::alloc_bit_block(unsigned factor)
{
    allocated += factor;
    return base_type::alloc_bit_block(factor);
}

void recording_allocator::free_bit_block(bm::word_t* block, size_t factor)
{
    freed += unsigned(factor);
    base_type::free_bit_block(block, factor);
}

/** Exercise equivalent bit/GAP/full storage, holes, empty sets and repeated calls. */
int main()
{
    test_hash_exception<exception_hash<false, true, true>>(1);
    test_hash_exception<exception_hash<true, false, true>>(2);
    test_hash_exception<exception_hash<true, true, false>>(3);
    using BV = bm::bvector<>;
    recording_hash hash;
    BV bv;
    bm::bvector_hash<BV, recording_hash> traversal(hash, bv.get_allocator());
    unsigned char digest = 0;
    traversal.compute(bv, &digest);
    assert(hash.bytes.size() == 15);
    assert(std::memcmp(hash.bytes.data(), "BMBLOCKSHA256V1", 15) == 0);
    bv.set_range(65536, 131071);
    traversal.compute(bv, &digest);
    const auto full = hash.bytes;
    assert(full.size() == 15 + 8192 + 8);
    for (unsigned i = 15; i < 15 + 8192; ++i)
    {
        assert(full[i] == 255);
    } // for i
    assert(full[15 + 8192] == 1);
    for (unsigned i = 15 + 8192 + 1; i < full.size(); ++i)
    {
        assert(full[i] == 0);
    } // for i
    bv.optimize();
    traversal.compute(bv, &digest);
    assert(hash.bytes == full);
    BV gap(bm::BM_GAP);
    gap.set(3);
    gap.set(196609);
    traversal.compute(gap, &digest);
    const auto sparse = hash.bytes;
    assert(sparse.size() == 15 + 2 * (8192 + 8));
    assert(sparse[15] == 8);
    assert(sparse[15 + 8192 + 8] == 2);
    assert(sparse[15 + 8192 + 8 + 8192] == 3);
    BV bit(bm::BM_BIT);
    bit.set(3);
    bit.set(196609);
    traversal.compute(bit, &digest);
    assert(hash.bytes == sparse);
    bit.set(65539);
    bit.clear_bit(65539);
    traversal.compute(bit, &digest);
    assert(hash.bytes == sparse);
    bit.clear();
    traversal.compute(bit, &digest);
    assert(hash.bytes.size() == 15);
    assert(hash.resets == 7);
#ifdef BM64ADDR
    BV high;
    high.set((BV::size_type(1) << 40) + 7);
    traversal.compute(high, &digest);
    assert(hash.bytes.size() == 15 + 8192 + 8);
    assert(hash.bytes[15] == 128);
    assert(hash.bytes[15 + 8192 + 3] == 1);
#endif
    BV parent_full;
    parent_full.set_range(0, BV::size_type(256) * 65536 - 1);
    traversal.compute(parent_full, &digest);
    assert(hash.bytes.size() == 15 + 256 * (8192 + 8));
    const auto parent_bytes = hash.bytes;
    parent_full.optimize();
    traversal.compute(parent_full, &digest);
    assert(hash.bytes == parent_bytes);
    exception_hash<true, true, true> no_throw_hash;
    NoThrowHash no_throw_traversal(no_throw_hash, parent_full.get_allocator());
    no_throw_traversal.compute(parent_full, &digest);
    using CustomBV = bm::bvector<recording_allocator>;
    CustomBV custom;
    custom.set(9);
    const unsigned before_alloc = recording_allocator::allocated;
    const unsigned before_free = recording_allocator::freed;
    {
        bm::bvector_hash<CustomBV, recording_hash> custom_hash(hash, custom.get_allocator());
        assert(recording_allocator::allocated == before_alloc + 2);
        custom_hash.compute(custom, &digest);
        custom_hash.compute(custom, &digest);
        assert(recording_allocator::allocated == before_alloc + 2);
    }
    assert(recording_allocator::freed == before_free + 2);
    return 0;
}
