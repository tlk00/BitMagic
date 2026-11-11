#ifndef BMBVHASH__H__INCLUDED__
#define BMBVHASH__H__INCLUDED__
/*
Copyright(c) 2026 Anatoliy Kuznetsov(anatoliy_kuznetsov at yahoo.com)

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


#ifndef BM__H__INCLUDED__
# error missing include (bm.h or bm64.h)
#endif
#include "bmdef.h"
#include "bmfunc.h"

namespace bm
{

/** Canonical block traversal with a caller-provided streaming hash functor.
    @tparam BV Bit-vector type supplying its allocator and block manager.
    @tparam Hash Functor implementing reset(), update(const void*, size_t),
        final(unsigned char*), and a positive static digest_size constant.
    The functor is borrowed and must outlive this non-copyable object. Calls on
    one instance are sequential. Independent instances may run concurrently.
    Each compute resets the functor, hashes ASCII BMBLOCKSHA256V1, then every
    nonzero logical block as 8192 LE bytes plus its absolute uint64 LE index.
    Allocated zero blocks and trailing logical size do not affect set content.
    GAP blocks expand into aligned scratch. FULL blocks supply all-one bytes.
    Two allocator-owned blocks are retained for GAP and byte-order scratch.
    @ingroup setalgo
*/
template<class BV, class Hash>
class bvector_hash
{
    static constexpr bool is_nothrow_hash_ =
        noexcept(static_cast<Hash*>(nullptr)->reset()) &&
        noexcept(static_cast<Hash*>(nullptr)->update(static_cast<const void*>(nullptr), size_t(0))) &&
        noexcept(static_cast<Hash*>(nullptr)->final(static_cast<unsigned char*>(nullptr)));
        ///< Exception guarantee of the exact functor operations used by compute.
public:
    typedef typename BV::allocator_type allocator_type; ///< Propagated allocator.
    /** Borrow hash and copy allocator. Allocate scratch or propagate failure.
        @param hash Mutable hash functor, exclusively borrowed during compute.
        @param allocator Allocator for both retained scratch blocks. Supply the
            source vector's get_allocator() to preserve its allocation policy.
    */
    explicit bvector_hash(Hash& hash, const allocator_type& allocator = allocator_type());
    /** Release retained scratch through the original allocator. */
    ~bvector_hash();
    bvector_hash(const bvector_hash&) = delete;
    bvector_hash& operator=(const bvector_hash&) = delete;
    /** Reset hash state while preserving allocator and scratch. Exceptions from
        the functor propagate; another init/compute starts a new hash stream.
        This operation is noexcept exactly when Hash::reset is noexcept.
    */
    void init() noexcept(noexcept(static_cast<Hash*>(nullptr)->reset()));
    /** Hash bv's logical set and write Hash::digest_size bytes to digest.
        @param bv Read-only source; no concurrent mutation is allowed.
        @param digest Caller-owned non-null output of at least digest_size bytes.
        Conditionally noexcept when reset, update(const void*, size_t) and final
        are noexcept. Otherwise functor exceptions propagate and scratch remains
        valid for the next call. Traversal performs no allocation.
    */
    void compute(const BV& bv, unsigned char* digest) noexcept(is_nothrow_hash_);
private:
    Hash&          hash_;       ///< Borrowed streaming hash functor.
    allocator_type allocator_;  ///< Scratch allocation policy.
    bm::word_t*    scratch_;    ///< Two aligned allocator blocks, retained until destruction.
};

template<class BV, class Hash>
bvector_hash<BV, Hash>::bvector_hash(Hash& hash, const allocator_type& allocator)
    : hash_(hash), allocator_(allocator), scratch_(allocator_.alloc_bit_block(2))
{
    static_assert(Hash::digest_size > 0, "Hash digest must have positive size");
    static_assert(sizeof(bm::word_t) == 4, "Canonical format requires 32-bit words");
}

template<class BV, class Hash>
bvector_hash<BV, Hash>::~bvector_hash()
{
    allocator_.free_bit_block(scratch_, 2);
}

template<class BV, class Hash>
void bvector_hash<BV, Hash>::init() noexcept(noexcept(static_cast<Hash*>(nullptr)->reset()))
{
    hash_.reset();
}

template<class BV, class Hash>
void bvector_hash<BV, Hash>::compute(const BV& bv, unsigned char* digest) noexcept(is_nothrow_hash_)
{
    BM_ASSERT(digest);
    init();
    const char domain[] = "BMBLOCKSHA256V1";
    hash_.update(static_cast<const void*>(domain), sizeof(domain) - 1);
    const auto& manager = bv.get_blocks_manager();
    if (manager.is_init())
    {
        const unsigned endian = 1;
        const bool little_endian = *reinterpret_cast<const unsigned char*>(&endian) == 1;
        unsigned char* canonical = reinterpret_cast<unsigned char*>(scratch_ + bm::set_block_size);
        for (unsigned top = 0; top < manager.top_block_size(); ++top)
        {
            const auto* blocks = manager.get_topblock(top);
            if (!blocks)
            {
                continue;
            }
            const bool full_parent = reinterpret_cast<const bm::word_t*>(blocks) == FULL_BLOCK_FAKE_ADDR;
            for (unsigned sub = 0; sub < bm::set_sub_array_size; ++sub)
            {
                const bm::word_t* block = full_parent ? FULL_BLOCK_REAL_ADDR : blocks[sub];
                if (!block)
                {
                    continue;
                }
                const bm::id64_t number = bm::id64_t(top) * bm::set_sub_array_size + sub;
                if (block == FULL_BLOCK_FAKE_ADDR)
                {
                    block = FULL_BLOCK_REAL_ADDR;
                }
                else if (BM_IS_GAP(block))
                {
                    bm::gap_convert_to_bitset(scratch_, BMGAP_PTR(block));
                    block = scratch_;
                }
                if (bm::bit_is_all_zero(block))
                {
                    continue;
                }
                if (little_endian)
                {
                    hash_.update(static_cast<const void*>(block), bm::set_block_size * sizeof(bm::word_t));
                }
                else
                {
                    for (unsigned word = 0; word < bm::set_block_size; ++word)
                    {
                        for (unsigned byte = 0; byte < 4; ++byte)
                        {
                            canonical[word * 4 + byte] = static_cast<unsigned char>(block[word] >> (byte * 8));
                        } // for byte
                    } // for word
                    hash_.update(static_cast<const void*>(canonical), bm::set_block_size * sizeof(bm::word_t));
                }
                unsigned char postfix[8];
                for (unsigned byte = 0; byte < 8; ++byte)
                {
                    postfix[byte] = static_cast<unsigned char>(number >> (byte * 8));
                } // for byte
                hash_.update(static_cast<const void*>(postfix), sizeof(postfix));
            } // for sub
        } // for top
    }
    hash_.final(digest);
}

} // namespace bm
#include "bmundef.h"
#endif
