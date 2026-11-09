#ifndef BM_TEST_NEON__H_INCLUDED__
#define BM_TEST_NEON__H_INCLUDED__
// Shared by the 32- and 64-bit stress drivers. The oracle never calls BM's
// counting, shifting or Boolean block implementations.
#ifdef BMNEONOPT
#include <algorithm>
#include <cstring>
#include <cstdint>
#include <iostream>
#include <random>
#include <vector>
#include "bmserial.h"
#include "bmalgo.h"
#include "bmbvimport.h"
#include "bmbuffer.h"

static void NeonCheck(bool ok, const char* what)
{
    if (!ok) { std::cerr << "NEON check failed: " << what << std::endl; std::abort(); }
}

static unsigned NeonScalarCount(unsigned w)
{
    unsigned n = 0;
    for (; w; w &= w-1) ++n;
    return n;
}

// Exercise the allocator contract independently of SIMD load requirements.
static void TestNEONAllocators()
{
    auto check_alignment = [](const void* ptr) {
        NeonCheck(ptr && reinterpret_cast<std::uintptr_t>(ptr) % 128 == 0,
                  "allocator 128-byte alignment");
    };
    for (unsigned n : {1u, 3u, 15u, 16u, 17u, 31u, 32u, 33u, 127u, 129u, 2048u})
    {
        bm::word_t* block = bm::block_allocator::allocate(n, 0);
        check_alignment(block);
        std::memset(block, 0xA5, n * sizeof(*block));
        bm::block_allocator::deallocate(block, n);

        void* ptr = bm::ptr_allocator::allocate(n, 0);
        check_alignment(ptr);
        std::memset(ptr, 0x5A, n * sizeof(void*));
        bm::ptr_allocator::deallocate(ptr, n);

        void* bytes = bm::aligned_new_malloc(n);
        check_alignment(bytes);
        std::memset(bytes, 0xA5, n);
        bm::aligned_free(bytes);
    } // for n

    bm::standard_allocator alloc;
    for (unsigned level = 0; level < bm::gap_levels; ++level)
    {
        bm::gap_word_t* gap = alloc.alloc_gap_block(level, bm::gap_len_table<true>::_len);
        check_alignment(gap);
        gap[0] = bm::gap_word_t(level << 1);
        alloc.free_gap_block(gap, bm::gap_len_table<true>::_len);
    } // for level

    bm::word_t* chunk = alloc.alloc_bit_block(3);
    for (unsigned i = 0; i < 3; ++i)
        check_alignment(chunk + i * bm::set_block_size);
    alloc.free_bit_block(chunk, 3);

    bm::standard_alloc_pool pool;
    alloc.set_pool(&pool);
    bm::word_t* block = alloc.alloc_bit_block();
    check_alignment(block);
    alloc.free_bit_block(block);
    bm::word_t* reused = alloc.alloc_bit_block();
    NeonCheck(reused == block, "allocator pool reuse");
    check_alignment(reused);
    alloc.free_bit_block(reused);
    alloc.set_pool(0);

    bm::byte_buffer<bm::standard_allocator> buffer;
    for (unsigned n : {1u, 129u, 8193u, 32769u})
    {
        buffer.resize(n);
        check_alignment(buffer.data());
        std::memset(buffer.data(), 0x5A, n);
    } // for n
    bm::byte_buffer<bm::standard_allocator> copy(buffer);
    check_alignment(copy.data());
    NeonCheck(std::memcmp(copy.data(), buffer.data(), buffer.size()) == 0,
              "aligned buffer copy");
}

static void TestNEONUtils()
{
    TestNEONAllocators();
    using namespace bm;
    const unsigned N = bm::set_block_size, W = bm::set_block_digest_wave_size;
    std::mt19937 rng(0x4e454f4eu);
    std::vector<unsigned> a(N), b(N), c(N), e(N), dst(N), ref(N);
    auto check_block = [&](const char* label) { NeonCheck(dst == ref, label); };
    std::cout << "NEON kernels: adversarial + random (seed 0x4e454f4e)" << std::endl;
    for (unsigned trial = 0; trial < 160; ++trial)
    {
        for (unsigned i = 0; i < N; ++i)
        {
            a[i] = unsigned(rng()); b[i] = unsigned(rng());
            c[i] = unsigned(rng()); e[i] = unsigned(rng());
            switch (trial)
            {
            case 0: a[i] = b[i] = c[i] = e[i] = 0; break;
            case 1: a[i] = b[i] = c[i] = e[i] = ~0u; break;
            case 2: a[i] = 0x55555555u; b[i] = ~a[i]; break;
            case 3: a[i] = i & 1u ? 0u : ~0u; b[i] = a[i]; break;
            case 4: a[i] = 0; b[i] = ~0u; break;
            case 5: a[i] = ~0u; b[i] = 0; break;
            case 6: a[i] = 1u; b[i] = 0x80000000u; break;
            case 7: a[i] = (i % W == W-1) ? 0x80000000u : 0u;
                    b[i] = (i % W == 0) ? 1u : 0u; break;
            default: break;
            }
        }
        unsigned ca = 0, cand = 0, cor = 0, cx = 0, cs = 0;
        unsigned gc = 1, gx = 1, prev = a[0] & 1u, px = (a[0]^b[0]) & 1u;
        for (unsigned i = 0; i < N; ++i)
        {
            ca += NeonScalarCount(a[i]); cand += NeonScalarCount(a[i]&b[i]);
            cor += NeonScalarCount(a[i]|b[i]); cx += NeonScalarCount(a[i]^b[i]);
            cs += NeonScalarCount(a[i]&~b[i]);
            for (unsigned bit = 0; bit < 32; ++bit)
            {
                unsigned v = (a[i] >> bit) & 1u, x = ((a[i]^b[i]) >> bit) & 1u;
                gc += v != prev; gx += x != px; prev = v; px = x;
            }
        }
        NeonCheck(VECT_BITCOUNT(a.data(), a.data()+N) == ca, "count");
        NeonCheck(VECT_BITCOUNT_AND(a.data(), a.data()+N, b.data()) == cand, "count AND");
        NeonCheck(VECT_BITCOUNT_OR(a.data(), a.data()+N, b.data()) == cor, "count OR");
        NeonCheck(VECT_BITCOUNT_XOR(a.data(), a.data()+N, b.data()) == cx, "count XOR");
        NeonCheck(VECT_BITCOUNT_SUB(a.data(), a.data()+N, b.data()) == cs, "count SUB");
        NeonCheck(VECT_BLOCK_CHANGE(a.data(), N) == gc, "runs");
        unsigned r, bits;
        VECT_BLOCK_CHANGE_BC(a.data(), &r, &bits);
        NeonCheck(r == gc && bits == ca, "runs + count");
        VECT_BLOCK_XOR_CHANGE(a.data(), b.data(), N, &r, &bits);
        NeonCheck(r == gx && bits == cx, "XOR runs + count");
        VECT_BLOCK_XOR_CHANGE(b.data(), a.data(), N, &r, &bits);
        NeonCheck(r == gx && bits == cx, "XOR operand order");
        for (unsigned length : {1u, 3u, 4u, 15u, 16u, 17u, W})
        {
            unsigned count = 0, runs = 1, last = (a[1]^b[1]) & 1u;
            for (unsigned i = 1; i <= length; ++i)
                for (unsigned k = 0; k < 32; ++k)
                {
                    unsigned bit = ((a[i]^b[i]) >> k) & 1u;
                    count += bit; runs += bit != last; last = bit;
                }
            NeonCheck(VECT_BITCOUNT_XOR(a.data()+1, a.data()+1+length, b.data()+1) == count, "count tail");
            VECT_BLOCK_XOR_CHANGE(a.data()+1, b.data()+1, length, &r, &bits);
            NeonCheck(r == runs && bits == count, "change tail");
        }
        dst = ref = a;
        for (unsigned i = 0; i < N; ++i) ref[i] = (a[i] & b[i]);
        {
            auto result = VECT_AND_BLOCK(dst.data(), b.data());
            check_block("VECT_AND_BLOCK");
            unsigned acc = 0;
            for (unsigned i = 0; i < N; ++i) acc |= ref[i];
            NeonCheck(result == acc, "VECT_AND_BLOCK return");
        }
        dst = ref = a;
        for (unsigned i = 0; i < W; ++i) ref[i] = (a[i] & b[i]);
        {
            auto result = VECT_AND_DIGEST(dst.data(), b.data());
            check_block("VECT_AND_DIGEST");
            unsigned acc = 0;
            for (unsigned i = 0; i < W; ++i) acc |= ref[i];
            NeonCheck(result == (acc == 0), "VECT_AND_DIGEST return");
        }
        dst = ref = a;
        for (unsigned i = 0; i < W; ++i) ref[i] = (b[i] & c[i]);
        {
            auto result = VECT_AND_DIGEST_2WAY(dst.data(), b.data(), c.data());
            check_block("VECT_AND_DIGEST_2WAY");
            unsigned acc = 0;
            for (unsigned i = 0; i < W; ++i) acc |= ref[i];
            NeonCheck(result == (acc == 0), "VECT_AND_DIGEST_2WAY return");
        }
        dst = ref = a;
        for (unsigned i = 0; i < W; ++i) ref[i] = ((a[i] & b[i]) & c[i]);
        {
            auto result = VECT_AND_DIGEST_3WAY(dst.data(), b.data(), c.data());
            check_block("VECT_AND_DIGEST_3WAY");
            unsigned acc = 0;
            for (unsigned i = 0; i < W; ++i) acc |= ref[i];
            NeonCheck(result == (acc == 0), "VECT_AND_DIGEST_3WAY return");
        }
        dst = ref = a;
        for (unsigned i = 0; i < W; ++i) ref[i] = ((((a[i] & b[i]) & c[i]) & e[i]) & a[i]);
        {
            auto result = VECT_AND_DIGEST_5WAY(dst.data(), b.data(), c.data(), e.data(), a.data());
            check_block("VECT_AND_DIGEST_5WAY");
            unsigned acc = 0;
            for (unsigned i = 0; i < W; ++i) acc |= ref[i];
            NeonCheck(result == (acc == 0), "VECT_AND_DIGEST_5WAY return");
        }
        dst = ref = a;
        for (unsigned i = 0; i < N; ++i) ref[i] = (a[i] | b[i]);
        {
            auto result = VECT_OR_BLOCK(dst.data(), b.data());
            check_block("VECT_OR_BLOCK");
            unsigned acc = ~0u;
            for (unsigned i = 0; i < N; ++i) acc &= ref[i];
            NeonCheck(result == (acc == ~0u), "VECT_OR_BLOCK return");
        }
        dst = ref = a;
        for (unsigned i = 0; i < N; ++i) ref[i] = (b[i] | c[i]);
        {
            auto result = VECT_OR_BLOCK_2WAY(dst.data(), b.data(), c.data());
            check_block("VECT_OR_BLOCK_2WAY");
            unsigned acc = ~0u;
            for (unsigned i = 0; i < N; ++i) acc &= ref[i];
            NeonCheck(result == (acc == ~0u), "VECT_OR_BLOCK_2WAY return");
        }
        dst = ref = a;
        for (unsigned i = 0; i < N; ++i) ref[i] = ((a[i] | b[i]) | c[i]);
        {
            auto result = VECT_OR_BLOCK_3WAY(dst.data(), b.data(), c.data());
            check_block("VECT_OR_BLOCK_3WAY");
            unsigned acc = ~0u;
            for (unsigned i = 0; i < N; ++i) acc &= ref[i];
            NeonCheck(result == (acc == ~0u), "VECT_OR_BLOCK_3WAY return");
        }
        dst = ref = a;
        for (unsigned i = 0; i < N; ++i) ref[i] = ((((a[i] | b[i]) | c[i]) | e[i]) | a[i]);
        {
            auto result = VECT_OR_BLOCK_5WAY(dst.data(), b.data(), c.data(), e.data(), a.data());
            check_block("VECT_OR_BLOCK_5WAY");
            unsigned acc = ~0u;
            for (unsigned i = 0; i < N; ++i) acc &= ref[i];
            NeonCheck(result == (acc == ~0u), "VECT_OR_BLOCK_5WAY return");
        }
        dst = ref = a;
        for (unsigned i = 0; i < N; ++i) ref[i] = (a[i] ^ b[i]);
        {
            auto result = VECT_XOR_BLOCK(dst.data(), b.data());
            check_block("VECT_XOR_BLOCK");
            unsigned acc = 0;
            for (unsigned i = 0; i < N; ++i) acc |= ref[i];
            NeonCheck(result == acc, "VECT_XOR_BLOCK return");
        }
        dst = ref = a;
        for (unsigned i = 0; i < N; ++i) ref[i] = (b[i] ^ c[i]);
        {
            auto result = VECT_XOR_BLOCK_2WAY(dst.data(), b.data(), c.data());
            check_block("VECT_XOR_BLOCK_2WAY");
            unsigned acc = 0;
            for (unsigned i = 0; i < N; ++i) acc |= ref[i];
            NeonCheck(result == acc, "VECT_XOR_BLOCK_2WAY return");
        }
        dst = ref = a;
        for (unsigned i = 0; i < N; ++i) ref[i] = (a[i] &~ b[i]);
        {
            auto result = VECT_SUB_BLOCK(dst.data(), b.data());
            check_block("VECT_SUB_BLOCK");
            unsigned acc = 0;
            for (unsigned i = 0; i < N; ++i) acc |= ref[i];
            NeonCheck(result == acc, "VECT_SUB_BLOCK return");
        }
        dst = ref = a;
        for (unsigned i = 0; i < W; ++i) ref[i] = (a[i] &~ b[i]);
        {
            auto result = VECT_SUB_DIGEST(dst.data(), b.data());
            check_block("VECT_SUB_DIGEST");
            unsigned acc = 0;
            for (unsigned i = 0; i < W; ++i) acc |= ref[i];
            NeonCheck(result == (acc == 0), "VECT_SUB_DIGEST return");
        }
        dst = ref = a;
        for (unsigned i = 0; i < W; ++i) ref[i] = (b[i] &~ c[i]);
        {
            auto result = VECT_SUB_DIGEST_2WAY(dst.data(), b.data(), c.data());
            check_block("VECT_SUB_DIGEST_2WAY");
            unsigned acc = 0;
            for (unsigned i = 0; i < W; ++i) acc |= ref[i];
            NeonCheck(result == (acc == 0), "VECT_SUB_DIGEST_2WAY return");
        }
        dst = ref = a;
        for (unsigned i = 0; i < W; ++i) ref[i] = ((a[i] &~ b[i]) &~ c[i]);
        {
            auto result = VECT_SUB_DIGEST_3WAY(dst.data(), b.data(), c.data());
            check_block("VECT_SUB_DIGEST_3WAY");
            unsigned acc = 0;
            for (unsigned i = 0; i < W; ++i) acc |= ref[i];
            NeonCheck(result == (acc == 0), "VECT_SUB_DIGEST_3WAY return");
        }
        dst = ref = a;
        for (unsigned i = 0; i < W; ++i) ref[i] = ((((a[i] &~ b[i]) &~ c[i]) &~ e[i]) &~ a[i]);
        {
            auto result = VECT_SUB_DIGEST_5WAY(dst.data(), b.data(), c.data(), e.data(), a.data());
            check_block("VECT_SUB_DIGEST_5WAY");
            unsigned acc = 0;
            for (unsigned i = 0; i < W; ++i) acc |= ref[i];
            NeonCheck(result == (acc == 0), "VECT_SUB_DIGEST_5WAY return");
        }
        dst = ref = a;
        unsigned intersection = 0;
        for (unsigned i = 0; i < W; ++i) { ref[i] |= b[i]&c[i]; intersection |= b[i]&c[i]; }
        NeonCheck(VECT_AND_OR_DIGEST_2WAY(dst.data(), b.data(), c.data()) == (intersection == 0), "AND-OR return");
        check_block("AND-OR");
        for (unsigned i = 0; i < N; ++i) ref[i] = a[i] ^ 0xa5a5a5a5u;
        VECT_XOR_ARR_2_MASK(dst.data(), a.data(), a.data()+N, 0xa5a5a5a5u); check_block("mask XOR");
        for (unsigned i = 0; i < N; ++i) ref[i] = ~a[i] & 0xa5a5a5a5u;
        VECT_ANDNOT_ARR_2_MASK(dst.data(), a.data(), a.data()+N, 0xa5a5a5a5u); check_block("mask SUB");
        ref = a; VECT_COPY_BLOCK(dst.data(), a.data()); check_block("copy");
        VECT_STREAM_BLOCK(dst.data(), a.data()); check_block("stream");
        std::vector<unsigned char> bytes(N*4+1);
        std::memcpy(bytes.data()+1, a.data(), N*4);
        VECT_COPY_BLOCK_UNALIGN(dst.data(), bytes.data()+1); check_block("unaligned copy");
        VECT_STREAM_BLOCK_UNALIGN(dst.data(), bytes.data()+1); check_block("unaligned stream");
        for (unsigned i = 0; i < N; ++i) ref[i] = ~a[i];
        VECT_INVERT_BLOCK(dst.data()); check_block("invert");
        NeonCheck(VECT_IS_ZERO_BLOCK(a.data()) == (ca == 0), "zero");
        NeonCheck(VECT_IS_ONE_BLOCK(a.data()) == (ca == N*32), "ones");
        VECT_SET_BLOCK(dst.data(), 7u); ref.assign(N, 7u); check_block("set");
        VECT_BLOCK_SET_DIGEST(dst.data(), 0u); std::fill(ref.begin(), ref.begin()+W, 0u); check_block("set digest");
        NeonCheck(VECT_IS_DIGEST_ZERO(dst.data()), "digest zero");

        bm::id64_t digest = (bm::id64_t(rng()) << 32) | rng();
        if (trial == 0) digest = 0;
        if (trial == 1) digest = ~bm::id64_t(0);
        if (trial == 2) digest = bm::id64_t(1) << 63;
        unsigned dc = 0;
        for (unsigned i = 0; i < N; ++i)
        {
            bool selected = (digest >> (i/W)) & 1u;
            ref[i] = a[i] ^ (selected ? b[i] : 0u);
            if (selected) dc += NeonScalarCount(a[i]);
        }
        NeonCheck(VECT_BIT_COUNT_DIGEST(a.data(), digest) == dc, "digest count");
        VECT_BIT_BLOCK_XOR(dst.data(), a.data(), b.data(), digest); check_block("partial XOR");
        dst = a; VECT_BIT_BLOCK_XOR_2WAY(dst.data(), b.data(), digest); check_block("inplace partial XOR");
        for (unsigned carry_in = 0; carry_in < 2; ++carry_in)
        {
            for (unsigned direction = 0; direction < 2; ++direction)
            {
                unsigned carry = carry_in, acc = 0;
                for (unsigned k = 0; k < N; ++k)
                {
                    unsigned i = direction ? N-1-k : k;
                    unsigned v = a[i];
                    ref[i] = direction ? ((v >> 1) | (carry << 31)) : ((v << 1) | carry);
                    carry = direction ? v & 1u : v >> 31;
                    acc |= ref[i];
                }
                dst = a;
                unsigned nonzero = 99;
                bool co = direction ? VECT_SHIFT_L1(dst.data(), &nonzero, carry_in)
                                    : VECT_SHIFT_R1(dst.data(), &nonzero, carry_in);
                check_block("shift");
                NeonCheck(co == bool(carry) && bool(nonzero) == bool(acc), "shift carry/nonzero");
            }
            dst = a;
            // Force empty waves between nonempty waves.
            for (unsigned i = 0; i < N; ++i) if ((i/W)%3 == 1) dst[i] = 0;
            bm::id64_t d = 0, expected = 0;
            unsigned carry = carry_in;
            for (unsigned i = 0; i < N; ++i)
            {
                if (dst[i]) d |= bm::id64_t(1) << (i/W);
                ref[i] = ((dst[i] << 1) | carry) & b[i];
                carry = dst[i] >> 31;
                if (ref[i]) expected |= bm::id64_t(1) << (i/W);
            }
            bool co = VECT_SHIFT_R1_AND(dst.data(), carry_in, b.data(), &d);
            check_block("shift AND");
            NeonCheck(co == bool(carry) && d == expected, "shift AND carry/digest");
        }
        unsigned pos = 123, expected = 0;
        while (expected < N*32 && !((a[expected/32] >> (expected%32)) & 1u)) ++expected;
        bool found = VECT_BIT_FIND_FIRST(a.data(), 0, &pos);
        NeonCheck(found == (expected < N*32) && (!found || pos == expected), "find first");
        expected = 0;
        while (expected < N*32 && !(((a[expected/32]^b[expected/32]) >> (expected%32)) & 1u)) ++expected;
        found = VECT_BIT_FIND_DIFF(a.data(), b.data(), &pos);
        NeonCheck(found == (expected < N*32) && (!found || pos == expected), "find diff");
        for (unsigned i = 0; i < N; ++i) ref[i] = b[i] | a[i];
        dst = b;
        bool full = neon_decode_arr<true>(dst.data(), bytes.data()+1, N);
        check_block("decode OR unaligned");
        NeonCheck(full == std::all_of(ref.begin(), ref.end(), [](unsigned v) { return v == ~0u; }), "decode OR return");
        for (unsigned i = 0; i < N; ++i) ref[i] = b[i] & a[i];
        dst = b; neon_decode_arr<false>(dst.data(), bytes.data()+1, N); check_block("decode AND unaligned");
    }
    // Single-bit positions exercise every vector lane and word/digest/block edge.
    for (unsigned bit : {0u, 1u, 31u, 32u, 63u, 64u, 127u, 128u, W*32-1, W*32, N*32-1})
    {
        a.assign(N, 0u); a[bit/32] = 1u << (bit%32);
        unsigned pos;
        NeonCheck(VECT_BIT_FIND_FIRST(a.data(), bit/32, &pos) && pos == bit, "offset find");
        NeonCheck(VECT_BITCOUNT(a.data(), a.data()+N) == 1, "single-bit count");
    }
    // Exact-sized arrays keep tail reads visible to AddressSanitizer.
    for (unsigned size = 1; size <= 130; ++size)
    {
        std::vector<unsigned> arr(size);
        for (unsigned i = 0; i < size; ++i) arr[i] = 0x7fffffc0u+i;
        for (unsigned target : {0u, 0x7fffffffu, 0x80000000u, ~0u})
        {
            unsigned expected = unsigned(std::lower_bound(arr.begin(), arr.end(), target)-arr.begin());
            NeonCheck(VECT_LOWER_BOUND_SCAN_U32(arr.data(), target, 0, size-1) == expected, "unsigned lower bound");
        }
        std::vector<bm::gap_word_t> gap(size+1);
        gap[0] = bm::gap_word_t((size << 3) | (size & 1u));
        for (unsigned i = 1; i < size; ++i) gap[i] = bm::gap_word_t(i*500);
        gap[size] = 65535;
        for (unsigned target : {0u, 1u, 499u, 500u, 32767u, 32768u, 65000u, 65535u})
        {
            unsigned expected = unsigned(std::lower_bound(gap.begin()+1, gap.end(), target)-gap.begin());
            unsigned is_set;
            NeonCheck(VECT_GAP_BFIND(gap.data(), target, &is_set) == expected, "GAP lower bound");
            NeonCheck(is_set == ((size&1u)^((expected-1)&1u)), "GAP parity");
            NeonCheck(VECT_GAP_TEST(gap.data(), target) == is_set, "GAP test");
        }
        NeonCheck(bm::gap_bit_count_unr(gap.data()) == bm::gap_bit_count(gap.data()), "GAP sum dispatch");
    }
    for (unsigned trial = 0; trial < 32; ++trial)
    {
        a.assign(N, 0u);
        // Long runs remain within GAP header capacity.
        for (unsigned i = 0; i < N; ++i) a[i] = ((i+trial) / (trial+1)) & 1u ? ~0u : 0u;
        std::vector<bm::gap_word_t> expected(8192), actual(8192);
        unsigned len = 1, prev = a[0]&1u;
        for (unsigned bit = 1; bit < N*32; ++bit)
        {
            unsigned v = (a[bit/32] >> (bit%32)) & 1u;
            if (v != prev) expected[len++] = bm::gap_word_t(bit-1);
            prev = v;
        }
        expected[len] = 65535; expected[0] = bm::gap_word_t((len<<3)|(a[0]&1u));
        unsigned got = VECT_BIT_TO_GAP(actual.data(), a.data(), unsigned(actual.size()));
        NeonCheck(got == len && std::equal(expected.begin(), expected.begin()+len+1, actual.begin()), "bit to GAP");
    }
    std::vector<unsigned> indices = {65536, 65536, 65537, 65567, 65568, 65599, 65663, 131071, 131072};
    NeonCheck(VECT_ARR_BLOCK_LOOKUP(indices.data(), unsigned(indices.size()), 1u, 0u) == 8, "index lookup");
    dst.assign(N, 0u); ref = dst;
    VECT_SET_BLOCK_BITS(dst.data(), indices.data(), 0u, 8u);
    for (unsigned i = 0; i < 8; ++i) { unsigned n = indices[i]&65535u; ref[n/32] |= 1u << (n%32); }
    check_block("indexed scatter duplicate words");
    // More than one batch, maximal lane values and a short tail.
    std::vector<unsigned> large(N*40+3, ~0u);
    NeonCheck(VECT_BITCOUNT(large.data(), large.data()+large.size()) == large.size()*32,
              "wide count accumulation and batch tail");
    for (unsigned size : {1u, 3u, 4u, 7u, 16u, 31u})
    {
        std::vector<unsigned char> raw(size*4+1);
        std::vector<unsigned> source(size), target(size);
        for (unsigned i = 0; i < size; ++i) source[i] = unsigned(rng());
        std::memcpy(raw.data()+1, source.data(), size*4);
        bm::decoder or_decoder(raw.data()+1);
        or_decoder.get_32_OR(target.data(), size);
        NeonCheck(target == source, "decoder OR tail and offset");
        std::fill(target.begin(), target.end(), ~0u);
        bm::decoder and_decoder(raw.data()+1);
        and_decoder.get_32_AND(target.data(), size);
        NeonCheck(target == source, "decoder AND tail and offset");
    }
    // Exercise actual dispatch through public operations and serialization.
    for (unsigned trial = 0; trial < 12; ++trial)
    {
        std::vector<unsigned> wa(N*4), wb(wa.size()), wx(wa.size());
        unsigned expected_count = 0, expected_xor = 0;
        for (unsigned i = 0; i < wa.size(); ++i)
        {
            wa[i] = unsigned(rng()); wb[i] = unsigned(rng());
            if (trial & 1u) // mix dense, zero and interval-compressible blocks
            {
                if (i/N == 1) wa[i] = wb[i] = 0;
                if (i/N == 2) { wa[i] = (i/8)&1u ? ~0u : 0u; wb[i] = ~wa[i]; }
            }
            wx[i] = wa[i]^wb[i];
            expected_count += NeonScalarCount(wa[i]);
            expected_xor += NeonScalarCount(wx[i]);
        }
        bm::bvector<> av, bv, expected;
        bm::bit_import_u32(av, wa.data(), unsigned(wa.size()), true);
        bm::bit_import_u32(bv, wb.data(), unsigned(wb.size()), true);
        bm::bit_import_u32(expected, wx.data(), unsigned(wx.size()), false);
        NeonCheck(av.count() == expected_count, "public count");
        NeonCheck(bm::count_xor(av, bv) == expected_xor, "public XOR count");
        for (unsigned bit = 0; bit < wa.size()*32; ++bit)
            NeonCheck(av.test(bit) == bool((wa[bit/32] >> (bit%32))&1u), "import/test oracle");
        bm::bvector<> product;
        product.bit_xor(av, bv);
        NeonCheck(product == expected, "public XOR materialization");
        bm::serializer<bm::bvector<> > serializer;
        bm::serializer<bm::bvector<> >::buffer buffer;
        for (unsigned level = 0; level <= 6; ++level)
        {
            serializer.set_compression_level(level);
            serializer.serialize(av, buffer);
            bm::bvector<> restored;
            bm::deserialize(restored, buffer.buf());
            NeonCheck(restored == av, "serialization round trip");
            bm::operation_deserializer<bm::bvector<> > od;
            restored = bv; product.bit_and(av, bv);
            od.deserialize(restored, buffer.buf(), nullptr, bm::set_AND);
            NeonCheck(restored == product, "serialized AND");
            restored = bv; product.bit_or(av, bv);
            od.deserialize(restored, buffer.buf(), nullptr, bm::set_OR);
            NeonCheck(restored == product, "serialized OR");
        }
    }
    std::cout << "NEON kernels DONE OK" << std::endl;
}
#endif // BMNEONOPT
#endif
