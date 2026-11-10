#ifndef BM_BWT__H__INCLUDED__
#define BM_BWT__H__INCLUDED__
/*
Copyright(c) 2002-2026 Anatoliy Kuznetsov(anatoliy_kuznetsov at yahoo.com)

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

/*! \file bmbwt.h
    \brief Binary cyclic Burrows-Wheeler transforms for serialization (internal).

    Include after bm.h or bm64.h, as for the other serialization utilities.

    BM applies BWT to the alphabet {0, 1}, one logical bit per symbol. Logical
    position i is bit (i & 31) of word (i >> 5), least-significant bit first.
    Sorting packed bytes or numeric words would therefore give a different
    transform. The transformed last column is packed back into the same bit
    layout, preserving both the number of bits and the population count.

    These are cyclic transforms: no end-of-string sentinel is appended to the
    binary alphabet. The primary index identifies the sorted row representing
    the unrotated input and restores its starting position during inversion.
    Periodic inputs may have equal rotations; reconstruction does not require
    every sorted row to represent a distinct bit sequence.

    The whole-block path sorts 65,536 rotations by prefix doubling with linear
    scratch storage. The wave path independently transforms each of 64 waves
    of 1,024 bits using an explicit byte-per-bit rotation matrix. Its matrix is
    temporary sorting storage, not the packed output or a serialized payload.
    A whole block needs one 16-bit primary index; each wave needs 10 bits, so
    all 64 wave indexes occupy 80 bytes. The serializer supplies the framing.

    The direct BIC payload encoder/decoder implement the internal whole-block
    bm::set_block_bwt_unpacked format (opcode 96). Callers retain their scratch
    buffers and supply or consume opcode framing.

    BWT alone does not compress: it rearranges bits so a subsequent BM block
    codec may encode runs or positions more efficiently. Complete encoded
    record sizes, including primary-index framing, determine codec selection.
*/

#include <string.h>
#include "bmbuffer.h"
#include "encoding.h"

namespace bm
{

const unsigned bwt_wave_bits = 1024; ///< Binary symbols per independently transformed wave.
const unsigned bwt_matrix_bytes = bm::bwt_wave_bits * bm::bwt_wave_bits; ///< Byte-per-bit wave rotation matrix size.
const unsigned bwt_index_bytes = 80; ///< Packed size of 64 ten-bit wave primary indexes.

/**
    \brief Serializer controls for optional binary BWT trials.
    \ingroup bvserial
    \internal
*/
struct bwt_options
{
    unsigned sample = 1;       ///< Evaluate one in every N eligible BIT blocks; must be positive.
    bool whole_block = false; ///< Transform one 65,536-bit block instead of its independent waves.
    bool diagnostics = false; ///< Collect per-block measurements for whole-block trials.
};

/**
    \brief Complete-record and transformed-payload measurements for one block.
    \ingroup bvserial
    \internal

    Sizes are bytes. The candidate includes BWT framing and its primary index;
    payload excludes that outer framing. A smaller payload need not produce a
    smaller candidate. Gap counts measure BM bit-block changes, not a change
    in population count: binary BWT preserves the number of set bits.
*/
struct bwt_block_observation
{
    id64_t vector = 0;           ///< Input vector ordinal supplied by the serializer.
    id64_t block = 0;            ///< Physical bit-block index within the input vector.
    unsigned baseline = 0;      ///< Ordinary-codec record size before the BWT trial.
    unsigned payload = 0;       ///< Encoded transformed block without outer BWT framing.
    unsigned candidate = 0;     ///< Complete BWT candidate record size.
    unsigned selected = 0;      ///< Smaller of the ordinary and BWT candidate record sizes.
    unsigned popcount = 0;      ///< Number of set bits in either original or transformed block.
    unsigned gaps_before = 0;   ///< BM bit-block change count before transformation.
    unsigned gaps_after = 0;    ///< BM bit-block change count after transformation.
    unsigned primary = 0;       ///< Whole-block primary row index, in 0..65,535.
    unsigned baseline_codec = 0; ///< Leading opcode of the ordinary encoded record.
    unsigned payload_codec = 0; ///< Leading opcode of the encoded transformed block.
    bool baseline_bic = false;  ///< Ordinary record uses a recognized BIC block opcode.
    bool payload_bic = false;   ///< Transformed payload uses a recognized BIC block opcode.
};

/**
    \brief Aggregate binary BWT trial statistics and optional block observations.
    \ingroup bvserial
    \internal
    \tparam Alloc Application bvector allocator (BV::allocator_type).

    Serializers pass their application's allocator explicitly. The default
    supports standalone helper use. Observation storage is allocated lazily;
    callers may reserve a known block count before collecting diagnostics.
    Payload savings exclude outer BWT framing; accepted-record savings include
    it. Timing fields are nanoseconds and are collected only when requested.
*/
template<class Alloc = bm::standard_allocator>
struct bwt_statistics
{
    bm::heap_vector<bwt_block_observation, Alloc, true> blocks; ///< Owned per-block diagnostic observations.
    id64_t transform_ns = 0;             ///< Accumulated whole-block transform time.
    id64_t payload_ns = 0;               ///< Accumulated transformed-payload encoding time.
    id64_t eligible = 0;                 ///< Eligible BWT trial opportunities before sampling.
    id64_t evaluated = 0;                ///< Sampled opportunities actually transformed.
    id64_t accepted = 0;                 ///< Complete BWT records smaller than the ordinary records.
    id64_t baseline_bytes = 0;           ///< Sum of ordinary record sizes for evaluated trials.
    id64_t transformed_payload_bytes = 0; ///< Sum of transformed payload sizes without BWT framing.
    id64_t saved_bytes = 0;              ///< Complete-record bytes saved by accepted trials.
    id64_t gap_before = 0;               ///< Sum of original BM bit-block change counts.
    id64_t gap_after = 0;                ///< Sum of transformed BM bit-block change counts.
    id64_t gap_improved = 0;             ///< Trials reducing the bit-block change count.
    id64_t gap_equal = 0;                ///< Trials preserving the bit-block change count.
    id64_t gap_regressed = 0;            ///< Trials increasing the bit-block change count.
    id64_t gap_reduction = 0;            ///< Total change-count reduction in improving trials.
    id64_t gap_increase = 0;             ///< Total change-count increase in regressing trials.
    id64_t payload_improved = 0;         ///< Transformed payloads smaller than ordinary records.
    id64_t payload_saved_bytes = 0;      ///< Payload-only byte savings in improving trials.
    id64_t bic_baseline = 0;             ///< Evaluated ordinary records using recognized BIC opcodes.
    id64_t bic_transformed = 0;          ///< Transformed payloads using recognized BIC opcodes.
    id64_t bic_to_bic_accepted = 0;       ///< Accepted BWT trials with BIC ordinary and transformed records.
    id64_t bic_to_bic_saved_bytes = 0;    ///< Complete-record savings for accepted BIC-to-BIC trials.
};

/**
    \brief Reusable allocator-aware storage for whole-block binary BWT.
    \ingroup bvserial
    \internal
    \tparam Alloc Application bvector allocator (BV::allocator_type).

    Forward transformation uses prefix doubling: classes initially describe
    single bits, then prefixes of length 2, 4, and so on. Stable counting sort
    orders pairs of class numbers without constructing the rotation matrix.
    Time is O(n log n); five unsigned arrays require O(n) storage. At n=65,536,
    their logical size is 1.25 MiB rather than a 4 GiB byte-per-bit matrix.

    Inversion reuses rank for occurrence counts and order for packed restored
    bits. Inverse-only calls allocate these two arrays lazily (264 KiB for a
    whole block); shifted, next_rank, and counts remain unallocated. The direct
    BIC payload decoder additionally retains an 8 KiB last-column staging block.
    Storage is retained across calls, including calls with smaller n.
    Reserving the maximum n with prepare() before a sequence of transforms
    avoids allocations inside that sequence. One workspace must not be shared
    by overlapping transform calls or overlap either bit buffer.
*/
template<class Alloc = bm::standard_allocator>
struct bwt_full_workspace
{
    bm::heap_vector<unsigned, Alloc, true> order;     ///< Sorted rotation starts; packed output scratch during inversion.
    bm::heap_vector<unsigned, Alloc, true> shifted;   ///< Rotation starts shifted backward by the current prefix length.
    bm::heap_vector<unsigned, Alloc, true> rank;      ///< Prefix equivalence classes; occurrence ranks during inversion.
    bm::heap_vector<unsigned, Alloc, true> next_rank; ///< Equivalence classes for doubled prefixes.
    bm::heap_vector<unsigned, Alloc, true> counts;   ///< Class histogram, then counting-sort write positions.
    bm::heap_vector<bm::word_t, Alloc, true> last_column; ///< Lazy payload-decoder staging block; no stack-sized block.

    /**
        \brief Size the five forward-transform arrays, retaining their capacities.
        \param n Positive number of binary symbols; defaults to a full BM block.
        \note Elements need not be initialized here: each algorithm writes the
            scratch entries it reads. Growing beyond capacity may allocate.
    */
    void prepare(unsigned n = bm::set_block_size * 32);
};


/**
    \brief Transform N packed binary symbols using cyclic prefix-doubling BWT.
    \ingroup bvserial
    \internal
    \tparam N Power-of-two symbol count, at least 8; symbols are bits, not bytes.
    \tparam Alloc Application allocator used by the reusable workspace.
    \param src Read-only input, at least ceil(N / 32) words, LSB-first per word.
    \param dst Output, at least ceil(N / 32) words; must not overlap src or ws.
    \param ws Reusable scratch arrays, resized to N entries without shrinking capacity.
    \return Primary row whose rotation starts at input position zero.
    \note Output row r contains the bit preceding the start of sorted rotation r.
        No sentinel is stored. Equal rotations of periodic input may receive
        different primary indexes from another BWT implementation, yet have the
        same last column and reconstruct the same bit sequence. Unused high bits
        in the last output word are zero; high input bits outside N are ignored.
*/
template<unsigned N, class Alloc>
inline
unsigned bwt_transform_bits(const word_t* src, word_t* dst,
                            bwt_full_workspace<Alloc>& ws);

/**
    \brief Invert N packed binary symbols in place using binary LF mapping.
    \ingroup bvserial
    \internal
    \tparam N Power-of-two symbol count, at least 8.
    \tparam Alloc Application allocator used by the reusable workspace.
    \param block Input last column, replaced by the original ceil(N / 32) words.
    \param primary Primary row from the matching forward transform; less than N.
    \param ws Reusable scratch arrays; must not overlap block.
    \note The first sorted column consists of all zeros followed by all ones.
        LF(r) = C(bit[r]) + Occ(bit[r], r), where C(0)=0, C(1)=zero_count,
        and Occ counts equal symbols strictly before row r in the last column.
        Starting at primary, LF traversal restores the original bits backward.
        The two-symbol alphabet needs only zero/one counters, not byte histograms.
        Unused high bits in the last output word are zero.
*/
template<unsigned N, class Alloc>
inline
void bwt_inverse_bits(word_t* block, unsigned primary,
                      bwt_full_workspace<Alloc>& ws);

/**
    \brief Transform all 65,536 bits of one BM block as a single cyclic sequence.
    \ingroup bvserial
    \internal
    \tparam Alloc Application allocator used by the reusable workspace.
    \param src Read-only input of bm::set_block_size words.
    \param dst Output of bm::set_block_size words; must not overlap src or ws.
    \param ws Reusable whole-block scratch arrays.
    \return Primary row in 0..65,535, represented by an unsigned short.
*/
template<class Alloc>
inline
unsigned short bwt_transform_full(const word_t* src, word_t* dst,
                                  bwt_full_workspace<Alloc>& ws);

/**
    \brief Restore one whole BM block in place from its binary BWT last column.
    \ingroup bvserial
    \internal
    \tparam Alloc Application allocator used by the reusable workspace.
    \param block Input/output array of bm::set_block_size words.
    \param primary Primary row of the original unrotated 65,536-bit sequence.
    \param ws Reusable whole-block scratch arrays; must not overlap block.
*/
template<class Alloc>
inline
void bwt_inverse_full(word_t* block, unsigned short primary,
                      bwt_full_workspace<Alloc>& ws);

/**
    \brief Sort fixed-size wave rotation indexes without STL algorithms.
    \ingroup bvserial
    \internal
    \param order Input/output permutation of 0..bm::bwt_wave_bits-1, sorted by
        rotation contents and then by original start index.
    \param matrix Read-only byte-per-bit rotation matrix of bm::bwt_matrix_bytes
        bytes; row i starts at matrix + i*bm::bwt_wave_bits. Must not overlap order.
    \note Iterative merge sort performs O(n log n) rotation comparisons with
        n=bm::bwt_wave_bits. It uses one n-element unsigned-short scratch array
        (2 KiB for a 1,024-bit wave), with no heap allocation or recursion.
        The start-index tie break preserves the previous total rotation order,
        including primary indexes for periodic input with equal rotations.
*/
inline
void bwt_sort_wave_rotations(
                        unsigned short (&order)[bm::bwt_wave_bits],
                        const unsigned char* matrix);

/**
    \brief Apply independent binary BWT to every 1,024-bit wave of a BM block.
    \ingroup bvserial
    \internal
    \param src Read-only input of bm::set_block_size words, LSB-first per word.
    \param dst Packed last columns in wave order; must not overlap src or matrix.
    \param matrix Reusable scratch of at least bm::bwt_matrix_bytes bytes, disjoint
        from src, dst, and primary. Each byte temporarily holds one bit symbol.
    \param primary Output row indexes, one per wave, each in 0..1,023.
    \note Rotations are sorted by binary symbols with 0 before 1. Equal rotations
        are ordered by their input start index to make wave indexes deterministic.
        No rotation crosses a wave boundary. The same 1 MiB matrix is reused for
        each wave, and the caller owns its allocation policy. Sorting uses another
        2 KiB of bounded automatic scratch storage.
*/
inline
void bwt_transform_block(const word_t* src, word_t* dst,
                         unsigned char* matrix,
                         unsigned short (&primary)[bm::block_waves]);

/**
    \brief Restore all independently transformed waves of one BM block in place.
    \ingroup bvserial
    \internal
    \param block Input/output array of bm::set_block_size words containing wave last columns.
    \param primary One primary row per wave, each less than bm::bwt_wave_bits.
    \note Zero/one occurrence ranks are local to each wave. Repeating a shorter
        LF cycle is valid for periodic waves; exactly 1,024 symbols are restored.
        Scratch arrays are bounded automatic storage; no heap allocation occurs.
*/
inline
void bwt_inverse_block(word_t* block,
                       const unsigned short (&primary)[bm::block_waves]);

/**
    \brief Pack the 64 wave primary indexes into 80 bytes.
    \ingroup bvserial
    \internal
    \param primary Input indexes, each in 0..1,023 (ten significant bits).
    \param packed Output bytes, disjoint from primary. Bit b of primary[i] is
        written at packed bit 10*i+b, LSB-first within each byte. No host-endian
        word representation is serialized.
*/
inline
void bwt_pack_indexes(const unsigned short (&primary)[bm::block_waves],
                             unsigned char (&packed)[bm::bwt_index_bytes]);

/**
    \brief Unpack the 80-byte wave-index representation.
    \ingroup bvserial
    \internal
    \param packed Input bytes in the format written by bm::bwt_pack_indexes().
    \param primary Output indexes, disjoint from packed, each in 0..1,023.
*/
inline
void bwt_unpack_indexes(const unsigned char (&packed)[bm::bwt_index_bytes],
                        unsigned short (&primary)[bm::block_waves]);


/**
    \brief Encode the direct-position payload of a whole-block binary BWT record.
    \ingroup bvserial
    \internal
    \tparam Alloc Application bvector allocator (BV::allocator_type).
        Serializers and deserializers pass their application's allocator
        explicitly; the default supports internal codec tests.

    Format: bm::set_block_bwt_unpacked (opcode 96). The writer receives the
    already transformed 65,536-bit last column and its primary row. It encodes
    the ascending set-bit positions directly with center-minimal binary
    interpolative coding (BIC), without a nested BM block-codec record or a
    periodic-chain formula. The reader rebuilds that last column and applies
    binary inverse BWT to recover the original block.

    The retained payload layout contains a zero polarity bit, the 16-bit
    primary row, zero occupancy/map fields of nine and one bits, a gamma-coded
    position count plus one, the BIC positions, and gamma(1) for no chain formula.
    Zero positions represent an empty last column. Reserved fields remain in
    the stream to preserve its existing layout; older unreleased polarity and
    formula variants are rejected. The caller emits or consumes the opcode.
    The writer's position buffer and inverse-BWT scratch use Alloc. The reader
    decodes BIC directly into the packed last column without a position buffer.
*/
template<class Alloc = bm::standard_allocator>
struct bwt_bic_encoder
{
    /**
        \brief Write a direct BIC payload for an already transformed binary block.
        \ingroup bvserial
        \internal
        \tparam ENC Encoder supported by bm::bit_out (bit-packed output).
        \tparam PositionsVect Allocator-aware vector of bm::gap_word_t.
        \param enc Destination positioned immediately after the caller's
            bm::set_block_bwt_unpacked opcode (96). Must have sufficient capacity
            when using a bounded memory encoder; transformed must remain valid
            and unchanged throughout the call.
        \param transformed Read-only BWT last column, containing
            bm::set_block_size words in BM's LSB-first bit order.
        \param primary Row identifying the unrotated input in the sorted cyclic
            rotations, in 0..65,535, returned by bm::bwt_transform_full().
        \param positions Exclusive caller-owned scratch; resized to 65,536 entries.
            Contents are replaced and capacity retained across calls.
        \note Writes the payload only and flushes its bit accumulator. It does
            not compute the forward transform or write the opcode. The position
            buffer reserves space for every possible set bit before collection.
            The caller retains the scratch between blocks; no local owner is
            constructed by this helper.
    */
    template<class ENC, class PositionsVect>
    static void encode(ENC& enc, const bm::word_t* transformed, unsigned primary,
                       PositionsVect& positions);
};

/**
    \brief Restore the direct BIC payload of a whole-block binary BWT record.
    \ingroup bvserial
    \internal
    \tparam Alloc Application bvector allocator (BV::allocator_type).

    Consumes the payload of bm::set_block_bwt_unpacked (opcode 96), whose layout is
    documented by bm::bwt_bic_encoder. The caller consumes the opcode first.
    Reconstructs the last column directly from BIC positions, then applies
    binary inverse BWT. Output is restored by assignment; the enclosing
    deserializer applies logical operations such as OR to the destination.
    All inverse scratch belongs to the caller and retains its capacity.
*/
template<class Alloc = bm::standard_allocator>
struct bwt_bic_decoder
{
    /**
        \brief Consume a direct BIC payload and optionally overwrite the output block.
        \ingroup bvserial
        \internal
        \tparam DEC Decoder supported by bm::bit_in (bit-packed input).
        \param dec Source positioned immediately after an already consumed
            bm::set_block_bwt_unpacked opcode (96). Advances past the payload
            even when reconstruction is skipped.
        \param block SIMD-aligned output array of bm::set_block_size words,
            overwritten with the original pre-BWT bits. This is block restoration, not an OR
            operation: the enclosing deserializer merges this output with the
            destination bvector. A null or other non-address BM sentinel
            (IS_VALID_ADDR(block) is false) requests payload consumption and
            validation without reconstructing bits or running inverse BWT.
        \param workspace Exclusive caller-owned inverse-transform scratch. Its
            capacity is retained between calls; dry traversal does not allocate it.
        \note Checks the supported polarity, reserved fields, position-count
            bound, and absent chain formula through bm::serialization_require().
            Malformed/truncated input handling additionally depends on DEC.
            Rebuilding uses direct BIC-to-bitset decoding. Skip mode uses dry
            BIC decoding to advance the source without storing positions.
    */
    template<class DEC>
    static void decode(DEC& dec, bm::word_t* block,
                       bm::bwt_full_workspace<Alloc>& workspace);
};


// -----------------------------------------------------------------------
// Implementation
// -----------------------------------------------------------------------

template<class Alloc>
inline
void bwt_full_workspace<Alloc>::prepare(unsigned n)
{
    BM_ASSERT(n);
    order.resize(n);
    shifted.resize(n);
    rank.resize(n);
    next_rank.resize(n);
    counts.resize(n);
}


template<unsigned N, class Alloc>
inline
unsigned bwt_transform_bits(const word_t* src, word_t* dst,
                            bwt_full_workspace<Alloc>& ws)
{
    static_assert(N >= 8 && !(N & (N - 1)), "BWT size must be a power of two >= 8");
    BM_ASSERT(src && dst && src != dst);
    const unsigned n = N;
    const unsigned mask = n - 1;
    ws.prepare(n);

    // The initial alphabet has just two classes. Extract logical bits in BM
    // order, rather than comparing the numeric values of their storage words.
    unsigned zeros = 0;
    for (unsigned i = 0; i < n; ++i)
    {
        const unsigned bit = (src[i >> 5] >> (i & 31u)) & 1u;
        ws.rank[i] = bit;
        zeros += !bit;
    } // for i
    if (!zeros || zeros == n)
    {
        // Every rotation of a uniform sequence is identical. No sorting is
        // needed; row zero is a valid primary for either binary symbol.
        if (n < 32)
            dst[0] = src[0] & (~word_t(0) >> ((32 - (n & 31)) & 31));
        else
            ::memcpy(dst, src, n / 8);
        return 0;
    }
    unsigned z = 0;
    unsigned o = zeros;
    for (unsigned i = 0; i < n; ++i)
    {
        ws.order[ws.rank[i] ? o++ : z++] = i;
    } // for i
    BM_ASSERT(z == zeros && o == n);

    unsigned classes = 2;
    for (unsigned len = 1; len < n && classes < n; len <<= 1)
    {
        // order already sorts prefixes of length len. Shifting those starts
        // backward by len orders the second half of each doubled prefix.
        // Stable counting sort on the first-half class then orders the pair.
        BM_ASSERT(classes <= n);
        ::memset(ws.counts.data(), 0, size_t(classes) * sizeof(unsigned));
        for (unsigned i = 0; i < n; ++i)
        {
            const unsigned p = (ws.order[i] - len) & mask;
            ws.shifted[i] = p;
            BM_ASSERT(ws.rank[p] < classes);
            ++ws.counts[ws.rank[p]];
        } // for i
        unsigned sum = 0;
        for (unsigned c = 0; c < classes; ++c)
        {
            const unsigned count = ws.counts[c];
            ws.counts[c] = sum;
            sum += count;
        } // for c
        BM_ASSERT(sum == n);
        for (unsigned i = 0; i < n; ++i)
        {
            const unsigned p = ws.shifted[i];
            BM_ASSERT(ws.counts[ws.rank[p]] < n);
            ws.order[ws.counts[ws.rank[p]]++] = p;
        } // for i
        ws.next_rank[ws.order[0]] = 0;
        unsigned next_classes = 1;
        for (unsigned i = 1; i < n; ++i)
        {
            const unsigned a = ws.order[i-1];
            const unsigned b = ws.order[i];
            // The power-of-two mask implements cyclic wraparound. Class
            // numbers describe bit strings, not packed-word boundaries.
            next_classes += ws.rank[a] != ws.rank[b] ||
                ws.rank[(a + len) & mask] != ws.rank[(b + len) & mask];
            ws.next_rank[b] = next_classes - 1;
        } // for i
        ws.rank.swap(ws.next_rank);
        classes = next_classes;
    } // for len

    ::memset(dst, 0, ((n + 31) / 32) * sizeof(word_t));
    unsigned primary = 0;
    for (unsigned row = 0; row < n; ++row)
    {
        const unsigned p = ws.order[row];
        BM_ASSERT(p < n);
        if (!p)
            primary = row;
        // The last column contains the symbol immediately before each
        // rotation start, including input[n-1] for the unrotated row.
        const unsigned prev = (p - 1) & mask;
        if ((src[prev >> 5] >> (prev & 31u)) & 1u)
            dst[row >> 5] |= word_t(1) << (row & 31u);
    } // for row
    return primary;
}


template<unsigned N, class Alloc>
inline
void bwt_inverse_bits(word_t* block, unsigned primary,
                      bwt_full_workspace<Alloc>& ws)
{
    static_assert(N >= 8 && !(N & (N - 1)), "BWT size must be a power of two >= 8");
    const unsigned n = N;
    BM_ASSERT(block && primary < n);
    // Inversion needs only ranks and packed output. Lazily retain their peak
    // capacities; callers processing whole blocks may reserve 64K ranks upfront.
    ws.rank.resize(n);
    ws.order.resize((n + 31) / 32);
    unsigned zeros = 0;
    unsigned ones = 0;
    for (unsigned i = 0; i < n; ++i)
    {
        const unsigned bit = (block[i >> 5] >> (i & 31u)) & 1u;
        // Exclusive occurrence rank within the row's own binary symbol.
        ws.rank[i] = bit ? ones++ : zeros++;
    } // for i
    BM_ASSERT(zeros + ones == n);
    ::memset(ws.order.data(), 0, ((n + 31) / 32) * sizeof(unsigned));
    unsigned row = primary;
    for (unsigned pos = n; pos; )
    {
        --pos;
        BM_ASSERT(row < n);
        const unsigned bit = (block[row >> 5] >> (row & 31u)) & 1u;
        if (bit)
            ws.order[pos >> 5] |= 1u << (pos & 31u);
        // The first column is 0^zeros followed by 1^ones. LF matches a last-
        // column occurrence to the same occurrence in that sorted column.
        row = (bit ? zeros : 0u) + ws.rank[row];
    } // for pos
    if (n < 32)
        block[0] = ws.order[0];
    else
        ::memcpy(block, ws.order.data(), n / 8);
}


template<class Alloc>
inline
unsigned short bwt_transform_full(const word_t* src, word_t* dst,
                                  bwt_full_workspace<Alloc>& ws)
{
    static_assert(bm::set_block_size * 32 == 65536, "Whole-block BWT requires a 16-bit primary index");
    return static_cast<unsigned short>(
        bm::bwt_transform_bits<bm::set_block_size * 32>(src, dst, ws));
}


template<class Alloc>
inline
void bwt_inverse_full(word_t* block, unsigned short primary,
                      bwt_full_workspace<Alloc>& ws)
{
    bm::bwt_inverse_bits<bm::set_block_size * 32>(block, primary, ws);
}


inline
void bwt_sort_wave_rotations(
                        unsigned short (&order)[bm::bwt_wave_bits],
                        const unsigned char* matrix)
{
    static_assert(bm::bwt_wave_bits >= 2 && !(bm::bwt_wave_bits & (bm::bwt_wave_bits - 1)),
                  "Binary BWT wave width must be a power of two");
    BM_ASSERT(matrix);
    unsigned short scratch[bm::bwt_wave_bits];
    unsigned short* src = order;
    unsigned short* dst = scratch;
    for (unsigned width = 1; width < bm::bwt_wave_bits; width <<= 1)
    {
        // Every run has exactly width entries: the wave size is a power of
        // two. Merge adjacent runs and alternate the source/destination arrays.
        for (unsigned first = 0; first < bm::bwt_wave_bits; first += width * 2)
        {
            const unsigned middle = first + width;
            const unsigned end = middle + width;
            unsigned left = first;
            unsigned right = middle;
            unsigned out = first;
            while (left < middle && right < end)
            {
                const unsigned a = src[left];
                const unsigned b = src[right];
                BM_ASSERT(a < bm::bwt_wave_bits && b < bm::bwt_wave_bits);
                const int cmp = ::memcmp(matrix + a * bm::bwt_wave_bits,
                                         matrix + b * bm::bwt_wave_bits,
                                         bm::bwt_wave_bits);
                // Equal binary rotations retain ascending start indexes,
                // exactly as with the original wave comparator.
                if (cmp < 0 || (cmp == 0 && a <= b))
                    dst[out++] = src[left++];
                else
                    dst[out++] = src[right++];
            } // while left, right
            while (left < middle)
            {
                dst[out++] = src[left++];
            } // while left
            while (right < end)
            {
                dst[out++] = src[right++];
            } // while right
            BM_ASSERT(out == end);
        } // for first
        unsigned short* previous = src;
        src = dst;
        dst = previous;
    } // for width
    if (src != order)
        ::memcpy(order, src, sizeof(order));
}


inline
void bwt_transform_block(const word_t* src, word_t* dst,
                         unsigned char* matrix,
                         unsigned short (&primary)[bm::block_waves])
{
    BM_ASSERT(src && dst && matrix && src != dst);
    unsigned short order[bm::bwt_wave_bits];
    for (unsigned wave = 0; wave < bm::block_waves; ++wave)
    {
        const word_t* in = src + wave * bm::set_block_digest_wave_size;
        word_t* out = dst + wave * bm::set_block_digest_wave_size;
        unsigned char bits[bm::bwt_wave_bits];
        unsigned ones = 0;
        for (unsigned i = 0; i < bm::bwt_wave_bits; ++i)
        {
            // Expanding each bit to a byte makes memcmp compare binary
            // symbols in logical order, independent of packed bit layout.
            bits[i] = static_cast<unsigned char>((in[i >> 5] >> (i & 31u)) & 1u);
            ones += bits[i];
        } // for i
        if (ones == 0 || ones == bm::bwt_wave_bits)
        {
            ::memcpy(out, in, bm::bwt_wave_bits / 8);
            primary[wave] = 0;
            continue;
        }
        for (unsigned r = 0; r < bm::bwt_wave_bits; ++r)
        {
            unsigned char* row = matrix + r * bm::bwt_wave_bits;
            ::memcpy(row, bits + r, bm::bwt_wave_bits - r);
            ::memcpy(row + bm::bwt_wave_bits - r, bits, r);
            order[r] = static_cast<unsigned short>(r);
        } // for r
        // The explicit start-index tie break makes the order total even for
        // periodic waves, whose distinct starts can describe equal rotations.
        bm::bwt_sort_wave_rotations(order, matrix);
        for (unsigned i = 0; i < bm::set_block_digest_wave_size; ++i)
        {
            out[i] = 0;
        } // for i
        for (unsigned row = 0; row < bm::bwt_wave_bits; ++row)
        {
            const unsigned r = order[row];
            BM_ASSERT(r < bm::bwt_wave_bits);
            if (r == 0)
                primary[wave] = static_cast<unsigned short>(row);
            if (bits[(r + bm::bwt_wave_bits - 1) & (bm::bwt_wave_bits - 1)])
                out[row >> 5] |= word_t(1) << (row & 31u);
        } // for row
    } // for wave
}


inline
void bwt_inverse_block(word_t* block,
                       const unsigned short (&primary)[bm::block_waves])
{
    BM_ASSERT(block);
    for (unsigned wave = 0; wave < bm::block_waves; ++wave)
    {
        BM_ASSERT(primary[wave] < bm::bwt_wave_bits);
        word_t* out = block + wave * bm::set_block_digest_wave_size;
        unsigned char last[bm::bwt_wave_bits];
        unsigned short zeros_before[bm::bwt_wave_bits];
        unsigned short ones_before[bm::bwt_wave_bits];
        unsigned zeros = 0;
        unsigned ones = 0;
        for (unsigned i = 0; i < bm::bwt_wave_bits; ++i)
        {
            last[i] = static_cast<unsigned char>((out[i >> 5] >> (i & 31u)) & 1u);
            zeros_before[i] = static_cast<unsigned short>(zeros);
            ones_before[i] = static_cast<unsigned short>(ones);
            if (last[i])
                ++ones;
            else
                ++zeros;
        } // for i
        BM_ASSERT(zeros + ones == bm::bwt_wave_bits);
        word_t restored[bm::set_block_digest_wave_size] = {};
        unsigned row = primary[wave];
        for (unsigned pos = bm::bwt_wave_bits; pos; )
        {
            --pos;
            BM_ASSERT(row < bm::bwt_wave_bits);
            const unsigned bit = last[row];
            if (bit)
                restored[pos >> 5] |= word_t(1) << (pos & 31u);
            row = bit ? zeros + ones_before[row] : zeros_before[row];
        } // for pos
        ::memcpy(out, restored, sizeof(restored));
    } // for wave
}


inline
void bwt_pack_indexes(const unsigned short (&primary)[bm::block_waves],
                             unsigned char (&packed)[bm::bwt_index_bytes])
{
    static_assert(bm::bwt_index_bytes * 8 == bm::block_waves * 10, "Each binary wave index occupies ten bits");
    ::memset(packed, 0, sizeof(packed));
    for (unsigned i = 0; i < bm::block_waves; ++i)
    {
        BM_ASSERT(primary[i] < bm::bwt_wave_bits);
        for (unsigned b = 0; b < 10; ++b)
        {
            if (primary[i] & (1u << b))
            {
                const unsigned bit = i * 10 + b;
                packed[bit >> 3] |= static_cast<unsigned char>(1u << (bit & 7u));
            }
        } // for b
    } // for i
}


inline
void bwt_unpack_indexes(const unsigned char (&packed)[bm::bwt_index_bytes],
                        unsigned short (&primary)[bm::block_waves])
{
    for (unsigned i = 0; i < bm::block_waves; ++i)
    {
        unsigned p = 0;
        for (unsigned b = 0; b < 10; ++b)
        {
            const unsigned bit = i * 10 + b;
            p |= ((packed[bit >> 3] >> (bit & 7u)) & 1u) << b;
        } // for b
        primary[i] = static_cast<unsigned short>(p);
    } // for i
}

template<class Alloc>
template<class ENC, class PositionsVect>
inline
void bwt_bic_encoder<Alloc>::encode(ENC& enc,
                                    const bm::word_t* transformed,
                                    unsigned primary, PositionsVect& positions)
{
    BM_ASSERT(transformed && primary < bm::gap_max_bits);
    positions.resize(bm::gap_max_bits);
    unsigned count = 0;
    bm::bit_out<ENC> out(enc);
    out.put_bit(0); // original polarity
    out.put_bits(primary, 16);
    out.put_bits(0, 9); // retained reserved occupancy field
    out.put_bit(0);     // retained reserved map field
    for (unsigned i = 0; i < bm::set_block_size; ++i)
    {
        bm::word_t word = transformed[i];
        while (word)
        {
            positions[count++] = bm::gap_word_t(i * 32 + bm::bit_scan_forward32(word));
            word &= word - 1;
        } // while word
    } // for i
    out.gamma(count + 1);
    if (count)
        out.bic_encode_u16_cm(positions.data(), count, 0, 65535);
    out.gamma(1); // no chain formula
    out.flush();
}


template<class Alloc>
template<class DEC>
inline
void bwt_bic_decoder<Alloc>::decode(DEC& dec, bm::word_t* block,
                                   bm::bwt_full_workspace<Alloc>& workspace)
{
    bm::bit_in<DEC> in(dec);
    bm::serialization_require(in.get_bit() == 0); // original polarity only
    const unsigned primary = in.get_bits(16);
    const unsigned reserved = in.get_bits(9);
    const unsigned reserved_map = in.get_bit();
    bm::serialization_require(!reserved && !reserved_map);
    const unsigned active = 65536;
    const unsigned residual_count_plus = in.gamma();
    bm::serialization_require(residual_count_plus && residual_count_plus <= active + 1);
    const unsigned residual_count = residual_count_plus - 1;
    bm::word_t* packed = nullptr;
    if (IS_VALID_ADDR(block))
    {
        BM_ASSERT(bm::is_aligned(block));
        // Retain the staging block with the caller's workspace. Dry traversal
        // never sizes it; validation completes before modifying the output.
        workspace.last_column.resize(bm::set_block_size);
        packed = workspace.last_column.data();
        BM_ASSERT(block != packed);
        bm::bit_block_set(packed, 0);
        in.bic_decode_u16_bitset(packed, residual_count,
                                0, bm::gap_word_t(active - 1));
    }
    else
    {
        in.bic_decode_u16_dry(residual_count, 0, bm::gap_word_t(active - 1));
    }
    bm::serialization_require(in.gamma() == 1); // no chain formula
    if (!IS_VALID_ADDR(block))
        return;
    bm::bit_block_copy(block, packed);
    bm::bwt_inverse_full(block, bm::short_t(primary), workspace);
}

} // namespace bm
#endif
