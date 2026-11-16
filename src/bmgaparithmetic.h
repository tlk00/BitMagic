#ifndef BM_GAP_ARITHMETIC__H__INCLUDED__
#define BM_GAP_ARITHMETIC__H__INCLUDED__
/* Copyright(c) 2026 Anatoliy Kuznetsov. Apache License 2.0. */
/** @file bmgaparithmetic.h
    Internal magnitude GAP codec for opcode 110 (bm::set_block_gap_arithmetic).
    Includes after bm.h. The serializer supplies the opcode and selection policy.
    An 18-bit logical body extent precedes flags, gamma run count/corrections,
    an 18-bit arithmetic extent, magnitude symbols, literal remainders, and
    unchanged BIC set/clear and global XOR exception sections. Separate zero/one
    16-symbol models reset at every block. Only nonterminal run lengths are coded;
    the final run follows from the 65,536-bit total. No physical lookahead guard
    or model table is transmitted. Shared XOR metadata resumes in the same word.
*/
#include "bmarithmetic.h"
#include "bmbuffer.h"
#include "bmfunc.h"

namespace bm
{
/** @internal Tested original-chain gates and complete-record savings policy.
    @ingroup bvserial
*/
struct gap_arithmetic_policy
{
    /// Return the minimum original GAP words at levels 7 through 10; zero excludes the level.
    static unsigned gate(unsigned level) BMNOEXCEPT
        { return (level == 7 || level == 8) ? 606u : (level == 9 || level == 10) ? 234u : 0u; }
    /// Return the minimum complete-record byte saving; zero excludes the level.
    static unsigned minimum_gain(unsigned level) BMNOEXCEPT
        { return level == 7 ? 24u : (level >= 8 && level <= 10) ? 1u : 0u; }
    /** Admit a strictly smaller complete record after the original-chain gate.
        @param level Compression level, eligible at levels 7 through 10.
        @param words Original GAP words including header, before cleanup.
        @param incumbent Complete incumbent score, including shared metadata/padding.
        @param candidate Complete candidate score under the same layout.
        @return True only when both tested gate and inclusive saving cutoff pass.
    */
    static bool accepts(unsigned level, unsigned words, size_t incumbent, size_t candidate) BMNOEXCEPT
    {
        return bm::gap_arithmetic_policy::gate(level) && words >= bm::gap_arithmetic_policy::gate(level) &&
            candidate < incumbent && incumbent - candidate >= bm::gap_arithmetic_policy::minimum_gain(level);
    }
};

/** @internal Borrowed native-order bit section, excluding physical word padding.
    @ingroup bvserial
*/
struct gap_arithmetic_bits
{
    const unsigned char* data = nullptr; ///< Borrowed word-aligned section origin.
    unsigned offset = 0; ///< Logical start relative to data, in bits.
    unsigned count = 0; ///< Meaningful bits; zero denotes an absent section.
};

/** @internal Retained producer scratch; propagates the serializer allocator.
    @ingroup bvserial
*/
template<class Alloc>
struct gap_arithmetic_encoder_workspace
{
    bm::heap_vector<bm::gap_word_t, Alloc, true> base; ///< Restored cleaned endpoints.
    bm::heap_vector<bm::gap_word_t, Alloc, true> positions; ///< Native exception decode scratch.
    bm::heap_vector<unsigned char, Alloc, true> trial; ///< Complete trial record.
    bm::heap_vector<unsigned char, Alloc, true> best; ///< Smallest complete record.
    bm::gap_arithmetic_bits inner; ///< Unchanged set/clear chain.
    bm::gap_arithmetic_bits outer; ///< Unchanged global XOR list.
    unsigned correction0 = 0; ///< First-run minimum correction.
    unsigned correction1 = 0; ///< Alternate-run minimum correction.
    /// Reset per-record metadata, preserving lazy bounded buffer capacity.
    void init();
};

template<class Alloc>
void gap_arithmetic_encoder_workspace<Alloc>::init()
{
    inner = bm::gap_arithmetic_bits();
    outer = bm::gap_arithmetic_bits();
    correction0 = correction1 = 0;
    base.resize(8192);
    positions.resize(65536);
    // Each extent is limited to 2^18 bits; one extra word covers final padding.
    trial.resize(32776);
    best.resize(32776);
}

/** @internal Retained decoder scratch, independent of producer state.
    @ingroup bvserial
*/
template<class Alloc>
struct gap_arithmetic_decoder_workspace
{
    bm::heap_vector<bm::gap_word_t, Alloc, true> gap; ///< Decoded cleaned GAP endpoints.
    bm::heap_vector<unsigned char, Alloc, true> widths; ///< Nonterminal magnitude widths.
    bm::heap_vector<bm::gap_word_t, Alloc, true> positions; ///< Exception position scratch.
    /// Reset logical sizes for runs, retaining bounded capacity; runs is 1..8191.
    void init(unsigned runs);
};

template<class Alloc>
void gap_arithmetic_decoder_workspace<Alloc>::init(unsigned runs)
{
    BM_ASSERT(runs && runs < 8192);
    gap.resize(runs + 1);
    widths.resize(runs - 1);
}

/** @internal Two-context magnitude producer; exception bits are borrowed unchanged.
    @ingroup bvserial
*/
struct gap_arithmetic_encoder
{
    /** Encode a body frame after the caller's opcode byte.
        @param enc Writable word channel with at least 32776 spare bytes.
        @param gap Restored cleaned GAP, including header and terminal 65535.
        @param c0 First-run correction; zero means absent.
        @param c1 Alternate-run correction; zero means absent.
        @param inner Borrowed set/clear exception bits, or an empty view.
        @param outer Borrowed global XOR exception bits, or an empty view.
        @return Logical body bits, excluding its 18-bit extent and final padding.
        Caller retains source views until this call finishes. No allocation.
    */
    template<class Encoder>
    static unsigned encode(Encoder& enc, const bm::gap_word_t* gap,
        unsigned c0, unsigned c1, bm::gap_arithmetic_bits inner, bm::gap_arithmetic_bits outer);
private:
    /// Copy a borrowed meaningful bit section into the current accumulator.
    template<class BitOut>
    static void copy(BitOut& bits, bm::gap_arithmetic_bits section);
    /// Patch width bits in native-order words; storage must include the final word.
    static void patch(unsigned char* data, unsigned offset, unsigned width, unsigned value);
};

template<class BitOut>
void gap_arithmetic_encoder::copy(BitOut& bits, bm::gap_arithmetic_bits section)
{
    if (!section.count)
    {
        return;
    }
    bm::decoder source(section.data + (section.offset / 32) * 4);
    bm::bit_in<bm::decoder> input(source);
    if (section.offset % 32)
    {
        input.get_bits(section.offset % 32);
    }
    for (unsigned left = section.count; left;)
    {
        const unsigned count = left < 31 ? left : 31;
        bits.put_bits(input.get_bits(count), count);
        left -= count;
    } // for left
}

inline
void gap_arithmetic_encoder::patch(unsigned char* data, unsigned offset, unsigned width, unsigned value)
{
    for (unsigned i = 0; i < width; ++i)
    {
        const unsigned pos = offset + i;
        unsigned word;
        ::memcpy(&word, data + (pos / 32) * 4, sizeof(word));
        const unsigned mask = 1u << (pos % 32);
        word = (word & ~mask) | (((value >> i) & 1u) << (pos % 32));
        ::memcpy(data + (pos / 32) * 4, &word, sizeof(word));
    } // for i
}

template<class Encoder>
unsigned gap_arithmetic_encoder::encode(Encoder& enc, const bm::gap_word_t* gap,
    unsigned c0, unsigned c1, bm::gap_arithmetic_bits inner, bm::gap_arithmetic_bits outer)
{
    const unsigned runs = gap[0] >> 3;
    BM_ASSERT(runs && runs < 8192 && gap[runs] == 65535);
    unsigned char* frame = enc.get_pos();
    const size_t start_size = enc.size();
    unsigned body_count, arithmetic_field, count = 0;
    {
        bm::bit_out<Encoder> bits(enc);
        bits.put_bits(0, 18);
        bits.put_bits((gap[0] & 1u) | (unsigned(c0 != 0) << 3) | (unsigned(c1 != 0) << 4) |
                      (unsigned(inner.count != 0) << 5) | (unsigned(outer.count != 0) << 6), 8);
        bits.gamma(runs);
        if (c0)
        {
            bits.gamma(c0);
        }
        if (c1)
        {
            bits.gamma(c1);
        }
        arithmetic_field = unsigned((enc.size() - start_size) * 8) + bits.pending_bits();
        bits.put_bits(0, 18);
        bm::arithmetic_model<16> models[2];
        models[0].init();
        models[1].init();
        bm::arithmetic_symbol_encoder<decltype(bits)> coder(bits);
        for (unsigned i = 0; i + 1 < runs; ++i)
        {
            unsigned value = i ? unsigned(gap[i + 1] - gap[i] - 1) : gap[1];
            const unsigned correction = (i & 1u) ? c1 : c0;
            BM_ASSERT(value >= correction);
            value -= correction;
            const unsigned symbol = bm::bit_scan_reverse32(value + 1);
            BM_ASSERT(symbol < 16);
            auto& model = models[(gap[0] & 1u) ^ (i & 1u)];
            coder.put(model, symbol);
        } // for i
        count = coder.finish();
        for (unsigned i = 0; i + 1 < runs; ++i)
        {
            unsigned value = i ? unsigned(gap[i + 1] - gap[i] - 1) : gap[1];
            value -= (i & 1u) ? c1 : c0;
            const unsigned k = bm::bit_scan_reverse32(value + 1);
            if (k)
            {
                bits.put_bits(value + 1 - (1u << k), k);
            }
        } // for i
        bm::gap_arithmetic_encoder::copy(bits, inner);
        bm::gap_arithmetic_encoder::copy(bits, outer);
        body_count = unsigned((enc.size() - start_size) * 8) + bits.pending_bits() - 18;
    }
    bm::serialization_require(body_count < (1u << 18) && count < (1u << 18));
    bm::gap_arithmetic_encoder::patch(frame, 0, 18, body_count);
    bm::gap_arithmetic_encoder::patch(frame, arithmetic_field, 18, count);
    return body_count;
}

/** @internal Magnitude GAP consumer with length-only traversal.
    @ingroup bvserial
*/
struct gap_arithmetic_decoder
{
    /** Restore or skip one frame after opcode 110 (bm::set_block_gap_arithmetic).
        @param dec Borrowed byte channel.
        @param block Exclusive aligned BIT output; null skips without allocating models/scratch.
        @param workspace Allocator-owned scratch, reused across blocks.
        @param flags Scratch aligned BIT block for existing BIC window flags.
        @param shared Optional enclosing accumulator; following metadata is preserved.
        Validates logical extents, supported mode, run sums and exception ordering.
    */
    template<class Decoder, class Alloc>
    static void decode(Decoder& dec, bm::word_t* block,
        bm::gap_arithmetic_decoder_workspace<Alloc>& workspace, bm::word_t* flags,
        bm::bit_in<Decoder>* shared = nullptr);
};

template<class Decoder, class Alloc>
void gap_arithmetic_decoder::decode(Decoder& dec, bm::word_t* block,
    bm::gap_arithmetic_decoder_workspace<Alloc>& workspace, bm::word_t* flags,
    bm::bit_in<Decoder>* shared)
{
    bm::bit_in<Decoder> local(dec);
    auto& bits = shared ? *shared : local;
    const unsigned extent = bits.get_bits(18);
    bm::serialization_require(extent >= 29);
    const auto origin = dec.get_pos();
    const unsigned pending = bits.pending_bits();
    auto consumed = [&]()
    {
        return size_t(dec.get_pos() - origin) * 8 + pending - bits.pending_bits();
    };
    auto check_extent = [&]()
    {
        bm::serialization_require(consumed() <= extent);
    };
    if (!block)
    {
        for (unsigned left = extent; left;)
        {
            const unsigned count = left < 31 ? left : 31;
            bits.get_bits(count);
            left -= count;
        } // for left
    }
    else
    {
        const unsigned control = bits.get_bits(8);
        bm::serialization_require(!(control & 0x86u)); // Mode zero; bit 7 remains reserved.
        const unsigned runs = bits.gamma();
        bm::serialization_require(runs && runs < 8192);
        const unsigned c0 = (control & 8u) ? bits.gamma() : 0;
        const unsigned c1 = (control & 16u) ? bits.gamma() : 0;
        bm::serialization_require(c0 < 65535 && c1 < 65535);
        const unsigned count = bits.get_bits(18);
        check_extent();
        bm::serialization_require(count >= 2 && count <= extent - consumed());
        workspace.init(runs);
        bm::arithmetic_virtual_reader<bm::bit_in<Decoder> > input(bits, count);
        bm::arithmetic_model<16, true> models[2];
        models[0].init();
        models[1].init();
        bm::arithmetic_symbol_decoder<decltype(input)> coder(input);
        for (unsigned i = 0; i + 1 < runs; ++i)
        {
            auto& model = models[(control & 1u) ^ (i & 1u)];
            const unsigned symbol = coder.get(model);
            workspace.widths[i] = static_cast<unsigned char>(symbol);
        } // for i
        bm::serialization_require(!input.remaining() && input.virtual_count() == 30);
        auto* gap = workspace.gap.data();
        gap[0] = bm::gap_word_t((runs << 3) | (control & 1u));
        unsigned sum = 0;
        for (unsigned i = 0; i + 1 < runs; ++i)
        {
            const unsigned k = workspace.widths[i];
            bm::serialization_require(k < 16 && consumed() + k <= extent);
            const unsigned remainder = k ? bits.get_bits(k) : 0;
            sum += (1u << k) + remainder + ((i & 1u) ? c1 : c0);
            bm::serialization_require(sum < 65536);
            gap[i + 1] = bm::gap_word_t(sum - 1);
        } // for i
        gap[runs] = 65535;
        bm::gap_convert_to_bitset(block, gap);
        if (control & 0x60u)
        {
            workspace.positions.resize(65536);
        }
        if (control & 32u)
        {
            unsigned previous_role = 2;
            for (unsigned chain = 0; chain < 2; ++chain)
            {
                unsigned n = 0;
                const unsigned header = bits.decode_array(workspace.positions.data(), flags, &n);
                check_extent();
                bm::serialization_require(n && n < 65536);
                const unsigned role = bool(header & bm::h3f_ex_arr_1);
                bm::serialization_require(role != previous_role);
                for (unsigned i = 0; i < n; ++i)
                {
                    const unsigned pos = workspace.positions[i];
                    bm::serialization_require(!i || pos > workspace.positions[i - 1]);
                    if (role)
                    {
                        bm::set_bit(block, pos);
                    }
                    else
                    {
                        bm::clear_bit(block, pos);
                    }
                } // for i
                if (header & bm::h3f_ex_arr_ex_EOC)
                {
                    break;
                }
                bm::serialization_require(!chain);
                previous_role = role;
            } // for chain
        }
        if (control & 64u)
        {
            const unsigned n = bits.gamma();
            const unsigned minimum = bits.gamma();
            bm::serialization_require(n && n <= 4095 && minimum && minimum <= 65536);
            const unsigned correction = minimum - 1;
            bm::serialization_require(bm::id64_t(n - 1) * correction <= 65535);
            bits.bic_decode_u16(workspace.positions.data(), n, 0,
                               bm::gap_word_t(65535 - (n - 1) * correction));
            check_extent();
            unsigned previous = 0;
            for (unsigned i = 0; i < n; ++i)
            {
                const unsigned pos = workspace.positions[i] + i * correction;
                bm::serialization_require(pos < 65536 && (!i || pos > previous));
                block[pos >> 5] ^= bm::word_t(1) << (pos & 31);
                previous = pos;
            } // for i
        }
    }
    bm::serialization_require(consumed() == extent);
    if (!shared)
    {
        const unsigned padding = (32 - ((18 + extent) % 32)) % 32;
        bm::serialization_require(!padding || !bits.get_bits(padding));
    }
}
} // namespace bm
#endif
