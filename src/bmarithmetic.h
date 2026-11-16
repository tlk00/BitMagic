#ifndef BM_ARITHMETIC__H__INCLUDED__
#define BM_ARITHMETIC__H__INCLUDED__
/* Copyright(c) 2026 Anatoliy Kuznetsov. Apache License 2.0. */
/** \file bmarithmetic.h
    Internal zero-order adaptive arithmetic block codecs.
    Include after bm.h. Opcode 106 (bm::set_block_arithmetic) framing is
    supplied by the serializer.
    The body starts with a patched 19-bit field: bit 18 selects bytes (otherwise
    pooled low/high nibbles); low 18 bits count arithmetic bits, including a
    32-bit physical zero lookahead tail but excluding final word padding.
    Counts start at one and increase by one, independently for each block.
    Opcode 107 (bm::set_block_bit_0runs_arithmetic) stores the 0runs map before
    a variable-length copied-word stream. Its count excludes physical guard
    zeros; the decoder supplies bounded virtual zero lookahead. No rescaling,
    polarity trials, transform, model tables or early size-budget cutoff is used.

    The coder maintains an inclusive 32-bit interval [low, high], subdivided
    in proportion to symbol frequencies. Symbols come from successive low-to-high
    fields of each 32-bit word; all nibble positions share one model. After each
    symbol, encoder and decoder update identical counts and expand the interval
    by removing resolved leading bits or a middle-half underflow condition.
    Arithmetic code bits are MSB-first; the adapters reverse batched chunks for
    BM's LSB-first bit accumulator. The enclosing serializer selects candidates
    using complete record sizes; this header only encodes or consumes payloads.
*/
#include "encoding.h"

namespace bm
{
/** Arithmetic alphabets eligible for serializer candidate selection.
    @ingroup bvserial
*/
enum class arithmetic_coding_mode
{
    nibble, ///< Try the pooled 16-symbol nibble model.
    byte,   ///< Try the 256-symbol byte model.
    all     ///< Compare both complete records, nibble first on ties.
};


/** Internal fixed-alphabet cumulative model; no allocation or shared state.
    @tparam Alphabet Symbol count, either 16 (nibbles) or 256 (bytes).
    Fixed bounds specialize initialization, Fenwick updates and symbol search.
    @tparam DecoderModel Select decoder-specific optimizations when true.
    Encoders use fixed-alphabet Fenwick loops. Nibble decoders use four search
    probes and a forced-inline switch for updates on Clang, GCC and MSVC.
    Frequencies and proper-node update order match the runtime-alphabet model.
    total_ owns the root sum; the root slot is neither queried nor updated.
    @ingroup bvserial
*/
template<unsigned Alphabet, bool DecoderModel = false>
class arithmetic_model
{
public:
    /// Reset all counts to one for the fixed 16- or 256-symbol alphabet.
    void init() BMNOEXCEPT;
    /** Read a prefix from the initialized model.
        @param symbol Symbol index below Alphabet.
        @return Cumulative frequency strictly below symbol.
    */
    unsigned prefix(unsigned symbol) const BMNOEXCEPT;
    /** Search the initialized cumulative model.
        @param value Cumulative position below total().
        @param prefix Receives the cumulative frequency below the result.
        @return Symbol containing value.
    */
    unsigned find(unsigned value, unsigned& prefix) const BMNOEXCEPT;
    /** Increment an observed symbol and its affected cumulative nodes.
        @param symbol Symbol index below Alphabet in an initialized model.
        Increments the symbol frequency and total by one.
    */
    void update(unsigned symbol) BMNOEXCEPT;
    /// Return the positive frequency of an in-range symbol.
    unsigned frequency(unsigned symbol) const BMNOEXCEPT { return counts_[symbol]; }
    /// Return the sum of all frequencies.
    unsigned total() const BMNOEXCEPT { return total_; }
private:
    static_assert((Alphabet == 16) || (Alphabet == 256), "Unsupported arithmetic alphabet");
    unsigned counts_[Alphabet]; ///< Symbol counts, reset for each operation.
    unsigned tree_[Alphabet + 1]; ///< Fenwick counts; root slot initialized but unused.
    unsigned total_ = 0; ///< Current count sum.
};

/** Internal arithmetic producer using the existing packed bit writer.
    @ingroup bvserial
*/
struct arithmetic_block_encoder
{
    /** Write the codec body, patch its header and retain final padding.
        @param enc Borrowed writable encoder with at least 32 KiB available.
        @param block Borrowed complete logical BIT block.
        @param bytes True selects 256 byte symbols; false selects 16 nibbles.
        @return Meaningful arithmetic payload bits, including lookahead.
        The caller supplies opcode 106 (bm::set_block_arithmetic);
        this method owns no allocation.
    */
    template<class Encoder>
    static unsigned encode(Encoder& enc, const bm::word_t* block, bool bytes);
};

/** @internal Arithmetic producer for copied 0runs words.
    @ingroup bvserial
*/
struct arithmetic_words_encoder
{
    /** Encode copied words with a physical or virtual zero guard.
        @param enc Borrowed destination after the run map.
        @param words Borrowed contiguous copied words.
        @param word_count Number of copied 32-bit words.
        @param bytes Select byte symbols when true, pooled nibbles otherwise.
        @param physical_guard Append 32 physical lookahead bits when true.
        @return Declared physical arithmetic payload bits.
    */
    template<class Encoder>
    static unsigned encode(Encoder& enc, const bm::word_t* words,
                           unsigned word_count, bool bytes, bool physical_guard);
private:
    /** Encode using a compile-time symbol width and model alphabet.
        @tparam Bytes True for bytes, false for nibbles.
        @param enc Borrowed writable destination with at least 32 KiB available.
        @param block Borrowed word_count contiguous words.
        @param word_count Number of copied words.
        @param physical_guard Append the physical zero lookahead when true.
        @return Meaningful arithmetic payload bits, excluding word padding.
    */
    template<bool Bytes, class Encoder>
    static unsigned encode_impl(Encoder& enc, const bm::word_t* block,
                                unsigned word_count, bool physical_guard);
};

/** @internal Reverse a chunk between arithmetic MSB-first and bit-I/O LSB-first order.
    @param value Chunk stored in the low count bits.
    @param count Number of meaningful bits, from 1 through 31.
    @return Reversed count-bit value, used by both encoder and decoder batching.
    The portable permutation can compile to a native bit-reversal instruction.
*/
BMFORCEINLINE
unsigned arithmetic_reverse_chunk(unsigned value, unsigned count) BMNOEXCEPT
{
    BM_ASSERT(count && (count < 32));
    // Swap adjacent 1-, 2-, 4-, 8- and 16-bit groups, then retain only the
    // requested reversed prefix. This changes bit order within a chunk only.
    value = ((value >> 1) & 0x55555555u) | ((value & 0x55555555u) << 1);
    value = ((value >> 2) & 0x33333333u) | ((value & 0x33333333u) << 2);
    value = ((value >> 4) & 0x0F0F0F0Fu) | ((value & 0x0F0F0F0Fu) << 4);
    value = ((value >> 8) & 0x00FF00FFu) | ((value & 0x00FF00FFu) << 8);
    value = (value >> 16) | (value << 16);
    return value >> (32 - count);
}

/** Internal bounded bit consumer; borrows and preserves the enclosing accumulator.
    @ingroup bvserial
*/
template<class BitIn>
class arithmetic_bit_reader
{
public:
    /// Borrow bits and limit arithmetic consumption to count logical bits.
    arithmetic_bit_reader(BitIn& bits, unsigned count) : bits_(bits), left_(count) {}
    /// Read one arithmetic bit, rejecting exhaustion before touching the parent.
    unsigned get()
    {
        bm::serialization_require(left_ != 0);
        --left_;
        return bits_.get_bit();
    }
    /** Read a bounded chunk in arithmetic code order.
        @param count Number of bits to consume, from 1 through 31.
        @return Bits with the first physical bit in position count - 1.
        Rejects exhaustion before touching the borrowed parent accumulator.
    */
    BMFORCEINLINE
    unsigned get_msb(unsigned count)
    {
        BM_ASSERT(count && (count < 32));
        bm::serialization_require(left_ >= count);
        left_ -= count;
        return bm::arithmetic_reverse_chunk(bits_.get_bits(count), count);
    }
    /// Return unconsumed bits within the declared payload.
    unsigned remaining() const BMNOEXCEPT { return left_; }
private:
    BitIn& bits_; ///< Shared enclosing bit reader.
    unsigned left_; ///< Remaining logical payload budget.
};

/** @internal Bound physical arithmetic bits and synthesize zero lookahead.
    @ingroup bvserial
*/
template<class BitIn>
class arithmetic_virtual_reader
{
public:
    /** Borrow a bit accumulator and its physical payload extent.
        @param bits Borrowed accumulator, retained for following metadata.
        @param count Number of physically encoded arithmetic bits.
    */
    arithmetic_virtual_reader(BitIn& bits, unsigned count)
        : bits_(bits), left_(count) {}
    /** Return one physical bit, or a virtual zero after the boundary. */
    unsigned get()
    {
        if (left_)
        {
            --left_;
            return bits_.get_bit();
        }
        ++virtual_count_;
        return 0;
    }
    /** Read a chunk, supplying virtual zeros after the physical boundary.
        @param count Number of bits to consume, from 1 through 31.
        @return Bits with the first physical bit in position count - 1.
        Updates physical and virtual consumption without reading past the payload.
    */
    BMFORCEINLINE
    unsigned get_msb(unsigned count)
    {
        BM_ASSERT(count && (count < 32));
        const unsigned physical = left_ < count ? left_ : count;
        unsigned value = 0;
        if (physical)
        {
            value = bits_.get_bits(physical);
            left_ -= physical;
        }
        virtual_count_ += count - physical;
        // Unread high bits of value are zero, so reversal places virtual zeros
        // after the physical bits in arithmetic (MSB-first) order.
        return bm::arithmetic_reverse_chunk(value, count);
    }
    /// Return physical bits left unread.
    unsigned remaining() const BMNOEXCEPT { return left_; }
    /// Return zero bits supplied without touching the accumulator.
    unsigned virtual_count() const BMNOEXCEPT { return virtual_count_; }
private:
    BitIn& bits_; ///< Borrowed enclosing bit accumulator.
    unsigned left_; ///< Remaining physical payload bits.
    unsigned virtual_count_ = 0; ///< Synthesized zero lookahead bits.
};

/** @internal Adaptive arithmetic symbol producer for a borrowed bit accumulator.
    @ingroup bvserial
*/
template<class BitOut>
class arithmetic_symbol_encoder
{
public:
    /// Start a fresh interval at the current position of the borrowed output.
    explicit arithmetic_symbol_encoder(BitOut& output) : bits(output) { bits.flush_if_full(); }
    /// Encode one in-range symbol using the supplied initialized model, then update it.
    template<class Model> void put(Model& model, unsigned symbol);
    /// Terminate once and return the number of emitted payload bits including termination.
    unsigned finish();
private:
    /// Emit the resolved bit followed by all deferred underflow complements.
    void emit(unsigned bit);
    BitOut& bits; ///< Borrowed enclosing accumulator.
    unsigned low = 0; ///< Inclusive lower endpoint.
    unsigned high = ~0u; ///< Inclusive upper endpoint.
    unsigned pending = 0; ///< Deferred underflow bits.
    unsigned count = 0; ///< Emitted bits in this interval stream.
};

template<class BitOut>
void arithmetic_symbol_encoder<BitOut>::emit(unsigned bit)
{
    bits.put_bit(bit);
    ++count;
    while (pending)
    {
        bits.put_bit(bit ^ 1u);
        ++count;
        --pending;
    } // while pending
}

template<class BitOut>
unsigned arithmetic_symbol_encoder<BitOut>::finish()
{
    ++pending;
    emit(low < 0x40000000u ? 0 : 1);
    return count;
}

template<class BitOut>
template<class Model>
void arithmetic_symbol_encoder<BitOut>::put(Model& model, unsigned symbol)
{
    const unsigned left = model.prefix(symbol);
    const unsigned right = left + model.frequency(symbol);
    const bm::id64_t range = bm::id64_t(high) - low + 1;
    const unsigned base = low;
    const unsigned total = model.total();
    // Subdivide the inclusive interval using cumulative frequencies.
    // The 64-bit range represents the initial width of 2^32 exactly.
    high = base + unsigned(range * right / total - 1);
    low = base + unsigned(range * left / total);
    BM_ASSERT(low <= high);
    for (;;)
    {
        const unsigned difference = low ^ high;
        if (!(difference & 0x80000000u))
        {
            // Collapse the common MSB prefix into one output chunk.
            // Cap at 31 to keep all shifts and bit_out counts defined.
            const unsigned shift = difference ?
                31 - bm::bit_scan_reverse32(difference) : 31;
            const unsigned prefix = bm::arithmetic_reverse_chunk(
                low >> (32 - shift), shift);
            if (pending)
            {
                // Underflow complements follow the first resolved bit,
                // before the remaining common prefix bits.
                emit(prefix & 1u);
                if (shift > 1)
                {
                    bits.put_bits(prefix >> 1, shift - 1);
                    count += shift - 1;
                }
            }
            else
            {
                bits.put_bits(prefix, shift);
                count += shift;
            }
            // put_bits may retain a full word; emit/guard put_bit
            // calls require space. Flush that word without padding.
            bits.flush_if_full();
            low <<= shift;
            high = (high << shift) | (~0u >> (32 - shift));
            continue;
        }
        else if ((low >= 0x40000000u) && (high < 0xC0000000u))
        {
            // E3 underflow: the interval straddles the midpoint but
            // fits in the middle half. Defer output and expand it.
            ++pending;
            low -= 0x40000000u;
            high -= 0x40000000u;
        }
        else
        {
            break;
        }
        low <<= 1;
        high = (high << 1) | 1u;
    } // for renormalization
    model.update(symbol);
}

/** @internal Adaptive arithmetic symbol consumer over a bounded payload reader.
    @ingroup bvserial
*/
template<class Reader>
class arithmetic_symbol_decoder
{
public:
    /// Seed a fresh code register from the borrowed physical or virtual reader.
    explicit arithmetic_symbol_decoder(Reader& reader);
    /// Decode and return one symbol, advancing the interval and supplied model.
    template<class Model> unsigned get(Model& model);
private:
    Reader& input; ///< Borrowed payload reader with caller-selected lookahead.
    unsigned low = 0; ///< Inclusive lower endpoint.
    unsigned high = ~0u; ///< Inclusive upper endpoint.
    unsigned code = 0; ///< Current code register.
};

template<class Reader>
arithmetic_symbol_decoder<Reader>::arithmetic_symbol_decoder(Reader& reader) : input(reader)
{
    for (unsigned i = 0; i < 32; ++i)
    {
        code = (code << 1) | input.get();
    } // for i
}

template<class Reader>
template<class Model>
unsigned arithmetic_symbol_decoder<Reader>::get(Model& model)
{
    bm::serialization_require((code >= low) && (code <= high));
    const bm::id64_t range = bm::id64_t(high) - low + 1;
    const unsigned total = model.total();
    // Map the code point back to a cumulative-frequency position. The
    // +1/-1 adjustments invert the encoder's inclusive interval rounding.
    const unsigned value = unsigned(((bm::id64_t(code) - low + 1) * total - 1) / range);
    // The validated code interval and bounded 64-bit numerator imply
    // value < total. Keep input validation on the interval above.
    BM_ASSERT(value < total);
    unsigned left;
    const unsigned symbol = model.find(value, left);
    const unsigned right = left + model.frequency(symbol);
    const unsigned base = low;
    high = base + unsigned(range * right / total - 1);
    low = base + unsigned(range * left / total);
    for (;;)
    {
        const unsigned difference = low ^ high;
        if (!(difference & 0x80000000u))
        {
            // Collapse consecutive E1/E2 steps. Shifting unsigned values
            // discards the common prefix, including all E2 corrections.
            // Cap at 31 for defined shifts and bit_in's chunk interface.
            const unsigned shift = difference ?
                31 - bm::bit_scan_reverse32(difference) : 31;
            low <<= shift;
            high = (high << shift) | (~0u >> (32 - shift));
            code = (code << shift) | input.get_msb(shift);
            continue;
        }
        else if ((low >= 0x40000000u) && (high < 0xC0000000u))
        {
            // Mirror E3 underflow in both the interval and code register.
            low -= 0x40000000u;
            high -= 0x40000000u;
            code -= 0x40000000u;
        }
        else
        {
            break;
        }
        low <<= 1;
        high = (high << 1) | 1u;
        code = (code << 1) | input.get();
    } // for renormalization
    model.update(symbol);
    return symbol;
}

/** Internal arithmetic consumer supporting reconstruction and traversal-only use.
    @ingroup bvserial
*/
struct arithmetic_block_decoder
{
    /** Read one codec body and leave shared positioned at following metadata.
        @param dec Borrowed byte decoder, used to create a standalone accumulator.
        @param block Exclusive complete block output, or null for bounded skipping.
        @param shared Optional enclosing bit reader; its partial word is preserved.
        Validates the declared bit extent before reading, and never allows
        arithmetic lookahead to consume a following logical metadata bit.
    */
    template<class Decoder>
    static void decode(Decoder& dec, bm::word_t* block,
                       bm::bit_in<Decoder>* shared = nullptr);
private:
    /** Decode a block payload using a fixed symbol width and model alphabet.
        @tparam Bytes True for bytes, false for nibbles.
        @param bits Borrowed accumulator positioned after the codec header.
        @param count Validated physical payload extent in bits.
        @param block Borrowed aligned complete block, overwritten on success.
    */
    template<bool Bytes, class Decoder>
    static void decode_impl(bm::bit_in<Decoder>& bits, unsigned count,
                            bm::word_t* block);
};

/** @internal Decoder for opcode 107 (bm::set_block_bit_0runs_arithmetic).
    @ingroup bvserial
*/
struct arithmetic_0runs_decoder
{
    /** Restore copied words into their positions from a preceding run map.
        @param dec Borrowed byte decoder after opcode 107 (bm::set_block_bit_0runs_arithmetic).
        @param block Exclusive destination block, or null for traversal.
        @param positions Retained allocator-backed word-position workspace.
        @param shared Optional enclosing bit accumulator for shared XOR records.
    */
    template<class Decoder, class Alloc>
    static void decode(Decoder& dec, bm::word_t* block,
                       bm::heap_vector<unsigned, Alloc, true>& positions,
                       bm::bit_in<Decoder>* shared = nullptr);
private:
    /** Decode copied words using a fixed symbol width and model alphabet.
        @tparam Bytes True for bytes, false for nibbles.
        @param bits Borrowed accumulator positioned after the codec header.
        @param count Validated physical payload extent in bits.
        @param block Borrowed destination with zeroed uncopied words.
        @param copied Number of copied words to reconstruct.
        @param positions Borrowed retained scatter positions for copied words.
    */
    template<bool Bytes, class Decoder, class Alloc>
    static void decode_impl(bm::bit_in<Decoder>& bits, unsigned count,
        bm::word_t* block, unsigned copied,
        const bm::heap_vector<unsigned, Alloc, true>& positions);
};

template<unsigned Alphabet, bool DecoderModel>
inline
void arithmetic_model<Alphabet, DecoderModel>::init() BMNOEXCEPT
{
    total_ = Alphabet;
    tree_[0] = 0;
    // With unit initial frequencies, each Fenwick node equals the number of
    // symbols it spans: the value of its least significant set bit.
    for (unsigned i = 0; i < Alphabet; ++i)
    {
        counts_[i] = 1;
        tree_[i + 1] = (i + 1) & (0u - (i + 1));
    } // for i
}

template<unsigned Alphabet, bool DecoderModel>
inline
unsigned arithmetic_model<Alphabet, DecoderModel>::prefix(unsigned symbol) const BMNOEXCEPT
{
    BM_ASSERT(symbol < Alphabet);
    unsigned sum = 0;
    // Clear one set bit at a time to sum disjoint ranges preceding symbol.
    for (unsigned i = symbol; i; i &= i - 1)
    {
        sum += tree_[i];
    } // for i
    return sum;
}

template<unsigned Alphabet, bool DecoderModel>
inline
unsigned arithmetic_model<Alphabet, DecoderModel>::find(unsigned value, unsigned& sum) const BMNOEXCEPT
{
    BM_ASSERT(value < total_);
    // Find the largest prefix whose sum is <= value. Its length is the symbol
    // index; the returned sum is the lower boundary of that symbol's interval.
    if constexpr (DecoderModel && (Alphabet == 16))
    {
        // total_ owns the root sum. Search only the four proper-tree levels;
        // the probes at 8, 4, 2 and 1 stay within [1, 15].
        unsigned index = sum = 0;
        {
            const unsigned next = index + 8;
            if (sum + tree_[next] <= value)
            {
                index = next;
                sum += tree_[next];
            }
        }
        {
            const unsigned next = index + 4;
            if (sum + tree_[next] <= value)
            {
                index = next;
                sum += tree_[next];
            }
        }
        {
            const unsigned next = index + 2;
            if (sum + tree_[next] <= value)
            {
                index = next;
                sum += tree_[next];
            }
        }
        {
            const unsigned next = index + 1;
            if (sum + tree_[next] <= value)
            {
                index = next;
                sum += tree_[next];
            }
        }
        return index;
    }
    else
    {
        unsigned index = sum = 0;
        // The full-alphabet sum is total_, so a root probe cannot succeed.
        // Start at half the alphabet; every subsequent probe stays below it.
        for (unsigned step = Alphabet / 2; step; step >>= 1)
        {
            const unsigned next = index + step;
            if (sum + tree_[next] <= value)
            {
                index = next;
                sum += tree_[next];
            }
        } // for step
        BM_ASSERT(index < Alphabet);
        return index;
    }
}

template<unsigned Alphabet, bool DecoderModel>
inline
void arithmetic_model<Alphabet, DecoderModel>::update(unsigned symbol) BMNOEXCEPT
{
    BM_ASSERT(symbol < Alphabet);
    ++counts_[symbol];
    ++total_;
    // total_ owns the full-alphabet sum; prefix/search never read the root.
    for (unsigned i = symbol + 1; i < Alphabet; i += i & (0u - i))
    {
        ++tree_[i];
    } // for i
}

/** @internal Increment a decoder nibble's frequency, total and fixed Fenwick path.
    @param symbol Symbol index below 16 in an initialized model.
    BMFORCEINLINE keeps this switch inside the codec's symbol loop on
    supported compilers.
*/
template<>
BMFORCEINLINE
void arithmetic_model<16, true>::update(unsigned symbol) BMNOEXCEPT
{
    BM_ASSERT(symbol < 16);
    ++counts_[symbol];
    ++total_;
    // Only proper Fenwick ancestors need updates. total_ owns the root sum.
    // Keep the same node order as the generic update loop.
    switch (symbol)
    {
    case 0:
        ++tree_[1]; ++tree_[2]; ++tree_[4]; ++tree_[8];
        break;
    case 1:
        ++tree_[2]; ++tree_[4]; ++tree_[8];
        break;
    case 2:
        ++tree_[3]; ++tree_[4]; ++tree_[8];
        break;
    case 3:
        ++tree_[4]; ++tree_[8];
        break;
    case 4:
        ++tree_[5]; ++tree_[6]; ++tree_[8];
        break;
    case 5:
        ++tree_[6]; ++tree_[8];
        break;
    case 6:
        ++tree_[7]; ++tree_[8];
        break;
    case 7:
        ++tree_[8];
        break;
    case 8:
        ++tree_[9]; ++tree_[10]; ++tree_[12];
        break;
    case 9:
        ++tree_[10]; ++tree_[12];
        break;
    case 10:
        ++tree_[11]; ++tree_[12];
        break;
    case 11:
        ++tree_[12];
        break;
    case 12:
        ++tree_[13]; ++tree_[14];
        break;
    case 13:
        ++tree_[14];
        break;
    case 14:
        ++tree_[15];
        break;
    case 15:
        break;
    default:
        BM_ASSERT(false);
        break;
    } // switch symbol
}

template<class Encoder>
unsigned arithmetic_words_encoder::encode(Encoder& enc, const bm::word_t* block,
                                          unsigned word_count, bool bytes, bool physical_guard)
{
    if (bytes)
    {
        return bm::arithmetic_words_encoder::encode_impl<true>(enc, block, word_count, physical_guard);
    }
    return bm::arithmetic_words_encoder::encode_impl<false>(enc, block, word_count, physical_guard);
}

template<bool Bytes, class Encoder>
unsigned arithmetic_words_encoder::encode_impl(Encoder& enc,
    const bm::word_t* block, unsigned word_count, bool physical_guard)
{
    const auto header = enc.get_pos();
    unsigned count = 0;
    {
        bm::bit_out<Encoder> bits(enc);
        bits.put_bits(0, 19);
        bm::arithmetic_model<Bytes ? 256 : 16> model;
        model.init();
        bm::arithmetic_symbol_encoder<decltype(bits)> coder(bits);
        constexpr unsigned width = Bytes ? 8 : 4;
        const unsigned symbols = word_count * 32 / width;
        const unsigned mask = (1u << width) - 1;
        // Visit low-to-high byte/nibble fields, sharing one adaptive model
        // across all fields. No symbol table is transmitted.
        for (unsigned i = 0; i < symbols; ++i)
        {
            const unsigned pos = i * width;
            const unsigned symbol = (block[pos >> 5] >> (pos & 31)) & mask;
            coder.put(model, symbol);
        } // for i
        // Select a terminating interval prefix and resolve all deferred bits.
        count = coder.finish();
        // Whole-block records store decoder lookahead; copied-word records
        // omit it and let the virtual reader supply zeros at their boundary.
        for (unsigned i = 0; physical_guard && (i < 32); ++i)
        {
            bits.put_bit(0);
            ++count;
        } // for i
    }
    BM_ASSERT((count >= (physical_guard ? 34u : 2u)) && (count < (1u << 18)));
    const unsigned packed = count | (unsigned(Bytes) << 18);
    // The third byte also contains five payload bits. Patch only the 19
    // reserved header bits after the bit writer has flushed its last word.
    header[0] = static_cast<unsigned char>(packed);
    header[1] = static_cast<unsigned char>(packed >> 8);
    header[2] = static_cast<unsigned char>((header[2] & 0xF8u) | (packed >> 16));
    return count;
}

template<class Encoder>
unsigned arithmetic_block_encoder::encode(Encoder& enc, const bm::word_t* block, bool bytes)
{
    // Reuse the word-stream coder, with a physical guard for the full block.
    return bm::arithmetic_words_encoder::encode(enc, block, bm::set_block_size,
                                                 bytes, true);
}

template<class Decoder>
void arithmetic_block_decoder::decode(Decoder& dec, bm::word_t* block,
                                      bm::bit_in<Decoder>* shared)
{
    bm::bit_in<Decoder> local(dec);
    auto& bits = shared ? *shared : local;
    const unsigned header = bits.get_bits(19);
    unsigned count = header & 0x3FFFFu;
    bm::serialization_require((count >= 34) && (count <= 131072));
    // Traversal needs only the declared extent, without model updates or output.
    // Consume through the selected accumulator to preserve following metadata.
    if (!block)
    {
        while (count)
        {
            const unsigned n = count > 16 ? 16 : count;
            bits.get_bits(n);
            count -= n;
        } // while count
        return;
    }
    if (header >> 18)
    {
        bm::arithmetic_block_decoder::decode_impl<true>(bits, count, block);
    }
    else
    {
        bm::arithmetic_block_decoder::decode_impl<false>(bits, count, block);
    }
}

template<bool Bytes, class Decoder>
void arithmetic_block_decoder::decode_impl(bm::bit_in<Decoder>& bits, unsigned count,
    bm::word_t* block)
{
    bm::arithmetic_bit_reader<bm::bit_in<Decoder> > input(bits, count);
    bm::arithmetic_model<Bytes ? 256 : 16, true> model;
    model.init();
    bm::arithmetic_symbol_decoder<decltype(input)> coder(input);
    bm::bit_block_set(block, 0);
    constexpr unsigned width = Bytes ? 8 : 4;
    // Decode the known number of symbols; there is no end-of-stream symbol.
    for (unsigned i = 0; i < 65536 / width; ++i)
    {
        const unsigned symbol = coder.get(model);
        const unsigned pos = i * width;
        block[pos >> 5] |= bm::word_t(symbol) << (pos & 31);
    } // for i
    // Finalization emits two more bits than the normalization count, so the
    // 32-bit initialization/lookahead consumes all but two physical guard bits.
    bm::serialization_require(input.remaining() == 2);
    // Consume both guard bits only after validating the remaining extent.
    bm::serialization_require(input.get_msb(2) == 0);
}

template<class Decoder, class Alloc>
void arithmetic_0runs_decoder::decode(Decoder& dec, bm::word_t* block,
    bm::heap_vector<unsigned, Alloc, true>& positions,
    bm::bit_in<Decoder>* shared)
{
    const unsigned start = shared ? shared->get_bits(8) : dec.get_8();
    bm::serialization_require(start <= 1);
    unsigned kind = start, word = 0, copied = 0;
    if (block)
    {
        bm::bit_block_set(block, 0);
        positions.resize(0);
    }
    // The map alternates zero runs and literal runs until all block words are
    // covered. Record only literal destinations; zero runs need no symbols.
    while (word < bm::set_block_size)
    {
        const unsigned length = shared ? shared->get_bits(16) : dec.get_16();
        bm::serialization_require(length && (length <= bm::set_block_size - word));
        if (kind)
        {
            if (block)
            {
                positions.resize(copied + length);
                // Flatten literal runs into scatter positions for the packed
                // copied-word stream. Traversal avoids this workspace entirely.
                for (unsigned i = 0; i < length; ++i)
                {
                    positions[copied + i] = word + i;
                } // for i
            }
            copied += length;
        }
        word += length;
        kind ^= 1;
    } // while word
    bm::bit_in<Decoder> local(dec);
    auto& bits = shared ? *shared : local;
    const unsigned header = bits.get_bits(19);
    unsigned count = header & 0x3FFFFu;
    bm::serialization_require((count >= 2) && (count < (1u << 18)));
    // Skip exactly the physical stream; virtual lookahead has no stored bits.
    if (!block)
    {
        while (count)
        {
            const unsigned n = count > 16 ? 16 : count;
            bits.get_bits(n);
            count -= n;
        } // while count
        return;
    }
    if (header >> 18)
    {
        bm::arithmetic_0runs_decoder::decode_impl<true>(bits, count, block, copied, positions);
    }
    else
    {
        bm::arithmetic_0runs_decoder::decode_impl<false>(bits, count, block, copied, positions);
    }
}

template<bool Bytes, class Decoder, class Alloc>
void arithmetic_0runs_decoder::decode_impl(bm::bit_in<Decoder>& bits, unsigned count,
    bm::word_t* block, unsigned copied,
    const bm::heap_vector<unsigned, Alloc, true>& positions)
{
    bm::arithmetic_virtual_reader<bm::bit_in<Decoder> > input(bits, count);
    bm::arithmetic_model<Bytes ? 256 : 16, true> model;
    model.init();
    bm::arithmetic_symbol_decoder<decltype(input)> coder(input);
    constexpr unsigned width = Bytes ? 8 : 4;
    // Run lengths determine the exact symbol count, including an empty stream.
    for (unsigned i = 0; i < copied * 32 / width; ++i)
    {
        const unsigned symbol = coder.get(model);
        const unsigned pos = i * width;
        // Symbols are contiguous in the coded stream, but their destination
        // words follow the literal positions recovered from the run map.
        block[positions[pos >> 5]] |= bm::word_t(symbol) << (pos & 31);
    } // for i
    // Termination contributes two bits; the decoder's 32-bit lookahead must
    // therefore use exactly 30 virtual zeros and exhaust all physical bits.
    bm::serialization_require((input.remaining() == 0) && (input.virtual_count() == 30));
}

} // namespace bm
#endif
