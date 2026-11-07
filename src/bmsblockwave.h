#ifndef BMSBLOCKWAVE__H__INCLUDED__
#define BMSBLOCKWAVE__H__INCLUDED__

/** \file bmsblockwave.h
    \brief Internal fixed-width and adaptive wave-compacted superblock codecs.
    \ingroup bvserial
    \internal

    Implements bm::set_sblock_bienc_wave1024 (90). A superblock contains
    16,384 waves of 1,024 bits. Removing unoccupied waves produces a compact
    position domain for gap correction and BIC encoding. The map can describe
    occupied waves or their complement; this polarity concerns the wave map,
    not the values of the source bits. Explicit endpoint bounds are optional.

    Cost records and position-mapping helpers are also used by the adaptive
    wave codec, bm::set_sblock_bienc_wave_adaptive (97).
*/

#include "bmbuffer.h"
#include "encoding.h"

namespace bm
{

/** \brief Internal encoded-record bit accounting shared by wave codecs.
    \ingroup bvserial
    \internal
*/
struct sblock_wave_cost
{
    bm::id64_t bytes = 0;   ///< Complete encoded record size in bytes.
    bm::id64_t header = 0;  ///< Metadata bits, including opcode and explicit bounds.
    bm::id64_t map = 0;     ///< BIC wave-map payload bits.
    bm::id64_t bic = 0;     ///< BIC corrected-position payload bits.
    bm::id64_t padding = 0; ///< Terminal alignment bits.

    /// Return whether the component bits sum to the complete record size.
    bool reconciles() const;
};

/** \brief Internal eligibility and layout costs for one wave-compacted record.
    \ingroup bvserial
    \internal
*/
struct sblock_wave_record
{
    bool eligible = false; ///< At least one wave is unoccupied.
    unsigned best = 0;     ///< Index of the smallest complete encoded layout.
    bm::sblock_wave_cost costs[4]; ///< Occupied/derived, occupied/explicit,
                                  ///< empty/derived, empty/explicit layouts.
};

/** \brief Remove unoccupied waves from sorted superblock positions.
    \ingroup bvserial
    \internal
    \tparam OccupiedVect Sorted occupied-wave index container.
    \tparam ValuesVect Resizable unsigned-position scratch container.
    \param original Strictly increasing source positions, valid for n elements.
    \param n Number of source positions.
    \param width Number of bits per wave; must be nonzero.
    \param occupied Sorted map containing every wave referenced by original.
    \param values Receives n compacted positions; retains allocated capacity.
*/
template<class OccupiedVect, class ValuesVect>
void sblock_wave_compact(const unsigned* original, unsigned n, unsigned width,
                        const OccupiedVect& occupied, ValuesVect& values);

/** \brief Restore one compacted position to its original superblock position.
    \ingroup bvserial
    \internal
    \tparam OccupiedVect Sorted occupied-wave index container.
    \param packed Position in the compacted wave domain.
    \param width Number of bits per wave; must be nonzero.
    \param occupied Sorted occupied-wave map containing the packed wave rank.
    \return Original superblock-relative bit position.
*/
template<class OccupiedVect>
unsigned sblock_wave_restore_position(unsigned packed, unsigned width,
                                     const OccupiedVect& occupied);

/** \brief Internal allocator-aware fixed 1,024-bit wave codec and workspace.
    \ingroup bvserial
    \internal
    \tparam Alloc Application bvector allocator (BV::allocator_type).
        Serializers and decoders pass their application allocator explicitly.
    \note Encoding retains map and position capacity across preparations.
        Production decoding uses the overload accepting reusable map scratch.
*/
template<class Alloc = bm::standard_allocator>
struct sblock_wave1024_model
{
    enum
    {
        wave_size = 1024,  ///< Number of bits in each fixed wave.
        wave_count = 16384 ///< Number of waves in one superblock.
    };
    bm::heap_vector<unsigned, Alloc, true> occupied; ///< Occupied wave indexes.
    bm::heap_vector<unsigned, Alloc, true> empty;    ///< Unoccupied wave indexes.
    bm::heap_vector<unsigned, Alloc, true> values;   ///< Gap-corrected compacted positions.
    unsigned correction = 0; ///< Uniform gap reduction per preceding position.
    unsigned upper = 0;      ///< Inclusive upper bound in the corrected domain.

    /** \brief Reset all model state while retaining allocated scratch capacity.
        \note Called by prepare() before every encoding trial. The caller's
            layout record is reset separately by prepare().
    */
    void init() BMNOEXCEPT;

    /** \brief Build wave maps and corrected positions for encoding trials.
        \param original Strictly increasing superblock-relative positions.
        \param n Number of positions; at least two, at most 65,535 for this codec.
        \param r Receives reset layout costs and eligibility status.
        \return True if at least one wave is unoccupied and compaction is eligible.
    */
    bool prepare(const unsigned* original, unsigned n, bm::sblock_wave_record& r);

    /** \brief Encode one wave-map and endpoint-layout combination.
        \tparam BitOut Bit accumulator supporting gamma, BIC and fixed-width fields.
        \tparam Position Callable returning the current total bit position.
        \param out Destination accumulator, shared across all payload fields.
        \param position Counts bits from the record start, including its opcode.
        \param sb Absolute superblock index.
        \param mode Layout 0 through 3: bit 1 selects the empty-wave map;
            bit 0 selects explicit corrected-position endpoints.
        \param cost Receives header, map and position payload bit counts.
            The caller supplies final byte size and padding after flushing out.
        \note Requires a successful prepare(). The caller writes the opcode
            and owns the accumulator's final flush.
    */
    template<class BitOut, class Position>
    void encode(BitOut& out, Position position, unsigned sb, unsigned mode,
                bm::sblock_wave_cost& cost) const;

    /** \brief Decode using local allocator-aware wave-map scratch.
        \tparam BitIn Bit input reader.
        \param in Source positioned after the opcode.
        \param output Destination with capacity for at least 65,535 positions.
        \param sb Receives the absolute superblock index.
        \return Number of restored superblock-relative positions.
        \note This convenience overload allocates scratch per call. Repeated
            decoding should use the overload with reusable map containers.
    */
    template<class BitIn>
    static unsigned decode(BitIn& in, unsigned* output, unsigned& sb);

    /** \brief Decode using reusable wave maps.
        \tparam BitIn Bit input reader.
        \tparam MapVect Allocator-aware unsigned vector.
        \param in Source positioned after the opcode.
        \param output Destination with capacity for at least 65,535 positions.
        \param sb Receives the absolute superblock index.
        \param map Exclusive encoded-map scratch, capacity retained.
        \param occupied Exclusive occupied-wave scratch, capacity retained.
        \return Number of restored superblock-relative positions.
        \note Malformed payloads are rejected by bm::serialization_require().
    */
    template<class BitIn, class MapVect>
    static unsigned decode(BitIn& in, unsigned* output, unsigned& sb,
                           MapVect& map, MapVect& occupied);
};

// ----------------------------------------------------------------------

inline
bool sblock_wave_cost::reconciles() const
{
    return bytes * 8 == header + map + bic + padding;
}



template<class OccupiedVect, class ValuesVect>
inline
void sblock_wave_compact(const unsigned* original, unsigned n, unsigned width,
                        const OccupiedVect& occupied, ValuesVect& values)
{
    BM_ASSERT(width && (!n || (original && !occupied.empty())));
    values.resize(n);
    for (unsigned i = 0, rank = 0; i < n; ++i)
    {
        const unsigned wave = original[i] / width;
        while (rank < occupied.size() && occupied[rank] != wave)
        {
            ++rank;
        }
        BM_ASSERT(rank < occupied.size());
        values[i] = rank * width + original[i] % width;
    } // for i
}

template<class OccupiedVect>
inline
unsigned sblock_wave_restore_position(unsigned packed, unsigned width,
                                     const OccupiedVect& occupied)
{
    BM_ASSERT(width && packed / width < occupied.size());
    return occupied[packed / width] * width + packed % width;
}

template<class Alloc>
inline
void sblock_wave1024_model<Alloc>::init() BMNOEXCEPT
{
    occupied.clear();
    empty.clear();
    values.clear();
    correction = 0;
    upper = 0;
}

template<class Alloc>
inline
bool sblock_wave1024_model<Alloc>::prepare(const unsigned* original, unsigned n,
                                          bm::sblock_wave_record& r)
{
    init();
    r = bm::sblock_wave_record();
    values.resize(n);
    BM_ASSERT(n >= 2);
    occupied.reserve(bm::sblock_wave1024_model<Alloc>::wave_count);
    empty.reserve(bm::sblock_wave1024_model<Alloc>::wave_count);
    for (unsigned i = 0; i < n; ++i)
    {
        unsigned wave = original[i] / bm::sblock_wave1024_model<Alloc>::wave_size;
        if (occupied.empty() || occupied.back() != wave)
        {
            occupied.push_back(wave);
        }
    } // for i
    bm::sblock_wave_compact(original, n,
        bm::sblock_wave1024_model<Alloc>::wave_size, occupied, values);
    if (occupied.size() == bm::sblock_wave1024_model<Alloc>::wave_count)
    {
        return false;
    }
    for (unsigned wave = 0, j = 0; wave < bm::sblock_wave1024_model<Alloc>::wave_count; ++wave)
    {
        if (j < occupied.size() && occupied[j] == wave)
        {
            ++j;
        }
        else
        {
            empty.push_back(wave);
        }
    } // for wave
    r.eligible = true;
    correction = values[1] - values[0] - 1;
    for (unsigned i = 2; i < n; ++i)
    {
        if (values[i] - values[i - 1] - 1 < correction)
        {
            correction = values[i] - values[i - 1] - 1;
        }
    } // for i
    upper = unsigned(occupied.size()) * bm::sblock_wave1024_model<Alloc>::wave_size -
            1 - (n - 1) * correction;
    for (unsigned i = 0; i < n; ++i)
    {
        values[i] -= i * correction;
    } // for i
    BM_ASSERT(values.back() <= upper);
    return true;
}

template<class Alloc>
template<class BitOut, class Position>
inline
void sblock_wave1024_model<Alloc>::encode(BitOut& out, Position position,
                                        unsigned sb, unsigned mode,
                                        bm::sblock_wave_cost& cost) const
{
    const unsigned n = unsigned(values.size());
    const auto& map_values = (mode & 2) ? empty : occupied;
    out.put_bits(mode, 2); // polarity and explicit-boundary flag
    unsigned sb_mode = sb < 8 ? 0 : sb <= 255 ? 1 : sb <= 65535 ? 2 : 3;
    out.put_bits(sb_mode, 2);
    switch (sb_mode)
    {
    case 0:
        out.gamma(sb + 1);
        break;
    case 1:
        out.put_bits(sb, 8);
        break;
    case 2:
        out.put_16_no(sb);
        break;
    default:
        out.put_32_no(sb);
        break;
    } // switch sb_mode
    out.gamma(n);
    out.gamma(unsigned(map_values.size()));
    out.gamma(correction + 1);
    cost.header = position(); // includes placeholder opcode
    out.bic_encode_u32_cm(map_values.data(), unsigned(map_values.size()), 0,
                         bm::sblock_wave1024_model<Alloc>::wave_count - 1);
    cost.map = position() - cost.header;
    if (mode & 1)
    {
        const bm::id64_t before = position();
        out.gamma(values.front() + 1);
        out.gamma(upper - values.back() + 1);
        cost.header += position() - before;
    }
    const bm::id64_t before = position();
    if (mode & 1)
    {
        if (n > 2)
        {
            out.bic_encode_u32_cm(values.data() + 1, n - 2,
                                 values.front() + 1, values.back() - 1);
        }
    }
    else
    {
        out.bic_encode_u32_cm(values.data(), n, 0, upper);
    }
    cost.bic = position() - before;
}

template<class Alloc>
template<class BitIn>
inline
unsigned sblock_wave1024_model<Alloc>::decode(BitIn& in, unsigned* output,
                                            unsigned& sb)
{
    bm::heap_vector<unsigned, Alloc, true> map;
    bm::heap_vector<unsigned, Alloc, true> occupied;
    return bm::sblock_wave1024_model<Alloc>::decode(in, output, sb, map, occupied);
}

template<class Alloc>
template<class BitIn, class MapVect>
inline
unsigned sblock_wave1024_model<Alloc>::decode(BitIn& in, unsigned* output,
                                            unsigned& sb, MapVect& map,
                                            MapVect& occupied)
{
    map.clear();
    occupied.clear();
    const unsigned mode = in.get_bits(2);
    const unsigned sb_mode = in.get_bits(2);
    switch (sb_mode)
    {
    case 0:
        sb = in.gamma() - 1;
        break;
    case 1:
        sb = in.get_bits(8);
        break;
    case 2:
        sb = in.get_16_no();
        break;
    default:
        sb = in.get_32_no();
        break;
    } // switch sb_mode
    const unsigned n = in.gamma();
    const unsigned count = in.gamma();
    const unsigned correction_code = in.gamma();
    bm::serialization_require(n >= 2 && n <= 65535 && count >= 1 &&
        count < bm::sblock_wave1024_model<Alloc>::wave_count && correction_code >= 1);
    const unsigned d = correction_code - 1;
    map.resize(count);
    in.bic_decode_u32_cm(map.data(), count, 0, bm::sblock_wave1024_model<Alloc>::wave_count - 1);
    for (unsigned i = 0; i < count; ++i)
    {
        bm::serialization_require(map[i] < bm::sblock_wave1024_model<Alloc>::wave_count &&
                                  (!i || map[i] > map[i - 1]));
    } // for i
    if (mode & 2)
    {
        occupied.reserve(bm::sblock_wave1024_model<Alloc>::wave_count - count);
        for (unsigned wave = 0, j = 0; wave < bm::sblock_wave1024_model<Alloc>::wave_count; ++wave)
        {
            if (j < count && map[j] == wave)
            {
                ++j;
            }
            else
            {
                occupied.push_back(wave);
            }
        } // for wave
    }
    else
    {
        occupied.swap(map);
    }
    const unsigned k = unsigned(occupied.size());
    const bm::id64_t packed_limit = bm::id64_t(k) * bm::sblock_wave1024_model<Alloc>::wave_size;
    const bm::id64_t correction_total = bm::id64_t(n - 1) * d;
    bm::serialization_require(k && n <= packed_limit && correction_total < packed_limit);
    const unsigned upper = unsigned(packed_limit - 1 - correction_total);
    if (mode & 1)
    {
        const unsigned first_code = in.gamma();
        const unsigned last_code = in.gamma();
        bm::serialization_require(first_code && last_code && first_code <= upper + 1 &&
                                  last_code <= upper + 1);
        output[0] = first_code - 1;
        output[n - 1] = upper - (last_code - 1);
        bm::serialization_require(output[0] < output[n - 1] &&
                                  output[n - 1] - output[0] >= n - 1);
        if (n > 2)
        {
            in.bic_decode_u32_cm(output + 1, n - 2,
                                 output[0] + 1, output[n - 1] - 1);
        }
    }
    else
    {
        in.bic_decode_u32_cm(output, n, 0, upper);
    }
    unsigned previous = 0;
    unsigned prior_corrected = 0;
    for (unsigned i = 0; i < n; ++i)
    {
        bm::serialization_require(output[i] <= upper &&
                                  (!i || output[i] > prior_corrected));
        prior_corrected = output[i];
        const bm::id64_t packed = bm::id64_t(output[i]) + bm::id64_t(i) * d;
        bm::serialization_require(packed < packed_limit);
        const unsigned original = bm::sblock_wave_restore_position(
            unsigned(packed), bm::sblock_wave1024_model<Alloc>::wave_size, occupied);
        bm::serialization_require(!i || original > previous);
        output[i] = original;
        previous = original;
    } // for i
    return n;
}


// ----------------------------------------------------------------------
// Adaptive-width wave codec: bm::set_sblock_bienc_wave_adaptive (97).
// ----------------------------------------------------------------------

/** \brief Internal bit sink for exact wave-codec size trials without output.
    \ingroup bvserial
    \internal
*/
struct sblock_wave_bit_counter
{
    bm::id64_t bits = 0; ///< Accumulated payload bits, excluding final padding.

    /** Count a fixed-width field.
        \param value Ignored field value.
        \param n Number of bits in the field.
    */
    void put_bits(unsigned value, unsigned n);
    /// Count a 16-bit field; value is ignored.
    void put_16_no(unsigned value);
    /// Count a 32-bit field; value is ignored.
    void put_32_no(unsigned value);
    /** Count a gamma-coded positive integer.
        \param v Nonzero integer to encode.
    */
    void gamma(unsigned v);
    /** Count the same recursive BIC decisions as the production writer.
        \param p Strictly increasing unsigned values.
        \param n Number of values.
        \param lo Inclusive domain lower bound.
        \param hi Inclusive domain upper bound.
    */
    void bic_encode_u32_cm(const unsigned* p, unsigned n, unsigned lo, unsigned hi);
};

/** \brief Internal allocator-aware adaptive-width wave codec and workspace.
    \ingroup bvserial
    \internal
    \tparam Alloc Application bvector allocator (BV::allocator_type).
    \note build() starts a new input and caches occupied 128-bit waves.
        prepare() reuses that cache across width trials. All vectors retain
        capacity; serializers own and reuse this model across inputs.
*/
template<class Alloc = bm::standard_allocator>
struct sblock_wave_adaptive_model
{
    unsigned wave_size = 128;   ///< Selected wave width, a power of two from 128 to 4,096.
    unsigned wave_count = 131072; ///< Number of waves at the selected width.
    bm::heap_vector<unsigned, Alloc, true> smallest; ///< Cached occupied 128-bit waves.
    bm::heap_vector<unsigned, Alloc, true> occupied; ///< Occupied waves at the trial width.
    bm::heap_vector<unsigned, Alloc, true> empty;    ///< Unoccupied waves at the trial width.
    bm::heap_vector<unsigned, Alloc, true> values;   ///< Gap-corrected compacted positions.
    unsigned correction = 0; ///< Uniform gap reduction per preceding position.
    unsigned upper = 0;      ///< Inclusive corrected-position upper bound.

    /** Reset input and trial state without releasing vector capacity.
        \note Preserves wave_size and wave_count configuration. Called by build().
    */
    void init() BMNOEXCEPT;
    /** Reset maps, positions and scalar trial results, preserving input cache.
        \note Keeps smallest, wave_size and wave_count for the next width trial.
            Called by prepare() before rebuilding trial state.
    */
    void init_trial() BMNOEXCEPT;
    /** Cache the occupied 128-bit waves for a new source list.
        \param original Strictly increasing superblock-relative positions.
        \param n Number of source positions, at least two for encoding.
    */
    void build(const unsigned* original, unsigned n);
    /** Build a corrected position domain at the configured width.
        \param original Same source list supplied to the preceding build().
        \param n Number of source positions, from two through 65,535.
        \param r Receives reset layout costs and eligibility.
        \return True when at least one wave is unoccupied.
        \note Requires build() for the current source and a supported width.
    */
    bool prepare(const unsigned* original, unsigned n, bm::sblock_wave_record& r);
    /** Encode the common wave payload after the adaptive width selector.
        \tparam BitOut Bit accumulator supporting gamma, BIC and fixed-width fields.
        \tparam Position Callable reporting current total record bit position.
        \param out Destination bit accumulator, shared across all payload fields.
        \param position Counts bits from the record start, including its opcode.
        \param sb Absolute superblock index.
        \param mode Layout 0 through 3: bit 1 selects the empty-wave map;
            bit 0 selects explicit corrected-position endpoints.
        \param cost Receives header, map and position payload bit counts.
        \note Requires successful prepare(). The caller writes the opcode and
            3-bit width selector, flushes out, and fills final bytes and padding.
    */
    template<class BitOut, class Position>
    void encode(BitOut& out, Position position, unsigned sb, unsigned mode,
                bm::sblock_wave_cost& cost) const;
    /** Decode using local allocator-aware wave-map scratch.
        \tparam BitIn Bit input reader.
        \param in Source positioned after the opcode, before the width selector.
        \param output Destination with capacity for at least 65,535 positions.
        \param sb Receives the absolute superblock index.
        \return Number of restored superblock-relative positions.
        \note Repeated decoding should use the reusable-scratch overload.
    */
    template<class BitIn>
    static unsigned decode(BitIn& in, unsigned* output, unsigned& sb);
    /** Decode using retained, exclusive wave-map scratch.
        \tparam BitIn Bit input reader.
        \tparam MapVect Allocator-aware unsigned-position container.
        \param in Source positioned after the opcode, before the width selector.
        \param output Destination with capacity for at least 65,535 positions.
        \param sb Receives the absolute superblock index.
        \param map Encoded-map scratch whose capacity is retained.
        \param occupied Expanded occupied-map scratch whose capacity is retained.
        \return Number of restored superblock-relative positions.
        \note Rejects malformed payloads through bm::serialization_require().
    */
    template<class BitIn, class MapVect>
    static unsigned decode(BitIn& in, unsigned* output, unsigned& sb,
                           MapVect& map, MapVect& occupied);
};

/** \brief Internal winning adaptive width and layout, or the incumbent size.
    \ingroup bvserial
    \internal
*/
struct sblock_wave_adaptive_choice
{
    unsigned width = 0;    ///< Wave width code: wave size equals 128 shifted by width.
    unsigned mode = 0;     ///< Wave-map polarity and endpoint layout, 0 through 3.
    bm::id64_t bytes = 0;  ///< Complete record size, including opcode and padding.
};

/** \brief Find a strictly smaller adaptive record across six widths and four layouts.
    \ingroup bvserial
    \internal
    \tparam Alloc Application bvector allocator.
    \param model Retained model; build() resets it for the supplied input.
    \param original Strictly increasing superblock-relative source positions.
    \param n Number of positions, from two through 65,535.
    \param sb Absolute superblock index.
    \param incumbent Existing encoded-record byte count to beat.
    \return Best width/layout and byte count; ties preserve the incumbent.
    \note The caller prepares the winning width again before encoding: the
        workspace otherwise contains results from the final trial width.
*/
template<class Alloc>
bm::sblock_wave_adaptive_choice sblock_wave_adaptive_select(
    bm::sblock_wave_adaptive_model<Alloc>& model, const unsigned* original,
    unsigned n, unsigned sb, bm::id64_t incumbent);

// ----------------------------------------------------------------------

inline
void sblock_wave_bit_counter::put_bits(unsigned /* value */, unsigned n)
{
    bits += n;
}

inline
void sblock_wave_bit_counter::put_16_no(unsigned /* value */)
{
    bits += 16;
}

inline
void sblock_wave_bit_counter::put_32_no(unsigned /* value */)
{
    bits += 32;
}

inline
void sblock_wave_bit_counter::gamma(unsigned v)
{
    BM_ASSERT(v);
    bits += 2 * bm::bit_scan_reverse32(v) + 1;
}

inline
void sblock_wave_bit_counter::bic_encode_u32_cm(const unsigned* p, unsigned n,
                                              unsigned lo, unsigned hi)
{
    bm::bic_encode_u32_cm_sink(*this, p, n, lo, hi);
}

template<class Alloc>
inline
void sblock_wave_adaptive_model<Alloc>::init() BMNOEXCEPT
{
    init_trial();
    smallest.clear();
}

template<class Alloc>
inline
void sblock_wave_adaptive_model<Alloc>::init_trial() BMNOEXCEPT
{
    occupied.clear();
    empty.clear();
    values.clear();
    correction = 0;
    upper = 0;
}

template<class Alloc>
inline
void sblock_wave_adaptive_model<Alloc>::build(const unsigned* original, unsigned n)
{
    init();
    smallest.reserve(n);
    for (unsigned i = 0; i < n; ++i)
    {
        unsigned w = original[i] >> 7;
        if (smallest.empty() || smallest.back() != w)
        {
            smallest.push_back(w);
        }
    } // for i
}

template<class Alloc>
inline
bool sblock_wave_adaptive_model<Alloc>::prepare(const unsigned* original,
                                               unsigned n, bm::sblock_wave_record& r)
{
    init_trial();
    r = bm::sblock_wave_record();
    values.resize(n);
    BM_ASSERT(n >= 2);
    occupied.reserve(wave_count);
    empty.reserve(wave_count);
    unsigned shift = 0;
    while ((128u << shift) < wave_size)
    {
        ++shift;
    }
    for (unsigned w : smallest)
    {
        w >>= shift;
        if (occupied.empty() || occupied.back() != w)
        {
            occupied.push_back(w);
        }
    } // for w
    bm::sblock_wave_compact(original, n, wave_size, occupied, values);
    if (occupied.size() == wave_count)
    {
        return false;
    }
    for (unsigned wave = 0, j = 0; wave < wave_count; ++wave)
    {
        if (j < occupied.size() && occupied[j] == wave)
        {
            ++j;
        }
        else
        {
            empty.push_back(wave);
        }
    } // for wave
    r.eligible = true;
    correction = values[1] - values[0] - 1;
    for (unsigned i = 2; i < n; ++i)
    {
        if (values[i] - values[i - 1] - 1 < correction)
        {
            correction = values[i] - values[i - 1] - 1;
        }
    } // for i
    upper = unsigned(occupied.size()) * wave_size - 1 - (n - 1) * correction;
    for (unsigned i = 0; i < n; ++i)
    {
        values[i] -= i * correction;
    } // for i
    BM_ASSERT(values.back() <= upper);
    return true;
}


template<class Alloc>
template<class BitOut, class Position>
inline
void sblock_wave_adaptive_model<Alloc>::encode(BitOut& out, Position position,
                                             unsigned sb, unsigned mode,
                                             bm::sblock_wave_cost& cost) const
{
    const unsigned n = unsigned(values.size());
    const auto& map_values = (mode & 2) ? empty : occupied;
    out.put_bits(mode, 2); // polarity and explicit-boundary flag
    unsigned sb_mode = sb < 8 ? 0 : sb <= 255 ? 1 : sb <= 65535 ? 2 : 3;
    out.put_bits(sb_mode, 2);
    switch (sb_mode)
    {
    case 0:
        out.gamma(sb + 1);
        break;
    case 1:
        out.put_bits(sb, 8);
        break;
    case 2:
        out.put_16_no(sb);
        break;
    default:
        out.put_32_no(sb);
        break;
    } // switch sb_mode
    out.gamma(n);
    out.gamma(unsigned(map_values.size()));
    out.gamma(correction + 1);
    cost.header = position(); // includes placeholder opcode
    out.bic_encode_u32_cm(map_values.data(), unsigned(map_values.size()), 0,
                         wave_count - 1);
    cost.map = position() - cost.header;
    if (mode & 1)
    {
        const bm::id64_t before = position();
        out.gamma(values.front() + 1);
        out.gamma(upper - values.back() + 1);
        cost.header += position() - before;
    }
    const bm::id64_t before = position();
    if (mode & 1)
    {
        if (n > 2)
        {
            out.bic_encode_u32_cm(values.data() + 1, n - 2,
                                 values.front() + 1, values.back() - 1);
        }
    }
    else
    {
        out.bic_encode_u32_cm(values.data(), n, 0, upper);
    }
    cost.bic = position() - before;
}

template<class Alloc>
template<class BitIn>
inline
unsigned sblock_wave_adaptive_model<Alloc>::decode(BitIn& in, unsigned* output,
                                                 unsigned& sb)
{
    bm::heap_vector<unsigned, Alloc, true> map;
    bm::heap_vector<unsigned, Alloc, true> occupied;
    return bm::sblock_wave_adaptive_model<Alloc>::decode(in, output, sb, map, occupied);
}

template<class Alloc>
template<class BitIn, class MapVect>
inline
unsigned sblock_wave_adaptive_model<Alloc>::decode(BitIn& in, unsigned* output,
                                                 unsigned& sb, MapVect& map,
                                                 MapVect& occupied)
{
    map.clear();
    occupied.clear();
    const unsigned width_code = in.get_bits(3);
    bm::serialization_require(width_code < 6);
    const unsigned wave_size = 128u << width_code;
    const unsigned wave_count = (1u << 24) / wave_size;
    const unsigned mode = in.get_bits(2);
    const unsigned sb_mode = in.get_bits(2);
    switch (sb_mode)
    {
    case 0:
        sb = in.gamma() - 1;
        break;
    case 1:
        sb = in.get_bits(8);
        break;
    case 2:
        sb = in.get_16_no();
        break;
    default:
        sb = in.get_32_no();
        break;
    } // switch sb_mode
    const unsigned n = in.gamma();
    const unsigned count = in.gamma();
    const unsigned correction_code = in.gamma();
    bm::serialization_require(n >= 2 && n <= 65535 && count >= 1 &&
        count < wave_count && correction_code >= 1);
    const unsigned d = correction_code - 1;
    map.resize(count);
    in.bic_decode_u32_cm(map.data(), count, 0, wave_count - 1);
    for (unsigned i = 0; i < count; ++i)
    {
        bm::serialization_require(map[i] < wave_count &&
                                  (!i || map[i] > map[i - 1]));
    } // for i
    if (mode & 2)
    {
        occupied.reserve(wave_count - count);
        for (unsigned wave = 0, j = 0; wave < wave_count; ++wave)
        {
            if (j < count && map[j] == wave)
            {
                ++j;
            }
            else
            {
                occupied.push_back(wave);
            }
        } // for wave
    }
    else
    {
        occupied.swap(map);
    }
    const unsigned k = unsigned(occupied.size());
    const bm::id64_t packed_limit = bm::id64_t(k) * wave_size;
    const bm::id64_t correction_total = bm::id64_t(n - 1) * d;
    bm::serialization_require(k && n <= packed_limit && correction_total < packed_limit);
    const unsigned upper = unsigned(packed_limit - 1 - correction_total);
    if (mode & 1)
    {
        const unsigned first_code = in.gamma();
        const unsigned last_code = in.gamma();
        bm::serialization_require(first_code && last_code && first_code <= upper + 1 &&
                                  last_code <= upper + 1);
        output[0] = first_code - 1;
        output[n - 1] = upper - (last_code - 1);
        bm::serialization_require(output[0] < output[n - 1] &&
                                  output[n - 1] - output[0] >= n - 1);
        if (n > 2)
        {
            in.bic_decode_u32_cm(output + 1, n - 2,
                                 output[0] + 1, output[n - 1] - 1);
        }
    }
    else
    {
        in.bic_decode_u32_cm(output, n, 0, upper);
    }
    unsigned previous = 0;
    unsigned prior_corrected = 0;
    for (unsigned i = 0; i < n; ++i)
    {
        bm::serialization_require(output[i] <= upper &&
                                  (!i || output[i] > prior_corrected));
        prior_corrected = output[i];
        const bm::id64_t packed = bm::id64_t(output[i]) + bm::id64_t(i) * d;
        bm::serialization_require(packed < packed_limit);
        const unsigned original = bm::sblock_wave_restore_position(
            unsigned(packed), wave_size, occupied);
        bm::serialization_require(!i || original > previous);
        output[i] = original;
        previous = original;
    } // for i
    return n;
}


template<class Alloc>
inline
bm::sblock_wave_adaptive_choice sblock_wave_adaptive_select(
    bm::sblock_wave_adaptive_model<Alloc>& model, const unsigned* original,
    unsigned n, unsigned sb, bm::id64_t incumbent)
{
    model.build(original, n);
    bm::sblock_wave_adaptive_choice best;
    best.bytes = incumbent;
    for (unsigned w = 0; w < 6; ++w)
    {
        model.wave_size = 128u << w;
        model.wave_count = (1u << 24) / model.wave_size;
        bm::sblock_wave_record record;
        if (!model.prepare(original, n, record))
        {
            continue;
        }
        // Count each map and payload once, then combine the four layouts.
        bm::id64_t maps[2];
        bm::id64_t payload[2];
        bm::id64_t bounds[2] = {0, 0};
        for (unsigned polarity = 0; polarity < 2; ++polarity)
        {
            const auto& map = polarity ? model.empty : model.occupied;
            bm::sblock_wave_bit_counter c;
            c.bic_encode_u32_cm(map.data(), unsigned(map.size()), 0, model.wave_count - 1);
            maps[polarity] = c.bits;
        } // for polarity
        bm::sblock_wave_bit_counter c;
        c.bic_encode_u32_cm(model.values.data(), n, 0, model.upper);
        payload[0] = c.bits;
        c.bits = 0;
        c.gamma(model.values.front() + 1);
        c.gamma(model.upper - model.values.back() + 1);
        bounds[1] = c.bits;
        c.bits = 0;
        if (n > 2)
        {
            c.bic_encode_u32_cm(model.values.data() + 1, n - 2,
                               model.values.front() + 1, model.values.back() - 1);
        }
        payload[1] = c.bits;
        for (unsigned mode = 0; mode < 4; ++mode)
        {
            bm::sblock_wave_bit_counter h;
            h.put_bits(0, 3 + 4); // width, polarity/bounds and superblock-index modes
            unsigned sm = sb < 8 ? 0 : sb <= 255 ? 1 : sb <= 65535 ? 2 : 3;
            if (!sm)
            {
                h.gamma(sb + 1);
            }
            else
            {
                h.put_bits(sb, sm == 1 ? 8 : sm == 2 ? 16 : 32);
            }
            h.gamma(n);
            h.gamma(unsigned((mode & 2) ? model.empty.size() : model.occupied.size()));
            h.gamma(model.correction + 1);
            bm::id64_t bits = h.bits + maps[mode >> 1] + bounds[mode & 1] + payload[mode & 1];
            bm::id64_t bytes = 1 + ((bits + 31) / 32) * 4;
            if (bytes < best.bytes)
            {
                best.bytes = bytes;
                best.width = w;
                best.mode = mode;
            }
        } // for mode
    } // for w
    return best;
}

} // namespace bm

#endif // BMSBLOCKWAVE__H__INCLUDED__
