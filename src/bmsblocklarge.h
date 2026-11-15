#ifndef BMSBLOCKLARGE__H__INCLUDED__
#define BMSBLOCKLARGE__H__INCLUDED__

/** \file bmsblocklarge.h
    \brief Internal high-count superblock BIC payload encoder and decoder.
    \ingroup bvserial
    \internal

    Implements bm::set_sblock_bienc_large (98). Positive records store set-bit
    positions; inverse records store clear-bit positions within a superblock.
    The caller writes the opcode and applies the decoded polarity. The legacy
    bm::set_sblock_bienc_v3 (68) format is independent of this payload.
    Payload codecs are independent of the serializer; whole-vector candidate
    selection is implemented in bmserial.h.
*/

#include "encoding.h"

namespace bm
{
/** \brief Internal wire-format constants for high-count superblocks.
    \ingroup bvserial
    \internal
*/
struct sblock_large_format
{
    static constexpr unsigned minimum = 65536; ///< Minimum positive-list cardinality.
    static constexpr unsigned maximum = 196607; ///< Maximum positive-list cardinality.

    /** \brief Internal flags following the leading polarity bit.
        \ingroup bvserial
        \internal
        The low two bits select the superblock-index representation. Endpoint
        flags select the first position and the distance from the last position
        to the superblock end. Positive lists use range3 to select one of two
        16-bit cardinality residual ranges; inverse lists use it to select a
        16-bit or 24-bit literal cardinality.
    */
    enum
    {
        sb16 = 1,       ///< Superblock index occupies 16 bits.
        sb32 = 2,       ///< Superblock index occupies 32 bits.
        sbgamma = 3,    ///< Superblock index plus one is gamma-coded.
        min16 = 4,      ///< First corrected position occupies 16 bits.
        min24 = 8,      ///< First position uses gamma high bits and 16 low bits.
        range3 = 16,    ///< Upper positive count range or 24-bit inverse count.
        max16 = 32,     ///< Corrected last-position tail occupies 16 bits.
        max24 = 64,     ///< Corrected last-position tail occupies 24 bits.
        corrected = 128 ///< Positions have a uniform minimum-gap correction.
    };
};

/** \brief Internal high-count superblock BIC payload encoder.
    \ingroup bvserial
    \internal
    \tparam OUT Bit-output type supporting gamma, BIC and fixed-width fields.
*/
template<class OUT>
struct sblock_large_encoder
{
    /** \brief Encode one positive or inverse high-count superblock payload.
        \ingroup bvserial
        \internal
        \tparam ValuesVect Resizable unsigned-position scratch container, using
            the caller's application allocator.
        \param out Destination bit accumulator; the caller owns its final flush.
        \param original Strictly increasing superblock-relative positions.
        \param n Number of positions; positive lists require minimum through
            maximum, while inverse lists may be empty but cannot fill a superblock.
        \param sb Absolute superblock index.
        \param correction Uniform gap reduction: corrected position i equals
            original[i] minus i times correction. The result must remain strictly
            increasing; correction must fit in 16 bits and be zero for n below two.
        \param values Reusable scratch receiving corrected positions. Its logical
            size becomes n; retained capacity avoids allocations on subsequent calls.
        \param inverse Encode clear-bit positions when true, set-bit positions otherwise.
        \note The first and last corrected positions are explicit; BIC encodes only
            the interior positions. Inverse empty lists retain endpoint placeholders.
    */
    template<class ValuesVect>
    static void encode(OUT& out, const unsigned* original, unsigned n,
                       unsigned sb, unsigned correction, ValuesVect& values,
                       bool inverse = false);
};

/** \brief Internal high-count superblock BIC payload decoder.
    \ingroup bvserial
    \internal
    \tparam IN Bit-input type supporting gamma, BIC and fixed-width fields.
*/
template<class IN>
struct sblock_large_decoder
{
    /** \brief Decode positions and metadata from one high-count superblock payload.
        \ingroup bvserial
        \internal
        \tparam ValuesVect Resizable unsigned-position container using the caller's
            application allocator.
        \param in Source bit accumulator positioned after the opcode.
        \param values Reusable output scratch. The first returned-count elements
            contain strictly increasing, restored superblock-relative positions.
            Its logical size is at least 65,536 because subsequent legacy records
            may reuse the same storage; capacity is retained across calls.
        \param sb Receives the absolute superblock index.
        \param inverse Receives true for clear-bit positions, false for set-bit positions.
        \return Number of valid positions in values, which may be zero.
        \note Format checks use bm::serialization_require(). BIC recursion has
            logarithmic depth; the position array is stored in the supplied container.
    */
    template<class ValuesVect>
    static unsigned decode(IN& in, ValuesVect& values, unsigned& sb, bool& inverse);
};

// ----------------------------------------------------------------------

template<class OUT>
template<class ValuesVect>
inline
void sblock_large_encoder<OUT>::encode(OUT& out, const unsigned* original,
                                     unsigned n, unsigned sb, unsigned correction,
                                     ValuesVect& values, bool inverse)
{
    bm::serialization_require(inverse ? n < bm::set_sub_total_bits :
        n >= bm::sblock_large_format::minimum && n <= bm::sblock_large_format::maximum);
    values.resize(n);
    if (n)
    {
        ::memcpy(values.data(), original, size_t(n) * sizeof(unsigned));
    }
    for (unsigned i = 0; i < n; ++i)
    {
        values[i] -= i * correction;
    } // for i

    unsigned lo = n ? values[0] : 0;
    unsigned tail = n ? bm::set_sub_total_bits - values[n - 1] : 1;
    unsigned flags = ((inverse ? n >= bm::sblock_large_format::minimum :
                      n >= 2 * bm::sblock_large_format::minimum) ?
                      bm::sblock_large_format::range3 : 0) |
                     (correction ? bm::sblock_large_format::corrected : 0);
    flags |= sb < 8 ? bm::sblock_large_format::sbgamma :
             sb > 65535 ? bm::sblock_large_format::sb32 :
             sb > 255 ? bm::sblock_large_format::sb16 : 0;
    flags |= lo > 65535 ? bm::sblock_large_format::min24 :
             lo > 255 ? bm::sblock_large_format::min16 : 0;
    flags |= tail > 65535 ? bm::sblock_large_format::max24 :
             tail > 255 ? bm::sblock_large_format::max16 : 0;

    out.put_bits(inverse ? 1u : 0u, 1);
    out.put_bits(flags, 8);
    if (inverse)
    {
        if (flags & bm::sblock_large_format::range3)
        {
            out.put_24_no(n);
        }
        else
        {
            out.put_16_no(static_cast<unsigned short>(n));
        }
    }
    else
    {
        out.put_16_no(static_cast<unsigned short>(n -
            (n >= 2 * bm::sblock_large_format::minimum ?
             2 * bm::sblock_large_format::minimum : bm::sblock_large_format::minimum)));
    }
    if (flags & bm::sblock_large_format::min24)
    {
        out.gamma(lo >> 16);
        out.put_16_no(static_cast<unsigned short>(lo & 65535));
    }
    else if (flags & bm::sblock_large_format::min16)
    {
        out.put_16_no(static_cast<unsigned short>(lo));
    }
    else
    {
        out.put_bits(lo, 8);
    }

    if (flags & bm::sblock_large_format::max24)
    {
        out.put_24_no(tail);
    }
    else if (flags & bm::sblock_large_format::max16)
    {
        out.put_16_no(static_cast<unsigned short>(tail));
    }
    else
    {
        out.put_bits(tail, 8);
    }

    if (correction)
    {
        unsigned code = correction < 16 ? 1 : correction < 256 ? 2 : 3;
        out.gamma(code);
        if (code == 1)
        {
            out.gamma(correction);
        }
        else if (code == 2)
        {
            out.put_bits(correction, 8);
        }
        else
        {
            out.put_16_no(static_cast<unsigned short>(correction));
        }
    }
    if ((flags & 3) == bm::sblock_large_format::sbgamma)
    {
        out.gamma(sb + 1);
    }
    else if (flags & bm::sblock_large_format::sb32)
    {
        out.put_32_no(sb);
    }
    else if (flags & bm::sblock_large_format::sb16)
    {
        out.put_16_no(static_cast<unsigned short>(sb));
    }
    else
    {
        out.put_bits(sb, 8);
    }

    if (n > 2)
    {
        out.bic_encode_u32_cm(values.data() + 1, n - 2, lo + 1, values[n - 1] - 1);
    }
}

template<class IN>
template<class ValuesVect>
inline
unsigned sblock_large_decoder<IN>::decode(IN& in, ValuesVect& values,
                                         unsigned& sb, bool& inverse)
{
    inverse = bool(in.get_bits(1));
    unsigned flags = in.get_bits(8);
    unsigned n = inverse ? ((flags & bm::sblock_large_format::range3) ?
                           in.get_24_no() : in.get_16_no()) :
        bm::sblock_large_format::minimum + ((flags & bm::sblock_large_format::range3) ?
                                    bm::sblock_large_format::minimum : 0) + in.get_16_no();
    bm::serialization_require(
        (flags & (bm::sblock_large_format::min16 | bm::sblock_large_format::min24)) !=
                 (bm::sblock_large_format::min16 | bm::sblock_large_format::min24) &&
        (flags & (bm::sblock_large_format::max16 | bm::sblock_large_format::max24)) !=
                 (bm::sblock_large_format::max16 | bm::sblock_large_format::max24));
    unsigned lo;
    if (flags & bm::sblock_large_format::min24)
    {
        unsigned block = in.gamma();
        bm::serialization_require(block < 256);
        lo = (block << 16) + in.get_16_no();
    }
    else
    {
        lo = (flags & bm::sblock_large_format::min16) ? in.get_16_no() : in.get_bits(8);
    }

    unsigned tail = (flags & bm::sblock_large_format::max24) ? in.get_24_no() :
                   (flags & bm::sblock_large_format::max16) ? in.get_16_no() : in.get_bits(8);
    bm::serialization_require(tail && tail <= bm::set_sub_total_bits);
    unsigned hi = bm::set_sub_total_bits - tail;
    unsigned c = 0;
    if (flags & bm::sblock_large_format::corrected)
    {
        unsigned code = in.gamma();
        bm::serialization_require(code >= 1 && code <= 3);
        c = code == 1 ? in.gamma() : code == 2 ? in.get_bits(8) : in.get_16_no();
        bm::serialization_require(c && c <= 65535);
    }
    if ((flags & 3) == bm::sblock_large_format::sbgamma)
    {
        unsigned v = in.gamma();
        bm::serialization_require(v);
        sb = v - 1;
    }
    else
    {
        sb = (flags & bm::sblock_large_format::sb32) ? in.get_32_no() :
             (flags & bm::sblock_large_format::sb16) ? in.get_16_no() : in.get_bits(8);
    }

    bm::serialization_require((n == 0 || (lo <= hi && hi - lo >= n - 1)) &&
        (n >= 2 || !c) &&
        bm::id64_t(hi) + bm::id64_t(n ? n - 1 : 0) * c < bm::set_sub_total_bits &&
        bm::id64_t(sb) <= bm::id64_t(bm::id_max) / bm::set_sub_total_bits);

    // Retain space for subsequent legacy records which reuse this pointer.
    // BIC recursion is log2(n); scratch storage stays on the heap.
    values.resize(n > 65536u ? n : 65536u);
    if (n)
    {
        values[0] = lo;
        values[n - 1] = hi;
    }
    if (n == 1)
    {
        bm::serialization_require(lo == hi);
    }
    if (n > 2)
    {
        in.bic_decode_u32_cm(values.data() + 1, n - 2, lo + 1, hi - 1);
    }

    for (unsigned i = 0; i < n; ++i)
    {
        values[i] += i * c;
        bm::serialization_require(values[i] < bm::set_sub_total_bits &&
                                  (!i || values[i] > values[i - 1]));
    } // for i
    return n;
}

// ----------------------------------------------------------------------



} // namespace bm

#endif // BMSBLOCKLARGE__H__INCLUDED__
