#ifndef BM_ADAPTIVE_DELTA_H_INCLUDED
#define BM_ADAPTIVE_DELTA_H_INCLUDED
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

/*! \file bmadaptivedelta.h
    \brief Adaptive minimum-delta correction for serialization (internal).
*/

#include "encoding.h"

namespace bm
{

/**
    \brief Windowed minimum-delta model for strictly increasing coordinates.
    \ingroup bvserial
    \internal

    Each window supplies a correction smaller than its minimum coordinate gap.
    Subtracting cumulative corrections preserves strict ordering while reducing
    the range for interpolative encoding. Window zero uses the caller's global
    correction; later windows use their minimum gap minus one. Position i uses
    window i / step, including the gap across the window boundary. Position zero
    is never corrected.

    \tparam MinimumEntries Minimum length eligible for codec selection. Keep
        separate limits for superblocks and other uses of the transform.

    Correction containers hold unsigned values and provide resize(), size(),
    and operator[] (for example, bm::heap_vector). The model owns no storage.
    Callers may reserve ceil(maximum coordinate count / minimum step) entries
    before repeated calls to avoid reallocating the correction storage.
    Coordinate and output types must represent the corresponding unsigned
    coordinates without truncation; intermediate arithmetic uses unsigned.
*/
template<unsigned MinimumEntries>
struct adaptive_delta_model
{
    static constexpr unsigned minimum_entries = MinimumEntries; ///< Selection threshold.

    /**
        \brief Test whether the codec may consider this window size.
        \param n Number of input coordinates.
        \param step Positive number of coordinates per window.
        \return True if n meets the threshold and contains two complete windows.
    */
    static bool eligible(unsigned n, unsigned step) BMNOEXCEPT;

    /**
        \brief Compute window corrections and changed-correction coding costs.
        \tparam Value Input coordinate type.
        \tparam CorrectionsVect Resizable container of unsigned corrections.
        \param values Non-null, strictly increasing array of n coordinates.
        \param n Positive coordinate count; eligibility is checked separately.
        \param step Positive number of coordinates per window.
        \param global Global correction, at most every input gap minus one.
        \param corrections Output, resized to ceil(n / step); element zero is global.
        \param scores Output bit costs: [0] absolute gamma codes, [1] reserved
            (zero), [2] sign plus gamma-coded difference from previous window.
            Costs exclude flags and the separately encoded global correction.
        \return Number of corrections differing from the previous window.
    */
    template<class Value, class CorrectionsVect>
    static unsigned prepare(const Value* values, unsigned n, unsigned step,
                            unsigned global, CorrectionsVect& corrections,
                            bm::id64_t (&scores)[3]);

    /**
        \brief Subtract cumulative window corrections from coordinates.
        \tparam Value Input coordinate type.
        \tparam CorrectionsVect Container of unsigned corrections.
        \tparam Output Destination coordinate type.
        \param values Non-null, strictly increasing array of n coordinates.
        \param n Positive number of coordinates.
        \param step Positive number of coordinates per window.
        \param corrections Exactly ceil(n / step) corrections from prepare().
        \param out Non-null destination for n coordinates. Separate buffers or
            exact in-place operation are supported; partial overlap is not.
        \return Total correction subtracted from the last coordinate.
    */
    template<class Value, class CorrectionsVect, class Output>
    static unsigned transform(const Value* values, unsigned n, unsigned step,
                              const CorrectionsVect& corrections, Output* out);

    /**
        \brief Restore coordinates by adding cumulative window corrections.
        \tparam Value Corrected input coordinate type.
        \tparam CorrectionsVect Container of unsigned corrections.
        \tparam Output Restored coordinate type.
        \param values Non-null, strictly increasing array of n corrected coordinates.
        \param n Positive number of coordinates.
        \param step Positive number of coordinates per window.
        \param corrections Exactly ceil(n / step) corrections used by transform().
        \param out Non-null destination for n coordinates. Separate buffers or
            exact in-place operation are supported; partial overlap is not.
        \return Total correction added to the last coordinate.
        \note Validate untrusted coordinates and corrections before calling.
    */
    template<class Value, class CorrectionsVect, class Output>
    static unsigned restore(const Value* values, unsigned n, unsigned step,
                            const CorrectionsVect& corrections, Output* out);
};

/// Internal superblock selection policy; other uses must not relax its limit.
/// \ingroup bvserial
/// \internal
using superblock_delta_model = adaptive_delta_model<64>;
static_assert(superblock_delta_model::minimum_entries == 64,
              "superblock minimum must remain 64");

/**
    \brief Encode adaptive minimum-delta correction flags and values.
    \ingroup bvserial
    \internal

    Writes correction metadata only; coordinate transformation and payload
    encoding belong to the caller. The two stages preserve grouped flags
    followed by changed values and permit separate coding-cost measurements.
    This stateless helper borrows caller-owned storage without allocating.
*/
struct adaptive_delta_encoder
{
    /**
        \brief Write grouped correction-change flags, excluding window zero.
        \ingroup bvserial
        \internal
        \tparam BitOut Bit writer providing put_bit().
        \tparam CorrectionsVect Container with size() and operator[].
        \param out Destination bit stream.
        \param c Nonempty unsigned correction array; size must fit unsigned.
            Eight consecutive flags share an all-changed prefix bit; the final
            incomplete group uses individual flags without a prefix.
    */
    template<class BitOut, class CorrectionsVect>
    static void encode_flags(BitOut& out, const CorrectionsVect& c);

    /**
        \brief Write changed corrections after their grouped flags.
        \ingroup bvserial
        \internal
        \tparam BitOut Bit writer providing put_bit() and gamma().
        \tparam CorrectionsVect Container with size() and operator[].
        \param out Destination bit stream, positioned after the change flags.
        \param c Nonempty unsigned correction array; all entries are at least c[0].
            Absolute-mode entries must be less than the maximum unsigned value.
        \param mode 0 for gamma(correction + 1), 2 for a sign bit followed by the
            gamma-coded magnitude of change from the previous window.
    */
    template<class BitOut, class CorrectionsVect>
    static void encode_values(BitOut& out, const CorrectionsVect& c,
                              unsigned mode);
};

/**
    \brief Decode grouped adaptive minimum-delta correction metadata.
    \ingroup bvserial
    \internal

    Reads both flags and changed values into caller-owned correction storage.
    Storage capacity is retained between calls. Coordinate restoration and
    format-specific coordinate validation belong to the caller.
*/
struct adaptive_delta_decoder
{
    /**
        \brief Read grouped flags and changed corrections, validating stream bounds.
        \ingroup bvserial
        \internal
        \tparam BitIn Bit reader providing get_bit() and gamma().
        \tparam CorrectionsVect Resizable container of unsigned corrections.
        \param in Source bit stream, positioned at the grouped change flags.
        \param windows Positive number of windows, including window zero.
        \param global Previously decoded correction for window zero.
        \param mode 0 for absolute gamma, 2 for signed previous-relative gamma.
        \param maximum Inclusive upper bound on corrections and decoded gamma codes.
            Thus absolute mode permits changed corrections only up to maximum - 1.
        \param c Output, resized to windows; element zero is global. Its contents
            are unspecified if stream validation fails.
            Pre-reserve the maximum window count when reusing this container.
        \note Bounds and mode checks remain active when BM_ASSERT is disabled.
    */
    template<class BitIn, class CorrectionsVect>
    static void decode(BitIn& in, unsigned windows, unsigned global,
                       unsigned mode, unsigned maximum, CorrectionsVect& c);
};


// -----------------------------------------------------------------------
// Implementation
// -----------------------------------------------------------------------

template<unsigned MinimumEntries>
inline
bool adaptive_delta_model<MinimumEntries>::eligible(
                                    unsigned n, unsigned step) BMNOEXCEPT
{
    BM_ASSERT(step);
    return step && n >= MinimumEntries && step <= n / 2;
}

template<unsigned MinimumEntries>
template<class Value, class CorrectionsVect>
inline
unsigned adaptive_delta_model<MinimumEntries>::prepare(
                    const Value* values, unsigned n, unsigned step,
                    unsigned global, CorrectionsVect& corrections,
                    bm::id64_t (&scores)[3])
{
    BM_ASSERT(values && n && step);
    BM_ASSERT(global < ~0u);
    const unsigned windows = 1 + (n - 1) / step;
    corrections.resize(windows);
    corrections[0] = global;
    scores[0] = scores[1] = scores[2] = 0;
    unsigned changes = 0;
    for (unsigned w = 1; w < windows; ++w)
    {
        unsigned gap = ~0u;
        const unsigned begin = w * step;
        const unsigned remaining = n - begin;
        const unsigned end = begin + (remaining < step ? remaining : step);
        for (unsigned i = begin; i < end; ++i)
        {
            BM_ASSERT(values[i] > values[i-1]);
            const unsigned delta = unsigned(values[i] - values[i-1]);
            BM_ASSERT(Value(delta) == values[i] - values[i-1]);
            gap = delta < gap ? delta : gap;
        } // for i
        BM_ASSERT(gap && global < gap);
        const unsigned c = corrections[w] = gap - 1;
        const unsigned prev = corrections[w-1];
        if (c == prev)
            continue;
        ++changes;
        scores[0] += bm::gamma_size(c + 1);
        scores[2] += 1 + bm::gamma_size(c > prev ? c-prev : prev-c);
    } // for w
    return changes;
}

template<unsigned MinimumEntries>
template<class Value, class CorrectionsVect, class Output>
inline
unsigned adaptive_delta_model<MinimumEntries>::transform(
                    const Value* values, unsigned n, unsigned step,
                    const CorrectionsVect& corrections, Output* out)
{
    BM_ASSERT(values && out && n && step);
    BM_ASSERT(corrections.size() == 1 + (n - 1) / step);
    unsigned sum = 0;
    const unsigned first = unsigned(values[0]);
    BM_ASSERT(bm::id64_t(values[0]) == first);
    out[0] = Output(first);
    BM_ASSERT(bm::id64_t(out[0]) == first);
    for (unsigned i = 1; i < n; ++i)
    {
        const unsigned value = unsigned(values[i]);
        BM_ASSERT(bm::id64_t(values[i]) == value);
        const unsigned c = corrections[i / step];
        BM_ASSERT(c <= ~0u - sum);
        sum += c;
        BM_ASSERT(sum <= value);
        out[i] = Output(value - sum);
        BM_ASSERT(bm::id64_t(out[i]) == value - sum);
        BM_ASSERT(out[i] > out[i-1]);
    } // for i
    return sum;
}

template<unsigned MinimumEntries>
template<class Value, class CorrectionsVect, class Output>
inline
unsigned adaptive_delta_model<MinimumEntries>::restore(
                    const Value* values, unsigned n, unsigned step,
                    const CorrectionsVect& corrections, Output* out)
{
    BM_ASSERT(values && out && n && step);
    BM_ASSERT(corrections.size() == 1 + (n - 1) / step);
    unsigned sum = 0;
    const unsigned first = unsigned(values[0]);
    BM_ASSERT(bm::id64_t(values[0]) == first);
    out[0] = Output(first);
    BM_ASSERT(bm::id64_t(out[0]) == first);
    for (unsigned i = 1; i < n; ++i)
    {
        const unsigned value = unsigned(values[i]);
        BM_ASSERT(bm::id64_t(values[i]) == value);
        const unsigned c = corrections[i / step];
        BM_ASSERT(c <= ~0u - sum);
        sum += c;
        BM_ASSERT(value <= ~0u - sum);
        out[i] = Output(value + sum);
        BM_ASSERT(bm::id64_t(out[i]) == value + sum);
        BM_ASSERT(out[i] > out[i-1]);
    } // for i
    return sum;
}

template<class BitOut, class CorrectionsVect>
inline
void adaptive_delta_encoder::encode_flags(BitOut& out, const CorrectionsVect& c)
{
    BM_ASSERT(c.size() && c.size() <= ~0u);
    unsigned base = 1, n = unsigned(c.size());
    for (; n - base >= 8; base += 8)
    {
        bool ones = true;
        for (unsigned j = 0; j < 8; ++j)
        {
            ones &= c[base+j] != c[base+j-1];
        } // for j
        out.put_bit(ones);
        if (!ones)
        {
            for (unsigned j = 0; j < 8; ++j)
            {
                out.put_bit(c[base+j] != c[base+j-1]);
            } // for j
        }
    } // for base
    for (; base < n; ++base)
    {
        out.put_bit(c[base] != c[base-1]);
    } // for base
}

template<class BitOut, class CorrectionsVect>
inline
void adaptive_delta_encoder::encode_values(BitOut& out, const CorrectionsVect& c,
                                           unsigned mode)
{
    BM_ASSERT(c.size() && c.size() <= ~0u);
    BM_ASSERT(mode == 0 || mode == 2);
    const unsigned windows = unsigned(c.size());
    for (unsigned w = 1; w < windows; ++w)
    {
        const unsigned value = c[w], prev = c[w-1];
        BM_ASSERT(value >= c[0]);
        if (value == prev)
            continue;
        if (mode == 2)
        {
            out.put_bit(value < prev);
            out.gamma(value > prev ? value-prev : prev-value);
        }
        else
        {
            BM_ASSERT(value < ~0u);
            out.gamma(value + 1);
        }
    } // for w
}

template<class BitIn, class CorrectionsVect>
inline
void adaptive_delta_decoder::decode(BitIn& in, unsigned windows, unsigned global,
                                    unsigned mode, unsigned maximum,
                                    CorrectionsVect& c)
{
    bm::serialization_require(windows && global <= maximum && (mode == 0 || mode == 2));
    c.resize(windows);
    c[0] = global;
    unsigned base = 1;
    for (; windows - base >= 8; base += 8)
    {
        if (in.get_bit())
        {
            for (unsigned j = 0; j < 8; ++j)
            {
                c[base+j] = 1;
            } // for j
        }
        else
        {
            for (unsigned j = 0; j < 8; ++j)
            {
                c[base+j] = in.get_bit();
            } // for j
        }
    } // for base
    for (; base < windows; ++base)
    {
        c[base] = in.get_bit();
    } // for base
    for (unsigned w = 1; w < windows; ++w)
    {
        const unsigned prev = c[w-1];
        if (!c[w])
        {
            c[w] = prev;
            continue;
        }
        const bool negative = mode == 2 && in.get_bit();
        const unsigned code = in.gamma();
        bm::serialization_require(code && code <= maximum);
        if (mode == 2)
        {
            bm::serialization_require(negative ? code <= prev : code <= maximum-prev);
            c[w] = negative ? prev-code : prev+code;
        }
        else
            c[w] = code - 1;
        bm::serialization_require(c[w] >= global && c[w] <= maximum);
    } // for w
}

} // namespace bm
#endif
