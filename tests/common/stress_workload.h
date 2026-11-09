#ifndef BM_TEST_STRESS_WORKLOAD_H
#define BM_TEST_STRESS_WORKLOAD_H

#include <bmbvimport.h>

/** @file stress_workload.h
    @brief Repetition budgets shared by the 32-bit and 64-bit stress programs.

    Test-only helpers; include after bm.h and the definition of is_fast. Apply only to
    audited stress trials and query sweeps, never fixture sizes.
*/

/** Select a repetition budget without increasing the original workload.
    @param normal_count Original number of trials (also used in pedantic mode).
    @param fast_count Maximum trials in fast mode; callers must preserve cases
        selected by trial parity or other trial-dependent state transitions.
    @return The selected count. Does not change mode or random-generator state.
*/
inline
unsigned StressRepetitions(unsigned normal_count, unsigned fast_count)
{
    return is_fast && fast_count < normal_count ? fast_count : normal_count;
}

/** Advance an independent stress-input sweep across its original value range.
    @tparam Size Unsigned input type used by the caller.
    @param value Current input, no greater than last.
    @param last Inclusive upper input bound; last + 1 must be representable.
    @param samples Approximate number of evenly spaced fast samples; positive.
    @param edge_values Number of dense neighboring values at each end, in
        addition to the endpoint; defaults to 32. Must be positive.
    @return Next input, or last + 1 after the final case. Normal mode advances
        by one. Fast mode also retains dense endpoint values (33 by default) and neighbors of
        every power of two. Never use to skip fixture construction or cumulative
        mutations unless the caller restores the equivalent pre-trial state.
        Does not consume random numbers.
*/
template<typename Size>
inline
Size NextStressTestValue(Size value, Size last, unsigned samples = 256,
                        unsigned edge_values = 32)
{
    BM_ASSERT(value <= last && samples && edge_values);
    if (!is_fast || value < edge_values || last - value <= edge_values)
    {
        return value + 1;
    }
    const Size stride = last / samples ? last / samples : 1;
    const Size remaining = last - value;
    Size next = value + (stride < remaining ? stride : remaining);
    if (next > last - edge_values)
    {
        next = last - edge_values;
    }
    for (Size boundary = 1; boundary <= last;)
    {
        const Size candidates[] = { Size(boundary - 1), boundary, Size(boundary + 1) };
        for (Size candidate : candidates)
        {
            if (candidate > value && candidate < next)
            {
                next = candidate;
            }
        } // for candidate
        if (boundary > last / 2)
        {
            break;
        }
        boundary *= 2;
    } // for boundary
    return next;
}

/** Validate every result bit against an independent flat reference bitmap.
    @tparam ReferenceVect Reference bitmap type exposing size() and a borrowed u32 buffer.
    @tparam ResultVect Tested bvector type, also used for the temporary imported bitmap.
    @param reference Borrowed independent bitmap; capacity must be positive.
    @param actual Borrowed result vector; must match every reference bit.
    Fast-only validation for dense stress fixtures. Full bitmap equality retains
    all content checks; redundant per-bit iterator/range probes are sampled.
    Small fixtures use the original exhaustive validator in the owning driver.
    Allocates a temporary bitmap; does not modify inputs or consume randomness.
*/
template<class ReferenceVect, class ResultVect>
inline
void CheckFastDenseVector(const ReferenceVect& reference, const ResultVect& actual)
{
    BM_ASSERT(is_fast && reference.size());
    const auto words = reference.size() / 32 + 1;
    ResultVect imported;
    bm::bit_import_u32(imported, reference.get_buf(), words, true);
    BM_ASSERT(imported.equal(actual));

    // Find the last reference bit independently from the raw word buffer.
    auto remaining = words;
    while (remaining && !reference.get_buf()[remaining - 1])
    {
        --remaining;
    } // while remaining
    typename ResultVect::size_type last = 0;
    const bool found = actual.find_reverse(last);
    BM_ASSERT(found == bool(remaining));
    if (remaining)
    {
        const unsigned word = reference.get_buf()[remaining - 1];
        unsigned bit = 31;
        while (!(word & (1u << bit)))
        {
            --bit;
        } // while bit
        BM_ASSERT(last == (remaining - 1) * 32 + bit);
    }

    auto end = reference.size() - 1;
    for (auto probe = decltype(end)(0); probe <= end;
         probe = NextStressTestValue(probe, end, 16, 1))
    {
        BM_ASSERT(actual.test(probe) == bool(reference.is_bit_true(probe)));
        auto en = actual.get_enumerator(probe);
        auto expected = imported.get_enumerator(probe);
        for (unsigned trial = 0; trial < 3; ++trial)
        {
            BM_ASSERT(en.valid() == expected.valid());
            if (!en.valid())
            {
                break;
            }
            BM_ASSERT(*en == *expected);
            ++en;
            ++expected;
        } // for trial
        const auto to = end - probe < 65 ? end : probe + 65;
        unsigned count = 0;
        for (auto position = probe; position <= to; ++position)
        {
            count += bool(reference.is_bit_true(position));
        } // for position
        BM_ASSERT(actual.count_range(probe, to) == count);
        BM_ASSERT(actual.is_all_one_range(probe, to) == (count == to - probe + 1));
    } // for probe
}

/** Select the next spacing for fast interval-enumerator stress trials.
    @param increment Current positive spacing, less than limit.
    @param limit Exclusive spacing limit from the original test.
    @return Next spacing or limit after the final trial. Includes contiguous and
        alternating bits, word/block boundary neighbors, and limit - 1. The
        caller uses this only in fast mode; no random numbers are consumed.
*/
inline
unsigned NextFastIntervalIncrement(unsigned increment, unsigned limit)
{
    BM_ASSERT(increment && increment < limit);
    const unsigned increments[] =
        { 1, 2, 3, 31, 32, 33, 63, 64, 65, 255, 256, 257,
          65535, 65536, 65537 };
    for (unsigned next : increments)
    {
        if (next > increment && next < limit)
        {
            return next;
        }
    } // for next
    return increment < limit - 1 ? limit - 1 : limit;
}

/** Map reduced interpolation trials across the original spacing-factor range.
    @param trial Zero-based trial index, less than count.
    @param count Selected trial count, at least four and no greater than original.
    @param original Original trial count.
    @return Original spacing-factor index. Keep the special first two factors,
        then include both ends of the remaining range. Does not consume random
        numbers or reduce generator buffer sizes. With count == original the
        mapping is the identity.
*/
inline
unsigned StressRepetitionIndex(unsigned trial, unsigned count, unsigned original)
{
    BM_ASSERT(count >= 4 && count <= original && trial < count);
    if (trial < 2)
    {
        return trial;
    }
    return 2u + unsigned(static_cast<unsigned long long>(trial - 2u) *
                         (original - 3u) / (count - 3u));
}

/** Advance an audited SELECT query sweep without reducing its fixture.
    @tparam Size Unsigned position type used by the caller.
    @param position Current query position, no greater than last.
    @param last Inclusive end of the query range; last + 1 must be representable.
    @return Next position, or last + 1 after the final query. Normal mode advances
        by one. Fast mode keeps the first and last 66 positions and queries on
        either side of every 4096-position boundary. Does not consume randomness.
        Descending sweeps can apply this to the offset from their starting point.
*/
template<typename Size>
inline
Size NextSelectTestPosition(Size position, Size last)
{
    BM_ASSERT(position <= last);
    if (!is_fast || position < 65 || last - position <= 65)
    {
        return position + 1;
    }
    const Size offset = position % 4096;
    const Size step = (offset == 0 || offset == 4095) ? 1 : 4095 - offset;
    // Enter the dense tail even when the next grid boundary is outside it.
    const Size tail = last - 65;
    return step < tail - position ? position + step : tail;
}

#endif // BM_TEST_STRESS_WORKLOAD_H
