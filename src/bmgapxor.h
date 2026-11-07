#ifndef BM_GAP_XOR_H_INCLUDED
#define BM_GAP_XOR_H_INCLUDED
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

For more information please visit:  http://bitmagic.io
*/

/*! \file bmgapxor.h
    \brief GAP cleanup and XOR-exception transforms for serialization (internal).

    Include after bm.h or bm64.h. Helpers prepare candidates for
    bm::set_block_gap_xor_framed (opcode 94); they do not read or write its
    wire format. A cleaned GAP base XOR the odd-parity exception positions
    reconstructs the original block. Uniform bases are mathematically valid
    even when the enclosing codec cannot encode them as nested GAP records.
*/

#include "bmfunc.h"

namespace bm
{

/**
    \brief Reduce arbitrary XOR edits to sorted positions with odd multiplicity.
    \ingroup bvserial
    \internal
    \tparam PositionsVect Resizable contiguous container of integer positions
        providing data(), size(), resize(), and operator[].
    \param positions Input edits, possibly empty, unordered, or duplicated.
        Replaced by strictly increasing positions occurring an odd number of
        times in the input. Even multiplicities cancel under XOR.
    \note In-place heapsort takes O(n log n) time and constant scratch space,
        including for ordered runs and duplicate-heavy inputs. Compaction
        shrinks the logical size; BM vectors retain their allocated capacity.
*/
template<class PositionsVect>
inline
void gap_xor_cancel(PositionsVect& positions);

/**
    \brief Split a source GAP block into a cleaned base and XOR exceptions.
    \ingroup bvserial
    \internal
    \tparam CleanedVect Resizable contiguous container of bm::gap_word_t.
    \tparam PositionsVect Resizable contiguous container of bm::gap_word_t.
    \param gap Read-only, valid GAP block including its header and terminal 65535.
    \param split Candidate policy in 0..3: 0 removes single-bit zero runs,
        1 removes single-bit one runs, 2 tries zero runs before one runs,
        and 3 tries one runs before zero runs.
    \param base Writable GAP scratch for bm::gap_split(), disjoint from gap,
        exception arrays, and vector storage; must accommodate the complete GAP.
    \param ex0 Writable scratch for zero-run exceptions, disjoint from all other
        buffers; must accommodate every exception produced by bm::gap_split().
    \param ex1 Writable scratch for one-run exceptions, with the same capacity
        and non-overlap requirements as ex0.
    \param cleaned Output copy of the cleaned GAP, including header and terminator.
    \param positions Output sorted, unique exception positions after XOR cancellation.
        Both output vectors are exclusive caller-owned storage and must be distinct.
    \note Scratch arrays have no capacity arguments: the caller must size them
        before the call. Reserving the largest expected cleaned GAP and combined
        exception count in the output vectors avoids repeated growth. A uniform
        cleaned base and an empty exception list are valid transform results.
*/
template<class CleanedVect, class PositionsVect>
inline
void gap_xor_prepare_candidate(const bm::gap_word_t* gap, unsigned split,
                               bm::gap_word_t* base, bm::gap_word_t* ex0,
                               bm::gap_word_t* ex1, CleanedVect& cleaned,
                               PositionsVect& positions);

/**
    \brief Subtract a uniform minimum-gap correction from exception coordinates.
    \ingroup bvserial
    \internal
    \param positions Strictly increasing array of n block-local bit indexes.
    \param n Coordinate count in 2..bm::gap_max_bits.
    \param corrected Destination for n corrected coordinates. Exact in-place
        operation is supported; otherwise the two arrays must be disjoint.
    \return Minimum adjacent input gap minus one. Original coordinate i is
        restored as corrected[i] + i * correction.
    \note Preserves strict ordering and the first coordinate while reducing
        the range for BIC encoding. No allocation is performed.
*/
inline
unsigned gap_xor_correct(const bm::gap_word_t* positions, unsigned n,
                         bm::gap_word_t* corrected);


// -----------------------------------------------------------------------
// Implementation
// -----------------------------------------------------------------------

template<class PositionsVect>
inline
void gap_xor_cancel(PositionsVect& positions)
{
    const size_t count = positions.size();
    bm::heap_sort(positions.data(), count);
    size_t out = 0;
    for (size_t i = 0; i < count;)
    {
        size_t j = i + 1;
        while (j < count && positions[j] == positions[i])
        {
            ++j;
        } // while j
        if ((j - i) & 1u)
        {
            positions[out++] = positions[i];
        }
        i = j;
    } // for i
    positions.resize(out);
}

template<class CleanedVect, class PositionsVect>
inline
void gap_xor_prepare_candidate(const bm::gap_word_t* gap, unsigned split,
                               bm::gap_word_t* base, bm::gap_word_t* ex0,
                               bm::gap_word_t* ex1, CleanedVect& cleaned,
                               PositionsVect& positions)
{
    BM_ASSERT(gap && base && ex0 && ex1 && split < 4);
    BM_ASSERT(gap != base && gap != ex0 && gap != ex1);
    BM_ASSERT(base != ex0 && base != ex1 && ex0 != ex1);
    BM_ASSERT(bm::gap_length(gap) >= 2);
    BM_ASSERT(gap[bm::gap_length(gap) - 1] == bm::gap_max_bits - 1);
    const unsigned hist0[2] = {0, split == 1 ? 0u : split == 2 ? 2u : 1u};
    const unsigned hist1[2] = {0, split == 0 ? 0u : split == 3 ? 2u : 1u};
    unsigned n0 = 0;
    unsigned n1 = 0;
    bm::gap_split(gap, bm::gap_length(gap) - 1, 2, hist0, hist1,
                  base, ex0, ex1, n0, n1);
    cleaned.resize(bm::gap_length(base));
    ::memcpy(cleaned.data(), base, cleaned.size() * sizeof(bm::gap_word_t));
    positions.resize(n0 + n1);
    if (n0)
    {
        ::memcpy(positions.data(), ex0, n0 * sizeof(bm::gap_word_t));
    }
    if (n1)
    {
        ::memcpy(positions.data() + n0, ex1, n1 * sizeof(bm::gap_word_t));
    }
    bm::gap_xor_cancel(positions);
}

inline
unsigned gap_xor_correct(const bm::gap_word_t* positions, unsigned n,
                         bm::gap_word_t* corrected)
{
    BM_ASSERT(positions && corrected && n >= 2 && n <= bm::gap_max_bits);
    unsigned min_delta = bm::gap_max_bits;
    for (unsigned i = 1; i < n; ++i)
    {
        BM_ASSERT(positions[i] > positions[i-1]);
        const unsigned delta = unsigned(positions[i] - positions[i-1]);
        min_delta = delta < min_delta ? delta : min_delta;
    } // for i
    const unsigned correction = min_delta - 1;
    for (unsigned i = 0; i < n; ++i)
    {
        BM_ASSERT(i * correction <= unsigned(positions[i]));
        corrected[i] = bm::gap_word_t(positions[i] - i * correction);
        BM_ASSERT(!i || corrected[i] > corrected[i-1]);
    } // for i
    return correction;
}

} // namespace bm
#endif
