#ifndef BMGAPROUTE__H__INCLUDED__
#define BMGAPROUTE__H__INCLUDED__
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

/*! \file bmgaproute.h
    \brief Optional accounting of emitted GAP serialization routes (internal).

    Include after bm.h or bm64.h, as for the other serialization utilities.
    Enable collection through serializer::set_gap_route_diagnostics(). The
    counters describe emitted routes; they do not affect codec selection.
    Text presentation is provided separately by bmdbg.h.
*/

namespace bm
{

/**
    \brief Counts and sizes for GAP blocks taking one serialization route.
    \ingroup bvserial
    \internal

    Word counts describe the source GAP representation in 16-bit words,
    not its compressed wire size. Exact bytes are collected for direct
    records only; grouped records share a superblock byte total.
*/
struct gap_route_bucket
{
    bm::id64_t blocks = 0; ///< Number of source GAP blocks in this bucket.
    bm::id64_t words = 0;  ///< Sum of source GAP lengths in 16-bit words.
    bm::id64_t bytes = 0;  ///< Exact direct-record bytes; unused for grouped buckets.
};

/**
    \brief Accumulated source-block and emitted-record route measurements.
    \ingroup bvserial
    \internal

    Direct/grouped arrays use the source GAP allocation level (0..3), not
    the serializer compression level. The level-zero length buckets refine
    that level into 1..4, 5..8, 9..16, 17..32, 33..64, and 65+ GAP words.
    They overlap the level-zero totals and must not be added to them.
    Storage is fixed-size, zero-initialized, and requires no allocator.
*/
struct gap_route_statistics
{
    bm::id64_t superblocks = 0;      ///< Sparse superblock records accounted for by this collector.
    bm::id64_t superblock_bytes = 0; ///< Exact combined bytes of those superblock records.
    bm::id64_t grouped_bit = 0;      ///< Source BIT blocks absorbed into those superblocks.
    bm::id64_t grouped_full = 0;     ///< Source full blocks absorbed into those superblocks.
    bm::id64_t grouped_empty = 0;    ///< Source empty blocks absorbed into those superblocks.
    bm::id64_t direct_bit = 0;       ///< Source BIT blocks emitted through the direct block codec.
    bm::id64_t direct_bit_bytes = 0; ///< Exact bytes emitted for direct BIT block records.
    bm::id64_t folded_gap = 0;       ///< Visited GAP blocks folded into zero/full representations.
    bm::gap_route_bucket direct[4];     ///< Direct GAP records by source GAP allocation level.
    bm::gap_route_bucket grouped[4];    ///< Superblock-contained GAP blocks by source allocation level.
    bm::gap_route_bucket direct_l0[6];  ///< Direct level-zero GAP records by source word-count bucket.
    bm::gap_route_bucket grouped_l0[6]; ///< Grouped level-zero GAP blocks by source word-count bucket.
};

/**
    \brief Classify a GAP word count into one of six diagnostic buckets.
    \ingroup bvserial
    \internal
    \param len Source GAP length in 16-bit words. Zero is accepted and maps
        to bucket zero, although normal GAP observations have positive lengths.
    \return Index 0..5 for lengths <=4, <=8, <=16, <=32, <=64, and >64.
*/
inline
unsigned gap_route_length_bucket(unsigned len) BMNOEXCEPT;

/**
    \brief Accumulate GAP-route counters from another serializer or worker.
    \ingroup bvserial
    \internal
    \param d Destination accumulator; existing counts are preserved and added to.
    \param s Source counters. Self-addition is supported and doubles counts.
    \note Grouped-bucket bytes are unused and are not merged; their wire size
        is tracked by superblock_bytes. Counter sums must fit bm::id64_t.
        No allocation is performed.
*/
inline
void gap_route_add(bm::gap_route_statistics& d,
                   const bm::gap_route_statistics& s) BMNOEXCEPT;


// -----------------------------------------------------------------------
// Implementation
// -----------------------------------------------------------------------

inline
unsigned gap_route_length_bucket(unsigned len) BMNOEXCEPT
{
    if (len <= 4) return 0;
    if (len <= 8) return 1;
    if (len <= 16) return 2;
    if (len <= 32) return 3;
    if (len <= 64) return 4;
    return 5;
}

inline
void gap_route_add(bm::gap_route_statistics& d,
                   const bm::gap_route_statistics& s) BMNOEXCEPT
{
    d.superblocks += s.superblocks;
    d.superblock_bytes += s.superblock_bytes;
    d.grouped_bit += s.grouped_bit;
    d.grouped_full += s.grouped_full;
    d.grouped_empty += s.grouped_empty;
    d.direct_bit += s.direct_bit;
    d.direct_bit_bytes += s.direct_bit_bytes;
    d.folded_gap += s.folded_gap;
    for (unsigned i = 0; i < 4; ++i)
    {
        d.direct[i].blocks += s.direct[i].blocks;
        d.direct[i].words += s.direct[i].words;
        d.direct[i].bytes += s.direct[i].bytes;
        d.grouped[i].blocks += s.grouped[i].blocks;
        d.grouped[i].words += s.grouped[i].words;
    } // for i
    for (unsigned i = 0; i < 6; ++i)
    {
        d.direct_l0[i].blocks += s.direct_l0[i].blocks;
        d.direct_l0[i].words += s.direct_l0[i].words;
        d.direct_l0[i].bytes += s.direct_l0[i].bytes;
        d.grouped_l0[i].blocks += s.grouped_l0[i].blocks;
        d.grouped_l0[i].words += s.grouped_l0[i].words;
    } // for i
}

} // namespace bm
#endif
