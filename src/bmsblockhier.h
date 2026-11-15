#ifndef BMSBLOCKHIER__H__INCLUDED__
#define BMSBLOCKHIER__H__INCLUDED__
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

/** @file bmsblockhier.h
    @brief Internal hierarchical wave-compacted superblock codec and wave selection.
    @internal @ingroup bvserial

    Implements opcode 109 (bm::set_sblock_bienc_wave_hierarchical). Occupied
    blocks are compacted first, then occupied waves within those blocks. Global
    gap correction reduces the remaining position domain. Maps and positions
    share one bit accumulator; the caller writes the opcode and final padding.
    Six wave widths, 128 through 4096 bits, compete with four position layouts.

    An optional FULL-block map separates complete blocks from residual positions.
    Decoding returns residual positions and FULL block indexes separately; the
    deserializer applies both to the destination using its normal OR semantics.
    The shared selector compares this codec with opcode 108
    (bm::set_sblock_bienc_wave_runs) using complete padded record sizes.
*/

#include "bmfunc.h"
#include "encoding.h"
#include "bmalloc.h"
#include "bmbuffer.h"
#include "bmsblockwave.h"

namespace bm
{

/** @brief Allocator-aware hierarchical wave workspace for encoding and decoding.
    @internal @ingroup bvserial
    @tparam Alloc Serializer allocator propagated to retained heap buffers.

    build() resets source coordinates and the FULL map; pairs() resets pair
    groups and dictionary state for prepared wave values. Capacity is retained.
    Only dictionary and frequency entries below entries are valid. Decoding uses
    caller-owned output and map buffers rather than this instance's source state.
    The serializer owns this model directly for its entire lifetime, including
    across vectors. Buffers grow lazily when capacity is insufficient; resetting
    or processing smaller inputs does not release their allocations. The
    deserializer likewise retains its separate map and position buffers.
*/
template<class Alloc>
struct sblock_hier_model
{
    bm::heap_vector<unsigned, Alloc, true> full_blocks; ///< Original FULL block indexes, if present.
    bm::heap_vector<unsigned, Alloc, true> blocks; ///< Original occupied block indexes.
    bm::heap_vector<unsigned, Alloc, true> source; ///< Block-compacted positions.
    bm::heap_vector<unsigned, Alloc, true> complement; ///< Complement map scratch.
    bm::heap_vector<unsigned, Alloc, true> runs; ///< Triples of run, gap and repeat.
    bm::heap_vector<bm::id64_t, Alloc, true> keys; ///< Sorted pair keys for frequency counting.
    bm::sblock_wave_adaptive_model<Alloc> wave; ///< Wave transform workspace.
    bm::id64_t dictionary[8]{}; ///< Pair dictionary in frequency rank order.
    unsigned frequencies[8]{}; ///< Frequencies used during dictionary selection.
    unsigned entries = 0; ///< Actual dictionary count.
    unsigned final_run = 0; ///< Terminal contiguous run length.

    /** Rebuild residual block coordinates and reset the FULL map.
        @param ids Borrowed array of strictly increasing superblock-relative IDs.
        @param n Number of residual positions, at least two.
        @pre All IDs are below 2^24 and exclude FULL blocks.
        @post blocks and source describe the residual input; wave.build() resets
            the wave state. The caller must prepare the selected wave width.
        @note For encoding, n is at most 65535 without FULL blocks, or 131071
            with a FULL map. Populate full_blocks after build() when needed;
            maps must be disjoint and total cardinality at most 196607.
    */
    void build(const unsigned* ids, unsigned n);

    /** Group run/gap pairs and rank up to eight dictionary entries.
        @pre wave.values contains nonempty, strictly increasing corrected IDs.
        @post runs, keys, dictionary, frequencies, entries and final_run describe
            the current values; prior pair state is discarded, retaining capacity.
        @note Frequency counts groups, not their expanded repeat counts.
    */
    void pairs();

    /** Write or dry-count a full-domain occupancy map.
        @tparam Out Shared bit accumulator or exact bit counter.
        @tparam Vect Container of unsigned coordinates.
        @param out Receives map fields without flushing or final padding.
        @param selected Borrowed nonempty, strictly increasing coordinates.
        @param domain Coordinate count, at least selected.size().
        @pre Every selected coordinate is below domain.
        @return 0=all occupied, 1=gamma runs, 2=BIC occupied, 3=BIC empty.
        @note Reuses complement scratch. Gamma wins equal bit costs; occupied
            BIC wins a tie with empty BIC. BIC stores count and polarity; gamma
            runs terminate at domain and derive the count during decoding.
    */
    template<class Out, class Vect>
    unsigned map(Out& out, const Vect& selected, unsigned domain);

    /** Encode the current record into one shared bit accumulator.
        @tparam Out Shared bit accumulator or exact bit counter.
        @param out Receives fields after the opcode, excluding final padding.
        @param sb Absolute superblock number.
        @param width Width code from 0 through 5, selecting 128 << width bits.
        @param mode 0=bounded-domain BIC, 1=BIC with explicit endpoints,
            2=run pairs, 3=run pairs with at most eight dictionary entries.
        @pre build(), wave.prepare() and pairs() describe the current input;
            width matches the prepared wave size. Mode 3 requires entries > 0.
        @note The caller writes opcode 109
            (bm::set_sblock_bienc_wave_hierarchical) and final accumulator padding.
    */
    template<class Out>
    void encode(Out& out, unsigned sb, unsigned width, unsigned mode);

    /** Write the current corrected positions without maps or headers.
        @tparam Out Shared bit accumulator or exact bit counter.
        @param out Receives payload fields without final padding.
        @param mode Position layout, 0 through 3 as described by encode().
        @pre The same prepared values and pair state required by encode().
    */
    template<class Out>
    void positions(Out& out, unsigned mode);

    /** Decode and validate a nonempty occupancy map.
        @tparam In Shared bit input accumulator.
        @tparam Vect Resizable container of unsigned coordinates.
        @param in Input positioned at the map selector; advances past the map.
        @param dst Receives sorted occupied coordinates, replacing prior contents.
        @param domain Known, nonzero coordinate count.
        @note Inverse BIC expands within dst, retaining capacity for the full
            domain. No temporary heap buffer is created. Invalid map bounds
            or ordering raise a serialization error.
    */
    template<class In, class Vect>
    static void read_map(In& in, Vect& dst, unsigned domain);

    /** Decode and validate an opcode 109
        (bm::set_sblock_bienc_wave_hierarchical) payload.
        @tparam In Shared bit input accumulator.
        @tparam Vect Container of unsigned map coordinates.
        @param in Input positioned after the opcode byte; advances past payload.
        @param output Caller-owned storage for at least 65535 residual positions
            without FULL blocks, or 131071 when the FULL map is present.
        @param sb Receives the absolute superblock number.
        @param blocks Receives sorted original residual block indexes.
        @param waves Receives occupied wave indexes in compacted block space.
        @param full Receives original FULL block indexes, empty if absent.
        @pre Output storage and the three map buffers do not alias each other.
        @return Number of strictly increasing residual positions written.
        @note FULL ranges are not expanded into output; the caller restores
            them separately. Invalid fields raise a serialization error and
            may leave output and scratch partially modified.
    */
    template<class In, class Vect>
    static unsigned decode(In& in, unsigned* output, unsigned& sb,
                           Vect& blocks, Vect& waves, Vect& full);
};

template<class Alloc>
void sblock_hier_model<Alloc>::build(const unsigned* ids, unsigned n)
{
    full_blocks.clear();
    blocks.clear();
    source.resize(n);
    // Replace each original block index by its rank among occupied blocks.
    // The low 16 bits stay unchanged, preserving offsets within each block.
    unsigned prior = ~0u;
    for (unsigned i = 0; i < n; ++i)
    {
        unsigned b = ids[i] >> 16;
        if (b != prior)
        {
            blocks.push_back(b);
        }
        prior = b;
        source[i] = unsigned(blocks.size() - 1) * 65536 + (ids[i] & 65535);
    } // for i
    wave.build(source.data(), n);
}

template<class Alloc>
void sblock_hier_model<Alloc>::pairs()
{
    runs.clear();
    keys.clear();
    entries = 0;
    // Each triple is (contiguous set run, following zero gap, repeat count).
    // Consecutive identical pairs share one triple; the final run has no gap.
    unsigned ones = 1;
    for (unsigned i = 1; i < wave.values.size(); ++i)
    {
        unsigned gap = wave.values[i] - wave.values[i - 1] - 1;
        if (!gap)
        {
            ++ones;
            continue;
        }
        unsigned len = unsigned(runs.size());
        if (len && runs[len - 3] == ones && runs[len - 2] == gap)
        {
            ++runs[len - 1];
        }
        else
        {
            runs.push_back(ones);
            runs.push_back(gap);
            runs.push_back(1);
        }
        ones = 1;
    } // for i
    final_run = ones;
    for (unsigned i = 0; i < runs.size(); i += 3)
    {
        keys.push_back((bm::id64_t(runs[i]) << 32) | runs[i + 1]);
    } // for i
    // Count one key per group, regardless of repeat count: each group emits
    // one dictionary reference. Sort to collect equal keys in a single scan.
    bm::heap_sort(keys.data(), keys.size());
    for (unsigned i = 0; i < keys.size();)
    {
        unsigned end = i + 1;
        while (end < keys.size() && keys[end] == keys[i])
        {
            ++end;
        } // while end
        // Keep the eight most frequent keys. Smaller keys break frequency
        // ties, making dictionary order and encoded references deterministic.
        unsigned freq = end - i, rank = 0;
        while (rank < entries &&
               (frequencies[rank] > freq ||
                (frequencies[rank] == freq && dictionary[rank] < keys[i])))
        {
            ++rank;
        } // while rank
        if (rank < 8)
        {
            unsigned last = entries < 8 ? entries++ : 7;
            for (unsigned j = last; j > rank; --j)
            {
                frequencies[j] = frequencies[j - 1];
                dictionary[j] = dictionary[j - 1];
            } // for j
            dictionary[rank] = keys[i];
            frequencies[rank] = freq;
        }
        i = end;
    } // for i
}

template<class Alloc>
template<class Out, class Vect>
unsigned sblock_hier_model<Alloc>::map(Out& out, const Vect& selected, unsigned domain)
{
    // A fully occupied domain needs only its selector bit.
    bool all = selected.size() == domain;
    out.put_bits(all, 1);
    if (all)
    {
        return 0;
    }
    complement.clear();
    for (unsigned i = 0, j = 0; i < domain; ++i)
    {
        if (j < selected.size() && selected[j] == i)
        {
            ++j;
        }
        else
        {
            complement.push_back(i);
        }
    } // for i
    // Common selector bits cancel. Include the polarity bit and cardinality
    // for each BIC alternative; gamma runs carry their own starting polarity.
    bm::sblock_wave_bit_counter run, a, b;
    bm::sblock_wave_gamma_map(run, selected, domain);
    a.gamma(unsigned(selected.size()));
    a.put_bits(0, 1);
    a.bic_encode_u32_cm(selected.data(), unsigned(selected.size()), 0, domain - 1);
    b.gamma(unsigned(complement.size()));
    b.put_bits(0, 1);
    b.bic_encode_u32_cm(complement.data(), unsigned(complement.size()), 0, domain - 1);
    bool gamma = run.bits <= a.bits && run.bits <= b.bits;
    out.put_bits(gamma, 1);
    if (gamma)
    {
        bm::sblock_wave_gamma_map(out, selected, domain);
        return 1;
    }
    else
    {
        bool inverse = b.bits < a.bits;
        out.put_bits(inverse, 1);
        if (inverse)
        {
            out.gamma(unsigned(complement.size()));
            out.bic_encode_u32_cm(complement.data(), unsigned(complement.size()), 0, domain - 1);
            return 3;
        }
        else
        {
            out.gamma(unsigned(selected.size()));
            out.bic_encode_u32_cm(selected.data(), unsigned(selected.size()), 0, domain - 1);
            return 2;
        }
    }
}

template<class Alloc>
template<class Out>
void sblock_hier_model<Alloc>::encode(Out& out, unsigned sb, unsigned width, unsigned mode)
{
    // Fixed prefix: width, reserved-zero extension, superblock index layout,
    // position layout and optional FULL-map flag. All fields share out.
    out.put_bits(width, 3);
    out.put_bits(0, 1);
    // Gamma codes small indexes; larger indexes use 8, 16 or 32 bits.
    unsigned sm = sb < 8 ? 0 : sb <= 255 ? 1 : sb <= 65535 ? 2 : 3;
    out.put_bits(sm, 2);
    out.put_bits(mode, 2);
    out.put_bits(!full_blocks.empty(), 1);
    if (!sm)
    {
        out.gamma(sb + 1);
    }
    else
    {
        out.put_bits(sb, sm == 1 ? 8 : sm == 2 ? 16 : 32);
    }
    // Count only residual positions. FULL ranges have their own root map.
    out.gamma(unsigned(source.size()));
    if (!full_blocks.empty())
    {
        map(out, full_blocks, 256);
    }
    map(out, blocks, 256);
    out.gamma(wave.correction + 1);
    map(out, wave.occupied, wave.wave_count);
    positions(out, mode);
}

template<class Alloc>
template<class Out>
void sblock_hier_model<Alloc>::positions(Out& out, unsigned mode)
{
    if (mode == 0)
    {
        out.bic_encode_u32_cm(wave.values.data(), unsigned(source.size()), 0, wave.upper);
    }
    else if (mode == 1)
    {
        // Explicit endpoints narrow the BIC domain for interior positions.
        out.gamma(wave.values.front() + 1);
        out.gamma(wave.upper - wave.values.back() + 1);
        if (source.size() > 2)
        {
            out.bic_encode_u32_cm(wave.values.data() + 1,
                                 unsigned(source.size() - 2),
                                 wave.values.front() + 1, wave.values.back() - 1);
        }
    }
    else if (mode == 2)
    {
        bm::sblock_wave_pairs(out, wave.values);
    }
    else
    {
        // Dictionary entries store (run length, gap length). A miss writes
        // the literal pair; both paths then emit the group repeat count.
        out.gamma(entries);
        for (unsigned i = 0; i < entries; ++i)
        {
            out.gamma(unsigned(dictionary[i] >> 32));
            out.gamma(unsigned(dictionary[i]));
        } // for i
        out.gamma(wave.values.front() + 1);
        out.gamma(unsigned(runs.size() / 3) + 1);
        for (unsigned i = 0; i < runs.size(); i += 3)
        {
            bm::id64_t key = (bm::id64_t(runs[i]) << 32) | runs[i + 1];
            unsigned rank = 0;
            while (rank < entries && dictionary[rank] != key)
            {
                ++rank;
            } // while rank
            out.put_bits(rank < entries, 1);
            if (rank < entries)
            {
                out.gamma(rank + 1);
            }
            else
            {
                out.gamma(runs[i]);
                out.gamma(runs[i + 1]);
            }
            out.gamma(runs[i + 2]);
        } // for i
        out.gamma(final_run);
    }
}

template<class Alloc>
template<class In, class Vect>
void sblock_hier_model<Alloc>::read_map(In& in, Vect& dst, unsigned domain)
{
    dst.clear();
    if (in.get_bits(1))
    {
        for (unsigned i = 0; i < domain; ++i)
        {
            dst.push_back(i);
        } // for i
        return;
    }
    if (in.get_bits(1))
    {
        // Alternating zero/one runs must cover exactly the known domain.
        unsigned bit = in.get_bits(1), pos = 0;
        while (pos < domain)
        {
            unsigned run = in.gamma();
            bm::serialization_require(run && run <= domain - pos);
            if (bit)
            {
                for (unsigned i = 0; i < run; ++i)
                {
                    dst.push_back(pos + i);
                } // for i
            }
            pos += run;
            bit ^= 1;
        } // while pos
    }
    else
    {
        bool inv = in.get_bits(1);
        unsigned n = in.gamma();
        bm::serialization_require(n && n < domain);
        dst.resize(n);
        in.bic_decode_u32_cm(dst.data(), n, 0, domain - 1);
        for (unsigned i = 0; i < n; ++i)
        {
            bm::serialization_require(dst[i] < domain && (!i || dst[i] > dst[i - 1]));
        } // for i
        if (inv)
        {
            // Keep empty coordinates in the prefix while writing occupied
            // coordinates into the suffix, from right to left. The write
            // cursor never crosses unread empty coordinates. Retain the full
            // domain capacity for later maps instead of swapping in a temporary.
            dst.resize(domain);
            unsigned remaining = n;
            unsigned write = domain;
            for (unsigned i = domain; i; --i)
            {
                if (remaining && dst[remaining - 1] == i - 1)
                {
                    --remaining;
                }
                else
                {
                    BM_ASSERT(write > remaining);
                    dst[--write] = i - 1;
                }
            } // for i
            BM_ASSERT(!remaining && write == n);
            // The occupied suffix starts at n. Move it to the prefix in
            // forward order so overlapping source and destination are safe.
            const unsigned occupied = domain - n;
            for (unsigned i = 0; i < occupied; ++i)
            {
                dst[i] = dst[n + i];
            } // for i
            dst.resize(occupied);
        }
    }
    bm::serialization_require(!dst.empty());
}

template<class Alloc>
template<class In, class Vect>
unsigned sblock_hier_model<Alloc>::decode(In& in, unsigned* output, unsigned& sb,
                                        Vect& blocks, Vect& waves, Vect& full)
{
    // Reject unsupported widths/extensions before using them in domain sizes.
    unsigned width = in.get_bits(3);
    bm::serialization_require(width < 6 && !in.get_bits(1));
    unsigned sm = in.get_bits(2), mode = in.get_bits(2);
    const bool has_full = bool(in.get_bits(1));
    sb = sm ? in.get_bits(sm == 1 ? 8 : sm == 2 ? 16 : 32) : in.gamma() - 1;
    unsigned n = in.gamma();
    bm::serialization_require(n >= 2 && n <= (has_full ? 131071u : 65535u));
    full.clear();
    if (has_full)
    {
        read_map(in, full, 256);
    }
    read_map(in, blocks, 256);
    bm::serialization_require(n + full.size() * 65536u <= 196607u);
    // The two root maps must be disjoint. Both are sorted, so advance the
    // smaller index until one map is exhausted.
    for (unsigned i = 0, j = 0; i < blocks.size() && j < full.size();)
    {
        bm::serialization_require(blocks[i] != full[j]);
        if (blocks[i] < full[j])
        {
            ++i;
        }
        else
        {
            ++j;
        }
    } // for i
    unsigned d = in.gamma() - 1, wsize = 128u << width;
    read_map(in, waves, unsigned(blocks.size()) * (65536 / wsize));
    // Gap correction subtracts i*d from the ith packed position. Validate
    // its total span in 64 bits before constructing the corrected BIC bound.
    bm::id64_t limit = bm::id64_t(waves.size()) * wsize;
    bm::id64_t delta = bm::id64_t(n - 1) * d;
    bm::serialization_require(delta < limit && n <= limit);
    unsigned upper = unsigned(limit - 1 - delta);
    if (mode < 2)
    {
        if (!mode)
        {
            in.bic_decode_u32_cm(output, n, 0, upper);
        }
        else
        {
            unsigned first = in.gamma(), tail = in.gamma();
            bm::serialization_require(first && tail && first <= upper + 1 && tail <= upper + 1);
            output[0] = first - 1;
            output[n - 1] = upper - tail + 1;
            bm::serialization_require(output[0] < output[n - 1] &&
                                      output[n - 1] - output[0] >= n - 1);
            if (n > 2)
            {
                in.bic_decode_u32_cm(output + 1, n - 2, output[0] + 1, output[n - 1] - 1);
            }
        }
    }
    else
    {
        bm::id64_t dict[8] = {};
        unsigned count = 0;
        if (mode == 3)
        {
            count = in.gamma();
            bm::serialization_require(count && count <= 8);
            for (unsigned i = 0; i < count; ++i)
            {
                unsigned a = in.gamma(), b = in.gamma();
                dict[i] = (bm::id64_t(a) << 32) | b;
            } // for i
        }
        // groups includes the terminal run, which is emitted without a gap.
        unsigned first = in.gamma(), groups = in.gamma(), used = 0;
        bm::serialization_require(first && first <= upper + 1 && groups && groups <= n);
        unsigned pos = first - 1;
        // Append a set run only after checking both output capacity and the
        // corrected coordinate domain. pos then points just past that run.
        auto append = [&](unsigned ones)
        {
            bm::serialization_require(ones && ones <= n - used &&
                                      pos <= upper && ones - 1 <= upper - pos);
            for (unsigned j = 0; j < ones; ++j)
            {
                output[used++] = pos++;
            } // for j
        };
        for (unsigned i = 1; i < groups; ++i)
        {
            unsigned ones, gap;
            if (mode == 3 && in.get_bits(1))
            {
                unsigned rank = in.gamma();
                bm::serialization_require(rank && rank <= count);
                ones = unsigned(dict[rank - 1] >> 32);
                gap = unsigned(dict[rank - 1]);
            }
            else
            {
                ones = in.gamma();
                gap = in.gamma();
            }
            unsigned repeat = in.gamma();
            bm::serialization_require(ones && gap && repeat && repeat <= n - used);
            for (unsigned j = 0; j < repeat; ++j)
            {
                append(ones);
                bm::serialization_require(pos <= upper && gap <= upper - pos);
                pos += gap;
            } // for j
        } // for i
        append(in.gamma());
        bm::serialization_require(used == n);
    }
    // Undo transforms in reverse: restore i*d, expand occupied-wave ranks,
    // then expand occupied-block ranks. Validate corrected values before
    // overwriting them; FULL blocks remain in the separate full map.
    unsigned prior = 0;
    for (unsigned i = 0; i < n; ++i)
    {
        bm::serialization_require(output[i] <= upper && (!i || output[i] > prior));
        prior = output[i];
        unsigned packed = unsigned(output[i] + bm::id64_t(i) * d);
        unsigned compact = waves[packed / wsize] * wsize + packed % wsize;
        output[i] = blocks[compact >> 16] * 65536 + (compact & 65535);
        bm::serialization_require(!i || output[i] > output[i - 1]);
    } // for i
    return n;
}


/** @brief Exact unpadded component costs of one complete wave record.
    @internal @ingroup bvserial
*/
struct sblock_wave_score
{
    bm::id64_t bytes = ~bm::id64_t(0); ///< Complete bytes, or unavailable sentinel.
    bm::id64_t bits = 0; ///< All fields excluding opcode and final padding.
    bm::id64_t map_bits = 0; ///< Occupancy maps including count and selectors.
    bm::id64_t payload_bits = 0; ///< Position payload including explicit endpoints.
    unsigned width = 0; ///< Width code, 128 shifted by this value.
    unsigned mode = 0; ///< Format-specific payload/map flags.
    unsigned occupied = 0; ///< Occupied waves at this width.
    unsigned map_mode = 0; ///< Hierarchical root coder times four plus wave coder; map() codes 0 through 3.

    /** Retain a strict complete-byte improvement, preserving the first tie.
        @param total All fields excluding opcode and padding.
        @param maps Occupancy map fields included in total.
        @param payload Position payload fields included in total.
        @param w Width code.
        @param m Format-specific mode.
        @post On improvement, replaces cost, width and mode fields; occupied
            and map_mode remain the caller's responsibility.
    */
    void consider(bm::id64_t total, bm::id64_t maps, bm::id64_t payload,
                  unsigned w, unsigned m);
};

inline
void sblock_wave_score::consider(bm::id64_t total, bm::id64_t maps,
                                 bm::id64_t payload, unsigned w, unsigned m)
{
    // One opcode byte plus the accumulator rounded up to complete 32-bit words.
    // Compare final bytes: saving unpadded bits alone need not shrink a record.
    const bm::id64_t size = 1 + 4 * ((total + 31) / 32);
    if (size < bytes)
    {
        bytes = size;
        bits = total;
        map_bits = maps;
        payload_bits = payload;
        width = w;
        mode = m;
    }
}

/** @brief Optional search observations, independent of emitted opcode statistics.
    @internal @ingroup bvserial
    Counts trials before the higher-count selector can replace their records.
*/
struct sblock_wave_score_observation
{
    unsigned n = 0; ///< Original set-bit count.
    unsigned blocks = 0; ///< Occupied source blocks; FULL blocks cannot occur.
    unsigned sb = 0; ///< Absolute superblock index.
    bm::id64_t incumbent = 0; ///< Best ordinary/adaptive size to beat.
    bm::sblock_wave_score flat[6]; ///< Best opcode 108 (bm::set_sblock_bienc_wave_runs) layout per width.
    /// Best opcode 109 (bm::set_sblock_bienc_wave_hierarchical) layout per width.
    bm::sblock_wave_score hier[6];
    bm::sblock_wave_score best_flat; ///< Best flat layout across widths.
    bm::sblock_wave_score best_hier; ///< Best hierarchical layout across widths.
    unsigned tree_leaf = 0; ///< Cheap selector leaf, numbered 0 through 6.
    /// Predicted opcode 108 (bm::set_sblock_bienc_wave_runs) or
    /// opcode 109 (bm::set_sblock_bienc_wave_hierarchical).
    unsigned tree_format = 108;
    bm::id64_t tree_regret = 0; ///< Exact bytes lost after incumbent fallback.
};

/** @brief Source metadata for the inexpensive wave-format decision tree.
    @internal @ingroup bvserial
*/
struct sblock_wave_tree_features
{
    unsigned gap_level0 = 0; ///< Number of level-zero GAP blocks.
    unsigned gap_words = 0; ///< Sum of GAP header lengths.
    unsigned occupied_runs = 0; ///< Runs of nonempty source blocks.
    unsigned occupied_span = 0; ///< Inclusive span of occupied source blocks.
};

/** @brief Select a wave format using only source metadata.
    @internal @ingroup bvserial
    @param features Source block survey and occupancy geometry.
    @param n Number of residual set positions.
    @param leaf Receives the stable leaf number for regret calibration.
    @return Opcode 108 (bm::set_sblock_bienc_wave_runs) or opcode 109
        (bm::set_sblock_bienc_wave_hierarchical); prediction is heuristic.
*/
inline
unsigned sblock_wave_tree_predict(const bm::sblock_wave_tree_features& features,
                                  unsigned n, unsigned& leaf)
{
    // Fixed source-metadata thresholds avoid preparing both formats. Leaves
    // remain stable for observation/regret reporting; this is not a size proof.
    if (features.gap_level0 < 8)
    {
        if (features.gap_words < 64)
        {
            leaf = n >= 512 ? 1 : 0;
        }
        else
        {
            leaf = n >= 64 ? 3 : 2;
        }
    }
    else
    {
        leaf = features.occupied_runs < 4 ? 4 :
               features.occupied_span < 32 ? 5 : 6;
    }
    return leaf == 1 || leaf == 3 || leaf == 4 || leaf == 5 ? 109 : 108;
}

/** @brief Share wave transforms and exact position costs for opcode 108
    (bm::set_sblock_bienc_wave_runs) and opcode 109
    (bm::set_sblock_bienc_wave_hierarchical).
    @internal @ingroup bvserial
    @tparam Alloc Serializer allocator.
    @param flat Retained primary workspace; rebuilt in original coordinates
        when both formats participate, otherwise in the chosen coordinate space.
    @param hierarchy Retained compacted-block workspace.
    @param original Borrowed sorted source IDs.
    @param n Source count, 2 through 65535.
    @param sb Absolute superblock index.
    @param rival Best non-wave record size.
    @param policy 1=both formats, 2=hierarchical only, 3=hierarchical only
        when at most 16 blocks are occupied (otherwise both), 4=tree choice;
        internal value 5=flat only.
    @param features Optional borrowed source metadata for prediction and
        observations; required for policy 4.
    @param hierarchical Receives whether 109 wins the final competition.
    @param observer Optional synchronous callback; null disables reporting.
        Its score argument is borrowed and valid only during the callback.
    @param context Borrowed callback state valid throughout this call.
    @pre original points to n strictly increasing IDs below 2^24, without
        FULL blocks; flat and hierarchy are distinct retained workspaces.
    @post Trial state is scratch, not a prepared winning record. The caller
        must prepare the chosen workspace at the returned width before encoding.
    @return Winning width, mode and bytes; ties retain the incumbent, then
        opcode 108 (bm::set_sblock_bienc_wave_runs). Width and mode are relevant
        only when bytes is strictly less than rival.
*/
template<class Alloc>
bm::sblock_wave_adaptive_choice sblock_wave_shared_select(
    bm::sblock_wave_adaptive_model<Alloc>& flat,
    bm::sblock_hier_model<Alloc>& hierarchy,
    const unsigned* original, unsigned n, unsigned sb, bm::id64_t rival,
    unsigned policy, bool& hierarchical,
    void (*observer)(void*, const bm::sblock_wave_score_observation&),
    void* context, const bm::sblock_wave_tree_features* features = nullptr)
{
    // Normalize policy once. Internal policy 5 means flat-only; policy 2
    // prepares compacted coordinates directly and avoids a second transform.
    hierarchy.build(original, n);
    unsigned tree_leaf = 0;
    const unsigned tree_format = features ?
        bm::sblock_wave_tree_predict(*features, n, tree_leaf) : 108;
    if (policy == 3)
    {
        policy = hierarchy.blocks.size() <= 16 ? 2u : 1u;
    }
    if (policy == 4)
    {
        bm::serialization_require(features);
        policy = tree_format == 109 ? 2u : 5u;
    }
    const unsigned* source = policy == 2 ? hierarchy.source.data() : original;
    flat.build(source, n);
    // Scores are fixed-size local state. Detailed observations escape only
    // through the optional synchronous callback below.
    bm::sblock_wave_score_observation scores;
    scores.n = n;
    scores.blocks = unsigned(hierarchy.blocks.size());
    scores.sb = sb;
    scores.incumbent = rival;
    scores.tree_leaf = tree_leaf;
    scores.tree_format = tree_format;
    bm::sblock_wave_bit_counter root;
    const unsigned root_mode = policy == 5 ? 0 :
        hierarchy.map(root, hierarchy.blocks, 256);
    bm::sblock_wave_bit_counter common;
    const unsigned sm = sb < 8 ? 0 : sb <= 255 ? 1 : sb <= 65535 ? 2 : 3;
    if (!sm)
    {
        common.gamma(sb + 1);
    }
    else
    {
        common.put_bits(sb, sm == 1 ? 8 : sm == 2 ? 16 : 32);
    }
    common.gamma(n);
    // At each width, removing empty blocks preserves occupied-wave order
    // and offsets. Both formats therefore share correction and payload costs;
    // only occupancy maps and fixed header bits differ.
    for (unsigned w = 0; w < 6; ++w)
    {
        flat.wave_size = 128u << w;
        flat.wave_count = (policy == 2 ?
            unsigned(hierarchy.blocks.size()) * 65536 : (1u << 24)) / flat.wave_size;
        bm::sblock_wave_record record;
        flat.prepare(source, n, record, true);
        auto& wave = hierarchy.wave;
        wave.wave_size = flat.wave_size;
        const unsigned per_block = 65536 / wave.wave_size;
        wave.wave_count = scores.blocks * per_block;
        wave.correction = flat.correction;
        wave.upper = flat.upper;
        wave.occupied.clear();
        unsigned rank = 0;
        for (unsigned index : flat.occupied)
        {
            if (policy == 2)
            {
                wave.occupied.push_back(index);
                continue;
            }
            const unsigned block = index / per_block;
            // Merge sorted wave block IDs with occupied source blocks. rank
            // only advances, and every wave belongs to one of these blocks.
            while (hierarchy.blocks[rank] != block)
            {
                ++rank;
            } // while rank
            wave.occupied.push_back(rank * per_block + index % per_block);
        } // for index
        // The two transformations preserve occupied-wave order and offsets.
        // Borrow the one computed value array throughout this width trial.
        wave.values.swap(flat.values);
        if (policy != 5)
        {
            hierarchy.pairs();
        }
        // Count each position layout once, including endpoint/dictionary
        // overhead. These are unpadded bits, rounded only for complete records.
        bm::id64_t payload[4] = {};
        bm::sblock_wave_bit_counter counter;
        counter.bic_encode_u32_cm(wave.values.data(), n, 0, wave.upper);
        payload[0] = counter.bits;
        counter.bits = 0;
        counter.gamma(wave.values.front() + 1);
        counter.gamma(wave.upper - wave.values.back() + 1);
        if (n > 2)
        {
            counter.bic_encode_u32_cm(wave.values.data() + 1, n - 2,
                                     wave.values.front() + 1, wave.values.back() - 1);
        }
        payload[1] = counter.bits;
        counter.bits = 0;
        bm::sblock_wave_pairs(counter, wave.values);
        payload[2] = counter.bits;
        if (policy != 5 && hierarchy.entries)
        {
            counter.bits = 0;
            hierarchy.positions(counter, 3);
            payload[3] = counter.bits;
        }
        counter.bits = common.bits;
        counter.gamma(wave.correction + 1);
        const bm::id64_t prefix = counter.bits;
        counter.bits = 0;
        const unsigned wave_mode = policy == 5 ? 0 :
            hierarchy.map(counter, wave.occupied, wave.wave_count);
        const bm::id64_t hierarchy_maps = root.bits + counter.bits;
        // Hierarchical fixed prefix is 9 bits; root and wave maps each carry
        // their own selectors. Dictionary mode is invalid without entries.
        for (unsigned mode = 0; policy != 5 && mode < 4; ++mode)
        {
            if (mode == 3 && !hierarchy.entries)
            {
                continue;
            }
            scores.hier[w].consider(9 + prefix + hierarchy_maps + payload[mode],
                                   hierarchy_maps, payload[mode], w, mode);
        } // for mode
        if (policy != 2 && flat.occupied.size() < flat.wave_count)
        {
            bm::id64_t bic[2], gamma[2], counts[2];
            for (unsigned polarity = 0; polarity < 2; ++polarity)
            {
                const auto& map = polarity ? flat.empty : flat.occupied;
                counter.bits = 0;
                counter.gamma(unsigned(map.size()));
                counts[polarity] = counter.bits;
                counter.bits = 0;
                counter.bic_encode_u32_cm(map.data(), unsigned(map.size()), 0,
                                         flat.wave_count - 1);
                bic[polarity] = counter.bits;
                counter.bits = 0;
                bm::sblock_wave_gamma_map(counter, map, flat.wave_count);
                gamma[polarity] = counter.bits;
            } // for polarity
            // Flat mode bits: bit 0=explicit endpoints, bit 1=inverse map,
            // bit 2=gamma map, bit 3=run-pair payload. Endpoints and pairs cannot
            // coexist. The flat fixed prefix contributes 10 bits.
            for (unsigned mode = 0; mode < 16; ++mode)
            {
                if ((mode & 8) && (mode & 1))
                {
                    continue;
                }
                const unsigned polarity = (mode >> 1) & 1;
                const bm::id64_t maps = counts[polarity] +
                    ((mode & 4) ? gamma[polarity] : bic[polarity]);
                const bm::id64_t values = payload[(mode & 8) ? 2 : (mode & 1)];
                scores.flat[w].consider(10 + prefix + maps + values, maps, values, w, mode);
            } // for mode
        }
        scores.flat[w].occupied = unsigned(flat.occupied.size());
        scores.hier[w].occupied = unsigned(flat.occupied.size());
        scores.hier[w].map_mode = root_mode * 4 + wave_mode;
        if (scores.flat[w].bytes < scores.best_flat.bytes)
        {
            scores.best_flat = scores.flat[w];
        }
        if (scores.hier[w].bytes < scores.best_hier.bytes)
        {
            scores.best_hier = scores.hier[w];
        }
        // Return borrowed values to their owning workspace before reuse.
        wave.values.swap(flat.values);
    } // for w
    // Exact dual-format trials can measure the heuristic loss after applying
    // the same ordinary-record fallback. Other policies do not score both.
    if (features && policy == 1)
    {
        const bm::id64_t predicted = tree_format == 109 ?
            scores.best_hier.bytes : scores.best_flat.bytes;
        const bm::id64_t exact = bm::min_value(scores.best_hier.bytes,
                                              scores.best_flat.bytes);
        scores.tree_regret = bm::min_value(rival, predicted) -
                             bm::min_value(rival, exact);
    }
    if (observer)
    {
        observer(context, scores);
    }
    if (policy == 5)
    {
        scores.best_hier.bytes = ~bm::id64_t(0);
    }
    // Strict comparisons preserve ordinary encoding on a tie, then flat over
    // hierarchical. Per-format ties retain the first width/mode visited.
    const auto& best = scores.best_hier.bytes < scores.best_flat.bytes ?
        scores.best_hier : scores.best_flat;
    bm::sblock_wave_adaptive_choice choice;
    choice.bytes = rival;
    hierarchical = false;
    if (best.bytes < rival)
    {
        choice.bytes = best.bytes;
        choice.width = best.width;
        choice.mode = best.mode;
        hierarchical = scores.best_hier.bytes < scores.best_flat.bytes;
    }
    return choice;
}

} // namespace bm
#endif
