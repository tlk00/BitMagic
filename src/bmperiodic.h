/*! \file bmperiodic.h
    Internal arithmetic-progression prediction for bvector BIT-block serialization.
    Supports independent and ordered multi-group searches and formula codecs.
    Include after bm.h or bm64.h, as for other serialization utilities.
    \ingroup bvserial
    \internal
    The residual is the XOR of the normalized block and all recorded chains.
*/
#ifndef BM_PERIODIC_H_INCLUDED
#define BM_PERIODIC_H_INCLUDED

#include "bmfunc.h"
#include "bmbuffer.h"
#include "encoding.h"

namespace bm
{

/** \brief Bounds and policy for periodic block search.
    \ingroup bvserial
    \internal
*/
struct periodic_options
{
    unsigned min_step = 1; ///< first arithmetic-progression step
    unsigned max_step = 16; ///< last arithmetic-progression step
    unsigned max_variants = 4; ///< actual residual serialization trials
    unsigned min_length = 4; ///< Minimum chain length.
    unsigned max_holes = 1; ///< 0 or 1 missing bit per chain
    unsigned max_chains = 256; ///< bounds search and formula size
    unsigned sample = 1; ///< deterministic 1-in-N BIT-block sample
    unsigned order = 0; ///< 0 legacy; 1 independent; 2 dense; 3 sparse; 4 middle
    unsigned trial_budget = 2304; ///< residual encodings per sampled block
    unsigned formula_chains = 256; ///< total serialized chains, independent of pool cap
    unsigned min_step_gain = 0; ///< bytes per step before counting it as marginal
    unsigned marginal_streak = 0; ///< consecutive marginal steps before stopping; 0 disables
    unsigned min_stop_step = 1; ///< do not stop for a marginal streak before this step
    bool diagnostics = false; ///< measure candidates regardless of score
};

/** \brief Diagnostic measurements for one searched block.
    \ingroup bvserial
    \internal
*/
struct periodic_order_block
{
    bm::id64_t vector = 0; ///< Source vector index.
    unsigned block = 0; ///< Source block index.
    unsigned baseline = 0; ///< Ordinary encoded bytes.
    unsigned encoded = 0; ///< Selected encoded bytes.
    unsigned trials = 0; ///< Number of trials.
    unsigned groups = 0; ///< Number of groups.
    unsigned chains = 0; ///< Number of chains.
    unsigned popcount = 0; ///< Block population count.
    unsigned gaps = 0; ///< Block transition count.
    unsigned last_step = 0; ///< Last search step.
    unsigned stop_reason = 0; ///< Search termination reason.
    unsigned formula_bytes = 0; ///< Formula size in bytes.
    unsigned residual_bytes = 0; ///< Residual size in bytes.
    bm::id64_t step_mask = 0; ///< Mask of searched steps.
};

/** \brief Internal allocator-aware serialization utility.
    \ingroup bvserial
    \internal
    \tparam Alloc Application bvector allocator (BV::allocator_type).
        The default supports standalone codec use; serializers and decoders
        pass their application allocator explicitly.
*/
template<class Alloc = bm::standard_allocator>
struct periodic_order_statistics
{
    bm::heap_vector<periodic_order_block, Alloc, true> blocks; ///< Owned per-block diagnostics.
    bm::id64_t discovered[44] = {}; ///< Chains discovered by length bucket.
    bm::id64_t retained[44] = {}; ///< Chains retained by length bucket.
    bm::id64_t selected[44] = {}; ///< Chains selected by length bucket.
    bm::id64_t trial_counts[257] = {}; ///< Trial formulas by chain count.
    bm::id64_t selected_counts[257] = {}; ///< Selected formulas by chain count.
    bm::id64_t group_counts[129] = {}; ///< Selected formulas by group count.
    bm::id64_t step_visits[129] = {}; ///< Search visits by step.
    bm::id64_t groups_step[129] = {}; ///< Accepted groups by step.
    bm::id64_t gain_step[129] = {}; ///< Saved bytes by step.
    // Header, starts, lengths, padding: meaningful bits and terminal padding.
    bm::id64_t trial_bits[4] = {}; ///< Trial header, start, length and padding bits.
    bm::id64_t selected_bits[4] = {}; ///< Selected header, start, length and padding bits.
    // Explicit-first gamma lengths or fixed interleaving.
    bm::id64_t trial_formula_modes[2] = {}; ///< Trials by formula layout.
    bm::id64_t selected_formula_modes[2] = {}; ///< Selected formulas by layout.
    bm::id64_t exhausted = 0; ///< Searches exhausting the trial budget.
    bm::id64_t multi_group = 0; ///< Selected records with multiple groups.
    bm::id64_t additions = 0; ///< Groups accepted during search.
};

/** \brief Select a chain-length diagnostics bucket.
    \ingroup bvserial
    \internal
    \param length Chain length in bits.
    \return Bucket index.
*/
inline
unsigned periodic_length_bucket(unsigned length) BMNOEXCEPT
{
    if (length <= 32) return length;
    unsigned b = 33;
    for (unsigned hi = 64; length > hi; hi *= 2)
    { ++b; } // for hi
    return b;
}

/** \brief Internal allocator-aware serialization utility.
    \ingroup bvserial
    \internal
    \tparam Alloc Application bvector allocator (BV::allocator_type).
        The default supports standalone codec use; serializers and decoders
        pass their application allocator explicitly.
*/
template<class Alloc = bm::standard_allocator>
struct periodic_statistics
{
    periodic_order_statistics<Alloc> order_stats; ///< Ordered search statistics.
    bm::id64_t eligible = 0; ///< Blocks eligible for periodic encoding.
    bm::id64_t scanned = 0; ///< Candidate scans performed.
    bm::id64_t trials = 0; ///< Number of trials.
    bm::id64_t accepted = 0; ///< Blocks accepting a periodic record.
    bm::id64_t baseline_bytes = 0; ///< Ordinary record bytes for sampled blocks.
    bm::id64_t saved_bytes = 0; ///< Bytes saved by accepted records.
    bm::id64_t residual_bytes = 0; ///< Residual size in bytes.
    bm::id64_t formula_bytes = 0; ///< Formula size in bytes.
    bm::id64_t chains_accepted = 0; ///< chains accepted.
    bm::id64_t exceptions_accepted = 0; ///< exceptions accepted.
    bm::id64_t inverted_accepted = 0; ///< inverted accepted.
    bm::id64_t approximate_accepted = 0; ///< approximate accepted.
    bm::id64_t chain_cap = 0; ///< Scans reaching the chain cap.
    bm::id64_t no_candidates = 0; ///< Blocks without candidates.
    bm::id64_t sampled = 0; ///< Blocks selected for search.
    bm::id64_t sampled_baseline_bytes = 0; ///< sampled baseline bytes.
    bm::id64_t raw_candidates = 0; ///< Candidates discovered.
    bm::id64_t positive_candidates = 0; ///< Candidates with positive heuristic score.
    bm::id64_t candidate_blocks = 0; ///< candidate blocks.
    bm::id64_t positive_blocks = 0; ///< positive blocks.
    bm::id64_t positive_trials = 0; ///< positive trials.
    bm::id64_t nonpositive_trials = 0; ///< nonpositive trials.
    bm::id64_t positive_wins = 0; ///< positive wins.
    bm::id64_t nonpositive_wins = 0; ///< nonpositive wins.
    bm::id64_t payload_improved = 0; ///< payload improved.
    bm::id64_t near_miss_32 = 0; ///< near miss 32.
    bm::id64_t near_miss_128 = 0; ///< near miss 128.
    bm::id64_t best_miss_bytes = 0; ///< Bytes by which best candidates missed the baseline.
    bm::id64_t best_candidate_blocks = 0; ///< best candidate blocks.
    bm::id64_t best_baseline_bytes = 0; ///< best baseline bytes.
    bm::id64_t best_formula_bytes = 0; ///< best formula bytes.
    bm::id64_t best_residual_bytes = 0; ///< best residual bytes.
    bm::id64_t best_chains = 0; ///< best chains.
    bm::id64_t best_covered = 0; ///< best covered.
    bm::id64_t best_exceptions = 0; ///< best exceptions.
    bm::id64_t sampled_density[4] = {}; ///< sampled density.
    bm::id64_t candidate_density[4] = {}; ///< candidate density.
    bm::id64_t accepted_density[4] = {}; ///< accepted density.
    bm::id64_t trial_score_bucket[6] = {}; ///< trial score bucket.
    bm::id64_t win_score_bucket[6] = {}; ///< win score bucket.
    bm::id64_t trial_step[129] = {}; ///< trial step.
    bm::id64_t best_step[129] = {}; ///< best step.
    bm::id64_t accepted_step[129] = {}; ///< accepted step.
    bm::id64_t saved_step[129] = {}; ///< saved step.
};

/** \brief Select a candidate-score diagnostics bucket.
    \ingroup bvserial
    \internal
    \param score Heuristic candidate score.
    \return Bucket index.
*/
inline
unsigned periodic_score_bucket(int score) BMNOEXCEPT
{
    if (score < -256) return 0;
    if (score < -64) return 1;
    if (score <= 0) return 2;
    if (score < 64) return 3;
    if (score < 256) return 4;
    return 5;
}

/** \brief Select a block-density diagnostics bucket.
    \ingroup bvserial
    \internal
    \param count Number of set bits.
    \return Bucket index.
*/
inline
unsigned periodic_density_bucket(unsigned count) BMNOEXCEPT
{
    if (count < 4096) return 0;
    if (count < 16384) return 1;
    if (count < 24576) return 2;
    return 3;
}

// Stable under worker scheduling: the same vector/block coordinates are
// sampled whether serialization runs sequentially or through the thread pool.
/** \brief Choose a deterministic block sample independent of worker scheduling.
    \ingroup bvserial
    \internal
    \param vector_index Source vector index.
    \param block_index Block index within the vector.
    \param sample Sampling divisor; must be nonzero.
    \return True when this block is selected.
*/
inline
bool periodic_take_sample(bm::id64_t vector_index, unsigned block_index,
                                 unsigned sample) BMNOEXCEPT
{
    BM_ASSERT(sample);
    if (sample == 1) return true;
    bm::id64_t x = vector_index * 0x9e3779b97f4a7c15ULL +
                   bm::id64_t(block_index) * 0xbf58476d1ce4e5b9ULL;
    x ^= x >> 30; x *= 0xbf58476d1ce4e5b9ULL;
    x ^= x >> 27; x *= 0x94d049bb133111ebULL;
    x ^= x >> 31;
    return x % sample == 0;
}

/** \brief Arithmetic progression of predicted bit positions.
    \ingroup bvserial
    \internal
*/
struct periodic_chain
{
    bm::gap_word_t start; ///< First position in the chain.
    unsigned length; ///< Number of positions in the chain.
    unsigned char holes; ///< Missing predicted bits in the chain.
};

/** \brief Internal allocator-aware serialization utility.
    \ingroup bvserial
    \internal
    \tparam Alloc Application bvector allocator (BV::allocator_type).
        The default supports standalone codec use; serializers and decoders
        pass their application allocator explicitly.
*/
template<class Alloc = bm::standard_allocator>
struct periodic_candidate
{
    unsigned step = 0; ///< Spacing between chain positions.
    unsigned max_holes = 0; ///< Maximum missing positions per chain.
    unsigned exceptions = 0; ///< Total missing predicted positions.
    unsigned covered = 0; ///< Total predicted positions.
    bool truncated = false; ///< Whether chain count reached its cap.
    int score = 0; ///< Heuristic ranking score.
    bm::heap_vector<periodic_chain, Alloc, true> chains; ///< Owned chain descriptions.
};

/** \brief Read a bit using the selected polarity.
    \ingroup bvserial
    \internal
    \param block Source bit block.
    \param pos Position within the block.
    \param inverted Read complemented polarity when true.
    \return Normalized bit value.
*/
inline
bool periodic_bit(const bm::word_t* block, unsigned pos,
                         bool inverted) BMNOEXCEPT
{ return bool(bm::test_bit(block, pos)) != inverted; }

/** \brief XOR one predicted position in a bit block.
    \ingroup bvserial
    \internal
    \param block Writable bit block.
    \param pos Position within the block.
*/
inline
void periodic_flip(bm::word_t* block, unsigned pos) BMNOEXCEPT
{ block[pos >> 5] ^= bm::word_t(1u << (pos & 31)); }

// Enumerate disjoint maximal chains separately for every residue modulo step.
// A chain begins and ends on a one; a single internal zero is an XOR exception.
/** \brief Find bounded disjoint arithmetic-progressions in a BIT block.
    \ingroup bvserial
    \internal
    \tparam Alloc Application allocator.
    \param block Source bit block.
    \param inverted Search complemented polarity when true.
    \param step Distance between successive positions.
    \param holes Maximum missing bits in a chain.
    \param min_length Minimum chain length.
    \param max_chains Maximum retained chains; must be nonzero.
    \param out Candidate reused for the result.
    \param diagnostics Optional chain discovery counters.
    \return True if at least one chain was retained.
*/
template<class Alloc>
inline
bool periodic_find(const bm::word_t* block, bool inverted,
                          unsigned step, unsigned holes, unsigned min_length,
                          unsigned max_chains, periodic_candidate<Alloc>& out,
                          periodic_order_statistics<Alloc>* diagnostics = nullptr)
{
    BM_ASSERT(block && step && max_chains && min_length);
    out.step = out.max_holes = out.exceptions = out.covered = 0;
    out.truncated = false; out.score = 0;
    out.chains.clear();
    out.chains.reserve(max_chains);
    out.step = step; out.max_holes = holes;
    const auto better = [](const periodic_chain& a, const periodic_chain& b)
        { return a.length - a.holes > b.length - b.holes; };
    for (unsigned residue = 0; residue < step; ++residue)
    {
        unsigned pos = residue;
        while (pos < 65536)
        {
            if (!bm::periodic_bit(block, pos, inverted)) { pos += step; continue; }
            const unsigned first = pos;
            unsigned last = pos, seen_holes = 0, used_holes = 0;
            for (unsigned next = pos + step; next < 65536; next += step)
            {
                if (bm::periodic_bit(block, next, inverted))
                { last = next; used_holes = seen_holes; }
                else if (++seen_holes > holes) break;
            } // for next
            const unsigned length = (last - first) / step + 1;
            if (length >= min_length)
            {
                if (diagnostics) ++diagnostics->discovered[bm::periodic_length_bucket(length)];
                const periodic_chain chain = {
                    bm::gap_word_t(first), length,
                    static_cast<unsigned char>(used_holes)};
                if (out.chains.size() < max_chains)
                {
                    out.chains.push_back(chain);
                    out.covered += length;
                    out.exceptions += used_holes;
                    if (out.chains.size() == max_chains)
                        bm::heap_make(out.chains.data(), out.chains.size(), better);
                }
                else
                {
                    out.truncated = true;
                    if (length - used_holes >
                        out.chains.front().length - out.chains.front().holes)
                    {
                        auto& worst = out.chains.front();
                        out.covered -= worst.length;
                        out.covered += length;
                        out.exceptions -= worst.holes;
                        out.exceptions += used_holes;
                        worst = chain;
                        // Only the root changed; the child heaps remain valid.
                        bm::heap_sift_down(out.chains.data(), 0, out.chains.size(), better);
                    }
                }
                pos = last + step;
            }
            else pos = first + step;
        }
    } // for residue
    if (out.chains.empty()) return false;
    bm::heap_sort(out.chains.data(), out.chains.size(),
              [](const periodic_chain& a, const periodic_chain& b)
              { return a.start < b.start; });
    if (diagnostics)
        for (const auto& c : out.chains)
        {
            ++diagnostics->retained[bm::periodic_length_bucket(c.length)];
        } // for c
    // This is only a ranking estimate; complete encoded bytes decide acceptance.
    out.score = int(out.covered) - int(2 * out.exceptions) -
                int(16 * out.chains.size());
    return true;
}

/** \brief XOR predicted chain positions into a block.
    \ingroup bvserial
    \internal
    \tparam Alloc Application allocator.
    \param block Writable residual or normalized source block.
    \param candidate Chains and their common step.
*/
template<class Alloc>
inline
void periodic_apply(bm::word_t* block,
                           const periodic_candidate<Alloc>& candidate)
{
    for (const auto& chain : candidate.chains)
    {
        for (unsigned j = 0, p = chain.start; j < chain.length;
             ++j, p += candidate.step)
        {
            bm::periodic_flip(block, p);
        } // for j
    } // for chain
}

// Formula: polarity(8), step(8), chain count(16), BIC starts, gamma lengths.
// The bit stream is padded to a 32-bit boundary by bit_out::flush.
/** \brief Encode the single-group formula following set_block_periodic_post (75).
    \ingroup bvserial
    \internal
    \tparam ENC Encoder type.
    \tparam Alloc Application allocator.
    \param enc Destination encoder.
    \param c Candidate chains to encode.
    \param inverted Original block polarity.
*/
template<class ENC, class Alloc>
inline
void periodic_write_formula(ENC& enc, const periodic_candidate<Alloc>& c,
                                   bool inverted)
{
    enc.put_8(static_cast<unsigned char>(inverted));
    enc.put_8(static_cast<unsigned char>(c.step));
    enc.put_16(bm::short_t(c.chains.size()));
    BM_ASSERT(c.chains.size() <= 256);
    bm::gap_word_t starts[256];
    for (size_t i = 0; i < c.chains.size(); ++i)
    {
        starts[i] = c.chains[i].start;
    } // for i
    bm::bit_out<ENC> bout(enc);
    bout.bic_encode_u16_cm(starts, unsigned(c.chains.size()), 0, 65535);
    for (const auto& chain : c.chains)
    {
        bout.gamma(chain.length);
    } // for chain
    bout.flush();
}

/** \brief Read and apply set_block_periodic_post (75).
    \ingroup bvserial
    \internal
    \tparam DEC Decoder type.
    \param dec Source decoder.
    \param block SIMD-aligned writable residual, or an invalid address to consume only.
*/
template<class DEC>
inline
void periodic_read_formula(DEC& dec, bm::word_t* block)
{
    const unsigned inverted = dec.get_8();
    const unsigned step = dec.get_8();
    const unsigned count = dec.get_16();
    bm::serialization_require(inverted <= 1 && step && step <= 64 &&
                           count && count <= 256);
    bm::gap_word_t starts[256];
    bm::bit_in<DEC> bin(dec);
    bin.bic_decode_u16_cm(starts, count, 0, 65535);
    for (unsigned i = 0; i < count; ++i)
    {
        const unsigned length = bin.gamma();
        bm::serialization_require(length && length <= 65536u &&
              (i == 0 || starts[i-1] < starts[i]) &&
              unsigned(starts[i]) + bm::id64_t(length - 1u) * step < 65536u);
        if (IS_VALID_ADDR(block))
            for (unsigned j = 0, p = starts[i]; j < length; ++j, p += step)
            { bm::periodic_flip(block, p); } // for j
    } // for i
    if (IS_VALID_ADDR(block) && inverted)
    {
        BM_ASSERT(bm::is_aligned(block));
        bm::bit_invert(block);
    }
}


// Ordered multi-group prediction and formula coding.

/** \brief Name an ordered-search policy.
    \ingroup bvserial
    \internal
    \param order Search order identifier.
    \return Static policy name, or invalid.
*/
inline
const char* periodic_order_name(unsigned order)
{
    const char* names[] = {"legacy", "independent", "dense", "sparse", "middle"};
    return order <= 4 ? names[order] : "invalid";
}

/** \brief Bit counts for the components of a grouped formula.
    \ingroup bvserial
    \internal
*/
struct periodic_formula_bits
{
    unsigned header = 0; ///< Formula header bits.
    unsigned starts = 0; ///< Formula start-position bits.
    unsigned lengths = 0; ///< Formula length bits.
    unsigned padding = 0; ///< Terminal padding bits.
};

/** \brief Count bits in a gamma-coded positive integer.
    \ingroup bvserial
    \internal
    \param value Positive integer.
    \return Encoded bit count.
*/
inline
unsigned periodic_gamma_bits(unsigned value)
{
    BM_ASSERT(value);
    return 2 * bm::bit_scan_reverse32(value) + 1;
}

// Count the exact center-minimal BIC bits, without restarting a real bit writer.
/** \brief Calculate center-minimal BIC bit cost without writing.
    \ingroup bvserial
    \internal
    \param a Sorted positions.
    \param n Number of positions.
    \param lo Inclusive lower bound.
    \param hi Inclusive upper bound.
    \return Encoded bit count.
*/
inline
unsigned periodic_start_bits(const bm::gap_word_t* a, unsigned n,
                                    unsigned lo = 0, unsigned hi = 65535)
{
    unsigned bits = 0;
    while (n)
    {
        unsigned m = n / 2, v = a[m], r = hi - lo - n + 1;
        if (r)
        {
            unsigned range = r + 1, w = bm::bit_scan_reverse32(range);
            unsigned c = unsigned((1ull << (w + 1)) - range);
            int64_t l = int64_t(r / 2) - c / 2 - (range & 1u);
            unsigned h = r / 2 + c / 2, x = v - lo - m;
            bits += w + (int64_t(x) <= l || x > h);
        }
        if (m) bits += bm::periodic_start_bits(a, m, lo, v - 1);
        a += m + 1; n -= m + 1; lo = v + 1;
    }
    return bits;
}

/** \brief Reduced step-1 start and end positions for formula encoding.
    \ingroup bvserial
    \internal
*/
struct periodic_interleaved_data
{
    bm::gap_word_t values[512]; ///< Interleaved reduced starts and ends.
    unsigned first; ///< Explicit first start plus one.
    unsigned hi; ///< Upper bound for encoded values.
    unsigned start_bits; ///< Encoded first-start bits.
    unsigned lengths_bits; ///< Encoded interleaved payload bits.
    bool literal; ///< Whether the first start uses a literal.
};

/** \brief Prepare reduced interleaved step-1 chain endpoints.
    \ingroup bvserial
    \internal
    \tparam Alloc Application allocator.
    \param g Step-1 group with separated chains.
    \return Values and exact bit-cost components.
*/
template<class Alloc>
inline
periodic_interleaved_data periodic_prepare_interleaved(
            const periodic_candidate<Alloc>& g)
{
    bm::serialization_require(g.step == 1 &&
                           !g.chains.empty() && g.chains.size() <= 256);
    periodic_interleaved_data data;
    const unsigned n = unsigned(g.chains.size());
    unsigned previous_end = 0;
    for (unsigned i = 0; i < n; ++i)
    {
        const auto& c = g.chains[i];
        const unsigned start = c.start;
        const unsigned end = start + unsigned(c.length) - 1;
        bm::serialization_require(c.length >= 4 && end < 65536 &&
                               (!i || start >= previous_end + 2));
        previous_end = end;
    } // for i
    data.first = unsigned(g.chains[0].start) + 1;
    const unsigned first_gamma = bm::periodic_gamma_bits(data.first);
    data.literal = first_gamma >= 16;
    data.start_bits = 1 + (data.literal ? 16 : first_gamma);
    const unsigned correction = 3;
    data.hi = 65535 - (n - 1) * correction - 2;
    bm::serialization_require(data.first <= data.hi);
    for (unsigned i = 0; i < n; ++i)
    {
        const auto& c = g.chains[i];
        const unsigned offset = i * correction;
        const unsigned start = c.start;
        const unsigned end = start + unsigned(c.length) - 1;
        bm::serialization_require(start >= offset &&
                               end >= offset + 2);
        data.values[2 * i] = bm::gap_word_t(start - offset);
        data.values[2 * i + 1] = bm::gap_word_t(end - offset - 2);
        bm::serialization_require(data.values[2 * i] < data.values[2 * i + 1] &&
                               (!i || data.values[2 * i - 1] < data.values[2 * i]));
    } // for i
    data.lengths_bits = bm::periodic_start_bits(data.values + 1, 2 * n - 1,
                                            data.first, data.hi);
    return data;
}

// Exact bit cost of the explicit-first gamma or interleaved formula. The
// caller adds the common residual payload, post opcode, and frame header.
/** \brief Calculate the exact grouped formula bit cost.
    \ingroup bvserial
    \internal
    \tparam GroupsVect Candidate group container.
    \param groups Candidate groups.
    \param mode Zero for explicit starts, one for interleaved step-1 endpoints.
    \return Header, start, length and padding bit counts.
*/
template<class GroupsVect>
inline
periodic_formula_bits periodic_score_explicit_groups(
        const GroupsVect& groups, unsigned mode)
{
    bm::serialization_require(mode <= 1 && !groups.empty() && groups.size() <= 128);
    periodic_formula_bits bits;
    bits.header = 8;
    for (const auto& g : groups)
    {
        const unsigned n = unsigned(g.chains.size());
        bm::serialization_require(g.step && g.step <= 128 && n && n <= 256);
        bits.header += n > 32 ? 16 : 13;
        if (mode && g.step == 1)
        {
            auto data = bm::periodic_prepare_interleaved(g);
            bits.starts += data.start_bits;
            bits.lengths += data.lengths_bits;
            continue;
        }
        bm::gap_word_t starts[256];
        const bool step1 = g.step == 1;
        const unsigned hi = step1 ? 65536 - 4 * n : 65535;
        for (unsigned i = 0; i < n; ++i)
        {
            const unsigned start = g.chains[i].start;
            bm::serialization_require(!step1 || start >= 4 * i);
            starts[i] = bm::gap_word_t(start - (step1 ? 4 * i : 0));
        } // for i
        if (step1)
        {
            const unsigned first = unsigned(starts[0]) + 1;
            const unsigned gamma = bm::periodic_gamma_bits(first);
            bits.starts += 1 + (gamma >= 16 ? 16 : gamma);
            if (n > 1)
                bits.starts += bm::periodic_start_bits(starts + 1, n - 1,
                                                   first, hi);
        }
        else
            bits.starts += bm::periodic_start_bits(starts, n, 0, hi);
        for (const auto& c : g.chains)
        {
            bm::serialization_require(c.length >= 4);
            bits.lengths += bm::periodic_gamma_bits(unsigned(c.length) - 3);
        } // for c
    } // for g
    const unsigned meaningful = bits.header + bits.starts + bits.lengths;
    bits.padding = (32 - meaningful % 32) % 32;
    return bits;
}

// Grouped periodic post markers: set_block_periodic_groups (76),
// set_block_periodic_groups_reduced (78), set_block_periodic_groups_first (80),
// and set_block_periodic_interleaved_fixed (84). The corresponding BWT
// compositions use set_block_bwt_periodic_post (77),
// set_block_bwt_periodic_reduced (79), set_block_bwt_periodic_first (81),
// and set_block_bwt_interleaved_fixed (85).
// v1: polarity(1), groups-1(7); per group: step-1(7), count mode(1),
// count-1(5 or 8), BIC starts, gamma(length-3). Exactly one terminal flush.
// The opt-in v2 marker reduces step-1 starts by four positions per preceding
// chain. Both ends know the chain count, so the BIC ceiling shrinks as well.
// The opt-in v3 marker writes the first reduced start explicitly (gamma or
// 16-bit literal, selected with one bit) and uses it as the BIC lower boundary
// for the remaining starts. All bits share the same accumulator and flush.
// Fixed interleaving uses three positions within each chain and two between
// chains as minimum distances. The first start is explicit.
/** \brief Encode grouped periodic formulas and preserve one shared bit accumulator.
    \ingroup bvserial
    \internal
    \tparam ENC Encoder type.
    \tparam GroupsVect Candidate group container.
    \param enc Destination encoder.
    \param groups Candidate groups.
    \param inverted Original block polarity.
    \param reduce_step1_starts Reduce step-1 start positions.
    \param explicit_first_start Encode the first reduced start explicitly.
    \param interleaved_mode Use interleaved step-1 endpoints when one.
    \return Header, start, length and padding bit counts.
*/
template<class ENC, class GroupsVect>
periodic_formula_bits periodic_write_groups(ENC& enc,
                     const GroupsVect& groups, bool inverted,
                     bool reduce_step1_starts = false,
                     bool explicit_first_start = false,
                     unsigned interleaved_mode = 0)
{
    bm::serialization_require(!groups.empty() && groups.size() <= 128);
    bm::serialization_require(!explicit_first_start || reduce_step1_starts);
    bm::serialization_require(interleaved_mode <= 1);
    periodic_formula_bits bits;
    bits.header = 8;
    bm::bit_out<ENC> out(enc);
    out.put_bits(unsigned(inverted), 1);
    out.put_bits(unsigned(groups.size()) - 1, 7);
    unsigned total = 0;
    for (const auto& g : groups)
    {
        unsigned n = unsigned(g.chains.size());
        total += n;
        bm::serialization_require(g.step && g.step <= 128 && n && total <= 256);
        out.put_bits(g.step - 1, 7);
        out.put_bits(unsigned(n > 32), 1);
        out.put_bits(n - 1, n > 32 ? 8 : 5);
        bits.header += n > 32 ? 16 : 13;
        if (interleaved_mode && g.step == 1)
        {
            auto data = bm::periodic_prepare_interleaved(g);
            out.put_bits(unsigned(data.literal), 1);
            if (data.literal) out.put_bits(data.first - 1, 16);
            else out.gamma(data.first);
            bits.starts += data.start_bits;
            bits.lengths += data.lengths_bits;
            if (n * 2 > 1)
                out.bic_encode_u16_cm(data.values + 1, 2 * n - 1,
                                      bm::gap_word_t(data.first),
                                      bm::gap_word_t(data.hi));
            continue;
        }
        bm::gap_word_t starts[256];
        const bool reduced = reduce_step1_starts && g.step == 1;
        const unsigned hi = reduced ? 65536 - 4 * n : 65535;
        for (unsigned i = 0; i < n; ++i)
        {
            const unsigned start = g.chains[i].start;
            bm::serialization_require(!reduced ||
                (start >= 4 * i && start <= 65532 &&
                 (!i || start >= unsigned(g.chains[i-1].start) + 5)));
            starts[i] = bm::gap_word_t(start - (reduced ? 4 * i : 0));
        } // for i
        if (reduced && explicit_first_start)
        {
            const unsigned first = unsigned(starts[0]) + 1;
            const unsigned gamma_bits = 2 * bm::bit_scan_reverse32(first) + 1;
            const bool literal = gamma_bits >= 16;
            out.put_bits(unsigned(literal), 1);
            if (literal) out.put_bits(first - 1, 16);
            else out.gamma(first);
            bits.starts += 1 + (literal ? 16 : gamma_bits);
            if (n > 1)
            {
                bits.starts += bm::periodic_start_bits(starts + 1, n - 1,
                                                   first, hi);
                out.bic_encode_u16_cm(starts + 1, n - 1,
                                      bm::gap_word_t(first), bm::gap_word_t(hi));
            }
        }
        else
        {
            bits.starts += bm::periodic_start_bits(starts, n, 0, hi);
            out.bic_encode_u16_cm(starts, n, 0, bm::gap_word_t(hi));
        }
        for (const auto& c : g.chains)
        {
            bm::serialization_require(c.length >= 4);
            unsigned v = unsigned(c.length) - 3;
            bits.lengths += 2 * bm::bit_scan_reverse32(v) + 1;
            if (v == 1) out.put_bits(1, 1); else out.gamma(v);
        } // for c
    } // for g
    unsigned meaningful = bits.header + bits.starts + bits.lengths;
    bits.padding = (32 - meaningful % 32) % 32;
    out.flush();
    return bits;
}

/** \brief Decode and apply a grouped periodic formula.
    \ingroup bvserial
    \internal
    \tparam DEC Decoder type.
    \param dec Source decoder.
    \param block SIMD-aligned writable residual, or an invalid address to consume only.
    \param reduce_step1_starts Decode reduced step-1 starts.
    \param explicit_first_start Decode the explicit first reduced start.
    \param interleaved_mode Decode interleaved step-1 endpoints when one.
*/
template<class DEC>
void periodic_read_groups(DEC& dec, bm::word_t* block,
                          bool reduce_step1_starts = false,
                          bool explicit_first_start = false,
                          unsigned interleaved_mode = 0)
{
    bm::serialization_require(!explicit_first_start || reduce_step1_starts);
    bm::serialization_require(interleaved_mode <= 1);
    bm::bit_in<DEC> in(dec);
    unsigned inverted = in.get_bits(1), ng = in.get_bits(7) + 1, total = 0;
    for (unsigned g = 0; g < ng; ++g)
    {
        unsigned step = in.get_bits(7) + 1;
        unsigned width = in.get_bits(1) ? 8 : 5;
        unsigned count = in.get_bits(width) + 1;
        total += count;
        bm::serialization_require(total <= 256);
        if (interleaved_mode && step == 1)
        {
            const bool literal = in.get_bits(1);
            const unsigned first = literal ? in.get_bits(16) + 1 : in.gamma();
            bm::serialization_require(first);
            const bm::id64_t reduction = (count - 1) * 3 + 2;
            bm::serialization_require(reduction <= 65535);
            const unsigned hi = 65535 - unsigned(reduction);
            bm::serialization_require(first <= hi &&
                                   first + 2 * count - 2 <= hi);
            bm::gap_word_t values[512];
            values[0] = bm::gap_word_t(first - 1);
            in.bic_decode_u16_cm(values + 1, 2 * count - 1,
                                 bm::gap_word_t(first), bm::gap_word_t(hi));
            unsigned previous_end = 0;
            for (unsigned i = 0; i < count; ++i)
            {
                const unsigned offset = i * 3;
                const unsigned start = unsigned(values[2 * i]) + offset;
                const unsigned end = unsigned(values[2 * i + 1]) + offset + 2;
                bm::serialization_require(start < 65536 && end < 65536 &&
                    end >= start + 3 && (!i || start >= previous_end + 2));
                previous_end = end;
                if (IS_VALID_ADDR(block))
                    for (unsigned p = start; p <= end; ++p)
                    { bm::periodic_flip(block, p); } // for p
            } // for i
            continue;
        }
        bm::gap_word_t starts[256];
        const bool reduced = reduce_step1_starts && step == 1;
        const unsigned hi = reduced ? 65536 - 4 * count : 65535;
        if (reduced && explicit_first_start)
        {
            const bool literal = in.get_bits(1);
            const unsigned first = literal ? in.get_bits(16) + 1 : in.gamma();
            bm::serialization_require(first && first - 1 <= hi &&
                                   first + count - 1 <= hi + 1);
            starts[0] = bm::gap_word_t(first - 1);
            if (count > 1)
                in.bic_decode_u16_cm(starts + 1, count - 1,
                                     bm::gap_word_t(first), bm::gap_word_t(hi));
        }
        else
            in.bic_decode_u16_cm(starts, count, 0, bm::gap_word_t(hi));
        for (unsigned i = 0; i < count; ++i)
        {
            if (reduced)
            {
                const unsigned start = unsigned(starts[i]) + 4 * i;
                bm::serialization_require(start < 65536 &&
                    (!i || start >= unsigned(starts[i-1]) + 5));
                starts[i] = bm::gap_word_t(start);
            }
        } // for i
        for (unsigned i = 0; i < count; ++i)
        {
            unsigned v = in.gamma();
            bm::serialization_require(v && v <= 65533);
            unsigned length = v + 3;
            bm::serialization_require((!i || starts[i-1] < starts[i]) &&
                unsigned(starts[i]) + bm::id64_t(length - 1) * step < 65536);
            if (IS_VALID_ADDR(block))
                for (unsigned j = 0, p = starts[i]; j < length; ++j, p += step)
                { bm::periodic_flip(block, p); } // for j
        } // for i
    } // for g
    if (IS_VALID_ADDR(block) && inverted)
    {
        BM_ASSERT(bm::is_aligned(block));
        bm::bit_invert(block);
    }
}

/** \brief Accumulate formula component bit counts.
    \ingroup bvserial
    \internal
    \param d Destination array of four counters.
    \param b Source component counts.
*/
inline
void periodic_accumulate_bits(bm::id64_t* d, const periodic_formula_bits& b)
{ d[0] += b.header; d[1] += b.starts; d[2] += b.lengths; d[3] += b.padding; }

/** \brief Termination details for an ordered search.
    \ingroup bvserial
    \internal
*/
struct periodic_search_trace
{
    unsigned last_step = 0; ///< Last search step.
    // 0: completed the range, 1: marginal streak, 2: trial budget,
    // 3: formula chain cap.
    unsigned stop_reason = 0; ///< Search termination reason.
};

/** \brief Reusable candidate slots; inactive slots retain their chain buffers.
    \ingroup bvserial
    \internal
    \tparam Alloc Application allocator.
    Search bounds limit the slot count and each slot's chain capacity. clear()
    resets only the active count, so memory depends on peak search size rather
    than the number of blocks processed.
*/
template<class Alloc>
class periodic_candidate_slots
{
    bm::heap_vector<bm::periodic_candidate<Alloc>, Alloc, false> slots_; ///< Retained slots.
    unsigned used_ = 0; ///< Active prefix length.
public:
    /// Reserve slots without changing the active prefix. \param count Required slot count.
    void reserve(unsigned count);
    /// Reset the active prefix without releasing nested storage.
    void clear() BMNOEXCEPT;
    /// Return the active slot count.
    unsigned size() const BMNOEXCEPT;
    /// Return whether there are no active candidates.
    bool empty() const BMNOEXCEPT;
    /// Return the first active candidate address.
    bm::periodic_candidate<Alloc>* begin() BMNOEXCEPT;
    const bm::periodic_candidate<Alloc>* begin() const BMNOEXCEPT;
    /// Return one past the active prefix.
    bm::periodic_candidate<Alloc>* end() BMNOEXCEPT;
    const bm::periodic_candidate<Alloc>* end() const BMNOEXCEPT;
    /// Access an active slot. \param i Index below size().
    bm::periodic_candidate<Alloc>& operator[](unsigned i) BMNOEXCEPT;
    const bm::periodic_candidate<Alloc>& operator[](unsigned i) const BMNOEXCEPT;
    /// Access the first active slot.
    bm::periodic_candidate<Alloc>& front() BMNOEXCEPT;
    const bm::periodic_candidate<Alloc>& front() const BMNOEXCEPT;
    /// Copy into the next retained slot. \param candidate Source candidate, not a slot invalidated by growth.
    void push_back(const bm::periodic_candidate<Alloc>& candidate);
    /// Deactivate the last slot, retaining its storage.
    void pop_back() BMNOEXCEPT;
    /// Copy active candidates, retaining capacity. \param source Distinct source slots.
    void assign(const periodic_candidate_slots& source);
};

template<class Alloc>
inline
void periodic_candidate_slots<Alloc>::reserve(unsigned count)
{
    if (slots_.size() < count)
    {
        slots_.reserve(count); // relocate owning candidates by construction
        slots_.resize(count);
    }
}
template<class Alloc>
inline
void periodic_candidate_slots<Alloc>::clear() BMNOEXCEPT { used_ = 0; }
template<class Alloc>
inline
unsigned periodic_candidate_slots<Alloc>::size() const BMNOEXCEPT { return used_; }
template<class Alloc>
inline
bool periodic_candidate_slots<Alloc>::empty() const BMNOEXCEPT { return !used_; }
template<class Alloc>
inline
bm::periodic_candidate<Alloc>* periodic_candidate_slots<Alloc>::begin() BMNOEXCEPT
{ return slots_.data(); }
template<class Alloc>
inline
const bm::periodic_candidate<Alloc>* periodic_candidate_slots<Alloc>::begin() const BMNOEXCEPT
{ return slots_.data(); }
template<class Alloc>
inline
bm::periodic_candidate<Alloc>* periodic_candidate_slots<Alloc>::end() BMNOEXCEPT
{ return used_ ? begin() + used_ : begin(); }
template<class Alloc>
inline
const bm::periodic_candidate<Alloc>* periodic_candidate_slots<Alloc>::end() const BMNOEXCEPT
{ return used_ ? begin() + used_ : begin(); }
template<class Alloc>
inline
bm::periodic_candidate<Alloc>& periodic_candidate_slots<Alloc>::operator[](unsigned i) BMNOEXCEPT
{ BM_ASSERT(i < used_); return begin()[i]; }
template<class Alloc>
inline
const bm::periodic_candidate<Alloc>& periodic_candidate_slots<Alloc>::operator[](unsigned i) const BMNOEXCEPT
{ BM_ASSERT(i < used_); return begin()[i]; }
template<class Alloc>
inline
bm::periodic_candidate<Alloc>& periodic_candidate_slots<Alloc>::front() BMNOEXCEPT
{ return (*this)[0]; }
template<class Alloc>
inline
const bm::periodic_candidate<Alloc>& periodic_candidate_slots<Alloc>::front() const BMNOEXCEPT
{ return (*this)[0]; }
template<class Alloc>
inline
void periodic_candidate_slots<Alloc>::push_back(const bm::periodic_candidate<Alloc>& candidate)
{
    BM_ASSERT(used_ < slots_.size());
    slots_[used_] = candidate;
    ++used_;
}
template<class Alloc>
inline
void periodic_candidate_slots<Alloc>::pop_back() BMNOEXCEPT
{ BM_ASSERT(used_); --used_; }
template<class Alloc>
inline
void periodic_candidate_slots<Alloc>::assign(const periodic_candidate_slots& source)
{
    BM_ASSERT(this != &source);
    reserve(source.size()); clear();
    for (const auto& candidate : source) { push_back(candidate); } // for candidate
}

/** \brief Exclusive, bounded scratch for periodic candidate search.
    \ingroup bvserial
    \internal
    \tparam Alloc Application allocator.
*/
template<class Alloc>
struct periodic_search_workspace
{
    bm::heap_vector<bm::word_t, Alloc, true> blocks; ///< Two aligned blocks, allocated through the application block allocator.
    bm::periodic_candidate_slots<Alloc> groups; ///< Active search groups.
    bm::periodic_candidate_slots<Alloc> result; ///< Best group combination.
    bm::periodic_candidate_slots<Alloc> candidates; ///< Legacy search candidates.
    bm::periodic_candidate<Alloc> pool; ///< Reusable discovered chains.
    bm::periodic_candidate<Alloc> trial; ///< Reusable trial prefix.
    bm::periodic_candidate<Alloc> best; ///< Best addition at the current step.
    bm::heap_vector<unsigned, Alloc, true> steps; ///< Ordered search steps.
    bm::heap_vector<unsigned, Alloc, true> counts; ///< Chain prefix lengths.
    bm::heap_vector<unsigned, Alloc, true> candidate_order; ///< Retained score-ranked candidate indexes; candidates stay in place.
    bm::heap_vector<unsigned, Alloc, true> selected; ///< Legacy selected candidate indexes.

    /** \brief Prepare two retained, SIMD-aligned scratch blocks.
        \note Allocation is lazy and bounded. A serializer may prewarm this
        storage when periodic encoding is expected on most input blocks.
    */
    void prepare_blocks();
};

template<class Alloc>
inline
void periodic_search_workspace<Alloc>::prepare_blocks()
{
    if (blocks.size() < 2u * bm::set_block_size)
    {
        blocks.resize(2u * bm::set_block_size);
    }
    BM_ASSERT(bm::is_aligned(blocks.data()));
    BM_ASSERT(bm::is_aligned(blocks.data() + bm::set_block_size));
}

/** \brief Search periodic groups using retained caller-owned scratch.
    \ingroup bvserial
    \internal
    \tparam EVAL Complete-record evaluation callback type.
    \tparam SAVE Winning-record retention callback type.
    \tparam Alloc Application allocator.
    \param original Read-only source block.
    \param inverted Normalize by complementing the source when true.
    \param plain_size Ordinary encoded size to improve upon.
    \param o Validated search bounds and ordering policy.
    \param stats Aggregate search counters.
    \param evaluate Encode residual and groups; return complete record size.
    \param save Retain the last evaluated winning record of the supplied size.
    \param workspace Exclusive scratch, reused across blocks without cumulative growth.
    \param trace Optional last-step and termination diagnostics.
    \return Best groups, valid until the workspace is reused or destroyed.
    Each step commits only its best strict improvement over the current record.
*/
template<class EVAL, class SAVE, class Alloc>
const bm::periodic_candidate_slots<Alloc>& periodic_search_order(
        const bm::word_t* original, bool inverted, size_t plain_size,
        const periodic_options& o, periodic_statistics<Alloc>& stats,
        EVAL evaluate, SAVE save, bm::periodic_search_workspace<Alloc>& workspace,
        periodic_search_trace* trace = nullptr)
{
    workspace.prepare_blocks();
    bm::word_t* current = workspace.blocks.data();
    bm::word_t* trial = current + bm::set_block_size;
    bm::bit_block_copy(current, original);
    if (inverted)
    {
        bm::bit_invert(current);
    }
    auto& groups = workspace.groups;
    auto& result = workspace.result;
    groups.clear(); result.clear();
    groups.reserve(o.max_step - o.min_step + 2);
    size_t best_size = plain_size;
    unsigned used = 0, total_chains = 0;
    bool exhausted = false;
    unsigned marginal = 0;
    auto& ds = stats.order_stats;
    auto& steps = workspace.steps;
    steps.clear(); steps.reserve(o.max_step - o.min_step + 1);
    if (o.order == 4)
    {
        unsigned mid = (o.min_step + o.max_step) / 2;
        steps.push_back(mid);
        for (unsigned d = 1; steps.size() < o.max_step - o.min_step + 1; ++d)
        {
            if (mid + d <= o.max_step) steps.push_back(mid + d);
            if (mid >= o.min_step + d) steps.push_back(mid - d);
        } // for d
    }
    else
    {
        for (unsigned step = o.min_step; step <= o.max_step; ++step)
        { steps.push_back(step); } // for step
        if (o.order == 3)
        {
            for (size_t left = 0, right = steps.size(); left < right && left < --right; ++left)
            {
                bm::simple_swap(&steps[left], &steps[right]);
            } // for left
        }
    }
    for (unsigned step : steps)
    {
        if (used == o.trial_budget)
        { exhausted = true; if (trace) trace->stop_reason = 2; break; }
        if (total_chains == o.formula_chains)
        { if (trace) trace->stop_reason = 3; break; }
        if (trace) trace->last_step = step;
        ++ds.step_visits[step];
        const size_t before = best_size;
        auto& best_addition = workspace.best;
        for (unsigned holes = 0; holes <= o.max_holes; ++holes)
        {
            auto& pool = workspace.pool;
            ++stats.scanned;
            bool found = bm::periodic_find(current, false, step, holes, o.min_length,
                                       o.max_chains, pool, &ds);
            stats.chain_cap += pool.truncated;
            if (!found) continue;
            ++stats.raw_candidates;
            stats.positive_candidates += pool.score > 0;
            // Deterministic ranking; equal benefits prefer shorter descriptions,
            // then smaller starts. Prefix choice still uses actual record bytes.
            bm::heap_sort(pool.chains.data(), pool.chains.size(),
                [](const periodic_chain& a, const periodic_chain& b)
                {
                    int x = int(a.length) - 2 * a.holes;
                    int y = int(b.length) - 2 * b.holes;
                    if (x != y) return x > y;
                    if (a.length != b.length) return a.length < b.length;
                    return a.start < b.start;
                });
            const unsigned chain_count = unsigned(pool.chains.size());
            const unsigned available = o.formula_chains - total_chains;
            unsigned limit = chain_count < available ? chain_count : available;
            auto& counts = workspace.counts; counts.clear();
            counts.reserve(10);
            for (unsigned n = 1; n < limit; n *= 2)
            { counts.push_back(n); } // for n
            counts.push_back(limit);
            for (unsigned n : counts)
            {
                if (used == o.trial_budget) { exhausted = true; break; }
                auto& c = workspace.trial;
                c.covered = c.exceptions = 0; c.score = 0; c.truncated = false;
                c.step = step; c.max_holes = holes;
                c.chains.resize(n);
                ::memcpy(c.chains.data(), pool.chains.data(), n * sizeof(periodic_chain));
                bm::heap_sort(c.chains.data(), c.chains.size(),
                    [](const periodic_chain& a, const periodic_chain& b)
                    { return a.start < b.start; });
                for (const auto& chain : c.chains)
                { c.covered += chain.length; c.exceptions += chain.holes; } // for chain
                bm::bit_block_copy(trial, current);
                bm::periodic_apply(trial, c);
                groups.push_back(c);
                ++used; ++stats.trials; ++stats.trial_step[step];
                ++ds.trial_counts[total_chains + n];
                size_t size = evaluate(trial, groups);
                if (size < best_size)
                {
                    best_size = size;
                    best_addition = c;
                    result.assign(groups);
                    save(size);
                }
                groups.pop_back();
            } // for n
        } // for holes
        if (best_size < before && o.order != 1)
        {
            bm::periodic_apply(current, best_addition);
            total_chains += unsigned(best_addition.chains.size());
            groups.push_back(best_addition);
            ++ds.additions; ++ds.groups_step[step];
            ds.gain_step[step] += before - best_size;
        }
        if (o.marginal_streak)
        {
            marginal = before - best_size >= o.min_step_gain ? 0 : marginal + 1;
            if (step >= o.min_stop_step && marginal >= o.marginal_streak)
            { if (trace) trace->stop_reason = 1; break; }
        }
    } // for step
    // Independent trials always use the original normalized block and one group.
    if (o.order == 1 && !result.empty())
    {
        ++ds.additions; ++ds.groups_step[result[0].step];
        ds.gain_step[result[0].step] += plain_size - best_size;
    }
    ds.exhausted += exhausted;
    return result;
}

/** \brief Merge ordered-search counters and observations.
    \ingroup bvserial
    \internal
    \tparam Alloc Application allocator.
    \param d Destination statistics.
    \param s Source statistics.
*/
template<class Alloc>
inline
void periodic_order_add(periodic_order_statistics<Alloc>& d,
                               const periodic_order_statistics<Alloc>& s)
{
    for (unsigned i = 0; i < 44; ++i)
    { d.discovered[i] += s.discovered[i]; d.retained[i] += s.retained[i];
      d.selected[i] += s.selected[i]; } // for i
    for (unsigned i = 0; i < 257; ++i)
    { d.trial_counts[i] += s.trial_counts[i]; d.selected_counts[i] += s.selected_counts[i]; } // for i
    for (unsigned i = 0; i < 129; ++i)
    { d.group_counts[i] += s.group_counts[i]; d.step_visits[i] += s.step_visits[i];
      d.groups_step[i] += s.groups_step[i]; d.gain_step[i] += s.gain_step[i]; } // for i
    for (unsigned i = 0; i < 4; ++i)
    { d.trial_bits[i] += s.trial_bits[i]; d.selected_bits[i] += s.selected_bits[i]; } // for i
    for (unsigned i = 0; i < 2; ++i)
    {
        d.trial_formula_modes[i] += s.trial_formula_modes[i];
        d.selected_formula_modes[i] += s.selected_formula_modes[i];
    } // for i
    const size_t offset = d.blocks.size(), count = s.blocks.size();
    d.blocks.resize(offset + count);
    for (size_t i = 0; i < count; ++i)
    {
        d.blocks[offset + i] = s.blocks[i];
    } // for i
    d.exhausted += s.exhausted; d.multi_group += s.multi_group; d.additions += s.additions;
}



} // namespace bm
#endif // BM_PERIODIC_H_INCLUDED
