// Allocator propagation regression shared by the 32- and 64-bit stress suites.
#ifndef BM_SERIALIZATION_ALLOCATOR_TESTS_H__
#define BM_SERIALIZATION_ALLOCATOR_TESTS_H__
#include <cassert>
#include "bm.h"
#include "bmserial.h"
#include <type_traits>
#include <utility>
#include <iostream>

namespace serialization_allocator_test
{
struct tracking_block_allocator
{
    static inline size_t allocations = 0; ///< Number of allocations requested.
    static inline size_t fail_after = size_t(-1); ///< Allocation-failure injection countdown.
    static inline size_t live = 0;        ///< Outstanding block allocations.

    static bm::word_t* allocate(size_t n, const void* hint)
    {
        if (!fail_after) throw std::bad_alloc();
        if (fail_after != size_t(-1)) --fail_after;
        bm::word_t* p = bm::block_allocator::allocate(n, hint);
        ++allocations;
        ++live;
        return p;
    }
    static void deallocate(bm::word_t* p, size_t n) noexcept
    {
        assert(live);
        --live;
        bm::block_allocator::deallocate(p, n);
    }
};
using Alloc = bm::mem_alloc<tracking_block_allocator, bm::ptr_allocator,
    bm::alloc_pool<tracking_block_allocator, bm::ptr_allocator> >;
using BV = bm::bvector<Alloc>;
using Serializer = bm::serializer<BV>;
using Stats = Serializer::block_selection_statistics_type;
static_assert(std::is_same<Stats::allocator_type, Alloc>::value,
              "Selection observations must use the application allocator");

struct reader : bm::deserializer<BV, bm::decoder>
{
    using bm::deserializer<BV, bm::decoder>::sblock_adaptive_corrections_;
    using bm::deserializer<BV, bm::decoder>::sblock_adaptive_flags_;
};
struct iterator : bm::serial_stream_iterator<bm::decoder, BV::block_idx_type, Alloc>
{
    using base = bm::serial_stream_iterator<bm::decoder, BV::block_idx_type, Alloc>;
    using base::base;
    using base::sblock_adaptive_corrections_;
    using base::sblock_adaptive_flags_;
};

template<class Reader> void check_scratch(Reader& r)
{
    const size_t before = tracking_block_allocator::allocations;
    r.sblock_adaptive_corrections_.reserve(bm::sblock_adaptive_max_windows);
    r.sblock_adaptive_flags_.reserve(bm::sblock_adaptive_max_windows);
    assert(tracking_block_allocator::allocations == before + 2);
}

// Independent transform and ownership checks use the same custom allocator
// as the bvector, including nested periodic chain descriptions.
inline void CheckCodecWorkspaces()
{
    {
        bm::periodic_search_workspace<Alloc> periodic;
        const size_t initial = tracking_block_allocator::allocations;
        periodic.prepare_blocks();
        assert(tracking_block_allocator::allocations == initial + 1);
        bm::word_t* block = periodic.blocks.data();
        bm::word_t* guard = block + bm::set_block_size;
        assert(bm::is_aligned(block) && bm::is_aligned(guard));
        periodic.prepare_blocks();
        assert(periodic.blocks.data() == block);
        assert(tracking_block_allocator::allocations == initial + 1);

        Serializer::periodic_candidate_type chain;
        chain.step = 1;
        chain.chains.push_back({7, 5, 0});
        Serializer::periodic_candidates_type groups;
        groups.push_back(chain);
        Serializer::buffer bytes;
        bytes.resize(4096);
        for (unsigned grouped = 0; grouped < 2; ++grouped)
        {
            bm::encoder enc(bytes.data(), bytes.size());
            if (grouped) bm::periodic_write_groups(enc, groups, true);
            else bm::periodic_write_formula(enc, chain, true);
            bm::bit_block_set(block, 0);
            // Set each word explicitly: the scalar block fill repeats bytes.
            for (unsigned i = 0; i < bm::set_block_size; ++i)
            {
                guard[i] = 0x12345678u;
            } // for i
            bm::decoder dec(bytes.data());
            if (grouped) bm::periodic_read_groups(dec, block);
            else bm::periodic_read_formula(dec, block);
            for (unsigned i = 0; i < bm::set_block_size; ++i)
            {
                const bm::word_t mask = i == 0 ? 0xF80u : 0u;
                assert(block[i] == ~mask);
                assert(guard[i] == 0x12345678u);
            } // for i
            assert(dec.get_pos() == bytes.data() + enc.size());
        } // for grouped
    }
    bm::bwt_full_workspace<Alloc> workspace;
    size_t before = tracking_block_allocator::allocations;
    workspace.prepare();
    assert(tracking_block_allocator::allocations == before + 5);
    before = tracking_block_allocator::allocations;
    workspace.prepare();
    assert(tracking_block_allocator::allocations == before);
    BM_DECLARE_TEMP_BLOCK(source)
    BM_DECLARE_TEMP_BLOCK(transformed)
    bm::bit_block_set(source, 0);
    for (unsigned i = 0; i < 65536; ++i)
    {
        if (i % 101 < 13) bm::set_bit(source, i);
    } // for i
    const auto primary = bm::bwt_transform_full(source, transformed, workspace);
    bm::bwt_inverse_full(transformed, primary, workspace);
    assert(!::memcmp(source, transformed, sizeof(bm::word_t) * bm::set_block_size));
    assert(tracking_block_allocator::allocations == before);

    Serializer::periodic_candidate_type candidate;
    Serializer::periodic_statistics_type statistics;
    assert(bm::periodic_find(source, false, 1, 0, 4, 256,
                            candidate, &statistics.order_stats));
    Serializer::periodic_candidates_type groups;
    groups.reserve(4);
    groups.push_back(candidate);
    groups.push_back(candidate);
    // Fail the outer allocation and each nested chain copy in turn.
    const size_t live_before_copy = tracking_block_allocator::live;
    for (size_t failure = 0; failure < 3; ++failure)
    {
        tracking_block_allocator::fail_after = failure;
        bool rejected = false;
        try { auto failed_copy = groups; }
        catch (const std::bad_alloc&) { rejected = true; }
        tracking_block_allocator::fail_after = size_t(-1);
        assert(rejected && tracking_block_allocator::live == live_before_copy);
    } // for failure
    const size_t old_capacity = groups.capacity();
    tracking_block_allocator::fail_after = 2;
    bool reserve_failed = false;
    try { groups.reserve(old_capacity + 1); }
    catch (const std::bad_alloc&) { reserve_failed = true; }
    tracking_block_allocator::fail_after = size_t(-1);
    assert(reserve_failed && groups.capacity() == old_capacity);
    assert(groups.size() == 2 && tracking_block_allocator::live == live_before_copy);
    auto copy = groups;
    const auto original_start = groups[0].chains[0].start;
    copy[0].chains[0].start = bm::gap_word_t(original_start + 1);
    assert(groups[0].chains[0].start == original_start);
    const size_t live_before_assignment = tracking_block_allocator::live;
    tracking_block_allocator::fail_after = 2;
    bool assignment_failed = false;
    try { groups = copy; }
    catch (const std::bad_alloc&) { assignment_failed = true; }
    tracking_block_allocator::fail_after = size_t(-1);
    assert(assignment_failed && tracking_block_allocator::live == live_before_assignment);
    assert(groups[0].chains[0].start == original_start);
    groups = copy;
    copy.clear();
    assert(groups[0].chains[0].start == original_start + 1);
    groups.pop_back();
    assert(groups.size() == 1);
    before = tracking_block_allocator::allocations;
    auto moved = std::move(groups);
    assert(groups.empty() && moved.size() == 1);
    groups = std::move(moved);
    assert(moved.empty() && groups.size() == 1);
    assert(tracking_block_allocator::allocations == before);

    bm::heap_vector<unsigned, Alloc, true> values;
    values.resize(256);
    for (unsigned i = 0; i < values.size(); ++i)
    {
        values[i] = 7 + i * 67;
    } // for i
    bm::sblock_wave_record record;
    bm::sblock_wave1024_model<Alloc> fixed;
    bm::sblock_wave_adaptive_model<Alloc> adaptive;
    assert(fixed.prepare(values.data(), unsigned(values.size()), record));
    adaptive.build(values.data(), unsigned(values.size()));
    assert(adaptive.prepare(values.data(), unsigned(values.size()), record));
    before = tracking_block_allocator::allocations;
    assert(fixed.prepare(values.data(), unsigned(values.size()), record));
    adaptive.build(values.data(), unsigned(values.size()));
    assert(adaptive.prepare(values.data(), unsigned(values.size()), record));
    assert(tracking_block_allocator::allocations == before);

    Serializer::buffer bytes;
    bytes.resize(65536);
    for (unsigned mode = 0; mode < 4; ++mode)
    {
        bm::encoder enc(bytes.data(), bytes.size());
        {
            bm::bit_out<bm::encoder> out(enc);
            bm::sblock_wave_cost cost;
            fixed.encode(out, [&]() { return bm::id64_t(0); },
                         3, mode, cost);
        }
        bm::decoder dec(bytes.data());
        bm::bit_in<bm::decoder> in(dec);
        unsigned restored[256], sb = 0;
        before = tracking_block_allocator::allocations;
        assert(bm::sblock_wave1024_model<Alloc>::decode(in, restored, sb) == 256);
        assert(tracking_block_allocator::allocations > before);
        assert(sb == 3 && !::memcmp(restored, values.data(), sizeof(restored)));
        bm::heap_vector<unsigned, Alloc, true> map, occupied;
        for (unsigned pass = 0; pass < 4; ++pass)
        {
            bm::decoder repeated(bytes.data());
            bm::bit_in<bm::decoder> bits(repeated);
            before = tracking_block_allocator::allocations;
            assert(bm::sblock_wave1024_model<Alloc>::decode(bits, restored, sb,
                                                         map, occupied) == 256);
            // Non-inverse mode swaps the two retained buffers; warm both sides.
            if (pass > 1) assert(tracking_block_allocator::allocations == before);
            assert(sb == 3 && !::memcmp(restored, values.data(), sizeof(restored)));
        } // for pass
    } // for mode
}

/** \brief Check fixed-wave reset, capacity retention and serializer reuse.
    Exercises eligible, smaller, ineligible and eligible inputs without
    reconstructing the model or serializer between records.
*/
inline
void CheckWaveStateReuse()
{
    using Model = bm::sblock_wave1024_model<Alloc>;
    Model model;
    bm::sblock_wave_record record;
    bm::heap_vector<unsigned, Alloc, true> positions;
    positions.resize(32768);
    for (unsigned pass = 0; pass < 3; ++pass)
    {
        for (unsigned kind = 0; kind < 4; ++kind)
        {
            const unsigned n = kind == 0 ? 32768u : kind == 2 ?
                               unsigned(Model::wave_count) : 100u;
            positions.resize(n);
            for (unsigned i = 0; i < n; ++i)
            {
                positions[i] = kind == 2 ? i * Model::wave_size :
                               kind == 0 ? 17 + i * 67 : 1 + i * 163841;
            } // for i
            record.eligible = true;
            record.best = 3;
            for (unsigned mode = 0; mode < 4; ++mode)
            {
                record.costs[mode].bytes = 7;
                record.costs[mode].header = 8;
                record.costs[mode].map = 9;
                record.costs[mode].bic = 10;
                record.costs[mode].padding = 11;
            } // for mode
            const size_t before = tracking_block_allocator::allocations;
            const bool eligible = model.prepare(positions.data(), n, record);
            if (pass)
            {
                assert(tracking_block_allocator::allocations == before);
            }
            assert(eligible == (kind != 2) && record.eligible == eligible);
            assert(record.best == 0);
            for (unsigned mode = 0; mode < 4; ++mode)
            {
                const auto& cost = record.costs[mode];
                assert(!cost.bytes && !cost.header && !cost.map &&
                       !cost.bic && !cost.padding);
            } // for mode
            if (!eligible)
            {
                assert(model.correction == 0 && model.upper == 0);
                assert(model.empty.empty() && model.occupied.size() == Model::wave_count);
            }
            else
            {
                // An independently prepared model must have identical state.
                Model fresh;
                bm::sblock_wave_record fresh_record;
                assert(fresh.prepare(positions.data(), n, fresh_record));
                assert(model.values.equal(fresh.values));
                assert(model.occupied.equal(fresh.occupied));
                assert(model.empty.equal(fresh.empty));
                assert(model.correction == fresh.correction && model.upper == fresh.upper);
            }
        } // for kind
    } // for pass

    const unsigned* occupied_data = model.occupied.data();
    const unsigned* empty_data = model.empty.data();
    const unsigned* values_data = model.values.data();
    const size_t occupied_capacity = model.occupied.capacity();
    const size_t empty_capacity = model.empty.capacity();
    const size_t values_capacity = model.values.capacity();
    const size_t before_init = tracking_block_allocator::allocations;
    model.init();
    assert(model.occupied.empty() && model.empty.empty() && model.values.empty());
    assert(model.correction == 0 && model.upper == 0);
    assert(model.occupied.data() == occupied_data && model.empty.data() == empty_data &&
           model.values.data() == values_data);
    assert(model.occupied.capacity() == occupied_capacity &&
           model.empty.capacity() == empty_capacity && model.values.capacity() == values_capacity);
    assert(tracking_block_allocator::allocations == before_init);

    // Adaptive width trials preserve the input cache while resetting trial state.
    bm::sblock_wave_adaptive_model<Alloc> adaptive;
    adaptive.build(positions.data(), unsigned(positions.size()));
    const unsigned* cached_data = adaptive.smallest.data();
    const size_t cached_size = adaptive.smallest.size();
    for (unsigned pass = 0; pass < 2; ++pass)
    {
        for (unsigned width = 0; width < 6; ++width)
        {
            adaptive.wave_size = 128u << width;
            adaptive.wave_count = bm::set_sub_total_bits / adaptive.wave_size;
            const size_t before = tracking_block_allocator::allocations;
            assert(adaptive.prepare(positions.data(), unsigned(positions.size()), record));
            if (pass)
            {
                assert(tracking_block_allocator::allocations == before);
            }
            assert(adaptive.smallest.data() == cached_data &&
                   adaptive.smallest.size() == cached_size);
            adaptive.init_trial();
            assert(adaptive.smallest.data() == cached_data &&
                   adaptive.smallest.size() == cached_size);
            assert(adaptive.occupied.empty() && adaptive.empty.empty() &&
                   adaptive.values.empty() && !adaptive.correction && !adaptive.upper);
            assert(adaptive.wave_size == (128u << width));
        } // for width
    } // for pass
    const size_t before_adaptive_init = tracking_block_allocator::allocations;
    adaptive.init();
    assert(adaptive.smallest.empty() && adaptive.smallest.data() == cached_data);
    assert(adaptive.occupied.empty() && adaptive.empty.empty() && adaptive.values.empty());
    assert(!adaptive.correction && !adaptive.upper && adaptive.wave_size == 4096);
    assert(tracking_block_allocator::allocations == before_adaptive_init);
    positions.resize(Model::wave_count);
    for (unsigned i = 0; i < positions.size(); ++i)
    {
        positions[i] = i * Model::wave_size;
    } // for i
    adaptive.build(positions.data(), unsigned(positions.size()));
    adaptive.wave_size = Model::wave_size;
    adaptive.wave_count = Model::wave_count;
    assert(!adaptive.prepare(positions.data(), unsigned(positions.size()), record));
    assert(!record.eligible && !record.best && !adaptive.correction && !adaptive.upper);

    BV sources[4];
    for (unsigned kind = 0; kind < 4; ++kind)
    {
        const unsigned n = kind == 2 ? unsigned(Model::wave_count) : kind == 1 ? 40u : 100u;
        for (unsigned i = 0; i < n; ++i)
        {
            sources[kind].set(kind == 2 ? i * Model::wave_size : 1 + i * 163841 + kind);
        } // for i
        sources[kind].optimize();
    } // for kind
    Serializer serializer;
    serializer.set_compression_level(6);
    serializer.set_sblock_wave1024_options(true);
    Serializer::buffer bytes;
    for (unsigned pass = 0; pass < 3; ++pass)
    {
        for (unsigned kind = 0; kind < 4; ++kind)
        {
            const size_t before = tracking_block_allocator::allocations;
            serializer.serialize(sources[kind], bytes);
            if (pass)
            {
                assert(tracking_block_allocator::allocations == before);
            }
            assert(bool(serializer.get_compression_stat()[bm::set_sblock_bienc_wave1024]) ==
                   (kind != 2));
            BV restored;
            bm::deserialize(restored, bytes.data());
            assert(restored == sources[kind]);
        } // for kind
    } // for pass
}

inline void CheckIteratorFailure(const unsigned char* bytes)
{
    const size_t live = tracking_block_allocator::live;
    for (size_t failure = 0; failure < 3; ++failure)
    {
        tracking_block_allocator::fail_after = failure;
        bool rejected = false;
        try { iterator stream(bytes); }
        catch (const std::bad_alloc&) { rejected = true; }
        tracking_block_allocator::fail_after = size_t(-1);
        assert(rejected && tracking_block_allocator::live == live);
    } // for failure
}

// Repeated codec work must retain scratch capacity rather than allocate per block.
inline void CheckRetainedScratch()
{
    BM_DECLARE_TEMP_BLOCK(column)
    BM_DECLARE_TEMP_BLOCK(restored)
    bm::bit_block_set(column, 0);
    for (unsigned i = 0; i < bm::gap_max_bits; ++i)
    { if (i % 101 < 13) bm::set_bit(column, i); } // for i
    Serializer::buffer payload;
    payload.resize(256 * 1024);
    bm::heap_vector<bm::gap_word_t, Alloc, true> positions;
    bm::bwt_full_workspace<Alloc> inverse;
    for (unsigned pass = 0; pass < 5; ++pass)
    {
        const size_t before = tracking_block_allocator::allocations;
        bm::encoder enc(payload.data(), payload.size());
        bm::bwt_bic_encoder<Alloc>::encode(enc, column, 0, positions);
        bm::decoder dec(payload.data());
        bm::bwt_bic_decoder<Alloc>::decode(dec, restored, inverse);
        assert(dec.size() == enc.size());
        if (pass) assert(tracking_block_allocator::allocations == before);
    } // for pass
    assert(inverse.rank.size() == bm::gap_max_bits);
    assert(inverse.order.size() == bm::set_block_size);
    assert(inverse.shifted.empty() && inverse.next_rank.empty() && inverse.counts.empty());

    BV source;
    for (unsigned i = 0; i < bm::gap_max_bits; ++i)
    { if (i % 101 < 13) source.set(i); } // for i
    Serializer serializer;
    serializer.set_compression_level(10);
    // Isolate block scratch from whole-vector index/output data structures.
    serializer.set_sblock_large_options(false);
    Serializer::buffer bytes;
    for (unsigned pass = 0; pass < 5; ++pass)
    {
        const size_t before = tracking_block_allocator::allocations;
        serializer.serialize(source, bytes);
#if !BM_SERIALIZATION_DIAGNOSTICS
        assert(serializer.get_bwt_statistics().blocks.empty());
        assert(serializer.get_periodic_statistics().order_stats.blocks.empty());
        assert(serializer.get_block_selection_statistics().blocks.empty());
        if (pass) assert(tracking_block_allocator::allocations == before);
#else
        (void)before; // Explicitly enabled development histories may grow.
        assert(!serializer.get_bwt_statistics().blocks.empty());
        assert(!serializer.get_periodic_statistics().order_stats.blocks.empty());
        assert(!serializer.get_block_selection_statistics().blocks.empty());
#endif
    } // for pass
    bm::serial_iterator_workspace<Alloc> iterator_workspace;
    for (unsigned pass = 0; pass < 5; ++pass)
    {
        const size_t before = tracking_block_allocator::allocations;
        iterator cursor(bytes.data(), &iterator_workspace);
        if (pass) assert(tracking_block_allocator::allocations == before);
    } // for pass

    // FULL OR a BWT record consumes the payload with a cold inverse workspace.
    {
        bm::encoder enc(payload.data(), payload.size());
        enc.put_8(bm::BM_HM_NO_BO | bm::BM_HM_NO_GAPL);
        enc.put_8(bm::set_block_bwt_unpacked);
        bm::bwt_bic_encoder<Alloc>::encode(enc, column, 0, positions);
        enc.put_8(bm::set_block_end);
        struct cold_reader : bm::deserializer<BV, bm::decoder>
        {
            bool inverse_empty() const
            {
                return this->bwt_decode_storage_.rank.empty() &&
                       this->bwt_decode_storage_.order.empty() &&
                       this->bwt_decode_storage_.last_column.empty();
            }
        };
        BV full; full.set_range(0, bm::gap_max_bits-1); full.optimize();
        full.set(bm::gap_max_bits+7);
        cold_reader reader;
        assert(reader.inverse_empty());
        const size_t before = tracking_block_allocator::allocations;
        assert(reader.deserialize(full, payload.data()) == enc.size());
        assert(reader.inverse_empty());
        assert(tracking_block_allocator::allocations == before);
        assert(IS_FULL_BLOCK(full.get_blocks_manager().get_block_ptr(0, 0)));
        bm::operation_deserializer<BV> operation;
        size_t allocated = tracking_block_allocator::allocations;
        operation.deserialize(full, payload.data(), bm::set_OR);
        // Only the three fixed iterator buffers are initialized, no BWT scratch.
        assert(tracking_block_allocator::allocations == allocated + 3);
        allocated = tracking_block_allocator::allocations;
        operation.deserialize(full, payload.data(), bm::set_OR);
        assert(tracking_block_allocator::allocations == allocated);
        assert(full.count() == bm::gap_max_bits+1 && full.test(bm::gap_max_bits+7));
    }

    // Index codecs retain corrected-value scratch and transactional output.
    bm::deserialization_index<BV> index, decoded_index;
    bm::deserialization_index_serializer<BV> index_writer;
    bm::deserialization_index_deserializer<BV> index_reader;
    Serializer::buffer index_bytes;
    for (unsigned pass = 0; pass < 5; ++pass)
    {
        const size_t before = tracking_block_allocator::allocations;
        for (unsigned count : {512u, 16u, 512u})
        {
            index.clear();
            index.set_serialized_size(20000);
            for (unsigned i = 0; i < count; ++i)
            { index.add_marker_offset(i * 31 + i % 3); } // for i
            index_writer.serialize(index, index_bytes);
            assert(index_reader.deserialize(decoded_index, index_bytes.data(),
                                            index_bytes.size()) == index_bytes.size());
            assert(decoded_index.size() == index.size());
            assert(decoded_index.serialized_size() == index.serialized_size());
        } // for count
        if (pass > 1) assert(tracking_block_allocator::allocations == before);
    } // for pass

    // Resetting a periodic candidate must retain its nested chain capacity.
    Serializer::periodic_candidate_type candidate;
    bm::periodic_find(column, false, 1, 0, 4, 256, candidate);
    const size_t before = tracking_block_allocator::allocations;
    bm::periodic_find(column, false, 1, 0, 4, 256, candidate);
    assert(tracking_block_allocator::allocations == before);
}

/** Validate legacy periodic ranking and retained storage through repeated encoding. */
inline void CheckPeriodicIndexRanking()
{
    BV source;
    for (unsigned i = 0; i < bm::gap_max_bits; ++i)
    {
        if (i % 37 < 18 || i % 137 < 6)
        {
            source.set(i);
        }
    } // for i
    for (unsigned diagnostics = 0; diagnostics < 2; ++diagnostics)
    {
        Serializer serializer;
        bm::periodic_options options;
        options.order = 0;
        options.max_step = 8;
        options.max_chains = 16;
        options.max_variants = 4;
        options.diagnostics = diagnostics != 0;
        serializer.set_periodic_options(options);
        Serializer::buffer bytes;
        Serializer::buffer expected;
        for (unsigned pass = 0; pass < 3; ++pass)
        {
            const size_t before = tracking_block_allocator::allocations;
            serializer.serialize(source, bytes);
            if (pass)
            {
                assert(tracking_block_allocator::allocations == before);
                assert(bytes.size() == expected.size());
                assert(::memcmp(bytes.data(), expected.data(), bytes.size()) == 0);
            }
            else
            {
                expected.resize(bytes.size());
                ::memcpy(expected.data(), bytes.data(), bytes.size());
            }
            BV restored;
            bm::deserialize(restored, bytes.data());
            assert(restored == source);
        } // for pass
        assert(serializer.get_periodic_statistics().trials > 0);
    } // for diagnostics
}

inline void Run()
{
    {
        CheckCodecWorkspaces();
        CheckWaveStateReuse();
        CheckRetainedScratch();
        CheckPeriodicIndexRanking();
        Stats stats;
        size_t before = tracking_block_allocator::allocations;
        bm::block_selection_observation observation;
        stats.add_block(observation);
        assert(tracking_block_allocator::allocations == before + 1);
        before = tracking_block_allocator::allocations;
        Stats copy(stats);
        assert(tracking_block_allocator::allocations == before + 1);
        bm::block_selection_add(copy, stats);
        assert(copy.blocks.size() == 2);
        stats = Stats();
        assert(stats.blocks.empty());
        reader r;
        check_scratch(r);

        BV bv;
        for (unsigned i = 0; i < 4096; ++i)
        {
            bv.set(7 + i * 67);
        } // for i
        for (unsigned level = 6; level <= 10; ++level)
        {
            Serializer serializer;
            serializer.set_compression_level(level);
            Serializer::buffer bytes;
            serializer.serialize(bv, bytes);
            BV restored;
            bm::deserialize(restored, bytes.data());
            assert(restored == bv);
            bm::operation_deserializer<BV> operation;
            restored.clear();
            operation.deserialize(restored, bytes.data(), bm::set_OR);
            assert(restored == bv);
            CheckIteratorFailure(bytes.data());
            const size_t before_iterator = tracking_block_allocator::allocations;
            iterator stream(bytes.data());
            assert(tracking_block_allocator::allocations == before_iterator + 3);
            check_scratch(stream);
        } // for level
    }
    assert(tracking_block_allocator::live == 0);
    std::cout << "Serialization allocator regression PASS\n";
}

} // namespace serialization_allocator_test
#endif
