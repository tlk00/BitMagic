// Opt-in tests; does not enable this experiment in library stress tests.
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <cstring>
#include <random>
#ifdef TEST_BM64
#include "bm64.h"
#else
#include "bm.h"
#endif
#include "bmserial.h"
#include "bmfio.h"

static void require(bool value, const char* message)
{ if (!value) throw std::runtime_error(message); }

static void formulas()
{
    std::mt19937 rng(98317);
    for (unsigned n : {1u, 2u, 31u, 32u, 33u, 128u, 256u})
        for (unsigned ng : {1u, 2u, 128u})
        {
            if (n * ng > 256) continue;
            for (bool inverse : {false, true})
            for (unsigned start_mode : {0u, 1u, 2u, 3u})
            {
                const bool reduced = start_mode != 0;
                const bool first = start_mode >= 2;
                const unsigned interleaved = start_mode >= 3 ? start_mode - 2 : 0;
                std::vector<bm::periodic_candidate<>> groups(ng);
                BM_DECLARE_TEMP_BLOCK(expected)
                BM_DECLARE_TEMP_BLOCK(restored)
                for (unsigned i = 0; i < bm::set_block_size; ++i)
                    restored[i] = expected[i] = rng();
                for (unsigned g = 0; g < ng; ++g)
                {
                    auto& c = groups[g]; c.step = g % 2 ? 1 : 128;
                    for (unsigned j = 0; j < n; ++j)
                        c.chains.push_back({bm::gap_word_t(j * 200 + g),
                                            bm::gap_word_t(4 + rng() % 16), 0});
                    bm::periodic_apply(expected, c);
                }
                if (inverse)
                    for (unsigned i = 0; i < bm::set_block_size; ++i) expected[i] = ~expected[i];
                unsigned char buf[8192]; bm::encoder enc(buf, sizeof(buf));
                auto bits = bm::periodic_write_groups(enc, groups, inverse,
                                             reduced, first, interleaved);
                require(enc.size() * 8 == bits.header + bits.starts + bits.lengths + bits.padding,
                        "formula bit accounting");
                if (first)
                {
                    const auto scored = bm::periodic_score_explicit_groups(groups,
                                                                           interleaved);
                    require(scored.header == bits.header &&
                            scored.starts == bits.starts &&
                            scored.lengths == bits.lengths &&
                            scored.padding == bits.padding,
                            "scored formula differs from emitted formula");
                }
                bm::bounded_decoder dec; dec.bind(buf, enc.size());
                bm::periodic_read_groups(dec, restored, reduced, first,
                                         interleaved);
                require(dec.size() == enc.size(), "formula exact consumption");
                require(!std::memcmp(expected, restored, sizeof(bm::word_t)*bm::set_block_size),
                        "overlapping formula reconstruction");
                bm::bounded_decoder skip; skip.bind(buf, enc.size());
                bm::periodic_read_groups(skip, static_cast<bm::word_t*>(nullptr),
                                         reduced, first, interleaved);
                require(skip.size() == enc.size(), "formula skip consumption");
                for (size_t size = 0; size < enc.size(); ++size)
                {
                    bm::bounded_decoder truncated; truncated.bind(buf, size);
                    bool rejected = false;
                    try { bm::periodic_read_groups(truncated, static_cast<bm::word_t*>(nullptr),
                                                   reduced, first,
                                                   interleaved); }
                    catch (...) { rejected = true; }
                    require(rejected, "truncated formula accepted");
                }
            }
        }
}

static void step1_boundaries()
{
    for (unsigned n : {1u, 2u, 32u, 256u})
    for (unsigned offset : {0u, 127u, 128u, 60000u, 65535u})
    {
        bm::periodic_candidate<> group; group.step = 1;
        // Five-bit spacing is the tightest possible for distinct length-4 chains.
        const unsigned first = std::min(offset, 65536 - (5 * n - 1));
        for (unsigned i = 0; i < n; ++i)
            group.chains.push_back({bm::gap_word_t(first + 5 * i), 4, 0});
        BM_DECLARE_TEMP_BLOCK(expected)
        BM_DECLARE_TEMP_BLOCK(restored)
        bm::bit_block_set(expected, 0);
        bm::bit_block_set(restored, 0);
        bm::periodic_apply(expected, group);
        for (unsigned mode = 0; mode <= 1; ++mode)
        {
            bm::bit_block_set(restored, 0);
            unsigned char buf[8192]; bm::encoder enc(buf, sizeof(buf));
            bm::serializer<bm::bvector<> >::periodic_candidates_type groups;
            groups.push_back(group);
            auto bits = bm::periodic_write_groups(enc, groups, false, true,
                                                   true, mode);
            require(enc.size() * 8 == bits.header + bits.starts + bits.lengths + bits.padding,
                    "boundary formula accounting");
            const auto scored = bm::periodic_score_explicit_groups(groups, mode);
            require(scored.header == bits.header && scored.starts == bits.starts &&
                    scored.lengths == bits.lengths && scored.padding == bits.padding,
                    "boundary formula scoring");
            bm::bounded_decoder dec; dec.bind(buf, enc.size());
            bm::periodic_read_groups(dec, restored, true, true, mode);
            require(dec.size() == enc.size() &&
                    !std::memcmp(expected, restored, sizeof(bm::word_t)*bm::set_block_size),
                    "step-1 boundary reconstruction");
        }
    }
}

static void schedules()
{
    for (unsigned lo : {1u, 17u, 127u, 128u})
        for (unsigned mode = 1; mode <= 4; ++mode)
        {
            bm::periodic_options options; options.order = mode;
            options.min_step = lo; options.max_step = 128;
            auto steps = bm::periodic_schedule<bm::bvector<>::allocator_type>(options);
            require(steps.size() == 129 - lo, "schedule length");
            if (mode == 2) require(steps.front() == lo, "dense direction");
            if (mode == 3) require(steps.front() == 128, "sparse direction");
            if (mode == 4) require(steps.front() == (lo + 128) / 2, "middle first");
            std::sort(steps.begin(), steps.end());
            for (unsigned i = 0; i < steps.size(); ++i)
                require(steps[i] == lo + i, "schedule duplicate or missing step");
        }
}

static void serialization()
{
    using BV = bm::bvector<>;
    unsigned multiple = 0, inverted = 0;
    for (unsigned pattern = 0; pattern < 4; ++pattern)
    {
        BV source;
        for (unsigned p = 0; p < 18000; p += 5) source.set(p);
        for (unsigned p = 24001; p < 60000; p += 127) source.set(p);
        if (pattern & 1)
            for (unsigned p = 900; p < 17000; p += 395) source.set(p, false);
        if (pattern & 2) { source.resize(65536); source.invert(); }
        bm::serializer<BV>::buffer plain;
        bm::serializer<BV> baseline; baseline.set_compression_level(6);
        baseline.serialize(source, plain);
        for (unsigned mode = 1; mode <= 4; ++mode)
            for (unsigned budget : {1u, 2304u})
            {
                bm::periodic_options options; options.order = mode;
                options.max_step = 128; options.trial_budget = budget;
                bm::serializer<BV> ser; ser.set_compression_level(6);
                ser.set_periodic_options(options);
                bm::serializer<BV>::buffer blob; ser.serialize(source, blob);
                BV restored;
                auto consumed = bm::deserialize(restored, blob.buf());
                require(restored == source && consumed == blob.size(), "whole-vector round trip");
                require(blob.size() <= plain.size(), "fallback larger than baseline");
                const auto& s = ser.get_periodic_statistics();
                require(s.trials <= s.sampled * budget, "trial budget exceeded");
                require(s.saved_bytes == plain.size() - blob.size(), "saved-byte accounting");
                unsigned gain = 0;
                for (unsigned step = 1; step <= 128; ++step)
                    gain += unsigned(s.order_stats.gain_step[step]);
                require(gain == s.saved_bytes, "incremental gains mismatch");
                multiple += unsigned(s.order_stats.multi_group);
                inverted += unsigned(s.inverted_accepted);
            }
        bm::serializer<BV> reduced; reduced.set_compression_level(6);
        reduced.set_block_selection_step2_options();
        reduced.set_periodic_step1_start_reduction(true);
        bm::serializer<BV>::buffer reduced_blob;
        reduced.serialize(source, reduced_blob);
        BV reduced_restored;
        auto used = bm::deserialize(reduced_restored, reduced_blob.buf());
        require(reduced_restored == source && used == reduced_blob.size(),
                "reduced selection round trip");
        bm::serializer<BV> first; first.set_compression_level(6);
        first.set_block_selection_step2_options();
        first.set_periodic_step1_first_start(true);
        bm::serializer<BV>::buffer first_blob;
        first.serialize(source, first_blob);
        BV first_restored;
        used = bm::deserialize(first_restored, first_blob.buf());
        require(first_restored == source && used == first_blob.size(),
                "explicit-first selection round trip");
        for (unsigned mode : {1u, 4u})
        {
            bm::serializer<BV> interleaved; interleaved.set_compression_level(6);
            interleaved.set_block_selection_step2_options();
            interleaved.set_periodic_step1_interleaved(mode);
            bm::serializer<BV>::buffer packed;
            interleaved.serialize(source, packed);
            BV decoded;
            used = bm::deserialize(decoded, packed.buf());
            require(decoded == source && used == packed.size(),
                    "interleaved selection round trip");
        }
    }
    require(multiple > 0, "multi-group acceptance not exercised");
    require(inverted > 0, "polarity restoration not exercised");
    std::cout << "multi_group=" << multiple << " inverted=" << inverted << '\n';
}

int main()
{
    try { formulas(); step1_boundaries(); schedules(); serialization(); }
    catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
    std::cout << "Periodic order tests PASS\n";
}
