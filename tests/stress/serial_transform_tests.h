#ifndef BM_STRESS_SERIAL_TRANSFORM_TESTS_H
#define BM_STRESS_SERIAL_TRANSFORM_TESTS_H

#include <set>
#include "bmserial_profile.h"

// Transform oracles and direct internal block-codec tests.
namespace serial_transform_test {
static const unsigned seed = 0x71b34a29u;
static unsigned context = 0;
static const char* family = "unset";
static unsigned case_pattern = 0, case_split = 0, case_length = 0, case_step = 0, case_width = 0;
#define TRANSFORM_CHECK(expr) do { if (!(expr)) { \
    std::cerr << "Transform invariant failed: " #expr << " seed=" << seed \
              << " family=" << family << " pattern=" << case_pattern \
              << " split=" << case_split << " length=" << case_length \
              << " step=" << case_step << " width=" << case_width \
              << " case=" << context << " line=" << __LINE__ << std::endl; \
    std::abort(); } } while (0)

static unsigned ReferenceGammaBits(unsigned value)
{
    unsigned bits = 1;
    while (value >>= 1) bits += 2;
    return bits;
}

static void CheckGapCorrection(const std::vector<bm::gap_word_t>& positions)
{
    unsigned expected = 65535;
    for (size_t i = 1; i < positions.size(); ++i)
        expected = std::min(expected, unsigned(positions[i] - positions[i-1] - 1));
    std::vector<bm::gap_word_t> out(positions.size());
    const auto before = positions;
    TRANSFORM_CHECK(bm::gap_xor_correct(positions.data(), unsigned(positions.size()), out.data()) == expected);
    TRANSFORM_CHECK(positions == before);
    for (size_t i = 0; i < positions.size(); ++i) {
        TRANSFORM_CHECK(out[i] == unsigned(positions[i]) - unsigned(i)*expected);
        if (i) TRANSFORM_CHECK(out[i] > out[i-1]);
    }
    auto inplace = positions;
    TRANSFORM_CHECK(bm::gap_xor_correct(inplace.data(), unsigned(inplace.size()), inplace.data()) == expected);
    TRANSFORM_CHECK(inplace == out);
}

template<unsigned N>
static void CheckBWTOracle(bm::bwt_full_workspace<>& ws, unsigned trials)
{
    std::mt19937 random(seed + N);
    for (unsigned trial = 0; trial < trials; ++trial)
    {
        case_length = N; case_pattern = trial;
        std::vector<bm::word_t> input((N+31)/32, 0), output((N+31)/32, ~bm::word_t(0));
        if (N == 16) input[0] = trial | 0xa5a50000u;
        else for (auto& word : input) word = random();
        const auto original = input;
        std::vector<unsigned> rows(N);
        for (unsigned i = 0; i < N; ++i) rows[i] = i;
        auto bit = [&](unsigned i) { return (input[i/32] >> (i%32)) & 1u; };
        std::sort(rows.begin(), rows.end(), [&](unsigned a, unsigned b) {
            for (unsigned i = 0; i < N; ++i) {
                unsigned x = bit((a+i)%N), y = bit((b+i)%N);
                if (x != y) return x < y;
            }
            return a < b;
        });
        unsigned primary = bm::bwt_transform_bits<N>(input.data(), output.data(), ws);
        TRANSFORM_CHECK(input == original && primary < N);
        for (unsigned i = 0; i < N; ++i) {
            TRANSFORM_CHECK(((output[i/32] >> (i%32)) & 1u) == bit((rows[i]+N-1)%N));
            TRANSFORM_CHECK(bit((rows[primary]+i)%N) == bit(i));
        }
        // Decode an independently selected valid primary row, including ties.
        unsigned reference_primary = unsigned(std::find(rows.begin(), rows.end(), 0u) - rows.begin());
        bm::bwt_inverse_bits<N>(output.data(), reference_primary, ws);
        for (unsigned i = 0; i < N; ++i)
            TRANSFORM_CHECK(((output[i/32] >> (i%32)) & 1u) == bit(i));
    }
}

static void TestAdaptiveDeltaFlags()
{
    family = "AdaptiveDeltaFlags";
    std::cout << "TestAdaptiveDeltaFlags seed=" << seed << std::endl;
    bm::heap_vector<unsigned, bm::bvector<>::allocator_type, true> restored;
    restored.reserve(26); // largest flags+9 scratch size below
    const unsigned* restored_data = restored.data();
    for (unsigned flags : {0u, 1u, 7u, 8u, 9u, 15u, 16u, 17u})
    for (unsigned pattern = 0; pattern < 3; ++pattern)
    for (unsigned mode : {0u, 2u})
    for (unsigned prefix : {0u, 1u, 31u})
    {
        case_length = flags; case_pattern = pattern; case_split = mode; case_step = prefix;
        bm::heap_vector<unsigned, bm::bvector<>::allocator_type, true> corrections;
        corrections.resize(flags+1);
        corrections[0] = 3;
        for (unsigned w = 1; w <= flags; ++w)
            corrections[w] = pattern == 0 ? 3u : pattern == 1 ? (w%2 ? 17u : 3u) :
                             (w%3 ? corrections[w-1] : (corrections[w-1] == 3 ? 17u : 3u));
        unsigned expected_bits = 0;
        for (unsigned base = 1; base <= flags;) {
            if (base+7 <= flags) {
                bool all_change = true;
                for (unsigned j = 0; j < 8; ++j)
                    all_change &= corrections[base+j] != corrections[base+j-1];
                expected_bits += all_change ? 1u : 9u;
                base += 8;
            } else { expected_bits += flags-base+1; break; }
        }
        for (unsigned w = 1; w <= flags; ++w) {
            if (corrections[w] == corrections[w-1]) continue;
            unsigned delta = corrections[w] > corrections[w-1] ?
                corrections[w]-corrections[w-1] : corrections[w-1]-corrections[w];
            expected_bits += mode == 0 ? ReferenceGammaBits(corrections[w]+1) : 1+ReferenceGammaBits(delta);
        }
        unsigned char storage[512] = {};
        bm::encoder enc(storage, sizeof(storage));
        {
            bm::bit_out<bm::encoder> out(enc);
            if (prefix) out.put_bits(0, prefix);
            bm::adaptive_delta_encoder::encode_flags(out, corrections);
            bm::adaptive_delta_encoder::encode_values(out, corrections, mode);
            TRANSFORM_CHECK(enc.size()*8 + out.pending_bits() == prefix+expected_bits);
            out.put_bits(21, 5); // ensure the reader stops at the exact bit boundary
        }
        bm::decoder dec(storage);
        bm::bit_in<bm::decoder> in(dec);
        if (prefix) TRANSFORM_CHECK(in.get_bits(prefix) == 0);
        restored.resize(flags+9);
        for (unsigned w = 0; w < restored.size(); ++w)
        {
            restored[w] = ~0u; // heap_vector does not zero-initialize trivial values
        } // for w
        bm::adaptive_delta_decoder::decode(in, flags+1, 3, mode, 65535, restored);
        TRANSFORM_CHECK(restored.equal(corrections));
        TRANSFORM_CHECK(restored.data() == restored_data);
        TRANSFORM_CHECK(in.get_bits(5) == 21);
    }
}

static void CheckBWTWaves()
{
    const unsigned n = bm::bwt_wave_bits;
    bm::byte_buffer<bm::bvector<>::allocator_type> matrix(bm::bwt_matrix_bytes);
    std::vector<bm::word_t> source(bm::set_block_size), output(bm::set_block_size);
    std::mt19937 random(seed + n);
    for (unsigned pattern = 0; pattern < 10; ++pattern)
    {
        case_length = n;
        case_pattern = pattern;
        for (unsigned i = 0; i < bm::set_block_size; ++i)
        {
            source[i] = pattern == 0 ? 0u : pattern == 1 ? ~0u :
                pattern == 2 ? 0x55555555u : pattern == 3 ? 0x00010001u :
                pattern == 4 ? (i % 32 == 0 ? 1u : 0u) :
                pattern == 5 ? (i % 32 == 0 ? ~1u : ~0u) :
                pattern == 6 ? (i % 32 < 15 ? 0u : ~0u) : random();
        } // for i
        const auto original = source;
        std::vector<unsigned> rows(n);
        for (unsigned i = 0; i < n; ++i)
        {
            rows[i] = i;
        } // for i
        // Sort actual cyclic bit sequences independently of the production
        // byte matrix and merge sort. Periodic ties must retain the same primary.
        const auto bit = [&](unsigned i) { return (source[i/32] >> (i%32)) & 1u; };
        std::sort(rows.begin(), rows.end(), [&](unsigned a, unsigned b)
        {
            for (unsigned i = 0; i < n; ++i)
            {
                const unsigned x = bit((a+i)%n), y = bit((b+i)%n);
                if (x != y) return x < y;
            } // for i
            return a < b;
        });
        unsigned short primary[bm::block_waves], unpacked[bm::block_waves];
        bm::bwt_transform_block(source.data(), output.data(), matrix.data(), primary);
        TRANSFORM_CHECK(source == original);
        for (unsigned row = 0; row < n; ++row)
        {
            TRANSFORM_CHECK(((output[row/32] >> (row%32)) & 1u) == bit((rows[row]+n-1)%n));
        } // for row
        TRANSFORM_CHECK(primary[0] == unsigned(std::find(rows.begin(), rows.end(), 0u) - rows.begin()));
        unsigned char packed[bm::bwt_index_bytes];
        bm::bwt_pack_indexes(primary, packed);
        bm::bwt_unpack_indexes(packed, unpacked);
        TRANSFORM_CHECK(!::memcmp(primary, unpacked, sizeof(primary)));
        bm::bwt_inverse_block(output.data(), unpacked);
        TRANSFORM_CHECK(output == original);
    } // for pattern
}

static void TestBWTTransforms()
{
    family = "BWTTransforms"; case_pattern = case_split = case_length = case_step = case_width = 0;
    std::cout << "TestBWTTransforms seed=" << seed << std::endl;
    CheckBWTWaves();
    bm::bwt_full_workspace<> ws, other;
    // A rotation-matrix oracle for every 8-bit string, including equal rows.
    for (unsigned bits = 0; bits < 256; ++bits)
    {
        context = bits;
        std::vector<unsigned> rotations(8);
        for (unsigned i = 0; i < 8; ++i) rotations[i] = i;
        std::sort(rotations.begin(), rotations.end(), [bits](unsigned a, unsigned b) {
            for (unsigned j = 0; j < 8; ++j)
            {
                unsigned x = (bits >> ((a+j)%8)) & 1u;
                unsigned y = (bits >> ((b+j)%8)) & 1u;
                if (x != y) return x < y;
            }
            return a < b;
        });
        unsigned expected = 0;
        for (unsigned row = 0; row < 8; ++row)
            expected |= ((bits >> ((rotations[row]+7)%8)) & 1u) << row;
        bm::word_t source = bits, output = 0, repeat = 0;
        unsigned primary = bm::bwt_transform_bits<8>(&source, &output, ws);
        unsigned p2 = bm::bwt_transform_bits<8>(&source, &repeat, other);
        TRANSFORM_CHECK(source == bits && (output & 255u) == expected);
        TRANSFORM_CHECK((repeat & 255u) == expected && primary == p2);
        bm::bwt_inverse_bits<8>(&output, primary, ws);
        TRANSFORM_CHECK((output & 255u) == bits);
        for (unsigned shift = 1; shift < 8; ++shift)
        {
            bm::word_t rotated = ((bits >> shift) | (bits << (8-shift))) & 255u;
            bm::word_t col = 0;
            bm::bwt_transform_bits<8>(&rotated, &col, other);
            TRANSFORM_CHECK((col & 255u) == expected);
        }
    }
    CheckBWTOracle<16>(ws, 65536);
    CheckBWTOracle<32>(ws, 64);
    CheckBWTOracle<64>(ws, 64);
    CheckBWTOracle<128>(ws, 64);
    CheckBWTOracle<16>(ws, 32); // workspace shrink after larger transforms
    std::mt19937 rng(seed);
    const unsigned words = bm::set_block_size, n = words * 32;
    std::vector<bm::word_t> source(words + 2), dst(words + 2), copy(words), rotated(words);
    for (unsigned pattern = 0; pattern < 18; ++pattern)
    {
        context = 1000 + pattern;
        std::fill(source.begin(), source.end(), 0);
        source.front() = source.back() = 0xa5a55a5au;
        for (unsigned i = 0; i < n; ++i)
        {
            bool bit = pattern == 0 ? false : pattern == 1 ? true :
                pattern == 2 ? i == 0 : pattern == 3 ? i == n-1 :
                pattern == 4 ? i != 0 : pattern == 5 ? i != n-1 :
                pattern == 6 ? bool(i & 1u) : pattern == 7 ? (i >= 123 && i < 60000) :
                pattern == 8 ? (i % 8 < 3) : pattern == 9 ? (i % 16 < 11) :
                bool(rng() & 1u);
            if (bit) source[1 + i/32] |= bm::word_t(1) << (i%32);
        }
        std::copy(source.begin()+1, source.end()-1, copy.begin());
        std::fill(dst.begin(), dst.end(), 0x5a5aa5a5u);
        unsigned short primary = bm::bwt_transform_full(source.data()+1, dst.data()+1, ws);
        TRANSFORM_CHECK(source.front() == 0xa5a55a5au && source.back() == 0xa5a55a5au);
        TRANSFORM_CHECK(dst.front() == 0x5a5aa5a5u && dst.back() == 0x5a5aa5a5u);
        TRANSFORM_CHECK(std::equal(copy.begin(), copy.end(), source.begin()+1));
        unsigned in_count = 0, out_count = 0;
        for (unsigned i = 0; i < n; ++i) {
            in_count += (copy[i/32] >> (i%32)) & 1u;
            out_count += (dst[1+i/32] >> (i%32)) & 1u;
        }
        TRANSFORM_CHECK(in_count == out_count);
        std::vector<bm::word_t> transformed(dst.begin()+1, dst.end()-1), independent(words);
        unsigned short p2 = bm::bwt_transform_full(copy.data(), independent.data(), other);
        TRANSFORM_CHECK(independent == transformed && p2 == primary);
        for (unsigned i = 0; i < n; ++i)
            if ((copy[((i+37)%n)/32] >> ((i+37)%32)) & 1u)
                rotated[i/32] |= bm::word_t(1) << (i%32);
        bm::bwt_transform_full(rotated.data(), independent.data(), other);
        TRANSFORM_CHECK(independent == transformed);
        std::fill(rotated.begin(), rotated.end(), 0);
        bm::bwt_inverse_full(dst.data()+1, primary, ws);
        TRANSFORM_CHECK(std::equal(copy.begin(), copy.end(), dst.begin()+1));
        TRANSFORM_CHECK(dst.front() == 0x5a5aa5a5u && dst.back() == 0x5a5aa5a5u);
    }
}

/**
    \brief Test direct BWT payload decoding independently of bvector serialization.
    \ingroup bvserial
    \internal
    Checks materialization, dry traversal, payload rejection, and record boundaries
    for bm::set_block_bwt_unpacked (opcode 96). Inverse-BWT correctness itself is
    checked against independent rotation oracles by TestBWTTransforms().
*/
static void TestBWTUnpackedDecoding()
{
    family = "BWTUnpackedDecoding";
    case_pattern = case_split = case_length = case_step = case_width = 0;
    std::cout << "TestBWTUnpackedDecoding seed=" << seed << std::endl;
    const unsigned n = bm::gap_max_bits;
    const bm::word_t guard = 0xa5a55a5au, tail = 0x917bac32u;
    std::vector<bm::word_t> column(bm::set_block_size), expected(column.size());
    std::vector<bm::word_t> actual(column.size() + 2);
    std::vector<unsigned char> bytes(256 * 1024);
    std::mt19937 rng(seed);
    bm::bwt_full_workspace<> ws, decode_workspace;
    bm::heap_vector<bm::gap_word_t, bm::bvector<>::allocator_type, true> position_scratch;
    for (unsigned pattern = 0; pattern < 12; ++pattern)
    {
        case_pattern = pattern;
        std::fill(column.begin(), column.end(), 0);
        for (unsigned i = 0; i < n; ++i)
        {
            const bool bit = pattern == 0 ? false : pattern == 1 ? true :
                pattern == 2 ? i == 0 : pattern == 3 ? i == n-1 :
                pattern == 4 ? i != 0 : pattern == 5 ? i != n-1 :
                pattern == 6 ? (i == 31 || i == 32 || i == 63 || i == 64 || i == n-1) :
                pattern == 7 ? bool(i & 1u) :
                pattern == 8 ? (i >= 123 && i < 60000) :
                rng() % 100 < (pattern == 9 ? 1u : pattern == 10 ? 50u : 99u);
            if (bit)
                column[i / 32] |= bm::word_t(1) << (i % 32);
        } // for i
        const auto original = column;
        for (unsigned primary : {0u, n-1})
        {
            case_split = primary;
            expected = column;
            bm::bwt_inverse_full(expected.data(), bm::short_t(primary), ws);
            bm::encoder enc(bytes.data(), bytes.size());
            bm::bwt_bic_encoder<>::encode(enc, column.data(), primary, position_scratch);
            const size_t payload_size = enc.size();
            enc.put_32(tail);
            const auto encoded = bytes;
            TRANSFORM_CHECK(column == original);
            for (unsigned mode = 0; mode < 3; ++mode)
            {
                case_width = mode;
                std::fill(actual.begin(), actual.end(), guard);
                bm::word_t* target = mode == 0 ? actual.data() + 1 :
                    mode == 1 ? nullptr : FULL_BLOCK_FAKE_ADDR;
                bm::decoder dec(bytes.data());
                bm::bwt_bic_decoder<>::decode(dec, target, decode_workspace);
                TRANSFORM_CHECK(dec.size() == payload_size);
                TRANSFORM_CHECK(dec.get_32() == tail);
                TRANSFORM_CHECK(actual.front() == guard && actual.back() == guard);
                if (!mode)
                {
                    TRANSFORM_CHECK(std::equal(expected.begin(), expected.end(), actual.begin()+1));
                }
                else
                {
                    TRANSFORM_CHECK(std::all_of(actual.begin(), actual.end(),
                        [](bm::word_t word) { return word == guard; }));
                }
                TRANSFORM_CHECK(bytes == encoded);
            } // for mode
        } // for primary
    } // for pattern

    // Rejected headers/counts and retired formula variants must not overwrite
    // the caller's block, including when rejection follows BIC materialization.
    for (unsigned variant = 0; variant < 6; ++variant)
    {
        case_pattern = variant;
        bm::encoder enc(bytes.data(), bytes.size());
        {
            bm::bit_out<bm::encoder> out(enc);
            out.put_bit(variant == 0);
            out.put_bits(65535, 16);
            out.put_bits(variant == 1, 9);
            out.put_bit(variant == 2);
            const unsigned count = variant == 3 ? n+1 : variant == 5 ? 3u : 0u;
            out.gamma(count + 1);
            if (variant == 5)
            {
                const bm::gap_word_t positions[] = {0, 32, 65535};
                out.bic_encode_u16_cm(positions, 3, 0, 65535);
            }
            out.gamma(variant >= 4 ? 2 : 1);
            out.flush();
        }
        const size_t payload_size = enc.size();
        enc.put_32(tail);
        for (unsigned mode = 0; mode < 3; ++mode)
        {
            case_width = mode;
            std::fill(actual.begin(), actual.end(), guard);
            bm::decoder dec(bytes.data());
            bool rejected = false;
            try
            {
                bm::bwt_bic_decoder<>::decode(dec, mode == 0 ? actual.data()+1 :
                    mode == 1 ? nullptr : FULL_BLOCK_FAKE_ADDR, decode_workspace);
            }
            catch (const std::logic_error&)
            {
                rejected = true;
            }
            TRANSFORM_CHECK(rejected);
            TRANSFORM_CHECK(dec.size() <= payload_size);
            TRANSFORM_CHECK(std::all_of(actual.begin(), actual.end(),
                [](bm::word_t word) { return word == guard; }));
        } // for mode
        std::vector<unsigned char> malformed(payload_size + 3);
        malformed[0] = bm::BM_HM_NO_BO | bm::BM_HM_NO_GAPL;
        malformed[1] = bm::set_block_bwt_unpacked;
        ::memcpy(malformed.data()+2, bytes.data(), payload_size);
        malformed.back() = bm::set_block_end;
        for (unsigned api = 0; api < 2; ++api)
        {
            bm::bvector<> full; full.set_range(0, bm::gap_max_bits-1); full.optimize();
            const auto original = full;
            bool rejected = false;
            try
            {
                if (!api)
                {
                    bm::deserializer<bm::bvector<>, bm::decoder> reader;
                    reader.deserialize(full, malformed.data());
                }
                else
                {
                    bm::operation_deserializer<bm::bvector<>> reader;
                    reader.deserialize(full, malformed.data(), bm::set_OR);
                }
            }
            catch (const std::logic_error&) { rejected = true; }
            TRANSFORM_CHECK(rejected && full == original);
        } // for api
    } // for variant
}

/**
    \brief Verify OR merge independently of source serialization policy.
    \ingroup bvserial
    \internal
    \param bytes Complete forced-opcode bvector BLOB.
    \param source Independently constructed original source bits.
    \param opcode Expected top-level opcode; asserted by the profile collector.
    \param superblock True when the record covers an entire superblock.
*/
static void CheckOpcodeORMerge(const std::vector<unsigned char>& bytes,
                               const bm::bvector<>& source,
                               unsigned opcode, bool superblock)
{
    using BV = bm::bvector<>;
    bm::serialization_opcode_collector<BV> collector;
    bm::serialization_opcode_profile profile;
    collector.add_blob(bytes.data(), bytes.size(), profile);
    TRANSFORM_CHECK(profile.reconciles() && profile.opcode[opcode].count);
    const unsigned extent = superblock ? bm::set_sub_total_bits : bm::gap_max_bits;
    const std::string blob(reinterpret_cast<const char*>(bytes.data()), bytes.size());
    for (unsigned state = 0; state < 6; ++state)
    {
        case_pattern = opcode; case_width = state;
        BV initial;
        if (state == 2 || state == 3)
        {
            for (unsigned i = 0; i < bm::gap_max_bits; ++i)
            {
                if (state == 2 ? i % 257 == 0 : i % 3 != 0) initial.set(i);
            } // for i
        }
        if (state == 4) initial.set_range(100, 200);
        if (state == 5) initial.set_range(0, extent - 1);
        initial.optimize();
        if (state >= 1 && state <= 3)
            initial.get_blocks_manager().deoptimize_block(0, 0, true);
        const bm::word_t* block = initial.get_blocks_manager().get_block_ptr(0, 0);
        TRANSFORM_CHECK(state != 0 || !block);
        TRANSFORM_CHECK(state != 4 || BM_IS_GAP(block));
        TRANSFORM_CHECK(state != 5 || IS_FULL_BLOCK(block));
        if (state >= 1 && state <= 3)
            TRANSFORM_CHECK(IS_VALID_ADDR(block) && !BM_IS_GAP(block));
        // A bit beyond this record detects clearing or overwriting adjacent data.
        BV::size_type last = 0;
        source.find_reverse(last);
        const BV::size_type guard_base = (std::max)(BV::size_type(extent),
            ((last >> bm::set_block_shift) + 1) * bm::gap_max_bits);
        initial.set(guard_base + bm::gap_max_bits + 7);
        BV expected(initial); expected |= source;
        for (unsigned api = 0; api < 3; ++api)
        {
            case_split = api;
            BV actual(initial);
            // Copying may normalize an empty block; restore the intended state.
            if (state >= 1 && state <= 3)
                actual.get_blocks_manager().deoptimize_block(0, 0, true);
            bm::deserializer<BV, bm::decoder> reader;
            bm::operation_deserializer<BV> operation;
            bm::streams_deserializer<BV> streaming;
            for (unsigned repeat = 0; repeat < 2; ++repeat)
            {
                case_length = repeat;
                if (api == 0)
                    TRANSFORM_CHECK(reader.deserialize(actual, bytes.data()) == bytes.size());
                else if (api == 1)
                    operation.deserialize(actual, bytes.data(), bm::set_OR);
                else
                {
                    std::stringstream input(blob + "TRAILER");
                    bm::streams_decoder channel(input);
                    TRANSFORM_CHECK(streaming.deserialize(actual, channel));
                    TRANSFORM_CHECK(size_t(input.tellg()) == bytes.size() && input.get() == 'T');
                }
                TRANSFORM_CHECK(actual == expected); // exact union and OR idempotence
            } // for repeat
        } // for api
    } // for state
}

/** \brief Forced-record OR tests for transformed blocks and superblocks.
    \ingroup bvserial
    \internal
    Covers absent, allocated-zero, sparse/dense BIT, GAP, and FULL targets.
*/
static void TestNewOpcodeORMerge()
{
    family = "NewOpcodeORMerge";
    std::cout << "TestNewOpcodeORMerge seed=" << seed << std::endl;
    using BV = bm::bvector<>;
    std::vector<unsigned char> bytes(1024 * 1024);
    bm::heap_vector<bm::gap_word_t, BV::allocator_type, true> positions;
    bm::bwt_full_workspace<> workspace;
    bm::heap_vector<bm::word_t, BV::allocator_type, true> block, transformed;
    block.resize(bm::set_block_size); transformed.resize(bm::set_block_size);
    bm::heap_vector<bm::word_t, BV::allocator_type, true> target, scratch;
    target.resize(bm::set_block_size); scratch.resize(bm::set_block_size);
    bm::serial_iterator_workspace<BV::allocator_type> iterator_workspace;
    for (unsigned pattern = 0; pattern < 3; ++pattern)
    {
        BV source; bm::bit_block_set(block.data(), 0);
        for (unsigned i = 0; i < bm::gap_max_bits; ++i)
        {
            if (pattern == 1 || (pattern == 2 && i % 101 < 13))
            { source.set(i); bm::set_bit(block.data(), i); }
        } // for i
        for (unsigned opcode : {unsigned(bm::set_block_bwt), unsigned(bm::set_block_bwt_unpacked)})
        {
            bm::encoder enc(bytes.data(), bytes.size());
            enc.put_8(bm::BM_HM_NO_BO | bm::BM_HM_NO_GAPL);
            enc.put_8(static_cast<unsigned char>(opcode));

            const unsigned short primary = bm::bwt_transform_full(block.data(), transformed.data(), workspace);
            if (opcode == bm::set_block_bwt_unpacked)
                bm::bwt_bic_encoder<>::encode(enc, transformed.data(), primary, positions);
            else
            {
                enc.put_16(8193); enc.put_8(bm::set_block_bit);
                enc.put_32(transformed.data(), bm::set_block_size);
                enc.put_8(bm::set_block_bwt_full_post); enc.put_16(primary);
            }
            enc.put_8(bm::set_block_end);
            // Exercise the internal block OR dispatch directly on nonzero data.
            for (unsigned i = 0; i < bm::set_block_size; ++i)
            { target[i] = 0x01234567u ^ i; } // for i
            bm::serial_stream_iterator<bm::decoder, BV::block_idx_type, BV::allocator_type>
                cursor(bytes.data(), &iterator_workspace);
            cursor.next();
            TRANSFORM_CHECK(cursor.get_block_type() == opcode);
            cursor.get_bit_block(target.data(), scratch.data(), bm::set_OR);
            for (unsigned i = 0; i < bm::set_block_size; ++i)
            { TRANSFORM_CHECK(target[i] == ((0x01234567u ^ i) | block[i])); } // for i
            TRANSFORM_CHECK(cursor.decoder().get_pos() == bytes.data()+enc.size()-1);
            TRANSFORM_CHECK(cursor.decoder().get_8() == bm::set_block_end);
            CheckOpcodeORMerge(std::vector<unsigned char>(bytes.data(), bytes.data()+enc.size()), source, opcode, false);
        } // for opcode
    } // for pattern
    {
        BV source; source.set(5); source.set(7); source.set(10);
        bm::encoder enc(bytes.data(), bytes.size());
        enc.put_8(bm::BM_HM_NO_BO | bm::BM_HM_NO_GAPL);
        enc.put_8(bm::set_block_gap_xor_framed);
        enc.put_8(bm::set_block_bit_1bit); enc.put_16(5);
        {
            bm::bit_out<bm::encoder> out(enc); out.gamma(2); out.gamma(3);
            const bm::gap_word_t ids[] = {7, 8};
            out.bic_encode_u16(ids, 2, 0, 65533); out.flush();
        }
        enc.put_8(bm::set_block_end);
        CheckOpcodeORMerge(std::vector<unsigned char>(bytes.data(), bytes.data()+enc.size()),
                           source, bm::set_block_gap_xor_framed, false);
    }
    std::vector<unsigned> ids(96);
    BV wave_source;
    for (unsigned i = 0; i < ids.size(); ++i)
    { ids[i] = 17 + i * 1031; wave_source.set(ids[i]); } // for i
    for (unsigned width : {0u, 128u, 1024u, 4096u})
    for (unsigned mode = 0; mode < 4; ++mode)
    {
        const unsigned opcode = width ? bm::set_sblock_bienc_wave_adaptive : bm::set_sblock_bienc_wave1024;
        bm::encoder enc(bytes.data(), bytes.size());
        enc.put_8(bm::BM_HM_NO_BO | bm::BM_HM_NO_GAPL | bm::BM_HM_SPARSE);
        enc.put_8(static_cast<unsigned char>(opcode));
        bm::sblock_wave_record record; bm::sblock_wave_cost cost;
        {
            bm::bit_out<bm::encoder> out(enc);
            auto position = [&]() { return bm::id64_t(enc.size()) * 8 + out.pending_bits(); };
            if (width)
            {
                bm::sblock_wave_adaptive_model<> model;
                model.build(ids.data(), unsigned(ids.size()));
                model.wave_size = width; model.wave_count = bm::set_sub_total_bits / width;
                TRANSFORM_CHECK(model.prepare(ids.data(), unsigned(ids.size()), record));
                unsigned code = 0;
                while ((128u << code) != width) { ++code; } // while code
                out.put_bits(code, 3); model.encode(out, position, 0, mode, cost);
            }
            else
            {
                bm::sblock_wave1024_model<> model;
                TRANSFORM_CHECK(model.prepare(ids.data(), unsigned(ids.size()), record));
                model.encode(out, position, 0, mode, cost);
            }
        }
        enc.put_8(bm::set_block_end);
        CheckOpcodeORMerge(std::vector<unsigned char>(bytes.data(), bytes.data()+enc.size()), wave_source, opcode, true);
    } // for width, mode
    for (bool inverse : {false, true})
    {
        std::vector<unsigned> large_ids(65536), corrected;
        BV source;
        if (inverse) source.set_range(0, bm::set_sub_total_bits-1);
        for (unsigned i = 0; i < large_ids.size(); ++i)
        {
            large_ids[i] = 17 + i;
            if (inverse) source.clear_bit(large_ids[i]); else source.set(large_ids[i]);
        } // for i
        bm::encoder enc(bytes.data(), bytes.size());
        enc.put_8(bm::BM_HM_NO_BO | bm::BM_HM_NO_GAPL | bm::BM_HM_SPARSE);
        enc.put_8(bm::set_sblock_bienc_large);
        {
            bm::bit_out<bm::encoder> out(enc);
            bm::sblock_large_encoder<decltype(out)>::encode(out, large_ids.data(), unsigned(large_ids.size()), 0, 0, corrected, inverse);
        }
        enc.put_8(bm::set_block_end);
        CheckOpcodeORMerge(std::vector<unsigned char>(bytes.data(), bytes.data()+enc.size()), source,
                           bm::set_sblock_bienc_large, true);
    } // for inverse
    BV adaptive_source;
    unsigned at = 10; adaptive_source.set(at);
    for (unsigned i = 1; i < 96; ++i)
    { at += i >= 32 && i < 64 ? 2048u : 2u; adaptive_source.set(at); } // for i
    adaptive_source.optimize();
    bm::serializer<BV> writer; writer.set_compression_level(6);
    writer.set_sblock_adaptive_options(true);
    bm::serializer<BV>::buffer blob; writer.serialize(adaptive_source, blob);
    CheckOpcodeORMerge(std::vector<unsigned char>(blob.data(), blob.data()+blob.size()),
                       adaptive_source, bm::set_sblock_bienc_adaptive, true);
    std::cout << "TestNewOpcodeORMerge OK" << std::endl;
}

// Compare cancellation with a bitwise parity oracle, independent of sorting.
static void CheckGapXorCancellation()
{
    bm::heap_vector<bm::gap_word_t, bm::bvector<>::allocator_type, true> edits;
    std::vector<bm::word_t> parity(bm::set_block_size);
    std::mt19937 random(seed + 94);
    const unsigned lengths[] = {0, 1, 2, 3, 4, 7, 8, 9, 32, 1024, 4096, 65536, 65537};
    for (unsigned pattern = 0; pattern < 6; ++pattern)
    {
        for (unsigned n : lengths)
        {
            case_pattern = pattern;
            case_length = n;
            edits.resize(n);
            std::fill(parity.begin(), parity.end(), 0);
            for (unsigned i = 0; i < n; ++i)
            {
                const unsigned pos = pattern == 0 ? i % bm::gap_max_bits :
                    pattern == 1 ? (n - i - 1) % bm::gap_max_bits :
                    pattern == 2 ? 17u :
                    pattern == 3 ? (i & 1u ? 0u : bm::gap_max_bits - 1) :
                    pattern == 4 ? unsigned(random() % 32) :
                                   unsigned(random() % bm::gap_max_bits);
                edits[i] = bm::gap_word_t(pos);
                parity[pos / 32] ^= bm::word_t(1) << (pos % 32);
            } // for i
            const auto* storage = edits.data();
            const size_t capacity = edits.capacity();
            bm::gap_xor_cancel(edits);
            TRANSFORM_CHECK(edits.data() == storage && edits.capacity() == capacity);
            size_t index = 0;
            for (unsigned pos = 0; pos < bm::gap_max_bits; ++pos)
            {
                if (parity[pos / 32] & (bm::word_t(1) << (pos % 32)))
                {
                    TRANSFORM_CHECK(index < edits.size());
                    TRANSFORM_CHECK(edits[index++] == pos);
                }
            } // for pos
            TRANSFORM_CHECK(index == edits.size());
        } // for n
    } // for pattern
}

static void TestGapXorTransforms()
{
    family = "GapXorTransforms"; case_pattern = case_split = case_length = case_step = case_width = 0;
    std::cout << "TestGapXorTransforms seed=" << seed << std::endl;
    CheckGapXorCancellation();
    for (const auto& positions : std::vector<std::vector<bm::gap_word_t>>{
             {0, 65535}, {0, 1, 2, 3}, {100, 200, 300, 400}, {0, 100, 101, 65535}})
        CheckGapCorrection(positions);
    std::vector<bm::gap_word_t> edits = {65535, 0, 4, 4, 5, 5, 5, 0, 2};
    bm::gap_xor_cancel(edits);
    TRANSFORM_CHECK((edits == std::vector<bm::gap_word_t>{2, 5, 65535}));
    std::vector<bm::word_t> bits(bm::set_block_size);
    std::vector<bm::gap_word_t> gap(65540);
    // Match the serializer scratch allocation, with canaries on both ends.
    const size_t scratch_words = bm::set_block_size * bm::tmp_buff_alloc_factor *
                                 sizeof(bm::word_t) / sizeof(bm::gap_word_t);
    std::vector<bm::gap_word_t> base(scratch_words+2, 0x5a5a),
        ex0(scratch_words+2, 0x5a5a), ex1(scratch_words+2, 0x5a5a);
    std::vector<bm::gap_word_t> cleaned, positions, corrected;
    bool empty_seen = false, single_seen = false, uniform_seen = false;
    bool first_seen = false, last_seen = false;
    size_t most_exceptions = 0;
    std::mt19937 rng(seed);
    for (unsigned pattern = 0; pattern < 70; ++pattern)
    {
        context = 2000 + pattern; case_pattern = pattern;
        std::fill(bits.begin(), bits.end(), pattern & 1 ? ~bm::word_t(0) : 0);
        const unsigned count = pattern < 4 ? pattern :
                               pattern >= 68 ? 3000u : 40 + pattern * 9;
        for (unsigned j = 0; j < count; ++j)
        {
            unsigned pos = j == 0 ? 0 : j == 1 ? 65535 :
                           j == 2 ? 1 : (rng() % 65536);
            bits[pos/32] ^= bm::word_t(1) << (pos%32);
        }
        unsigned len = bm::bit_to_gap(gap.data(), bits.data(), bm::gap_max_bits);
        TRANSFORM_CHECK(len != 0);
        std::vector<bm::gap_word_t> original(gap.begin(), gap.begin() + len + 1);
        for (unsigned split = 0; split < 4; ++split)
        {
            context = 2000 + pattern * 4 + split; case_split = split; case_length = len;
            bm::gap_xor_prepare_candidate(gap.data(), split, base.data()+1, ex0.data()+1, ex1.data()+1, cleaned, positions);
            TRANSFORM_CHECK(base.front() == 0x5a5a && base.back() == 0x5a5a);
            TRANSFORM_CHECK(ex0.front() == 0x5a5a && ex0.back() == 0x5a5a);
            TRANSFORM_CHECK(ex1.front() == 0x5a5a && ex1.back() == 0x5a5a);
            empty_seen |= positions.empty();
            single_seen |= positions.size() == 1;
            uniform_seen |= cleaned.size() < 3;
            first_seen |= !positions.empty() && positions.front() == 0;
            last_seen |= !positions.empty() && positions.back() == 65535;
            most_exceptions = std::max(most_exceptions, positions.size());
            TRANSFORM_CHECK(std::equal(original.begin(), original.end(), gap.begin()));
            TRANSFORM_CHECK(std::adjacent_find(positions.begin(), positions.end(), std::greater_equal<bm::gap_word_t>()) == positions.end());
            std::vector<bm::word_t> actual(bm::set_block_size);
            bm::gap_convert_to_bitset(actual.data(), cleaned.data());
            for (auto p : positions) actual[p/32] ^= bm::word_t(1) << (p%32);
            TRANSFORM_CHECK(actual == bits);
            if (positions.size() >= 2)
            {
                CheckGapCorrection(positions);
                corrected.resize(positions.size());
                unsigned c = bm::gap_xor_correct(positions.data(), unsigned(positions.size()), corrected.data());
                for (size_t i = 0; i < positions.size(); ++i)
                    TRANSFORM_CHECK(unsigned(corrected[i]) + unsigned(i) * c == positions[i]);
                std::vector<bm::gap_word_t> inplace(positions);
                TRANSFORM_CHECK(bm::gap_xor_correct(inplace.data(), unsigned(inplace.size()), inplace.data()) == c);
                TRANSFORM_CHECK(inplace == corrected);
            }
        }
    }
    TRANSFORM_CHECK(empty_seen && single_seen && uniform_seen && first_seen && last_seen);
    TRANSFORM_CHECK(most_exceptions >= 1000);
}

static void TestAdaptiveDeltaTransforms()
{
    family = "AdaptiveDeltaTransforms"; case_pattern = case_split = case_length = case_step = case_width = 0;
    std::cout << "TestAdaptiveDeltaTransforms seed=" << seed << std::endl;
    typedef bm::superblock_delta_model model;
    TRANSFORM_CHECK(!model::eligible(~0u, (~0u / 2) + 1));
    bm::heap_vector<unsigned, bm::bvector<>::allocator_type, true> corrections;
    corrections.reserve(17); // ceil(257 / 16), covering all cases below
    const unsigned* corrections_data = corrections.data();
    for (unsigned step : {16u, 32u, 64u, 128u})
    for (unsigned n : {1u, 2u, 15u, 16u, 17u, 31u, 32u, 33u, 63u, 64u, 65u, 127u, 128u, 129u, 255u, 256u, 257u})
    for (unsigned pattern = 0; pattern < 5; ++pattern)
    for (unsigned global : {0u, 1u})
    {
        context = 3000 + step + n*4 + pattern;
        case_step = step; case_length = n; case_pattern = pattern; case_split = global;
        TRANSFORM_CHECK(model::eligible(n, step) == (n >= 64 && n >= 2*step));
        std::vector<unsigned> values(n), corrected(n), restored(n), inplace;
        values[0] = pattern == 3 ? 0xffff0000u : 7u;
        for (unsigned i = 1; i < n; ++i) {
            unsigned w = i / step;
            unsigned c = pattern == 0 ? 3u : pattern == 1 ? w+1 :
                         pattern == 2 ? 20u - w : (w%2 ? 20u : 1u);
            values[i] = values[i-1] + (pattern == 4 ? (i%step == 0 ? 2u : 23u) : c + 1 + (i%11 == 0 ? 2u : 0u));
        }
        bm::id64_t scores[3];
        unsigned changes = model::prepare(values.data(), n, step, global, corrections, scores);
        TRANSFORM_CHECK(corrections.data() == corrections_data);
        TRANSFORM_CHECK(corrections[0] == global);
        unsigned expected_changes = 0, total = 0;
        bm::id64_t absolute_score = 0, relative_score = 0;
        for (unsigned w = 1; w < corrections.size(); ++w) {
            unsigned min_gap = ~0u;
            for (unsigned i = w*step; i < std::min((w+1)*step,n); ++i)
                min_gap = std::min(min_gap, values[i]-values[i-1]);
            TRANSFORM_CHECK(corrections[w] == min_gap-1);
            expected_changes += corrections[w] != corrections[w-1];
            if (corrections[w] != corrections[w-1]) {
                unsigned delta = corrections[w] > corrections[w-1] ?
                    corrections[w] - corrections[w-1] : corrections[w-1] - corrections[w];
                absolute_score += ReferenceGammaBits(corrections[w] + 1);
                relative_score += 1 + ReferenceGammaBits(delta);
            }
        }
        TRANSFORM_CHECK(changes == expected_changes);
        TRANSFORM_CHECK(scores[0] == absolute_score && scores[1] == 0 && scores[2] == relative_score);
        for (unsigned i = 1; i < n; ++i) total += corrections[i/step];
        TRANSFORM_CHECK(model::transform(values.data(), n, step, corrections, corrected.data()) == total);
        unsigned prefix = 0;
        TRANSFORM_CHECK(corrected[0] == values[0]);
        for (unsigned i = 1; i < n; ++i) {
            prefix += corrections[i/step];
            TRANSFORM_CHECK(corrected[i] == values[i] - prefix);
            TRANSFORM_CHECK(corrected[i] > corrected[i-1]);
        }
        TRANSFORM_CHECK(model::restore(corrected.data(), n, step, corrections, restored.data()) == total);
        TRANSFORM_CHECK(restored == values);
        inplace = values;
        model::transform(inplace.data(), n, step, corrections, inplace.data());
        TRANSFORM_CHECK(inplace == corrected);
        model::restore(inplace.data(), n, step, corrections, inplace.data());
        TRANSFORM_CHECK(inplace == values);
    }

    // A small record reserves the full decoder workspace. Growing and shrinking
    // records must reuse it, including the largest legal window count.
    struct scratch_reader : bm::deseriaizer_base<bm::decoder, unsigned>
    {
        using bm::deseriaizer_base<bm::decoder, unsigned>::read_bic_sb_arr;
        using bm::deseriaizer_base<bm::decoder, unsigned>::sblock_adaptive_corrections_;
        using bm::deseriaizer_base<bm::decoder, unsigned>::sblock_adaptive_flags_;
    } reader;
    corrections.reserve(bm::sblock_adaptive_max_windows);
    const unsigned* correction_storage = nullptr;
    const unsigned char* flag_storage = nullptr;
    for (unsigned n : {64u, 65535u, 257u, 96u, 64u})
    for (unsigned si : {0u, 1u})
    {
        const unsigned step = bm::sblock_adaptive_min_step << si;
        std::vector<unsigned> values(n), corrected(n), restored(n);
        values[0] = 7;
        for (unsigned i = 1; i < n; ++i)
        {
            values[i] = values[i-1] + ((i / step) % 2 ? 4u : 2u);
        } // for i
        bm::id64_t scores[3];
        model::prepare(values.data(), n, step, 1, corrections, scores);
        model::transform(values.data(), n, step, corrections, corrected.data());
        const unsigned mode = si ? 0u : 2u;
        std::vector<unsigned char> bytes(n*8 + 1024);
        bm::encoder enc(bytes.data(), bytes.size());
        {
            bm::bit_out<bm::encoder> out(enc);
            out.put_32_no(7); out.gamma(1);
            out.put_bits(si, 2); out.put_bits(mode, 2); out.gamma(n);
            out.put_24_no(corrected[0]); out.put_24_no(corrected[n-1]);
            out.gamma(2); // global correction + 1
            bm::adaptive_delta_encoder::encode_flags(out, corrections);
            bm::adaptive_delta_encoder::encode_values(out, corrections, mode);
            out.flush_if_full();
            out.bic_encode_u32_cm(corrected.data()+1, n-2, corrected[0]+1, corrected[n-1]-1);
        }
        bm::decoder dec(bytes.data());
        unsigned sb = 0;
        TRANSFORM_CHECK(reader.read_bic_sb_arr(dec, bm::set_sblock_bienc_adaptive,
                                             restored.data(), &sb) == n);
        TRANSFORM_CHECK(sb == 7 && restored == values);
        if (!correction_storage)
        {
            correction_storage = reader.sblock_adaptive_corrections_.data();
            flag_storage = reader.sblock_adaptive_flags_.data();
        }
        TRANSFORM_CHECK(reader.sblock_adaptive_corrections_.data() == correction_storage);
        TRANSFORM_CHECK(reader.sblock_adaptive_flags_.data() == flag_storage);
        TRANSFORM_CHECK(reader.sblock_adaptive_corrections_.capacity() >= bm::sblock_adaptive_max_windows);
        TRANSFORM_CHECK(reader.sblock_adaptive_flags_.capacity() >= bm::sblock_adaptive_max_windows);
    } // for si, n
}

static void TestWaveCompaction()
{
    family = "WaveCompaction"; case_pattern = case_split = case_length = case_step = case_width = 0;
    std::cout << "TestWaveCompaction seed=" << seed << std::endl;
    bm::sblock_wave_adaptive_model<> adaptive;
    bm::sblock_wave1024_model<> fixed;
    std::mt19937 rng(seed);
    for (unsigned width = 128; width <= 4096; width *= 2)
    for (unsigned pattern = 0; pattern < 12; ++pattern)
    {
        context = 4000 + width + pattern; case_width = width; case_pattern = pattern;
        std::set<unsigned> sample;
        if (pattern == 0) { sample.insert(1); sample.insert(width-1); }
        else if (pattern == 1) {
            for (unsigned i = 0; i < 40; ++i) {
                unsigned p = i * width;
                if (p) sample.insert(p-1);
                sample.insert(p); sample.insert(p+1);
            }
            sample.insert((1u<<24)-1);
        } else if (pattern < 5) {
            unsigned count = (1u<<24)/width;
            unsigned missing = pattern == 2 ? 0 : pattern == 3 ? count/2 : count-1;
            for (unsigned w = 0; w < count; ++w) if (w != missing) sample.insert(w*width);
        } else for (unsigned i = 0; i < 300; ++i) sample.insert(rng() % (1u<<24));
        std::vector<unsigned> original(sample.begin(), sample.end()), packed, occupied;
        std::set<unsigned> waves;
        for (unsigned x : original) waves.insert(x/width);
        occupied.assign(waves.begin(), waves.end());
        for (unsigned x : original) {
            unsigned rank = unsigned(std::lower_bound(occupied.begin(), occupied.end(), x/width) - occupied.begin());
            packed.push_back(rank*width + x%width);
        }
        unsigned correction = ~0u;
        for (size_t i = 1; i < packed.size(); ++i)
            correction = std::min(correction, packed[i]-packed[i-1]-1);
        adaptive.wave_size = width; adaptive.wave_count = (1u<<24)/width;
        adaptive.build(original.data(), unsigned(original.size()));
        bm::sblock_wave_record record;
        TRANSFORM_CHECK(adaptive.prepare(original.data(), unsigned(original.size()), record));
        TRANSFORM_CHECK(record.eligible);
        std::vector<unsigned> expected_empty;
        for (unsigned w = 0; w < adaptive.wave_count; ++w)
            if (!waves.count(w)) expected_empty.push_back(w);
        TRANSFORM_CHECK(std::equal(adaptive.empty.begin(), adaptive.empty.end(), expected_empty.begin(), expected_empty.end()));
        TRANSFORM_CHECK(std::equal(adaptive.occupied.begin(), adaptive.occupied.end(), occupied.begin(), occupied.end()) && adaptive.correction == correction);
        TRANSFORM_CHECK(adaptive.upper == occupied.size()*width-1-(original.size()-1)*correction);
        for (size_t i = 0; i < original.size(); ++i) {
            TRANSFORM_CHECK(adaptive.values[i] == packed[i] - i*correction);
            TRANSFORM_CHECK(bm::sblock_wave_restore_position(packed[i], width, occupied) == original[i]);
            if (i) TRANSFORM_CHECK(adaptive.values[i] > adaptive.values[i-1]);
        }
        if (width == 1024) {
            bm::sblock_wave_record fixed_record;
            TRANSFORM_CHECK(fixed.prepare(original.data(), unsigned(original.size()), fixed_record));
            TRANSFORM_CHECK(fixed.occupied.equal(adaptive.occupied) && fixed.empty.equal(adaptive.empty));
            TRANSFORM_CHECK(fixed.values.equal(adaptive.values) && fixed.correction == adaptive.correction && fixed.upper == adaptive.upper);
        }
        std::vector<unsigned> all(adaptive.wave_count);
        for (unsigned i = 0; i < adaptive.wave_count; ++i) all[i] = i*width;
        adaptive.build(all.data(), unsigned(all.size()));
        TRANSFORM_CHECK(!adaptive.prepare(all.data(), unsigned(all.size()), record));
        TRANSFORM_CHECK(!record.eligible);
        if (width == 1024) {
            bm::sblock_wave_record fixed_record;
            TRANSFORM_CHECK(!fixed.prepare(all.data(), unsigned(all.size()), fixed_record));
        }
    }
}

static void TestPeriodicBoundary()
{
    family = "PeriodicBoundary"; case_pattern = case_split = case_length = case_step = case_width = 0;
    context = 5000;
    std::vector<bm::word_t> block(bm::set_block_size, ~bm::word_t(0));
    bm::periodic_candidate<> candidate;
    TRANSFORM_CHECK(bm::periodic_find(block.data(), false, 1, 0, 1, 4, candidate));
    TRANSFORM_CHECK(candidate.chains.size() == 1 && candidate.chains[0].length == 65536 && candidate.covered == 65536);
    bm::periodic_apply(block.data(), candidate);
    TRANSFORM_CHECK(std::count(block.begin(), block.end(), bm::word_t(0)) == bm::set_block_size);
    unsigned char bytes[64] = {};
    bm::encoder enc(bytes, sizeof(bytes));
    bm::periodic_write_formula(enc, candidate, false);
    bm::decoder dec(bytes);
    bm::periodic_read_formula(dec, block.data());
    TRANSFORM_CHECK(std::count(block.begin(), block.end(), ~bm::word_t(0)) == bm::set_block_size);
}
static void TestSuperblockEndIndexGather()
{
    family = "SuperblockEndIndexGather";
    std::cout << "TestSuperblockEndIndexGather" << std::endl;
    using BV = bm::bvector<>;
    using Index = bm::deserialization_index<BV>;
    const unsigned codes[] = {68, 90, 92, 97, 98};
    const BV::size_type groups[] = {0, 2, bm::id_max / bm::set_sub_total_bits};
    for (unsigned kind = 0; kind < 5; ++kind)
    for (bool bookmarks : {false, true})
    {
        // Opcode 98 selection explicitly excludes bookmarked serialization.
        if (kind == 4 && bookmarks) continue;
        case_pattern = kind; case_split = unsigned(bookmarks);
        BV source;
        std::vector<BV::size_type> sample;
        const unsigned count = kind == 4 ? 70000u : kind == 2 ? 4093u :
                               (kind == 1 || kind == 3) ? 100u : 2000u;
        for (auto group : groups)
        {
            unsigned at = 10;
            for (unsigned i = 0; i < count; ++i)
            {
                const unsigned local = kind == 0 ? 1+i*5003 :
                    (kind == 1 || kind == 3) ? 1+i*163841 : kind == 2 ? at : 1+i*3;
                const BV::size_type pos = group * bm::set_sub_total_bits + local;
                source.set(pos);
                if (i == 0 || i == count/2 || i == count-1) sample.push_back(pos);
                at += ((i+1)/32)%2 ? 96u : 2u;
            }
        }
        // A dense intervening block makes bookmark distances exceed the
        // writer's 512-byte threshold, so bookmarked cases contain real jumps.
        std::mt19937 noise(seed);
        for (unsigned i = 0; i < bm::gap_max_bits; ++i)
            if (noise() & 1u) source.set(BV::size_type(bm::set_sub_total_bits) + i);
        source.optimize();
        bm::serializer<BV> writer;
        writer.set_compression_level(kind ? 10 : 6);
        if (kind < 4) writer.set_sblock_large_options(false);
        if (kind == 1 || kind == 2) {
            writer.set_sblock_wave_adaptive_options(false);
            writer.set_sblock_wave1024_options(kind == 1);
            writer.set_sblock_adaptive_options(kind == 2);
        }
        writer.set_bookmarks(bookmarks, 4);
        bm::serializer<BV>::buffer bytes;
        writer.serialize(source, bytes);
        bm::serialization_opcode_collector<BV> collector;
        bm::serialization_opcode_profile profile;
        collector.set_collect_ids();
        unsigned penultimate = 0, last = 0;
        collector.add_blob(bytes.data(), bytes.size(), profile,
            [&](const unsigned char* record, size_t) { penultimate = last; last = *record; });
        std::cout << "  level=" << (kind ? 10 : 6) << " bookmarks=" << bookmarks << " opcodes:";
        for (unsigned op = 0; op < 256; ++op)
            if (profile.opcode[op].count) std::cout << " " << op << ":" << profile.opcode[op].count;
        std::cout << std::endl;
        TRANSFORM_CHECK(profile.reconciles() && profile.cardinalities_reconcile());
        TRANSFORM_CHECK(profile.ids == source.count());
        TRANSFORM_CHECK(profile.opcode[codes[kind]].count >= 3);
        TRANSFORM_CHECK(penultimate == codes[kind] && last == bm::set_block_end);
        TRANSFORM_CHECK(profile.opcode[bm::set_block_end].count == 1 && profile.opcode[bm::set_block_end].bytes == 1);

        Index ram_index, stream_index, persisted;
        BV dry;
        bm::deserializer<BV, bm::decoder>::block_starts_vect_type ram_starts, stream_starts;
        bm::deserializer<BV, bm::decoder> ram;
        ram.set_deserialization_index_construct(&ram_index);
        ram.set_record_block_starts(&ram_starts);
        TRANSFORM_CHECK(ram.deserialize(dry, bytes.data()) == bytes.size());
        TRANSFORM_CHECK(dry.none());
        const std::string blob(reinterpret_cast<const char*>(bytes.data()), bytes.size());
        std::stringstream stream(blob + "TRAILER");
        bm::streams_decoder channel(stream);
        bm::streams_deserializer<BV> streaming;
        streaming.set_deserialization_index_construct(&stream_index);
        streaming.set_record_block_starts(&stream_starts);
        TRANSFORM_CHECK(streaming.deserialize(dry, channel));
        TRANSFORM_CHECK(dry.none() && size_t(stream.tellg()) == bytes.size());
        TRANSFORM_CHECK(stream.get() == 'T');
        TRANSFORM_CHECK(ram_index.equal(stream_index));
        TRANSFORM_CHECK(ram_index.serialized_size() == bytes.size());
        TRANSFORM_CHECK(ram_starts.equal(stream_starts) && ram_starts.size() == ram_index.size());
        TRANSFORM_CHECK(ram_starts.back() == bm::set_total_blocks);
        TRANSFORM_CHECK(ram_index.has_bookmarks() == bookmarks);
        if (bookmarks) {
            BV::block_idx_type origin, target;
            size_t target_offset;
            TRANSFORM_CHECK(ram_index.find_bookmark_target(
                bm::set_total_blocks-1, origin, target, target_offset));
            TRANSFORM_CHECK(target > origin && target_offset < bytes.size());
        }
        size_t cursor = 0, ordinal = 0, offset = 0;
        while (ram_index.find_next_marker_offset(cursor, ordinal, offset)) cursor = offset;
        TRANSFORM_CHECK(cursor == bytes.size()-1);
        bm::deserialization_index_serializer<BV> index_writer;
        bm::deserialization_index_serializer<BV>::buffer index_bytes;
        index_writer.serialize(ram_index, index_bytes);
        bm::deserialization_index_deserializer<BV> index_reader;
        TRANSFORM_CHECK(index_reader.deserialize(persisted, index_bytes.data(), index_bytes.size()) == index_bytes.size());
        TRANSFORM_CHECK(persisted.equal(ram_index));

        // First/middle/last positions, scattered selection, empty selection,
        // and a requested block in an entirely empty intervening superblock.
        for (unsigned selection = 0; selection < sample.size()+3; ++selection)
        {
            case_length = selection;
            BV mask, digest, expected(source), actual, streamed;
            if (selection < sample.size()) mask.set(sample[selection]);
            else if (selection == sample.size()) for (auto pos : sample) mask.set(pos);
            else if (selection == sample.size()+2) mask.set(BV::size_type(3)*bm::set_sub_total_bits+7);
            mask.build_block_digest(digest);
            expected &= mask;
            bm::deserializer<BV, bm::decoder> reader;
            reader.set_deserialization_index_use(&persisted);
            reader.set_block_digest_vector_use(&digest);
            TRANSFORM_CHECK(reader.deserialize(actual, bytes.data()) == bytes.size());
            BV extra(actual); extra -= source; TRANSFORM_CHECK(extra.none());
            // Whole-superblock materialization is allowed; gather applies its mask.
            actual &= mask; TRANSFORM_CHECK(actual == expected);
            std::stringstream gathered_input(blob + "TRAILER");
            bm::streams_decoder gathered_channel(gathered_input);
            bm::streams_deserializer<BV> gathered_reader;
            gathered_reader.set_deserialization_index_use(&persisted);
            gathered_reader.set_block_digest_vector_use(&digest);
            TRANSFORM_CHECK(gathered_reader.deserialize(streamed, gathered_channel));
            extra = streamed; extra -= source; TRANSFORM_CHECK(extra.none());
            streamed &= mask; TRANSFORM_CHECK(streamed == expected);
            TRANSFORM_CHECK(size_t(gathered_input.tellg()) == bytes.size() && gathered_input.get() == 'T');
        }
    }
    std::cout << "TestSuperblockEndIndexGather OK" << std::endl;
}


// Mixed single-block transforms separated by zero runs and a superblock
// boundary. The profiles are assertions: fallback codecs cannot silently turn
// this into a test of ordinary serialization only.
static void TestTransformedBlockIndexGather()
{
    family = "TransformedBlockIndexGather";
    std::cout << "TestTransformedBlockIndexGather" << std::endl;
    using BV = bm::bvector<>;
    using Index = bm::deserialization_index<BV>;
    const unsigned blocks[] = {0, 1, 2, 7, 8, 9, 256, 257, 258, 4, 260};
    BV source;
    for (unsigned n = 0; n < 11; ++n)
    {
        const BV::size_type base = BV::size_type(blocks[n]) * bm::gap_max_bits;
        const unsigned kind = n < 9 ? n % 3 : 3;
        // 72: periodic BIT block; 94: irregular runs with isolated XOR
        // exceptions; 96: the same run model with a larger residual.
        std::mt19937 rng(kind == 1 ? 123u : kind == 2 ? 135u : seed);
        if (kind == 1 || kind == 2)
        {
            for (unsigned i = 0; i < bm::gap_max_bits;)
            {
                unsigned length = 100 + rng()%3000;
                if (rng() & 1u)
                    source.set_range(base+i, base+(std::min)(65535u, i+length));
                i += length+1;
            }
            for (unsigned i = 0; i < (kind == 1 ? 50u : 650u); ++i)
                source.flip(base+rng()%65536);
        }
        else
            for (unsigned i = 0; i < bm::gap_max_bits; ++i)
                if (kind == 0 ? i%17 < 8 : (rng() & 1u)) source.set(base+i);
    }
    source.optimize();
    for (bool bookmarks : {false, true})
    {
        case_split = unsigned(bookmarks);
        bm::serializer<BV> writer;
        writer.set_compression_level(10);
        writer.set_sblock_large_options(false);
        writer.set_bookmarks(bookmarks, 4);
        bm::serializer<BV>::buffer bytes;
        writer.serialize(source, bytes);
        bm::serialization_opcode_collector<BV> collector;
        bm::serialization_opcode_profile profile;
        collector.set_collect_ids();
        collector.add_blob(bytes.data(), bytes.size(), profile);
        std::cout << "  level=10 bookmarks=" << bookmarks << " opcodes:";
        for (unsigned op = 0; op < 256; ++op)
            if (profile.opcode[op].count) std::cout << " " << op << ":" << profile.opcode[op].count;
        std::cout << std::endl;
        TRANSFORM_CHECK(profile.reconciles() && profile.cardinalities_reconcile());
        TRANSFORM_CHECK(profile.ids == source.count());
        for (unsigned op : {72u, 94u, 96u}) TRANSFORM_CHECK(profile.opcode[op].count == 3);
        if (!bookmarks)
            CheckOpcodeORMerge(std::vector<unsigned char>(bytes.data(), bytes.data()+bytes.size()),
                               source, bm::set_block_bwt_unpacked, false);

        const std::string blob(reinterpret_cast<const char*>(bytes.data()), bytes.size());
        BV full, streamed_full;
        bm::deserializer<BV, bm::decoder> full_reader;
        TRANSFORM_CHECK(full_reader.deserialize(full, bytes.data()) == bytes.size());
        TRANSFORM_CHECK(full == source);
        std::stringstream full_input(blob + "TRAILER");
        bm::streams_decoder full_channel(full_input);
        bm::streams_deserializer<BV> full_stream_reader;
        TRANSFORM_CHECK(full_stream_reader.deserialize(streamed_full, full_channel));
        TRANSFORM_CHECK(streamed_full == source);
        TRANSFORM_CHECK(size_t(full_input.tellg()) == bytes.size() && full_input.get() == 'T');

        Index ram_index, stream_index, persisted;
        BV dry;
        bm::deserializer<BV, bm::decoder>::block_starts_vect_type ram_starts, stream_starts;
        bm::deserializer<BV, bm::decoder> ram;
        ram.set_deserialization_index_construct(&ram_index);
        ram.set_record_block_starts(&ram_starts);
        TRANSFORM_CHECK(ram.deserialize(dry, bytes.data()) == bytes.size() && dry.none());
        std::stringstream index_input(blob + "TRAILER");
        bm::streams_decoder index_channel(index_input);
        bm::streams_deserializer<BV> index_reader;
        index_reader.set_deserialization_index_construct(&stream_index);
        index_reader.set_record_block_starts(&stream_starts);
        TRANSFORM_CHECK(index_reader.deserialize(dry, index_channel) && dry.none());
        TRANSFORM_CHECK(size_t(index_input.tellg()) == bytes.size() && index_input.get() == 'T');
        TRANSFORM_CHECK(ram_index.equal(stream_index));
        TRANSFORM_CHECK(ram_index.serialized_size() == bytes.size());
        TRANSFORM_CHECK(ram_starts.equal(stream_starts) && ram_starts.size() == ram_index.size());
        for (unsigned block : blocks)
            TRANSFORM_CHECK(std::find(ram_starts.begin(), ram_starts.end(), block) != ram_starts.end());
        TRANSFORM_CHECK(ram_index.has_bookmarks() == bookmarks);
        if (bookmarks)
        {
            BV::block_idx_type origin, target;
            size_t offset;
            TRANSFORM_CHECK(ram_index.find_bookmark_target(260, origin, target, offset));
            TRANSFORM_CHECK(target > origin && offset < bytes.size());
        }
        bm::deserialization_index_serializer<BV> persist_writer;
        bm::deserialization_index_serializer<BV>::buffer index_bytes;
        persist_writer.serialize(ram_index, index_bytes);
        bm::deserialization_index_deserializer<BV> persist_reader;
        TRANSFORM_CHECK(persist_reader.deserialize(persisted, index_bytes.data(), index_bytes.size()) == index_bytes.size());
        TRANSFORM_CHECK(persisted.equal(ram_index));

        // Each block alone exercises both skipping and decoding every codec.
        // Also select all transformed blocks, empty blocks, no blocks, and a
        // pair straddling the superblock boundary. Exact block-level equality
        // detects accidental materialization of unrequested single-block data.
        for (unsigned selection = 0; selection < 15; ++selection)
        {
            case_length = selection;
            BV mask, digest, expected(source);
            auto select_block = [&](unsigned b) {
                const BV::size_type base = BV::size_type(b)*bm::gap_max_bits;
                mask.set_range(base, base+bm::gap_max_bits-1);
            };
            if (selection < 11) select_block(blocks[selection]);
            else if (selection == 11) for (unsigned b : blocks) select_block(b);
            else if (selection == 12) { select_block(3); select_block(259); }
            else if (selection == 14) { select_block(255); select_block(256); }
            mask.build_block_digest(digest);
            expected &= mask;
            // Exercise both freshly constructed and persisted indexes.
            for (bool reloaded : {false, true})
            {
                case_pattern = unsigned(reloaded);
                BV actual, streamed;
                bm::deserializer<BV, bm::decoder> reader;
                reader.set_deserialization_index_use(reloaded ? &persisted : &ram_index);
                reader.set_block_digest_vector_use(&digest);
                TRANSFORM_CHECK(reader.deserialize(actual, bytes.data()) == bytes.size());
                TRANSFORM_CHECK(actual == expected);
                std::stringstream input(blob + "TRAILER");
                bm::streams_decoder channel(input);
                bm::streams_deserializer<BV> stream_reader;
                stream_reader.set_deserialization_index_use(reloaded ? &persisted : &ram_index);
                stream_reader.set_block_digest_vector_use(&digest);
                TRANSFORM_CHECK(stream_reader.deserialize(streamed, channel));
                TRANSFORM_CHECK(streamed == expected);
                TRANSFORM_CHECK(size_t(input.tellg()) == bytes.size() && input.get() == 'T');
            }
        }
    }
    std::cout << "TestTransformedBlockIndexGather OK" << std::endl;
}


// Cutoff zero must disable grouping without bypassing the block-statistics
// survey needed by ordinary BIT and GAP encoders.
static void TestZeroSparseCutoff()
{
    family = "ZeroSparseCutoff";
    using BV = bm::bvector<>;
    for (bool optimized : {false, true})
    for (unsigned cutoff : {0u, 256u})
    {
        case_pattern = unsigned(optimized); case_length = cutoff;
        BV source;
        for (unsigned block : {0u, 1u, 256u, 257u})
            for (unsigned i = 0; i < 100; i += 2)
                source.set(BV::size_type(block)*bm::gap_max_bits+i);
        if (optimized) source.optimize();
        bm::serializer<BV> writer;
        writer.set_compression_level(6);
        writer.set_sparse_cutoff(cutoff);
        bm::serializer<BV>::buffer bytes;
        writer.serialize(source, bytes);
        bm::serialization_opcode_collector<BV> collector;
        bm::serialization_opcode_profile profile;
        collector.add_blob(bytes.data(), bytes.size(), profile);
        TRANSFORM_CHECK(profile.reconciles());
        TRANSFORM_CHECK(profile.opcode[bm::set_sblock_bienc_v3].count == (cutoff ? 2u : 0u));
        BV actual;
        bm::deserializer<BV, bm::decoder> reader;
        TRANSFORM_CHECK(reader.deserialize(actual, bytes.data()) == bytes.size());
        TRANSFORM_CHECK(actual == source);
    }
    std::cout << "TestZeroSparseCutoff OK" << std::endl;
}

// Profile real bit-plane BLOBs from an integer sparse vector, then exercise
// the sparse-vector (rather than bare-bvector) gather and index path.
static void TestSparseVectorOpcodeIndexGather(bool xor_refs = false)
{
    family = xor_refs ? "SparseVectorXorOpcodeIndexGather" : "SparseVectorOpcodeIndexGather";
    std::cout << family << std::endl;
    using BV = bm::bvector<>;
    using SV = bm::sparse_vector<unsigned, BV>;
    using Index = bm::sparse_vector_deserialization_index<BV>;
    SV source(bm::no_null);
    source.resize(SV::size_type(3) * bm::set_sub_total_bits);
    BV plane_bits[4];
    BV* planes[4] = {&plane_bits[0], &plane_bits[1], &plane_bits[2], &plane_bits[3]};
    const unsigned blocks[3][3] = {{0, 7, 256}, {1, 8, 257}, {2, 9, 258}};
    std::vector<SV::size_type> requested;
    for (unsigned kind = 0; kind < 3; ++kind)
    for (unsigned n = 0; n < 3; ++n)
    {
        BV& plane = *planes[kind];
        const SV::size_type base = SV::size_type(blocks[kind][n]) * bm::gap_max_bits;
        std::mt19937 rng(kind == 1 ? 123u : kind == 2 ? 135u : seed);
        if (kind == 2 && xor_refs)
        {
            rng.seed(231);
            for (unsigned i = 0; i < bm::gap_max_bits; ++i)
                if (rng()%100 < 11) plane.set(base+i);
        }
        else if (kind == 1 || kind == 2)
        {
            for (unsigned i = 0; i < bm::gap_max_bits;)
            {
                unsigned length = 100 + rng()%3000;
                if (rng() & 1u)
                    plane.set_range(base+i, base+(std::min)(65535u, i+length));
                i += length+1;
            }
            for (unsigned i = 0; i < (kind == 1 ? 50u : 650u); ++i)
                plane.flip(base+rng()%65536);
        }
        else
            for (unsigned i = 0; i < bm::gap_max_bits; ++i)
                if (i%17 < 8) plane.set(base+i);
        auto e = plane.get_enumerator(base);
        TRANSFORM_CHECK(e.valid() && *e < base+bm::gap_max_bits);
        requested.push_back(*e);
    }
    // Independent sparse plane spanning two superblocks, interleaved with
    // the three transformed planes above.
    for (unsigned group : {0u, 2u})
        for (unsigned i = 0; i < 2000; ++i)
            planes[3]->set(SV::size_type(group)*bm::set_sub_total_bits + 1+i*5003);
    // Keep the BWT plane above the sparse-superblock population threshold
    // in both occupied groups, so its individual block records stay visible.
    std::mt19937 filler(seed);
    for (unsigned block : {4u, 260u})
        for (unsigned i = 0; i < bm::gap_max_bits; ++i)
            if (filler()%100 < 85)
                planes[0]->set(SV::size_type(block)*bm::gap_max_bits+i);
    requested.push_back(1);
    requested.push_back(SV::size_type(2)*bm::set_sub_total_bits+1);
    if (xor_refs) {
        // An inverted reference must use a one-block full residual, never a
        // run marker extending into subsequent blocks.
        planes[0]->set_range(SV::size_type(6)*bm::gap_max_bits,
                            SV::size_type(7)*bm::gap_max_bits-1);
        requested.push_back(SV::size_type(6)*bm::gap_max_bits+17);
    }
    // Assemble through the sparse-vector API so its active-slice mask stays
    // coherent with the plane contents.
    BV occupied;
    for (const BV& plane : plane_bits) occupied |= plane;
    for (auto e = occupied.first(); e.valid(); ++e)
    {
        unsigned value = 0;
        for (unsigned p = 0; p < 4; ++p)
            value |= unsigned(planes[p]->test(*e)) << p;
        source.set(*e, value);
    }
    if (xor_refs)
    {
        BV blocks_digest, transformed_positions;
        for (unsigned p = 0; p < 3; ++p) transformed_positions |= plane_bits[p];
        transformed_positions.build_block_digest(blocks_digest);
        std::mt19937 reference(seed+1);
        // A common high plane makes low planes XOR-related. Keep sparse
        // plane 3 unchanged to verify that attached references disable groups.
        for (auto e = blocks_digest.first(); e.valid(); ++e)
            for (unsigned i = 0; i < bm::gap_max_bits; ++i)
            {
                auto pos = SV::size_type(*e)*bm::gap_max_bits+i;
                unsigned value = source.get(pos);
                if (reference() & 1u) value = (value ^ 7u) | (1u << 31);
                source.set(pos, value);
            }
    }
    source.optimize();

    for (bool bookmarks : {false, true})
    {
        case_split = unsigned(bookmarks);
        bm::sparse_vector_serializer<SV> writer;
        writer.set_xor_ref(xor_refs); // Explicit policy, independent of BMXORCOMP.
        writer.get_bv_serializer().set_compression_level(10);
        if (!xor_refs) writer.get_bv_serializer().set_sblock_large_options(false);
        writer.set_bookmarks(bookmarks, 4);
        bm::sparse_vector_serial_layout<SV> layout;
        writer.serialize(source, layout);
        TRANSFORM_CHECK(layout.size() > 0);
        Index index, stream_index, persisted;
        bm::sparse_vector_deserializer<SV> index_builder;
        index_builder.construct_deserialization_index(index, layout.data());
        TRANSFORM_CHECK(index.has_plane_offsets());
        bm::serialization_opcode_collector<BV> collector;
        bm::serializer<BV>::bv_ref_vector_type profile_refs;
        if (xor_refs) {
            for (unsigned p = 0; p < source.get_bmatrix().rows(); ++p)
                if (const BV* plane = source.get_bmatrix().get_row(p)) profile_refs.add(plane, p);
            collector.set_ref_vectors(&profile_refs);
        }
        collector.set_collect_ids();
        for (unsigned pi = 0; pi < (xor_refs ? 5u : 4u); ++pi)
        {
            const unsigned p = pi == 4 ? 31 : pi;
            const size_t offset = index.plane_offset(p);
            const auto* row = index.get_row(p);
            TRANSFORM_CHECK(offset && row && row->serialized_size());
            TRANSFORM_CHECK(offset + row->serialized_size() <= index.digest_offset());
            TRANSFORM_CHECK(layout.get_plane(p) == layout.data()+offset);
            bm::serialization_opcode_profile profile;
            unsigned residuals[256] = {}, pending = 0;
            collector.add_blob(layout.data()+offset, row->serialized_size(), profile,
                [&](const unsigned char* record, size_t) {
                    unsigned op = *record;
                    if (pending) ++residuals[op];
                    pending = (op >= bm::set_block_xor_ref8 && op <= bm::set_block_xor_chain) ||
                              (op >= bm::set_block_xor_ref8_um && op <= bm::set_block_xor_ref32_um);
                });
            TRANSFORM_CHECK(profile.reconciles());
            TRANSFORM_CHECK(profile.cardinalities_reconcile());
            TRANSFORM_CHECK(profile.ids == source.get_bmatrix().get_row(p)->count());
            std::cout << "  plane=" << p << " bookmarks=" << bookmarks << " opcodes:";
            for (unsigned op = 0; op < 256; ++op)
                if (profile.opcode[op].count) std::cout << " " << op << ":" << profile.opcode[op].count;
            std::cout << std::endl;
            if (p < 3) {
                unsigned code = p == 0 ? 72 : p == 1 ? 94 : 96;
                TRANSFORM_CHECK(profile.opcode[code].count >= 3);
                if (xor_refs) TRANSFORM_CHECK(residuals[code] >= 3);
                if (xor_refs && p == 0) TRANSFORM_CHECK(residuals[bm::set_block_1one] > 0);
                if (xor_refs && p == 1) TRANSFORM_CHECK(profile.opcode[bm::set_block_ref_eq].count > 0);
            }
            if (xor_refs) {
                // References disable group routing. Compare the same plane
                // with the full level-10 block policy and no references.
                bm::serializer<BV> ordinary;
                ordinary.set_compression_level(10);
                ordinary.set_sparse_cutoff(0);
                ordinary.set_sblock_large_options(false);
                ordinary.set_bookmarks(bookmarks, 4);
                bm::serializer<BV>::buffer plain;
                ordinary.serialize(*source.get_bmatrix().get_row(p), plain);
                TRANSFORM_CHECK(row->serialized_size() <= plain.size());
                for (unsigned code : {56u, 68u, 90u, 92u, 97u, 98u})
                    TRANSFORM_CHECK(profile.opcode[code].count == 0);
            }
            else if (p == 3)
                TRANSFORM_CHECK(profile.opcode[68].count + profile.opcode[90].count +
                                profile.opcode[92].count + profile.opcode[97].count == 2);
        }
        bm::sparse_vector_deserialization_index_serializer<BV> index_writer;
        bm::sparse_vector_deserialization_index_serializer<BV>::buffer index_bytes;
        TRANSFORM_CHECK(index_writer.serialize(index, index_bytes) > 0);
        bm::sparse_vector_deserialization_index_deserializer<BV> index_reader;
        TRANSFORM_CHECK(index_reader.deserialize(persisted, index_bytes.data(), index_bytes.size()) == index_bytes.size());
        TRANSFORM_CHECK(index.equal(persisted));
        TRANSFORM_CHECK(index.count_offsets() > 0);
        if (bookmarks) TRANSFORM_CHECK(index.count_bookmarks() > 0);

        SV full(bm::no_null);
        bm::sparse_vector_deserializer<SV> reader;
        reader.deserialize(full, layout.data());
        TRANSFORM_CHECK(source.equal(full));
        std::stringstream stream(std::ios::in | std::ios::out | std::ios::binary);
        bm::streams_encoder stream_writer(stream);
        TRANSFORM_CHECK(writer.serialize(source, stream_writer) && stream_writer.finish());
        const std::string blob = stream.str();
        TRANSFORM_CHECK(blob.size() == layout.size());
        TRANSFORM_CHECK(std::memcmp(blob.data(), layout.data(), blob.size()) == 0);
        stream.write("TRAILER", 7);
        stream.seekg(0);
        bm::streams_sparse_vector_deserializer<SV> stream_reader;
        TRANSFORM_CHECK(stream_reader.construct_deserialization_index(stream_index, stream));
        TRANSFORM_CHECK(stream_index.equal(index));
        TRANSFORM_CHECK(size_t(stream.tellg()) == blob.size() && stream.get() == 'T');
        stream.seekg(0);
        SV full_stream(bm::no_null);
        TRANSFORM_CHECK(stream_reader.deserialize(full_stream, stream));
        TRANSFORM_CHECK(source.equal(full_stream));
        TRANSFORM_CHECK(size_t(stream.tellg()) == blob.size() && stream.get() == 'T');

        for (unsigned selection = 0; selection < requested.size()+4; ++selection)
        {
            case_length = selection;
            BV mask;
            if (selection < requested.size()) mask.set(requested[selection]);
            else if (selection == requested.size())
                for (auto pos : requested) mask.set(pos);
            else if (selection == requested.size()+2)
                mask.set(SV::size_type(1)*bm::set_sub_total_bits+17);
            else if (selection == requested.size()+3) {
                mask.set(65536); mask.set(SV::size_type(2)*bm::set_sub_total_bits+65536);
            }
            SV expected(bm::no_null), actual(bm::no_null), streamed(bm::no_null);
            bm::sparse_vector_deserializer<SV> baseline;
            baseline.deserialize(expected, layout.data(), mask);
            for (auto e = mask.first(); e.valid(); ++e)
                if (expected.get(*e) != source.get(*e)) {
                    std::cerr << " gather pos=" << *e << " expected=" << expected.get(*e)
                              << " source=" << source.get(*e) << std::endl;
                    TRANSFORM_CHECK(false);
                }
            for (auto pos : requested)
                if (!mask.test(pos)) TRANSFORM_CHECK(expected.get(pos) == 0);
            reader.set_deserialization_index(&persisted);
            reader.set_deserialization_index_use();
            reader.deserialize(actual, layout.data(), mask);
            TRANSFORM_CHECK(expected.equal(actual));
            stream_reader.set_deserialization_index(&stream_index);
            stream_reader.set_deserialization_index_use();
            stream.seekg(0);
            TRANSFORM_CHECK(stream_reader.deserialize(streamed, stream, &mask));
            TRANSFORM_CHECK(expected.equal(streamed));
            TRANSFORM_CHECK(size_t(stream.tellg()) == blob.size() && stream.get() == 'T');
        }
    }
    std::cout << "TestSparseVectorOpcodeIndexGather OK" << std::endl;
}

#undef TRANSFORM_CHECK
} // namespace serial_transform_test
#endif
