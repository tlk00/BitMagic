#ifndef BM_TEST_SPARSE_SERIAL_LEVELS_H
#define BM_TEST_SPARSE_SERIAL_LEVELS_H

#include <type_traits>
#include <cstdint>
#include <limits>
#include <bmserial_profile.h>
#include "stress_workload.h"

// Additional serialization contracts shared by the new fixtures. The existing
// 32-/64-bit generators, reference comparisons and bookmark policies stay local.
namespace sparse_serial_test
{
template<class Check>
void ForLevels(bool pedantic, Check check)
{
    static size_t ordinal = 0;
    const unsigned selected = unsigned(ordinal++ % (bm::set_compression_max+1));
    for (unsigned level = 0; level <= bm::set_compression_max; ++level)
        if (pedantic || level == selected) check(level);
}

inline std::vector<unsigned> Levels(bool pedantic)
{
    std::vector<unsigned> levels;
    ForLevels(pedantic, [&](unsigned level) { levels.push_back(level); });
    return levels;
}

// Explicit float gather partitions are disjoint: 0..5, 6, and 7..10.
// Normal mode chooses one representative per requested partition so a shard
// never silently becomes an empty job because of rotating level selection.
inline
std::vector<unsigned> FloatGatherLevels(bool pedantic, unsigned partition)
{
    BM_ASSERT(partition <= 3);
    if (!partition)
        return Levels(pedantic);
    const unsigned first = partition == 1 ? 0u : partition == 2 ? 6u : 7u;
    const unsigned last = partition == 1 ? 5u : partition == 2 ? 6u : bm::set_compression_max;
    const unsigned selected = partition == 1 ? 0u : partition == 2 ? 6u : bm::set_compression_max;
    std::vector<unsigned> levels;
    for (unsigned level = first; level <= last; ++level)
    {
        if (pedantic || level == selected)
        {
            levels.push_back(level);
        }
    } // for level
    return levels;
}

// Bound repeated single-point gathers while retaining both ends of the source
// and evenly spaced interior samples. Block-island masks still visit every block.
template<class PointsVect>
void LimitGatherPoints(PointsVect& points, size_t limit)
{
    BM_ASSERT(limit >= 2);
    const size_t count = points.size();
    if (count <= limit)
        return;
    for (size_t i = 0; i < limit; ++i)
    {
        points[i] = points[i * (count - 1) / (limit - 1)];
    } // for i
    points.resize(limit);
}

template<class SV, bool Float = std::is_same<typename SV::value_type, float>::value>
struct Codecs
{
    using writer = bm::sparse_vector_serializer<SV>;
    using layout = bm::sparse_vector_serial_layout<SV>;
    using reader = bm::sparse_vector_deserializer<SV>;
    using stream_reader = bm::streams_sparse_vector_deserializer<SV>;
    using index_writer = bm::sparse_vector_deserialization_index_serializer<typename SV::bvector_type>;
    using index_reader = bm::sparse_vector_deserialization_index_deserializer<typename SV::bvector_type>;
};
template<class SV> struct Codecs<SV, true>
{
    using writer = bm::sparse_vector_float_serializer<SV>;
    using layout = bm::sparse_vector_float_serial_layout<SV>;
    using reader = bm::sparse_vector_float_deserializer<SV>;
    using stream_reader = bm::streams_sparse_vector_float_deserializer<SV>;
    using index_writer = bm::sparse_vector_float_deserialization_index_serializer<SV>;
    using index_reader = bm::sparse_vector_float_deserialization_index_deserializer<SV>;
};

inline uint32_t FloatBits(float value)
{
    uint32_t bits;
    static_assert(sizeof(bits) == sizeof(value), "32-bit float required");
    std::memcpy(&bits, &value, sizeof(bits));
    return bits;
}

/** Supplement equality checks with exact float bit probes.
    @tparam SV Float sparse-vector type; source and restored share the type.
    @param source Expected vector, borrowed for this call.
    @param restored Vector to validate, borrowed for this call.
    @param requested Optional borrowed gather mask; NULL checks all positions.
    @param fast Driver fast-mode flag. Large non-nullable vectors retain exact
        plane equality and sample redundant point probes. Small masks and NULL
        maps retain their full point enumeration.
    @pre Include after defining is_fast; fast may be true only in that mode.
    Does not modify vectors or random-generator state.
*/
template<class SV> void CheckFloatBits(const SV& source, const SV& restored,
    const typename SV::bvector_type* requested = nullptr, bool fast = false)
{
    if constexpr (std::is_same<typename SV::value_type, float>::value)
    {
        assert(source.size() == restored.size());
        const auto* assigned = source.get_null_bvector();
        const auto* other = restored.get_null_bvector();
        assert(bool(assigned) == bool(other));
        if (assigned)
        {
            assert(*assigned == *other);
            typename SV::bvector_type selected(*assigned, bm::finalization::READWRITE);
            if (requested) selected &= *requested;
            for (auto en = selected.first(); en.valid(); ++en)
                assert(FloatBits(source.get(*en)) == FloatBits(restored.get(*en)));
        }
        else if (fast && (requested ? requested->count() : source.size()) > 65536)
        {
            // equal() compares all sign/exponent/mantissa planes, including
            // zero gaps. Bound only the supplementary per-position probes.
            assert(source.equal(restored));
            for (typename SV::size_type i = 0; i < source.size();
                 i = NextStressTestValue(i, source.size() - 1, 32))
            {
                if (!requested || requested->test(i))
                {
                    assert(FloatBits(source.get(i)) == FloatBits(restored.get(i)));
                }
            } // for i
        }
        else if (requested)
        {
            for (auto en = requested->first(); en.valid(); ++en)
                if (*en < source.size())
                    assert(FloatBits(source.get(*en)) == FloatBits(restored.get(*en)));
        }
        else
            for (typename SV::size_type i = 0; i < source.size(); ++i)
                assert(FloatBits(source.get(i)) == FloatBits(restored.get(i)));
    }
}

template<class SV> void Sync(SV& value)
{
    // RSC-backed float children share a NULL map; their sizes are synchronized
    // from that map, as in the existing float stress helpers. Ordinary RSC
    // vectors retain their explicit logical size, including trailing NULLs.
    value.sync(true, std::is_same<typename SV::value_type, float>::value);
}

template<class SV> void Compare(const SV& source, const SV& restored,
    const typename SV::bvector_type* requested = nullptr)
{
    assert(source.size() == restored.size());
    assert(source.equal(restored));
    CheckFloatBits(source, restored, requested);
}

// Profile complete plane records using validated offsets, not a byte scan.
// No ID collection is needed: index-only traversal can also classify XOR records
// without borrowing reconstructed reference planes from the composite object.
template<class SV>
void ProfileSparse(const unsigned char* bytes, size_t size,
                   bm::serialization_opcode_profile& total)
{
    using BV = typename SV::bvector_type;
    typename bm::sparse_vector_deserializer<SV>::deserialization_index_type index;
    bm::sparse_vector_deserializer<SV> reader;
    reader.construct_deserialization_index(index, bytes);
    bm::serialization_opcode_collector<BV> collector;
    for (unsigned p = 0; p < index.rows(); ++p)
    {
        const auto* row = index.get_row(p);
        if (!row || !row->serialized_size()) continue;
        const size_t offset = index.plane_offset(p), length = row->serialized_size();
        assert(offset && offset <= size && length <= size-offset);
        collector.add_blob(bytes+offset, length, total);
    }
}

template<class SV>
bm::serialization_opcode_profile Profile(const unsigned char* bytes, size_t size)
{
    bm::serialization_opcode_profile total;
    if constexpr (std::is_same<typename SV::value_type, float>::value)
    {
        const size_t header = 3 + 3*sizeof(size_t);
        assert(size >= header && !std::memcmp(bytes, "bf0", 3));
        size_t lengths[3]; std::memcpy(lengths, bytes+3, sizeof(lengths));
        size_t offset = header;
        for (unsigned part = 0; part < 3; ++part)
        {
            assert(offset <= size && lengths[part] <= size-offset);
            if (!part)
            {
                bm::serialization_opcode_collector<typename SV::bvector_type> collector;
                collector.add_blob(bytes+offset, lengths[part], total);
            }
            else ProfileSparse<typename SV::sparse_vector_u>(bytes+offset, lengths[part], total);
            offset += lengths[part];
        }
        assert(offset == size);
    }
    else ProfileSparse<SV>(bytes, size, total);
    assert(total.reconciles());
    return total;
}

template<class SV>
bm::serialization_opcode_profile CheckOne(const SV& source, unsigned level,
    bool xor_refs, bool gather, bool profile = false)
{
    using C = Codecs<SV>;
    using BV = typename SV::bvector_type;
    std::cout << "Sparse serialization ";
    if constexpr (std::is_same<typename SV::value_type, float>::value)
        std::cout << "float ";
    else if constexpr (SV::is_remap_support::value)
        std::cout << "string ";
    else
        std::cout << (SV::is_rsc_support::value ? "RSC " : "integer ")
                  << (std::is_signed<typename SV::value_type>::value ? "signed " : "unsigned ")
                  << sizeof(typename SV::value_type)*8 << "-bit ";
    std::cout << "level=" << level << " XOR=" << xor_refs
              << " size=" << source.size() << std::endl;
    typename C::writer writer;
    writer.set_compression_level(level);
    assert(writer.get_compression_level() == level);
    writer.set_compression_level(bm::set_compression_max+1);
    assert(writer.get_compression_level() == level);
    writer.set_xor_ref(xor_refs);
    typename C::layout ram;
    writer.serialize(source, ram);
    if (profile) std::cout << "  RAM serialized" << std::endl;
    typename C::reader reader;
    SV restored(source.get_null_support());
    reader.deserialize(restored, ram.data()); Sync(restored);
    Compare(source, restored);
    if (profile) std::cout << "  RAM restored" << std::endl;
    std::stringstream io(std::ios::in | std::ios::out | std::ios::binary);
    io.write("prefix", 6);
    bm::streams_encoder channel(io);
    // Repeated output and reader reuse also exercise retained child state.
    for (unsigned repeat = 0; repeat < 2; ++repeat)
        assert(writer.serialize(source, channel));
    assert(channel.finish());
    io.write("TRAILER", 7);
    const std::string bytes = io.str();
    assert(bytes.size() == 6+2*ram.size()+7);
    for (unsigned repeat = 0; repeat < 2; ++repeat)
        assert(!std::memcmp(bytes.data()+6+repeat*ram.size(), ram.data(), ram.size()));
    typename C::stream_reader streamed;
    io.seekg(6);
    for (unsigned repeat = 0; repeat < 2; ++repeat)
    {
        assert(streamed.deserialize(restored, io)); Sync(restored);
        Compare(source, restored);
        assert(size_t(io.tellg()) == 6+(repeat+1)*ram.size());
    }
    assert(io.get() == 'T');
    if (profile) std::cout << "  stream round trips complete" << std::endl;
    if (gather)
    {
        typename C::reader::deserialization_index_type ram_index, stream_index, persisted;
        reader.construct_deserialization_index(ram_index, ram.data());
        io.clear(); io.seekg(6);
        assert(streamed.construct_deserialization_index(stream_index, io));
        assert(size_t(io.tellg()) == 6+ram.size());
        assert(ram_index.equal(stream_index));
        typename C::index_writer index_writer;
        typename C::index_writer::buffer index_bytes;
        index_writer.serialize(ram_index, index_bytes);
        typename C::index_reader index_reader;
        assert(index_reader.deserialize(persisted, index_bytes.data(), index_bytes.size()) == index_bytes.size());
        assert(ram_index.equal(persisted));
        for (unsigned pattern = 0; pattern < 3; ++pattern)
        {
            BV mask;
            if (pattern && source.size())
            {
                mask.set(0); mask.set(source.size()-1);
                for (auto i : {65535u, 65536u, 131071u, 131072u})
                    if (i < source.size()) mask.set(i);
                if (const BV* assigned = source.get_null_bvector())
                {
                    unsigned n = 0;
                    for (auto en = assigned->first(); en.valid() && n < 8; ++en, ++n)
                        mask.set(*en);
                }
                if (pattern == 2 && source.size() > 100)
                    mask.set_range(100, std::min(typename SV::size_type(300), source.size()-1));
            }
            SV expected(source.get_null_support());
            reader.deserialize(expected, ram.data(), mask); Sync(expected);
            for (auto* index : {&ram_index, &stream_index, &persisted})
            {
                typename C::reader indexed;
                indexed.set_deserialization_index(index);
                indexed.set_deserialization_index_use(true);
                indexed.deserialize(restored, ram.data(), mask); Sync(restored);
                Compare(expected, restored, &mask);
                streamed.set_deserialization_index(index);
                streamed.set_deserialization_index_use(true);
                io.clear(); io.seekg(6);
                assert(streamed.deserialize(restored, io, &mask)); Sync(restored);
                Compare(expected, restored, &mask);
                assert(size_t(io.tellg()) == 6+ram.size());
            }
        }
        if (source.size())
        {
            const auto from = source.size() > 65542 ? source.size()-65542 : 0;
            const auto to = std::min(typename SV::size_type(from+1024), source.size()-1);
            BV range; range.set_range(from, to);
            SV expected(source.get_null_support());
            reader.deserialize_range(expected, ram.data(), from, to); Sync(expected);
            streamed.set_deserialization_index(&persisted);
            io.clear(); io.seekg(6);
            assert(streamed.deserialize(restored, io, &range)); Sync(restored);
            Compare(expected, restored, &range);
            assert(size_t(io.tellg()) == 6+ram.size());
        }
    }
    if (profile) std::cout << "  indexed gathers complete" << std::endl;
    return profile ? Profile<SV>(ram.data(), ram.size()) : bm::serialization_opcode_profile();
}

template<class SV>
void CheckLevels(const SV& source, bool pedantic, bool gather = false, bool xor_refs = true)
{
    ForLevels(pedantic, [&](unsigned level) { CheckOne(source, level, xor_refs, gather); });
}


// Three copies of each transformed block keep each value plane above the
// sparse-superblock cutoff, while keeping the rank-space fixture manageable.
template<class BV> std::vector<unsigned> OpcodeValues()
{
    const unsigned count = 9*bm::gap_max_bits;
    std::vector<unsigned> values(count);
    for (unsigned kind = 0; kind < 3; ++kind)
    {
        BV plane;
        std::mt19937 rng(kind == 1 ? 123u : 135u);
        if (!kind)
        {
            for (unsigned i = 0; i < bm::gap_max_bits; ++i)
                if (i%17 < 8) plane.set(i);
        }
        else
        {
            for (unsigned pos = 0; pos < bm::gap_max_bits;)
            {
                const unsigned length = 100 + rng()%3000;
                if (rng() & 1u) plane.set_range(pos, std::min(65535u, pos+length));
                pos += length+1;
            }
            for (unsigned i = 0; i < (kind == 1 ? 50u : 650u); ++i)
                plane.flip(rng()%bm::gap_max_bits);
        }
        for (unsigned copy = 0; copy < 3; ++copy)
            for (auto en = plane.first(); en.valid(); ++en)
                values[(copy*3+kind)*bm::gap_max_bits+unsigned(*en)] |= 1u << kind;
    }
    return values;
}

template<class SV> void CheckOpcodeSource(SV& source)
{
    Sync(source); source.optimize();
    const auto profile = CheckOne(source, 10, false, true, true);
    std::cout << "  transformed planes: BWT=" << profile.opcode[72].count
              << " GAP-XOR=" << profile.opcode[94].count
              << " BWT-GAP-XOR=" << profile.opcode[96].count << std::endl;
    assert(profile.opcode[72].count >= 3);
    assert(profile.opcode[94].count >= 3);
    assert(profile.opcode[96].count >= 3);
    // Internal XOR can replace equal planes and transform related residuals.
    // Full and indexed gather checks remain identical, without pinning its choice.
    const auto xor_profile = CheckOne(source, 10, true, true, true);
    (void)xor_profile;
}

template<class BV> void RunOpcodeFixtures()
{
    using U = bm::sparse_vector<unsigned, BV>;
    using I = bm::sparse_vector<int, BV>;
    using R = bm::rsc_sparse_vector<unsigned, U>;
    using S = bm::str_sparse_vector<char, BV, 16>;
    using F = bm::sparse_vector_float<U>;
    using RF = bm::sparse_vector_float<R>;
    const auto values = OpcodeValues<BV>();
    const typename BV::size_type value_count = 9 * bm::gap_max_bits;
    assert(values.size() == value_count);
    const typename BV::size_type base = sizeof(typename BV::size_type) > 4
        ? typename BV::size_type(1ull << 32) : 0;
    std::cout << "Sparse typed opcode fixtures (signed/string/RSC/float/float-RSC)" << std::endl;
    {
        I source(bm::use_null);
        for (unsigned i = 0; i < values.size(); ++i) source.set(base+i, int(values[i]));
        source.set(base+value_count+17, std::numeric_limits<int>::min());
        CheckOpcodeSource(source);
    }
    {
        S source(bm::use_null);
        for (unsigned i = 0; i < values.size(); ++i)
        { char text[] = {char('@' | values[i]), 0}; source.set(i, text); }
        CheckOpcodeSource(source);
        source.remap();
        // Remapping changes plane distributions, so only require decoded results.
        CheckOne(source, 10, false, true, true);
        CheckOne(source, 10, true, true, true);
    }
    {
        U original(bm::use_null);
        for (unsigned i = 0; i < values.size(); ++i) original.set(base+2*i, values[i]);
        original.resize(base+2*value_count+101); // trailing NULLs
        R source; source.load_from(original); Sync(source);
        CheckOpcodeSource(source);
    }
    {
        F source(bm::use_null);
        for (unsigned i = 0; i < values.size(); ++i)
        {
            uint32_t bits = 0x3f800000u | values[i]; float value;
            std::memcpy(&value, &bits, sizeof(value)); source.set(base+i, value);
        }
        source.set(base+value_count+17, -3.5f);
        CheckOpcodeSource(source);
    }
    {
        RF source;
        for (unsigned i = 0; i < values.size(); ++i)
        {
            uint32_t bits = 0x3f800000u | values[i]; float value;
            std::memcpy(&value, &bits, sizeof(value)); source.set(base+2*i, value);
        }
        Sync(source);
        CheckOpcodeSource(source);
    }
}

template<class SV> void CheckSuperblockSource(SV& source)
{
    Sync(source); source.optimize();
    const auto plain = CheckOne(source, 8, false, true, true);
    std::cout << "  large superblock records=" << plain.opcode[98].count << std::endl;
    assert(plain.opcode[98].count);
    const auto refs = CheckOne(source, 8, true, true, true);
    // Float signs in these fixtures are empty. All populated data planes use
    // references, so no superblock opcode should be selected for those planes.
    for (unsigned code : {56u, 68u, 90u, 92u, 97u, 98u})
        assert(!refs.opcode[code].count);
}

template<class BV> void RunSuperblockFixtures()
{
    using U = bm::sparse_vector<unsigned, BV>;
    using I = bm::sparse_vector<int, BV>;
    using R = bm::rsc_sparse_vector<unsigned, U>;
    using S = bm::str_sparse_vector<char, BV, 16>;
    using F = bm::sparse_vector_float<U>;
    using RF = bm::sparse_vector_float<R>;
    const typename BV::size_type base = sizeof(typename BV::size_type) > 4
        ? typename BV::size_type(1ull << 32) : 0;
    std::cout << "Sparse typed superblock fixtures" << std::endl;
    {
        I source(bm::use_null);
        for (unsigned i = 0; i < 210001; ++i) source.set(base+i, i%3 == 1 ? 1 : 0);
        CheckSuperblockSource(source);
    }
    {
        S source(bm::use_null);
        for (unsigned i = 0; i < 210001; ++i) source.set(i, i%3 == 1 ? "A" : "@");
        CheckSuperblockSource(source);
    }
    {
        U original(bm::use_null);
        for (unsigned i = 0; i < 210001; ++i) original.set(base+2*i, i%3 == 1 ? 1 : 0);
        CheckSuperblockSource(original);
        R source; source.load_from(original);
        CheckSuperblockSource(source);
    }
    {
        F source(bm::use_null);
        for (unsigned i = 0; i < 210001; ++i)
        {
            uint32_t bits = 0x3f800000u | unsigned(i%3 == 1); float value;
            std::memcpy(&value, &bits, sizeof(value)); source.set(base+i, value);
        }
        CheckSuperblockSource(source);
    }
    {
        RF source;
        for (unsigned i = 0; i < 210001; ++i)
        {
            uint32_t bits = 0x3f800000u | unsigned(i%3 == 1); float value;
            std::memcpy(&value, &bits, sizeof(value)); source.set(base+2*i, value);
        }
        CheckSuperblockSource(source);
    }
    if constexpr (sizeof(typename BV::size_type) > 4)
    {
        // A separate bounded case reaches superblock index 65536, requiring
        // the wide index representation, without repeating that huge logical
        // prefix for every level/type fixture.
        I wide(bm::no_null);
        const auto high = typename BV::size_type(1) << 40;
        for (unsigned i = 0; i < 70000; ++i) wide.set(high+1+i*3, 1);
        const auto profile = CheckOne(wide, 8, false, true, true);
        assert(profile.opcode[98].count);
    }

}

template<class BV> void RunLevelFixtures(bool pedantic)
{
    std::cout << "Sparse typed level fixtures" << std::endl;
    using U = bm::sparse_vector<unsigned, BV>;
    using I = bm::sparse_vector<int, BV>;
    using R = bm::rsc_sparse_vector<unsigned, U>;
    using RI = bm::rsc_sparse_vector<int, I>;
    using S = bm::str_sparse_vector<char, BV, 32>;
    using F = bm::sparse_vector_float<U>;
    using RF = bm::sparse_vector_float<R>;
    const typename BV::size_type base = sizeof(typename BV::size_type) > 4
        ? typename BV::size_type(1ull << 32) : 0;
    auto check = [&](auto& source)
    {
        Sync(source);
        CheckLevels(source, pedantic, true, false);
        source.optimize();
        CheckLevels(source, pedantic, true, true);
    };
    U empty(bm::use_null); check(empty);
    U zero; zero.resize(131073); check(zero);
    R all_null; all_null.resize(131073); check(all_null);
    F empty_float(bm::use_null); check(empty_float);
    U u(bm::use_null); I si(bm::use_null);
    for (unsigned i = 0; i < 1000; ++i)
    {
        auto pos = base + (i*2654435761u)%200000u;
        u.set(pos, (i*2654435761u) ^ (i >> 3));
        si.set(pos, (i&1) ? -int(i)-1 : int(i));
    }
    u.set(base, 0); u.set(base+65536, std::numeric_limits<unsigned>::max());
    si.set(base, 0); si.set(base+65535, std::numeric_limits<int>::min());
    si.set(base+65536, std::numeric_limits<int>::max());
    u.resize(base+210000); si.resize(base+210000);
    check(u); check(si);
    bm::sparse_vector<unsigned long long, BV> u64(bm::use_null);
    bm::sparse_vector<long long, BV> i64(bm::use_null);
    for (unsigned i = 0; i < 64; ++i)
    {
        u64.set(base+65530+i, 1ull << i);
        i64.set(base+131060+i, i&1 ? std::numeric_limits<long long>::min()+i : std::numeric_limits<long long>::max()-i);
    }
    check(u64); check(i64);
    R r; r.load_from(u); check(r);
    RI ri; ri.load_from(si); check(ri);
    S strings(bm::use_null);
    strings.set(0, ""); strings.set(65535, "alpha");
    strings.set(65536, "alphabet"); strings.set(131072, "z");
    strings.resize(210000); check(strings);
    strings.remap(); check(strings);
    F floats(bm::use_null); RF rsc_floats;
    const uint32_t bits[] = {0, 0x80000000u, 1, 0x007fffffu, 0x00800000u,
        0x3f800000u, 0xbf800000u, 0x7f7fffffu, 0xff7fffffu,
        0x7f800000u, 0xff800000u, 0x7fc12345u, 0xffc54321u};
    for (unsigned i = 0; i < sizeof(bits)/sizeof(bits[0]); ++i)
    {
        float value; std::memcpy(&value, &bits[i], sizeof(value));
        const auto pos = base + (i ? 65530+i*65537 : 0);
        floats.set(pos, value); rsc_floats.set(pos, value);
        assert(FloatBits(floats.get(pos)) == bits[i]);
    }
    Sync(rsc_floats);
    for (unsigned i = 0; i < sizeof(bits)/sizeof(bits[0]); ++i)
    {
        const auto pos = base + (i ? 65530+i*65537 : 0);
        assert(FloatBits(rsc_floats.get(pos)) == bits[i]);
    }
    check(floats); check(rsc_floats);
    F no_null;
    for (unsigned i = 0; i < 100; ++i) no_null.push_back(float(int(i)-50)*0.125f);
    check(no_null);
    bm::sparse_vector_float_serializer<F> reused;
    for (unsigned level : {10u, 0u, 8u, 6u, 7u, 9u, 3u})
    for (const F* source : {&floats, &no_null, &floats})
    {
        reused.set_compression_level(level); reused.set_xor_ref(false);
        bm::sparse_vector_float_serializer<F> fresh;
        fresh.set_compression_level(level); fresh.set_xor_ref(false);
        bm::sparse_vector_float_serial_layout<F> a, b;
        reused.serialize(*source, a); fresh.serialize(*source, b);
        assert(a.size() == b.size() && !std::memcmp(a.data(), b.data(), a.size()));
    }
    std::cout << "Sparse typed level fixtures OK" << std::endl;
}

} // namespace sparse_serial_test
#endif
