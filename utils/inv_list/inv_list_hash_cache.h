#ifndef INV_LIST_HASH_CACHE_H_INCLUDED
#define INV_LIST_HASH_CACHE_H_INCLUDED
/* Copyright(c) 2026 Anatoliy Kuznetsov. Apache License, Version 2.0. */
/** @file inv_list_hash_cache.h
    @brief Bounded cache, SHA-256 and in-memory round-trip verification for inv_list.
    Application-only disk protocol and pipeline; canonical traversal is bmbvhash.h.
    Included after inv_list_parallel.h to reuse checked collection indices.
*/
#include <algorithm>
#include <array>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <limits>
#include <atomic>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <future>
#include <memory>
#include <string>
#include <system_error>
#include <vector>
#include "inv_list_sha256.h"
#include "bmbvhash.h"

namespace inv_hash
{
using Digest = std::array<unsigned char, 32>;
using BV = bm::bvector<>;
using Clock = std::chrono::steady_clock;

/** @internal SHA functor with a value-returning adapter for cache provenance. */
class Sha256 : public inv_list_sha256
{
public:
    /** Finalize once and return the 32 digest bytes. */
    Digest finish();
};

inline
Digest Sha256::finish()
{
    Digest result{};
    final(result.data());
    return result;
}

/** Encode value into the first width bytes of destination, little endian. */
static void put_le(unsigned char* destination, uint64_t value, unsigned width)
{
    for (unsigned i = 0; i < width; ++i)
    {
        destination[i] = static_cast<unsigned char>(value >> (i * 8));
    } // for i
}

/** Decode width little-endian bytes borrowed from source. */
static uint64_t get_le(const unsigned char* source, unsigned width)
{
    uint64_t result = 0;
    for (unsigned i = 0; i < width; ++i)
    {
        result |= uint64_t(source[i]) << (i * 8);
    } // for i
    return result;
}

/** Return checked file size of path; filesystem exceptions propagate. */
static uint64_t file_size(const std::string& path)
{
    return std::filesystem::file_size(path);
}

/** Read exactly size bytes into caller-owned data; throw with context on failure. */
static void read_exact(std::istream& stream, void* data, size_t size)
{
    stream.read(static_cast<char*>(data), std::streamsize(size));
    if (!stream)
    {
        throw std::runtime_error("Truncated input");
    }
}

/** Reuse inv_list's checked index reader; return owned offsets for the cache pipeline. */
static std::vector<uint64_t> offsets(const std::string& path, const std::string& index)
{
    const auto source = read_inv_blob_offsets(index, bm::id64_t(file_size(path)));
    return std::vector<uint64_t>(source.begin(), source.end());
}

/** Decode bv from borrowed bytes of size using worker-owned reader; demand exact consumption. */
static void decode(BV& bv, const unsigned char* bytes, size_t size,
                   bm::streams_deserializer<BV>& reader)
{
    bm::serialization_profile_streambuf view(bytes, size);
    std::istream input(&view);
    bm::streams_decoder channel(input);
    if (!reader.deserialize(bv, channel) || channel.tell() != size)
    {
        throw std::runtime_error("BLOB decode/length mismatch");
    }
}

/** @internal Cache pipeline options; all paths are caller supplied. */
struct Options
{
    std::string mode; ///< cache-bv, cache-raw, cached or roundtrip.
    std::string candidate; ///< Candidate .bv input.
    std::string candidate_index; ///< Candidate BLOB index.
    std::string raw; ///< LE length-prefixed raw uint32 vectors.
    std::string cache; ///< Digest cache input or fresh output.
    unsigned threads = 4; ///< Background worker count.
    uint64_t budget = 512ull * 1024 * 1024; ///< Combined staged input byte limit.
    unsigned batch_vectors = 512; ///< Maximum vectors per storage slot.
    uint64_t max_vectors = 0; ///< Prefix vector limit, zero means whole input.
    unsigned start_level = 1; ///< First in-memory round-trip compression level.
    unsigned end_level = 10; ///< Last round-trip level, inclusive.
    bool silent = false; ///< Suppress per-batch progress.
    bool mapped = false; ///< Unix mapped BLOB windows.
};

/** Hash exactly length bytes from file start with constant scratch storage. */
static Digest file_digest(const std::string& path, uint64_t length)
{
    std::ifstream input(path, std::ios::binary);
    Sha256 hash;
    std::vector<unsigned char> buffer(65536);
    while (length)
    {
        const size_t take = size_t(std::min<uint64_t>(length, buffer.size()));
        read_exact(input, buffer.data(), take);
        hash.update(buffer.data(), take);
        length -= take;
    } // while length
    return hash.finish();
}

/** @internal Per-vector input slices and deterministic output record, owned by a slot. */
struct Job
{
    uint64_t ordinal = 0; ///< Absolute vector ordinal.
    uint64_t begin = 0; ///< Candidate BLOB start or raw payload offset.
    uint64_t end = 0; ///< Candidate BLOB end or raw payload end.
    std::vector<unsigned char> raw; ///< Raw IDs bounded by the staging budget.
    std::array<unsigned char, 40> record{}; ///< SHA-256 plus LE uint64 cardinality.
};

/** @internal One owned batch; workers borrow it until all futures are joined. */
struct Batch
{
    std::unique_ptr<inv_list_read_buffer> candidate; ///< Candidate storage/mapping owner.
    std::vector<Job> jobs; ///< Input and output records.
    size_t staged = 0; ///< Charged bytes across both inputs.
    uint64_t begin = 0; ///< First candidate BLOB offset.
};

/** @internal Per-worker accumulated CPU wall durations and logical counts. */
struct Metrics
{
    double decode = 0; ///< Candidate decode/import and cardinality seconds.
    double check = 0; ///< Equality or hashing seconds.
    uint64_t ids = 0; ///< Candidate set bits.
    std::array<uint64_t, 11> vectors{}; ///< Successful round trips per level.
    std::array<uint64_t, 11> bytes{}; ///< Rebuilt encoded bytes per level.
    std::array<double, 11> encode_seconds{}; ///< Serialization worker seconds.
    std::array<double, 11> restore_seconds{}; ///< Rebuilt deserialization worker seconds.
    std::array<double, 11> compare_seconds{}; ///< Exact comparison worker seconds.
};

/** Return seconds since start using the monotonic clock. */
static double elapsed(Clock::time_point start)
{
    return std::chrono::duration<double>(Clock::now() - start).count();
}

/** Round-trip original through each selected serializer preset entirely in RAM.
    @tparam Serializer Encoder compatible with bm::serializer<BV>.
    @param original Worker-owned golden vector, optimized once without changing bits.
    @param options Level range, total buffer budget and worker count.
    @param ordinal Input vector index used by serializer policies and diagnostics.
    @param metric Worker-local successful level counts, sizes and durations.
    @param serializer Reusable worker-owned encoder, reset to each preset.
    @param reader Reusable worker-owned decoder for rebuilt vectors.
    @param buffer Worker-owned bounded encoded storage, grown without overlapping allocations.
    @param capacity Current buffer capacity, retained across vectors.
    @pre options.threads is positive and 1 <= start_level <= end_level <= 10.
    Throws on conservative buffer-bound overflow, encoding, decoding or equality failure.
*/
template<class Serializer>
static void roundtrip_vector(BV& original, const Options& options, uint64_t ordinal,
                             Metrics& metric, Serializer& serializer,
                             bm::streams_deserializer<BV>& reader,
                             std::unique_ptr<unsigned char[]>& buffer, size_t& capacity)
{
    {
        BV optimized(original);
        optimized.optimize();
        if (!original.equal(optimized))
        {
            throw std::runtime_error("Roundtrip optimization changed golden contents");
        }
        original.swap(optimized);
    }
    BV::statistics stats;
    original.calc_stat(&stats);
    const uint64_t buffer_limit = options.budget / 2 / options.threads;
    if (stats.max_serialize_mem > buffer_limit)
    {
        throw std::runtime_error("Roundtrip serialization bound exceeds worker budget; increase -read-memory-mib or reduce -threads");
    }
    if (stats.max_serialize_mem > capacity)
    {
        buffer.reset();
        capacity = 0;
        buffer.reset(new unsigned char[stats.max_serialize_mem]);
        capacity = stats.max_serialize_mem;
    }
    for (unsigned level = options.start_level; level <= options.end_level; ++level)
    {
        try
        {
            auto phase = Clock::now();
            auto settings = inv_serializer_options();
            settings.level = level;
            configure_inv_serializer(serializer, settings, plain_compression_policy::configure);
            serializer.set_periodic_vector_index(ordinal);
            serializer.set_bwt_vector_index(ordinal);
            const size_t bytes = size_t(serializer.serialize(original, buffer.get(), capacity));
            if (!bytes || bytes > capacity)
            {
                throw std::runtime_error("Invalid rebuilt BLOB size");
            }
            metric.encode_seconds[level] += elapsed(phase);
            phase = Clock::now();
            BV restored;
            decode(restored, buffer.get(), bytes, reader);
            metric.restore_seconds[level] += elapsed(phase);
            phase = Clock::now();
            if (!original.equal(restored))
            {
                throw std::runtime_error("Rebuilt vector comparison failed");
            }
            metric.compare_seconds[level] += elapsed(phase);
            ++metric.vectors[level];
            metric.bytes[level] += bytes;
        }
        catch (const std::exception& error)
        {
            throw std::runtime_error("roundtrip level=" + std::to_string(level) + " " + error.what());
        }
    } // for level
}

/** @internal Exclusive cache writer lock; held through publication or failure cleanup. */
class CacheLock
{
public:
    /** Acquire output's lock directory atomically or throw if it is held.
        @param output Fresh cache output path; parent directory must exist.
    */
    explicit CacheLock(const std::string& output);
    /** Release this writer's empty lock directory without throwing. */
    ~CacheLock();
    CacheLock(const CacheLock&) = delete;
    CacheLock& operator=(const CacheLock&) = delete;
private:
    std::filesystem::path path_; ///< Owned exclusive lock directory.
};

inline
CacheLock::CacheLock(const std::string& output) : path_(output + ".lock")
{
    if (!std::filesystem::create_directory(path_))
    {
        throw std::runtime_error("Hash cache writer lock exists: " + path_.string());
    }
}

inline
CacheLock::~CacheLock()
{
    std::error_code error;
    std::filesystem::remove(path_, error);
}

/** Execute one run. Futures are joined before their owning slot can be released.
    @param options Validated options, borrowed for the whole run.
    @return Zero after complete validation or cache publication; throws on any failure.
*/
static int run(const Options& options)
{
    const bool raw = options.mode == "cache-raw";
    const bool build = raw || options.mode == "cache-bv";
    const bool roundtrip = options.mode == "roundtrip";
    const bool verify_cache = options.mode == "cached" || (roundtrip && !options.cache.empty());
    const bool hash_needed = build || verify_cache;
    const uint64_t staging_budget = roundtrip ? options.budget / 2 : options.budget;
    if (hash_needed && !inv_list_sha256::available())
    {
        throw std::runtime_error("SHA-256 support disabled in this build");
    }
    const auto start = Clock::now();
    std::vector<uint64_t> candidate_offsets;
    if (!raw)
    {
        candidate_offsets = offsets(options.candidate, options.candidate_index);
    }
    const std::string source_path = raw ? options.raw : options.candidate;
    const uint64_t source_size = file_size(source_path);
    std::ifstream candidate_file(source_path, std::ios::binary);
    std::unique_ptr<CacheLock> writer_lock;
    std::ifstream cache_input;
    std::ofstream cache_output;
    std::string partial = options.cache + ".partial";
    uint64_t count = raw ? std::numeric_limits<uint64_t>::max() : candidate_offsets.size() - 1;
    if (options.max_vectors)
    {
        count = std::min(count, options.max_vectors);
    }
    if (build)
    {
        writer_lock.reset(new CacheLock(options.cache));
        if (std::filesystem::exists(options.cache) || std::filesystem::exists(partial))
        {
            throw std::runtime_error("Refusing to overwrite cache or partial: " + options.cache);
        }
        cache_output.open(partial, std::ios::binary);
        std::array<unsigned char, 64> header{};
        cache_output.write(reinterpret_cast<const char*>(header.data()), header.size());
        if (!cache_output)
        {
            throw std::runtime_error("Cannot create hash cache output: " + partial);
        }
    }
    if (verify_cache)
    {
        const uint64_t bytes = file_size(options.cache);
        cache_input.open(options.cache, std::ios::binary);
        std::array<unsigned char, 64> header{};
        read_exact(cache_input, header.data(), header.size());
        const uint64_t cached_count = get_le(header.data() + 8, 8);
        if (std::memcmp(header.data(), "BMSHA001", 8) || cached_count != count ||
            bytes < 96 || cached_count != (bytes - 96) / 40 || (bytes - 96) % 40)
        {
            throw std::runtime_error("Cache format/count/length mismatch");
        }
        Digest stored{};
        cache_input.seekg(std::streamoff(bytes - stored.size()));
        read_exact(cache_input, stored.data(), stored.size());
        if (stored != file_digest(options.cache, bytes - stored.size()))
        {
            throw std::runtime_error("Cache integrity mismatch");
        }
        cache_input.seekg(64);
    }
    std::unique_ptr<inv_list_read_buffer> candidate_loader;
    if (!raw)
    {
        candidate_loader.reset(new inv_list_read_buffer(options.candidate, options.mapped));
    }
    std::unique_ptr<Sha256> source_hash;
    if (build)
    {
        source_hash.reset(new Sha256);
    }
    std::unique_ptr<Batch> active;
    // Declaration order makes futures join before active storage unwinds.
    std::atomic<size_t> next_job{0};
    std::atomic<uint64_t> roundtrip_done{0}; ///< Vectors fully checked at all selected levels.
    auto next_roundtrip_report = Clock::now() + std::chrono::seconds(30);
    std::vector<std::future<Metrics>> workers;
    Metrics totals;
    uint64_t prepared = 0, completed = 0, raw_position = 0, input_bytes = 0;
    size_t peak_staged = 0;
    unsigned pending_length = 0;
    bool pending_raw = false;
    /** Drain every worker before propagating the first exception, then publish records. */
    auto drain = [&]()
    {
        std::exception_ptr error;
        for (auto& worker : workers)
        {
            try
            {
                if (roundtrip && !options.silent)
                {
                    while (worker.wait_for(std::chrono::seconds(1)) != std::future_status::ready)
                    {
                        if (Clock::now() >= next_roundtrip_report)
                        {
                            std::cerr << "Progress roundtrip completed=" << roundtrip_done.load()
                                      << " total=" << count << " elapsed_s=" << elapsed(start) << std::endl;
                            next_roundtrip_report = Clock::now() + std::chrono::seconds(30);
                        }
                    }
                }
                const auto metric = worker.get();
                totals.decode += metric.decode;
                totals.check += metric.check;
                totals.ids += metric.ids;
                for (unsigned level = options.start_level; roundtrip && level <= options.end_level; ++level)
                {
                    totals.vectors[level] += metric.vectors[level];
                    totals.bytes[level] += metric.bytes[level];
                    totals.encode_seconds[level] += metric.encode_seconds[level];
                    totals.restore_seconds[level] += metric.restore_seconds[level];
                    totals.compare_seconds[level] += metric.compare_seconds[level];
                } // for level
            }
            catch (...)
            {
                if (!error)
                {
                    error = std::current_exception();
                }
            }
        } // for worker
        workers.clear();
        if (error)
        {
            inv_list_check_worker(error, 1, std::string(roundtrip ? "Roundtrip worker failed: input=" : "Hash worker failed: input=") + source_path + " cache=" + options.cache);
        }
        if (active)
        {
            if (build)
            {
                for (const auto& job : active->jobs)
                {
                    cache_output.write(reinterpret_cast<const char*>(job.record.data()), std::streamsize(job.record.size()));
                } // for job
                if (!cache_output)
                {
                    throw std::runtime_error("Cache record write failed");
                }
            }
            completed += active->jobs.size();
            if (!options.silent)
            {
                if (roundtrip)
                {
                    std::cerr << "Progress roundtrip completed=" << completed
                              << " total=" << count << " elapsed_s=" << elapsed(start) << std::endl;
                }
                else
                {
                    std::cerr << "Progress verify vectors=" << completed << " elapsed_s=" << elapsed(start) << '\n';
                }
            }
            active.reset();
        }
    };
    while (prepared < count && (!raw || pending_raw || raw_position < source_size))
    {
        auto next = std::make_unique<Batch>();
        next->begin = raw ? 0 : candidate_offsets[size_t(prepared)];
        uint64_t limit = staging_budget / 2;
        while (prepared < count && next->jobs.size() < options.batch_vectors)
        {
            if (raw && !pending_raw)
            {
                if (raw_position == source_size)
                {
                    break;
                }
                unsigned char length[4];
                read_exact(candidate_file, length, 4);
                source_hash->update(length, 4);
                pending_length = unsigned(get_le(length, 4));
                raw_position += 4;
                pending_raw = true;
            }
            const uint64_t payload = raw ? uint64_t(pending_length) * 4 :
                candidate_offsets[size_t(prepared + 1)] - candidate_offsets[size_t(prepared)];
            const uint64_t charge = payload + 40;
            if (charge > staging_budget)
            {
                throw std::runtime_error("Vector exceeds staging budget at " + std::to_string(prepared));
            }
            if (next->jobs.empty())
            {
                if (charge > staging_budget - (active ? active->staged : 0))
                {
                    drain();
                }
                limit = std::min(std::max(limit, charge), staging_budget - (active ? active->staged : 0));
            }
            if (charge > limit - next->staged)
            {
                break;
            }
            Job job;
            job.ordinal = prepared;
            if (raw)
            {
                job.raw.resize(size_t(payload));
                read_exact(candidate_file, job.raw.data(), job.raw.size());
                source_hash->update(job.raw.data(), job.raw.size());
                raw_position += payload;
                pending_raw = false;
            }
            else
            {
                job.begin = candidate_offsets[size_t(prepared)];
                job.end = candidate_offsets[size_t(prepared + 1)];
            }
            if (verify_cache)
            {
                read_exact(cache_input, job.record.data(), job.record.size());
            }
            next->jobs.push_back(std::move(job));
            next->staged += size_t(charge);
            ++prepared;
        } // while prepared
        peak_staged = std::max(peak_staged, next->staged + (active ? active->staged : 0));
        if (!raw)
        {
            next->candidate.reset(new inv_list_read_buffer(options.candidate, options.mapped));
            const uint64_t bytes = next->jobs.back().end - next->begin;
            candidate_loader->load(candidate_file, next->begin, bytes);
            candidate_loader->swap_storage(*next->candidate);
            if (build)
            {
                source_hash->update(next->candidate->data(), size_t(bytes));
            }
            input_bytes += bytes;
        }
        drain();
        active = std::move(next);
        next_job = 0;
        for (unsigned worker = 0; worker < options.threads; ++worker)
        {
            workers.push_back(std::async(std::launch::async, [&]()
            {
                Metrics metric;
                std::unique_ptr<Sha256> vector_hash;
                std::unique_ptr<bm::bvector_hash<BV, Sha256>> traversal;
                if (hash_needed)
                {
                    vector_hash.reset(new Sha256);
                    traversal.reset(new bm::bvector_hash<BV, Sha256>(*vector_hash));
                }
                bm::streams_deserializer<BV> source_reader;
                std::unique_ptr<bm::streams_deserializer<BV>> restored_reader;
                std::unique_ptr<bm::serializer<BV>> serializer;
                if (roundtrip)
                {
                    serializer.reset(new bm::serializer<BV>);
                    restored_reader.reset(new bm::streams_deserializer<BV>);
                }
                std::unique_ptr<unsigned char[]> rebuilt;
                size_t rebuilt_capacity = 0;
                size_t index;
                while ((index = next_job.fetch_add(1)) < active->jobs.size())
                {
                    auto& job = active->jobs[index];
                    try
                    {
                        auto phase = Clock::now();
                        BV candidate;
                        if (raw)
                        {
                            std::vector<unsigned> ids(4096);
                            uint64_t previous = 0;
                            bool first = true;
                            for (size_t offset = 0; offset < job.raw.size(); )
                            {
                                const size_t take = std::min(ids.size(), (job.raw.size() - offset) / 4);
                                for (size_t i = 0; i < take; ++i)
                                {
                                    const unsigned id = unsigned(get_le(job.raw.data() + offset + i * 4, 4));
                                    if ((!first && id <= previous) || id >= bm::id_max)
                                    {
                                        throw std::runtime_error("Raw IDs must be strictly sorted, unique and in range");
                                    }
                                    ids[i] = id;
                                    previous = id;
                                    first = false;
                                } // for i
                                candidate.set(ids.data(), BV::size_type(take), bm::BM_SORTED);
                                offset += take * 4;
                            } // for offset
                        }
                        else
                        {
                            decode(candidate, active->candidate->data() + size_t(job.begin - active->begin),
                                   size_t(job.end - job.begin), source_reader);
                        }
                        const uint64_t cardinality = candidate.count();
                        metric.decode += elapsed(phase);
                        phase = Clock::now();
                        if (verify_cache && cardinality != get_le(job.record.data() + 32, 8))
                        {
                            throw std::runtime_error("Cached population count mismatch");
                        }
                        if (hash_needed)
                        {
                            Digest actual{};
                            traversal->compute(candidate, actual.data());
                            if (verify_cache && !std::equal(actual.begin(), actual.end(), job.record.begin()))
                            {
                                throw std::runtime_error("Cached SHA-256 mismatch");
                            }
                            std::copy(actual.begin(), actual.end(), job.record.begin());
                            put_le(job.record.data() + 32, cardinality, 8);
                        }
                        metric.check += elapsed(phase);
                        if (roundtrip)
                        {
                            roundtrip_vector(candidate, options, job.ordinal, metric,
                                             *serializer, *restored_reader, rebuilt, rebuilt_capacity);
                            ++roundtrip_done;
                        }
                        metric.ids += cardinality;
                    }
                    catch (const std::exception& error)
                    {
                        throw std::runtime_error("vector=" + std::to_string(job.ordinal) + " blob_begin=" + std::to_string(job.begin) + " blob_end=" + std::to_string(job.end) + " " + error.what());
                    }
                } // while index
                return metric;
            }));
        } // for worker
    } // while prepared
    drain();
    if (!completed)
    {
        throw std::runtime_error("Empty collection");
    }
    if (build)
    {
        std::array<unsigned char, 64> header{};
        std::memcpy(header.data(), "BMSHA001", 8);
        put_le(header.data() + 8, completed, 8);
        put_le(header.data() + 16, raw ? 1 : 2, 8);
        put_le(header.data() + 24, raw ? raw_position : input_bytes, 8);
        const auto source = source_hash->finish();
        std::copy(source.begin(), source.end(), header.begin() + 32);
        cache_output.seekp(0);
        cache_output.write(reinterpret_cast<const char*>(header.data()), header.size());
        cache_output.close();
        if (!cache_output)
        {
            throw std::runtime_error("Cache header/close failed");
        }
        const auto checksum = file_digest(partial, file_size(partial));
        std::ofstream footer(partial, std::ios::binary | std::ios::app);
        footer.write(reinterpret_cast<const char*>(checksum.data()), checksum.size());
        footer.close();
        if (!footer)
        {
            throw std::runtime_error("Cache footer failed");
        }
        std::filesystem::rename(partial, options.cache);
    }
    std::cout << "PASS mode=" << options.mode << " vectors=" << completed << " ids=" << totals.ids
              << " wall_s=" << elapsed(start) << " decode_worker_s=" << totals.decode
              << " check_worker_s=" << totals.check << " input_bytes=" << (raw ? raw_position : input_bytes)
              << " peak_staged_bytes=" << peak_staged << " budget_bytes=" << options.budget << '\n';
    if (roundtrip)
    {
        for (unsigned level = options.start_level; level <= options.end_level; ++level)
        {
            if (totals.vectors[level] != completed)
            {
                throw std::logic_error("Roundtrip level accounting mismatch");
            }
            std::cout << "Roundtrip level=" << level << " vectors=" << totals.vectors[level]
                      << " bytes=" << totals.bytes[level]
                      << " encode_worker_s=" << totals.encode_seconds[level]
                      << " decode_worker_s=" << totals.restore_seconds[level]
                      << " compare_worker_s=" << totals.compare_seconds[level] << '\n';
        } // for level
        std::cout << "Roundtrip complete." << '\n'
                  << "Total round trips=" << completed * (options.end_level - options.start_level + 1) << '\n'
                  << "Roundtrip buffer budget per worker=" << options.budget / 2 / options.threads << '\n';
    }
    std::cout << "Verification complete." << '\n'
              << "Total vectors=" << completed << '\n'
              << "Total ints=" << totals.ids << '\n';
    return 0;
}


} // namespace inv_hash
#endif
