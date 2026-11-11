/*
Copyright(c) 2026 Anatoliy Kuznetsov
Licensed under the Apache License, Version 2.0.
https://www.apache.org/licenses/LICENSE-2.0
*/
/* Experimental canonical-block SHA-256 validation for inv_list collections.
   Uses OpenSSL EVP and application-internal block access. Never changes inputs.
   Build and usage are documented in this directory's README.md. */
#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <future>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <openssl/evp.h>
#include "bm.h"
#include "bmserial.h"
#include "bmfio.h"
#include "bmserial_profile.h"
#include "inv_list_read_buffer.h"

using Digest = std::array<unsigned char, 32>;
using BV = bm::bvector<>;
using Clock = std::chrono::steady_clock;

/** @internal Owning streaming SHA-256 context. Each worker owns its contexts. */
class Sha256
{
    std::unique_ptr<EVP_MD_CTX, decltype(&EVP_MD_CTX_free)> context_{EVP_MD_CTX_new(), EVP_MD_CTX_free}; ///< EVP owner.
public:
    /** Initialize a fresh SHA-256 stream or throw on allocation/provider failure. */
    Sha256();
    /** Append bytes at data, borrowed only for this call, of length size. */
    void update(const void* data, size_t size);
    /** Finish this context once and return its digest. */
    Digest finish();
};

Sha256::Sha256()
{
    if (!context_ || EVP_DigestInit_ex(context_.get(), EVP_sha256(), nullptr) != 1)
    {
        throw std::runtime_error("SHA-256 initialization failed");
    }
}

void Sha256::update(const void* data, size_t size)
{
    if (EVP_DigestUpdate(context_.get(), data, size) != 1)
    {
        throw std::runtime_error("SHA-256 update failed");
    }
}

Digest Sha256::finish()
{
    Digest result{};
    unsigned length = 0;
    if (EVP_DigestFinal_ex(context_.get(), result.data(), &length) != 1 || length != result.size())
    {
        throw std::runtime_error("SHA-256 finalization failed");
    }
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

/** Hash logical nonzero blocks of bv, independent of allocation and GAP choices.
    Canonical stream is ASCII BMBLOCKSHA256V1 followed by records of 8192
    little-endian block bytes and a uint64 LE absolute block number postfix.
    Fixed heap scratch is reused for GAP expansion and portable byte encoding.
    The caller owns bv; no enumerated-ID array or mutation is used.
*/
static Digest block_digest(const BV& bv)
{
    Sha256 hash;
    const char domain[] = "BMBLOCKSHA256V1";
    hash.update(domain, sizeof(domain) - 1);
    BV::allocator_type allocator;
    const auto release = [&](bm::word_t* block) { allocator.free_bit_block(block); };
    std::unique_ptr<bm::word_t, decltype(release)> scratch(allocator.alloc_bit_block(), release);
    std::vector<unsigned char> bytes(bm::set_block_size * sizeof(bm::word_t));
    const auto& manager = bv.get_blocks_manager();
    if (!manager.is_init())
    {
        return hash.finish();
    }
    const uint64_t limit = uint64_t(manager.top_block_size()) * bm::set_sub_array_size;
    for (uint64_t number = 0; number < limit; ++number)
    {
        int finished = 0;
        const bm::word_t* block = manager.get_block(BV::block_idx_type(number), &finished);
        if (!block)
        {
            continue;
        }
        if (block == FULL_BLOCK_FAKE_ADDR)
        {
            block = FULL_BLOCK_REAL_ADDR;
        }
        else if (BM_IS_GAP(block))
        {
            bm::gap_convert_to_bitset(scratch.get(), BMGAP_PTR(block));
            block = scratch.get();
        }
        if (bm::bit_is_all_zero(block))
        {
            continue;
        }
        const uint32_t endian = 1;
        if (*reinterpret_cast<const unsigned char*>(&endian) == 1)
        {
            hash.update(block, bytes.size());
        }
        else
        {
            for (unsigned word = 0; word < bm::set_block_size; ++word)
            {
                put_le(bytes.data() + word * 4, block[word], 4);
            } // for word
            hash.update(bytes.data(), bytes.size());
        }
        unsigned char postfix[8];
        put_le(postfix, number, 8);
        hash.update(postfix, sizeof(postfix));
    } // for number
    return hash.finish();
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

/** Parse an existing inv_list index and validate exact coverage of path's bytes. */
static std::vector<uint64_t> offsets(const std::string& path, const std::string& index)
{
    std::ifstream input(index);
    std::string magic, extra;
    uint64_t bytes = 0, count = 0;
    if (!(input >> magic >> bytes >> count) || magic != "BM_INV_BLOB_INDEX_V1" ||
        bytes != file_size(path) || !count || count > 100000000)
    {
        throw std::runtime_error("Invalid index header: " + index);
    }
    std::vector<uint64_t> result;
    result.reserve(size_t(count + 1));
    for (uint64_t i = 0; i <= count; ++i)
    {
        uint64_t position = 0;
        if (!(input >> position) || position > bytes || (!i && position) ||
            (i && position <= result.back()))
        {
            throw std::runtime_error("Invalid index offset: " + index);
        }
        result.push_back(position);
    } // for i
    if (result.back() != bytes || input >> extra)
    {
        throw std::runtime_error("Index coverage mismatch: " + index);
    }
    return result;
}

/** Decode a fresh vector from borrowed bytes and demand exact indexed consumption. */
static void decode(BV& bv, const unsigned char* bytes, size_t size)
{
    bm::serialization_profile_streambuf view(bytes, size);
    std::istream input(&view);
    bm::streams_decoder channel(input);
    bm::streams_deserializer<BV> reader;
    if (!reader.deserialize(bv, channel) || channel.tell() != size)
    {
        throw std::runtime_error("BLOB decode/length mismatch");
    }
}

/** @internal Command-line experiment options; all paths are caller supplied. */
struct Options
{
    std::string mode; ///< exact, pair-sha, cache-bv, cache-raw, or cached.
    std::string candidate; ///< Candidate .bv input.
    std::string candidate_index; ///< Candidate BLOB index.
    std::string golden; ///< Golden .bv input.
    std::string golden_index; ///< Golden BLOB index.
    std::string raw; ///< LE length-prefixed raw uint32 vectors.
    std::string cache; ///< Digest cache input or fresh output.
    unsigned threads = 4; ///< Background worker count.
    uint64_t budget = 512ull * 1024 * 1024; ///< Combined staged input byte limit.
    unsigned batch_vectors = 512; ///< Maximum vectors per storage slot.
    uint64_t max_vectors = 0; ///< Prefix experiment limit, zero means whole input.
    bool mapped = false; ///< Unix mapped BLOB windows.
};

/** Parse argc/argv or throw for missing, unknown or inconsistent arguments. */
static Options parse(int argc, char** argv)
{
    Options result;
    for (int i = 1; i < argc; ++i)
    {
        const std::string key = argv[i];
        if (key == "--mmap")
        {
            result.mapped = true;
            continue;
        }
        if (++i == argc)
        {
            throw std::runtime_error("Missing value for " + key);
        }
        const std::string value = argv[i];
        if (key == "--mode")
        {
            result.mode = value;
        }
        else if (key == "--candidate")
        {
            result.candidate = value;
        }
        else if (key == "--candidate-index")
        {
            result.candidate_index = value;
        }
        else if (key == "--golden")
        {
            result.golden = value;
        }
        else if (key == "--golden-index")
        {
            result.golden_index = value;
        }
        else if (key == "--raw")
        {
            result.raw = value;
        }
        else if (key == "--cache")
        {
            result.cache = value;
        }
        else
        {
            size_t consumed = 0;
            if (value.empty() || value.find_first_not_of("0123456789") != std::string::npos)
            {
                throw std::runtime_error("Invalid numeric value: " + value);
            }
            const uint64_t number = std::stoull(value, &consumed);
            if (key == "--threads" && number && number <= 256)
            {
                result.threads = unsigned(number);
            }
            else if (key == "--batch-vectors" && number && number <= 1000000)
            {
                result.batch_vectors = unsigned(number);
            }
            else if (key == "--max-vectors" && number)
            {
                result.max_vectors = number;
            }
            else if (key == "--memory-mib" && number && number <= uint64_t(std::numeric_limits<size_t>::max()) / (1024 * 1024))
            {
                result.budget = number * 1024 * 1024;
            }
            else
            {
                throw std::runtime_error("Invalid option/value: " + key);
            }
        }
    } // for i
    const bool pair = result.mode == "exact" || result.mode == "pair-sha";
    const bool raw = result.mode == "cache-raw";
    const bool build = raw || result.mode == "cache-bv";
    if ((!pair && !build && result.mode != "cached") ||
        (!raw && (result.candidate.empty() || result.candidate_index.empty())) ||
        (pair && (result.golden.empty() || result.golden_index.empty())) ||
        ((build || result.mode == "cached") && result.cache.empty()) ||
        (raw && result.raw.empty()))
    {
        throw std::runtime_error("Modes: exact/pair-sha --candidate BV --candidate-index IDX --golden BV --golden-index IDX; cache-bv --candidate BV --candidate-index IDX --cache NEW; cache-raw --raw RAW --cache NEW; cached --candidate BV --candidate-index IDX --cache EXISTING. Common: --threads N --memory-mib N --batch-vectors N --mmap --max-vectors N (prefix experiment limit)");
    }
    return result;
}

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
    uint64_t golden_begin = 0; ///< Golden BLOB offset.
    uint64_t golden_end = 0; ///< Golden BLOB end.
    std::vector<unsigned char> raw; ///< Raw IDs bounded by the staging budget.
    std::array<unsigned char, 40> record{}; ///< SHA-256 plus LE uint64 cardinality.
};

/** @internal One owned batch; workers borrow it until all futures are joined. */
struct Batch
{
    std::unique_ptr<inv_list_read_buffer> candidate; ///< Candidate storage/mapping owner.
    std::unique_ptr<inv_list_read_buffer> golden; ///< Golden storage/mapping owner.
    std::vector<Job> jobs; ///< Input and output records.
    size_t staged = 0; ///< Charged bytes across both inputs.
    uint64_t begin = 0; ///< First candidate BLOB offset.
    uint64_t golden_begin = 0; ///< First golden BLOB offset.
};

/** @internal Per-worker accumulated CPU wall durations and logical counts. */
struct Metrics
{
    double decode = 0; ///< Candidate plus golden decode/import seconds.
    double check = 0; ///< Equality or hashing seconds.
    uint64_t ids = 0; ///< Candidate set bits.
};

/** Return seconds since start using the monotonic clock. */
static double elapsed(Clock::time_point start)
{
    return std::chrono::duration<double>(Clock::now() - start).count();
}

/** Execute one run. Futures are joined before their owning slot can be released.
    @param options Validated options, borrowed for the whole run.
    @return Zero after complete validation or cache publication; throws on any failure.
*/
static int run(const Options& options)
{
    const bool raw = options.mode == "cache-raw";
    const bool build = raw || options.mode == "cache-bv";
    const bool pair = options.mode == "exact" || options.mode == "pair-sha";
    const auto start = Clock::now();
    std::vector<uint64_t> candidate_offsets, golden_offsets;
    if (!raw)
    {
        candidate_offsets = offsets(options.candidate, options.candidate_index);
    }
    if (pair)
    {
        golden_offsets = offsets(options.golden, options.golden_index);
    }
    if (pair && golden_offsets.size() != candidate_offsets.size())
    {
        throw std::runtime_error("Collection vector count mismatch");
    }
    const std::string source_path = raw ? options.raw : options.candidate;
    const uint64_t source_size = file_size(source_path);
    std::ifstream candidate_file(source_path, std::ios::binary);
    std::ifstream golden_file(options.golden, std::ios::binary);
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
        if (std::filesystem::exists(options.cache) || std::filesystem::exists(partial))
        {
            throw std::runtime_error("Refusing to overwrite cache or partial: " + options.cache);
        }
        cache_output.open(partial, std::ios::binary);
        std::array<unsigned char, 64> header{};
        cache_output.write(reinterpret_cast<const char*>(header.data()), header.size());
    }
    if (options.mode == "cached")
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
    std::unique_ptr<inv_list_read_buffer> candidate_loader, golden_loader;
    if (!raw)
    {
        candidate_loader.reset(new inv_list_read_buffer(options.candidate, options.mapped));
    }
    if (pair)
    {
        golden_loader.reset(new inv_list_read_buffer(options.golden, options.mapped));
    }
    Sha256 source_hash;
    std::unique_ptr<Batch> active;
    // Declaration order makes futures join before active storage unwinds.
    std::atomic<size_t> next_job{0};
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
                const auto metric = worker.get();
                totals.decode += metric.decode;
                totals.check += metric.check;
                totals.ids += metric.ids;
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
            std::rethrow_exception(error);
        }
        if (active)
        {
            if (build)
            {
                for (const auto& job : active->jobs)
                {
                    cache_output.write(reinterpret_cast<const char*>(job.record.data()), job.record.size());
                } // for job
                if (!cache_output)
                {
                    throw std::runtime_error("Cache record write failed");
                }
            }
            completed += active->jobs.size();
            std::cerr << "progress vectors=" << completed << " elapsed_s=" << elapsed(start) << '\n';
            active.reset();
        }
    };
    while (prepared < count && (!raw || pending_raw || raw_position < source_size))
    {
        auto next = std::make_unique<Batch>();
        next->begin = raw ? 0 : candidate_offsets[size_t(prepared)];
        next->golden_begin = pair ? golden_offsets[size_t(prepared)] : 0;
        uint64_t limit = options.budget / 2;
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
                source_hash.update(length, 4);
                pending_length = unsigned(get_le(length, 4));
                raw_position += 4;
                pending_raw = true;
            }
            const uint64_t payload = raw ? uint64_t(pending_length) * 4 :
                candidate_offsets[size_t(prepared + 1)] - candidate_offsets[size_t(prepared)];
            const uint64_t golden_bytes = pair ? golden_offsets[size_t(prepared + 1)] - golden_offsets[size_t(prepared)] : 0;
            const uint64_t charge = payload + golden_bytes + 40;
            if (charge > options.budget)
            {
                throw std::runtime_error("Vector exceeds staging budget at " + std::to_string(prepared));
            }
            if (next->jobs.empty())
            {
                if (charge > options.budget - (active ? active->staged : 0))
                {
                    drain();
                }
                limit = std::min(std::max(limit, charge), options.budget - (active ? active->staged : 0));
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
                source_hash.update(job.raw.data(), job.raw.size());
                raw_position += payload;
                pending_raw = false;
            }
            else
            {
                job.begin = candidate_offsets[size_t(prepared)];
                job.end = candidate_offsets[size_t(prepared + 1)];
            }
            if (pair)
            {
                job.golden_begin = golden_offsets[size_t(prepared)];
                job.golden_end = golden_offsets[size_t(prepared + 1)];
            }
            if (options.mode == "cached")
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
                source_hash.update(next->candidate->data(), size_t(bytes));
            }
            input_bytes += bytes;
        }
        if (pair)
        {
            next->golden.reset(new inv_list_read_buffer(options.golden, options.mapped));
            const uint64_t bytes = next->jobs.back().golden_end - next->golden_begin;
            golden_loader->load(golden_file, next->golden_begin, bytes);
            golden_loader->swap_storage(*next->golden);
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
                            decode(candidate, active->candidate->data() + size_t(job.begin - active->begin), size_t(job.end - job.begin));
                        }
                        BV golden;
                        if (pair)
                        {
                            decode(golden, active->golden->data() + size_t(job.golden_begin - active->golden_begin), size_t(job.golden_end - job.golden_begin));
                        }
                        const uint64_t cardinality = candidate.count();
                        metric.decode += elapsed(phase);
                        phase = Clock::now();
                        if (options.mode == "exact")
                        {
                            if (!candidate.equal(golden))
                            {
                                throw std::runtime_error("Exact comparison mismatch");
                            }
                        }
                        else
                        {
                            const Digest actual = block_digest(candidate);
                            if (options.mode == "pair-sha" && actual != block_digest(golden))
                            {
                                throw std::runtime_error("Pair SHA-256 mismatch");
                            }
                            if (options.mode == "cached" &&
                                (!std::equal(actual.begin(), actual.end(), job.record.begin()) || cardinality != get_le(job.record.data() + 32, 8))) throw std::runtime_error("Cached SHA-256/count mismatch");
                            std::copy(actual.begin(), actual.end(), job.record.begin());
                            put_le(job.record.data() + 32, cardinality, 8);
                        }
                        metric.check += elapsed(phase);
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
        const auto source = source_hash.finish();
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
    return 0;
}

/** CLI entry point. Report failures with nonzero exit status and preserve partials. */
int main(int argc, char** argv)
{
    try
    {
        return run(parse(argc, argv));
    }
    catch (const std::exception& error)
    {
        std::cerr << "FAIL " << error.what() << '\n';
        return 1;
    }
}
