/*
Copyright(c) 2026 Anatoliy Kuznetsov

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy at http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing, software distributed
under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
CONDITIONS OF ANY KIND, either express or implied.
*/

/** \example sample27.cpp
    Serialize bvectors directly to a seekable binary file using streams_encoder.
    Sample a second vector with random_subset and reuse one serializer to write
    the vectors to separate files in sequence. Delete both files on exit.
    Restore each vector through a buffered file reader, then gather selected
    positions from both the file and a RAM BLOB. Compare both gather results
    with the source AND the selection mask. Optional verification checks byte
    compatibility between file and RAM serialization.
    \sa bm::streams_encoder
    \sa bm::streams_decoder
    \sa bm::streams_deserializer
    \sa bm::deserialization_index
    \sa bm::serializer
    \sa bm::random_subset

    \par Related serialization examples
    - \ref sample4.cpp "bvsample04": RAM serialization, buffer ownership, and ordinary deserialization.
    - \ref sample14.cpp "bvsample14": set algebra and count operations on serialized RAM BLOBs.
    - \ref sample22.cpp "bvsample22": bookmarks and selective range deserialization.
    - \ref bvsample01_64.cpp "bvsample01_64": large-address mode basics.
*/
/*! \file sample27.cpp
    \brief Buffered file I/O and indexed gather deserialization from file and RAM.
*/

#include <chrono>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "bm.h"
#include "bmserial.h"
#include "bmfio.h"
#include "bmrandom.h"
#include "bmalgo.h"

using bvector_type = bm::bvector<>;
using serializer_type = bm::serializer<bvector_type>;
using file_reader_type = bm::streams_deserializer<bvector_type>;
using dindex_type = bm::deserialization_index<bvector_type>;
using dindex_serializer_type = bm::deserialization_index_serializer<bvector_type>;
using dindex_deserializer_type = bm::deserialization_index_deserializer<bvector_type>;
namespace fs = std::filesystem;

// Deterministic source data with ordinary blocks, full runs and sparse bits.
// Compile with -DBM64ADDR to put the same patterns beyond the 32-bit range.
static void fill_vector(bvector_type& first)
{
#ifdef BM64ADDR
    const bvector_type::size_type base = bvector_type::size_type(1) << 40;
#else
    const bvector_type::size_type base = 0;
#endif
    unsigned rnd = 1;
    for (unsigned block = 0; block < 48; ++block)
        for (unsigned j = 0; j < 512; ++j)
        {
            rnd = rnd * 1664525u + 1013904223u;
            first.set(base + block * 65536 + (rnd >> 16));
        }
    first.set_range(base + 65536 * 50, base + 65536 * 53 - 1);
    for (unsigned i = 0; i < 4096; ++i)
        first.set(base + bvector_type::size_type(i) * 8191);
    first.optimize();
}

// A private directory prevents overwriting or deleting an existing user file.
// Destruction removes demo files on both normal completion and exceptions.
class demo_files
{
public:
    explicit demo_files(const fs::path& parent)
    {
        const auto stamp = std::chrono::high_resolution_clock::now()
                               .time_since_epoch().count();
        for (unsigned attempt = 0; attempt < 100; ++attempt)
        {
            directory = parent / ("sample27-" + std::to_string(stamp) +
                                   "-" + std::to_string(attempt));
            if (fs::create_directory(directory))
                return;
        }
        throw std::runtime_error("Cannot create a private demo directory");
    }
    ~demo_files()
    {
        std::error_code error;
        fs::remove_all(directory, error);
        if (error)
            std::cerr << "Cannot remove " << directory << ": " << error.message() << '\n';
    }
    demo_files(const demo_files&) = delete;
    demo_files& operator=(const demo_files&) = delete;
    fs::path directory;
};

/** Restore one complete bvector using a bounded input buffer.
    @param[out] restored Mutable destination, cleared before deserialization.
    @param path Binary file containing one serialized bvector.
    @param reader Reusable reader retaining its working memory.
    @note I/O failures are reported by exception in this demo.
*/
static
void deserialize_file(bvector_type& restored, const fs::path& path,
                      file_reader_type& reader)
{
    restored.clear();
    std::ifstream file(path, std::ios::binary);
    bm::streams_decoder input(file);
    if (!reader.deserialize(restored, input))
        throw std::runtime_error("File deserialization failed");
}

/** Build reusable positional metadata without restoring the complete vector.
    @param path File whose exact BLOB will be used for subsequent gathers.
    @param reader Reusable buffered reader, initially without index attachments.
    @param index Output marker/bookmark offsets belonging to this BLOB only.
    @note This initial scan has an I/O and CPU cost; reuse the index across queries.
*/
static
void build_file_index(const fs::path& path, file_reader_type& reader,
                      dindex_type& index)
{
    std::ifstream file(path, std::ios::binary);
    bm::streams_decoder input(file);
    bvector_type scratch;
    reader.set_deserialization_index_construct(&index);
    const bool ok = reader.deserialize(scratch, input);
    reader.unset_deserialization_index();
    if (!ok)
        throw std::runtime_error("File index construction failed");
}

/** Serialize a deserialization index in RAM and save it to a file.
    @param path Destination file for the standalone index record.
    @param index Index constructed for one exact serialized bvector BLOB.
    @return Number of bytes written to the file.
    @note The index serialization API is memory based. This helper writes its
          small BM-managed buffer to disk in one operation.
*/
static
std::size_t save_dindex(const fs::path& path, const dindex_type& index)
{
    dindex_serializer_type serializer;
    dindex_serializer_type::buffer buffer;
    const std::size_t size = serializer.serialize(index, buffer);

    std::ofstream file(path, std::ios::binary | std::ios::trunc);
    if (!file)
        throw std::runtime_error("Cannot open deserialization-index output file");
    file.write(reinterpret_cast<const char*>(buffer.data()),
               std::streamsize(size));
    file.close();
    if (!file)
        throw std::runtime_error("Cannot write deserialization-index file");
    return size;
}

/** Restore a deserialization index from a standalone file record.
    @param[out] index Destination index, replaced by the restored record.
    @param path File previously written by save_dindex().
    @note The complete index is intentionally buffered in RAM because indexes
          are small and the current persistence API is memory based.
*/
static
void load_dindex(dindex_type& index, const fs::path& path)
{
    std::ifstream file(path, std::ios::binary | std::ios::ate);
    if (!file)
        throw std::runtime_error("Cannot open deserialization-index input file");
    const std::streampos end = file.tellg();
    if (end <= std::streampos(0))
        throw std::runtime_error("Invalid deserialization-index file size");
    const std::size_t size = std::size_t(std::streamoff(end));
    dindex_serializer_type::buffer buffer;
    buffer.resize(size, false);
    file.seekg(0);
    file.read(reinterpret_cast<char*>(buffer.data()), std::streamsize(size));
    if (!file)
        throw std::runtime_error("Cannot read deserialization-index file");

    dindex_deserializer_type deserializer;
    const std::size_t consumed = deserializer.deserialize(index, buffer.data(), size);
    if (consumed != size)
        throw std::runtime_error("Deserialization-index file has trailing data");
}

/** Gather exact mask positions from a file using its previously built index.
    @param[out] gathered Mutable destination, cleared first; must not alias mask.
    @param path Binary file matching index.
    @param mask Logical bit positions requested by the application.
    @param index Marker/bookmark offsets for this exact BLOB.
    @param reader Reusable reader whose index/digest attachments are unset on return.
    @note Produces source AND mask without materializing all unrequested blocks.
*/
static
void gather_file(bvector_type& gathered, const fs::path& path,
                 const bvector_type& mask,
                 const dindex_type& index, file_reader_type& reader)
{
    gathered.clear(true);
    bvector_type digest;
    
    // construct a technical bvector it is charting what is approximately is needed
    // (digest gives block size precision, block 64K elements)
    mask.build_block_digest(digest); // Translate bit positions into requested blocks.
    std::ifstream file(path, std::ios::binary);
    bm::streams_decoder input(file);
    reader.set_deserialization_index_use(&index);
    reader.set_block_digest_vector_use(&digest);
    const bool ok = reader.deserialize(gathered, input);
    reader.unset_block_digest_vector();
    reader.unset_deserialization_index();
    if (!ok)
        throw std::runtime_error("File gather deserialization failed");
    // Selection during decoding is at block granularity. Apply the exact mask
    // afterward to remove other bits from the selected blocks.
    gathered &= mask;
}

/** Gather the same positions from an already resident RAM serialization.
    @param[out] gathered Mutable destination, cleared first; must not alias mask.
    @param blob RAM BLOB with exactly the same bytes as the indexed file.
    @param mask Logical positions to retain.
    @param index Reusable offsets into this exact serialized representation.
    @note Produces source AND mask using the existing RAM deserializer.
*/
static
void gather_ram(bvector_type& gathered, const serializer_type::buffer& blob,
                const bvector_type& mask, const dindex_type& index)
{
    gathered.clear();
    bvector_type digest;

    mask.build_block_digest(digest); // digest provides approximate block sized deserialization request
    bm::deserializer<bvector_type, bm::decoder> reader;
    
    reader.set_deserialization_index_use(&index);
    reader.set_block_digest_vector_use(&digest);
    reader.deserialize(gathered, blob.data());
    
    gathered &= mask; // Same exact-mask step as the file path.
}

/** Create disjoint query intervals, including positions where source bits are zero.
    @param[out] mask Mutable destination, cleared before constructing the mask.
    @param source Nonempty vector whose first block locates the sample address range.
    @note Produces a selection mask spanning nearby and distant blocks in either
          address mode. The output must not alias source.
*/
static
void make_gather_mask(bvector_type& mask, const bvector_type& source)
{
    mask.clear();
    bvector_type::size_type first = 0, last = 0;
    if (!source.find_range(first, last))
        throw std::runtime_error("Expected a nonempty sample vector");
    const auto base = (first >> bm::set_block_shift) << bm::set_block_shift;
    mask.set_range(base + 2 * 65536 + 128, base + 3 * 65536 + 255);
    mask.set_range(base + 50 * 65536 + 64, base + 50 * 65536 + 4095);
    mask.set_range(base + 400 * 65536 + 32, base + 401 * 65536 + 128);
}

// This full-file allocation is ONLY for optional validation, not file output.
static
void verify_file(const fs::path& path, const bvector_type& original,
                 size_t length, serializer_type& serializer)
{
    std::ifstream file;
    file.exceptions(std::ios::badbit | std::ios::failbit);
    file.open(path, std::ios::binary | std::ios::ate);
    if (file.tellg() != std::streampos(std::streamoff(length)))
        throw std::runtime_error("Unexpected file length");
    file.seekg(0);
    std::vector<unsigned char> bytes(length);
    file.read(reinterpret_cast<char*>(bytes.data()), std::streamsize(bytes.size()));

    serializer_type::buffer ram;
    serializer.serialize(original, ram); // unchanged RAM API/settings
    if (ram.size() != length || std::memcmp(ram.data(), bytes.data(), length))
        throw std::runtime_error("File bytes differ from RAM serialization");

    bvector_type restored;
    bm::deserialize(restored, bytes.data());
    if (restored != original)
        throw std::runtime_error("Restored vector differs from the input");
}

int main(int argc, char* argv[])
{
    try
    {
        fs::path parent = ".";
        bool verify = false, have_path = false;
        for (int i = 1; i < argc; ++i)
        {
            const std::string arg = argv[i];
            if (arg == "--help" || arg == "-h")
            {
                std::cout << "Usage: " << argv[0] << " [--verify] [directory]\n"
                          << "Creates two temporary demo files; both are deleted before exit.\n";
                return 0;
            }
            if (arg == "--verify")
                verify = true;
            else if (!have_path && !arg.empty() && arg[0] != '-')
            { parent = arg; have_path = true; }
            else
                throw std::runtime_error("Invalid arguments; use --help");
        }

        bvector_type first, second;
        fill_vector(first);
        bm::random_subset<bvector_type> sampler;
        const auto sample_count = first.count() / 4;
        sampler.sample(second, first, sample_count);
        if (second.count() != sample_count || bm::count_sub(second, first))
            throw std::runtime_error("Invalid random subset");
        std::cout << "Source bits: " << first.count()
                  << "; random subset bits: " << second.count() << '\n';

        demo_files files(parent);
        const fs::path paths[2] = {files.directory / "first.bv",
                                   files.directory / "subset.bv"};
        const fs::path index_paths[2] = {files.directory / "first.didx",
                                        files.directory / "subset.didx"};
        const bvector_type* vectors[2] = {&first, &second};
        size_t lengths[2] = {0, 0};

        // Construct exactly one serializer, before opening either file.
        serializer_type serializer;
        // Keep the default compression level for forward compatibility.
        serializer.set_bookmarks(true, 16);
        for (unsigned i = 0; i < 2; ++i)
        {
            // Each file has its own channel. The serializer and its retained
            // allocator-typed working buffer are reused for the second file.
            // Binary, seekable, and NOT ios::app: headers/bookmarks need patches.
            std::ofstream file(paths[i], std::ios::binary | std::ios::trunc);
            if (!file)
                throw std::runtime_error("Cannot open output file");
            bm::streams_encoder output(file); // borrows this file stream

            // Success means accepted; safe flushes release working memory.
            if (!serializer.serialize(*vectors[i], output))
                throw std::runtime_error("Serialization: write failed");
            lengths[i] = output.size();
            // Complete this channel before closing its file and opening another.
            // finish() does not close the stream or promise disk durability.
            if (!output.finish())
                throw std::runtime_error("Output finalization failed");
            file.close();
            if (!file)
                throw std::runtime_error("Closing output file failed");
            std::cout << "Wrote " << paths[i] << ": " << lengths[i] << " bytes\n";
        }
        file_reader_type reader; // Reuse input working memory across both files.
        bvector_type mask;
        make_gather_mask(mask, first);
        for (unsigned i = 0; i < 2; ++i)
        {
            // The full restore is a separate demonstration, not a prerequisite
            // for file gathering or for constructing the deserialization index.
            {
                bvector_type restored;
                deserialize_file(restored, paths[i], reader);
                if (!restored.equal(*vectors[i]))
                    throw std::runtime_error("Stream restore differs from source");
                std::cout << paths[i].filename() << ": stream restore equal()=true\n";
            }
            dindex_type index;
            build_file_index(paths[i], reader, index);
            const std::size_t index_size = save_dindex(index_paths[i], index);
            dindex_type restored_index;
            load_dindex(restored_index, index_paths[i]);
            if (!index.equal(restored_index))
                throw std::runtime_error("Restored deserialization index differs");
            std::cout << "  Saved and restored " << index_paths[i].filename()
                      << ": " << index_size << " bytes; equal()=true\n";
            bvector_type from_file;
            gather_file(from_file, paths[i], mask, restored_index, reader);

            // RAM gathering is an independent demo: here a complete RAM BLOB
            // is deliberately allocated. File gathering above needs no such BLOB.
            serializer_type::buffer ram;
            serializer.serialize(*vectors[i], ram);
            bvector_type from_ram;
            gather_ram(from_ram, ram, mask, restored_index);
            bvector_type expected(*vectors[i], bm::finalization::READWRITE);
            expected &= mask;
            const bool file_equal = expected.equal(from_file);
            const bool ram_equal = expected.equal(from_ram);
            std::cout << "  Gather bits: " << expected.count()
                      << "; source AND mask: file equal()=" << std::boolalpha << file_equal
                      << ", RAM equal()=" << ram_equal << '\n';
            if (!file_equal || !ram_equal)
                throw std::runtime_error("Gather result differs from source AND mask");
        }
        if (verify)
        {
            for (unsigned i = 0; i < 2; ++i)
                verify_file(paths[i], *vectors[i], lengths[i], serializer);
            std::cout << "Verification OK: exact RAM bytes and both vectors restored\n";
        }
        // Check cleanup errors here; the guard also attempts cleanup on failure,
        // after all file streams have been closed during exception unwinding.
        fs::remove_all(files.directory);
        std::cout << "Deleted all demo files and their temporary directory\n";

    }
    catch (const std::exception& ex)
    {
        // streams_encoder also propagates exceptions if enabled on the stream.
        // The demo-file guard removes partial output; failed I/O is not retried.
        std::cerr << "sample27: " << ex.what() << '\n';
        return 1;
    }
    return 0;
}
