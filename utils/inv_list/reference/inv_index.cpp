/** @file inv_index.cpp
    Development utility for compressed sidecar inspection and limited indexed reads.
    Build separately; does not modify BLOB collections or publish reference files.
*/
#include <algorithm>
#include <iostream>
#include "../inv_list_index.h"

/** Parse a decimal unsigned argument with no sign or trailing characters.
    @param text Null-terminated decimal argument.
    @return Parsed 64-bit value; throws on invalid input or overflow.
*/
static bm::id64_t number(const char* text)
{
    std::string value(text);
    size_t end = 0;
    if (value.empty() || value.find_first_not_of("0123456789") != std::string::npos)
    {
        throw std::runtime_error("Expected unsigned decimal argument");
    }
    const auto result = std::stoull(value, &end);
    if (end != value.size())
    {
        throw std::runtime_error("Invalid numeric argument");
    }
    return result;
}

/** Validate samples spread across a collection using exact indexed slices.
    @param path Borrowed BLOB filename, opened read-only.
    @param offsets Validated N+1 boundaries, borrowed for this call.
    @return Number of distinct BLOBs successfully decoded.
    @throws std::runtime_error on I/O, decoding, or consumed-size mismatch.
*/
static size_t check_reads(const std::string& path,
                          const std::vector<bm::id64_t>& offsets)
{
    std::ifstream input(path, std::ios::binary);
    const size_t count = offsets.size() - 1;
    std::vector<size_t> samples = {0, count - 1};
    for (size_t i = 1; i < 8; ++i)
    {
        samples.push_back((count - 1) * i / 8);
    } // for i
    if (count > 1)
    {
        samples.push_back(1);
    }
    std::sort(samples.begin(), samples.end());
    samples.erase(std::unique(samples.begin(), samples.end()), samples.end());
    for (size_t i : samples)
    {
        const bm::id64_t length = offsets[i + 1] - offsets[i];
        // Bound development-tool memory for malformed or unexpectedly large BLOBs.
        if (length > 512ull * 1024 * 1024)
        {
            throw std::runtime_error("Sample exceeds 512 MiB validation limit");
        }
        input.seekg(std::streamoff(offsets[i]));
        std::string blob(size_t(length), '\0');
        input.read(&blob[0], std::streamsize(length));
        if (!input)
        {
            throw std::runtime_error("Cannot read indexed sample");
        }
        std::istringstream stream(blob, std::ios::binary);
        bm::streams_decoder channel(stream);
        bm::streams_deserializer<bm::bvector<> > reader;
        bm::bvector<> bv;
        if (!reader.deserialize(bv, channel) ||
            stream.tellg() != std::streampos(std::streamoff(length)))
        {
            throw std::runtime_error("Indexed sample decode or boundary mismatch");
        }
    } // for i
    return samples.size();
}

/** Run write, dump, or check without modifying BLOB collections.
    @param argc Argument count.
    @param argv write OUT LEVEL (decimal boundaries on stdin); dump IN FILE_BYTES;
                check INDEX BLOB (samples indexed BLOBs).
    @return Zero on success, one on invalid input or failed validation.
*/
int main(int argc, char** argv)
{
    try
    {
        const std::string action = argc > 1 ? argv[1] : "";
        if (action == "write" && argc == 4)
        {
            const auto level = number(argv[3]);
            if (level > 10)
            {
                throw std::runtime_error("Compression level must be 0..10");
            }
            std::vector<bm::id64_t> offsets;
            std::string value;
            while (std::cin >> value)
            {
                offsets.push_back(number(value.c_str()));
            } // while value
            if (!std::cin.eof())
            {
                throw std::runtime_error("Cannot read boundary input");
            }
            write_inv_blob_offsets(argv[2], offsets, unsigned(level));
            if (read_inv_blob_offsets(argv[2], offsets.back()) != offsets)
            {
                throw std::runtime_error("Written offset mismatch");
            }
            std::cout << "Compared offsets=" << offsets.size() << '\n';
        }
        else if (action == "dump" && argc == 4)
        {
            const auto bytes = number(argv[3]);
            const auto offsets = read_inv_blob_offsets(argv[2], bytes);
            for (auto offset : offsets)
            {
                std::cout << offset << '\n';
            } // for offset
        }
        else if (action == "check" && argc == 4)
        {
            std::ifstream blob(argv[3], std::ios::binary | std::ios::ate);
            const auto end = blob.tellg();
            if (end < 0)
            {
                throw std::runtime_error("Cannot size BLOB file");
            }
            const auto offsets = read_inv_blob_offsets(argv[2], bm::id64_t(end));
            const auto samples = check_reads(argv[3], offsets);
            std::cout << "Compared offsets=" << offsets.size()
                      << " indexed_reads=" << samples << '\n';
        }
        else
        {
            throw std::runtime_error("Usage: inv_index write OUT LEVEL (boundaries on stdin) | dump IN FILE_BYTES | check INDEX BLOB");
        }
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
