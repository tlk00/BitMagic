#ifndef INV_LIST_INDEX_H_INCLUDED
#define INV_LIST_INDEX_H_INCLUDED

/** @file inv_list_index.h
    Compressed BLOB sidecar I/O for inv_list and its inspection utility.
    Stores V2 BM sparse-vector lengths and rejects other formats. Uses
    64-bit values even in a 32-bit BM build; vector counts remain bounded.
    This sidecar is independent of BM's internal deserialization index.
*/
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include "bm.h"
#include "bmserial.h"
#include "bmsparsevec.h"
#include "bmsparsevec_serial.h"
#include "bmfio.h"

/** @internal Sparse vector of unsigned 64-bit BLOB lengths. */
using inv_blob_lengths = bm::sparse_vector<bm::id64_t, bm::bvector<> >;

/** Read and validate every boundary, returning owned offsets including zero.
    @param path Compressed V2 sidecar to read.
    @param file_bytes Exact associated BLOB file size, checked against the header.
    @return N+1 strictly increasing boundaries for N nonempty BLOBs.
    @throws std::runtime_error on invalid, truncated, stale or trailing data.
*/
inline
std::vector<bm::id64_t> read_inv_blob_offsets(
    const std::string& path, bm::id64_t file_bytes)
{
    std::ifstream input(path, std::ios::binary);
    std::string header, magic, extra;
    std::getline(input, header);
    std::istringstream fields(header);
    bm::id64_t indexed_bytes = 0, count = 0;
    unsigned level = 0;
    if (!(fields >> magic >> indexed_bytes >> count >> level) ||
        magic != "BM_INV_BLOB_INDEX_V2" || indexed_bytes != file_bytes ||
        !count || count > 100000000 || level > 10 || (fields >> extra))
    {
        throw std::runtime_error("Invalid or stale compressed BLOB index header");
    }
    inv_blob_lengths lengths;
    bm::streams_sparse_vector_deserializer<inv_blob_lengths> reader;
    if (!reader.deserialize(lengths, input) || lengths.size() != count ||
        input.peek() != std::char_traits<char>::eof())
    {
        throw std::runtime_error("Invalid compressed BLOB index payload");
    }
    std::vector<bm::id64_t> offsets;
    offsets.reserve(size_t(count + 1));
    offsets.push_back(0);
    for (bm::id64_t i = 0; i < count; ++i)
    {
        const bm::id64_t length = lengths.get(inv_blob_lengths::size_type(i));
        if (!length || length > file_bytes - offsets.back())
        {
            throw std::runtime_error("Invalid BLOB index length");
        }
        offsets.push_back(offsets.back() + length);
    } // for i
    if (offsets.back() != file_bytes)
    {
        throw std::runtime_error("BLOB index does not cover exact file");
    }
    return offsets;
}

/** Write a V2 sidecar with collection-level BM compression and 64-bit lengths.
    @param path Output file, replaced by this call; caller handles publication.
    @param offsets Validated N+1 boundaries including zero, borrowed for this call.
    @param level Collection compression level, in [0,10].
    @throws std::runtime_error on invalid boundaries or write failure.
*/
inline
void write_inv_blob_offsets(const std::string& path,
                           const std::vector<bm::id64_t>& offsets, unsigned level)
{
    if (offsets.size() < 2 || offsets.size() > 100000001 || offsets.front() ||
        level > 10)
    {
        throw std::runtime_error("Invalid BLOB index write arguments");
    }
    inv_blob_lengths lengths;
    auto inserter = lengths.get_back_inserter();
    for (size_t i = 1; i < offsets.size(); ++i)
    {
        if (offsets[i] <= offsets[i - 1])
        {
            throw std::runtime_error("Invalid BLOB index write offset");
        }
        inserter = offsets[i] - offsets[i - 1];
    } // for i
    inserter.flush();
    lengths.optimize();
    bm::sparse_vector_serializer<inv_blob_lengths> writer;
    writer.set_compression_level(level);
    std::ofstream output(path, std::ios::binary | std::ios::trunc);
    output << "BM_INV_BLOB_INDEX_V2 " << offsets.back() << ' '
           << offsets.size() - 1 << ' ' << level << '\n';
    bm::streams_encoder channel(output);
    if (!writer.serialize(lengths, channel) || !channel.finish())
    {
        throw std::runtime_error("Cannot serialize BLOB index");
    }
    output.close();
    if (!output)
    {
        throw std::runtime_error("Error writing BLOB index");
    }
}

#endif
