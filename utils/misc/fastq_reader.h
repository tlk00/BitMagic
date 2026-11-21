#ifndef BITMAGIC_UTILS_FASTQ_READER_H
#define BITMAGIC_UTILS_FASTQ_READER_H

/*
Copyright(c) 2026 Anatoliy Kuznetsov(anatoliy_kuznetsov at yahoo.com)

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

For more information please visit:  http://bitmagic.io
*/

/// @file fastq_reader.h Streaming reader for unwrapped, four-line FASTQ records.
#include <cstdint>
#include <fstream>
#include <stdexcept>
#include <string>

namespace seqio
{

/// One FASTQ record. Sequence bases are normalized to uppercase A/C/G/T/N.
struct FastqRecord
{
    std::string defline;  ///< Header, including the leading '@'.
    std::string sequence; ///< Normalized sequence.
    std::string separator; ///< Separator, including the leading '+'.
    std::string quality;  ///< Original ASCII quality characters.
};

/// Read one FASTQ record at a time. The caller can reuse a record's buffers.
class FastqReader
{
public:
    /// Open an input stream owned by this reader. Throw if opening fails.
    /// @param filename Path to a plain four-line FASTQ file.
    explicit FastqReader(const std::string& filename);

    /// Read the next record and advance the record number. Throw on malformed input.
    /// @param record Caller-owned buffers overwritten during the read.
    /// @return True for a complete record, false at clean EOF.
    bool next(FastqRecord& record);

    /// Return the current one-based record number, or zero before the first read.
    std::uint64_t record_number() const { return record_number_; }

private:
    /// Read a line and remove a trailing carriage return from CRLF files.
    /// @param line Destination buffer.
    /// @return True on success, false at clean EOF.
    bool read_line(std::string& line);

    /// Throw a failure identifying the file and current record.
    /// @param reason Description of the input error.
    [[noreturn]] void fail(const char* reason) const;

    std::string filename_;        ///< Filename for diagnostics.
    std::ifstream stream_;        ///< Input stream.
    std::uint64_t record_number_ = 0; ///< Number of headers consumed.
};

inline
FastqReader::FastqReader(const std::string& filename)
    : filename_(filename), stream_(filename, std::ios::binary)
{
    if (!stream_)
    {
        throw std::runtime_error("Could not open FASTQ: " + filename_);
    }
}

inline
bool FastqReader::read_line(std::string& line)
{
    if (!std::getline(stream_, line))
    {
        if (stream_.bad() || !stream_.eof())
        {
            fail("input read failed");
        }
        return false;
    }
    if (!line.empty() && line.back() == '\r')
    {
        line.pop_back();
    }
    return true;
}

inline
void FastqReader::fail(const char* reason) const
{
    throw std::runtime_error(filename_ + ": record " + std::to_string(record_number_) + ": " + reason);
}

inline
bool FastqReader::next(FastqRecord& record)
{
    if (!read_line(record.defline))
    {
        return false;
    }
    ++record_number_;
    if (record.defline.empty() || record.defline.front() != '@')
    {
        fail("expected an '@' defline");
    }
    if (!read_line(record.sequence) || !read_line(record.separator) || !read_line(record.quality))
    {
        fail("incomplete four-line FASTQ record");
    }
    if (record.separator.empty() || record.separator.front() != '+')
    {
        fail("expected a '+' separator after one sequence line");
    }
    if (record.sequence.size() != record.quality.size())
    {
        fail("sequence and quality lengths differ");
    }

    for (char& base : record.sequence)
    {
        switch (base)
        {
        case 'a': base = 'A'; break;
        case 'c': base = 'C'; break;
        case 'g': base = 'G'; break;
        case 't': base = 'T'; break;
        case 'n': base = 'N'; break;
        case 'A': case 'C': case 'G': case 'T': case 'N': break;
        default: fail("sequence contains a symbol outside A/C/G/T/N");
        } // switch base
    } // for base
    return true;
}

} // namespace seqio

#endif // BITMAGIC_UTILS_FASTQ_READER_H
