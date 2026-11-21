/*
Copyright(c) 2025 Anatoliy Kuznetsov(anatoliy_kuznetsov at yahoo.com)

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

#ifndef BITMAGIC_UTILS_FASTA_READER_H
#define BITMAGIC_UTILS_FASTA_READER_H

/// @file fasta_reader.h Streaming wrapped FASTA input for examples and utilities.
/// Sequence lines are uppercased without restricting the alphabet.
#include <cctype>
#include <cstddef>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace seqio
{

/// Stream nonempty FASTA lines without allocating a chromosome-sized record.
class FastaReader
{
public:
    /// Open filename for text input. Throw on failure.
    /// @param filename Path to a plain FASTA file.
    /// @param strip_cr Remove a trailing carriage return after text-mode translation.
    explicit FastaReader(const std::string& filename, bool strip_cr = true);

    /// Advance to the next nonempty defline (including '>') or uppercase sequence line.
    /// @param line Caller-owned reusable buffer, valid until modified by the caller.
    /// @return True for a line, false at clean EOF. Throw on an I/O error.
    bool next_line(std::string& line);

private:
    std::string filename_; ///< Input path for diagnostics.
    std::ifstream input_;  ///< Owned text input stream.
    bool strip_cr_;        ///< Strip CR when reading CRLF files in Unix text mode.
};

inline
FastaReader::FastaReader(const std::string& filename, bool strip_cr)
    : filename_(filename), input_(filename), strip_cr_(strip_cr)
{
    if (!input_)
    {
        throw std::runtime_error("Could not open FASTA input: " + filename_);
    }
}

inline
bool FastaReader::next_line(std::string& line)
{
    while (std::getline(input_, line))
    {
        if (strip_cr_ && !line.empty() && line.back() == '\r')
        {
            line.pop_back();
        }
        if (line.empty())
        {
            continue;
        }
        if (line.front() != '>')
        {
            for (char& base : line)
            {
                base = static_cast<char>(std::toupper(static_cast<unsigned char>(base)));
            } // for base
        }
        return true;
    } // while getline
    if (input_.bad() || !input_.eof())
    {
        throw std::runtime_error("Could not read FASTA input: " + filename_);
    }
    return false;
}

/// Concatenate sequence lines, accepting headerless input and skipping empty records.
/// @param reader Input consumed from its current position to EOF.
/// @param sequence Replaced with uppercase bases.
/// @param starts Optional output replaced with offsets of nonempty records after the first.
inline
void read_fasta(FastaReader& reader, std::vector<char>& sequence,
                std::vector<std::size_t>* starts = nullptr)
{
    sequence.clear();
    if (starts)
    {
        starts->clear();
    }
    bool new_record = false;
    std::string line;
    while (reader.next_line(line))
    {
        if (line.front() == '>')
        {
            new_record = true;
            continue;
        }
        if (starts && new_record && !sequence.empty())
        {
            starts->push_back(sequence.size());
        }
        new_record = false;
        sequence.insert(sequence.end(), line.begin(), line.end());
    } // while reader.next_line
}

} // namespace seqio

#endif // BITMAGIC_UTILS_FASTA_READER_H
