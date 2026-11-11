#ifndef INV_LIST_READ_BUFFER_H_INCLUDED
#define INV_LIST_READ_BUFFER_H_INCLUDED
/*
Copyright(c) 2026 Anatoliy Kuznetsov
Licensed under the Apache License, Version 2.0.
https://www.apache.org/licenses/LICENSE-2.0
*/
/** @file inv_list_read_buffer.h
    @brief Bounded sequential BLOB batches for indexed inv_list verification/decode.
    @note Application-internal. Input files must remain unchanged while reading.
    Unix mapping covers one batch at a time; hints are advisory, not RAM reservations.
*/
#include <fstream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include "bm.h"
#if defined(__unix__) || defined(__APPLE__)
#include <fcntl.h>
#include <sys/mman.h>
#include <unistd.h>
#endif

/** @brief Own a sequential input batch or a read-only Unix mapping window.
    @internal
    Borrowed data stays valid until clear(), load() or destruction. The caller
    must join every worker before releasing the batch. Instances cannot be copied.
*/
class inv_list_read_buffer
{
public:
    /** Open the optional mapping handle.
        @param path Stable source file, owned by the caller.
        @param mapped Request Unix windowed mapping instead of buffered reads.
        @throws std::runtime_error If the requested backend cannot open the file.
    */
    inv_list_read_buffer(const std::string& path, bool mapped);
    /** Release the current batch and optional file handle after workers join. */
    ~inv_list_read_buffer();
    inv_list_read_buffer(const inv_list_read_buffer&) = delete;
    inv_list_read_buffer& operator=(const inv_list_read_buffer&) = delete;
    /** Load the next consecutive byte range, replacing the previous batch.
        @param source Buffered input positioned at begin, used by buffered mode.
        @param begin Absolute file offset, equal to the previous range end.
        @param bytes Positive range size, representable as size_t and streamsize.
        @pre No worker retains a reference to the previous batch.
        @throws std::runtime_error On invalid sizes, mapping or short input.
    */
    void load(std::ifstream& source, bm::id64_t begin, bm::id64_t bytes);
    /** Release storage after workers join; preserve sequential position. */
    void clear();
    /** @return Borrowed batch start, or nullptr when no batch is loaded. */
    const unsigned char* data() const { return data_; }
private:
    bool mapped_; ///< Selected backend.
    bm::id64_t next_offset_ = 0; ///< Required start of the next batch.
    std::unique_ptr<unsigned char[]> buffer_; ///< Buffered batch owner.
    const unsigned char* data_ = nullptr; ///< Borrowed usable range start.
#if defined(__unix__) || defined(__APPLE__)
    int fd_ = -1; ///< Read-only mapping handle, retained across batches.
    void* mapping_ = MAP_FAILED; ///< Page-aligned mapped range owner.
    size_t mapped_bytes_ = 0; ///< Mapped length including page alignment prefix.
#endif
};

inline
inv_list_read_buffer::inv_list_read_buffer(const std::string& path, bool mapped)
    : mapped_(mapped)
{
#if defined(__unix__) || defined(__APPLE__)
    if (mapped_)
    {
        fd_ = ::open(path.c_str(), O_RDONLY);
        if (fd_ < 0)
        {
            throw std::runtime_error("Cannot open mapped BV input");
        }
    }
#else
    (void)path;
    if (mapped_)
    {
        throw std::runtime_error("Mapped BV input requires Unix");
    }
#endif
}

inline
inv_list_read_buffer::~inv_list_read_buffer()
{
    clear();
#if defined(__unix__) || defined(__APPLE__)
    if (fd_ >= 0)
    {
        ::close(fd_);
    }
#endif
}

inline
void inv_list_read_buffer::clear()
{
    data_ = nullptr;
    buffer_.reset();
#if defined(__unix__) || defined(__APPLE__)
    if (mapping_ != MAP_FAILED)
    {
        ::munmap(mapping_, mapped_bytes_);
        mapping_ = MAP_FAILED;
        mapped_bytes_ = 0;
    }
#endif
}

inline
void inv_list_read_buffer::load(std::ifstream& source, bm::id64_t begin,
                                bm::id64_t bytes)
{
    clear();
    if (!bytes || begin != next_offset_ ||
        bytes > (std::numeric_limits<size_t>::max)() ||
        bytes > bm::id64_t((std::numeric_limits<std::streamsize>::max)()))
    {
        throw std::runtime_error("Invalid sequential BV batch range");
    }
#if defined(__unix__) || defined(__APPLE__)
    if (mapped_)
    {
        const long page_size = ::sysconf(_SC_PAGESIZE);
        if (page_size <= 0)
        {
            throw std::runtime_error("Cannot determine mapping page size");
        }
        const bm::id64_t aligned = begin - begin % bm::id64_t(page_size);
        const size_t prefix = size_t(begin - aligned);
        if (bytes > (std::numeric_limits<size_t>::max)() - prefix ||
            aligned > bm::id64_t((std::numeric_limits<off_t>::max)()))
        {
            throw std::runtime_error("BV mapping range is too large");
        }
        mapped_bytes_ = size_t(bytes) + prefix;
        mapping_ = ::mmap(nullptr, mapped_bytes_, PROT_READ, MAP_PRIVATE, fd_, off_t(aligned));
        if (mapping_ == MAP_FAILED)
        {
            throw std::runtime_error("Cannot map BV batch");
        }
#ifdef MADV_SEQUENTIAL
        ::madvise(mapping_, mapped_bytes_, MADV_SEQUENTIAL);
#endif
#ifdef MADV_WILLNEED
        ::madvise(mapping_, mapped_bytes_, MADV_WILLNEED);
#endif
        data_ = static_cast<const unsigned char*>(mapping_) + prefix;
    }
    else
#endif
    {
        buffer_.reset(new unsigned char[size_t(bytes)]);
        source.read(reinterpret_cast<char*>(buffer_.get()), std::streamsize(bytes));
        if (!source)
        {
            throw std::runtime_error("Truncated sequential BV batch");
        }
        data_ = buffer_.get();
    }
    next_offset_ = begin + bytes;
}
#endif
