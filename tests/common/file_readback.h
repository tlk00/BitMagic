// Checked fixture I/O for the 32- and 64-bit serialization stress tests.
// The including stress suite supplies its portable filesystem alias bmtest_fs.
#ifndef BM_TEST_FILE_READBACK_H__
#define BM_TEST_FILE_READBACK_H__

#include <cerrno>
#include <cstdlib>
#include <iostream>
#include <string>
#include <system_error>

namespace file_readback
{
/**
    Log a failed fixture operation and abort, leaving the fixture on disk.
    @param file Failed input stream, with its original state still set.
    @param path Fixture name used to query the on-disk length.
    @param offset Intended absolute byte offset.
    @param count Requested read length in bytes, or zero for a seek.
    @param io_errno Value of errno immediately after the failed operation.
*/
[[noreturn]] inline
void report_failure(const std::istream& file, const std::string& path,
                    std::streamoff offset, std::streamsize count, int io_errno)
{
    std::cerr << "CheckFileSerialization " << (count ? "read" : "seek")
              << " failed: path=" << path << " file_offset=" << offset;
    if (count)
    {
        std::cerr << " requested=" << count << " gcount=" << file.gcount();
    }
    std::cerr << " eof=" << file.eof() << " fail=" << file.fail()
              << " bad=" << file.bad() << " errno=" << io_errno;
    std::error_code ec;
    const auto size = bmtest_fs::file_size(path, ec);
    if (ec)
    {
        std::cerr << " file_size unavailable (" << ec.message() << ")";
    }
    else
    {
        std::cerr << " file_size=" << size;
    }
    std::cerr << std::endl;
    std::abort();
}

/** Position file at absolute byte offset or abort with diagnostics for path. */
inline
void seek(std::istream& file, const std::string& path, std::streamoff offset)
{
    errno = 0;
    if (!file.seekg(offset))
    {
        report_failure(file, path, offset, 0, errno);
    }
}

/**
    Read count bytes into dest or abort with diagnostics for path.
    The caller has positioned file at the given absolute byte offset.
*/
inline
void read(std::istream& file, const std::string& path, std::streamoff offset,
          char* dest, std::streamsize count)
{
    errno = 0;
    if (!file.read(dest, count))
    {
        report_failure(file, path, offset, count, errno);
    }
}
} // namespace file_readback

#endif
