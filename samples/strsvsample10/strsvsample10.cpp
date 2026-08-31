/*
Copyright(c) 2002-2026 Anatoliy Kuznetsov(anatoliy_kuznetsov at yahoo.com)

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

/** \example strsvsample10.cpp
  Example of mmap-backed gather deserialization of bm::str_sparse_vector<>.

  \sa bm::str_sparse_vector
  \sa bm::sparse_vector_serializer
  \sa bm::sparse_vector_deserializer
  \sa bm::sparse_vector_deserialization_index
*/

/*! \file strsvsample10.cpp
    \brief Example: str_sparse_vector<> mmap-backed gather deserialization
*/

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <iomanip>
#include <cstdio>
#include <cstring>
#include <stdexcept>

#if defined(__unix__) || defined(__APPLE__)
# define BM_SAMPLE_HAS_MMAP 1
# include <fcntl.h>
# include <sys/mman.h>
# include <sys/stat.h>
# include <unistd.h>
#else
# define BM_SAMPLE_HAS_MMAP 0
#endif

#include "bm.h"
#include "bmstrsparsevec.h"
#include "bmsparsevec_serial.h"

#include "bmundef.h" /* clear the pre-proc defines from BM */

using namespace std;

typedef bm::bvector<> bvector_type;
typedef bm::str_sparse_vector<char, bvector_type, 16> str_sv_type;
typedef bm::sparse_vector_serializer<str_sv_type> serializer_type;
typedef bm::sparse_vector_deserializer<str_sv_type> deserializer_type;
typedef bm::sparse_vector_serial_layout<str_sv_type> layout_type;
typedef deserializer_type::deserialization_index_type deserialization_index_type;

#if BM_SAMPLE_HAS_MMAP

static const char* sample_blob_name = "strsvsample10_data.bm";
static const str_sv_type::size_type sample_vector_size = 4u * 1024u * 1024u;
static const unsigned gather_count = 10;
static const unsigned island_count = 3;
static const unsigned island_len = 192;


struct memory_profile
{
    size_t plain_strings_estimate = 0;
    size_t sparse_vector_memory = 0;
    size_t serialized_blob_size = 0;
    size_t deserialization_index_memory = 0;
    size_t average_gathered_memory = 0;
};


static
double bytes_to_mb(size_t size)
{
    return double(size) / (1024.0 * 1024.0);
}


static
double percent_of(size_t size, size_t total)
{
    return total ? (100.0 * double(size) / double(total)) : 0.0;
}


static
void print_mb(const char* label, size_t size)
{
    cout << "  " << label << " = "
         << fixed << setprecision(2) << bytes_to_mb(size) << " MB"
         << defaultfloat << endl;
}


static
void make_value(str_sv_type::size_type idx, string& value)
{
    // The value pattern intentionally repeats prefixes and categories. This is
    // representative for dictionary-like strings and lets str_sparse_vector<>
    // demonstrate remapping and bit-plane compression.
    static const char* groups[] = {
        "acct", "order", "event", "metric", "token", "region", "device", "user"
    };

    value = groups[idx & 7u];
    value += "-";
    value += to_string(unsigned((idx >> 8) & 0xFFFFu));
    value += "-";
    value += to_string(unsigned((idx * 2654435761u) >> 20));
}


static
void build_source_vector(str_sv_type& str_sv, memory_profile& profile)
{
    cout << "Stage 1: constructing source string sparse vector" << endl;

    auto bi = str_sv.get_back_inserter();
    string value;
    profile.plain_strings_estimate = sizeof(vector<string>) +
        size_t(sample_vector_size) * sizeof(string);
    for (str_sv_type::size_type i = 0; i < sample_vector_size; ++i)
    {
        make_value(i, value);
        profile.plain_strings_estimate += value.size() + 1;
        bi = value.c_str();
    }
    bi.flush();

    // Remapping analyzes character frequencies per string plane and replaces
    // common symbols with smaller codes. This reduces both in-memory succinct
    // vector size and the serialized BLOB size.
    str_sv.remap();

    // Optimization compresses bit-vector planes before serialization. A real
    // preparation job would normally do this once, then save the BLOB for later
    // serving from disk.
    BM_DECLARE_TEMP_BLOCK(tb)
    str_sv.optimize(tb);

    str_sv_type::statistics st;
    str_sv.calc_stat(&st);
    profile.sparse_vector_memory = st.memory_used;

    cout << "  source vector size = " << str_sv.size() << " elements" << endl;
    print_mb("plain std::vector<std::string> estimate", profile.plain_strings_estimate);
    print_mb("optimized str_sparse_vector<> memory", profile.sparse_vector_memory);
}


static
void save_serialized_blob(const char* fname, const unsigned char* buf, size_t size)
{
    ofstream fout(fname, ios::out | ios::binary | ios::trunc);
    if (!fout.good())
        throw runtime_error("Cannot create serialized BLOB file");

    fout.write(reinterpret_cast<const char*>(buf), streamsize(size));
    if (!fout.good())
        throw runtime_error("Cannot write serialized BLOB file");
}


static
memory_profile serialize_to_file(const char* fname)
{
    memory_profile profile;
    str_sv_type str_sv(bm::no_null);
    build_source_vector(str_sv, profile);

    cout << "Stage 1: serializing vector with bookmarks" << endl;

    serializer_type serializer;

    // Bookmarks are serialized skip points. They add small metadata overhead to
    // the BLOB, but later allow the deserializer and its index to avoid long
    // sequential scans when only a small fraction of elements is requested.
    serializer.set_bookmarks(true, 16);

    layout_type layout;
    serializer.serialize(str_sv, layout);
    profile.serialized_blob_size = layout.size();

    save_serialized_blob(fname, layout.buf(), layout.size());
    cout << "  serialized file = " << fname << endl;
    print_mb("serialized BLOB size", profile.serialized_blob_size);
    cout << "  serialized BLOB is "
         << fixed << setprecision(1)
         << percent_of(profile.serialized_blob_size,
                       profile.sparse_vector_memory)
         << "% of optimized str_sparse_vector<> memory"
         << defaultfloat << endl;

    // When this function exits, the source vector and temporary serialized
    // layout are released. Stage 2 below opens only the file, mirroring a
    // lower-RAM retrieval process that does not keep the full vector resident.
    return profile;
}


class mapped_file
{
public:
    mapped_file(const char* fname)
        : fd_(-1), addr_(0), size_(0)
    {
        fd_ = ::open(fname, O_RDONLY);
        if (fd_ < 0)
            throw runtime_error("Cannot open serialized BLOB file");

        struct stat st;
        if (::fstat(fd_, &st) != 0)
        {
            close_fd();
            throw runtime_error("Cannot stat serialized BLOB file");
        }
        if (st.st_size <= 0)
        {
            close_fd();
            throw runtime_error("Serialized BLOB file is empty");
        }
        size_ = size_t(st.st_size);

        addr_ = ::mmap(0, size_, PROT_READ, MAP_PRIVATE, fd_, 0);
        if (addr_ == MAP_FAILED)
        {
            addr_ = 0;
            close_fd();
            throw runtime_error("mmap() failed for serialized BLOB file");
        }

#ifdef MADV_RANDOM
        // Gather requests touch a few separated regions. The random-access hint
        // helps platforms that use madvise() to tune paging decisions.
        ::madvise(addr_, size_, MADV_RANDOM);
#endif
    }

    ~mapped_file()
    {
        if (addr_)
            ::munmap(addr_, size_);
        close_fd();
    }

    const unsigned char* data() const
        { return static_cast<const unsigned char*>(addr_); }

    size_t size() const { return size_; }

private:
    void close_fd()
    {
        if (fd_ >= 0)
        {
            ::close(fd_);
            fd_ = -1;
        }
    }

private:
    int     fd_;
    void*   addr_;
    size_t  size_;
};


static
void build_gather_mask(bvector_type& mask_bv, unsigned gather_idx)
{
    mask_bv.clear();

    // Gather deserialization uses a bit-vector mask: 1 means "materialize this
    // logical element from the serialized BLOB". This sample asks for three
    // small islands per request, which is a typical shape for sparse retrieval.
    for (unsigned i = 0; i < island_count; ++i)
    {
        str_sv_type::size_type base =
            str_sv_type::size_type((gather_idx + 1) * (i + 3)) * 7919u;
        base += str_sv_type::size_type(i) * (sample_vector_size / 4);
        base %= sample_vector_size - island_len - 1;
        mask_bv.set_range(base, base + island_len - 1);
    }
}


static
void print_gather_ranges(const bvector_type& mask_bv)
{
    bvector_type::enumerator en = mask_bv.first();
    while (en.valid())
    {
        bvector_type::size_type from = *en;
        bvector_type::size_type to = from;
        for (++en; en.valid() && *en == to + 1; ++en)
            to = *en;
        cout << " [" << from << ".." << to << "]";
    }
}


static
size_t build_deserialization_index(deserialization_index_type& dindex,
                                   const unsigned char* mapped_buf,
                                   size_t mapped_size)
{
    cout << "Stage 2: building deserialization index" << endl;

    // The deserialization index is a map of the serialized stream. It records
    // marker offsets and bookmark landing points so index-assisted gather
    // deserialization can jump over unrelated compressed blocks. This reduces
    // disk paging and CPU decompression work when requests are small.
    deserializer_type index_builder;
    index_builder.construct_deserialization_index(dindex, mapped_buf);
    dindex.optimize();

    size_t index_memory = dindex.memory_used();
    cout << "  index resident memory = "
         << fixed << setprecision(2) << bytes_to_mb(index_memory) << " MB"
         << " (" << setprecision(3) << percent_of(index_memory, mapped_size)
         << "% of mapped BLOB)" << defaultfloat << endl;
    return index_memory;
}


static
memory_profile gather_from_mapped_blob(const char* fname)
{
    memory_profile profile;

    cout << "Stage 2: opening serialized BLOB with mmap" << endl;

    mapped_file mapped(fname);
    print_mb("mapped BLOB size", mapped.size());

    deserialization_index_type dindex;
    profile.deserialization_index_memory =
        build_deserialization_index(dindex, mapped.data(), mapped.size());

    cout << "Stage 2: gather deserialization from mapped BLOB" << endl;

    size_t gathered_memory_total = 0;
    for (unsigned i = 0; i < gather_count; ++i)
    {
        bvector_type mask_bv;
        build_gather_mask(mask_bv, i);

        str_sv_type str_sv_out(bm::no_null);
        deserializer_type deserializer;
        deserializer.set_deserialization_index(&dindex);
        deserializer.set_deserialization_index_use(true);

        // Only masked positions are reconstructed into str_sv_out. The full
        // source vector remains absent from memory.
        deserializer.deserialize(str_sv_out, mapped.data(), mask_bv);

        str_sv_type::statistics st;
        str_sv_out.calc_stat(&st);
        gathered_memory_total += st.memory_used;

        cout << "  gather " << (i + 1)
             << ": requested " << mask_bv.count() << " elements in";
        print_gather_ranges(mask_bv);
        cout << endl;
    }

    cout << "Stage 2: resident memory used by selected results" << endl;
    profile.average_gathered_memory = gathered_memory_total / gather_count;
    print_mb("average gathered str_sparse_vector<> memory",
             profile.average_gathered_memory);
    return profile;
}

#endif


int main(void)
{
#if BM_SAMPLE_HAS_MMAP
    try
    {
        memory_profile profile = serialize_to_file(sample_blob_name);
        memory_profile retrieval_profile = gather_from_mapped_blob(sample_blob_name);
        profile.deserialization_index_memory =
            retrieval_profile.deserialization_index_memory;
        profile.average_gathered_memory =
            retrieval_profile.average_gathered_memory;

        cout << "Summary: storage and resident-memory progression" << endl;
        print_mb("plain std::vector<std::string> estimate", profile.plain_strings_estimate);
        print_mb("optimized str_sparse_vector<> memory", profile.sparse_vector_memory);
        print_mb("serialized BLOB on disk", profile.serialized_blob_size);
        print_mb("resident deserialization index", profile.deserialization_index_memory);
        print_mb("average gathered output vector", profile.average_gathered_memory);

        if (std::remove(sample_blob_name) == 0)
            cout << "Removed generated file: " << sample_blob_name << endl;
        else
            cout << "Generated file left on disk: " << sample_blob_name << endl;
    }
    catch (std::exception& ex)
    {
        cerr << ex.what() << endl;
        return 1;
    }
#else
    cout << "strsvsample10 requires POSIX mmap APIs; this platform is not supported by the sample." << endl;
#endif

    return 0;
}
