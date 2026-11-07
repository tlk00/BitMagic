#ifndef BM_SERIAL_PROFILE_H_INCLUDED
#define BM_SERIAL_PROFILE_H_INCLUDED

/** \file bmserial_profile.h
    \brief Internal diagnostics for opcodes in completed serialization streams.
    \ingroup bvserial
    \internal

    This optional, development utility inspects final bvector
    BLOBs. It does not participate in compressor trials or winner selection,
    and ordinary serialization does not invoke it. Tests use it to verify
    opcode coverage; inv_list uses it for opcode reports and surveys. Sparse
    vector tests profile individual serialized bit-plane BLOBs with this API.

    A deserialization-index pass identifies and validates top-level record
    spans without scanning arbitrary payload bytes for opcode values. The
    default mode collects counts and byte sizes. Optional original-ID metrics
    perform an additional full reconstruction pass and count restored set bits
    per record, including transformed records. Reference-XOR reconstruction
    requires the original reference vectors supplied through set_ref_vectors().

    Reuse the collector to retain its decoder, index and scratch allocations.
    The input view borrows the BLOB without copying it. Human-readable labels
    and report formatting are provided separately in bmdbg.h.

    Example for an existing serialized BLOB and application bvector type BV:
    \code
    bm::serialization_opcode_collector<BV> collector;
    bm::serialization_opcode_profile profile;
    collector.set_collect_ids(); // Optional: reconstruct original set bits.
    collector.add_blob(blob.data(), blob.size(), profile);
    assert(profile.reconciles() && profile.cardinalities_reconcile());
    \endcode
*/

#include "bmfio.h"
#include <cstdint>
#include <istream>
#include <limits>
#include <stdexcept>
#include <streambuf>

namespace bm
{

/** \brief Internal aggregate measurements of emitted serialization records.
    \ingroup bvserial
    \internal
    \note ID-list BLOBs have no opcode record markers and are accounted for
        separately. Cardinalities describe original set bits, not residuals.
*/
struct serialization_opcode_profile
{
    /** \brief Internal measurements for one opcode value.
        \ingroup bvserial
        \internal
    */
    struct entry
    {
        std::uint64_t count = 0; ///< Number of emitted records with this opcode.
        std::uint64_t bytes = 0; ///< Complete record bytes, including opcode.
        std::uint64_t ids = 0;   ///< Original set bits covered by these records.
    };
    entry opcode[256]; ///< Measurements indexed by opcode byte.
    std::uint64_t blobs = 0; ///< Number of profiled BLOBs.
    std::uint64_t bytes = 0; ///< Total profiled BLOB bytes.
    std::uint64_t header_bytes = 0; ///< Headers preceding opcode records.
    std::uint64_t id_list_bytes = 0; ///< Complete marker-free ID-list BLOB bytes.
    std::uint64_t records = 0; ///< Total top-level opcode records.
    std::uint64_t ids = 0; ///< Original set bits from cardinality-enabled BLOBs.
    std::uint64_t id_list_ids = 0; ///< Set bits from marker-free ID-list BLOBs.
    std::uint64_t cardinality_blobs = 0; ///< BLOBs profiled with ID collection enabled.

    /** Accumulate another profile without changing it.
        \param r Profile whose counters are added to this aggregate.
    */
    void add(const bm::serialization_opcode_profile& r) noexcept;
    /** Verify that all BLOBs have ID metrics and their record totals agree.
        \return True when every BLOB has cardinalities and they sum to ids.
    */
    bool cardinalities_reconcile() const noexcept;
    /** Sum complete opcode-record sizes, excluding headers and ID-list BLOBs.
        \return Total bytes attributed to opcode records.
    */
    std::uint64_t record_bytes() const noexcept;
    /** Check record counts, minimum record sizes and complete byte accounting.
        \return True when record counters and byte totals reconcile.
    */
    bool reconciles() const noexcept;
};

/** \brief Internal checked, read-only stream view of one borrowed BLOB.
    \ingroup bvserial
    \internal
    \note The backing bytes must remain valid and unchanged while the view is
        in use. The get area borrows their address; no input-byte copy is made.
*/
class serialization_profile_streambuf : public std::streambuf
{
public:
    /** Construct an input view over a bounded memory range.
        \param data Non-null address of serialized bytes.
        \param size Number of accessible bytes, representable as streamoff.
        \throws std::invalid_argument For a null address or unsupported size.
    */
    serialization_profile_streambuf(const unsigned char* data, size_t size);
protected:
    /** Seek within the borrowed input range.
        \param off Signed byte displacement from dir.
        \param dir Beginning, current position or end of input.
        \param which Requested stream mode, which must include input.
        \return New absolute position, or position -1 for an invalid request.
    */
    pos_type seekoff(off_type off, std::ios_base::seekdir dir,
                     std::ios_base::openmode which) override;
    /** Seek to an absolute input position.
        \param pos Requested position from the input beginning.
        \param which Requested stream mode, which must include input.
        \return New absolute position, or position -1 for an invalid request.
    */
    pos_type seekpos(pos_type pos, std::ios_base::openmode which) override;
};

/** \brief Internal reusable collector for completed bvector BLOBs.
    \ingroup bvserial
    \internal
    \tparam BV Application bvector type, including its allocator.
    \note Owns retained decoding scratch. A collector is exclusive mutable
        state; use separate instances for concurrent profiling. ID metrics are
        disabled by default to avoid the additional full reconstruction pass.
*/
template<class BV>
class serialization_opcode_collector
{
    bm::streams_deserializer<BV> reader_; ///< Retained stream decoder and codec scratch.
    bm::deserialization_index<BV> index_; ///< Validated top-level record offsets.
    BV scratch_; ///< Reconstructed bits when collecting original-ID metrics.
    bool collect_ids_ = false; ///< Enable reconstruction and cardinality accounting.
    typename bm::deserializer<BV, bm::decoder>::block_starts_vect_type starts_; ///< Block starts.
    bm::heap_vector<std::uint64_t, typename BV::allocator_type, true> record_ids_; ///< Original bit counts.
public:
    /** Borrow reference planes for reference-XOR reconstruction.
        \param refs Reference vectors, or null to unset them. Borrowed vectors
            must outlive their use by this collector.
        \note Required for reference-XOR BLOBs when collecting ID metrics.
    */
    void set_ref_vectors(const typename bm::serializer<BV>::bv_ref_vector_type* refs);
    /** Enable or disable exact original-ID cardinalities.
        \param enabled True adds a full reconstruction pass for each BLOB.
    */
    void set_collect_ids(bool enabled = true);
    /** Validate and accumulate measurements for one complete BLOB.
        \param data Non-null address of the BLOB; borrowed for this call.
        \param size Exact, nonzero BLOB size in bytes.
        \param profile Aggregate receiving counts, bytes and optional ID metrics.
        \throws std::exception For invalid input or decoding/allocation failure.
        \note Rejects truncated or extra input through exact-consumption checks.
    */
    void add_blob(const unsigned char* data, size_t size,
                  bm::serialization_opcode_profile& profile);
    /** Validate, accumulate and observe final top-level record spans.
        \tparam Observer Callable accepting a record address and byte length.
        \param data Non-null address of the complete BLOB, borrowed for this call.
        \param size Exact, nonzero BLOB size in bytes.
        \param profile Aggregate receiving counts, bytes and optional ID metrics.
        \param observe Called once per opcode record after span validation;
            pointers are borrowed from data. ID-list BLOBs have no such records.
        \throws std::exception For invalid input or decoding/allocation failure.
        \note The observer must not modify the BLOB or reenter this collector.
            If it throws, profile may already contain partial accumulation.
    */
    template<class Observer>
    void add_blob(const unsigned char* data, size_t size,
                  bm::serialization_opcode_profile& profile, Observer observe);
};

// ----------------------------------------------------------------------

inline
void serialization_opcode_profile::add(const bm::serialization_opcode_profile& r) noexcept
{
    for (unsigned i = 0; i < 256; ++i)
    {
        opcode[i].count += r.opcode[i].count;
        opcode[i].bytes += r.opcode[i].bytes;
        opcode[i].ids += r.opcode[i].ids;
    } // for i
    blobs += r.blobs;
    bytes += r.bytes;
    header_bytes += r.header_bytes;
    id_list_bytes += r.id_list_bytes;
    records += r.records;
    ids += r.ids;
    id_list_ids += r.id_list_ids;
    cardinality_blobs += r.cardinality_blobs;
}

inline
bool serialization_opcode_profile::cardinalities_reconcile() const noexcept
{
    std::uint64_t total = id_list_ids;
    for (const auto& e : opcode)
    {
        total += e.ids;
    } // for e
    return cardinality_blobs == blobs && total == ids;
}

inline
std::uint64_t serialization_opcode_profile::record_bytes() const noexcept
{
    std::uint64_t total = 0;
    for (const auto& e : opcode)
    {
        total += e.bytes;
    } // for e
    return total;
}

inline
bool serialization_opcode_profile::reconciles() const noexcept
{
    std::uint64_t count = 0;
    for (const auto& e : opcode)
    {
        count += e.count;
        if (e.bytes < e.count)
        {
            return false;
        }
    } // for e
    return count == records && header_bytes + id_list_bytes + record_bytes() == bytes;
}

inline
serialization_profile_streambuf::serialization_profile_streambuf(
    const unsigned char* data, size_t size)
{
    if (!data || size > size_t((std::numeric_limits<std::streamoff>::max)()))
    {
        throw std::invalid_argument("Invalid serialization profile input");
    }
    char* begin = reinterpret_cast<char*>(const_cast<unsigned char*>(data));
    setg(begin, begin, begin + size);
}

inline
serialization_profile_streambuf::pos_type serialization_profile_streambuf::seekoff(
    off_type off, std::ios_base::seekdir dir, std::ios_base::openmode which)
{
    if (!(which & std::ios_base::in))
    {
        return pos_type(off_type(-1));
    }
    const off_type size = egptr() - eback();
    off_type base = dir == std::ios_base::beg ? 0 :
                    dir == std::ios_base::cur ? gptr() - eback() :
                    dir == std::ios_base::end ? size : -1;
    if (base < 0 || off < -base || off > size - base)
    {
        return pos_type(off_type(-1));
    }
    const off_type pos = base + off;
    setg(eback(), eback() + pos, egptr());
    return pos_type(pos);
}

inline
serialization_profile_streambuf::pos_type serialization_profile_streambuf::seekpos(
    pos_type pos, std::ios_base::openmode which)
{
    return seekoff(off_type(pos), std::ios_base::beg, which);
}

template<class BV>
inline
void serialization_opcode_collector<BV>::set_ref_vectors(
    const typename bm::serializer<BV>::bv_ref_vector_type* refs)
{
    reader_.set_ref_vectors(refs);
}

template<class BV>
inline
void serialization_opcode_collector<BV>::set_collect_ids(bool enabled)
{
    collect_ids_ = enabled;
}

template<class BV>
inline
void serialization_opcode_collector<BV>::add_blob(
    const unsigned char* data, size_t size, bm::serialization_opcode_profile& profile)
{
    add_blob(data, size, profile, [](const unsigned char*, size_t) {});
}

template<class BV>
template<class Observer>
inline
void serialization_opcode_collector<BV>::add_blob(
    const unsigned char* data, size_t size, bm::serialization_opcode_profile& profile,
    Observer observe)
{
    if (!data || !size)
    {
        throw std::invalid_argument("Empty serialization profile BLOB");
    }
    bm::serialization_profile_streambuf view(data, size);
    std::istream stream(&view);
    bm::streams_decoder input(stream);
    reader_.unset_deserialization_index();
    reader_.unset_block_digest_vector();
    reader_.set_deserialization_index_construct(&index_);
    scratch_.clear();
    if (!reader_.deserialize(scratch_, input) || input.tell() != size ||
        index_.serialized_size() != size)
    {
        throw std::logic_error("Serialization profile BLOB did not decode completely");
    }
    reader_.unset_deserialization_index();
    std::uint64_t blob_ids = 0;
    if (collect_ids_)
    {
        bm::serialization_profile_streambuf decoded_view(data, size);
        std::istream decoded_stream(&decoded_view);
        bm::streams_decoder decoded_input(decoded_stream);
        scratch_.clear();
        reader_.set_record_block_starts(&starts_);
        bool ok;
        try
        {
            ok = reader_.deserialize(scratch_, decoded_input);
        }
        catch (...)
        {
            reader_.set_record_block_starts(nullptr);
            throw;
        }
        reader_.set_record_block_starts(nullptr);
        if (!ok || decoded_input.tell() != size || starts_.size() != index_.size())
        {
            throw std::logic_error("Opcode cardinality traversal differs from index");
        }
        blob_ids = scratch_.count();
        record_ids_.resize(starts_.size());
        if (!record_ids_.empty())
        {
            ::memset(record_ids_.data(), 0, record_ids_.size() * sizeof(std::uint64_t));
        }
        std::uint64_t assigned = 0;
        for (size_t i = 0; i < starts_.size(); ++i)
        {
            const std::uint64_t first = std::uint64_t(starts_[i]) << bm::set_block_shift;
            const std::uint64_t end = i + 1 < starts_.size() ?
                std::uint64_t(starts_[i + 1]) << bm::set_block_shift :
                std::uint64_t(scratch_.size());
            if (end < first && i + 1 < starts_.size())
            {
                throw std::logic_error("Non-monotonic opcode block ranges");
            }
            if (end > first && first < scratch_.size())
            {
                record_ids_[i] = scratch_.count_range(typename BV::size_type(first),
                    typename BV::size_type((std::min)(end,
                        std::uint64_t(scratch_.size())) - 1));
            }
            assigned += record_ids_[i];
        } // for i
        if (!(data[0] & bm::BM_HM_ID_LIST) && assigned != blob_ids)
        {
            throw std::logic_error("Opcode original-ID counts do not reconcile");
        }
    }
    if (data[0] & bm::BM_HM_ID_LIST)
    {
        if (!index_.empty())
        {
            throw std::logic_error("ID-list profile has block markers");
        }
        profile.bytes += size;
        ++profile.blobs;
        profile.id_list_bytes += size;
        if (collect_ids_)
        {
            ++profile.cardinality_blobs;
            profile.ids += blob_ids;
            profile.id_list_ids += blob_ids;
        }
        return;
    }
    size_t marker_index = 0, marker = 0;
    if (!index_.find_next_marker_offset(0, marker_index, marker) || !marker || marker >= size)
    {
        throw std::logic_error("Missing serialization profile marker");
    }
    size_t local_bytes = marker;
    const size_t header_bytes = marker;
    while (marker < size)
    {
        size_t next = size;
        if (!index_.find_next_marker_offset(marker, marker_index, next))
        {
            next = size;
        }
        if (next <= marker || next > size)
        {
            throw std::logic_error("Invalid serialization profile marker order");
        }
        local_bytes += next - marker;
        marker = next;
    }
    if (local_bytes != size)
    {
        throw std::logic_error("Serialization profile bytes do not reconcile");
    }
    profile.bytes += size;
    ++profile.blobs;
    profile.header_bytes += header_bytes;
    if (collect_ids_)
    {
        ++profile.cardinality_blobs;
        profile.ids += blob_ids;
    }
    size_t record_index = 0;
    marker = header_bytes;
    marker_index = 0;
    while (marker < size)
    {
        size_t next = size;
        if (!index_.find_next_marker_offset(marker, marker_index, next))
        {
            next = size;
        }
        auto& entry = profile.opcode[data[marker]];
        if (collect_ids_)
        {
            entry.ids += record_ids_[record_index++];
        }
        ++entry.count;
        entry.bytes += next - marker;
        ++profile.records;
        observe(data + marker, next - marker);
        marker = next;
    }
}

} // namespace bm

#endif // BM_SERIAL_PROFILE_H_INCLUDED
