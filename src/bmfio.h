#ifndef BMFIO__H__INCLUDED__
#define BMFIO__H__INCLUDED__
/*
Copyright(c) 2026 Anatoliy Kuznetsov
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy at http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing, software distributed
under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
CONDITIONS OF ANY KIND, either express or implied.
*/

/** @file bmfio.h
    @brief Buffered bvector serialization and deserialization over C++ streams.
*/

#include <ostream>
#include <istream>
#include <limits>
#include "bmdef.h"
#include "bmconst.h"
#include "bmsparsevec_serial.h"
#include "bmsparsevec_float_serial.h"


namespace bm
{
/** Buffered serialization to a seekable binary C++ output stream.

    Open files in binary output mode, without ios::app. The stream is borrowed
    and must not be repositioned or written independently while this channel is
    in use. Offsets are relative to its initial position. No framing is added.

    @ingroup bvserial
    @sa streams_decoder

    The serializer owns working memory; this channel owns its RAM encoder.
    flush() releases submitted memory for immediate reuse. patch() also retains
    no caller memory. An asynchronous implementation of this protocol must copy
    bytes before returning, preserve write/patch order, and join pending work in
    finish(). unbind_buffer() must be non-throwing, including during unwinding.

    finish() is caller-owned, repeatable and non-closing. Here it synchronizes
    the ostream buffer to report delayed stream errors; it does not imply disk
    durability. Stream exceptions propagate. Otherwise errors are sticky and
    is_good() reports them; stream().rdstate() provides stream-specific details.
*/
class streams_encoder
{
public:
    /** Bind a seekable output stream and capture its current position as origin.
        @param stream Borrowed binary output stream, opened without append mode.
        @note Stream exceptions propagate; otherwise is_good() reports initialization failure.
    */
    explicit streams_encoder(std::ostream& stream);

    /** Inspect channel and stream state.
        @return True if no sticky channel error exists and the stream is good.
    */
    bool is_good() const BMNOEXCEPT;

    /** Access the borrowed stream for channel-specific error interrogation.
        @return The underlying output stream; do not write or seek it independently during use.
    */
    std::ostream& stream() BMNOEXCEPT;

    /** Measure all encoded bytes relative to the channel origin.
        @return Submitted bytes plus bytes still pending in the RAM encoder.
    */
    size_t size() const BMNOEXCEPT;

    /** Access the RAM encoder used by serialization helpers.
        @return Stable encoder reference; its buffer binding may change after reserve().
    */
    bm::encoder& get_encoder() BMNOEXCEPT;

    /** Attach serializer-owned working memory without writing to the stream.
        @param data Writable buffer retained until unbind_buffer().
        @param capacity Available buffer size in bytes.
        @pre No buffer is currently bound.
    */
    void bind_buffer(unsigned char* data, size_t capacity) BMNOEXCEPT;

    /** Release the borrowed buffer and reset the RAM encoder without I/O.
        @note Pending bytes must already have been submitted unless abandoning a failed call.
    */
    void unbind_buffer() BMNOEXCEPT;

    /** Write pending encoded bytes and reset the RAM cursor.
        @return True if all pending bytes were accepted, or there were none.
        @note Success permits immediate buffer reuse. Stream exceptions propagate;
        otherwise failures are sticky and return false.
    */
    bool flush();

    /** Replace previously encoded bytes while preserving the append position.
        @param offset First byte to replace, relative to the channel origin.
        @param data Replacement bytes, borrowed only for this call.
        @param count Number of replacement bytes.
        @return True if the replacement succeeded; false on a range or I/O error.
        @note Patches in pending memory avoid I/O. Other patches flush, seek, write,
        and restore the append position. Stream exceptions propagate.
    */
    bool patch(size_t offset, const unsigned char* data, size_t count);

    /** Copy construction is disabled because the stream and buffer are borrowed.
        @param other Channel that must not be copied.
    */
    streams_encoder(const streams_encoder& other) = delete;

    /** Copy assignment is disabled to preserve unique channel state.
        @param other Channel that must not be assigned.
        @return No value can be returned; this operation is deleted.
    */
    streams_encoder& operator=(const streams_encoder& other) = delete;

    /** Grow the serializer-owned buffer after submitting pending bytes.
        @tparam Buffer Byte buffer providing resize(), data(), and size().
        @param buffer Owner of the currently bound working memory.
        @param capacity Minimum required capacity in bytes.
        @return True if the bound buffer has sufficient capacity; false on I/O failure.
        @note Allocation and stream exceptions propagate. The old borrow is released
        before allocation so that unwinding cannot retain a dangling pointer.
    */
    template<class Buffer>
    bool reserve(Buffer& buffer, size_t capacity);

    /** Complete caller-owned output synchronization without closing the stream.
        @return True if pending bytes and the stream buffer were flushed successfully.
        @note Repeatable and non-closing; does not promise disk durability. Stream
        exceptions propagate. Additional serialization calls may follow success.
    */
    bool finish();

private:
    /** Check conversion to the stream offset and transfer-count types.
        @param n Unsigned byte count or offset.
        @return True if both streamoff and streamsize can represent n.
    */
    static bool fits(size_t n) BMNOEXCEPT;

    /** Check whether an origin-relative position fits the stream offset type.
        @param n Byte offset relative to origin_.
        @return True if the nonnegative absolute position can be represented.
    */
    bool position_fits(size_t n) const BMNOEXCEPT;

    std::ostream& stream_;       ///< Borrowed output stream.
    std::streampos origin_;       ///< Initial physical output position.
    bm::encoder enc_;             ///< RAM encoder bound to serializer memory.
    unsigned char* data_;         ///< Borrowed working buffer, or null.
    size_t capacity_;             ///< Bound buffer capacity in bytes.
    size_t written_;              ///< Bytes submitted since origin_.
    bool good_;                   ///< Sticky channel success state.
};

/** @internal Bounds-checked RAM decoding; never performs I/O.
    Scalar and array coding delegates to the original RAM decoder.
 */
class bounded_decoder
{
public:
    /** Signal an attempted read or position outside the bound RAM window.
        @note This internal exception carries no stream I/O state.
    */
    struct exhausted {};

    /** Construct an unbound, empty RAM decoder without allocating memory.
    */
    bounded_decoder() BMNOEXCEPT;

    /** Bind a readable RAM window and reset its cursor.
        @param p Borrowed input memory; may be null for an empty window.
        @param n Number of valid bytes in the window.
    */
    void bind(const unsigned char* p, size_t n) BMNOEXCEPT;

    /** Inspect bytes consumed from the current window.
        @return Cursor offset relative to the bound memory.
    */
    size_t size() const BMNOEXCEPT;

    /** Inspect the readable suffix of the window.
        @return Number of bytes not yet consumed.
    */
    size_t remaining() const BMNOEXCEPT;

    /** Access the current RAM position.
        @return Pointer to the next byte, or the end of the bound window.
        @pre A memory window is bound.
    */
    const unsigned char* get_pos() const BMNOEXCEPT;

    /** Move the RAM cursor without reading or performing I/O.
        @param p Position within the bound window, including its end.
        @throws exhausted If p is outside the window.
    */
    void set_pos(const unsigned char* p);

    /** Move the RAM cursor by a signed byte displacement.
        @param delta Relative displacement; negative values move backward.
        @throws exhausted If the resulting position is outside the window.
    */
    void seek(int delta);

    /** Read one 8-bit unsigned value and advance the RAM cursor.
        @return Decoded value using the existing RAM byte-order convention.
        @throws exhausted If fewer than 1 readable bytes remain.
    */
    unsigned char get_8();

    /** Read one 16-bit unsigned value and advance the RAM cursor.
        @return Decoded value using the existing RAM byte-order convention.
        @throws exhausted If fewer than 2 readable bytes remain.
    */
    bm::short_t get_16();

    /** Read one 24-bit unsigned value and advance the RAM cursor.
        @return Decoded value using the existing RAM byte-order convention.
        @throws exhausted If fewer than 3 readable bytes remain.
    */
    bm::word_t get_24();

    /** Read one 32-bit unsigned value and advance the RAM cursor.
        @return Decoded value using the existing RAM byte-order convention.
        @throws exhausted If fewer than 4 readable bytes remain.
    */
    bm::word_t get_32();

    /** Read one 48-bit unsigned value and advance the RAM cursor.
        @return Decoded value using the existing RAM byte-order convention.
        @throws exhausted If fewer than 6 readable bytes remain.
    */
    bm::id64_t get_48();

    /** Read one 64-bit unsigned value and advance the RAM cursor.
        @return Decoded value using the existing RAM byte-order convention.
        @throws exhausted If fewer than 8 readable bytes remain.
    */
    bm::id64_t get_64();

    /** Read the mask-prefixed h64 representation from the RAM window.
        @return Decoded 64-bit unsigned value.
        @throws exhausted If the mask or a selected value byte is unavailable.
    */
    bm::id64_t get_h64();

    /** Read an array of 16-bit values using the original RAM decoder.
        @param dst Destination array, or null to skip the encoded values.
        @param count Number of array elements to consume; must be positive.
        @throws exhausted If the complete encoded array does not fit in the window.
    */
    void get_16(bm::short_t* dst, unsigned count);

    /** Read an array of 32-bit values using the original RAM decoder.
        @param dst Destination array, or null to skip the encoded values.
        @param count Number of array elements to consume; must be positive.
        @throws exhausted If the complete encoded array does not fit in the window.
    */
    void get_32(bm::word_t* dst, unsigned count);

    /** Decode 32-bit words and combine them into an existing array with OR.
        @param dst Writable destination containing count words, or null to skip input.
        @param count Number of words to read and combine, normally bm::set_block_size.
        @pre Count and destination alignment satisfy the original RAM decoder contract.
        @return True if every resulting word is all ones; false when dst is null.
        @throws exhausted If the complete encoded array is unavailable.
    */
    bool get_32_OR(bm::word_t* dst, unsigned count);

    /** Decode 32-bit words and combine them into an existing array with AND.
        @param dst Writable destination containing count words, or null to skip input.
        @param count Number of words to read and combine, normally bm::set_block_size.
        @pre Count and destination alignment satisfy the original RAM decoder contract.
        @throws exhausted If the complete encoded array is unavailable.
    */
    void get_32_AND(bm::word_t* dst, unsigned count);

    /** Copy raw bytes and advance the RAM cursor.
        @param dst Destination buffer, or null to consume without copying.
        @param count Number of bytes to consume.
        @throws exhausted If count exceeds the readable suffix.
    */
    void memcpy(unsigned char* dst, size_t count);

private:
    /** Validate a byte-count request without advancing the cursor.
        @param n Required readable bytes.
        @throws exhausted If n exceeds remaining().
    */
    void require(size_t n) const;

    /** Validate array size without overflowing a count-times-width calculation.
        @param n Required element count.
        @param width Encoded bytes per element; must be nonzero.
        @throws exhausted If the array exceeds the readable suffix.
    */
    void require_array(unsigned n, unsigned width) const;

    const unsigned char* data_;  ///< Start of the borrowed RAM window.
    size_t count_;                ///< Valid bytes in the bound window.
    size_t pos_;                  ///< Current byte offset within the window.
};

/** @internal Allow checked RAM exhaustion to propagate through bit_in helpers.
    This specialization leaves the legacy RAM decoder's noexcept behavior intact.
*/
template<> struct decoder_noexcept<bounded_decoder>
{
    enum { value = false }; ///< Checked decoding may throw exhausted.
};

/** Buffered input from a seekable binary C++ stream.
    @ingroup bvserial
    @sa streams_deserializer

    No framing is added to the bvector format. The deserializer lends its retained
    working buffer for each call; this channel owns the RAM decoder. prepare()
    and seek() are the only reading/repositioning boundaries. RAM primitives and
    bit_in never perform I/O. A successful complete() restores the actual stream
    position to the byte following the BLOB, discarding any buffered read-ahead.

    Offsets are relative to the BLOB origin captured by begin(). Between completed
    deserialize calls the caller may read a custom payload or reposition stream().
    Do not access the stream independently during a deserialize call.

    Stream exceptions propagate. Other failures return false and remain sticky;
    error() distinguishes invalid/truncated input from stream I/O failure. Partial
    deserialization may have modified the destination vector on failure.
*/
class streams_decoder
{
public:
    enum { streaming = true }; ///< Select exact physical BLOB-end accounting.

    /** Persistent channel error categories, independent of stream exception masks. */
    enum error_code
    {
        no_error,      ///< No channel error has been recorded.
        io_error,      ///< Input, positioning, or stream-state failure.
        invalid_input  ///< Invalid request or exhausted encoded input.
    };

    /** Borrow an input stream without reading or capturing a BLOB origin yet.
        @param stream Seekable binary stream retained for the channel lifetime.
    */
    explicit streams_decoder(std::istream& stream) BMNOEXCEPT;

    /** Inspect both sticky channel status and the current stream state.
        @return True if no channel error exists and the stream is good.
    */
    bool is_good() const BMNOEXCEPT;

    /** Inspect the first recorded channel error.
        @return Persistent error category; inspect stream().rdstate() for stream details.
    */
    error_code error() const BMNOEXCEPT;

    /** Access the borrowed stream for errors or between-record application payloads.
        @return Underlying input stream; independent access is forbidden during decoding.
    */
    std::istream& stream() BMNOEXCEPT;

    /** Access the checked RAM decoder used by the shared traversal.
        @return Stable decoder reference whose memory binding may change at I/O boundaries.
    */
    bounded_decoder& get_decoder() BMNOEXCEPT;

    /** Inspect the logical input cursor independently of physical read-ahead.
        @return Byte offset relative to the BLOB origin captured by begin().
    */
    size_t tell() const BMNOEXCEPT;

    /** Start a BLOB at the stream current position and borrow working memory.
        @param data Writable input buffer owned by the deserializer.
        @param capacity Buffer capacity in bytes.
        @return True if origin capture and binding succeeded; false on stream failure.
        @note Resets the window and relative cursor, but never clears a sticky error.
        Stream exceptions propagate.
    */
    bool begin(unsigned char* data, size_t capacity);

    /** Release borrowed memory and empty the RAM window without performing I/O.
        @note Safe during exception unwinding; does not clear channel errors.
    */
    void unbind_buffer() BMNOEXCEPT;

    /** Record a channel failure while preserving the first error.
        @param e Error category to record if no earlier error exists.
        @return Always false, for direct propagation from a failed operation.
    */
    bool fail(error_code e = invalid_input) BMNOEXCEPT;

    /** Make an upcoming encoded record available in the RAM window.
        @param n Desired readable bytes; must not exceed the bound capacity.
        @return True after a successful refill or at normal EOF; false on channel failure.
        @note Preserves unread bytes and fills the remaining buffer. A final partial
        window is valid: checked RAM reads detect exhaustion if the record needs more
        bytes. EOF/fail exceptions caused only by this short read are handled locally;
        other stream exceptions propagate.
    */
    bool prepare(size_t n);

    /** Move the logical cursor to a BLOB-relative position.
        @param pos Destination byte offset relative to the current BLOB origin.
        @return True on success; false if positioning fails.
        @note An in-window seek adjusts only the RAM cursor. Other seeks reposition the
        stream and invalidate the window. Stream exceptions propagate.
    */
    bool seek(size_t pos);

    /** Restore the underlying stream to the logical cursor after decoding.
        @return True if the final physical seek succeeded; false on stream failure.
        @note Discards read-ahead so the caller can read a trailing application payload.
        The traversal must already have located the exact BLOB end. Does not close the
        stream; stream exceptions propagate.
    */
    bool complete();

    /** Copy construction is disabled because the stream and buffer are borrowed.
        @param other Channel that must not be copied.
    */
    streams_decoder(const streams_decoder& other) = delete;

    /** Copy assignment is disabled to preserve unique channel state.
        @param other Channel that must not be assigned.
        @return No value can be returned; this operation is deleted.
    */
    streams_decoder& operator=(const streams_decoder& other) = delete;

    /** Grow the owner buffer while preserving unread input and its logical offset.
        @tparam Buffer Byte buffer providing resize(), data(), and size().
        @param buffer Owner of the currently bound memory.
        @param capacity Minimum required capacity in bytes.
        @return True if capacity is sufficient after rebinding; false on channel failure.
        @note Does not perform I/O. Allocation exceptions propagate after releasing the
        old borrow, preventing a dangling pointer during unwinding.
    */
    template<class Buffer>
    bool reserve(Buffer& buffer, size_t capacity);

private:
    /** Seek the physical stream using a checked origin-relative offset.
        @param pos Destination byte offset relative to origin_.
        @return True if seekg succeeded; false on offset overflow or I/O failure.
        @note Marks I/O failure before the potentially throwing stream operation.
    */
    bool seek_stream(size_t pos);

    std::istream& stream_;                           ///< Borrowed input stream.
    std::streampos origin_ = std::streampos(0);        ///< Current BLOB physical origin.
    bounded_decoder dec_;                            ///< Checked decoder over borrowed memory.
    unsigned char* data_ = 0;                         ///< Deserializer-owned buffer, or null.
    size_t capacity_ = 0;                            ///< Bound capacity in bytes.
    size_t offset_ = 0;                              ///< BLOB-relative start of the RAM window.
    bool eof_ = false;                               ///< Last refill reached physical EOF.
    error_code error_ = no_error;                    ///< First recorded channel failure.
};

/** Streaming bvector deserialization, including range and indexed gather.
    @ingroup bvserial
    @sa streams_decoder, deserializer
    @tparam BV Bit-vector type supplying its allocator and decoding operations.

    Uses the same traversal and coding as the RAM deserializer. Settings and
    reference/index attachment APIs are inherited. deserialize() returns true
    only after the underlying stream has been restored to the exact BLOB end.
    Gather and range calls parse any remaining tail without materializing it.
    An index built by this interface records the exact physical BLOB size.

    The BM allocator-typed scratch buffer is retained across calls. Its
    initial size is 65 KiB. It grows to 257 KiB for legacy arrays/gamma codes
    and sparse superblocks (65536 unsigned entries), independently of the
    bvector's logical size.
    No sparse-vector format or implementation is changed by this interface.

    Example:
    @code
        bm::streams_decoder input(file); // Borrow seekable binary input.
        bm::streams_deserializer<BV> reader; // Retain reusable working memory.
        bool ok = reader.deserialize(bv, input); // Decode and restore the BLOB end.
        // On success, file can now read an application-specific trailing payload.
    @endcode
*/
template<class BV>
class streams_deserializer : public bm::deserializer<BV, bm::bounded_decoder>
{
    /** Shared RAM/stream traversal with checked primitive decoding. */
    typedef bm::deserializer<BV, bm::bounded_decoder> parent_type;
public:
    /** Working-memory owner using the bit-vector allocator contract. */
    typedef bm::byte_buffer<typename BV::allocator_type> buffer_type;

    /** Deserialize one BLOB and restore the underlying stream to its exact end.
        @tparam IN Buffered input channel implementing the streams_decoder protocol.
        @param bv Destination vector; decoding has the inherited OR semantics.
        @param input Channel borrowing the stream and this reader's working buffer.
        @return True only if decoding and final positioning both succeeded.
        @note Exhausted RAM input records failure and returns false. Other
        exceptions propagate after recording failure. The destination may be
        partially modified on failure; its allocation strategy is restored and
        the buffer borrow released on every exit. Range and gather settings are
        inherited; final completion also parses any unrequested tail.
    */
    template<class IN> bool deserialize(BV& bv, IN& input);

private:
    /** @internal Restore destination strategy and release memory on every exit.
        @tparam IN Input channel whose buffer must be unbound without throwing.
    */
    template<class IN> struct release
    {
        IN& in;                 ///< Channel borrowing the retained working memory.
        BV& target;             ///< Destination whose original strategy is restored.
        bm::strategy strategy;  ///< Allocation strategy saved at entry.

        /** Restore the target strategy and unbind channel memory without I/O. */
        ~release();
    };

    /** @internal Adapt an input channel and buffer owner to the shared traversal.
        @tparam IN Input channel with checked RAM decoding and explicit I/O boundaries.
    */
    template<class IN> struct source
    {
        enum { streaming = true }; ///< Request complete physical BLOB consumption.

        /** Bind references to the channel and retained working buffer.
            @param input Channel used for reads, positioning, and status.
            @param memory Owner of the channel's working memory.
        */
        source(IN& input, buffer_type& memory);

        /** Access the checked RAM decoder without reading the stream.
            @return Decoder owned by the input channel.
        */
        bm::bounded_decoder& get_decoder();

        /** Ensure buffer capacity, then prepare the next RAM decoding window.
            @param n Maximum bytes needed for the upcoming encoded record.
            @return True if reserve() and prepare() succeed; false on channel failure.
        */
        bool prepare(size_t n);

        /** Forward a BLOB-relative cursor move to the channel.
            @param pos Destination offset relative to the current BLOB origin.
            @return True if channel positioning succeeded.
        */
        bool seek(size_t pos);

        /** Read the logical cursor without physical stream I/O.
            @return Byte offset relative to the current BLOB origin.
        */
        size_t tell() const;

        IN& in;               ///< Borrowed channel used by the traversal.
        buffer_type& buffer;  ///< Retained memory owner used when capacity grows.
    };
    buffer_type buffer_;      ///< Allocator-typed memory retained across BLOBs.
};

/** Buffered restore of sparse vectors from seekable binary streams.
    Supports sparse_vector, rsc_sparse_vector and str_sparse_vector in the current digest-table
    format and native byte order. Successful calls leave the stream immediately
    after the object. Older inline-offset headers are not supported here.
    Inherited index attachment and XOR-reference settings apply to this reader.
    I/O errors return false unless the stream is configured to throw. Invalid
    metadata returns false; other decoding exceptions propagate. Targets may
    be partially modified on failure.
    @tparam SV Sparse-vector type.
*/
template<class SV>
class streams_sparse_vector_deserializer : public sparse_vector_deserializer<SV>
{
    typedef sparse_vector_deserializer<SV> base_type; ///< Shared sparse-vector state.
public:
    typedef typename SV::bvector_type bvector_type; ///< Plane and mask type.
    typedef typename base_type::deserialization_index_type deserialization_index_type; ///< Plane index collection.
    /** Restore an object or gather selected logical addresses.
        @param sv Mutable destination with matching NULL support.
        @param stream Seekable input positioned at the object start.
        @param mask Optional address mask, distinct from destination planes.
        @return True on successful decode and final positioning.
    */
    bool deserialize(SV& sv, std::istream& stream, const bvector_type* mask = 0);
    /** Construct reusable plane indexes without restoring the full vector.
        @param index Output indexes for this exact BLOB.
        @param stream Seekable input positioned at the object start.
        @return True on successful index construction and final positioning.
    */
    bool construct_deserialization_index(deserialization_index_type& index,
                                         std::istream& stream);
private:
    /** Read metadata, then decode planes in dependency order.
        @param sv Destination or index-construction scratch vector.
        @param stream Borrowed input stream.
        @param mask Optional logical selection.
        @param build Optional index construction target.
        @return True if the complete operation succeeded.
    */
    bool read(SV& sv, std::istream& stream, const bvector_type* mask,
              deserialization_index_type* build);
    /** Seek to a checked object-relative offset.
        @param stream Borrowed stream.
        @param origin Physical object start.
        @param offset Object-relative byte offset.
        @return True if seeking succeeded.
    */
    static bool seek(std::istream& stream, std::streampos origin, size_t offset);
    streams_deserializer<bvector_type> reader_; ///< Retained plane decoder buffer.
    typename serializer<bvector_type>::buffer metadata_; ///< Small offset-table buffer.
};

/** Buffered restore of legacy bf0 float BLOBs on the native architecture.
    Successful calls leave the stream after the complete composite object.
    Failed calls return false or propagate configured stream/decoder exceptions;
    the destination may be partially modified. Shared NULL planes are restored
    before mantissas, including rank-compressed storage.
    @tparam SV Float sparse-vector specialization.
*/
template<class SV>
class streams_sparse_vector_float_deserializer
{
public:
    typedef typename SV::bvector_type bvector_type; ///< Address mask type.
    typedef sparse_vector_float_deserialization_index<SV> deserialization_index_type; ///< Composite index.
    typedef typename serializer<bvector_type>::bv_ref_vector_type bv_ref_vector_type; ///< XOR references.
    /** Restore all values or a logical gather.
        @param sv Mutable target with matching NULL support.
        @param stream Seekable input at the BLOB start.
        @param mask Optional logical address mask.
        @return True after successful restoration and final positioning.
    */
    bool deserialize(SV& sv, std::istream& stream, const bvector_type* mask = 0);
    /** Build indexes for the two integer components.
        @param index Output index for this exact BLOB.
        @param stream Seekable input at the BLOB start.
        @return True after successful construction and final positioning.
    */
    bool construct_deserialization_index(deserialization_index_type& index, std::istream& stream);
    /** Attach a borrowed index. @param index Matching index or null to detach. */
    void set_deserialization_index(deserialization_index_type* index);
    /** Control index-assisted reads. @param enable Whether to use attached indexes. */
    void set_deserialization_index_use(bool enable = true);
    /** Attach external XOR references. @param refs References or null to detach. */
    void set_xor_ref(bv_ref_vector_type* refs);
private:
    /** Validate the native header and calculate absolute component boundaries.
        @param stream Input at the header.
        @param positions Output sign, exponent, mantissa and end positions.
        @return True for valid sizes within the seekable input extent.
    */
    static bool header(std::istream& stream, std::streampos (&positions)[4]);
    streams_deserializer<bvector_type> signs_; ///< Retained sign decoder.
    streams_sparse_vector_deserializer<typename SV::sparse_vector_u> exponents_; ///< Exponent decoder.
    streams_sparse_vector_deserializer<typename SV::sparse_vector_u> mantissas_; ///< Mantissa decoder.
};

namespace sv_di_detail
{
/** @internal Stream output adapter for composite index codecs. */
struct stream_output
{
    std::ostream& stream; ///< Borrowed stream.
    std::streampos origin; ///< Physical record origin.
    size_t pos = 0; ///< Bytes written relative to origin.
    /** Capture stream origin. @param s Seekable output. */
    explicit stream_output(std::ostream& s);
    /** Write bytes. @param p Source. @param n Length. @return I/O success. */
    bool write(const unsigned char* p, size_t n);
    /** Patch framing and restore position. @param at Relative offset. @param p Source.
        @param n Length. @return I/O success. */
    bool patch(size_t at, const unsigned char* p, size_t n);
};
/** @internal Stream input adapter with a checked physical extent. */
struct stream_input
{
    std::istream& stream; ///< Borrowed stream.
    size_t remaining = 0; ///< Remaining bytes within current record/physical extent.
    size_t pos = 0; ///< Bytes read relative to origin.
    /** Capture input extent and restore the initial position. @param s Seekable input. */
    explicit stream_input(std::istream& s);
    /** Read within the current bound. @param p Destination. @param n Length. @return I/O success. */
    bool read(unsigned char* p, size_t n);
};
} // namespace sv_di_detail

// --------------------------------------------
// Implementations for: streams_encoder
// --------------------------------------------

inline
streams_encoder::streams_encoder(std::ostream& stream)
        : stream_(stream), origin_(stream.tellp()), enc_(0, 0),
          data_(0), capacity_(0), written_(0), good_(false)
{
    good_ = stream_.good() && origin_ != std::streampos(-1);
}

inline
bool streams_encoder::is_good() const BMNOEXCEPT
{
    return good_ && stream_.good();
}

inline
std::ostream& streams_encoder::stream() BMNOEXCEPT
{
    return stream_;
}

inline
size_t streams_encoder::size() const BMNOEXCEPT
{
    return written_ + (data_ ? enc_.size() : 0);
}

inline
bm::encoder& streams_encoder::get_encoder() BMNOEXCEPT
{
    return enc_;
}

inline
void streams_encoder::bind_buffer(unsigned char* data, size_t capacity) BMNOEXCEPT
{
    BM_ASSERT(!data_);
    data_ = data; capacity_ = capacity;
    enc_ = bm::encoder(data, capacity);
}

inline
void streams_encoder::unbind_buffer() BMNOEXCEPT
{
    data_ = 0; capacity_ = 0;
    enc_ = bm::encoder(0, 0);
}

inline
bool streams_encoder::flush()
{
    if (!is_good()) return false;
    const size_t count = data_ ? enc_.size() : 0;
    if (!count) return true;
    if (!fits(count) || count > size_t(-1) - written_ ||
        !position_fits(written_ + count))
        return good_ = false;
    // Mark failure before potentially throwing stream operations.
    good_ = false;
    stream_.write(reinterpret_cast<const char*>(data_),
                  std::streamsize(count));
    if (!stream_.good()) return false;
    written_ += count;
    enc_.reset();
    return good_ = true;
}

inline
bool streams_encoder::patch(size_t offset, const unsigned char* data, size_t count)
{
    if (!is_good()) return false;
    const size_t total = size();
    if (offset > total || count > total - offset)
        return good_ = false;
    // A patch may target this region before it is submitted.
    if (offset >= written_)
    {
        if (count) ::memcpy(data_ + offset - written_, data, count);
        return true;
    }
    if (!flush() || !fits(offset) || !fits(written_) || !fits(count))
        return good_ = false;
    good_ = false;
    stream_.seekp(origin_ + std::streamoff(offset));
    if (!stream_.good()) return false;
    stream_.write(reinterpret_cast<const char*>(data), std::streamsize(count));
    if (!stream_.good()) return false;
    stream_.seekp(origin_ + std::streamoff(written_));
    return good_ = stream_.good();
}



inline
bool streams_encoder::finish()
{
    if (!flush()) return false;
    good_ = false;
    stream_.flush();
    return good_ = stream_.good();
}

inline
bool streams_encoder::fits(size_t n) BMNOEXCEPT
{
    return n <= size_t((std::numeric_limits<std::streamoff>::max)()) &&
           n <= size_t((std::numeric_limits<std::streamsize>::max)());
}

inline
bool streams_encoder::position_fits(size_t n) const BMNOEXCEPT
{
    const std::streamoff base = std::streamoff(origin_);
    return base >= 0 && n <= size_t(
        (std::numeric_limits<std::streamoff>::max)() - base);
}


template<class Buffer>
bool streams_encoder::reserve(Buffer& buffer, size_t capacity)
{
    if (capacity <= capacity_) return is_good();
    if (!flush()) return false;
    // Clear the borrowed pointer before allocation (which may throw).
    unbind_buffer();
    buffer.resize(capacity, false);
    bind_buffer(buffer.data(), buffer.size());
    return true;
}

// --------------------------------------------
// Implementations for: bounded_decoder
// --------------------------------------------

inline
bounded_decoder::bounded_decoder() BMNOEXCEPT : data_(0), count_(0), pos_(0)
{
    
}

inline
void bounded_decoder::bind(const unsigned char* p, size_t n) BMNOEXCEPT
{
    data_ = p;
    count_ = n;
    pos_ = 0;
}

inline
size_t bounded_decoder::size() const BMNOEXCEPT
{
    return pos_;
}

inline
size_t bounded_decoder::remaining() const BMNOEXCEPT
{
    return count_ - pos_;
}

inline
const unsigned char* bounded_decoder::get_pos() const BMNOEXCEPT
{
    return data_ + pos_;
}

inline
void bounded_decoder::set_pos(const unsigned char* p)
{
    if (p < data_ || p > data_ + count_) throw exhausted();
    pos_ = size_t(p - data_);
}

inline
void bounded_decoder::seek(int delta)
{
    if (delta < 0)
    {
        size_t n = size_t(-static_cast<long long>(delta));
        if (n > pos_) throw exhausted();
        pos_ -= n;
    }
    else
    {
        require(size_t(delta));
        pos_ += size_t(delta);
    }
}

inline
unsigned char bounded_decoder::get_8()
{
    require(1);
    return data_[pos_++];
}

inline
bm::short_t bounded_decoder::get_16()
{
    require(2);
    bm::decoder d(data_ + pos_);
    pos_ += 2;
    return d.get_16();
}

inline
bm::word_t bounded_decoder::get_24()
{
    require(3);
    bm::decoder d(data_ + pos_);
    pos_ += 3;
    return d.get_24();
}

inline
bm::word_t bounded_decoder::get_32()
{
    require(4);
    bm::decoder d(data_ + pos_);
    pos_ += 4;
    return d.get_32();
}

inline
bm::id64_t bounded_decoder::get_48()
{
    require(6);
    bm::decoder d(data_ + pos_);
    pos_ += 6;
    return d.get_48();
}

inline
bm::id64_t bounded_decoder::get_64()
{
    require(8);
    bm::decoder d(data_ + pos_);
    pos_ += 8;
    return d.get_64();
}

inline
bm::id64_t bounded_decoder::get_h64()
{
    bm::id64_t value = 0;
    unsigned mask = get_8();
    for (unsigned i = 0; i < 8; ++i)
        if (mask & (1u << i)) value |= bm::id64_t(get_8()) << (i * 8);
    return value;
}

inline
void bounded_decoder::get_16(bm::short_t* dst, unsigned count)
{
    require_array(count, 2);
    bm::decoder d(data_ + pos_); d.get_16(dst, count); pos_ += size_t(count) * 2;
}

inline
void bounded_decoder::get_32(bm::word_t* dst, unsigned count)
{
    require_array(count, 4);
    bm::decoder d(data_ + pos_); d.get_32(dst, count); pos_ += size_t(count) * 4;
}

inline
bool bounded_decoder::get_32_OR(bm::word_t* dst, unsigned count)
{
    require_array(count, 4);
    bm::decoder d(data_ + pos_);
    bool full = d.get_32_OR(dst, count); pos_ += size_t(count) * 4; return full;
}

inline
void bounded_decoder::get_32_AND(bm::word_t* dst, unsigned count)
{
    require_array(count, 4);
    bm::decoder d(data_ + pos_); d.get_32_AND(dst, count); pos_ += size_t(count) * 4;
}

inline
void bounded_decoder::memcpy(unsigned char* dst, size_t count)
{
    require(count);
    if (dst && count) ::memcpy(dst, data_ + pos_, count);
    pos_ += count;
}

inline
void bounded_decoder::require(size_t n) const
{
    if (n > remaining()) throw exhausted();
}

inline
void bounded_decoder::require_array(unsigned n, unsigned width) const
{
    if (size_t(n) > remaining() / width) throw exhausted();
}


// --------------------------------------------
// Implementations for: streams_decoder
// --------------------------------------------

inline
streams_decoder::streams_decoder(std::istream& stream) BMNOEXCEPT
        : stream_(stream)
{
    
}

inline
bool streams_decoder::is_good() const BMNOEXCEPT
{
    return error_ == no_error && stream_.good();
}

inline
streams_decoder::error_code streams_decoder::error() const BMNOEXCEPT
{
    return error_;
}

inline
std::istream& streams_decoder::stream() BMNOEXCEPT
{
    return stream_;
}

inline
bounded_decoder& streams_decoder::get_decoder() BMNOEXCEPT
{
    return dec_;
}

inline
size_t streams_decoder::tell() const BMNOEXCEPT
{
    return offset_ + dec_.size();
}

inline
bool streams_decoder::begin(unsigned char* data, size_t capacity)
{
    if (!is_good()) return fail(io_error);
    error_ = io_error;
    origin_ = stream_.tellg();
    if (!stream_.good() || origin_ == std::streampos(-1)) return false;
    data_ = data; capacity_ = capacity; offset_ = 0; eof_ = false;
    dec_.bind(data, 0);
    error_ = no_error;
    return true;
}

inline
void streams_decoder::unbind_buffer() BMNOEXCEPT
{
    data_ = 0;
    capacity_ = 0;
    dec_.bind(0, 0);
}

inline
bool streams_decoder::fail(error_code e) BMNOEXCEPT
{
    if (error_ == no_error) error_ = e;
    return false;
}

inline
bool streams_decoder::prepare(size_t n)
{
    if (!is_good()) return false;
    if (n > capacity_) return fail();
    if (dec_.remaining() >= n || eof_) return true;
    const size_t pos = tell(), left = dec_.remaining();
    if (left) ::memmove(data_, dec_.get_pos(), left);
    offset_ = pos;
    dec_.bind(data_, left);
    error_ = io_error;
    try
    {
        stream_.read(reinterpret_cast<char*>(data_ + left),
                     std::streamsize(capacity_ - left));
    }
    catch (const std::ios_base::failure&)
    {
        // Reading the final partial window is normal, even when the caller
        // enables failbit/eofbit exceptions. Bounds checks detect truncation.
        if (!stream_.eof() || stream_.bad()) throw;
    }
    if (stream_.bad() || (stream_.fail() && !stream_.eof())) return false;
    size_t got = size_t(stream_.gcount());
    eof_ = stream_.eof();
    if (eof_) stream_.clear();
    dec_.bind(data_, left + got);
    error_ = no_error;
    return true;
}

inline
bool streams_decoder::seek(size_t pos)
{
    if (!is_good()) return false;
    const size_t extent = dec_.size() + dec_.remaining();
    if (pos >= offset_ && pos - offset_ <= extent)
    {
        dec_.set_pos(data_ + (pos - offset_));
        return true;
    }
    if (!seek_stream(pos)) return false;
    offset_ = pos; eof_ = false; dec_.bind(data_, 0);
    return true;
}

inline
bool streams_decoder::complete()
{
    if (!is_good()) return false;
    const size_t pos = tell();
    if (!seek_stream(pos)) return false;
    offset_ = pos; eof_ = false; dec_.bind(data_, 0);
    return true;
}



inline
bool streams_decoder::seek_stream(size_t pos)
{
    const std::streamoff base = std::streamoff(origin_);
    if (base < 0 || pos > size_t((std::numeric_limits<std::streamoff>::max)() - base))
        return fail(io_error);
    error_ = io_error;
    stream_.seekg(origin_ + std::streamoff(pos));
    if (!stream_.good()) return false;
    error_ = no_error;
    return true;
}


template<class Buffer>
bool streams_decoder::reserve(Buffer& buffer, size_t capacity)
{
    if (capacity <= capacity_) return is_good();
    if (!is_good()) return false;
    const size_t pos = tell(), left = dec_.remaining();
    if (left) ::memmove(data_, dec_.get_pos(), left);
    unbind_buffer(); // no dangling borrow if allocation throws
    buffer.resize(capacity, true);
    data_ = buffer.data(); capacity_ = buffer.size(); offset_ = pos;
    dec_.bind(data_, left);
    return true;
}

// --------------------------------------------
// Implementations for: streams_deserializer<BV>
// --------------------------------------------

template<class BV>
template<class IN>
bool streams_deserializer<BV>::deserialize(BV& bv, IN& input)
{
    const size_t capacity = 8 * bm::set_block_size * sizeof(bm::word_t) + 1024;
    if (!buffer_.size()) buffer_.resize(capacity, false);
    release<IN> guard = {input, bv, bv.get_new_blocks_strat()};
    try
    {
        if (!input.begin(buffer_.data(), buffer_.size())) return false;
        source<IN> view(input, buffer_);
        const size_t count = this->deserialize_from(bv, view);
        if (count == size_t(-1)) return input.fail();
        return input.complete();
    }
    catch (const bm::bounded_decoder::exhausted&)
    {
        return input.fail();
    }
    catch (...)
    {
        input.fail();
        throw;
    }
}

template<class BV>
template<class IN>
streams_deserializer<BV>::release<IN>::~release()
{
    target.set_new_blocks_strat(strategy);
    in.unbind_buffer();
}

template<class BV>
template<class IN>
streams_deserializer<BV>::source<IN>::source(IN& input, buffer_type& memory)
    : in(input), buffer(memory)
{}

template<class BV>
template<class IN>
bm::bounded_decoder& streams_deserializer<BV>::source<IN>::get_decoder()
{
    return in.get_decoder();
}

template<class BV>
template<class IN>
bool streams_deserializer<BV>::source<IN>::prepare(size_t n)
{
    return in.reserve(buffer, n) && in.prepare(n);
}

template<class BV>
template<class IN>
bool streams_deserializer<BV>::source<IN>::seek(size_t pos)
{
    return in.seek(pos);
}

template<class BV>
template<class IN>
size_t streams_deserializer<BV>::source<IN>::tell() const
{
    return in.tell();
}

// --------------------------------------------
// Implementations for: streams_sparse_vector_deserializer<SV>
// --------------------------------------------

template<class SV>
bool streams_sparse_vector_deserializer<SV>::seek(
    std::istream& stream, std::streampos origin, size_t offset)
{
    if (offset > size_t((std::numeric_limits<std::streamoff>::max)())) return false;
    stream.seekg(origin + std::streamoff(offset));
    return bool(stream);
}

template<class SV>
bool streams_sparse_vector_deserializer<SV>::deserialize(
    SV& sv, std::istream& stream, const bvector_type* mask)
{
    return read(sv, stream, mask, 0);
}

template<class SV>
bool streams_sparse_vector_deserializer<SV>::construct_deserialization_index(
    deserialization_index_type& index, std::istream& stream)
{
    SV scratch;
    return read(scratch, stream, 0, &index);
}

template<class SV>
bool streams_sparse_vector_deserializer<SV>::read(
    SV& sv, std::istream& stream, const bvector_type* mask,
    deserialization_index_type* build)
{
    struct guard_type
    {
        streams_deserializer<bvector_type>& reader;
        typename base_type::bv_ref_vector_type& refs;
        ~guard_type()
        {
            reader.unset_deserialization_index(); reader.unset_block_digest_vector();
            reader.set_ref_vectors(0); refs.reset();
        }
    } guard = {reader_, this->bv_ref_};
    reader_.unset_deserialization_index(); reader_.unset_block_digest_vector();
    reader_.set_ref_vectors(0); this->bv_ref_.reset();
    if (!stream) return false;
    const std::streampos origin = stream.tellg();
    if (origin == std::streampos(-1)) return false;
    unsigned char header[33] = {};
    stream.read(reinterpret_cast<char*>(header), 2);
    if (!stream || header[0] != 'B') return false;
    // A subsequent BLOB may switch from remapped strings to plain strings.
    if (sv.is_null_external()) sv.clear_all_preserve_null(true, 0);
    else sv.clear_all(true, 0);
    if (header[1] == 'Z')
    {
        if (build) build->reset(0);
        return true;
    }
    if (header[1] != (sv.is_compressed() ? 'C' : 'M')) return false;
    stream.read(reinterpret_cast<char*>(header)+2, 31);
    if (!stream || header[2] != (unsigned char)globals<true>::byte_order() ||
        header[3] != 0 || (header[4] != 1 && header[4] != 2)) return false;
    bm::decoder hdr(header);
    unsigned char version = 0;
    const unsigned planes = this->load_header(hdr, sv, version);
    if (!this->digest_offset_ || this->digest_offset_ < sizeof(header) ||
        this->digest_offset_ > (std::numeric_limits<size_t>::max)()) return false;
    if (!seek(stream, origin, size_t(this->digest_offset_))) return false;
    this->plane_digest_bv_.clear();
    streams_decoder digest_input(stream);
    if (!reader_.deserialize(this->plane_digest_bv_, digest_input)) return false;
    typename bvector_type::size_type last_plane = 0;
    if (this->plane_digest_bv_.find_reverse(last_plane) && last_plane >= planes) return false;
    const unsigned count = unsigned(this->plane_digest_bv_.count());
    this->off_vect_.resize(planes);
    for (unsigned i = 0; i < planes; ++i) this->off_vect_[i] = 0;
    metadata_.resize(size_t(planes)*8 + 128, false);
    streams_decoder table(stream);
    if (!table.begin(metadata_.data(), metadata_.size())) return false;
    struct table_guard_type
    {
        streams_decoder& input;
        ~table_guard_type() { input.unbind_buffer(); }
    } table_guard = {table};
    try
    {
        if (!table.prepare(metadata_.size())) return false;
        bounded_decoder& dec = table.get_decoder();
        const unsigned kind = dec.get_8();
        if (kind == '6')
        {
            for (unsigned i = 0; i < planes; ++i)
                if (this->plane_digest_bv_.test(i))
                {
                    bm::id64_t offset = dec.get_64();
                    if (offset > (std::numeric_limits<size_t>::max)()) return false;
                    this->off_vect_[i] = size_t(offset);
                }
        }
        else if (kind == '3' && count >= 4)
        {
            this->off32_vect_.resize(count);
            unsigned lo = dec.get_32(), hi = dec.get_32();
            if (hi <= lo || hi-lo < count-1) return false;
            this->off32_vect_[0] = lo; this->off32_vect_[count-1] = hi;
            bm::bit_in<bounded_decoder> bits(dec);
            bits.bic_decode_u32_cm(this->off32_vect_.data()+1, count-2, lo, hi);
            unsigned j = 0;
            for (unsigned i = 0; i < planes; ++i)
                if (this->plane_digest_bv_.test(i)) this->off_vect_[i] = this->off32_vect_[j++];
        }
        else return false;
        if (!table.complete()) return false;
    }
    catch (const bounded_decoder::exhausted&) { return false; }
    const std::streampos end = stream.tellg();
    if (end == std::streampos(-1)) return false;
    table.unbind_buffer(); // metadata buffer may be reused for remapping below
    size_t previous = 0;
    for (unsigned i = 0; i < planes; ++i)
    {
        if (this->plane_digest_bv_.test(i) && !this->off_vect_[i]) return false;
        if (size_t offset = this->off_vect_[i])
        {
            if (offset < sizeof(header) || offset <= previous || offset >= this->digest_offset_)
                return false;
            previous = offset;
        }
    }
    if (build)
    {
        build->reset(planes);
        build->set_plane_offsets(this->off_vect_, size_t(this->digest_offset_));
    }
    this->resize_stream_target(sv, typename SV::size_type(this->sv_size_));
    sv.get_bmatrix().allocate_rows(planes);
    reader_.set_ref_vectors(this->bv_ref_ptr_ ? this->bv_ref_ptr_ : &this->bv_ref_);
    const bvector_type* selection = mask;
    bvector_type digest;
    size_t remap_offset = sizeof(header); // valid even when no planes are present
    // RSC masks use logical addresses until the complete NULL plane is loaded.
    for (int row = int(planes)-1; row >= 0; --row)
    {
        const unsigned i = unsigned(row);
        const bool null_plane = sv.is_nullable() && i == planes-1;
        const size_t offset = this->off_vect_[i];
        if (!offset)
        {
            if (null_plane && !sv.is_null_external() && !build) return false;
            if (null_plane && SV::is_rsc_support::value && mask && !build)
            {
                const bvector_type* nulls = sv.get_null_bvector();
                this->not_null_mask_bv_.bit_and(*nulls, *mask, bvector_type::opt_compress);
                this->rsc_mask_bv_.clear();
                this->rsc_compressor_.compress(this->rsc_mask_bv_, *nulls, this->not_null_mask_bv_);
                selection = &this->rsc_mask_bv_;
            }
            continue;
        }
        // A reused dynamic string target can retain more rows than this BLOB.
        // Its NULL plane remains at the target's current last row.
        bvector_type* bv = null_plane
            ? sv.get_bmatrix().get_row(unsigned(sv.get_bmatrix().get_null_idx()))
            : sv.get_create_slice(i);
        if (null_plane && sv.is_null_external())
        {
            this->null_decode_scratch_bv_.clear();
            bv = &this->null_decode_scratch_bv_;
        }
        if (!this->bv_ref_ptr_) this->bv_ref_.add(bv, i);
        reader_.unset_deserialization_index(); reader_.unset_block_digest_vector();
        if (build) reader_.set_deserialization_index_construct(build->construct_row(i));
        else if (selection && !(null_plane && SV::is_rsc_support::value) &&
                 this->deserialization_index_ && this->deserialization_index_use_)
        {
            const auto* index_row = this->deserialization_index_->get_row(i);
            if (index_row)
            {
                selection->build_block_digest(digest);
                reader_.set_deserialization_index_use(index_row);
                reader_.set_block_digest_vector_use(&digest);
            }
        }
        if (!seek(stream, origin, offset)) return false;
        streams_decoder input(stream);
        if (!reader_.deserialize(*bv, input)) return false;
        if constexpr (SV::is_remap_support::value)
            if (offset == previous)
            {
                const std::streampos pos = stream.tellg();
                if (pos == std::streampos(-1) || pos < origin) return false;
                remap_offset = size_t(pos-origin);
            }
        if (build) bv->clear(true);
        else if (null_plane && SV::is_rsc_support::value && mask)
        {
            const bvector_type* nulls = sv.get_null_bvector();
            this->not_null_mask_bv_.bit_and(*nulls, *mask, bvector_type::opt_compress);
            this->rsc_mask_bv_.clear();
            this->rsc_compressor_.compress(this->rsc_mask_bv_, *nulls, this->not_null_mask_bv_);
            selection = &this->rsc_mask_bv_;
        }
    }
    if constexpr (SV::is_remap_support::value)
    {
        if (remap_offset >= this->digest_offset_) return false;
        if (build) build->set_remap_offset(remap_offset);
        if (!seek(stream, origin, remap_offset)) return false;
        streams_decoder remap(stream);
        if (!remap.begin(metadata_.data(), metadata_.size())) return false;
        table_guard_type remap_guard = {remap};
        try
        {
            if (!this->load_remap_stream(sv, remap, metadata_) || !remap.complete()) return false;
        }
        catch (const bounded_decoder::exhausted&) { return false; }
        if (stream.tellg() != origin + std::streamoff(this->digest_offset_)) return false;
    }
    if (!build)
    {
        // Delay AND until XOR dependencies no longer need the original blocks.
        for (unsigned i = 0; i < planes; ++i)
        {
            if (i == planes-1 && sv.is_nullable() &&
                (SV::is_rsc_support::value || sv.is_null_external())) continue;
            bvector_type* bv = sv.get_bmatrix().get_row(i);
            if (!bv) continue;
            if (selection) bv->bit_and(*selection, bvector_type::opt_compress);
            else if (!SV::is_rsc_support::value && this->sv_size_)
                bv->keep_range(0, typename SV::size_type(this->sv_size_-1));
            if (this->is_final_ == bm::finalization::READONLY)
                bv->optimize_freeze(this->temp_block_, bvector_type::opt_compress);
        }
        if (selection && sv.is_nullable() && !SV::is_rsc_support::value &&
            !sv.is_null_external() && sv.get_bmatrix().get_null_idx() >= planes)
        {
            bvector_type* nulls = sv.get_bmatrix().get_row(unsigned(sv.get_bmatrix().get_null_idx()));
            if (nulls) nulls->bit_and(*selection, bvector_type::opt_compress);
        }
        if (sv.max_vector_size == 1 && !sv.is_null_external())
            if (sv.get_bmatrix().get_row(sv.sv_value_slices)) sv.mark_null_idx(sv.sv_value_slices);
        sv.sync(true, true);
    }
    stream.seekg(end);
    return bool(stream);
}

// --------------------------------------------
// Implementations for: streams_sparse_vector_float_deserializer
// --------------------------------------------

template<class SV>
bool streams_sparse_vector_float_deserializer<SV>::header(
    std::istream& stream, std::streampos (&positions)[4])
{
    char signature[3]; size_t sizes[3];
    if (!stream.read(signature, 3) || std::memcmp(signature, "bf0", 3) ||
        !stream.read(reinterpret_cast<char*>(sizes), sizeof(sizes))) return false;
    positions[0] = stream.tellg();
    if (positions[0] == std::streampos(-1)) return false;
    stream.seekg(0, std::ios::end);
    const std::streampos end = stream.tellg();
    if (!stream || end < positions[0]) return false;
    for (unsigned i = 0; i < 3; ++i)
    {
        if (!sizes[i] || sizes[i] > uintmax_t(end - positions[i])) return false;
        positions[i+1] = positions[i] + std::streamoff(sizes[i]);
    }
    stream.seekg(positions[0]);
    return bool(stream);
}

template<class SV>
void streams_sparse_vector_float_deserializer<SV>::set_deserialization_index(
    deserialization_index_type* index)
{
    exponents_.set_deserialization_index(index ? &index->exponent_index_ : 0);
    mantissas_.set_deserialization_index(index ? &index->mantissa_index_ : 0);
}

template<class SV>
void streams_sparse_vector_float_deserializer<SV>::set_deserialization_index_use(bool enable)
{
    exponents_.set_deserialization_index_use(enable);
    mantissas_.set_deserialization_index_use(enable);
}

template<class SV>
void streams_sparse_vector_float_deserializer<SV>::set_xor_ref(bv_ref_vector_type* refs)
{
    exponents_.set_xor_ref(refs);
    mantissas_.set_xor_ref(refs);
}

template<class SV>
bool streams_sparse_vector_float_deserializer<SV>::deserialize(
    SV& sv, std::istream& stream, const bvector_type* mask)
{
    std::streampos positions[4];
    if (!header(stream, positions)) return false;
    sv.clear();
    streams_decoder input(stream);
    if (!signs_.deserialize(sv.signs_, input) || stream.tellg() != positions[1]) return false;
    if (mask) sv.signs_ &= *mask;
    if (!exponents_.deserialize(sv.exponents_, stream, mask) ||
        stream.tellg() != positions[2]) return false;
    if constexpr (is_rsc_sparse_vector<typename SV::sparse_vector_u>::value)
        sv.exponents_.sync(false, false);
    sv.attach_mantissa_null_plane_();
    if (!mantissas_.deserialize(sv.mantissas_, stream, mask) ||
        stream.tellg() != positions[3]) return false;
    return bool(stream);
}

template<class SV>
bool streams_sparse_vector_float_deserializer<SV>::construct_deserialization_index(
    deserialization_index_type& index, std::istream& stream)
{
    std::streampos positions[4];
    if (!header(stream, positions)) return false;
    index.reset();
    stream.seekg(positions[1]);
    if (!exponents_.construct_deserialization_index(index.exponent_index_, stream) ||
        stream.tellg() != positions[2]) return false;
    if (!mantissas_.construct_deserialization_index(index.mantissa_index_, stream) ||
        stream.tellg() != positions[3]) return false;
    return bool(stream);
}

// --------------------------------------------
// Implementations for: sv_di_detail::stream_output
// --------------------------------------------
inline sv_di_detail::stream_output::stream_output(std::ostream& s)
    : stream(s), origin(s.tellp())
{
    if (origin == std::streampos(-1)) stream.setstate(std::ios::failbit);
}
inline bool sv_di_detail::stream_output::write(const unsigned char* p, size_t n)
{
    require(n <= size_t((std::numeric_limits<std::streamsize>::max)()));
    const size_t next = add(pos, n);
    if (!stream.write(reinterpret_cast<const char*>(p), std::streamsize(n))) return false;
    pos = next; return true;
}
inline bool sv_di_detail::stream_output::patch(size_t at, const unsigned char* p, size_t n)
{
    require(at <= pos && n <= pos - at && pos <= size_t((std::numeric_limits<std::streamoff>::max)()));
    stream.seekp(origin + std::streamoff(at));
    stream.write(reinterpret_cast<const char*>(p), std::streamsize(n));
    stream.seekp(origin + std::streamoff(pos)); return bool(stream);
}
// --------------------------------------------
// Implementations for: sv_di_detail::stream_input
// --------------------------------------------
inline sv_di_detail::stream_input::stream_input(std::istream& s) : stream(s)
{
    const std::streampos start = stream.tellg();
    if (start == std::streampos(-1)) { stream.setstate(std::ios::failbit); return; }
    stream.seekg(0, std::ios::end);
    const std::streampos end = stream.tellg();
    if (!stream || end < start) { stream.setstate(std::ios::failbit); return; }
    require(bm::id64_t(end - start) <= bm::id64_t(size_t(-1)));
    remaining = size_t(end - start); stream.seekg(start);
}
inline bool sv_di_detail::stream_input::read(unsigned char* p, size_t n)
{
    require(n <= remaining && n <= size_t((std::numeric_limits<std::streamsize>::max)()));
    if (!stream.read(reinterpret_cast<char*>(p), std::streamsize(n))) return false;
    remaining -= n; pos += n; return true;
}

// --------------------------------------------
// Stream implementations for: sparse_vector_deserialization_index_serializer
// --------------------------------------------
template<class BV>
size_t sparse_vector_deserialization_index_serializer<BV>::serialize(const index_type& index, std::ostream& stream)
{
    sv_di_detail::stream_output out(stream);
    if (!stream) return 0;
    return serialize_to(index, out);
}
// --------------------------------------------
// Stream implementations for: sparse_vector_deserialization_index_deserializer
// --------------------------------------------
template<class BV>
size_t sparse_vector_deserialization_index_deserializer<BV>::deserialize(index_type& index, std::istream& stream)
{
    sv_di_detail::stream_input in(stream);
    if (!stream) return 0;
    return deserialize_from(index, in);
}
// --------------------------------------------
// Stream implementations for: sparse_vector_float_deserialization_index_serializer
// --------------------------------------------
template<class SV>
size_t sparse_vector_float_deserialization_index_serializer<SV>::serialize(const index_type& index, std::ostream& stream)
{
    sv_di_detail::stream_output out(stream);
    if (!stream) return 0;
    return serialize_to(index, out);
}
// --------------------------------------------
// Stream implementations for: sparse_vector_float_deserialization_index_deserializer
// --------------------------------------------
template<class SV>
size_t sparse_vector_float_deserialization_index_deserializer<SV>::deserialize(index_type& index, std::istream& stream)
{
    sv_di_detail::stream_input in(stream);
    if (!stream) return 0;
    return deserialize_from(index, in);
}
} // namespace bm
#endif
