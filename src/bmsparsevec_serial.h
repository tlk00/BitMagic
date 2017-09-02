#ifndef BMSPARSEVEC_SERIAL__H__INCLUDED__
#define BMSPARSEVEC_SERIAL__H__INCLUDED__
/*
Copyright(c) 2017 Anatoliy Kuznetsov(anatoliy_kuznetsov at yahoo.com)

Permission is hereby granted, free of charge, to any person 
obtaining a copy of this software and associated documentation 
files (the "Software"), to deal in the Software without restriction, 
including without limitation the rights to use, copy, modify, merge, 
publish, distribute, sublicense, and/or sell copies of the Software, 
and to permit persons to whom the Software is furnished to do so, 
subject to the following conditions:

The above copyright notice and this permission notice shall be included 
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, 
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, 
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR 
OTHER DEALINGS IN THE SOFTWARE.

For more information please visit:  http://bmagic.sourceforge.net

*/

#include "bmdef.h"
#include "bmsparsevec.h"

namespace bm
{

/*!
    \brief layout class for serialization buffer structure
    
    Class keeps a memory block sized for the target sparse vector BLOB.
    This class also provides acess to bit-plane memory, so it becomes possible
    to use parallel storage methods to save bit-plains into
    different storage shards.
 
*/
template<class SV>
struct sparse_vector_serial_layout
{
    typedef typename SV::value_type value_type;

    sparse_vector_serial_layout()
    : buffer_(0), capacity_(0), serialized_size_(0)
    {
    }
    
    ~sparse_vector_serial_layout()
    {
        if (buffer_)
            ::free(buffer_);
    }
    
    /*!
        \brief resize capacity
        \param capacity - new capacity
        \return new buffer or 0 if failed
    */
    unsigned char* reserve(size_t capacity)
    {
        if (capacity == 0)
        {
            freemem();
            return 0;
        }
        if (capacity <= capacity_)  // buffer reduction - avoid
        {
            serialized_size_ = 0;
            return buffer_;
        }
        // buffer growth
        freemem();
        buffer_ = (unsigned char*) ::malloc(capacity);
        if (buffer_)
        {
            capacity_ = capacity;
        }
        return buffer_;
    }
    
    /// return current serialized size
    size_t  size() const { return serialized_size_; }
    
    /// Set new serialized size
    void resize(size_t ssize) { serialized_size_ = ssize; }
    
    /// return serialization buffer capacity
    size_t  capacity() const { return capacity_; }
    
    /// free memory
    void freemem()
    {
        if (buffer_)
        {
            ::free(buffer_);
            buffer_ = 0; capacity_ = serialized_size_ = 0;
        }
    }
    
    /// Set plain output pointer and size
    void set_plain(unsigned i, unsigned char* ptr, unsigned buf_size)
    {
        plain_ptrs_[i] = ptr;
        plane_size_[i] = buf_size;
    }
    
    /// Get plain pointer
    const unsigned char* get_plain(unsigned i) const
    {
        return plain_ptrs_[i];
    }
    
    /// Return serializatio buffer pointer
    const unsigned char* buf() const { return buffer_; }
    
private:
    sparse_vector_serial_layout(const sparse_vector_serial_layout&);
    void operator=(const sparse_vector_serial_layout&);
protected:
    unsigned char* buffer_;                       ///< serialization buffer
    size_t         capacity_;                      ///< buffer capacity
    size_t         serialized_size_;               ///< serialized size

    unsigned char* plain_ptrs_[sizeof(value_type)*8]; ///< pointers on serialized bit-palins
    unsigned plane_size_[sizeof(value_type)*8];       ///< serialized plain size
};

// -------------------------------------------------------------------------


/*!
    \brief Serialize sparse vector into a buffer(s) structure
 
 Serialization format:
 <pre>

 | HEADER | BITVECTRORS |

 Header structure:
   BYTE+BYTE: Magic-signature 'BM'
   BYTE : Byte order ( 0 - Big Endian, 1 - Little Endian)
   BYTE : Number of Bit-vector plains (total)
   INT64: Vector size
   INT64: Offset of plain 0 from the header start (value 0 means plain is empty)
   INT64: Offset of plain 1 from
   ...
   INT32: reserved

 </pre>
 
    \param sv         - sparse vector to serialize
    \param sv_layout  - buffer structure to keep the result
    \param bv_serialization_flags - bit-vector serialization flags
    as defined in bm::serialization_flags    
    
    \return "0" - success, "-1" memory allocation error
    
    @sa serialization_flags
*/
template<class SV>
int sparse_vector_serialize(
                const SV&                        sv,
                sparse_vector_serial_layout<SV>& sv_layout,
                unsigned                         bv_serialization_flags = 0)
{
    typename SV::statistics sv_stat;
    sv.calc_stat(&sv_stat);
    
    unsigned char* buf = sv_layout.reserve(sv_stat.max_serialize_mem);
    if (!buf) // memory allocation error
    {
        return -1;
    }
    bm::encoder enc(buf, (unsigned)sv_layout.capacity());
    unsigned plains = sv.plains();

    
    // calculate header size in bytes
    unsigned h_size = 1 + 1 + 1 + 1 + 8 + (8 * plains) + 4;

    // ptr where bit-vectors start
    unsigned char* buf_ptr = buf + h_size;

    unsigned i;
    for (i = 0; i < plains; ++i)
    {
        const typename SV::bvector_type_ptr bv = sv.plain(i);
        if (!bv)  // empty plain
        {
            sv_layout.set_plain(i, 0, 0);
            continue;
        }
        unsigned buf_size =
            bm::serialize(*bv, buf_ptr, 0, bv_serialization_flags);
        sv_layout.set_plain(i, buf_ptr, buf_size);
        buf_ptr += buf_size;
        
    } // for i
    
    sv_layout.resize(buf_ptr - buf);
    
    
    // save header
    ByteOrder bo = globals<true>::byte_order();
    
    enc.put_8('B');
    enc.put_8('M');
    enc.put_8((unsigned char)bo);
    enc.put_8((unsigned char)plains);
    enc.put_64(sv.size());
    
    for (i = 0; i < plains; ++i)
    {
        const unsigned char* p = sv_layout.get_plain(i);
        if (!p)
        {
            enc.put_64(0);
            continue;
        }
        size_t offset = p - buf;
        enc.put_64(offset);
    }

    return 0;
}

// -------------------------------------------------------------------------


template<class SV>
int sparse_vector_deserialize(SV& sv,
                              const unsigned char* buf,
                              bm::word_t* temp_block=0)
{
    typedef typename SV::bvector_type   bvector_type;

    // TODO: implement correct processing of byte-order corect deserialization
//    ByteOrder bo_current = globals<true>::byte_order();

    bm::decoder dec(buf);
    unsigned char h1 = dec.get_8();
    unsigned char h2 = dec.get_8();

    BM_ASSERT(h1 == 'B' && h2 == 'M');
    if (h1 != 'B' && h2 != 'M')  // no magic header? issue...
    {
        return -1;
    }
    
    //unsigned char bv_bo =
        dec.get_8();
    unsigned plains = dec.get_8();
    
    if (!plains || plains > sv.plains())
    {
        return -2; // incorrect number of plains for the target svector
    }
    
    sv.clear();
    bm::id64_t sv_size = dec.get_64();
    if (sv_size == 0)
    {
        return 0;  // empty vector
    }
    sv.resize((unsigned)sv_size);
    
    unsigned i = 0;
    for (i = 0; i < plains; ++i)
    {
        size_t offset = (size_t) dec.get_64();
        if (offset == 0) // null vector
        {
            continue;
        }
        const unsigned char* bv_buf_ptr = buf + offset;
        bvector_type*  bv = sv.get_plain(i);
        BM_ASSERT(bv);
        
        bm::deserialize(*bv, bv_buf_ptr, temp_block);
        if (!temp_block)
        {
            typename bvector_type::blocks_manager_type& bv_bm =
                                                bv->get_blocks_manager();
            temp_block = bv_bm.check_allocate_tempblock();
        }
    } // for i
    return 0;
}

// -------------------------------------------------------------------------


    
} // namespace bm

#include "bmundef.h"

#endif
