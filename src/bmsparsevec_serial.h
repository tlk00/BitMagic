#ifndef BMSPARSEVEC_SERIAL__H__INCLUDED__
#define BMSPARSEVEC_SERIAL__H__INCLUDED__
/*
Copyright(c) 2002-2017 Anatoliy Kuznetsov(anatoliy_kuznetsov at yahoo.com)

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

/*! \file bmsparsevec_serial.h
    \brief Serialization for sparse_vector<>
*/

#include <vector>

#include "bmsparsevec.h"
#include "bmserial.h"
#include "bmdef.h"

namespace bm
{

/** \defgroup svserial Sparse vector serialization
    Sparse vector serialization
    \ingroup svector
 */


/*!
    \brief layout class for serialization buffer structure
    
    Class keeps a memory block sized for the target sparse vector BLOB.
    This class also provides acess to bit-plane memory, so it becomes possible
    to use parallel storage methods to save bit-plains into
    different storage shards.
    
    \ingroup svserial
 
*/
template<class SV>
struct sparse_vector_serial_layout
{
    typedef typename SV::value_type   value_type;
    typedef typename SV::bvector_type bvector_type;
    typedef typename serializer<bvector_type>::buffer buffer_type;

    sparse_vector_serial_layout()
    {
    }
    
    ~sparse_vector_serial_layout()
    {
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
        buf_.reinit(capacity);
        return buf_.data();
    }
    
    /// return current serialized size
    size_t  size() const { return buf_.size();  }
    
    /// Set new serialized size
    void resize(size_t ssize) { buf_.resize(ssize);  }
    
    /// return serialization buffer capacity
    size_t  capacity() const { return buf_.capacity(); }
    
    /// free memory
    void freemem()
    {
        buf_.release();
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
    
    /// Return serialization buffer pointer
    const unsigned char* buf() const { return buf_.buf(); /*return buffer_;*/ }
    
private:
    sparse_vector_serial_layout(const sparse_vector_serial_layout&);
    void operator=(const sparse_vector_serial_layout&);
protected:
    buffer_type    buf_;                       ///< serialization buffer
    unsigned char* plain_ptrs_[SV::sv_plains]; ///< pointers on serialized bit-palins
    unsigned plane_size_[SV::sv_plains];       ///< serialized plain size
};

// -------------------------------------------------------------------------


/*!
    \brief Serialize sparse vector into a memory buffer(s) structure
 
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
    \param temp_block - temporary buffer 
                        (allocate with BM_DECLARE_TEMP_BLOCK(x) for speed)
    \param bv_serialization_flags - bit-vector serialization flags
    as defined in bm::serialization_flags    
    
    \ingroup svserial
    
    @sa serialization_flags
*/
template<class SV>
void sparse_vector_serialize(
                const SV&                        sv,
                sparse_vector_serial_layout<SV>& sv_layout,
                bm::word_t*                      temp_block = 0)
{
    typename SV::statistics sv_stat;
    sv.calc_stat(&sv_stat);
    
    unsigned char* buf = sv_layout.reserve(sv_stat.max_serialize_mem);
    bm::encoder enc(buf, (unsigned)sv_layout.capacity());
    unsigned plains = sv.stored_plains();

    // calculate header size in bytes
    unsigned h_size = 1 + 1 + 1 + 1 + 8 + (8 * plains) + 4;

    // ptr where bit-plains start
    unsigned char* buf_ptr = buf + h_size;

    bm::serializer<typename SV::bvector_type > bvs(temp_block);
    bvs.gap_length_serialization(false);
    bvs.set_compression_level(4);
    
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
            bvs.serialize(*bv, buf_ptr, sv_stat.max_serialize_mem);
        
        sv_layout.set_plain(i, buf_ptr, buf_size);
        buf_ptr += buf_size;
        sv_stat.max_serialize_mem -= buf_size;
        
    } // for i
    
    sv_layout.resize(size_t(buf_ptr - buf));
    
    
    // save the header
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
        size_t offset = size_t(p - buf);
        enc.put_64(offset);
    }
}

// -------------------------------------------------------------------------

/*!
    \brief Deserialize svector<>
    \param sv         - target sparse vector
    \param buf        - source memory buffer
    \param temp_block - temporary block buffer to avoid re-allocations
 
    \return error non-zero codes means failure
 
    \ingroup svector
*/
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
    
    if (!plains || plains > sv.stored_plains())
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

/**
    Seriaizer for compressed collections
*/
template<class CBC>
class compressed_collection_serializer
{
public:
    typedef CBC                                  compressed_collection_type;
    typedef typename CBC::bvector_type           bvector_type;
    typedef typename CBC::buffer_type            buffer_type;
    typedef typename CBC::statistics             statistics_type;
    typedef typename CBC::address_resolver_type  address_resolver_type;
    
public:
    void serialize(const CBC&    buffer_coll,
                   buffer_type&  buf,
                   bm::word_t*   temp_block = 0);
};

/**
    Deseriaizer for compressed collections
*/
template<class CBC>
class compressed_collection_deserializer
{
public:
    typedef CBC                                  compressed_collection_type;
    typedef typename CBC::bvector_type           bvector_type;
    typedef typename CBC::buffer_type            buffer_type;
    typedef typename CBC::statistics             statistics_type;
    typedef typename CBC::address_resolver_type  address_resolver_type;
    typedef typename CBC::container_type         container_type;

public:
    int deserialize(CBC&                 buffer_coll,
                    const unsigned char* buf,
                    bm::word_t*          temp_block=0);
};


// -------------------------------------------------------------------------

/**
    \brief Serialize compressed collection into memory buffer

Serialization format:


<pre>
 | MAGIC_HEADER | ADDRESS_BITVECTROR | LIST_OF_BUFFER_SIZES | BUFFER(s)
 
   MAGIC_HEADER:
   BYTE+BYTE: Magic-signature 'BM'
   BYTE : Byte order ( 0 - Big Endian, 1 - Little Endian)
 
   ADDRESS_BITVECTROR:
   INT64: address bit-vector size
   <memblock>: serialized address bit-vector
 
   LIST_OF_BUFFER_SIZES:
   INT64 - buffer sizes count
   INT32 - buffer size 0
   INT32 - buffer size 1
   ...
 
   BUFFERS:
   <memblock>: block0
   <memblock>: block1
   ...
 
</pre>
*/

template<class CBC>
void compressed_collection_serializer<CBC>::serialize(const CBC&    buffer_coll,
                                                      buffer_type&  buf,
                                                      bm::word_t*   temp_block)
{
    statistics_type st;
    buffer_coll.calc_stat(&st);
    
    buf.resize(st.max_serialize_mem);
    
    // ptr where bit-plains start
    unsigned char* buf_ptr = buf.data();
    
    bm::encoder enc(buf.data(), buf.capacity());
    ByteOrder bo = globals<true>::byte_order();
    enc.put_8('B');
    enc.put_8('C');
    enc.put_8((unsigned char)bo);
    
    unsigned char* mbuf1 = enc.get_pos(); // bookmark position
    enc.put_64(0);  // address vector size (reservation)

    buf_ptr = enc.get_pos();

    const address_resolver_type& addr_res = buffer_coll.resolver();
    const bvector_type& bv = addr_res.get_bvector();
    {
        bm::serializer<bvector_type > bvs(temp_block);
        bvs.gap_length_serialization(false);
        bvs.set_compression_level(4);

        size_t addr_bv_size = bvs.serialize(bv, buf_ptr, buf.size());
        buf_ptr += addr_bv_size;

        enc.set_pos(mbuf1); // rewind to bookmark
        enc.put_64(addr_bv_size); // save the address vector size
    }
    enc.set_pos(buf_ptr); // restore stream position
    size_t coll_size = buffer_coll.size();
    
    enc.put_64(coll_size);
    
    // pass 1 (save buffer sizes)
    {
        for (unsigned i = 0; i < buffer_coll.size(); ++i)
        {
            const buffer_type& cbuf = buffer_coll.get(i);
            unsigned sz = (unsigned)cbuf.size();
            enc.put_32(sz);
        } // for i
    }
    // pass 2 (save buffers)
    {
        for (unsigned i = 0; i < buffer_coll.size(); ++i)
        {
            const buffer_type& cbuf = buffer_coll.get(i);
            unsigned sz = (unsigned)cbuf.size();
            enc.memcpy(cbuf.buf(), sz);
        } // for i
    }
    buf.resize(enc.size());
    
}

// -------------------------------------------------------------------------
template<class CBC>
int compressed_collection_deserializer<CBC>::deserialize(
                                CBC&                 buffer_coll,
                                const unsigned char* buf,
                                bm::word_t*          temp_block)
{
    // TODO: implement correct processing of byte-order corect deserialization
    //    ByteOrder bo_current = globals<true>::byte_order();
    
    bm::decoder dec(buf);
    unsigned char h1 = dec.get_8();
    unsigned char h2 = dec.get_8();

    BM_ASSERT(h1 == 'B' && h2 == 'C');
    if (h1 != 'B' && h2 != 'C')  // no magic header? issue...
    {
        return -1;
    }
    //unsigned char bv_bo =
        dec.get_8();

    // -----------------------------------------
    // restore address resolver
    //
    bm::id64_t addr_bv_size = dec.get_64();
    
    const unsigned char* bv_buf_ptr = dec.get_pos();
    
    address_resolver_type& addr_res = buffer_coll.resolver();
    bvector_type& bv = addr_res.get_bvector();
    bv.clear();
    
    bm::deserialize(bv, bv_buf_ptr, temp_block);
    addr_res.sync();
    
    unsigned addr_cnt = bv.count();
    
    dec.seek((int)addr_bv_size);
    
    // -----------------------------------------
    // read buffer sizes
    //
    bm::id64_t coll_size = dec.get_64();
    if (coll_size != addr_cnt)
    {
        return -2; // buffer size collection does not match address vector
    }
    
    std::vector<unsigned> buf_size_vec(coll_size);
    {
        for (unsigned i = 0; i < coll_size; ++i)
        {
            unsigned sz = dec.get_32();
            buf_size_vec[i] = sz;
        } // for i
    }

    {
        container_type& buf_vect = buffer_coll.container();
        buf_vect.resize(coll_size);
        for (unsigned i = 0; i < coll_size; ++i)
        {
            unsigned sz = buf_size_vec[i];
            buffer_type& b = buf_vect.at(i);
            b.resize(sz);
            dec.memcpy(b.data(), sz);
        } // for i
    }
    
    buffer_coll.sync();
    
    return 0;
}


} // namespace bm

#include "bmundef.h"

#endif
