#ifndef BMBUFFER__H__INCLUDED__
#define BMBUFFER__H__INCLUDED__
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

#include "bmdef.h"

namespace bm
{

/**
    Byte buffer pointer
    @internal
*/
class byte_buffer_ptr
{
public:
    byte_buffer_ptr()
        : byte_buf_(0), size_(0)
    {}
    
    /// construct byte buffer pointer
    ///
    byte_buffer_ptr(unsigned char* in_buf, size_t in_size)
        : byte_buf_(in_buf), size_(in_size)
    {}
    
    /// Set buffer pointer
    void set_buf(unsigned char* in_buf, size_t in_size)
    {
        byte_buf_ = in_buf; size_= in_size;
    }

    /// Get buffer size
    size_t size() const { return size_; }
    
    /// Get read access to buffer memory
    const unsigned char* buf() const { return byte_buf_; }

    /// Get write access to buffer memory
    unsigned char* data() { return byte_buf_; }

    bool operator==(const byte_buffer_ptr& lhs) const
    {
        if (this == &lhs)
            return true;
        if (size_ != lhs.size_)
            return false;
        int cmp = ::memcmp(byte_buf_, lhs.byte_buf_, size_);
        return (cmp  == 0);
    }

protected:
    unsigned char* byte_buf_;     ///< byte buffer pointer to hold data
    size_t         size_;         ///< current buffer size
};

/**
    Byte buffer template, extention of byte_buffer_ptr memory management
    \internal
*/
template<typename BVAlloc>
class byte_buffer : public byte_buffer_ptr
{
public:
    typedef BVAlloc                                          bv_allocator_type;
    typedef typename bv_allocator_type::block_allocator_type allocator_type;
    
public:
    byte_buffer() : capacity_(0), alloc_factor_(0)
    {}
    
    byte_buffer(size_t in_capacity)
        : capacity_(0), alloc_factor_(0)
    {
        allocate(in_capacity);
    }
    
    byte_buffer(const byte_buffer& lhs)
    {
        byte_buf_ = 0;
        size_ = capacity_ = alloc_factor_ = 0;
        if (lhs.byte_buf_)
        {
            copy_from(lhs.byte_buf_, lhs.size_);
        }
    }
    
#ifndef BM_NO_CXX11
    /// Move constructor
    byte_buffer(byte_buffer&& in_buf) BMNOEXEPT
    {
        byte_buf_ = in_buf.byte_buf_;
        in_buf.byte_buf_ = 0;
        this->size_ = in_buf.size_;
        capacity_ = in_buf.capacity_;
        in_buf.size_ = in_buf.capacity_ = 0;
        alloc_factor_ = in_buf.alloc_factor_;
    }
    
    /// Move assignment operator
    byte_buffer& operator=(byte_buffer&& lhs) BMNOEXEPT
    {
        move_from(lhs);
        return *this;
    }
#endif

    byte_buffer& operator=(const byte_buffer& lhs)
    {
        if (this == &lhs)
            return *this;

        copy_from(lhs.buf(), lhs.size());
        return *this;
    }
    
    ~byte_buffer()
    {
        free_buffer();
    }
    
    /// swap content with another buffer
    void swap(byte_buffer& other) BMNOEXEPT
    {
        if (this == &other)
            return;
        unsigned char* btmp = byte_buf_;
        byte_buf_ = other.byte_buf_;
        other.byte_buf_ = btmp;
        
        bm::xor_swap(this->size_, other.size_);
        bm::xor_swap(capacity_, other.capacity_);
        bm::xor_swap(alloc_factor_, other.alloc_factor_);
    }
    
    /// take/move content from another buffer
    void move_from(byte_buffer& other) BMNOEXEPT
    {
        if (this == &other)
            return;
        free_buffer();
        this->byte_buf_ = other.byte_buf_; other.byte_buf_ = 0;
        this->size_ = other.size_;
        this->capacity_ = other.capacity_;
        this->alloc_factor_ = other.alloc_factor_;
    }


    /// Free underlying memory
    void release()
    {
        free_buffer();
        this->size_ = capacity_ = 0;
    }

    /// copy data from an external buffer
    ///
    void copy_from(const unsigned char* in_buf, size_t in_size)
    {
        if (in_size)
        {
            allocate(in_size);
            ::memcpy(byte_buf_, in_buf, in_size);
        }
        this->size_ = in_size;
    }
    
    
    /// Get buffer capacity
    size_t capacity() const { return capacity_; }

    /// adjust current size (buffer content preserved)
    void resize(size_t new_size)
    {
        if (new_size <= capacity_)
        {
            this->size_ = new_size;
            return;
        }
        byte_buffer tmp_buffer(new_size); // temp with new capacity
        tmp_buffer = *this;
        this->swap(tmp_buffer);
        
        this->size_ = new_size;
    }
    
    /// reserve new capacity (buffer content preserved)
    void reserve(size_t new_capacity)
    {
        if (new_capacity <= capacity_)
            return;
        
        byte_buffer tmp_buffer(new_capacity);
        tmp_buffer = *this;
        this->swap(tmp_buffer);
    }
    
    /// reserve new capacity (buffer content NOT preserved, size set to 0)
    void reinit(size_t new_capacity)
    {
        allocate(new_capacity);
        this->size_ = 0;
    }
    
    /// reserve new capacity (buffer content NOT preserved, size set to 0)
    /// @sa reinit
    void reallocate(size_t new_capacity)
    {
        reinit(new_capacity);
    }

    /// try to shrink the capacity to more optimal size
    void optimize()
    {
        if (!byte_buf_)
            return;
        size_t words = compute_words(size_);
        if (words < alloc_factor_) // possible to shrink
        {
            byte_buffer tmp_buffer(*this);
            this->swap(tmp_buffer);
        }
    }
    
    /// return memory consumtion
    size_t mem_usage() const
    {
        return sizeof(capacity_) + sizeof(alloc_factor_) +
               capacity();
    }
    
private:
    /// Override from the base class
    void set_buf(unsigned char* buf, size_t size);

    /// compute number of words for the desired capacity
    static size_t compute_words(size_t capacity)
    {
        size_t words = (capacity / sizeof(bm::word_t))+1;
        return words;
    }
    
    void allocate(size_t new_capacity)
    {
        if (byte_buf_ && new_capacity <= capacity_)
            return;
        
        free_buffer();
    
        size_t words = compute_words(new_capacity);
        
        bm::word_t* p = allocator_type::allocate(words, 0);
        this->byte_buf_ = (unsigned char*) p;
        this->size_ = 0;
        alloc_factor_ = (unsigned)words;
        capacity_ = alloc_factor_ * sizeof(bm::word_t);
    }

    void free_buffer()
    {
        if (byte_buf_)
        {
            allocator_type::deallocate((bm::word_t*)byte_buf_, alloc_factor_);
            this->byte_buf_ = 0;
        }
    }
private:
    size_t         capacity_;     ///< current capacity
    size_t         alloc_factor_; ///< number of blocks allocated for buffer
};

/**
    Heap allocated scalar-type matrix

    This class can be seen as a matrix (row-column)
    access adaptor on top of the heap allocated buffer.

    @internal
*/
template<typename Val, unsigned ROWS, unsigned COLS, typename BVAlloc>
class heap_matrix
{
public:
    typedef BVAlloc                                          bv_allocator_type;
    typedef bm::byte_buffer<bv_allocator_type>               buffer_type;
    typedef Val                                              value_type;
    typedef unsigned                                         size_type;

    enum params
    {
        n_rows = ROWS,
        n_columns = COLS,
        size_in_bytes = sizeof(value_type) * COLS * ROWS,
        row_size_in_bytes = sizeof(value_type) * COLS
    };

    static unsigned rows() { return ROWS; }
    static unsigned cols() { return COLS; }

    /**
        By default object is constructed NOT allocated.
    */
    heap_matrix()
        : buffer_()
    {}

    heap_matrix(bool init)
        : buffer_(init ? size_in_bytes : 0)
    {
        if (init)
            buffer_.resize(size_in_bytes);
    }

    /**
        Post construction allocation, initialization
    */
    void init()
    {
        buffer_.resize(size_in_bytes);
    }
    
    value_type get(unsigned row_idx, unsigned col_idx) const
    {
        BM_ASSERT(row_idx < ROWS);
        BM_ASSERT(col_idx < COLS);
        BM_ASSERT(buffer_.size());
        const unsigned char* buf = buffer_.buf() + row_idx * row_size_in_bytes;
        return ((const value_type*)buf)[col_idx];
    }

    const value_type* row(unsigned row_idx) const
    {
        BM_ASSERT(row_idx < ROWS);
        BM_ASSERT(buffer_.size());
        const unsigned char* buf = buffer_.buf() + row_idx * row_size_in_bytes;
        return (const value_type*) buf;
    }

    value_type* row(unsigned row_idx)
    {
        BM_ASSERT(row_idx < ROWS);
        BM_ASSERT(buffer_.size());

        unsigned char* buf = buffer_.data() + row_idx * row_size_in_bytes;
        return (value_type*)buf;
    }

    /** memset all buffer to all zeroes */
    void set_zero()
    {
        ::memset(buffer_.data(), 0, size_in_bytes);
    }
    
    /*!  swap content
    */
    void swap(heap_matrix& other) BMNOEXEPT
    {
        buffer_.swap(other.buffer_);
    }
    
    /*!  move content from another matrix
    */
    void move_from(heap_matrix& other) BMNOEXEPT
    {
        buffer_.move_from(other.buffer_);
    }

    /** Get low-level buffer access */
    buffer_type& get_buffer() { return buffer_; }
    /** Get low-level buffer access */
    const buffer_type& get_buffer() const { return buffer_; }

    /*! remapping: vect[idx] = matrix[idx, vect[idx] ]
    */
    template<typename VECT_TYPE>
    void remap(VECT_TYPE* vect, size_type size) const
    {
        BM_ASSERT(size < ROWS);
        const unsigned char* buf = buffer_.buf();
        for (size_type i = 0; i < size; ++i)
        {
            const value_type* row = buf + i * row_size_in_bytes;
            VECT_TYPE v0 = vect[i];
            BM_ASSERT(size_type(v0) < COLS);
            value_type remap_v = row[unsigned(v0)];
            vect[i] = VECT_TYPE(remap_v);
        } // for i
    }

protected:
    buffer_type     buffer_;
};

} // namespace bm

#include "bmundef.h"

#endif
