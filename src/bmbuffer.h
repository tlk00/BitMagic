#ifndef BMBUFFER__H__INCLUDED__
#define BMBUFFER__H__INCLUDED__
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

You have to explicitly mention BitMagic project in any derivative product,
its WEB Site, published materials, articles or any other work derived from this
project or based on our code or know-how.

For more information please visit:  http://bitmagic.io

*/


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
    byte_buffer_ptr(unsigned char* buf, size_t size)
        : byte_buf_(buf), size_(size)
    {}
    
    /// Set buffer pointer
    void set_buf(unsigned char* buf, size_t size)
    {
        byte_buf_ = buf; size_= size;
    }

    /// Get buffer size
    size_t size() const { return size_; }
    
    /// Get read access to buffer memory
    const unsigned char* buf() const { return byte_buf_; }

    /// Get write access to buffer memory
    unsigned char* data() { return byte_buf_; }

    bool operator==(const byte_buffer_ptr& buf) const
    {
        if (this == &buf)
            return true;
        if (size_ != buf.size_)
            return false;
        int cmp = ::memcmp(byte_buf_, buf.byte_buf_, size_);
        return (cmp  == 0);
    }

protected:
    unsigned char* byte_buf_;     ///< byte buffer pointer to hold data
    size_t         size_;         ///< current buffer size
};

} // namespace bm


#endif
