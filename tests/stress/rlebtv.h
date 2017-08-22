/*
Copyright (c) 2002,2003 Anatoliy Kuznetsov.

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
*/

#ifndef RLEBTV__H__INCLUDED__
#define RLEBTV__H__INCLUDED__

#include <cassert>

/*
    GAP vector is designed to be very small.
    The size of buffer is fixed 128 short values.
    First word is a start flag.
    Used for debugging purposes.
*/

class gap_vector
{
public:
    gap_vector(gap_word_t init = 0)
    {
        m_buf[0] = 1 << 3;  // initial length = 1;

        if (init)
        {
            m_buf[0] ^= 1;
        }
        m_buf[1] = gap_max_bits - 1;
    }

    /// Checks if bit pos 1 or 0. Returns 0 if 0 and non zero otherwise.
    int is_bit_true(unsigned pos) const;
    int test(unsigned pos) const;

    /// Sets bit number pos to 1
    /// returns true if operation sucessful, false if GAP threashold reached and
    /// we need to convert block to true bitblock
    bool set_bit(unsigned pos);

    /// Sets bit number pos to 0
    bool clear_bit(unsigned pos);

    /// Counts number of bits ON 
    unsigned bit_count() const;
    
    unsigned count_range(unsigned left, 
                         unsigned right,
                         unsigned* d=0) const;

    void control() const;
    
    void convert_to_bitset(unsigned* dest) const;

    int combine_and(const gap_word_t* other);

    const gap_word_t* get_buf() const { return m_buf; }

    void invert() { *m_buf = 1 - *m_buf; }

    void temp_invert() const;

    void combine_or(const gap_vector& vect);
    void combine_sub(const gap_vector& vect);
    void combine_xor(const gap_vector& vect);

    gap_word_t* get_buf() { return m_buf; }

    int compare(const gap_vector& vect);

private:
    gap_word_t   m_buf[bm::gap_max_buff_len+3];    
};


inline void gap_vector::combine_or(const gap_vector& vect)
{

    gap_word_t tmp_buf[bm::gap_max_buff_len * 3] = {0,}; // temporary result
    unsigned len, dsize;
    gap_word_t* res = bm::gap_operation_or(m_buf, vect.get_buf(), tmp_buf, dsize);
    len = bm::gap_length(res);

    if (res == tmp_buf)
    {
        assert(len);
        ::memcpy(m_buf, tmp_buf, (len+1)*sizeof(gap_word_t));
    }
    else
    {
        assert(res == m_buf);
    }

}

inline void gap_vector::combine_xor(const gap_vector& vect)
{

    gap_word_t tmp_buf[gap_max_buff_len * 3] = {0,}; // temporary result
    unsigned len, dsize;
    gap_word_t* res = 
        bm::gap_operation_xor(m_buf, vect.get_buf(), tmp_buf, dsize);
    len = bm::gap_length(res);

    if (res == tmp_buf)
    {
        assert(len);
        ::memcpy(m_buf, tmp_buf, (len+1)*sizeof(gap_word_t));
    }
    else
    {
        assert(res == m_buf);
    }

}


inline void gap_vector::combine_sub(const gap_vector& vect)
{
    unsigned dsize;
    gap_word_t tmp_buf[bm::gap_max_buff_len * 3] = {0,}; // temporary result
    gap_word_t* res = 
        bm::gap_operation_sub(m_buf, vect.get_buf(), tmp_buf, dsize);
    unsigned len = bm::gap_length(res);

    if (res == tmp_buf)
    {
        assert(len);
        ::memcpy(m_buf, tmp_buf, (len+1) * sizeof(gap_word_t));
    }
    else
    {
        assert(res == m_buf);
    }

}



inline void gap_vector::temp_invert() const
{
    gap_word_t* buf = const_cast<gap_word_t*>(m_buf);
    *buf ^= 1;
}

inline void gap_vector::control() const
{
    unsigned sum = bm::gap_control_sum(m_buf);
    if(sum != bm::gap_max_bits-1)
    {
        assert(0);
        cout << "GAP Control failed." << endl;
        exit(1);
    }
    unsigned len = bm::gap_length(m_buf);
    gap_word_t prev = m_buf[1];
    for (unsigned i = 2; i < len; ++i)
    {
        if (i != len-1)
        {
            if (m_buf[i] <= prev)
            {
                assert(0);
                cout << "GAP sequence control failed." << endl;
                exit(1);
            }
        }
        else
        {
            assert(m_buf[i] == bm::gap_max_bits-1);
        }
        prev = m_buf[i];
    }
}

inline int gap_vector::is_bit_true(unsigned pos) const
{
    return bm::gap_test(m_buf, pos);
}

inline int gap_vector::test(unsigned pos) const
{
    return bm::gap_test(m_buf, pos);
}


unsigned gap_vector::bit_count() const
{
    return bm::gap_bit_count(m_buf);
}

unsigned gap_vector::count_range(unsigned left, 
                                 unsigned right,
                                 unsigned*) const
{
    return bm::gap_bit_count_range(m_buf, (gap_word_t)left, (gap_word_t)right);
}



// returns destination size


int gap_vector::combine_and(const gap_word_t* other)
{

    gap_word_t tmp_buf[gap_max_buff_len * 3] = {0,}; // temporary result
    unsigned dsize;

    gap_word_t* res = bm::gap_operation_and(m_buf, other, tmp_buf, dsize);
    unsigned len = bm::gap_length(res);

    if (res == tmp_buf)
    {
        assert(len);
        ::memcpy(m_buf, tmp_buf, (len+1) * sizeof(gap_word_t));
    }
    else
    {
        assert(res == m_buf);
    }
    return 0;

}


inline bool gap_vector::set_bit(unsigned pos)
{
    unsigned is_set;
    bm::gap_set_value(1, m_buf, pos, &is_set);
    return is_set!=0;
}

inline bool gap_vector::clear_bit(unsigned pos)
{
    unsigned is_set;
    bm::gap_set_value(0, m_buf, pos, &is_set);
    return is_set!=0;
}

inline void gap_vector::convert_to_bitset(unsigned* dest) const
{
    bm::gap_convert_to_bitset(dest, m_buf);
}

inline int gap_vector::compare(const gap_vector& vect)
{
    return bm::gapcmp(m_buf, vect.get_buf());
}


#endif


