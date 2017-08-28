#ifndef BMSPARSEVEC_H__INCLUDED__
#define BMSPARSEVEC_H__INCLUDED__
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


For more information please visit:   http://bmagic.sourceforge.net

*/

#include <memory.h>
#include <stdexcept>

#include "bmdef.h"
#include "bm.h"
#include "bmfunc.h"


namespace bm
{

template<class Val, class BV>
class sparse_vector
{
public:
    typedef Val                                      value_type;
    typedef bm::id_t                                 size_type;
    typedef BV                                       bvector_type;
    typedef bvector_type*                            bvector_type_ptr;
	typedef const value_type&                        const_reference;
    typedef typename BV::allocator_type              allocator_type;
    typedef typename bvector_type::allocation_policy allocation_policy_type;
    
    /*! Statistical information about  memory allocation details. */
    struct statistics : public bv_statistics
    {};

public:
    sparse_vector(allocation_policy_type ap = allocation_policy_type(),
                  size_type bv_size = bm::id_max,
                  const allocator_type&   alloc  = allocator_type());
    sparse_vector(const sparse_vector<Val, BV>& sv);
    
    sparse_vector& operator = (const sparse_vector<Val, BV>& sv)
    {
        clear();
        resize(sv.size());
        bv_size_ = sv.bv_size_;
        alloc_ = sv.alloc_;
    
        for (size_type i = 0; i < sizeof(Val)*8; ++i)
        {
            const bvector_type* bv = sv.plains_[i];
            if (bv)
                plains_[i] = new bvector_type(*bv);
        } // for i
        return *this;
    }
    
    
    ~sparse_vector();
    
    
    value_type operator[](size_type idx) const { return this->get(idx); }
    
    /*!
        \brief Import list of integers from a C style array
        \param arr  - source array
        \param size - source size
        \param offset - target index in the sparse vector
    */
    void import(const value_type* arr, size_type size, size_type offset = 0);
    
    
    /*! \brief return size of the vector
        \return size of sparse vector
    */
    size_type size() const;
    
    
    /*! \brief resize vector
        \param sz - new size
    */
    void resize(size_type sz);
    
    /*! \brief resize to zero, free memory
    */
    void clear();
    
    /*!
        \brief access specified element with bounds checking
        \param i - element index
        \return value of the element
    */
    value_type at(size_type i) const;
    
    /*!
        \brief get specified element without bounds checking
        \param i - element index
        \return value of the element
    */
    value_type get(bm::id_t idx) const;
    
    /*!
        \brief set specified element with bounds checking and automatic resize
        \param idx - element index
        \param v   - element value
    */
    void set(size_type idx, value_type v);

    /*!
        \brief push value back into vector
        \param v   - element value
    */
    void push_back(value_type v);
    
    /*!
        \brief check if another sparse vector has the same content and size
        \return true, if it is the same
    */
    bool equal(const sparse_vector<Val, BV>& sv) const;


    /*!
        \brief run memory optimization for all vector plains
        \param temp_block - pre-allocated memory block to avoid unnecessary re-allocs
        \param opt_mode - requested compression depth
        \param stat - memory allocation statistics after optimization
    */
    void optimize(bm::word_t* temp_block = 0,
                  typename bvector_type::optmode opt_mode = bvector_type::opt_compress,
                  typename sparse_vector<Val, BV>::statistics* stat = 0);
    
    /*!
        \brief join all with another sparse vector using OR operation
        \param sv - argument vector to join with
        \return slf reference
    */
    sparse_vector<Val, BV>& join(const sparse_vector<Val, BV>& sv);


    /*!
       @brief Calculates memory statistics.

       @param st - pointer on statistics structure to be filled in. 

       Function fills statistics structure containing information about how 
       this vector uses memory and estimation of max. amount of memory 
       bvector needs to serialize itself.

       @sa statistics
    */
    void calc_stat(struct sparse_vector<Val, BV>::statistics* st) const;
    

    /*!
        \brief get access to bit-plain, function checks and creates a plain
    */
    bvector_type_ptr get_plain(unsigned i);
    
    /*!
        \brief get total number of bit-plains in the vector
    */
    static unsigned plains() { return unsigned(sizeof(Val)*8); }
    
    /*!
        \brief get access to bit-plain as is (can return NULL)
    */
    bvector_type_ptr plain(unsigned i) { return plains_[i]; }
    const bvector_type_ptr plain(unsigned i) const { return plains_[i]; }
    
    /*!
        \brief free memory in bit-plain
    */
    void free_plain(unsigned i);
    
    /*!
        \brief clear range (assign bit 0 for all plains)
        \param left  - interval start
        \param right - interval end (closed interval)
    */
    sparse_vector<Val, BV>& clear_range(size_type left, size_type right);
private:

    /*! \brief free all internal vectors
    */
    void free_vectors();

    /*! \brief set value without checking boundaries
    */
protected:
    void set_value(size_type idx, value_type v);
private:
    
    size_type                bv_size_;
    allocator_type           alloc_;
    allocation_policy_type   ap_;
    
    bvector_type_ptr    plains_[sizeof(Val)*8];
    size_type           size_;
};

//---------------------------------------------------------------------


template<class Val, class BV>
sparse_vector<Val, BV>::sparse_vector(
        allocation_policy_type  ap,
        size_type               bv_size,
        const allocator_type&   alloc)
: bv_size_(bv_size),
  alloc_(alloc),
  ap_(ap),
  size_(0)
{
    ::memset(plains_, 0, sizeof(plains_));
}

//---------------------------------------------------------------------

template<class Val, class BV>
sparse_vector<Val, BV>::sparse_vector(const sparse_vector<Val, BV>& sv)
: bv_size_ (sv.bv_size_),
  alloc_(sv.alloc_),
  ap_(sv.ap_),
  size_(sv.size_)
{
    for (size_type i = 0; i < sizeof(Val)*8; ++i)
    {
        const bvector_type* bv = sv.plains_[i];
        if (bv)
            plains_[i] = new bvector_type(*bv);
        else
            plains_[i] = 0;
    } // for i
    
}

//---------------------------------------------------------------------


template<class Val, class BV>
sparse_vector<Val, BV>::~sparse_vector()
{
    free_vectors();
}



//---------------------------------------------------------------------


template<class Val, class BV>
void sparse_vector<Val, BV>::import(const value_type* arr,
                                    size_type         size,
                                    size_type         offset)
{
    if (size == 0)
        throw std::range_error("sparse vector range error");
    
    // cleap all plains in the range to provide corrrect import of 0 values
    this->clear_range(offset, offset + size - 1);
    
    size_type i;
    for (i = 0; i < size; ++i)
    {
        value_type v = arr[i];
        if (!v)
            continue;
        
        unsigned b_list[sizeof(Val)*8];
        unsigned bcnt = bm::bit_list_4(v, b_list);

        for (unsigned j = 0; j < bcnt; ++j)
        {
            unsigned p = b_list[j];
            bvector_type* bv = get_plain(p);
            bv->set_bit(i + offset);
        } // for j
    } // for i
    if (i > size_)
        size_ = i;
}

//---------------------------------------------------------------------

template<class Val, class BV>
typename sparse_vector<Val, BV>::size_type
sparse_vector<Val, BV>::size() const
{
    return size_;
}

//---------------------------------------------------------------------

template<class Val, class BV>
void sparse_vector<Val, BV>::resize(typename sparse_vector<Val, BV>::size_type sz)
{
    if (sz == size_)  // nothing to do
        return;
    
    if (!sz) // resize to zero is an equivalent of non-destructive deallocation
    {
        clear();
        return;
    }
    
    if (sz < size_) // vector shrink
    {
        // clear the tails
        this->clear_range(sz, size_-1);
    }
    
    size_ = sz;
}

//---------------------------------------------------------------------


template<class Val, class BV>
typename sparse_vector<Val, BV>::bvector_type_ptr
   sparse_vector<Val, BV>::get_plain(unsigned i)
{
    BM_ASSERT(i < (sizeof(Val)*8));

    bvector_type_ptr bv = plains_[i];
    if (!bv)
    {
        bv = new bvector_type(ap_.strat, ap_.glevel_len,
                              bv_size_,
                              alloc_);
        plains_[i] = bv;
    }
    return bv;
}

//---------------------------------------------------------------------

template<class Val, class BV>
typename sparse_vector<Val, BV>::value_type
sparse_vector<Val, BV>::at(typename sparse_vector<Val, BV>::size_type i) const
{
    if (i >= size_)
        throw std::range_error("sparse vector range error");
    return this->get(i);
}

//---------------------------------------------------------------------


template<class Val, class BV>
typename sparse_vector<Val, BV>::value_type
sparse_vector<Val, BV>::get(bm::id_t i) const
{
    BM_ASSERT(i < size_);
    
    value_type v = 0;
    for (unsigned j = 0; j < sizeof(Val)*8; ++j)
    {
        const bvector_type* bv = this->plains_[j];
        if (bv)
        {
            unsigned b = bv->test(i);
            v |= (b << j);
        }
    }
    return v;
}

//---------------------------------------------------------------------

template<class Val, class BV>
void sparse_vector<Val, BV>::set(size_type idx, value_type v)
{ 
    if (idx >= size_)
    {
        size_ = idx+1;
    }
    set_value(idx, v);
}

//---------------------------------------------------------------------

template<class Val, class BV>
void sparse_vector<Val, BV>::push_back(value_type v)
{
    set_value(size_, v);
    ++size_;
}

//---------------------------------------------------------------------

template<class Val, class BV>
void sparse_vector<Val, BV>::set_value(size_type idx, value_type v)
{
    // TODO: optimize to clear and set in just one pass

    // clear the plains
    for (unsigned i = 0; i < sizeof(Val) * 8; ++i)
    {
        bvector_type* bv = plains_[i];
        if (bv)
            bv->clear_bit(idx);
    }

    // set bits in plains
    unsigned b_list[sizeof(Val) * 8];
    unsigned bcnt = bm::bit_list_4(v, b_list);

    for (unsigned j = 0; j < bcnt; ++j)
    {
        unsigned p = b_list[j];
        bvector_type* bv = get_plain(p);
        bv->set_bit(idx);
    } // for j
}



//---------------------------------------------------------------------

template<class Val, class BV>
void sparse_vector<Val, BV>::clear()
{
    free_vectors();
    size_ = 0;
    ::memset(plains_, 0, sizeof(plains_));
}

//---------------------------------------------------------------------


template<class Val, class BV>
void sparse_vector<Val, BV>::free_vectors()
{
    for (size_type i = 0; i < sizeof(Val)*8; ++i)
        delete plains_[i];
}

//---------------------------------------------------------------------


template<class Val, class BV>
void sparse_vector<Val, BV>::free_plain(unsigned i)
{
    BM_ASSERT(i < sizeof(Val)*8);
    bvector_type* bv = plains_[i];
    delete bv;
    plains_[i] = 0;
}

//---------------------------------------------------------------------

template<class Val, class BV>
sparse_vector<Val, BV>&
sparse_vector<Val, BV>::clear_range(
    typename sparse_vector<Val, BV>::size_type left,
    typename sparse_vector<Val, BV>::size_type right)
{
    if (right < left)
    {
        return clear_range(right, left);
    }
    
    for (unsigned i = 0; i < sizeof(Val) * 8; ++i)
    {
        bvector_type* bv = plains_[i];
        if (bv)
        {
            bv->set_range(left, right, false);
        }
    } // for i
    
    return *this;
}

//---------------------------------------------------------------------

template<class Val, class BV>
void sparse_vector<Val, BV>::calc_stat(
     struct sparse_vector<Val, BV>::statistics* st) const
{
    BM_ASSERT(st);
    
	st->bit_blocks = st->gap_blocks = 0; 
	st->max_serialize_mem = st->memory_used = 0;
 
    for (unsigned j = 0; j < sizeof(Val)*8; ++j)
    {
        const bvector_type* bv = this->plains_[j];
        if (bv)
        {
            typename bvector_type::statistics stbv;
            bv->calc_stat(&stbv);
            
            st->bit_blocks += stbv.bit_blocks;
            st->gap_blocks += stbv.gap_blocks;
            st->max_serialize_mem += stbv.max_serialize_mem + 8;
            st->memory_used += stbv.memory_used;
            
        }
    } // for j
    // header accounting
    st->max_serialize_mem += 1 + 1 + 1 + 1 + 8 + (8 * sizeof(Val) * 8);

}

//---------------------------------------------------------------------

template<class Val, class BV>
void sparse_vector<Val, BV>::optimize(
    bm::word_t*                                  temp_block, 
    typename bvector_type::optmode               opt_mode,
    typename sparse_vector<Val, BV>::statistics* st)
{
    for (unsigned j = 0; j < sizeof(Val) * 8; ++j)
    {
        bvector_type* bv = this->plains_[j];
        if (bv)
        {
            if (!bv->any())  // empty vector?
            {
                delete this->plains_[j];
                this->plains_[j] = 0;
                continue;
            }
            
            typename bvector_type::statistics stbv;
            
            bv->optimize(temp_block, opt_mode, &stbv);
            
            
            if (st)
            {
                st->bit_blocks += stbv.bit_blocks;
                st->gap_blocks += stbv.gap_blocks;
                st->max_serialize_mem += stbv.max_serialize_mem + 8;
                st->memory_used += stbv.memory_used;
            }

            if (!temp_block)
            {
                typename bvector_type::blocks_manager_type& bv_bm =
                                                    bv->get_blocks_manager();
                temp_block = bv_bm.check_allocate_tempblock();
            }
        }
    } // for j

}

//---------------------------------------------------------------------

template<class Val, class BV>
sparse_vector<Val, BV>&
sparse_vector<Val, BV>::join(const sparse_vector<Val, BV>& sv)
{
    size_type arg_size = sv.size();
    if (size_ < arg_size)
    {
        resize(arg_size);
    }
    for (unsigned j = 0; j < sizeof(Val) * 8; ++j)
    {
        bvector_type* arg_bv = sv.plains_[j];
        if (arg_bv)
        {
            bvector_type* bv = this->plains_[j];
            if (!bv)
            {
                bv = get_plain(j);
            }
            *bv |= *arg_bv;
        }
    } // for j
    return *this;
}


//---------------------------------------------------------------------

template<class Val, class BV>
bool sparse_vector<Val, BV>::equal(const sparse_vector<Val, BV>& sv) const
{
    size_type arg_size = sv.size();
    if (size_ != arg_size)
    {
        return false;
    }
    for (unsigned j = 0; j < sizeof(Val) * 8; ++j)
    {
        const bvector_type* bv = plains_[j];
        const bvector_type* arg_bv = sv.plains_[j];
        if (!bv && !arg_bv)
            continue;
        // check if any not NULL and not empty
        if (!bv && arg_bv)
        {
            if (arg_bv->any())
                return false;
        }
        if (bv && !arg_bv)
        {
            if (bv->any())
                return false;
        }
        // both not NULL
        int cmp = bv->compare(*arg_bv);
        if (cmp != 0)
            return false;
    } // for j
    return true;
}

//---------------------------------------------------------------------


} // namespace bm

#include "bmundef.h"


#endif
