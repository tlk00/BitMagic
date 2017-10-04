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
#include "bmtrans.h"
#include "bmalgo_impl.h"


namespace bm
{

/** \defgroup svector Sparse vector
    Sparse vector for integer types using bit transposition transform
    \ingroup bmagic
 */


/*!
   \brief sparse vector with runtime compression using bit transposition method
   \ingroup svector
*/
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
    /*!
        \brief Sparse vector constructor
        \param ap - allocation strategy for underlying bit-vectors
        Default allocation policy uses BM_BIT setting (fastest access)
        \param bv_max_size - maximum possible size of underlying bit-vectors
        Please note, this is NOT size of svector itself, it is dynamic upper limit
        which should be used very carefully if we surely know the ultimate size
        \param alloc - allocator for bit-vectors
        
        \sa bvector<>
        \sa bm::bvector<>::allocation_policy
        \sa bm::startegy
    */
    sparse_vector(allocation_policy_type ap = allocation_policy_type(),
                  size_type bv_max_size = bm::id_max,
                  const allocator_type&   alloc  = allocator_type());
    
    /*! copy-ctor */
    sparse_vector(const sparse_vector<Val, BV>& sv);


#ifndef BMNOCXX11
    /*! move-ctor */
    sparse_vector(sparse_vector<Val, BV>&& sv) BMNOEXEPT;


    /*! move assignmment operator */
    sparse_vector<Val,BV>& operator = (sparse_vector<Val, BV>&& sv) BMNOEXEPT
    {
        if (this != &sv)
        {
            clear();
            swap(sv);
        }
        return *this;
    }
#endif

    /*! copy assignmment operator */
    sparse_vector<Val,BV>& operator = (const sparse_vector<Val, BV>& sv)
    {
        if (this != &sv)
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
        }
        return *this;
    }
    
    
    ~sparse_vector() BMNOEXEPT;
    
    /*! \brief content exchange
    */
    void swap(sparse_vector<Val, BV>& sv) BMNOEXEPT;
    
    /*!
        \brief get specified element without bounds checking
        \param idx - element index
        \return value of the element
    */
    value_type operator[](size_type idx) const { return this->get(idx); }
    
    /*!
        \brief Import list of elements from a C-style array
        \param arr  - source array
        \param size - source size
        \param offset - target index in the sparse vector
    */
    void import(const value_type* arr, size_type size, size_type offset = 0);

    /*!
        \brief Bulk export list of elements to a C-style array
        
        For efficiency, this is left as a low level function, 
        it does not do any bounds checking on the target array, it will 
        override memory and crash if you are not careful with allocation
        
        This method uses cache-miss optimized algorithm which is by far faster 
        than random element access functions.
     
        \param arr  - dest array
        \param size - dest size
        \param offset - target index in the sparse vector to export from
        \param zero_mem - set to false if target array is pre-initialized 
                          with 0s to avoid performance penalty
     
        \return number of exported elements
    */
    size_type extract(value_type* arr,
                      size_type size,
                      size_type offset = 0,
                      bool      zero_mem = true);

    
    
    /*! \brief return size of the vector
        \return size of sparse vector
    */
    size_type size() const;
    
    
    /*! \brief return true if vector is empty
        \return true if empty
    */
    bool empty() const;
    
    
    /*! \brief resize vector
        \param sz - new size
    */
    void resize(size_type sz);
    
    /*! \brief resize to zero, free memory
    */
    void clear() BMNOEXEPT;
    
    /*!
        \brief access specified element with bounds checking
        \param idx - element index
        \return value of the element
    */
    value_type at(size_type idx) const;
    
    /*!
        \brief get specified element without bounds checking
        \param idx - element index
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
    void free_vectors() BMNOEXEPT;

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
        size_type               bv_max_size,
        const allocator_type&   alloc)
: bv_size_(bv_max_size),
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
    if (this != &sv)
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
}

//---------------------------------------------------------------------
#ifndef BMNOCXX11

template<class Val, class BV>
sparse_vector<Val, BV>::sparse_vector(sparse_vector<Val, BV>&& sv) BMNOEXEPT
{
    if (this != &sv)
    {
        bv_size_ = 0;
        alloc_ = sv.alloc_;
        ap_ = sv.ap_;
        size_ = sv.size_;
        
        for (size_type i = 0; i < sizeof(Val)*8; ++i)
        {
            plains_[i] = sv.plains_[i];
            sv.plains_[i] = 0;
        }
        sv.size_ = 0;
    }
}

#endif



//---------------------------------------------------------------------

template<class Val, class BV>
sparse_vector<Val, BV>::~sparse_vector() BMNOEXEPT
{
    free_vectors();
}

//---------------------------------------------------------------------

template<class Val, class BV>
void sparse_vector<Val, BV>::swap(sparse_vector<Val, BV>& sv) BMNOEXEPT
{
    if (this != &sv)
    {
        bm::xor_swap(bv_size_, sv.bv_size_);
        
        allocator_type alloc_tmp = alloc_;
        alloc_ = sv.alloc_;
        sv.alloc_ = alloc_tmp;
        
        allocation_policy_type ap_tmp = ap_;
        ap_ = sv.ap_;
        sv.ap_ = ap_tmp;
        
        for (size_type i = 0; i < sizeof(Val)*8; ++i)
        {
            bvector_type* bv_tmp = plains_[i];
            plains_[i] = sv.plains_[i];
            sv.plains_[i] = bv_tmp;
        } // for i
        
        bm::xor_swap(size_, sv.size_);        
    }
}


//---------------------------------------------------------------------

template<class Val, class BV>
void sparse_vector<Val, BV>::import(const value_type* arr,
                                    size_type         size,
                                    size_type         offset)
{
    unsigned b_list[sizeof(Val)*8];
    unsigned row_len[sizeof(Val)*8] = {0, };
    
    const unsigned transpose_window = 256;
    bm::tmatrix<bm::id_t, sizeof(Val)*8, transpose_window> tm; // matrix accumulator
    
    if (size == 0)
        throw std::range_error("sparse vector range error");
    
    // clear all plains in the range to provide corrrect import of 0 values
    this->clear_range(offset, offset + size - 1);
    
    // transposition algorithm uses bitscen to find index bits and store it
    // in temporary matrix (list for each bit plain), matrix here works
    // when array gets to big - the list gets loaded into bit-vector using
    // bulk load algorithm, which is faster than single bit access
    //
    
    size_type i;
    for (i = 0; i < size; ++i)
    {
        unsigned bcnt = bm::bitscan_popcnt(arr[i], b_list);
        const unsigned bit_idx = i + offset;
        
        for (unsigned j = 0; j < bcnt; ++j)
        {
            unsigned p = b_list[j];
            unsigned rl = row_len[p];
            tm.row(p)[rl] = bit_idx;
            row_len[p] = ++rl;
            
            if (rl == transpose_window)
            {
                bvector_type* bv = get_plain(p);
                const bm::id_t* r = tm.row(p);
                bm::combine_or(*bv, r, r + rl);
                row_len[p] = 0;
                tm.row(p)[0] = 0;
            }
        } // for j
    } // for i
    
    // process incomplete transposition lines
    //
    for (unsigned k = 0; k < tm.rows(); ++k)
    {
        unsigned rl = row_len[k];
        if (rl)
        {
            bvector_type* bv = get_plain(k);
            const bm::id_t* r = tm.row(k);
            bm::combine_or(*bv, r, r + rl);
        }
    } // for k
    
    
    if (i > size_)
        size_ = i;
}


//---------------------------------------------------------------------

template<class Val, class BV>
typename sparse_vector<Val, BV>::size_type
sparse_vector<Val, BV>::extract(value_type* arr,
                                size_type   size,
                                size_type   offset,
                                bool        zero_mem)
{
    if (size == 0)
        return 0;
    if (zero_mem)
        ::memset(arr, 0, sizeof(value_type)*size);
    
    size_type start = offset;
    size_type end = start + size;
    if (end > size_)
        end = size_;
    
	bool masked_scan = !(offset == 0 && size == this->size());

    if (masked_scan)
    {
        bvector_type bv_mask;
        for (size_type i = 0; i < sizeof(Val)*8; ++i)
        {
            const bvector_type* bv = plains_[i];
            if (bv)
            {
                value_type mask = (1 << i);            
                bv_mask.set_range(offset, end - 1);
                bv_mask.bit_and(*bv);
                for (typename BV::enumerator en(&bv_mask, 0); en.valid(); ++en)
                {
                    size_type idx = *en - offset;
                    BM_ASSERT(idx < size);
                    arr[idx] |= mask;
                } // for
                bv_mask.clear();
            }
        } // for i
    }
    else
    {
        for (size_type i = 0; i < sizeof(Val)*8; ++i)
        {
            const bvector_type* bv = plains_[i];
            if (bv)
            {
                value_type mask = (1 << i);            
                for (typename BV::enumerator en(bv, 0); en.valid(); ++en)
                {
                    size_type idx = *en;
                    BM_ASSERT(idx < size);
                    arr[idx] |= mask;
                } // for
            }
        } // for i
    }

    return 0;
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
bool sparse_vector<Val, BV>::empty() const
{
    return (size_ == 0);
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
sparse_vector<Val, BV>::at(typename sparse_vector<Val, BV>::size_type idx) const
{
    if (idx >= size_)
        throw std::range_error("sparse vector range error");
    return this->get(idx);
}

//---------------------------------------------------------------------


template<class Val, class BV>
typename sparse_vector<Val, BV>::value_type
sparse_vector<Val, BV>::get(bm::id_t i) const
{
    BM_ASSERT(i < size_);
    
    value_type v = 0;
    const bvector_type* bv;
    for (unsigned j = 0; j < sizeof(Val)*8; ++j)
    {
        if ((bv = this->plains_[j])!=0)   v |= ((bv->test(i))<<j);
        if ((bv = this->plains_[++j])!=0) v |= ((bv->test(i))<<j);
        if ((bv = this->plains_[++j])!=0) v |= ((bv->test(i))<<j);
        if ((bv = this->plains_[++j])!=0) v |= ((bv->test(i))<<j);
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
void sparse_vector<Val, BV>::clear() BMNOEXEPT
{
    free_vectors();
    size_ = 0;
    ::memset(plains_, 0, sizeof(plains_));
}

//---------------------------------------------------------------------


template<class Val, class BV>
void sparse_vector<Val, BV>::free_vectors() BMNOEXEPT
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
