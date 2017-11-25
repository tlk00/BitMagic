/*
Copyright (c) 2003 Anatoliy Kuznetsov.
*/
#ifndef BM_COMPRESSED_VECTOR__H__INCLUDED__
#define BM_COMPRESSED_VECTOR__H__INCLUDED__

#include "bmmalloc.h"


namespace bm
{


template<class V, class S, class A=bm::malloc_allocator<V> >
class compressed_scalar_vector
{
public:
    typedef S                       set_type;
    typedef S                       bvector_type;
    typedef V                       value_type;
    typedef bm::id_t                key_type;

public:
    compressed_scalar_vector(S* set=0);
    compressed_scalar_vector(unsigned size, S* set=0);
    compressed_scalar_vector(const compressed_scalar_vector<V, S, A>& vect);    
    compressed_scalar_vector<V, S, A>& operator=(const compressed_scalar_vector<V, S, A>& vect);
    
    ~compressed_scalar_vector();
    
    unsigned size() const     { return size_; }
    unsigned capacity() const { return capacity_; }
    
    void set(unsigned idx, const V& value);
    
    void resize(unsigned new_size);
private:
    void reserve_new(unsigned capacity);    
private:
    S*         set_;
    V*         data_;
    unsigned   size_;
    unsigned   capacity_;
    unsigned   max_idx_;
    unsigned*  set_block_count_cache_;    
};

// ---------------------------------------------------------------------------

template<class V, class S, class A>
compressed_scalar_vector<V, S, A>::compressed_scalar_vector(S* set)
: set_(set ? set : new S()),
  data_(0),
  size_(0),
  capacity_(0),
  max_idx_(0),
  set_block_count_cache_(0)
{
    assert(!set_->any());
}

template<class V, class S, class A>
compressed_scalar_vector<V, S, A>::compressed_scalar_vector(unsigned size, S* set)
: set_(set ? set : new S()),
  max_idx_(0),
  set_block_count_cache_(0)
{
    assert(!set_->any());
    data_ = size ? A::allocate(size, 0) : 0;    
    size_ = capacity_ = size_;
}

template<class V, class S, class A>
compressed_scalar_vector<V, S, A>::compressed_scalar_vector
                               (const compressed_scalar_vector<V, S, A>& vect)
: set_(new S(vect.set_),
  max_idx_(vect.max_idx_),
  set_block_count_cache_(0)
{
    capacity_ = size_ = vect.size();
    if (size_)
    {
        data_ = A::allocate(size_, 0) : 0;
        ::memcpy(data_, vect.data_, sizeof(V) * size_);
    }
    else
    {
        data_ = 0;
    }
}

template<class V, class S, class A>
compressed_scalar_vector<V, S, A>::~compressed_scalar_vector()
{
    delete set_;
    if (data_)
        A::deallocate(data_, capacity_);
    delete set_block_count_cache_;
}


template<class V, class S, class A>
compressed_scalar_vector<V, S, A>& operator=(const compressed_scalar_vector<V, S, A>& vect)
{
    *set_ = *(vect.set_);
    
    if (capacity_ < vect.size())
    {
        A::deallocate(data_, capacity_);
        capacity_ = size_ = vect.size();
        data_ = A::allocate(size_, 0);        
    }
    else
    {
        size_ = vect.size();
    }
    ::memcpy(data_, vect.data_, sizeof(V) * size_);
}


template<class V, class S, class A>
void compressed_scalar_vector<V, S, A>::set(unsigned idx, const V& value)
{
    assert(idx < size_);
    
    unsigned i;
    if (idx >= max_idx_)
    {
        i = set_->count();
        max_idx_ = idx;
    }
    else
    {
        i = set_->count_range(0, idx);
    }
    
    if (set_->test(idx))
    {
        --i;
    }
    else
    {
        set_->set(idx);
        if (i >= size())
        {
            resize(i+1);
        }
    }
    data_[i] = value;
}


template<class V, class S, class A>
void compressed_scalar_vector<V, S, A>::resize(unsigned new_size)
{
    if (new_size > capacity)
    {
        reserve_new(new_size);
    }
    size_ = new_size;
}

template<class V, class S, class A>
void compressed_scalar_vector<V, S, A>::reserve_new(unsigned capacity)
{
    V* data_new = A::allocate(capacity, 0);
    ::memcpy(data_, data_new, sizeof(V) * size_);
    capacity_ = capacity;
}


} // namespace bm

#endif

