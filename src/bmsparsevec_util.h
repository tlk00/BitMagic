#ifndef BMSPARSEVEC_UTIL_H__INCLUDED__
#define BMSPARSEVEC_UTIL_H__INCLUDED__
/*
Copyright(c) 2002-2017 Anatoliy Kuznetsov(anatoliy_kuznetsov at yahoo.com)

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

#include <memory.h>
#include <stdexcept>

#include "bmdef.h"


namespace bm
{

/*!
    \brief Bit-bector prefix sum address resolver using bit-vector prefix sum
    as a space compactor.
 
    @internal
*/
template<class BV>
class bvps_addr_resolver
{
public:
    typedef bm::id_t                                 size_type;
    typedef BV                                       bvector_type;
    typedef typename bvector_type::blocks_count      bvector_blocks_psum_type;
public:
    bvps_addr_resolver();
    bvps_addr_resolver(const bvps_addr_resolver& addr_res);
    
    /*!
        \brief Resolve id to integer id (address)
     
        Resolve id to address with full sync and range checking
     
        \param id_from - input id to resolve
        \param id_to   - output id
     
        \return true if id is known and resolved successfully
    */
    bool resolve(bm::id_t id_from, bm::id_t* id_to) const;
    
    /*!
        \brief Resolve id to integer id (address) without sync check
     
        If prefix sum table is not sync - unexpected behavior
     
        \param id_from - input id to resolve
        \param id_to   - output id
     
        \return true if id is known and resolved successfully
    */
    bool get(bm::id_t id_from, bm::id_t* id_to) const;
    
    /*!
        \brief Set id (bit) to address resolver
    */
    void set(bm::id_t id_from);
    
    
    /*!
        \brief Re-calculate prefix sum table
        \param force - force recalculation even if it is already recalculated
    */
    void calc_prefix_sum(bool force = false);
    
    /*!
        \brief Re-calculate prefix sum table (same as calc_prefix_sum)
        \param force - force recalculation even if it is already recalculated
        @sa calc_prefix_sum
    */
    void sync(bool force = false) { calc_prefix_sum(force); }
    
    /*!
        \brief returns true if prefix sum table is in sync with the vector
    */
    bool in_sync() const { return in_sync_; }
    
    /*!
        \brief Unsync the prefix sum table
    */
    void unsync() { in_sync_ = false; }
    
    /*!
        \brief Get const reference to the underlying bit-vector
    */
    const bvector_type& get_bvector() const { return addr_bv_; }

    /*!
        \brief Get writable reference to the underlying bit-vector
     
        Access to underlying vector assumes modification and
        loss of prefix sum acceleration structures.
        Need to call sync() at the end of transaction.
    */
    bvector_type& get_bvector()  { unsync(); return addr_bv_; }
    
    /*!
        \brief optimize underlying bit-vector
    */
    void optimize(bm::word_t* temp_block = 0);

    
private:
    bvector_type              addr_bv_;   ///< bit-vector for id translation
    bvector_blocks_psum_type  bv_blocks_; ///< prefix sum for fast translation
    bool                      in_sync_;   ///< flag if prefix sum is in-sync with vector
};


/*!
    \brief sparse vector based address resolver
    (no space compactor, just bit-plane compressors provided by sparse_vector)
 
    @internal
*/
template<class SV>
class sv_addr_resolver
{
public:
    typedef bm::id_t                              size_type;
    typedef SV                                    sparse_vector_type;
    typedef typename SV::bvector_type             bvector_type;
public:
    sv_addr_resolver();
    sv_addr_resolver(const sv_addr_resolver& addr_res);
    
    /*!
        \brief Resolve id to integer id (address)
     
        \param id_from - input id to resolve
        \param id_to   - output id
     
        \return true if id is known and resolved successfully
    */
    bool resolve(bm::id_t id_from, bm::id_t* id_to) const;
    
    /*!
        \brief Resolve id to integer id (address)
     
        \param id_from - input id to resolve
        \param id_to   - output id
     
        \return true if id is known and resolved successfully
    */
    bool get(bm::id_t id_from, bm::id_t* id_to) const;
    
    /*!
        \brief Set id (bit) to address resolver
    */
    void set(bm::id_t id_from);
    
    /*!
        \brief optimize underlying sparse vectors
    */
    void optimize(bm::word_t* temp_block = 0);

    /*!
    \brief Get const reference to the underlying bit-vector of set values
    */
    const bvector_type& get_bvector() const { return set_flags_bv_; }

private:
    bvector_type              set_flags_bv_;   ///< bit-vector of set flags
    sparse_vector_type        addr_sv_;     ///< sparse vector for address map
    unsigned                  max_addr_;    ///< maximum allocated address/index
};


/**
    \brief Compressed (sparse collection of objects)
    @internal
*/
template<class Value, class BV>
class compressed_collection
{
public:
    typedef bm::id_t                             key_type;
    typedef bm::id_t                             address_type;
    typedef Value                                value_type;
    typedef Value                                mapped_type;
    typedef BV                                   bvector_type;
    typedef std::vector<value_type>              container_type;
    typedef bm::bvps_addr_resolver<bvector_type> address_resolver_type;
    
public:
    compressed_collection();
    
    /**
        Add new value to compressed collection.
        Must be added in sorted key order (growing).
     
        Unsorted will not be added!
     
        \return true if value was added (does not violate sorted key order)
    */
    bool push_back(key_type key, const value_type& val);
    
    /**
        find and return const associated value (with bounds/presense checking)
    */
    const value_type& at(key_type key) const;

    /**
        find and return associated value (with bounds/presense checking)
    */
    value_type& at(key_type key);

    /** Checkpoint method to prepare collection for reading
    */
    void sync();
    
    /** perform memory optimizations/compression
    */
    void optimize();

    /** Resolve key address (index) in the dense vector
    */
    bool resolve(key_type key, address_type* addr) const;
    
    /** Get access to associated value by resolved address
    */
    const value_type& get(address_type addr) const;
protected:
    void throw_range_error(const char* err_msg) const;
    
private:
    address_resolver_type             addr_res_;    ///< address resolver
    container_type                    dense_vect_;  ///< compressed space container
    key_type                          last_add_;    ///< last added element
};

/**
    \brief Compressed (sparse collection of objects)
    @internal
*/
template<class BV>
class compressed_buffer_collection :
               public compressed_collection<typename serializer<BV>::buffer, BV>
{
public:
    typedef typename serializer<BV>::buffer     buffer_type;
    typedef compressed_collection<typename serializer<BV>::buffer, BV> parent_type;
public:

    bool move_buffer(typename parent_type::key_type key, buffer_type& buffer)
    {
        bool added = push_back(key, buffer_type());
        if (!added)
            return added;
        buffer_type& buf = this->at(key);
        buf.swap(buffer);
        return added;
    }

};


//---------------------------------------------------------------------

template<class BV>
bvps_addr_resolver<BV>::bvps_addr_resolver()
: in_sync_(false)
{
    addr_bv_.init();
}

//---------------------------------------------------------------------

template<class BV>
bvps_addr_resolver<BV>::bvps_addr_resolver(const bvps_addr_resolver& addr_res)
: addr_bv_(addr_res.addr_bv_),
  in_sync_(addr_res.in_sync_)
{
    if (in_sync_)
    {
        bv_blocks_.copy_from(addr_res.bv_blocks_);
    }
    addr_bv_.init();
}

//---------------------------------------------------------------------

template<class BV>
bool bvps_addr_resolver<BV>::resolve(bm::id_t id_from, bm::id_t* id_to) const
{
    BM_ASSERT(id_to);

    if (in_sync_)
    {
        *id_to = addr_bv_.count_to_test(id_from, bv_blocks_);
    }
    else  // slow access
    {
        bool found = addr_bv_.test(id_from);
        if (!found)
        {
            *id_to = 0;
        }
        else
        {
            *id_to = addr_bv_.count_range(0, id_from);
        }
    }
    return (bool) *id_to;
}

//---------------------------------------------------------------------

template<class BV>
bool bvps_addr_resolver<BV>::get(bm::id_t id_from, bm::id_t* id_to) const
{
    BM_ASSERT(id_to);
    BM_ASSERT(in_sync_);

    *id_to = addr_bv_.count_to_test(id_from, bv_blocks_);
    return (bool)*id_to;
}

//---------------------------------------------------------------------

template<class BV>
void bvps_addr_resolver<BV>::set(bm::id_t id_from)
{
    BM_ASSERT(!addr_bv_.test(id_from));
    
    in_sync_ = false;
    addr_bv_.set(id_from);
}


//---------------------------------------------------------------------


template<class BV>
void bvps_addr_resolver<BV>::calc_prefix_sum(bool force)
{
    if (in_sync_ && force == false)
        return;  // nothing to do
    
    addr_bv_.running_count_blocks(&bv_blocks_); // compute popcount prefix list
    in_sync_ = true;
}

//---------------------------------------------------------------------

template<class BV>
void bvps_addr_resolver<BV>::optimize(bm::word_t* temp_block)
{
    addr_bv_.optimize(temp_block);
}

//---------------------------------------------------------------------




template<class SV>
sv_addr_resolver<SV>::sv_addr_resolver()
: max_addr_(0)
{
    set_flags_bv_.init();
}

//---------------------------------------------------------------------

template<class SV>
sv_addr_resolver<SV>::sv_addr_resolver(const sv_addr_resolver& addr_res)
: set_flags_bv_(addr_res.set_flags_bv_),
  addr_sv_(addr_res.addr_sv_),
  max_addr_(addr_res.max_addr_)
{
}

//---------------------------------------------------------------------

template<class SV>
bool sv_addr_resolver<SV>::resolve(bm::id_t id_from, bm::id_t* id_to) const
{
    BM_ASSERT(id_to);
    
    bool found = set_flags_bv_.test(id_from);
    *id_to = found ? addr_sv_.at(id_from) : 0;
    return found;
}

//---------------------------------------------------------------------

template<class SV>
void sv_addr_resolver<SV>::set(bm::id_t id_from)
{
    bool found = set_flags_bv_.test(id_from);
    if (!found)
    {
        set_flags_bv_.set(id_from);
        ++max_addr_;
        addr_sv_.set(id_from, max_addr_);
    }
}

//---------------------------------------------------------------------

template<class SV>
void sv_addr_resolver<SV>::optimize(bm::word_t* temp_block)
{
    set_flags_bv_.optimize(temp_block);
    addr_sv_.optimize(temp_block);
}

//---------------------------------------------------------------------



template<class Value, class BV>
compressed_collection<Value, BV>::compressed_collection()
: last_add_(0)
{
}


//---------------------------------------------------------------------

template<class Value, class BV>
bool compressed_collection<Value, BV>::push_back(key_type key, const value_type& val)
{
    if (dense_vect_.empty()) // adding first one
    {
    }
    else
    if (key <= last_add_)
    {
        BM_ASSERT(0);
        return false;
    }
    
    addr_res_.set(key);
    last_add_ = key;
    dense_vect_.push_back(val);
    return true;
}

//---------------------------------------------------------------------

template<class Value, class BV>
void compressed_collection<Value, BV>::sync()
{
    addr_res_.sync();
}

//---------------------------------------------------------------------

template<class Value, class BV>
void compressed_collection<Value, BV>::optimize()
{
    addr_res_.optimize();
}

//---------------------------------------------------------------------

template<class Value, class BV>
bool compressed_collection<Value, BV>::resolve(key_type key,
                                               address_type* addr) const
{
    bool found = addr_res_.resolve(key, addr);
    *addr -= found;
    return found;
}

//---------------------------------------------------------------------


template<class Value, class BV>
const typename compressed_collection<Value, BV>::value_type&
compressed_collection<Value, BV>::get(address_type addr) const
{
    return dense_vect_.at(addr);
}

//---------------------------------------------------------------------

template<class Value, class BV>
const typename compressed_collection<Value, BV>::value_type&
compressed_collection<Value, BV>::at(key_type key) const
{
    address_type idx;
    bool found = addr_res_.resolve(key, &idx);
    if (!found)
    {
        throw_range_error("compressed collection item not found");
    }
    return get(idx-1);
}

//---------------------------------------------------------------------

template<class Value, class BV>
typename compressed_collection<Value, BV>::value_type&
compressed_collection<Value, BV>::at(key_type key)
{
    address_type idx;
    bool found = addr_res_.resolve(key, &idx);
    if (!found)
    {
        throw_range_error("compressed collection item not found");
    }
    return dense_vect_.at(idx-1);
}


//---------------------------------------------------------------------

template<class Value, class BV>
void compressed_collection<Value, BV>::throw_range_error(const char* err_msg) const
{
#ifndef BM_NO_STL
    throw std::range_error(err_msg);
#else
    BM_ASSERT_THROW(false, BM_ERR_RANGE);
#endif
}

//---------------------------------------------------------------------

} // namespace bm

#include "bmundef.h"


#endif
