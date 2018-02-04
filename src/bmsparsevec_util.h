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
 
    \internal
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

//---------------------------------------------------------------------

template<class BV>
bvps_addr_resolver<BV>::bvps_addr_resolver()
: in_sync_(false)
{
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
}

//---------------------------------------------------------------------

template<class BV>
bool bvps_addr_resolver<BV>::resolve(bm::id_t id_from, bm::id_t* id_to) const
{
    BM_ASSERT(id_to);

    bool found = addr_bv_.test(id_from);
    if (!found)
    {
        *id_to = ~0u;
        return false;
    }

    if (in_sync_)
    {
        *id_to = addr_bv_.count_to(id_from, bv_blocks_);
    }
    else  // slow access
    {
        *id_to = addr_bv_.count_range(0, id_from);
    }
    return found;
}

//---------------------------------------------------------------------

template<class BV>
bool bvps_addr_resolver<BV>::get(bm::id_t id_from, bm::id_t* id_to) const
{
    BM_ASSERT(id_to);
    BM_ASSERT(in_sync_);

    bool found = addr_bv_.test(id_from);
    if (!found)
    {
        *id_to = ~0u;
        return found;
    }

    *id_to = addr_bv_.count_to(id_from, bv_blocks_);
    return found;
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


} // namespace bm

#include "bmundef.h"


#endif
