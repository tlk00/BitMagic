#ifndef BMALGO__H__INCLUDED__
#define BMALGO__H__INCLUDED__
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

#include "bm.h"
#include "bmfunc.h"
#include "bmdef.h"

#include "bmalgo_impl.h"



namespace bm
{

/**
    @brief bit-vector visitor scanner to traverse each 1 bit using C++ visitor
 
    @param bv - bit vector to scan
    @param bit_functor (should support add_bits() and add_range() methods
 
    \ingroup setalgo
*/
template<class BV, class Func>
void for_each_bit(const BV&    bv,
                  Func&        bit_functor)
{
    const typename BV::blocks_manager_type& bman = bv.get_blocks_manager();
    bm::word_t*** blk_root = bman.top_blocks_root();
    
    unsigned tsize = bman.top_block_size();
    for (unsigned i = 0; i < tsize; ++i)
    {
        bm::word_t** blk_blk = blk_root[i];
        if (!blk_blk)
        {
            continue;
        }
        unsigned r = i * bm::set_array_size;
        for (unsigned j = 0; j < bm::set_array_size; ++j)
        {
            const bm::word_t* block = blk_blk[j];
            if (block)
            {
                unsigned nb = r + j;
                if (BM_IS_GAP(block))
                {
                    bm::for_each_gap_blk(BMGAP_PTR(block),
                                         nb * bm::bits_in_block,
                                         bit_functor);
                }
                else
                {
                    // TODO: optimize FULL BLOCK ADDRESS scan
                    block = BLOCK_ADDR_SAN(block);
                    
                    bm::for_each_bit_blk(block, nb * bm::bits_in_block,
                                         bit_functor);
                }
            }
        } // for j
    }  // for i

}


/**
    @brief bit-vector visitor scanner to traverse each 1 bit using C callback
 
    @param bv - bit vector to scan
    @param handle_ptr - handle to private memory used by callback
    @param callback_ptr - callback function
 
    \ingroup setalgo
 
    @sa bit_visitor_callback_type
*/
template<class BV>
void visit_each_bit(const BV&                 bv,
                    void*                     handle_ptr,
                    bit_visitor_callback_type callback_ptr)
{
    // private adaptor for C-style callbacks
    struct callback_adaptor
    {
        callback_adaptor(void* h, bit_visitor_callback_type cb_func)
        : handle_(h), func_(cb_func)
        {}
        
        void add_bits(bm::id_t offset, const unsigned char* bits, unsigned size)
        {
            for (unsigned i = 0; i < size; ++i)
                func_(handle_, offset + bits[i]);
        }
        void add_range(bm::id_t offset, unsigned size)
        {
            for (unsigned i = 0; i < size; ++i)
                func_(handle_, offset + i);
        }
        
        void* handle_;
        bit_visitor_callback_type func_;
    };
    
    callback_adaptor func(handle_ptr, callback_ptr);
    bm::for_each_bit(bv, func);
}



} // bm

#include "bmundef.h"

#endif
