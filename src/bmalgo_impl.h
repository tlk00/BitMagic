#ifndef BMALGO_IMPL__H__INCLUDED__
#define BMALGO_IMPL__H__INCLUDED__
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

/*! \file bmalgo_impl.h
    \brief Algorithms for bvector<>
*/

#ifdef _MSC_VER
#pragma warning( push )
#pragma warning( disable : 4311 4312 4127)
#endif

#include "bmdef.h"
#include "bmutil.h"

namespace bm
{

/*! 
   @defgroup setalgo bvector<> algorithms 
 
   Set algebra algorithms using bit-vectors, arrays. Binary metrics and distances. Random sub-sets.
   @ingroup bvector
 */

/*! 
    @defgroup distance Binary-distance metrics
 
    Distance metrics and algorithms to compute binary distances
    @ingroup setalgo
 */


/*! 
    \brief    Distance metrics codes defined for vectors A and B
    \ingroup  distance
*/
enum distance_metric
{
    COUNT_AND = set_COUNT_AND,          //!< (A & B).count()
    COUNT_XOR = set_COUNT_XOR,          //!< (A ^ B).count()
    COUNT_OR  = set_COUNT_OR,           //!< (A | B).count()
    COUNT_SUB_AB = set_COUNT_SUB_AB,    //!< (A - B).count()
    COUNT_SUB_BA = set_COUNT_SUB_BA,    //!< (B - A).count()
    COUNT_A      = set_COUNT_A,         //!< A.count()
    COUNT_B      = set_COUNT_B          //!< B.count()
};

/**
    Convert set operation into compatible distance metric
    \ingroup  distance
*/
inline
distance_metric operation2metric(set_operation op)
{
    BM_ASSERT(is_const_set_operation(op));
    if (op == set_COUNT) op = set_COUNT_B;
    // distance metric is created as a set operation sub-class
    // simple cast will work to convert
    return (distance_metric) op;
}

/*! 
    \brief Distance metric descriptor, holds metric code and result.
    \sa distance_operation
*/

struct distance_metric_descriptor
{
     distance_metric   metric;
     bm::id_t          result;
     
     distance_metric_descriptor(distance_metric m)
     : metric(m),
       result(0)
    {}
    distance_metric_descriptor()
    : metric(bm::COUNT_XOR),
      result(0)
    {}
    
    /*! 
        \brief Sets metric result to 0
    */
    void reset()
    {
        result = 0;
    }
};



/*!
    \brief Internal function computes different distance metrics.
    \internal 
    \ingroup  distance
     
*/
inline
void combine_count_operation_with_block(const bm::word_t*           blk,
                                        const bm::word_t*           arg_blk,
                                        distance_metric_descriptor* dmit,
                                        distance_metric_descriptor* dmit_end)
                                            
{     
     gap_word_t* g1 = BMGAP_PTR(blk);
     gap_word_t* g2 = BMGAP_PTR(arg_blk);

     unsigned gap = BM_IS_GAP(blk);
     unsigned arg_gap = BM_IS_GAP(arg_blk);
     
     if (gap) // first block GAP-type
     {
         if (arg_gap)  // both blocks GAP-type
         {
             for (distance_metric_descriptor* it = dmit;it < dmit_end; ++it)
             {
                 distance_metric_descriptor& dmd = *it;
                 
                 switch (dmd.metric)
                 {
                 case bm::COUNT_AND:
                     dmd.result += gap_count_and(g1, g2);
                     break;
                 case bm::COUNT_OR:
                     dmd.result += gap_count_or(g1, g2);
                     break;
                 case bm::COUNT_SUB_AB:
                     dmd.result += gap_count_sub(g1, g2);
                     break;
                 case bm::COUNT_SUB_BA:
                     dmd.result += gap_count_sub(g2, g1);
                     break;
                 case bm::COUNT_XOR:
                     dmd.result += gap_count_xor(g1, g2);
                    break;
                 case bm::COUNT_A:
                    dmd.result += gap_bit_count_unr(g1);
                    break;
                 case bm::COUNT_B:
                    dmd.result += gap_bit_count_unr(g2);
                    break;
                 default:
                     BM_ASSERT(0);
                 } // switch
                                     
             } // for it
             
             return;

         }
         else // first block - GAP, argument - BITSET
         {
             for (distance_metric_descriptor* it = dmit;it < dmit_end; ++it)
             {
                 distance_metric_descriptor& dmd = *it;
                 
                 switch (dmd.metric)
                 {
                 case bm::COUNT_AND:
                     if (arg_blk)
                        dmd.result += gap_bitset_and_count(arg_blk, g1);
                     break;
                 case bm::COUNT_OR:
                     if (!arg_blk)
                        dmd.result += gap_bit_count_unr(g1);
                     else
                        dmd.result += gap_bitset_or_count(arg_blk, g1); 
                     break;
                 case bm::COUNT_SUB_AB:
                     {
                     bm::word_t  BM_VECT_ALIGN temp_bit_blk[bm::set_block_size] BM_VECT_ALIGN_ATTR;

                     gap_convert_to_bitset((bm::word_t*) temp_bit_blk, g1);
                     dmd.result += 
                       bit_operation_sub_count((bm::word_t*)temp_bit_blk, arg_blk);
                     }
                     break;
                 case bm::COUNT_SUB_BA:
                     dmd.metric = bm::COUNT_SUB_AB; // recursive call to SUB_AB
                     combine_count_operation_with_block(arg_blk,
                                                        blk,
                                                        it, it+1);
                     dmd.metric = bm::COUNT_SUB_BA; // restore status quo
                     break;
                 case bm::COUNT_XOR:
                     if (!arg_blk)
                        dmd.result += gap_bit_count_unr(g1);
                     else
                        dmd.result += gap_bitset_xor_count(arg_blk, g1);
                     break;
                 case bm::COUNT_A:
                    if (g1)
                        dmd.result += gap_bit_count_unr(g1);
                    break;
                 case bm::COUNT_B:
                    if (arg_blk)
                    {
                        dmd.result += 
                          bit_block_count(arg_blk);
                    }
                    break;
                 default:
                     BM_ASSERT(0);
                 } // switch
                                     
             } // for it
             
             return;
         
         }
     } 
     else // first block is BITSET-type
     {     
         if (arg_gap) // second argument block is GAP-type
         {
             for (distance_metric_descriptor* it = dmit;it < dmit_end; ++it)
             {
                 distance_metric_descriptor& dmd = *it;
                 
                 switch (dmd.metric)
                 {
                 case bm::COUNT_AND:
                     if (blk) 
                        dmd.result += gap_bitset_and_count(blk, g2);                         
                     break;
                 case bm::COUNT_OR:
                     if (!blk)
                        dmd.result += gap_bit_count_unr(g2);
                     else
                        dmd.result += gap_bitset_or_count(blk, g2);
                     break;
                 case bm::COUNT_SUB_AB:
                     if (blk)
                        dmd.result += gap_bitset_sub_count(blk, g2);
                     break;
                 case bm::COUNT_SUB_BA:
                     dmd.metric = bm::COUNT_SUB_AB; // recursive call to SUB_AB
                     combine_count_operation_with_block(arg_blk,
                                                        //arg_gap,
                                                        blk,
                                                        //gap,
                                                        it, it+1);
                     dmd.metric = bm::COUNT_SUB_BA; // restore status quo
                     break;
                 case bm::COUNT_XOR:
                     if (!blk)
                        dmd.result += gap_bit_count_unr(g2);
                     else
                        dmd.result += gap_bitset_xor_count(blk, g2); 
                    break;
                 case bm::COUNT_A:
                    if (blk)
                    {
                        dmd.result += 
                            bit_block_count(blk);
                    }
                    break;
                 case bm::COUNT_B:
                    if (g2)
                        dmd.result += gap_bit_count_unr(g2);
                    break;
                 default:
                     BM_ASSERT(0);
                 } // switch
                                     
             } // for it
             
             return;
         }
     }

     // --------------------------------------------
     //
     // Here we combine two plain bitblocks 

     for (distance_metric_descriptor* it = dmit; it < dmit_end; ++it)
     {
         distance_metric_descriptor& dmd = *it;
        bit_operation_count_func_type gfunc = 
            operation_functions<true>::bit_operation_count(dmd.metric);
        if (gfunc)
        {
            dmd.result += (*gfunc)(blk, arg_blk);
        }
        else
        {
            switch (dmd.metric)
            {
            case bm::COUNT_A:
                if (blk)
                    dmd.result += bm::bit_block_count(blk);
                break;
            case bm::COUNT_B:
                if (arg_blk)
                    dmd.result += bm::bit_block_count(arg_blk);
                break;
            case bm::COUNT_AND:
            case bm::COUNT_XOR:
            case bm::COUNT_OR:
            case bm::COUNT_SUB_AB:
            case bm::COUNT_SUB_BA:
            default:
                BM_ASSERT(0);
            } // switch
        }

     } // for it
}

/*!
\brief Internal function computes AND distance.
\internal 
\ingroup  distance
*/
inline
unsigned combine_count_and_operation_with_block(const bm::word_t* blk,
                                                const bm::word_t* arg_blk)
{
    unsigned gap = BM_IS_GAP(blk);
    unsigned arg_gap = BM_IS_GAP(arg_blk);
    if (gap) // first block GAP-type
    {
        if (arg_gap)  // both blocks GAP-type
        {
            return gap_count_and(BMGAP_PTR(blk), BMGAP_PTR(arg_blk));
        }
        else // first block - GAP, argument - BITSET
        {
            return gap_bitset_and_count(arg_blk, BMGAP_PTR(blk));
        }
    } 
    else // first block is BITSET-type
    {     
        if (arg_gap) // second argument block is GAP-type
        {
            return gap_bitset_and_count(blk, BMGAP_PTR(arg_blk));
        }
    }

    // --------------------------------------------
    // Here we combine two plain bitblocks 

    return bit_operation_and_count(blk, arg_blk);
}


/*!
    \brief Internal function computes different existense of distance metric.
    \internal 
    \ingroup  distance
*/
inline
void combine_any_operation_with_block(const bm::word_t* blk,
                                      unsigned gap,
                                      const bm::word_t* arg_blk,
                                      unsigned arg_gap,
                                      distance_metric_descriptor* dmit,
                                      distance_metric_descriptor* dmit_end)
                                            
{
     gap_word_t* res=0;
     
     gap_word_t* g1 = BMGAP_PTR(blk);
     gap_word_t* g2 = BMGAP_PTR(arg_blk);
     
     if (gap) // first block GAP-type
     {
         if (arg_gap)  // both blocks GAP-type
         {
             gap_word_t tmp_buf[bm::gap_max_buff_len * 3]; // temporary result
             
             for (distance_metric_descriptor* it = dmit;it < dmit_end; ++it)
             {
                 distance_metric_descriptor& dmd = *it;
                 if (dmd.result)
                 {
                     continue;
                 }
                 res = 0;
                 unsigned dsize = 0;
                 switch (dmd.metric)
                 {
                 case bm::COUNT_AND:
                     dmd.result += gap_operation_any_and(g1, g2);
                     break;
                 case bm::COUNT_OR:
                     res = gap_operation_or(g1, g2, tmp_buf, dsize);
                     break;
                 case bm::COUNT_SUB_AB:
                     dmd.result += gap_operation_any_sub(g1, g2); 
                     break;
                 case bm::COUNT_SUB_BA:
                     dmd.result += gap_operation_any_sub(g2, g1); 
                     break;
                 case bm::COUNT_XOR:
                    dmd.result += gap_operation_any_xor(g1, g2); 
                    break;
                 case bm::COUNT_A:
                    res = g1;
                    break;
                 case bm::COUNT_B:
                    res = g2;
                    break;
                 default:
                     BM_ASSERT(0);
                 } // switch
                if (res)
                    dmd.result += !gap_is_all_zero(res);
                                     
             } // for it
             
             return;

         }
         else // first block - GAP, argument - BITSET
         {
             for (distance_metric_descriptor* it = dmit;it < dmit_end; ++it)
             {
                 distance_metric_descriptor& dmd = *it;
                 if (dmd.result)
                 {
                     continue;
                 }
                 
                 switch (dmd.metric)
                 {
                 case bm::COUNT_AND:
                     if (arg_blk)
                        dmd.result += gap_bitset_and_any(arg_blk, g1);
                     break;
                 case bm::COUNT_OR:
                     if (!arg_blk)
                        dmd.result += !gap_is_all_zero(g1);
                     else
                        dmd.result += gap_bitset_or_any(arg_blk, g1); 
                     break;
                 case bm::COUNT_SUB_AB:
                     {
                     bm::word_t  BM_VECT_ALIGN temp_blk[bm::set_block_size] BM_VECT_ALIGN_ATTR;
                     gap_convert_to_bitset((bm::word_t*) temp_blk, g1);
                     dmd.result += 
                       bit_operation_sub_any((bm::word_t*)temp_blk, arg_blk);
                     }
                     break;
                 case bm::COUNT_SUB_BA:
                     dmd.metric = bm::COUNT_SUB_AB; // recursive call to SUB_AB
                     combine_any_operation_with_block(arg_blk,
                                                      arg_gap,
                                                      blk,
                                                      gap,
                                                      it, it+1);
                     dmd.metric = bm::COUNT_SUB_BA; // restore status quo
                     break;
                 case bm::COUNT_XOR:
                     if (!arg_blk)
                        dmd.result += !gap_is_all_zero(g1);
                     else
                        dmd.result += gap_bitset_xor_any(arg_blk, g1);
                     break;
                 case bm::COUNT_A:
                    if (g1)
                        dmd.result += !gap_is_all_zero(g1);
                    break;
                 case bm::COUNT_B:
                    if (arg_blk)
                    {
                        dmd.result += 
                          !bit_is_all_zero(arg_blk);
                    }
                    break;
                 default:
                     BM_ASSERT(0);
                 } // switch
                                     
             } // for it
             
             return;
         
         }
     } 
     else // first block is BITSET-type
     {     
         if (arg_gap) // second argument block is GAP-type
         {
             for (distance_metric_descriptor* it = dmit;it < dmit_end; ++it)
             {
                 distance_metric_descriptor& dmd = *it;
                 if (dmd.result)
                 {
                     continue;
                 }
                 
                 switch (dmd.metric)
                 {
                 case bm::COUNT_AND:
                     if (blk) 
                        dmd.result += gap_bitset_and_any(blk, g2);                         
                     break;
                 case bm::COUNT_OR:
                     if (!blk)
                        dmd.result += !gap_is_all_zero(g2);
                     else
                        dmd.result += gap_bitset_or_any(blk, g2);
                     break;
                 case bm::COUNT_SUB_AB:
                     if (blk)
                        dmd.result += gap_bitset_sub_any(blk, g2);
                     break;
                 case bm::COUNT_SUB_BA:
                     dmd.metric = bm::COUNT_SUB_AB; // recursive call to SUB_AB
                     combine_any_operation_with_block(arg_blk,
                                                      arg_gap,
                                                      blk,
                                                      gap,
                                                      it, it+1);
                     dmd.metric = bm::COUNT_SUB_BA; // restore status quo
                     break;
                 case bm::COUNT_XOR:
                     if (!blk)
                        dmd.result += !gap_is_all_zero(g2);
                     else
                        dmd.result += gap_bitset_xor_any(blk, g2); 
                    break;
                 case bm::COUNT_A:
                    if (blk)
                    {
                        dmd.result+=
                            !bm::bit_is_all_zero(blk);
                    }
                    break;
                 case bm::COUNT_B:
                    if (g2)
                        dmd.result += !gap_is_all_zero(g2);
                    break;
                 default:
                     BM_ASSERT(0);
                 } // switch
                                     
             } // for it
             
             return;
         }
     }

     // --------------------------------------------
     //
     // Here we combine two plain bitblocks 

     for (distance_metric_descriptor* it = dmit; it < dmit_end; ++it)
     {
        distance_metric_descriptor& dmd = *it;
        if (dmd.result)
        {
            continue;
        }

        switch (dmd.metric)
        {
        case bm::COUNT_AND:
            dmd.result += 
            bit_operation_and_any(blk, arg_blk);
            break;
        case bm::COUNT_OR:
            dmd.result += 
            bit_operation_or_any(blk, arg_blk);
            break;
        case bm::COUNT_SUB_AB:
            dmd.result += 
            bit_operation_sub_any(blk, arg_blk);
            break;
        case bm::COUNT_SUB_BA:
            dmd.result += 
            bit_operation_sub_any(arg_blk, blk);
            break;
        case bm::COUNT_XOR:
            dmd.result += 
            bit_operation_xor_any(blk, arg_blk);
            break;
        case bm::COUNT_A:
            if (blk)
                dmd.result += !bit_is_all_zero(blk);
            break;
        case bm::COUNT_B:
            if (arg_blk)
                dmd.result += !bit_is_all_zero(arg_blk);
            break;
        default:
            BM_ASSERT(0);
        } // switch

     } // for it
}



/*!
    Convenience internal function to compute combine count for one metric
    \internal
    \ingroup  distance
*/
inline
unsigned combine_count_operation_with_block(const bm::word_t* blk,
                                            const bm::word_t* arg_blk,
                                            distance_metric metric)
{
    distance_metric_descriptor dmd(metric);
    combine_count_operation_with_block(blk, //gap, 
                                       arg_blk, //arg_gap, 
                                       &dmd, &dmd+1);
    return dmd.result;
}


/*!
    Convenience internal function to compute combine any for one metric
    \internal
    \ingroup  distance
*/
inline
unsigned combine_any_operation_with_block(const bm::word_t* blk,
                                          unsigned gap,
                                          const bm::word_t* arg_blk,
                                          unsigned arg_gap,
                                          distance_metric metric)
{
    distance_metric_descriptor dmd(metric);
    combine_any_operation_with_block(blk, gap, 
                                     arg_blk, arg_gap, 
                                     &dmd, &dmd+1);
    return dmd.result;
}

/*!
    \brief Staging function for distance operation

    \return temp block allocated (or NULL)

    \internal
*/
inline
void distance_stage(const distance_metric_descriptor* dmit,
                    const distance_metric_descriptor* dmit_end,
                    bool*                             is_all_and)
{
    for (const distance_metric_descriptor* it = dmit; it < dmit_end; ++it)
    {
        if (it->metric != bm::COUNT_AND)
        {
            *is_all_and = false;
        } 
    } // for
}

/*!
    \brief Distance computing template function.

    Function receives two bitvectors and an array of distance metrics
    (metrics pipeline). Function computes all metrics saves result into
    corresponding pipeline results (distance_metric_descriptor::result)
    An important detail is that function reuses metric descriptors, 
    incrementing received values. It allows you to accumulate results 
    from different calls in the pipeline.
    
    \param bv1      - argument bitvector 1 (A)
    \param bv2      - argument bitvector 2 (B)
    \param dmit     - pointer to first element of metric descriptors array
                      Input-Output parameter, receives metric code as input,
                      computation is added to "result" field
    \param dmit_end - pointer to (last+1) element of metric descriptors array
    \ingroup  distance
    
*/
template<class BV>
void distance_operation(const BV& bv1, 
                        const BV& bv2, 
                        distance_metric_descriptor* dmit,
                        distance_metric_descriptor* dmit_end)
{
    const typename BV::blocks_manager_type& bman1 = bv1.get_blocks_manager();
    const typename BV::blocks_manager_type& bman2 = bv2.get_blocks_manager();

    bool is_all_and = true; // flag is distance operation is just COUNT_AND
    distance_stage(dmit, dmit_end, &is_all_and);

    bm::word_t*** blk_root = bman1.top_blocks_root();
    unsigned block_idx = 0;
    unsigned i, j;
    
    const bm::word_t* blk;
    const bm::word_t* arg_blk;

    BM_SET_MMX_GUARD

    unsigned top_block_size = bman1.top_block_size();
    unsigned ebs2 = bman2.top_block_size();
    if (ebs2 > top_block_size)
        top_block_size = ebs2;

    for (i = 0; i < top_block_size; ++i)
    {
        bm::word_t** blk_blk = blk_root ? blk_root[i] : 0;

        if (blk_blk == 0) // not allocated
        {
            // AND operation requested - we can skip this portion here 
            if (is_all_and)
            {
                block_idx += bm::set_array_size;
                continue;
            }
            const bm::word_t* const* bvbb = bman2.get_topblock(i);
            if (bvbb == 0) 
            {
                block_idx += bm::set_array_size;
                continue;
            }

            blk = 0;
            for (j = 0; j < bm::set_array_size; ++j,++block_idx)
            {                
                arg_blk = bman2.get_block(i, j);
                if (!arg_blk) 
                    continue;
                combine_count_operation_with_block(blk, 
                                                   arg_blk, 
                                                   dmit, dmit_end);
            } // for j
            continue;
        }

        for (j = 0; j < bm::set_array_size; ++j, ++block_idx)
        {
            blk = bman1.get_block(i, j);
            if (blk == 0 && is_all_and)
                continue;

            arg_blk = bman2.get_block(i, j);

            if (!blk & !arg_blk)
                continue;
                
            combine_count_operation_with_block(blk, 
                                               arg_blk, 
                                               dmit, dmit_end);
        } // for j

    } // for i
}


/*!
\brief Distance AND computing template function.


\param bv1      - argument bitvector 1 (A)
\param bv2      - argument bitvector 2 (B)
\ingroup  distance

*/
template<class BV>
unsigned distance_and_operation(const BV& bv1, 
                                const BV& bv2)
{
    const typename BV::blocks_manager_type& bman1 = bv1.get_blocks_manager();
    const typename BV::blocks_manager_type& bman2 = bv2.get_blocks_manager();
    
    if (!bman1.is_init() || !bman2.is_init())
        return 0;

    bm::word_t*** blk_root     = bman1.top_blocks_root();
    bm::word_t*** blk_root_arg = bman2.top_blocks_root();
    unsigned count = 0;

    unsigned top_block_size =
        bm::min_value(bman1.top_block_size(),bman2.top_block_size());

    for (unsigned i = 0; i < top_block_size; ++i)
    {
        bm::word_t** blk_blk;
        bm::word_t** blk_blk_arg;
        if ((blk_blk = blk_root[i]) == 0 || (blk_blk_arg= blk_root_arg[i]) == 0)
        {
            continue;
        }
        for (unsigned j = 0; j < bm::set_array_size; j+=4)
        {
            (blk_blk[j] && blk_blk_arg[j]) ? 
                count += combine_count_and_operation_with_block(BLOCK_ADDR_SAN(blk_blk[j]), BLOCK_ADDR_SAN(blk_blk_arg[j]))
                :0;
            (blk_blk[j+1] && blk_blk_arg[j+1]) ? 
                count += combine_count_and_operation_with_block(BLOCK_ADDR_SAN(blk_blk[j+1]), BLOCK_ADDR_SAN(blk_blk_arg[j+1]))
                :0;
            (blk_blk[j+2] && blk_blk_arg[j+2]) ? 
                count += combine_count_and_operation_with_block(BLOCK_ADDR_SAN(blk_blk[j+2]), BLOCK_ADDR_SAN(blk_blk_arg[j+2]))
                :0;
            (blk_blk[j+3] && blk_blk_arg[j+3]) ? 
                count += combine_count_and_operation_with_block(BLOCK_ADDR_SAN(blk_blk[j+3]), BLOCK_ADDR_SAN(blk_blk_arg[j+3]))
                :0;
        } // for j

    } // for i
    return count;
}




/*!
    \brief Distance screening template function.

    Function receives two bitvectors and an array of distance metrics
    (metrics pipeline). Function computes possybility of a metric(any bit) 
    saves result into corresponding pipeline results 
    (distance_metric_descriptor::result)
    An important detail is that function reuses metric descriptors, 
    incrementing received values. It allows you to accumulate results 
    from different calls in the pipeline.
    
    \param bv1      - argument bitvector 1 (A)
    \param bv2      - argument bitvector 2 (B)
    \param dmit     - pointer to first element of metric descriptors array
                      Input-Output parameter, receives metric code as input,
                      computation is added to "result" field
    \param dmit_end - pointer to (last+1) element of metric descriptors array
    \ingroup  distance
*/
template<class BV>
void distance_operation_any(const BV& bv1, 
                            const BV& bv2, 
                            distance_metric_descriptor* dmit,
                            distance_metric_descriptor* dmit_end)
{
    const typename BV::blocks_manager_type& bman1 = bv1.get_blocks_manager();
    const typename BV::blocks_manager_type& bman2 = bv2.get_blocks_manager();
    
    bool is_all_and = true; // flag is distance operation is just COUNT_AND
    distance_stage(dmit, dmit_end, &is_all_and);
  
    bm::word_t*** blk_root = bman1.top_blocks_root();
    unsigned block_idx = 0;
    unsigned i, j;
    
    const bm::word_t* blk;
    const bm::word_t* arg_blk;
    bool  blk_gap;
    bool  arg_gap;

    unsigned top_block_size = bman1.top_block_size();
    unsigned ebs2 = bman2.top_block_size();
    if (ebs2 > top_block_size)
        top_block_size = ebs2;

    for (i = 0; i < top_block_size; ++i)
    {
        bm::word_t** blk_blk = blk_root ? blk_root[i] : 0;

        if (blk_blk == 0) // not allocated
        {
            // AND operation requested - we can skip this portion here 
            if (is_all_and)
            {
                block_idx += bm::set_array_size;
                continue;
            }

            const bm::word_t* const* bvbb = bman2.get_topblock(i);
            if (bvbb == 0) 
            {
                block_idx += bm::set_array_size;
                continue;
            }

            blk = 0;
            blk_gap = false;

            for (j = 0; j < bm::set_array_size; ++j,++block_idx)
            {                
                arg_blk = bman2.get_block(i, j);
                if (!arg_blk) 
                    continue;
                arg_gap = BM_IS_GAP(arg_blk);
                
                combine_any_operation_with_block(blk, blk_gap,
                                                 arg_blk, arg_gap,
                                                 dmit, dmit_end);

                // check if all distance requests alredy resolved
                bool all_resolved = false;
                distance_metric_descriptor* it=dmit;
                do
                {
                    if (!it->result)
                    {
                        all_resolved = false;
                        break;
                    }
                    ++it;
                } while (it < dmit_end);
                if (all_resolved)
                    return;
            } // for j

            continue;
        }

        for (j = 0; j < bm::set_array_size; ++j, ++block_idx)
        {
            blk = bman1.get_block(i, j);
            //BLOCK_ADDR_SAN(blk_blk[j]);
            if (blk == 0 && is_all_and)
                continue;

            arg_blk = bman2.get_block(i, j);

            if (blk == 0 && arg_blk == 0)
                continue;
                
            arg_gap = BM_IS_GAP(arg_blk);
            blk_gap = BM_IS_GAP(blk);
            
            combine_any_operation_with_block(blk, blk_gap,
                                             arg_blk, arg_gap,
                                             dmit, dmit_end);
            
            // check if all distance requests alredy resolved
            bool all_resolved = false;
            distance_metric_descriptor* it=dmit;
            do
            {
                if (!it->result)
                {
                    all_resolved = false;
                    break;
                }
                ++it;
            } while (it < dmit_end);
            if (all_resolved)
                return;

        } // for j

    } // for i
}



/*!
   \brief Computes bitcount of AND operation of two bitsets
   \param bv1 - Argument bit-vector.
   \param bv2 - Argument bit-vector.
   \return bitcount of the result
   \ingroup  distance
*/
template<class BV>
bm::id_t count_and(const BV& bv1, const BV& bv2)
{
    return distance_and_operation(bv1, bv2);
}

/*!
   \brief Computes if there is any bit in AND operation of two bitsets
   \param bv1 - Argument bit-vector.
   \param bv2 - Argument bit-vector.
   \return non zero value if there is any bit
   \ingroup  distance
*/
template<class BV>
bm::id_t any_and(const BV& bv1, const BV& bv2)
{
    distance_metric_descriptor dmd(bm::COUNT_AND);
    
    distance_operation_any(bv1, bv2, &dmd, &dmd+1);
    return dmd.result;
}



/*!
   \brief Computes bitcount of XOR operation of two bitsets
   \param bv1 - Argument bit-vector.
   \param bv2 - Argument bit-vector.
   \return bitcount of the result
   \ingroup  distance
*/
template<class BV>
bm::id_t count_xor(const BV& bv1, const BV& bv2)
{
    distance_metric_descriptor dmd(bm::COUNT_XOR);
    
    distance_operation(bv1, bv2, &dmd, &dmd+1);
    return dmd.result;
}

/*!
   \brief Computes if there is any bit in XOR operation of two bitsets
   \param bv1 - Argument bit-vector.
   \param bv2 - Argument bit-vector.
   \return non zero value if there is any bit
   \ingroup  distance
*/
template<class BV>
bm::id_t any_xor(const BV& bv1, const BV& bv2)
{
    distance_metric_descriptor dmd(bm::COUNT_XOR);
    
    distance_operation_any(bv1, bv2, &dmd, &dmd+1);
    return dmd.result;
}



/*!
   \brief Computes bitcount of SUB operation of two bitsets
   \param bv1 - Argument bit-vector.
   \param bv2 - Argument bit-vector.
   \return bitcount of the result
   \ingroup  distance
*/
template<class BV>
bm::id_t count_sub(const BV& bv1, const BV& bv2)
{
    distance_metric_descriptor dmd(bm::COUNT_SUB_AB);
    
    distance_operation(bv1, bv2, &dmd, &dmd+1);
    return dmd.result;
}


/*!
   \brief Computes if there is any bit in SUB operation of two bitsets
   \param bv1 - Argument bit-vector.
   \param bv2 - Argument bit-vector.
   \return non zero value if there is any bit
   \ingroup  distance
*/
template<class BV>
bm::id_t any_sub(const BV& bv1, const BV& bv2)
{
    distance_metric_descriptor dmd(bm::COUNT_SUB_AB);
    
    distance_operation_any(bv1, bv2, &dmd, &dmd+1);
    return dmd.result;
}


/*!    
   \brief Computes bitcount of OR operation of two bitsets
   \param bv1 - Argument bit-vector.
   \param bv2 - Argument bit-vector.
   \return bitcount of the result
   \ingroup  distance
*/
template<class BV>
bm::id_t count_or(const BV& bv1, const BV& bv2)
{
    distance_metric_descriptor dmd(bm::COUNT_OR);
    
    distance_operation(bv1, bv2, &dmd, &dmd+1);
    return dmd.result;
}

/*!
   \brief Computes if there is any bit in OR operation of two bitsets
   \param bv1 - Argument bit-vector.
   \param bv2 - Argument bit-vector.
   \return non zero value if there is any bit
   \ingroup  distance
*/
template<class BV>
bm::id_t any_or(const BV& bv1, const BV& bv2)
{
    distance_metric_descriptor dmd(bm::COUNT_OR);
    
    distance_operation_any(bv1, bv2, &dmd, &dmd+1);
    return dmd.result;
}


/**
    \brief Internal algorithms scans the input for the block range limit
    \internal
*/
template<class It>
It block_range_scan(It  first, It last, unsigned nblock, unsigned* max_id)
{
    It right;
    for (right = first; right != last; ++right)
    {
        unsigned v = unsigned(*right);
        BM_ASSERT(v < bm::id_max);
        if (v >= *max_id)
            *max_id = v;
        unsigned nb = v >> bm::set_block_shift;
        if (nb != nblock)
            break;
    }
    return right;
}

/**
    \brief OR Combine bitvector and the iterable sequence 

    Algorithm combines bvector with sequence of integers 
    (represents another set). When the agrument set is sorted 
    this method can give serious increase in performance.
    
    \param bv     - destination bitvector
    \param first  - first element of the iterated sequence
    \param last   - last element of the iterated sequence
    
    \ingroup setalgo
*/
template<class BV, class It>
void combine_or(BV& bv, It  first, It last)
{
    typename BV::blocks_manager_type& bman = bv.get_blocks_manager();
    if (!bman.is_init())
        bman.init_tree();
    
    unsigned max_id = 0;

    while (first < last)
    {
        unsigned nblock = unsigned((*first) >> bm::set_block_shift);     
        It right = bm::block_range_scan(first, last, nblock, &max_id);

        if (max_id >= bv.size())
        {
            BM_ASSERT(max_id < bm::id_max);
            bv.resize(max_id + 1);
        }

        // now we have one in-block array of bits to set
        
        label1:
        
        int block_type;
        bm::word_t* blk =
            bman.check_allocate_block(nblock, 
                                      true, 
                                      bv.get_new_blocks_strat(), 
                                      &block_type);
        if (!blk) 
            continue;
                        
        if (block_type == 1) // gap
        {            
            bm::gap_word_t* gap_blk = BMGAP_PTR(blk);
            unsigned threshold = bm::gap_limit(gap_blk, bman.glen());
            
            for (; first < right; ++first)
            {
                unsigned is_set;
                unsigned nbit   = (*first) & bm::set_block_mask; 
                
                unsigned new_block_len =
                    bm::gap_set_value(true, gap_blk, nbit, &is_set);
                if (new_block_len > threshold) 
                {
                    bman.extend_gap_block(nblock, gap_blk);
                    ++first;
                    goto label1;  // block pointer changed - goto reset
                }
            }
        }
        else // bit block
        {
            for (; first < right; ++first)
            {
                unsigned nbit   = unsigned(*first & bm::set_block_mask); 
                unsigned nword  = unsigned(nbit >> bm::set_word_shift); 
                nbit &= bm::set_word_mask;
                blk[nword] |= (1u << nbit);
            } // for
        }
    } // while
    
    bv.forget_count();
}


/**
    \brief XOR Combine bitvector and the iterable sequence 

    Algorithm combines bvector with sequence of integers 
    (represents another set). When the agrument set is sorted 
    this method can give serious increase in performance.
    
    \param bv     - destination bitvector
    \param first  - first element of the iterated sequence
    \param last   - last element of the iterated sequence
    
    \ingroup setalgo
*/
template<class BV, class It>
void combine_xor(BV& bv, It  first, It last)
{
    typename BV::blocks_manager_type& bman = bv.get_blocks_manager();
    if (!bman.is_init())
        bman.init_tree();
    
    unsigned max_id = 0;

    while (first < last)
    {
        unsigned nblock = unsigned((*first) >> bm::set_block_shift);     
        It right = block_range_scan(first, last, nblock, &max_id);

        if (max_id >= bv.size())
        {
            BM_ASSERT(max_id < bm::id_max);
            bv.resize(max_id + 1);
        }

        // now we have one in-block array of bits to set
        
        label1:
        
        int block_type;
        bm::word_t* blk =
            bman.check_allocate_block(nblock, 
                                      true, 
                                      bv.get_new_blocks_strat(), 
                                      &block_type,
                                      false /* no null return */);
        BM_ASSERT(blk); 
                        
        if (block_type == 1) // gap
        {
            bm::gap_word_t* gap_blk = BMGAP_PTR(blk);
            unsigned threshold = bm::gap_limit(gap_blk, bman.glen());
            
            for (; first < right; ++first)
            {
                unsigned is_set;
                unsigned nbit   = (*first) & bm::set_block_mask; 
                
                is_set = bm::gap_test_unr(gap_blk, nbit);
                BM_ASSERT(is_set <= 1);
                is_set ^= 1; 
                
                unsigned new_block_len =
                    gap_set_value(is_set, gap_blk, nbit, &is_set);
                if (new_block_len > threshold) 
                {
                    bman.extend_gap_block(nblock, gap_blk);
                    ++first;
                    goto label1;  // block pointer changed - goto reset
                }
            }
        }
        else // bit block
        {
            for (; first < right; ++first)
            {
                unsigned nbit   = unsigned(*first & bm::set_block_mask); 
                unsigned nword  = unsigned(nbit >> bm::set_word_shift); 
                nbit &= bm::set_word_mask;
                blk[nword] ^= (1u << nbit);
            } // for
        }
    } // while
    
    bv.forget_count();
}



/**
    \brief SUB Combine bitvector and the iterable sequence 

    Algorithm combines bvector with sequence of integers 
    (represents another set). When the agrument set is sorted 
    this method can give serious increase in performance.
    
    \param bv     - destination bitvector
    \param first  - first element of the iterated sequence
    \param last   - last element of the iterated sequence
    
    \ingroup setalgo
*/
template<class BV, class It>
void combine_sub(BV& bv, It  first, It last)
{
    typename BV::blocks_manager_type& bman = bv.get_blocks_manager();
    if (!bman.is_init())
        bman.init_tree();
    
    unsigned max_id = 0;

    while (first < last)
    {
        unsigned nblock = unsigned((*first) >> bm::set_block_shift);     
        It right = block_range_scan(first, last, nblock, &max_id);

        if (max_id >= bv.size())
        {
            BM_ASSERT(max_id < bm::id_max);
            bv.resize(max_id + 1);
        }

        // now we have one in-block array of bits to set
        
        label1:
        
        int block_type;
        bm::word_t* blk =
            bman.check_allocate_block(nblock, 
                                      false, 
                                      bv.get_new_blocks_strat(), 
                                      &block_type);

        if (!blk)
            continue;
                        
        if (block_type == 1) // gap
        {
            bm::gap_word_t* gap_blk = BMGAP_PTR(blk);
            unsigned threshold = bm::gap_limit(gap_blk, bman.glen());
            
            for (; first < right; ++first)
            {
                unsigned is_set;
                unsigned nbit   = (*first) & bm::set_block_mask; 
                
                is_set = bm::gap_test_unr(gap_blk, nbit);
                if (!is_set)
                    continue;
                
                unsigned new_block_len =
                    gap_set_value(false, gap_blk, nbit, &is_set);
                if (new_block_len > threshold) 
                {
                    bman.extend_gap_block(nblock, gap_blk);
                    ++first;
                    goto label1;  // block pointer changed - goto reset
                }
            }
        }
        else // bit block
        {
            for (; first < right; ++first)
            {
                unsigned nbit   = unsigned(*first & bm::set_block_mask); 
                unsigned nword  = unsigned(nbit >> bm::set_word_shift); 
                nbit &= bm::set_word_mask;
                blk[nword] &= ~((bm::word_t)1 << nbit);
            } // for
        }
    } // while
    
    bv.forget_count();
}

/**
    \brief AND Combine bitvector and the iterable sequence 

    Algorithm combines bvector with sorted sequence of integers 
    (represents another set).
    
    \param bv     - destination bitvector
    \param first  - first element of the iterated sequence
    \param last   - last element of the iterated sequence
    
    \ingroup setalgo
*/
template<class BV, class It>
void combine_and_sorted(BV& bv, It  first, It last)
{
    bm::id_t prev = 0;
    for ( ;first < last; ++first)
    {
        bm::id_t id = *first;
        BM_ASSERT(id >= prev); // make sure it's sorted
        bv.set_bit_and(id, true);
        if (++prev < id) 
        {
            bv.set_range(prev, id-1, false);
        }
        prev = id;
    }
}


/**
    \brief AND Combine bitvector and the iterable sequence 

    Algorithm combines bvector with sequence of integers 
    (represents another set). When the agrument set is sorted 
    this method can give serious increase in performance.
    
    \param bv     - destination bitvector
    \param first  - first element of the iterated sequence
    \param last   - last element of the iterated sequence
    
    \ingroup setalgo
    \sa combine_and_sorted
*/
template<class BV, class It>
void combine_and(BV& bv, It  first, It last)
{
    BV bv_tmp;
    combine_or(bv_tmp, first, last);
    bv &= bv_tmp;
}

/*!
    \brief Compute number of bit intervals (GAPs) in the bitvector
    
    Algorithm traverses bit vector and count number of uninterrupted
    intervals of 1s and 0s.
    <pre>
    For example: 
      empty vector   - 1 interval
      00001111100000 - gives us 3 intervals
      10001111100000 - 4 intervals
      00000000000000 - 1 interval
      11111111111111 - 1 interval
    </pre>
    \return Number of intervals
    \ingroup setalgo
*/
template<class BV>
bm::id_t count_intervals(const BV& bv)
{
    const typename BV::blocks_manager_type& bman = bv.get_blocks_manager();
    
    if (!bman.is_init())
        return 1;

    bm::word_t*** blk_root = bman.top_blocks_root();
    typename BV::blocks_manager_type::block_count_change_func func(bman);
    for_each_block(blk_root, bman.top_block_size(), func);

    return func.count();        
}

/*!
    \brief Export bitset from an array of binary data representing
    the bit vector. 

    Input array can be array of ints, chars or other native C types.
    Method works with iterators, which makes it compatible with 
    STL containers and C arrays

    \param bv     - destination bitvector
    \param first  - first element of the iterated sequence
    \param last   - last element of the iterated sequence

    \ingroup setalgo
*/
template<typename BV, class It>
void export_array(BV& bv, It first, It last)
{
    typename BV::blocks_manager_type& bman = bv.get_blocks_manager();
    if (!bman.is_init())
        bman.init_tree();
    
    unsigned inp_word_size = sizeof(*first);
    size_t array_size = size_t(last - first);
    size_t bit_cnt = array_size * inp_word_size * 8;
    int block_type;
    bm::word_t tmp;
    bm::word_t b1, b2, b3, b4;

    if (bit_cnt >= bv.size())
        bv.resize((bm::id_t)bit_cnt + 1);
    else 
        bv.set_range((bm::id_t)bit_cnt, bv.size() - 1, false);
    
    switch (inp_word_size)
    {
    case 1:
        {
            size_t word_cnt = array_size / 4;
            for (unsigned i = 0; i < bm::set_total_blocks; ++i)
            {
                bm::word_t* blk =
                    bman.check_allocate_block(i, 
                                              false, 
                                              BM_BIT, 
                                              &block_type,
                                              false /*no NULL ret*/);
                if (block_type == 1) // gap
                {
                    blk = bman.deoptimize_block(i);
                }
                
                bm::word_t* wrd_ptr = blk;
                if (word_cnt >= bm::set_block_size) {
                    bm::word_t* wrd_end = blk + bm::set_block_size;
                    do {
                        b1 = bm::word_t(*first++); b2 = bm::word_t(*first++);
                        b3 = bm::word_t(*first++); b4 = bm::word_t(*first++);
                        tmp = b1 | (b2 << 8u) | (b3 << 16u) | (b4 << 24u);
                        *wrd_ptr++ = tmp;
                    } while (wrd_ptr < wrd_end);
                    word_cnt -= bm::set_block_size;
                } 
                else 
                {
                    size_t to_convert = size_t(last - first);
                    for (size_t j = 0; j < to_convert / 4; ++j)
                    {
                        b1 = bm::word_t(*first++); b2 = bm::word_t(*first++);
                        b3 = bm::word_t(*first++); b4 = bm::word_t(*first++);
                        tmp = b1 | (b2 << 8u) | (b3 << 16u) | (b4 << 24u);
                        *wrd_ptr++ = tmp;
                    }
                    while (1)
                    {
                        if (first == last) break;
                        *wrd_ptr = bm::word_t(*first++);
                        if (first == last) break;
                        *wrd_ptr |= bm::word_t(*first++) << 8u;
                        if (first == last) break;
                        *wrd_ptr |= bm::word_t(*first++) << 16u;
                        if (first == last) break;
                        *wrd_ptr |= bm::word_t(*first++) << 24u;
                        ++wrd_ptr;
                    }
                }
                if (first == last) break;
            } // for
        }
        break;
    case 2:
        {
            size_t word_cnt = array_size / 2;
            for (unsigned i = 0; i < bm::set_total_blocks; ++i)
            {
                bm::word_t* blk =
                    bman.check_allocate_block(i, 
                                              false, 
                                              BM_BIT, 
                                              &block_type,
                                              false /*no NULL ret*/);
                if (block_type) // gap
                    blk = bman.deoptimize_block(i);
                
                bm::word_t* wrd_ptr = blk;
                if (word_cnt >= bm::set_block_size) {
                    bm::word_t* wrd_end = blk + bm::set_block_size;
                    do {
                        b1 = bm::word_t(*first++); b2 = bm::word_t(*first++);
                        tmp = b1 | (b2 << 16);
                        *wrd_ptr++ = tmp;
                    } while (wrd_ptr < wrd_end);
                    word_cnt -= bm::set_block_size;
                } 
                else 
                {
                    size_t to_convert = size_t(last - first);
                    for (unsigned j = 0; j < to_convert / 2; ++j)
                    {
                        b1 = bm::word_t(*first++); b2 = bm::word_t(*first++);
                        tmp = b1 | (b2 << 16u);
                        *wrd_ptr++ = tmp;
                    }
                    while (1)
                    {
                        if (first == last) break;
                        *wrd_ptr = bm::word_t(*first++);
                        if (first == last) break;
                        *wrd_ptr |= bm::word_t((*first++) << 16u);
                        ++wrd_ptr;
                    }
                }
                if (first == last) break;
            } // for
        }
        break;
    case 4:
        {
            size_t word_cnt = array_size;
            for (unsigned i = 0; i < bm::set_total_blocks; ++i)
            {
                bm::word_t* blk =
                    bman.check_allocate_block(i, 
                                              false, 
                                              BM_BIT, 
                                              &block_type,
                                              false /*no NULL ret*/);
                if (block_type == 1) // gap
                    blk = bman.deoptimize_block(i);
                
                bm::word_t* wrd_ptr = blk;
                if (word_cnt >= bm::set_block_size) {
                    bm::word_t* wrd_end = blk + bm::set_block_size;
                    do
                    {
                        *wrd_ptr++ = bm::word_t(*first++);
                    } while (wrd_ptr < wrd_end);
                    word_cnt -= bm::set_block_size;
                } 
                else 
                {
                    while (1)
                    {
                        if (first == last) break;
                        *wrd_ptr = bm::word_t(*first++);
                        ++wrd_ptr;
                    }
                }
                if (first == last) break;
            } // for
        }
        break;
    default:
        BM_ASSERT(0);
    } // switch

}


/*!
   \brief for-each visitor, calls a special visitor functor for each 1 bit group
 
   \param block - bit block buffer pointer
   \param offset - global block offset (number of bits)
   \param bit_functor - functor must support .add_bits(offset, bits_ptr, size)
 
   @ingroup bitfunc
   @internal
*/
template<class Func>
void for_each_bit_blk(const bm::word_t* block, bm::id_t offset,
                      Func&  bit_functor)
{
    unsigned char bits[bm::set_bitscan_wave_size*32];

    unsigned offs = offset;
    unsigned cnt;
    const word_t* block_end = block + bm::set_block_size;
    do
    {
        cnt = bm::bitscan_wave(block, bits);
        if (cnt)
            bit_functor.add_bits(offs, bits, cnt);
        offs += bm::set_bitscan_wave_size * 32;
        block += bm::set_bitscan_wave_size;
    } while (block < block_end);
}


/*!
   \brief for-each visitor, calls a special visitor functor for each 1 bit range
 
   \param block - bit block buffer pointer
   \param offset - global block offset (number of bits)
   \param bit_functor - functor must support .add_range(offset, bits_ptr, size)
 
   @ingroup gapfunc
   @internal
*/
template<typename T, typename Func>
void for_each_gap_blk(const T* buf, bm::id_t offset,
                      Func&  bit_functor)
{
    const T* pcurr = buf + 1;
    const T* pend = buf + (*buf >> 3);

    if (*buf & 1)
    {
        bit_functor.add_range(offset, *pcurr + 1);
        ++pcurr;
    }
    for (++pcurr; pcurr <= pend; pcurr += 2)
    {
        T prev = *(pcurr-1);
        bit_functor.add_range(offset + prev + 1, *pcurr - prev);
    }
}


} // namespace bm

#ifdef _MSC_VER
#pragma warning( pop )
#endif

#endif
