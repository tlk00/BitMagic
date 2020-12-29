#ifndef BMXORFUNC__H__INCLUDED__
#define BMXORFUNC__H__INCLUDED__
/*
Copyright(c) 2002-2019 Anatoliy Kuznetsov(anatoliy_kuznetsov at yahoo.com)

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

/*! \file bmxor.h
    \brief Functions and utilities for XOR filters (internal)
*/

#include "bmdef.h"
#include "bmutil.h"


namespace bm
{

/**
    XOR complementarity type between 2 blocks
    @internal
 */
enum xor_complement_match
{
    e_no_xor_match = 0,
    e_xor_match_GC,
    e_xor_match_BC,
    e_xor_match_iBC,
    e_xor_match_EQ
};

/*!
    Function calculates basic complexity statistics on XOR product of
    two blocks (b1 XOR b2)
    @ingroup bitfunc
    @internal
*/
inline
void bit_block_xor_change32(const bm::word_t* BMRESTRICT block,
                            const bm::word_t* BMRESTRICT xor_block,
                            unsigned size,
                            unsigned* BMRESTRICT gc,
                            unsigned* BMRESTRICT bc) BMNOEXCEPT
{
    BM_ASSERT(gc && bc);

    unsigned gap_count = 1;
    unsigned bit_count = 0;

    bm::word_t  w, w0, w_prev, w_l;
    w = w0 = *block ^ *xor_block;
    bit_count += word_bitcount(w);

    const int w_shift = int(sizeof(w) * 8 - 1);
    w ^= (w >> 1);
    gap_count += bm::word_bitcount(w);
    gap_count -= (w_prev = (w0 >> w_shift)); // negative value correction

    const bm::word_t* block_end = block + size; // bm::set_block_size;
    for (++block, ++xor_block; block < block_end; ++block, ++xor_block)
    {
        w = w0 = *block ^ *xor_block;
        bit_count += bm::word_bitcount(w);
        ++gap_count;
        if (!w)
        {
            gap_count -= !w_prev;
            w_prev = 0;
        }
        else
        {
            w ^= (w >> 1);
            gap_count += bm::word_bitcount(w);

            w_l = w0 & 1;
            gap_count -= (w0 >> w_shift);  // negative value correction
            gap_count -= !(w_prev ^ w_l);  // word border correction

            w_prev = (w0 >> w_shift);
        }
    } // for

    *gc = gap_count;
    *bc = bit_count;
}

/*!
    Function calculates number of times when bit value changed
    @internal
*/
inline
void bit_block_xor_change(const bm::word_t* BMRESTRICT block,
                              const bm::word_t* BMRESTRICT xor_block,
                              unsigned size,
                              unsigned* BMRESTRICT gc,
                              unsigned* BMRESTRICT bc) BMNOEXCEPT
{
#ifdef VECT_BLOCK_XOR_CHANGE
    VECT_BLOCK_XOR_CHANGE(block, xor_block, size, gc, bc);
#else
    bm::bit_block_xor_change32(block, xor_block, size, gc, bc);
#endif
}

/**
    Structure to compute XOR gap-count profile by sub-block waves
    @ingroup bitfunc
    @internal
*/
struct block_waves_xor_descr
{
    // stats for s-block
    unsigned short sb_gc[bm::block_waves]; ///< GAP counts
    unsigned short sb_bc[bm::block_waves]; ///< BIT counts

    // stats ref.block XOR mask
    unsigned short sb_xor_gc[bm::block_waves]; ///< XOR-mask GAP count
    unsigned short sb_xor_bc[bm::block_waves]; ///< XOR-mask GAP count
};

/**
    Capture the XOR filter results (xor block against ref.block)
    @internal
 */
struct block_xor_match_descr
{
    typedef bvector_size_type   size_type;

    bm::xor_complement_match  match_type; ///< match type
    unsigned                  block_gain; ///< XOR filter improvement (best)
    size_type                 ref_idx;    ///< reference vector index
    bm::id64_t                xor_d64;    ///< recorded digest

    // sub-block waves masks for metrics
    bm::id64_t   gc_d64;
    bm::id64_t   bc_d64;
    bm::id64_t   ibc_d64;

    // recorded gains for metrics
    unsigned     gc_gain;
    unsigned     bc_gain;
    unsigned     ibc_gain;

    block_xor_match_descr() : match_type(e_no_xor_match) {}
};

/**
    XOR match pair
    @internal
 */
struct match_pair
{
    bvector_size_type         ref_idx;    ///< reference vector index
    bm::id64_t                xor_d64;    ///< recorded digest

    match_pair(){}
    match_pair(bvector_size_type idx, bm::id64_t d64)
    : ref_idx(idx), xor_d64(d64)
    {}

};

/**
    XOR match chain
    @internal
 */
template<typename BLOCK_IDX>
struct block_match_chain
{
    BLOCK_IDX   nb;
    unsigned    chain_size;
    unsigned    ref_idx[64];
    bm::id64_t  xor_d64[64];
    bm::xor_complement_match match;
};

/**
    Greedy algorithm to find additional matches
    improving the inital best match block on its match type

    @param match_pairs_vect - [out] target vector of best match pairs
    @param match_vect - [in/out] vector of all found match descriptors

    @return number of new finds (if any)
 */
template<typename PVT, typename VT>
typename VT::size_type
greedy_refine_match_vector(PVT&                      match_pairs_vect,
                           VT&                       match_vect,
                           typename VT::size_type    best_ref_idx,
                           bm::id64_t                d64,
                           bm::xor_complement_match  match_type)
{
    BM_ASSERT(match_type && d64);
    match_pairs_vect.resize(0);

    bm::id64_t  d64_acc(d64);

    // pass 1 (exact match)
    //
    typename VT::size_type sz = match_vect.size();
    for (typename VT::size_type i = 0; (i < sz) && (d64_acc != ~0ull); ++i)
    {
        block_xor_match_descr& xmd = match_vect[i];
        if (xmd.ref_idx == best_ref_idx) // self hit
            continue;
        if (xmd.match_type == match_type) // best compatible match types
        {
            bm::id64_t d64_new = ~d64_acc & xmd.xor_d64;
            if (d64_new)
            {
                d64_acc |= d64_new;    // add mask to accum
                match_pairs_vect.push_back(bm::match_pair(xmd.ref_idx, d64_new));
            }
            xmd.match_type = e_no_xor_match; // mark it to exclude from pass 2
        }
    } // for i

    // pass 2 (extended match)
    //
    const unsigned min_gain_cut_off = 50;
    for (typename VT::size_type i = 0; (i < sz) && (d64_acc != ~0ull); ++i)
    {
        block_xor_match_descr& xmd = match_vect[i];
        if (xmd.ref_idx == best_ref_idx || !xmd.match_type) // self hit or none
            continue;
        BM_ASSERT(xmd.match_type != match_type);

        bm::id64_t d64_new = 0;
        switch (match_type)
        {
        case e_xor_match_GC:
            if (xmd.gc_gain > min_gain_cut_off)
                d64_new = ~d64_acc & xmd.gc_d64;
            break;
        case e_xor_match_BC:
            if (xmd.bc_gain > min_gain_cut_off)
                d64_new = ~d64_acc & xmd.bc_d64;
            break;
        case e_xor_match_iBC:
            if (xmd.ibc_gain > min_gain_cut_off)
                d64_new = ~d64_acc & xmd.ibc_d64;
            break;
        default:
            break;
        } // switch

        if (d64_new) // some improvement found
        {
            d64_acc |= d64_new;    // add mask to accum
            match_pairs_vect.push_back(bm::match_pair(xmd.ref_idx, d64_new));
            xmd.match_type = e_no_xor_match;
        }
    } // for

    return match_pairs_vect.size();
}

/**
    Check effective bit-rate for the XOR encode vector
    @return 1 - < 256 (8bit), 2 - < 65536 (16-bit) or 0 - 32-bit
    @internal
 */
template<typename PVT>
unsigned char check_pair_vect_vbr(const PVT& match_pairs_vect,
                             typename PVT::size_type  ref_idx)
{
    typename PVT::size_type max_idx = 0;
    if (ref_idx > max_idx)
        max_idx = ref_idx;
    for (typename PVT::size_type i = 0; i < match_pairs_vect.size(); ++i)
    {
        const match_pair& mp = match_pairs_vect[i];
        if (mp.ref_idx > max_idx)
            max_idx = mp.ref_idx;
    } // for i
    if (max_idx < 256)
        return 1;
    if (max_idx < 65536)
        return 2;
    return 0;
}


/**
    Compute reference (non-XOR) 64-dim complexity descriptor for the
    s-block.
    Phase 1 of the XOR filtering process is to establish the base metric

    @internal
*/
inline
void compute_s_block_descr(const bm::word_t* BMRESTRICT block,
                        block_waves_xor_descr& BMRESTRICT x_descr) BMNOEXCEPT
{
    // TODO: SIMD (for loop can go inside VECT to minimize LUT re-inits)
    for (unsigned i = 0; i < bm::block_waves; ++i)
    {
        unsigned off = (i * bm::set_block_digest_wave_size);
        const bm::word_t* sub_block = block + off;
        unsigned gc, bc;
        // TODO: optimize to compute GC and BC in a single pass
        #if defined(VECT_BLOCK_CHANGE)
            gc = VECT_BLOCK_CHANGE(sub_block, bm::set_block_digest_wave_size);
        #else
            gc = bm::bit_block_change32(sub_block, bm::set_block_digest_wave_size);
        #endif
        x_descr.sb_gc[i] = (unsigned short) gc;
        bc = bm::bit_count_min_unroll(
                    sub_block, sub_block + bm::set_block_digest_wave_size);
        x_descr.sb_bc[i] = (unsigned short) bc;
    } // for i
    // TODO: comute and return d64 - use it later
}


/**
    Compute reference complexity descriptor based on XOR vector.
    Returns the digest of sub-blocks where XOR filtering improved the metric
    (function needs reference to estimate the improvement).

    part of Phase 2 of the XOR filtering process

    @sa compute_sub_block_complexity_descr

    @internal
*/
inline
void compute_xor_complexity_descr(
                        const bm::word_t* BMRESTRICT block,
                        const bm::word_t* BMRESTRICT xor_block,
                        bm::block_waves_xor_descr& BMRESTRICT x_descr,
                        bm::block_xor_match_descr& BMRESTRICT xmd)
{
    // digest of all ZERO sub-blocks
    bm::id64_t d0 = ~bm::calc_block_digest0(block);

    // Pass 1: compute XOR descriptors
    //
    for (unsigned i = 0; i < bm::block_waves; ++i)
    {
        unsigned off = (i * bm::set_block_digest_wave_size);
        const bm::word_t* sub_block = block + off;
        const bm::word_t* xor_sub_block = xor_block + off;

        unsigned xor_gc, xor_bc;
        bm::bit_block_xor_change(sub_block, xor_sub_block,
                                 bm::set_block_digest_wave_size,
                                 &xor_gc, &xor_bc);
        x_descr.sb_xor_gc[i] = (unsigned short)xor_gc;
        x_descr.sb_xor_bc[i] = (unsigned short)xor_bc;
    } // for i

    // Pass 2: find the best match
    //
    unsigned block_gc_gain(0), block_bc_gain(0), block_ibc_gain(0);
    bm::id64_t gc_digest(0), bc_digest(0), ibc_digest(0);
    const unsigned wave_max_bits = bm::set_block_digest_wave_size * 32;

    for (unsigned i = 0; i < bm::block_waves; ++i)
    {
        bm::id64_t dmask = (1ull << i);
        if (d0 & dmask)
            continue;

        unsigned xor_gc = x_descr.sb_xor_gc[i];
        if (xor_gc <= 1)
        {
            gc_digest |= dmask;
            block_gc_gain += x_descr.sb_gc[i]; // all previous GAPs are gone
        }
        else if (xor_gc < x_descr.sb_gc[i]) // some improvement in GAPs
        {
            gc_digest |= dmask;
            block_gc_gain += (x_descr.sb_gc[i] - xor_gc);
        }
        unsigned xor_bc = x_descr.sb_xor_bc[i];
        if (xor_bc < x_descr.sb_bc[i]) // improvement in BITS
        {
            bc_digest |= dmask;
            block_bc_gain += (x_descr.sb_bc[i] - xor_bc);
        }
        unsigned xor_ibc = wave_max_bits - xor_bc;
        unsigned wave_ibc = wave_max_bits - x_descr.sb_bc[i];
        if (xor_ibc < wave_ibc) // improvement in 0 BITS
        {
            ibc_digest |= dmask;
            block_ibc_gain += (wave_ibc - xor_ibc);
        }

    } // for i

    // Save metric0))) results into XOR match descriptor
    //

    xmd.gc_d64 = gc_digest;
    xmd.bc_d64 = bc_digest;
    xmd.ibc_d64 = ibc_digest;

    xmd.gc_gain = block_gc_gain;
    xmd.bc_gain = block_bc_gain;
    xmd.ibc_gain = block_ibc_gain;


    // Find the winning metric and its digest mask
    //

    if (!(block_gc_gain | block_bc_gain | block_ibc_gain)) // match not found
    {
        // this is to check if xor filter has
        // canceled out a whole sub-block wave
        //
        bm::id64_t d0_x = ~bm::calc_block_digest0(xor_block);
        if (d0 == d0_x)
        {
            xmd.match_type = bm::e_xor_match_GC;
            xmd.block_gain = bm::block_waves;
            xmd.xor_d64 = d0;
            return;
        }

        xmd.match_type = bm::e_no_xor_match; xmd.block_gain = 0; xmd.xor_d64 = 0;
        return;
    }

    if (block_gc_gain > block_bc_gain)
    {
        if (block_gc_gain > block_ibc_gain)
        {
            xmd.match_type = bm::e_xor_match_GC;
            xmd.block_gain = block_gc_gain;
            xmd.xor_d64 = gc_digest;
            return;
        }
    }
    else // GC_gain <= BC_gain
    {
        if (block_bc_gain > block_ibc_gain)
        {
            xmd.match_type = bm::e_xor_match_BC; xmd.block_gain = block_bc_gain;
            xmd.xor_d64 = bc_digest;
            return;
        }
    }

    xmd.match_type = bm::e_xor_match_iBC;
    xmd.block_gain = block_ibc_gain;
    xmd.xor_d64 = ibc_digest;
    return;
}

/**
    Build partial XOR product of 2 bit-blocks using digest mask

    @param target_block - target := block ^ xor_block
    @param block - arg1
    @param xor_block - arg2
    @param digest - mask for each block wave to XOR (1) or just copy (0)

    @internal
*/
inline
void bit_block_xor(bm::word_t*  target_block,
                   const bm::word_t*  block, const bm::word_t* xor_block,
                   bm::id64_t digest) BMNOEXCEPT
{
    BM_ASSERT(target_block);
    BM_ASSERT(block);
    BM_ASSERT(xor_block);
    BM_ASSERT(digest);

#ifdef VECT_BIT_BLOCK_XOR
    VECT_BIT_BLOCK_XOR(target_block, block, xor_block, digest);
#else
    for (unsigned i = 0; i < bm::block_waves; ++i)
    {
        const bm::id64_t mask = (1ull << i);
        unsigned off = (i * bm::set_block_digest_wave_size);
        const bm::word_t* sub_block = block + off;
        bm::word_t* t_sub_block = target_block + off;

        const bm::word_t* sub_block_end = sub_block + bm::set_block_digest_wave_size;

        if (digest & mask) // XOR filtered sub-block
        {
            const bm::word_t* xor_sub_block = xor_block + off;
            // TODO: unroll
            for (; sub_block < sub_block_end; )
                *t_sub_block++ = *sub_block++ ^ *xor_sub_block++;
        }
        else // just copy source
        {
            for (; sub_block < sub_block_end;)
                *t_sub_block++ = *sub_block++;
        }
    } // for i
#endif
}


/**
    Build partial XOR product of 2 bit-blocks using digest mask

    @param target_block - target := target ^ xor_block
    @param xor_block - arg
    @param digest - mask for each block wave to XOR (1)

    @internal
*/
inline
void bit_block_xor(bm::word_t* target_block, const bm::word_t*  xor_block,
                   bm::id64_t digest) BMNOEXCEPT
{
    BM_ASSERT(target_block);
    BM_ASSERT(xor_block);
    BM_ASSERT(digest);

    // TODO: SIMD
    //

    while (digest)
    {
        bm::id64_t t = bm::bmi_blsi_u64(digest); // d & -d;

        unsigned wave = bm::word_bitcount64(t - 1);
        unsigned off = wave * bm::set_block_digest_wave_size;

        const bm::word_t* xor_sub_block = xor_block + off;
        bm::word_t* t_sub_block = target_block + off;
        const bm::word_t* t_sub_block_end = t_sub_block + bm::set_block_digest_wave_size;

        // TODO: unroll
        for (; t_sub_block < t_sub_block_end; )
            *t_sub_block++ ^= *xor_sub_block++;

        digest = bm::bmi_bslr_u64(digest); // d &= d - 1;
    } // while
}

/**
    Dry XOR operation on two GAP blocks to compute BC/GC metrics
 */
template<typename T>
void gap_operation_dry_xor(const T*  vect1, const T*  vect2,
                           unsigned& gc,    unsigned& bc) BMNOEXCEPT
{
    const T*  cur1 = vect1;
    const T*  cur2 = vect2;

    unsigned gc_cnt(1), bc_cnt(0);

    T bitval1 = (T)(*cur1++ & 1);
    T bitval2 = (T)(*cur2++ & 1);
    T bitval = bitval1 ^ bitval2;
    T bitval_prev = bitval;

    T c1(*cur1), c2(*cur2);
    T res(0), res_prev(0);

    while (1)
    {
        bitval = bitval1 ^ bitval2;
        if (bitval != bitval_prev)
        {
            ++gc_cnt;
            bitval_prev = bitval;
            res_prev = res;
        }

        if (c1 < c2)
        {
            res = c1;
            if (bitval)
            {
                bc_cnt += res - res_prev;
                res_prev = res;
            }
            ++cur1; c1 = *cur1; bitval1 ^= 1;
        }
        else // >=
        {
            res = c2;
            if (bitval)
            {
                bc_cnt += res - res_prev;
                res_prev = res;
            }
            if (c2 < c1)
            {
                bitval2 ^= 1;
            }
            else  // equal
            {
                if (c2 == (bm::gap_max_bits - 1))
                    break;
                ++cur1; c1 = *cur1;
                bitval1 ^= 1; bitval2 ^= 1;
            }
            ++cur2; c2 = *cur2;
        }
    } // while

    gc = gc_cnt;
    bc = bc_cnt;
}

/**
    List of reference bit-vectors with their true index associations

    Each referece vector would have two alternative indexes:
     - index(position) in the reference list
     - index(row) in the external bit-matrix (plane index)

    @internal
*/
template<typename BV>
class bv_ref_vector
{
public:
    typedef BV                                          bvector_type;
    typedef typename bvector_type::size_type            size_type;
    typedef bvector_type*                               bvector_type_ptr;
    typedef const bvector_type*                         bvector_type_const_ptr;
    typedef typename bvector_type::allocator_type       bv_allocator_type;


    typedef bm::block_match_chain<size_type>            block_match_chain_type;
    typedef
    bm::dynamic_heap_matrix<block_match_chain_type, bv_allocator_type>
                                                             matrix_chain_type;

public:

    /// reset the collection (resize(0))
    void reset()
    {
        rows_acc_ = 0;
        ref_bvects_.resize(0); ref_bvects_rows_.resize(0);
    }

    /**
        Add reference vector
        @param bv - bvector pointer
        @param ref_idx - reference (row) index
    */
    void add(const bvector_type* bv, size_type ref_idx)
    {
        BM_ASSERT(bv);
        ref_bvects_.push_back(bv);
        ref_bvects_rows_.push_back(ref_idx);
    }

    /// Get reference list size
    size_type size() const BMNOEXCEPT { return (size_type)ref_bvects_.size(); }

    /// Get reference vector by the index in this ref-vector
    const bvector_type* get_bv(size_type idx) const BMNOEXCEPT
                                        { return ref_bvects_[idx]; }

    /// Get reference row index by the index in this ref-vector
    size_type get_row_idx(size_type idx) const BMNOEXCEPT
                        { return (size_type)ref_bvects_rows_[idx]; }

    /// not-found value for find methods
    static
    size_type not_found() BMNOEXCEPT { return ~(size_type(0)); }

    /// Find vector index by the reference index
    /// @return ~0 if not found
    size_type find(std::size_t ref_idx) const BMNOEXCEPT
    {
        size_type sz = size();
        for (size_type i = 0; i < sz; ++i) // TODO: optimization
            if (ref_idx == ref_bvects_rows_[i])
                return i;
        return not_found();
    }

    /// Find vector index by the pointer
    /// @return ~0 if not found
    size_type find_bv(const bvector_type* bv) const BMNOEXCEPT
    {
        size_type sz = size();
        for (size_type i = 0; i < sz; ++i)
            if (bv == ref_bvects_[i])
                return i;
        return not_found();
    }

    /// Fill block allocation digest for all vectors in the reference collection
    ///   @param bv_blocks - [out] bvector of blocks statistics
    ///
    void fill_alloc_digest(bvector_type& bv_blocks) const
    {
        size_type sz = size();
        if (sz)
        {
            for (size_type i = 0; i < sz; ++i)
                ref_bvects_[i]->fill_alloc_digest(bv_blocks);
            BM_DECLARE_TEMP_BLOCK(tb);
            bv_blocks.optimize(tb);
        }
    }

    /// Reset and build vector of references from a basic bit-matrix
    ///  all NULL rows are skipped, not added to the ref.vector
    /// @sa add_vectors
    ///
    template<class BMATR>
    void build(const BMATR& bmatr)
    {
        reset();
        add_vectors(bmatr);
    }

    /// Append basic bit-matrix to the list of reference vectors
    /// @sa build
    /// @sa add_sparse_vector
    template<typename BMATR>
    void add_vectors(const BMATR& bmatr)
    {
        size_type rows = bmatr.rows();
        for (size_type r = 0; r < rows; ++r)
        {
            bvector_type_const_ptr bv = bmatr.get_row(r);
            if (bv)
                add(bv, rows_acc_ + r);
        } // for r
        rows_acc_ += unsigned(rows);
    }

    /// Add bit-transposed sparse vector as a bit-matrix
    /// @sa add_vectors
    ///
    template<class SV>
    void add_sparse_vector(const SV& sv)
    {
        add_vectors(sv.get_bmatrix());
    }

    /** Utility function to resize matrix based on number of vectors and blocks
    */
    void resize_xor_matrix(matrix_chain_type& matr,
                           size_type total_blocks) const
    {
        BM_ASSERT(total_blocks);
        matr.resize(ref_bvects_.size(), total_blocks);
    }

    /** Calculate blocks digest and resize XOR distance matrix
        based on total number of available blocks
     */
    void build_nb_digest_and_xor_matrix(matrix_chain_type& matr,
                                        bvector_type& bv_blocks) const
    {
        fill_alloc_digest(bv_blocks);
        size_type cnt = bv_blocks.count();
        resize_xor_matrix(matr, cnt);
    }

protected:
    typedef bm::heap_vector<bvector_type_const_ptr, bv_allocator_type, true> bvptr_vector_type;
    typedef bm::heap_vector<std::size_t, bv_allocator_type, true> bv_plane_vector_type;

protected:
    unsigned                 rows_acc_ = 0;     ///< total rows accumulator
    bvptr_vector_type        ref_bvects_;       ///< reference vector pointers
    bv_plane_vector_type     ref_bvects_rows_;  ///< reference vector row idxs
};

// --------------------------------------------------------------------------
//
// --------------------------------------------------------------------------

/**
    XOR similarity model

    @internal
*/

template<typename BV>
struct xor_sim_model
{
public:
    typedef BV                                          bvector_type;
    typedef typename bvector_type::size_type            size_type;
    typedef typename bvector_type::allocator_type       bv_allocator_type;


    typedef bm::block_match_chain<size_type>            block_match_chain_type;
    typedef
    bm::dynamic_heap_matrix<block_match_chain_type, bv_allocator_type>
                                                             matrix_chain_type;

public:
    matrix_chain_type  matr; ///< model matrix
    bvector_type       bv_blocks;  ///< blocks digest
};

// --------------------------------------------------------------------------
//
// --------------------------------------------------------------------------

/**
    XOR scanner to search for complement-similarities in
    collections of bit-vectors

    @internal
*/
template<typename BV>
class xor_scanner
{
public:
    typedef bm::bv_ref_vector<BV>                 bv_ref_vector_type;
    typedef BV                                    bvector_type;
    typedef typename bvector_type::allocator_type bv_allocator_type;
    typedef typename bvector_type::size_type      size_type;

    typedef bm::heap_vector<bm::block_xor_match_descr, bv_allocator_type, true>
            xor_matches_vector_type;
    typedef bm::heap_vector<bm::match_pair, bv_allocator_type, true>
            match_pairs_vector_type;
    typedef typename bv_ref_vector_type::matrix_chain_type
                                                    matrix_chain_type;
    typedef bm::heap_vector<bm::word_t*, bv_allocator_type, true>
                                                    bv_blocks_vector_type;

public:

    xor_scanner() {}
    ~xor_scanner()
    {
        free_blocks();
    }


    void set_ref_vector(const bv_ref_vector_type* ref_vect) BMNOEXCEPT
    { ref_vect_ = ref_vect; }

    const bv_ref_vector_type& get_ref_vector() const BMNOEXCEPT
    { return *ref_vect_; }

    /** Compute statistics for the anchor search block
        @param block - bit-block target
    */
    void compute_s_block_stats(const bm::word_t* block) BMNOEXCEPT;

    /** Scan for all candidate bit-blocks to find mask or match
        @return XOR referenece match type
    */
    bm::xor_complement_match
    search_best_xor_mask(const bm::word_t* s_block,
                              size_type ri,
                              size_type ridx_from,
                              size_type ridx_to,
                              unsigned i, unsigned j,
                              bm::word_t* tx_block);
    /**
        Run a search to add possible XOR match chain additions
     */
    size_type refine_match_chain();

    /** Scan all candidate gap-blocks to find best XOR match
    */
    bool search_best_xor_gap(const bm::word_t* block,
                             size_type         ridx_from,
                             size_type         ridx_to,
                             unsigned i, unsigned j);

    /**
        XOR all match blocks to target using their digest masks
     */
     void apply_xor_match_vector(bm::word_t* target_xor_block,
                           const bm::word_t* block,
                           const bm::word_t* ref_block,
                            bm::id64_t        d64,
                            const match_pairs_vector_type& pm_vect,
                            unsigned i, unsigned j) const BMNOEXCEPT;

    /**
        Calculate matrix of best XOR match metrics per block
        for the attached collection of bit-vectors
     */
    void compute_sim_model(const bv_ref_vector_type& ref_vect,
                           xor_sim_model<BV> &sim_model);

    /**
        Calculate matrix of best XOR match metrics per block
        for the attached collection of bit-vectors
     */
    void compute_sim_model(xor_sim_model<BV> &sim_model);

    /**
        Check if XOR transform simplified block enough for
        compressibility objective
     */
    bool validate_xor(const bm::word_t* xor_block) const BMNOEXCEPT;

    size_type found_ridx() const BMNOEXCEPT { return found_ridx_; }

    /// Return best match type of a found block
    ///
    bm::xor_complement_match get_best_match_type() const BMNOEXCEPT
        {   return x_block_mtype_; }

    const bm::word_t* get_found_block() const BMNOEXCEPT
        { return found_block_xor_; }
    unsigned get_x_best_metric() const BMNOEXCEPT { return x_best_metric_; }
    bm::id64_t get_xor_digest() const BMNOEXCEPT { return x_d64_; }

    unsigned get_s_bc() const BMNOEXCEPT { return s_bc_; }
    unsigned get_s_gc() const BMNOEXCEPT { return s_gc_; }
    unsigned get_s_block_best() const BMNOEXCEPT
                    { return s_block_best_metric_; }


    bm::block_waves_xor_descr& get_descr() BMNOEXCEPT { return x_descr_; }

    static
    bm::xor_complement_match best_metric(unsigned bc, unsigned gc,
                                        unsigned* best_metric);

    xor_matches_vector_type& get_match_vector() BMNOEXCEPT
        { return match_vect_; }

    match_pairs_vector_type& get_match_pairs() BMNOEXCEPT
        { return chain_match_vect_; }

    /// Return block from the reference vector [vect_idx, block_i, block_j]
    ///
    const bm::word_t* get_ref_block(size_type ri,
                                    unsigned i, unsigned j) const BMNOEXCEPT
    { return ref_vect_->get_bv(ri)->get_blocks_manager().get_block_ptr(i, j); }

protected:

    /// Deoptimize vertical slice of GAP blocks
    /// @param nb - block number
    ///
    void deoptimize_gap_blocks(size_type nb);

    /// Free the collection of temp blocks
    void free_blocks() BMNOEXCEPT;

    /// Sync TEMP vector size
    void sync_nb_vect();

private:
    xor_scanner(const xor_scanner&) = delete;
    xor_scanner& operator=(const xor_scanner&) = delete;

private:
    const bv_ref_vector_type*   ref_vect_ = 0; ///< ref.vect for XOR filter
    bv_blocks_vector_type       nb_blocks_vect_;   ///< pointers to temp blocks
    bv_allocator_type           alloc_;            ///< allocator to produce blocks

    bm::block_waves_xor_descr        x_descr_;  ///< XOR desriptor

    // S-block statistics
    //
    unsigned                      s_bc_;     ///< bitcount
    unsigned                      s_gc_;     ///< gap count
    unsigned                      x_best_metric_; ///< min(gc, bc, ibc)

    unsigned                      s_block_best_metric_; ///< s-block metric
    bm::xor_complement_match      x_block_mtype_;       ///< metric type

    // scan related metrics
    bm::id64_t                    x_d64_;        ///< search digest
    size_type                     found_ridx_;   ///< match vector (in references)
    const bm::word_t*             found_block_xor_;

    // match chain members:
    //
    xor_matches_vector_type       match_vect_; ///< vector of match descr
    match_pairs_vector_type       chain_match_vect_; ///< refined match pairs
};

// --------------------------------------------------------------------------
//
// --------------------------------------------------------------------------
//
//  Naming conventions and glossary:
//
//  s_block - source serialization block to be XOR filtered (anchor block)
//  best_ref_block - best reference block picked for XOR transform
//  xt_block - s_block XOR best_ref_block
//



template<typename BV>
void xor_scanner<BV>::compute_s_block_stats(const bm::word_t* block) BMNOEXCEPT
{
    BM_ASSERT(IS_VALID_ADDR(block));
    BM_ASSERT(!BM_IS_GAP(block));
    BM_ASSERT(ref_vect_->size() > 0);

    // TODO: one pass compute or use approx GC, BC based on pre-computed sub-waves
    bm::compute_s_block_descr(block, x_descr_);
    bm::bit_block_change_bc(block, &s_gc_, &s_bc_);

    x_block_mtype_ = best_metric(s_bc_, s_gc_, &s_block_best_metric_);
    x_best_metric_ = s_block_best_metric_;
}

// --------------------------------------------------------------------------

template<typename BV>
bm::xor_complement_match
xor_scanner<BV>::search_best_xor_mask(const bm::word_t* s_block,
                                       size_type s_ri,
                                       size_type ridx_from,
                                       size_type ridx_to,
                                       unsigned i, unsigned j,
                                       bm::word_t* tx_block)
{
    BM_ASSERT(ridx_from <= ridx_to);
    BM_ASSERT(IS_VALID_ADDR(s_block));
    BM_ASSERT(tx_block);

    if (ridx_to > ref_vect_->size())
        ridx_to = ref_vect_->size();

    bm::xor_complement_match rb_found = e_no_xor_match;
    bm::id64_t d64 = 0;
    found_block_xor_ = 0;

    unsigned best_block_gain = 0;
    int best_ri = -1;

    match_vect_.resize(0);

    unsigned s_gc(0);
    bool s_gap = BM_IS_GAP(s_block);
    if (s_gap)
    {
        const bm::gap_word_t* gap_s_block = BMGAP_PTR(s_block);
        s_gc = bm::gap_length(gap_s_block);
        s_block = nb_blocks_vect_.at(s_ri);
        BM_ASSERT(s_block);
    }

    // scan pass: over all reference vectors
    //
    for (size_type ri = ridx_from; ri < ridx_to; ++ri)
    {
        const bm::word_t* ref_block = get_ref_block(ri, i, j);
        if (BM_IS_GAP(ref_block))
        {
            if (nb_blocks_vect_.size() > ri)
                ref_block = nb_blocks_vect_[ri];
        }
        if (!IS_VALID_ADDR(ref_block))
            continue;

        BM_ASSERT(s_block != ref_block);

        bm::block_xor_match_descr xmd;
        bm::compute_xor_complexity_descr(s_block, ref_block, x_descr_, xmd);
        if (xmd.xor_d64) // candidate XOR block found
        {
            xmd.ref_idx = ri;
            BM_ASSERT(xmd.match_type);
            if (xmd.block_gain > best_block_gain)
            {
                best_block_gain = xmd.block_gain;
                best_ri = int(ri);
                d64 = xmd.xor_d64;
                if (xmd.block_gain >= bm::gap_max_bits)
                    break;
            }
            match_vect_.push_back(xmd); // place into vector of matches
        }
    } // for ri

    found_ridx_ = size_type(best_ri);
    x_d64_ = d64;

    if (best_ri != -1) // found some gain, validate it now
    {
        // assumed that if XOR compression c_level is at the highest
        const float bie_bits_per_int = 3.0f; // c_level_ < 6 ? 3.75f : 3.0f;
        const unsigned bie_limit =
                unsigned(float(bm::gap_max_bits) / bie_bits_per_int);

        unsigned xor_bc, xor_gc, xor_ibc;
        const bm::word_t* ref_block = get_ref_block(size_type(best_ri), i, j);
        bool r_gap = BM_IS_GAP(ref_block);
        if (r_gap)
            ref_block = nb_blocks_vect_[size_type(best_ri)];
        found_block_xor_ = ref_block;

        // TODO: one pass operation?
        bm::bit_block_xor(tx_block, s_block, ref_block, d64);
        bm::bit_block_change_bc(tx_block, &xor_gc, &xor_bc);

        if (!xor_bc) // check if completely identical block?
        {
            x_best_metric_ = xor_bc;
            found_ridx_ = size_type(best_ri);
            x_block_mtype_ = rb_found = e_xor_match_BC;

            unsigned block_pos;
            bool found = bm::block_find_first_diff(s_block, ref_block, &block_pos);
            if (!found)
                x_block_mtype_ = rb_found = e_xor_match_EQ;
        }
        else // find the best matching metric (GC, BC, iBC, ..)
        {
            if (xor_gc < x_best_metric_)
            {
                if (xor_gc < bie_limit)
                    rb_found = e_xor_match_GC;
                x_block_mtype_ = e_xor_match_GC;
                x_best_metric_ = xor_gc;
            }
            if (xor_bc < x_best_metric_)
            {
                if (xor_bc < bie_limit)
                    rb_found = e_xor_match_BC;
                x_block_mtype_ = e_xor_match_BC;
                x_best_metric_ = xor_bc;
            }
            else
            {
                xor_ibc = bm::gap_max_bits - xor_bc;
                if (xor_ibc < x_best_metric_)
                {
                    if (xor_ibc < bie_limit)
                        rb_found = e_xor_match_iBC;
                    x_block_mtype_ = e_xor_match_iBC;
                    x_best_metric_ = xor_ibc;
                }
            }

            // double check if XOR improves compression
            // with accounted serialization overhead
            //
            if (rb_found)
            {
                BM_ASSERT(get_s_block_best() > x_best_metric_);
                unsigned gain = get_s_block_best() - x_best_metric_;
                gain *= 3; // use bit estimate (speculative)
                // gain should be greater than overhead for storing
                // reference data: xor token, digest-64, block idx
                unsigned gain_min = unsigned(sizeof(char) + sizeof(unsigned));
                if (d64 != ~0ull) // if mask is all 1s - it is not used
                    gain_min += (unsigned)sizeof(bm::id64_t);
                gain_min *= 8; // in bits
                if (gain > gain_min)
                    return rb_found;

                if ((r_gap & s_gap) && (d64 != ~0ULL)) // both blocks are GAPs
                {
                    bm::bit_block_xor_2way(tx_block, s_block, ref_block);
                    bm::bit_block_change_bc(tx_block, &xor_gc, &xor_bc);
                    if (xor_gc < s_gc)
                    {
                        x_d64_ = ~0ULL;
                        return e_xor_match_GC;
                    }
                }
                return e_no_xor_match;
            }
        }
    }
    return rb_found;
}

// --------------------------------------------------------------------------

template<typename BV>
typename xor_scanner<BV>::size_type xor_scanner<BV>::refine_match_chain()
{
    size_type match_size = 0;
    if (x_d64_ == ~0ull || !x_d64_)
        return match_size;
    bm::xor_complement_match mtype = get_best_match_type();
    match_size = (size_type)
        bm::greedy_refine_match_vector(
            chain_match_vect_, match_vect_, found_ridx_, x_d64_, mtype);
    return match_size;
}

// --------------------------------------------------------------------------

template<typename BV>
void xor_scanner<BV>::compute_sim_model(const bv_ref_vector_type& ref_vect,
                                        xor_sim_model<BV> &sim_model)
{
    const bv_ref_vector_type* ref_vect_curr = this->ref_vect_; // save ref-vect

    ref_vect_ = &ref_vect;
    compute_sim_model(sim_model);

    ref_vect_ = ref_vect_curr; // restore state
}

template<typename BV>
void xor_scanner<BV>::compute_sim_model(bm::xor_sim_model<BV>& sim_model)
{
    BM_ASSERT(ref_vect_);

    sim_model.bv_blocks.clear(false);
    size_type rsize = ref_vect_->size();
    ref_vect_->build_nb_digest_and_xor_matrix(sim_model.matr,
                                              sim_model.bv_blocks);

    sync_nb_vect();

    BM_DECLARE_TEMP_BLOCK(xor_tmp_block);

    typename bvector_type::enumerator en(sim_model.bv_blocks);
    for (size_type col = 0; en.valid(); ++en, ++col)
    {
        size_type nb = *en;

        deoptimize_gap_blocks(nb);

        unsigned i0, j0;
        bm::get_block_coord(nb, i0, j0);

        for (size_type ri=0; true; ++ri)
        {
            bm::block_match_chain<size_type>* m_row = sim_model.matr.row(ri);
            BM_ASSERT(col < sim_model.matr.cols());
            bm::block_match_chain<size_type>& bmc = m_row[col];
            bmc.nb = nb;
            bmc.chain_size = 0;
            bmc.match = e_no_xor_match;

            if (ri == rsize-1)
                break;

            const bm::word_t* s_block = get_ref_block(ri, i0, j0);
            if (!IS_VALID_ADDR(s_block))
                continue;

            const bm::word_t* s_block_stat=s_block;
            if (BM_IS_GAP(s_block))
            {
                s_block_stat = nb_blocks_vect_.at(ri);
                BM_ASSERT(s_block_stat);
            }
            compute_s_block_stats(s_block_stat);
            bmc.match = search_best_xor_mask(s_block, ri, ri+1, rsize,
                                              i0, j0, xor_tmp_block);

            // take index in the ref-vector (not translated to a plain number)
            //
            size_type ridx = found_ridx();
            switch (bmc.match)
            {
            case e_no_xor_match:
                continue;
            case e_xor_match_EQ:
                bmc.chain_size++;
                bmc.ref_idx[0] = unsigned(ridx);
                break;
            case e_xor_match_GC:
                bmc.chain_size++;
                bmc.ref_idx[0] = unsigned(ridx);
                bmc.xor_d64[0] = ~0ULL;
                break;
            default:
                bmc.chain_size++;
                bmc.ref_idx[0] = ridx;
                bmc.xor_d64[0] = get_xor_digest();
            } // switch

        } // for ri
    } // for en
}

// --------------------------------------------------------------------------

template<typename BV>
bool xor_scanner<BV>::search_best_xor_gap(const bm::word_t* block,
                                          size_type ridx_from,
                                          size_type ridx_to,
                                          unsigned i, unsigned j)
{
    BM_ASSERT(ridx_from <= ridx_to);
    BM_ASSERT(BM_IS_GAP(block));

    if (ridx_to > ref_vect_->size())
        ridx_to = ref_vect_->size();

    x_d64_ = 0;

    const bm::gap_word_t* gap_block = BMGAP_PTR(block);
    unsigned gap_len = bm::gap_length(gap_block);
    if (gap_len <= 3)
        return false;
    unsigned bc = bm::gap_bit_count_unr(gap_block);

    bool kb_found = false;
    unsigned best_gap_metric = gap_len;
    if (bc < best_gap_metric)
        best_gap_metric = bc;

    for (size_type ri = ridx_from; ri < ridx_to; ++ri)
    {
        const bvector_type* bv = ref_vect_->get_bv(ri);
        BM_ASSERT(bv);
        const typename bvector_type::blocks_manager_type& bman = bv->get_blocks_manager();
        const bm::word_t* block_xor = bman.get_block_ptr(i, j);
        if (!IS_VALID_ADDR(block_xor) || !BM_IS_GAP(block_xor))
            continue;

        const bm::gap_word_t* gap_xor_block = BMGAP_PTR(block_xor);
        unsigned gap_xor_len = bm::gap_length(gap_block);
        if (gap_xor_len <= 3)
            continue;

        BM_ASSERT(block != block_xor);

        unsigned res_gc, res_bc;
        bm::gap_operation_dry_xor(gap_block, gap_xor_block, res_gc, res_bc);
/*
        unsigned res_len;
        bm::gap_operation_xor(gap_block, gap_xor_block, tmp_buf, res_len);
        unsigned glen = bm::gap_length(tmp_buf);
        if (res_len > glen) // size overflow
            continue;
        unsigned res_bc = bm::gap_bit_count_unr(tmp_buf); */
        if (!res_bc) // identical block
        {
            best_gap_metric = res_bc;
            kb_found = true;
            found_ridx_ = ri;
            found_block_xor_ = (const bm::word_t*)gap_xor_block;
            x_block_mtype_ = e_xor_match_BC;
        }

/*
        unsigned res_len;
        bool f = bm::gap_operation_dry_xor(gap_block, gap_xor_block, res_len, best_gap_len); */
        if ((res_gc < best_gap_metric))
        {
            unsigned gain = best_gap_metric - res_gc;
            if (gain > 2)
            {
                best_gap_metric = res_gc;
                kb_found = true;
                found_ridx_ = ri;
                found_block_xor_ = (const bm::word_t*)gap_xor_block;
                x_block_mtype_ = e_xor_match_GC;
            }
        }
        if (res_bc < best_gap_metric)
        {
            unsigned gain = best_gap_metric - res_bc;
            if (gain > 2)
            {
                best_gap_metric = res_bc;
                kb_found = true;
                found_ridx_ = ri;
                found_block_xor_ = (const bm::word_t*)gap_xor_block;
                x_block_mtype_ = e_xor_match_BC;
            }
        }
        unsigned res_ibc = bm::gap_max_bits - res_bc;
        if (res_ibc < best_gap_metric)
        {
            unsigned gain = best_gap_metric - res_ibc;
            if (gain > 2)
            {
                best_gap_metric = res_ibc;
                kb_found = true;
                found_ridx_ = ri;
                found_block_xor_ = (const bm::word_t*)gap_xor_block;
                x_block_mtype_ = e_xor_match_iBC;
            }
        }

        if (best_gap_metric <= 1)
            break;
    } // for ri

    return kb_found;
}

// --------------------------------------------------------------------------

template<typename BV>
void xor_scanner<BV>::apply_xor_match_vector(
                       bm::word_t* target_xor_block,
                       const bm::word_t* block,
                       const bm::word_t* ref_block,
                       bm::id64_t        d64,
                       const match_pairs_vector_type& pm_vect,
                       unsigned i, unsigned j) const BMNOEXCEPT
{
    bm::bit_block_xor(target_xor_block, block, ref_block, d64);
    auto sz = pm_vect.size();
    for (typename match_pairs_vector_type::size_type k = 0; k < sz; ++k)
    {
        const bm::match_pair& mp = pm_vect[k];
        const bm::word_t* block_ref = get_ref_block(mp.ref_idx, i, j);
        bm::bit_block_xor(target_xor_block, block_ref, mp.xor_d64);
    } // for k

}

// --------------------------------------------------------------------------

template<typename BV>
bool xor_scanner<BV>::validate_xor(const bm::word_t* xor_block) const BMNOEXCEPT
{
    const float bie_bits_per_int = 3.0f;
    const unsigned bie_limit =
            unsigned(float(bm::gap_max_bits) / bie_bits_per_int);

    unsigned bc, gc;
    bm::bit_block_change_bc(xor_block, &gc, &bc);
    unsigned xor_best_metric;
    bm::xor_complement_match mtype = best_metric(bc, gc, &xor_best_metric);
    if (mtype && (xor_best_metric < get_s_block_best()))
    {
        unsigned gain = get_s_block_best() - xor_best_metric;
        gain *= 3; // use bit estimate (speculative)
        // gain should be greater than overhead for storing
        // reference data: xor token, digest-64, block idx
        unsigned gain_min =
           unsigned (sizeof(char) + sizeof(bm::id64_t) + sizeof(unsigned));
        gain_min *= 8; // in bits
        if (!bc || ((gain > gain_min) && (xor_best_metric < bie_limit)))
            return true;
    }
    return false;
}

// --------------------------------------------------------------------------


template<typename BV>
bm::xor_complement_match
xor_scanner<BV>::best_metric(unsigned bc, unsigned gc, unsigned* best_metric)
{
    BM_ASSERT(best_metric);
    unsigned ibc = bm::gap_max_bits - bc;
    if (!ibc)
    {
        *best_metric = gc;
        return e_xor_match_GC;
    }
    if (gc < bc) // GC < GC
    {
        if (gc < ibc)
        {
            *best_metric = gc;
            return e_xor_match_GC;
        }
    }
    else // GC >= BC
    {
        if (bc < ibc)
        {
            *best_metric = bc;
            return e_xor_match_BC;
        }
    }
    *best_metric = ibc;
    return e_xor_match_iBC;
}

// --------------------------------------------------------------------------


template<typename BV>
void xor_scanner<BV>::free_blocks() BMNOEXCEPT
{
    size_t sz = nb_blocks_vect_.size();
    for (size_t i = 0; i < sz; ++i)
    {
        bm::word_t* blk = nb_blocks_vect_[i];
        if (blk)
            alloc_.free_bit_block(blk);
    }
    nb_blocks_vect_.resize(0);
}

// --------------------------------------------------------------------------

template<typename BV>
void xor_scanner<BV>::deoptimize_gap_blocks(size_type nb)
{
    size_type rsize = ref_vect_->size();
    BM_ASSERT(nb_blocks_vect_.size() == rsize);
    unsigned i0, j0;
    bm::get_block_coord(nb, i0, j0);

    for (size_type ri=0; ri < rsize; ++ri)
    {
        const bm::word_t* block = get_ref_block(ri, i0, j0);
        if (BM_IS_GAP(block))
        {
            bm::word_t* t_block = nb_blocks_vect_.at(ri);
            if (!t_block)
            {
                t_block = alloc_.alloc_bit_block();
                nb_blocks_vect_[ri] = t_block;
            }
            bm::gap_convert_to_bitset(t_block, BMGAP_PTR(block));
        }
    } // for ri

}

// --------------------------------------------------------------------------

template<typename BV>
void xor_scanner<BV>::sync_nb_vect()
{
    size_type rsize = ref_vect_->size();
    if (nb_blocks_vect_.size() == rsize)
        return;
    free_blocks();
    nb_blocks_vect_.resize(rsize);
    bm::word_t** vect_data = nb_blocks_vect_.data();
    for (size_type i = 0; i < rsize; ++i)
        vect_data[i] = 0;
}

// --------------------------------------------------------------------------

} // namespace bm

#endif
