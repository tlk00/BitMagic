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

/*!
    Function calculates number of times when bit value changed
    @internal
*/
inline
unsigned bit_block_xor_change32(const bm::word_t* BMRESTRICT block,
                                const bm::word_t* BMRESTRICT xor_block,
                                unsigned size) BMNOEXCEPT
{
    unsigned gap_count = 1;

    bm::word_t  w, w0, w_prev, w_l;
    w = w0 = *block ^ *xor_block;

    const int w_shift = int(sizeof(w) * 8 - 1);
    w ^= (w >> 1);
    BM_INCWORD_BITCOUNT(gap_count, w);
    gap_count -= (w_prev = (w0 >> w_shift)); // negative value correction

    const bm::word_t* block_end = block + size; // bm::set_block_size;
    for (++block, ++xor_block; block < block_end; ++block, ++xor_block)
    {
        w = w0 = *block ^ *xor_block;
        ++gap_count;
        if (!w)
        {
            gap_count -= !w_prev;
            w_prev = 0;
        }
        else
        {
            w ^= (w >> 1);
            BM_INCWORD_BITCOUNT(gap_count, w);

            w_l = w0 & 1;
            gap_count -= (w0 >> w_shift);  // negative value correction
            gap_count -= !(w_prev ^ w_l);  // word border correction

            w_prev = (w0 >> w_shift);
        }
    } // for
    return gap_count;
}

/*!
    Function calculates number of times when bit value changed
    @internal
*/
inline
unsigned bit_block_xor_change(const bm::word_t* BMRESTRICT block,
                              const bm::word_t* BMRESTRICT xor_block,
                              unsigned size) BMNOEXCEPT
{
#ifdef VECT_BLOCK_XOR_CHANGE
    return VECT_BLOCK_XOR_CHANGE(block, xor_block, size);
#else
    return bit_block_xor_change32(block, xor_block, size);
#endif
}

/**
    Structure to compute XOR gap-count profile by sub-block waves

    @ingroup bitfunc
    @internal
*/
struct block_waves_xor_descr
{
    unsigned short sb_change[bm::block_waves];
    unsigned short sb_xor_change[bm::block_waves];
};

/**
    Compute reference (non-XOR) 64-dim complexity descriptor for the
    target block.

    Phase 1 of the XOR filtering process

    @internal
*/
inline
void compute_complexity_descr(
                        const bm::word_t* BMRESTRICT block,
                        block_waves_xor_descr& BMRESTRICT x_descr) BMNOEXCEPT
{
    for (unsigned i = 0; i < bm::block_waves; ++i)
    {
        unsigned off = (i * bm::set_block_digest_wave_size);
        const bm::word_t* sub_block = block + off;
        #if defined(VECT_BLOCK_CHANGE)
            unsigned change VECT_BLOCK_CHANGE(sub_block,
                                              bm::set_block_digest_wave_size);
        #else
            unsigned change = bm::bit_block_change32(sub_block,
                                                bm::set_block_digest_wave_size);
        #endif
        x_descr.sb_change[i] = (unsigned short) change;
    } // for i
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
bm::id64_t compute_xor_complexity_descr(
                        const bm::word_t* BMRESTRICT block,
                        const bm::word_t* BMRESTRICT xor_block,
                        block_waves_xor_descr& BMRESTRICT x_descr,
                        unsigned& BMRESTRICT block_gain) BMNOEXCEPT
{
    block_gain = 0; // approximate block gain (sum of sub-waves)
    bm::id64_t digest = 0;
    for (unsigned i = 0; i < bm::block_waves; ++i)
    {
        unsigned off = (i * bm::set_block_digest_wave_size);
        const bm::word_t* sub_block = block + off;
        const bm::word_t* xor_sub_block = xor_block + off;

        unsigned xor_change =
                bm::bit_block_xor_change(sub_block, xor_sub_block,
                                         bm::set_block_digest_wave_size);

        x_descr.sb_xor_change[i] = (unsigned short)xor_change;
        if (xor_change <= 1)
        {
            digest |= (1ull << i);
            block_gain += x_descr.sb_change[i];
        }
        else
        {
            if (xor_change < x_descr.sb_change[i]) // detected improvement
            {
                digest |= (1ull << i);
                block_gain += (x_descr.sb_change[i] - xor_change);
            }
        }
    } // for i
    return digest;
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
                   const bm::word_t*  block, const bm::word_t*  xor_block,
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
    List of reference bit-vectors with their true index associations

    Each referece vector would have two alternative indexes:
     - index(position) in the reference list
     - index(row) in the external bit-matrix (plain index)

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
public:

    /// reset the collection (resize(0))
    void reset()
    {
        ref_bvects_.resize(0);
        ref_bvects_rows_.resize(0);
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
        for (size_type i = 0; i < sz; ++i)
            if (ref_idx == ref_bvects_rows_[i])
                return i;
        return not_found();
    }

    /// build vector of references from a basic bit-matrix
    ///  all NULL rows are skipped, not added to the ref.vector
    template<class BMATR>
    void build(const BMATR& bmatr)
    {
        reset();
        size_type rows = bmatr.rows();
        for (size_type r = 0; r < rows; ++r)
        {
            bvector_type_const_ptr bv = bmatr.get_row(r);
            if (bv)
                add(bv, r);
        } // for r
    }

protected:
    typedef bm::heap_vector<bvector_type_const_ptr, bv_allocator_type, true> bvptr_vector_type;
    typedef bm::heap_vector<std::size_t, bv_allocator_type, true> bv_plain_vector_type;

protected:
    bvptr_vector_type        ref_bvects_;       ///< reference vector pointers
    bv_plain_vector_type     ref_bvects_rows_;  ///< reference vector row idxs
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
    typedef bm::bv_ref_vector<BV>                bv_ref_vector_type;
    typedef BV                                   bvector_type;
    typedef typename bvector_type::size_type     size_type;

public:
    void set_ref_vector(const bv_ref_vector_type* ref_vect) BMNOEXCEPT
    { ref_vect_ = ref_vect; }

    const bv_ref_vector_type& get_ref_vector() const BMNOEXCEPT
    { return *ref_vect_; }

    /** Compute statistics for the anchor search vector
        @param block - bit-block target
    */
    void compute_x_block_stats(const bm::word_t* block) BMNOEXCEPT;

    /** Scan for all candidate bit-blocks to find mask or match
        @return true if XOR complement or matching vector found
    */
    bool search_best_xor_mask(const bm::word_t* block,
                              size_type ridx_from,
                              size_type ridx_to,
                              unsigned i, unsigned j,
                              bm::word_t* tb);

    /** Scan all candidate gap-blocks to find best XOR match
    */
    bool search_best_xor_gap(const bm::word_t* block,
                             size_type         ridx_from,
                             size_type         ridx_to,
                             unsigned i, unsigned j);

    /**
        Validate serialization target
    */
    bool validate_found(bm::word_t* xor_block,
                        const bm::word_t* block) const BMNOEXCEPT;

    size_type found_ridx() const BMNOEXCEPT { return found_ridx_; }
    const bm::word_t* get_found_block() const BMNOEXCEPT
    { return found_block_xor_; }
    unsigned get_x_best_metric() const BMNOEXCEPT { return x_best_metric_; }
    bm::id64_t get_xor_digest() const BMNOEXCEPT { return x_d64_; }

    /// true if completely identical vector found
    bool is_eq_found() const BMNOEXCEPT { return !x_best_metric_; }


    unsigned get_x_bc() const BMNOEXCEPT { return x_bc_; }
    unsigned get_x_gc() const BMNOEXCEPT { return x_gc_; }
    unsigned get_x_block_best() const BMNOEXCEPT
                    { return x_block_best_metric_; }


    bm::block_waves_xor_descr& get_descr() BMNOEXCEPT { return x_descr_; }

private:
    const bv_ref_vector_type*        ref_vect_ = 0; ///< ref.vect for XOR filter

    // target x-block related statistics
    //
    bm::block_waves_xor_descr     x_descr_;  ///< XOR desriptor
    unsigned                      x_bc_;     ///< bitcount
    unsigned                      x_gc_;     ///< gap count
    unsigned                      x_best_metric_; /// dynamic min(gc, bc)
    unsigned                      x_block_best_metric_; ///< min(gc, bc) initial

    // scan related metrics
    bm::id64_t                    x_d64_;        ///< search digest
    size_type                     found_ridx_;   ///< match vector (in references)
    const bm::word_t*             found_block_xor_;
};

// --------------------------------------------------------------------------
//
// --------------------------------------------------------------------------

template<typename BV>
void xor_scanner<BV>::compute_x_block_stats(const bm::word_t* block) BMNOEXCEPT
{
    BM_ASSERT(IS_VALID_ADDR(block));
    BM_ASSERT(!BM_IS_GAP(block));
    BM_ASSERT(ref_vect_->size() > 0);

    bm::compute_complexity_descr(block, x_descr_);
    bm::bit_block_change_bc(block, &x_gc_, &x_bc_);
    x_block_best_metric_ = x_best_metric_ = x_gc_ < x_bc_ ? x_gc_ : x_bc_;
}

// --------------------------------------------------------------------------

template<typename BV>
bool xor_scanner<BV>::search_best_xor_mask(const bm::word_t* block,
                                           size_type ridx_from,
                                           size_type ridx_to,
                                           unsigned i, unsigned j,
                                           bm::word_t* tb)
{
    BM_ASSERT(ridx_from <= ridx_to);
    BM_ASSERT(IS_VALID_ADDR(block));
    BM_ASSERT(!BM_IS_GAP(block));
    BM_ASSERT(tb);

    if (ridx_to > ref_vect_->size())
        ridx_to = ref_vect_->size();

    bool kb_found = false;
    bm::id64_t d64 = 0;
    found_block_xor_ = 0;

    unsigned best_block_gain = 0;
    int best_ri = -1;

    for (size_type ri = ridx_from; ri < ridx_to; ++ri)
    {
        const bvector_type* bv = ref_vect_->get_bv(ri);
        BM_ASSERT(bv);
        const typename bvector_type::blocks_manager_type& bman =
                                                bv->get_blocks_manager();
        const bm::word_t* block_xor = bman.get_block_ptr(i, j);
        if (!IS_VALID_ADDR(block_xor) || BM_IS_GAP(block_xor))
            continue;

        BM_ASSERT(block != block_xor);

        unsigned block_gain = 0;

        bm::id64_t xor_d64 =
            bm::compute_xor_complexity_descr(block, block_xor, x_descr_, block_gain);
        if (xor_d64) // candidate XOR block
        {
            if (block_gain > best_block_gain)
            {
                best_block_gain = block_gain;
                best_ri = int(ri);
                d64 = xor_d64;
            }
        }
    } // for ri

    if (best_ri != -1) // found some gain
    {
        unsigned xor_bc, xor_gc;
        const bvector_type* bv = ref_vect_->get_bv(size_type(best_ri));
        const typename bvector_type::blocks_manager_type& bman = bv->get_blocks_manager();
        const bm::word_t* block_xor = bman.get_block_ptr(i, j);

        bm::bit_block_xor(tb, block, block_xor, d64);
        bm::bit_block_change_bc(tb, &xor_gc, &xor_bc);

        if (xor_gc < x_best_metric_ && xor_gc < bm::bie_cut_off)
        {
            x_best_metric_ = xor_gc;
            kb_found = true;
            found_ridx_ = size_type(best_ri);
            found_block_xor_ = block_xor;
        }
        if (xor_bc < x_best_metric_ && xor_bc < bm::bie_cut_off)
        {
            x_best_metric_ = xor_bc;
            kb_found = true;
            found_ridx_ = size_type(best_ri);
            found_block_xor_ = block_xor;
            if (!xor_bc) // completely identical block?
            {
                unsigned pos;
                bool f = bm::bit_find_first_diff(block, block_xor, &pos);
                x_best_metric_ += f;
            }
        }
    }

    x_d64_ = d64;
    return kb_found;
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

    const bm::gap_word_t* gap_block = BMGAP_PTR(block);
    unsigned gap_len = bm::gap_length(gap_block);
    if (gap_len <= 3)
        return false;
    unsigned best_gap_len = gap_len;
    bool kb_found = false;

    for (size_type ri = ridx_from; ri < ridx_to; ++ri)
    {
        const bvector_type* bv = ref_vect_->get_bv(ri);
        BM_ASSERT(bv);
        const typename bvector_type::blocks_manager_type& bman = bv->get_blocks_manager();
        const bm::word_t* block_xor = bman.get_block_ptr(i, j);
        if (!IS_VALID_ADDR(block_xor))
            continue;
        if (!BM_IS_GAP(block_xor))
            continue;

        const bm::gap_word_t* gap_xor_block = BMGAP_PTR(block_xor);
        unsigned gap_xor_len = bm::gap_length(gap_block);
        if (gap_xor_len <= 3)
            continue;

        BM_ASSERT(block != block_xor);

        unsigned res_len;
        bool f = bm::gap_operation_dry_xor(gap_block, gap_xor_block, res_len, best_gap_len);
        if (f && (res_len < best_gap_len))
        {
            best_gap_len = res_len;
            kb_found = true;
            found_ridx_ = ri;
            found_block_xor_ = (const bm::word_t*)gap_xor_block;
        }

    } // for ri

    return kb_found;
}

// --------------------------------------------------------------------------

template<typename BV>
bool xor_scanner<BV>::validate_found(bm::word_t* xor_block,
                                     const bm::word_t* block) const BMNOEXCEPT
{
    bm::id64_t d64 = get_xor_digest();
    BM_ASSERT(d64);
    const bm::word_t* key_block = get_found_block();
    bm::bit_block_xor(xor_block, block, key_block, d64);

    unsigned bc, gc;
    bm::bit_block_change_bc(xor_block, &gc, &bc);
    unsigned xor_best_metric = gc < bc ? gc : bc;
    if (xor_best_metric < get_x_block_best())
    {
        unsigned gain = get_x_block_best() - xor_best_metric;
        gain *= 4; // take 4 as bit estimate for BIC
        // gain should be greater than overhead for storing
        // reference data: xor token, digest-64, block idx
        unsigned gain_min = unsigned (sizeof(char) + sizeof(bm::id64_t) + sizeof(unsigned));
        gain_min *= 8; // in bits
        if (gain > gain_min)
        {
            BM_ASSERT(xor_best_metric < bm::bie_cut_off);
            return true;
        }
    }
    return false;
}

// --------------------------------------------------------------------------

} // namespace bm

#endif
