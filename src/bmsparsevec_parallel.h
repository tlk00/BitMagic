#ifndef BMSPARSEVEC_PARALLEL__H__INCLUDED__
#define BMSPARSEVEC_PARALLEL__H__INCLUDED__
/*
Copyright(c) 2020 Anatoliy Kuznetsov(anatoliy_kuznetsov at yahoo.com)

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

/*! \file bmsparsevec_parallel.h
    \brief Parallel planner for operations with sparse vectors
*/

#include "bmsparsevec_serial.h"

namespace bm
{

/**
    Builder class to prepare a batch of tasks for parallel optimization of
    a sparse vector
    @ingroup bmtasks
 */
template<typename SVect, typename Lock>
class optimize_plan_builder
{
public:
    typedef SVect                                     sparse_vector_type;
    typedef Lock                                      lock_type;

    typedef typename sparse_vector_type::bvector_type bvector_type;
    typedef typename bvector_type::allocator_type     allocator_type;
    typedef typename bvector_type::optmode            optmode_type;
    typedef typename sparse_vector_type::statistics   sv_statistics_type;

    struct task_batch : public bm::task_batch<allocator_type>
    {
        typedef bm::task_batch<allocator_type>          parent_type;
        typedef typename parent_type::task_vector_type  task_vector_type;
    };

    void build_plan(task_batch& batch,
                    sparse_vector_type& sv,
                    typename bvector_type::optmode opt_mode,
                    typename sparse_vector_type::statistics* st)
    {
        typename task_batch::task_vector_type& tv = batch.get_task_vector();
        auto rsize = sv.get_bmatrix().rows();
        tv.reserve(rsize);
        for (unsigned k = 0; k < rsize; ++k)
        {
            bvector_type* bv = sv.get_bmatrix().get_row(k);
            if (bv)
            {
                bm::task_description& tdescr = tv.add();
                tdescr.init(task_run, (void*)&tdescr,
                            (void*)bv, (void*)st, opt_mode);
            }
        } // for
    }


protected:
    /// Task execution Entry Point
    /// @internal
    static void* task_run(void* argp)
    {
        if (!argp)
            return 0;
        bm::task_description* tdescr = static_cast<bm::task_description*>(argp);

        bvector_type* bv = static_cast<bvector_type*>(tdescr->ctx0);
        sv_statistics_type* st = static_cast<sv_statistics_type*>(tdescr->ctx1);
        optmode_type opt_mode = static_cast<optmode_type>(tdescr->param0);

        typename bvector_type::statistics stbv;
        stbv.reset();
        BM_DECLARE_TEMP_BLOCK(tb)
        bv->optimize(tb, opt_mode, &stbv);

        if (st)
        {
            static lock_type lk;
            bm::lock_guard<lock_type> lg(lk);
            st->add(stbv);
        }
        return 0;
    }
};

/**
    Parallel plan builder for the XOR filter scanner
    @ingroup bmtasks
 */
template<typename BV>
class compute_sim_matrix_plan_builder
{
public:
    typedef BV                                       bvector_type;
    typedef typename BV::size_type                   size_type;
    typedef typename bvector_type::allocator_type    allocator_type;
    typedef bm::bv_ref_vector<BV>                    bv_ref_vector_type;


    struct task_batch : public bm::task_batch<allocator_type>
    {
        typedef bm::task_batch<allocator_type>          parent_type;
        typedef typename parent_type::task_vector_type  task_vector_type;
    };

    void build_plan(task_batch& batch,
                    bm::xor_sim_model<BV>& sim_model,
                    const bv_ref_vector_type& ref_vect,
                    const bm::xor_sim_params& xs_params)
    {
        sim_model.bv_blocks.clear(true);
        ref_vect.build_nb_digest_and_xor_matrix(sim_model.matr,
                                                 sim_model.bv_blocks);

        typename bvector_type::size_type nb_count = sim_model.bv_blocks.count();
        typename task_batch::task_vector_type& tv = batch.get_task_vector();
        tv.reserve(nb_count);

        typename bvector_type::enumerator en(sim_model.bv_blocks);
        for (size_type col = 0; en.valid(); ++en, ++col)
        {
            size_type nb = *en;
            bm::task_description& tdescr = tv.add();
            tdescr.init(task_run, (void*)&tdescr,
                        (void*)&sim_model, (void*)&ref_vect, nb);
            tdescr.payload0.u64 = col; // rank
            tdescr.payload1.void_ptr = (void*)&xs_params;
        } // for en
    }

protected:

    /// Task execution Entry Point
    /// @internal
    static void* task_run(void* argp)
    {
        thread_local bm::xor_scanner<BV> xor_scan;

        if (!argp)
            return 0;
        bm::task_description* tdescr =
            static_cast<bm::task_description*>(argp);
        bm::xor_sim_model<BV>* sim_model =
            static_cast<bm::xor_sim_model<BV>*>(tdescr->ctx0);
        const bv_ref_vector_type* bv_ref =
            static_cast<const bv_ref_vector_type*>(tdescr->ctx1);

        xor_scan.set_ref_vector(bv_ref);
        size_type nb   = (size_type) tdescr->param0;
        size_type rank = (size_type) tdescr->payload0.u64;
        const bm::xor_sim_params* params =
            static_cast<const bm::xor_sim_params*>(tdescr->payload1.void_ptr);

        xor_scan.sync_nb_vect();
        xor_scan.compute_sim_model(*sim_model, nb, rank, *params);
        return 0;
    }
};


/**
    Parallel plan builder for succinct sparse vector serialization

    @sa sparse_vector_serializer
    @ingroup bmtasks
 */
template<typename SV>
class sv_serialization_plan_builder
{
public:
    typedef SV                                       sparse_vector_type;
    typedef typename SV::bvector_type                bvector_type;
    typedef typename SV::size_type                   size_type;
    typedef typename bvector_type::allocator_type    allocator_type;
    typedef bm::bv_ref_vector<bvector_type>          bv_ref_vector_type;
    typedef bm::xor_sim_model<bvector_type>          xor_sim_model_type;


    struct serialization_params
    {
        serialization_params()
        : sb_bookmarks_(false),
          sb_range_(0),
          compression_level_(bm::set_compression_default),
          bv_ref_ptr_(0), sim_model_ptr_(0)
        {}

        bool            sb_bookmarks_; ///< Bookmarks flag
        unsigned        sb_range_;     ///< Desired bookmarks interval
        unsigned        compression_level_;

        const bv_ref_vector_type*   bv_ref_ptr_;
        const xor_sim_model_type*   sim_model_ptr_;
    };

    struct task_batch : public bm::task_batch<allocator_type>
    {
        typedef bm::task_batch<allocator_type>          parent_type;
        typedef typename parent_type::task_vector_type  task_vector_type;

        serialization_params s_params;
    };

public:
    sv_serialization_plan_builder()
    {}

    void set_bookmarks(bool enable, unsigned bm_interval = 256) BMNOEXCEPT
        { s_params_.sb_bookmarks_ = enable; s_params_.sb_range_ = bm_interval; }

    void set_xor_ref(const bv_ref_vector_type* bv_ref_ptr) BMNOEXCEPT
        {  s_params_.bv_ref_ptr_ = bv_ref_ptr; }

    void set_sim_model(const xor_sim_model_type* sim_model) BMNOEXCEPT
        { s_params_.sim_model_ptr_ = sim_model; }


    void build_plan(task_batch& batch,
                    sparse_vector_serial_layout<SV>&  sv_layout,
                    const sparse_vector_type& sv)
    {
        typename task_batch::task_vector_type& tv = batch.get_task_vector();
        unsigned planes = sv.stored_planes();
        tv.reserve(planes + 1); // +1 for finalization task

        batch.s_params = s_params_;

        for (unsigned i = 0; i < planes; ++i)
        {
            typename SV::bvector_type_const_ptr bv = sv.get_plane(i);
            if (!bv)  // empty plane
            {
                sv_layout.set_plane(i, 0, 0);
                continue;
            }

            bm::task_description& tdescr = tv.add();
            tdescr.init(task_run, (void*)&tdescr,
                        (void*)bv, (void*)&batch.s_params, 0);
            tdescr.ret = (void*)&sv_layout;

            if (s_params_.bv_ref_ptr_)
            {
                BM_ASSERT(batch.s_params.sim_model_ptr_);
                tdescr.payload0.u32 =
                    (unsigned)s_params_.bv_ref_ptr_->find_bv(bv);
                BM_ASSERT(tdescr.payload0.u32
                          != s_params_.bv_ref_ptr_->not_found());
            }
            else
            {
                // ref vector not set: see set_xor_ref()
                BM_ASSERT(!batch.s_params.sim_model_ptr_);
            }

        } // for i

        // Add barrier task at the end to finalize the compression

        bm::task_description& tdescr = tv.add();
        tdescr.init(task_run_final, (void*)&tdescr,
                    (void*)0, (void*)&batch.s_params, 0);
        tdescr.flags = bm::task_description::barrier_ok;
        tdescr.ret = (void*)&sv_layout;
    }
protected:
    /// Task execution Entry Point
    /// @internal
    static void* task_run(void* argp)
    {
        if (!argp)
            return 0;
        //TODO: full implementation
        return 0;
    }

    static void* task_run_final(void* argp)
    {
        if (!argp)
            return 0;
        //TODO: full implementation
        return 0;
    }
protected:
    serialization_params s_params_;
};

} // namespace bm

#endif
