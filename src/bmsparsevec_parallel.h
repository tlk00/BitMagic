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

namespace bm
{

/**
    Builder class to prepare a batch of tasks for parallel optimization of
    a sparse vector
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

    class task_batch : public bm::task_batch<allocator_type>
    {
    };

    void build_plan(task_batch& batch,
                    sparse_vector_type& sv,
                    typename bvector_type::optmode opt_mode,
                    typename sparse_vector_type::statistics* st)
    {
        auto& tv = batch.get_task_vector();
        auto rsize = sv.get_bmatrix().rows();
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
        bm::task_description* tdescr = (bm::task_description*) argp;

        bvector_type* bv = static_cast<bvector_type*>(tdescr->ctx0);
        sv_statistics_type* st = static_cast<sv_statistics_type*>(tdescr->ctx1);
        optmode_type opt_mode = static_cast<optmode_type>(tdescr->param0);

        typename bvector_type::statistics stbv;
        stbv.reset();
        BM_DECLARE_TEMP_BLOCK(tb);
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

} // namespace bm

#endif
