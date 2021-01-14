#ifndef BMTASK__H__INCLUDED__
#define BMTASK__H__INCLUDED__
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

/*! \file bmtask.h
    \brief Task definitions for parallel programming with BitMagic

    The design intent is to make tasks reasonably compatible with
    different threading and execution models, possibly with different
    languages and non-C++ runtimes.

*/

#include "bmbuffer.h"

namespace bm
{

/** Typedef for a call-back function pointer (pthread conformant signature)
*/
typedef void* (*task_func_type)(void*);

/**
    Structure describing a task as a callback function, argument and result.

    Structure contains pointers but being a descriptor does not own
    their scope or nature.
*/
struct task_description
{
    enum task_flags
    {
        no_flag = 0,     ///< no flag specified
        barrier_ok = 1u,  ///< barrier waits all prev.tasks done without error
        barrier_any = (1u << 1), ///< barrier waits all prev tasks done (success or not)
        barrier_ok_delayed = (1u << 2)
    };

    task_func_type  func;      ///< pthread-like function callback
    void*           argp;      ///< arg pointer
    void*           ret;       ///< ret pointer

    void*           ctx0;      ///< reserved
    void*           ctx1;      ///< reserved
    bm::id64_t      param0;    ///< reserved

    /// Union to add extra flexible payload to tasks
    union
    {
        int                i32;;
        unsigned           u32;
        unsigned long long u64;
        float              fp32;
        double             fp64;
        void*              void_ptr;
    } payload0, payload1;

    bm::id64_t      flags;     ///< task flags to designate barriers
    unsigned        err_code;  ///< error code
    unsigned        done;      ///< 0 - pending

    task_description() BMNOEXCEPT {}

    task_description(task_func_type  f, void* argptr = 0) BMNOEXCEPT
    {
        init(f, argptr, 0, 0, 0);
    }

    void init(task_func_type  f, void* argptr,
              void* c0, void* c1, bm::id64_t p0) BMNOEXCEPT
    {
        func = f; argp = argptr;
        ret = 0; ctx0 = c0; ctx1 = c1;
        param0 = p0; flags = 0; err_code = done = 0;
        payload0.u64 = payload1.u64 = 0;
    }
};

/**
    Interface definition for a group of tasks (batch)
 */
class task_batch_base
{
public:
    typedef unsigned size_type;
public:
    virtual ~task_batch_base() {}

    /// Return size of batch
    virtual size_type size() = 0;

    /// Get task by number
    virtual task_description* get_task(size_type i) = 0;
};

/**
    Basic collection of tasks for parallel execution
 */
template<typename BVAlloc>
class task_batch : public task_batch_base
{
public:
    typedef BVAlloc                     bv_allocator_type;
    typedef task_batch_base::size_type  size_type;
    typedef
    bm::heap_vector<bm::task_description, bv_allocator_type, true> task_vector_type;

public:
    virtual size_type size() { return (size_type) task_vect_.size(); }
    virtual bm::task_description* get_task(size_type i)
        { return &task_vect_.at(i); }

    /// Get access to internal task vector
    ///
    task_vector_type& get_task_vector()  BMNOEXCEPT { return task_vect_; }
    const task_vector_type& get_task_vector() const  BMNOEXCEPT
    { return task_vect_; }

protected:
    task_vector_type      task_vect_; ///< list of tasks

};


/**
    Run task batch sequentially

    Function is used for testing and debugging purposes or as a reference
    to implement custom parallel executors

    @param tasks - collection of tasks to run
 */
inline
void run_task_batch(task_batch_base & tasks)
{
    task_batch_base::size_type batch_size = tasks.size();
    for (task_batch_base::size_type i = 0; i < batch_size; ++i)
    {
        bm::task_description* tdescr = tasks.get_task(i);
        tdescr->argp = tdescr; // restore the self referenece
        tdescr->ret = tdescr->func(tdescr->argp);
        tdescr->done = 1;
    } // for
}


/**
    Simple scoped lock guard
    @internal
 */
template<typename Lock> class lock_guard
{
public:
    lock_guard(Lock& lk) : lk_(lk) { lk_.lock(); }
    ~lock_guard() { lk_.unlock(); }
private:
    lock_guard(const lock_guard<Lock>&) = delete;
    lock_guard<Lock>& operator=(const lock_guard<Lock>&) = delete;
private:
    Lock& lk_;
};


} // namespace bm

#endif

