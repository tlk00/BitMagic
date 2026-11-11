/** @file test_worker_error.cpp
    @brief Induce captured and pool-level failures on actual background workers.
*/
#include <atomic>
#include <exception>
#include <iostream>
#include <stdexcept>
#include <string>
#include "bm.h"
#include "bmthreadpool.h"
#include "inv_list_progress.h"
#include "inv_list_worker_error.h"

/** Run failure cases and require every task to finish before releasing storage.
    @return Zero on success, one for a failed propagation assertion.
*/
int main()
{
    try
    {
        using Pool = bm::thread_pool<bm::task_descr*, bm::spin_lock<bm::pad0_struct> >;
        Pool pool;
        bm::task_batch<bm::bvector<>::allocator_type> batch;
        std::exception_ptr errors[5];
        batch.get_task_vector().reserve(5);
        for (unsigned i = 0; i < 5; ++i)
        {
            bm::task_function_t task([&, i](void*) -> int
            {
                if (i == 3)
                {
                    return 19;
                }
                if (i == 4)
                {
                    throw std::runtime_error("uncaught task failure");
                }
                try
                {
                    if (i == 0)
                    {
                        throw std::runtime_error("induced deserialize error");
                    }
                    if (i == 1)
                    {
                        throw std::bad_alloc();
                    }
                    throw 7;
                }
                catch (...)
                {
                    errors[i] = std::current_exception();
                    return 1;
                }
            });
            batch.add(task, 0);
        } // for i
        pool.start(4);
        for (unsigned i = 0; i < 5; ++i)
        {
            pool.get_job_queue().push(batch.get_task(i));
        } // for i
        pool.set_stop_mode(Pool::stop_when_done);
        pool.join();
        const std::string expected[] = {"induced deserialize error", std::bad_alloc().what(), "non-standard worker exception",
                                  "task error code=19", "task error code=-1"};
        for (unsigned i = 0; i < 5; ++i)
        {
            inv_list_wait_for_task(pool, batch, i, false, []() {});
            bool caught = false;
            try
            {
                inv_list_check_worker(errors[i], batch.get_task(i)->err_code, "fixture worker");
            }
            catch (const std::runtime_error& e)
            {
                const std::string text(e.what());
                caught = text.find("fixture worker") != std::string::npos &&
                         text.find(expected[i]) != std::string::npos;
            }
            if (!caught)
            {
                throw std::runtime_error("Worker failure was lost");
            }
        } // for i
        inv_list_check_worker({}, 0, "successful task");
        std::cout << "Captured standard/allocation/non-standard exceptions and pool error codes PASS\n";
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
    return 0;
}
