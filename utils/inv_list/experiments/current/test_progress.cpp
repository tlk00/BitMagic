// Focused application progress checks; no dataset or serialized stream required.
#include "../../inv_list_progress.h"
#include <cassert>
#include <iostream>

/** Run deterministic checks; any argument also enables a real 31-second wait. */
int main(int argc, char**)
{
    using clock_type = inv_list_progress_timer::clock_type;
    const auto start = clock_type::time_point{};
    inv_list_progress_timer timer(true, start);
    assert(!timer.due(false, start + std::chrono::seconds(29)));
    assert(timer.due(false, start + std::chrono::seconds(30)));
    timer.reported(start + std::chrono::seconds(30));
    assert(!timer.due(false, start + std::chrono::seconds(59)));
    assert(timer.due(false, start + std::chrono::seconds(60)));
    assert(timer.due(true, start + std::chrono::seconds(31)));
    timer.reported(start + std::chrono::seconds(31));
    assert(!timer.due(false, start + std::chrono::seconds(60)));
    assert(timer.due(false, start + std::chrono::seconds(61)));
    inv_list_progress_timer disabled(false, start);
    assert(!disabled.due(true, start + std::chrono::hours(1)));

    using pool_type = bm::thread_pool<bm::task_descr*, bm::spin_lock<bm::pad0_struct>>;
    pool_type pool;
    bm::task_batch<bm::bvector<>::allocator_type> batch;
    bm::task_function_t task([](void*) { return 0; });
    batch.add(task, nullptr);
    batch.add(task, nullptr);
    assert(inv_list_completed_tasks(batch, 12) == 10);
    std::thread worker([&]()
    {
        batch.get_task(1)->done.store(true, std::memory_order_release);
        std::this_thread::sleep_for(std::chrono::milliseconds(250));
        batch.get_task(0)->done.store(true, std::memory_order_release);
    });
    unsigned checks = 0;
    bool saw_later_completion = false;
    inv_list_wait_for_task(pool, batch, 0, true, [&]()
    {
        ++checks;
        saw_later_completion |= inv_list_completed_tasks(batch, 12) == 11;
    });
    worker.join();
    assert(checks >= 2 && saw_later_completion);
    assert(inv_list_completed_tasks(batch, 12) == 12);
    inv_list_wait_for_task(pool, batch, 0, false, [&]() { assert(false); });
    if (argc > 1)
    {
        batch.get_task(0)->done.store(false, std::memory_order_release);
        inv_list_progress_timer realtime(true);
        const auto before = clock_type::now();
        unsigned reports = 0;
        std::thread slow_worker([&]()
        {
            std::this_thread::sleep_for(std::chrono::seconds(31));
            batch.get_task(0)->done.store(true, std::memory_order_release);
        });
        inv_list_wait_for_task(pool, batch, 0, true, [&]()
        {
            if (realtime.due(false))
            {
                const double seconds = std::chrono::duration<double>(clock_type::now() - before).count();
                assert(seconds >= 30);
                ++reports;
                realtime.reported();
                std::cout << "Real deadline report during worker wait: " << seconds << " seconds" << std::endl;
            }
        });
        slow_worker.join();
        assert(reports >= 1);
    }
    std::cout << "30-second deadline/count reset, disabled reporting, ordered-wait progress and completion snapshots PASS\n";
}
