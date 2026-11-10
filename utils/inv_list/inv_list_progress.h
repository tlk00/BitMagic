#ifndef INV_LIST_PROGRESS_H_INCLUDED
#define INV_LIST_PROGRESS_H_INCLUDED

/** @file inv_list_progress.h
    @brief Internal main-thread progress timing and bounded worker waits.
    Reports are best effort: blocking filesystem operations cannot be interrupted.
*/
#include <chrono>
#include <thread>
#include <stdexcept>
#include "bm.h"
#include "bmthreadpool.h"

/** @internal Main-thread timer; vector-count triggers supplement a 30-second deadline. */
class inv_list_progress_timer
{
public:
    using clock_type = std::chrono::steady_clock;
    using time_point = clock_type::time_point;

    /** Initialize reporting state.
        @param enabled Whether progress was requested by the user.
        @param now Initial monotonic time; explicit values permit deterministic tests.
    */
    explicit inv_list_progress_timer(bool enabled, time_point now = clock_type::now());

    /** Test whether a report is due without changing the timer.
        @param count_due A vector-count or final-report trigger.
        @param now Current monotonic time.
        @return True when enabled and either the trigger or 30-second deadline is met.
    */
    bool due(bool count_due, time_point now = clock_type::now()) const;

    /** Record an emitted report, restarting its deadline.
        @param now Time at which the report was emitted.
    */
    void reported(time_point now = clock_type::now());

private:
    bool enabled_;        ///< False preserves explicitly disabled progress.
    time_point last_;     ///< Time of the most recent report, or phase start.
};

/** @internal Count completed tasks without reading worker-owned payloads.
    @param batch Current batch, whose task descriptors remain alive during the call.
    @param prepared Total dispatched tasks, including previous drained batches.
    @return Previous batches plus tasks whose completion flag has been published.
*/
template<class Batch>
bm::id64_t inv_list_completed_tasks(Batch& batch, bm::id64_t prepared);

/** @internal Wait for an ordered task while allowing main-thread progress.
    @param pool Worker pool; the caller owns it and joins workers on failure.
    @param batch Current batch, which must outlive the wait.
    @param index Index of the task to wait for.
    @param enabled Enable periodic progress checks; otherwise retain the ordinary wait.
    @param report Main-thread callback; it decides whether the 30-second report is due.
*/
template<class Pool, class Batch, class Report>
void inv_list_wait_for_task(Pool& pool, Batch& batch, unsigned index,
                            bool enabled, Report report);

inline
inv_list_progress_timer::inv_list_progress_timer(bool enabled, time_point now)
    : enabled_(enabled), last_(now)
{}

inline
bool inv_list_progress_timer::due(bool count_due, time_point now) const
{
    return enabled_ && (count_due || now - last_ >= std::chrono::seconds(30));
}

inline
void inv_list_progress_timer::reported(time_point now)
{
    last_ = now;
}

template<class Batch>
bm::id64_t inv_list_completed_tasks(Batch& batch, bm::id64_t prepared)
{
    if (prepared < batch.size())
    {
        throw std::logic_error("Progress batch exceeds dispatched task count");
    }
    bm::id64_t completed = prepared - batch.size();
    for (unsigned i = 0; i < batch.size(); ++i)
    {
        completed += batch.get_task(i)->done.load(std::memory_order_acquire);
    } // for i
    return completed;
}

template<class Pool, class Batch, class Report>
void inv_list_wait_for_task(Pool& pool, Batch& batch, unsigned index,
                            bool enabled, Report report)
{
    if (!enabled)
    {
        bm::thread_pool_executor<Pool>::wait_for_batch_done(pool, batch, index, index);
        return;
    }
    const bm::task_descr* task = batch.get_task(index);
    auto next_check = std::chrono::steady_clock::now();
    unsigned spins = 0;
    while (!task->done.load(std::memory_order_acquire))
    {
        if (pool.is_stopped() == Pool::stop_now)
        {
            throw std::runtime_error("Worker pool stopped before task completion");
        }
        // Preserve the executor's yielding wait and prompt completion response.
        // Sample the clock sparsely and invoke the reporting callback at most
        // ten times per second, independent of how fast the spin loop runs.
        if (!(spins++ & 63u))
        {
            const auto now = std::chrono::steady_clock::now();
            if (now >= next_check)
            {
                report();
                next_check = now + std::chrono::milliseconds(100);
            }
        }
        std::this_thread::yield();
    } // while task pending
}

#endif
