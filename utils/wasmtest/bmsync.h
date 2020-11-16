#pragma once

#include <queue>
#include <mutex>
#include <condition_variable>

namespace bm
{

/// Pad 60 bytes so that the final  ocupiles 64 bytes (1 cache line)
struct pad60_struct { char c[60]; };
/// Empty padding
struct pad0_struct {  };

/**
    Spin-lock with two-phase acquire (read + cas)
    padding parameter optionally adds a buffer to avoid CPU cache
    line contention.
    TODO: test if padding realy helps in our case
 */
template<class Pad = bm::pad0_struct>
class spin_lock
{
public:
    explicit spin_lock() {}

    /// Lock the lock
    void lock() noexcept
    {
        while(1) // spin loop
        {
            unsigned locked = locked_.load(std::memory_order_relaxed);
            if (!locked &&
                 locked_.compare_exchange_weak(locked, true,
                                               std::memory_order_acquire))
                break;
#if defined(BMSSE2OPT) || defined(BMSSE42OPT) || defined(BMAVX2OPT) || defined(BMAVX512OPT)
            _mm_pause();
#endif
        } // while
    }

    /// Unlock the lock
    void unlock() noexcept
    {
        locked_.store(false, std::memory_order_release);
    }
private:
    spin_lock(const spin_lock&)=delete;
    spin_lock& operator=(const spin_lock&)=delete;

private:
    std::atomic<unsigned> locked_{false};
    Pad p_;
};

/// Wait for multiple threads to exit
///
template<typename TCont>
void join_multiple_threads(TCont& tcont)
{
    typename TCont::iterator it_end = tcont.end();
    for (typename TCont::iterator it = tcont.begin(); it != it_end; ++it)
    {
        if (it->joinable())
            it->join();
    } // for it
}



/// Typedef for a call-back function pointer
///
typedef void* (*func_call_type)(void*);

/// Task descriptor, encapsulates pthread-like callbacks with argument
/// and return
///
struct task_descr
{
    func_call_type  func;   ///< pointer to pthread-like function callback
    void*           argp;   ///< arg pointer
    void*           ret;    ///< ret pointer

    task_descr(){}
    explicit task_descr(func_call_type  f, void* argptr = 0)
        : func(f), argp(argptr), ret(0)
    {}
};


template<typename QValue, typename Lock> class thread_pool;

/**
    Thread sync queue

    TODO: use data structures which can guarantee "noexcept" for trivial ops
 */
template<typename Value, typename Lock>
class queue_sync
{
public:
    typedef Value     value_type;
    typedef Lock      lock_type;

    /// constructor
    ///
    queue_sync(){}

    /// Push value to the back of the queue
    /// @param v - value to put in the queue
    ///
    /// @sa push_no_lock
    ///
    void push(value_type v)
    {
        {
            std::lock_guard<lock_type> lg(dq_lock_);
            data_queue_.push(v);
        }
        queue_push_cond_.notify_one();
    }

    /// Push value to the back of the queue without lock protection
    /// It is assumed that processing did not start and we are just staging
    /// the batch
    ///
    /// @param v - value to put in the queue
    /// @sa push
    ///
    void push_no_lock(value_type v) { data_queue_.push(v); }


    /// Extract value
    /// @param v - [out] value returned
    /// @return true if extracted
    ///
    bool try_pop(value_type& v)
    {
        std::lock_guard<lock_type> guard(dq_lock_);
        if (data_queue_.empty())
            return false;
        v = data_queue_.front();
        data_queue_.pop();
        return true;
    }

    /// @return true if empty
    bool empty() const
    {
        std::lock_guard<lock_type> guard(dq_lock_);
        return data_queue_.empty();
    }

    /// lock the queue access
    /// @sa push_no_lock, unlock
    void lock() { dq_lock_.lock(); }

    /// unlock the queue access
    /// @sa push_no_lock, lock
    void unlock()
    {
        dq_lock_.unlock();
        // lock-unlock is done to protect bulk push, need to wake up
        // all waiting workers
        queue_push_cond_.notify_all();
    }

    template<typename QV, typename L> friend class bm::thread_pool;
private:
    queue_sync(const queue_sync&) = delete;
    queue_sync& operator=(const queue_sync&) = delete;
private:
    std::queue<value_type>    data_queue_; ///< queue object
    mutable lock_type         dq_lock_;    ///< lock for queue

    // signal structure for wait on empty queue
protected:
    mutable std::mutex        signal_mut_; ///< signal mutex for q submissions
    std::condition_variable   queue_push_cond_;   ///< mutex paired conditional


};


/**
    Thread pool with custom (thread safe) queue
    QValue - task queue value parameter
    Lock   - locking protection type (like std::mutex or spinlock)
*/
template<typename QValue, typename Lock>
class thread_pool
{
public:
    typedef QValue    value_type;
    typedef Lock      lock_type;
    typedef bm::queue_sync<QValue, lock_type> queue_type;

    /**
        Stop modes for threads:
        0 - keep running/waiting for jobs
        1 - wait for empty task queue then stop threads
        2 - stop threads now even if there are pending tasks
    */
    enum stop_mode
    {
        no_stop = 0,          ///< keep spinning on busy-wait
        stop_when_done = 1,   ///< stop if tsak queue is empty
        stop_now = 2          ///< stop right now
    };

public:
    thread_pool(stop_mode sm = no_stop)
        : stop_flag_(sm)
    {}

    ~thread_pool()
    {
        int is_stop = stop_flag_;
        if (!is_stop) // finish the outstanding jobs and close threads
            set_stop_mode(stop_when_done);
        join();
    }

    /** Setup the criteria for threads shutdown
        Also notifies all threads on a new directive
        @param sm - stop mode
    */
    void set_stop_mode(stop_mode sm)
    {
        stop_flag_ = sm;
        job_queue_.queue_push_cond_.notify_all();
    }

    /**
        Start thread pool worker threads.
        @param tcount - number of threads to start
     */
    void start(unsigned tcount)
    {
        int is_stop = stop_flag_.load(std::memory_order_relaxed);
        if (is_stop == stop_now) // immediate stop requested
            return;
        // TODO: consider lock protect of thread_vect_ member
        for(unsigned i = 0;i < tcount; ++i)
        {
            thread_vect_.emplace_back(
                            std::thread(&thread_pool::worker_func,this));
        } // for
    }

    /**
        Wait for threads to stop
    */
    void join()
    {
        bm::join_multiple_threads(thread_vect_);
        thread_vect_.resize(0);
    }

    /**
        Conditional spin-wait for the queue to empty
     */
     void wait_empty_queue()
     {
        const std::chrono::duration<int, std::milli> wait_duration(10);
        while(1)
        {
            if (job_queue_.empty())
                break;
            std::cv_status wait_res;
            {
                std::unique_lock<std::mutex> lk(task_done_mut_);
                wait_res = task_done_cond_.wait_for(lk, wait_duration);
            }
            if (wait_res == std::cv_status::timeout)
            {
#if defined(BMSSE2OPT) || defined(BMSSE42OPT) || defined(BMAVX2OPT) || defined(BMAVX512OPT)
                _mm_pause();
#else
                std::this_thread::yield();
#endif
            }
        } // while
     }

    /// Get access to the job submission queue
    ///
    queue_type& get_job_queue() noexcept { return job_queue_; }

protected:

    /// Internal worker wrapper with busy-wait spin loop
    /// making pthread-like call for tasks
    ///
    void worker_func()
    {
        const std::chrono::duration<int, std::milli> wait_duration(10);
        while(1)
        {
            int is_stop = stop_flag_.load(std::memory_order_relaxed);
            if (is_stop == stop_now) // immediate stop requested
                break;

            task_descr task_descr;
            if (job_queue_.try_pop(task_descr))
            {
                // TODO: consider try-catch here
                void* ret = task_descr.func(task_descr.argp);
                (void)ret; // TODO: add results processing

                task_done_cond_.notify_one();
                continue;
            }
            // queue appears to be empty, check if requested to stop
            //
            is_stop = stop_flag_.load(std::memory_order_relaxed);
            if (is_stop)
                return;

            // enter a temporal condition wait
            //   notifications are treated as unreliable re-verified
            //   via spin over the poll of the queue
            std::cv_status wait_res;
            {
                std::unique_lock<std::mutex> lk(job_queue_.signal_mut_);
                wait_res =
                    job_queue_.queue_push_cond_.wait_for(lk, wait_duration);
            }
            if (wait_res == std::cv_status::timeout)
            {
                is_stop = stop_flag_.load(std::memory_order_relaxed);
                if (is_stop == stop_now) // immediate stop requested
                    return;

                std::this_thread::yield();
            }
        } // while
    }


private:
    thread_pool(const thread_pool&)=delete;
    thread_pool& operator=(const thread_pool&)=delete;

private:
    queue_type                job_queue_;    ///< queue (thread sync)
    std::vector<std::thread>  thread_vect_;  ///< threads servicing queue
    std::atomic<int>          stop_flag_{0}; ///< stop flag to all threads

    // notification channel for results wait
    mutable std::mutex        task_done_mut_; ///< signal mutex for task done
    std::condition_variable   task_done_cond_;///< mutex paired conditional

};



} // bm
