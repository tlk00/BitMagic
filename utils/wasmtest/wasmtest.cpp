#include <iostream>
#include <thread>
#include <chrono>
#include <vector>
#include <list>

using namespace std;

#include <unistd.h>
#include <pthread.h>
#include <assert.h>

#include "bmsync.h"

#ifdef __EMSCRIPTEN__
#include <emscripten/threading.h>
#endif



static atomic<unsigned> c;


extern "C"
{
    void* myPThreadFunc(void* /*argp*/)
    {
        ++c;
        //std::this_thread::sleep_for(std::chrono::milliseconds(100));
        return NULL;
    }
}

using namespace std;

bm::spin_lock<> lk;

struct Two_Spins
{
    bm::spin_lock<bm::pad60_struct> l1;
    bm::spin_lock<bm::pad60_struct> l2;
};

int main(void)
{
    c = 0;

    cout << "Sizeof spin_lock: " << sizeof(bm::spin_lock<>) << endl;
    cout << "Sizeof 2 spin_locks: " << sizeof(Two_Spins) << endl;


    bm::queue_sync<unsigned, std::mutex> qu_sync;
    qu_sync.push(1);
    assert(!qu_sync.empty());
    unsigned v(0);
    bool b = qu_sync.try_pop(v);
    assert(b);
    assert(v == 1);
    assert(qu_sync.empty());

    std::atomic<bool> locked{false};

    bool was_locked = locked.load(std::memory_order_relaxed);

    std::cout << "Testing threads.." << endl;

    #ifdef __EMSCRIPTEN__
        if (emscripten_has_threading_support())
        {
            auto cores = emscripten_num_logical_cores();
            cout << "Threading enabled. Cores = " << cores << endl;
        }
        else
        {
            return 1;
        }
    #endif


    {
        locked.compare_exchange_weak(was_locked,
                                     true, std::memory_order_acquire);

//        typedef
//        bm::thread_pool<bm::task_descr, bm::spin_lock<bm::pad60_struct> > pool_type;

        typedef
        bm::thread_pool<bm::task_descr, bm::spin_lock<bm::pad0_struct> > pool_type;

//        typedef
//        bm::thread_pool<bm::task_descr, std::mutex> pool_type;

        auto time_before = std::chrono::high_resolution_clock::now();

        pool_type tpool;  // our thread pool here (no threads created yet)
        pool_type::queue_type& qu = tpool.get_job_queue();

        // push task descriptors in non blocking mode because threads are
        // not running yet
        //
        for (unsigned i = 0; i < 1000; ++i)
            qu.push_no_lock(bm::task_descr(myPThreadFunc));

        tpool.start(3); // start the processing

        tpool.wait_empty_queue();
        cout << "after wait_empty_queue()=" << c.load(std::memory_order_relaxed) << endl;

        // keep pushing jobs
        //
        for (unsigned i = 0; i < 1000; ++i)
            qu.push(bm::task_descr(myPThreadFunc));

        tpool.wait_empty_queue();

        // submit next batch at once then unlock processing using scope lock
        {
            std::lock_guard<pool_type::queue_type> scope_lck(qu);
            for (unsigned i = 0; i < 1000; ++i)
                qu.push_no_lock(bm::task_descr(myPThreadFunc));
        }

        tpool.set_stop_mode(pool_type::stop_when_done);
        tpool.join();

        auto time_now = std::chrono::high_resolution_clock::now();
        auto diff = time_now - time_before;
        auto ms = std::chrono::duration <double, std::milli> (diff).count();
        std::cout << "Total: " << ms << " ms" << std::endl;


         cout << "after join=" << c.load(std::memory_order_relaxed) << endl;
    }

    cout << was_locked << endl;
    cout << "end" << endl;
    return 0;
}
