#include <iostream>
#include <thread>
#include <chrono>
#include <vector>
#include <string>
#include <memory>
#include <list>
#include <cassert>

using namespace std;

#include <unistd.h>
#include <pthread.h>
#include <assert.h>
#include <math.h>

#include "bm.h"
#include "bmthreadpool.h"

#ifdef __EMSCRIPTEN__
#include <emscripten/threading.h>
#include <emscripten/emscripten.h>
#include <emscripten/bind.h>

using namespace emscripten;
#endif
using namespace std;





#if 0
extern "C"
{
    void* myPThreadFunc(void* /*argp*/)
    {
        ++c;
        //std::this_thread::sleep_for(std::chrono::milliseconds(100));
        return NULL;
    }
}
#endif


/// Test method for simple Fibonacci
///

extern "C"
{
EMSCRIPTEN_KEEPALIVE
int my_fib(int x)
{
  if (x < 1)
    return 0;
  if (x == 1)
    return 1;
  return my_fib(x-1) + my_fib(x-2);
}
}


int my_sq(int x) noexcept
{
    return x * x;
}

class TestVector_proxy;

/// Test class to try various bindings
///
class TestVector
{
public:
    TestVector() noexcept {}

    unsigned get_value() const noexcept
    { return 26; }

    void set_str(const string& s)
    {
        str_ = s;
    }

    std::string get_hello() const noexcept
    { return str_; }

    TestVector_proxy* createProxy();
private:
    std::string str_;
};

class TestVector_proxy
{
public:
    TestVector_proxy(TestVector* vec)
    : vec_ptr_(vec)
    {
        vect_buf_.resize(5);
        for (size_t i = 0; i < vect_buf_.size(); ++i)
            vect_buf_[i] = 100;
    }
    void set(const std::string& s) { vec_ptr_->set_str(s); }

    int getBufPtr() noexcept
    {
        std::vector<unsigned>::value_type* buf_ptr = vect_buf_.data();
        return (int) buf_ptr;
    }

    int getBufSize() noexcept
    {
        return (int) vect_buf_.size();
    }

    unsigned sum() const noexcept
    {
        unsigned s = 0;
        for (size_t i = 0; i < vect_buf_.size(); ++i)
            s += vect_buf_[i];
        return s;
    }

private:
    TestVector* vec_ptr_;
    std::vector<unsigned> vect_buf_;
};


TestVector_proxy* TestVector::createProxy()
{
    return new TestVector_proxy(this);
}

// ------------------------------------------------------------------
// Thread management
//

typedef
bm::thread_pool<bm::task_description*, bm::spin_lock<bm::pad0_struct> > thread_pool_type;
typedef bm::thread_pool_executor<thread_pool_type> thread_pool_exec_type;

thread_pool_type  g_tpool; ///< global thread pool


extern "C"
{

/**
    Initialize the thread pool
 */
EMSCRIPTEN_KEEPALIVE
void init_thread_pool(int pool_size)
{
    g_tpool.start(pool_size);
}

EMSCRIPTEN_KEEPALIVE
void stop_thread_pool()
{
    g_tpool.stop();
}


} // extern C


typedef int JS_heap_ptr_type; ///< JavaScript heap offset


/// Collection of floating point vectors (memory buffers) for numerical compute
/// available from the JS side
///
class FloatVectorStore
{
public:
    typedef std::vector<float>                          float_vector_type;
    typedef std::vector<unique_ptr<float_vector_type> > float_vector_ptr_type;
public:

    /**
        Reset storage
    */
    void reset() { fp_vectors_.resize(0); }

    /**
        Add a new vector to storage, return JS ready 32-bit heap pointer
        @param vect_size - target vector size
        @return pointer to created bugger
     */
    JS_heap_ptr_type add_vector_js(size_t vect_size)
    {
        fp_vectors_.emplace_back(new float_vector_type(vect_size));
        float* p = fp_vectors_.back()->data();
        return (JS_heap_ptr_type) p;
    }

    JS_heap_ptr_type get_ptr(size_t idx)
    {
        float* p = fp_vectors_.at(idx)->data();
        return (JS_heap_ptr_type) p;
    }

    size_t get_size(size_t idx)
    {
        return fp_vectors_.at(idx)->size();
    }


    int size() const noexcept
    {
        return fp_vectors_.size();
    }

    // test function to check sanity
    float sum(size_t idx) const
    {
        float s = 0.0;
        const float_vector_type& fp_vect = *(fp_vectors_[idx].get());
        for (size_t i = 0; i < fp_vect.size(); ++i)
            s += fp_vect[i];
        return s;
    }

protected:
    float_vector_ptr_type   fp_vectors_; ///< floating point vectors
};


extern "C"
{

/**
    Initialize the thread pool
 */
EMSCRIPTEN_KEEPALIVE
FloatVectorStore* create_FloatStore()
{
    return new FloatVectorStore();
}

/// Calculate mean by JS heap pointer
///
EMSCRIPTEN_KEEPALIVE
float fp_mean(JS_heap_ptr_type js_arr, unsigned size) noexcept
{
    const float* arr = reinterpret_cast<float*>(js_arr);
    float s(0.0);
    for (size_t i = 0; i < size; ++i)
        s += arr[i];
    return s / size;
}

EMSCRIPTEN_KEEPALIVE
float fp_pearson(JS_heap_ptr_type js_arr_x,
                 JS_heap_ptr_type js_arr_y,  unsigned size) noexcept
{
float res;
//for (unsigned pass = 0; pass < 10; ++pass)
{
    const float* arr_x = reinterpret_cast<float*>(js_arr_x);
    const float* arr_y = reinterpret_cast<float*>(js_arr_y);

    unsigned n = 0;
    double xMean(0.0), yMean(0.0), numerator(0.0), d(0.0), e(0.0);
    for (size_t i = 0; i < size; ++i)
    {
        float x = arr_x[i];
        float y = arr_y[i];
        ++n;
        xMean = (n * xMean + x) / (n + 1);
        yMean = (n * yMean + y) / (n + 1);
        numerator = numerator + (x - xMean) * (y - yMean);
        d = d + (x - xMean) * (x - xMean);
        e = e + (y - yMean) * (y - yMean);
    } // for i

    res = (n==0) ? 0.0 : numerator / (sqrt(d) * sqrt(e));
}
    return res;
}

}

// ----------------------------------------------------------------

/// Batch builder for parallel MT processing
///
///
class PearsonsParallelBuilder
{
public:
    typedef std::mutex                                lock_type;
    typedef typename bm::bvector<>                    bvector_type;
    typedef typename bvector_type::allocator_type     allocator_type;
public:
    PearsonsParallelBuilder() {}


    struct task_batch : public bm::task_batch<allocator_type>
    {
        typedef bm::task_batch<allocator_type>          parent_type;
        typedef typename parent_type::task_vector_type  task_vector_type;
    };

    void add_to_batch(task_batch& batch,
                      JS_heap_ptr_type js_arr_x,
                      JS_heap_ptr_type js_arr_y,  unsigned arr_size)
    {
        typename task_batch::task_vector_type& tv = batch.get_task_vector();
        bm::task_description& tdescr = tv.add();
        tdescr.init(task_run_func, (void*)&tdescr,
                    (void*)js_arr_x, (void*)js_arr_y, arr_size);
    }

    /// Add job to a batch embedded with the parallel builder
    ///
    void add_to_default_batch(JS_heap_ptr_type js_arr_x,
                              JS_heap_ptr_type js_arr_y,
                              unsigned arr_size)
    {
        add_to_batch(default_batch_, js_arr_x, js_arr_y, arr_size);
    }

    /// Reset the deafult batch
    void reset_default_batch()
    {
        default_batch_.get_task_vector().resize(0);
    }

    /// Run batch and wait for its  execution on thread pool
    ///
    void sync_run_default_batch()
    {
        thread_pool_exec_type exec;
        exec.run(g_tpool, default_batch_, true);
    }

    float get_pcoeff_for_task(int idx) const
    {
        const typename task_batch::task_vector_type& tv =
                                    default_batch_.get_task_vector();
        const bm::task_description& tdescr = tv.at(idx);
        return tdescr.payload.fp32;
    }


protected:

    /// Task execution Entry Point Function
    /// @internal
    static void* task_run_func(void* argp)
    {
        if (!argp)
            return 0;


  //std::thread::id this_id = std::this_thread::get_id();
  //std::cout << "thread: " << this_id << endl;;
  //printf("%i\n", this_id);

        bm::task_description* tdescr = (bm::task_description*) argp;
        JS_heap_ptr_type js_arr_x = (JS_heap_ptr_type) tdescr->ctx0;
        JS_heap_ptr_type js_arr_y = (JS_heap_ptr_type) tdescr->ctx1;

        float pcoeff =
                fp_pearson(js_arr_x, js_arr_y, unsigned(tdescr->param0));
        tdescr->payload.fp32 = pcoeff;
        return nullptr;
    }

private:
    task_batch              default_batch_; // default batch
};



// ------------------------------------------------------------------
// Bindings
//


EMSCRIPTEN_BINDINGS(BitMagic) {
    emscripten::function("init_thread_pool", &init_thread_pool);
    emscripten::function("stop_thread_pool", &stop_thread_pool);
    emscripten::function("create_FloatStore", &create_FloatStore, allow_raw_pointers());
    emscripten::function("my_sq", &my_sq);
    emscripten::function("fp_mean", &fp_mean, allow_raw_pointers());
    emscripten::function("fp_pearson", &fp_pearson, allow_raw_pointers());

//    emscripten::function("my_fib", &my_fib);

emscripten::class_<PearsonsParallelBuilder>("PearsonsParallelBuilder")
	.constructor<>()
    .function("add_to_default_batch", &PearsonsParallelBuilder::add_to_default_batch)
    .function("reset_default_batch", &PearsonsParallelBuilder::reset_default_batch)
    .function("sync_run_default_batch", &PearsonsParallelBuilder::sync_run_default_batch)
    .function("get_pcoeff_for_task", &PearsonsParallelBuilder::get_pcoeff_for_task)
;

emscripten::class_<FloatVectorStore>("FloatVectorStore")
    .function("add_vector_js", &FloatVectorStore::add_vector_js)
    .function("reset", &FloatVectorStore::reset)
    .function("size", &FloatVectorStore::size)
    .function("get_ptr", &FloatVectorStore::get_ptr)
    .function("get_size", &FloatVectorStore::get_size)
    .function("sum", &FloatVectorStore::sum)
;


emscripten::class_<TestVector_proxy>("TestVector_proxy")
    .function("set", &TestVector_proxy::set)
    .function("get_buf_ptr", &TestVector_proxy::getBufPtr)
    .function("get_buf_size", &TestVector_proxy::getBufSize)
    .function("sum", &TestVector_proxy::sum)
    ;


emscripten::class_<TestVector>("TestVector")
	.constructor<>()
	.function("get_value", &TestVector::get_value)
    .function("set_str", &TestVector::set_str)
    .function("get_hello", &TestVector::get_hello)
    .function("createProxy", &TestVector::createProxy, allow_raw_pointers())
    ;
}




