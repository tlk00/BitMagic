#ifndef BM_SIMD_BENCH__H_INCLUDED__
#define BM_SIMD_BENCH__H_INCLUDED__
// Focused backend comparisons: build this same driver against native and
// translated headers, with identical compiler flags and seed.
#include <chrono>
#include <cstdlib>
#include <iomanip>
#include "bmbvimport.h"
#include "bmxor.h"

struct simd_bench_options
{
    bool kernels = false, dense = false, gap = false;
    unsigned repeats = 200, blocks = 64, seed = 1337;
    bool selected() const { return kernels || dense || gap; }
};

static volatile unsigned long long simd_bench_sink = 0;

static void SIMDBenchCheck(bool value, const char* message)
{
    if (!value) { std::cerr << "SIMD benchmark check failed: " << message << '\n'; std::exit(1); }
}

template<class Func>
static void MeasureSIMD(const char* label, unsigned repeats, unsigned units,
                        bool timings, bool report, Func fn)
{
    const unsigned long long expected = fn(); // warmup and checksum reference
    const unsigned samples = report ? 5u : 1u;
    const unsigned iterations = timings ? repeats : 1u;
    double ns[5] = {};
    for (unsigned sample = 0; sample < samples; ++sample)
    {
        unsigned long long sum = 0;
        auto start = std::chrono::steady_clock::now();
        for (unsigned r = 0; r < iterations; ++r)
        {
#if defined(__GNUC__) || defined(__clang__)
            // Prevent invariant input loads / complete sweeps being hoisted.
            __asm__ __volatile__("" : : : "memory");
#endif
            sum += fn();
        }
        auto end = std::chrono::steady_clock::now();
        simd_bench_sink = sum;
        SIMDBenchCheck(sum == expected * iterations, label);
        ns[sample] = std::chrono::duration<double, std::nano>(end-start).count()
                   / (double(iterations)*units);
    }
    if (timings)
    {
        std::sort(ns, ns+samples);
        if (report)
        {
            // Match the perf driver's established "name; duration" output.
            const double elapsed_ms = ns[2] * double(iterations) * units / 1.0e6;
            bm::chrono_taker<>::print_duration(std::cout, label, elapsed_ms);
        }
    }
}

static void RunSIMDBenchmarks(const simd_bench_options& opt, bool timings,
                              bool report = true)
{
    const unsigned N = bm::set_block_size;
    std::mt19937 rng(opt.seed);
    // bit_block_t supplies the alignment required by the translated baseline.
    std::vector<bm::bit_block_t> a(opt.blocks), b(opt.blocks);
    for (unsigned k = 0; k < opt.blocks; ++k)
        for (unsigned i = 0; i < N; ++i)
        {
            a[k].begin()[i] = unsigned(rng());
            b[k].begin()[i] = unsigned(rng());
        }
    if (opt.kernels)
    {
        auto measure = [&](const char* label, auto operation)
        {
            MeasureSIMD(label, opt.repeats, opt.blocks, timings, report, [&]()
            {
                unsigned long long sum = 0;
                for (unsigned k = 0; k < opt.blocks; ++k)
                    sum += operation(a[k].begin(), b[k].begin());
                return sum;
            });
        };
        measure("kernel/count", [](const unsigned* x, const unsigned*) { return bm::bit_block_count(x); });
        measure("kernel/count-and", [](const unsigned* x, const unsigned* y) { return bm::bit_block_and_count(x,y); });
        measure("kernel/count-or", [](const unsigned* x, const unsigned* y) { return bm::bit_block_or_count(x,y); });
        measure("kernel/count-xor", [](const unsigned* x, const unsigned* y) { return bm::bit_block_xor_count(x,y); });
        measure("kernel/count-sub", [](const unsigned* x, const unsigned* y) { return bm::bit_block_sub_count(x,y); });
        measure("kernel/runs", [](const unsigned* x, const unsigned*) { return bm::bit_block_calc_change(x); });
        measure("kernel/runs-count", [](const unsigned* x, const unsigned*)
        {
            unsigned gc, bc; bm::bit_block_change_bc(x, &gc, &bc);
            return (static_cast<unsigned long long>(gc) << 32) | bc;
        });
        measure("kernel/xor-runs-count", [](const unsigned* x, const unsigned* y)
        {
            unsigned gc, bc; bm::bit_block_xor_change(x, y, N, &gc, &bc);
            return (static_cast<unsigned long long>(gc) << 32) | bc;
        });
    }
    if (opt.dense)
    {
        std::vector<unsigned> words_a(size_t(opt.blocks)*N), words_b(words_a.size());
        for (unsigned k = 0; k < opt.blocks; ++k)
        {
            std::copy(a[k].begin(), a[k].end(), words_a.data()+size_t(k)*N);
            std::copy(b[k].begin(), b[k].end(), words_b.data()+size_t(k)*N);
        }
        bm::bvector<> av, bv;
        bm::bit_import_u32(av, words_a.data(), unsigned(words_a.size()), false);
        bm::bit_import_u32(bv, words_b.data(), unsigned(words_b.size()), false);
        bm::bvector<>::statistics st_a, st_b;
        av.calc_stat(&st_a); bv.calc_stat(&st_b);
        SIMDBenchCheck(st_a.gap_blocks == 0 && st_b.gap_blocks == 0 &&
                       st_a.bit_blocks == opt.blocks && st_b.bit_blocks == opt.blocks,
                       "dense fixtures must consist of real bit blocks");
        if (timings && report)
            std::cout << "Dense fixture: bit_blocks=" << st_a.bit_blocks
                      << ", gap_blocks=0 (per vector)\n";
        MeasureSIMD("dense/count", opt.repeats, opt.blocks, timings, report, [&]() { return av.count(); });
        MeasureSIMD("dense/count-and", opt.repeats, opt.blocks, timings, report, [&]() { return bm::count_and(av,bv); });
        MeasureSIMD("dense/count-or", opt.repeats, opt.blocks, timings, report, [&]() { return bm::count_or(av,bv); });
        MeasureSIMD("dense/count-xor", opt.repeats, opt.blocks, timings, report, [&]() { return bm::count_xor(av,bv); });
        MeasureSIMD("dense/count-sub", opt.repeats, opt.blocks, timings, report, [&]() { return bm::count_sub(av,bv); });
        MeasureSIMD("dense/intervals", opt.repeats, opt.blocks, timings, report, [&]() { return bm::count_intervals(av); });
        bm::bvector<> tmp;
        MeasureSIMD("dense/materialize-xor", opt.repeats, opt.blocks, timings, report, [&]()
        {
            tmp.bit_xor(av, bv, bm::bvector<>::opt_none);
            return tmp.count();
        });
    }
    if (opt.gap)
    {
        // Long regular runs: conversion is valid within the GAP capacity.
        std::vector<std::vector<bm::gap_word_t> > gaps(opt.blocks, std::vector<bm::gap_word_t>(1024));
        for (unsigned k = 0; k < opt.blocks; ++k)
        {
            for (unsigned i = 0; i < N; ++i) a[k].begin()[i] = ((i+k)/8)&1u ? ~0u : 0u;
            bm::bit_to_gap(gaps[k].data(), a[k].begin(), 1024);
        }
        MeasureSIMD("gap/count", opt.repeats, opt.blocks, timings, report, [&]()
        {
            unsigned long long sum = 0;
            for (unsigned k = 0; k < opt.blocks; ++k) sum += bm::gap_bit_count_unr(gaps[k].data());
            return sum;
        });
        bm::gap_word_t buffer[1024];
        MeasureSIMD("gap/convert", opt.repeats, opt.blocks, timings, report, [&]()
        {
            unsigned long long sum = 0;
            for (unsigned k = 0; k < opt.blocks; ++k)
            {
                unsigned len = bm::bit_to_gap(buffer, a[k].begin(), 1024);
                sum += len + buffer[0] + buffer[len];
            }
            return sum;
        });
    }
}
#endif
