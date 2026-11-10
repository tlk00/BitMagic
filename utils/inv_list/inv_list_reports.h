#ifndef INV_LIST_REPORTS_H_INCLUDED
#define INV_LIST_REPORTS_H_INCLUDED
/* Copyright(c) 2026 Anatoliy Kuznetsov
   Licensed under the Apache License, Version 2.0.
   https://www.apache.org/licenses/LICENSE-2.0 */

/** @file inv_list_reports.h
    @brief Application presentation helpers for inv_list collection progress.

    Formats caller-owned snapshots without collecting statistics, reading CLI
    globals, or controlling workers. Research-policy renderers accept statistics
    by const reference; no input references are retained. Output field names are
    preserved for the existing benchmark reporting scripts.
*/
#include <cstddef>
#include <iomanip>
#include <ostream>
#include <vector>
#include "bm.h"

/** @brief Immutable-by-convention snapshot supplied by the ordered main thread.
    @internal
*/
struct inv_progress_values
{
    bm::id64_t vectors; ///< Records written or checked in order.
    bm::id64_t prepared; ///< Dispatched records.
    bm::id64_t completed; ///< Completed records, including later finished jobs.
    std::streamoff input_bytes; ///< Raw input consumed in order.
    bm::id64_t encoded_bytes; ///< Encoded output or input consumed in order.
    double elapsed; ///< Wall-clock seconds since the operation began.
    double last_seconds; ///< Last emitted vector's compression time.
    std::size_t last_bytes; ///< Last emitted vector's encoded byte count.
};

/** Format compression progress without a newline.
    @param out Caller-owned output stream; formatting flags may change.
    @param values Snapshot valid for this synchronous call.
*/
inline
void inv_report_compress_progress(std::ostream& out, const inv_progress_values& values);

/** Format verification/decode progress without a newline.
    @param out Caller-owned output stream; formatting flags may change.
    @param values Snapshot valid for this synchronous call.
    @param verify True for raw-ID verification, false for independent decode.
*/
inline
void inv_report_read_progress(std::ostream& out, const inv_progress_values& values,
                              bool verify);

/** Format optional research-policy progress fields.
    @tparam Statistics Policy statistics exposing the referenced fields.
    @param out Caller-owned output stream; formatting flags may change.
    @param s Read-only statistics snapshot, valid throughout this call.
*/
template<class Statistics>
void inv_report_block_selection_step2_policy(std::ostream& out, const Statistics& s);

/** Format optional research-policy progress fields.
    @tparam Statistics Policy statistics exposing the referenced fields.
    @param out Caller-owned output stream; formatting flags may change.
    @param s Read-only statistics snapshot, valid throughout this call.
*/
template<class Statistics>
void inv_report_gap_route_compression_policy(std::ostream& out, const Statistics& s);

/** Format compression totals after stream closure and reconciliation.
    @param out Caller-owned output stream; formatting flags may change.
    @param diagnostics Read-only snapshot value for this report.
    @param prepared Read-only snapshot value for this report.
    @param emitted Read-only snapshot value for this report.
    @param peak_input_bytes Read-only snapshot value for this report.
    @param worker_seconds Read-only snapshot value for this report.
    @param head_wait_seconds Read-only snapshot value for this report.
    @param input_consumed Read-only snapshot value for this report.
    @param fsize Read-only snapshot value for this report.
    @param total_ints Read-only snapshot value for this report.
    @param total_min_delta Read-only snapshot value for this report.
    @param output_bytes Read-only snapshot value for this report.
*/
inline
void inv_report_compression_summary(std::ostream& out,
    bool diagnostics,
    bm::id64_t prepared,
    bm::id64_t emitted,
    std::size_t peak_input_bytes,
    double worker_seconds,
    double head_wait_seconds,
    std::streamoff input_consumed,
    std::streamoff fsize,
    bm::id64_t total_ints,
    bm::id64_t total_min_delta,
    bm::id64_t output_bytes);

/** Format read totals and optional superblock length distribution; retain no references.
    @param out Caller-owned output stream; formatting flags may change.
    @param verify Read-only snapshot value for this report.
    @param prepared Read-only snapshot value for this report.
    @param emitted Read-only snapshot value for this report.
    @param worker_seconds Read-only snapshot value for this report.
    @param peak_staged_bytes Read-only snapshot value for this report.
    @param show_lengths Read-only snapshot value for this report.
    @param sblock_exact Read-only snapshot value for this report.
    @param sblock_bin_count Read-only snapshot value for this report.
    @param sblock_bin_bytes Read-only snapshot value for this report.
    @param sblock_count Read-only snapshot value for this report.
    @param sblock_bytes Read-only snapshot value for this report.
    @param sblock_length_sum Read-only snapshot value for this report.
    @param sblock_min Read-only snapshot value for this report.
    @param sblock_max Read-only snapshot value for this report.
    @param input_consumed Read-only snapshot value for this report.
    @param input_size Read-only snapshot value for this report.
    @param total_ints Read-only snapshot value for this report.
*/
inline
void inv_report_read_summary(std::ostream& out,
    bool verify,
    bm::id64_t prepared,
    bm::id64_t emitted,
    double worker_seconds,
    std::size_t peak_staged_bytes,
    bool show_lengths,
    const std::vector<bm::id64_t>& sblock_exact,
    const bm::id64_t (&sblock_bin_count)[64],
    const bm::id64_t (&sblock_bin_bytes)[64],
    bm::id64_t sblock_count,
    bm::id64_t sblock_bytes,
    bm::id64_t sblock_length_sum,
    unsigned sblock_min,
    unsigned sblock_max,
    std::streamoff input_consumed,
    std::streamoff input_size,
    bm::id64_t total_ints);

inline
void inv_report_compress_progress(std::ostream& out, const inv_progress_values& v)
{
    out << "Progress compress vectors=" << v.vectors
        << " prepared=" << v.prepared << " completed=" << v.completed
        << " input_bytes=" << v.input_bytes << " output_bytes=" << v.encoded_bytes
        << " elapsed_s=" << std::fixed << std::setprecision(2) << v.elapsed
        << " vectors_per_s=" << (v.elapsed > 0 ? double(v.vectors) / v.elapsed : 0)
        << " last_vector_s=" << v.last_seconds << " last_vector_bytes=" << v.last_bytes;
}

inline
void inv_report_read_progress(std::ostream& out, const inv_progress_values& v, bool verify)
{
    out << "Progress " << (verify ? "verify" : "decode")
        << " vectors=" << v.vectors << " prepared=" << v.prepared
        << " completed=" << v.completed;
    if (verify)
    {
        out << " input_bytes=" << v.input_bytes;
    }
    out << " encoded_bytes=" << v.encoded_bytes
        << " elapsed_s=" << std::fixed << std::setprecision(2) << v.elapsed
        << " vectors_per_s=" << (v.elapsed > 0 ? double(v.vectors) / v.elapsed : 0);
}

template<class Statistics>
void inv_report_block_selection_step2_policy(std::ostream& out, const Statistics& s)
{
    out << " evaluated=" << s.blocks.size() << " saved_bytes=" << s.saved
        << " extra_saved_bytes=" << (static_cast<long long>(s.fourway_bytes) - static_cast<long long>(s.selected_bytes))
        << " step2_wins=" << s.chosen[4]
        << " step2_trials=" << s.step2_trials;
    if (s.bwt_min_gain)
    {
        out << " bwt_rejected=" << s.threshold_rejected;
    }
}

template<class Statistics>
void inv_report_gap_route_compression_policy(std::ostream& out, const Statistics& s)
{
    bm::id64_t direct = 0, grouped = 0;
    for (unsigned i = 0; i < 4; ++i)
    {
        direct += s.direct[i].blocks;
        grouped += s.grouped[i].blocks;
    } // for i
    out << " superblocks=" << s.superblocks
        << " grouped_gap=" << grouped << " direct_gap=" << direct;
}

inline
void inv_report_compression_summary(std::ostream& out,
    bool diagnostics,
    bm::id64_t prepared,
    bm::id64_t emitted,
    std::size_t peak_input_bytes,
    double worker_seconds,
    double head_wait_seconds,
    std::streamoff input_consumed,
    std::streamoff fsize,
    bm::id64_t total_ints,
    bm::id64_t total_min_delta,
    bm::id64_t output_bytes)
{
    if (diagnostics)
    {
        out << "Parallel prepared=" << prepared
         << " emitted=" << emitted
         << " peak_batch_input_bytes=" << peak_input_bytes
         << " worker_seconds=" << worker_seconds
         << " ordered_wait_seconds=" << head_wait_seconds << std::endl;
    }
    out << "\nTotal vectors=" << emitted << std::endl;
    out << "Input consumed bytes=" << input_consumed << std::endl;
    out << "Input file bytes=" << fsize << std::endl;
    out << "Total ints=" << total_ints << std::endl;
    out << "  min-deltas = " << total_min_delta << std::endl;
    out << "  min delta ratio = " << std::setprecision(3)
         << double(total_min_delta) / double(total_ints) << std::endl;
    out << "BV size = " << output_bytes << std::endl;
    out << "BV Bits per/int = " << std::setprecision(3)
         << double(output_bytes * 8ull) / double(total_ints) << std::endl;
}

inline
void inv_report_read_summary(std::ostream& out,
    bool verify,
    bm::id64_t prepared,
    bm::id64_t emitted,
    double worker_seconds,
    std::size_t peak_staged_bytes,
    bool show_lengths,
    const std::vector<bm::id64_t>& sblock_exact,
    const bm::id64_t (&sblock_bin_count)[64],
    const bm::id64_t (&sblock_bin_bytes)[64],
    bm::id64_t sblock_count,
    bm::id64_t sblock_bytes,
    bm::id64_t sblock_length_sum,
    unsigned sblock_min,
    unsigned sblock_max,
    std::streamoff input_consumed,
    std::streamoff input_size,
    bm::id64_t total_ints)
{
    out << (verify ? "Verification complete." : "Decode complete.") << std::endl;
    out << "Parallel read prepared=" << prepared << " emitted=" << emitted
         << " worker_seconds=" << worker_seconds
         << " peak_staged_input_bytes=" << peak_staged_bytes << std::endl;
    out << "Total vectors=" << emitted << std::endl;
    if (show_lengths)
    {
        bm::id64_t cumulative = 0;
        const bm::id64_t rank75 = (sblock_count * 3 + 3) / 4;
        unsigned p75 = 0;
        for (unsigned len = 1; len < sblock_exact.size(); ++len)
        {
            cumulative += sblock_exact[len];
            if (cumulative >= rank75)
            {
                p75 = len;
                break;
            }
        } // for len
        out << "SBlock BIC summary records=" << sblock_count
             << " bytes=" << sblock_bytes
             << " min=" << (sblock_count ? sblock_min : 0)
             << " max=" << sblock_max
             << " mean=" << std::fixed << std::setprecision(3)
             << (sblock_count ? double(sblock_length_sum) / double(sblock_count) : 0.0)
             << " p75=" << p75 << std::endl;
        for (unsigned bin = 0; bin < 64; ++bin)
        {
            unsigned lo, hi;
            if (bin < 16)
            {
                lo = hi = bin + 1;
            }
            else
            {
                const unsigned e = 4 + (bin - 16) / 4;
                const unsigned q = (bin - 16) % 4;
                const unsigned base = 1u << e;
                lo = base + q * base / 4;
                hi = base + (q + 1) * base / 4 - 1;
                if (bin == 16)
                {
                    lo = 17;
                }
            }
            out << "SBlock BIC bin=" << bin
                 << " lo=" << lo << " hi=" << hi
                 << " records=" << sblock_bin_count[bin]
                 << " bytes=" << sblock_bin_bytes[bin] << std::endl;
        } // for bin
    }
    if (verify)
    {
        out << "Input consumed bytes=" << input_consumed << std::endl;
        out << "Input file bytes=" << input_size << std::endl;
        out << "Total ints=" << total_ints << std::endl;
    }
}

#endif
