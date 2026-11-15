#ifndef INV_LIST_WAVE_SCORES_H_INCLUDED
#define INV_LIST_WAVE_SCORES_H_INCLUDED
/** @file inv_list_wave_scores.h
    @brief Optional exact wave-selection TSV export, separate from codec logic.
    @internal Models are borrowed observations. A mutex protects the presentation
    sink shared by compression workers. No per-record history is retained in RAM.
*/
#include <algorithm>
#include <fstream>
#include <mutex>
#include <stdexcept>
#include <string>
#include "bmsblockhier.h"

/** @brief Thread-safe optional search-score presentation sink. */
class inv_wave_score_report
{
    std::ofstream output_; ///< Buffered output owned by this report.
    std::mutex mutex_; ///< Serializes complete seven-row groups.
    unsigned long long groups_ = 0; ///< Stable group IDs within this export.
public:
    /** Open a fresh report destination.
        @param path Output file path, exclusively owned for this process.
    */
    void open(const std::string& path);
    /// Flush buffered output and report write failures before successful exit.
    void finish();
    /** Export one search observation synchronously.
        @param context Borrowed inv_wave_score_report pointer.
        @param data Borrowed observation valid only for this call.
    */
    static void observe(void* context, const bm::sblock_wave_score_observation& data);
};

inline
void inv_wave_score_report::open(const std::string& path)
{
    output_.open(path);
    if (!output_)
    {
        throw std::runtime_error("Cannot open wave-score report: " + path);
    }
    output_ << "group\tscope\tsb\tn\tblocks\trival\tflat_bytes\thier_bytes"
        "\tflat_width\thier_width\tflat_mode\thier_mode\tflat_map_bits"
        "\thier_map_bits\tflat_header_bits\thier_header_bits"
        "\tflat_payload_bits\thier_payload_bits\tflat_padding_bits\thier_padding_bits\tflat_waves\thier_waves\thier_map_mode\n";
}

inline
void inv_wave_score_report::finish()
{
    if (output_.is_open())
    {
        output_.flush();
        if (!output_)
        {
            throw std::runtime_error("Error flushing wave-score report");
        }
    }
}

inline
void inv_wave_score_report::observe(void* context,
                                    const bm::sblock_wave_score_observation& data)
{
    auto& report = *static_cast<inv_wave_score_report*>(context);
    std::lock_guard<std::mutex> lock(report.mutex_);
    auto row = [&](unsigned scope, const bm::sblock_wave_score& flat,
                   const bm::sblock_wave_score& hier)
    {
        const bool valid = flat.bytes != ~bm::id64_t(0);
        report.output_ << report.groups_ << '\t' << scope << '\t' << data.sb
            << '\t' << data.n << '\t' << data.blocks << '\t' << data.incumbent
            << '\t' << (valid ? flat.bytes : 0) << '\t' << hier.bytes
            << '\t' << (128u << flat.width) << '\t' << (128u << hier.width)
            << '\t' << flat.mode << '\t' << hier.mode
            << '\t' << flat.map_bits << '\t' << hier.map_bits
            << '\t' << flat.bits - flat.map_bits - flat.payload_bits
            << '\t' << hier.bits - hier.map_bits - hier.payload_bits
            << '\t' << flat.payload_bits << '\t' << hier.payload_bits
            << '\t' << (valid ? (flat.bytes-1)*8-flat.bits : 0)
            << '\t' << (hier.bytes-1)*8-hier.bits
            << '\t' << flat.occupied << '\t' << hier.occupied << '\t' << hier.map_mode << '\n';
    };
    for (unsigned w = 0; w < 6; ++w)
    {
        row(w, data.flat[w], data.hier[w]);
    } // for w
    row(6, data.best_flat, data.best_hier);
    ++report.groups_;
    if (!report.output_)
    {
        throw std::runtime_error("Error writing wave-score report");
    }
}
/** @brief Aggregate exact size regret by inexpensive selector leaf. */
class inv_wave_leaf_report
{
    /** @brief Aggregate observations and regret for one selector leaf. */
    struct bucket
    {
        unsigned long long groups = 0; ///< Number of observed groups.
        unsigned long long losses = 0; ///< Groups with nonzero regret.
        unsigned long long bytes = 0; ///< Total regret in bytes.
        unsigned long long maximum = 0; ///< Largest regret in bytes.
    };
    bucket leaves_[7]; ///< Aggregates indexed by selector leaf.
    std::mutex mutex_; ///< Protects aggregates from compression workers.
    std::string path_; ///< Output destination, or empty when disabled.
public:
    /** Set the output destination to be written by finish().
        @param path Output file path, exclusively owned for this process.
    */
    void open(const std::string& path) { path_ = path; }
    /** Aggregate one synchronous search observation.
        @param context Borrowed inv_wave_leaf_report pointer.
        @param row Borrowed observation valid only for this call.
    */
    static void observe(void* context, const bm::sblock_wave_score_observation& row);
    /// Write the aggregate report and report output failures before successful exit.
    void finish();
};

inline
void inv_wave_leaf_report::observe(void* context,
                                   const bm::sblock_wave_score_observation& row)
{
    auto& report = *static_cast<inv_wave_leaf_report*>(context);
    std::lock_guard<std::mutex> lock(report.mutex_);
    auto& b = report.leaves_[row.tree_leaf];
    ++b.groups;
    b.losses += row.tree_regret != 0;
    b.bytes += row.tree_regret;
    b.maximum = std::max(b.maximum, static_cast<unsigned long long>(row.tree_regret));
}

inline
void inv_wave_leaf_report::finish()
{
    if (path_.empty())
    {
        return;
    }
    std::ofstream output(path_);
    if (!output)
    {
        throw std::runtime_error("Cannot open leaf report: " + path_);
    }
    output << "leaf\tformat\tgroups\tlosses\tregret_bytes\tmax_regret\n";
    for (unsigned leaf = 0; leaf < 7; ++leaf)
    {
        const auto& b = leaves_[leaf];
        output << leaf << '\t' << ((leaf == 1 || leaf == 3 || leaf == 4 || leaf == 5) ? 109 : 108)
               << '\t' << b.groups << '\t' << b.losses << '\t' << b.bytes
               << '\t' << b.maximum << '\n';
    } // for leaf
    if (!output)
    {
        throw std::runtime_error("Error writing leaf report: " + path_);
    }
}

#endif
