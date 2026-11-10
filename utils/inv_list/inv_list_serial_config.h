#ifndef INV_LIST_SERIAL_CONFIG_H_INCLUDED
#define INV_LIST_SERIAL_CONFIG_H_INCLUDED
/*
Copyright(c) 2026 Anatoliy Kuznetsov
Licensed under the Apache License, Version 2.0.
https://www.apache.org/licenses/LICENSE-2.0
*/

/** @file inv_list_serial_config.h
    @brief Application serializer settings shared by serial and threaded inv_list.

    Retained codec experiments supply a configuration callback. Settings are
    passed explicitly; this header has no dependency on CLI globals or reporting.
    This is an application utility, not a public BitMagic library interface.
*/

#include "bm.h"
#include "bmserial.h"

/** @brief Snapshot of preset settings and explicit serializer overrides.
    @internal
    @ingroup bvserial
*/
struct inv_list_serial_options
{
    unsigned level = bm::set_compression_default; ///< Library compression preset.
    unsigned gap_tolerance = 100; ///< Sparse-superblock GAP grouping tolerance.
    unsigned inverse_limit = 196607; ///< Maximum inverse-list search cardinality.
    unsigned bwt_min_gain = 0; ///< Explicit minimum complete-record BWT saving.
    bool bwt_gain_specified = false; ///< Apply bwt_min_gain to BWT presets.
    bool block_endpoints = false; ///< Enable complete-size block endpoint selection.
    bool sblock_correction = false; ///< Enable global superblock correction selection.
    bool sblock_large = false; ///< Enable large-count superblocks explicitly.
    bool wave_fixed = false; ///< Enable fixed-width superblock wave compaction.
    bool wave_adaptive = false; ///< Enable adaptive-width superblock wave compaction.
    bool bwt_shared = false; ///< Share BWT across competing block encodings.
};

/** @brief Apply a preset, experiment policy, then explicit serializer overrides.
    @internal
    @ingroup bvserial
    @tparam ConfigurePolicy Callable accepting the destination serializer by reference.
    @param serializer Exclusively owned destination; existing settings are reset by the preset.
    @param options Read-only configuration snapshot; no reference is retained.
    @param configure_policy Synchronous callback applying experiment settings after the preset.
    @pre Options and callback satisfy the serializer setters' documented constraints.
    @note The callback and any referenced state must remain valid through this call.
*/
template<class ConfigurePolicy>
void configure_inv_serializer(bm::serializer<bm::bvector<> >& serializer,
                              const inv_list_serial_options& options,
                              ConfigurePolicy configure_policy);

template<class ConfigurePolicy>
void configure_inv_serializer(bm::serializer<bm::bvector<> >& serializer,
                              const inv_list_serial_options& options,
                              ConfigurePolicy configure_policy)
{
    serializer.byte_order_serialization(false);
    serializer.gap_length_serialization(false);
    serializer.set_compression_level(options.level);
    serializer.set_sblock_gap_tolerance(options.gap_tolerance);
    configure_policy(serializer);
    if (options.level >= 9 && options.bwt_gain_specified)
    {
        serializer.set_block_selection_bwt_min_gain(options.bwt_min_gain);
    }
    if (options.block_endpoints)
    {
        serializer.set_block_endpoint_selection(true);
    }
    if (options.sblock_correction)
    {
        serializer.set_sblock_correction_selection(true);
    }
    if (options.sblock_large)
    {
        serializer.set_sblock_large_options(true);
    }
    serializer.set_sblock_large_inverse_limit(options.inverse_limit);
    if (options.wave_fixed)
    {
        serializer.set_sblock_wave1024_options(true);
    }
    if (options.wave_adaptive)
    {
        serializer.set_sblock_wave_adaptive_options(true);
    }
    if (options.bwt_shared)
    {
        serializer.set_bwt_shared_options(true);
    }
}

#endif
