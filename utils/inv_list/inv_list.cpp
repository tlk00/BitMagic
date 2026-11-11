/*
Copyright(c) 2002-2026 Anatoliy Kuznetsov(anatoliy_kuznetsov at yahoo.com)

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

For more information please visit:  http://bitmagic.io
*/

/** \example inv_list.cpp
  Utility to compress test sets of inverted lists (Gov2 corpus)
  and study compression characteristics and different comression levels
*/

/*! \file inv_list.cpp
    \brief Utility to compress test sets of inverted lists
*/

// This development tool consumes per-block codec histories and timings.
#ifndef BM_SERIALIZATION_DIAGNOSTICS
# define BM_SERIALIZATION_DIAGNOSTICS 1
#endif

#include <iostream>
#include <chrono>
#include <thread>
#include <time.h>
#include <stdio.h>
#include <cstdlib>
#include <algorithm>
#include <cmath>
#include <limits>
#include <iomanip>
#include <exception>
#include <sstream>
#include <array>
#include <queue>
#include <cstdint>
#include <cstring>
#include <random>
#include <atomic>
#include <functional>


#ifdef _MSC_VER
#pragma warning( push )
#pragma warning( disable : 4996)
#endif

#include <vector>
#include <map>
#include <memory>

#include "bm.h"
#include "bmalgo.h"
#include "bmcomplexity.h"
#include "bmserial.h"
#include "bmfio.h"
#include "bmsparsevec.h"
#include "bmsparsevec_compr.h"
#include "bmsparsevec_algo.h"
#include "bmsparsevec_serial.h"
#include "bmalgo_similarity.h"
#include "bmthreadpool.h"


#include "bmdbg.h"
#include "bmtimer.h"
#include "inv_list_serial_config.h"
#include "bmundef.h" /* clear the pre-proc defines from BM */

using namespace std;

static
void show_help()
{
    std::cout
        << "BitMagic Inverted List Compression Test (c) 2019" << std::endl
        << "-u32in u32-input-file       -- raw 32-bit unsigned int file" << std::endl
        << "-bvout bvect-output-file    -- bit-vector compressed out file" << std::endl
        << "-bvin bvect-input-file      -- bit-vector compressed in file" << std::endl
        << std::endl
        << "-max-vectors N              -- process at most N input vectors in compression/verification" << std::endl
        << "-progress-vectors N         -- progress every N vectors or 30 seconds (default off)" << std::endl
        << "-bv-index-out FILE          -- write BLOB offsets during threaded compression" << std::endl
        << "-bv-index-in FILE           -- use BLOB offsets for threaded verify/decode" << std::endl
        << "-block-select-step2         -- also try independent step 2 on whole-block BWT" << std::endl
        << "-sblock-wave-adaptive       -- select adaptive wave widths using exact dry costs" << std::endl
        << "-bwt-shared                 -- share one original-polarity BWT across 72/96 selection" << std::endl
        << "-gap-xor                    -- experiment: GAP cleanup plus XOR residual" << std::endl
        << "-gap-route-stats            -- measure sparse-superblock and direct GAP record routes" << std::endl
        << "-sblock-len-stats           -- histogram BIC superblock lengths during indexed decode" << std::endl
        << "-opcode-survey FILE         -- verified per-vector opcode CSV; use -bvout /dev/null" << std::endl
        << "-profile-opcodes             -- measure emitted top-level opcode and payload bytes" << std::endl
        << "-sblock-wave                -- emit fixed 1024-bit wave superblocks when smaller" << std::endl
        << "-sblock-large-inverse-limit N -- maximum zero positions trialed (default 196607)" << std::endl
        << "-sblock-large               -- optional positive/inverse superblocks with exact fallback" << std::endl
        << "-sblock-correction-select   -- compare both opcode-68 global correction choices" << std::endl
        << "-block-endpoint-select      -- trial both existing opcode-65 endpoint layouts" << std::endl
        << "-compression-diagnostics         -- print worker timing diagnostics" << std::endl
        << "-sblock-adaptive            -- select adaptive 16/32/64/128-position delta BIC superblocks" << std::endl
        << "-sblock-gap-tolerance N     -- GAP-heavy superblock tolerance (14..100; default 100, 14 is legacy)" << std::endl
        << "-bwt-min-gain N             -- require N complete-record bytes saved by BWT (default 16; 0 disables)" << std::endl
        << "-periodic-step1-interleaved N -- BIC start/end pairs: 1 fixed, 4 gamma-vs-fixed formula choice" << std::endl
        << "-threads N                  -- parallel compression pool (1..256; omitted = sequential)" << std::endl
        << "-read-memory-mib N          -- indexed read batch budget, raw + encoded (1..8192, default 512)" << std::endl
        << "-read-mmap                  -- Unix read-only windows with read-ahead hints" << std::endl
        << "-batch-vectors N            -- MT input window (default 512; maximum 512)" << std::endl
        << "-level N                    -- compression level 0..10 (default 6)"     << std::endl
        << "-silent (-s)                -- suppress legacy status; explicit progress remains enabled"          << std::endl
        << "--silent                    -- suppress all progress, including explicit/timed reports" << std::endl
        << "-verify                     -- verify compressed version "             << std::endl
        << "Cross-bvector (xbv) research modes:" << std::endl
        << "-pair-search -bvin FILE -csv FILE -- analyze XOR pairs" << std::endl
        << "-pair-index -bvin FILE -pair-index-out FILE -- build a size-sorted shard index" << std::endl
        << "-pair-index-in FILE         -- use a precomputed index for search or bounded XOR build/verify" << std::endl
        << "-rs-bench -bvin FILE -pair-index-in FILE -rs-out PREFIX -- RS screening experiment" << std::endl
        << "-rs-set shard|stratified -rs-size N -rs-seed N -rs-shard ID (default auto)" << std::endl
        << "-rs-thresholds LIST         -- normalized L1 cutoffs (default 0.1,0.25,0.5,1)" << std::endl
        << "-rs-repeats N -rs-min-bytes N -rs-byte-coverage P (defaults 2,4096,90)" << std::endl
        << "-rs-cluster-plan -bvin FILE -pair-index-in FILE -rs-cluster-out FILE -- plan RS clusters" << std::endl
        << "-rs-cluster-feature regions3|blocks8 (default regions3); -rs-cluster-seed N (default 42)" << std::endl
        << "-rs-cluster-profiles-in PREFIX -- replay clustering from a completed RS benchmark" << std::endl
        << "-rs-cluster-vectors N (default 7000); -rs-cluster-memory N (default 8 GiB)" << std::endl
        << "-pair-cluster-in FILE       -- search clusters from a completed RS manifest" << std::endl
        << "-pair-cluster-parallel      -- evaluate one RS cluster per BM thread-pool task" << std::endl
        << "-pair-cluster-memory N      -- total estimated active-cluster memory bytes (default 24 GiB)" << std::endl
        << "-pair-limit N              -- first N vectors (0 = all; default all)" << std::endl
        << "-pair-min-gain N           -- minimum structural gain (default 256)" << std::endl
        << "-pair-min-pct P            -- minimum target gain percent (default 10)" << std::endl
        << "-pair-max-ratio R          -- maximum serialized size ratio (0 = off; otherwise >= 1)" << std::endl
        << "-pair-threads N            -- pair evaluation threads (default 1; maximum 256)" << std::endl
        << "-pair-shard-vectors N      -- maximum vectors per greedy shard (default 15000)" << std::endl
        << "-pair-shard-bytes N        -- maximum summed serialized bytes per shard (default 3 GiB)" << std::endl
        << "-pair-shard-tail-pct P     -- merge a final shard below P percent (default 5)" << std::endl
        << "-pair-inspect              -- serialize best candidates for inspection" << std::endl
        << "-ref-eval -bvin FILE -pairs-in FILE -ref-csv FILE -- compare plain XOR and reference-table serialization" << std::endl
        << "-ref-top N                 -- reference candidates per target (default 4; maximum 64)" << std::endl
        << "-xor-build -bvin FILE -pairs-in FILE -xorout FILE -- build an XOR collection" << std::endl
        << "-xor-verify -bvin FILE -xorin FILE -- verify an XOR collection against FILE" << std::endl
        << "-xor-refs FILE             -- reference metadata file (required for build/verify)" << std::endl
        << "-xor-index FILE            -- payload offset index file (required for build/verify)" << std::endl
        << "-xor-memory N              -- indexed build/verify working-memory limit (default 8 GiB)" << std::endl
        << "-decode                     -- run decode test (in-memory)"            << std::endl
        << "-timing (-t)                -- evaluate timing/duration of operations" << std::endl
      ;
}




// Arguments
//
std::string  bv_in_file;
std::string  bv_out_file;

std::string  u32_in_file;

bool         is_timing = false;
bool         is_verify = false;
bool         is_silent = false;
bool         progress_silent = false; ///< Suppress explicitly requested structured progress.
bool         is_decode = false;

bool pair_search = false, pair_inspect = false;
bool pair_cluster_parallel = false;
bm::id64_t pair_cluster_memory = 24ull * 1024 * 1024 * 1024;
std::string pair_csv;
bm::id64_t pair_limit = 0, pair_min_gain = 256;
double pair_min_pct = 10.0, pair_max_ratio = 0.0;
unsigned pair_threads = 1;
bool pair_index = false;
std::string pair_index_in, pair_index_out;
bm::id64_t pair_shard_vectors = 15000;
bm::id64_t pair_shard_bytes = 3ull * 1024 * 1024 * 1024;
double pair_shard_tail_pct = 5.0;

bool rs_bench = false;
std::string rs_out, rs_set = "shard", rs_thresholds = "0.1,0.25,0.5,1";
bm::id64_t rs_size = 7000, rs_seed = 42, rs_min_bytes = 4096, rs_repeats = 2;
bm::id64_t rs_shard = (std::numeric_limits<bm::id64_t>::max)();
double rs_byte_coverage = 90;
bool rs_cluster_plan = false;
std::string rs_cluster_out, rs_cluster_in, rs_cluster_profiles_in, rs_cluster_feature = "regions3";
bm::id64_t rs_cluster_seed = 42, rs_cluster_vectors = 7000;
bm::id64_t rs_cluster_memory = 8ull * 1024 * 1024 * 1024;

bool ref_eval = false;
std::string ref_pairs_in, ref_csv;
unsigned ref_top = 4;

bool xor_build = false, xor_verify = false;
std::string xor_out_file, xor_in_file, xor_refs_file, xor_index_file;
bm::id64_t xor_memory = 8ull * 1024 * 1024 * 1024;

bm::id64_t max_vectors = 0;      // zero means the complete input collection
bm::id64_t progress_vectors = 0; // zero suppresses periodic progress
std::string bv_index_out_file, bv_index_in_file;
bool block_selection_step2 = false;
bool gap_xor_search = false;
bool compression_diagnostics = false;
bool profile_opcodes = false;
bool block_endpoint_select = false;
bool sblock_correction_select = false;
bool sblock_large = false;
unsigned sblock_large_inverse_limit = 196607;
std::string opcode_survey_file;
bool bwt_shared = false;
bool sblock_wave1024_encode = false;
bool sblock_wave_adaptive = false;
bool gap_route_stats = false;
bool sblock_len_stats = false;
bool sblock_adaptive_search = false;
unsigned sblock_gap_tolerance = 100;
unsigned bwt_min_gain = 0;
bool bwt_min_gain_specified = false;
unsigned periodic_step1_interleaved = 0;
unsigned compression_threads = 0, compression_batch_vectors = 512;
bool compression_batch_specified = false;
bm::id64_t read_memory_bytes = 512ull * 1024 * 1024; ///< Combined raw/encoded batch budget.
bool read_mmap = false; ///< Use bounded Unix read-only mapping windows.
bool read_options_specified = false; ///< Reject read options outside indexed reads.
unsigned     c_level = bm::set_compression_default;


static
int parse_args(int argc, char *argv[])
{
    for (int i = 1; i < argc; ++i)
    {
        std::string arg = argv[i];
        if ((arg == "-h") || (arg == "--help"))
        {
            show_help();
            return 0;
        }
        
        if (arg == "-max-vectors" || arg == "-progress-vectors")
        {
            if (++i == argc) throw std::runtime_error(arg + " requires a value");
            const std::string value(argv[i]);
            if (value.empty() || value.find_first_not_of("0123456789") != std::string::npos)
                throw std::runtime_error("Invalid positive integer for " + arg);
            const bm::id64_t n = std::stoull(value);
            if (!n) throw std::runtime_error(arg + " must be positive");
            if (arg == "-max-vectors") max_vectors = n;
            else progress_vectors = n;
            continue;
        }
        if (arg == "-bv-index-out" || arg == "-bv-index-in")
        {
            if (++i == argc) throw std::runtime_error(arg + " requires a file");
            if (arg == "-bv-index-out") bv_index_out_file = argv[i];
            else bv_index_in_file = argv[i];
            continue;
        }
        if (arg == "-block-select-step2")
            { block_selection_step2 = true; continue; }
        if (arg == "-opcode-survey") {
            if (++i == argc) throw std::runtime_error("-opcode-survey requires a CSV file");
            opcode_survey_file=argv[i]; profile_opcodes=true; continue;
        }
        if (arg == "-profile-opcodes") { profile_opcodes = true; continue; }
        if (arg == "-sblock-wave-adaptive") { sblock_wave_adaptive=true;sblock_wave1024_encode=true;continue; }
        if (arg == "-sblock-wave") { sblock_wave1024_encode = true; continue; }
        if (arg == "-sblock-large-inverse-limit") {
            if(++i==argc)throw std::runtime_error(arg+" requires an integer");
            auto n=std::stoul(argv[i]);
            if(n>=bm::set_sub_total_bits/2)throw std::runtime_error("inverse limit must be below half a superblock");
            sblock_large_inverse_limit=unsigned(n);continue;
        }
        if (arg == "-sblock-large") { sblock_large = true; continue; }
        if (arg == "-sblock-correction-select") { sblock_correction_select = true; continue; }
        if (arg == "-block-endpoint-select") { block_endpoint_select = true; continue; }
        if (arg == "-bwt-shared") { bwt_shared = true; continue; }
        if (arg == "-compression-diagnostics") { compression_diagnostics = true; continue; }
        if (arg == "-gap-xor") { gap_xor_search = true; continue; }
        if (arg == "-gap-route-stats") { gap_route_stats = true; continue; }
        if (arg == "-sblock-len-stats") { sblock_len_stats = true; continue; }
        if (arg.find("-sblock-delta-") == 0)
            throw std::runtime_error("Retired delta experiment option; use -sblock-adaptive (fully adaptive)");
        if (arg == "-sblock-adaptive") { sblock_adaptive_search = true; continue; }
        if (arg == "-sblock-gap-tolerance")
        {
            if (++i == argc) throw std::runtime_error(arg + " requires a value");
            const std::string value(argv[i]);
            if (value.empty() || value.find_first_not_of("0123456789") != std::string::npos)
                throw std::runtime_error("Invalid integer for " + arg);
            const unsigned long parsed = std::stoul(value);
            if (parsed < 14 || parsed > 100)
                throw std::runtime_error("-sblock-gap-tolerance must be 14..100");
            sblock_gap_tolerance = unsigned(parsed);
            continue;
        }
        if (arg == "-periodic-step1-interleaved")
        {
            if (++i == argc) throw std::runtime_error(arg + " requires a mode");
            const std::string value(argv[i]);
            if (value != "1" && value != "4")
                throw std::runtime_error("Invalid interleaved mode (expected 1 or 4)");
            periodic_step1_interleaved = unsigned(value[0] - '0');
            continue;
        }
        if (arg == "-bwt-min-gain")
        {
            if (++i == argc) throw std::runtime_error(arg + " requires a value");
            const std::string value(argv[i]);
            if (value.empty() || value.find_first_not_of("0123456789") != std::string::npos)
                throw std::runtime_error("Invalid integer for " + arg);
            const auto n = std::stoull(value);
            if (n > 65535) throw std::runtime_error("Invalid range for " + arg);
            bwt_min_gain = unsigned(n); bwt_min_gain_specified = true; continue;
        }
        if (arg == "-read-memory-mib")
        {
            if (++i == argc)
            {
                throw std::runtime_error(arg + " requires a value");
            }
            const std::string value(argv[i]);
            if (value.empty() || value.find_first_not_of("0123456789") != std::string::npos)
            {
                throw std::runtime_error("Invalid integer for " + arg);
            }
            const auto mib = std::stoull(value);
            if (!mib || mib > 8192)
            {
                throw std::runtime_error("-read-memory-mib must be 1..8192");
            }
            read_memory_bytes = mib * 1024 * 1024;
            read_options_specified = true;
            continue;
        }
        if (arg == "-read-mmap")
        {
#if defined(__unix__) || defined(__APPLE__)
            read_mmap = true;
            read_options_specified = true;
            continue;
#else
            throw std::runtime_error("-read-mmap requires Unix");
#endif
        }
        if (arg == "-threads" || arg == "-batch-vectors")
        {
            if (++i == argc) throw std::runtime_error(arg + " requires a value");
            const std::string value(argv[i]);
            if (value.empty() || value.find_first_not_of("0123456789") != std::string::npos)
                throw std::runtime_error("Invalid integer for " + arg);
            const auto n = std::stoull(value);
            const unsigned limit = arg == "-threads" ? 256u : 512u;
            if (!n || n > limit) throw std::runtime_error("Invalid range for " + arg);
            if (arg == "-threads") compression_threads = unsigned(n);
            else { compression_batch_vectors = unsigned(n); compression_batch_specified = true; }
            continue;
        }
        if (arg == "-pair-search") { pair_search = true; continue; }
        if (arg == "-pair-cluster-parallel") { pair_cluster_parallel = true; continue; }
        if (arg == "-pair-index") { pair_index = true; continue; }
        if (arg == "-rs-bench") { rs_bench = true; continue; }
        if (arg == "-rs-cluster-plan") { rs_cluster_plan = true; continue; }
        if (arg == "-rs-cluster-out" || arg == "-rs-cluster-feature" ||
            arg == "-rs-cluster-profiles-in" ||
            arg == "-pair-cluster-in" || arg == "-pair-cluster-memory" ||
            arg == "-xor-memory" || arg == "-rs-cluster-seed" ||
            arg == "-rs-cluster-vectors" || arg == "-rs-cluster-memory")
        {
            if (++i == argc) throw std::runtime_error(arg + " requires a value");
            const std::string value(argv[i]);
            if (arg == "-rs-cluster-out") rs_cluster_out = value;
            else if (arg == "-rs-cluster-feature") rs_cluster_feature = value;
            else if (arg == "-rs-cluster-profiles-in") rs_cluster_profiles_in = value;
            else if (arg == "-pair-cluster-in") rs_cluster_in = value;
            else
            {
                if (value.empty() || value.find_first_not_of("0123456789") != std::string::npos)
                    throw std::runtime_error("Invalid integer for " + arg);
                const auto n = std::stoull(value);
                if (arg == "-rs-cluster-seed") rs_cluster_seed = n;
                else if (arg == "-rs-cluster-vectors") rs_cluster_vectors = n;
                else if (arg == "-rs-cluster-memory") rs_cluster_memory = n;
                else if (arg == "-pair-cluster-memory") pair_cluster_memory = n;
                else xor_memory = n;
            }
            continue;
        }
        if (arg == "-rs-out" || arg == "-rs-set" || arg == "-rs-thresholds" ||
            arg == "-rs-size" || arg == "-rs-seed" || arg == "-rs-min-bytes" ||
            arg == "-rs-repeats" || arg == "-rs-shard" || arg == "-rs-byte-coverage")
        {
            if (++i == argc) throw std::runtime_error(arg + " requires a value");
            const std::string value(argv[i]);
            if (arg == "-rs-out") rs_out = value;
            else if (arg == "-rs-set") rs_set = value;
            else if (arg == "-rs-thresholds") rs_thresholds = value;
            else if (arg == "-rs-byte-coverage")
            {
                std::size_t pos = 0; rs_byte_coverage = std::stod(value, &pos);
                if (pos != value.size() || !std::isfinite(rs_byte_coverage) ||
                    rs_byte_coverage <= 0 || rs_byte_coverage > 100)
                    throw std::runtime_error("RS byte coverage must be in (0,100]");
            }
            else
            {
                if (value.empty() || value.find_first_not_of("0123456789") != std::string::npos)
                    throw std::runtime_error("Invalid integer for " + arg);
                const auto n = std::stoull(value);
                if (arg == "-rs-size") rs_size = n;
                else if (arg == "-rs-seed") rs_seed = n;
                else if (arg == "-rs-min-bytes") rs_min_bytes = n;
                else if (arg == "-rs-repeats") rs_repeats = n;
                else rs_shard = n;
            }
            continue;
        }
        if (arg == "-pair-inspect") { pair_inspect = true; continue; }
        if (arg == "-ref-eval") { ref_eval = true; continue; }
        if (arg == "-xor-build") { xor_build = true; continue; }
        if (arg == "-xor-verify") { xor_verify = true; continue; }
        if (arg == "-csv" || arg == "-pair-limit" ||
            arg == "-pair-min-gain" || arg == "-pair-min-pct" ||
            arg == "-pair-max-ratio" || arg == "-pair-threads" ||
            arg == "-pair-index-in" || arg == "-pair-index-out" ||
            arg == "-pair-shard-vectors" || arg == "-pair-shard-bytes" ||
            arg == "-pair-shard-tail-pct" ||
            arg == "-pairs-in" || arg == "-ref-csv" || arg == "-ref-top" ||
            arg == "-xorout" || arg == "-xorin" ||
            arg == "-xor-refs" || arg == "-xor-index")
        {
            if (++i == argc)
                throw std::runtime_error(arg + " requires a value");
            const std::string value(argv[i]);
            if (arg == "-csv") { pair_csv = value; continue; }
            if (arg == "-pair-index-in") { pair_index_in = value; continue; }
            if (arg == "-pair-index-out") { pair_index_out = value; continue; }
            if (arg == "-pairs-in") { ref_pairs_in = value; continue; }
            if (arg == "-ref-csv") { ref_csv = value; continue; }
            if (arg == "-xorout") { xor_out_file = value; continue; }
            if (arg == "-xorin") { xor_in_file = value; continue; }
            if (arg == "-xor-refs") { xor_refs_file = value; continue; }
            if (arg == "-xor-index") { xor_index_file = value; continue; }
            std::size_t pos = 0;
            if (arg == "-pair-max-ratio")
            {
                pair_max_ratio = std::stod(value, &pos);
                if (pos != value.size() || !std::isfinite(pair_max_ratio) ||
                    (pair_max_ratio != 0 && pair_max_ratio < 1))
                    throw std::runtime_error("Pair maximum ratio must be 0 (off) or at least 1");
            }
            else if (arg == "-pair-min-pct" || arg == "-pair-shard-tail-pct")
            {
                const double pct = std::stod(value, &pos);
                if (pos != value.size() || !std::isfinite(pct) || pct < 0 || pct > 100)
                    throw std::runtime_error("Invalid percentage for " + arg);
                if (arg == "-pair-min-pct") pair_min_pct = pct;
                else pair_shard_tail_pct = pct;
            }
            else
            {
                if (value.empty() || value.find_first_not_of("0123456789") != std::string::npos)
                    throw std::runtime_error("Invalid nonnegative integer for " + arg);
                auto n = std::stoull(value, &pos);
                if (arg == "-pair-limit") pair_limit = n;
                else if (arg == "-pair-shard-vectors")
                {
                    if (!n) throw std::runtime_error("Pair shard vector limit must be positive");
                    pair_shard_vectors = n;
                }
                else if (arg == "-pair-shard-bytes")
                {
                    if (!n) throw std::runtime_error("Pair shard byte limit must be positive");
                    pair_shard_bytes = n;
                }
                else if (arg == "-pair-threads")
                {
                    if (!n || n > 256)
                        throw std::runtime_error("Pair thread count must be from 1 to 256");
                    pair_threads = unsigned(n);
                }
                else if (arg == "-ref-top")
                {
                    if (!n || n > 64)
                        throw std::runtime_error("Reference candidate count must be from 1 to 64");
                    ref_top = unsigned(n);
                }
                else pair_min_gain = n;
            }
            continue;
        }
        if (arg == "-v" || arg == "-verify")
        {
            is_verify = true;
            continue;
        }
        if (arg == "-decode")
        {
            is_decode = true;
            continue;
        }
        if (arg == "-l" || arg == "-level")
        {
            if (i + 1 < argc)
            {
                const char* lvl = argv[++i];
                char *end;
                errno = 0;
                const unsigned long parsed_level = std::strtoul(lvl, &end, 10);
                if (errno == ERANGE || end == lvl || *end || parsed_level > bm::set_compression_max)
                {
                    std::cerr << "Error parsing -level: range error for "
                              << lvl<< std::endl;
                    return 1;
                }
                c_level = unsigned(parsed_level);
            }
            else
            {
                std::cerr << "Error: -level requires compression level number" << std::endl;
                return 1;
            }
            continue;
        }

        if (arg == "-bvout" || arg == "--bvout")
        {
            if (i + 1 < argc)
            {
                bv_out_file = argv[++i];
            }
            else
            {
                std::cerr << "Error: -bvout requires file name" << std::endl;
                return 1;
            }
            continue;
        }
        if (arg == "-bvin" || arg == "--bvin")
        {
            if (i + 1 < argc)
            {
                bv_in_file = argv[++i];
            }
            else
            {
                std::cerr << "Error: -bvin requires file name" << std::endl;
                return 1;
            }
            continue;
        }

        if (arg == "-svin" || arg == "--svin" ||
            arg == "-svout" || arg == "--svout")
        {
            throw std::runtime_error("Retired SV option; use experiments/sparse_vector/svexp.cpp");
        }

        if (arg == "-u32in" || arg == "--u32in")
        {
            if (i + 1 < argc)
            {
                u32_in_file = argv[++i];
            }
            else
            {
                std::cerr << "Error: -u32in requires file name" << std::endl;
                return 1;
            }
            continue;
        }
        
        if (arg == "--silent")
        {
            is_silent = true;
            progress_silent = true;
            continue;
        }
        if (arg == "-silent" || arg == "-s" || arg == "--s")
        {
            is_silent = true;
            continue;
        }
        if (arg == "-timing" || arg == "--timing" || arg == "-t" || arg == "--t")
        {
            is_timing = true;
            continue;
        }
        throw std::runtime_error("Unknown option: " + arg);
    } // for i
    return 0;
}


// Globals
//

typedef bm::sparse_vector<unsigned, bm::bvector<> > sparse_vector_u32;
typedef bm::rsc_sparse_vector<unsigned, sparse_vector_u32> rsc_sparse_vector_u32;

bm::chrono_taker<>::duration_map_type  timing_map;


/// Read 32-bit vector size-prefix format (length:0, 1, 2, 3, ....)
///
template<class VT>
int io_read_u32_coll(std::ifstream& fin, VT& vec)
{
    typedef typename VT::value_type  value_type;
    vec.resize(0);
    if (!fin.good())
        return -1;
    value_type len;
    fin.read((char*) &len, std::streamsize(sizeof(len)));
    if (!fin.good())
        return -1;
    if (!len)
        return -2; // 0-len detected (broken file)
    
    vec.resize(len);
    fin.read((char*) &vec[0], std::streamsize(len*sizeof(value_type)));
    if (!fin.good())
        return -1;

    return 0;
}

/// Check if input vector is monotonously sorted (true inverted list)
/// along the way in computes a minimal delta between values
///
template<typename VT>
int validate_inp_vec(const VT& vec,
                     typename VT::value_type& min_delta,
                     typename VT::value_type& min_delta_cnt
                     )
{
    typename VT::value_type md_cnt = min_delta_cnt = 0;
    auto sz = vec.size();
    if (!sz)
        return -1;
    auto i_prev = vec[0];
    typename VT::value_type md = ~(typename VT::value_type(0)); // max possible uint
    for (typename VT::size_type i = 1; i < sz; ++i)
    {
        auto val = vec[i];
        if (val <= i_prev)
            return -2;
        typename VT::value_type d1 = val - i_prev;
        if (d1 < md)
        {
            md_cnt = 0; md = d1;
        }
        md_cnt += (d1 == md);
        i_prev = val;
    }
    min_delta = md;
    min_delta_cnt = md_cnt;
    return 0;
}

/** Compare a decoded vector logically with an independently imported ID list.
    @tparam ValuesVect Contiguous container of sorted, unique unsigned IDs.
    @tparam BV Bit-vector type used by the reader.
    @param vec Read-only expected IDs; storage remains owned by the caller.
    @param bv Read-only decoded vector.
    @return Zero for exact logical equality, minus one for any difference.
    @pre Input IDs are sorted and unique, as required by the collection format.
*/
template<typename ValuesVect, typename BV>
int compare_vect(const ValuesVect& vec, const BV& bv)
{
    if (vec.size() != bv.count())
    {
        return -1;
    }
    BV expected;
    if (!vec.empty())
    {
        expected.set(vec.data(), typename BV::size_type(vec.size()), bm::BM_SORTED);
    }
    return bv.equal(expected) ? 0 : -1;
}

/// Debug utility to detect super sparse bit-vectors
/// which probably get bad compression rate
///
template<typename BV>
bool is_super_sparse(const BV& bv)
{
    typename BV::statistics st;
    bv.calc_stat(&st);
    auto bc = bv.count();
    auto blocks_count = st.gap_blocks + st.bit_blocks;
    if (bc <= blocks_count)
        return true;
    auto bc_parity = blocks_count * 6;
    return (bc <= bc_parity);
}


/** Convert a sorted nonempty ID list into an optimized bvector and append it.
    @tparam ValuesVect Contiguous container of sorted, unique unsigned IDs.
    @param bv_output Destination collection stream, owned by the caller.
    @param vec Nonempty input list; remains unchanged.
    @param bvs Exclusively owned serializer, reused for encoding.
    @throws std::runtime_error If serialization fails.
*/
template<typename ValuesVect>
void write_as_bvector(bm::streams_encoder& bv_output,
                     const ValuesVect& vec,
                     bm::serializer<bm::bvector<> >& bvs)
{
    BM_DECLARE_TEMP_BLOCK(tb)
    bm::bvector<> bv;
    bv.set(&vec[0], bm::bvector<>::size_type(vec.size()), bm::BM_SORTED);

    bv.optimize(tb);
    if (!bvs.serialize(bv, bv_output))
    {
        throw std::runtime_error("Error write to bvect out file");
    }
}

/** Capture CLI serializer settings for synchronous configuration.
    @return Independent settings value; no references to CLI globals are retained.
    @note CLI configuration is read-only during compression worker execution.
*/
static inv_list_serial_options inv_serializer_options()
{
    inv_list_serial_options options;
    options.level = c_level;
    options.gap_tolerance = sblock_gap_tolerance;
    options.inverse_limit = sblock_large_inverse_limit;
    options.bwt_min_gain = bwt_min_gain;
    options.bwt_gain_specified = bwt_min_gain_specified;
    options.block_endpoints = block_endpoint_select;
    options.sblock_correction = sblock_correction_select;
    options.sblock_large = sblock_large;
    options.wave_fixed = sblock_wave1024_encode;
    options.wave_adaptive = sblock_wave_adaptive;
    options.bwt_shared = bwt_shared;
    return options;
}

/// read the input collection sequence, write using various compression schemes
///
static
void compress_inv_dump_file(const std::string& fname,
                            const std::string& bv_out_fname)
{
    bm::id64_t total_ints = 0;
    bm::id64_t min_delta_ints = 0;

    cout << "Reading input collection: " << fname << endl;
    if (!bv_out_fname.empty())
        cout << "Writing to BV collection: " << bv_out_fname << endl;
    else
        cout << "NO BV collection specified" << endl;
    cout << "Compression level: " << c_level << endl;

    bm::chrono_taker tt1(std::cout, "1. Convert collection", 1, &timing_map);

    vector<unsigned> vec;
    std::ifstream fin(fname.c_str(), std::ios::in | std::ios::binary);
    if (!fin.good())
    {
        throw std::runtime_error("Cannot open input file");
    }
    
    std::ofstream bv_file;
    std::unique_ptr<bm::streams_encoder> bv_output;
    if (!bv_out_fname.empty())
    {
        bv_file.open(bv_out_fname, std::ios::out | std::ios::binary);
        if (!bv_file.good())
            throw std::runtime_error("Cannot open bvect out file");
        bv_output.reset(new bm::streams_encoder(bv_file));
        if (!bv_output->is_good())
            throw std::runtime_error("Cannot initialize bvect out stream");
    }
    fin.seekg(0, std::ios::end);
    std::streamsize fsize = fin.tellg();

    fin.seekg(0, std::ios::beg);
    
    // initialize serializer
    //
    bm::serializer<bm::bvector<> > bvs;
    configure_inv_serializer(bvs, inv_serializer_options(), [](auto& bvs)
    {
        if (block_selection_step2)
        {
            bvs.set_block_selection_step2_options();
            if (bwt_min_gain_specified)
            {
                bvs.set_block_selection_bwt_min_gain(bwt_min_gain);
            }
            if (periodic_step1_interleaved)
            {
                bvs.set_periodic_step1_interleaved(periodic_step1_interleaved);
            }
        }
    });

    // main loop to read sample vectors
    //
    const auto progress_start = std::chrono::steady_clock::now();
    bm::id64_t processed_vectors = 0;
    bm::id64_t i;
    for (i = 0; true; ++i)
    {
        if (max_vectors && i >= max_vectors) break;
        if (block_selection_step2) bvs.set_periodic_vector_index(i);
        if (block_selection_step2) bvs.set_bwt_vector_index(i);
        int ret = io_read_u32_coll(fin, vec);
        if (ret != 0)
            throw std::runtime_error("Error reading input file");
        unsigned min_delta, min_delta_cnt;

        {
            ret = validate_inp_vec(vec, min_delta, min_delta_cnt); // check if we have sorted unique set
            if (ret != 0)
                throw std::runtime_error("Input vector validation failed");
            if (!min_delta || !min_delta_cnt)
                throw std::runtime_error("Input vector validation delta error");
            
            min_delta_ints += min_delta_cnt;
        }
        total_ints += vec.size(); // remember the total size of the collection
        
        // Serialize consecutive self-delimiting bit-vector BLOBs.
        //
        size_t bv_blob_size = 0;
        const auto vector_start = std::chrono::steady_clock::now();
        if (!bv_out_fname.empty())
        {
            const size_t bv_start = bv_output->size();
            write_as_bvector(*bv_output, vec, bvs);
            bv_blob_size = bv_output->size() - bv_start;
        }
        
        std::streamsize fpos_curr = fin.tellg();
        const bm::id64_t completed = ++processed_vectors;
        if (progress_vectors && !progress_silent &&
            (completed % progress_vectors == 0 ||
             (max_vectors && completed == max_vectors) || fpos_curr == fsize))
        {
            const auto now = std::chrono::steady_clock::now();
            const double elapsed = std::chrono::duration<double>(now - progress_start).count();
            const double last = std::chrono::duration<double>(now - vector_start).count();
            std::ostringstream progress;
            progress << "Progress compress vectors=" << completed
                     << " input_bytes=" << fpos_curr
                     << " output_bytes=" << (bv_output ? bv_output->size() : 0)
                     << " elapsed_s=" << std::fixed << std::setprecision(2) << elapsed
                     << " vectors_per_s=" << (elapsed > 0 ? double(completed) / elapsed : 0)
                     << " last_vector_s=" << last
                     << " last_vector_bytes=" << bv_blob_size;

            std::cerr << progress.str() << std::endl;
        }
        if (fpos_curr == fsize)
            break;
        
        if (!is_silent)
        {
            cout << "\r" << i << "   " << fpos_curr << " / " << fsize
                 << " ( size=" << vec.size() << " ) "
                 << flush;
        }

    } // for i

    if (bv_output && !bv_output->finish())
        throw std::runtime_error("Error finishing bvect out file");
    
    // print statistics about test set
    
    cout << endl;
    cout << "Total vectors=" << processed_vectors << endl;
    cout << "Input consumed bytes=" << fin.tellg() << endl;
    cout << "Input file bytes=" << fsize << endl;
    cout << "Total ints=" << total_ints << endl;
    cout << "  min-deltas = " << min_delta_ints << endl;
    {
        double min_delta_ratio = double(min_delta_ints) / double(total_ints);
        cout << "  min delta ratio = " << std::setprecision(3) << min_delta_ratio << endl;
    }
    if (!bv_out_fname.empty())
    {
        bm::id64_t bv_size = (bm::id64_t)bv_output->size();
        cout << "BV size = " << bv_size << endl;
        double bv_bits_per_int = double(bv_size * 8ull) / double(total_ints);
        cout << "BV Bits per/int = " << std::setprecision(3) << bv_bits_per_int << endl;
        
        bv_file.close();
        if (!bv_file.good())
            throw std::runtime_error("Error closing bvect out file");
    }
    

}

// ------------------------------------------------------------------------

/// read the input collection sequence and dump file, verify correctness
///
static
void verify_inv_dump_file(const std::string& fname,
                          const std::string& bv_in_fname)
{
    bm::id64_t total_ints = 0;


    bm::chrono_taker tt1(std::cout, "2. Verify collection", 1, &timing_map);


    cout << "Reading input collection: " << fname << endl;
    if (!bv_in_fname.empty())
        cout << "Reading BV collection: " << bv_in_fname << endl;
    else
        cout << "NO BV collection specified" << endl;


    vector<unsigned> vec;
    std::ifstream fin(fname.c_str(), std::ios::in | std::ios::binary);
    if (!fin.good())
    {
        throw std::runtime_error("Cannot open input file");
    }
    
    if (bv_in_fname.empty())
        throw std::runtime_error("No bvect dump file specified");
    std::ifstream bv_file(bv_in_fname.c_str(), std::ios::in | std::ios::binary);
    if (!bv_file.good())
        throw std::runtime_error("Cannot open bvect dump file");
    fin.seekg(0, std::ios::end);
    std::streamsize fsize = fin.tellg();
    fin.seekg(0, std::ios::beg);
    bv_file.seekg(0, std::ios::end);
    std::streamsize bv_fsize = bv_file.tellg();
    bv_file.seekg(0, std::ios::beg);
    if (!fin.good() || !bv_file.good() || fsize < 0 || bv_fsize < 0)
        throw std::runtime_error("Error positioning input files");

    bm::streams_decoder bv_input(bv_file);
    bm::streams_deserializer<bm::bvector<> > bv_reader;

    // main loop to read sample vectors
    //
    const auto progress_start = std::chrono::steady_clock::now();
    bm::id64_t processed_vectors = 0;
    bm::id64_t i;
    for (i = 0; true; ++i)
    {
        if (max_vectors && i >= max_vectors) break;
        int ret = io_read_u32_coll(fin, vec);
        if (ret != 0)
            throw std::runtime_error("Error reading input file");

        total_ints += vec.size(); // remember the total size of the collection
        
        bm::bvector<> bv;
        if (!bv_reader.deserialize(bv, bv_input))
            throw std::runtime_error("Error reading bvect dump file");
        if (compare_vect(vec, bv) != 0)
            throw std::runtime_error("Vector comparison failed");
        
        std::streamsize fpos_curr = fin.tellg();
        const bm::id64_t completed = ++processed_vectors;
        if (progress_vectors && !progress_silent &&
            (completed % progress_vectors == 0 ||
             (max_vectors && completed == max_vectors) || fpos_curr == fsize))
        {
            const double elapsed = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - progress_start).count();
            std::cerr << "Progress verify vectors=" << completed
                      << " input_bytes=" << fpos_curr
                      << " encoded_bytes=" << bv_file.tellg()
                      << " elapsed_s=" << std::fixed << std::setprecision(2) << elapsed
                      << " vectors_per_s=" << (elapsed > 0 ? double(completed) / elapsed : 0)
                      << std::endl;
        }
        if (fpos_curr == fsize)
            break;
        
        if (!is_silent)
        {
            cout << "\r" << fpos_curr << "/" << fsize
                 << " ( size=" << vec.size() << " )              "
                 << flush;
        }
    } // for i

    if (bv_file.tellg() != bv_fsize)
        throw std::runtime_error("Extra data in bvect dump file");
    
    cout << endl;
    cout << "Verification complete." << endl;
    cout << "Total vectors=" << processed_vectors << endl;
    cout << "Input consumed bytes=" << fin.tellg() << endl;
    cout << "Input file bytes=" << fsize << endl;
    cout << "Total ints=" << total_ints << endl;
}

/// read and decode the compressed dump file
///
static
void decode_test_dump_file(const std::string& bv_in_fname)
{
    bm::chrono_taker tt1(std::cout, "3. Decode collection", 1, &timing_map);

    if (bv_in_fname.empty())
        throw std::runtime_error("Cannot open bvect dump file");
    std::ifstream bv_file(bv_in_fname.c_str(), std::ios::in | std::ios::binary);
    if (!bv_file.good())
        throw std::runtime_error("Cannot open bvect dump file");
    bv_file.seekg(0, std::ios::end);
    std::streamsize fsize = bv_file.tellg();
    bv_file.seekg(0, std::ios::beg);
    if (!bv_file.good() || fsize < 0)
        throw std::runtime_error("Error positioning bvect dump file");
    if (!fsize)
        throw std::runtime_error("Empty bvect dump file");

    bm::streams_decoder bv_input(bv_file);
    bm::streams_deserializer<bm::bvector<> > bv_reader;

    // main loop to read sample vectors
    //
    const auto progress_start = std::chrono::steady_clock::now();
    bm::id64_t i = 0;
    for (;; ++i)
    {
        std::streamsize bv_start = bv_file.tellg();
        if (bv_start == fsize)
            break;
        if (bv_start < 0 || bv_start > fsize)
            throw std::runtime_error("Invalid position in bvect dump file");

        bm::bvector<> bv;
        if (!bv_reader.deserialize(bv, bv_input))
            throw std::runtime_error("Error reading bvect dump file");

        std::streamsize fpos_curr = bv_file.tellg();
        if (fpos_curr <= bv_start || fpos_curr > fsize)
            throw std::runtime_error("Invalid position in bvect dump file");

        if (progress_vectors && !progress_silent &&
            ((i + 1) % progress_vectors == 0 || fpos_curr == fsize))
        {
            const double elapsed = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - progress_start).count();
            std::cerr << "Progress decode vectors=" << (i + 1)
                      << " encoded_bytes=" << fpos_curr
                      << " elapsed_s=" << std::fixed << std::setprecision(2) << elapsed
                      << " vectors_per_s=" << (elapsed > 0 ? double(i + 1) / elapsed : 0)
                      << std::endl;
        }
        if (!is_silent)
        {
            cout << "\r" << fpos_curr << "/" << fsize
                 << flush;
        }
    } // for i

    cout << endl;
    cout << "Decode complete." << endl;
    cout << "Total vectors=" << i << endl;
}



#include "inv_list_xbv_pair_search.h"
#include "inv_list_xbv_csv.h"
#include "inv_list_xbv_pair_shard.h"
#include "inv_list_xbv_pair_cluster_parallel.h"
#include "inv_list_xbv_ref_eval.h"
#include "inv_list_xbv_xor_store.h"
#include "inv_list_xbv_rs_bench.h"
#include "inv_list_xbv_rs_cluster.h"
#include "inv_list_parallel.h"

int main(int argc, char *argv[])
{
    if (argc < 2)
    {
        show_help();
        return 1;
    }
    
    try
    {
        for (int i = 1; i < argc; ++i)
        {
            if (std::string(argv[i]) == "-h" || std::string(argv[i]) == "--help")
            {
                show_help();
                return 0;
            }
        } // for i
        auto ret = parse_args(argc, argv);
        if (ret != 0)
            return ret;
        
        const unsigned analysis_modes = unsigned(rs_bench) + unsigned(rs_cluster_plan) + unsigned(pair_search) + unsigned(pair_index) + unsigned(ref_eval) +
                                        unsigned(xor_build) + unsigned(xor_verify);
        if ((max_vectors || progress_vectors) && analysis_modes)
            throw std::runtime_error("Vector limit/progress applies only to compression, verification, or decode");
        if (max_vectors && is_decode && u32_in_file.empty())
            throw std::runtime_error("-max-vectors applies to input collection compression/verification; decode reads its complete BLOB");
        if (block_selection_step2 && (analysis_modes || pair_cluster_parallel ||
            is_verify || is_decode || bv_out_file.empty()))
            throw std::runtime_error("-block-select-step2 requires independent BV compression");
        if (gap_xor_search && (!block_selection_step2 || !compression_threads ||
                               c_level != 6 || u32_in_file.empty() ||
                               bv_out_file.empty() || is_verify || is_decode))
            throw std::runtime_error("-gap-xor requires threaded level-6 -block-select-step2 BV compression");
        if (gap_route_stats && (!block_selection_step2 || !compression_threads ||
                               c_level != 6 || u32_in_file.empty() ||
                               bv_out_file.empty() || is_verify || is_decode || gap_xor_search))
            throw std::runtime_error("-gap-route-stats requires threaded level-6 -block-select-step2 BV compression without -gap-xor");
        if (!opcode_survey_file.empty() &&
            (bv_out_file != "/dev/null" || !bv_index_out_file.empty() ||
             sblock_correction_select || sblock_large))
            throw std::runtime_error("-opcode-survey requires -bvout /dev/null, no index/SV output, and retained baseline codecs");
        if (profile_opcodes && (!compression_threads || u32_in_file.empty() ||
                                bv_out_file.empty() || is_verify || is_decode || analysis_modes))
            throw std::runtime_error("-profile-opcodes requires threaded BV compression");
        if ((block_endpoint_select || sblock_correction_select || sblock_large) && (!compression_threads || !block_selection_step2 ||
                                      c_level != 6 || bv_out_file.empty() ||
                                      u32_in_file.empty() || is_verify || is_decode))
            throw std::runtime_error("Endpoint/correction selection requires threaded level-6 BV compression");
        if (bwt_shared && (!compression_threads || !block_selection_step2 ||
                            c_level != 6 || bv_out_file.empty() ||
                            u32_in_file.empty() || is_verify || is_decode))
            throw std::runtime_error("Shared BWT requires threaded level-6 -block-select-step2 BV compression");
        if (sblock_wave1024_encode && (!compression_threads || !block_selection_step2 ||
                                   c_level != 6 || bv_out_file.empty() ||
                                   u32_in_file.empty() || is_verify || is_decode))
            throw std::runtime_error("-sblock-wave requires threaded level-6 BV compression");
        if (sblock_len_stats && (!is_decode || !compression_threads || bv_in_file.empty() ||
                                 bv_index_in_file.empty() || is_verify || !u32_in_file.empty()))
            throw std::runtime_error("-sblock-len-stats requires indexed threaded BV decode");
        if (sblock_adaptive_search && (!block_selection_step2 || !compression_threads ||
                                      c_level != 6 || u32_in_file.empty() ||
                                      bv_out_file.empty() || is_verify || is_decode ||
                                      gap_route_stats))
            throw std::runtime_error("-sblock-adaptive requires threaded level-6 -block-select-step2 BV compression");
        if (bwt_min_gain_specified && !block_selection_step2 && c_level < 9)
            throw std::runtime_error("-bwt-min-gain requires a block-selection mode");
        if (periodic_step1_interleaved && !block_selection_step2)
            throw std::runtime_error("-periodic-step1-interleaved requires -block-select-step2");
        if (compression_threads &&
            (analysis_modes ||
             (is_verify || is_decode ? bv_index_in_file.empty() || bv_in_file.empty() :
              (u32_in_file.empty() || bv_out_file.empty()))))
            throw std::runtime_error("-threads requires BV compression or indexed BV verify/decode");
        if (read_options_specified && (!compression_threads || (!is_verify && !is_decode)))
        {
            throw std::runtime_error("Read options require indexed threaded verify/decode");
        }
        if (compression_batch_specified && !compression_threads)
            throw std::runtime_error("-batch-vectors requires -threads");
        if (!bv_index_out_file.empty() &&
            (!compression_threads || is_verify || is_decode || bv_out_file.empty()))
            throw std::runtime_error("-bv-index-out requires threaded BV compression");
        if (!bv_index_in_file.empty() &&
            (!compression_threads || (!is_verify && !is_decode) || bv_in_file.empty()))
            throw std::runtime_error("-bv-index-in requires threaded BV verify or decode");

        if (analysis_modes > 1)
            throw std::runtime_error("RS benchmarking, cluster planning, pair indexing/search, reference evaluation, and XOR build/verification are separate modes");
        if (rs_bench)
        {
            if (is_verify || is_decode || !u32_in_file.empty() || !bv_out_file.empty())
                throw std::runtime_error("RS benchmark must run separately");
            run_rs_benchmark();
            return 0;
        }
        if (rs_cluster_plan)
        {
            if (is_verify || is_decode || !u32_in_file.empty() || !bv_out_file.empty())
                throw std::runtime_error("RS cluster planning must run separately");
            build_rs_cluster_plan();
            return 0;
        }
        if (pair_index)
        {
            if (is_verify || is_decode || !u32_in_file.empty() || !bv_out_file.empty())
                throw std::runtime_error("Pair indexing must run separately from compression/verification/decode");
            build_pair_shard_index();
            return 0;
        }
        if (pair_search)
        {
            if (is_verify || is_decode || !u32_in_file.empty() || !bv_out_file.empty())
                throw std::runtime_error("Pair search must run separately from compression/verification/decode");
            if (pair_cluster_parallel) analyze_pairs_cluster_parallel();
            else if (pair_index_in.empty())
            {
                const pair_search_options options{bv_in_file, pair_csv, pair_limit,
                    pair_min_gain, pair_min_pct, pair_max_ratio, pair_threads,
                    c_level, pair_inspect, is_silent};
                analyze_pairs(options);
            }
            else analyze_pairs_sharded();
            return 0;
        }
        if (ref_eval)
        {
            if (is_verify || is_decode || !u32_in_file.empty() || !bv_out_file.empty())
                throw std::runtime_error("Reference evaluation must run separately from compression/verification/decode");
            analyze_references();
            return 0;
        }
        if (xor_build)
        {
            if (is_verify || is_decode || !u32_in_file.empty() || !bv_out_file.empty())
                throw std::runtime_error("XOR build must run separately from compression/verification/decode");
            build_xor_store();
            return 0;
        }
        if (xor_verify)
        {
            if (is_verify || is_decode || !u32_in_file.empty() || !bv_out_file.empty())
                throw std::runtime_error("XOR verification must run separately from compression/verification/decode");
            verify_xor_store();
            return 0;
        }
        if (!u32_in_file.empty())
        {
            if (!is_verify)
            {
                cout << "Compression" << endl;
                if (compression_threads)
                {
                    if (sblock_adaptive_search && gap_xor_search)
                        compress_inv_dump_file_parallel<sblock_adaptive_gap_xor_policy>(u32_in_file, bv_out_file);
                    else if (sblock_adaptive_search)
                        compress_inv_dump_file_parallel<sblock_adaptive_compression_policy>(u32_in_file, bv_out_file);
                    else if (gap_route_stats)
                        compress_inv_dump_file_parallel<gap_route_compression_policy>(u32_in_file, bv_out_file);
                    else if (gap_xor_search)
                        compress_inv_dump_file_parallel<gap_xor_compression_policy>(u32_in_file, bv_out_file);
                    else if (block_selection_step2)
                        compress_inv_dump_file_parallel<block_selection_step2_policy>(u32_in_file, bv_out_file);
                    else
                        compress_inv_dump_file_parallel<plain_compression_policy>(u32_in_file, bv_out_file);
                }
                else compress_inv_dump_file(u32_in_file, bv_out_file);
            }
            else
            {
                if (is_verify)
                {
                    cout << "Verification." << endl;
                    if (compression_threads)
                        read_inv_dump_file_parallel(u32_in_file, bv_in_file, true);
                    else verify_inv_dump_file(u32_in_file, bv_in_file);
                }
            }
        }
        if (is_decode)
        {
            cout << "Decode test." << endl;
            if (compression_threads)
                read_inv_dump_file_parallel(std::string(), bv_in_file, false);
            else decode_test_dump_file(bv_in_file);
        }

                
        if (is_timing)  // print all collected timings
        {
            std::cout << std::endl << "Timings (ms):" << std::endl;
            bm::chrono_taker<>::print_duration_map(std::cout, timing_map);
        }
    }
    catch (std::exception& ex)
    {
        std::cerr << "Error:" << ex.what() << std::endl;
        return 1;
    }

    return 0;
}


#ifdef _MSC_VER
#pragma warning( pop )
#endif
