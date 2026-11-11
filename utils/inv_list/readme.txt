inv_list: inverted-list compression and validation utility

Build from utils/inv_list:
    make
On AArch64, select native ARM NEON with:
    make neon
Build and execution are separate. List supported CLI options with ./bminv --help.

INPUT AND OUTPUT
Raw input is a sequence of lists in the Gov2 flat format: a little-endian
32-bit length followed by that many sorted, unique unsigned 32-bit IDs.
Input datasets remain under /Volumes/USB_1TB/int_set_bench/Gov2Flat by default
in the current harness; local int_set_bench inputs are also preserved.

BV output is a sequence of self-delimiting BitMagic bvector BLOBs without
per-list length prefixes. An optional .idx sidecar stores BLOB boundaries for
threaded verification and decoding. Older length-prefixed collections need
to be regenerated. Retired -svin/-svout options fail explicitly; use the
separate experiments/sparse_vector utility.

CURRENT GOV2 WORKFLOW
See experiments/current/README.md for the full sorted/unsorted harness,
USB output requirements, profiling, and focused validation commands.
Compression presets 0..10 are selected with -level; the default is 6.
See ../../docs/compression.md for the level table and terminology.

Examples (run from utils/inv_list):
    ./bminv -s -level 10 -threads 10 -u32in INPUT -bvout OUTPUT \
        -bv-index-out OUTPUT.idx -progress-vectors 10000
    ./bminv -s -threads 10 -u32in INPUT -bvin OUTPUT \
        -bv-index-in OUTPUT.idx -verify
    ./bminv -s -threads 10 -bvin OUTPUT -bv-index-in OUTPUT.idx -decode

Use -max-vectors N to bound input compression/verification. Use -threads 1
for controlled single-worker timings. Omit -threads for the serial path.
Threaded progress updates after N vectors or 30 seconds when enabled; the
serial path currently uses vector-count boundaries.

RESEARCH AND DIAGNOSTICS
Cross-bvector (xbv) pair search, sharded search, reference evaluation, XOR stores,
and RS clustering remain available in bminv. Their run_gov2_xbv_*.sh wrappers
provide workflow-specific configuration. See experiments/RESEARCH_HISTORY.txt
for historical rationale and experiment records; use --help and current
scripts as the authority for supported arguments.

Generic opcode profiling and playlist survey/replay tools are documented in
experiments/README.md. Inputs stay outside cleanup operations; generated
collections and logs belong on the USB experiment drive.

WORKSPACE AND DEVELOPMENT DIAGNOSTICS
Serializers/deserializers retain scratch capacity; destroy them to release it.
BM_SERIALIZATION_DIAGNOSTICS is enabled by inv_list.cpp for codec histories.
Use a consistent macro value across translation units. Optional histories may
grow during an operation; use bounded workloads for detailed diagnostics.

RELEASE CLEANUP AND PORTABLE COLLECTION RUNNER
Use run_collection.sh for Gov2, ClueWeb, scalar or SIMD variants:
    bash run_collection.sh --bin ./bminv --input /data/gov2.sorted \
        --out /results/gov2-sorted --levels 6,7,8,9,10 --threads 10
The output directory must not exist; its parent must exist. Each level gets
its own stream, index, compression/verification/decode logs and exact commands.
Omit --threads for serial execution. Add --max-vectors 1000 for a smoke run.
Repeat with another --input/--out for unsorted or ClueWeb collections, and
another --bin/--out for SIMD comparisons. All streams are retained.
This replaces run_gov2_{sorted,unsorted}[_arm].sh and run_clueweb_*.sh.
The profiled full-Gov2 workflow in experiments/current remains available.

xbv means cross-bvector: headers use inv_list_xbv_*, workflow runners use
run_gov2_xbv_*, and Python analysis/plotting tools include xbv in their names.
The established -pair-*, -rs-*, -ref-* and -xor-* CLI interfaces and output
formats are unchanged. For the former run_complexity.sh example, use:
    ./bminv -pair-search -bvin INPUT.bv -csv OUTPUT.csv -pair-limit 10000 \
        -pair-min-pct 5 -pair-max-ratio 2 -pair-threads 4 -pair-inspect

Standalone wave/full BWT, periodic search, compose/multistep selection and
start-reduction/first-start ablations have been retired from the CLI. Their
switches now fail as unknown options; -diag/-d was ineffective and is removed.
Use presets for normal compression. Step-2, interleaved step-1 and targeted
codec controls remain for retained comparisons; library formats are unchanged.

SILENT PROGRESS
Use --silent to suppress all progress, including explicitly requested vector-count
and timed updates. Final statistics and errors remain visible. Legacy -s/-silent
retains structured progress requested by -progress-vectors for existing harnesses.
run_inv_test_pl.sh --silent hides console progress while retaining full phase logs.

INDEXED READ MEMORY AND UNIX MAPPING
Threaded verify/decode now reads adjacent indexed BLOBs in one sequential batch,
then workers decode bounded RAM slices. Full raw-ID comparison is retained.
-read-memory-mib N budgets raw IDs plus encoded bytes per batch (default 512 MiB,
maximum 8192 MiB). -batch-vectors also bounds records per batch (default 512).
The budget excludes the index, worker decoder/bitvector/comparison workspaces,
optional diagnostic records and filesystem cache. Leave RAM for those users.
A single vector exceeding the budget fails before allocating its raw payload;
increase the budget explicitly. All workers finish before storage is reused.
On Unix, -read-mmap selects read-only batch windows with advisory sequential
and WILLNEED hints when available. Windows uses the portable buffered reader.
Mapping covers the compressed input; raw IDs still use sequential buffered I/O.
Mapping assumes input files remain unchanged. It does not reserve physical RAM.
