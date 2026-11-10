> November 3 decoder cleanup: opcode 92 now accepts only the retained adaptive
> step descriptor and absolute (mode 0) or signed previous-relative (mode 2)
> corrections. The unreleased fixed-step header and global-relative mode 1
> decoder branches are removed. Current emitted bytes and all released opcodes
> are unchanged. Use [current validation](../current/README.md).

> Historical experiment notes. The exploratory measurement code and old
> runners described below were retired on 2026-11-01. Use
> [current highest-compression verification](../current/README.md).
> Retained codec correctness tests remain where applicable.

# Retired Elias delta correction experiment (2026-10-30)

Elias delta was not a productive change for superblock correction encoding.
The full Gov2 sorted/unsorted control and candidate streams each verified
and decoded 1,107,205 vectors and 5,979,715,441 IDs. The experiment compared
absolute corrections and signed previous-relative differences using both
gamma and Elias delta, with complete-record fallback.

| Collection | Gamma control bytes | Candidate bytes | Saved bytes | Elias trials | Selected (absolute / previous) |
| --- | ---: | ---: | ---: | ---: | ---: |
| sorted | 2,678,408,172 | 2,678,406,280 | 1,892 (0.000071%) | 1,386,166 | 473 (410 / 63) |
| unsorted | 5,441,839,210 | 5,441,839,206 | 4 | 4,757,373 | 1 (1 / 0) |

These gains do not justify another codec and search path. The user requested
retirement: experimental opcode 90, encoder/decoder, CLI, diagnostics, tests
and runner were removed. Retain gamma absolute/previous-relative selection,
grouped flags, and adaptive steps under opcode 92. Released opcodes are unchanged.
Single-pass timings did not establish a speed advantage. Reports/logs remain:
`/Volumes/USB_1TB/BitMagic-experiments/2026-10-30/runs/inv_list/gov2-elias-delta-full-20261030-091640/report.txt`.

# Retired gamma RLE flag experiment (2026-10-29)

Full Gov2 sorted and unsorted comparison: see
`/Volumes/USB_1TB/BitMagic-experiments/2026-10-29/runs/inv_list/gov2-rle-flags-full-20261029-190624/report.txt`.
Both control and candidate streams verified and decoded 1,107,205 vectors
and 5,979,715,441 IDs. Sorted: 2,678,408,172 control bytes versus
2,678,208,892 RLE bytes, saving 199,280 bytes (0.00744%) after 1,014,985
RLE trials and 9,054 selected records. Unsorted: identical 5,441,839,210
bytes; 746,356 RLE trials selected no records. The size gain does not justify
the extra search, opcode, and decoding path. Experimental opcode 90 and its
encoder/decoder and CLI/harness have been retired. Grouped flags (opcode 92)
remain the adaptive superblock format. Existing released opcodes are unchanged.

# Current fully adaptive superblock validation

The 20261028-184913 full sorted cleanup validation passed and exactly matched
C: 2,678,868,407 bytes, 14,132 adaptive superblocks, 1,977,776 saved bytes.
Routine output now contains the summary and four compact step rows. Pass
`-compression-detail` to bminv for bit accounting and histograms.
The shared delta template retains the 64-entry superblock minimum; the unsuccessful adaptive GAP-XOR experiment has been removed.

`-sblock-adaptive` now always uses configuration C: select 16/32/64/128
position windows and absolute or signed previous-relative gamma corrections.
It uses opcode 92 with the existing step/coder descriptor, grouped change
flags, one final bit-writer flush, and strict complete-record fallback to
ordinary BIC. Existing opcode-92 formats still decode. Forced step/coder
APIs, command-line options, gate diagnostics, global-relative emission and
A/B comparison runners have been removed.

From `utils/inv_list`, build separately with `make neon`, then run:

```bash
bash experiments/sblock_adaptive/run_gov2_adaptive_validate.sh \
  --bin ./bminv --input-dir /Volumes/USB_1TB/int_set_bench/Gov2Flat \
  --datasets sorted --threads 10 --progress 10000 \
  --out "/Volumes/USB_1TB/BitMagic-experiments/2026-10-28/runs/inv_list/gov2-sorted-adaptive-clean-$(date +%Y%m%d-%H%M%S)"
```

Defaults to the full collection. It uses BWT/step2, 16-byte cutoff, GAP-XOR,
and 100% grouping tolerance. USB mount and free space are checked. Every
stream is verified against raw IDs and decoded, then validated by the report
before its `.bv` is deleted (use `--keep-streams` to retain it). Reports,
indices and JSON remain. Raw input stays local.

Full sorted C reference (20261028-164917): 1,107,205 vectors,
5,979,715,441 IDs, 2,678,868,407 bytes, 14,132 selected superblocks,
1,977,776 adaptive bytes saved. A/B/C all passed. These are historical reference values, not an exact-size acceptance gate.

Focused tests: compile `test_sblock_adaptive.cpp` and `test_sblock_delta.cpp`
with `c++ -std=c++17 -O2 -pthread -I ../../src` and run them.

## Historical experiments

The material below records previous experiments; removed comparison commands
are not supported by the current binary.

# Adaptive superblock BIC experiment

The optional `-sblock-adaptive` mode uses opcode 92. Opcode 91 and its
interleaved change-flag layout were experimental and have been retired.
The grouped layout is tried only for
sparse-superblock BIC arrays of at least 64 positions with at least one change
in their 32-position window correction. Each window includes its incoming
gap. The first window uses the global minimum gap minus one; later windows
have change flags stored in eight-flag groups before the gamma-coded
`correction + 1` values. An all-one group uses one marker bit; other full
groups use a marker plus eight raw bits. Partial groups are raw. The encoder
skips the trial when these groups cannot reduce flag bits. This metadata
precedes the BIC array. The BIC upper bound is the original last
position minus the sum of all applied corrections. The decoder restores the
original positions after BIC. Released superblock opcodes remain decodable.

The serializer encodes the ordinary record first and accepts opcode 92 only
when its complete record is smaller. Ties retain the ordinary record. The
trial report compares the emitted collection with the verified GAP-route
survey's indexed control, checks the per-record savings against the complete
stream difference, and requires full raw-input verification and decode.

Build separately with `make neon` in `utils/inv_list`. Run
`bash experiments/sblock_adaptive/run_gov2_hflags.sh --help` for
options. The experiment writes streams and logs only to USB_1TB; raw input
and source stay in the BitMagic repository.

`--gap-xor` combines the adaptive superblock selector with individual GAP-XOR.
The sparse-superblock decision happens first. When it is selected, its blocks
are not visited individually. Otherwise each GAP block may use GAP-XOR. The
report checks that the full-stream gain equals the sum of both selected gains.

`run_gov2_sblock_eligibility.sh` tests whether relaxing the GAP-length
superblock rejection admits useful additional groups. The serializer now uses
100% by default; the legacy threshold is 14% and the older-control runners
set it explicitly. `--gap-tolerance 25` permits a smaller relaxation, while
100% admits every otherwise eligible group below the existing 65,536-position
limit. This is a routing experiment; it
does not add a codec or change the decoder. The report compares the full
verified stream with the combined adaptive + GAP-XOR control, because newly
grouped blocks no longer take the individual GAP-XOR route.

For Gov2.unsorted, `run_gov2_unsorted_eligibility_pair.sh` builds a matching
14% combined control, verifies and decodes it, then runs and verifies the
relaxed stream. The paired report isolates the tolerance change. Both streams
are kept on USB_1TB for inspection; source and raw input remain local.

## Completed window/coder comparison

The full sorted A/B/C run at `gov2-sorted-delta-compare-20261028-164917`
verified all configurations. A (32/absolute) emitted 2,679,110,703 bytes;
B (32/absolute-or-previous) emitted 2,679,043,023; C (adaptive windows and
coder) emitted 2,678,868,407. C is retained. Removed variants and their
command-line controls are not supported. The decoder retains all previously
supported opcode-92 layouts. Full logs remain on USB.
