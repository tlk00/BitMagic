# Opcode-68 global correction selection experiment

Opt-in `-sblock-correction-select` compares the historical opcode-68 record
against the opposite global-correction choice. Correction is the existing
minimum gap minus one, capped at 65535. The historical policy applies it only
when `(length - 1) * correction > 250`. Zero correction needs no second trial.
Both complete records use the existing writer, metadata and final accumulator
padding. Strictly smaller records win; ties preserve historical bytes. This
initial experiment writes trial records; it is not a dry-score performance
optimization. Old opcode definitions and decoders are unchanged.

The selected opcode-68 candidate then competes with retained 90/92/97 formats.
Counters describe pre-competition decisions: no correction possible, trial,
tie, historical strict win, newly applied correction, newly omitted correction,
and candidate savings. Final collection savings are measured separately and
must not exceed those candidate savings. No superblock splitting is included.
The option is disabled by default and isolated for later removal or promotion.

## Run

From utils/inv_list, build separately:

```bash
make neon
bash experiments/sblock_correction/run_gov2_correction.sh \
  --bin ./bminv --input-dir /Volumes/USB_1TB/int_set_bench/Gov2Flat \
  --threads 10 --progress 1000 --max-vectors 10000 \
  --out "/Volumes/USB_1TB/BitMagic-experiments/2026-11-02/runs/inv_list/gov2-correction-10k-$(date +%Y%m%d-%H%M%S)"
```

Runs baseline and selected modes on both sorted and unsorted prefixes with all
retained highest-compression selectors. Each of the four streams receives raw-ID
verification, independent threaded decode, opcode byte/cardinality profiling and
reconciliation. Generated streams are deleted after successful child reporting;
logs, hashes, indices and JSON remain. Build/source/input files stay local.
Report: `report.txt`; machine-readable comparison: `comparison.json`.
Each arm also has its complete opcode structure report under
`MODE/DATASET/profile/report.txt`.

The comparison shows prefix opcode bytes, share and bits/original-ID against
available full references. Full unsorted per-opcode ID counts are available;
the older sorted full reference has byte shares and collection bits/ID only.
Unavailable per-opcode efficiencies show n/a. The reference paths are recorded
in comparison.json; unavailable USB references do not fabricate a comparison.
A 10K prefix can contain a large fraction of IDs because vector lengths vary.
These reports diagnose skew; they cannot establish that a prefix represents
all omitted vectors. No claim of population-wide gain should precede full runs.

Focused test from repository root:

```bash
c++ -std=c++17 -O1 -g -DBMNEONOPT -fsanitize=address,undefined -Isrc \
  utils/inv_list/experiments/sblock_correction/test_correction.cpp \
  -o /tmp/test_sblock_correction
/tmp/test_sblock_correction
```

Validation completed: NEON build; ASan/UBSan direct min-cost/tie/legacy-decode
checks over 240 generated fixtures, including four strict improvements; whole
serializer tests with and without rival codecs; opcode/cardinality reconciliation;
and baseline/selected two-vector sorted and unsorted encode/raw-verify/decode
runs with report generation. Those Gov2 smoke prefixes had no superblock
candidates; selector coverage comes from the synthetic tests. The disabled-mode
unsorted smoke stream matches the earlier retained profiler baseline byte-for-byte.
This is focused validation; the 10K comparison remains for the user run.

## Distributed follow-up (preferred next experiment)

The 10K-prefix comparison at
`/Volumes/USB_1TB/BitMagic-experiments/2026-11-02/runs/inv_list/gov2-correction-10k-20261102-102412`
verified all four streams. Sorted saved 12 bytes; unsorted saved none.
Sorted had 26 positive-correction trials among 7,806 candidates, four omission
wins (16 pre-competition bytes), 21 ties and one historical strict win.
All 7,355 unsorted candidates had no positive global correction.
The prefix covered 80.586% of IDs but only 0.274%/0.328% of final opcode-68
records in sorted/unsorted. It contains every final unsorted opcode-65 record.
These findings do not settle correction selection for the sparse-vector tail.

The follow-up samples 10,000 **whole vectors** using zero-based indices
`floor((2*i+1)*total_vectors/(2*sample_count))`, one midpoint per equal index
stratum. Sorted and unsorted use identical indices. It does not sample blocks
or IDs within vectors, and does not force the largest endpoint vectors into
the sample. No inverse-probability or population-size projection is claimed.
The actual sample counts, opcode shares and efficiencies are reported.

```bash
make neon
bash experiments/sblock_correction/run_gov2_correction_sample.sh \
  --bin ./bminv --input-dir /Volumes/USB_1TB/int_set_bench/Gov2Flat \
  --sample-vectors 10000 --threads 10 --progress 1000 \
  --out "/Volumes/USB_1TB/BitMagic-experiments/2026-11-02/runs/inv_list/gov2-correction-sample10k-$(date +%Y%m%d-%H%M%S)"
```

Source Gov2 files remain local. Sample copies are generated under `OUT/sample`
on USB. They and `sample_manifest.json` are retained for future experiments
(including the later four-quarter superblock trial). The manifest records
source sizes/mtimes, selected indices, original offsets, per-vector counts and
sample SHA-256 hashes. Source framing is scanned completely; copied record
bytes are preserved; paired counts and sample verification totals must agree.

The four current-codec runs are under `OUT/comparison`; the final combined
table is also copied to `OUT/report.txt`. The table includes per-opcode total
bytes and bits/ID, full-reference comparisons and opcode-68 record/byte
coverage. Sampled input bytes are additional disk usage; generated serialized
streams are deleted only after successful verification/reporting. No compression
codec change is introduced by this sampling follow-up.

Sampling validation: synthetic paired-record tests cover exact midpoint indices,
byte-preserving copies, all-vector sampling, invalid counts and truncated input.
A real ten-vector distributed sample scanned both full Gov2 inputs (1,107,205
records each), selected matching 7,176-ID samples, and passed all four
baseline/selected sorted/unsorted raw verification, independent decode and
manifest/opcode report reconciliation runs. This is a harness smoke test, not
the requested 10K experiment.
