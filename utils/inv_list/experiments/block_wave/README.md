> Historical experiment notes. The exploratory measurement code and old
> runners described below were retired on 2026-11-01. Use
> [current highest-compression verification](../current/README.md).
> Retained codec correctness tests remain where applicable.

# Opcode-65 empty-wave measurements

Build separately in `utils/inv_list` with `make neon`. Run:

```bash
bash experiments/block_wave/run_gov2_block_wave.sh --bin ./bminv \
  --threads 10 \
  --out /Volumes/USB_1TB/BitMagic-experiments/2026-10-30/runs/inv_list/gov2-unsorted-block-wave-$(date +%Y%m%d-%H%M%S)
```

No vector limit means the full collection. `--max-vectors 10000` limits a trial
run to the first 10K vectors. This is a prefix, not an evenly spaced sample.
The input stays in `/Volumes/USB_1TB/int_set_bench/Gov2Flat/gov2.unsorted`.
The harness checks the USB mount, enables current BWT/step2 (16-byte cutoff),
GAP-XOR and superblock 68/92/90 selection, verifies raw IDs, independently
decodes the collection, then prints `block-wave-report.txt`.
Verified `.bv` streams are deleted unless `--keep-streams` is specified.
Logs, indices, profiles and reports remain on USB.

The new `-block-wave-trial` option is purely observational. A per-worker index
pass locates FINAL opcode-65 records, including their exact record byte spans.
It does not measure speculative serializer candidates or embedded records.
The report reconciles its count and baseline bytes against the independent
main-thread opcode profiler. No existing opcode writer or decoder is changed.

For each block, try 32/64/128/256/512/1024-bit waves. Skip a wave size only when
all its waves are occupied. Pack occupied waves, retaining bit offsets within
each wave, then subtract `i * (minimum adjacent gap - 1)` from position i.
For every size compare occupied/empty BIC maps and derived/explicit bounds.

Hypothetical record (scratch memory only):

- One opcode byte (placeholder 255; never emitted and no opcode allocated).
- Position count via delta16s, two mode bits, gamma(correction + 1).
- Gamma map count and BIC wave indexes over `[0, 65536 / wave_size - 1]`.
- Optional gamma(first + 1) and gamma(derived_upper - last + 1).
- BIC u16 positions, or interior positions with explicit endpoints.
- One shared bit accumulator, flushed only at the end.

The derived upper bound is `occupied_count * wave_size - 1 - (n - 1) * correction`.
The correction-only control omits the map and uses one boundary-mode bit and
`65535 - (n - 1) * correction` as its bound. Its header follows the proposed
format, not opcode 65's existing endpoint heuristic/metadata.

Fixed-size scores have no size identifier. Adaptive wave selection adds THREE
bits to each candidate before rounding the accumulator to complete 32-bit
words, then chooses the smallest candidate with ordinary-byte fallback.
Strict byte savings only; ties keep ordinary encoding. The correction-only
control is reported separately and is not included in adaptive wave selection.

All savings are projected, including actual word padding. Component totals are
for the smallest candidate per eligible block, including losers; savings totals
count winners only. Worker seconds are summed steady-clock elapsed durations
of evaluation under parallel load, not isolated CPU benchmarks. The probe also
adds an index pass and decodes opcode-65 arrays; these costs are included in
collection wall time but excluded from individual candidate timings.

Focused test (includes independent scratch decoding and no-mutation checks):

```bash
c++ -std=c++17 -O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer \
  -I ../../src experiments/block_wave/test_block_wave.cpp -o /tmp/test_block_wave
/tmp/test_block_wave
```

These are measurement helpers isolated under `experiments/block_wave` for easy
retirement. The reusable record observer belongs to `bmserial_profile.h`.

## Deferred BWT experiment

Measure empty-wave correction on the transformed bit-block inside emitted BWT
records. Compare a fixed 128-bit wave with an adaptive wave-size choice, using
the original BWT record as the baseline. Include the wave map, size selector,
BWT primary index, framing, and final padding in complete-record costs. First
profile eligible blocks and projected wins; only add an encoder/decoder if the
verified gain justifies the added decode work. Keep the existing 16-byte BWT
selection cutoff when comparing collection-level results. The opcode-65 wave
probe above does not answer this question because it sees different blocks.
