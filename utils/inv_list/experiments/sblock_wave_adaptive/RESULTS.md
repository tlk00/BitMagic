# Retained adaptive wave widths: full Gov2 results (2026-11-01)

Keep **90** (fixed 1024-bit waves) and **97** (adaptive 128..4096-bit waves).
Reserve the exhaustive adaptive-width search for the **highest, size-first
compression levels**. It remains opt-in; this is not a balanced default or a
completed redesign of numeric compression levels. Current maximum is level 6,
with the optional selectors explicitly enabled by the highest-mode harness.

Evidence: /Volumes/USB_1TB/BitMagic-experiments/2026-11-01/runs/inv_list/gov2-adaptive-waves-20261101-094920

Each dataset: 1,107,205 vectors, 5,979,715,441 IDs; all four streams verified
against raw IDs and independently decoded. Opcode accounting and actual byte
reductions matched the exploratory projections exactly.

| Dataset | Baseline bytes | Adaptive bytes | Saved bytes | Saved % | Baseline encode s | Adaptive encode s |
|---|---:|---:|---:|---:|---:|---:|
| sorted | 2,662,267,916 | 2,658,874,844 | 3,393,072 | 0.127450 | 746.763 | 1066.829 |
| unsorted | 5,440,422,788 | 5,440,396,028 | 26,760 | 0.000492 | 1056.271 | 1405.147 |

Sorted weighted bits/ID: 3.561732 -> 3.557192. Unsorted: 7.278504 -> 7.278468.
Encoding overhead was 42.9% sorted and 33.0% unsorted in this one paired run.
Sorted decoding was 49.891 -> 64.854 s; unsorted 67.858 -> 58.554 s. Do not
claim a stable decoding effect from these inconsistent single-run timings.

97 selected 512,954 sorted and 5,011 unsorted records. Sorted selected widths
128/256/512/2048/4096: 176871/106781/78330/82976/67996. Unsorted:
380/63/37/1517/3014. All five alternative widths contribute. 1024 has no
97 wins because 90 supplies the same transformation without width bits.
Fixed 90 still selected 543,503 sorted records (1,031 unsorted), so retain it.

Future compression-level design: exhaustive width search is a size-first
option; investigate cheaper occupancy-guided selection and safe pruning before
considering a balanced default. Keep ordinary BIC maps; hierarchy/min-delta
maps were rejected and gamma runs deferred for their marginal extra savings.

Cleanup retires per-step/histogram and map/wave/endpoint measurement probes,
97 width/savings instrumentation, and the old paired measurement runners.
Keep wire codecs, exact dry scoring, strict fallback, score/write consistency,
focused codec tests and opcode profiling. Remaining compact codec/worker
summaries require -compression-diagnostics (--diagnostics in the harness).
Use ../current/run_gov2_highest.sh for post-cleanup paired-dataset validation.

## 2026-11-13 shared 97/108 selection

Single-worker L10 NEON -O2, sequential before/after runs on 10,000 distributed
vectors from each collection (every other ordinal in the earlier 20K sample).
Each dataset represents 55,980,282 IDs.

| Dataset | Encode before / after, seconds | Reduction | Decode before / after, seconds |
|---|---:|---:|---:|
| sorted | 91.714 / 70.760 | 22.8% | 1.162 / 1.150 |
| unsorted | 120.459 / 96.810 | 19.6% | 4.233 / 4.147 |

Streams, indexes and opcode profiles matched exactly. Raw-ID verification,
independent decoding, ASan/UBSan and 64-bit NEON focused tests passed.
Timings are one-pass observations. Decode differences do not establish a speedup.

[Benchmark report](/Volumes/USB_1TB/BitMagic-experiments/gov2-wave-shared-sample10k-20261113/comparison.txt)

The preceding full L10 trial emitted sorted 90/97/108 counts of
529,345 / 483,574 / 51,007 and unsorted counts of 658 / 3,555 / 2,128.
Old formats remain active candidates. Format 108 can represent their transforms,
but its width and codec fields add framing cost. Consolidation requires measuring
that cost and deciding whether to retain a compact fixed-width/BIC subformat.

These results precede the added reserved extension bit in opcode 108.
The revised prefix is opcode, three width bits, reserved zero bit, four layout
bits. Readers reject a nonzero reserved bit until an extension is defined.
Earlier unreleased 108 streams need their original decoder or regeneration.

## 2026-11-13 full consolidation and retirement

Both configurations validated 1,107,205 vectors / 5,979,715,441 IDs per dataset,
with 10 workers, raw-ID verification, independent decode and opcode accounting.
The comparison included opcode 108's reserved-zero extension bit.

| Dataset | Control bytes | Consolidated bytes | Penalty | Encode control/new, s | Decode control/new, s |
|---|---:|---:|---:|---:|---:|
| sorted | 2,642,196,054 | 2,642,702,642 | 506,588 B (0.019173%) | 1105.009 / 1114.541 | 23.029 / 24.506 |
| unsorted | 5,331,824,192 | 5,331,826,516 | 2,324 B (0.000044%) | 1510.480 / 1450.282 | 66.458 / 66.269 |

The consolidated streams emitted no 90/97 records. Opcode 108 counts were
1,032,542 sorted and 5,880 unsorted. Timings are single-pass observations.
The modest size penalty was accepted to remove two unreleased formats and
reduce implementation and selection complexity. Opcodes 90 and 97 are now
vacant; their writers, readers, routing and obsolete options are retired.
The retained 108 selector and transform helpers implement all six widths.

[Completed comparison](/Volumes/USB_1TB/BitMagic-experiments/gov2-wave-consolidation-20261113-170909/comparison.txt)

A new full levels-1–10 reference matrix is prepared for
`/Volumes/USB_1TB/gov2-reference-Nov-2026.v4`. Its runtime validation is pending.
Earlier reference directories remain historical artifacts and require their
original decoder for streams containing retired experimental formats.
