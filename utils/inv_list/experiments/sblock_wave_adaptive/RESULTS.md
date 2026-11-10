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
