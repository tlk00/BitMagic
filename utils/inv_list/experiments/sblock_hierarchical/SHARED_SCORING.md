# Shared exact scoring, Gov2.sorted L8, 2026-11-14

The joint selector prepares corrected positions once per wave width, scores
ordinary position payloads once, and compares both complete map/header layouts.
Only the selected record is written. The cap-8 dictionary uses the same grouped
pairs. The standalone 108 path remains available without hierarchy overhead;
the 109-only path uses compacted coordinates and avoids unused flat maps.

All three full runs validated 1,107,205 vectors / 5,979,715,441 IDs. Before,
after and diagnostics streams are byte-identical by SHA-256 and direct cmp:
2,664,247,433 bytes (3.564380 bits/ID). Opcode counts and byte spans also match.

| Run | Encode s | Verify s | Decode s |
|---|---:|---:|---:|
| Before, independent scoring | 633.261 | 43.778 | 22.457 |
| After, shared exact scoring | 537.395 | 42.936 | 22.421 |
| After, TSV diagnostics enabled | 540.504 | 43.316 | 22.968 |

Observed encoding time decreased by 95.866 s (15.14%). Native ARM NEON,
ten workers, one pass per run. The baseline overlapped focused builds/tests;
cache and scheduling effects also limit interpretation of these timings.

## Strict comparisons

Search observations cover 1,757,456 groups before final emission selection.
Among all groups, 108 is strictly smaller in 1,033,601, 109 in 225,556, and
498,299 tie on complete padded bytes. The incumbent can beat both wave formats.

Among the 1,038,866 groups where a wave format strictly beats the ordinary/
adaptive incumbent, the results are:

| Outcome | Groups |
|---|---:|
| 108 strictly smaller | 739,107 |
| 109 strictly smaller | 34,819 |
| Equal padded size | 264,940 |

The emitted 1,004,047 opcode-108 records include those 264,940 ties.
All raw strict-108 comparisons collectively cost 109 another 60,964,071 map
bits. Its position payload is 21,179,858 bits cheaper on that same population,
but the map cost and padding outweigh that improvement. In 578,157 raw
strict-108 groups, both selected the same width and equal position-payload
costs, isolating the difference to maps, framing and padding.

## Fast predictor assessment

Observed occupied-block counts range from 2 to 174. The no-NULL condition
(256 occupied blocks) has no examples here. FULL blocks cannot occur in these
source lists because their cardinality is below 65,536.

No tested occupied-block threshold that skips any work is lossless.

| Hypothetical gate | Groups affected | Wrong choices costing bytes | Total byte cost | Largest group loss |
|---|---:|---:|---:|---:|
| Skip 108 when occupied blocks <= 4 | 70,607 (4.02%) | 7,036 | 28,144 | 4 |
| Skip 108 when occupied blocks <= 16 | 636,724 (36.23%) | 71,221 | 285,052 | 8 |
| Skip 109 when occupied blocks >= 64 | 310,821 (17.69%) | 4,811 | 39,364 | 2,768 |
| Skip 109 when occupied blocks >= 128 | 27,960 (1.59%) | 2,692 | 17,988 | 1,180 |

These are dry gate costs against exact joint selection with ordinary/adaptive
fallback, before any later selection. They measure byte loss, not runtime
saved. They are dataset observations and are not universal loss bounds.
The <=16 gate costs 0.010699% of this collection. It is a candidate for a
separate speed-mode experiment. Exact shared scoring remains the policy for maximum
compression. No heuristic gate has been enabled.

[Complete diagnostics and all occupied-block bins](/Volumes/USB_1TB/BitMagic-experiments/gov2-l8-wave-score-diagnostics-20261114/score-report.txt)

## Reproduction

Build separately with make neon in utils/inv_list. For the before/after runs,
use run_gov2_highest.sh with --datasets sorted --level 8 --threads 10
--wave-mode both --keep-streams, and fresh USB output paths. The after-run
can use --reference BEFORE_ROOT to require matching stream hashes.

For TSV diagnostics, also pass --wave-score-out
OUT/sorted/profile/wave-scores.tsv. The reporter accepts that path through
experiments/sblock_hierarchical/report_scores.py. Scope 0..5 denotes each width;
scope 6 compares each format's best width. The optional callback belongs to
the serialization support code. The thread-safe TSV sink belongs to inv_list.

The independent scorer oracle regression passed 368 decisions under
ASan/UBSan, covering both policies, tie/rival handling and near-limit counts.
The exporter smoke test reconciled each component sum to the complete record.
