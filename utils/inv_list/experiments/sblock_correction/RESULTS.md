# Opcode-68 correction selection: prefix and worst-ratio playlists

## Unsplit metadata control, 2026-11-02

Run: `/Volumes/USB_1TB/BitMagic-experiments/2026-11-02/runs/inv_list/gov2-op68-unsplit-control-20261102-171452`.
The same gamma counts/endpoints/correction fields and CM BIC were used without
partitioning, omitting quarter occupancy metadata. All whole-stream baseline
sizes now exactly match bminv (byte-order/GAP-length headers disabled). Record
counts, IDs and bytes reconcile with the frozen survey; quarter results exactly
reproduce the preceding trial. All baseline vectors and both layouts of every
final opcode-68 record round-tripped with exact trial-record consumption.
The 76-case focused test passed normally and under ASan/UBSan.

| Dataset | Unsplit-only saving with fallback | Quarter-only saving with fallback | Joint saving with fallback | Quarters add beyond ordinary + unsplit |
|---|---:|---:|---:|---:|
| sorted | 16 B (4 records) | 120,276 B | 120,284 B | 120,268 B |
| unsorted | 16 B (4 records) | 969,792 B | 969,796 B | 969,780 B |

Forced unsplit encoding is larger by 32,524 / 48,180 bytes than ordinary 68.
Thus this specific metadata replacement alone is unproductive. The modest
quarter gains arise from partitioning and its associated local boundary ranges,
rather than simply replacing ordinary metadata with these gamma fields.
Local min correction remains negligible. This does not rule out other metadata
encodings or partitioning designs. No new production opcode is justified yet.
Single-worker exploration times were 121.275 / 321.701 seconds including all
baseline work and verification; these are not production performance benchmarks.

Future higher-count routing is recorded separately in
[SUPERBLOCK_ELIGIBILITY.md](SUPERBLOCK_ELIGIBILITY.md), including verified current
65,535-ID cap and the user's opcode-implied count-base / 16-bit residual idea.
It requires a new population survey, not these currently eligible playlists.

## Four-quarter trial, 2026-11-02

Run: `/Volumes/USB_1TB/BitMagic-experiments/2026-11-02/runs/inv_list/gov2-op68-quarters-20261102-170008`.
Method and measurement format: [QUARTERS.md](QUARTERS.md).
Both largest-byte 10K playlists completed. The offline format stores quarter
occupancy and per-quarter counts/endpoints/corrections with CM BIC interiors,
sharing one bit accumulator. No production opcode or decoder was changed.

| Dataset | Final opcode-68 records | Strict-fallback saving | Target opcode saving | Playlist stream saving | Local-correction saving before fallback |
|---|---:|---:|---:|---:|---:|
| sorted | 11,550 | 120,276 B | 0.074206% | 0.038959% | 0 B |
| unsorted | 18,419 | 969,792 B | 0.148527% | 0.120601% | 264 B |

All baseline vectors round-tripped; all 29,969 candidates round-tripped with
exact record consumption, including losers. Target record counts, original IDs
and bytes exactly match the frozen reference. Histograms reconcile. Focused
76-case boundary tests also passed normally and under ASan/UBSan.
The standalone writer leaves byte-order headers enabled (one extra byte/vector);
projected stream percentages use the original bminv stream sizes.

Global correction is zero in every final record. Positive local corrections
occur in 1,577 sorted quarters (20,019 IDs) and 4,481 unsorted quarters
(2,647,550 IDs), but correction metadata consumes almost all their potential.
Only 66 unsorted candidate records improve from local correction, by 264 bytes.
Other trial savings combine partitioning and changed boundary metadata; an
unsplit control with identical metadata was not measured.

Forced splitting increases total target bytes by 173,408 sorted / 199,896
unsorted. Strict fallback wins on 3,654 / 9,104 records. These are modest savings
on selected byte-heavy workloads, not whole-collection estimates.

Mean final record lengths are 15,361.27 / 30,475.43 IDs, mean sizes
14,033.17 / 35,449.21 bytes. Length-histogram 75th-percentile upper bounds are
32,768 / 65,536 IDs (quantized, not exact percentiles). This workload does
exercise substantial superblocks. Single-worker exploration times were
119.333 / 316.823 seconds, including baseline encoding and all trial validation;
they do not estimate production encoder/decoder overhead.

Recommendation: do not integrate a split opcode based on these results.
If continuing, first compare identical metadata without splitting to separate
header improvements from partitioning. Four-quarter local min correction
itself has not demonstrated a useful gain here.

## Largest opcode-68 byte contributors, 2026-11-02

Run: `/Volumes/USB_1TB/BitMagic-experiments/2026-11-02/runs/inv_list/gov2-op68-high-bytes-20261102-165053`.
Independent 10K playlists rank total opcode-68 bytes per vector. All four
streams passed raw verification, independent decode and byte/ID accounting;
baseline target totals matched the frozen survey.

| Dataset | Stream bytes, both arms | Opcode-68 bytes | Full opcode-68 byte coverage | Savings |
|---|---:|---:|---:|---:|
| sorted | 308,727,389 | 162,083,062 | 35.7169% | 0 |
| unsorted | 804,135,317 | 652,939,027 | 29.5247% | 0 |

Sorted has 461,977,696 whole-vector IDs and 11,550 final opcode-68 records;
unsorted has 696,441,508 IDs and 18,419 final records. Opcode-68 weighted
bits/ID are 7.3083 and 9.3056. Before competing superblock formats, sorted
has 19,340 candidates: 19,339 with zero legal positive global correction,
one trial tied. All 18,516 unsorted candidates have zero correction.
Thus this experiment finds essentially no global correction to select on
the byte-heavy workload. It does not measure local correction opportunities.
Per-record length/size distributions were not collected by this run.

Encode baseline/selected: sorted 20.353/18.619 seconds, unsorted
47.590/47.129 seconds. These profile-inclusive single passes do not establish
a performance improvement. Equal sizes are verified; byte identity was not
independently compared after generated streams were removed by the harness.

Next proposed measurement: local minimum gaps and complete-record costs for
four address-space quarters on these retained playlists, including framing
and final competition with existing formats. No split codec is implemented.

## Full survey and reusable playlists, 2026-11-02

Persistent survey:
`/Volumes/USB_1TB/BitMagic-experiments/surveys/gov2-op68-survey-20261102-122955`

Replay:
`/Volumes/USB_1TB/BitMagic-experiments/2026-11-02/runs/inv_list/gov2-op68-playlist-eval-20261102-132317`

Both full surveys round-tripped every vector against its optimized original.
Full byte equivalents match the retained baseline: sorted 2,658,874,844 bytes,
unsorted 5,440,396,028 bytes, each containing 5,979,715,441 original IDs.
The sorted full opcode-68 efficiency is now known: 9.758833 bits/ID,
453,798,745 bytes and 692,757 records. Unsorted is 12.051954 bits/ID,
2,211,498,269 bytes and 2,163,601 records.

Ranked independently by exact target bytes/IDs, filtered at 256 target bytes
per vector. Sorted has 446,937 eligible vectors; unsorted has 1,057,462.
Each playlist contains the highest-ratio 10,000 vectors.

| Dataset | Baseline stream bytes | Selected bytes | Saved | Stream saving | Opcode-68 bits/ID before -> after |
|---|---:|---:|---:|---:|---:|
| sorted | 3,082,076 | 3,081,628 | 448 | 0.014536% | 18.547862 -> 18.545140 |
| unsorted | 2,745,924 | 2,690,656 | 55,268 | 2.012729% | 20.341521 -> 19.929097 |

All four streams passed raw-ID verification, independent decode and opcode
byte/cardinality accounting. Replay baseline opcode-68 bytes, IDs and record
counts exactly reproduce the frozen survey selections. Encode times were
3.513/3.595 seconds sorted and 3.448/3.442 seconds unsorted; these one-pass,
profile-inclusive timings do not establish a performance effect.

Sorted: 20,000 superblock candidates, 18,166 without positive global correction,
1,834 trials, 1,583 ties, 112 omission wins and 139 historical strict wins.
Unsorted: 20,000 candidates, eight without positive correction, 19,992 trials,
6,172 ties, 13,817 omission wins and three historical strict wins.
Neither dataset had a win from newly applying correction below the old gate.
All candidate savings survived final rival-codec selection. Omission savings
average four bytes per winning record, consistent with accumulator-word costs;
header/BIC/padding contributions have not been independently isolated.

The difficult playlists contain little overall byte mass: sorted covers 0.673%
of full opcode-68 bytes, unsorted 0.123%. Their whole vectors average only
132.214 and 107.206 IDs respectively. They are useful fast development cases,
not an estimate of population-wide gains. High bits/ID still favors relatively
small vectors despite the 256-byte filter. A large-byte-contribution playlist
would complement these before judging broader effects or four-quarter layouts.

Decision: keep correction selection optional pending broader measurement.
The evidence supports testing whether correction metadata pays for itself;
it does not support simply removing the historical gate and always correcting.
The survey and playlists are retained for additional stratifications. After
relocation, both 10K playlists were restored and their selected raw-record hashes
checked successfully. Frozen survey artifacts are preserved; replay references
and the library catalog point at the persistent directory.

## Earlier 10K prefix

Prefix run:
`/Volumes/USB_1TB/BitMagic-experiments/2026-11-02/runs/inv_list/gov2-correction-10k-20261102-102412`

Sorted saved 12 bytes (16 pre-competition); unsorted saved none. The prefix
covered 80.586% of IDs but only 0.274%/0.328% of full final opcode-68 records.
The unsorted prefix offered no positive global corrections. Its lack of gains
therefore did not characterize the sparse tail, as the playlist run confirms.
