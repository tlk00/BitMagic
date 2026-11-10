# Gov2 sparse-vector collection experiment

Build separately from the run:

```bash
cd /Users/anatoliykuznetsov/dev/BitMagic/utils/inv_list
make neon svexp-neon
```

The runner processes the complete sorted collection by default:

```bash
bash experiments/sparse_vector/run_gov2_sv.sh \
  --input /Volumes/USB_1TB/int_set_bench/Gov2Flat/gov2.sorted \
  --bv-reference-dir /Volumes/USB_1TB/BitMagic-experiments/2026-10-23/runs/inv_list/gov2-sorted-step1-gamma-fixed-full-20261025-233445 \
  --threads 10 --batch 512 --progress 1000
```

For unsorted Gov2, replace `gov2.sorted` with `gov2.unsorted` and omit
`--bv-reference-dir`; the runner builds, verifies, and decodes a matching
BV reference first. Use `--max-vectors 1000` for a preliminary bounded run.

Each input posting list becomes either a sparse vector of raw document IDs or
a sparse vector containing the minimum gap, first ID, and successive gaps
minus that minimum. For each representation, the experiment emits plain,
intra-vector plane-XOR, and BWT/periodic-selection streams. XOR and BWT/periodic
selection are separate modes. Selection uses compression level 6, the
gamma-versus-fixed step-1 formula, independent BWT residual step 2, and a
16-byte minimum BWT gain measured on complete block records. There is no
inter-vector reference search.

Each `*.sv` is a concatenation of self-delimiting sparse-vector records.
The accompanying `*.idx` identifies record boundaries for parallel verify
and decode. The report gives both stream bits per original ID and stored bits
per original ID including the index. Both metrics divide total collection
bytes by the total original IDs, so vectors are weighted by their sizes.
The final table includes absolute bytes, both savings comparisons, and
encode/verify/decode elapsed times. Every sparse
stream must reconstruct the original posting lists and decode successfully
before the report prints `Correctness: PASS`.

The runner checks that the USB_1TB USB is mounted and stops starting a new
stream when fewer than 4 GiB remain. Source files, scripts, and raw Gov2 input
stay in the BitMagic checkout; generated streams and logs go to the USB.

## Full Gov2.sorted result (2026-10-26)

All six SV streams verified against the input and decoded all 1,107,205
posting lists (5,979,715,441 original IDs). The reused, independently
verified BV reference is the full sorted collection under
`gov2-sorted-step1-gamma-fixed-full-20261025-233445`. Results and phase logs:
`/Volumes/USB_1TB/BitMagic-experiments/2026-10-23/runs/inv_list/gov2-sorted-svexp-20261026-113739`.
Bits/ID below divides stream bytes by original IDs; the index adds about
0.0163 bits/ID to every variant.

| Variant | Stream bytes | Bits/ID | Encode seconds |
| --- | ---: | ---: | ---: |
| BV reference | 2,690,277,862 | 3.5992 | 637.53 |
| SV raw, plain | 5,543,268,989 | 7.4161 | 77.41 |
| SV raw, XOR | 4,762,811,379 | 6.3720 | 139.14 |
| SV raw, selection | 4,598,312,369 | 6.1519 | 2227.74 |
| SV gaps, plain | 3,192,289,608 | 4.2708 | 38.75 |
| SV gaps, XOR | 3,173,991,938 | 4.2463 | 79.22 |
| SV gaps, selection | 3,190,840,643 | 4.2689 | 884.76 |

Minimum-gap reduction is the dominant SV improvement. With that representation,
plane XOR saves 18,297,670 bytes (0.57%) versus plain SV; block selection saves
1,448,965 bytes (0.045%) at much higher encode cost. The smallest SV stream,
gaps+XOR, is 483,714,076 bytes (17.98%) larger than the BV reference. Raw SV
benefits more from the codecs, but remains substantially larger than gaps SV.
The BV timing came from its earlier run, so treat cross-run timing comparisons
as indicative. Keep gaps+XOR as the sorted-SV control for further experiments;
these results do not justify further selector tuning on sorted SV alone.
