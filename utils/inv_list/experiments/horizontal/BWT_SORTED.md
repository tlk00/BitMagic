# Sorted Gov2 whole-block BWT measurements

Build separately, from the main BitMagic repository:

```bash
cd /Users/anatoliykuznetsov/dev/BitMagic/utils/inv_list
make
```

Then run the first 10,000 vectors with eight workers:

```bash
caffeinate -i bash experiments/horizontal/run_gov2_sorted_bwt_full.sh
```

The wrapper reads local `/Volumes/USB_1TB/int_set_bench/Gov2Flat/gov2.sorted` and creates a timestamped output directory under `/Volumes/USB_1TB/BitMagic-experiments/2026-10-23/runs/inv_list/`. It checks the USB mount at startup and the main harness checks output-volume availability before each phase. Existing output directories are refused. Source and scripts stay in the main repo. The harness does not build or resume a failed run.

Defaults: 10,000 vectors, progress every 25, level 6, eight compression workers, batch window 512, and every eligible BIT block evaluated. GAP blocks retain their normal path. Both baseline and BWT output are retained, verified against the original input, and decoded. Override limits for a smoke run with `--max-vectors 3 --threads 2 --batch-vectors 2 --progress 1`.

The BWT candidate uses one cyclic transform of the entire 65,536-bit block. It passes the result to the existing level-6 `encode_block()` selector, which includes BIC codecs. This experiment measures the existing path, without forcing BIC, adding a new codec, XOR search, or periodic search. The whole BWT record includes six bytes for its wrapper and primary index; ordinary encoding wins ties. `-bwt-diag` adds measurements only and is disabled by default.

Outputs:

- `report.txt`: actual stream savings, verification results, BIC-to-BIC savings, GAP statistics, and separate aggregate worker elapsed times for BWT construction and transformed-payload encoding. These worker timings are sums across workers, not end-to-end wall time.
- `metrics.tsv`: sizes, phase wall times and aggregate statistics.
- `blocks.tsv`: vector/block coordinates; baseline, transformed payload, framed candidate and selected sizes; actual saved bytes; popcount; GAP counts before/after; primary index; baseline/payload opcode and BIC flags. This can be joined with periodic logs by vector/block coordinates. A future best-of-codecs calculation is not a serialized hybrid.
- `codec_transitions.tsv`: baseline/payload opcode pairs, accepted/rejected counts. Opcode constants are in `src/bmserial.h`.
- `run.txt`, per-phase logs/timings, and the two serialized `.bv` collections.

The report checks per-block sums against serializer aggregates and the whole-stream byte delta. It distinguishes payload improvements lost to the six-byte wrapper from actual wins, and GAP reductions that do not yield byte savings.

Validation performed during preparation: standalone make build; three sorted vectors encoded, verified and decoded in both modes; all 791 BWT observations reconciled with aggregate counts and byte savings; diagnostic-on/off output bytes identical; single-thread and multithread BWT bytes and coordinate-bearing observations identical. This is a focused check, not the requested 10K benchmark or a full stress run.
