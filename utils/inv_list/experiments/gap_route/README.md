# GAP route survey

`-gap-route-stats` adds accounting only. It emits the same level-6 stream as
`-block-select-step2 -periodic-step1-interleaved 4 -bwt-min-gain 16`, with
the diagnostic disabled by default. The survey counts optimized GAP blocks
absorbed by sparse-superblock BIC records and GAP blocks written individually.
It records the exact emitted bytes for each sparse-superblock record and for
each individual GAP/BIT record. Sparse-superblock bytes are a group total:
they cannot be apportioned reliably among its constituent blocks.

The final report contains route counts and bytes, four GAP levels, and six
level-0 GAP-length buckets (lengths in 16-bit words). The grouped GAP counts
plus direct GAP counts can be reconciled with the optimized-block scanner.
`Other framing/spans` includes stream headers, empty/full-block spans, and
any record bytes outside the three measured routes.

Build and execution are separate. From `utils/inv_list`, run `make neon`,
then run `experiments/gap_route/run_gov2_gap_route.sh --help`. The harness
requires an attached USB_1TB drive, writes the stream and logs there,
verifies every requested vector against the raw input, decodes every BLOB,
and produces `report.txt`. Omit `--max-vectors` for the full collection.

For the BIC superblock length survey, use the already written `survey.bv` and
`survey.idx` with `bminv -decode -sblock-len-stats -threads 10 -bv-index-in`.
The decoder prints 64 length intervals, the exact minimum, maximum, mean,
and 75th percentile. `report_sblock_lengths.py RUN_DIR` validates the decoded
record count and byte sum against `report.txt`, then writes
`sblock-length.csv` and a short summary. Length is the number of set-bit
positions given to the BIC coder, not the number of resident GAP blocks.
