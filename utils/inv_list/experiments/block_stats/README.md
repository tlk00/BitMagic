# Gov2 optimized bvector block composition

Build separately from the scan:

```bash
cd /Users/anatoliykuznetsov/dev/BitMagic/utils/inv_list
make gov2-block-stats-neon
```

Run the full sorted collection (or set `--max-vectors 10000` for a prefix):

```bash
bash experiments/block_stats/run_gov2_block_stats.sh
```

The scanner reads the raw length-prefixed u32 lists, validates sorted unique
IDs, constructs each `bvector<>` with `BM_SORTED`, calls `optimize(tb)`, then
collects `calc_stat()` and exact GAP lengths from the resident blocks. This
matches the reconstruction path used by `bminv`; it does not re-optimize
previously serialized bvector streams. Generated statistics go to the mounted
USB_1TB USB, while source and raw input remain in the BitMagic checkout.

`report.txt` gives BIT/GAP/full/zero counts, fractions, the four GAP levels,
and per-vector percentiles. `gap_lengths.csv` is an exact GAP-word-length by
level histogram; `per_vector.csv` records each list's counts. Implicit full
and zero counts are measured within the block span from zero through each
vector's last ID. BIT and GAP percentages use both materialized and logical
denominators. GAP word counts and BIT bytes describe resident structure,
not serialized output bytes.
