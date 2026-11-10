# Experiment cleanup — 2026-11-07

Removed historical horizontal/BWT comparison runners and their coupled report
scripts, the completed opcode-retirement comparison harness, one-off editing
scripts, source snapshots, Finder metadata, and generated Python caches.

Retained findings, codec fixtures, current preset validation, generic profiling,
playlists, sparse-vector tools and targeted diagnostics. No BM codec or bminv CLI
implementation was changed. The shared progress renderer moved from horizontal/
to current/; active consumers were updated. Historical commands in research notes
remain as records rather than supported instructions.

Removed source/script artifacts:

- `bwt_wave128/report_bwt_retired.py`
- `bwt_wave128/run_gov2_bwt_retired.sh`
- `horizontal/move_gov2_experiments_to_usb.py`
- `horizontal/report_bwt_blocks.py`
- `horizontal/report_gov2_block_select.py`
- `horizontal/report_gov2_bwt_full.py`
- `horizontal/report_gov2_bwt_multistep.py`
- `horizontal/report_gov2_bwt_periodic.py`
- `horizontal/report_gov2_bwt_step2.py`
- `horizontal/report_gov2_bwt_threshold.py`
- `horizontal/report_gov2_bwt_xor.py`
- `horizontal/report_gov2_periodic.py`
- `horizontal/report_gov2_periodic_orders.py`
- `horizontal/report_gov2_periodic_pair.py`
- `horizontal/report_gov2_step1_first_start.py`
- `horizontal/report_gov2_step1_interleaved.py`
- `horizontal/report_gov2_step1_start_reduce.py`
- `horizontal/run_gov2_block_select.sh`
- `horizontal/run_gov2_bwt.sh`
- `horizontal/run_gov2_bwt_decode_compare.py`
- `horizontal/run_gov2_bwt_full.sh`
- `horizontal/run_gov2_bwt_multistep.sh`
- `horizontal/run_gov2_bwt_periodic.sh`
- `horizontal/run_gov2_bwt_step2.sh`
- `horizontal/run_gov2_bwt_threshold.sh`
- `horizontal/run_gov2_periodic.sh`
- `horizontal/run_gov2_periodic_orders.sh`
- `horizontal/run_gov2_periodic_orders_usb.sh`
- `horizontal/run_gov2_periodic_pair.sh`
- `horizontal/run_gov2_sorted_bwt_full.sh`
- `horizontal/run_gov2_sorted_periodic_full.sh`
- `horizontal/run_gov2_step1_first_start.sh`
- `horizontal/run_gov2_step1_interleaved.sh`
- `horizontal/run_gov2_step1_start_reduce.sh`
- `horizontal/scratch_sources/add_bwt_switches.py`
- `horizontal/scratch_sources/before/bmserial.h`
- `horizontal/scratch_sources/before/inv_list.cpp`
- `horizontal/scratch_sources/edit_decode.py`
- `horizontal/scratch_sources/edit_inv.py`
- `horizontal/scratch_sources/edit_serial.py`
- `horizontal/scratch_sources/fail-verify.sh`
- `horizontal/scratch_sources/fix_bwt_switches.py`
- `horizontal/scratch_sources/fixtures.py`
- `horizontal/scratch_sources/mock_bminv.py`
- `horizontal/scratch_sources/progress_edit.py`
- `horizontal/scratch_sources/progress_harness_edit.py`
- `horizontal/scratch_sources/test_bwt_core.cpp`
- `horizontal/scratch_sources/test_bwt_full.cpp`

## Release preparation — 2026-11-10

Retired stale CLI modes and policy statistics listed in RETIREMENT_CANDIDATES.md.
Renamed cross-bvector headers, scripts and pair-search parity test with `xbv`;
CLI options and research output schemas remain compatible. Consolidated the six
collection runners into `../run_collection.sh`; moved the complexity example to
`../readme.txt`. No datasets, codec implementations or codec tests were deleted.

Removed completed/obsolete trial sources and runners:

- `sblock_correction/run_large.py`
- `sblock_correction/run_polarity.py`
- `sblock_correction/quarter_trial.cpp`
- `sblock_correction/report_quarters.py`
- `horizontal/progress_console.py` (temporary relocation shim; no active runners
  found when removed)

Large-count/polarity runners called a no-longer-supported harness flag and their
old baselines predated presets. Findings and OP98_POSITIVE_BASELINE.json remain
historical evidence. Direct correction/large-count codec tests are retained.
