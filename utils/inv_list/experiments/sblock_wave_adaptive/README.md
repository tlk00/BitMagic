# Adaptive wave-compacted superblocks

Opcode 108 (`set_sblock_bienc_wave_runs`) is the retained wave format at levels
8–10. It selects 128, 256, 512, 1024, 2048 or 4096-bit windows, compacts empty
windows, and applies global min-gap correction. Wave maps use BIC or gamma
runs; corrected positions use BIC or repeated run-length pairs. Map polarity
and BIC endpoint layouts are selected by complete-record dry scores.

The bit after the three-bit width selector is reserved and written as zero.
Readers reject one until an extension layout is defined. All payloads share
one accumulator and final padding. Other superblock codecs remain fallbacks.

Unreleased opcodes 90 and 97 and the paired consolidation runner were retired
following the full-collection comparison in [RESULTS.md](RESULTS.md). Their
numeric IDs are vacant. Historical streams containing them require the older
binary. Reports remain on USB_1TB.

The optional API `set_sblock_wave_options(true)` and CLI `-sblock-wave` enable
108 explicitly; normal callers select a compression level. The old fixed-wave,
adaptive-only, and replacement-trial options have been removed.

Current collection validation uses `experiments/current/run_gov2_highest.sh`.
Reference reconstruction uses `reference/run_gov2_reference.sh`. The retained
focused source `test_sblock_wave_adaptive.cpp` covers dry-score parity against
complete encoding, all supported layouts, reconstruction, ties and rejection
of an undefined extension bit. It is development code, excluded from bminv.
