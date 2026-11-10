# CLI retirement status — 2026-11-10

Retired from bminv: `-bwt`/`-bwt-sample`, `-bwt-full`/`-bwt-diag`, standalone
`-periodic` and its order/search/budget controls, `-block-select`,
`-block-select-compose`, `-block-select-multistep` and multistep controls,
`-periodic-step1-start-reduce`, `-periodic-step1-first-start`, and ineffective
`-diag`/`-d`. Parser state, validation, policies and their reporting were removed.
Old options fail before output is opened. Library codecs/decoders were not removed.

Retained intentionally: `-block-select-step2`, `-periodic-step1-interleaved`,
codec selectors, `-gap-route-stats` and `-sblock-len-stats`. Opcode playlist surveys,
GAP-XOR/routing and sparse-vector comparison harnesses still use these controls.
Do not replace their historical explicit policy with a preset without checking
complete-record parity. Normal runs use levels 0–10 and optional BWT gain.

Cross-bvector research is preserved under `xbv` filenames. Codec fixtures,
including the independent wave dry-score oracle, remain available. Historical
notes describe past commands, not the current supported CLI.
