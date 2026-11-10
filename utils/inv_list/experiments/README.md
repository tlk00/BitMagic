# inv_list experiments and validation

Build utilities separately; harnesses consume an existing binary.

| Directory | Retained purpose |
|---|---|
| [current](current/README.md) | Current compression presets, full Gov2 verification, shared progress renderer and focused CLI/progress checks. |
| opcode_profile | Generic opcode size and bits-per-ID reporting. |
| opcode_playlist | Survey, selection and replay of representative vector playlists. |
| block_stats, gap_route | Collection composition and routing diagnostics. |
| sparse_vector | Separate sparse-vector collection experiments. |
| block_endpoint, gap_xor, sblock_adaptive, sblock_correction, sblock_wave_adaptive | Retained codec tests, targeted diagnostics and comparison tools. |
| horizontal, bwt_wave128 | Historical findings and retained C++ fixtures; collection runners retired. |
| block_small_wave, block_wave, gap_bwt, sblock_map, sblock_probe | Historical findings. |

Historical Markdown files can contain commands for removed tools and old format
variants. Use `current/README.md` for the supported full-collection workflow.
Generated collections and run logs belong on the USB experiment drive, outside
this source directory.

See [the cleanup record](CLEANUP.md) for the removed artifacts.
