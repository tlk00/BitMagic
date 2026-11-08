# TAgent

TAgent is a remote, unattended BitMagic build and correctness-test supervisor.
Implementation languages are Bash and Python (standard library only). Supported
initial targets are x86-64 Linux/WSL and Apple Silicon macOS, matching the current
`build_all.sh` scripts. Intel macOS and ARM Linux are not supported yet.

## Agreed requirements

- Start from an operator-cloned BitMagic repository. Watch its `origin/master`
  by default. Fetch periodically (default five minutes); any new commit,
  including documentation changes, triggers a complete cycle. Test the current
  tip on every manual startup by default, even if it was tested earlier. With
  `--startup changed`, test only when the commit or coverage settings differ from
  the last successful run. An automatic self-restart resumes its pending cycle.
  During perf/stress execution, check the remote at that
  interval without changing the checkout. If a new tip appears, finish the
  currently running perf processes or active stress groups, then fetch,
  fast-forward and start a fresh build/test cycle. Queued perf processes and
  stress groups are skipped while active jobs drain; remaining old-cycle builds
  and groups are skipped and intermediate commits may be skipped. No partial cycle is recorded
  as successful. Checks run at test process heartbeats (up to 10 seconds later).
- Update the containing checkout with fast-forward-only Git operations. Refuse
  tracked local modifications, divergent history, and conflicting untracked
  files; never reset, stash, or discard operator work.
- Put editable settings prominently at the top of `run_tagent.sh`; command-line
  options override them. Validate Python, tools/versions, platform, usable CPU
  features, writable storage and minimum free space before work. Recheck free
  space during operation. One instance may own a checkout at a time.
- Build through `tests/perf/build_all.sh` and `tests/stress/build_all.sh`. Build
  failures print a brief diagnostic and log location, signal distress and exit
  nonzero immediately.
- Default coverage is perf release variants and every debug/release stress
  variant of `tests/stress/t.cpp`. Exclude stress64. `--pedantic` enables stress
  serialization sweeps and debug/release perf SIMD variants as well.
- Run perf first, with a shared worker limit; timing comparisons are not an
  acceptance criterion. Then run every stress build using the existing parallel
  runner, one build at a time. All scheduled tests must succeed.
- Default workers: max(1, available logical CPUs minus two). CPU availability
  respects affinity and Linux cgroup quotas where detectable. Allow an explicit
  worker override. This limits processes, rather than pinning/reserving cores.
- Run binaries only when supported by the host CPU/OS. Log skipped instruction
  sets. AVX512 is opt-in. Record expected outputs explicitly; stale binaries
  cannot stand in for a missing build. Isolate test files per worker.
- On test failure, stop all remaining jobs, preserve logs, print the failing
  build/command/status, signal distress, and wait for a newer commit. Fetch with
  progressive delays (1, 2, 4, 8, then 10 minutes). Rebuild when one arrives.
  Quit nonzero if no fix arrives within the configurable recovery window
  (60 minutes by default). A successful fetch without new code does not reset
  that deadline. A newer commit starts a new attempt.
- Network failures use the same configurable recovery window and backoff. A
  recovered connection resumes normal operation. Other operational failures
  (storage, permissions, tools) print diagnostics and exit nonzero.
- Detect changes under `tests/tagent` during updates and automatically restart
  through the launcher, preserving arguments and persistent state. If restart
  cannot proceed, exit with operator instructions. Local editable launcher
  settings should instead be supplied through environment/CLI when unattended
  self-updates are wanted, since tracked edits prevent safe updates.
- After success, persist the tested commit and coverage settings and wait for
  changes. Periodic terminal status overwrites one line and shows last successful
  fetch, attempt count, current commit, next action and free disk. Redirected
  output uses timestamped lines. Print progress during builds/tests too.
- Sound plays on the test machine; `--no-sound` disables it. Prefer local audio
  playback, with a local console-speaker fallback. Missing sound never prevents
  testing. Events: startup beep, build-success cha-cha, workload diesel pulse,
  failure SOS three times, idle disappointment. Initial synthesized sounds are
  deliberately simple and can be refined later. Ctrl-C/termination stops all
  child process groups, preserves logs and releases the checkout lock.

## Operation

Requires Bash 3.2+, Python 3.9+, Git, make, and the C++ toolchain used by BitMagic.
No pip packages are required. Run on a dedicated clean checkout:

```sh
bash tests/tagent/run_tagent.sh --help
bash tests/tagent/run_tagent.sh --jobs 6
bash tests/tagent/run_tagent.sh --startup changed
bash tests/tagent/run_tagent.sh --pedantic --retry-minutes 60 --no-sound
```

Defaults are at the top of the launcher; environment variables of the same names
can override them. Logs, persistent state, generated sounds and isolated test
working directories default to `TAgent-logs/` at repository root (untracked; add
it to your local Git exclude file if desired). Each cycle has a timestamp/commit directory,
`manifest.json`, separate build/perf/stress logs, and stress-group logs.
`agent.log` records lifecycle events. Agent-created cycle folders older than seven
days are removed periodically, including failed and superseded cycles. Change
`LOG_RETENTION_DAYS` in the launcher or use `--log-retention-days N`; use a longer
retention if failure evidence must be kept. The active cycle, state, sound files,
other folders, and older `.tagent/` directories are never removed by cleanup.

`STARTUP` at the top of the launcher defaults to `build`; `--startup changed`
overrides it to reuse a previous successful result. `--force` still forces a
cycle even with `--startup changed`. `--once` performs the selected startup
policy and exits; with `changed` it can exit without testing. It still uses
network recovery but exits immediately on test failure. `--check` validates
local prerequisites without fetching, building or testing. A failed initial
fetch retries within the configured window and exits if recovery fails. Failure exit codes: 1 operational/build, 2 invalid
arguments, 3 test/recovery exhausted, 130 interrupted. Automatic self-update
restarts retain options. The running agent code is not hot-reloaded.

## Focused validation

```sh
python3 -m unittest discover -s tests/tagent -p 'test_*.py'
bash -n tests/tagent/run_tagent.sh tests/perf/build_all.sh
```

The supervisor tests use temporary local Git repositories, synthetic build/test
programs and the existing stress runner. They check recovery deadlines, saved
success state, fast-forward updates/self-restart dispatch, failure propagation,
worker isolation, process timeouts and Ctrl-C descendant cleanup. They do not
constitute a full BitMagic regression, native Linux/WSL validation, or an audio
quality check. First deployment should use `--check`, then an attended `--once`
cycle before unattended operation.

## Steampunk and dieselpunk sound modes

The default `dieselpunk` theme uses crankshaft rumble, cylinder knocks and an air
horn. During steampunk build warm-up, a mechanical spark telegraph sends
**"New BitMagic updates, new BitMagic updates"** in Morse, including the comma,
over an independent steam hiss that gradually swells in volume. The dry metal
key track sounds only at contact closure/release; held marks have no sustained
tone or keyed hiss. Dot/dash lengths are the time between the two impacts.
Steam and key use independent synthesis and are mixed without shared compression,
so Morse never modulates the steam. Key strikes have brief inharmonic metal
resonances, with the main strike around 2 kHz and softer upper partials. The dot
lasts 105 ms; the full message lasts about 38.4 seconds, with playback waits allowing it to finish. This plays at
the start of each perf/stress build script, including the initial cycle; it is
not a separate notification for each fetched commit.

`steampunk` uses exhaust chuffs, steam, mechanical bells and a steam-whistle
chord. Both include a build-up, build-success cha-cha, descending idle sigh,
victory horn and Morse SOS repeated three times. All effects and the combined
demo are 44.1 kHz mono PCM and need no extra Python packages. The dispatch
message is a bundled mono voice/steam/wheel recording.

```sh
# Demo only: no Git checks, fetching, builds or tests. Ctrl-C stops playback.
bash tests/tagent/run_tagent.sh --sound-demo --sound-theme steampunk
bash tests/tagent/run_tagent.sh --sound-demo --sound-theme dieselpunk

# One sound, or normal operation with a chosen theme and volume (0..1).
bash tests/tagent/run_tagent.sh --sound-preview depart --sound-theme steampunk
bash tests/tagent/run_tagent.sh --sound-preview depart-perf --sound-theme steampunk
bash tests/tagent/run_tagent.sh --sound-preview depart-stress --sound-theme steampunk
bash tests/tagent/run_tagent.sh --sound-preview dispatch --sound-theme steampunk
bash tests/tagent/run_tagent.sh --sound-theme dieselpunk --sound-volume 0.25
```

`SOUND_THEME` and `SOUND_VOLUME` are editable launcher defaults with environment
and CLI overrides. `--no-sound` disables all playback. Audio targets the test
machine, trying afplay, PipeWire, PulseAudio and ALSA, then the local console bell.
Unavailable audio does not prevent testing. Generated WAVs remain in the state
directory for inspection.

In steampunk mode, the perf and stress phases each get a prerecorded female
station announcement, with quiet radio static, followed by 25 seconds of Vasily
Agapkin's march "Farewell of Slavianka" through a station-speaker filter and the
existing 23-second departure: platform bell, steam whistle, then accelerating
chuffs and bogie knocks as the train pulls away. Each sequence lasts about 54
seconds and needs no TTS engine. Testing starts without waiting for playback;
the announcements are played once per phase, not per worker or build.
The mono music excerpt is from the [United States Coast Guard Band recording](https://commons.wikimedia.org/wiki/File:Farewell_of_Slavianka_(USCGB).ogg),
identified there as a U.S. government public-domain performance, with the
composition also listed as public domain in the United States. The source audio
was truncated, resampled, filtered, faded and attenuated for this effect. The
extended departure doubles bell decay, whistle and acceleration timing while
preserving tonal pitch. During active tests, every five minutes TAgent
chooses independently with equal probability between the passing freight train
and the male dispatch message with quiet steam and wheels. The first plays about
five minutes after tests start (heartbeat granularity up to ten seconds).
The train clock spans perf and stress builds and stops when the cycle ends; no
passing trains play during builds, idle or recovery. Other sound events do not
reset that clock. Diesel work clips retain their 15-second interval.
Steampunk clips last 10.5 seconds: a heavy freight train approaches, passes the listener, and recedes. Prominent layered cylinder chuffs and steam breaths lead the slower paired bogie knocks
and low rolling rumble of the carriages. Their cadence is separate from the
exhaust; mild pitch/cadence change and distance envelopes suggest the passing
train. This is a procedural impression, not a recording or physical simulation.
Diesel clips last three seconds. Their
pulse rate follows active perf processes or the stress runner's START/PASS/FAIL
counts relative to the job budget; it does not measure CPU utilization. Clips
play in background threads and are queued to avoid overlapping sounds. On failure,
queued non-emergency sounds are discarded, and the distress signal follows the
clip currently playing. The demo includes both low and full
worker occupancy; use `--jobs 8` if demonstrating on a single-worker machine.

The steampunk demo includes both phase departures, two passing trains, the dispatch
message, clear/caution build reports, and successful/failed build summaries. The random choice applies to periodic test sounds, not to demo previews.


## Build dispatcher

TAgent builds perf and stress concurrently when the worker budget is at least
two; with one worker it builds them sequentially. Each family keeps its own log.
If either build fails, TAgent waits for the other build to finish, reports the
combined diagnostics and exits without starting tests.

While builds run, TAgent prints a status report after 60 seconds, then after
another 70 seconds, 80 seconds, and so on. In steampunk mode the dispatcher
announces it too. `BUILD_REPORT_SECONDS` in the launcher or
`--build-report-seconds N` changes the initial interval; each later interval
still increases by 10 seconds, and a new build cycle resets it. Checks run at
subprocess heartbeats, so reports can be up to about 10 seconds later. Each
report shows the ready and pending expected perf and stress build targets;
ready means the executable has been produced. The report also counts complete
GCC/Clang-style `warning:` diagnostic lines across the current perf and stress
build logs so far; repeated warnings in separate variants count separately.
Compiler summary lines such as `3 warnings generated` are not counted. Other
warning formats may not be recognized; the currently unfinished line is counted
when it completes, or when its build ends. This is a count of log diagnostic
lines, not unique issues or a full code analysis.

With eSpeak NG, the dispatcher announces the exact warning and build-target
counts. The voice is selected from `DISPATCHER_VOICE` (default
`en-us+f3,en-us+f4,en-us+m3`; set one voice or a comma-separated list via the
launcher or `--dispatcher-voice`). At zero warnings, randomized reports mention
clear tracks, weather, snow, avalanches or signals; above zero the dispatcher
calls caution on the line.
The wording is railway-themed fiction, not a real dispatcher instruction.
At the end of the build, success or failure, TAgent prints the recognized
compiler/linker `error:` (including `fatal error:`) and `warning:` line totals
across both build logs. A failed build remains failed even if its failure produced
no recognized diagnostic lines. In steampunk mode, the dispatcher announces the
same totals. With eSpeak NG they are spoken naturally; without TTS, prerecorded
fragments announce exact numbers digit by digit. Unrecognized formats are not
included. Build failures still propagate their original exit status.
When TTS is unavailable or disabled, bundled female-voice WAVs report either no
warnings or a generic caution; the fallback never invents a warning count. All
voice playback remains optional and follows the normal no-audio fallback. The
prerecorded progress fallback gives generic clear/caution audio rather than
inventing dynamic counts; the full counts remain visible on the console.

Preview either case without building (the warning preview uses a simulated 3):

```sh
bash tests/tagent/run_tagent.sh --sound-theme steampunk --sound-preview build-clear
bash tests/tagent/run_tagent.sh --sound-theme steampunk --sound-preview build-warning
bash tests/tagent/run_tagent.sh --sound-theme steampunk --sound-preview build-final-ok
bash tests/tagent/run_tagent.sh --sound-theme steampunk --sound-preview build-final-failed
```

## Optional local text-to-speech

`run_tagent.sh` checks for `espeak-ng` at startup. `TTS_MODE=auto` (default)
uses it only for the steampunk dispatch voice; `TTS_MODE=off` always uses the
bundled WAV. `TTS_VOICE=en-us+m3` selects the eSpeak voice. Both settings can be
edited at the top of the launcher, set in the environment, or overridden by
`--tts auto|off` and `--tts-voice NAME`. Synthesis runs in the audio thread, has
a 15-second timeout, and needs no network. Missing/failing TTS is reported and
falls back to the WAV; missing/failing audio output is reported and tests continue
without further playback attempts. Restart TAgent after connecting an audio device.
Demos and previews still generate WAV files when no output device is available.
`--no-sound` disables both synthesized and prerecorded sound.

On Ubuntu, install `espeak-ng` and an audio player such as `pulseaudio-utils`,
then try a simulated six-worker announcement without starting tests:

```sh
bash tests/tagent/run_tagent.sh --sound-theme steampunk --jobs 6 --sound-preview dispatch
bash tests/tagent/run_tagent.sh --sound-theme steampunk --tts off --sound-preview dispatch
```

The preview count is the configured worker limit; during a test cycle it is the
number of perf processes or stress groups last observed by the runner. Like the
other periodic steampunk effect, live dispatch is selected by a coin toss every
five minutes. Neither missing eSpeak nor an unavailable sound device blocks tests.

The demo also exports `demo-steampunk.wav` or `demo-dieselpunk.wav` in the state
directory, with short gaps between events. This allows listening later when the
current host has no audio device. Demo files include SOS at the end.
