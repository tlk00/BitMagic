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
- Default coverage runs all produced perf release and debug variants, and every
  debug/release stress variant of `tests/stress/t.cpp`. Exclude stress64. Perf
  debug variants always use `--fast`; `--pedantic` builds additional x86 perf
  debug SIMD variants and enables release stress serialization sweeps.
- Run perf first, with a shared worker limit; timing comparisons are not an
  acceptance criterion. Then run every stress build using the existing parallel
  runner, one build at a time. All scheduled tests must succeed. Stress variants
  share learned heavy-first timings across cycles in `TAgent-logs/stress-order.tsv`
  (or `stress-order.tsv` under a custom state directory); individual test logs
  remain in their cycle directories.
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
  testing. Events: startup beep, build-success cha-cha, passing freight train,
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

## Steampunk sounds

During steampunk build warm-up, a mechanical spark telegraph sends
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
chord. The sounds include a build-up, build-success cha-cha, descending idle sigh,
victory horn and Morse SOS repeated three times. All effects and the combined
demo are 44.1 kHz mono PCM and need no extra Python packages. The dispatch
message is a bundled mono voice/steam/wheel recording.

```sh
# Demo only: no Git checks, fetching, builds or tests. Ctrl-C stops playback.
bash tests/tagent/run_tagent.sh --sound-demo --sound-theme steampunk
# All effects plus the full depot and five-minute train backgrounds.
bash tests/tagent/run_tagent.sh --sound-preview

# One sound, or normal operation with a chosen theme and volume (0..1).
bash tests/tagent/run_tagent.sh --sound-preview depart --sound-theme steampunk
bash tests/tagent/run_tagent.sh --sound-preview depot --sound-theme steampunk
bash tests/tagent/run_tagent.sh --sound-preview train-ambient --sound-theme steampunk
bash tests/tagent/run_tagent.sh --sound-preview depart-perf --sound-theme steampunk
bash tests/tagent/run_tagent.sh --sound-preview depart-stress --sound-theme steampunk
bash tests/tagent/run_tagent.sh --sound-preview dispatch --sound-theme steampunk
bash tests/tagent/run_tagent.sh --sound-volume 0.25
```

`SOUND_THEME` and `SOUND_VOLUME` are editable launcher defaults with environment
and CLI overrides. `--no-sound` disables all playback. Audio targets the test
machine, trying afplay, PipeWire, PulseAudio and ALSA, then the local console bell.
Unavailable audio does not prevent testing. Generated WAVs remain in the state
directory for inspection.

During normal steampunk operation, a quiet depot ambience loops in the background.
When perf or stress tests start, it transitions to a five-minute mono moving-train
loop with recorded wheel sounds and occasional quiet tonal Morse for BITMAGIC.
At the end of the tests, the background returns to the depot.
It fades to 10% over about two seconds before each sound event, then returns
over about three seconds. On Linux, both backgrounds use a continuous PCM stream
to one player, including while foreground announcements play. Linux tries PipeWire, PulseAudio and ALSA in
that order. On macOS, playback uses one-second WAV chunks with short tapers at
the boundaries. The background is independent of builds and tests;
an unavailable player disables the ambience without stopping work. Set
`DEPOT_SOUND=no` in the launcher or use `--no-depot` to disable both backgrounds;
`--no-sound` disables all audio. The normal sound demo and individual previews
play their effects without the continuous depot layer.

The depot recording includes modified excerpts from [Sherman tank, rijden en
manoeuvreren](https://commons.wikimedia.org/wiki/File:Sherman_tank,_rijden_en_manoeuvreren_-_SoundCloud_-_Beeld_en_Geluid.ogg),
by Beeld en Geluid, licensed under [CC BY-SA 3.0](https://creativecommons.org/licenses/by-sa/3.0/).
The continuous excerpt was resampled, normalized and mixed with quiet workshop sounds;
the resulting `assets/depot.wav` is distributed under the same license. To
reproduce it, convert the source to 22050 Hz, 16-bit mono WAV and run
`python3 tests/tagent/generate_depot_ambient.py SHERMAN.wav`.
The moving-train loop combines and modifies [Стук колёс поезда](https://commons.wikimedia.org/wiki/File:%D0%A1%D1%82%D1%83%D0%BA_%D0%BA%D0%BE%D0%BB%D1%91%D1%81_%D0%BF%D0%BE%D0%B5%D0%B7%D0%B4%D0%B0.ogg)
by Ural-66 and [Passenger Train on the run](https://freesound.org/people/ceejay/sounds/76078/)
by ceejay, plus [Railway Line & Signal Lights](https://freesound.org/people/dbache/sounds/351425/)
by dbache, all [CC0](https://creativecommons.org/publicdomain/zero/1.0/).
The Freesound sources are their publicly available compressed previews.
The sources were mixed to mono, resampled, normalized and crossfaded, with quiet
Morse added. To reproduce the mix, convert both source recordings to 22050 Hz,
16-bit mono WAV and run the mixer, which uses only Python's standard library:

```sh
python3 tests/tagent/generate_train_ambient.py WHEELS.wav BETWEEN_WAGONS.wav FREIGHT.wav
```

In steampunk mode, the perf and stress phases each get a prerecorded female
station announcement, with quiet radio static, followed by 25 seconds of Vasily
Agapkin's march "Farewell of Slavianka" through a station-speaker filter and the
existing 23-second departure: platform bell, steam whistle, then accelerating
chuffs and bogie knocks as the train pulls away. Each sequence lasts about 54
seconds and needs no TTS engine. The march fades over its last three seconds;
the departure fades over its last five, and a passing train fades over its last
two seconds. Testing starts without waiting for playback;
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
reset that clock.
Steampunk clips last 10.5 seconds: a heavy freight train approaches, passes the listener, and recedes. Prominent layered cylinder chuffs and steam breaths lead the slower paired bogie knocks
and low rolling rumble of the carriages. Their cadence is separate from the
exhaust; mild pitch/cadence change and distance envelopes suggest the passing
train. This is a procedural impression, not a recording or physical simulation.
Clips
play in background threads and are queued to avoid overlapping sounds. On failure,
queued non-emergency sounds are discarded, and the distress signal follows the
clip currently playing. The demo includes both low and full
worker occupancy; use `--jobs 8` if demonstrating on a single-worker machine.

The steampunk demo includes both phase departures, two passing trains, the dispatch
message, remote update/no-update sounds, clear/caution build reports, and
successful/failed build summaries. Each successful fetch or remote check prints
the result. During tests a remote check detects commits without fetching or
changing the checkout until a safe boundary. Optional TTS announces new updates
or a clear line; without TTS, different short sound effects indicate each case.
After fetched changes are applied and a build is due, the dispatcher announces
that the freight train is diverted to the build depot. A check during active
tests only announces that updates were detected.
The random choice applies to periodic test sounds, not to demo previews.

```sh
bash tests/tagent/run_tagent.sh --sound-preview repo-new
bash tests/tagent/run_tagent.sh --sound-preview repo-build
bash tests/tagent/run_tagent.sh --sound-preview repo-unchanged
```

While perf variants run, every ten minutes TAgent prints the progress of every
active variant and announces one variant in rotation using the separate Piper
conductor voice or a female eSpeak fallback. Progress counts complete timing-result lines
ending in `ms`, `sec` or `min`; blank lines, headers, informational lines and
incomplete lines are excluded. These are observed benchmark results, not a
prediction of the total number of tests. Without TTS, a short neutral sound
plays and the exact count remains on the console. Preview a simulated report:

```sh
bash tests/tagent/run_tagent.sh --sound-preview perf-progress
```


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

With Piper or eSpeak NG, the dispatcher announces the exact warning and build-target
counts. For eSpeak, the voice is selected from `DISPATCHER_VOICE` (default
`en-us+f3,en-us+f4,en-us+m3`; set one voice or a comma-separated list via the
launcher or `--dispatcher-voice`). At zero warnings, randomized reports mention
clear tracks, weather, snow, avalanches or signals; above zero the dispatcher
calls caution on the line.
The wording is railway-themed fiction, not a real dispatcher instruction.
At the end of the build, success or failure, TAgent prints the recognized
compiler/linker `error:` (including `fatal error:`) and `warning:` line totals
across both build logs. A failed build remains failed even if its failure produced
no recognized diagnostic lines. In steampunk mode, the dispatcher announces the
same totals. With TTS they are spoken dynamically; without TTS, prerecorded
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

`TTS_MODE=auto` (default) prefers local Piper, falls back to `espeak-ng`, then
to prerecorded audio or non-verbal effects if neither works. `TTS_MODE=off`
disables both engines. By default Piper uses `~/piper-demo/.venv/bin/python`
and keeps the dispatcher voice under `LOG_DIR/piper-voices/en_US-lessac-medium.onnx`
with its `.onnx.json` alongside. The conductor and engineer models are also
stored in `LOG_DIR/piper-voices`; configure them with `PIPER_CONDUCTOR_MODEL` /
`--piper-conductor-model` and `PIPER_ENGINEER_MODEL` / `--piper-engineer-model`.
`LOG_DIR` means the selected `--log-dir` or the default `TAgent-logs`.
The conductor voice is prepared at startup, in the background during normal
operation. TAgent reports when the files are present or downloaded, and logs
which engine actually speaks the first conductor report (and any switch to the
fallback). A sound preview waits for voice preparation before playing.
Change `PIPER_PYTHON` and `PIPER_MODEL` in the launcher or environment, or use
`--piper-python` and `--piper-model`, to select another installation or voice.
When needed, TAgent creates the default user-local virtual environment and
installs `piper-tts` in the background; it does not use sudo or modify system
Python. An existing configured interpreter is only modified if it belongs to a
virtual environment. A missing named voice is downloaded into staging and
validated by generating a short WAV before installation. TAgent uses eSpeak or
the existing sound fallback until Piper is ready. Setup is bounded and failures
do not stop tests. A standalone sound preview waits for the setup attempt.
Set `PIPER_INSTALL=no` or pass `--no-piper-install` to disable engine installation;
set `PIPER_DOWNLOAD=no` or pass `--no-piper-download` to disable voice downloads.
`TTS_VOICE=en-us+m3` selects the
eSpeak voice. These settings can be
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

The demo also exports `demo-steampunk.wav` in the state
directory, with short gaps between events. This allows listening later when the
current host has no audio device. Demo files include SOS at the end.

### Engineer Joe

During perf and stress testing, Joe offers a short original railway, office, or
BitMagic joke every randomly selected 10–20 minutes. His separate Piper voice is
`en_US-danny-low`, automatically provisioned under `LOG_DIR/piper-voices`. Missing
Piper falls back to the English male eSpeak voice; missing audio does not stop tests.
The usual ambient ducking applies. `--no-joker` disables the role; `--no-sound`
and `--tts off` also silence it.

The editable collection in `joker.txt` contains 240 original lines, including mild
adult double entendres, obscure C++ humor, standards committee jokes, and Git commits. All lines are shuffled without replacement across cycles
within one agent process: at least 40 hours before the deck repeats. Restarting
TAgent reshuffles it. Joe identifies himself as train engineer at each report.
Each test phase starts a fresh 10–20 minute timer; jokes wait until audio is idle.
Preview with `./run_tagent.sh --sound-preview joker`.

Perf progress includes an approximate percentage of timed test elements, based on
189 entries in a complete v10.0.1 log. Variants and future test changes can alter
the total; this is not an estimate of elapsed time. Active runs are capped at 99%.
If a run exceeds the baseline, the conductor announces a delay and its actual
completed count instead of a percentage.

Spoken update reports identify dispatch; a TAgent self-update announces that the
launcher is restarting before it restarts. The prerecorded departure announcements
remain station announcements.

On a parallel stress failure, TAgent prints the variant, failed group and exit
status, diagnostic log path, its last 12 lines, and quoted copy commands. A compact
`TAgent-logs/last-test-failure.txt` preserves the report; recovery exhaustion points
to it so the original failure remains easy to find.

### Log storage on another drive

Use `--log-dir PATH` (also accepted as `--state-dir`) or set `LOG_DIR` at the
top of `run_tagent.sh`. The command-line option overrides the environment.
For Windows drive D under WSL, use the Linux mount path:

```sh
bash tests/tagent/run_tagent.sh --log-dir /mnt/d/TAgent-logs
```

TAgent creates missing directories and checks actual file creation, writing,
removal, and available space before starting work. `--min-free-gb` sets the reserve
(default 5 GB); periodic checks continue during operation. Startup prints the
resolved directory. Logs, test working files, saved success state and learned
stress timing history all use this location. Existing logs/state are not moved;
stop TAgent and copy them yourself if needed. The option survives self-restarts.

Stress release variants run before all debug variants. Debug variants use the
runner `--fast` mode, including during `--pedantic` TAgent runs; release variants
use `--pedantic` when requested. The runner cannot combine these two modes.

After a played Joe announcement, Dispatch may reply after a random 3–7 second
pause, choosing from five responses. A failure, shutdown, or another queued
announcement cancels the response. Preview the reply with `--sound-preview
joker-reply`; `--sound-preview joker` plays the exchange when local audio is available.

When the SOS signal plays, available TTS follows it with an engineer mayday and
a dispatcher confirmation that the emergency response team is on the way. Without
TTS the SOS remains available. Preview the sequence with `--sound-preview fail`.

SOS and both distress voices use 1.25 times the configured sound volume, capped
at 1.0. The radio exchange has short static bursts at transmission boundaries and
a brief pause before each call.

Synthesized crew radio messages and prerecorded build-report fallbacks use short
static bursts. Normal messages retain their configured volume; distress audio uses
1.25 times that volume, capped at 1.0.

`--restart-on-update` interrupts active build or test processes after a remote
check detects a new commit, then fetches, fast-forwards and restarts the launcher
for a fresh build and test cycle. Without it, active tests finish at their usual
safe boundary. In-flight logs are retained. Checks during active work occur at
subprocess heartbeats; Git/network recovery rules still apply.

TAgent also recompiles the sample programs selected by `samples/Makefile` during
the build phase. `xsample04b` is excluded. Samples are built but not executed.
The shared CPU budget reserves up to two workers for perf and stress builds and
uses the remaining workers for parallel sample compilation (at least one when
samples run). Output is saved in `build-samples.log`; sample warnings are included
in the combined diagnostics and live spoken reports. Any sample build failure
stops the cycle before tests.

`--fast` runs every perf and stress variant in fast mode. Perf uses one timed
repetition and smaller large fixtures; stress uses reduced checks. `--fast` and `--pedantic` are mutually exclusive; TAgent rejects
the combination at startup.
Fast runs have distinct saved success state so `--startup changed` does not
confuse them with normal runs.

At the start of each testing phase and at random 10–20 minute intervals while
tests are running, the conductor announces the next destination:
Performance Test Station or Stress Test Station. With `--fast`, each
announcement describes express service and the reduced workload for its phase.
Preview with `--sound-preview station-perf` or `--sound-preview station-stress`
(add `--fast` for the express announcement). Normal stress announces the regular
schedule; `--pedantic` announces extended inspections and potential delays. Reports
wait for foreground audio to finish and can be delayed by subprocess heartbeats.
When TTS is unavailable the phase
departure sound still plays.

Foreground audio is serialized. Repeated routine status messages coalesce while
queued; phase-specific messages are discarded if their test phase has ended,
and pending build-status reports are discarded when the build finishes. Repeated
failure signals do not interrupt an already queued or active SOS and distress
exchange. Before normal shutdown or launcher restart, TAgent closes the audio
queue to new messages and drains it; Ctrl-C still stops child processes promptly.

During foreground announcements the background ambience fades down to 10% of
its normal level, then fades back up. Linux keeps its continuous audio stream.

Every produced, supported perf debug binary runs with `--fast`, including in
normal and `--pedantic` cycles. Release perf uses `--fast` only when the TAgent
`--fast` option is supplied.

At startup TAgent logs its PID, parent PID, process group, session and terminal.
SIGINT, SIGTERM and SIGHUP are reported by name and number with current process
context and cycle before child cleanup. Python's standard signal handler does not
provide the sender identity; these fields do not identify who sent a signal.

On a catchable interrupt, TAgent stops child jobs first, then, if TTS and audio
are available, the engineer reports that emergency brakes are applied and
Dispatch acknowledges. This announces a stop; a signal does not itself restart
TAgent. Preview with `--sound-preview interrupt`. Without audio or TTS, shutdown
continues without the exchange.

The engineer's Piper voice is prepared even with `--no-joker`, because emergency
announcements also use it. The interrupt preview waits for voice setup. The first
engineer utterance logs the engine and model actually used, or the eSpeak fallback.

TAgent stores its generated audio, cycle logs, scratch files, and default Piper
voice models under the configured log directory. Child processes receive a
`TMPDIR`, `TMP`, and `TEMP` pointing to its `tmp` subdirectory. Voice models
are retained when old cycle logs are removed. Explicit Piper model paths remain
supported and may point elsewhere. The Piper virtual environment itself remains
user-local by default so WSL log drives need not support executable venvs. Existing
voices are not moved automatically; TAgent downloads missing models into the new
location, or you can copy the `.onnx` and `.onnx.json` files there.
