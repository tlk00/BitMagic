#!/bin/bash
# TAgent launcher. Environment and CLI overrides survive self-restarts.
STARTUP=${STARTUP:-build}        # build | changed
POLL_MINUTES=${POLL_MINUTES:-5}
RETRY_MINUTES=${RETRY_MINUTES:-60}
LOG_DIR=${LOG_DIR:-}            # empty = repo/TAgent-logs; WSL D: /mnt/d/TAgent-logs
LOG_RETENTION_DAYS=${LOG_RETENTION_DAYS:-7}
BUILD_REPORT_SECONDS=${BUILD_REPORT_SECONDS:-60}
RESERVE_CPUS=${RESERVE_CPUS:-2}
JOBS=${JOBS:-0}                  # 0 = automatic
MIN_FREE_GB=${MIN_FREE_GB:-5}
INV_READ_MEMORY_MIB=${INV_READ_MEMORY_MIB:-} # empty = auto (Linux: 70% available RAM; Mac: 8192)
INV_READ_MMAP=${INV_READ_MMAP:-yes}            # yes | no; Unix mapping windows
INV_CACHE_DIR=${INV_CACHE_DIR:-}  # empty = auto /mnt/c/TAgent-cache on Linux, if present
REMOTE=${REMOTE:-origin}
BRANCH=${BRANCH:-master}
SOUND=${SOUND:-yes}
DEPOT_SOUND=${DEPOT_SOUND:-yes}  # yes | no (steampunk ambience)
SOUND_THEME=${SOUND_THEME:-steampunk}
TTS_MODE=${TTS_MODE:-auto}      # auto | off
PIPER_PYTHON=${PIPER_PYTHON:-$HOME/piper-demo/.venv/bin/python}
PIPER_MODEL=${PIPER_MODEL:-}  # empty = LOG_DIR/piper-voices/en_US-lessac-medium.onnx
PIPER_CONDUCTOR_MODEL=${PIPER_CONDUCTOR_MODEL:-}  # empty = LOG_DIR/piper-voices
PIPER_ENGINEER_MODEL=${PIPER_ENGINEER_MODEL:-}    # empty = LOG_DIR/piper-voices
PIPER_FIREMAN_MODEL=${PIPER_FIREMAN_MODEL:-}      # empty = LOG_DIR/piper-voices/en_US-ryan-medium.onnx
PIPER_FOREMAN_MODEL=${PIPER_FOREMAN_MODEL:-}      # empty = LOG_DIR/piper-voices/en_AU-librivox-medium.onnx
PIPER_INSTALL=${PIPER_INSTALL:-yes}    # yes | no: install Piper in its user venv
PIPER_DOWNLOAD=${PIPER_DOWNLOAD:-yes}  # yes | no: download a missing voice
TTS_VOICE=${TTS_VOICE:-en-us+m3}
DISPATCHER_VOICE=${DISPATCHER_VOICE:-en-us+f3,en-us+f4,en-us+m3}  # comma-separated
SOUND_VOLUME=${SOUND_VOLUME:-0.35}      # 0.0 .. 1.0
PYTHON=${PYTHON:-python3}

set -eu
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
command -v "$PYTHON" >/dev/null || { echo 'Python 3.9+ is required' >&2; exit 1; }
"$PYTHON" -c 'import sys; sys.exit(0 if sys.version_info >= (3,9) else 1)' || {
    echo 'Python 3.9+ is required' >&2; exit 1;
}
# TTS is optional; preserve an empty result so Python knows it was checked.
TAGENT_TTS_ENGINE="$(command -v espeak-ng || true)"
TAGENT_TTS_CHECKED=1
export TAGENT_TTS_ENGINE TAGENT_TTS_CHECKED TTS_MODE TTS_VOICE DISPATCHER_VOICE BUILD_REPORT_SECONDS
export PIPER_PYTHON PIPER_MODEL PIPER_CONDUCTOR_MODEL PIPER_ENGINEER_MODEL PIPER_FIREMAN_MODEL PIPER_FOREMAN_MODEL PIPER_INSTALL PIPER_DOWNLOAD
export SOUND_THEME SOUND_VOLUME DEPOT_SOUND
export STARTUP LOG_RETENTION_DAYS LOG_DIR INV_READ_MEMORY_MIB INV_READ_MMAP INV_CACHE_DIR
export POLL_MINUTES RETRY_MINUTES RESERVE_CPUS JOBS MIN_FREE_GB REMOTE BRANCH SOUND PYTHON
exec "$PYTHON" "$SCRIPT_DIR/tagent.py" "$@"
