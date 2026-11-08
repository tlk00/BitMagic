#!/bin/bash
# TAgent launcher. Environment and CLI overrides survive self-restarts.
POLL_MINUTES=${POLL_MINUTES:-5}
RETRY_MINUTES=${RETRY_MINUTES:-60}
RESERVE_CPUS=${RESERVE_CPUS:-2}
JOBS=${JOBS:-0}                  # 0 = automatic
MIN_FREE_GB=${MIN_FREE_GB:-5}
REMOTE=${REMOTE:-origin}
BRANCH=${BRANCH:-master}
SOUND=${SOUND:-yes}
SOUND_THEME=${SOUND_THEME:-dieselpunk}  # steampunk | dieselpunk
SOUND_VOLUME=${SOUND_VOLUME:-0.35}      # 0.0 .. 1.0
PYTHON=${PYTHON:-python3}

set -eu
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
command -v "$PYTHON" >/dev/null || { echo 'Python 3.9+ is required' >&2; exit 1; }
"$PYTHON" -c 'import sys; sys.exit(0 if sys.version_info >= (3,9) else 1)' || {
    echo 'Python 3.9+ is required' >&2; exit 1;
}
export SOUND_THEME SOUND_VOLUME
export POLL_MINUTES RETRY_MINUTES RESERVE_CPUS JOBS MIN_FREE_GB REMOTE BRANCH SOUND PYTHON
exec "$PYTHON" "$SCRIPT_DIR/tagent.py" "$@"
