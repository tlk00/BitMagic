#!/bin/sh
#
# parallel exucion of tests in different terminal windows for MacOS

open -a Terminal "./stress64_debug"
open -a Terminal "./stress64_release"
open -a Terminal "./stress64_release_sse2"
open -a Terminal "./stress64_release_sse42"
open -a Terminal "./stress64_release_avx2"

