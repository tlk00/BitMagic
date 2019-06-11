#!/bin/sh
#
# parallel execion of tests in different terminal windows for MacOS

open -a Terminal "./stress64_debug"
open -a Terminal "./stress64_debug_sse2"
open -a Terminal "./stress64_debug_sse42"
open -a Terminal "./stress64_debug_avx2"

