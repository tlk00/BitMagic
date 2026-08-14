# bvsample23: Interval enumeration

This example interprets runs of set bits as closed intervals and traverses
them with `bm::interval_enumerator<>` from `bmintervals.h`.

It covers the direct `valid()`/`advance()` interface, STL-style iteration,
access through both `start()`/`end()` and pair values, and construction or
repositioning at an arbitrary bit index.

The `extend_left` option controls whether positioning inside an interval
returns that interval's true left boundary or starts at the requested
position. Positioning in a gap advances to the next interval, while moving
beyond the final interval invalidates the enumerator.
