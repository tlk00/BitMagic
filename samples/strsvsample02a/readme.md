# strsvsample02a

This example contributed by Andrei Shkeda (NCBI).


Example shows how to use std::sort on BitMagic succinct vector and generate sorted succinct vector.

Notable is that example uses succinct vector specific fast comparator, caching the last value, which is often reused in many sort implementations.

Sorted vector later remapped to use optimal codes, example shows how to collect runtime memory consumption of succinct containers.

This example demonstrates how sorted vector uses less memory due to compressive memory techniques of BitMagic Library. 
