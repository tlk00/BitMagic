# strsvsample02a

This example is contributed by Andrei Shkeda (NCBI).


Example shows how to use `std::sort()` on BitMagic succinct vector for applications.

This case is useful where raw data does not fit in RAM and otherwise you would have to use disk sort. Operations in compressed memory space makes it possible without touching I/O.

Notable is that example uses succinct vector specific fast comparator, caching the last value, which is often reused in many sort implementations.

Sorted vector later remapped to use optimal codes, example shows how to collect runtime memory consumption of succinct containers.

This example demonstrates how sorted vector uses less memory due to compressive memory techniques of BitMagic Library. 


Expected output:

> ./strsvsample02a 
> 
> 1.std::sort() of index: ; 1.323 sec
> 
> Sort validation Ok.
> 
> Memory footprint statistics:
> 
> SV sorted(before remap) memory_used : 828740
> 
> SV unsorted memory_used             : 1011268
> 
> SV sorted(after remap) memory_used  : 344300
