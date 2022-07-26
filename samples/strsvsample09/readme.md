# strsvsample09

This example shows how to implement quick sort algorithm in compressive memory.

Example generates a simple synthetic random data, loads it into succinct sparse vector, benchmarks memory consumption and shows a few ways to implement quick sort and insertion sort.

# Quick Findings

1. Compressive sort uses signifiacntly less memory. In a real application memory footprint will be even worse due to heap fragmentation. 
2. Quick sort wins with a huge margine for compressive search over the insertion sort.
3. Sorting in non-compressive form (std::sort() ) is faster but it can be potencially compensated using parallel sort (not covered in this example).



## expected output

	Loading 950573 elements...
	std::vector<string> mem    = 30545041
	Succinct vector vector mem = 4261404
	Quick Sort... 1
	Quick Sort... 2
	Insertion Sort... 
	std::sort...
	Sort validation Ok.
	1. quick sort (succint); 23.86 sec
	2. quick sort 2 (succint); 13.33 sec
	3. insertion sort (succint); 4.845 min
	4. std::sort(); 523.9 ms


