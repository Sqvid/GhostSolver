#!/bin/env gnuplot
reset session

numFiles=2

file(ft, n) = sprintf("test%d.%s", n, ft)

do for [i=1:numFiles] {
	set term gif size 800,800 animate delay 1
	set output file("gif", i)
	stats file("dat", i) u 1 nooutput
	nBlocks = STATS_blocks

	do for [j=1:nBlocks] {
		plot file("dat", i) index j-1 u 1:2 with lines,\
			"" index j-1 u 1:3 with lines,\
			"" index j-1 u 1:4 with lines
	}
}
