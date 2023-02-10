#!/bin/env gnuplot
reset session

numFiles=5

file(ft, n) = sprintf("test%d.%s", n, ft)

do for [i=1:numFiles] {
	set term gif size 800,800 animate delay 10
	set output file("gif", i)

	stats file("dat", i) u 1 nooutput
	nBlocks = STATS_blocks

	do for [j=1:nBlocks] {
		set multiplot layout 3,1

		set title sprintf("Dataset %d: density", i)
		plot file("dat", i) index j-1 u 1:2 with lines linetype rgb "blue"

		set title sprintf("Dataset %d: velocity", i)
		plot "" index j-1 u 1:3 with lines linetype rgb "red"

		set title sprintf("Dataset %d: pressure", i)
		plot "" index j-1 u 1:4 with lines linetype rgb "violet"

		unset multiplot
	}
}
