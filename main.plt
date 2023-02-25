#!/bin/env gnuplot
reset session

numFiles = 5

file(ft, n) = sprintf("test%d.%s", n, ft)

do for [i=1:numFiles] {
	set term gif size 800,800 animate delay 10
	set output file("gif", i)

	stats file("dat", i) u 1 nooutput
	nBlocks = STATS_blocks

	do for [j=1:nBlocks] {
		set multiplot layout 2,2

		set title sprintf("Dataset %d: density", i)
		splot file("dat", i) index j-1 u 1:2:3 w l lt rgb "blue"

		set title sprintf("Dataset %d: x-velocity", i)
		splot "" index j-1 u 1:2:4 w l lt rgb "red"

		set title sprintf("Dataset %d: y-velocity", i)
		splot "" index j-1 u 1:2:5 w l lt rgb "green"

		set title sprintf("Dataset %d: pressure", i)
		splot "" index j-1 u 1:2:6 w l lt rgb "violet"

		unset multiplot
	}
}
