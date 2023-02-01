#!/bin/env gnuplot

reset session
set term gif size 800,800 animate delay 1

set output "test1.gif"

stats "test1.dat" u 1 nooutput
N = STATS_blocks

do for [i=1:N] {
	plot "test1.dat" index i-1 u 1:2 with lines,\
		"" index i-1 u 1:3 with lines,\
		"" index i-1 u 1:4 with lines
}
