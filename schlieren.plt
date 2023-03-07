#!/bin/env gnuplot

stats "test1.dat" u 1 nooutput
last = STATS_blocks - 1

# Create GIF.
set term gif size 800,800 animate delay 10
set palette gray
set output "schlieren.gif"

set view map

do for [block=1:last] {
	splot "test1.dat" i block u 1:2:7 w pm3d
}

# Create PNG.
set term png size 800, 800
set output "schlieren.png"

set view map

# Plot final frame as a PNG.
splot "test1.dat" i last u 1:2:7 w pm3d

