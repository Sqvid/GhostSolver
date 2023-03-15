#!/bin/env gnuplot

# Sod's test in X
stats "test1X.dat" nooutput
last = STATS_blocks - 1

set term png size 800,800
set output "test1XDensity.png"
splot "test1X.dat" i last u 1:2:3 w l

set output "test1XVelocityX.png"
splot "test1X.dat" i last u 1:2:4 w l

set output "test1XVelocityY.png"
splot "test1X.dat" i last u 1:2:5 w l

set output "test1XPressure.png"
splot "test1X.dat" i last u 1:2:6 w l

# Sod's test in Y
stats "test1Y.dat" nooutput
last = STATS_blocks - 1

set output "test1YDensity.png"
splot "test1Y.dat" i last u 1:2:3 w l

set output "test1YVelocityX.png"
splot "test1Y.dat" i last u 1:2:4 w l

set output "test1YVelocityY.png"
splot "test1Y.dat" i last u 1:2:5 w l

set output "test1YPressure.png"
splot "test1Y.dat" i last u 1:2:6 w l

# Cylindrical explosion test
stats "cylExpl.dat" nooutput
last = STATS_blocks - 1

set output "cylExplDensity.png"
splot "cylExpl.dat" i last u 1:2:3 w l

set output "cylExplVelX.png"
splot "cylExpl.dat" i last u 1:2:4 w l

set output "cylExplVelY.png"
splot "cylExpl.dat" i last u 1:2:5 w l

set output "cylExplPressure.png"
splot "cylExpl.dat" i last u 1:2:6 w l

# Rigid tests
stats "single.dat" nooutput
last = STATS_blocks - 1

set output "single.png"
splot "single.dat" i last u 1:2:3 w l

stats "separate.dat" nooutput
last = STATS_blocks - 1

set output "separate.png"
splot "separate.dat" i last u 1:2:3 w l
set output "sepSchlieren.png"
set view map
set palette gray
splot "separate.dat" i last u 1:2:7 w pm3d

stats "overlap.dat" nooutput
last = STATS_blocks - 1

set output "overlap.png"
splot "overlap.dat" i last u 1:2:3 w l
