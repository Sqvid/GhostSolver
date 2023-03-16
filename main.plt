#!/bin/env gnuplot

set term png size 800,800

# Sod's test in X.
filename = "test1X.dat"
stats filename nooutput
last = STATS_blocks - 1

set key on
set view 60,30

set output "test1XDensity.png"
splot filename i last u 1:2:3 w l

set output "test1XVelocityX.png"
splot filename i last u 1:2:4 w l

set output "test1XVelocityY.png"
splot filename i last u 1:2:5 w l

set output "test1XPressure.png"
splot filename i last u 1:2:6 w l

# Sod's test in Y.
filename = "test1Y.dat"
stats filename nooutput
last = STATS_blocks - 1

set key on
set view 60,30

set output "test1YDensity.png"
splot filename i last u 1:2:3 w l

set output "test1YVelocityX.png"
splot filename i last u 1:2:4 w l

set output "test1YVelocityY.png"
splot filename i last u 1:2:5 w l

set output "test1YPressure.png"
splot filename i last u 1:2:6 w l

# Cylindrical explosion test.
filename = "cylExpl.dat"
stats filename nooutput
last = STATS_blocks - 1

set key on
set view 60,30

set output "cylExplDensity.png"
splot filename i last u 1:2:3 w l

set output "cylExplVelX.png"
splot filename i last u 1:2:4 w l

set output "cylExplVelY.png"
splot filename i last u 1:2:5 w l

set output "cylExplPressure.png"
splot filename i last u 1:2:6 w l

# Rigid tests.
# Single circle.
filename = "single.dat"
stats filename nooutput
last = STATS_blocks - 1

set key on
set view 60,30

set output "singleDensity.png"
splot filename i last u 1:2:3 w l

set output "singleVelocityX.png"
splot filename i last u 1:2:4 w l

set output "singleVelocityY.png"
splot filename i last u 1:2:5 w l

set output "singlePressure.png"
splot filename i last u 1:2:6 w l

set output "singleSchlieren.png"
set view map
set palette gray
splot filename i last u 1:2:7 w pm3d


# Separated circles.
filename = "separate.dat"
stats filename nooutput
last = STATS_blocks - 1

set key on
set view 60,30

set output "separate.png"
splot filename i last u 1:2:3 w l

set output "sepVelocityX.png"
splot filename i last u 1:2:4 w l

set output "sepVelocityY.png"
splot filename i last u 1:2:5 w l

set output "sepPressure.png"
splot filename i last u 1:2:6 w l

set output "sepSchlieren.png"
set view map
set palette gray
splot filename i last u 1:2:7 w pm3d

# Overlapping cirlces.
filename = "overlap.dat"
stats filename nooutput
last = STATS_blocks - 1

set key on
set view 60,30

set output "overlap.png"
splot filename i last u 1:2:3 w l

set output "overlapVelocityX.png"
splot filename i last u 1:2:4 w l

set output "overlapVelocityY.png"
splot filename i last u 1:2:5 w l

set output "overlapPressure.png"
splot filename i last u 1:2:6 w l

set output "overlapSchlieren.png"
set view map
set palette gray
splot filename i last u 1:2:7 w pm3d
