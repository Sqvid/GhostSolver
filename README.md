# A Riemann Ghost-Fluid Solver for the Euler Equations

This repository contains the code for a two-dimensional ghost-fluid solver for
the [Euler Equations](https://en.wikipedia.org/wiki/Euler_equations_(fluid_dynamics)).

The majority of the computation can be found within the implementation of the
Simulation class. This solver supports a variety of slope-limiters and has been
tested against several standard problems from the literature (See Riemann
Solvers and Numerical Methods for Fluid Dynamics: A Practical Introduction, by
E. F. Toro).
