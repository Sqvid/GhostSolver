#ifndef GHOSTSOLVER_TORO2DTESTS_HPP
#define GHOSTSOLVER_TORO2DTESTS_HPP

// Initial distribution of primitive variables for Toro's tests.
//
//        rho_L   v_L     p_L     rho_R    v_R     p_R     time
// _____________________________________________________________
// 1 |    1       0       1       .125     0       .1      .25
// 2 |    1       -2      .4      1        2       .4      .15
// 3 |    1       0       1000    1        0       .01     .012
// 4 |    1       0       .01     1        0       100     .035
// 5 |    5.99924 19.5975 460.894 5.99242 -6.19633 46.0950 .035
//
// Each test is run on a x-domain [0, 1].
// L subscript for x <= 0.5, R subscript for x > 0.5
// gamma = 1.4 for all cases.
//
// Source for test initial conditions:
// Riemann Solvers and Numerical Methods for Fluid Dynamics, E.F. Toro,
// ISBN : 978-3-662-03492-7

double test1DensityX(double x, double y);
double test1DensityY(double x, double y);
double test1VelocityX(double x, double y);
double test1VelocityY(double x, double y);
double test1PressureX(double x, double y);
double test1PressureY(double x, double y);
double constantVY(double x, double y);
double constantVX(double x, double y);

// Cylindrical explosion test also from Toro.
double cylExplDensity(double x, double y);
double cylExplVelocityX(double x, double y);
double cylExplVelocityY(double x, double y);
double cylExplPressure(double x, double y);

// Rigid-body tests, initial conditions.
double rigidTestDensity(double x, double y);
double rigidTestVelocityX(double x, double y);
double rigidTestVelocityY(double x, double y);
double rigidTestPressure(double x, double y);

#endif
