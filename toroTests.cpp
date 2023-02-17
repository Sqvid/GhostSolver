#include "toroTests.hpp"

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
// Source for test initail conditions:
// Riemann Solvers and Numerical Methods for Fluid Dynamics, E.F. Toro,
// ISBN : 978-3-662-03492-7

double test1Density(double x) {
	return x <= 0.5 ? 1 : 0.125;
}

double test1Velocity(double x) {
	x = 0;
	return x;
}

double test1Pressure(double x) {
	return x <= 0.5 ? 1 : 0.1;
}

double test2Density(double x) {
	x = 1;
	return x;
}

double test2Velocity(double x) {
	return x <= 0.5 ? -2 : 2;
}

double test2Pressure(double x) {
	x = 0.4;
	return x;
}

double test3Density(double x) {
	x = 1;
	return x;
}

double test3Velocity(double x) {
	x = 0;
	return x;
}

double test3Pressure(double x) {
	return x <= 0.5 ? 1000 : 0.01;
}

double test4Density(double x) {
	x = 1;
	return x;
}

double test4Velocity(double x) {
	x = 0;
	return x;
}

double test4Pressure(double x) {
	return x <= 0.5 ? 0.01 : 100;
}

double test5Density(double x) {
	return x <= 0.5 ? 5.99924 : 5.99242;
}

double test5Velocity(double x) {
	return x <= 0.5 ? 19.5975 : -6.19633;
}

double test5Pressure(double x) {
	return x <= 0.5 ? 460.894 : 46.0950;
}

double constantVY(double y) {
	y = 0;

	return y;
}
