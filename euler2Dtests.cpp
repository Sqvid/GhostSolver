#include "euler2Dtests.hpp"

double test1DensityX(double x, __attribute__((unused)) double y) {
	return x <= 0.5 ? 1 : 0.125;
}

double test1DensityY(__attribute__((unused)) double x, double y) {
	return y <= 0.5 ? 1 : 0.125;
}

double test1VelocityX(__attribute__((unused)) double x, __attribute__((unused)) double y) {
	return 0;
}

double test1VelocityY(__attribute__((unused)) double x, __attribute__((unused)) double y) {
	return 0;
}

double test1PressureX(double x, __attribute__((unused)) double y) {
	return x <= 0.5 ? 1 : 0.1;
}

double test1PressureY(__attribute__((unused)) double x, double y) {
	return y <= 0.5 ? 1 : 0.1;
}

double constantVX(__attribute__((unused)) double x, __attribute__((unused)) double y) {
	return 0;
}

double constantVY(__attribute__((unused)) double x, __attribute__((unused)) double y) {
	return 0;
}

double explLS(double x, double y) {
	double r = 0.4;

	return r*r - x*x - y*y;
}

double cylExplDensity(double x, double y) {
	double phi = explLS(x, y);

	return phi >= 0 ? 1 : 0.125;
}

double cylExplVelocityX(__attribute__((unused)) double x, __attribute__((unused)) double y) {
	return 0;
}

double cylExplVelocityY(__attribute__((unused)) double x, __attribute__((unused)) double y) {
	return 0;
}

double cylExplPressure(double x, double y) {
	double phi = explLS(x, y);

	return phi >= 0 ? 1 : 0.1;
}

double rigidTestDensity(double x, __attribute__((unused)) double y) {
	return x <= 0.2 ? 1.3764 : 1;
}

double rigidTestVelocityX(double x, __attribute__((unused)) double y) {
	return x <= 0.2 ? 0.394 : 0;
}

double rigidTestVelocityY(__attribute__((unused)) double x, __attribute__((unused)) double y) {
	return 0;
}

double rigidTestPressure(double x, __attribute__((unused)) double y) {
	return x <= 0.2 ? 1.5698 : 1;
}
