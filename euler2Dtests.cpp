#include "euler2Dtests.hpp"

double test1DensityX(double x, double y) {
	return x <= 0.5 ? 1 : 0.125;
}

double test1DensityY(double x, double y) {
	return y <= 0.5 ? 1 : 0.125;
}

double test1VelocityX(double x, double y) {
	x = 0;
	return x;
}

double test1VelocityY(double x, double y) {
	y = 0;
	return y;
}

double test1PressureX(double x, double y) {
	return x <= 0.5 ? 1 : 0.1;
}

double test1PressureY(double x, double y) {
	return y <= 0.5 ? 1 : 0.1;
}

double constantVX(double x, double y) {
	x = 0;

	return y;
}

double constantVY(double x, double y) {
	y = 0;

	return y;
}

double explLS(double x, double y) {
	double r = 0.4;

	return x*x + y*y - r*r;
}

double cylExplDensity(double x, double y) {
	double phi = explLS(x, y);

	return phi <= 0 ? 1 : 0.125;
}

double cylExplVelocityX(double x, double y) {
	x = 0;
	y = 0;

	return x + y;
}

double cylExplVelocityY(double x, double y) {
	x = 0;
	y = 0;

	return x + y;
}

double cylExplPressure(double x, double y) {
	double phi = explLS(x, y);

	return phi <= 0 ? 1 : 0.1;
}

double rigidTestDensity(double x, double y) {
	y = 0;

	return x + y <= 0.2 ? 1.3764 : 1;
}

double rigidTestVelocityX(double x, double y) {
	y = 0;

	return x + y <= 0.2 ? 0.394 : 0;
}

double rigidTestVelocityY(double x, double y) {
	x = 0;
	y = 0;

	return x + y;
}

double rigidTestPressure(double x, double y) {
	y = 0;

	return x + y <= 0.2 ? 1.5698 : 1;
}
