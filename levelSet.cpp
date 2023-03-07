#include <cmath>

#include "levelSet.hpp"

double circleLS(double x, double y, double t) {
	auto x0 {0.6}, y0 {0.5}, r {0.2}, vX {0.0};

	// Centre displacement.
	auto dX {vX*t};

	return r - std::sqrt((x - x0 - dX)*(x - x0 - dX) + (y - y0)*(y - y0));
}

double separateCirclesLS(double x, double y, double t) {
	auto x1 {0.6}, y1 {0.25}, r1 {0.2};
	auto x2 {0.6}, y2 {0.75}, r2 {0.2};

	auto dX {0*t};

	// Centre displacement.

	auto phi1 =  r1 - std::sqrt((x - x1 - dX)*(x - x1 - dX) + (y - y1)*(y - y1));
	auto phi2 =  r2 - std::sqrt((x - x2 - dX)*(x - x2 - dX) + (y - y2)*(y - y2));

	// Inside circle 1.
	if (phi1 >= 0) {
		return phi1;

	// Inside circle 2.
	} else if (phi2 >= 0) {
		return phi2;
	}

	// Not inside either circle, return distance to closest surface.
	return phi1 > phi2 ? phi1 : phi2;
}

TwoVector findNormal(std::function<double (double, double, double)> levelSet,
		double x, double y, double t) {

	TwoVector normal;
	double dx {1e-6};
	double dy = dx;

	normal.x = (levelSet(x + dx, y, t) - levelSet(x - dx, y, t)) / (2 * dx);
	normal.y = (levelSet(x, y + dy, t) - levelSet(x, y - dy, t)) / (2 * dy);
	auto mag = normal.mag();

	// Normalise vector.
	normal.x /= mag;
	normal.y /= mag;

	return normal;
}
