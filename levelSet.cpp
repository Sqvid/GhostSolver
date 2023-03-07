#include <cmath>

#include "levelSet.hpp"

double circleLS(double x, double y, double t) {
	auto x0 {0.6}, y0 {0.6}, r {0.2}, vX {0.0};

	// Centre displacement.
	auto dX {vX*t};

	return r - std::sqrt((x - x0 - dX)*(x - x0 - dX) + (y - y0)*(y - y0));
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
