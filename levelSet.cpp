#include <cmath>

#include "levelSet.hpp"

double circleLS(double x, double y) {
	auto x0 {0.6}, y0 {0.6}, r {0.2};

	return r - std::sqrt((x - x0)*(x - x0) + (y - y0)*(y - y0));
}

TwoVector findNormal(std::function<double (double, double)> lsFunc,
		double x, double y) {

	TwoVector normal;
	double dx {1e-6};
	double dy = dx;

	normal.x = (lsFunc(x + dx, y) - lsFunc(x - dx, y)) / (2 * dx);
	normal.y = (lsFunc(x, y + dy) - lsFunc(x, y - dy)) / (2 * dy);
	auto mag = normal.mag();

	// Normalise vector.
	normal.x /= mag;
	normal.y /= mag;

	return normal;
}
