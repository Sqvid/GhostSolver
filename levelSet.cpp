#include <cmath>

#include "levelSet.hpp"

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

// Helper function to calculate distance between (x1, y1) and (x2, y2).
double distBetween(double x1, double y1, double x2, double y2) {
	return std::sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2));
}

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

double overlapCirclesLS(double x, double y, double t) {
	auto x0 {0.6}, y1 {0.35}, y2 {0.65}, r {0.2};

	auto dX {0*t};

	// Centre displacement.

	auto phi1 =  r - std::sqrt((x - x0 - dX)*(x - x0 - dX) + (y - y1)*(y - y1));
	auto phi2 =  r - std::sqrt((x - x0 - dX)*(x - x0 - dX) + (y - y2)*(y - y2));

	// Outside both circles.
	if (phi1 < 0 && phi2 < 0) {
		return phi1;

	// Inside circle 1.
	} else if (phi1 >= 0 && phi2 < 0) {
		return phi1;

	// Inside circle 2.
	} else if (phi1 < 0 && phi2 >= 0) {
		return phi2;
	}

	// Inside the intersection of both circles. Return closest intersection
	// point.
	auto yDiff = y2 - y1;
	auto halfChord = std::sqrt(r*r - yDiff * yDiff * 0.5 * 0.5);

	// (xI, yI) is the coordinate of the nearest intersection.
	auto xI {0.0};
	auto yI = y1 + 0.5*yDiff;
	double distToInterface;

	if (x < x0) {
		xI = x0 - halfChord;

	} else {
		xI = x0 + halfChord;
	}

	return distToInterface = std::sqrt((x - xI)*(x - xI) + (y - yI)*(y - yI));
}

double squareLS(double x, double y, double t) {
	// (xC, yC) is the centre of the square. hSide is half the side length.
	auto xC {0.6}, yC {0.5}, hSide {0.2};

	// Bounds of the square. L, R, U, D = left, right, up, down.
	auto xL = xC - hSide;
	auto xR = xC + hSide;
	auto yU = yC + hSide;
	auto yD = yC - hSide;

	// Inside the square.
	if (x >= xL && x <= xR && y >= yD && y <= yU) {
		auto lDist = x - xL;
		auto rDist = xR - x;
		// The closest x-bound.
		auto xDist = lDist <= rDist ? lDist : rDist;

		auto uDist = yU - y;
		auto dDist = y - yD;
		// The closest y-bound.
		auto yDist = uDist <= dDist ? uDist : dDist;

		// The closest side.
		return xDist <= yDist ? xDist : yDist;
	}


	// Outside the square. By convention the return value here is negative.
	// On the left.
	if (x < xL) {
		// Top-left corner.
		if (y > yU) {
			return -distBetween(x, y, xL, yU);

		// Bottom-left corner.
		} else if (y < yD) {
			return -distBetween(x, y, xL, yD);

		// Directly left of the square.
		} else {
			return x - xL;
		}

	// On the right.
	} else if (x > xR) {
		// Top-right corner.
		if (y > yU) {
			return -distBetween(x, y, xR, yU);

		// Bottom-right corner.
		} else if (y < yD) {
			return -distBetween(x, y, xR, yD);

		// Directly right of the square.
		} else {
			return xR - x;
		}

	// Directly above or below the square.
	} else {
		// Above.
		if (y > yU) {
			return yU - y;

		// Below
		} else {
			return y - yD;
		}
	}
}
