#ifndef GHOSTSOLVER_TWOVECTOR_HPP
#define GHOSTSOLVER_TWOVECTOR_HPP

#include <cmath>

// A simple two-dimensional linear algebra vector.
struct TwoVector {
	double x;
	double y;

	// Public member functions.
	// Constructor.
	TwoVector(double xVal = 0, double yVal = 0);

	double mag();

	// Operator overloads:
	// Dot-product.
	double operator*(const TwoVector& v);

	// Multiplication by a scalar.
	TwoVector operator*(const double& a);
	friend TwoVector operator*(const double& a, const TwoVector& v);

	// Vector addition/subtraction.
	TwoVector operator+(const TwoVector& v);
	TwoVector operator-(const TwoVector& v);
};

#endif
