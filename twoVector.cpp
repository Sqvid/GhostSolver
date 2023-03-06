#include "twoVector.hpp"

// Constructor.
TwoVector::TwoVector(double xVal, double yVal) {
	x = xVal;
	y = yVal;
}

double TwoVector::mag() {
	return std::sqrt(x*x + y*y);
}

// Operator overloads:
// Dot-product.
double TwoVector::operator*(const TwoVector& v) {
	return x * v.x + y * v.y;
}

// Multiplication by a scalar.
TwoVector TwoVector::operator*(const double& a) {
	TwoVector result;

	x = a * x;
	y = a * y;

	return result;
}

TwoVector operator*(const double& a, const TwoVector& v) {
	TwoVector result;

	result.x = a * v.x;
	result.y = a * v.y;

	return result;
}

// Vector addition/subtraction.
TwoVector TwoVector::operator+(const TwoVector& v) {
	TwoVector result;

	x = x + v.x;
	y = y + v.y;

	return result;
}

TwoVector TwoVector::operator-(const TwoVector& v) {
	TwoVector result;

	x = x - v.x;
	y = y - v.y;

	return result;
}
