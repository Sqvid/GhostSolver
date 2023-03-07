#ifndef GHOSTSOLVER_LEVELSET_HPP
#define GHOSTSOLVER_LEVELSET_HPP

#include <functional>

#include "twoVector.hpp"

double circleLS(double x, double y, double t);
double separateCirclesLS(double x, double y, double t);

TwoVector findNormal(std::function<double (double, double, double)> levelSet,
		double x, double y, double t);

#endif
