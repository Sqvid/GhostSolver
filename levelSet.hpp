#ifndef GHOSTSOLVER_LEVELSET_HPP
#define GHOSTSOLVER_LEVELSET_HPP

#include <functional>

#include "twoVector.hpp"

TwoVector findNormal(std::function<double (double, double, double)> levelSet,
		double x, double y, double t);

double circleLS(double x, double y, double t);
double squareLS(double x, double y, double t);
double separateCirclesLS(double x, double y, double t);
double overlapCirclesLS(double x, double y, double t);

#endif
