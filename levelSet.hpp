#ifndef GHOSTSOLVER_LEVELSET_HPP
#define GHOSTSOLVER_LEVELSET_HPP

#include <functional>

#include "twoVector.hpp"

double circleLS(double x, double y);
TwoVector findNormal(std::function<double (double, double)> lsFunc,
		double x, double y, double dx , double dy);

#endif
