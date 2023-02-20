#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <stdexcept>

#include "eulerData.hpp"
#include "simulation.hpp"
#include "slopeLimiter.hpp"

namespace fvm::internal {

	double xiRight(double r) {
		return 2 / (1 + r);
	}

	double minbee(double r) {
		double limited;

		if (r <= 0) {
			limited = 0;

		} else if (r <= 1) {
			limited = r;

		} else {
			limited = std::min(1.0, xiRight(r));
		}

		return limited;
	}

	double superbee(double r) {
		double limited;

		if (r <= 0) {
			limited = 0;

		} else if (r <= 0.5) {
			limited = 2*r;

		} else if (r <= 1) {
			limited = 1;

		} else {
			limited = std::min({r, xiRight(r), 2.0},
					[](const double& a, const double& b) {
						return a < b;
					});
		}

		return limited;
	}

	double vanAlbada(double r) {
		double limited;

		if (r <= 0) {
			limited = 0;

		} else {
			limited = std::min((r*(1 + r))/(1 + r*r), xiRight(r));
		}

		return limited;
	}

	double vanLeer(double r) {
		double limited;

		if (r <= 0) {
			limited = 0;

		} else {
			double xiR = xiRight(r);
			limited = std::min(xiR * r, xiR);
		}

		return limited;
	}

	double limit(double r, SlopeLimiter slType) {
		double limited = r;

		switch (slType) {
			case SlopeLimiter::minbee:
				limited = minbee(r);
				break;
			case SlopeLimiter::superbee:
				limited = superbee(r);
				break;
			case SlopeLimiter::vanAlbada:
				limited = vanAlbada(r);
				break;
			case SlopeLimiter::vanLeer:
				limited = vanLeer(r);
				break;
			default:
				break;
		}

		return limited;
	}
}
