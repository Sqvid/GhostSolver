#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <stdexcept>

#include "eulerData.hpp"
#include "simulation.hpp"
#include "slopeLimiter.hpp"

namespace fvm {
	const double slopeTolerence_ = 0.00001;

	double xiRight(double r) {
		return 2 / (1 + r);
	}

	double minbee(double r) {
		if (r <= 0) {
			return 0;

		} else if (r > 0 && r <= 1) {
			return r;

		} else {
			return std::min(1.0, xiRight(r));
		}
	}

	double superbee(double r) {
		if (r <= 0) {
			return 0;

		} if (r > 0 && r <= 0.5) {
			return 2*r;

		} else if (r > 0.5 && r <= 1) {
			return 1;

		} else {
			return std::min({r, xiRight(r), 2.0},
					[](const double& a, const double& b) {
						return a < b;
					});
		}
	}

	double vanAlbada(double r) {
		if (r <= 0) {
			return 0;

		} else {
			return std::min((r*(1 + r))/(1 + r*r), xiRight(r));
		}
	}

	double vanLeer(double r) {
		if (r <= 0) {
			return 0;

		} else {
			auto xiR = xiRight(r);
			return std::min(xiR * r, xiR);
		}
	}

	double limit(double r, SlopeLimiter slType) {
		double result = r;

		switch (slType) {
			case SlopeLimiter::minbee:
				result = minbee(r);
				break;
			case SlopeLimiter::superbee:
				result = superbee(r);
				break;
			case SlopeLimiter::vanAlbada:
				result = vanAlbada(r);
				break;
			case SlopeLimiter::vanLeer:
				result = vanLeer(r);
				break;
			default:
				break;
		}

		return result;
	}

	void linearReconst(EulerData& eulerData, CellVector& lIfaces, CellVector& rIfaces, SlopeLimiter slType) {
		// Alias the numerical data.
		eulerData.setMode(EulerDataMode::conserved);
		CellVector& u = eulerData.data();

		// Index for energy. This is the quantity we are going to limit
		// on.
		int eIndex = static_cast<int>(ConservedQuant::energy);

		// Calculate reconstructed interface values.
		for (size_t i = 1; i < u.size() - 1; i++) {
			QuantArray deltaLeft = u[i] - u[i - 1];
			QuantArray deltaRight = u[i + 1] - u[i];

			QuantArray delta = 0.5 * (deltaLeft + deltaRight);

			double r;

			if (std::fabs(deltaRight[eIndex]) < slopeTolerence_) {
				r = 0;
			} else {
				r = deltaLeft[eIndex] / deltaRight[eIndex];
			}

			QuantArray uLeft = u[i] - 0.5 * delta * limit(r, slType);
			QuantArray uRight = u[i] + 0.5 * delta * limit(r, slType);

			lIfaces[i - 1] = uLeft;
			rIfaces[i - 1] = uRight;
		}
	}
}
