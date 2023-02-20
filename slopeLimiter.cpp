#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <stdexcept>

#include "eulerData.hpp"
#include "simulation.hpp"
#include "slopeLimiter.hpp"

namespace fvm {
	const double slopeTolerence = 0.00001;

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

	// TODO: Move this function into Simulation class.
	void internal::linearReconst(EulerData& eulerData, Grid& lIfaces, Grid& rIfaces, SlopeLimiter slType) {
		// Alias the numerical data.
		eulerData.setMode(EulerDataMode::conserved);
		Grid& u = eulerData.data();

		// Index for energy. This is the quantity we are going to limit
		// on.
		int eIndex = static_cast<int>(ConservedQuant::energy);

		// Calculate reconstructed interface values.
		for (size_t i = 1; i < u.size() - 1; ++i) {
			for (size_t j = 1; j < u[0].size() - 1; ++j) {
				Cell deltaLeft = u[i][j] - u[i - 1][j];
				Cell deltaRight = u[i + 1][j] - u[i][j];

				Cell delta = 0.5 * (deltaLeft + deltaRight);

				double r;

				if (std::fabs(deltaRight[eIndex]) < slopeTolerence) {
					r = 0;
				} else {
					r = deltaLeft[eIndex] / deltaRight[eIndex];
				}

				Cell uLeft = u[i][j] - 0.5 * delta * limit(r, slType);
				Cell uRight = u[i][j] + 0.5 * delta * limit(r, slType);

				lIfaces[i][j] = uLeft;
				rIfaces[i][j] = uRight;
			}
		}
	}
}
