#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <stdexcept>

#include "eulerData.hpp"
#include "slopeLimiter.hpp"

namespace fvm {
	SlopeLimiter::SlopeLimiter(SlopeLimiterType slType) {
		slopeTolerence_ = 0.00001;
		slType_ = slType;

		switch (slType) {
			case SlopeLimiterType::none:
				limit_ = doNothing_;
				break;

			case SlopeLimiterType::minbee:
				limit_ = minbee_;
				break;

			case SlopeLimiterType::superbee:
				limit_ = superbee_;
				break;

			case SlopeLimiterType::vanAlbada:
				limit_ = superbee_;
				break;

			case SlopeLimiterType::vanLeer:
				limit_ = superbee_;
				break;
		}
	}

	bool compare(const double& a, const double& b) {
		return a < b;
	}

	// TODO: Is this needed?
	constexpr double SlopeLimiter::xiLeft_(double r) {
		return (2*r) / (1 + r);
	}

	constexpr double SlopeLimiter::xiRight_(double r) {
		return 2 / (1 + r);
	}

	constexpr double SlopeLimiter::minbee_(double r) {
		if (r <= 0) {
			return 0;

		} else if (r > 0 && r <= 1) {
			return r;

		} else {
			return std::min(1.0, xiRight_(r));
		}
	}

	constexpr double SlopeLimiter::superbee_(double r) {
		if (r <= 0) {
			return 0;

		} if (r > 0 && r <= 0.5) {
			return 2*r;

		} else if (r > 0.5 && r <= 1) {
			return 1;

		} else {
			return std::min({r, xiRight_(r), 2.0},
					[](const double& a, const double& b) {
						return a < b;
					});
		}
	}

	constexpr double SlopeLimiter::vanAlbada_(double r) {
		if (r <= 0) {
			return 0;

		} else {
			return std::min((r*(1 + r))/(1 + r*r), xiRight_(r));
		}
	}

	constexpr double SlopeLimiter::vanLeer_(double r) {
		if (r <= 0) {
			return 0;

		} else {
			return std::min(xiLeft_(r), xiRight_(r));
		}
	}

	// TODO: There must be a better way than requiring such a function.
	constexpr double SlopeLimiter::doNothing_(double r) {
		return r;
	}

	void SlopeLimiter::linearReconst(EulerData& eulerData, CellVector& lIfaces, CellVector& rIfaces) {
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

			// This is probably a terrible hack to avoid dividing by
			// zero.
			//if (std::fabs(deltaLeft[eIndex]) < slopeTolerence_) {
			//	deltaLeft[eIndex] = 0;
			//}

			if (std::fabs(deltaRight[eIndex]) < slopeTolerence_) {
				r = 0;
			} else {
				r = deltaLeft[eIndex] / deltaRight[eIndex];
			}

			QuantArray uLeft = u[i] - 0.5 * delta * limit_(r);
			QuantArray uRight = u[i] + 0.5 * delta * limit_(r);

			lIfaces[i - 1] = uLeft;
			rIfaces[i - 1] = uRight;
		}
	}
}
