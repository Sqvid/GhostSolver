#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <stdexcept>

#include "eulerData.hpp"
#include "slopeLimiter.hpp"

namespace fvm {
	const double slopeTolerence = 0.0001;

	SlopeLimiter::SlopeLimiter(SlopeLimiterType slType) {
		slType_ = slType;

		switch (slType) {
			case SlopeLimiterType::minbee:
				limit_ = minbee_;
				break;

			case SlopeLimiterType::none:
				limit_ = doNothing_;
				break;
		}
	}

	constexpr double SlopeLimiter::xiLeft_(double r) {
		return (2*r) / (1 + r);
	}

	constexpr double SlopeLimiter::xiRight_(double r) {
		return 2 / (1 + r);
	}

	constexpr double SlopeLimiter::minbee_(double r) {
		if (r <= 0) {
			return 0;
		} if (r > 0 && r <= 1) {
			return r;
		} else {
			return std::min(1.0, xiRight_(r));
		}
	}

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
			if (std::fabs(deltaLeft[eIndex]) < slopeTolerence) {
				deltaLeft[eIndex] = 0;
			}

			if (std::fabs(deltaRight[eIndex]) < slopeTolerence) {
				r = 10;
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
