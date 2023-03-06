#ifndef GHOSTSOLVER_SLOPELIMITER_HPP
#define GHOSTSOLVER_SLOPELIMITER_HPP

#include "eulerData.hpp"

namespace fvm {
	// Slope limiter options.
	enum class SlopeLimiter {
		none = -1,
		minbee,
		superbee,
		vanAlbada,
		vanLeer
	};

	// Not intended for the end-user.
	namespace internal {
		const double slopeTolerence = 1e-6;

		double minbee(double r);
		double superbee(double r);
		double vanAlbada(double r);
		double vanLeer(double r);

		double limit(double r, SlopeLimiter slType);
	}
}

#endif
