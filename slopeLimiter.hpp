#ifndef GHOSTSOLVER_SLOPELIMITER_HPP
#define GHOSTSOLVER_SLOPELIMITER_HPP

#include <functional>

#include "eulerData.hpp"

using std::size_t;

namespace fvm {
	// Slope limiter options.
	enum class SlopeLimiter {
		none = -1,
		minbee,
		superbee,
		vanAlbada,
		vanLeer
	};

	void linearReconst(EulerData& eulerData, CellVector& lIfaces, CellVector& rIfaces, SlopeLimiter slType);
}

#endif
