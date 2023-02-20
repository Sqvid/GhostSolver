#ifndef GHOSTSOLVER_SLOPELIMITER_HPP
#define GHOSTSOLVER_SLOPELIMITER_HPP

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

	// Not intended for the end-user.
	namespace internal {
		/** @brief Finds a linear reconstruction of #EulerData, and fills #lIfaces
		 * and #rIfaces with the interface values.
		 *
		 * @param eulerData Original piecewise volume-averaged data.
		 * @param lIfaces Reference to a CellVector in which to store left-interface
		 * values.
		 * @param rIfaces Reference to a CellVector in which to store
		 * right-interface values.
		 * @param the type of slope-limiter to use.
		 */
		void linearReconst(EulerData& eulerData, Grid& lIfaces, Grid& rIfaces, SlopeLimiter slType);
	}
}

#endif
