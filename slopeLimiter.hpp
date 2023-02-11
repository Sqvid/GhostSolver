#ifndef GHOSTSOLVER_SLOPELIMITER_HPP
#define GHOSTSOLVER_SLOPELIMITER_HPP

#include <functional>

#include "eulerData.hpp"

using std::size_t;

namespace fvm {
	// Slope limiter options.
	enum class SlopeLimiterType {
		none = -1,
		minbee,
		superbee,
		vanAlbada,
		vanLeer
	};

	class SlopeLimiter {
		public:
			// Constructor:
			SlopeLimiter(SlopeLimiterType slType = SlopeLimiterType::none);

			void linearReconst(EulerData& eulerData,
					CellVector& lIfaces, CellVector& rIfaces);

		private:
			// Private member data:
			SlopeLimiterType slType_;
			double slopeTolerence_;

			// Private member functions:
			std::function<double (double)> limit_;
			static double xiRight_(double r);
			static double minbee_(double r);
			static double superbee_(double r);
			static double vanAlbada_(double r);
			static double vanLeer_(double r);
			static double doNothing_(double r);
	};
}

#endif
