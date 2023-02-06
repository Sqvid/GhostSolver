#ifndef GHOSTSOLVER_SLOPELIMITER_HPP
#define GHOSTSOLVER_SLOPELIMITER_HPP

#include <functional>

#include "fvmData.hpp"

using std::size_t;

namespace fvm {
	// Slope limiter options.
	enum class SlopeLimiterType {
		none = -1,
		minbee
	};

	class SlopeLimiter {
		public:
			// Constructor:
			SlopeLimiter(SlopeLimiterType slType = SlopeLimiterType::none);

			void linearReconstruct(EulerData& eulerData, CellVector& interfaces);

		private:
			// Private member data:
			SlopeLimiterType slType_;

			// Private member functions:
			std::function<double (double)> limit_;

			constexpr static double xiLeft_(double r);
			constexpr static double xiRight_(double r);
			constexpr static double minbee_(double r);
			constexpr static double doNothing_(double r);
	};
}

#endif
