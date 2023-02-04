#ifndef GHOSTSOLVER_FVMSIMULATION_HPP
#define GHOSTSOLVER_FVMSIMULATION_HPP

#include <array>
#include <functional>
#include <vector>
#include <fstream>

#include "fvmData.hpp"

using std::size_t;

namespace fvm {
	// Flux scheme options.
	enum class FluxScheme {
		laxFriedrichs,
		richtmyer,
		force,
	};

	// flux scheme.
	enum class SlopeLimiter {
		noLimiter = -1,
		minbee
	};

	// The simulation contains the parameters of the simulation as well as the
	// associated scheme, and flux expression being used.
	class Simulation {
		public:
			// Constructor
			Simulation(unsigned int nCells, double xStart, double xEnd,
					double tStart, double tEnd, double cfl,
					double gamma,
					std::function<double (double)> densityDist,
					std::function<double (double)> velocityDist,
					std::function<double (double)> pressureDist,
					FluxScheme fluxScheme,
					SlopeLimiter slopeLimiter = SlopeLimiter::noLimiter);

			// Accessors
			// Getters
			unsigned int nCells() { return nCells_; }
			double xStart() { return xStart_; }
			double xEnd() { return xEnd_; }
			double tStart() { return tStart_; }
			double tEnd() { return tEnd_; }
			double tNow() { return tNow_; };
			double cfl() { return cfl_; }
			double dx() { return dx_; }
			double dt() { return dt_; }
			double gamma() { return gamma_; }
			EulerData data() { return eulerData_; }

			// Wrapper around EulerData qunatity getter.
			QuantArray getQuantity(size_t i) { return eulerData_.getQuantity(i); }

			// Public member functions
			void step();
			void saveToFile(std::ofstream& output);

		private:
			// Private member data
			unsigned int nCells_;
			double xStart_;
			double xEnd_;
			double tStart_;
			double tEnd_;
			double tNow_;
			double cfl_;
			double dx_;
			double dt_;
			double gamma_;
			FluxScheme fluxScheme_;
			std::function<double (double)> densityDist_;
			std::function<double (double)> velocityDist_;
			std::function<double (double)> pressureDist_;
			EulerData eulerData_;
			TripletVector flux_;

			// Private member functions
			double calcTimeStep_();
			QuantArray lfFlux_(size_t i);
			QuantArray richtmyerFlux_(size_t i);
			QuantArray calcFlux_(size_t i);
			QuantArray fluxExpr_(QuantArray u);
			// Wrappers around EulerData conversion functions.
			void makeConserved_() { eulerData_.makeConserved(); }
			void makePrimitive_() { eulerData_.makePrimitive(); }

	};
}

#endif
