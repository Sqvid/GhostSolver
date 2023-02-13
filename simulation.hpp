#ifndef GHOSTSOLVER_FVMSIMULATION_HPP
#define GHOSTSOLVER_FVMSIMULATION_HPP

#include <array>
#include <functional>
#include <vector>
#include <fstream>

#include "eulerData.hpp"
#include "slopeLimiter.hpp"

using std::size_t;

namespace fvm {
	// Flux scheme options.
	enum class FluxScheme {
		laxFriedrichs,
		richtmyer,
		force,
		hllc
	};

	/** The main class that holds the input parameters of the simulation.
	 * The simulation contains the parameters of the simulation as well as
	 * the associated scheme, and flux expression being used.
	 *
	 * @brief Holds the initial conditions and numerical scheme options.
	 * @param nCells Number of cells in the simulation.
	 * @param xStart Minimum x-value.
	 * @param xEnd Maximum x-value.
	 * @param tStart Start time of the simulation.
	 * @param tEnd End time of the simulation.
	 * @param cfl Courant number to use.
	 * @param gamma Adiabatic constant of the medium.
	 * @param fluxScheme Flux scheme to use.
	 * @param slType Type of slope-limiting to use, if any.
	 */
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
					SlopeLimiterType slType = SlopeLimiterType::none);

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

			// Wrapper around EulerData quantity getter.
			//QuantArray getQuantity(size_t i) { return eulerData_[i]; }

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
			SlopeLimiterType slType_;
			SlopeLimiter sLimiter_;
			std::function<double (double)> densityDist_;
			std::function<double (double)> velocityDist_;
			std::function<double (double)> pressureDist_;
			EulerData eulerData_;
			CellVector flux_;
			// Left and right reconstructed interface values.
			CellVector lSlopeIfaces_;
			CellVector rSlopeIfaces_;

			// Private member functions
			double calcTimeStep_();
			QuantArray lfFlux_(const QuantArray& uLeft, const QuantArray& uRight);
			QuantArray richtmyerFlux_(const QuantArray& uLeft, const QuantArray& uRight);
			QuantArray forceFlux_(const QuantArray& uLeft, const QuantArray& uRight);
			QuantArray hllcFlux_(const QuantArray& uLeft, const QuantArray& uRight);
			QuantArray calcFlux_(const QuantArray& uLeft, const QuantArray& uRight);
			QuantArray fluxExpr_(QuantArray u);
			// Wrappers around EulerData mode conversion function.
			void setMode(EulerDataMode want) { eulerData_.setMode(want); }
	};
}

#endif
