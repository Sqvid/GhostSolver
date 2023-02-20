#ifndef GHOSTSOLVER_FVMSIMULATION_HPP
#define GHOSTSOLVER_FVMSIMULATION_HPP

#include <array>
#include <cstddef>
#include <fstream>
#include <functional>
#include <vector>

#include "eulerData.hpp"
#include "slopeLimiter.hpp"

using std::size_t;

namespace fvm {
	// The spatial dimensions in which the simulation is being carried out..
	enum class Axis {
		x,
		y
	};

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
			Simulation(int nCells, double xStart, double xEnd, double tStart,
					double tEnd, double cfl, double gamma,
					std::function<double (double)> densityDist,
					std::function<double (double)> velocityDistX,
					std::function<double (double)> velocityDistY,
					std::function<double (double)> pressureDist,
					FluxScheme fluxScheme,
					SlopeLimiter slType = SlopeLimiter::none);

			// Accessors
			// Getters
			/// @brief Get number of cells.
			int nCells() { return nCells_; }
			/// @brief Get the minimum x value.
			double xStart() { return xStart_; }
			/// @brief Get the maximum x value.
			double xEnd() { return xEnd_; }
			/// @brief Get the start time.
			double tStart() { return tStart_; }
			/// @brief Get the end time.
			double tEnd() { return tEnd_; }
			/// @brief Get the current simulation time.
			double tNow() { return tNow_; };
			/// @brief Get the Courant number being used.
			double cfl() { return cfl_; }
			/// @brief Get the x-cell length.
			double dx() { return dx_; }
			/// @brief Get the current timestep.
			double dt() { return dt_; }
			/// @brief Get the adiabatic constant.
			double gamma() { return gamma_; }
			/// @brief Get the underlying data.
			const EulerData& data() { return eulerData_; }

			// Public member functions.
			void step();

			// Operator overloads.
			const CellVector& operator[](size_t i) { return eulerData_[i]; }
			friend std::ofstream& operator<<(std::ofstream& output, Simulation& sim);

		private:
			// Private member data
			int nCells_;
			double xStart_;
			double xEnd_;
			double yStart_;
			double yEnd_;
			double tStart_;
			double tEnd_;
			double tNow_;
			double cfl_;
			double dx_;
			double dy_;
			double dt_;
			double gamma_;
			FluxScheme fluxScheme_;
			SlopeLimiter slType_;
			EulerData eulerData_;
			Grid flux_;
			// Left and right reconstructed interface values.
			//CellVector lSlopeIfaces_;
			//CellVector rSlopeIfaces_;

			// Private member functions
			double calcTimeStep_();
			Cell lfFlux_(const Cell& uLeft, const Cell& uRight);
			Cell richtmyerFlux_(const Cell& uLeft, const Cell& uRight);
			Cell forceFlux_(const Cell& uLeft, const Cell& uRight);
			Cell hllcFlux_(const Cell& uLeft, const Cell& uRight);
			Cell calcFlux_(const Cell& uLeft, const Cell& uRight);
			Cell fluxExpr_(Cell u);
			// Wrappers around EulerData mode conversion function.
			void setMode(EulerDataMode want) { eulerData_.setMode(want); }
	};
}

#endif
