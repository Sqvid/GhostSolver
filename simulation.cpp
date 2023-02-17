#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "simulation.hpp"
#include "slopeLimiter.hpp"

namespace fvm {
	// Public member function definitions:

	// Constructor:
	Simulation::Simulation(int nCells, double xStart, double xEnd,
			double tStart, double tEnd, double cfl, double gamma,
			std::function<double (double)> densityDist,
			std::function<double (double)> velocityDistX,
			std::function<double (double)> velocityDistY,
			std::function<double (double)> pressureDist,
			FluxScheme fluxScheme,
			SlopeLimiter slType)
			// Initialiser list
			: eulerData_(gamma) {

		// Check for sane input values.
		if (nCells <= 0) {
			throw std::invalid_argument("nCells must be > 0.");
		} else if (tStart < 0) {
			throw std::invalid_argument("tStart must be >= 0");
		} else if (tEnd <= tStart) {
			throw std::invalid_argument("tEnd must be > tStart");
		}

		nCells_ = nCells;
		xStart_ = xStart;
		xEnd_ = xEnd;
		tStart_ = tStart;
		tEnd_ = tEnd;
		tNow_ = tStart;
		cfl_ = cfl;
		dx_ = (xEnd_ - xStart_) / nCells_;
		dt_ = 0;
		gamma_ = gamma;
		fluxScheme_ = fluxScheme;
		slType_ = slType;

		eulerData_.data().resize(nCells_ + 2);
		flux_.resize(nCells_ + 1);
		lSlopeIfaces_.resize(nCells_);
		rSlopeIfaces_.resize(nCells_);

		// Indices for primitive quantities.
		int dIndex = static_cast<int>(PrimitiveQuant::density);
		int vIndexX = static_cast<int>(PrimitiveQuant::velocityX);
		int vIndexY = static_cast<int>(PrimitiveQuant::velocityY);
		int pIndex = static_cast<int>(PrimitiveQuant::pressure);

		// Populate the solution space with the initail function.
		for (size_t i = 0; i < eulerData_.data().size(); ++i) {
			// ith cell centre
			double x = xStart_ + (i - 0.5) * dx_;

			Cell cellValues;
			cellValues[dIndex] = densityDist(x);
			cellValues[vIndexX] = velocityDistX(x);
			cellValues[vIndexY] = velocityDistY(x);
			cellValues[pIndex] = pressureDist(x);

			eulerData_.setQuantity(i, cellValues);
		}

		// Apply transmissive boundary conditions.
		eulerData_.setQuantity(0, eulerData_[1]);
		eulerData_.setQuantity(nCells_ + 1, eulerData_[nCells_]);
	}

	// Evolve the simulation one timestep.
	void Simulation::step() {
		// All calculations must be done in conserved mode.
		eulerData_.setMode(EulerDataMode::conserved);

		dt_ = calcTimeStep_();
		tNow_ += dt_;

		// If slope limiting has been requested.
		if (slType_ != SlopeLimiter::none) {
			// FIXME: Calculate these from reconstructed boundary
			// cells.
			flux_[0] = calcFlux_(eulerData_[0], eulerData_[1]);
			flux_[nCells_] = calcFlux_(eulerData_[nCells_], eulerData_[nCells_ + 1]);

			internal::linearReconst(eulerData_, lSlopeIfaces_, rSlopeIfaces_, slType_);

			// Half-timestep evolution.
			for (std::size_t i = 0; i < nCells_; ++i) {
				Cell& uLeft = lSlopeIfaces_[i];
				Cell& uRight = rSlopeIfaces_[i];

				Cell cellChange = 0.5 * (dt_/dx_) * (fluxExpr_(uRight) - fluxExpr_(uLeft));

				uLeft = uLeft - cellChange;
				uRight = uRight - cellChange;
			}

			for (size_t i = 0; i < nCells_ - 1; ++i) {
				Cell uRight = rSlopeIfaces_[i];
				Cell uNextLeft = lSlopeIfaces_[i + 1];

				flux_[i + 1] = calcFlux_(uRight, uNextLeft);
			}

		} else {
			// Compute flux vector.
			for (size_t i = 0; i < flux_.size(); ++i) {
				flux_[i] = calcFlux_(eulerData_[i], eulerData_[i + 1]);
			}
		}

		for (unsigned int i = 1; i <= nCells_; ++i) {
			Cell cellValues = eulerData_[i];
			Cell newValues = cellValues - (dt_/dx_) * (flux_[i] - flux_[i - 1]);

			eulerData_.setQuantity(i, newValues);
		}

		// Apply boundary conditions.
		eulerData_.setQuantity(0, eulerData_[1]);
		eulerData_.setQuantity(nCells_ + 1, eulerData_[nCells_]);
	}

	// Output simulation data in Gnuplot format.
	std::ofstream& operator<<(std::ofstream& output, Simulation& sim) {
		// Output data in primitive form.
		sim.eulerData_.setMode(EulerDataMode::primitive);

		for (size_t i = 0; i < sim.nCells_; ++i) {
			double x = sim.xStart_ + i * sim.dx_;

			// Indices of primitive quantities.
			int dIndex = static_cast<int>(PrimitiveQuant::density);
			int vIndexX = static_cast<int>(PrimitiveQuant::velocityX);
			int pIndex = static_cast<int>(PrimitiveQuant::pressure);

			// Values of primitive quantities.
			Cell cellValues = sim[i];

			output << x << " "
				<< cellValues[dIndex] << " "
				<< cellValues[vIndexX] << " "
				<< cellValues[pIndex] << "\n";
		}

		return output;
	}

	// Private member function definitions:

	// Heuristic to find a timestep to keep cell updates stable.
	// Timestep is computed such that it is stable for the fastest
	// information wave.
	double Simulation::calcTimeStep_() {
		double vMax = 0;

		// Indices for conserved quantities.
		int dIndex = static_cast<int>(ConservedQuant::density);
		int moIndexX = static_cast<int>(ConservedQuant::momentumX);
		int moIndexY = static_cast<int>(ConservedQuant::momentumY);
		int eIndex = static_cast<int>(ConservedQuant::energy);

		for (size_t i = 1; i < eulerData_.size(); ++i) {
			Cell cellValues = eulerData_[i];
			double rho = cellValues[dIndex];
			double rhoVX = cellValues[moIndexX];
			double rhoVY = cellValues[moIndexY];
			double e = cellValues[eIndex];

			// v = velocity + sound speed
			auto cSound = std::sqrt((gamma_*(gamma_ - 1)
						*(e - 0.5*((rhoVX*rhoVX + rhoVY*rhoVY) / rho))) / rho);

			auto vX = std::fabs(rhoVX / rho) + cSound;

			if (vX > vMax) {
				vMax = vX;
			}
		}

		return cfl_ * (dx_ / vMax);
	}

	// Return the fluxes for the conserved quantities.
	Cell Simulation::fluxExpr_(Cell u) {
		Cell flux;

		int dIndex = static_cast<int>(ConservedQuant::density);
		int moIndexX = static_cast<int>(ConservedQuant::momentumX);
		int moIndexY = static_cast<int>(ConservedQuant::momentumY);
		int eIndex = static_cast<int>(ConservedQuant::energy);

		double rho = u[dIndex];
		double rhoVX = u[moIndexX];
		double rhoVY = u[moIndexY];
		double e = u[eIndex];
		double p = (gamma_ - 1) * (e - (rhoVX*rhoVX)/(2*rho));

		flux[dIndex] = rhoVX;
		flux[moIndexX] = (rhoVX * rhoVX)/rho + p;
		flux[moIndexY] = (rhoVX * rhoVY)/rho;
		flux[eIndex] = (e + p) * (rhoVX / rho);

		return flux;
	}

	// Lax-Friedrichs flux function.
	Cell Simulation::lfFlux_(const Cell& uLeft, const Cell& uRight) {
		return 0.5 * (dx_/dt_) * (uLeft - uRight)
			+ 0.5 * (fluxExpr_(uRight) + fluxExpr_(uLeft));
	}

	// Richtmeyer flux function.
	Cell Simulation::richtmyerFlux_(const Cell& uLeft, const Cell& uRight) {
		Cell halfStepUpdate =
			0.5 * (uLeft + uRight)
			- 0.5 * (dt_/dx_) * (fluxExpr_(uRight) - fluxExpr_(uLeft));

		return fluxExpr_(halfStepUpdate);
	}

	Cell Simulation::forceFlux_(const Cell& uLeft, const Cell& uRight) {
		return 0.5 * (lfFlux_(uLeft, uRight) + richtmyerFlux_(uLeft, uRight));
	}

	// @brief Helper function to convert QuantArray variables to primitive form.
	// @note Useful for Riemann-based schemes.
	// @warn Does not check original state of variables. Use with caution!
	Cell makePrimQuants(const Cell& u, const double gamma) {
		Cell prim;
		prim = u;

		int dIndex = static_cast<size_t>(ConservedQuant::density);
		int moIndex = static_cast<int>(ConservedQuant::momentumX);
		int eIndex = static_cast<int>(ConservedQuant::energy);

		// Convert energy to pressure.
		prim[eIndex] = (gamma - 1) * (u[eIndex]
				- 0.5 * ((u[moIndex] * u[moIndex]) / u[dIndex]));

		// Convert momentum to velocity.
		prim[moIndex] = u[moIndex] / u[dIndex];

		return prim;
	}

	Cell Simulation::hllcFlux_(const Cell& uLeft, const Cell& uRight) {
		Cell flux;

		// Make primitive version copies of the given cell values.
		Cell primL = makePrimQuants(uLeft, gamma_);
		Cell primR = makePrimQuants(uRight, gamma_);

		// Indices for primitive quantities.
		int dIndex = static_cast<int>(PrimitiveQuant::density);
		int vIndex = static_cast<int>(PrimitiveQuant::velocityX);
		int pIndex = static_cast<int>(PrimitiveQuant::pressure);

		// Index for energy.
		int eIndex = static_cast<int>(ConservedQuant::energy);

		// Alias quantities for convenience.
		double dL = primL[dIndex];
		double vL = primL[vIndex];
		double pL = primL[pIndex];

		double dR = primR[dIndex];
		double vR = primR[vIndex];
		double pR = primR[pIndex];

		double cSoundL = std::sqrt((gamma_ * pL) / dL);
		double cSoundR = std::sqrt((gamma_ * pR) / dR);

		// Find approximate left and right sound speeds.
		double sPlus = std::max(std::fabs(vL) + cSoundL, std::fabs(vR) + cSoundR);

		double sL = -sPlus;
		double sR = sPlus;

		// Approximate contact sound speed.
		double sStar = (pR - pL + dL*vL*(sL - vL) - dR*vR*(sR - vR))
					/ (dL*(sL - vL) - dR*(sR - vR));

		if ( sL >= 0 ) {
			flux = fluxExpr_(uLeft);

		} else if (sStar >= 0) {
			Cell hllcL = dL * ((sL - vL) / (sL - sStar))
					* Cell({1, sStar,
					uLeft[eIndex]/dL + (sStar - vL)*(sStar + pL/(dL * (sL - vL)))});

			flux = fluxExpr_(uLeft) + sL * (hllcL - uLeft);

		} else if (sR >= 0) {
			Cell hllcR = dR * ((sR - vR) / (sR - sStar))
					* Cell({1, sStar,
					uRight[eIndex]/dR + (sStar - vR)*(sStar + pR/(dR * (sR - vR)))});

			flux = fluxExpr_(uRight) + sR * (hllcR - uRight);

		} else {
			flux = fluxExpr_(uRight);
		}

		return flux;
	}

	// Return the appropriate value for the given cell, flux scheme, and flux expression.
	Cell Simulation::calcFlux_(const Cell& uLeft, const Cell& uRight) {
		Cell flux;

		switch (fluxScheme_) {
			case FluxScheme::laxFriedrichs:
				flux = lfFlux_(uLeft, uRight);
				break;

			case FluxScheme::richtmyer:
				flux = richtmyerFlux_(uLeft, uRight);
				break;

			case FluxScheme::force:
				flux = forceFlux_(uLeft, uRight);
				break;

			case FluxScheme::hllc:
				flux = hllcFlux_(uLeft, uRight);
				break;
		}

		return flux;
	}
}
