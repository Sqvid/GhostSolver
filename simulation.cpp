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
		yStart_ = xStart_;
		yEnd_ = xEnd_;
		tStart_ = tStart;
		tEnd_ = tEnd;
		tNow_ = tStart;
		cfl_ = cfl;
		dx_ = (xEnd_ - xStart_) / nCells_;
		dy_ = dx_;
		dt_ = 0;
		gamma_ = gamma;
		fluxScheme_ = fluxScheme;
		slType_ = slType;

		// Resize data and flux grids in x.
		eulerData_.data().resize(nCells_ + 2);
		flux_.resize(nCells_ + 1);

		// Resize data and flux grids in y.
		for (int i = 0; i < nCells_ + 2; ++i) {
			if (i < nCells_ + 1) {
				flux_[i].resize(nCells_ + 1);
			}

			eulerData_.data()[i].resize(nCells_ + 2);
		}

		//lSlopeIfaces_.resize(nCells_);
		//rSlopeIfaces_.resize(nCells_);

		// Indices for primitive quantities.
		int dIndex = static_cast<int>(PrimitiveQuant::density);
		int vIndexX = static_cast<int>(PrimitiveQuant::velocityX);
		int vIndexY = static_cast<int>(PrimitiveQuant::velocityY);
		int pIndex = static_cast<int>(PrimitiveQuant::pressure);

		// Populate the solution space with the initail function.
		for (size_t i = 0; i < eulerData_.data().size(); ++i) {
			// ith cell centre x-position.
			double x = xStart_ + (i - 0.5) * dx_;

			for (size_t j = 0; j < eulerData_.data()[0].size(); ++j) {
				// ith cell centre y-position.
				double y = yStart_ + (j - 0.5) * dy_;

				Cell cellValues;
				cellValues[dIndex] = densityDist(x);
				cellValues[vIndexX] = velocityDistX(x);
				cellValues[vIndexY] = velocityDistY(y);
				cellValues[pIndex] = pressureDist(x);

				eulerData_.setQuantity(i, j, cellValues);
			}
		}

		// Apply transmissive boundary conditions.
		for (int j = 1; j < nCells; ++j) {
			// x-boundaries.
			eulerData_.setQuantity(0, j, eulerData_[1][j]);
			eulerData_.setQuantity(nCells_ + 1, j, eulerData_[nCells_][j]);
		}
	}

	// Evolve the simulation one timestep.
	void Simulation::step() {
		// All calculations must be done in conserved mode.
		eulerData_.setMode(EulerDataMode::conserved);

		dt_ = calcTimeStep_();
		tNow_ += dt_;

		// If slope limiting has been requested.
		/*if (slType_ != SlopeLimiter::none) {
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
		*/
			// Compute flux vector.
		for (size_t i = 0; i < flux_.size(); ++i) {
			for (size_t j = 0; j < flux_[0].size(); ++j) {
				// x-fluxes.
				flux_[i][j] = calcFlux_(eulerData_[i][j], eulerData_[i + 1][j]);
			}
		}
		//}

		// Apply finite difference calculations.
		for (int i = 1; i < nCells_ + 1; ++i) {
			for (int j = 1; j < nCells_ + 1; ++j) {
				Cell cell = eulerData_[i][j];
				Cell newCell = cell - (dt_/dx_) * (flux_[i][j] - flux_[i - 1][j]);

				eulerData_.setQuantity(i, j, newCell);
			}
		}

		// Apply boundary conditions in x.
		for (int j = 1; j < nCells_ + 1; ++j) {
			eulerData_.setQuantity(0, j, eulerData_[1][j]);
			eulerData_.setQuantity(nCells_ + 1, j, eulerData_[nCells_][j]);
		}
	}

	// Output simulation data in Gnuplot format.
	std::ofstream& operator<<(std::ofstream& output, Simulation& sim) {
		// Output data in primitive form.
		sim.eulerData_.setMode(EulerDataMode::primitive);

		for (int i = 1; i < sim.nCells_ + 1; ++i) {
			double x = sim.xStart_ + (i - 1) * sim.dx_;

			for (int j = 1; j < sim.nCells_ + 1; ++j) {
				double y = sim.yStart_ + (j - 1) * sim.dy_;

				// Indices of primitive quantities.
				int dIndex = static_cast<int>(PrimitiveQuant::density);
				int vIndexX = static_cast<int>(PrimitiveQuant::velocityX);
				int vIndexY = static_cast<int>(PrimitiveQuant::velocityY);
				int pIndex = static_cast<int>(PrimitiveQuant::pressure);

				// Values of primitive quantities.
				Cell cell = sim[i][j];

				output << x << " "
					<< y << " "
					<< cell[dIndex] << " "
					<< cell[vIndexX] << " "
					<< cell[vIndexY] << " "
					<< cell[pIndex] << "\n";
			}

			output << "\n";
		}

		return output;
	}

	// Private member function definitions:

	// Heuristic to find a timestep to keep cell updatas stable.
	// Timestep is computed such that it is stable for the fastest
	// information wave.
	double Simulation::calcTimeStep_() {
		double vMax = 0;

		// Indices for conserved quantities.
		int dIndex = static_cast<int>(ConservedQuant::density);
		int moIndexX = static_cast<int>(ConservedQuant::momentumX);
		int moIndexY = static_cast<int>(ConservedQuant::momentumY);
		int eIndex = static_cast<int>(ConservedQuant::energy);

		for (int i = 1; i < nCells_ + 1; ++i) {
			for (int j = 1; j < nCells_ + 1; ++j) {
				Cell cell = eulerData_[i][j];

				double rho = cell[dIndex];
				double rhoVX = cell[moIndexX];
				double rhoVY = cell[moIndexY];
				double e = cell[eIndex];

				// v = velocity + sound speed
				auto cSound = std::sqrt((gamma_*(gamma_ - 1)
							*(e - 0.5*((rhoVX*rhoVX + rhoVY*rhoVY) / rho))) / rho);

				auto vX = std::fabs(rhoVX / rho) + cSound;

				if (vX > vMax) {
					vMax = vX;
				}
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
		//Cell halfStepUpdate = (0.5 * (uLeft + uRight)) - (0.5 * (dt_/dx_) * (fluxExpr_(uRight) - fluxExpr_(uLeft)));
		Cell half1 = 0.5 * (uLeft + uRight);
		Cell half2 = 0.5 * (dt_/dx_) * (fluxExpr_(uRight) - fluxExpr_(uLeft));
		Cell halfStepUpdate = half1 - half2;

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

	// TODO: Make axis-specific.
	Cell Simulation::hllcFlux_(const Cell& uLeft, const Cell& uRight) {
		Cell flux;

		// Make primitive version copies of the given cell values.
		Cell primL = makePrimQuants(uLeft, gamma_);
		Cell primR = makePrimQuants(uRight, gamma_);

		// Indices for primitive quantities.
		int dIndex = static_cast<int>(PrimitiveQuant::density);
		int vIndexX = static_cast<int>(PrimitiveQuant::velocityX);
		int vIndexY = static_cast<int>(PrimitiveQuant::velocityY);
		int pIndex = static_cast<int>(PrimitiveQuant::pressure);

		// Index for energy.
		int eIndex = static_cast<int>(ConservedQuant::energy);

		// Alias quantities for convenience.
		double dL = primL[dIndex];
		double vxL = primL[vIndexX];
		double vyL = primL[vIndexY];
		double pL = primL[pIndex];

		double dR = primR[dIndex];
		double vxR = primR[vIndexX];
		double vyR = primR[vIndexY];
		double pR = primR[pIndex];

		double cSoundL = std::sqrt((gamma_ * pL) / dL);
		double cSoundR = std::sqrt((gamma_ * pR) / dR);

		// Find approximate left and right sound speeds.
		double sPlus = std::max(std::fabs(vxL) + cSoundL, std::fabs(vxR) + cSoundR);

		double sL = -sPlus;
		double sR = sPlus;

		// Approximate contact sound speed.
		double sStar = (pR - pL + dL*vxL*(sL - vxL) - dR*vxR*(sR - vxR))
					/ (dL*(sL - vxL) - dR*(sR - vxR));

		if ( sL >= 0 ) {
			flux = fluxExpr_(uLeft);

		} else if (sStar >= 0) {
			Cell hllcL = dL * ((sL - vxL) / (sL - sStar))
					* Cell({1, sStar, vyL,
					uLeft[eIndex]/dL + (sStar - vxL)*(sStar + pL/(dL * (sL - vxL)))});

			flux = fluxExpr_(uLeft) + sL * (hllcL - uLeft);

		} else if (sR >= 0) {
			Cell hllcR = dR * ((sR - vxR) / (sR - sStar))
					* Cell({1, sStar, vyR,
					uRight[eIndex]/dR + (sStar - vxR)*(sStar + pR/(dR * (sR - vxR)))});

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
