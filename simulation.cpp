#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <iostream>
#include <stdexcept>

#include "simulation.hpp"

using std::size_t;

namespace fvm {
	// Public member function definitions:

	// Constructor:
	Simulation::Simulation(unsigned int nCells, double xStart, double xEnd,
			double tStart, double tEnd, double cfl, double gamma,
			std::function<double (double)> densityDist,
			std::function<double (double)> velocityDist,
			std::function<double (double)> pressureDist,
			FluxScheme fluxScheme,
			SlopeLimiterType slType)
			// Initialiser list
			: sLimiter_(slType), eulerData_(gamma) {

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
		eulerData_.data().resize(nCells_ + 2);
		flux_.resize(nCells_ + 1);
		fluxScheme_ = fluxScheme;
		lSlopeIfaces_.resize(nCells_);
		rSlopeIfaces_.resize(nCells_);
		slType_ = slType;
		densityDist_ = densityDist;
		velocityDist_ = velocityDist;
		pressureDist_ = pressureDist;

		// Indices for primitive quantities.
		int dIndex = static_cast<int>(PrimitiveQuant::density);
		int vIndex = static_cast<int>(PrimitiveQuant::velocity);
		int pIndex = static_cast<int>(PrimitiveQuant::pressure);

		// Populate the solution space with the initail function.
		for (size_t i = 0; i < eulerData_.data().size(); ++i) {
			// ith cell centre
			double x = xStart_ + (i - 0.5) * dx_;
			QuantArray cellValues;
			cellValues[dIndex] = densityDist_(x);
			cellValues[vIndex] = velocityDist_(x);
			cellValues[pIndex] = pressureDist_(x);

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
		if (slType_ != SlopeLimiterType::none) {
			// FIXME: Calculate these from reconstructed boundary
			// cells.
			flux_[0] = calcFlux_(eulerData_[0], eulerData_[1]);
			flux_[nCells_] = calcFlux_(eulerData_[nCells_], eulerData_[nCells_ + 1]);

			sLimiter_.linearReconst(eulerData_, lSlopeIfaces_, rSlopeIfaces_);

			// Half-timestep evolution.
			for (std::size_t i = 0; i < nCells_; ++i) {
				QuantArray& uLeft = lSlopeIfaces_[i];
				QuantArray& uRight = rSlopeIfaces_[i];

				QuantArray cellChange = 0.5 * (dt_/dx_)
					* (fluxExpr_(uRight) - fluxExpr_(uLeft));

				uLeft = uLeft - cellChange;
				uRight = uRight - cellChange;
			}

			for (size_t i= 0; i < nCells_ - 1; ++i) {
				QuantArray uRight = rSlopeIfaces_[i];
				QuantArray uNextLeft = lSlopeIfaces_[i + 1];

				flux_[i + 1] = calcFlux_(uRight, uNextLeft);
			}

		} else {
			// Compute flux vector.
			for (size_t i = 0; i < flux_.size(); ++i) {
				flux_[i] = calcFlux_(eulerData_[i], eulerData_[i + 1]);
			}
		}

		for (unsigned int i = 1; i <= nCells_; ++i) {
			QuantArray cellValues = eulerData_[i];
			QuantArray newValues = cellValues - (dt_/dx_) * (flux_[i] - flux_[i - 1]);

			eulerData_.setQuantity(i, newValues);
		}

		// Apply boundary conditions.
		eulerData_.setQuantity(0, eulerData_[1]);
		eulerData_.setQuantity(nCells_ + 1, eulerData_[nCells_]);
	}

	// Output simulation data in Gnuplot format.
	void Simulation::saveToFile(std::ofstream& output) {
		// Output data in primitive form.
		eulerData_.setMode(EulerDataMode::primitive);

		for (size_t i = 0; i < nCells_; ++i) {
			double x = xStart_ + i * dx_;

			// Indices of primitive quantities.
			int dIndex = static_cast<int>(fvm::PrimitiveQuant::density);
			int vIndex = static_cast<int>(fvm::PrimitiveQuant::velocity);
			int pIndex = static_cast<int>(fvm::PrimitiveQuant::pressure);

			// Values of primitive quantities.
			QuantArray cellValues = eulerData_[i];

			output << x << " "
				<< cellValues[dIndex] << " "
				<< cellValues[vIndex] << " "
				<< cellValues[pIndex] << "\n";
		}
	}

	// Private member function definitions:

	// Heuristic to find a timestep to keep cell updates stable.
	// Timestep is computed such that it is stable for the fastest
	// information wave.
	double Simulation::calcTimeStep_() {
		double vMax = 0;

		// Indices for conserved quantities.
		int dIndex = static_cast<int>(ConservedQuant::density);
		int moIndex = static_cast<int>(ConservedQuant::momentum);
		int eIndex = static_cast<int>(ConservedQuant::energy);

		for (size_t i = 1; i < eulerData_.size(); ++i) {
			QuantArray cellValues = eulerData_[i];
			double rho = cellValues[dIndex];
			double rhoV = cellValues[moIndex];
			double e = cellValues[eIndex];

			// v = velocity + sound speed
			double v = std::fabs(rhoV / rho)
				+ std::sqrt((gamma_*(gamma_ - 1)*(e - 0.5*((rhoV*rhoV) / rho))) / rho);

			if (v > vMax) {
				vMax = v;
			}
		}

		return cfl_ * (dx_ / vMax);
	}

	// Return the fluxes for the conserved quantities.
	QuantArray Simulation::fluxExpr_(QuantArray u) {
		QuantArray flux;

		int dIndex = static_cast<int>(ConservedQuant::density);
		int moIndex = static_cast<int>(ConservedQuant::momentum);
		int eIndex = static_cast<int>(ConservedQuant::energy);

		double rho = u[dIndex];
		double rhoV = u[moIndex];
		double e = u[eIndex];
		double p = (gamma_ - 1) * (e - (rhoV*rhoV)/(2*rho));

		flux[dIndex] = rhoV;
		flux[moIndex] = (rhoV*rhoV)/rho + p;
		flux[eIndex] = (e + p) * (rhoV / rho);

		return flux;
	}

	// Lax-Friedrichs flux function.
	QuantArray Simulation::lfFlux_(const QuantArray& uLeft, const QuantArray& uRight) {
		return 0.5 * (dx_/dt_) * (uLeft - uRight)
			+ 0.5 * (fluxExpr_(uRight) + fluxExpr_(uLeft));
	}

	// Richtmeyer flux function.
	QuantArray Simulation::richtmyerFlux_(const QuantArray& uLeft, const QuantArray& uRight) {
		QuantArray halfStepUpdate =
			0.5 * (uLeft + uRight)
			- 0.5 * (dt_/dx_) * (fluxExpr_(uRight) - fluxExpr_(uLeft));

		return fluxExpr_(halfStepUpdate);
	}

	QuantArray Simulation::forceFlux_(const QuantArray& uLeft, const QuantArray& uRight) {
		return 0.5 * (lfFlux_(uLeft, uRight) + richtmyerFlux_(uLeft, uRight));
	}

	QuantArray Simulation::hllcFlux_(const QuantArray& uLeft, const QuantArray& uRight) {
		QuantArray flux;

		// Make primitive version copies of the given cell values.
		QuantArray primL = makePrimQuants(uLeft, gamma_);
		QuantArray primR = makePrimQuants(uRight, gamma_);

		// Indices for primitive quantities.
		int dIndex = static_cast<int>(PrimitiveQuant::density);
		int vIndex = static_cast<int>(PrimitiveQuant::velocity);
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

		} else if (sL < 0 && sStar >= 0) {
			QuantArray hllcL = dL * ((sL - vL) / (sL - sStar))
					* QuantArray({1, sStar,
					uLeft[eIndex]/dL + (sStar - vL)*(sStar + pL/(dL * (sL - vL)))});

			flux = fluxExpr_(uLeft) + sL * (hllcL - uLeft);

		} else if (sStar < 0 && sR >= 0) {
			QuantArray hllcR = dR * ((sR - vR) / (sR - sStar))
					* QuantArray({1, sStar,
					uRight[eIndex]/dR + (sStar - vR)*(sStar + pR/(dR * (sR - vR)))});

			flux = fluxExpr_(uRight) + sL * (hllcR - uRight);

		} else if (sR < 0) {
			flux = fluxExpr_(uRight);
		}

		return flux;
	}

	// Return the appropriate value for the given cell, flux scheme, and flux expression.
	QuantArray Simulation::calcFlux_(const QuantArray& uLeft, const QuantArray& uRight) {
		QuantArray flux;

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
