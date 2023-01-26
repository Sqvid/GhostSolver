#include <cmath>
#include <cstddef>
#include <functional>
#include <iostream>
#include <vector>

#include "FvmSimulation.hpp"

using std::size_t;

namespace fvm {
	void EulerData::convertToConserved(double gamma) {
		int d = static_cast<int>(PrimitiveQuant::density);
		int v = static_cast<int>(PrimitiveQuant::velocity);
		int p = static_cast<int>(PrimitiveQuant::pressure);

		for (size_t i = 0; i < data_.size(); ++i) {
			// Convert velocity to momentum.
			data_[i][v] *= data_[i][d];

			// Convert pressure to energy.
			data_[i][p] = (data_[i][p] / (gamma - 1))
				+ 0.5 * ((data_[i][v]*data_[i][v]) / data_[i][d]);
		}

		mode_ = EulerDataMode::conserved;
	}

	void EulerData::convertToPrimitive(double gamma) {
		int d = static_cast<size_t>(ConservedQuant::density);
		int mo = static_cast<int>(ConservedQuant::momentum);
		int e = static_cast<int>(ConservedQuant::energy);

		for (size_t i = 0; i < data_.size(); ++i) {
			// Convert energy to pressure.
			data_[i][e] = (gamma - 1) * (data_[i][e] - 0.5 * ((data_[i][mo]*data_[i][mo]) / data_[i][d]));

			// Convert momentum to velocity.
			data_[i][mo] /= data_[i][d];
		}

		mode_ = EulerDataMode::primitive;
	}

	double EulerData::getQuantity(size_t tripletIndex, size_t quantIndex) {
		size_t i = tripletIndex;
		size_t q = quantIndex;

		return data_[i][q];
	}

	void EulerData::setQuantity(size_t tripletIndex, size_t quantIndex, double value) {
		size_t i = tripletIndex;
		size_t q = quantIndex;

		data_[i][q] = value;
	}

	// Constructor declaration
	Simulation::Simulation(unsigned int nCells, double xStart, double xEnd,
			double tStart, double tEnd, double cfl, double gamma,
			std::function<double (double)> densityDist,
			std::function<double (double)> velocityDist,
			std::function<double (double)> pressureDist,
			FluxScheme fluxScheme) {

		// TODO: Add checks to enforce sane values. E.g. Positive nCells.
		nCells_ = nCells;
		xStart_ = xStart;
		xEnd_ = xEnd;
		tStart_ = tStart;
		tEnd_ = tEnd;
		tNow_ = tStart;
		cfl_ = cfl;
		dx_ = (xEnd_ - xStart_) / nCells_;
		gamma_ = gamma;
		eulerData_.data().resize(nCells_ + 2);
		flux_.resize(nCells_ + 1);
		fluxScheme_ = fluxScheme;
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


			eulerData_.setQuantity(i, dIndex, densityDist_(x));
			eulerData_.setQuantity(i, vIndex, velocityDist_(x));
			eulerData_.setQuantity(i, pIndex, pressureDist_(x));
		}

		// TODO: Support other boundary conditions.
		// Apply transmissive boundary conditions.
		eulerData_.setQuantity(0, vIndex, eulerData_.getQuantity(1, vIndex));
		eulerData_.setQuantity(nCells_ + 1, vIndex, eulerData_.getQuantity(nCells_, vIndex));

		eulerData_.setQuantity(0, dIndex, eulerData_.getQuantity(1, dIndex));
		eulerData_.setQuantity(nCells_ + 1, dIndex, eulerData_.getQuantity(nCells_, dIndex));

		eulerData_.setQuantity(0, dIndex, eulerData_.getQuantity(1, dIndex));
		eulerData_.setQuantity(nCells_ + 1, dIndex, eulerData_.getQuantity(nCells_, dIndex));

		dt_ = calcTimeStep_();
	}

	// Heuristic to find a timestep to keep cell updates stable.
	double Simulation::calcTimeStep_() {
		double vMax = 0;

		// Indices for primitive quantities.
		int vIndex = static_cast<int>(PrimitiveQuant::velocity);
		int pIndex = static_cast<int>(PrimitiveQuant::pressure);

		for (size_t i = 1; i < eulerData_.size(); ++i) {
			double velocity = eulerData_.getQuantity(i, vIndex);
			double pressure = eulerData_.getQuantity(i, pIndex);

			// v = velocity + sound speed
			double v = std::fabs(velocity)
				+ std::sqrt((gamma_ * pressure) / pressure);

			if (v > vMax) {
				vMax = v;
			}
		}

		return cfl_ * (dx_ / vMax);
	}

	// Evolve the simulation one timestep.
	void Simulation::step() {
		dt_ = calcTimeStep_();
		tNow_ += dt_;

		int d = static_cast<int>(ConservedQuant::density);
		int mo = static_cast<int>(ConservedQuant::momentum);
		int e = static_cast<int>(ConservedQuant::energy);

		// Compute flux vector.
		for (size_t i = 0; i < flux_.size(); ++i) {
			flux_[i][d] = calcFlux_(i, ConservedQuant::density);
			flux_[i][mo] = calcFlux_(i, ConservedQuant::momentum);
			flux_[i][e] = calcFlux_(i, ConservedQuant::energy);
		}

		for (unsigned int i = 1; i <= nCells_; ++i) {
			double newDensity = eulerData_.getQuantity(i, d) - (dt_/dx_) * (flux_[i][d] - flux_[i - 1][d]);
			double newMomentum = eulerData_.getQuantity(i, mo) - (dt_/dx_) * (flux_[i][mo] - flux_[i - 1][mo]);
			double newEnergy = eulerData_.getQuantity(i, e) - (dt_/dx_) * (flux_[i][e] - flux_[i - 1][e]);

			eulerData_.setQuantity(i, d, newDensity);
			eulerData_.setQuantity(i, mo, newMomentum);
			eulerData_.setQuantity(i, e, newEnergy);
		}

		// Apply boundary conditions.
		eulerData_.setQuantity(0, 0, eulerData_.getQuantity(1, 0));
		eulerData_.setQuantity(0, 1, eulerData_.getQuantity(1, 1));
		eulerData_.setQuantity(0, 2, eulerData_.getQuantity(1, 2));

		eulerData_.setQuantity(nCells_ + 1, 0, eulerData_.getQuantity(nCells_, 0));
		eulerData_.setQuantity(nCells_ + 1, 1, eulerData_.getQuantity(nCells_, 1));
		eulerData_.setQuantity(nCells_ + 1, 2, eulerData_.getQuantity(nCells_, 2));
	}

	// Return the fluxes for the conserved quantities.
	double Simulation::fluxExpr_(size_t i, ConservedQuant quant) {
		double rhoV = eulerData_.getQuantity(i, static_cast<int>(ConservedQuant::momentum));
		double rho = eulerData_.getQuantity(i, static_cast<int>(ConservedQuant::density));
		double e = eulerData_.getQuantity(i, static_cast<int>(ConservedQuant::energy));
		double p = (gamma_ - 1) * (e - (rhoV*rhoV)/(2*rho));

		if (quant == ConservedQuant::density) {
			return rhoV;

		} else if (quant == ConservedQuant::momentum) {
			return (rhoV*rhoV)/rho + p;

		} else if (quant == ConservedQuant::energy) {
			return (e + p) * (rhoV / rho);
		}

		return 0;
	}

	double Simulation::lfFlux_(size_t i, ConservedQuant quant) {
		int q = static_cast<int>(quant);

		return 0.5 * (dx_/dt_) * (eulerData_.getQuantity(i, q) - eulerData_.getQuantity(i + 1, q))
			+ 0.5 * (fluxExpr_(i + 1, quant)
			+ fluxExpr_(i, quant));
	}

	double Simulation::richtmyerFlux_(size_t i, ConservedQuant quant) {
		int q = static_cast<int>(quant);

		double halfStepUpdate =
			0.5 * (eulerData_.getQuantity(i, q) + eulerData_.getQuantity(i + 1, q))
			- 0.5 * (dt_/dx_)
			* (fluxExpr_(i + 1, quant) - fluxExpr_(i, quant));

		return fluxExpr_(halfStepUpdate, quant);
	}

	// Return the appropriate value for the given cell, flux scheme, and flux expression.
	double Simulation::calcFlux_(size_t i, ConservedQuant quant) {
		if (fluxScheme_ == FluxScheme::laxFriedrich) {
			return lfFlux_(i, quant);

		} else if (fluxScheme_ == FluxScheme::richtmyer) {
			return richtmyerFlux_(i, quant);

		} else if (fluxScheme_ == FluxScheme::force) {
			return 0.5 * (lfFlux_(i, quant) + richtmyerFlux_(i, quant));
		}

		return 0;
	}
}
