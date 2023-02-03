#include <cmath>
#include <cstddef>
#include <cstdio>
#include <functional>
#include <iostream>
#include <iterator>
#include <utility>
#include <vector>

#include "FvmSimulation.hpp"

using std::size_t;

namespace fvm {
	// Operator overloads for QuantArray.
	QuantArray operator+(QuantArray a, QuantArray b) {
		QuantArray ans;

		for (size_t i = 0; i < a.size(); ++i) {
			ans[i] = a[i] + b[i];
		}

		return ans;
	}

	QuantArray operator-(QuantArray a, QuantArray b) {
		QuantArray ans;

		for (size_t i = 0; i < a.size(); ++i) {
			ans[i] = a[i] - b[i];
		}

		return ans;
	}

	QuantArray operator*(double a, QuantArray u) {
		QuantArray ans;

		for (size_t i = 0; i < u.size(); ++i) {
			ans[i] = a * u[i];
		}

		return ans;
	}

	QuantArray operator*(QuantArray u, double a) {
		QuantArray ans;

		for (size_t i = 0; i < u.size(); ++i) {
			ans[i] = a * u[i];
		}

		return ans;
	}

	QuantArray operator/(QuantArray u, double a) {
		QuantArray ans;

		for (size_t i = 0; i < u.size(); ++i) {
			ans[i] = u[i] / a;
		}

		return ans;
	}

	void EulerData::makeConserved(double gamma) {
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

	void EulerData::makePrimitive(double gamma) {
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

	QuantArray EulerData::getQuantity(size_t tripletIndex) {
		size_t i = tripletIndex;

		return data_[i];
	}

	void EulerData::setQuantity(size_t tripletIndex, QuantArray value) {
		size_t i = tripletIndex;

		data_[i] = value;
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
			QuantArray cellValues;
			cellValues[dIndex] = densityDist_(x);
			cellValues[vIndex] = velocityDist_(x);
			cellValues[pIndex] = pressureDist_(x);

			eulerData_.setQuantity(i, cellValues);
		}

		// TODO: Support other boundary conditions.
		// Apply transmissive boundary conditions.
		eulerData_.setQuantity(0, eulerData_.getQuantity(1));
		eulerData_.setQuantity(nCells_ + 1, eulerData_.getQuantity(nCells_));

		dt_ = calcTimeStep_();
	}

	// Heuristic to find a timestep to keep cell updates stable.
	double Simulation::calcTimeStep_() {
		double vMax = 0;

		// Indices for primitive quantities.
		int vIndex = static_cast<int>(PrimitiveQuant::velocity);
		int pIndex = static_cast<int>(PrimitiveQuant::pressure);

		for (size_t i = 1; i < eulerData_.size(); ++i) {
			QuantArray cellValues = eulerData_.getQuantity(i);
			double velocity = cellValues[vIndex];
			double pressure = cellValues[pIndex];

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

		// Compute flux vector.
		for (size_t i = 0; i < flux_.size(); ++i) {
			flux_[i] = calcFlux_(i);
		}

		for (unsigned int i = 1; i <= nCells_; ++i) {
			QuantArray cellValues = eulerData_.getQuantity(i);
			QuantArray newValues = cellValues - (dt_/dx_) * (flux_[i] - flux_[i - 1]);

			eulerData_.setQuantity(i, newValues);
		}

		// Apply boundary conditions.
		eulerData_.setQuantity(0, eulerData_.getQuantity(1));
		eulerData_.setQuantity(nCells_ + 1, eulerData_.getQuantity(nCells_));
	}

	// Output simulation data in Gnuplot format.
	void Simulation::saveToFile(std::ofstream& output) {
		for (size_t i = 0; i < nCells_; ++i) {
			double x = xStart_ + i * dx_;

			// Indices of primitive quantities.
			int dIndex = static_cast<int>(fvm::PrimitiveQuant::density);
			int vIndex = static_cast<int>(fvm::PrimitiveQuant::velocity);
			int pIndex = static_cast<int>(fvm::PrimitiveQuant::pressure);

			// Values of primitive quantities.
			QuantArray cellValues = eulerData_.getQuantity(i);

			output << x << " "
				<< cellValues[dIndex] << " "
				<< cellValues[vIndex] << " "
				<< cellValues[pIndex] << "\n";
		}

		// Delimit Gnuplot code block with two blank lines.
		output << "\n\n";
	}

	// Return the flux for the conserved quantity in cell i.
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
	QuantArray Simulation::lfFlux_(size_t i) {
		QuantArray u = eulerData_.getQuantity(i);
		QuantArray uNext = eulerData_.getQuantity(i + 1);

		return 0.5 * (dx_/dt_) * (u - uNext)
			+ 0.5 * (fluxExpr_(uNext) + fluxExpr_(u));
	}

	// Richtmeyer flux function.
	QuantArray Simulation::richtmyerFlux_(size_t i) {
		QuantArray u = eulerData_.getQuantity(i);
		QuantArray uNext = eulerData_.getQuantity(i + 1);

		QuantArray halfStepUpdate =
			0.5 * (u + uNext)
			- 0.5 * (dt_/dx_) * (fluxExpr_(uNext) - fluxExpr_(u));

		// FIXME: This doesn't make any sense.
		return fluxExpr_(halfStepUpdate);
	}

	// Return the appropriate value for the given cell, flux scheme, and flux expression.
	QuantArray Simulation::calcFlux_(size_t i) {
		QuantArray flux;

		switch (fluxScheme_) {
			case FluxScheme::laxFriedrich:
				flux = lfFlux_(i);
				break;

			case FluxScheme::richtmyer:
				flux = richtmyerFlux_(i);
				break;

			case FluxScheme::force:
				flux = 0.5 * (lfFlux_(i) + richtmyerFlux_(i));
				break;
		}

		return flux;
	}
}
