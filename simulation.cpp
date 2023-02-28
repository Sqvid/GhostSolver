#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "simulation.hpp"
#include "slopeLimiter.hpp"

namespace fvm {
	// Public member function definitions:

	// Constructor:
	Simulation::Simulation(int nCells, double xStart, double xEnd,
			double tStart, double tEnd, double cfl, double gamma,
			std::function<double (double, double)> densityDist,
			std::function<double (double, double)> velocityDistX,
			std::function<double (double, double)> velocityDistY,
			std::function<double (double, double)> pressureDist,
			FluxScheme fluxScheme,
			SlopeLimiter slType)
			// Initialiser list
			: eulerData_(gamma) {

		// Check for sane input values.
		if (nCells <= 0) {
			throw std::invalid_argument("nCells must be > 0.");

		} else if (tStart < 0) {
			throw std::invalid_argument("tStart must be >= 0.");

		} else if (tEnd <= tStart) {
			throw std::invalid_argument("tEnd must be > tStart.");
		}

		nCells_ = nCells;
		nBoundary_ = 2;
		nTotal_ = nCells_ + 2*nBoundary_;
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

		// Resize grids in x.
		eulerData_.data().resize(nTotal_);
		flux_.resize(nTotal_);
		lSlopeIfaces_.resize(nTotal_);
		rSlopeIfaces_.resize(nTotal_);

		// Resize grids in y.
		for (int i = 0; i < nTotal_; ++i) {
			eulerData_.data()[i].resize(nTotal_);
			flux_[i].resize(nTotal_);
			lSlopeIfaces_[i].resize(nTotal_);
			rSlopeIfaces_[i].resize(nTotal_);
		}

		// Indices for primitive quantities.
		constexpr int dIndex = static_cast<int>(PrimitiveQuant::density);
		constexpr int vIndexX = static_cast<int>(PrimitiveQuant::velocityX);
		constexpr int vIndexY = static_cast<int>(PrimitiveQuant::velocityY);
		constexpr int pIndex = static_cast<int>(PrimitiveQuant::pressure);

		// Populate the solution space with the initial function.
		for (size_t i = 0; i < eulerData_.xSize(); ++i) {
			// ith cell centre x-position.
			double x = xStart_ + (i - nBoundary_ + 0.5) * dx_;

			for (size_t j = 0; j < eulerData_.ySize(); ++j) {
				// ith cell centre y-position.
				double y = yStart_ + (j - nBoundary_ + 0.5) * dy_;

				Cell cell;
				cell[dIndex] = densityDist(x, y);
				cell[vIndexX] = velocityDistX(x, y);
				cell[vIndexY] = velocityDistY(x, y);
				cell[pIndex] = pressureDist(x, y);

				eulerData_.setCell(i, j, cell);
			}
		}
	}

	// Evolve the simulation one timestep.
	void Simulation::step() {
		// All calculations must be done in conserved mode.
		eulerData_.setMode(EulerDataMode::conserved);

		dt_ = calcTimeStep_();
		tNow_ += dt_;

		// Apply boundary conditions in x.
		for (int i = 0; i < nBoundary_; ++i) {
			for (int j = nBoundary_; j < nCells_ + nBoundary_; ++j) {
				eulerData_.setCell(i, j, eulerData_[nBoundary_][j]);
				eulerData_.setCell(nCells_ + nBoundary_ + i, j, eulerData_[nCells_ + nBoundary_ - 1][j]);
			}
		}

		calcFluxGrid_(Axis::x);

		// Apply finite difference calculations in x.
		for (int i = nBoundary_; i < nCells_ + nBoundary_; ++i) {
			for (int j = nBoundary_; j < nCells_ + nBoundary_; ++j) {
				Cell diff = (dt_/dx_) * (flux_[i][j] - flux_[i - 1][j]);
				Cell newCell = eulerData_[i][j] - diff;

				eulerData_.setCell(i, j, newCell);
			}
		}

		// Apply boundary conditions in y.
		for (int i = nBoundary_; i < nCells_ + nBoundary_; ++i) {
			for (int j = 0; j < nBoundary_; ++j) {
				eulerData_.setCell(i, j, eulerData_[i][nBoundary_]);
				eulerData_.setCell(i, nCells_ + nBoundary_ + j, eulerData_[i][nCells_ + nBoundary_ - 1]);
			}
		}

		calcFluxGrid_(Axis::y);

		// Apply finite difference calculations in y.
		for (int i = nBoundary_; i < nCells_ + nBoundary_; ++i) {
			for (int j = nBoundary_; j < nCells_ + nBoundary_; ++j) {
				Cell diff = (dt_/dy_) * (flux_[i][j] - flux_[i][j - 1]);
				Cell newCell = eulerData_[i][j] - diff;

				eulerData_.setCell(i, j, newCell);
			}
		}
	}

	// Output simulation data in Gnuplot format.
	std::ofstream& operator<<(std::ofstream& output, Simulation& sim) {
		// Output data in primitive form.
		sim.eulerData_.setMode(EulerDataMode::primitive);

		for (int i = sim.nBoundary_; i < sim.nCells_ + sim.nBoundary_; ++i) {
			double x = sim.xStart_ + (i - sim.nBoundary_) * sim.dx_;

			for (int j = sim.nBoundary_; j < sim.nCells_ + sim.nBoundary_; ++j) {
				double y = sim.yStart_ + (j - sim.nBoundary_) * sim.dy_;

				// Indices of primitive quantities.
				constexpr int dIndex = static_cast<int>(PrimitiveQuant::density);
				constexpr int vIndexX = static_cast<int>(PrimitiveQuant::velocityX);
				constexpr int vIndexY = static_cast<int>(PrimitiveQuant::velocityY);
				constexpr int pIndex = static_cast<int>(PrimitiveQuant::pressure);

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

	void Simulation::calcFluxGrid_(Axis ax) {
		if (slType_ != SlopeLimiter::none) {
			// Linearly reconstruct the data and store interface values in
			// lSlopeIfaces_, and rSlopeIfaces_.
			linearReconst_(ax);

			// Half-timestep evolution.
			for (int i = nBoundary_; i < nCells_ + nBoundary_; ++i) {
				for (int j = nBoundary_; j < nCells_ + nBoundary_; ++j) {
					Cell& uLeft = lSlopeIfaces_[i][j];
					Cell& uRight = rSlopeIfaces_[i][j];

					// FIXME: This will fail if dx != dy.
					Cell cellChange = 0.5 * (dt_/dx_) * (fluxExpr_(uRight, ax) - fluxExpr_(uLeft, ax));

					uLeft = uLeft - cellChange;
					uRight = uRight - cellChange;
				}
			}

			// Calculate fluxes with the half-evolved interface values.
			switch (ax) {
				case Axis::x:
					for (int i = nBoundary_ - 1; i < nCells_ + nBoundary_; ++i) {
						for (int j = nBoundary_ - 1; j < nCells_ + nBoundary_; ++j) {
							Cell uRight = rSlopeIfaces_[i][j];
							Cell uNextLeft = lSlopeIfaces_[i + 1][j];

							flux_[i][j] = calcFlux_(uRight, uNextLeft, ax);
						}
					}

					break;

				case Axis::y:
					for (int i = nBoundary_ - 1; i < nCells_ + nBoundary_; ++i) {
						for (int j = nBoundary_ - 1; j < nCells_ + nBoundary_; ++j) {
							Cell uRight = rSlopeIfaces_[i][j];
							Cell uNextLeft = lSlopeIfaces_[i][j + 1];

							flux_[i][j] = calcFlux_(uRight, uNextLeft, ax);
						}
					}

					break;
			}

		} else {
			switch (ax) {
				case Axis::x:
					for (int i = nBoundary_ - 1; i < nCells_ + nBoundary_; ++i) {
						for (int j = nBoundary_ - 1; j < nCells_ + nBoundary_; ++j) {
							flux_[i][j] = calcFlux_(eulerData_[i][j], eulerData_[i + 1][j], ax);
						}
					}

					break;

				case Axis::y:
					for (int i = nBoundary_ - 1; i < nCells_ + nBoundary_; ++i) {
						for (int j = nBoundary_ - 1; j < nCells_ + nBoundary_; ++j) {
							flux_[i][j] = calcFlux_(eulerData_[i][j], eulerData_[i][j + 1], ax);
						}
					}

					break;
			}
		}
	}

	// Private member function definitions:
	// ====================================================================== //

	// Heuristic to find a timestep to keep cell updatas stable.
	// Timestep is computed such that it is stable for the fastest
	// information wave.
	double Simulation::calcTimeStep_() {
		// Indices for conserved quantities.
		constexpr int dIndex = static_cast<int>(ConservedQuant::density);
		constexpr int moIndexX = static_cast<int>(ConservedQuant::momentumX);
		constexpr int moIndexY = static_cast<int>(ConservedQuant::momentumY);
		constexpr int eIndex = static_cast<int>(ConservedQuant::energy);

		double max = 0;
		for (int i = nBoundary_; i < nCells_ + nBoundary_; ++i) {
			for (int j = nBoundary_; j < nCells_ + nBoundary_; ++j) {
				Cell cell = eulerData_[i][j];

				double rho = cell[dIndex];
				double rhoVX = cell[moIndexX];
				double rhoVY = cell[moIndexY];
				double e = cell[eIndex];

				// v = velocity + sound speed
				auto cSound = std::sqrt((gamma_*(gamma_ - 1)
							*(e - 0.5*((rhoVX*rhoVX + rhoVY*rhoVY) / rho))) / rho);

				auto sX = std::fabs(rhoVX / rho) + cSound;
				auto sY = std::fabs(rhoVY / rho) + cSound;

				auto r = std::max(sX / dx_, sY / dy_);

				if (r > max) {
					max = r;
				}
			}
		}

		return cfl_ / max;
	}

	// FIXME: Adjust boundary when reconstruction is fixed with extra ghost
	// cells.
	void Simulation::linearReconst_(Axis ax) {
		// Alias the numerical data.
		eulerData_.setMode(EulerDataMode::conserved);
		Grid& u = eulerData_.data();

		// Index for energy. This is the quantity we are going to limit
		// on.
		constexpr int eIndex = static_cast<int>(ConservedQuant::energy);

		// Calculate reconstructed interface values.
		for (int i = nBoundary_ - 1; i < nCells_ + nBoundary_ + 1; ++i) {
			for (int j = nBoundary_ - 1; j < nCells_ + nBoundary_ + 1; ++j) {
				Cell deltaLeft {};
				Cell deltaRight {};

				switch (ax) {
					case Axis::x:
						deltaLeft = u[i][j] - u[i - 1][j];
						deltaRight = u[i + 1][j] - u[i][j];
						break;

					case Axis::y:
						deltaLeft = u[i][j] - u[i][j - 1];
						deltaRight = u[i][j + 1] - u[i][j];
						break;
				}

				Cell delta = 0.5 * (deltaLeft + deltaRight);

				double r;

				if (std::fabs(deltaRight[eIndex]) < internal::slopeTolerence) {
					r = 0;
				} else {
					r = deltaLeft[eIndex] / deltaRight[eIndex];
				}

				Cell uLeft = u[i][j] - 0.5 * delta * internal::limit(r, slType_);
				Cell uRight = u[i][j] + 0.5 * delta * internal::limit(r, slType_);

				lSlopeIfaces_[i][j] = uLeft;
				rSlopeIfaces_[i][j] = uRight;
			}
		}
	}

	// Return the fluxes for the conserved quantities.
	Cell Simulation::fluxExpr_(Cell u, Axis ax) {
		Cell flux;

		constexpr int dIndex = static_cast<int>(ConservedQuant::density);
		constexpr int moIndexX = static_cast<int>(ConservedQuant::momentumX);
		constexpr int moIndexY = static_cast<int>(ConservedQuant::momentumY);
		constexpr int eIndex = static_cast<int>(ConservedQuant::energy);

		double rho = u[dIndex];

		if (rho == 0) {
			throw std::logic_error("Density in flux expression must never be zero!");
		}

		double rhoVX = u[moIndexX];
		double rhoVY = u[moIndexY];
		double e = u[eIndex];
		double p = (gamma_ - 1) * (e - (rhoVX*rhoVX + rhoVY*rhoVY)/(2*rho));

		switch (ax) {
			case Axis::x:
				flux[dIndex] = rhoVX;
				flux[moIndexX] = (rhoVX * rhoVX)/rho + p;
				flux[moIndexY] = (rhoVX * rhoVY)/rho;
				flux[eIndex] = (e + p) * (rhoVX / rho);
				break;

			case Axis::y:
				flux[dIndex] = rhoVY;
				flux[moIndexX] = (rhoVX * rhoVY)/rho;
				flux[moIndexY] = (rhoVY * rhoVY)/rho + p;
				flux[eIndex] = (e + p) * (rhoVY / rho);
				break;
		}

		return flux;
	}

	// Lax-Friedrichs flux function.
	Cell Simulation::lfFlux_(const Cell& uLeft, const Cell& uRight, Axis ax) {
		return 0.5 * (dx_/dt_) * (uLeft - uRight)
			+ 0.5 * (fluxExpr_(uRight, ax) + fluxExpr_(uLeft, ax));
	}

	// Richtmeyer flux function.
	Cell Simulation::richtmyerFlux_(const Cell& uLeft, const Cell& uRight, Axis ax) {
		//Cell halfStepUpdate = (0.5 * (uLeft + uRight)) - (0.5 * (dt_/dx_) * (fluxExpr_(uRight) - fluxExpr_(uLeft)));
		Cell half1 = 0.5 * (uLeft + uRight);
		Cell half2 = 0.5 * (dt_/dx_) * (fluxExpr_(uRight, ax) - fluxExpr_(uLeft, ax));
		Cell halfStepUpdate = half1 - half2;

		return fluxExpr_(halfStepUpdate, ax);
	}

	Cell Simulation::forceFlux_(const Cell& uLeft, const Cell& uRight, Axis ax) {
		return 0.5 * (lfFlux_(uLeft, uRight, ax) + richtmyerFlux_(uLeft, uRight, ax));
	}

	// @brief Helper function to convert Cell variables to primitive form.
	// @note Useful for Riemann-based schemes.
	// @warn Does not check original state of variables. Use with caution!
	Cell makePrimQuants(const Cell& u, const double gamma) {
		Cell prim;
		prim = u;

		constexpr int dIndex = static_cast<size_t>(ConservedQuant::density);
		constexpr int moIndexX = static_cast<int>(ConservedQuant::momentumX);
		constexpr int moIndexY = static_cast<int>(ConservedQuant::momentumY);
		constexpr int eIndex = static_cast<int>(ConservedQuant::energy);

		auto rho = u[dIndex];
		auto rhoVX = u[moIndexX];
		auto rhoVY = u[moIndexY];
		auto e = u[eIndex];

		// Convert energy to pressure.
		prim[eIndex] = (gamma - 1) * (e - 0.5 * ((rhoVX*rhoVX + rhoVY*rhoVY) / rho));

		// Convert momentum to velocity.
		prim[moIndexX] = rhoVX / rho;
		prim[moIndexY] = rhoVY / rho;

		return prim;
	}

	Cell Simulation::hllcFlux_(const Cell& uLeft, const Cell& uRight, Axis ax) {
		Cell flux;

		// Make primitive version copies of the given cell values.
		Cell primL = makePrimQuants(uLeft, gamma_);
		Cell primR = makePrimQuants(uRight, gamma_);

		// Indices for primitive quantities.
		constexpr int dIndex = static_cast<int>(PrimitiveQuant::density);
		constexpr int vIndexX = static_cast<int>(PrimitiveQuant::velocityX);
		constexpr int vIndexY = static_cast<int>(PrimitiveQuant::velocityY);
		constexpr int pIndex = static_cast<int>(PrimitiveQuant::pressure);

		// Index for energy.
		constexpr int eIndex = static_cast<int>(ConservedQuant::energy);

		// Alias quantities for convenience.
		double rhoL = primL[dIndex];
		double vxL = primL[vIndexX];
		double vyL = primL[vIndexY];
		double pL = primL[pIndex];

		double rhoR = primR[dIndex];
		double vxR = primR[vIndexX];
		double vyR = primR[vIndexY];
		double pR = primR[pIndex];

		// Sound speed estimates.
		double cSoundL = std::sqrt((gamma_ * pL) / rhoL);
		double cSoundR = std::sqrt((gamma_ * pR) / rhoR);

		double sPlus {}, sL {}, sR {}, sStar {};
		Cell hllcL {}, hllcR {};
		switch (ax) {
			case Axis::x:
				// Find approximate left and right sound speeds.
				sPlus = std::max(std::fabs(vxL) + cSoundL, std::fabs(vxR) + cSoundR);
				sL = -sPlus;
				sR = sPlus;

				// Approximate contact velocity in x.
				sStar = (pR - pL + rhoL*vxL*(sL - vxL) - rhoR*vxR*(sR - vxR))
					/ (rhoL*(sL - vxL) - rhoR*(sR - vxR));

				// Approximate intermediate x-states.
				hllcL = rhoL * ((sL - vxL) / (sL - sStar))
						* Cell({1, sStar, vyL,
						uLeft[eIndex]/rhoL + (sStar - vxL)*(sStar + pL/(rhoL * (sL - vxL)))});

				hllcR = rhoR * ((sR - vxR) / (sR - sStar))
					* Cell({1, sStar, vyR,
					uRight[eIndex]/rhoR + (sStar - vxR)*(sStar + pR/(rhoR * (sR - vxR)))});

				break;

			case Axis::y:
				// Find approximate left and right sound speeds.
				sPlus = std::max(std::fabs(vyL) + cSoundL, std::fabs(vyR) + cSoundR);
				sL = -sPlus;
				sR = sPlus;

				// Approximate contact velocity in y.
				sStar = (pR - pL + rhoL*vyL*(sL - vyL) - rhoR*vyR*(sR - vyR))
					/ (rhoL*(sL - vyL) - rhoR*(sR - vyR));

				// Approximate intermediate y-states.
				hllcL = rhoL * ((sL - vyL) / (sL - sStar))
						* Cell({1, vxL, sStar,
						uLeft[eIndex]/rhoL + (sStar - vyL)*(sStar + pL/(rhoL * (sL - vyL)))});

				hllcR = rhoR * ((sR - vyR) / (sR - sStar))
					* Cell({1, vxR, sStar,
					uRight[eIndex]/rhoR + (sStar - vyR)*(sStar + pR/(rhoR * (sR - vyR)))});


				break;
		}

		if ( sL >= 0 ) {
			flux = fluxExpr_(uLeft, ax);

		} else if (sStar >= 0) {
			flux = fluxExpr_(uLeft, ax) + sL * (hllcL - uLeft);

		} else if (sR >= 0) {
			flux = fluxExpr_(uRight, ax) + sR * (hllcR - uRight);

		} else {
			flux = fluxExpr_(uRight, ax);
		}

		return flux;
	}

	// Return the appropriate value for the given cell, flux scheme, and flux expression.
	Cell Simulation::calcFlux_(const Cell& uLeft, const Cell& uRight, Axis ax) {
		Cell flux;

		switch (fluxScheme_) {
			case FluxScheme::laxFriedrichs:
				flux = lfFlux_(uLeft, uRight, ax);
				break;

			case FluxScheme::richtmyer:
				flux = richtmyerFlux_(uLeft, uRight, ax);
				break;

			case FluxScheme::force:
				flux = forceFlux_(uLeft, uRight, ax);
				break;

			case FluxScheme::hllc:
				flux = hllcFlux_(uLeft, uRight, ax);
				break;
		}

		return flux;
	}

	void Simulation::findInterface() {
		// Map of positions inside and outside level-set function;
		std::vector<std::vector<int>> lsMap;

		// Will track the positons that are inside (1) and outside (0) our
		// level-set function.
		lsMap.resize(nTotal_);
		for (auto& col : lsMap) {
			col.resize(nTotal_);
		}

		// Map of interface cells. Interfaces are represented as 1's; all other
		// cells are 0.
		auto ifMap = lsMap;

		auto circleLS = [](double x, double y) {
			auto x0 {0.0}, y0 {0.0}, r {.4};

			return r*r - (x - x0)*(x - x0) - (y - y0)*(y - y0);
		};

		// Find all cells on or inside the boundary.
		for (int i = 0; i < nTotal_; ++i) {
			auto x = xStart_ + (i - nBoundary_ + 0.5) * dx_;

			for (int j = 0; j < nTotal_; ++j) {
				auto y = yStart_ + (j - nBoundary_ + 0.5) * dy_;

				lsMap[i][j] = (circleLS(x, y) >= 0);
			}
		}

		// Find interface cells.
		// This works by finding cells that have a different sign from one or
		// more of their neighbours.
		for (int i = 1; i < nTotal_ - 1; ++i) {
			for (int j = 1; j < nTotal_ - 1; ++j) {
				// Since lsMap was populated with 0's and 1's, we can use some
				// boolean logic here.
				if (lsMap[i][j] && !((lsMap[i][j + 1] && lsMap[i][j - 1])
							&& (lsMap[i + 1][j] && lsMap[i - 1][j]))) {

					ifMap[i][j] = 1;
				}
			}
		}
	}
}
