#include <algorithm>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "levelSet.hpp"
#include "simulation.hpp"
#include "slopeLimiter.hpp"
#include "twoVector.hpp"

using std::size_t;

namespace fvm {
	// Public member function definitions:

	// Constructor:
	Simulation::Simulation(int nCells, double xStart, double xEnd,
			double tStart, double tEnd, double cfl, double gamma,
			std::function<double (double, double)> densityDist,
			std::function<double (double, double)> velocityDistX,
			std::function<double (double, double)> velocityDistY,
			std::function<double (double, double)> pressureDist,
			std::function<double (double, double, double)> levelSet,
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
		levelSet_ = levelSet;

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

		// Populate the solution space with the initial function.
		for (int i = 0; i < eulerData_.xSize(); ++i) {
			// ith cell centre x-position.
			double x = xStart_ + (i - nBoundary_ + 0.5) * dx_;

			for (int j = 0; j < eulerData_.ySize(); ++j) {
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
		populateInterfaceCells_();
		populateGhostRegion_();

		// Visualise level-set zero-contour with ascii.
		//std::cout << "\n\n";
		//for (int i = 1; i < nTotal_ - 1; ++i) {
		//	for (int j = 1; j < nTotal_ - 1; ++j) {
		//		if (isInterfaceCell(i, j)) {
		//			std::cout << "-";

		//		} else {
		//			std::cout << "@";
		//		}

		//		std::cout << " ";
		//	}

		//	std::cout << "\n";
		//}

		// All calculations must be done in conserved mode.
		eulerData_.setMode(EulerDataMode::conserved);

		dt_ = calcTimeStep_();
		tNow_ += dt_;

		// Apply boundary conditions in x.
		for (int i = 0; i < nBoundary_; ++i) {
			for (int j = nBoundary_; j < nCells_ + nBoundary_; ++j) {
				eulerData_.setCell(i, j, eulerData_[nBoundary_][j]);
				eulerData_.setCell(nCells_ + nBoundary_ + i, j,
						eulerData_[nCells_ + nBoundary_ - 1][j]);
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
				eulerData_.setCell(i, nCells_ + nBoundary_ + j,
						eulerData_[i][nCells_ + nBoundary_ - 1]);
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

		// Get the max density. This is used to scale the schlieren plot.
		double maxGradRho = 0;
		std::vector<std::vector<double>> rhoGradients(sim.nTotal_, std::vector<double>(sim.nTotal_));

		for (int i = sim.nBoundary_; i < sim.nCells_ + sim.nBoundary_; ++i) {
			for (int j = sim.nBoundary_; j < sim.nCells_ + sim.nBoundary_; ++j) {
				// Values of primitive quantities.
				// Mock-schlieren variable calculation.
				// X-derivitive of density.
				auto dRhoX = (sim[i + 1][j][dIndex] - sim[i - 1][j][dIndex]) / (2 * sim.dx_);
				// Y-derivitive of density.
				auto dRhoY = (sim[i][j + 1][dIndex] - sim[i][j - 1][dIndex]) / (2 * sim.dy_);
				TwoVector gradRho(dRhoX, dRhoY);
				auto gradMag = gradRho.mag();

				if (gradMag > maxGradRho) {
					maxGradRho = gradMag;
				}

				rhoGradients[i][j] = gradMag;
			}
		}

		for (int i = sim.nBoundary_; i < sim.nCells_ + sim.nBoundary_; ++i) {
			double x = sim.xStart_ + (i - sim.nBoundary_) * sim.dx_;

			for (int j = sim.nBoundary_; j < sim.nCells_ + sim.nBoundary_; ++j) {
				double y = sim.yStart_ + (j - sim.nBoundary_) * sim.dy_;

				// Values of primitive quantities.
				Cell cell = sim[i][j];

				// The mock-schlieren variable.
				double k = 15;
				double b = 0.8;
				auto schlieren = b * std::exp((-k * rhoGradients[i][j]) / maxGradRho);

				output << x << " "
					<< y << " "
					<< cell[dIndex] << " "
					<< cell[vIndexX] << " "
					<< cell[vIndexY] << " "
					<< cell[pIndex] << " "
					<< schlieren << "\n";
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
					Cell cellChange = 0.5 * (dt_/dx_) * (fluxExpr_(uRight, ax)
							- fluxExpr_(uLeft, ax));

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
		const Grid& u = eulerData_.data();

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
		Cell hllcL {}, hllcR {}; switch (ax) {
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

	 bool Simulation::isInterfaceCell_(int i, int j) {
		 if (i <= 0 || i >= nTotal_ - 1 || j <= 0 || j >= nTotal_ - 1) {
			 throw std::out_of_range("isInterfaceCell: Cell or neighbours are out of bounds!");
		 }

		 auto x = xStart_ + (i - nBoundary_ + 0.5) * dx_;
		 auto y = yStart_ + (j - nBoundary_ + 0.5) * dy_;

		 // the cell at i,j is an interface cell if the level set function at
		 // its cell-centre is >= 0 and at least one of its neighbours has a
		 // level set < 0.
		 return levelSet_(x, y, tNow_) >= 0 && !(levelSet_(x - dx_, y, tNow_) >= 0
				 && levelSet_(x + dx_, y, tNow_) >= 0 && levelSet_(x, y - dy_, tNow_) >= 0
				 && levelSet_(x, y + dy_, tNow_) >= 0);
	 }

	// Get cell values at point (x, y) via bilinear interpolation.
	Cell Simulation::blInterpolate_(TwoVector v) {
		auto x = v.x;
		auto y = v.y;

		if (x < xStart_ || y < yStart_) {
			throw std::logic_error("Interpolation point is outside domain!");
		}

		// Indices of lower-left cell-centre.
		int i, j;

		// x and y bounds of the lower-left cell-centre.
		double x1, x2, y1, y2;

		// Nearest cell wall, or centre.
		int xMesh = static_cast<int>(2 * ((x - xStart_) / dx_));
		int yMesh = static_cast<int>(2 * ((y - yStart_) / dy_));

		// xMesh is even # of dxs => it is a cell wall.
		if (xMesh % 2 == 0) {
			i = (xMesh / 2) + nBoundary_ - 1;

		// xMesh is the location of a cell centre.
		} else {
			i = ((xMesh + 1) / 2) + nBoundary_ - 1;
		}

		x1 = xStart_ + (i - nBoundary_ + 0.5) * dx_;
		x2 = x1 + dx_;

		if (yMesh % 2 == 0) { j = (yMesh / 2) + nBoundary_ - 1;

		// yMesh is the location of a cell centre.
		} else {
			j = ((yMesh + 1) / 2) + nBoundary_ - 1;
		}

		y1 = yStart_ + (j - nBoundary_ + 0.5) * dy_;
		y2 = y1 + dy_;

		// Lower-left, upper-left, lower-right, and upper-right cell centres
		// that bound our coordinate.
		Cell ll = eulerData_[i][j];
		Cell ul = eulerData_[i][j + 1];
		Cell lr = eulerData_[i + 1][j];
		Cell ur = eulerData_[i + 1][j + 1];

		// Linear interpolation in x.
		Cell lower = ((x2 - x) * ll + (x - x1) * lr) / (x2 - x1);
		Cell upper = ((x2 - x) * ul + (x - x1) * ur) / (x2 - x1);

		Cell interpolated =((y2 - y) * lower + (y - y1) * upper) / (y2 - y1);

		return interpolated;
	}

	void Simulation::populateInterfaceCells_() {
		eulerData_.setMode(EulerDataMode::primitive);

		for (int i = 1; i < nTotal_ - 1; ++i) {
			for (int j = 1; j < nTotal_ - 1; ++j) {
				// If this is not an interface cell there is nothing to do.
				if (!isInterfaceCell_(i, j)) {
					continue;
				}

				// Interface cell (x, y) position.
				auto xI = xStart_ + (i - nBoundary_ + 0.5) * dx_;
				auto yI = yStart_ + (j - nBoundary_ + 0.5) * dy_;

				// vI is an interface cell centre.
				TwoVector vI(xI, yI);
				// nI is the normal vector at point (xI, yI)
				auto nI = findNormal(levelSet_, xI, yI, tNow_);

				// vP is the closest point to vI that lies on the level set
				// zero-contour.
				auto vP = vI - levelSet_(xI, yI, tNow_)*nI;

				// These are the two "adjacent" cells for the Riemann GFM.
				// vP1 is position vector for the cell inside the ghost-region,
				// while vP2 is for the cell in the real-region.
				auto vP2 = vP - 1.5 * dx_ * nI;

				// Get real state via bilinear interpolation.
				auto realState = blInterpolate_(vP2);

				auto vRealX = realState[vIndexX];
				auto vRealY = realState[vIndexY];
				TwoVector v(vRealX, vRealY);

				// Normal velocity.
				double vNormalMag = v * nI;
				auto vTangent = v - (vNormalMag * nI);
				//auto vTangentMag = vTangent.mag();

				auto rotRealState = realState;
				rotRealState[vIndexX] = -vNormalMag;
				//rotRealState[vIndexY] = vTangentMag;

				// The ghost state is the same as the real state but with the normal
				// velocity flipped.
				auto rotGhostState = rotRealState;
				rotGhostState[vIndexX] = vNormalMag;

				double rhoL = rotGhostState[dIndex];
				double vxL = rotGhostState[vIndexX];
				//double vyL = rotGhostState[vIndexY];
				double pL = rotGhostState[pIndex];

				double rhoR = rotRealState[dIndex];
				double vxR = rotRealState[vIndexX];
				//double vyR = rotRealState[vIndexY];
				double pR = rotRealState[pIndex];

				// Sound speed estimates.
				double cSoundL = std::sqrt((gamma_ * pL) / rhoL);
				double cSoundR = std::sqrt((gamma_ * pR) / rhoR);

				double sPlus {}, sL {}, sR {}, sStar {};
				Cell hllcR {};
				// Find approximate left and right sound speeds.
				sPlus = std::max(std::fabs(vxL) + cSoundL, std::fabs(vxR) + cSoundR);
				sL = -sPlus;
				sR = sPlus;

				// Approximate contact velocity in x.
				sStar = (pR - pL + rhoL*vxL*(sL - vxL) - rhoR*vxR*(sR - vxR))
					/ (rhoL*(sL - vxL) - rhoR*(sR - vxR));

				hllcR = rhoR * ((sR - vxR) / (sR - sStar))
					* Cell({1, sStar, 0,
							rotRealState[eIndex]/rhoR + (sStar - vxR)*(sStar + pR/(rhoR * (sR - vxR)))});

				Cell rotInterState = hllcR;

				//if ( sL >= 0 ) {
				//	rotInterState = rotGhostState;

				//} else if (sStar >= 0) {
				//	rotInterState = hllcL;

				//} else if (sR >= 0) {
				//	rotInterState = hllcR;

				//} else {
				//	rotInterState = rotRealState;
				//}

				Cell interState = rotInterState;
				auto vInter = rotInterState[vIndexX]*nI + vTangent;
				interState[vIndexX] = vInter.x;
				interState[vIndexY] = vInter.y;

				// Assign the intermediate state to the current interface cell.
				eulerData_[i][j] = interState;
			}
		}
	}

	void Simulation::populateGhostRegion_() {
		double unknown = 1e100;
		double huge = 2 * unknown;

		// Initialise sweepGrid;
		Grid sweepGrid(nTotal_);
		for (auto& vec : sweepGrid) {
			vec.resize(nTotal_, Cell({-huge, -huge, -huge, -huge}));
		}

		const int sweepX = nTotal_ - 1;
		const int sweepY = nTotal_ - 1;

		// Initialise values inside interface.
		for (int i = 1; i < sweepX; ++i) {
			auto x = xStart_ + (i - nBoundary_ + 0.5) * dx_;

			for (int j = 1; j < sweepY; ++j) {
				auto y = yStart_ + (j - nBoundary_ + 0.5) * dy_;

				if (isInterfaceCell_(i, j)) {
					sweepGrid[i][j] = eulerData_[i][j];

				} else if (levelSet_(x, y, tNow_) >= 0) {
					sweepGrid[i][j] = Cell({huge, huge, huge, huge});
				}
			}
		}

		// Lambda that does constant extrapolation.
		auto doExtrapolation = [&](int i, int j, double x, double y) {
			auto& cell = sweepGrid[i][j];

			// Extrapolate every quantity via Eikonal equation.
			for (size_t q = 0; q < cell.size(); ++q) {

				auto quant = sweepGrid[i][j][q];

				// Cell is available for an update.
				if (levelSet_(x, y, tNow_) > 0 && !isInterfaceCell_(i, j)) {
					auto qX = std::min(sweepGrid[i - 1][j][q], sweepGrid[i + 1][j][q]);
					auto qY = std::min(sweepGrid[i][j - 1][q], sweepGrid[i][j + 1][q]);
					auto n = findNormal(levelSet_, x, y, tNow_);

					// Helper variable.
					auto beta = std::abs((n.y * dx_) / (n.x * dy_));

					auto qGhost = (qX + beta*qY) / (1 + beta);

					if (qGhost < quant) {
						cell[q] = qGhost;
					}
				}
			}
		};

		// Positive x-sweep.
		for (int i = 1; i < sweepX; ++i) {
			auto x = xStart_ + (i - nBoundary_ + 0.5) * dx_;

			for (int j = 1; j < sweepY; ++j) {
				auto y = yStart_ + (j - nBoundary_ + 0.5) * dy_;

				// Extrapolate every quantity via Eikonal equation.
				doExtrapolation(i, j, x, y);
			}
		}

		// Positive y-sweep.
		for (int j = 1; j < sweepY; ++j) {
			auto y = yStart_ + (j - nBoundary_ + 0.5) * dy_;

			for (int i = 1; i < sweepX; ++i) {
				auto x = xStart_ + (i - nBoundary_ + 0.5) * dx_;

				// Extrapolate every quantity via Eikonal equation.
				doExtrapolation(i, j, x, y);
			}
		}

		// Negative x-sweep.
		for (int i = sweepGrid.size() - 1; i > 0; --i) {
			auto x = xStart_ + (i - nBoundary_ + 0.5) * dx_;

			for (int j = sweepGrid[0].size() - 1; j > 0; --j) {
				auto y = yStart_ + (j - nBoundary_ + 0.5) * dy_;

				// Extrapolate every quantity via Eikonal equation.
				doExtrapolation(i, j, x, y);
			}
		}

		// Negative y-sweep.
		for (int j = sweepGrid[0].size() - 1; j > 0; --j) {
			auto y = yStart_ + (j - nBoundary_ + 0.5) * dy_;

			for (int i = sweepGrid.size() - 1; i > 0; --i) {
				auto x = xStart_ + (i - nBoundary_ + 0.5) * dx_;

				// Extrapolate every quantity via Eikonal equation.
				doExtrapolation(i, j, x, y);
			}
		}

		// Copy values over to the real mesh.
		for (int i = 1; i < nTotal_ - 1; ++i) {
			auto x = xStart_ + (i - nBoundary_ + 0.5) * dx_;

			for (int j = 1; j < nTotal_ - 1; ++j) {
				auto y = yStart_ + (j - nBoundary_ + 0.5) * dy_;

				if (levelSet_(x, y, tNow_) > 0 && !isInterfaceCell_(i, j)) {
					eulerData_[i][j] = sweepGrid[i][j];
				}
			}
		}
	}
}
