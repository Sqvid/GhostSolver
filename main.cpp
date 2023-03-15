#include <cstddef>
#include <exception>
#include <fstream>
#include <iostream>

#include "simulation.hpp"
#include "slopeLimiter.hpp"
#include "euler2Dtests.hpp"
#include "levelSet.hpp"

void runSimulation(fvm::Simulation& sim, std::ofstream& output) {
	output << sim << "\n\n";

	// Change variables to conserved for solving.
	for (double t = sim.tStart(); t < sim.tEnd(); t = sim.tNow()) {
		sim.step();
		output << sim << "\n\n";
	}

	output << sim;
}

int main(void) {
	try {
		fvm::Simulation sodsX(100, 0, 1, 0, 0.25, 0.9, 1.4,
				&test1DensityX, &test1VelocityX, &constantVY,
				&test1PressureX, fvm::FluxScheme::hllc,
				fvm::SlopeLimiter::superbee);

		fvm::Simulation sodsY(100, 0, 1, 0, 0.25, 0.9, 1.4,
				&test1DensityY, &constantVX, &test1VelocityY,
				&test1PressureY, fvm::FluxScheme::hllc,
				fvm::SlopeLimiter::superbee);

		fvm::Simulation cylExplosion(100, -1, 1, 0, 0.25, 0.9, 1.4,
				&cylExplDensity, &cylExplVelocityX, &cylExplVelocityY,
				&cylExplPressure, fvm::FluxScheme::hllc,
				fvm::SlopeLimiter::superbee);

		fvm::Simulation singleCircle(100, 0, 1, 0, 0.4, 0.9, 1.4,
				&rigidTestDensity, &rigidTestVelocityX, &rigidTestVelocityY,
				&rigidTestPressure, fvm::FluxScheme::hllc,
				fvm::SlopeLimiter::superbee, &circleLS);

		fvm::Simulation sepCircles(100, 0, 1, 0, 0.4, 0.9, 1.4,
				&rigidTestDensity, &rigidTestVelocityX, &rigidTestVelocityY,
				&rigidTestPressure, fvm::FluxScheme::hllc,
				fvm::SlopeLimiter::superbee, &separateCirclesLS);

		fvm::Simulation overlapCircles(100, 0, 1, 0, 0.4, 0.9, 1.4,
				&rigidTestDensity, &rigidTestVelocityX, &rigidTestVelocityY,
				&rigidTestPressure, fvm::FluxScheme::hllc,
				fvm::SlopeLimiter::superbee, &overlapCirclesLS);

		std::ofstream sodsXOut("test1X.dat");
		std::ofstream sodsYOut("test1Y.dat");
		std::ofstream cylOut("cylExpl.dat");
		std::ofstream circleOut("single.dat");
		std::ofstream sepOut("separate.dat");
		std::ofstream overlapOut("overlap.dat");

		if (!sodsXOut || !sodsYOut || !cylOut) {
			std::cerr << "Error: could not create output file.\n";
			return 1;
		}

		runSimulation(sodsX, sodsXOut);
		runSimulation(sodsY, sodsYOut);
		runSimulation(cylExplosion, cylOut);
		runSimulation(singleCircle, circleOut);
		runSimulation(sepCircles, sepOut);
		runSimulation(overlapCircles, overlapOut);

	} catch (std::exception const& ex) { //std::exception const& ex){
		std::cerr << "Error: " << ex.what() << "\n";
	}

	return 0;
}
