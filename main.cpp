#include <cstddef>
#include <exception>
#include <fstream>
#include <iostream>

#include "simulation.hpp"
#include "slopeLimiter.hpp"
#include "toro2DTests.hpp"

void runSimulation(fvm::Simulation& sim, std::ofstream& output) {
	output << sim << "\n\n";

	// Change variables to conserved for solving.
	int nFrame = 0;
	for (double t = sim.tStart(); t < sim.tEnd(); t = sim.tNow()) {
		sim.step();
		++nFrame;

		if (nFrame % 5 == 0) {
			// Change variables back to primitive for output.
			output << sim << "\n\n";
		}
	}

	output << sim;
}

int main(void) {
	try {
		fvm::Simulation test1(100, -1, 1, 0, 0.25, 0.9, 1.4,
				&cylExplDensity, &cylExplVelocityX, &cylExplVelocityY, &cylExplPressure,
				fvm::FluxScheme::hllc, fvm::SlopeLimiter::vanLeer);

		std::ofstream output1("test1.dat");

		if (!output1) {
			std::cerr << "Error: could not create output file.\n";
			return 1;
		}

		runSimulation(test1, output1);

	} catch (std::exception const& ex) { //std::exception const& ex){
		std::cerr << "Error: " << ex.what() << "\n";
	}

	return 0;
}
