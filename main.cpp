#include <cstddef>
#include <exception>
#include <fstream>
#include <iostream>

#include "simulation.hpp"
#include "slopeLimiter.hpp"
#include "toroTests.hpp"

void runSimulation(fvm::Simulation& sim, std::ofstream& output) {
	sim.saveToFile(output);
	output << "\n\n";

	// Change variables to conserved for solving.
	int nFrame = 0;
	for (double t = sim.tStart(); t < sim.tEnd(); t = sim.tNow()) {
		sim.step();
		++nFrame;

		if (nFrame % 3 == 0) {
			// Change variables back to primitive for output.
			sim.saveToFile(output);
			output << "\n\n";
		}
	}

	sim.saveToFile(output);
}

int main(void) {
	try {
		fvm::Simulation test1(150, 0, 1, 0, 0.25, 0.8, 1.4,
				&test1Density, &test1Velocity, &test1Pressure,
				fvm::FluxScheme::hllc, fvm::SlopeLimiterType::none);

		fvm::Simulation test2(150, 0, 1, 0, 0.150, 0.8, 1.4,
				&test2Density, &test2Velocity, &test2Pressure,
				fvm::FluxScheme::hllc, fvm::SlopeLimiterType::none);

		fvm::Simulation test3(150, 0, 1, 0, 0.012, 0.8, 1.4,
				&test3Density, &test3Velocity, &test3Pressure,
				fvm::FluxScheme::hllc, fvm::SlopeLimiterType::none);

		fvm::Simulation test4(150, 0, 1, 0, 0.035, 0.8, 1.4,
				&test4Density, &test4Velocity, &test4Pressure,
				fvm::FluxScheme::hllc, fvm::SlopeLimiterType::none);

		fvm::Simulation test5(150, 0, 1, 0, 0.035, 0.8, 1.4,
				&test5Density, &test5Velocity, &test5Pressure,
				fvm::FluxScheme::hllc, fvm::SlopeLimiterType::none);

		std::ofstream output1("test1.dat");
		std::ofstream output2("test2.dat");
		std::ofstream output3("test3.dat");
		std::ofstream output4("test4.dat");
		std::ofstream output5("test5.dat");

		if (!output1 || !output2 || !output3 || !output4 || !output5) {
			std::cerr << "Error: could not create output file.\n";
			return 1;
		}

		runSimulation(test1, output1);
		runSimulation(test2, output2);
		runSimulation(test3, output3);
		runSimulation(test4, output4);
		runSimulation(test5, output5);

	} catch (std::exception const& ex) { //std::exception const& ex){
		std::cerr << ex.what() << "\n";
	}

	return 0;
}
