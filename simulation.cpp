#include <fstream>
#include <iostream>

#include "fvmSimulation.hpp"
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
	fvm::Simulation test1(200, 0, 1, 0, 0.25, 0.8, 1.4, &test1Density, &test1Velocity,
			&test1Pressure, fvm::FluxScheme::force, fvm::SlopeLimiterType::minbee);
	fvm::Simulation test2(200, 0, 1, 0, 0.25, 0.8, 1.4, &test1Density, &test1Velocity,
			&test1Pressure, fvm::FluxScheme::force, fvm::SlopeLimiterType::none);
	//fvm::Simulation test2(200, 0, 1, 0, 0.15, 0.8, 1.4, &test2Density, &test2Velocity,
	//		&test2Pressure, fvm::FluxScheme::force);
	fvm::Simulation test3(200, 0, 1, 0, 0.012, 0.8, 1.4, &test3Density, &test3Velocity,
			&test3Pressure, fvm::FluxScheme::force);
	fvm::Simulation test4(200, 0, 1, 0, 0.035, 0.8, 1.4, &test4Density, &test4Velocity,
			&test4Pressure, fvm::FluxScheme::force);
	fvm::Simulation test5(200, 0, 1, 0, 0.035, 0.8, 1.4, &test5Density, &test5Velocity,
			&test5Pressure, fvm::FluxScheme::force);

	std::ofstream output1("test1.dat");
	std::ofstream output2("test2.dat");
	std::ofstream output3("test3.dat");
	std::ofstream output4("test4.dat");
	std::ofstream output5("test5.dat");

	if (!output1) {
		std::cerr << "Error: could not create output file.\n";
		return 1;
	}

	runSimulation(test1, output1);
	runSimulation(test2, output2);
	runSimulation(test3, output3);
	runSimulation(test4, output4);
	runSimulation(test5, output5);

	return 0;
}
