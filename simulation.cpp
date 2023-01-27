#include <fstream>
#include <iostream>

#include "FvmSimulation.hpp"
double densityDist(double x);
double velocityDist(double x);
double pressureDist(double x);

int main(void) {
	fvm::Simulation test1(800, 0, 1, 0, 0.25, 0.8, 1.4, &densityDist, &velocityDist,
			&pressureDist, fvm::FluxScheme::force);
	std::ofstream output("euler.dat");


	if (!output) {
		std::cerr << "Error: could not create output file.\n";
		return 1;
	}

	test1.saveToFile(output);

	// Change variables to conserved for solving.
	test1.convertToConserved();
	int nFrame = 0;
	for (double t = test1.tStart(); t < test1.tEnd(); t = test1.tNow()) {
		test1.step();
		++nFrame;

		if (nFrame % 3 == 0) {
			// Change variables back to primitive for output.
			test1.convertToPrimitive();
			test1.saveToFile(output);
			test1.convertToConserved();
		}
	}

	return 0;
}

// Initial distribution of primitive variables
double densityDist(double x) {
	return x <= 0.5 ? 1 : 0.125;
}

double velocityDist(double x) {
	x = 0;
	return x;
}

double pressureDist(double x) {
	return x <= 0.5 ? 1 : 0.1;
}
