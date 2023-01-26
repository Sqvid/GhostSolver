#include <fstream>
#include <iostream>

#include "FvmSimulation.hpp"
double densityDist(double x);
double velocityDist(double x);
double pressureDist(double x);
void outputSim(std::ofstream& output, fvm::Simulation& sim);

int main(void) {
	fvm::Simulation test1(800, 0, 1, 0, 0.25, 0.8, 1.4, &densityDist, &velocityDist,
			&pressureDist, fvm::FluxScheme::force);
	std::ofstream output("euler.dat");


	if (!output) {
		std::cerr << "Error: could not create output file.\n";
		return 1;
	}

	outputSim(output, test1);

	// Change variables to conserved for solving.
	test1.convertToConserved();
	int nFrame = 0;
	for (double t = test1.tStart(); t < test1.tEnd(); t = test1.tNow()) {
		test1.step();
		++nFrame;

		if (nFrame % 3 == 0) {
			// Change variables back to primitive for output.
			test1.convertToPrimitive();
			outputSim(output, test1);
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

void outputSim(std::ofstream& output, fvm::Simulation& sim) {
	for (size_t i = 0; i < sim.nCells(); ++i) {
		double x = sim.xStart() + i * sim.dx();

		// Indices of primitive quantities.
		int d = static_cast<int>(fvm::PrimitiveQuant::density);
		int v = static_cast<int>(fvm::PrimitiveQuant::velocity);
		int p = static_cast<int>(fvm::PrimitiveQuant::pressure);

		// Values of primitive quantities.
		double density = sim.getQuantity(i, d);
		double velocity = sim.getQuantity(i, v);
		double pressure = sim.getQuantity(i, p);

		output << x << " " << density << " " << velocity << " " << pressure << "\n";
	}

	// Delimit Gnuplot code block with two blank lines.
	output << "\n\n";
}
