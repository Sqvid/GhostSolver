#include <fstream>
#include <iostream>

#include "FvmSimulation.hpp"

double test1Density(double x);
double test1Velocity(double x);
double test1Pressure(double x);
double test2Density(double x);
double test2Velocity(double x);
double test2Pressure(double x);
double test3Density(double x);
double test3Velocity(double x);
double test3Pressure(double x);
double test4Density(double x);
double test4Velocity(double x);
double test4Pressure(double x);
double test5Density(double x);
double test5Velocity(double x);
double test5Pressure(double x);

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
			&test1Pressure, fvm::FluxScheme::force);
	fvm::Simulation test2(200, 0, 1, 0, 0.15, 0.8, 1.4, &test2Density, &test2Velocity,
			&test2Pressure, fvm::FluxScheme::force);
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
	//test1.saveToFile(output1);

	//// Change variables to conserved for solving.
	//test1.convertToConserved();
	//int nFrame = 0;
	//for (double t = test1.tStart(); t < test1.tEnd(); t = test1.tNow()) {
	//	test1.step();
	//	++nFrame;

	//	if (nFrame % 3 == 0) {
	//		// Change variables back to primitive for output.
	//		test1.convertToPrimitive();
	//		test1.saveToFile(output1);
	//		test1.convertToConserved();
	//	}
	//}
	//runSimulation(test2, output2);
	//runSimulation(test3, output3);
	//runSimulation(test4, output4);
	//runSimulation(test5, output5);

	return 0;
}

// Initial distribution of primitive variables for Toro's tests.
//
//        rho_L   v_L     p_L     rho_R    v_R     p_R     time
// _____________________________________________________________
// 1 |    1       0       1       .125     0       .1      .25
// 2 |    1       -2      .4      1        2       .4      .15
// 3 |    1       0       1000    1        0       .01     .012
// 4 |    1       0       .01     1        0       100     .035
// 5 |    5.99924 19.5975 460.894 5.99242 -6.19633 46.0950 .035
//
// Each test is run on a x-domain [0, 1].
// L subscript for x <= 0.5, R subscript for x > 0.5
// gamma = 1.4 for all cases.
//
// Source for test initail conditions:
// Riemann Solvers and Numerical Methods for Fluid Dynamics, E.F. Toro,
// ISBN : 978-3-662-03492-7

double test1Density(double x) {
	return x <= 0.5 ? 1 : 0.125;
}

double test1Velocity(double x) {
	x = 0;
	return x;
}

double test1Pressure(double x) {
	return x <= 0.5 ? 1 : 0.1;
}

double test2Density(double x) {
	x = 1;
	return x;
}

double test2Velocity(double x) {
	return x <= 0.5 ? -2 : 2;
}

double test2Pressure(double x) {
	x = 0.4;
	return x;
}

double test3Density(double x) {
	x = 1;
	return x;
}

double test3Velocity(double x) {
	x = 0;
	return x;
}

double test3Pressure(double x) {
	return x <= 0.5 ? 1000 : 0.01;
}

double test4Density(double x) {
	x = 1;
	return x;
}

double test4Velocity(double x) {
	x = 0;
	return x;
}

double test4Pressure(double x) {
	return x <= 0.5 ? 0.01 : 100;
}

double test5Density(double x) {
	return x <= 0.5 ? 5.99924 : 5.99242;
}

double test5Velocity(double x) {
	return x <= 0.5 ? 19.5975 : -6.19633;
}

double test5Pressure(double x) {
	return x <= 0.5 ? 460.894 : 46.0950;
}
