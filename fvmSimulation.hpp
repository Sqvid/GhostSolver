#ifndef GHOSTSOLVER_FVMSIMULATION_HPP
#define GHOSTSOLVER_FVMSIMULATION_HPP

#include <array>
#include <functional>
#include <vector>
#include <fstream>

using std::size_t;

namespace fvm {
	// A vector of arrays of length three (triplet).
	typedef std::array<double, 3> QuantArray;
	typedef std::vector<QuantArray> TripletVector;

	// Operator overloads for QuantArray.
	QuantArray operator+(QuantArray a, QuantArray b);
	QuantArray operator-(QuantArray a, QuantArray b);
	QuantArray operator*(double a, QuantArray u);
	QuantArray operator*(QuantArray u, double a);
	QuantArray operator/(QuantArray u, double a);

	enum class PrimitiveQuant {
		density = 0,
		velocity,
		pressure
	};

	enum class ConservedQuant {
		density = 0,
		momentum,
		energy
	};

	// The types of variables EulerData can hold
	enum class EulerDataMode {
		primitive,
		conserved
	};

	// The type of flux scheme
	enum class FluxScheme {
		laxFriedrichs,
		richtmyer,
		force,
	};

	class EulerData {
		public:
			// Constructor
			// Data is assumed to be in primitive form by default.
			EulerData(double gamma, EulerDataMode mode = EulerDataMode::primitive)
				// Initialiser list
				: gamma_(gamma), mode_(mode) {};

			void makeConserved();
			void makePrimitive();

			// EulerData accessors
			// Getters
			size_t size() { return data_.size(); }
			TripletVector& data() { return data_; }
			EulerDataMode mode() { return mode_; }
			QuantArray getQuantity(size_t i) { return data_[i]; }

			// Setters
			void setQuantity(size_t i, QuantArray newValues) { data_[i] = newValues; }

		private:
			// Private member data
			TripletVector data_;
			double gamma_;
			EulerDataMode mode_;
	};

	// The simulation contains the parameters of the simulation as well as the
	// associated scheme, and flux expression being used.
	class Simulation {
		public:
			// Constructor
			Simulation(unsigned int nCells, double xStart, double xEnd,
					double tStart, double tEnd, double cfl,
					double gamma,
					std::function<double (double)> densityDist,
					std::function<double (double)> velocityDist,
					std::function<double (double)> pressureDist,
					FluxScheme fluxScheme);

			// Accessors
			// Getters
			unsigned int nCells() { return nCells_; }
			double xStart() { return xStart_; }
			double xEnd() { return xEnd_; }
			double tStart() { return tStart_; }
			double tEnd() { return tEnd_; }
			double tNow() { return tNow_; };
			double cfl() { return cfl_; }
			double dx() { return dx_; }
			double dt() { return dt_; }
			double gamma() { return gamma_; }
			EulerData data() { return eulerData_; }

			// Wrapper around EulerData qunatity getter.
			QuantArray getQuantity(size_t i) { return eulerData_.getQuantity(i); }

			// Public member functions
			void step();
			void saveToFile(std::ofstream& output);

		private:
			// Private member data
			unsigned int nCells_;
			double xStart_;
			double xEnd_;
			double tStart_;
			double tEnd_;
			double tNow_;
			double cfl_;
			double dx_;
			double dt_;
			double gamma_;
			FluxScheme fluxScheme_;
			std::function<double (double)> densityDist_;
			std::function<double (double)> velocityDist_;
			std::function<double (double)> pressureDist_;
			EulerData eulerData_;
			TripletVector flux_;

			// Private member functions
			double calcTimeStep_();
			QuantArray lfFlux_(size_t i);
			QuantArray richtmyerFlux_(size_t i);
			QuantArray calcFlux_(size_t i);
			QuantArray fluxExpr_(QuantArray u);
			// Wrappers around EulerData conversion functions.
			void makeConserved_() { eulerData_.makeConserved(); }
			void makePrimitive_() { eulerData_.makePrimitive(); }

	};
}

#endif
