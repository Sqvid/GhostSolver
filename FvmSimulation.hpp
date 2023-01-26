#ifndef PRACTICAL4_FVMSIMULATION_HPP
#define PRACTICAL4_FVMSIMULATION_HPP

#include <array>
#include <cstddef>
#include <functional>
#include <vector>

namespace fvm {
	// A vector of arrays of length three (triplet)
	typedef std::vector<std::array<double, 3>> TripletVector;

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
		laxFriedrich,
		richtmyer,
		force,
		godunov
	};

	// TODO; Move forward declaration.
	class Simulation;

	class EulerData {
		public:
			// Constructor
			// Data is assumed to be in primitive form by default.
			EulerData(EulerDataMode mode = EulerDataMode::primitive)
				// Initialiser list
				: mode_(mode) {};

			void convertToConserved(double gamma);
			void convertToPrimitive(double gamma);

			// EulerData accessors
			// Getters
			size_t size() { return data_.size(); }
			TripletVector& data() { return data_; }
			EulerDataMode mode() { return mode_; }
			std::array<double, 3>& front() { return data_.front(); }
			std::array<double, 3>& back() { return data_.back(); }
			double& velocity(size_t i) { return data_[i][static_cast<int>(PrimitiveQuant::velocity)]; }
			double& density(size_t i) { return data_[i][static_cast<int>(PrimitiveQuant::density)]; }
			double& pressure(size_t i) { return data_[i][static_cast<int>(PrimitiveQuant::pressure)]; }

		private:
			// Private member data
			TripletVector data_;
			EulerDataMode mode_;

			// Private member functions

			// Friend classes and functions
			friend class Simulation;
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
			EulerData& data() { return eulerData_; }
			//Setters
			void convertToConserved() {eulerData_.convertToConserved(gamma_);}
			void convertToPrimitive() {eulerData_.convertToPrimitive(gamma_);}

			// Public member functions
			void step();
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
			double lfFlux_(size_t i, ConservedQuant q);
			double richtmyerFlux_(size_t i, ConservedQuant q);
			double calcFlux_(size_t i, ConservedQuant q);
			double fluxExpr_(size_t i, ConservedQuant q);
	};
}

#endif
