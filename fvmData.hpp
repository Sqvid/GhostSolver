#ifndef GHOSTSOLVER_FVMDATA_HPP
#define GHOSTSOLVER_FVMDATA_HPP

#include <array>
#include <vector>

using std::size_t;

namespace fvm {
	// A vector of arrays of length three (triplet).
	typedef std::array<double, 3> QuantArray;
	typedef std::vector<QuantArray> TripletVector;

	// Operator overloads for QuantArray.
	// Vector addition and subtraction.
	QuantArray operator+(QuantArray a, QuantArray b);
	QuantArray operator-(QuantArray a, QuantArray b);
	// vector multiplication/division by scalar. 
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

	// The types of variables EulerData can hold.
	enum class EulerDataMode {
		primitive,
		conserved
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
			double gamma_;
			TripletVector data_;
			EulerDataMode mode_;
	};
}

#endif
