#ifndef GHOSTSOLVER_FVMDATA_HPP
#define GHOSTSOLVER_FVMDATA_HPP

#include <array>
#include <cstddef>
#include <vector>

using std::size_t;

namespace fvm {
	/** Array that stores the physical quantities of a cell. These
	 * quantities might be primitive (such as density, velocity, and
	 * pressure), or conserved (such as density, momentum, and energy).
	 * @brief Measured physical quantities within a cell.
	 */
	typedef std::array<double, 3> QuantArray;
	// The cell vector stores the QuantArrays for every cell.
	/** Vector of cell values.
	 * @brief Vector of cell values.
	 */
	typedef std::vector<QuantArray> CellVector;

	// Operator overloads for QuantArray.
	// Vector addition and subtraction.
	QuantArray operator+(QuantArray a, QuantArray b);
	QuantArray operator-(QuantArray a, QuantArray b);
	// vector multiplication/division by scalar.
	QuantArray operator*(double a, QuantArray u);
	QuantArray operator*(QuantArray u, double a);
	QuantArray operator/(QuantArray u, double a);

	QuantArray makePrimQuants(const QuantArray& u, const double gamma);

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

	/** Holds the simulation data. This includes current cell values, and
	 * the type of quantities (e.g. primitive or conserved). Also provides
	 * helper functions to manipulate and retrieve the state of the data.
	 *
	 * @brief Holds the simulation data.
	 */
	class EulerData {
		public:
			// Constructor
			// Data is assumed to be in primitive form by default.
			EulerData(double gamma, EulerDataMode mode = EulerDataMode::primitive)
				// Initialiser list
				: gamma_(gamma), mode_(mode) {};

			// EulerData accessors
			// Getters
			QuantArray& operator[](size_t i) { return data_[i]; }
			size_t size() { return data_.size(); }
			CellVector& data() { return data_; }
			EulerDataMode mode() { return mode_; }

			// Setters
			void setQuantity(size_t i, QuantArray newValues) { data_[i] = newValues; }
			void setMode(EulerDataMode want);

		private:
			// Private member data
			double gamma_;
			CellVector data_;
			EulerDataMode mode_;
			void makeConserved_();
			void makePrimitive_();
	};
}

#endif
