#ifndef GHOSTSOLVER_FVMDATA_HPP
#define GHOSTSOLVER_FVMDATA_HPP

#include <array>
#include <cstddef>
#include <ostream>
#include <vector>

using std::size_t;

namespace fvm {
	/** Array that stores the physical quantities of a cell. These
	 * quantities might be primitive (such as density, velocity, and
	 * pressure), or conserved (such as density, momentum, and energy).
	 * @brief Measured physical quantities within a cell.
	 */
	typedef std::array<double, 4> Cell;
	/** Vector of cells.
	 * @brief Vector of cell values.
	 */
	typedef std::vector<Cell> CellVector;
	typedef std::vector<CellVector> Grid;

	// Operator overloads for Cell.
	// Vector addition and subtraction.
	Cell operator+(Cell a, Cell b);
	Cell operator-(Cell a, Cell b);
	// Vector multiplication/division by scalar.
	Cell operator*(double a, Cell u);
	Cell operator*(Cell u, double a);
	Cell operator/(Cell u, double a);

	enum class PrimitiveQuant {
		density = 0,
		velocityX,
		velocityY,
		pressure
	};

	// Primitive indices.
	constexpr auto dIndex = static_cast<int>(PrimitiveQuant::density);
	constexpr auto vIndexX = static_cast<int>(PrimitiveQuant::velocityX);
	constexpr auto vIndexY = static_cast<int>(PrimitiveQuant::velocityY);
	constexpr auto pIndex = static_cast<int>(PrimitiveQuant::pressure);

	enum class ConservedQuant {
		density = dIndex,
		momentumX = vIndexX,
		momentumY = vIndexY,
		energy = pIndex
	};

	// Conservative indices.
	constexpr auto moIndexX = static_cast<int>(ConservedQuant::momentumX);
	constexpr auto moIndexY = static_cast<int>(ConservedQuant::momentumY);
	constexpr auto eIndex = static_cast<int>(ConservedQuant::energy);

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
			explicit EulerData(double gamma, EulerDataMode mode = EulerDataMode::primitive)
				// Initialiser list
				: gamma_(gamma), mode_(mode) {};

			// EulerData accessors
			// Getters
			size_t xSize() { return data_.size(); }
			size_t ySize() { return data_[0].size(); }
			Grid& data() { return data_; }
			EulerDataMode mode() { return mode_; }

			// Setters
			void setCell(size_t i, size_t j, Cell newValues) { data_[i][j] = newValues; }
			void setMode(EulerDataMode want);

			// Operator overloads.
			CellVector& operator[](size_t i) { return data_[i]; }

		private:
			// Private member data
			double gamma_;
			Grid data_;
			EulerDataMode mode_;
			void makeConserved_();
			void makePrimitive_();
	};
}

#endif
