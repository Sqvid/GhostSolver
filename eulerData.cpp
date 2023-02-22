#include "eulerData.hpp"

namespace fvm {
	// Operator overloads for Cell.
	Cell operator+(Cell a, Cell b) {
		Cell ans;

		for (size_t i = 0; i < a.size(); ++i) {
			ans[i] = a[i] + b[i];
		}

		return ans;
	}

	Cell operator-(Cell a, Cell b) {
		Cell ans;

		for (size_t i = 0; i < a.size(); ++i) {
			ans[i] = a[i] - b[i];
		}

		return ans;
	}

	Cell operator*(double a, Cell u) {
		Cell ans;

		for (size_t i = 0; i < u.size(); ++i) {
			ans[i] = a * u[i];
		}

		return ans;
	}

	Cell operator*(Cell u, double a) {
		Cell ans;

		for (size_t i = 0; i < u.size(); ++i) {
			ans[i] = a * u[i];
		}

		return ans;
	}

	Cell operator/(Cell u, double a) {
		Cell ans;

		for (size_t i = 0; i < u.size(); ++i) {
			ans[i] = u[i] / a;
		}

		return ans;
	}

	std::ostream& operator<<(std::ostream& stream, Cell cell) {
		for (size_t i = 0; i < cell.size(); ++i) {
			stream << cell[i];
		}

		return stream;
	}

	// Put the data in the correct mode.
	void EulerData::setMode(EulerDataMode want) {
		// The mode is already correctly set; return early.
		if (mode_ == want) {
			return;
		}

		switch (want) {
			case EulerDataMode::primitive:
				makePrimitive_();
				break;
			case EulerDataMode::conserved:
				makeConserved_();
				break;
		}
	}

	void EulerData::makeConserved_() {
		int dIndex = static_cast<int>(PrimitiveQuant::density);
		int vIndexX = static_cast<int>(PrimitiveQuant::velocityX);
		int vIndexY = static_cast<int>(PrimitiveQuant::velocityY);
		int pIndex = static_cast<int>(PrimitiveQuant::pressure);

		for (size_t i = 0; i < data_.size(); ++i) {
			for (size_t j = 0; j < data_.size(); ++j) {
				auto rho = data_[i][j][dIndex];
				auto p = data_[i][j][pIndex];

				// Convert velocity to momentum.
				data_[i][j][vIndexX] *= rho;
				data_[i][j][vIndexY] *= rho;

				auto rhoVX = data_[i][j][vIndexX];
				auto rhoVY = data_[i][j][vIndexY];

				// Convert pressure to energy.
				data_[i][j][pIndex] = (p / (gamma_ - 1)) + 0.5 * ((rhoVX*rhoVX + rhoVY*rhoVY) / rho);
			}
		}

		mode_ = EulerDataMode::conserved;
	}

	void EulerData::makePrimitive_() {
		int dIndex = static_cast<int>(ConservedQuant::density);
		int moIndexX = static_cast<int>(ConservedQuant::momentumX);
		int moIndexY = static_cast<int>(ConservedQuant::momentumY);
		int eIndex = static_cast<int>(ConservedQuant::energy);

		for (size_t i = 0; i < this->xSize(); ++i) {
			for (size_t j = 0; j < this->ySize(); ++j) {
				auto rho = data_[i][j][dIndex];
				auto rhoVX = data_[i][j][moIndexX];
				auto rhoVY = data_[i][j][moIndexY];
				auto e = data_[i][j][eIndex];

				// Convert energy to pressure.
				data_[i][j][eIndex] = (gamma_ - 1) * (e - 0.5 * ((rhoVX*rhoVX + rhoVY*rhoVY) / rho));

				// Convert momentum to velocity.
				data_[i][j][moIndexX] /= rho;
				data_[i][j][moIndexY] /= rho;
			}
		}

		mode_ = EulerDataMode::primitive;
	}
}
