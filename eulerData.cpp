#include "eulerData.hpp"

namespace fvm {
	// Operator overloads for QuantArray.
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
			auto rho = data_[i][dIndex];
			auto p = data_[i][pIndex];

			// Convert velocity to momentum.
			data_[i][vIndexX] *= rho;
			data_[i][vIndexY] *= rho;

			auto rhoVX = data_[i][vIndexX];
			auto rhoVY = data_[i][vIndexY];

			// Convert pressure to energy.
			data_[i][pIndex] = (p / (gamma_ - 1)) + 0.5 * ((rhoVX*rhoVX + rhoVY*rhoVY) / rho);
		}

		mode_ = EulerDataMode::conserved;
	}

	void EulerData::makePrimitive_() {
		int dIndex = static_cast<int>(ConservedQuant::density);
		int moIndexX = static_cast<int>(ConservedQuant::momentumX);
		int moIndexY = static_cast<int>(ConservedQuant::momentumY);
		int eIndex = static_cast<int>(ConservedQuant::energy);

		for (size_t i = 0; i < data_.size(); ++i) {
			auto d = data_[i][dIndex];
			auto rhoVX = data_[i][moIndexX];
			auto rhoVY = data_[i][moIndexY];
			auto e = data_[i][eIndex];

			// Convert energy to pressure.
			data_[i][eIndex] = (gamma_ - 1) * (e - 0.5 * ((rhoVX*rhoVX + rhoVY*rhoVY) / d));

			// Convert momentum to velocity.
			data_[i][moIndexX] /= d;
			data_[i][moIndexY] /= d;
		}

		mode_ = EulerDataMode::primitive;
	}
}
