#include "eulerData.hpp"

namespace fvm {
	// Operator overloads for QuantArray.
	QuantArray operator+(QuantArray a, QuantArray b) {
		QuantArray ans;

		for (size_t i = 0; i < a.size(); ++i) {
			ans[i] = a[i] + b[i];
		}

		return ans;
	}

	QuantArray operator-(QuantArray a, QuantArray b) {
		QuantArray ans;

		for (size_t i = 0; i < a.size(); ++i) {
			ans[i] = a[i] - b[i];
		}

		return ans;
	}

	QuantArray operator*(double a, QuantArray u) {
		QuantArray ans;

		for (size_t i = 0; i < u.size(); ++i) {
			ans[i] = a * u[i];
		}

		return ans;
	}

	QuantArray operator*(QuantArray u, double a) {
		QuantArray ans;

		for (size_t i = 0; i < u.size(); ++i) {
			ans[i] = a * u[i];
		}

		return ans;
	}

	QuantArray operator/(QuantArray u, double a) {
		QuantArray ans;

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
				makePrimitive();
				break;
			case EulerDataMode::conserved:
				makeConserved();
				break;
		}
	}

	void EulerData::makeConserved() {
		// Already in conserved form. Return early.
		if (mode_ == EulerDataMode::conserved) {
			return;
		}

		int dIndex = static_cast<int>(PrimitiveQuant::density);
		int vIndex = static_cast<int>(PrimitiveQuant::velocity);
		int pIndex = static_cast<int>(PrimitiveQuant::pressure);

		for (size_t i = 0; i < data_.size(); ++i) {
			// Convert velocity to momentum.
			data_[i][vIndex] *= data_[i][dIndex];

			// Convert pressure to energy.
			data_[i][pIndex] = (data_[i][pIndex] / (gamma_ - 1))
				+ 0.5 * ((data_[i][vIndex] * data_[i][vIndex]) / data_[i][dIndex]);
		}

		mode_ = EulerDataMode::conserved;
	}

	void EulerData::makePrimitive() {
		// Already in primitive form. Return early.
		if (mode_ == EulerDataMode::primitive) {
			return;
		}

		int dIndex = static_cast<size_t>(ConservedQuant::density);
		int moIndex = static_cast<int>(ConservedQuant::momentum);
		int eIndex = static_cast<int>(ConservedQuant::energy);

		for (size_t i = 0; i < data_.size(); ++i) {
			// Convert energy to pressure.
			data_[i][eIndex] = (gamma_ - 1) * (data_[i][eIndex]
					- 0.5 * ((data_[i][moIndex] * data_[i][moIndex]) / data_[i][dIndex]));

			// Convert momentum to velocity.
			data_[i][moIndex] /= data_[i][dIndex];
		}

		mode_ = EulerDataMode::primitive;
	}
}
