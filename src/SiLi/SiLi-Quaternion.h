#pragma once

#include "SiLi.h"

namespace SiLi {
	/*
	 * Quaternion
	 *
	 * internal representation SiLi::Matrix<4, 1>(w, x, y, z)
	 */
	template <typename T = DefaultType>
	class Quaternion {
	private:
		Matrix<4, 1, T> mValue;

	public:
		auto real() -> T& {
			return mValue(0);
		}
		auto imag() -> typename Matrix<3, 1, T>::View {
			return mValue.template view<3, 1>(1, 0);
		}

		auto real() const -> T const& {
			return mValue(0);
		}
		auto imag() const -> typename Matrix<3, 1, T>::CView {
			return mValue.template view<3, 1>(1, 0);
		}

	public:
		Quaternion()
			: mValue {{T(1.), T(0.), T(0.), T(0.)}}
		{}
		// _angle in radian
		Quaternion(T _angle, Matrix<3, 1, T> const& _axis)
		{
			using std::cos;
			using std::sin;
			real() = cos(_angle * T(0.5));
			imag() = sin(_angle * T(0.5)) * _axis.normalized();
		}

		Quaternion(Quaternion const& _q)
			: mValue {_q.mValue}
		{}
		Quaternion(Matrix<4, 1, T> const& _matrix)
			: mValue(_matrix)
		{}

		Quaternion(Matrix<3, 3, T> const& mat) {
			fromMatrix(mat);
		}

		/**
		 * computes quaternion for minimal rotation from v1 to v2
		 */
		Quaternion(Matrix<3, 1, T> v1, Matrix<3, 1, T> v2) {
			v1 = v1.normalized();
			v2 = v2.normalized();

			real() = T(1.) + v1.t() * v2;
			imag() = cross(v1, v2);

			auto n = norm();
			using std::abs;
			if (abs(n) < T(1e-9)) {
				real() = T(1.);
				imag() = T(0.);
			}
			*this = normalized();
		}
	private:
		void fromMatrix(Matrix<3, 3, T> const& mat) {
			auto trace = sum(mat.diag());
			using std::sqrt;
			if (trace > T(0.)) {
				auto s = sqrt(trace+ T(1.)) * T(2.);
				mValue(0) = T(0.25) * s;
				mValue(1) = (mat(2, 1) - mat(1, 2)) / s;
				mValue(2) = (mat(0, 2) - mat(2, 0)) / s;
				mValue(3) = (mat(1, 0) - mat(0, 1)) / s;
			} else if (mat(0, 0) > mat(1, 1) and mat(0, 0) > mat(2, 2)) {
				auto s = sqrt(1. + mat(0, 0) - mat(1, 1) - mat(2, 2)) * T(2.);
				mValue(0) = (mat(2, 1) - mat(1, 2)) / s;
				mValue(1) = T(0.25) * s;
				mValue(2) = (mat(0, 1) + mat(1, 0)) / s;
				mValue(3) = (mat(0, 2) + mat(2, 0)) / s;
			} else if (mat(1, 1) > mat(2, 2)) {
				auto s = sqrt(T(1.) + mat(1, 1) - mat(0, 0) - mat(2, 2)) * T(2.);
				mValue(0) = (mat(0, 2) - mat(2, 0)) / s;
				mValue(1) = (mat(0, 1) + mat(1, 0)) / s;
				mValue(2) = T(0.25) * s;
				mValue(3) = (mat(1, 2) + mat(2, 1)) / s;
			} else {
				auto s = sqrt(T(1.) + mat(2, 2) - mat(0, 0) - mat(1, 1)) * T(2.);
				mValue(0) = (mat(1, 0) - mat(0, 1)) / s;
				mValue(1) = (mat(0, 2) + mat(2, 0)) / s;
				mValue(2) = (mat(1, 2) + mat(2, 1)) / s;
				mValue(3) = T(0.25) * s;
			}
			*this = normalized();
		}

	public:
		auto mat() const -> SiLi::Matrix<3, 3, T> {
			auto const& q = mValue;
			auto R = SiLi::make_eye<3, 3, double>();
			R(0, 0) = 1 - 2 * (q(2)*q(2) + q(3)*q(3));
			R(0, 1) = -2*q(0)*q(3) + 2*q(1)*q(2);
			R(0, 2) = 2*q(0)*q(2) + 2*q(1)*q(3);

			R(1, 0) = 2*q(0)*q(3) + 2*q(1)*q(2);
			R(1, 1) = 1 - 2*(q(1)*q(1) + q(3)*q(3));
			R(1, 2) = -2*q(0)*q(1) + 2*q(2)*q(3);

			R(2, 0) = -2*q(0)*q(2) + 2*q(1)*q(3);
			R(2, 1) = 2*q(0)*q(1) + 2*q(2)*q(3);
			R(2, 2) = 1-2*(q(1)*q(1) + q(2)*q(2));
			return R;
		}

		auto operator()(int e) -> T& {
			return mValue(e);
		}

		auto operator()(int e) const -> T const& {
			return mValue(e);
		}

		auto operator+(Quaternion const& _q) const -> Quaternion {
			return Quaternion (mValue + _q.mValue);
		}
		auto operator+=(Quaternion const& _q) -> Quaternion& {
			mValue += _q.mValue;
			return *this;
		}
		auto operator-(Quaternion const& _q) const -> Quaternion {
			return {mValue - _q.mValue};
		}
		auto operator-=(Quaternion const& _q) -> Quaternion& {
			mValue -= _q.mValue;
			return *this;
		}
		auto operator-() const -> Quaternion {
			return {-mValue};
		}
		auto operator*(T scalar) const -> Quaternion {
			return {mValue*scalar};
		}
		auto operator*=(T scalar) -> Quaternion& {
			mValue *= scalar;
			return *this;
		}

		operator SiLi::Matrix<3, 3, T>() const {
			return mat();
		}

		auto norm() const -> T {
			return mValue.norm();
		}

		auto normalized() const -> Quaternion {
			return {mValue * (T(1.) / norm())};
		}

		auto conjugate() const -> Quaternion {
			Quaternion q;
			q.real() =  real();
			q.imag() = -imag();
			return q;
		}

		auto dot(Quaternion q) const -> T {
			return mValue.t() * q.mValue;
		}

		auto operator*(Quaternion const& q) const -> Quaternion {
			Quaternion r;
			r.real() = real() * q.real() - imag().t() * q.imag();
			r.imag() = real() * q.imag() + q.real() * imag() + cross(imag(), q.imag());
			return r;
		}
		auto operator*=(Quaternion const& q) -> Quaternion& {
			*this = *this * q;
			return *this;
		}

		auto rotate(SiLi::Matrix<3, 1, T> const& v) const -> SiLi::Matrix<3, 1, T> {
			Quaternion q;
			q.real() = 0.;
			q.imag() = v;
			return (*this * q * this->conjugate()).imag();
		}

		auto slerp(Quaternion<T> v1, T factor) -> Quaternion {
			auto v0 = normalized();
			v1 = v1.normalized();

			auto dot = v0.dot(v1);
			if (dot > 1 - 1e-9) {
				return (v0 + (v1 - v0) * factor).normalized();
			}

			if (dot < 0.) {
				v1 = -v1;
				dot = -dot;
			}

			using std::max;
			using std::min;
			dot = min(T(1.), max(T(-1), dot));
			using std::acos;
			using std::sin;
			using std::cos;
			auto theta_0 = acos(dot);
			auto theta   = theta_0 * factor;

			auto v2 = (v1 - v0*dot).normalized();
			auto r =  v0 * cos(theta) + v2 * sin(theta);
			return r.normalized();
		}
	};

/*
 * ostream
 */
template<typename T>
std::ostream& operator<< (std::ostream& stream, Quaternion<T> const& q) {
	stream << q(0) << " " << q(1) << " " << q(2) << " " << q(3);
	return stream;
}

}
