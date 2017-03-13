#pragma once

#include "SiLi.h"

namespace SiLi {
	template <typename T = DefaultType>
	class Quaternion {
	private:
		Matrix<4, 1, T> mValue;

	public:
		auto real() -> T& {
			return mValue(0);
		}
		auto img() -> typename Matrix<3, 1, T>::View {
			return mValue.template view<3, 1>(1, 0);
		}

		auto real() const -> T const& {
			return mValue(0);
		}
		auto img() const -> typename Matrix<3, 1, T>::CView {
			return mValue.template view<3, 1>(1, 0);
		}

	public:
		Quaternion()
			: mValue {{1., 0., 0., 0.}}
		{}
		// _angle in radian
		Quaternion(T _angle, Matrix<3, 1, T> const& _axis)
		{
			using std::cos;
			using std::sin;
			real() = 0.5 * cos(_angle);
			img()  = 0.5 * sin(_angle) * _axis;
		}

		Quaternion(T _angle, T _x, T _y, T _z)
		{
			using std::cos;
			using std::sin;
			real() = cos(_angle * 0.5);
			img()  = sin(_angle * 0.5) * Matrix<3, 1, T>{{_x, _y, _z}};
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
		Quaternion(Matrix<4, 4, T> const& mat) {
			fromMatrix(mat.template view<3, 3>(0, 0));
		}
		/**
		 * computes quaternion for minimal rotation from v1 to v2
		 */
		Quaternion(Matrix<3, 1, T> v1, Matrix<3, 1, T> v2) {
			v1 = v1.normalize();
			v2 = v2.normalize();

			real() = 1. + v1.t() * v2;
			img()  = cross(v1, v2);

			auto n = norm();
			if (std::abs(n) < 1e-9) {
				real() = 1.;
				img() = 0;
			}
			*this = normalize();
		}
	private:
		void fromMatrix(Matrix<3, 3, T> const& mat) {
			auto trace = sum(mat.diag());
			using std::sqrt;
			if (trace > 0) {
				auto s = sqrt(trace+1.0) * 2.;
				mValue(0) = 0.25 * s;
				mValue(1) = (mat(2, 1) - mat(1, 2)) / s;
				mValue(2) = (mat(0, 2) - mat(2, 0)) / s;
				mValue(3) = (mat(1, 0) - mat(0, 1)) / s;
			} else if (mat(0, 0) > mat(1, 1) and mat(0, 0) > mat(2, 2)) {
				auto s = sqrt(1. + mat(0, 0) - mat(1, 1) - mat(2, 2)) * 2.;
				mValue(0) = (mat(2, 1) - mat(1, 2)) / s;
				mValue(1) = 0.25 * s;
				mValue(2) = (mat(0, 1) + mat(1, 0)) / s;
				mValue(3) = (mat(0, 2) + mat(2, 0)) / s;
			} else if (mat(1, 1) > mat(2, 2)) {
				auto s = sqrt(1. + mat(1, 1) - mat(0, 0) - mat(2, 2)) * 2.;
				mValue(0) = (mat(0, 2) - mat(2, 0)) / s;
				mValue(1) = (mat(0, 1) + mat(1, 0)) / s;
				mValue(2) = 0.25 * s;
				mValue(3) = (mat(1, 2) + mat(2, 1)) / s;
			} else {
				auto s = sqrt(1. + mat(2, 2) - mat(0, 0) - mat(1, 1)) * 2.;
				mValue(0) = (mat(1, 0) - mat(0, 1)) / s;
				mValue(1) = (mat(0, 2) + mat(2, 0)) / s;
				mValue(2) = (mat(1, 2) + mat(2, 1)) / s;
				mValue(3) = 0.25 * s;
			}
			*this = normalize();
		}

	public:
		auto mat() const -> SiLi::Matrix<4, 4, T> {
			auto const& q = mValue;
			auto R = SiLi::make_eye<4, 4, double>();
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
			return {mValue + _q.mValue};
		}
		auto operator-=(Quaternion const& _q) -> Quaternion& {
			mValue += _q.mValue;
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
		operator SiLi::Matrix<4, 4,T>() const {
			return mat();
		}
		operator SiLi::Matrix<3, 3, T>() const {
			return mat().view<3, 3>(0, 0);
		}
		auto norm() const -> T {
			return mValue.norm();
		}
		auto normalize() const -> Quaternion {
			return {mValue * (1. / mValue.norm())};
		}

		auto conjugate() const -> Quaternion {
			Quaternion q;
			q.real() = real();
			q.img()  = -img();
			return q;
		}
		auto dot(Quaternion q) const -> T {
			return mValue.t() * q.mValue;
		}

		auto operator*(Quaternion const& q) const -> Quaternion {
			Quaternion r;
			r.real() = real() * q.real() - img().t() * q.img();
			r.img()  = real() * q.img() + q.real() * img() + cross(img(), q.img());
			return r;
		}
		auto operator*=(Quaternion const& q) -> Quaternion& {
			*this = *this * q;
			return *this;
		}

		auto rotate(SiLi::Matrix<3, 1, T> const& v) const -> SiLi::Matrix<3, 1, T> {
			Quaternion q;
			q.real() = 0.;
			q.img()  = v;
			return (*this * q * this->conjugate()).img();
		}



		auto slerp(Quaternion<T> v1, T factor) -> Quaternion {
			auto v0 = normalize();
			v1 = v1.normalize();

			auto dot = v0.dot(v1);
			if (dot > 1 - 1e-9) {
				return (v0 + (v1 - v0) * factor).normalize();
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

			auto v2 = (v1 - v0*dot).normalize();
			auto r =  v0 * cos(theta) + v2 * sin(theta);
			return r.normalize();
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
