#include <SiLi/SiLi-Quaternion.h>
#include <catch2/catch.hpp>

TEST_CASE("quat_init", "") {
	SiLi::Quaternion<double> q1;
	CHECK(q1(0) == 1.);
	CHECK(q1(1) == 0.);
	CHECK(q1(2) == 0.);
	CHECK(q1(3) == 0.);

	CHECK(q1.real() == 1.);
	CHECK(q1.imag()(0) == 0.);
	CHECK(q1.imag()(1) == 0.);
	CHECK(q1.imag()(2) == 0.);
}

TEST_CASE("quat_init2", "") {
	SiLi::Quaternion<double> q1(M_PI*0.5, {{1., 0., 0.}});

	SiLi::Matrix<4, 1, double> e({1./std::sqrt(2), 1./std::sqrt(2.), 0., 0.});
	CHECK(q1(0) == Approx(e(0)).margin(1e-9));
	CHECK(q1(1) == Approx(e(1)).margin(1e-9));
	CHECK(q1(2) == Approx(e(2)).margin(1e-9));
	CHECK(q1(3) == Approx(e(3)).margin(1e-9));

	CHECK(q1.real() == Approx(  e(0)).margin(1e-9));
	CHECK(q1.imag()(0) == Approx(e(1)).margin(1e-9));
	CHECK(q1.imag()(1) == Approx(e(2)).margin(1e-9));
	CHECK(q1.imag()(2) == Approx(e(3)).margin(1e-9));
}

TEST_CASE("quat_mul1", "") {
	SiLi::Quaternion<double> q1(M_PI*0.5, {{ 1., 0., 0.}});
	SiLi::Quaternion<double> q2(M_PI*0.5, {{-1., 0., 0.}});

	auto q3 = q1 * q2;

	SiLi::Matrix<4, 1, double> e({1., 0., 0., 0.});
	CHECK(q3(0) == Approx(e(0)).margin(1e-9));
	CHECK(q3(1) == Approx(e(1)).margin(1e-9));
	CHECK(q3(2) == Approx(e(2)).margin(1e-9));
	CHECK(q3(3) == Approx(e(3)).margin(1e-9));
}

TEST_CASE("quat_mul2", "") {
	SiLi::Quaternion<double> q1(M_PI*0.5, {{1., 0., 0.}});
	auto q2 = q1.conjugate();
	auto q3 = q1 * q2;

	SiLi::Matrix<4, 1, double> e({1., 0., 0., 0.});
	CHECK(q3(0) == Approx(e(0)).margin(1e-9));
	CHECK(q3(1) == Approx(e(1)).margin(1e-9));
	CHECK(q3(2) == Approx(e(2)).margin(1e-9));
	CHECK(q3(3) == Approx(e(3)).margin(1e-9));
}

TEST_CASE("quat_slerp", "") {
	SiLi::Quaternion<double> q1(M_PI, {{1., 0., 0.}});
	SiLi::Quaternion<double> q2(0., {{1., 0., 0.}});
	auto q3 = q1.slerp(q2, 0.5);

	SiLi::Quaternion<double> e({M_PI*0.5, {{1., 0., 0.}}});
	CHECK(q3(0) == Approx(e(0)).margin(1e-9));
	CHECK(q3(1) == Approx(e(1)).margin(1e-9));
	CHECK(q3(2) == Approx(e(2)).margin(1e-9));
	CHECK(q3(3) == Approx(e(3)).margin(1e-9));
}

TEST_CASE("quat_to_mat", "") {
	SiLi::Quaternion<double> q1(M_PI * 0.25, {{1., 0., 0.}});

	SiLi::Matrix<3, 3, double> e ({{1., 0., 0.},
	                              {0., 1./std::sqrt(2.), -1./std::sqrt(2.)},
	                              {0., 1./std::sqrt(2.),  1./std::sqrt(2.)}});
	auto M = q1.mat();
	for (auto r : {0, 1, 2}) {
		for (auto c : {0, 1, 2}) {
			CHECK(e(r, c) == Approx(M(r, c)).margin(1e-9));
		}
	}
}

TEST_CASE("mat_to_quat", "") {
	SiLi::Matrix<3, 3, double> m ({{1., 0., 0.},
	                              {0., 1./std::sqrt(2.), -1./std::sqrt(2.)},
	                              {0., 1./std::sqrt(2.),  1./std::sqrt(2.)}});
	SiLi::Quaternion<double> q = m;
	SiLi::Quaternion<double> e(M_PI * 0.25, {{1., 0., 0.}});

	CHECK(e(0) == Approx(q(0)).margin(1e-9));
	CHECK(e(1) == Approx(q(1)).margin(1e-9));
	CHECK(e(2) == Approx(q(2)).margin(1e-9));
	CHECK(e(3) == Approx(q(3)).margin(1e-9));
}

TEST_CASE("quat_minimalrotation", "") {
	SiLi::Matrix<3, 1, double> v1 ({1., 0., 0.});
	SiLi::Matrix<3, 1, double> v2 ({0.5, 0.5, 0.});

	SiLi::Quaternion<double> q (v1, v2);

	auto v3 = q.rotate(v1);

	CHECK(double(v2.t() * v3) * (1./ (v2.norm() * v3.norm())) == Approx(1.).margin(1e-9));
}
