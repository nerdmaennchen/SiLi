#include <SiLi/SiLi-Quaternion.h>
#include <gtest/gtest.h>

#if __GNUC__ > 4

TEST(SiLi, quat_init) {
	SiLi::Quaternion<double> q1;
	EXPECT_EQ(q1(0), 1.);
	EXPECT_EQ(q1(1), 0.);
	EXPECT_EQ(q1(2), 0.);
	EXPECT_EQ(q1(3), 0.);

	EXPECT_EQ(q1.real(), 1.);
	EXPECT_EQ(q1.imag()(0), 0.);
	EXPECT_EQ(q1.imag()(1), 0.);
	EXPECT_EQ(q1.imag()(2), 0.);
}

TEST(SiLi, quat_init2) {
	SiLi::Quaternion<double> q1(M_PI*0.5, {{1., 0., 0.}});

	SiLi::Matrix<4, 1, double> e({1./std::sqrt(2), 1./std::sqrt(2.), 0., 0.});
	EXPECT_NEAR(q1(0), e(0), 1e-9);
	EXPECT_NEAR(q1(1), e(1), 1e-9);
	EXPECT_NEAR(q1(2), e(2), 1e-9);
	EXPECT_NEAR(q1(3), e(3), 1e-9);

	EXPECT_NEAR(q1.real(),   e(0), 1e-9);
	EXPECT_NEAR(q1.imag()(0), e(1), 1e-9);
	EXPECT_NEAR(q1.imag()(1), e(2), 1e-9);
	EXPECT_NEAR(q1.imag()(2), e(3), 1e-9);
}

TEST(SiLi, quat_mul1) {
	SiLi::Quaternion<double> q1(M_PI*0.5, {{ 1., 0., 0.}});
	SiLi::Quaternion<double> q2(M_PI*0.5, {{-1., 0., 0.}});

	auto q3 = q1 * q2;

	SiLi::Matrix<4, 1, double> e({1., 0., 0., 0.});
	EXPECT_NEAR(q3(0), e(0), 1e-9);
	EXPECT_NEAR(q3(1), e(1), 1e-9);
	EXPECT_NEAR(q3(2), e(2), 1e-9);
	EXPECT_NEAR(q3(3), e(3), 1e-9);
}

TEST(SiLi, quat_mul2) {
	SiLi::Quaternion<double> q1(M_PI*0.5, {{1., 0., 0.}});
	auto q2 = q1.conjugate();
	auto q3 = q1 * q2;

	SiLi::Matrix<4, 1, double> e({1., 0., 0., 0.});
	EXPECT_NEAR(q3(0), e(0), 1e-9);
	EXPECT_NEAR(q3(1), e(1), 1e-9);
	EXPECT_NEAR(q3(2), e(2), 1e-9);
	EXPECT_NEAR(q3(3), e(3), 1e-9);
}

TEST(SiLi, quat_slerp) {
	SiLi::Quaternion<double> q1(M_PI, {{1., 0., 0.}});
	SiLi::Quaternion<double> q2(0., {{1., 0., 0.}});
	auto q3 = q1.slerp(q2, 0.5);

	SiLi::Quaternion<double> e({M_PI*0.5, {{1., 0., 0.}}});
	EXPECT_NEAR(q3(0), e(0), 1e-9);
	EXPECT_NEAR(q3(1), e(1), 1e-9);
	EXPECT_NEAR(q3(2), e(2), 1e-9);
	EXPECT_NEAR(q3(3), e(3), 1e-9);
}

TEST(SiLi, quat_to_mat) {
	SiLi::Quaternion<double> q1(M_PI * 0.25, {{1., 0., 0.}});

	SiLi::Matrix<3, 3, double> e ({{1., 0., 0.},
	                              {0., 1./std::sqrt(2.), -1./std::sqrt(2.)},
	                              {0., 1./std::sqrt(2.),  1./std::sqrt(2.)}});
	auto M = q1.mat();
	for (auto r : {0, 1, 2}) {
		for (auto c : {0, 1, 2}) {
			EXPECT_NEAR(e(r, c), M(r, c), 1e-9);
		}
	}
}

TEST(SiLi, mat_to_quat) {
	SiLi::Matrix<3, 3, double> m ({{1., 0., 0.},
	                              {0., 1./std::sqrt(2.), -1./std::sqrt(2.)},
	                              {0., 1./std::sqrt(2.),  1./std::sqrt(2.)}});
	SiLi::Quaternion<double> q = m;
	SiLi::Quaternion<double> e(M_PI * 0.25, {{1., 0., 0.}});

	EXPECT_NEAR(e(0), q(0), 1e-9);
	EXPECT_NEAR(e(1), q(1), 1e-9);
	EXPECT_NEAR(e(2), q(2), 1e-9);
	EXPECT_NEAR(e(3), q(3), 1e-9);
}

TEST(SiLi, quat_minimalrotation) {
	SiLi::Matrix<3, 1, double> v1 ({1., 0., 0.});
	SiLi::Matrix<3, 1, double> v2 ({0.5, 0.5, 0.});

	SiLi::Quaternion<double> q (v1, v2);

	auto v3 = q.rotate(v1);

	EXPECT_NEAR(double(v2.t() * v3) * (1./ (v2.norm() * v3.norm())), 1., 1e-9);
}

#endif
