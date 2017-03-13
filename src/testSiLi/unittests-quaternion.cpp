#include <SiLi/SiLi-Quaternion.h>
#include <gtest/gtest.h>

TEST(SiLi, quat_init) {
	SiLi::Quaternion<double> q1;
	EXPECT_EQ(q1(0), 1.);
	EXPECT_EQ(q1(1), 0.);
	EXPECT_EQ(q1(2), 0.);
	EXPECT_EQ(q1(3), 0.);

	EXPECT_EQ(q1.real(), 1.);
	EXPECT_EQ(q1.img()(0), 0.);
	EXPECT_EQ(q1.img()(1), 0.);
	EXPECT_EQ(q1.img()(2), 0.);
}

TEST(SiLi, quat_init2) {
	SiLi::Quaternion<double> q1(M_PI*0.5, 1., 0., 0.);

	SiLi::Matrix<4, 1, double> e({1./std::sqrt(2), 1./std::sqrt(2.), 0., 0.});
	EXPECT_NEAR(q1(0), e(0), 1e-9);
	EXPECT_NEAR(q1(1), e(1), 1e-9);
	EXPECT_NEAR(q1(2), e(2), 1e-9);
	EXPECT_NEAR(q1(3), e(3), 1e-9);

	EXPECT_NEAR(q1.real(),   e(0), 1e-9);
	EXPECT_NEAR(q1.img()(0), e(1), 1e-9);
	EXPECT_NEAR(q1.img()(1), e(2), 1e-9);
	EXPECT_NEAR(q1.img()(2), e(3), 1e-9);
}

TEST(SiLi, quat_mul1) {
	SiLi::Quaternion<double> q1(M_PI*0.5, 1., 0., 0.);
	SiLi::Quaternion<double> q2(M_PI*0.5, -1., 0., 0.);

	auto q3 = q1 * q2;

	SiLi::Matrix<4, 1, double> e({1., 0., 0., 0.});
	EXPECT_NEAR(q3(0), e(0), 1e-9);
	EXPECT_NEAR(q3(1), e(1), 1e-9);
	EXPECT_NEAR(q3(2), e(2), 1e-9);
	EXPECT_NEAR(q3(3), e(3), 1e-9);
}

TEST(SiLi, quat_mul2) {
	SiLi::Quaternion<double> q1(M_PI*0.5, 1., 0., 0.);
	auto q2 = q1.conjugate();
	auto q3 = q1 * q2;

	SiLi::Matrix<4, 1, double> e({1., 0., 0., 0.});
	EXPECT_NEAR(q3(0), e(0), 1e-9);
	EXPECT_NEAR(q3(1), e(1), 1e-9);
	EXPECT_NEAR(q3(2), e(2), 1e-9);
	EXPECT_NEAR(q3(3), e(3), 1e-9);
}

TEST(SiLi, quat_mul3) {
	SiLi::Quaternion<double> q1(M_PI, 1., 0., 0.);
	SiLi::Quaternion<double> q2(1., 0., 0., 0.);
	auto q3 = q1.slerp(q2, 0.5);

	SiLi::Quaternion<double> e({M_PI*0.5, 1., 0., 0.});
	EXPECT_NEAR(q3(0), e(0), 1e-9);
	EXPECT_NEAR(q3(1), e(1), 1e-9);
	EXPECT_NEAR(q3(2), e(2), 1e-9);
	EXPECT_NEAR(q3(3), e(3), 1e-9);
}

TEST(SiLi, quat_to_mat) {
	SiLi::Quaternion<double> q1(M_PI * 0.25, 1., 0., 0.);

	SiLi::Matrix<4, 4, double> e ({{1., 0., 0., 0.},
	                              {0., 1./std::sqrt(2.), -1./std::sqrt(2.), 0.},
	                              {0., 1./std::sqrt(2.),  1./std::sqrt(2.), 0.},
	                              {0., 0., 0., 1.}});
	auto M = q1.mat();
	for (auto r : {0, 1, 2, 3}) {
		for (auto c : {0, 1, 2, 3}) {
			EXPECT_NEAR(e(r, c), M(r, c), 1e-9);
		}
	}
}

TEST(SiLi, mat_to_quat) {
	SiLi::Matrix<4, 4, double> m ({{1., 0., 0., 0.},
	                              {0., 1./std::sqrt(2.), -1./std::sqrt(2.), 0.},
	                              {0., 1./std::sqrt(2.),  1./std::sqrt(2.), 0.},
	                              {0., 0., 0., 1.}});
	SiLi::Quaternion<double> q = m;
	SiLi::Quaternion<double> e(M_PI * 0.25, 1., 0., 0.);

	EXPECT_NEAR(e(0), q(0), 1e-9);
	EXPECT_NEAR(e(1), q(1), 1e-9);
	EXPECT_NEAR(e(2), q(2), 1e-9);
	EXPECT_NEAR(e(3), q(3), 1e-9);

}









TEST(SiLi, quat_test) {
	SiLi::Quaternion<double> q1, q2;
	q1 += q2;
	auto q3 = q1 + q2;
	auto q4 = q1 - q3;
	q1 -= q4;
	auto q5 = -q1;
	auto q6 = q5 * 5.;
	q6 *= 3.;
	auto r = q2.norm();
	(void)r;
	auto q7 = q6.conjugate();

	auto q8 = q7 * q6;
	q8 *= q7;

	auto q9 = q8.slerp(q7, 0.5);
	auto a = q8.dot(q9);
	(void)a;

	SiLi::Quaternion<double> q10;

	auto mat = q8.mat();
	SiLi::Quaternion<double> q11(mat.view<3, 3>(0, 0));

	std::cout << "q10: " << q10 << "\n";
	std::cout << "q11: " << q11 << "\n";
	std::cout << "mat: " << mat << "\n";


}
