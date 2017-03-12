#include <SiLi/SiLi.h>
#include <gtest/gtest.h>

using namespace SiLi;
template <int rows, int cols>
using M = Matrix<rows, cols>;


TEST(SiLi, init0) {
	M< 0, 0>  m1({});
	M<10, 0>  m2({});
	M< 0, 10> m3({});

	EXPECT_EQ(0, m1.num_rows());
	EXPECT_EQ(0, m1.num_cols());

	EXPECT_EQ(10, m2.num_rows());
	EXPECT_EQ(0, m2.num_cols());

	EXPECT_EQ(0, m3.num_rows());
	EXPECT_EQ(10, m3.num_cols());
}

TEST(SiLi, init1) {
	M< 1, 1>  m1({});
	M<10, 1>  m2({});
	M< 1, 10> m3({});

	EXPECT_EQ(1, m1.num_rows());
	EXPECT_EQ(1, m1.num_cols());

	EXPECT_EQ(10, m2.num_rows());
	EXPECT_EQ(1, m2.num_cols());

	EXPECT_EQ(1, m3.num_rows());
	EXPECT_EQ(10, m3.num_cols());
}

TEST(SiLi, matrixView_access) {
	M<2, 2>  m({{11, 12}, {21, 22}});
	M<2, 2>::View view = m;

	EXPECT_EQ(2, m.num_rows());
	EXPECT_EQ(2, m.num_cols());

	EXPECT_EQ(2, view.num_rows());
	EXPECT_EQ(2, view.num_cols());

	EXPECT_EQ(11,  m(0, 0));
	EXPECT_EQ(21,  m(1, 0));
	EXPECT_EQ(12,  m(0, 1));
	EXPECT_EQ(22,  m(1, 1));


	EXPECT_EQ(view(0, 0),  m(0, 0));
	EXPECT_EQ(view(1, 0),  m(1, 0));
	EXPECT_EQ(view(0, 1),  m(0, 1));
	EXPECT_EQ(view(1, 1),  m(1, 1));

	EXPECT_EQ(&view(0, 0),  &m(0, 0));
	EXPECT_EQ(&view(1, 0),  &m(1, 0));
	EXPECT_EQ(&view(0, 1),  &m(0, 1));
	EXPECT_EQ(&view(1, 1),  &m(1, 1));
}

TEST(SiLi, matrixView_access_write) {
	M<2, 2>  m({{11, 12}, {21, 22}});
	M<2, 2>::View view = m;

	EXPECT_EQ(11,  m(0, 0));
	EXPECT_EQ(21,  m(1, 0));
	EXPECT_EQ(12,  m(0, 1));
	EXPECT_EQ(22,  m(1, 1));

	m(0, 0) *= 2;
	m(1, 0) *= 2;
	m(0, 1) *= 2;
	m(1, 1) *= 2;

	EXPECT_EQ(22,  view(0, 0));
	EXPECT_EQ(42,  view(1, 0));
	EXPECT_EQ(24,  view(0, 1));
	EXPECT_EQ(44,  view(1, 1));

	view(0, 0) /= 2;
	view(1, 0) /= 2;
	view(0, 1) /= 2;
	view(1, 1) /= 2;

	EXPECT_EQ(11,  m(0, 0));
	EXPECT_EQ(21,  m(1, 0));
	EXPECT_EQ(12,  m(0, 1));
	EXPECT_EQ(22,  m(1, 1));
}

TEST(SiLi, matrixView_const_view) {
	M<2, 2>  m({{11, 12}, {21, 22}});
	M<2, 2>::View const view = m;
	M<2, 2>::CView view2 = m;

	EXPECT_EQ(11,  view(0, 0));
	EXPECT_EQ(21,  view(1, 0));
	EXPECT_EQ(12,  view(0, 1));
	EXPECT_EQ(22,  view(1, 1));

	EXPECT_EQ(11,  view2(0, 0));
	EXPECT_EQ(21,  view2(1, 0));
	EXPECT_EQ(12,  view2(0, 1));
	EXPECT_EQ(22,  view2(1, 1));


	//view(0, 0) = 0; // not possible
	//view2(0, 0) = 0; // not possible

	EXPECT_EQ(11,  view(0, 0));
	EXPECT_EQ(21,  view(1, 0));
	EXPECT_EQ(12,  view(0, 1));
	EXPECT_EQ(22,  view(1, 1));
}

TEST(SiLi, matrixView_view_to_mat) {
	M<2, 2>  m({{11, 12}, {21, 22}});
	M<2, 2>::View const view = m;

	EXPECT_EQ(11,  view(0, 0));
	EXPECT_EQ(21,  view(1, 0));
	EXPECT_EQ(12,  view(0, 1));
	EXPECT_EQ(22,  view(1, 1));

	M<1, 1>  m2({{0}});
	m2 = view.mat<1, 1> (1, 1);

	auto m3 = view.mat<1, 1> (1, 1);

	EXPECT_EQ(1, m2.num_rows());
	EXPECT_EQ(1, m2.num_cols());

	EXPECT_EQ(22, m2(0, 0));

	EXPECT_EQ(1, m3.num_rows());
	EXPECT_EQ(1, m3.num_cols());

	EXPECT_EQ(22, m3(0, 0));
}

TEST(SiLi, matrixView_view) {
	M<3, 3>  m({{11, 12, 13},
	            {21, 22, 23},
	            {31, 32, 33}});
	auto const view1 = m.view<2, 1>(1, 2);

	EXPECT_EQ(23,  view1(0, 0));
	EXPECT_EQ(33,  view1(1, 0));

	auto view2 = m.view<2, 1>(1, 2);
	EXPECT_EQ(23,  view2(0, 0));
	EXPECT_EQ(33,  view2(1, 0));

	auto view3 = m.view<2, 1>(1, 2).view<2, 1>(0, 0);
	EXPECT_EQ(23,  view3(0, 0));
	EXPECT_EQ(33,  view3(1, 0));
}

TEST(SiLi, matrixView_iterator0) {
	std::vector<double> v;
	M<3, 3>  m({{11, 12, 13},
	            {21, 22, 23},
	            {31, 32, 33}});
	auto iter = m.begin();
	(void)iter;
}
TEST(SiLi, matrixView_iterator) {
	std::vector<double> v;
	M<3, 3>  m({{11, 12, 13},
	            {21, 22, 23},
	            {31, 32, 33}});

	for (auto x : m) {
		v.push_back(x);
	}
	ASSERT_EQ(9, v.size());
	EXPECT_EQ(11, v[0]);
	EXPECT_EQ(12, v[1]);
	EXPECT_EQ(13, v[2]);
	EXPECT_EQ(21, v[3]);
	EXPECT_EQ(22, v[4]);
	EXPECT_EQ(23, v[5]);
	EXPECT_EQ(31, v[6]);
	EXPECT_EQ(32, v[7]);
	EXPECT_EQ(33, v[8]);
}
TEST(SiLi, matrixView_iterator2) {
	std::vector<double> v;
	M<3, 3>  m({{11, 12, 13},
	            {21, 22, 23},
	            {31, 32, 33}});

	for (auto& x : m) {
		x += 1.;
	}

	for (auto x : m) {
		v.push_back(x);
	}
	ASSERT_EQ(9, v.size());
	EXPECT_EQ(12, v[0]);
	EXPECT_EQ(13, v[1]);
	EXPECT_EQ(14, v[2]);
	EXPECT_EQ(22, v[3]);
	EXPECT_EQ(23, v[4]);
	EXPECT_EQ(24, v[5]);
	EXPECT_EQ(32, v[6]);
	EXPECT_EQ(33, v[7]);
	EXPECT_EQ(34, v[8]);
}

#if __GNUC__ > 4
TEST(SiLi, matrixView_make_matrix) {
	auto m = SiLi::make_mat({{11, 12, 13},
	                         {21, 22, 23}});

	EXPECT_EQ(2, m.num_rows());
	EXPECT_EQ(3, m.num_cols());

	EXPECT_EQ(11,  m(0, 0));
	EXPECT_EQ(21,  m(1, 0));
	EXPECT_EQ(12,  m(0, 1));
	EXPECT_EQ(22,  m(1, 1));
	EXPECT_EQ(13,  m(0, 2));
	EXPECT_EQ(23,  m(1, 2));
}
#endif

TEST(SiLi, matrixView_det11) {
	auto m = SiLi::Matrix<1, 1, double>({{11.}});
	EXPECT_EQ(m.det(), 11.);
}

TEST(SiLi, matrixView_det22) {
	auto m = SiLi::make_mat<double, 2, 2>({{11., 12.},
	                         {21., 22.}});
	EXPECT_EQ(m.det(), -10);
}

TEST(SiLi, matrixView_det33) {
	auto m = SiLi::make_mat<double, 3, 3>({{1, 0, 0},
	                         {0, 1, 0},
	                         {0, 0, 1}});
	EXPECT_EQ(m.det(), 1);
}

TEST(SiLi, matrixView_inv33) {
	auto m = SiLi::make_mat<double, 3, 3>({{1., 0., 0.},
	                         {0., 1., 0.},
	                         {0., 0., 1.}});

	std::vector<double> v {1, 0, 0, 0, 1, 0, 0, 0, 1};
	std::vector<double> list;
	for (auto x : inv(m)) {
		list.push_back(x);
	}
	ASSERT_EQ(v.size(), list.size());
	for (size_t i(0); i < v.size(); ++i) {
		EXPECT_EQ(v[i], list[i]);
	}
}

TEST(SiLi, matrixView_inv44) {
	auto m = SiLi::make_mat<double, 4, 4>({{1., 2., 3., 4.},
	                         {6., 6., 7., 8.},
	                         {9., 10., 10., 12.},
	                         {13., 14., 15., 15.}});


	std::vector<double> v {-13./19, 17./19, -2./19, -4./19, 7./38, -53./38, 23./38, 4/19.,
	                       7./19, 4./19, -15./19, 8./19, 1./19, 6./19, 6./19, -7./19};
	std::vector<double> list;
	for (auto x : inv(m)) {
		list.push_back(x);
	}
	ASSERT_EQ(v.size(), list.size());
	for (size_t i(0); i < v.size(); ++i) {
		EXPECT_NEAR(v[i], list[i], 1e-9);
	}
}

TEST(SiLi, matrixView_svd) {
	auto m = SiLi::make_mat<double, 2, 2>({{1., 2.}, {2., 3.}});

	auto svd = m.svd();
	auto S = make_diag(svd.S);

	auto re = svd.U * S * svd.V.t();
	EXPECT_NEAR(std::abs(svd.U.det()), 1., 1e-9);
	EXPECT_NEAR(std::abs(svd.V.det()), 1., 1e-9);
	EXPECT_NEAR(m(0, 0), re(0, 0), 1e-9);
	EXPECT_NEAR(m(1, 0), re(1, 0), 1e-9);
	EXPECT_NEAR(m(0, 1), re(0, 1), 1e-9);
	EXPECT_NEAR(m(1, 1), re(1, 1), 1e-9);
}

TEST(SiLi, matrixView_svd2) {
	auto m = SiLi::make_mat<double, 3, 3>({{1., 2., 3.}, {2., 3., 4.}, {4., 1., 10.}});

	auto svd = m.svd();
	auto S = make_diag(svd.S);

	auto re = svd.U * S * svd.V.t();
	EXPECT_NEAR(std::abs(svd.U.det()), 1., 1e-9);
	EXPECT_NEAR(std::abs(svd.V.det()), 1., 1e-9);
	EXPECT_NEAR(m(0, 0), re(0, 0), 1e-9);
	EXPECT_NEAR(m(1, 0), re(1, 0), 1e-9);
	EXPECT_NEAR(m(2, 0), re(2, 0), 1e-9);
	EXPECT_NEAR(m(0, 1), re(0, 1), 1e-9);
	EXPECT_NEAR(m(1, 1), re(1, 1), 1e-9);
	EXPECT_NEAR(m(2, 1), re(2, 1), 1e-9);
	EXPECT_NEAR(m(0, 2), re(0, 2), 1e-9);
	EXPECT_NEAR(m(1, 2), re(1, 2), 1e-9);
	EXPECT_NEAR(m(2, 2), re(2, 2), 1e-9);
}

TEST(SiLi, matrixView_svd3) {
	auto m = SiLi::make_mat<double, 4, 4>({{1., 2., 3., 10.}, {2., 3., 4., 20.}, {4., -1., 10., 50.}, {10., -11., 12., -13.}});

	auto svd = m.svd();
	auto S = make_diag(svd.S);

	auto re = svd.U * S * svd.V.t();
	EXPECT_NEAR(std::abs(svd.U.det()), 1., 1e-9);
	EXPECT_NEAR(std::abs(svd.V.det()), 1., 1e-9);
	EXPECT_NEAR(m(0, 0), re(0, 0), 1e-9);
	EXPECT_NEAR(m(1, 0), re(1, 0), 1e-9);
	EXPECT_NEAR(m(2, 0), re(2, 0), 1e-9);
	EXPECT_NEAR(m(3, 0), re(3, 0), 1e-9);
	EXPECT_NEAR(m(0, 1), re(0, 1), 1e-9);
	EXPECT_NEAR(m(1, 1), re(1, 1), 1e-9);
	EXPECT_NEAR(m(2, 1), re(2, 1), 1e-9);
	EXPECT_NEAR(m(3, 1), re(3, 1), 1e-9);
	EXPECT_NEAR(m(0, 2), re(0, 2), 1e-9);
	EXPECT_NEAR(m(1, 2), re(1, 2), 1e-9);
	EXPECT_NEAR(m(2, 2), re(2, 2), 1e-9);
	EXPECT_NEAR(m(3, 2), re(3, 2), 1e-9);
	EXPECT_NEAR(m(0, 3), re(0, 3), 1e-9);
	EXPECT_NEAR(m(1, 3), re(1, 3), 1e-9);
	EXPECT_NEAR(m(2, 3), re(2, 3), 1e-9);
	EXPECT_NEAR(m(3, 3), re(3, 3), 1e-9);
}



TEST(SiLi, matrixView_transposed) {
	auto m = SiLi::make_mat<double, 2, 2>({{11., 12.},
	                         {21., 22.}});
	auto view = m.t();
	EXPECT_EQ(view.det(), -10);
	EXPECT_EQ(view(0, 0), 11.);
	EXPECT_EQ(view(1, 0), 12.);
	EXPECT_EQ(view(0, 1), 21.);
	EXPECT_EQ(view(1, 1), 22.);
}

TEST(SiLi, matrixView_transposed2) {
	auto m = SiLi::make_mat<double, 3, 2>({{11., 12.},
	                         {21., 22.},
	                         {31., 32.}});
	auto view = m.t();
	EXPECT_EQ(view.num_rows(), 2);
	EXPECT_EQ(view.num_cols(), 3);
	EXPECT_EQ(view(0, 0), 11.);
	EXPECT_EQ(view(1, 0), 12.);
	EXPECT_EQ(view(0, 1), 21.);
	EXPECT_EQ(view(1, 1), 22.);
	EXPECT_EQ(view(0, 2), 31.);
	EXPECT_EQ(view(1, 2), 32.);
}

TEST(SiLi, matrixView_add) {
	auto m1 = SiLi::make_mat<double, 2, 2>({{11., 12.},
	                          {21., 22.}});
	auto m2 = SiLi::make_mat<double, 2, 2>({{111., 112.},
	                          {121., 122.}});

	auto m3 = m1 + m2;
	EXPECT_EQ(m3(0, 0), 122.);
	EXPECT_EQ(m3(1, 0), 142.);
	EXPECT_EQ(m3(0, 1), 124.);
	EXPECT_EQ(m3(1, 1), 144.);
}

TEST(SiLi, matrixView_sub) {
	auto m1 = SiLi::make_mat<double, 2, 2>({{11., 12.},
	                          {21., 22.}});
	auto m2 = SiLi::make_mat<double, 2, 2>({{111., 112.},
	                          {121., 122.}});

	auto m3 = m2 - m1;
	EXPECT_EQ(m3(0, 0), 100.);
	EXPECT_EQ(m3(1, 0), 100.);
	EXPECT_EQ(m3(0, 1), 100.);
	EXPECT_EQ(m3(1, 1), 100.);
}

TEST(SiLi, matrixView_scale) {
	auto m1 = SiLi::make_mat<double, 2, 2>({{11., 12.},
	                          {21., 22.}});

	auto m3 = m1 * 2.;
	EXPECT_EQ(m3(0, 0), 22.);
	EXPECT_EQ(m3(1, 0), 42.);
	EXPECT_EQ(m3(0, 1), 24.);
	EXPECT_EQ(m3(1, 1), 44.);
}
TEST(SiLi, matrixView_neg) {
	auto m1 = SiLi::make_mat<double, 2, 2>({{11., 12.},
	                          {21., 22.}});

	auto m3 = -m1;
	EXPECT_EQ(m3(0, 0), -11.);
	EXPECT_EQ(m3(1, 0), -21.);
	EXPECT_EQ(m3(0, 1), -12.);
	EXPECT_EQ(m3(1, 1), -22.);
}
TEST(SiLi, matrixView_diag) {
	auto m = SiLi::make_mat<double, 3, 2>({{11., 12.},
	                         {21., 22.},
	                         {31., 32.}});

	EXPECT_EQ(11., m.diag()(0, 0));
	EXPECT_EQ(22., m.diag()(1, 0));

	EXPECT_EQ(11., m.diag().t()(0, 0));
	EXPECT_EQ(22., m.diag().t()(0, 1));

	EXPECT_EQ(11., m.diag().t().diag()(0, 0));
	EXPECT_EQ(1, m.diag().t().diag().num_rows());
	EXPECT_EQ(1, m.diag().t().diag().num_cols());

	EXPECT_EQ(11., m.t().diag()(0, 0));
	EXPECT_EQ(22., m.t().diag()(1, 0));
}

TEST(SiLi, matrixView_diag2) {
	auto m = SiLi::make_mat<double, 3, 2>({{11., 12.},
	                         {21., 22.},
	                         {31., 32.}});

	EXPECT_EQ(11., m.diag()(0));
	EXPECT_EQ(22., m.diag()(1));
	EXPECT_EQ(11., m.diag().t()(0));
	EXPECT_EQ(22., m.diag().t()(1));
}

TEST(SiLi, matrixView_make_diag) {
	auto v = SiLi::make_mat<double, 3, 1>({{11.},
	                         {21.},
	                         {31.}});
	auto m = make_diag(v);

	EXPECT_EQ(11., m(0, 0));
	EXPECT_EQ(0.,  m(0, 1));
	EXPECT_EQ(0.,  m(0, 2));
	EXPECT_EQ(0.,  m(1, 0));
	EXPECT_EQ(21., m(1, 1));
	EXPECT_EQ(0.,  m(1, 2));
	EXPECT_EQ(0.,  m(2, 0));
	EXPECT_EQ(0.,  m(2, 1));
	EXPECT_EQ(31., m(2, 2));
}

TEST(SiLi, matrixView_element_wise_mul) {
	auto v = SiLi::make_mat<double, 3, 1>({{11.},
	                         {21.},
	                         {31.}});
	auto m = make_eye<double, 3, 3>();
	m.diag() &= v;

	EXPECT_EQ(11., m(0, 0));
	EXPECT_EQ(0.,  m(0, 1));
	EXPECT_EQ(0.,  m(0, 2));
	EXPECT_EQ(0.,  m(1, 0));
	EXPECT_EQ(21., m(1, 1));
	EXPECT_EQ(0.,  m(1, 2));
	EXPECT_EQ(0.,  m(2, 0));
	EXPECT_EQ(0.,  m(2, 1));
	EXPECT_EQ(31., m(2, 2));
}
TEST(SiLi, matrixView_element_wise_mul2) {
	auto v = SiLi::make_mat<double, 3, 1>({{11.},
	                         {21.},
	                         {31.}});
	auto m = make_eye<double, 3, 3>();
	m.diag() = m.diag() & v;

	EXPECT_EQ(11., m(0, 0));
	EXPECT_EQ(0.,  m(0, 1));
	EXPECT_EQ(0.,  m(0, 2));
	EXPECT_EQ(0.,  m(1, 0));
	EXPECT_EQ(21., m(1, 1));
	EXPECT_EQ(0.,  m(1, 2));
	EXPECT_EQ(0.,  m(2, 0));
	EXPECT_EQ(0.,  m(2, 1));
	EXPECT_EQ(31., m(2, 2));
}

TEST(SiLi, matrixView_swap) {
	auto m = SiLi::make_mat<double, 3, 2>({{11., 12.},
	                         {21., 22.},
	                         {31., 32.}});

	auto view1 = m.view<1, 2>(0, 0);
	auto view2 = m.view<1, 2>(2, 0);

	swap(view1, view2);
	EXPECT_EQ(31., m(0, 0));
	EXPECT_EQ(32., m(0, 1));
	EXPECT_EQ(21., m(1, 0));
	EXPECT_EQ(22., m(1, 1));
	EXPECT_EQ(11., m(2, 0));
	EXPECT_EQ(12., m(2, 1));
}

TEST(SiLi, matrixView_rows) {
	auto m = SiLi::make_mat<double, 3, 2>({{11., 12.},
	                         {21., 22.},
	                         {31., 32.}});

	auto rows = m.rows();

	EXPECT_EQ(11., rows[0](0, 0));
	EXPECT_EQ(12., rows[0](0, 1));
	EXPECT_EQ(21., rows[1](0, 0));
	EXPECT_EQ(22., rows[1](0, 1));
	EXPECT_EQ(31., rows[2](0, 0));
	EXPECT_EQ(32,  rows[2](0, 1));
}

TEST(SiLi, matrixView_cols) {
	auto m = SiLi::make_mat<double, 3, 2>({{11., 12.},
	                         {21., 22.},
	                         {31., 32.}});

	auto cols = m.cols();

	EXPECT_EQ(11., cols[0](0, 0));
	EXPECT_EQ(12., cols[1](0, 0));
	EXPECT_EQ(21., cols[0](1, 0));
	EXPECT_EQ(22., cols[1](1, 0));
	EXPECT_EQ(31., cols[0](2, 0));
	EXPECT_EQ(32,  cols[1](2, 0));
}
TEST(SiLi, matrixView_single_element_to_double) {
	auto m = SiLi::make_mat<double, 3, 1>({{1.},
	                         {2.},
	                         {3.}});

	double s = m.t() * m;
	EXPECT_EQ(14., s);
}
