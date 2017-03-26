#include <SiLi/SiLi.h>
#include <gtest/gtest.h>

#if __GNUC__ > 4

using namespace SiLi;
template <int rows, int cols>
using M = Matrix<rows, cols>;


TEST(SiLi, init0) {
	M< 0, 0>  m1(0);
	M<10, 0>  m2(0);
	M< 0, 10> m3(0);

	EXPECT_EQ(0, m1.num_rows());
	EXPECT_EQ(0, m1.num_cols());

	EXPECT_EQ(10, m2.num_rows());
	EXPECT_EQ(0, m2.num_cols());

	EXPECT_EQ(0, m3.num_rows());
	EXPECT_EQ(10, m3.num_cols());
}

TEST(SiLi, init1) {
	M< 1, 1>  m1(0);
	M<10, 1>  m2(0);
	M< 1, 10> m3(0);

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

TEST(SiLi, matrixView_make_vec) {
	auto m = SiLi::make_vec({11, 21, 31});

	EXPECT_EQ(3, m.num_rows());
	EXPECT_EQ(1, m.num_cols());

	EXPECT_EQ(11,  m(0, 0));
	EXPECT_EQ(21,  m(1, 0));
	EXPECT_EQ(31,  m(2, 0));
}


TEST(SiLi, matrixView_det11) {
	auto m = SiLi::Matrix<1, 1, double>({{11.}});
	EXPECT_EQ(m.det(), 11.);
}

TEST(SiLi, matrixView_det22) {
	auto m = SiLi::make_mat<2, 2, double>({{11., 12.},
	                         {21., 22.}});
	EXPECT_NEAR(m.det(), -10, 1e-9);
}

TEST(SiLi, matrixView_det33) {
	auto m = SiLi::make_mat<3, 3, double>({{1, 0, 0},
	                         {0, 1, 0},
	                         {0, 0, 1}});
	EXPECT_EQ(m.det(), 1);
}

TEST(SiLi, matrixView_det33_b) {
	auto m = SiLi::make_mat<3, 3, double>({{1., 2., 3.}, {2., 3., 4.}, {4., 1., 10.}});

	EXPECT_NEAR(m.det(), -12., 1e-9);
}

TEST(SiLi, matrixView_det44) {
	auto m = SiLi::make_mat<4, 4, double>({{1., 2., 3., 10.}, {2., 3., 4., 20.}, {4., -1., 10., 50.}, {10., -11., 12., -13.}});

	EXPECT_NEAR(m.det(), 2248., 1e-9);
}



TEST(SiLi, matrixView_inv33) {
	auto m = SiLi::make_mat<3, 3, double>({{1., 0., 0.},
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
	auto m = SiLi::make_mat<4, 4, double>({{1., 2., 3., 4.},
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

TEST(SiLi, matrixView_svd2x2) {
	auto m = SiLi::make_mat<2, 2, double>({{1., 2.}, {2., 3.}});

	auto svd = m.svd();
	auto S = make_diag(svd.S);

	auto re = svd.U * S * svd.V.t_view();
	EXPECT_NEAR(std::abs(svd.U.det()), 1., 1e-9);
	EXPECT_NEAR(std::abs(svd.V.det()), 1., 1e-9);
	EXPECT_NEAR(m(0, 0), re(0, 0), 1e-9);
	EXPECT_NEAR(m(1, 0), re(1, 0), 1e-9);
	EXPECT_NEAR(m(0, 1), re(0, 1), 1e-9);
	EXPECT_NEAR(m(1, 1), re(1, 1), 1e-9);
}

TEST(SiLi, matrixView_svd3x3) {
	auto m = SiLi::make_mat<3, 3, double>({{1., 2., 3.}, {2., 3., 4.}, {4., 1., 10.}});

	auto svd = m.svd();
	auto S = make_diag(svd.S);

	auto re = svd.U * S * svd.V.t_view();
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

TEST(SiLi, matrixView_svd4x4) {
	auto m = SiLi::make_mat<4, 4, double>({{1., 2., 3., 10.}, {2., 3., 4., 20.}, {4., -1., 10., 50.}, {10., -11., 12., -13.}});

	auto svd = m.svd();
	auto S = make_diag(svd.S);

	auto re = svd.U * S * svd.V.t_view();
	EXPECT_NEAR(std::abs(svd.U.det()), 1., 1e-9);
	EXPECT_NEAR(std::abs(svd.V.det()), 1., 1e-9);
	for (int row(0); row<m.num_rows(); ++row) {
		for (int col(0); col<m.num_cols(); ++col) {
			EXPECT_NEAR(m(row, col), re(row, col), 1e-9);
		}
	}
}

TEST(SiLi, matrixView_svd5x5) {
	auto m = SiLi::make_mat<5, 5, double>({{1., 2., 3., 10., -100.}, {2., 3., 4., 20., 101.}, {4., -1., 10., 50., -102.}, {10., -11., 12., -13., -103.}, {104., -105., 106., 107., -108.}});

	auto svd = m.svd();
	auto S = make_diag(svd.S);

	auto re = svd.U * S * svd.V.t_view();
	EXPECT_NEAR(std::abs(svd.U.det()), 1., 1e-9);
	EXPECT_NEAR(std::abs(svd.V.det()), 1., 1e-9);
	for (int row(0); row<m.num_rows(); ++row) {
		for (int col(0); col<m.num_cols(); ++col) {
			EXPECT_NEAR(m(row, col), re(row, col), 1e-9);
		}
	}
}


TEST(SiLi, svd_random) {
	auto m = SiLi::make_mat<4, 5, double>({
			 {-0.050500,  -0.978799,   0.165052,   1.532855,  -0.848506},
			 {-1.569933,   0.264375,  -0.156239,  -0.590324,   0.318738},
			 {-1.046523,  -0.668761,  -0.304812,  -0.026862,   1.034121},
			 {-0.414202,   0.031721,  -1.523946,  -0.423158,  -0.925514}
	});

	auto U_true = SiLi::make_mat<4, 5, double>({
		{ 0.5820622,  -0.4900696,   0.6054104,   0.2334813, 0.},
		{-0.6432321,  -0.0995143,   0.2632539,   0.7120722, 0.},
		{-0.4258450,   0.0033398,   0.6534692,  -0.6257978, 0.},
		{-0.2571226,  -0.8659778,  -0.3703425,  -0.2163721, 0.}
	});
	auto S_true = SiLi::make_mat<5, 5, double>({
		   {2.41447,         0,         0,         0,         0},
		   {      0,   1.88627,         0,         0,         0},
		   {      0,         0,   1.80314,         0,         0},
		   {      0,         0,         0,   0.81464,         0},
		   {      0,         0,         0,         0,         0}

	});
	auto V_true = SiLi::make_mat<5, 5, double>({
		{ 0.634753,   0.284251,  -0.540357,  -0.472804,  -0.027580},
		{-0.191820,   0.224606,  -0.538916,   0.455870,   0.643824},
		{ 0.297461,   0.664458,   0.235140,   0.549661,  -0.335551},
		{ 0.576595,  -0.172883,   0.505653,   0.056356,   0.615458},
		{-0.373295,   0.630365,   0.326508,  -0.513163,   0.305529}
	});

	auto svd = m.svd();
	auto S = make_diag(svd.S);

	auto re = svd.U * S * svd.V.t_view();
	auto diff = m - re;;
	EXPECT_NEAR(0., diff.norm(), 1e-9);
	auto Udet = svd.U.view<4, 4>(0, 0).det();
	EXPECT_NEAR(1., std::abs(Udet), 1e-9);
	EXPECT_NEAR(1., std::abs(svd.V.det()), 1e-9);


	for (int row(0); row<U_true.num_rows(); ++row) {
		for (int col(0); col<U_true.num_cols(); ++col) {
			EXPECT_NEAR(std::abs(svd.U(row, col)), std::abs(U_true(row, col)), 1e-5);
		}
	}
	for (int row(0); row<S_true.num_rows(); ++row) {
		for (int col(0); col<S_true.num_cols(); ++col) {
			EXPECT_NEAR(std::abs(S(row, col)), std::abs(S_true(row, col)), 1e-5);
		}
	}
	for (int row(0); row<V_true.num_rows(); ++row) {
		for (int col(0); col<V_true.num_cols(); ++col) {
			EXPECT_NEAR(std::abs(svd.V(row, col)), std::abs(V_true(row, col)), 1e-5);
		}
	}
}


TEST(SiLi, matrixView_transposed) {
	auto m = SiLi::make_mat<2, 2, double>({{11., 12.},
	                         {21., 22.}});
	auto view = m.t_view();
	EXPECT_NEAR(view.det(), -10., 1e-9);
	EXPECT_EQ(view(0, 0), 11.);
	EXPECT_EQ(view(1, 0), 12.);
	EXPECT_EQ(view(0, 1), 21.);
	EXPECT_EQ(view(1, 1), 22.);
}

TEST(SiLi, matrixView_transposed2) {
	auto m = SiLi::make_mat<3, 2, double>({{11., 12.},
	                         {21., 22.},
	                         {31., 32.}});
	auto view = m.t_view();
	EXPECT_EQ(view.num_rows(), 2);
	EXPECT_EQ(view.num_cols(), 3);
	EXPECT_EQ(view(0, 0), 11.);
	EXPECT_EQ(view(1, 0), 12.);
	EXPECT_EQ(view(0, 1), 21.);
	EXPECT_EQ(view(1, 1), 22.);
	EXPECT_EQ(view(0, 2), 31.);
	EXPECT_EQ(view(1, 2), 32.);
}

TEST(SiLi, matrixView_transposed_transposed) {
	auto m = SiLi::make_mat<2, 2, double>({{11., 12.},
	                         {21., 22.}});
	auto view = m.t_view().t_view();
	EXPECT_NEAR(view.det(), -10., 1e-9);
	EXPECT_EQ(view(0, 0), 11.);
	EXPECT_EQ(view(1, 0), 21.);
	EXPECT_EQ(view(0, 1), 12.);
	EXPECT_EQ(view(1, 1), 22.);
}

TEST(SiLi, matrixView_implicit_cast_to_T) {
	auto m = SiLi::make_mat<1, 1, double>({{11.}});

	double x = m;
	EXPECT_EQ(x, 11.);
}

TEST(SiLi, matrixView_add) {
	auto m1 = SiLi::make_mat<2, 2, double>({{11., 12.},
	                          {21., 22.}});
	auto m2 = SiLi::make_mat<2, 2, double>({{111., 112.},
	                          {121., 122.}});

	auto m3 = m1 + m2;
	EXPECT_EQ(m3(0, 0), 122.);
	EXPECT_EQ(m3(1, 0), 142.);
	EXPECT_EQ(m3(0, 1), 124.);
	EXPECT_EQ(m3(1, 1), 144.);
}

TEST(SiLi, matrixView_sub) {
	auto m1 = SiLi::make_mat<2, 2, double>({{11., 12.},
	                          {21., 22.}});
	auto m2 = SiLi::make_mat<2, 2, double>({{111., 112.},
	                          {121., 122.}});

	auto m3 = m2 - m1;
	EXPECT_EQ(m3(0, 0), 100.);
	EXPECT_EQ(m3(1, 0), 100.);
	EXPECT_EQ(m3(0, 1), 100.);
	EXPECT_EQ(m3(1, 1), 100.);
}

TEST(SiLi, matrixView_scale) {
	auto m1 = SiLi::make_mat<2, 2, double>({{11., 12.},
	                          {21., 22.}});

	auto m3 = m1 * 2.;
	EXPECT_EQ(m3(0, 0), 22.);
	EXPECT_EQ(m3(1, 0), 42.);
	EXPECT_EQ(m3(0, 1), 24.);
	EXPECT_EQ(m3(1, 1), 44.);
}
TEST(SiLi, matrixView_neg) {
	auto m1 = SiLi::make_mat<2, 2, double>({{11., 12.},
	                          {21., 22.}});

	auto m3 = -m1;
	EXPECT_EQ(m3(0, 0), -11.);
	EXPECT_EQ(m3(1, 0), -21.);
	EXPECT_EQ(m3(0, 1), -12.);
	EXPECT_EQ(m3(1, 1), -22.);
}
TEST(SiLi, matrixView_diag) {
	auto m = SiLi::make_mat<3, 2, double>({{11., 12.},
	                         {21., 22.},
	                         {31., 32.}});

	EXPECT_EQ(11., m.diag()(0, 0));
	EXPECT_EQ(22., m.diag()(1, 0));

	EXPECT_EQ(11., m.diag().t_view()(0, 0));
	EXPECT_EQ(22., m.diag().t_view()(0, 1));

	EXPECT_EQ(11., m.diag().t_view().diag()(0, 0));
	EXPECT_EQ(1, m.diag().t_view().diag().num_rows());
	EXPECT_EQ(1, m.diag().t_view().diag().num_cols());

	EXPECT_EQ(11., m.t_view().diag()(0, 0));
	EXPECT_EQ(22., m.t_view().diag()(1, 0));
}

TEST(SiLi, matrixView_diag2) {
	auto m = SiLi::make_mat<3, 2, double>({{11., 12.},
	                         {21., 22.},
	                         {31., 32.}});

	EXPECT_EQ(11., m.diag()(0));
	EXPECT_EQ(22., m.diag()(1));
	EXPECT_EQ(11., m.diag().t_view()(0));
	EXPECT_EQ(22., m.diag().t_view()(1));
}

TEST(SiLi, matrixView_make_diag) {
	auto v = SiLi::make_mat<3, 1, double>({{11.},
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

TEST(SiLi, matrixView_make_diag2) {
	auto v = SiLi::make_mat<3, 1, double>({{11.},
	                         {21.},
	                         {31.}});
	auto m = make_diag<3,4>(v);

	EXPECT_EQ(11., m(0, 0));
	EXPECT_EQ(0.,  m(0, 1));
	EXPECT_EQ(0.,  m(0, 2));
	EXPECT_EQ(0.,  m(0, 3));
	EXPECT_EQ(0.,  m(1, 0));
	EXPECT_EQ(21., m(1, 1));
	EXPECT_EQ(0.,  m(1, 2));
	EXPECT_EQ(0.,  m(1, 3));
	EXPECT_EQ(0.,  m(2, 0));
	EXPECT_EQ(0.,  m(2, 1));
	EXPECT_EQ(31., m(2, 2));
	EXPECT_EQ(0.,  m(2, 3));
}

TEST(SiLi, matrixView_make_diag3) {
	auto v = SiLi::make_mat<3, 1, double>({{11.},
	                         {21.},
	                         {31.}});
	auto m = make_diag<4,3>(v);

	EXPECT_EQ(11., m(0, 0));
	EXPECT_EQ(0.,  m(0, 1));
	EXPECT_EQ(0.,  m(0, 2));
	EXPECT_EQ(0.,  m(1, 0));
	EXPECT_EQ(21., m(1, 1));
	EXPECT_EQ(0.,  m(1, 2));
	EXPECT_EQ(0.,  m(2, 0));
	EXPECT_EQ(0.,  m(2, 1));
	EXPECT_EQ(31., m(2, 2));
	EXPECT_EQ(0.,  m(3, 0));
	EXPECT_EQ(0.,  m(3, 1));
	EXPECT_EQ(0.,  m(3, 2));
}

TEST(SiLi, matrixView_element_wise_mul) {
	auto v = SiLi::make_mat<3, 1, double>({{11.},
	                         {21.},
	                         {31.}});
	auto m = make_eye<3, 3, double>();
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
	auto v = SiLi::make_mat<3, 1, double>({{11.},
	                         {21.},
	                         {31.}});
	auto m = make_eye<3, 3, double>();
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
	auto m = SiLi::make_mat<3, 2, double>({{11., 12.},
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
	auto m = SiLi::make_mat<3, 2, double>({{11., 12.},
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
	auto m = SiLi::make_mat<3, 2, double>({{11., 12.},
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
	auto m = SiLi::make_mat<3, 1, double>({{1.},
	                         {2.},
	                         {3.}});

	double s = m.t_view() * m;
	EXPECT_EQ(14., s);
}

TEST(SiLi, matrixView_abs) {
	auto m = SiLi::make_mat<3, 4, double>({{-1., 2., -3., 4.},
	                         {6., -6., 7., -8.},
	                         {9., 10., -10., 12.}});
	auto a = abs(m);
	EXPECT_EQ(a(0, 0), 1.);
	EXPECT_EQ(a(0, 1), 2.);
	EXPECT_EQ(a(0, 2), 3.);
	EXPECT_EQ(a(0, 3), 4.);
	EXPECT_EQ(a(1, 0), 6.);
	EXPECT_EQ(a(1, 1), 6.);
	EXPECT_EQ(a(1, 2), 7.);
	EXPECT_EQ(a(1, 3), 8.);
	EXPECT_EQ(a(2, 0), 9.);
	EXPECT_EQ(a(2, 1), 10.);
	EXPECT_EQ(a(2, 2), 10.);
	EXPECT_EQ(a(2, 3), 12.);
}

TEST(SiLi, matrixView_isfinite_true) {
	auto m = SiLi::make_mat<3, 4, double>({{-1., 2., -3., 4.},
	                         {6., -6., 7., -8.},
	                         {9., 10., -10., 12.}});
	EXPECT_EQ(isfinite(m), true);
}
TEST(SiLi, matrixView_isfinite_false) {
	auto m = SiLi::make_mat<3, 4, double>({{-1., 2., -3., 4.},
	                         {6., -6./0., 7., -8.},
	                         {9., 10., -10., 12.}});
	EXPECT_EQ(isfinite(m), false);
}

TEST(SiLi, matrixView_join_rows) {
	M<2, 1>  m1({{11.}, {21.}});
	M<2, 2>  m2({{111., 112.}, {121., 122.}});

	auto m3 = m1.join_rows(m2);
	EXPECT_EQ(11., m3(0, 0));
	EXPECT_EQ(21., m3(1, 0));
	EXPECT_EQ(111., m3(0, 1));
	EXPECT_EQ(112., m3(0, 2));
	EXPECT_EQ(121., m3(1, 1));
	EXPECT_EQ(122., m3(1, 2));
}

TEST(SiLi, matrixView_join_cols) {
	M<1, 2>  m1({{11., 21.}});
	M<2, 2>  m2({{111., 112.}, {121., 122.}});

	auto m3 = m1.join_cols(m2);
	EXPECT_EQ(11., m3(0, 0));
	EXPECT_EQ(21., m3(0, 1));
	EXPECT_EQ(111., m3(1, 0));
	EXPECT_EQ(112., m3(1, 1));
	EXPECT_EQ(121., m3(2, 0));
	EXPECT_EQ(122., m3(2, 1));
}


TEST(SiLi, matrixView_rowRange1) {
	M<4, 2>  m({
			{00., 01.},
			{10., 11.},
			{20., 21.},
			{30., 31.}
	});

	auto rows = m.rows<0, 3, 2>();
	EXPECT_EQ(rows.size(), 2);

	EXPECT_EQ(rows[0](0, 0), 00.);
	EXPECT_EQ(rows[0](0, 1), 01.);

	EXPECT_EQ(rows[1](0, 0), 20.);
	EXPECT_EQ(rows[1](0, 1), 21.);
}

TEST(SiLi, matrixView_rowRange2) {
	M<4, 2>  m({
			{00., 01.},
			{10., 11.},
			{20., 21.},
			{30., 31.}
	});

	auto rows = m.rows<3, 0, -2>();
	EXPECT_EQ(rows.size(), 2);

	EXPECT_EQ(rows[0](0, 0), 30.);
	EXPECT_EQ(rows[0](0, 1), 31.);

	EXPECT_EQ(rows[1](0, 0), 10.);
	EXPECT_EQ(rows[1](0, 1), 11.);

}


TEST(SiLi, matrixView_colRange1) {
	M<2, 4>  m({
			{00., 01. , 02., 03.},
			{10., 11. , 12., 13.}
	});

	auto cols = m.cols<0, 3, 2>();
	EXPECT_EQ(cols.size(), 2);

	EXPECT_EQ(cols[0](0, 0), 00.);
	EXPECT_EQ(cols[0](1, 0), 10.);

	EXPECT_EQ(cols[1](0, 0), 02.);
	EXPECT_EQ(cols[1](1, 0), 12.);
}

TEST(SiLi, matrixView_colRange2) {
	M<2, 4>  m({
			{00., 01. , 02., 03.},
			{10., 11. , 12., 13.}
	});

	auto cols = m.cols<3, 0, -2>();
	EXPECT_EQ(cols.size(), 2);

	EXPECT_EQ(cols[0](0, 0), 03.);
	EXPECT_EQ(cols[0](1, 0), 13.);

	EXPECT_EQ(cols[1](0, 0), 01.);
	EXPECT_EQ(cols[1](1, 0), 11.);
}
#else
	TEST(SiLi, gcc4) {
		EXPECT_TRUE(false); // gcc4 is not supported
	}

#endif
