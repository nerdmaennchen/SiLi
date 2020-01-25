#include <SiLi/SiLi.h>

#include <catch2/catch.hpp>

using namespace SiLi;
template <int rows, int cols>
using M = Matrix<rows, cols>;

TEST_CASE("testing init0","") {
	M< 0, 0>  m1(0);
	M<10, 0>  m2(0);
	M< 0, 10> m3(0);

	CHECK(0 == m1.num_rows());
	CHECK(0 == m1.num_cols());

	CHECK(10 == m2.num_rows());
	CHECK(0 == m2.num_cols());

	CHECK(0 == m3.num_rows());
	CHECK(10 == m3.num_cols());
}

TEST_CASE("testing init1","") {
	M< 1, 1>  m1(0);
	M<10, 1>  m2(0);
	M< 1, 10> m3(0);

	CHECK(1 == m1.num_rows());
	CHECK(1 == m1.num_cols());

	CHECK(10 == m2.num_rows());
	CHECK(1 == m2.num_cols());

	CHECK(1 == m3.num_rows());
	CHECK(10 == m3.num_cols());
}

TEST_CASE("dynamic initialization","") {
	SiLi::Matrix<-1, -1, double> m;
	m = SiLi::make_mat<3, 2, double>({{11., 12.},
	                         {21., 22.},
	                         {31., 32.}});
	CHECK(m(0, 0) == 11.);
	CHECK(m(0, 1) == 12.);

	CHECK(m(1, 0) == 21.);
	CHECK(m(1, 1) == 22.);

	CHECK(m(2, 0) == 31.);
	CHECK(m(2, 1) == 32.);


}

TEST_CASE("matrixView_access","") {
	M<2, 2>  m({{11, 12}, {21, 22}});
	M<2, 2>::View view { m.view<2, 2>(0, 0) };

	CHECK(2 == m.num_rows());
	CHECK(2 == m.num_cols());

	CHECK(2 == view.num_rows());
	CHECK(2 == view.num_cols());

	CHECK(11 ==  m(0, 0));
	CHECK(21 ==  m(1, 0));
	CHECK(12 ==  m(0, 1));
	CHECK(22 ==  m(1, 1));


	CHECK(view(0, 0) ==  m(0, 0));
	CHECK(view(1, 0) ==  m(1, 0));
	CHECK(view(0, 1) ==  m(0, 1));
	CHECK(view(1, 1) ==  m(1, 1));

	CHECK(&view(0, 0) ==  &m(0, 0));
	CHECK(&view(1, 0) ==  &m(1, 0));
	CHECK(&view(0, 1) ==  &m(0, 1));
	CHECK(&view(1, 1) ==  &m(1, 1));
}

TEST_CASE("matrixView_access_write","") {
	M<2, 2>  m({{11, 12}, {21, 22}});
	M<2, 2>::View view { m.view<2, 2>(0, 0) };

	CHECK(11 == m(0, 0));
	CHECK(21 == m(1, 0));
	CHECK(12 == m(0, 1));
	CHECK(22 == m(1, 1));

	m(0, 0) *= 2;
	m(1, 0) *= 2;
	m(0, 1) *= 2;
	m(1, 1) *= 2;

	CHECK(22 == view(0,  0));
	CHECK(42 == view(1,  0));
	CHECK(24 == view(0,  1));
	CHECK(44 == view(1,  1));

	view(0, 0) /= 2;
	view(1, 0) /= 2;
	view(0, 1) /= 2;
	view(1, 1) /= 2;

	CHECK(11 == m(0, 0));
	CHECK(21 == m(1, 0));
	CHECK(12 == m(0, 1));
	CHECK(22 == m(1, 1));
}

TEST_CASE("matrixView_const_view","") {
	M<2, 2>  m({{11, 12}, {21, 22}});
	M<2, 2>::View const view { m.view<2, 2>(0, 0) };
	M<2, 2>::CView view2 { m.view<2, 2>(0, 0) };

	CHECK(11 ==  view(0, 0));
	CHECK(21 ==  view(1, 0));
	CHECK(12 ==  view(0, 1));
	CHECK(22 ==  view(1, 1));

	CHECK(11 ==  view2(0, 0));
	CHECK(21 ==  view2(1, 0));
	CHECK(12 ==  view2(0, 1));
	CHECK(22 ==  view2(1, 1));


	//view(0, 0) = 0; // not possible
	//view2(0, 0) = 0; // not possible

	CHECK(11 ==  view(0, 0));
	CHECK(21 ==  view(1, 0));
	CHECK(12 ==  view(0, 1));
	CHECK(22 ==  view(1, 1));
}

TEST_CASE("matrixView_view_to_mat","") {
	M<2, 2>  m({{11, 12}, {21, 22}});
	M<2, 2>::View const view { m.view<2, 2>(0, 0) };

	CHECK(11 ==  view(0, 0));
	CHECK(21 ==  view(1, 0));
	CHECK(12 ==  view(0, 1));
	CHECK(22 ==  view(1, 1));

	M<1, 1>  m2({{0}});
	m2 = view.mat<1, 1> (1, 1);

	auto m3 = view.mat<1, 1> (1, 1);

	CHECK(1 == m2.num_rows());
	CHECK(1 == m2.num_cols());

	CHECK(22 == m2(0, 0));

	CHECK(1 == m3.num_rows());
	CHECK(1 == m3.num_cols());

	CHECK(22 == m3(0, 0));
}

TEST_CASE("matrixView_view","") {
	M<3, 3>  m({{11, 12, 13},
	            {21, 22, 23},
	            {31, 32, 33}});
	auto const view1 { m.view<2, 1>(1, 2) };

	CHECK(23 ==  view1(0, 0));
	CHECK(33 ==  view1(1, 0));

	auto view2 { m.view<2, 1>(1, 2) };
	CHECK(23 ==  view2(0, 0));
	CHECK(33 ==  view2(1, 0));

	auto view3 { m.view<2, 1>(1, 2).view<2, 1>(0, 0) };
	CHECK(23 ==  view3(0, 0));
	CHECK(33 ==  view3(1, 0));
}


TEST_CASE("matrixView_iterator0","") {
	std::vector<double> v;
	M<3, 3>  m({{11, 12, 13},
	            {21, 22, 23},
	            {31, 32, 33}});
	auto iter = m.begin();
	(void)iter;
}

TEST_CASE("matrixView_iterator1","") {
	std::vector<double> v;
	M<3, 3>  m({{11, 12, 13},
	            {21, 22, 23},
	            {31, 32, 33}});

	for (auto x : m) {
		v.push_back(x);
	}
	REQUIRE(9 == v.size());
	CHECK(11 == v[0]);
	CHECK(12 == v[1]);
	CHECK(13 == v[2]);
	CHECK(21 == v[3]);
	CHECK(22 == v[4]);
	CHECK(23 == v[5]);
	CHECK(31 == v[6]);
	CHECK(32 == v[7]);
	CHECK(33 == v[8]);
}

TEST_CASE("matrixView_iterator2","") {
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
	REQUIRE(9 == v.size());
	CHECK(12 == v[0]);
	CHECK(13 == v[1]);
	CHECK(14 == v[2]);
	CHECK(22 == v[3]);
	CHECK(23 == v[4]);
	CHECK(24 == v[5]);
	CHECK(32 == v[6]);
	CHECK(33 == v[7]);
	CHECK(34 == v[8]);
}

TEST_CASE("matrixView_make_matrix","") {
	auto m = SiLi::make_mat({{11, 12, 13},
	                         {21, 22, 23}});

	CHECK(2 == m.num_rows());
	CHECK(3 == m.num_cols());

	CHECK(11 ==  m(0, 0));
	CHECK(21 ==  m(1, 0));
	CHECK(12 ==  m(0, 1));
	CHECK(22 ==  m(1, 1));
	CHECK(13 ==  m(0, 2));
	CHECK(23 ==  m(1, 2));
}

TEST_CASE("matrixView_make_vec","") {
	auto m = SiLi::make_vec({11, 21, 31});

	CHECK(3 == m.num_rows());
	CHECK(1 == m.num_cols());

	CHECK(11 ==  m(0, 0));
	CHECK(21 ==  m(1, 0));
	CHECK(31 ==  m(2, 0));
}

TEST_CASE("matrixView_det11","") {
	auto m = SiLi::Matrix<1, 1, double>({{11.}});
	CHECK(m.det() == 11.);
}

TEST_CASE("matrixView_det22","") {
	auto m = SiLi::make_mat<2, 2, double>({{11., 12.},
	                         {21., 22.}});
	REQUIRE(m.det() == Approx(-10).margin(1e-9));
}

TEST_CASE("matrixView_det33","") {
	auto m = SiLi::make_mat<3, 3, double>({{1, 0, 0},
	                         {0, 1, 0},
	                         {0, 0, 1}});
	CHECK(m.det() == 1);
}

TEST_CASE("matrixView_det33_b","") {
	auto m = SiLi::make_mat<3, 3, double>({{1., 2., 3.}, {2., 3., 4.}, {4., 1., 10.}});

	CHECK(m.det() == Approx(-12.).margin(1e-9));
}

TEST_CASE("matrixView_det44","") {
	auto m = SiLi::make_mat<4, 4, double>({{1., 2., 3., 10.}, {2., 3., 4., 20.}, {4., -1., 10., 50.}, {10., -11., 12., -13.}});

	CHECK(m.det() == Approx(2248.).margin(1e-9));
}

TEST_CASE("matrixView_det44_dyn","") {
	SiLi::Matrix<-1, -1, double> m;
	m = SiLi::make_mat<4, 4, double>({{1., 2., 3., 10.}, {2., 3., 4., 20.}, {4., -1., 10., 50.}, {10., -11., 12., -13.}});

	CHECK(m.det() == Approx(2248.).margin(1e-9));
}

TEST_CASE("matrixView_det44_dyn2","") {
	auto m = SiLi::make_mat<double>({{1., 2., 3., 10.}, {2., 3., 4., 20.}, {4., -1., 10., 50.}, {10., -11., 12., -13.}});

	CHECK(m.det() == Approx(2248.).margin(1e-9));
}


TEST_CASE("matrixView_det44_inv33","") {
	auto m = SiLi::make_mat<3, 3, double>({{1., 0., 0.},
	                         {0., 1., 0.},
	                         {0., 0., 1.}});

	std::vector<double> v {1, 0, 0, 0, 1, 0, 0, 0, 1};
	std::vector<double> list;
	for (auto x : inv(m)) {
		list.push_back(x);
	}
	REQUIRE(v.size() == list.size());
	for (size_t i(0); i < v.size(); ++i) {
		CHECK(v[i] == list[i]);
	}
}

TEST_CASE("matrixView_inv44","") {
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
	REQUIRE(v.size() == list.size());
	for (size_t i(0); i < v.size(); ++i) {
		CHECK(v[i] == Approx(list[i]).margin(1e-9));
	}
}

TEST_CASE("matrixView_inv44_dyn","") {
	SiLi::Matrix<-1, -1, double> m;
	m = SiLi::make_mat<4, 4, double>({{1., 2., 3., 4.},
	                         {6., 6., 7., 8.},
	                         {9., 10., 10., 12.},
	                         {13., 14., 15., 15.}});


	std::vector<double> v {-13./19, 17./19, -2./19, -4./19, 7./38, -53./38, 23./38, 4/19.,
	                       7./19, 4./19, -15./19, 8./19, 1./19, 6./19, 6./19, -7./19};
	std::vector<double> list;
	for (auto x : inv(m)) {
		list.push_back(x);
	}
	REQUIRE(v.size() == list.size());
	for (size_t i(0); i < v.size(); ++i) {
		CHECK(v[i] == Approx(list[i]).margin(1e-9));
	}
}

TEST_CASE("matrixView_pinv44","") {
	auto m = SiLi::make_mat<4, 4, double>({{1., 2., 3., 4.},
	                         {6., 6., 7., 8.},
	                         {9., 10., 10., 12.},
	                         {13., 14., 15., 15.}});


	std::vector<double> v {-13./19, 17./19, -2./19, -4./19, 7./38, -53./38, 23./38, 4/19.,
	                       7./19, 4./19, -15./19, 8./19, 1./19, 6./19, 6./19, -7./19};
	std::vector<double> list;
	for (auto x : m.pinv()) {
		list.push_back(x);
	}
	REQUIRE(v.size() == list.size());
	for (size_t i(0); i < v.size(); ++i) {
		CHECK(v[i] == Approx(list[i]).margin(1e-9));
	}
}

TEST_CASE("matrixView_pinv44_dyn","") {
	auto m = SiLi::make_mat<double>({{1., 2., 3., 4.},
	                         {6., 6., 7., 8.},
	                         {9., 10., 10., 12.},
	                         {13., 14., 15., 15.}});


	std::vector<double> v {-13./19, 17./19, -2./19, -4./19, 7./38, -53./38, 23./38, 4/19.,
	                       7./19, 4./19, -15./19, 8./19, 1./19, 6./19, 6./19, -7./19};
	std::vector<double> list;
	for (auto x : m.pinv()) {
		list.push_back(x);
	}
	REQUIRE(v.size() == list.size());
	for (size_t i(0); i < v.size(); ++i) {
		CHECK(v[i] == Approx(list[i]).margin(1e-9));
	}
}

TEST_CASE("matrixView_fast_inv33","") {
	auto m = SiLi::make_mat<3, 3, double>({{1., 0., 0.},
	                         {0., 1., 0.},
	                         {0., 0., 1.}});

	std::vector<double> v {1, 0, 0, 0, 1, 0, 0, 0, 1};
	std::vector<double> list;
	for (auto x : fastInv(m)) {
		list.push_back(x);
	}
	REQUIRE(v.size() == list.size());
	for (size_t i(0); i < v.size(); ++i) {
		CHECK(v[i] == list[i]);
	}
}

TEST_CASE("matrixView_fast_inv44","") {
	auto m = SiLi::make_mat<4, 4, double>({{1., 2., 3., 4.},
	                         {6., 6., 7., 8.},
	                         {9., 10., 10., 12.},
	                         {13., 14., 15., 15.}});


	std::vector<double> v {-13./19, 17./19, -2./19, -4./19, 7./38, -53./38, 23./38, 4/19.,
	                       7./19, 4./19, -15./19, 8./19, 1./19, 6./19, 6./19, -7./19};
	std::vector<double> list;
	for (auto x : fastInv(m)) {
		list.push_back(x);
	}
	REQUIRE(v.size() == list.size());
	for (size_t i(0); i < v.size(); ++i) {
		CHECK(v[i] == Approx(list[i]).margin(1e-9));
	}
}

TEST_CASE("matrixView_fast_inv33_dyn","") {
	auto m = SiLi::make_mat<double>({{1., 0., 0.},
	                         {0., 1., 0.},
	                         {0., 0., 1.}});

	std::vector<double> v {1, 0, 0, 0, 1, 0, 0, 0, 1};
	std::vector<double> list;
	for (auto x : fastInv(m)) {
		list.push_back(x);
	}
	REQUIRE(v.size() == list.size());
	for (size_t i(0); i < v.size(); ++i) {
		CHECK(v[i] == list[i]);
	}
}

TEST_CASE("matrixView_fast_inv44_dyn","") {
	auto m = SiLi::make_mat<double>({{1., 2., 3., 4.},
	                         {6., 6., 7., 8.},
	                         {9., 10., 10., 12.},
	                         {13., 14., 15., 15.}});


	std::vector<double> v {-13./19, 17./19, -2./19, -4./19, 7./38, -53./38, 23./38, 4/19.,
	                       7./19, 4./19, -15./19, 8./19, 1./19, 6./19, 6./19, -7./19};
	std::vector<double> list;
	for (auto x : fastInv(m)) {
		list.push_back(x);
	}
	REQUIRE(v.size() == list.size());
	for (size_t i(0); i < v.size(); ++i) {
		CHECK(v[i] == Approx(list[i]).margin(1e-9));
	}
}


TEST_CASE("matrixView_svd2x2","") {
	auto m = SiLi::make_mat<2, 2, double>({{1., 2.}, {2., 3.}});

	auto svd = m.svd();
	auto S = make_diag(svd.S);

	auto re = svd.U * S * svd.V.t_view();
	CHECK(std::abs(svd.U.det()) == Approx(1.).margin(1e-9));
	CHECK(std::abs(svd.V.det()) == Approx(1.).margin(1e-9));
	CHECK(m(0, 0) == Approx(re(0, 0)).margin(1e-9));
	CHECK(m(1, 0) == Approx(re(1, 0)).margin(1e-9));
	CHECK(m(0, 1) == Approx(re(0, 1)).margin(1e-9));
	CHECK(m(1, 1) == Approx(re(1, 1)).margin(1e-9));
}

TEST_CASE("matrixView_svd2x2_dyn","") {
	auto m = SiLi::make_mat<double>({{1., 2.}, {2., 3.}});

	auto svd = m.svd();
	auto S = make_diag(svd.S);


	auto re = svd.U * S * svd.V.t_view();
	CHECK(std::abs(svd.U.det()) == Approx(1.).margin(1e-9));
	CHECK(std::abs(svd.V.det()) == Approx(1.).margin(1e-9));
	CHECK(m(0, 0) == Approx(re(0, 0)).margin(1e-9));
	CHECK(m(1, 0) == Approx(re(1, 0)).margin(1e-9));
	CHECK(m(0, 1) == Approx(re(0, 1)).margin(1e-9));
	CHECK(m(1, 1) == Approx(re(1, 1)).margin(1e-9));
}

TEST_CASE("matrixView_svd3x3","") {
	auto m = SiLi::make_mat<3, 3, double>({{1., 2., 3.}, {2., 3., 4.}, {4., 1., 10.}});

	auto svd = m.svd();
	auto S = make_diag(svd.S);

	auto re = svd.U * S * svd.V.t_view();
	CHECK(std::abs(svd.U.det()) == Approx(1.).margin(1e-9));
	CHECK(std::abs(svd.V.det()) == Approx(1.).margin(1e-9));
	CHECK(m(0, 0) == Approx(re(0, 0)).margin(1e-9));
	CHECK(m(1, 0) == Approx(re(1, 0)).margin(1e-9));
	CHECK(m(2, 0) == Approx(re(2, 0)).margin(1e-9));
	CHECK(m(0, 1) == Approx(re(0, 1)).margin(1e-9));
	CHECK(m(1, 1) == Approx(re(1, 1)).margin(1e-9));
	CHECK(m(2, 1) == Approx(re(2, 1)).margin(1e-9));
	CHECK(m(0, 2) == Approx(re(0, 2)).margin(1e-9));
	CHECK(m(1, 2) == Approx(re(1, 2)).margin(1e-9));
	CHECK(m(2, 2) == Approx(re(2, 2)).margin(1e-9));
}

TEST_CASE("matrixView_svd4x4", "") {
	auto m = SiLi::make_mat<4, 4, double>({{1., 2., 3., 10.}, {2., 3., 4., 20.}, {4., -1., 10., 50.}, {10., -11., 12., -13.}});

	auto svd = m.svd();
	auto S = make_diag(svd.S);

	auto re = svd.U * S * svd.V.t_view();
	CHECK(std::abs(svd.U.det()) == Approx(1.).margin(1e-9));
	CHECK(std::abs(svd.V.det()) == Approx(1.).margin(1e-9));
	for (int row(0); row<m.num_rows(); ++row) {
		for (int col(0); col<m.num_cols(); ++col) {
			CHECK(m(row, col) == Approx(re(row, col)).margin(1e-9));
		}
	}
}

TEST_CASE("matrixView_svd5x5","") {
	auto m = SiLi::make_mat<5, 5, double>({{1., 2., 3., 10., -100.}, {2., 3., 4., 20., 101.}, {4., -1., 10., 50., -102.}, {10., -11., 12., -13., -103.}, {104., -105., 106., 107., -108.}});

	auto svd = m.svd();
	auto S = make_diag(svd.S);

	auto re = svd.U * S * svd.V.t_view();
	CHECK(std::abs(svd.U.det()) == Approx(1.).margin(1e-9));
	CHECK(std::abs(svd.V.det()) == Approx(1.).margin(1e-9));
	for (int row(0); row<m.num_rows(); ++row) {
		for (int col(0); col<m.num_cols(); ++col) {
			CHECK(m(row, col) == Approx(re(row, col)).margin(1e-9));
		}
	}
}

TEST_CASE("svd_random","") {
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
	CAPTURE(m);
	CAPTURE(re);
	CAPTURE(diff);
	CHECK(0. == Approx(diff.norm()).margin(1e-9));
	auto Udet = svd.U.view<4, 4>(0, 0).det();
	CHECK(1. == Approx(std::abs(Udet)).margin(1e-9));
	CHECK(1. == Approx(std::abs(svd.V.det())).margin(1e-9));


	for (int row(0); row<U_true.num_rows(); ++row) {
		for (int col(0); col<U_true.num_cols(); ++col) {
			CHECK(std::abs(svd.U(row, col)) == Approx(std::abs(U_true(row, col))).margin(1e-5));
		}
	}
	for (int row(0); row<S_true.num_rows(); ++row) {
		for (int col(0); col<S_true.num_cols(); ++col) {
			CHECK(std::abs(S(row, col)) == Approx(std::abs(S_true(row, col))).margin(1e-5));
		}
	}
	for (int row(0); row<V_true.num_rows(); ++row) {
		for (int col(0); col<V_true.num_cols(); ++col) {
			CHECK(std::abs(svd.V(row, col)) == Approx(std::abs(V_true(row, col))).margin(1e-5));
		}
	}
}

TEST_CASE("svd_random_online_sqr","") {
	for (int i(0); i<1000; ++i) {
		auto v = []() { return (double(rand())/RAND_MAX) * 20. - 10.; };
		auto m = SiLi::make_mat<4, 4, double>({
				 {v(), v(), v(), v()},
				 {v(), v(), v(), v()},
				 {v(), v(), v(), v()},
				 {v(), v(), v(), v()}
		});

		auto svd = m.svd();
		auto S = make_diag(svd.S);

		auto re = svd.U * S * svd.V.t_view();
		auto diff = m - re;;
		CHECK(0. == Approx(diff.norm()).margin(1e-9));
		auto Udet = svd.U.view<4, 4>(0, 0).det();
		CHECK(1. == Approx(std::abs(Udet)).margin(1e-9));
		CHECK(1. == Approx(std::abs(svd.V.det())).margin(1e-9));
	}
}

TEST_CASE("svd_random_online","") {
	for (int i(0); i<100; ++i) {
		auto v = []() { return std::round(((double(rand())/RAND_MAX) * 20. - 10.)*100.)/100.; };
		auto m = SiLi::make_mat<4, 5, double>({
				 {v(), v(), v(), v(), v()},
				 {v(), v(), v(), v(), v()},
				 {v(), v(), v(), v(), v()},
				 {v(), v(), v(), v(), v()}
		}).t();

		auto svd = m.svd();
		auto S = make_diag(svd.S);

		auto re = svd.U * S * svd.V.t_view();
		auto diff = m - re;;
		CHECK(0. == Approx(diff.norm()).margin(1e-9));
		//!TODO
/*		auto Udet = svd.U.view<4, 4>(0, 0).det();
		CHECK(1. == Approx(std::abs(Udet)));.margin(1e-9)*/
		CHECK(1. == Approx(std::abs(svd.V.det())).margin(1e-9));
	}
}

TEST_CASE("matrxView_pow1","") {
	auto m = SiLi::make_mat<3, 3, double>({{1., 0., 0.},
	                         {0., 1., 0.},
	                         {0., 0., 1.}});

	std::vector<double> v {1, 0, 0, 0, 1, 0, 0, 0, 1};
	std::vector<double> list;
	for (auto x : (m^1.5)) {
		list.push_back(x);
	}
	REQUIRE(v.size() == list.size());
	for (size_t i(0); i < v.size(); ++i) {
		CHECK(v[i] == list[i]);
	}
}

TEST_CASE("matrxView_pow2","") {
	auto m = SiLi::make_mat<3, 3, double>({{1., 0., 0.},
	                         {0., 4., 0.},
	                         {0., 0., 1.}});

	std::vector<double> v {1, 0, 0, 0, 2, 0, 0, 0, 1};
	std::vector<double> list;
	for (auto x : (m^0.5)) {
		list.push_back(x);
	}
	REQUIRE(v.size() == list.size());
	for (size_t i(0); i < v.size(); ++i) {
		CHECK(v[i] == list[i]);
	}
}

TEST_CASE("matrxView_pow3","") {
	auto m = SiLi::make_mat<3, 3, double>({{1., 0., 0.},
	                         {0., 4., 0.},
	                         {0., 0., 1.}});

	std::vector<double> v {1, 0, 0, 0, 0.25, 0, 0, 0, 1};
	std::vector<double> list;
	for (auto x : (m^(-1.))) {
		list.push_back(x);
	}
	REQUIRE(v.size() == list.size());
	for (size_t i(0); i < v.size(); ++i) {
		CHECK(v[i] == list[i]);
	}
}

TEST_CASE("matrxView_rank1","") {
	auto m = SiLi::make_mat<3, 3, double>({{1., 0., 0.},
	                         {0., 4., 0.},
	                         {0., 0., 1.}});

	CHECK(m.rank() == 3);
}

TEST_CASE("matrxView_rank2","") {
	auto m = SiLi::make_mat<3, 3, double>({{1., 0., 0.},
	                         {0., 4., 0.},
	                         {0., 0., 0.}});

	CHECK(m.rank() == 2);
}


TEST_CASE("matrxView_transposed","") {
	auto m = SiLi::make_mat<2, 2, double>({{11., 12.},
	                         {21., 22.}});
	auto view = m.t_view();
	CHECK(view.det() == Approx(-10.).margin(1e-9));
	CHECK(view(0, 0) == 11.);
	CHECK(view(1, 0) == 12.);
	CHECK(view(0, 1) == 21.);
	CHECK(view(1, 1) == 22.);
}

TEST_CASE("matrxView_transposed2","") {
	auto m = SiLi::make_mat<3, 2, double>({{11., 12.},
	                         {21., 22.},
	                         {31., 32.}});
	auto view = m.t_view();
	CHECK(view.num_rows() == 2);
	CHECK(view.num_cols() == 3);
	CHECK(view(0, 0) == 11.);
	CHECK(view(1, 0) == 12.);
	CHECK(view(0, 1) == 21.);
	CHECK(view(1, 1) == 22.);
	CHECK(view(0, 2) == 31.);
	CHECK(view(1, 2) == 32.);
}

TEST_CASE("matrixView_transposed_transposed", "") {
	auto m = SiLi::make_mat<2, 2, double>({{11., 12.},
	                         {21., 22.}});
	auto view = m.t_view().t_view();
	CHECK(view.det() == Approx(-10.).margin(1e-9));
	CHECK(view(0, 0) == 11.);
	CHECK(view(1, 0) == 21.);
	CHECK(view(0, 1) == 12.);
	CHECK(view(1, 1) == 22.);
}

TEST_CASE("matrixView_implicit_cast_to_T", "") {
	auto m = SiLi::make_mat<1, 1, double>({{11.}});

	double x = m;
	CHECK(x == 11.);
}

TEST_CASE("matrixView_add", "") {
	auto m1 = SiLi::make_mat<2, 2, double>({{11., 12.},
	                          {21., 22.}});
	auto m2 = SiLi::make_mat<2, 2, double>({{111., 112.},
	                          {121., 122.}});

	auto m3 = m1 + m2;
	CHECK(m3(0, 0) == 122.);
	CHECK(m3(1, 0) == 142.);
	CHECK(m3(0, 1) == 124.);
	CHECK(m3(1, 1) == 144.);
}

TEST_CASE("matrixView_sub", "") {
	auto m1 = SiLi::make_mat<2, 2, double>({{11., 12.},
	                          {21., 22.}});
	auto m2 = SiLi::make_mat<2, 2, double>({{111., 112.},
	                          {121., 122.}});

	auto m3 = m2 - m1;
	CHECK(m3(0, 0) == 100.);
	CHECK(m3(1, 0) == 100.);
	CHECK(m3(0, 1) == 100.);
	CHECK(m3(1, 1) == 100.);
}

TEST_CASE("matrixView_scale", "") {
	auto m1 = SiLi::make_mat<2, 2, double>({{11., 12.},
	                          {21., 22.}});

	auto m3 = m1 * 2.;
	CHECK(m3(0, 0) == 22.);
	CHECK(m3(1, 0) == 42.);
	CHECK(m3(0, 1) == 24.);
	CHECK(m3(1, 1) == 44.);
}
TEST_CASE("matrixView_neg", "") {
	auto m1 = SiLi::make_mat<2, 2, double>({{11., 12.},
	                          {21., 22.}});

	auto m3 = -m1;
	CHECK(m3(0, 0) == -11.);
	CHECK(m3(1, 0) == -21.);
	CHECK(m3(0, 1) == -12.);
	CHECK(m3(1, 1) == -22.);
}
TEST_CASE("matrixView_diag", "") {
	auto m = SiLi::make_mat<3, 2, double>({{11., 12.},
	                         {21., 22.},
	                         {31., 32.}});

	CHECK(11. == m.diag()(0, 0));
	CHECK(22. == m.diag()(1, 0));

	CHECK(11. == m.diag().t_view()(0, 0));
	CHECK(22. == m.diag().t_view()(0, 1));

	CHECK(11. == m.diag().t_view().diag()(0, 0));
	CHECK( 1 == m.diag().t_view().diag().num_rows());
	CHECK( 1 == m.diag().t_view().diag().num_cols());

	CHECK(11. == m.t_view().diag()(0, 0));
	CHECK(22. == m.t_view().diag()(1, 0));
}

TEST_CASE("matrixView_diag_dyn", "") {
	SiLi::Matrix<-1, -1, double> m;
	m = SiLi::make_mat<3, 2, double>({{11., 12.},
	                         {21., 22.},
	                         {31., 32.}});

	CHECK(11. == m.diag()(0, 0));
	CHECK(22. == m.diag()(1, 0));

	CHECK(11. == m.diag().t_view()(0, 0));
	CHECK(22. == m.diag().t_view()(0, 1));

	CHECK(11. == m.diag().t_view().diag()(0, 0));
	CHECK( 1 == m.diag().t_view().diag().num_rows());
	CHECK( 1 == m.diag().t_view().diag().num_cols());

	CHECK(11. == m.t_view().diag()(0, 0));
	CHECK(22. == m.t_view().diag()(1, 0));
}

TEST_CASE("matrixView_diag_write_dyn", "") {
	SiLi::Matrix<-1, -1, double> m;
	m = SiLi::make_mat<3, 2, double>({{11., 12.},
	                         {21., 22.},
	                         {31., 32.}});

	m.diag()(0) = 1.;
	m.diag()(1) = 2.;

	CHECK(1. == m(0, 0));
	CHECK(2. == m(1, 1));

	m.diag() += 1.;

	CHECK(2. == m(0, 0));
	CHECK(3. == m(1, 1));

}



TEST_CASE("matrixView_diag2", "") {
	auto m = SiLi::make_mat<3, 2, double>({{11., 12.},
	                         {21., 22.},
	                         {31., 32.}});

	CHECK(11. == m.diag()(0));
	CHECK(22. == m.diag()(1));
	CHECK(11. == m.diag().t_view()(0));
	CHECK(22. == m.diag().t_view()(1));
}

TEST_CASE("matrixView_make_diag", "") {
	auto v = SiLi::make_mat<3, 1, double>({{11.},
	                         {21.},
	                         {31.}});
	auto m = make_diag(v);

	CHECK(11. == m(0, 0));
	CHECK(0. ==  m(0, 1));
	CHECK(0. ==  m(0, 2));
	CHECK(0. ==  m(1, 0));
	CHECK(21. == m(1, 1));
	CHECK(0. ==  m(1, 2));
	CHECK(0. ==  m(2, 0));
	CHECK(0. ==  m(2, 1));
	CHECK(31. == m(2, 2));
}

TEST_CASE("matrixView_make_diag2", "") {
	auto v = SiLi::make_mat<3, 1, double>({{11.},
	                         {21.},
	                         {31.}});
	auto m = make_diag<3, 4>(v);

	CHECK(11. == m(0, 0));
	CHECK(0. ==  m(0, 1));
	CHECK(0. ==  m(0, 2));
	CHECK(0. ==  m(0, 3));
	CHECK(0. ==  m(1, 0));
	CHECK(21. == m(1, 1));
	CHECK(0. ==  m(1, 2));
	CHECK(0. ==  m(1, 3));
	CHECK(0. ==  m(2, 0));
	CHECK(0. ==  m(2, 1));
	CHECK(31. == m(2, 2));
	CHECK(0. ==  m(2, 3));
}

TEST_CASE("matrixView_make_diag3", "") {
	auto v = SiLi::make_mat<3, 1, double>({{11.},
	                         {21.},
	                         {31.}});
	auto m = make_diag<4, 3>(v);

	CHECK(11. == m(0, 0));
	CHECK(0. ==  m(0, 1));
	CHECK(0. ==  m(0, 2));
	CHECK(0. ==  m(1, 0));
	CHECK(21. == m(1, 1));
	CHECK(0. ==  m(1, 2));
	CHECK(0. ==  m(2, 0));
	CHECK(0. ==  m(2, 1));
	CHECK(31. == m(2, 2));
	CHECK(0. ==  m(3, 0));
	CHECK(0. ==  m(3, 1));
	CHECK(0. ==  m(3, 2));
}

TEST_CASE("matrixView_element_wise_mul", "") {
	auto v = SiLi::make_mat<3, 1, double>({{11.},
	                         {21.},
	                         {31.}});
	auto m = make_eye<3, 3, double>();
	m.diag() &= v;

	CHECK(11. == m(0, 0));
	CHECK(0. ==  m(0, 1));
	CHECK(0. ==  m(0, 2));
	CHECK(0. ==  m(1, 0));
	CHECK(21. == m(1, 1));
	CHECK(0. ==  m(1, 2));
	CHECK(0. ==  m(2, 0));
	CHECK(0. ==  m(2, 1));
	CHECK(31. == m(2, 2));
}
TEST_CASE("matrixView_element_wise_mul2", "") {
	auto v = SiLi::make_mat<3, 1, double>({{11.},
	                         {21.},
	                         {31.}});
	auto m = make_eye<3, 3, double>();
	m.diag() = m.diag() & v;

	CHECK(11. == m(0, 0));
	CHECK(0. ==  m(0, 1));
	CHECK(0. ==  m(0, 2));
	CHECK(0. ==  m(1, 0));
	CHECK(21. == m(1, 1));
	CHECK(0. ==  m(1, 2));
	CHECK(0. ==  m(2, 0));
	CHECK(0. ==  m(2, 1));
	CHECK(31. == m(2, 2));
}

TEST_CASE("matrixView_swap", "") {
	auto m = SiLi::make_mat<3, 2, double>({{11., 12.},
	                         {21., 22.},
	                         {31., 32.}});

	auto view1 = m.view<1, 2>(0, 0);
	auto view2 = m.view<1, 2>(2, 0);

	swap(view1, view2);
	CHECK(31. == m(0, 0));
	CHECK(32. == m(0, 1));
	CHECK(21. == m(1, 0));
	CHECK(22. == m(1, 1));
	CHECK(11. == m(2, 0));
	CHECK(12. == m(2, 1));
}

TEST_CASE("matrixView_rows", "") {
	auto m = SiLi::make_mat<3, 2, double>({{11., 12.},
	                         {21., 22.},
	                         {31., 32.}});

	auto rows = m.rows();

	CHECK(11. == rows[0](0, 0));
	CHECK(12. == rows[0](0, 1));
	CHECK(21. == rows[1](0, 0));
	CHECK(22. == rows[1](0, 1));
	CHECK(31. == rows[2](0, 0));
	CHECK(32. == rows[2](0, 1));
}

TEST_CASE("matrixView_cols", "") {
	auto m = SiLi::make_mat<3, 2, double>({{11., 12.},
	                         {21., 22.},
	                         {31., 32.}});

	auto cols = m.cols();

	CHECK(11. == cols[0](0, 0));
	CHECK(12. == cols[1](0, 0));
	CHECK(21. == cols[0](1, 0));
	CHECK(22. == cols[1](1, 0));
	CHECK(31. == cols[0](2, 0));
	CHECK(32. == cols[1](2, 0));
}
TEST_CASE("matrixView_single_element_to_double", "") {
	auto m = SiLi::make_mat<3, 1, double>({{1.},
	                         {2.},
	                         {3.}});

	double s = m.t_view() * m;
	CHECK(14. == s);
}

TEST_CASE("matrixView_abs", "") {
	auto m = SiLi::make_mat<3, 4, double>({{-1., 2., -3., 4.},
	                         {6., -6., 7., -8.},
	                         {9., 10., -10., 12.}});
	auto a = abs(m);
	CHECK(a(0, 0) == 1.);
	CHECK(a(0, 1) == 2.);
	CHECK(a(0, 2) == 3.);
	CHECK(a(0, 3) == 4.);
	CHECK(a(1, 0) == 6.);
	CHECK(a(1, 1) == 6.);
	CHECK(a(1, 2) == 7.);
	CHECK(a(1, 3) == 8.);
	CHECK(a(2, 0) == 9.);
	CHECK(a(2, 1) == 10.);
	CHECK(a(2, 2) == 10.);
	CHECK(a(2, 3) == 12.);
}

TEST_CASE("matrixView_isfinite_true", "") {
	auto m = SiLi::make_mat<3, 4, double>({{-1., 2., -3., 4.},
	                         {6., -6., 7., -8.},
	                         {9., 10., -10., 12.}});
	CHECK(isfinite(m) == true);
}
TEST_CASE("matrixView_isfinite_false", "") {
	auto m = SiLi::make_mat<3, 4, double>({{-1., 2., -3., 4.},
	                         {6., -6./0., 7., -8.},
	                         {9., 10., -10., 12.}});
	CHECK(isfinite(m) == false);
}

TEST_CASE("matrixView_join_rows", "") {
	M<2, 1>  m1({{11.}, {21.}});
	M<2, 2>  m2({{111., 112.}, {121., 122.}});

	auto m3 = m1.join_rows(m2);
	CHECK(11. == m3(0, 0));
	CHECK(21. == m3(1, 0));
	CHECK(111. == m3(0, 1));
	CHECK(112. == m3(0, 2));
	CHECK(121. == m3(1, 1));
	CHECK(122. == m3(1, 2));
}

TEST_CASE("matrixView_join_cols", "") {
	M<1, 2>  m1({{11., 21.}});
	M<2, 2>  m2({{111., 112.}, {121., 122.}});

	auto m3 = m1.join_cols(m2);
	CHECK(11. == m3(0, 0));
	CHECK(21. == m3(0, 1));
	CHECK(111. == m3(1, 0));
	CHECK(112. == m3(1, 1));
	CHECK(121. == m3(2, 0));
	CHECK(122. == m3(2, 1));
}


TEST_CASE("matrixView_rowRange1", "") {
	M<4, 2>  m({
			{00., 01.},
			{10., 11.},
			{20., 21.},
			{30., 31.}
	});

	auto rows = m.rows<0, 3, 2>();
	CHECK(rows.size() == 2);

	CHECK(rows[0](0, 0) == 00.);
	CHECK(rows[0](0, 1) == 01.);

	CHECK(rows[1](0, 0) == 20.);
	CHECK(rows[1](0, 1) == 21.);
}

TEST_CASE("matrixView_rowRange2", "") {
	M<4, 2>  m({
			{00., 01.},
			{10., 11.},
			{20., 21.},
			{30., 31.}
	});

	auto rows = m.rows<3, 0, -2>();
	CHECK(rows.size() == 2);

	CHECK(rows[0](0, 0) == 30.);
	CHECK(rows[0](0, 1) == 31.);

	CHECK(rows[1](0, 0) == 10.);
	CHECK(rows[1](0, 1) == 11.);

}


TEST_CASE("matrixView_colRange1", "") {
	M<2, 4>  m({
			{00., 01. , 02., 03.},
			{10., 11. , 12., 13.}
	});

	auto cols = m.cols<0, 3, 2>();
	CHECK(cols.size() == 2);

	CHECK(cols[0](0, 0) == 00.);
	CHECK(cols[0](1, 0) == 10.);

	CHECK(cols[1](0, 0) == 02.);
	CHECK(cols[1](1, 0) == 12.);
}

TEST_CASE("matrixView_colRange2", "") {
	M<2, 4>  m({
			{00., 01. , 02., 03.},
			{10., 11. , 12., 13.}
	});

	auto cols = m.cols<3, 0, -2>();
	CHECK(cols.size() == 2);

	CHECK(cols[0](0, 0) == 03.);
	CHECK(cols[0](1, 0) == 13.);

	CHECK(cols[1](0, 0) == 01.);
	CHECK(cols[1](1, 0) == 11.);
}

TEST_CASE("inv10x10", "") {
	SiLi::Matrix<10, 10, double> mat{0};
	mat.diag() = 0.1;

	auto inv = mat.inv();

	for (int row(0); row < mat.num_rows(); ++row) {
		for (int col(0); col < mat.num_cols(); ++col) {
			if (col == row) {
				CHECK(inv(row, col) == Approx(10.).margin(1e-9));
			} else {
				CHECK(inv(row, col) == Approx(0.).margin(1e-9));
			}
		}
	}

/*
	std::vector<double> v {1, 0, 0, 0, 1, 0, 0, 0, 1};
	std::vector<double> list;
	for (auto x : inv(m)) {
		list.push_back(x);
	}
	REQUIRE(v.size() == list.size());
	for (size_t i(0); i < v.size(); ++i) {
		CHECK(v[i], list[i]);
	}*/
}
