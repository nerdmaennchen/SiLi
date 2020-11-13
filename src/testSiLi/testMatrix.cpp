#include <SiLi/SiLi.h>
#include <catch2/catch.hpp>

using namespace SiLi;

TEST_CASE("make_mat", "[init]") {
	SECTION("init 3x1") {
		constexpr auto m = SiLi::Matrix<3, 1, double>{2., 3., 4.}; // Critical
		static_assert(std::is_same_v<decltype(m)::value_t, double>);
		static_assert(3 == decltype(m)::Rows);
		static_assert(1 == decltype(m)::Cols);
	}

	SECTION("init 1x2") {
		constexpr auto m = SiLi::Matrix{{{2., 3.}}}; // Critical
		static_assert(std::is_same_v<decltype(m)::value_t, double>);
		static_assert(1 == decltype(m)::Rows);
		static_assert(2 == decltype(m)::Cols);
		static_assert(2. == m(0, 0));
		static_assert(3. == m(0, 1));
	}

	SECTION("init 2x2") {
		constexpr auto m = SiLi::Matrix{{{2., 3.},
		                                  {0., 4.}}}; // Critical
		static_assert(std::is_same_v<decltype(m)::value_t, double>);
		static_assert(2 == decltype(m)::Rows);
		static_assert(2 == decltype(m)::Cols);
		static_assert(2. == m(0, 0));
		static_assert(3. == m(0, 1));
		static_assert(0. == m(1, 0));
		static_assert(4. == m(1, 1));
	}

	SECTION("addition") {
		constexpr auto m1 = SiLi::Matrix{{{2., 3.}}};
		constexpr auto m2 = SiLi::Matrix{{{4., 5.}}};
		constexpr auto z  = m1 + m2; // Critical

		static_assert(std::is_same_v<decltype(z)::value_t, double>);
		static_assert(1 == decltype(z)::Rows);
		static_assert(2 == decltype(z)::Cols);
		static_assert(6. == z(0, 0));
		static_assert(8. == z(0, 1));
	}

	SECTION("addition with different type") {
		constexpr auto m1 = SiLi::Matrix{{{2., 3.}}};
		constexpr auto m2 = SiLi::Matrix{{{4, 5}}};
		constexpr auto z  = m1 + m2; // Critical
		static_assert(std::is_same_v<decltype(m1)::value_t, double>);
		static_assert(std::is_same_v<decltype(m2)::value_t, int>);

		static_assert(std::is_same_v<decltype(z)::value_t, double>);
		static_assert(1 == decltype(z)::Rows);
		static_assert(2 == decltype(z)::Cols);
		static_assert(6. == z(0, 0));
		static_assert(8. == z(0, 1));
	}


	SECTION("addition - inplace") {
		constexpr auto m1 = SiLi::Matrix{{{2., 3.}}};
		constexpr auto m2 = SiLi::Matrix{{{4., 5.}}};
		constexpr auto z  = [=](auto x) {
			return x += m2; // Critical
		}(m1);

		static_assert(std::is_same_v<decltype(z)::value_t, double>);
		static_assert(1 == decltype(z)::Rows);
		static_assert(2 == decltype(z)::Cols);
		static_assert(6. == z(0, 0));
		static_assert(8. == z(0, 1));
	}

	SECTION("addition - inplace with different type") {
		constexpr auto m1 = SiLi::Matrix{{{2., 3.}}};
		constexpr auto m2 = SiLi::Matrix{{{4, 5}}};
		constexpr auto z  = [=](auto x) {
			return x += m2; // Critical
		}(m1);

		static_assert(std::is_same_v<decltype(z)::value_t, double>);
		static_assert(1 == decltype(z)::Rows);
		static_assert(2 == decltype(z)::Cols);
		static_assert(6. == z(0, 0));
		static_assert(8. == z(0, 1));
	}

	SECTION("subsrtaction") {
		constexpr auto m1 = SiLi::Matrix{{{2., 3.}}};
		constexpr auto m2 = SiLi::Matrix{{{4., 5.}}};
		constexpr auto z  = m1 - m2; // Critical

		static_assert(std::is_same_v<decltype(z)::value_t, double>);
		static_assert(1 == decltype(z)::Rows);
		static_assert(2 == decltype(z)::Cols);
		static_assert(-2. == z(0, 0));
		static_assert(-2. == z(0, 1));
	}

	SECTION("substraction with type change") {
		constexpr auto m1 = SiLi::Matrix{{{2., 3.}}};
		constexpr auto m2 = SiLi::Matrix{{{4, 5}}};
		constexpr auto z  = m1 - m2; // Critical
		static_assert(std::is_same_v<decltype(m1)::value_t, double>);
		static_assert(std::is_same_v<decltype(m2)::value_t, int>);

		static_assert(std::is_same_v<decltype(z)::value_t, double>);
		static_assert(1 == decltype(z)::Rows);
		static_assert(2 == decltype(z)::Cols);
		static_assert(-2. == z(0, 0));
		static_assert(-2. == z(0, 1));
	}


	SECTION("substraction - inplace") {
		constexpr auto m1 = SiLi::Matrix{{{2., 3.}}};
		constexpr auto m2 = SiLi::Matrix{{{4., 5.}}};
		constexpr auto z  = [=](auto x) {
			return x -= m2; // Critical
		}(m1);

		static_assert(std::is_same_v<decltype(z)::value_t, double>);
		static_assert(1 == decltype(z)::Rows);
		static_assert(2 == decltype(z)::Cols);
		static_assert(-2. == z(0, 0));
		static_assert(-2. == z(0, 1));
	}

	SECTION("substraction - inplace with different type") {
		constexpr auto m1 = SiLi::Matrix{{{2., 3.}}};
		constexpr auto m2 = SiLi::Matrix{{{4, 5}}};
		constexpr auto z  = [=](auto x) {
			return x -= m2; // Critical
		}(m1);

		static_assert(std::is_same_v<decltype(z)::value_t, double>);
		static_assert(1 == decltype(z)::Rows);
		static_assert(2 == decltype(z)::Cols);
		static_assert(-2. == z(0, 0));
		static_assert(-2. == z(0, 1));
	}

	SECTION("multiplication") {
		constexpr auto m1 = SiLi::Matrix{{{2., 3.}}};
		constexpr auto m2 = SiLi::Matrix{{{4.}, {5.}}};
		constexpr auto z = m1 * m2; // Critical;

		static_assert(std::is_same_v<decltype(z), const SiLi::Matrix<1, 1, double>>);
		static_assert(23. == double{z});
	}

	SECTION("operator==") {
		static constexpr auto l = SiLi::Matrix{{{2., 3.},
		                                        {0., 4.}}};
		static constexpr auto r = SiLi::Matrix{{{7., 8.},
		                                        {1., 5.}}};
		static_assert(not (l == r));
	}
	SECTION("operator==") {
		static constexpr auto l = SiLi::Matrix{{{2., 3.},
		                                        {0., 4.}}};
		static constexpr auto r = SiLi::Matrix{{{2., 3.},
		                                        {0., 4.}}};
		static_assert(l == r);
	}



	SECTION("operator!=") {
		static constexpr auto l = SiLi::Matrix{{{2., 3.},
		                                        {0., 4.}}};
		static constexpr auto r = SiLi::Matrix{{{7., 8.},
		                                        {1., 5.}}};
		static_assert(l != r);
	}
	SECTION("operator!=") {
		static constexpr auto l = SiLi::Matrix{{{2., 3.},
		                                        {0., 4.}}};
		static constexpr auto r = SiLi::Matrix{{{2., 3.},
		                                        {0., 4.}}};
		static_assert(not (l != r));
	}



	SECTION("diag 2x2") {
		static constexpr auto m = SiLi::Matrix{{{2., 3.},
		                                        {0., 4.}}};
		static constexpr auto z = view_diag(m); // Critical
		static_assert(std::is_same_v<decltype(z)::value_t, const double>);
		static_assert(2 == decltype(z)::Rows);
		static_assert(1 == decltype(z)::Cols);
		static_assert(2. == z(0, 0));
		static_assert(4. == z(1, 0));
	}

	SECTION("diag 2x2 - multiplication") {
		static constexpr auto m = SiLi::Matrix{{{2., 3.},
		                                         {0., 4.}}};
		static constexpr auto z = view_diag(m) * 3.; // Critical
		static_assert(std::is_same_v<decltype(z)::value_t, double>);
		static_assert(2 == decltype(z)::Rows);
		static_assert(1 == decltype(z)::Cols);
		static_assert(6.  == z(0, 0));
		static_assert(12. == z(1, 0));
	}

	SECTION("diag 2x2") {
		static constexpr auto m = SiLi::Matrix{{{2., 3.},
		                                         {0., 4.}}};
		static constexpr auto z = diag(m); // Critical
		static_assert(std::is_same_v<decltype(z)::value_t, double>);
		static_assert(2 == decltype(z)::Rows);
		static_assert(1 == decltype(z)::Cols);
		static_assert(2. == z(0, 0));
		static_assert(4. == z(1, 0));
	}


	SECTION("diag 2x2 - view") {
		static constexpr auto m = SiLi::Matrix{{{2., 3.},
		                                         {0., 4.}}};
		static constexpr auto z = [](auto x) {
			view_diag(x) = 7.; // Critical
			return x;
		}(m);
		static_assert(std::is_same_v<decltype(z)::value_t, double>);
		static_assert(2 == decltype(z)::Rows);
		static_assert(2 == decltype(z)::Cols);
		static_assert(7. == z(0, 0));
		static_assert(3. == z(0, 1));
		static_assert(0. == z(1, 0));
		static_assert(7. == z(1, 1));
	}
	SECTION("view 3x3 ") {
		static constexpr auto m = SiLi::Matrix{{{11, 12, 13},
		                                         {21, 22, 23},
		                                         {31, 32, 33}}};

		SECTION("check 0,0-1,1") {
			static constexpr auto z = view<0, 0, 1, 1>(m); // Critical
			static_assert(std::is_same_v<decltype(z)::value_t, const int>);
			static_assert(1 == decltype(z)::Rows);
			static_assert(1 == decltype(z)::Cols);
			static_assert(11 == z(0, 0));
		}

		SECTION("check 1,0-2,1") {
			static constexpr auto z = view<1, 0, 2, 1>(m); // Critical
			static_assert(std::is_same_v<decltype(z)::value_t, const int>);
			static_assert(1 == decltype(z)::Rows);
			static_assert(1 == decltype(z)::Cols);
			static_assert(21 == z(0, 0));
		}
		SECTION("check 0,1-1,2") {
			static constexpr auto z = view<0, 1, 1, 2>(m); // Critical
			static_assert(std::is_same_v<decltype(z)::value_t, const int>);
			static_assert(1 == decltype(z)::Rows);
			static_assert(1 == decltype(z)::Cols);
			static_assert(12 == z(0, 0));
		}
		SECTION("check 0,0-3,3") {
			static constexpr auto z = view<0, 0, 3, 3>(m); // Critical
			static_assert(std::is_same_v<decltype(z)::value_t, const int>);
			static_assert(3 == decltype(z)::Rows);
			static_assert(3 == decltype(z)::Cols);
			static_assert(11 == z(0, 0));
			static_assert(21 == z(1, 0));
			static_assert(31 == z(2, 0));
			static_assert(12 == z(0, 1));
			static_assert(22 == z(1, 1));
			static_assert(32 == z(2, 1));
			static_assert(13 == z(0, 2));
			static_assert(23 == z(1, 2));
			static_assert(33 == z(2, 2));

			SECTION("check 1,0-3,2") {
				static constexpr auto z2 = view<1, 0, 3, 2>(z); // Critical
				static_assert(2 == decltype(z2)::Rows);
				static_assert(2 == decltype(z2)::Cols);
				static_assert(21 == z2(0, 0));
				static_assert(31 == z2(1, 0));
				static_assert(22 == z2(0, 1));
				static_assert(32 == z2(1, 1));
				SECTION("check 1,0-2,2") {
					static constexpr auto z3 = view<1, 0, 2, 2>(z2); // Critical
					static_assert(1 == decltype(z3)::Rows);
					static_assert(2 == decltype(z3)::Cols);
					static_assert(31 == z3(0, 0));
					static_assert(32 == z3(0, 1));
				}
			}
		}
	}
	SECTION("row 3x3") {
		static constexpr auto m = SiLi::Matrix{{{11, 12, 13},
		                                        {21, 22, 23},
		                                        {31, 32, 33}}};
		SECTION("row 0") {
			static constexpr auto z = view_row<0>(m); // Critical
			static_assert(z == Matrix{{{11, 12, 13}}});
		}
		SECTION("row 1") {
			static constexpr auto z = view_row<1>(m); // Critical
			static_assert(z == Matrix{{{21, 22, 23}}});
		}
		SECTION("row 2") {
			static constexpr auto z = view_row<2>(m); // Critical
			static_assert(z == Matrix{{{31, 32, 33}}});
		}
	}
	SECTION("col 3x3") {
		static constexpr auto m = SiLi::Matrix{{{11, 12, 13},
		                                        {21, 22, 23},
		                                        {31, 32, 33}}};
		SECTION("col 0") {
			static constexpr auto z = view_col<0>(m); // Critical	}
			static_assert(z == Matrix{{{11}, {21}, {31}}});
		}
		SECTION("col 1") {
			static constexpr auto z = view_col<1>(m); // Critical
			static_assert(z == Matrix{{{12}, {22}, {32}}});
		}
		SECTION("col 2") {
			static constexpr auto z = view_col<2>(m); // Critical
			static_assert(z == Matrix{{{13}, {23}, {33}}});
		}
	}

	SECTION("join_rows 2x2") {
		static constexpr auto m0 = SiLi::Matrix{{{22., 23.},
		                                         {20., 24.}}};
		static constexpr auto m1 = SiLi::Matrix{{{12., 13., 14.},
		                                         {10., 14., 17.}}};

		static constexpr auto z = join_rows(m0, m1); // Critical
		static_assert(z == SiLi::Matrix{{{22., 23., 12., 13., 14.},
		                                 {20., 24., 10., 14., 17.}}});
	}

	SECTION("join_cols 2x2") {
		static constexpr auto m0 = SiLi::Matrix{{{22., 23.},
		                                         {20., 24.},
		                                         {26., 27.}}};
		static constexpr auto m1 = SiLi::Matrix{{{12., 13.},
		                                         {10., 14.}}};

		static constexpr auto z = join_cols(m0, m1); // Critical
		static_assert(z == SiLi::Matrix{{{22., 23.},
		                                 {20., 24.},
		                                 {26., 27.},
		                                 {12., 13.},
		                                 {10., 14.}}});
	}
	SECTION("det 1x1") {
		static constexpr auto m = SiLi::Matrix{{{4.}}};
		static constexpr auto z = det(m); // Critical
		static_assert(4. == z);
	}

	SECTION("det 2x2") {
		static constexpr auto m = SiLi::Matrix{{{4., 5.},
		                                        {6., 8.}}};
		static constexpr auto z = det(m); // Critical
		static_assert(2. == z);
	}
	SECTION("det 2x2 - zero") {
		static constexpr auto m = SiLi::Matrix{{{4., 5.},
		                                        {4., 5.}}};
		static constexpr auto z = det(m); // Critical
		static_assert(0. == z);
	}
	SECTION("det 2x2 - negative") {
		static constexpr auto m = SiLi::Matrix{{{4.,  5.},
		                                        {6., -8.}}};
		static constexpr auto z = det(m); // Critical
		static_assert(-62. == z);
	}
	SECTION("det 3x3") {
		static constexpr auto m = SiLi::Matrix{{{4.,  5., 6.},
		                                        {6., -8., 9.},
		                                        {3.,  4., 1.}}};
		static constexpr auto z = det(m); // Critical
		static_assert(217. == z);
	}

	SECTION("det 3x3 - zero") {
		static constexpr auto m = SiLi::Matrix{{{4.,  5., 6.},
		                                        {3.,  4., 1.},
		                                        {3.,  4., 1.}}};
		static constexpr auto z = det(m); // Critical
		static_assert(0. == z);
	}


	SECTION("det 3x3 - negative") {
		static constexpr auto m = SiLi::Matrix{{{4., 5., 6.},
		                                        {6., 8., 9.},
		                                        {3., 4., 1.}}};
		static constexpr auto z = det(m); // Critical
		static_assert(-7. == z);
	}
	SECTION("det 4x4") {
		static constexpr auto m = SiLi::Matrix{{{ 1.,   2.,  3.,  10.},
		                                        { 2.,   3.,  4.,  20.},
		                                        { 4.,  -1., 10.,  50.},
		                                        {10., -11., 12., -13.}}};
		static constexpr auto z = det(m); // Critical
		static_assert(2248. == z);
	}

	SECTION("element multiplication") {
		constexpr auto m1 = SiLi::Matrix{{{2., 3.}}};
		constexpr auto m2 = SiLi::Matrix{{{4., 5.}}};
		constexpr auto z  = element_multi(m1, m2); // Critical

		static_assert(std::is_same_v<decltype(z)::value_t, double>);
		static_assert( 1  == decltype(z)::Rows);
		static_assert( 2  == decltype(z)::Cols);
		static_assert( 8. == z(0, 0));
		static_assert(15. == z(0, 1));
	}

	SECTION("element multiplication with different type") {
		constexpr auto m1 = SiLi::Matrix{{{2., 3.}}};
		constexpr auto m2 = SiLi::Matrix{{{4, 5}}};
		constexpr auto z  = element_multi(m1, m2); // Critical
		static_assert(std::is_same_v<decltype(m1)::value_t, double>);
		static_assert(std::is_same_v<decltype(m2)::value_t, int>);

		static_assert(std::is_same_v<decltype(z)::value_t, double>);
		static_assert( 1  == decltype(z)::Rows);
		static_assert( 2  == decltype(z)::Cols);
		static_assert( 8. == z(0, 0));
		static_assert(15. == z(0, 1));
	}

	SECTION("cross") {
		constexpr auto m1 = SiLi::Matrix{{{1.}, {0.}, {0.}}};
		constexpr auto m2 = SiLi::Matrix{{{0.}, {1.}, {0.}}};
		constexpr auto z   = cross(m1, m2); // Critical
		static_assert(z == SiLi::Matrix{{{0.}, {0.}, {1.}}});
	}

	SECTION("inv 2x2") {
		constexpr auto m = SiLi::Matrix{{{11., 12.},
		                                 {31., 32.}}};
		constexpr auto z = inv(m); // Critical
		static_assert(std::abs(std::get<0>(z)) > 0);

		static_assert(std::abs(std::get<1>(z)(0, 0) - (-1.60)) < 1.e-9);
		static_assert(std::abs(std::get<1>(z)(1, 0) - ( 1.55)) < 1.e-9);
		static_assert(std::abs(std::get<1>(z)(0, 1) - ( 0.60)) < 1.e-9);
		static_assert(std::abs(std::get<1>(z)(1, 1) - (-0.55)) < 1.e-9);
	}

	SECTION("inv 3x3") {
		constexpr static auto m = SiLi::Matrix{{{ 2., 1.,  3.},
                                                { 1., 3., -3.},
                                                {-2., 4.,  4.}}};
		constexpr static auto z = inv(m); // Critical
		static_assert(std::abs(std::get<0>(z)) > 0);
		static_assert(std::abs(std::get<1>(z)(0, 0) - ( 0.3000)) < 1.e-9);
		static_assert(std::abs(std::get<1>(z)(1, 0) - ( 0.0250)) < 1.e-9);
		static_assert(std::abs(std::get<1>(z)(2, 0) - ( 0.1250)) < 1.e-9);
		static_assert(std::abs(std::get<1>(z)(0, 1) - ( 0.1000)) < 1.e-9);
		static_assert(std::abs(std::get<1>(z)(1, 1) - ( 0.1750)) < 1.e-9);
		static_assert(std::abs(std::get<1>(z)(2, 1) - (-0.1250)) < 1.e-9);
		static_assert(std::abs(std::get<1>(z)(0, 2) - (-0.1500)) < 1.e-9);
		static_assert(std::abs(std::get<1>(z)(1, 2) - ( 0.1125)) < 1.e-9);
		static_assert(std::abs(std::get<1>(z)(2, 2) - ( 0.0625)) < 1.e-9);
	}

	SECTION("transposed 2x2 - view") {
		static constexpr auto m = SiLi::Matrix{{{2., 3.},
		                                        {0., 4.}}};
		static constexpr auto z = view_trans(m); // Critical
		static_assert(2 == decltype(z)::Rows);
		static_assert(2 == decltype(z)::Cols);
		static_assert(2. == z(0, 0));
		static_assert(0. == z(0, 1));
		static_assert(3. == z(1, 0));
		static_assert(4. == z(1, 1));
	}

	SECTION("transposed 2x2") {
		static constexpr auto m = SiLi::Matrix{{{2., 3.},
		                                        {0., 4.}}};
		static constexpr auto z = trans(m); // Critical
		static_assert(2 == decltype(z)::Rows);
		static_assert(2 == decltype(z)::Cols);
		static_assert(2. == z(0, 0));
		static_assert(0. == z(0, 1));
		static_assert(3. == z(1, 0));
		static_assert(4. == z(1, 1));
	}

	SECTION("transposed^2 2x2 - view") {
		static constexpr auto m = SiLi::Matrix{{{2., 3.},
		                                        {0., 4.}}};
		static constexpr auto z = view_trans(view_trans(m)); // Critical

		static_assert(2 == decltype(z)::Rows);
		static_assert(2 == decltype(z)::Cols);
		static_assert(2. == z(0, 0));
		static_assert(3. == z(0, 1));
		static_assert(0. == z(1, 0));
		static_assert(4. == z(1, 1));
	}

	SECTION("transposed 2x2 - view - diag") {
		static constexpr auto m = SiLi::Matrix{{{2., 3.},
		                                        {0., 4.}}};
		static constexpr auto z0 = view_trans(m); // Critical
		static constexpr auto z  = view_diag(z0); // Critical
		static_assert(2 == decltype(z)::Rows);
		static_assert(1 == decltype(z)::Cols);
		static_assert(2. == z(0));
		static_assert(4. == z(1));
	}
}
