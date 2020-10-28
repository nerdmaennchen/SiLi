#include <SiLi/SiLi.h>
#include <catch2/catch.hpp>

TEST_CASE("check example from docs", "[doc]") {
	SECTION("at (0))") {
		auto m = SiLi::Matrix{{{3, 4},
		                       {5, 6}}};
		at<0, 0>(m) = 5;
		at<0, 1>(m) = 6;
		CHECK(at<0, 0>(m) == 5);
		CHECK(at<0, 1>(m) == 6);
		CHECK(at<1, 0>(m) == 5);
		CHECK(at<1, 1>(m) == 6);
	}

	SECTION("at (0) (0))") {
		auto m = SiLi::Matrix{{{3, 4},
		                       {5, 6}}};
		at<0, 0>(view<0, 0, 2, 2>(m)) = 5;
		at<0, 1>(view<0, 0, 2, 2>(m)) = 6;
		CHECK(at<0, 0>(m) == 5);
		CHECK(at<0, 1>(m) == 6);
		CHECK(at<1, 0>(m) == 5);
		CHECK(at<1, 1>(m) == 6);
	}

	SECTION("at (0) (1))") {
		auto const m = SiLi::Matrix{{{3, 4},
		                             {5, 6}}};
		CHECK(at<0, 0>(view<0, 0, 2, 2>(m)) == 3);
		CHECK(at<0, 1>(view<0, 0, 2, 2>(m)) == 4);
		CHECK(at<1, 0>(view<0, 0, 2, 2>(m)) == 5);
		CHECK(at<1, 1>(view<0, 0, 2, 2>(m)) == 6);
	}



	SECTION("at (1))") {
		auto m = SiLi::Matrix{{{3},
		                       {5}}};
		at<0>(m) = 5;
		at<1>(m) = 6;

		CHECK(at<0>(m) == 5);
		CHECK(at<1>(m) == 6);
	}
	SECTION("rows") {
		auto m = SiLi::Matrix{{{3, 4, 7},
		                       {5, 6, 8}}};
		CHECK(rows(m) == 2);
	}
	SECTION("cols") {
		auto m = SiLi::Matrix{{{3, 4, 7},
		                       {5, 6, 8}}};
		CHECK(cols(m) == 3);
	}
	SECTION("+v") {
		auto m = SiLi::Matrix{{{3, 4, 7},
							   {5, 6, 8}}};
		auto m2 = +m;
		CHECK((m == m2));
	}
	SECTION("a+b") {
		auto a = SiLi::Matrix{{{3, 4, 7},
		                       {5, 6, 8}}};
		auto b = SiLi::Matrix{{{1, 2, 3},
		                       {4, 5, 6}}};
		auto c = a + b;
		auto z = SiLi::Matrix{{{4, 6, 10},
		                       {9, 11, 14}}};
		CHECK((c==z));
	}

	SECTION("a+=b") {
		auto a = SiLi::Matrix{{{3, 4, 7},
		                       {5, 6, 8}}};
		auto b = SiLi::Matrix{{{1, 2, 3},
		                       {4, 5, 6}}};
		a += b;
		auto z = SiLi::Matrix{{{4, 6, 10},
		                       {9, 11, 14}}};
		CHECK((a==z));
	}

	SECTION("-v") {
		auto m = SiLi::Matrix{{{3, 4, 7},
							   {5, 6, 8}}};
		auto m2 = -m;
		auto z = SiLi::Matrix{{{-3, -4, -7},
		                       {-5, -6, -8}}};
		CHECK((m2 == z));
	}
	SECTION("a-b") {
		auto a = SiLi::Matrix{{{3, 4, 7},
		                       {5, 6, 8}}};
		auto b = SiLi::Matrix{{{1, 2, 3},
		                       {4, 5, 6}}};
		auto c = a - b;
		auto z = SiLi::Matrix{{{2,  2, 4},
		                       {1,  1, 2}}};
		CHECK((c==z));
	}

	SECTION("a-=b") {
		auto a = SiLi::Matrix{{{3, 4, 7},
		                       {5, 6, 8}}};
		auto b = SiLi::Matrix{{{1, 2, 3},
		                       {4, 5, 6}}};
		a -= b;
		auto z = SiLi::Matrix{{{2,  2, 4},
		                       {1,  1, 2}}};
		CHECK((a==z));
	}

	SECTION("a*b") {
		auto a = SiLi::Matrix{{{3, 4, 7},
		                       {5, 6, 8}}};
		auto b = SiLi::Matrix{{{1, 2},
		                       {4, 5},
		                       {3, 6}}};
		auto c = a * b;
		auto z = SiLi::Matrix{{{40, 68},
		                       {53, 88}}};
		CHECK((c==z));
	}

	SECTION("a*s") {
		auto a = SiLi::Matrix{{{3, 4, 7},
		                       {5, 6, 8}}};
		auto c = a * 5;
		auto z = SiLi::Matrix{{{15, 20, 35},
		                       {25, 30, 40}}};
		CHECK((c==z));
	}

	SECTION("s*a") {
		auto a = SiLi::Matrix{{{3, 4, 7},
		                       {5, 6, 8}}};
		auto c = 5 * a;
		auto z = SiLi::Matrix{{{15, 20, 35},
		                       {25, 30, 40}}};
		CHECK((c==z));
	}

	SECTION("a*=5") {
		auto a = SiLi::Matrix{{{3, 4, 7},
		                       {5, 6, 8}}};
		a *= 5;
		auto z = SiLi::Matrix{{{15, 20, 35},
		                       {25, 30, 40}}};
		CHECK((a==z));
	}

	SECTION("a/s") {
		auto a = SiLi::Matrix{{{15, 20, 35},
		                       {25, 30, 40}}};

		auto c = a / 5;
		auto z = SiLi::Matrix{{{3, 4, 7},
		                       {5, 6, 8}}};
		CHECK((c==z));
	}

	SECTION("a/=s") {
		auto a = SiLi::Matrix{{{15, 20, 35},
		                       {25, 30, 40}}};

		a /= 5;
		auto z = SiLi::Matrix{{{3, 4, 7},
		                       {5, 6, 8}}};
		CHECK((a==z));
	}

	SECTION("a==b") {
		auto a = SiLi::Matrix{{{15, 20, 35},
		                       {25, 30, 40}}};
		auto b = SiLi::Matrix{{{15, 20, 35},
		                       {25, 30, 40}}};
		auto c = SiLi::Matrix{{{15, 20, 35},
		                       {25, 30, 43}}};
		CHECK((a == b));
		CHECK(not (a == c));
	}

	SECTION("a!=b") {
		auto a = SiLi::Matrix{{{15, 20, 35},
		                       {25, 30, 40}}};
		auto b = SiLi::Matrix{{{15, 20, 35},
		                       {25, 30, 40}}};
		auto c = SiLi::Matrix{{{15, 20, 35},
		                       {25, 30, 43}}};
		CHECK(not (a != b));
		CHECK((a != c));
	}

	SECTION("view_diag") {
		auto a = SiLi::Matrix{{{15, 20, 35},
		                       {25, 30, 40}}};
		auto v = view_diag(a);
		CHECK(rows(v) == 2);
		CHECK(cols(v) == 1);
		CHECK(v(0) == 15);
		CHECK(v(1) == 30);
		v(0) = 7;
		v(1) = 8;
		auto z = SiLi::Matrix{{{ 7, 20, 35},
		                       {25,  8, 40}}};
		CHECK((a == z));
	}

	SECTION("diag") {
		auto a = SiLi::Matrix{{{15, 20, 35},
		                       {25, 30, 40}}};
		auto v = diag(a);
		CHECK(rows(v) == 2);
		CHECK(cols(v) == 1);
		CHECK(v(0) == 15);
		CHECK(v(1) == 30);
	}

	SECTION("view_upper_diag") {
		auto a = SiLi::Matrix{{{15, 20, 35},
		                       {25, 30, 40}}};
		auto v = view_upper_diag(a);
		CHECK(rows(v) == 2);
		CHECK(cols(v) == 1);
		CHECK(v(0) == 20);
		CHECK(v(1) == 40);
		v(0) = 7;
		v(1) = 8;
		auto z = SiLi::Matrix{{{15,  7, 35},
		                       {25, 30,  8}}};
		CHECK((a == z));
	}

	SECTION("view (0)") {
		auto a = SiLi::Matrix{{{ 1,  2,  3,  4,  5},
		                       { 6,  7,  8,  9, 10},
		                       {11, 12, 13, 14, 15},
		                       {16, 17, 18, 19, 20}}};
		auto v0 = view<1, 1, 3, 3>(a);
		auto z0 = SiLi::Matrix{{{ 7,  8},
		                        {12, 13}}};
		CHECK((v0 == z0));

		auto v1 = view<2, 3, SiLi::End, SiLi::End>(a);
		auto z1 = SiLi::Matrix{{{14, 15},
		                        {19, 20}}};
		CHECK((v1 == z1));
	}

	SECTION("join_rows") {
		auto a = SiLi::Matrix{{{ 1,  2,  3},
		                       {11, 12, 13}}};
		auto b = SiLi::Matrix{{{ 4,  5,  6},
		                       {14, 15, 16}}};

		auto c = join_rows(a, b);

		auto z = SiLi::Matrix{{{ 1,  2,  3,  4,  5,  6},
		                       {11, 12, 13, 14, 15, 16}}};
		CHECK((c == z));
	}

	SECTION("join_cols") {
		auto a = SiLi::Matrix{{{ 1,  2,  3},
		                       {11, 12, 13}}};
		auto b = SiLi::Matrix{{{ 4,  5,  6},
		                       {14, 15, 16}}};

		auto c = join_cols(a, b);

		auto z = SiLi::Matrix{{{ 1,  2,  3},
		                       {11, 12, 13},
		                       { 4,  5,  6},
		                       {14, 15, 16}}};
		CHECK((c == z));
	}

	SECTION("det") {
		auto a = SiLi::Matrix{{{ 1,  2},
		                       {11, 12}}};
		auto v = det(a);
		CHECK(v == -10);
	}

	SECTION("cross") {
		auto a = SiLi::Matrix{{{ 1},
		                       { 2},
		                       { 3}}};

		auto b = SiLi::Matrix{{{ 4},
		                       { 5},
		                       { 6}}};

		auto v = cross(a, b);

		auto z = SiLi::Matrix{{{-3},
		                       { 6},
		                       {-3}}};
		CHECK((z == v));
	}

	SECTION("sum") {
		auto a = SiLi::Matrix{{{ 3,  4,  7},
		                       { 5,  6,  8}}};
		auto v = sum(a);
		CHECK(v == 33);
	}

	SECTION("sum_rows") {
		auto a = SiLi::Matrix{{{ 3,  4,  7},
		                       { 5,  6,  8}}};
		auto v = sum_rows(a);
		auto z = SiLi::Matrix{{{14},
		                       {19}}};
		CHECK((v == z));
	}

	SECTION("sum_cols") {
		auto a = SiLi::Matrix{{{ 3,  4,  7},
		                       { 5,  6,  8}}};
		auto v = sum_cols(a);
		auto z = SiLi::Matrix{{{8, 10, 15}}};
		CHECK((v == z));
	}

	SECTION("inv") {
		auto a = SiLi::Matrix{{{ 1.,  2.},
		                       {11., 12.}}};

		auto [v, m] = inv(a);
		auto z = SiLi::Matrix{{{-1.2,  0.2},
		                       { 1.1, -0.1}}};
		CHECK(v == -10);
		CHECK(sum(z - m) < 1e-9);
	}

	SECTION("view_trans") {
		auto a = SiLi::Matrix{{{ 1,  2,  3},
		                       {11, 12, 13}}};
		auto v = view_trans(a);
		auto z0 = SiLi::Matrix{{{ 1, 11},
		                        { 2, 12},
		                        { 3, 13}}};

		CHECK((v == z0));

		v(1, 0) = 99;
		auto z1 = SiLi::Matrix{{{ 1, 99,  3},
		                        {11, 12, 13}}};
		CHECK((a == z1));
	}

	SECTION("trans") {
		auto a = SiLi::Matrix{{{ 1,  2,  3},
		                       {11, 12, 13}}};
		auto v = trans(a);
		auto z = SiLi::Matrix{{{ 1, 11},
		                       { 2, 12},
		                       { 3, 13}}};

		CHECK((v == z));
	}

	SECTION("makeI (0)") {
		auto a = SiLi::makeI<3, 4, int>();
		auto z = SiLi::Matrix{{{1, 0, 0, 0},
		                       {0, 1, 0, 0},
		                       {0, 0, 1, 0}}};
		CHECK((a == z));
	}

	SECTION("makeI (1)") {
		auto a = SiLi::makeI<3, int>();
		auto z = SiLi::Matrix{{{1, 0, 0},
		                       {0, 1, 0},
		                       {0, 0, 1}}};
		CHECK((a == z));
	}

	SECTION("nom") {
		auto a = SiLi::Matrix{{{ 3.},
		                       { 4.},
		                       { 5.}}};
		auto v = norm(a);
		auto z = std::sqrt(50);
		CHECK(v == z);
	}

	SECTION("dot") {
		auto a = SiLi::Matrix{{{ 1},
		                       { 2},
		                       { 3}}};

		auto b = SiLi::Matrix{{{ 4},
		                       { 5},
		                       { 6}}};

		auto v = dot(a, b);

		CHECK(v == 32);
	}

	SECTION("outerProd") {
		auto a = SiLi::Matrix{{{ 1},
		                       { 2},
		                       { 3}}};

		auto b = SiLi::Matrix{{{ 4},
		                       { 5},
		                       { 6}}};

		auto v = outerProd(a, b);

		auto z = SiLi::Matrix{{{ 4,  5,  6},
		                       { 8, 10, 12},
		                       {12, 15, 18}}};
		CHECK((z == v));
	}
}
