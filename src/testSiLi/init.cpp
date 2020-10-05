#if 0
#include <SiLi/SiLi-Mat.h>
#include <catch2/catch.hpp>

template <size_t Ns>
void f(std::array<double, Ns> list) {
//	constexpr auto Rows = sizeof...(list);
}

template <typename T, int Ns>
using carray = T const (&)[Ns];

/*template<typename T, int N0, int ...Ns>
auto g(carray<T, N0>, carray<T, Ns>... list) noexcept {
	return SiLi_Mat::Matrix<sizeof...(list)+1, N0, T>{};
}

template <size_t ...Ns>
auto h(std::array<int, Ns>... list) {
	std::array a = {list...};
}*/



TEST_CASE("initialize matrix","[init]") {
	SECTION("matrices with one element") {
//		auto m = SiLi_Mat::Matrix<1, 3, int> {{0, 0, 0}};
		//auto m = SiLi_Mat::Matrix<2, 3> (
		auto m = SiLi_Mat::make_mat<2, 3> (
			0, 0, 0,
			0, 0, 0
		);
		auto m2 = SiLi_Mat::make_mat<3, 2> (
			0, 0, 0,
			0, 0, 0
		);
		auto m3 = m * m2;


		CHECK(m.rows() == 2);
		CHECK(m.cols() == 3);
	}
	std::array x = {0., 0.};
	f(x);
/*	g({0., 0.}, {0., 0.});
	std::array a1 = {0, 0};
	std::array a2 = {0, 0};
	h(a1, std::array{0, 0});*/
/*	auto x = h(
		Row{0, 0, 0},
		Row{0, 0, 0},
	);*/
}

#endif
