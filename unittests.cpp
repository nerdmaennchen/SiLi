#include "SiLi.h"
#include "SiLi-iostream.h"


template <int rows, int cols, typename P1, typename P2, typename T>
bool operator==(SiLi::MatrixView<rows, cols, P1, T> const& lhs,
                SiLi::MatrixView<rows, cols, P2, T> const& rhs) {

	for (int i(0); i < rows; ++i) {
		for (int j(0); j < cols; ++j) {
			if (lhs(i, j) != rhs(i, j)) return false;
		}
	}
	return true;
}
template <int rows, int cols, typename Prop, typename T>
auto asMat(SiLi::MatrixView<rows, cols, Prop, T> const& view) -> SiLi::Matrix<rows, cols> {
	return {view};
}
template <typename T1>
bool equal(T1 const& t1) {
	return t1 == asMat(t1);
}
template <typename T1, typename T2>
bool equal(T1 const& t1, T2 const& t2) {
	return t1 == t2;
}



using namespace SiLi;
template <int rows, int cols>
using M = Matrix<rows, cols>;

#define CHECK(a, b) \
	if ((a) != (b)) { \
		std::cout << "unittest fail in line: " << __LINE__ << "\n"; \
		std::cout << std::boolalpha << false << ": " << #a << " -> " << a << " != " << #b << "\n"; \
	}

/** Tests
*/
void unittests() {
	std::cout << "Running unittests: \n";
	// tests of 0x0 matrices
	{
		std::cout << "- Matrices of size 0\n";
		M< 0, 0>  m1({});
		M<10, 0>  m2({});
		M< 0, 10> m3({});

		{ auto m = m3 * m2; }
		{ auto m = m1 + m1; }
	}
	// tests of 1x1 matrices
	{
		std::cout << "- Matrices of size 1x1\n";
		M<1, 1> m1({1});
		M<1, 1> m2({2});
		CHECK(m1(0, 0), 1);
		CHECK(m2(0, 0), 2);
		CHECK((m1 + m2)(0, 0), 3);
	}

	// tests
	{
		SiLi::Matrix<4, 4> c({{1, 2, 3, 4}, {5, 6, 7, 8}, {9, 10, 11, 12}, {13, 14, 15, 16}});
		std::cout << c << std::endl;
		std::cout << c.t() << std::endl;
		std::cout << c.t().subview<3, 2>(0, 1) << std::endl;
		std::cout << c.t().subview<2, 3>(2, 1) << std::endl;
		std::cout << c.t().subview<2, 3>(2, 1).t() << std::endl;
		std::cout << c.t().subview<3, 2>(0, 1) + c.t().subview<2, 3>(2, 1).t() <<  std::endl;
	}
	// tests read only properties
	{
		SiLi::Matrix<4, 4> c({{1, 2, 3, 4}, {5, 6, 7, 8}, {9, 10, 11, 12}, {13, 14, 15, 16}});
		auto v1 { c.subview<3, 2>(0, 0) };
		auto v2 { v1 };
		auto const v3 { c.subview<3, 2>(0, 0) };
//		auto v4 { v3 };

		debugPrint(v2);
		debugPrint(v3);
		debugPrint(v3.t());
//		debugPrint(v4);
		v2 = v1;
//		v4 = v1;
	}

}

