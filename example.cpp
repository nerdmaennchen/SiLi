
#include "SiLi.h"
#include "SiLi-iostream.h"

void unittests();

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






int main() {
	SiLi::Matrix<3, 2> b({{1, 0}, {0, 1}, {1, 0}});
	SiLi::Matrix<2, 3> c({{1, 0, 0}, {0, 1, 0}});
	std::cout << b << std::endl;
	std::cout << c << std::endl;
	std::cout << b * c << std::endl;
	std::cout << (b * c).t() << std::endl;

	std::cout << SiLi::Matrix<1, 1>(.1f).inv() << std::endl;

	auto subview = b.subview<3, 1>(0, 0);
	subview(0,0) = -1;
	subview(1,0) = -2;
	subview(2,0) = -3;
	std::cout << b << std::endl;
	std::cout << subview << std::endl;
	subview *= 2;
	std::cout << subview << std::endl;
	std::cout << b << std::endl;

	{
		SiLi::Matrix<2, 2> b({{1, 2}, {3, 4}});
		SiLi::Matrix<2, 2> c({{1, 2}, {3, 4}});
		auto view = (b+c).submat<2, 2>(0, 0);
		view(2, 2) = 5;

		SiLi::Matrix<2, 2> d({{1, 2}, {3, 4}});
		std::cout << view << std::endl;
	}
	{
		SiLi::Matrix<2, 2> b({{1, 2}, {3, 4}});
		SiLi::Matrix<2, 2> c({{1, 2}, {3, 4}});
		b = c;
		b.subview<2, 2>(0, 0) = c;
		c = b.subview<2, 2>(0, 0);
		c.subview<2, 2>(0, 0) = b.subview<2, 2>(0, 0);
	}
	{
		SiLi::Matrix<2, 2> b({{1, 2}, {3, 4}});
		SiLi::Matrix<2, 2> c({{1, 2}, {3, 4}});
		b += c;
		b -= c;
		b + c;
		b - c;
		b * c;

		SiLi::Matrix<2, 2> d(b.subview<2, 2>(0, 0));

		(b+c).submat<1, 1>(1, 1);
		b.submat<1, 1>(1, 1);

		SiLi::Matrix<4, 4> e(0);

		b.subview<1, 1>(1, 1) = e.subview<1, 1>(0, 0);

		// lines that shouldn't compile
//		b.submat<1, 1>(1, 1) = e.subview<1, 1>(0, 0);
//		SiLi::Matrix<2, 2>(0).subview<1, 1>(0, 0);
	}
	{
		SiLi::Matrix<2, 2> b({{1, 2}, {3, 4}});
		std::cout << b << std::endl;
		std::cout << b.t() << std::endl;

		SiLi::Matrix<4, 4> c({{1, 2, 3, 4}, {5, 6, 7, 8}, {9, 10, 11, 12}, {13, 14, 15, 16}});

		std::cout << std::boolalpha;
		std::cout << equal(c) << std::endl;
		std::cout << equal(c.t()) << std::endl;
		std::cout << equal(c.subview<2, 2>(1, 1)) << std::endl;
		std::cout << equal(c.subview<2, 2>(1, 1).t()) << std::endl;
		std::cout << equal(c.subview<2, 2>(2, 1)) << std::endl;
		std::cout << equal(c.subview<2, 2>(2, 1).t()) << std::endl;
		std::cout << equal(c.subview<3, 3>(0, 1)) << std::endl;
		std::cout << equal(c.subview<3, 3>(0, 1).t()) << std::endl;
		std::cout << equal(c.subview<3, 3>(1, 0)) << std::endl;
		std::cout << equal(c.subview<3, 3>(1, 0).t()) << std::endl;
	}
	std::cout << "next section\n";
	{
		SiLi::Matrix<3, 2> b({{1, 0}, {0, 1}, {1, 0}});
		SiLi::Matrix<2, 3> c({{1, 0, 0}, {0, 1, 0}});
		SiLi::Matrix<3, 3> d({{1, 0, 0}, {0, 1, 0}, {1, 0, 0}});
		SiLi::Matrix<3, 3> e({{1, 0, 1}, {0, 1, 0}, {0, 0, 0}});

		std::cout << equal(b * c) << std::endl;
		std::cout << equal((b * c).t()) << std::endl;
		std::cout << equal(b * c, d) << std::endl;
		std::cout << equal((b * c).t(), e) << std::endl;

		//auto const& f = b;
		//SiLi::MatrixView<2, 3, 2, true, const float> g = f.t();
		//g = f.t();

	}
	unittests();
}
