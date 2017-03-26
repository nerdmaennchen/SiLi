
#include <SiLi/SiLi.h>
#include <SiLi/SiLi-simple.h>
#include <iostream>

using namespace SiLi;
int main() {
	{
		auto v = Matrix<3, 1, double>({{2., 0., 0.}});
		std::cout << "norm: " << v.norm() << std::endl;
		std::cout << "elementwise: " << v.applyElementWise([](double const& v) { return v/2.; }).t() << std::endl;
		std::cout << "norm: " << (v * 0.5).t() << std::endl;
		std::cout << "norm: " << (v * (1. / v.norm())).t() << std::endl;
//		std::cout << "norm: " << v.normalized() << std::endl;
		return 0;
	}
	{
		auto m = Matrix<-1, 2, double>();
		m = Matrix<3, 2, double>({{1., 0.}, {0., 1.}, {3., 4.}});
		std::cout << m.num_rows() << " " << m.num_cols() << std::endl;
		std::cout << m << std::endl;
	}

	{
		auto m = Matrix<-1, -1, double>();
		m = Matrix<3, 2, double>({{1., 0.}, {0., 1.}, {3., 4.}});
		std::cout << m.num_rows() << " " << m.num_cols() << std::endl;
		std::cout << m << std::endl;
	}

	{
		auto m1 = Matrix<-1, 2, double>();
		m1 = Matrix<3, 2, double>({{1., 0.}, {0., 1.}, {3., 4.}});

		auto m2 = Matrix<2, -1, double>();
		m2 = Matrix<3, 2, double>({{1., 0.}, {0., 1.}, {3., 4.}}).t();

		auto m3 = m1 * m2;
		auto m4 = m2 * m1;

		std::cout << "m1: " << std::endl;
		std::cout << m1.num_rows() << " " << m1.num_cols() << std::endl;
		std::cout << m1 << std::endl;

		std::cout << "m2: " << std::endl;
		std::cout << m2.num_rows() << " " << m2.num_cols() << std::endl;
		std::cout << m2 << std::endl;


		std::cout << "m3: " << std::endl;
		std::cout << m3.num_rows() << " " << m3.num_cols() << std::endl;
		std::cout << m3 << std::endl;

		std::cout << "m4: " << std::endl;
		std::cout << m4.num_rows() << " " << m4.num_cols() << std::endl;
		std::cout << m4 << std::endl;
	}

	{
		auto m = Matrix<1, -1, double>();
		m = Matrix<1, 4, double>({{1., 0., 0., 1.}});

		for (int x(0); x<5; ++x) {
			std::cout << "m(" << x << ")" << std::endl;
			std::cout << m.num_rows() << " " << m.num_cols() << std::endl;
			std::cout << m << std::endl;
			m = m.join_cols(m*2.);
		}
	}



/*	auto b = SiLi::make_mat<3, 2, double>({{1., 0.}, {0., 1.}, {1., 0.}}); // SiLi::Matrix<3, 2>
	auto c = SiLi::make_mat<2, 3, double>({{1., 0., 0.}, {0., 1., 0.}});   // SiLi::Matrix<2, 3>
	std::cout << b << std::endl;
	std::cout << c << std::endl;
	std::cout << b * c << std::endl;
	std::cout << (b * c).t() << std::endl;

	for (auto row : b.rows<1,2>()) {
		std::cout << row;
	}
	std::cout << std::endl;

	std::cout << SiLi::Matrix<1, 1>(.1f).inv() << std::endl;

	auto view = b.view<3, 1>(0, 0);
	std::cout << view << std::endl;
	std::cout << b.diag() << std::endl;*/


}
