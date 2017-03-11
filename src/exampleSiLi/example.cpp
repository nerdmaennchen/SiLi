
#include <SiLi/SiLi.h>
#include <SiLi/SiLi-iostream.h>

int main() {
	auto b = SiLi::make_mat<double, 3, 2>({{1., 0.}, {0., 1.}, {1., 0.}}); // SiLi::Matrix<3, 2>
	auto c = SiLi::make_mat<double, 2, 3>({{1., 0., 0.}, {0., 1., 0.}});   // SiLi::Matrix<2, 3>
	std::cout << b << std::endl;
	std::cout << c << std::endl;
	std::cout << b * c << std::endl;
	std::cout << (b * c).t() << std::endl;

	std::cout << SiLi::Matrix<1, 1>(.1f).inv() << std::endl;

	auto view = b.view<3, 1>(0, 0);
	std::cout << view << std::endl;
	std::cout << b.diag() << std::endl;


}
