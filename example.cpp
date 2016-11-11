
#include "SiLi.h"
#include <iostream>

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
}



