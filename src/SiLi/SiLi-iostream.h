#pragma once

#include <iostream>
#include <typeinfo>

namespace SiLi
{
template<int rows, int cols, typename Prop, typename T>
std::ostream& operator<< (std::ostream& stream, SiLi::MatrixView<rows, cols, Prop, T const> const& view) {
	for (int i(0); i < view.num_rows(); ++i) {
		for (int j(0); j < view.num_cols(); ++j) {
			stream << view(i, j) << "\t";
		}
		stream << "\n";
	}
	return stream;
}
template<int rows, int cols, typename Prop, typename T>
void _debugPrint(std::string const& _name, MatrixView<rows, cols, Prop, T const> const& view) {
	std::cout << _name << "\n";
	std::cout << view;
	std::cout << "type: " << typeid(T).name()
	          << ", transposed: " << std::boolalpha << Prop::transposed
//	          << ", read only: " << Prop::readOnly"
	          << "\n\n";
}
#define debugPrint(a) _debugPrint(#a, a)

}
