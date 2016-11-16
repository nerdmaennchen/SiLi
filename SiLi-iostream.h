#pragma once

#include <iostream>

namespace SiLi
{
template<int rows, int cols, int rowStride, typename T>
std::ostream& operator<< (std::ostream& stream, SiLi::MatrixView<rows, cols, rowStride, T> const& view) {
	for (int i(0); i < view.num_rows(); ++i) {
		for (int j(0); j < view.num_cols(); ++j) {
			stream << view(i, j) << "\t";
		}
		stream << "\n";
	}
	return stream;
}

}
