#pragma once

#include <iostream>

#include "concepts.h"

namespace SiLi {

template <Viewable V>
auto operator<<(std::ostream& os, V const& v) -> auto& {
	for (int row {0}; row < rows_v<V>; ++row) {
		for (int col {0}; col < cols_v<V>; ++col) {
			os << v(row, col) << " ";
		}
		os << "\n";
	}
	return os;
}

}
