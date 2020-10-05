#pragma once

#include "concepts.h"

#include <array>

namespace SiLi2 {

template<int _rows, int _cols, typename T> requires (_rows >= 0 and _cols >= 0)
class Matrix<_rows, _cols, T> {
	std::array<T, _cols*_rows> vals;

public:
	using value_t = T;

	static constexpr int Rows   = _rows;
	static constexpr int Cols   = _cols;
	static constexpr int Stride = _cols;

	constexpr Matrix() : vals{} {}

	constexpr Matrix(T const (&values)[Rows][Cols]) {
		for (int row{0}; row < Rows; ++row) {
			for (int col{0}; col < Cols; ++col) {
				this->operator()(row, col) = values[row][col];
			}
		}
	}

	template <Viewable L>
	constexpr Matrix(L const& view) requires (L::Rows == Rows and L::Cols == Cols) {
		for (int row{0}; row < Rows; ++row) {
			for (int col{0}; col < Cols; ++col) {
				this->operator()(row, col) = view(row, col);
			}
		}
	}


	constexpr auto operator()(int row, int col) -> T& {
		return vals[col + row * Cols];
	}
	constexpr auto operator()(int row, int col) const -> T const& {
		return vals[col + row * Cols];
	}

	constexpr auto operator()(int row) -> T& requires (Cols == 1 and Rows != 1) {
		return vals[row * Stride];
	}
	constexpr auto operator()(int row) const -> T const& requires (Cols == 1 and Rows != 1) {
		return vals[row * Stride];
	}
	constexpr auto operator()(int col) -> T& requires (Rows == 1) {
		return vals[col];
	}
	constexpr auto operator()(int col) const -> T const& requires (Rows == 1) {
		return vals[col];
	}



	constexpr auto view() {
		return MatrixView<_rows, _cols, _cols, T>{vals.data()};
	}
	constexpr auto view() const {
		return MatrixView<_rows, _cols, _cols, T const>{vals.data()};
	}
	constexpr auto data() -> T* {
		return vals.data();
	}
	constexpr auto data() const -> T const* {
		return vals.data();
	}

	constexpr auto operator=(T const& s) -> Matrix& {
		for (int row{0}; row < Rows; ++row) {
			for (int col{0}; col < Cols; ++col) {
				this->operator()(row, col) = s;
			}
		}
		return *this;
	}

	template <Viewable V>
	constexpr auto operator=(V const& v) -> Matrix& requires (Rows == V::Rows and Cols == V::Cols) {
		for (int row{0}; row < Rows; ++row) {
			for (int col{0}; col < Cols; ++col) {
				this->operator()(row, col) = v(row, col);
			}
		}
		return *this;
	}

};

template <int rows, int cols, typename T>
Matrix(T const (&)[rows][cols]) -> Matrix<rows, cols, T>;

template <int rows, int cols, int stride, typename T>
Matrix(MatrixView<rows, cols, stride, T> const&) -> Matrix<rows, cols, T>;


}
