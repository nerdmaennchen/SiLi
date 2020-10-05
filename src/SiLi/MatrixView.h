#pragma once

#include "concepts.h"

namespace SiLi2 {

template<int, int, typename, typename...>
class Matrix;

template<int, int, int, typename, typename...>
class MatrixView;

template<int _rows, int _cols, int _stride, typename T> requires (_rows >= 0 and _cols >= 0 and _stride >= 0)
class MatrixView<_rows, _cols, _stride, T> {
	T* mData;
public:
	using value_t = T;
	static constexpr int Rows   = _rows;
	static constexpr int Cols   = _cols;
	static constexpr int Stride = _stride;

	constexpr int num_rows() const { return _rows; }
	constexpr int num_cols() const { return _cols; }
	constexpr int stride() const { return _stride; }

	constexpr MatrixView(T* _data)
		: mData{_data}
	{};

	constexpr MatrixView(Matrix<_rows, _cols, T> matrix)
		: mData{matrix.data()}
	{};

	constexpr auto data() -> T* {
		return mData;
	}

	constexpr auto data() const -> T const* {
		return mData;
	}

	constexpr auto operator()(int row, int col) -> T& {
		return mData[col + row * Stride];
	}
	constexpr auto operator()(int row, int col) const -> T const& {
		return mData[col + row * Stride];
	}

	constexpr auto operator()(int row) -> T& requires (Cols == 1 and Rows != 1) {
		return mData[row * Stride];
	}
	constexpr auto operator()(int row) const -> T const& requires (Cols == 1 and Rows != 1) {
		return mData[row * Stride];
	}
	constexpr auto operator()(int col) -> T& requires (Rows == 1) {
		return mData[col];
	}
	constexpr auto operator()(int col) const -> T const& requires (Rows == 1) {
		return mData[col];
	}



	constexpr auto operator=(T const& s) -> MatrixView& {
		for (int row{0}; row < Rows; ++row) {
			for (int col{0}; col < Cols; ++col) {
				this->operator()(row, col) = s;
			}
		}
		return *this;
	}
	template <Viewable V>
	constexpr auto operator=(V const& v) -> MatrixView& requires (Rows == V::Rows and Cols == V::Cols) {
		for (int row{0}; row < Rows; ++row) {
			for (int col{0}; col < Cols; ++col) {
				this->operator()(row, col) = v(row, col);
			}
		}
		return *this;
	}


};

template <int _rows, int _cols, typename T>
MatrixView(Matrix<_rows, _cols, T>&) -> MatrixView<_rows, _cols, _cols, T>;

template <int _rows, int _cols, typename T>
MatrixView(Matrix<_rows, _cols, T> const&) -> MatrixView<_rows, _cols, _cols, T const>;


}
