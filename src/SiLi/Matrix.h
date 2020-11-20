#pragma once

#include "concepts.h"

#include <array>

namespace SiLi {

/*! Represents a matrix
 *
 * Fullfills the _concept::Matrix concept.
 *
  * If _rows or _cols is 1 it fullfills _concept::Vector.
 *
 * \caption Template Parameters
 * \param _rows number of rows of the matrix, must be larger or equal to zero
 * \param _cols number of columns of the matrix, must be larger or equal to zero
 * \param T     type of the elements
 *
 * \caption Methods
 * \param data() returns pointer to the underlying data structure
 * \param m(row,col) access element at ``row`` and ``col``
 */
template<int _rows, int _cols, typename T> requires (_rows >= 0 and _cols >= 0)
class Matrix<_rows, _cols, T> {
	std::array<T, _cols*_rows> vals;

public:
	using value_t = T;

	static constexpr int  Rows       = _rows;
	static constexpr int  Cols       = _cols;
	static constexpr int  Stride     = _cols;
	static constexpr bool Transposed = false;

	constexpr Matrix() : vals{} {}

	template <typename ...S>
	constexpr Matrix(S... _values)
		: vals{std::forward<S>(_values)...}
	{}

	constexpr Matrix(T const (&values)[Rows][Cols]) {
		for (int row{0}; row < Rows; ++row) {
			for (int col{0}; col < Cols; ++col) {
				this->operator()(row, col) = values[row][col];
			}
		}
	}
	template <_concept::Matrix V>
	constexpr Matrix(V const& view) requires (V::Rows == Rows and V::Cols == Cols) {
		*this = view;
	}

	constexpr auto data() -> T* {
		return vals.data();
	}
	constexpr auto data() const -> T const* {
		return vals.data();
	}


	explicit constexpr operator T() requires (Rows == 1 and Cols == 1) {
		return get(*this, 0, 0);
	}
	explicit constexpr operator T const() const requires (Rows == 1 and Cols == 1) {
		return get(*this, 0, 0);
	}


	constexpr auto operator()(int row, int col) -> T& {
		return get(*this, row, col);
	}
	constexpr auto operator()(int row, int col) const -> T const& {
		return get(*this, row, col);
	}

	constexpr auto operator()(int entry) -> T& requires (Cols == 1 or Rows == 1) {
		return get(*this, entry);
	}
	constexpr auto operator()(int entry) const -> T const& requires (Cols == 1 or Rows == 1) {
		return get(*this, entry);
	}
	constexpr auto operator[](int entry) -> T& requires (Cols == 1 or Rows == 1) {
		return get(*this, entry);
	}
	constexpr auto operator[](int entry) const -> T const& requires (Cols == 1 or Rows == 1) {
		return get(*this, entry);
	}

	template <int row, int col>
	constexpr auto at() -> T& {
		return get<row, col>(*this);
	}
	template <int row, int col>
	constexpr auto at() const -> T const& {
		return get<row, col>(*this);
	}

	template <int entry>
	constexpr auto at() -> T& requires (Cols == 1 or Rows == 1) {
		return get<entry>(*this);
	}
	template <int entry>
	constexpr auto at() const -> T const& requires (Cols == 1 or Rows == 1) {
		return get<entry>(*this);
	}





	constexpr auto view() {
		return View<_rows, _cols, _cols, T, false>{data()};
	}
	constexpr auto view() const {
		return View<_rows, _cols, _cols, T const, false>{data()};
	}

	constexpr auto operator=(T const& s) -> Matrix& {
		for_each_constexpr<Matrix>([&]<int row, int col>() constexpr {
			at<row, col>() = s;
		});
		return *this;
	}

	template <_concept::Matrix V>
	constexpr auto operator=(V const& v) -> Matrix& requires (Rows == V::Rows and Cols == V::Cols) {
		for_each_constexpr<Matrix>([&]<int row, int col>() constexpr {
			at<row, col>() = v.template at<row, col>();
		});
		return *this;
	}
};

template <int rows, int cols, typename T>
Matrix(T const (&)[rows][cols]) -> Matrix<rows, cols, T>;

template <int rows, int cols, int stride, typename T, bool _transposed>
Matrix(View<rows, cols, stride, T, _transposed> const&) -> Matrix<rows, cols, std::remove_const_t<T>>;

template<int rows, typename T>
using Vector = Matrix<rows, 1, T>;

}
