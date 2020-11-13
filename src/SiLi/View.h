#pragma once

#include "concepts.h"

namespace SiLi {

/*! Represents a view onto a matrix
 *
 * Fullfills the _concept::Matrix concept.
 *
 * If _rows or _cols is 1 it fullfills the _concept::Vector concept.
 *
 * \caption Template Parameters
 * \param _rows number of rows of the view, must be larger or equal to zero
 * \param _cols number of columns of the view, must be larger or equal to zero
 * \param _stride length of the stride
 * \param T     type of the elements
 * \param _transposed indicates if index access of _rows and _cols should be inverted to simulate transpose access
 *
 * * \caption Methods
 * \param data() returns pointer to the underlying data structure
 * \param m(row,col) access element at ``row`` and ``col``

 */

template<int _rows, int _cols, int _stride, typename T, bool _transposed> requires (_rows >= 0 and _cols >= 0 and _stride >= 0)
class View<_rows, _cols, _stride, T, _transposed> {
	T& mData;
public:
	using value_t = T;
	static constexpr int  Rows       = _rows;
	static constexpr int  Cols       = _cols;
	static constexpr int  Stride     = _stride;
	static constexpr bool Transposed = _transposed;

	constexpr View(T* _data)
		: mData{*_data}
	{};

	constexpr View(Matrix<_rows, _cols, T> matrix)
		: mData{*matrix.data()}
	{};

	View(View& v)
		: mData{v.mData}
	{}

	View(View const& v) = delete;

	constexpr auto data() -> T* {
		return &mData;
	}

	constexpr auto data() const -> T const* {
		return &mData;
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



	constexpr auto operator=(T const& s) -> View& {
		for_each_constexpr<View>([&]<int row, int col>() {
			this->operator()(row, col) = s;
		});
		return *this;
	}
	template <_concept::Matrix V>
	constexpr auto operator=(V const& v) -> View& requires (Rows == V::Rows and Cols == V::Cols) {
		for_each_constexpr<View>([&]<int row, int col>() {
			this->operator()(row, col) = v(row, col);
		});
		return *this;
	}


};

template <int _rows, int _cols, typename T>
View(Matrix<_rows, _cols, T>&) -> View<_rows, _cols, _cols, T, false>;

template <int _rows, int _cols, typename T>
View(Matrix<_rows, _cols, T> const&) -> View<_rows, _cols, _cols, T const, false>;

}
