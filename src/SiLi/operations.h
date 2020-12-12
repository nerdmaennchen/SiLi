#pragma once

#include "concepts.h"

#include <algorithm>
#include <cmath>
#include <tuple>

namespace SiLi {

struct End {};
template <typename T>
concept CEnd = std::is_same_v<End, T>;

template <typename T>
constexpr auto value() -> typename std::add_rvalue_reference<value_t<T>>::type;

/*! Access element
 * \shortexample at<Row, Col>(m)
 * \group Free Matrix Functions
 *
 * Checking accessibility at compile time.
 *
 * \caption Template Parameters
 * \param Row access row. Valid value [0-Row)
 * \param Col access column. Validvalue [0-Col)
 *
 * \caption Parameters
 * \param m   _concept::Matrix
 * \return    returns reference to underlying element
 *
 * \code
 *   auto m = SiLi::Matrix{{{3, 4},
 *                          {5, 6}}};
 *   at<0, 0>(m) = 5;
 *   at<0, 1>(m) = 6;
 *   std::cout << m << "\n"; // prints {{5, 6},
 *                                      {5, 6}}
 * \endcode
 */
template <int Row, int Col, _concept::Matrix M>
constexpr auto at(M&& m) -> auto& {
	static_assert(0 <= Row and Row < rows_v<M>, "accessing outside of valid range");
	static_assert(0 <= Col and Col < cols_v<M>, "accessing outside of valid range");
	return m.template at<Row, Col>();
}

/*! Access element
 * \shortexample at<Idx>(v);
 * \group Free Vector Functions
 *
 * Checking accessibility at compile time.
 *
 * \param Idx row access
 * \param m   _concept::Vector
 * \return    returns reference to underlying element
 *
 * \code
 *   auto m = SiLi::Matrix{{{3},
 *                          {5}}};
 *   at<0>(m) = 5;
 *   at<1>(m) = 6;
 *   std::cout << m << "\n"; // prints {{5},
 *                                      {6}}
 * \endcode
 */
template <int Idx, _concept::Vector V>
constexpr auto at(V&& v) -> auto& {
	if constexpr (cols_v<V> == 1) {
		return at<Idx, 0>(v);
	} else {
		return at<0, Idx>(v);
	}
}

/*! Number of rows
 * \shortexample rows(m)
 * \group Free Matrix Functions
 *
 * \param m _concept::Matrix
 * \return  number of rows
 *
 * \code
 *   auto m = SiLi::Matrix{{{3, 4, 7},
 *                          {5, 6, 8}}};
 *   std::cout << rows(m) << "\n"; // prints 2
 * \endcode
 */
template <_concept::Matrix M>
constexpr auto rows(M const& m) {
	return rows_v<M>;
}


/*! Number of columns
 * \shortexample cols(m)
 * \group Free Matrix Functions
 *
 * \param  m _concept::Matrix
 * \return   number of columns
 *
 * \code
 *   auto m = SiLi::Matrix{{{3, 4, 7},
 *                          {5, 6, 8}}};
 *   std::cout << cols(m) << "\n"; // prints 3
 * \endcode
 */
template <_concept::Matrix M>
constexpr auto cols(M const& m) {
	return cols_v<M>;
}

namespace details {
template <_concept::Matrix V, typename Operator>
constexpr auto apply(V const& v, Operator op) {
	using U = decltype(op(value<V>()));
	auto ret = Matrix<rows_v<V>, cols_v<V>, U>{};
	for_each_constexpr<V>([&]<auto row, auto col>() {
		at<row, col>(ret) = op(at<row, col>(v));
	});
	return ret;
}

template <_concept::Matrix L, _concept::Matrix R, typename Operator> requires (rows_v<L> == rows_v<R> and cols_v<L> == cols_v<R>)
constexpr auto apply(L const& l, R const& r, Operator op) {
	using U = decltype(op(value<L>(), value_t<R>()));
	auto ret = Matrix<rows_v<L>, cols_v<L>, U>{};
	for_each_constexpr<L>([&]<auto row, auto col>() {
		at<row, col>(ret) = op(at<row, col>(l), at<row, col>(r));
	});
	return ret;
}

template <_concept::Matrix V, typename Operator>
constexpr void self_assign_apply(V&& v, Operator op) {
	for_each_constexpr<V>([&]<auto row, auto col>() {
		op(at<row, col>(v));
	});
}

template <_concept::Matrix L, _concept::Matrix R, typename Operator> requires (rows_v<L> == rows_v<R> and cols_v<L> == cols_v<R>)
constexpr void self_assign_apply(L&& l, R const& r, Operator op) {
	for_each_constexpr<L>([&]<auto row, auto col>() {
		op(at<row, col>(l), at<row, col>(r));
	});
}
}

/*! Copy of m
 * \shortexample +m
 * \group Matrix Operations
 *
 * \param m _concept::Matrix
 * \return copy of m with ``+`` applied to every element
 *
 * \code
 *   auto m = SiLi::Matrix{{{3, 4, 7},
 *                          {5, 6, 8}}};
 *   auto m2 = +m;
 *   std::cout << m << "\n";
 *   std::cout << m2 << "\n";
 * \endcode
 */
template <_concept::Matrix M>
constexpr auto operator+(M const& m) {
	return details::apply(m, [](auto e) constexpr { return +e; });
}

/*! Elementwise addition
 * \shortexample l + r
 * \group Matrix Operations
 *
 * \param l left _concept::Matrix of addition
 * \param r right _concept::Matrix of addition
 * \return  Matrix with l and r
 *
 * l and r must have the same dimension.
 *
 * \code
 *   auto a = SiLi::Matrix{{{3, 4, 7},
 *                          {5, 6, 8}}};
 *   auto b = SiLi::Matrix{{{1, 2, 3},
 *                          {4, 5, 6}}};
 *   auto c = a + b;
 *   std::cout << c << "\n"; // prints {{4,  6, 10},
 *                                      {9, 11, 14}}
 * \endcode
 */
template <_concept::Matrix L, _concept::Matrix R>
constexpr auto operator+(L const& l, R const& r) {
	return details::apply(l, r, [](auto _l, auto _r) constexpr { return _l + _r; });
}

/*! Elementwise addition
 * \shortexample l += r
 * \group Matrix Operations
 *
 * \param l left _concept::Matrix
 * \param r right _concept::Matrix
 * \return  reference to l
 *
 * l and r must have the same dimension.
 *
 * \code
 *   auto a = SiLi::Matrix{{{3, 4, 7},
 *                          {5, 6, 8}}};
 *   auto b = SiLi::Matrix{{{1, 2, 3},
 *                          {4, 5, 6}}};
 *   auto a += b;
 *   std::cout << a << "\n"; // prints {{4,  6, 10},
 *                                      {9, 11, 14}}
 * \endcode
 */
template <_concept::Matrix L, _concept::Matrix R>
constexpr auto operator+=(L& l, R const& r) -> auto& {
//	l = l + r;
	details::self_assign_apply(l, r, [](auto& _l, auto _r) constexpr { _l += _r; });
	return l;
}

/*! Elementwise negation
 * \shortexample -m
 * \group Matrix Operations
 *
 * \param m _concept::Matrix
 * \return  Matrix with negated elements
 *
 * \code
 *   auto m = SiLi::Matrix{{{3, 4, 7},
 *                          {5, 6, 8}}};
 *   auto m2 = -m;
 *   std::cout << m << "\n";
 *   std::cout << m2 << "\n";
 * \endcode
 */
template <_concept::Matrix M>
constexpr auto operator-(M const& m) {
	return details::apply(m, [](auto e) constexpr { return -e; });
}

/*! Elementwise subtraction
 * \shortexample l - r
 * \group Matrix Operations
 *
 * \param l left _concept::Matrix
 * \param r right _concept::Matrix
 * \return  matrix with (l - r)
 *
 * l and r must have the same dimension.
 *
 * \code
 *   auto a = SiLi::Matrix{{{3, 4, 7},
 *                          {5, 6, 8}}};
 *   auto b = SiLi::Matrix{{{1, 2, 3},
 *                          {4, 5, 6}}};
 *   auto c = a - b;
 *   std::cout << c << "\n"; // prints {{2,  2, 4},
 *                                      {1,  1, 2}}
 * \endcode
 */
template <_concept::Matrix L, _concept::Matrix R>
constexpr auto operator-(L const& l, R const& r) {
	return details::apply(l, r, [](auto _l, auto _r) constexpr { return _l - _r; });
}

/*! Elementwise subtraction
 * \shortexample l += r
 * \group Matrix Operations
 *
 * \param l left _concept::Matrix
 * \param r right _concept::Matrix
 * \return  reference to l
 *
 * l and r must have the same dimension.
 *
 * \code
 *   auto a = SiLi::Matrix{{{3, 4, 7},
 *                          {5, 6, 8}}};
 *   auto b = SiLi::Matrix{{{1, 2, 3},
 *                          {4, 5, 6}}};
 *   auto a += b;
 *   std::cout << a << "\n"; // prints {{2,  2, 4},
 *                                      {1,  1, 2}}
 * \endcode
 */
template <_concept::Matrix L, _concept::Matrix R>
constexpr auto operator-=(L& l, R const& r) -> auto& {
	details::self_assign_apply(l, r, [](auto& _l, auto _r) constexpr { _l -= _r; });
	return l;
}

/*! Matrix multiplication
 * \shortexample l * r
 * \group Matrix Operations
 *
 * \param l _concept::Matrix
 * \param r _concept::Matrix
 * \return Matrix with result of l and r multiplication
 *
 * l must have the same number of columns as r has rows.
 *
 * \code
 *   auto a = SiLi::Matrix{{{3, 4, 7},
 *                          {5, 6, 8}}};
 *   auto b = SiLi::Matrix{{{1, 2},
 *                          {4, 5},
 *                          {3, 6}};
 *   auto c = a * b; // c is of type SiLi::Matrix<2, 2, int>
 *   std::cout << c << "\n"; // prints {{40, 68},
 *                                      {53, 88}}
 * \endcode
 */
template <_concept::Matrix L, _concept::Matrix R> requires (L::Cols == R::Rows)
constexpr auto operator*(L const& l, R const& r) {
	using U = decltype(std::declval<typename L::value_t>() * std::declval<typename R::value_t>());

	if constexpr (L::Rows == 1 and R::Cols == 1) {
		auto ret = U{};
		for_constexpr<0, L::Cols>([&]<int I>() {
			ret += at<I>(l) * at<I>(r);
		});
		return Matrix{{{ret}}};
	} else {
		auto ret = Matrix<L::Rows, R::Cols, U>{};
		for_constexpr<0, L::Rows>([&]<int row>() {
			for_constexpr<0, R::Cols>([&]<int col>() {
				at<row, col>(ret) = U{view_row<row>(l) * view_col<col>(r)};
			});
		});
		return ret;
	}
}

/*! Scalar multiplication
 * \shortexample l * s
 * \group Matrix Operations
 *
 * \param l _concept::Matrix
 * \param s scalar value
 * \return Matrix with each element multiplied with s
 *
 * \code
 *   auto a = SiLi::Matrix{{{3, 4, 7},
 *                          {5, 6, 8}}};
 *   auto c = a * 5;
 *   std::cout << c << "\n"; // prints {{15, 20, 35},
 *                                      {25, 30, 40}}
 * \endcode
 */
template <_concept::Matrix L>
constexpr auto operator*(L const& l, value_t<L> const& s) {
	return details::apply(l, [s](auto _l) { return _l * s; });
}

/*! Scalar multiplication
 * \shortexample s * r
 * \group Matrix Operations
 *
 * \param s scalar value
 * \param r _concept::Matrix
 * \return Matrix with each element multiplied with s
 *
 * \code
 *   auto a = SiLi::Matrix{{{3, 4, 7},
 *                          {5, 6, 8}}};
 *   auto c = 5 * a;
 *   std::cout << c << "\n"; // prints {{15, 20, 35},
 *                                      {25, 30, 40}}
 * \endcode
 */
template <_concept::Matrix R>
constexpr auto operator*(value_t<R> const& s, R&& r) {
	return std::forward<R>(r) * s;
}

/*! Scalar multiplication
 * \shortexample l *= s
 * \group Matrix Operations
 *
 * \param l _concept::Matrix
 * \param s scalar value
 * \return Referenc to l
 *
 * \code
 *   auto a = SiLi::Matrix{{{3, 4, 7},
 *                          {5, 6, 8}}};
 *   auto a *= 5;
 *   std::cout << a << "\n"; // prints {{15, 20, 35},
 *                                      {25, 30, 40}}
 * \endcode
 */
template <_concept::Matrix L>
constexpr auto operator*=(L&& l, value_t<L> const& s) -> auto& {
	details::self_assign_apply(l, [s](auto& _l) { _l *= s; });
	return l;
}

/*! Divide elements by a scalar
 * \shortexample l / s
 * \group Matrix Operations
 *
 * \param l _concept::Matrix
 * \param s scalar
 * \return Matrix same size and type as l
 *
 * \code
 *   auto a = SiLi::Matrix{{{15, 20, 35},
 *                          {25, 30, 40}}};
 *   auto c = a / 5;
 *   std::cout << c << "\n"; // prints {{3, 4, 7},
 *                                      {5, 6, 8}}
 * \endcode
 */
template <_concept::Matrix L>
constexpr auto operator/(L const& l, value_t<L> const& s) {
	return details::apply(l, [s](auto _l) { return _l / s; });
}


/*! Divide elementwise by a scalar
 * \shortexample l /= s
 * \group Matrix Operations
 *
 * \param l _concept::Matrix
 * \param s scalar
 * \return  reference to l
 *
 * \code
 *   auto a = SiLi::Matrix{{{15, 20, 35},
 *                          {25, 30, 40}}};
 *   a /= 5;
 *   std::cout << c << "\n"; // prints {{3, 4, 7},
 *                                      {5, 6, 8}}
 * \endcode
 */
template <_concept::Matrix L>
constexpr auto operator/=(L&& l, value_t<L> const& s) -> L& {
	details::self_assign_apply(l, [s](auto& _l) { _l /= s; });
	return l;
}

/*! Elementwise comparision
 * \shortexample l == r
 * \group Free Matrix Functions
 *
 * \param l _concept::Matrix
 * \param r _concept::Matrix
 * \return  true if all elements are equal, otherwise false
 *
 * l and r must have same dimensions
 *
 * \code
 *   auto a = SiLi::Matrix{{{15, 20, 35},
 *                          {25, 30, 40}}};
 *
 *   auto b = SiLi::Matrix{{{15, 20, 35},
 *                          {25, 30, 40}}};
 *   auto c = SiLi::Matrix{{{15, 20, 35},
 *                          {25, 30, 43}}};
 *   std::cout << (a == b) << "\n"; // prints 1
 *   std::cout << (a == c) << "\n"; // prints 0
 * \endcode
 */
template <_concept::Matrix L, _concept::Matrix R> requires (L::Rows == R::Rows and L::Cols == R::Cols)
constexpr auto operator==(L const& l, R const& r) -> bool {
	return for_each_constexpr<L>([&]<int row, int col>() {
		return at<row, col>(l) == at<row, col>(r);
	});
}

/*! Elementwise comparision
 * \shortexample l != r
 * \group Free Matrix Functions
 *
 * \param l _concept::Matrix
 * \param r _concept::Matrix
 * \return  false if all element are equal, otherwise true
 *
 * l and r must have same dimensions
 *
 * \code
 *   auto a = SiLi::Matrix{{{15, 20, 35},
 *                          {25, 30, 40}}};
 *
 *   auto b = SiLi::Matrix{{{15, 20, 35},
 *                          {25, 30, 40}}};
 *   auto c = SiLi::Matrix{{{15, 20, 35},
 *                          {25, 30, 43}}};
 *   std::cout << (a != b) << "\n"; // prints 0
 *   std::cout << (a != c) << "\n"; // prints 1
 * \endcode
 */
template <_concept::Matrix L, _concept::Matrix R> requires (L::Rows == R::Rows and L::Cols == R::Cols)
constexpr auto operator!=(L const& l, R const& r) -> bool {
	return not (l == r);
}

/*! Diagonal view
 * \shortexample view_diag(m)
 * \group Free Matrix Functions
 *
 * \param m _concept::Matrix
 * \return  View representing a column vector that gives access to the diagonal of m
 *
 * \code
 *   auto a = SiLi::Matrix{{{15, 20, 35},
 *                          {25, 30, 40}}};
 *   auto v = view_diag(a);
 *   std::cout << rows(v) << "x" << cols(v) << "\n"; // prints 2x1
 *   std::cout << v(0) << " " << v(1) << "\n"; // prints 15 and 30
 *   v(0) = 7;
 *   v(1) = 8;
 *   std::cout << a << "\n"; // prints {{ 7, 20, 35}
 *                                      {25,  8, 40}}
 * \endcode
 */
template <_concept::Matrix M>
constexpr auto view_diag(M&& m) {
	using U = value_t<M>;
	using std::min;
	constexpr auto length = min(cols_v<M>, rows_v<M>);
	constexpr auto stride = stride_v<M>+1;
	return View<length, 1, stride, U, false>{m.data()};
}

/*! Diagonal
 * \shortexample diag(m)
 * \group Free Matrix Functions
 *
 * \param m _concept::Matrix
 * \return  Matrix as a column vector with the diagonal values of m
 *
 * \code
 *   auto a = SiLi::Matrix{{{15, 20, 35},
 *                          {25, 30, 40}}};
 *   auto v = diag(a);
 *   std::cout << rows(v) << "x" << cols(v) << "\n"; // prints 2x1
 *   std::cout << v(0) << " " << v(1) << "\n"; // prints 15 and 30
 * \endcode
 */
template <_concept::Matrix M>
constexpr auto diag(M const& m) {
	return Matrix{view_diag(m)};
}

/*! Super diagonal view
 * \shortexample view_upper_diag(m)
 * \group Free Matrix Functions
 *
 * \param m _concept::Matrix
 * \return  View representing a column vector that gives access to the super diagonal of m
 *
 * The super diagonal is the diagonal above the diagonal
 *
 * \code
 *   auto a = SiLi::Matrix{{{15, 20, 35},
 *                          {25, 30, 40}}};
 *   auto v = view_upper_diag(a);
 *   std::cout << rows(v) << "x" << cols(v) << "\n"; // prints 2x1
 *   std::cout << v(0) << " " << v(1) << "\n"; // prints 20 and 40
 *   v(0) = 7;
 *   v(1) = 8;
 *   std::cout << a << "\n"; // prints {{15,  7, 35}
 *                                      {25, 30,  8}}
 * \endcode
 */
template <_concept::Matrix M>
constexpr auto view_upper_diag(M&& m) {
	return view_diag(view<0, 1, End, End>(m));
}

/*! View
 * \shortexample view<R0, C0, R1, C1>(m)
 * \group Free Matrix Functions
 *
 * \param start_row starting row
 * \param start_col starting column
 * \param end_row   end row (exclusive). Use SiLi::End to indicate rows(v).
 * \param end_col   end column (exclusive). Use SiLi::End to indicate cols(v).
 * \param v         _concept::Matrix
 * \return          View
 *
 *
 * Returns a view onto a matrix
 *
 * \code
 *   auto a = SiLi::Matrix{{{ 1,  2,  3,  4,  5},
 *                          { 6,  7,  8,  9, 10},
 *                          {11, 12, 13, 14, 15},
 *                          {16, 17, 18, 19, 20}}};
 *   auto v0 = view<1, 1, 3, 3>(a);
 *   std::cout << rows(v0) << "x" << cols(v0) << "\n"; // prints 2x2
 *   std::cout << v0 << "\n"; // prints {{ 7,  8},
 *                                        12, 13}
 *
 *   auto v1 = view<2, 3, SiLi::End, SiLi::End>(a);
 *   std::cout << rows(v1) << "x" << cols(v1) << "\n"; // prints 2x2
 *   std::cout << v1 << "\n"; // prints {{14, 15},
 *                                        19, 20}

 * \endcode
 */
template <int start_row, int start_col, int end_row, int end_col, _concept::Matrix M>
    requires (end_row <= rows_v<M> and end_col <= cols_v<M>)
constexpr auto view(M&& m) {
	using U = value_t<M>;

	constexpr auto rows       = end_row - start_row;
	constexpr auto cols       = end_col - start_col;
	constexpr auto stride     = stride_v<M>;
	constexpr auto transposed = transposed_v<M>;

	if constexpr (not transposed) {
		return View<rows, cols, stride, U, transposed>{m.data() + start_col + start_row * stride};
	} else {
		return View<rows, cols, stride, U, transposed>{m.data() + start_row + start_col * stride};
	}
}

template <int start_row, int start_col, CEnd end_row, int end_col, _concept::Matrix V>
constexpr auto view(V&& v) {
	return view<start_row, start_col, rows_v<V>, end_col>(std::forward<V>(v));
}
template <int start_row, int start_col, int end_row, CEnd end_col, _concept::Matrix V>
constexpr auto view(V&& v) {
	return view<start_row, start_col, end_row, cols_v<V>>(std::forward<V>(v));
}
template <int start_row, int start_col, CEnd end_row, CEnd end_col, _concept::Matrix V>
constexpr auto view(V&& v) {
	return view<start_row, start_col, rows_v<V>, cols_v<V>>(std::forward<V>(v));
}

template <int start_row, int end_row, _concept::Matrix V> requires (end_row <= rows_v<V> and cols_v<V> == 1)
constexpr auto view(V&& v) {
	return view<start_row, 0, end_row, 1>(v);
}

template <int start_col, int end_col, _concept::Matrix V> requires (end_col <= cols_v<V> and rows_v<V> == 1 and cols_v<V> != 1)
constexpr auto view(V&& v) {
	return view<0, start_col, 1, end_col>(v);
}

template <int start, CEnd end, _concept::Matrix V>
constexpr auto view(V&& v) {
	return view<start, length_v<V>>(std::forward<V>(v));
}


template <int _row, _concept::Matrix V>
constexpr auto view_row(V&& v) {
	return view<_row, 0, _row+1, cols_v<V>>(std::forward<V>(v));
}

template <int _col, _concept::Matrix V>
constexpr auto view_col(V&& v) {
	return view<0, _col, rows_v<V>, _col+1>(std::forward<V>(v));
}

/*! Joins rows
 * \shortexample join_rows(l, r)
 * \group Free Matrix Functions
 *
 * \param l _concept::Matrix
 * \param r _concept::Matrix
 * \return  Matrix of joined rows of l and r
 *
 * l and r must have same number of rows.
 *
 * \code
 *   auto a = SiLi::Matrix{{{ 1,  2,  3},
 *                          {11, 12, 13}}};
 *   auto b = SiLi::Matrix{{{ 4,  5,  6},
 *                          {14, 15, 16}}};
*    auto v = join_rows(a, b);
 *   std::cout << rows(v) << "x" << cols(v) << "\n"; // prints 2x6
 *   std::cout << v << "\n"; // prints {{ 1,  2,  3,  4,  5,  6},
 *                                      {11, 12, 13, 14, 15, 16}}
 * \endcode
 */
template <_concept::Matrix L, _concept::Matrix R> requires (L::Rows == R::Rows)
constexpr auto join_rows(L const& l, R const& r) {
	auto matrix = Matrix<L::Rows, L::Cols + R::Cols, typename L::value_t>{};
	view<0,       0, L::Rows,         L::Cols>(matrix) = l;
	view<0, L::Cols, R::Rows, L::Cols+R::Cols>(matrix) = r;
	return matrix;
}

/*! Joins columns
 * \shortexample join_cols(a, b)
 * \group Free Matrix Functions
 *
 * \param l Viewable
 * \param r Viewable
 * \return  Matrix of joined columns of l and r
 *
 * l and r must have same number of columns.
 *
 * \code
 *   auto a = SiLi::Matrix{{{ 1,  2,  3},
 *                          {11, 12, 13}}};
 *   auto b = SiLi::Matrix{{{ 4,  5,  6},
 *                          {14, 15, 16}}};
 *   auto v = join_cols(a, b);
 *   std::cout << rows(v) << "x" << cols(v) << "\n"; // prints 4x3
 *   std::cout << v << "\n"; // prints {{ 1,  2,  3},
 *                                      {11, 12, 13},
 *                                      { 4,  5,  6},
 *                                      {14, 15, 16}}
 * \endcode
 */
template <_concept::Matrix L, _concept::Matrix R> requires (L::Cols == R::Cols)
constexpr auto join_cols(L const& l, R const& r) {
	auto matrix = Matrix<L::Rows+R::Rows, L::Cols, typename L::value_t>{};
	view<      0, 0,         L::Rows, L::Cols>(matrix) = l;
	view<L::Rows, 0, L::Rows+R::Rows, R::Cols>(matrix) = r;
	return matrix;
}

// !TODO needs unit tests
// lu decomposition, returns L value
template<_concept::Matrix V> requires (rows_v<V> == cols_v<V>)
constexpr auto luDecomposition_L(V const& v) {
	using T = std::decay_t<value_t<V>>;

	constexpr auto N = rows_v<V>;
	auto L = Matrix<N, N, T>{};

	for_constexpr<0, N>([&]<auto K>() {
		auto kRow = view_row<K>(L);
		auto kCol = view_col<K>(L);
		for_constexpr<K, N>([&]<auto I>() {
			at<K, I>(L) = at<K, I>(v) - T{kRow * view_col<I>(L)};
		});
		for_constexpr<K+1, N>([&]<auto I>() {
			at<I, K>(L) = (at<I, K>(v) - T{view_row<I>(L) * kCol}) / at<K, K>(L);
		});
	});
	return L;
}


// compute 1x1 determinant
template <_concept::Matrix V> requires (V::Rows == 1 and V::Cols == 1)
constexpr auto det(V const& v) {
	return at<0, 0>(v);
}

// compute 2x2 determinant
template <_concept::Matrix V> requires (V::Rows == 2 and V::Cols == 2)
constexpr auto det(V const& v) {
	return at<0, 0>(v)*at<1, 1>(v) - at<0, 1>(v)*at<1, 0>(v);
}

// compute 3x3 determinant
template <_concept::Matrix V> requires (V::Rows == 3 and V::Cols == 3)
constexpr auto det(V const& v) {
	return   (at<0, 0>(v)*at<1, 1>(v)*at<2, 2>(v)
	        + at<0, 1>(v)*at<1, 2>(v)*at<2, 0>(v)
	        + at<0, 2>(v)*at<1, 0>(v)*at<2, 1>(v))
	       - (at<0, 2>(v)*at<1, 1>(v)*at<2, 0>(v)
	        + at<0, 1>(v)*at<1, 0>(v)*at<2, 2>(v)
	        + at<0, 0>(v)*at<1, 2>(v)*at<2, 1>(v));
}

/*! Compute determinant
 * \shortexample det(m)
 * \group Free Matrix Functions
 *
 * \param m _concept::Matrix must have same number of rows and columns
 * \return  Determinant of m
 *
 * \code
 *   auto a = SiLi::Matrix{{{ 1,  2},
 *                          {11, 12}}};
 *   auto v = det(a);
 *   std::cout << v << "\n"; // prints -10
 * \endcode
 */
template <_concept::Matrix M> requires (M::Rows == M::Cols and M::Rows > 3)
constexpr auto det(M const& m) {
	using T = typename M::value_t;
	auto retValue = T{1};

	auto L        = luDecomposition_L(m);

	for_constexpr<0, M::Rows>([&]<int i>() {
		retValue *= at<i, i>(L);
	});

	using std::isfinite;
	if (not isfinite(retValue)) {
		retValue = T{0};
	}
	return retValue;
}

//element wise product
template <_concept::Matrix L, _concept::Matrix R> requires (L::Rows == R::Rows and L::Cols == R::Cols)
constexpr auto element_multi(L const& l, R const& r) {
	return details::apply(l, r, [](auto _l, auto _r) constexpr { return _l * _r; });
}

/*! Cross product of two vectors.
 * \shortexample cross(l, r)
 * \group Free Matrix Functions
 *
 * \param l _concept::Vector of length 3.
 * \param r _concept::Vector of length.
 * \return  Matrix of size 3×1 with cross product of l and r.
 *
 * \code
 *   auto a = SiLi::Matrix{{{ 1},
 *                          { 2},
 *                          { 3}}};
 *
 *   auto b = SiLi::Matrix{{{ 4},
 *                          { 5},
 *                          { 6}}};

 *   auto v = cross(a, b);
 *   std::cout << v << "\n"; // prints {{-3},
 *                                      { 6},
 *                                      {-3}}
 * \endcode
 */
template <_concept::Vector L, _concept::Vector R> requires (length_v<R> == 3 and length_v<L> == 3)
constexpr auto cross(L const& l, R const& r) {
	using U = decltype(std::declval<typename L::value_t>() * std::declval<typename R::value_t>());
	auto ret = Matrix<3, 1, U>{};
	ret(0) = l(1) * r(2) - l(2) * r(1);
	ret(1) = l(2) * r(0) - l(0) * r(2);
	ret(2) = l(0) * r(1) - l(1) * r(0);
	return ret;
}

/*! Sum of all elements.
 * \shortexample sum(m)
 * \group Free Matrix Functions
 *
 * \param m _concept::Matrix
 * \return  Sum off all elements of m.
 *
 * \code
 *   auto a = SiLi::Matrix{{{3, 4, 7},
 *                          {5, 6, 8}}};
 *   auto c = sum(a)
 *   std::cout << c << "\n"; // prints 33
 * \endcode
 */
template <_concept::Matrix M>
constexpr auto sum(M const& m) {
	auto acc = value_t<M>{};
	for_each_constexpr<M>([&]<auto row, auto col>() {
		acc += at<row, col>(m);
	});
	return acc;
}

/*! Sum of each row.
 * \shortexample sum_rows(m)
 * \group Free Matrix Functions
 *
 * \param m _concept::Matrix
 * \return  Column vector with the sums of each row of m
 *
 * \code
 *   auto a = SiLi::Matrix{{{3, 4, 7},
 *                          {5, 6, 8}}};
 *   auto c = sum_rows(a)
 *   std::cout << c << "\n"; // prints {{14}
 *                                      {19}}
 * \endcode
 */
template <_concept::Matrix M>
constexpr auto sum_rows(M const& m) {
	using T = value_t<M>;
	auto c = Matrix<rows_v<M>, 1, T>{};
	for_constexpr<0, rows_v<M>>([&]<int row>() {
		at<row>(c) = sum(view_row<row>(m));
	});
	return c;
}

/*! Sum of each column.
 * \shortexample sum_cols(m)
 * \group Free Matrix Functions
 *
 * \param m _concept::Matrix
 * \return  Row vector with the sums of each column of m
 *
 * \code
 *   auto a = SiLi::Matrix{{{3, 4, 7},
 *                          {5, 6, 8}}};
 *   auto c = sum_cols(a)
 *   std::cout << c << "\n"; // prints {{8, 10, 15}}
 * \endcode
 */
template <_concept::Matrix M>
constexpr auto sum_cols(M const& m) {
	using T = value_t<M>;
	auto c = Matrix<1, cols_v<M>, T>{};
	for_constexpr<0, cols_v<M>>([&]<int col>() {
		at<col>(c) = sum(view_col<col>(m));
	});
	return c;
}



// inverse of 1x1
template <_concept::Matrix V> requires (V::Rows == V::Cols and V::Rows == 1)
constexpr auto inv(V v) -> std::tuple<typename V::value_t, V> {
	using T = typename V::value_t;
	using std::abs;
	auto d = det(v);
	if (abs(d) < 1.e-5) {
		return {at<0, 0>(v), v};
	}

	at<0, 0>(v) = T(1) / at<0, 0>(v);
	return {d, v};
}

// inverse of 2x2
template <_concept::Matrix V> requires (V::Rows == V::Cols and V::Rows == 2)
constexpr auto inv(V v) -> std::tuple<typename V::value_t, V> {
	using T = typename V::value_t;
	using std::abs;
	auto d = det(v);
	if (abs(d) < 1.e-5) {
		return {at<0, 0>(v), v};
	}
	auto c = T(1) / d;
	return {d, Matrix{{{ at<1, 1>(v)*c, -at<0, 1>(v)*c},
	                   {-at<1, 0>(v)*c,  at<0, 0>(v)*c}}}};

}

//inverse of 3x3
template <_concept::Matrix M> requires (M::Rows == M::Cols and M::Rows == 3)
constexpr auto inv(M const& m) -> std::tuple<typename M::value_t, M> {
	using T = typename M::value_t;
	using std::abs;
	auto d = det(m);
	if (abs(d) < 1.e-5) {
		return {at<0, 0>(m), m};
	}
	auto c = T(1) / d;
	auto ret = M{};

	for_each_constexpr<M>([&]<int row, int col>() {
		auto tl = at<(row+1)%3, (col+1)%3>(m);
		auto br = at<(row+2)%3, (col+2)%3>(m);
		auto tr = at<(row+1)%3, (col+2)%3>(m);
		auto bl = at<(row+2)%3, (col+1)%3>(m);
		at<col, row>(ret) = (tl*br-tr*bl) * c;
	});
	return {d, ret};
}

/*! Compute inverse
 * \shortexample inv(m)
 * \group Free Matrix Functions
 *
 * \param  m _concept::Matrix must have same number of rows as number of cols.
 * \return   Returns tuple of determinant and inverse, if detereminant is zero the inverse is invalid.
 *
 * \code
 *   auto a = SiLi::Matrix{{{ 1.,  2.},
 *                          {11., 12.}}};
 *   auto [v, m] = inv(a);
 *   std::cout << v << "\n"; // prints -10
 *   std::cout << m << "\n"; // prints {{-1.2,  0.2},
 *                                      { 1.1, -0.1}}
 * \endcode
 */
template <_concept::Matrix M> requires (M::Rows == M::Cols and M::Rows > 3)
constexpr auto inv(M m) -> std::tuple<typename M::value_t, M> {
	using T = typename M::value_t;
	constexpr int N = M::Rows;
	auto det        = T{1};

	bool failed{false};
	for_constexpr<0, N>([&]<int p>() -> bool {
		auto pivot = at<p, p>(m);
		det = det * pivot;
		if (fabs(pivot) < 1.e-5) {
			failed = true;
			return false;
		}
		view_col<p>(m) /= -pivot;
		for_each_constexpr<M>([&]<int i, int j>() {
			if (i == p or j == p) return;
			at<i, j>(m) += at<p, j>(m) * at<i, p>(m);
		});
		view_row<p>(m) /= pivot;
		at<p, p>(m) = T{1}/ pivot;
		return true;
	});
	if (failed) {
		return {0., m};
	}
	return {det, m};
}

/*! Transposed view
 * \shortexample view_trans(m)
 * \group Free Matrix Functions
 *
 * \param  m _concept::Matrix
 * \return   Matrix of the transposed of m
 *
 * \code
 *   auto a = SiLi::Matrix{{{ 1,  2,  3},
 *                          {11, 12,  4}}};
 *   auto v = view_trans(a);
 *   std::cout << m << "\n"; // prints {{1, 11},
 *                                      {2, 12},
 *                                      {3, 13}}
 * \endcode
 */
template <_concept::Matrix M>
constexpr auto view_trans(M&& m) {
	return View<cols_v<M>, rows_v<M>, stride_v<M>, value_t<M>, not transposed_v<M>>{m.data()};
}

/*! Transposed
 * \shortexample trans(m)
 * \group Free Matrix Functions
 *
 * \param  m _concept::Matrix
 * \return   Matrix of the transposed of m
 *
 * \code
 *   auto a = SiLi::Matrix{{{ 1,  2,  3},
 *                          {11, 12,  4}}};
 *   auto v = trans(a);
 *   std::cout << m << "\n"; // prints {{1, 11},
 *                                      {2, 12},
 *                                      {3, 13}}
 * \endcode
 */
template <_concept::Matrix M>
constexpr auto trans(M&& m) {
	return Matrix{view_trans(m)};
}

/*! Identity matrix
 * \shortexample SiLi::makeI<Rows, Cols, int>()
 * \group Free Matrix Functions
 *
 * \param Rows number of rows
 * \param Cols number of columns
 * \param T    type of each element
 * \return     matrix with ones on the diagonal
 *
 * \code
 *   auto a = SiLi::makeI<3, 4, int>();
 *   std::cout << a << "\n"; // prints {{1, 0, 0, 0},
 *                                      {0, 1, 0, 0},
 *                                      {0, 0, 1, 0}}
 * \endcode
 */
template <int Rows, int Cols, typename T>
constexpr auto makeI() {
	auto i = Matrix<Rows, Cols, T>{};
	view_diag(i) = T{1};
	return i;
}

/*! Identity matrix
 * \shortexample SiLi::makeI<N, int>()
 * \group Free Matrix Functions
 *
 * \param N size of the matrix
 * \param T type of each element
 * \return  matrix with ones on the diagonal
 *
 * \code
 *   auto a = SiLi::makeI<3, int>();
 *   std::cout << a << "\n"; // prints {{1, 0, 0},
 *                                      {0, 1, 0},
 *                                      {0, 0, 1}}
 * \endcode
 */
template <int N, typename T>
constexpr auto makeI() {
	return makeI<N, N, T>();
}


/*! Compute norm
 * \shortexample norm(v)
 * \group Free Vector Functions
 *
 * \param v _concept::Vector
 * \return  the norm (also considered as the length).
 *
 * \code
 *   auto a = SiLi::Matrix{{{ 3.},
 *                          { 4.},
 *                          { 5.}}};
 *
 *   auto v = norm(a);
 *   std::cout << v << "\n"; // prints 7.0711
 * \endcode
 */
template <_concept::Vector V>
constexpr auto norm(V const& v) {
	auto acc = value_t<V>{};
	for_each_constexpr<V>([&]<auto row, auto col>() {
		acc += v(row, col)*v(row, col);
	});
	using std::sqrt;
	return sqrt(acc);
}

/*! Compute abs
 * \shortexample abs(m)
 * \group Free Matrix Functions
 *
 * \param m _concept::Matrix
 * \return  the absolute values
 *
 * \code
 *   auto a = SiLi::Matrix{{{ -3.},
 *                          { -4.},
 *                          { -5.}}};
 *
 *   auto v = abs(a);
 *   std::cout << v << "\n"; // prints {{3},
 *                           // {4},
 *                           // {5}}
 * \endcode
 */
template <_concept::Matrix M>
constexpr auto abs(M const& m) {
	using std::abs;
	return details::apply(m, [](auto e) constexpr { return abs(e); });
}


/*! Dot product of two vectors
 * \shortexample dot(l, r)
 * \group Free Vector Functions
 *
 * \param l _concept::Vector
 * \param r _concept::Vector
 * \return  the dot product (also called scalar product) of l and r. l and r must be vectors of the same length
 *
 * \code
 *   auto a = SiLi::Matrix{{{ 1},
 *                          { 2},
 *                          { 3}}};
 *
 *   auto b = SiLi::Matrix{{{ 4},
 *                          { 5},
 *                          { 6}}};
 *
 *   auto v = dot(a, b);
 *   std::cout << v << "\n"; // prints 32
 * \endcode
 */
template <_concept::Vector L, _concept::Vector R> requires(length_v<L> == length_v<R>)
constexpr auto dot(L const& l, R const& r) {
	using U = decltype(value<L>() * value_t<R>());
	auto acc = U{};
	for_constexpr<0, length_v<L>>([&]<int I>() {
		acc += at<I>(l) * at<I>(r);
	});
	return acc;
}

/*! Outer product
 * \shortexample outerProd(l, r)
 * \group Free Vector Functions
 *
 * \param l _concept::Vector
 * \param r _concept::Vector
 * \return  Matrix with the outer product of l and r
 *
 * \code
 *   auto a = SiLi::Matrix{{{ 1},
 *                          { 2},
 *                          { 3}}};
 *
 *   auto b = SiLi::Matrix{{{ 4},
 *                          { 5},
 *                          { 6}}};
 *
 *   auto v = outerProd(a, b);
 *   std::cout << v << "\n"; // prints {{ 4,  5,  6},
 *                           //         { 8, 10, 12},
 *                           //         {12, 15, 18}}
 * \endcode
 */
template <_concept::Vector L, _concept::Vector R>
constexpr auto outerProd(L const& l, R const& r) {
	if constexpr (cols_v<L> != 1) {
		return outerProd(view_trans(l), r);
	} else if constexpr (rows_v<R> != 1) {
		return outerProd(l, view_trans(r));
	} else {
		return l * r;
	}
}

}
