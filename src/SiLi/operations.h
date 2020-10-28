/*!\mainpage My Personal Index Page
 *
 * \section intro_sec Introduction
 *
 * This is the introduction.
 *
 * \section install_sec Installation
 *
 * \subsection step1 Step 1: Opening the box
 */
/** \addtogroup freefunctions
 *
 *  Free functions
 */

#pragma once

#include "concepts.h"

#include <cmath>
#include <optional>

namespace SiLi {

struct End {};
template <typename T>
concept CEnd = std::is_same_v<End, T>;

template <typename T>
constexpr auto value() -> std::add_rvalue_reference<value_t<T>>::type;

/*!\brief Access element of matrix.
 * \ingroup freefunctions
 *
 * Accessing element of matrix. Checking accessibility at compile time.
 *
 * \param Row access
 * \param Col access
 * \param v   viewable
 * \return    returns reference to underlying element
 *
 * #### Example ####
 * \code
 *   auto m = SiLi::Matrix{{{3, 4},
 *                          {5, 6}}};
 *   at<0, 0>(m) = 5;
 *   at<0, 1>(m) = 6;
 *   std::cout << m << "\n"; // prints {{5, 6},
 *                                      {5, 6}}
 * \endcode
 */
template <int Row, int Col, Viewable V>
constexpr auto at(V&& v) -> auto& {
	static_assert(0 <= Row and Row < rows_v<V>, "accessing outside of valid range");
	static_assert(0 <= Col and Col < cols_v<V>, "accessing outside of valid range");
	return v(Row, Col);
}

/*!\brief Access element of matrix.
 * \ingroup freefunctions
 *
 * Accessing element of vector. Checking accessibility at compile time.
 *
 * \param Idx row access
 * \param v   viewable must have columns==1 or rows==1
 * \return    returns reference to underlying element
 *
 *  #### Example ####
 * \code
 *   auto m = SiLi::Matrix{{{3},
 *                          {5}}};
 *   at<0>(m) = 5;
 *   at<1>(m) = 6;
 *   std::cout << m << "\n"; // prints {{5},
 *                                      {6}}
 * \endcode
 */
template <int Idx, ViewableVector V>
constexpr auto at(V&& v) -> auto& {
	if constexpr (cols_v<V> == 1) {
		return at<Idx, 0>(v);
	} else {
		return at<0, Idx>(v);
	}
}

/*!\brief Retrive number of rows
 * \ingroup freefunctions
 *
 * \param v viewable
 * \return  numbers of rows
 *
 *  #### Example ####
 * \code
 *   auto m = SiLi::Matrix{{{3, 4, 7},
 *                          {5, 6, 8}}};
 *   std::cout << rows(m) << "\n"; // prints 2
 * \endcode
 */
template <Viewable V>
constexpr auto rows(V const& v) {
	return rows_v<V>;
}


/*!\brief Retrive numbers of columns
 * \ingroup freefunctions
 *
 * \param  v viewable
 * \return   number of columns
 *
 *  #### Example ####
 * \code
 *   auto m = SiLi::Matrix{{{3, 4, 7},
 *                          {5, 6, 8}}};
 *   std::cout << cols(m) << "\n"; // prints 3
 * \endcode
 */
template <Viewable V>
constexpr auto cols(V const& v) {
	return cols_v<V>;
}

namespace details {
template <Viewable V, typename Operator>
constexpr auto apply(V const& v, Operator op) {
	using U = decltype(op(value<V>()));
	auto ret = Matrix<rows_v<V>, cols_v<V>, U>{};
	for_each_constexpr<V>([&]<auto row, auto col>() {
		at<row, col>(ret) = op(at<row, col>(v));
	});
	return ret;
}

template <Viewable L, Viewable R, typename Operator> requires (rows_v<L> == rows_v<R> and cols_v<L> == cols_v<R>)
constexpr auto apply(L const& l, R const& r, Operator op) {
	using U = decltype(op(value<L>(), value_t<R>()));
	auto ret = Matrix<rows_v<L>, cols_v<L>, U>{};
	for_each_constexpr<L>([&]<auto row, auto col>() {
		at<row, col>(ret) = op(at<row, col>(l), at<row, col>(r));
	});
	return ret;
}

template <Viewable V, typename Operator>
constexpr auto self_assign_apply(V&& v, Operator op) -> auto& {
	for_each_constexpr<V>([&]<auto row, auto col>() {
		op(at<row, col>(v));
	});
	return v;
}

template <Viewable L, Viewable R, typename Operator> requires (rows_v<L> == rows_v<R> and cols_v<L> == cols_v<R>)
constexpr auto self_assign_apply(L&& l, R const& r, Operator op) -> auto& {
	for_each_constexpr<L>([&]<auto row, auto col>() {
		op(at<row, col>(l), at<row, col>(r));
	});
	return l;
}
}

/*!\brief return copy of v
 * \ingroup freefunctions
 *
 * \param v viewable
 * \return copy of v with `+` applied to every element
 *
 *  #### Example ####
 * \code
 *   auto m = SiLi::Matrix{{{3, 4, 7},
 *                          {5, 6, 8}}};
 *   auto m2 = +m;
 *   std::cout << m << "\n";
 *   std::cout << m2 << "\n";
 * \endcode
 */
template <Viewable V>
constexpr auto operator+(V const& v) {
	return details::apply(v, [](auto v) { return +v; });
}

/*!\brief Elementwise addition
 * \ingroup freefunctions
 *
 * \param l left matrix of addition
 * \param r right matrix of addition
 * \return  matrix with l and r
 *
 * l and r must have the same dimension.
 *
 *  #### Example ####
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
template <Viewable L, Viewable R>
constexpr auto operator+(L const& l, R const& r) {
	return details::apply(l, r, [](auto l, auto r) { return l + r; });
}

/*!\brief Elementwis addition
 * \ingroup freefunctions
 *
 * \param l left Viewable
 * \param r right Viewable
 * \return  reference to l
 *
 * l and r must have the same dimension.
 *
 *  #### Example ####
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
template <Viewable L, Viewable R>
constexpr auto operator+=(L&& l, R const& r) -> auto& {
	return details::self_assign_apply(l, r, [](auto& l, auto r) { l += r; });
}

/*!\brief Elementwise negation of elements
 * \ingroup freefunctions
 *
 * \param v Viewable
 * \return  matrix with negated elements
 *
 *  #### Example ####
 * \code
 *   auto m = SiLi::Matrix{{{3, 4, 7},
 *                          {5, 6, 8}}};
 *   auto m2 = -m;
 *   std::cout << m << "\n";
 *   std::cout << m2 << "\n";
 * \endcode
 */
template <Viewable V>
constexpr auto operator-(V const& v) {
	return details::apply(v, [](auto v) { return -v; });
}

/*!\brief Elementwise subtraction
 * \ingroup freefunctions
 *
 * \param l left Viewable
 * \param r right Viewable
 * \return  matrix with (l - r)
 *
 * l and r must have the same dimension.
 *
 *  #### Example ####
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
template <Viewable L, Viewable R>
constexpr auto operator-(L const& l, R const& r) {
	return details::apply(l, r, [](auto l, auto r) { return l - r; });
}

/*!\brief Elementwise subtraction
 * \ingroup freefunctions
 *
 * \param l left Viewable
 * \param r right Viewable
 * \return  reference to l
 *
 * l and r must have the same dimension.
 *
 *  #### Example ####
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
template <Viewable L, Viewable R>
constexpr auto operator-=(L& l, R const& r) -> auto& {
	return details::self_assign_apply(l, r, [](auto& l, auto r) { l -= r; });
}

/*!\brief Matrix multiplication.
 * \ingroup freefunctions
 *
 * \param l Viewable
 * \param r Viewable
 * \return Matrix with result of l and r multiplication
 *
 * l must have the same number of columns as r has rows.
 *
 *  #### Example ####
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
template <Viewable L, Viewable R> requires (L::Cols == R::Rows)
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
				ret(row, col) = U{view_row<row>(l) * view_col<col>(r)};
			});
		});
		return ret;
	}
}

/*!\brief Scalar multiplication
 * \ingroup freefunctions
 *
 * \param l Viewable
 * \param s scalar value
 * \return Matrix with each element multiplied with s
 *
 *  #### Example ####
 * \code
 *   auto a = SiLi::Matrix{{{3, 4, 7},
 *                          {5, 6, 8}}};
 *   auto c = a * 5;
 *   std::cout << c << "\n"; // prints {{15, 20, 35},
 *                                      {25, 30, 40}}
 * \endcode
 */
template <Viewable L>
constexpr auto operator*(L const& l, value_t<L> const& s) {
	return details::apply(l, [s](auto l) { return l * s; });
}

/*!\brief Scalar multiplication
 * \ingroup freefunctions
 *
 * \param s scalar value
 * \param r Viewable
 * \return Matrix with each element multiplied with s
 *
 *  #### Example ####
 * \code
 *   auto a = SiLi::Matrix{{{3, 4, 7},
 *                          {5, 6, 8}}};
 *   auto c = 5 * a;
 *   std::cout << c << "\n"; // prints {{15, 20, 35},
 *                                      {25, 30, 40}}
 * \endcode
 */
template <Viewable R>
constexpr auto operator*(value_t<R> const& s, R&& r) {
	return std::forward<R>(r) * s;
}

/*!\brief Scalar multiplication
 * \ingroup freefunctions
 *
 * \param l Viewable
 * \param s scalar value
 * \return Referenc to l
 *
 *  #### Example ####
 * \code
 *   auto a = SiLi::Matrix{{{3, 4, 7},
 *                          {5, 6, 8}}};
 *   auto a *= 5;
 *   std::cout << a << "\n"; // prints {{15, 20, 35},
 *                                      {25, 30, 40}}
 * \endcode
 */
template <Viewable L>
constexpr auto operator*=(L&& l, value_t<L> const& s) -> auto& {
	return details::self_assign_apply(l, [s](auto& l) { l *= s; });
}

/*!\brief Divide elementwise with scalar
 * \ingroup freefunctions
 *
 * \param l Viewable
 * \param s scalar
 * \return Matrix same size and type as l
 *
 *  #### Example ####
 * \code
 *   auto a = SiLi::Matrix{{{15, 20, 35},
 *                          {25, 30, 40}}};
 *   auto c = a / 5;
 *   std::cout << c << "\n"; // prints {{3, 4, 7},
 *                                      {5, 6, 8}}
 * \endcode
 */
template <Viewable L>
constexpr auto operator/(L const& l, value_t<L> const& s) {
	return details::apply(l, [s](auto l) { return l / s; });
}


/*!\brief Divide elementwise with scalar
 * \ingroup freefunctions
 *
 * \param l Viewable
 * \param s scalar
 * \return  reference to l
 *
 *  #### Example ####
 * \code
 *   auto a = SiLi::Matrix{{{15, 20, 35},
 *                          {25, 30, 40}}};
 *   a /= 5;
 *   std::cout << c << "\n"; // prints {{3, 4, 7},
 *                                      {5, 6, 8}}
 * \endcode
 */
template <Viewable L>
constexpr auto operator/=(L&& l, value_t<L> const& s) -> L& {
	return details::self_assign_apply(l, [s](auto& l) { l /= s; });
}

/*!\brief Elementwise comparision of two matrices
 * \ingroup freefunctions
 *
 * \param l Viewable
 * \param r Viewable
 * \return  true if all elements are equal, otherwise false
 *
 * l and r must have same dimensions
 *
 *  #### Example ####
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
template <Viewable L, Viewable R> requires (L::Rows == R::Rows and L::Cols == R::Cols)
constexpr auto operator==(L const& l, R const& r) -> bool {
	return for_each_constexpr<L>([&]<int row, int col>() {
		return at<row, col>(l) == at<row, col>(r);
	});
}

/*!\brief Elementwise comparision of two matrices
 * \ingroup freefunctions
 *
 * \param l Viewable
 * \param r Viewable
 * \return  false if all element are equal, otherwise true
 *
 * l and r must have same dimensions
 *
 *  #### Example ####
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
template <Viewable L, Viewable R> requires (L::Rows == R::Rows and L::Cols == R::Cols)
constexpr auto operator!=(L const& l, R const& r) -> bool {
	return not (l == r);
}

/*!\brief Create view of diagonal
 * \ingroup freefunctions
 *
 * \param v Viewable
 * \return  Viewable representing a column vector that gives access to the diagonal of v
 *
 *  #### Example ####
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
template <Viewable V>
constexpr auto view_diag(V&& v) {
	using U = value_t<V>;
	using std::min;
	constexpr auto length = min(cols_v<V>, rows_v<V>);
	constexpr auto stride = stride_v<V>+1;
	return View<length, 1, stride, U, false>{v.data()};
}

/*!\brief Create diagonal
 * \ingroup freefunctions
 *
 * \param v Viewable
 * \return  Matrix as a column vector with the diagonal values of v
 *
 *  #### Example ####
 * \code
 *   auto a = SiLi::Matrix{{{15, 20, 35},
 *                          {25, 30, 40}}};
 *   auto v = diag(a);
 *   std::cout << rows(v) << "x" << cols(v) << "\n"; // prints 2x1
 *   std::cout << v(0) << " " << v(1) << "\n"; // prints 15 and 30
 * \endcode
 */
template <Viewable V>
constexpr auto diag(V const& v) {
	return Matrix{view_diag(v)};
}

/*!\brief Create view of super diagonal
 * \ingroup freefunctions
 *
 * \param v Viewable
 * \return  Viewable representing a column vector that gives access to the super diagonal of v
 *
 * The super diagonal is the diagonal above the diagonal
 *
 *  #### Example ####
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
template <Viewable V>
constexpr auto view_upper_diag(V&& v) {
	return view_diag(view<0, 1, End, End>(v));
}

/*!\brief Create view of matrix
 * \ingroup freefunctions
 *
 * \param start_row starting row
 * \param start_col starting column
 * \param end_row   end row (exclusive). Use SiLi::End to indicate rows(v).
 * \param end_col   end column (exclusive). Use SiLi::End to indicate cols(v).
 * \param v         Viewable
 * \return          Viewable
 *
 *
 * Returns a view onto a matrix
 *
 *  #### Example ####
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
template <int start_row, int start_col, int end_row, int end_col, Viewable V>
    requires (end_row <= rows_v<V> and end_col <= cols_v<V>)
constexpr auto view(V&& v) {
	using U = value_t<V>;

	constexpr auto rows       = end_row - start_row;
	constexpr auto cols       = end_col - start_col;
	constexpr auto stride     = stride_v<V>;
	constexpr auto transposed = transposed_v<V>;

	if constexpr (not transposed) {
		return View<rows, cols, stride, U, transposed>{v.data() + start_col + start_row * stride};
	} else {
		return View<rows, cols, stride, U, transposed>{v.data() + start_row + start_col * stride};
	}
}

template <int start_row, int start_col, CEnd end_row, int end_col, Viewable V>
constexpr auto view(V&& v) {
	return view<start_row, start_col, rows_v<V>, end_col>(std::forward<V>(v));
}
template <int start_row, int start_col, int end_row, CEnd end_col, Viewable V>
constexpr auto view(V&& v) {
	return view<start_row, start_col, end_row, cols_v<V>>(std::forward<V>(v));
}
template <int start_row, int start_col, CEnd end_row, CEnd end_col, Viewable V>
constexpr auto view(V&& v) {
	return view<start_row, start_col, rows_v<V>, cols_v<V>>(std::forward<V>(v));
}

template <int start_row, int end_row, Viewable V> requires (end_row <= rows_v<V> and cols_v<V> == 1)
constexpr auto view(V&& v) {
	return view<start_row, 0, end_row, 1>(v);
}

template <int start_col, int end_col, Viewable V> requires (end_col <= cols_v<V> and rows_v<V> == 1 and cols_v<V> != 1)
constexpr auto view(V&& v) {
	return view<0, start_col, 1, end_col>(v);
}

template <int start, CEnd end, Viewable V>
constexpr auto view(V&& v) {
	return view<start, length_v<V>>(std::forward<V>(v));
}


template <int _row, Viewable V>
constexpr auto view_row(V&& v) {
	return view<_row, 0, _row+1, cols_v<V>>(std::forward<V>(v));
}

template <int _col, Viewable V>
constexpr auto view_col(V&& v) {
	return view<0, _col, rows_v<V>, _col+1>(std::forward<V>(v));
}

/*!\brief Joins rows of two Viewables.
 * \ingroup freefunctions
 *
 * \param l Viewable
 * \param r Viewable
 * \return  Matrix of joined rows of l and r
 *
 * l and r must have same number of rows.
 *
 *  #### Example ####
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
template <Viewable L, Viewable R> requires (L::Rows == R::Rows) 
constexpr auto join_rows(L const& l, R const& r) {
	auto matrix = Matrix<L::Rows, L::Cols + R::Cols, typename L::value_t>{};
	view<0,       0, L::Rows,         L::Cols>(matrix) = l;
	view<0, L::Cols, R::Rows, L::Cols+R::Cols>(matrix) = r;
	return matrix;
}

/*!\brief Joins columns of two Viewables.
 * \ingroup freefunctions
 *
 * \param l Viewable
 * \param r Viewable
 * \return  Matrix of joined columns of l and r
 *
 * l and r must have same number of columns.
 *
 *  #### Example ####
 * \code
 *   auto a = SiLi::Matrix{{{ 1,  2,  3},
 *                          {11, 12, 13}}};
 *   auto b = SiLi::Matrix{{{ 4,  5,  6},
 *                          {14, 15, 16}}};
 *   auto v = join_rows(a, b);
 *   std::cout << rows(v) << "x" << cols(v) << "\n"; // prints 4x3
 *   std::cout << v << "\n"; // prints {{ 1,  2,  3},
 *                                      {11, 12, 13},
 *                                      { 4,  5,  6},
 *                                      {14, 15, 16}}
 * \endcode
 */
template <Viewable L, Viewable R> requires (L::Cols == R::Cols)
constexpr auto join_cols(L const& l, R const& r) {
	auto matrix = Matrix<L::Rows+R::Rows, L::Cols, typename L::value_t>{};
	view<      0, 0,         L::Rows, L::Cols>(matrix) = l;
	view<L::Rows, 0, L::Rows+R::Rows, R::Cols>(matrix) = r;
	return matrix;
}

// !TODO needs unit tests
// lu decomposition, returns L value
template<Viewable V> requires (rows_v<V> == cols_v<V>)
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
template <Viewable V> requires (V::Rows == 1 and V::Cols == 1)
constexpr auto det(V const& v) {
	return v(0, 0);
}

// compute 2x2 determinant
template <Viewable V>
constexpr auto det(V const& v) requires (V::Rows == 2 and V::Cols == 2) {
	return v(0, 0)*v(1, 1) - v(0, 1)*v(1, 0);
}

// compute 3x3 determinant
template <Viewable V>
constexpr auto det(V const& v) requires (V::Rows == 3 and V::Cols == 3){
	return   (v(0, 0)*v(1, 1)*v(2, 2)
	        + v(0, 1)*v(1, 2)*v(2, 0)
	        + v(0, 2)*v(1, 0)*v(2, 1))
	       - (v(0, 2)*v(1, 1)*v(2, 0)
	        + v(0, 1)*v(1, 0)*v(2, 2)
	        + v(0, 0)*v(1, 2)*v(2, 1));
}

/*!\brief Determinant of Viewable
 * \ingroup freefunctions
 *
 * \param v Viewable must have same number of rows and columns
 * \return  Determinant of v
 *
 *  #### Example ####
 * \code
 *   auto a = SiLi::Matrix{{{ 1,  2},
 *                          {11, 12}}};
 *   auto v = det(a);
 *   std::cout << v << "\n"; // prints -10
 * \endcode
 */
template <Viewable V>
constexpr auto det(V const& v) requires (V::Rows == V::Cols and V::Rows > 3) {
	using T = typename V::value_t;
	auto retValue = T{1};

	auto L        = luDecomposition_L(v);

	for_constexpr<0, V::Rows>([&]<int i>() {
		retValue *= L(i, i);
	});

	using std::isfinite;
	if (not isfinite(retValue)) {
		retValue = T{0};
	}
	return retValue;
}

//element wise product
template <Viewable L, Viewable R> requires (L::Rows == R::Rows and L::Cols == R::Cols)
constexpr auto element_multi(L const& l, R const& r) {
	using U = decltype(std::declval<typename L::value_t>() * std::declval<typename R::value_t>());
	auto ret = Matrix<L::Rows, L::Cols, U>{};
	for_each_constexpr<L>([&]<auto row, auto col>() {
		ret(row, col) = l(row, col) * r(row, col);
	});
	return ret;
}

/*!\brief Cross product of two vectors.
 * \ingroup freefunctions
 *
 * \param l Viewable of size 3×1.
 * \param r Viewable of size 3×1.
 * \return  Matrix of size 3×1 with cross product of l and r.
 *
 *  #### Example ####
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
template <Viewable L, Viewable R> requires (L::Rows == R::Rows and L::Cols == R::Cols and L::Rows == 3 and L::Cols == 1)
constexpr auto cross(L const& l, R const& r) {
	using U = decltype(std::declval<typename L::value_t>() * std::declval<typename R::value_t>());
	auto ret = Matrix<3, 1, U>{};
	ret(0) = l(1) * r(2) - l(2) * r(1);
	ret(1) = l(2) * r(0) - l(0) * r(2);
	ret(2) = l(0) * r(1) - l(1) * r(0);
	return ret;
}

/*!\brief Sum of all elements.
 * \ingroup freefunctions
 *
 * \param v Viewable
 * \return  Sum off all elements of v.
 *
 *  #### Example ####
 * \code
 *   auto a = SiLi::Matrix{{{3, 4, 7},
 *                          {5, 6, 8}}};
 *   auto c = sum(a)
 *   std::cout << c << "\n"; // prints 33
 * \endcode
 */
template <Viewable V>
constexpr auto sum(V const& v) {
	auto acc = value_t<V>{};
	for_each_constexpr<V>([&]<auto row, auto col>() {
		acc += v(row, col);
	});
	return acc;
}

/*!\brief Sum of each row.
 * \ingroup freefunctions
 *
 * \param v Viewable
 * \return  Column vector with the sums of each row of v
 *
 *  #### Example ####
 * \code
 *   auto a = SiLi::Matrix{{{3, 4, 7},
 *                          {5, 6, 8}}};
 *   auto c = sum_rows(a)
 *   std::cout << c << "\n"; // prints {{14}
 *                                      {19}}
 * \endcode
 */
template <Viewable V>
constexpr auto sum_rows(V const& v) {
	using T = value_t<V>;
	auto c = Matrix<rows_v<V>, 1, T>{};
	for_constexpr<0, rows_v<V>>([&]<int row>() {
		at<row>(c) = sum(view_row<row>(v));
	});
	return c;
}

/*!\brief Sum of each column.
 * \ingroup freefunctions
 *
 * \param v Viewable
 * \return  Row vector with the sums of each column of v
 *
 *  #### Example ####
 * \code
 *   auto a = SiLi::Matrix{{{3, 4, 7},
 *                          {5, 6, 8}}};
 *   auto c = sum_cols(a)
 *   std::cout << c << "\n"; // prints {{8, 10, 15}}
 * \endcode
 */
template <Viewable V>
constexpr auto sum_cols(V const& v) {
	using T = value_t<V>;
	auto c = Matrix<1, cols_v<V>, T>{};
	for_constexpr<0, cols_v<V>>([&]<int col>() {
		at<col>(c) = sum(view_col<col>(v));
	});
	return c;
}



// inverse of 1x1
template <Viewable V> requires (V::Rows == V::Cols and V::Rows == 1)
constexpr auto inv(V v) -> std::tuple<typename V::value_t, V> {
	using T = typename V::value_t;
	using std::abs;
	auto d = det(v);
	if (abs(d) < 1.e-5) {
		return {v(0, 0), v};
	}

	v(0, 0) = T(1) / v(0, 0);
	return {d, v};
}

// inverse of 2x2
template <Viewable V> requires (V::Rows == V::Cols and V::Rows == 2)
constexpr auto inv(V v) -> std::tuple<typename V::value_t, V> {
	using T = typename V::value_t;
	using std::abs;
	auto d = det(v);
	if (abs(d) < 1.e-5) {
		return {v(0, 0), v};
	}
	auto c = T(1) / d;
	return {d, Matrix{{{ v(1, 1), -v(0, 1)},
	                   {-v(1, 0),  v(0, 0)}}}*c};
}

//inverse of 3x3
template <Viewable V> requires (V::Rows == V::Cols and V::Rows == 3)
constexpr auto inv(V v) -> std::tuple<typename V::value_t, V> {
	using T = typename V::value_t;
	using std::abs;
	auto d = det(v);
	if (abs(d) < 1.e-5) {
		return {v(0, 0), v};
	}
	auto c = T(1) / d;
	auto ret = V{};

	for_each_constexpr<V>([&]<int row, int col>() {
		auto tl = v((row+1)%3, (col+1)%3);
		auto br = v((row+2)%3, (col+2)%3);
		auto tr = v((row+1)%3, (col+2)%3);
		auto bl = v((row+2)%3, (col+1)%3);
		ret(col, row) = (tl*br-tr*bl) * c;
	});
	return {d, ret};
}

/*!\brief inverse of a matrix
 * \ingroup freefunctions
 *
 * \param  v Viewable must have same number of rows as number of cols.
 * \return   Returns tuple of determinant and inverse, if detereminant is zero the inverse is invalid.
 *
 *  #### Example ####
 * \code
 *   auto a = SiLi::Matrix{{{ 1.,  2.},
 *                          {11., 12.}}};
 *   auto [v, m] = inv(a);
 *   std::cout << v << "\n"; // prints -10
 *   std::cout << m << "\n"; // prints {{-1.2,  0.2},
 *                                      { 1.1, -0.1}}
 * \endcode
 */
template <Viewable V> requires (V::Rows == V::Cols and V::Rows > 3)
constexpr auto inv(V v) -> std::tuple<typename V::value_t, V> {
	using T = typename V::value_t;
	constexpr int N = V::Rows;
	auto det   = T{1};

	bool failed{false};
	for_constexpr<0, N>([&]<int p>() -> bool {
		auto pivot = v(p, p);
		det = det * pivot;
		if (fabs(pivot) < 1.e-5) {
			failed = true;
			return false;
		}
		view_col<p>(v) /= -pivot;
		for_each_constexpr<V>([&]<int i, int j>() {
			if (i == p or j == p) return;
			v(i, j) += v(p, j) * v(i, p);
		});
		view_row<p>(v) /= pivot;
		v(p, p) = T{1}/ pivot;
		return true;
	});
	if (failed) {
		return {0., v};
	}
	return {det, v};
}

/*!\brief view of the transposed
 * \ingroup freefunctions
 *
 * \param  v Viewable
 * \return   Viewable of tho transposed of v
 *
 *  #### Example ####
 * \code
 *   auto a = SiLi::Matrix{{{ 1,  2,  3},
 *                          {11, 12,  4}}};
 *   auto v = view_trans(a);
 *   std::cout << m << "\n"; // prints {{1, 11},
 *                                      {2, 12},
 *                                      {3, 13}}
 * \endcode
 */
template <typename V>
constexpr auto view_trans(V&& v) requires (Viewable<V>) {
	return View<cols_v<V>, rows_v<V>, stride_v<V>, value_t<V>, not transposed_v<V>>{v.data()};
}

/*!\brief transposed
 * \ingroup freefunctions
 *
 * \param  v Viewable
 * \return   Matrix of the transposed of v
 *
 *  #### Example ####
 * \code
 *   auto a = SiLi::Matrix{{{ 1,  2,  3},
 *                          {11, 12,  4}}};
 *   auto v = trans(a);
 *   std::cout << m << "\n"; // prints {{1, 11},
 *                                      {2, 12},
 *                                      {3, 13}}
 * \endcode
 */
template <Viewable V>
constexpr auto trans(V&& v) {
	return Matrix{view_trans(v)};
}

/*!\brief construct an identity matrix
 * \ingroup freefunctions
 *
 * \param Rows number of rows
 * \param Cols number of columns
 * \param T    type of each element
 * \return     matrix with ones on the diagonal
 *
 *  #### Example ####
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

/*!\brief construct a quadratic identity matrix
 * \ingroup freefunctions
 *
 * \param N size of the matrix
 * \param T type of each element
 * \erutrn matrix with ones on the diagonal
 *
 *  #### Example ####
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


/*!\brief computes the norm
 * \ingroup freefunctions
 *
 * \param v ViewableVector
 * \return the norm (also considered as the length).
 *
 *  #### Example ####
 * \code
 *   auto a = SiLi::Matrix{{{ 3.},
 *                          { 4.},
 *                          { 5.}}};
 *
 *   auto v = norm(a);
 *   std::cout << v << "\n"; // prints 7.0711
 * \endcode
 */
template <ViewableVector V>
constexpr auto norm(V const& v) {
	auto acc = value_t<V>{};
	for_each_constexpr<V>([&]<auto row, auto col>() {
		acc += v(row, col)*v(row, col);
	});
	using std::sqrt;
	return sqrt(acc);
}

/*!\brief dot product of two vectors
 * \ingroup freefunctions
 *
 * \param l ViewableVector
 * \param r ViewableVector
 * \return  the dot product (also called scalar product) of l and r. l and r must be vectors of the same length
 *
 *  #### Example ####
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
template <ViewableVector L, ViewableVector R> requires(length_v<L> == length_v<R>)
constexpr auto dot(L const& l, R const& r) {
	using U = decltype(value<L>() * value_t<R>());
	auto acc = U{};
	for_constexpr<0, length_v<L>>([&]<int I>() {
		acc += at<I>(l) * at<I>(r);
	});
	return acc;
}

/*!\brief compute outer product of two vectors
 * \ingroup freefunctions
 *
 * \param l ViewableVector
 * \param r ViewableVector
 * \return  Matrix with the outer product of l and r
 *
 *  #### Example ####
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
 *                                      { 8, 10, 12},
 *                                      {12, 15, 18}}
 * \endcode
 */
template <ViewableVector L, ViewableVector R>
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
