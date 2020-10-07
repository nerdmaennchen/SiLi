#pragma once

#include "concepts.h"

#include <cmath>
#include <optional>

namespace SiLi {

template <typename T>
constexpr auto value() -> std::add_rvalue_reference<value_t<T>>::type;

template <Viewable L, Viewable R> requires (L::Rows == R::Rows and L::Cols == R::Cols)
constexpr auto operator+(L const& l, R const& r) {
	using U = decltype(value<L>() + value_t<R>());
	auto ret = Matrix<L::Rows, L::Cols, U>{};
	for_each_constexpr<L>([&]<auto row, auto col>() {
		ret(row, col) = l(row, col) + r(row, col);
	});
	return ret;
}

template <Viewable L, Viewable R> requires (L::Rows == R::Rows and L::Cols == R::Cols)
constexpr auto operator+=(L& l, R const& r) -> auto& {
	for_each_constexpr<L>([&]<auto row, auto col>() {
		l(row, col) += r(row, col);
	});
	return l;
}

template <Viewable V>
constexpr auto operator-(V const& l) {
	using U = decltype(-value<V>());
	auto ret = Matrix<V::Rows, V::Cols, U>{};
	for_each_constexpr<V>([&]<auto row, auto col>() {
		ret(row, col) = -l(row, col);
	});
	return ret;
}

template <Viewable L, Viewable R> requires (L::Rows == R::Rows and L::Cols == R::Cols)
constexpr auto operator-(L const& l, R const& r) {
	using U = decltype(std::declval<typename L::value_t>() - std::declval<typename R::value_t>());
	auto ret = Matrix<L::Rows, L::Cols, U>{};
	for_each_constexpr<L>([&]<auto row, auto col>() {
		ret(row, col) = l(row, col) - r(row, col);
	});
	return ret;
}

template <Viewable L, Viewable R> requires (L::Rows == R::Rows and L::Cols == R::Cols)
constexpr auto operator-=(L& l, R const& r) -> auto& {
	for_each_constexpr<L>([&]<auto row, auto col>() {
		l(row, col) -= r(row, col);
	});
	return l;
}

template <Viewable L, Viewable R> requires (L::Cols == R::Rows)
constexpr auto operator*(L const& l, R const& r) {
	using U = decltype(std::declval<typename L::value_t>() * std::declval<typename R::value_t>());

	if constexpr (L::Rows == 1 and R::Cols == 1) {
		auto ret = U{};
		for_constexpr<0, L::Cols>([&]<int i>() {
			ret += l(i) * r(i);
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

template <Viewable L>
constexpr auto operator*(L&& l, value_t<L> const& s) {
	using U = decltype(value<L>() * value<L>());
	auto ret = Matrix<rows_v<L>, cols_v<L>, U>{};
	for_each_constexpr<L>([&]<auto row, auto col>() {
		ret(row, col) = l(row, col) * s;
	});

	return ret;
}
template <Viewable R>
constexpr auto operator*(value_t<R> const& s, R&& r) {
	return std::forward<R>(r) * s;
}


template <Viewable L>
constexpr auto operator*=(L& l, value_t<L> const& s) -> auto& {
	for_each_constexpr<L>([&]<auto row, auto col>() {
		l(row, col) *= s;
	});
	return l;
}


template <Viewable L, typename T>
constexpr auto operator/(L const& l, T const& s) {
	using U = decltype(std::declval<typename L::value_t>() * std::declval<T>());
	auto ret = Matrix<L::Rows, L::Cols, U>{};
	for_each_constexpr<L>([&]<auto row, auto col>() {
		ret(row, col) = l(row, col) / s;
	});

	return ret;
}

template <Viewable L, typename T>
constexpr auto operator/=(L&& l, T const& s) -> L& {
	for_each_constexpr<L>([&]<auto row, auto col>() {
		l(row, col) /= s;
	});
	return l;
}



template <Viewable L, Viewable R>
constexpr auto operator==(L const& l, R const& r) -> bool requires (L::Rows == R::Rows and L::Cols == R::Cols) {
	for (int row{0}; row < L::Rows; ++row) {
		for (int col{0}; col < L::Cols; ++col) {
			if (l(row, col) != r(row, col)) {
				return false;
			}
		}
	}
	return true;
}
template <Viewable L, Viewable R>
constexpr auto operator!=(L const& l, R const& r) -> bool requires (L::Rows == R::Rows and L::Cols == R::Cols) {
	for (int row{0}; row < L::Rows; ++row) {
		for (int col{0}; col < L::Cols; ++col) {
			if (l(row, col) != r(row, col)) {
				return true;
			}
		}
	}
	return false;
}


// diagonal view
template <Viewable V>
constexpr auto view_diag(V&& v) {
	using U = value_t<V>;
	using std::min;
	constexpr auto length = min(cols_v<V>, rows_v<V>);
	constexpr auto stride = stride_v<V>+1;
	return View<length, 1, stride, U, false>{v.data()};
}

template <Viewable V>
constexpr auto diag(V&& v) {
	return Matrix{view_diag(std::forward<V>(v))};
}

template <int start_row, int start_col, int end_row, int end_col, Viewable V>
constexpr auto view(V&& v) requires (end_row <= rows_v<V> and end_col <= cols_v<V>) {
	using U = value_t<V>;

	constexpr auto rows       = end_row - start_row;
	constexpr auto cols       = end_col - start_col;
	constexpr auto stride     = stride_v<V>;
	constexpr auto transposed = transposed_v<V>;

	return View<rows, cols, stride, U, transposed>{v.data() + start_col + start_row * stride};
}

template <int _row, Viewable V>
constexpr auto view_row(V&& v) {
	return view<_row, 0, _row+1, cols_v<V>>(std::forward<V>(v));
}

template <int _col, Viewable V>
constexpr auto view_col(V&& v) {
	return view<0, _col, rows_v<V>, _col+1>(std::forward<V>(v));
}

template <Viewable L, Viewable R>
constexpr auto join_rows(L const& l, R const& r) requires (L::Rows == R::Rows) {
	auto matrix = Matrix<L::Rows, L::Cols + R::Cols, typename L::value_t>{};
	view<0,       0, L::Rows,         L::Cols>(matrix) = l;
	view<0, L::Cols, R::Rows, L::Cols+R::Cols>(matrix) = r;
	return matrix;
}

template <Viewable L, Viewable R>
constexpr auto join_cols(L const& l, R const& r) requires (L::Cols == R::Cols) {
	auto matrix = Matrix<L::Rows+R::Rows, L::Cols, typename L::value_t>{};
	view<      0, 0,         L::Rows, L::Cols>(matrix) = l;
	view<L::Rows, 0, L::Rows+R::Cols, R::Cols>(matrix) = r;
	return matrix;
}

// !TODO needs unit tests
// lu decomposition, returns L value
template<Viewable V>
constexpr auto luDecomposition_L(V const& v) requires (rows_v<V> == cols_v<V>) {
	using T = std::decay_t<value_t<V>>;

	constexpr auto N = rows_v<V>;
	auto L = Matrix<N, N, T>{};

	for_constexpr<0, N>([&]<auto K>() {
		auto kRow = view_row<K>(L);
		auto kCol = view_col<K>(L);
		for_constexpr<K, N>([&]<auto I>() {
			L(K, I) = v(K, I) - T{kRow * view_col<I>(L)};
		});
		for_constexpr<K+1, N>([&]<auto I>() {
			L(I, K) = (v(I, K) - T{view_row<I>(L) * kCol}) / L(K, K);
		});
	});
	return L;
}


// compute 1x1 determinant
template <Viewable V>
constexpr auto det(V const& v) requires (V::Rows == 1 and V::Cols == 1) {
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

// compute determinante of any other matrix
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

// cross product for 3x1 and 3x1 matrices
template <Viewable L, Viewable R> requires (L::Rows == R::Rows and L::Cols == R::Cols and L::Rows == 3 and L::Cols == 1)
constexpr auto cross(L const& l, R const& r) {
	using U = decltype(std::declval<typename L::value_t>() * std::declval<typename R::value_t>());
	auto ret = Matrix<3, 1, U>{};
	ret(0) = l(1) * r(2) - l(2) * r(1);
	ret(1) = l(2) * r(0) - l(0) * r(2);
	ret(2) = l(0) * r(1) - l(1) * r(0);
	return ret;
}

// sum of all elements
template <Viewable V>
constexpr auto sum(V const& v) {
	auto acc = value_t<V>{};
	for_each_constexpr<V>([&]<auto row, auto col>() {
		acc += v(row, col);
	});
	return acc;
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

// transpose view
template <typename V>
constexpr auto view_trans(V&& v) requires (Viewable<V>) {
	return View<cols_v<V>, rows_v<V>, stride_v<V>, value_t<V>, not transposed_v<V>>{v.data()};
}

// transpose
template <Viewable V>
constexpr auto trans(V&& v) {
	return Matrix{view_trans(v)};
}

}
