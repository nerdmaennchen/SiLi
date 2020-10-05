#pragma once

#include "concepts.h"
#include <cmath>
#include <optional>

namespace SiLi2 {

// !TODO for_constexpr, is there a standard solution?
template <auto Iter, auto End, typename L>
constexpr void for_constexpr(L lambda) {
	if constexpr(Iter != End) {
		using R = decltype(lambda.template operator()<Iter>());
		if constexpr (std::is_same_v<bool, R>) {
			auto r = lambda.template operator()<Iter>();
			if (r) {
				for_constexpr<Iter+1, End>(lambda);
			}
		} else {
			lambda.template operator()<Iter>();
			for_constexpr<Iter+1, End>(lambda);
		}
	}
}

// !TODO for_each_constexpr, is there a standard solution?
template <Viewable V, typename L>
constexpr auto for_each(L const& lambda) {
	for_constexpr<0, V::Rows>([&]<auto row>() {
		for_constexpr<0, V::Cols>([&]<auto col>() {
			lambda.template operator()<row, col>();
		});
	});
}

template <Viewable L, Viewable R> requires (L::Rows == R::Rows and L::Cols == R::Cols)
constexpr auto operator+(L const& l, R const& r) {
	using U = decltype(std::declval<typename L::value_t>() + std::declval<typename R::value_t>());
	auto ret = Matrix<L::Rows, L::Cols, U>{};
	for_each<L>([&]<auto row, auto col>() {
		ret(row, col) = l(row, col) + r(row, col);
	});
	return ret;
}

template <Viewable L, Viewable R> requires (L::Rows == R::Rows and L::Cols == R::Cols)
constexpr auto operator+=(L& l, R const& r) {
	static_assert(not std::is_const_v<L>, "first parameter of operator += must be a non-const object");
	for_each<L>([&]<auto row, auto col>() {
		l(row, col) += r(row, col);
	});
	return l;
}

template <Viewable L>
constexpr auto operator-(L const& l) {
	using U = decltype(-std::declval<typename L::value_t>());
	auto ret = Matrix<L::Rows, L::Cols, U>{};
	for_each<L>([&]<auto row, auto col>() {
		ret(row, col) = -l(row, col);
	});
	return ret;
}

template <Viewable L, Viewable R> requires (L::Rows == R::Rows and L::Cols == R::Cols)
constexpr auto operator-(L const& l, R const& r) {
	using U = decltype(std::declval<typename L::value_t>() - std::declval<typename R::value_t>());
	auto ret = Matrix<L::Rows, L::Cols, U>{};
	for_each<L>([&]<auto row, auto col>() {
		ret(row, col) = l(row, col) - r(row, col);
	});
	return ret;
}

template <Viewable L, Viewable R> requires (L::Rows == R::Rows and L::Cols == R::Cols)
constexpr auto operator-=(L& l, R const& r) {
	static_assert(not std::is_const_v<L>, "first parameter of operator -= must be a non-const object");

	for_each<L>([&]<auto row, auto col>() {
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
		return ret;
	} else {
		auto ret = Matrix<L::Rows, R::Cols, U>{};
		for_constexpr<0, L::Rows>([&]<int row>() {
			for_constexpr<0, R::Cols>([&]<int col>() {
				ret(row, col) = view_row<row>(l) * view_col<col>(r);
			});
		});
		return ret;
	}
}

template <Viewable L, typename T>
constexpr auto operator*(L const& l, T const& s) {
	using U = decltype(std::declval<typename L::value_t>() * std::declval<T>());
	auto ret = Matrix<L::Rows, L::Cols, U>{};
	for_each<L>([&]<auto row, auto col>() {
		ret(row, col) = l(row, col) * s;
	});

	return ret;
}

template <Viewable L, typename T>
constexpr auto operator*=(L&& l, T const& s) -> L& {
	for_each<L>([&]<auto row, auto col>() {
		l(row, col) *= s;
	});
	return l;
}


template <Viewable L, typename T>
constexpr auto operator/(L const& l, T const& s) {
	using U = decltype(std::declval<typename L::value_t>() * std::declval<T>());
	auto ret = Matrix<L::Rows, L::Cols, U>{};
	for_each<L>([&]<auto row, auto col>() {
		ret(row, col) = l(row, col) / s;
	});

	return ret;
}

template <Viewable L, typename T>
constexpr auto operator/=(L&& l, T const& s) -> L& {
	for_each<L>([&]<auto row, auto col>() {
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


// read only diagonal view
template <Viewable L>
constexpr auto view_diag(L& l) {
	using T = typename L::value_t;
	using U = std::conditional_t<std::is_const_v<L>, const T, T>;
	using std::min;
	constexpr auto length = min(L::Cols, L::Rows);
	constexpr auto stride = L::Stride+1;
	return MatrixView<length, 1, stride, U>{l.data()};
}

template <Viewable L>
constexpr auto diag(L& l) {
	using U = typename L::value_t;
	using std::min;
	constexpr auto length = min(L::Cols, L::Rows);
	auto mat = Matrix<length, 1, U>{};
	for_constexpr<0, length>([&]<int row>() {
		mat(row, 0) = l(row, row);
	});
	return mat;
}

template <int start_row, int start_col, int end_row, int end_col, Viewable L>
constexpr auto view(L& l) {
	using T = typename L::value_t;
	using U = std::conditional_t<std::is_const_v<L>, const T, T>;

	constexpr int rows = end_row - start_row;
	constexpr int cols = end_col - start_col;
	static_assert(rows <= L::Rows, "new View must be smaller than old view");
	static_assert(cols <= L::Cols, "new View must be smaller than old view");

	constexpr int stride = L::Stride;

	return MatrixView<rows, cols, stride, U>{l.data() + start_col + start_row * stride};
}

template <int _row, Viewable V>
constexpr auto view_row(V& v) {
	return view<_row, 0, _row+1, V::Cols>(v);
}

template <int _col, Viewable V>
constexpr auto view_col(V& v) {
	return view<0, _col, V::Rows, _col+1>(v);
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
constexpr auto luDecomposition_L(V& v) requires (V::Rows == V::Cols) {
	constexpr auto N = V::Rows;
	auto L = Matrix<N, N, typename V::value_t>{};

	for_constexpr<0, N>([&]<auto K>() {
		auto kRow = view_row<K>(L);
		auto kCol = view_col<K>(L);
		for_constexpr<K, N>([&]<auto I>() {
			L(K, I) = v(K, I) - kRow * view_col<I>(L);
		});
		for_constexpr<K+1, N>([&]<auto I>() {
			L(I, K) = (v(I, K) - view_row<I>(L) * kCol) / L(K, K);
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
	for_each<L>([&]<auto row, auto col>() {
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
	using T = typename V::value_t;
	auto acc = T{};
	for_each<V>([&]<auto row, auto col>() {
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

	for_each<V>([&]<int row, int col>() {
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
		for_each<V>([&]<int i, int j>() {
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


}
