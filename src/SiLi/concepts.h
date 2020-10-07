#pragma once

#include <type_traits>

namespace SiLi {

template<int, int, typename, typename...>
class Matrix;

template<int, int, int, typename, bool, typename...>
class View;


template <typename T>
struct is_matrix : std::false_type {};

template<int _rows, int _cols, typename T>
struct is_matrix<Matrix<_rows, _cols, T>> : std::true_type {};
template<int _rows, int _cols, typename T>
struct is_matrix<Matrix<_rows, _cols, T>&> : std::true_type {};
template<int _rows, int _cols, typename T>
struct is_matrix<Matrix<_rows, _cols, T>&&> : std::true_type {};
template<int _rows, int _cols, typename T>
struct is_matrix<Matrix<_rows, _cols, T> const&> : std::true_type {};


template <typename T>
struct is_view : std::false_type {};

template<int _rows, int _cols, int _stride, typename T, bool _transposed>
struct is_view<View<_rows, _cols, _stride, T, _transposed>> : std::true_type {};
template<int _rows, int _cols, int _stride, typename T, bool _transposed>
struct is_view<View<_rows, _cols, _stride, T, _transposed>&> : std::true_type {};
template<int _rows, int _cols, int _stride, typename T, bool _transposed>
struct is_view<View<_rows, _cols, _stride, T, _transposed>&&> : std::true_type {};
template<int _rows, int _cols, int _stride, typename T, bool _transposed>
struct is_view<View<_rows, _cols, _stride, T, _transposed> const&> : std::true_type {};


template <typename T>
concept Viewable = is_matrix<T>::value or is_view<T>::value;

// value_t for finding the underlying value
namespace detail {
template <typename T>
struct value_t;

template <Viewable V>
struct value_t<V> : std::type_identity<typename V::value_t> {};
template <Viewable V>
struct value_t<V&> : std::type_identity<typename V::value_t> {};
template <Viewable V>
struct value_t<V&&> : std::type_identity<typename V::value_t> {};
template <Viewable V>
struct value_t<V const&> : std::type_identity<const typename V::value_t> {};
}

template <typename T>
using value_t = detail::value_t<T>::type;

// find the underlying rows, cols and stride
namespace detail {
template <typename T> struct rows;
template <Viewable V> struct rows<V> : std::integral_constant<int, V::Rows> {};
template <Viewable V> struct rows<V&> : std::integral_constant<int, V::Rows> {};
template <Viewable V> struct rows<V&&> : std::integral_constant<int, V::Rows> {};
template <Viewable V> struct rows<V const&> : std::integral_constant<int, V::Rows> {};

template <typename T> struct cols;
template <Viewable V> struct cols<V> : std::integral_constant<int, V::Cols> {};
template <Viewable V> struct cols<V&> : std::integral_constant<int, V::Cols> {};
template <Viewable V> struct cols<V&&> : std::integral_constant<int, V::Cols> {};
template <Viewable V> struct cols<V const&> : std::integral_constant<int, V::Cols> {};

template <typename T> struct stride;
template <Viewable V> struct stride<V> : std::integral_constant<int, V::Stride> {};
template <Viewable V> struct stride<V&> : std::integral_constant<int, V::Stride> {};
template <Viewable V> struct stride<V&&> : std::integral_constant<int, V::Stride> {};
template <Viewable V> struct stride<V const&> : std::integral_constant<int, V::Stride> {};

template <typename T> struct transposed;
template <Viewable V> struct transposed<V> : std::integral_constant<bool, V::Transposed> {};
template <Viewable V> struct transposed<V&> : std::integral_constant<bool, V::Transposed> {};
template <Viewable V> struct transposed<V&&> : std::integral_constant<bool, V::Transposed> {};
template <Viewable V> struct transposed<V const&> : std::integral_constant<bool, V::Transposed> {};
}

template <typename T> constexpr int  rows_v       = detail::rows<T>::value;;
template <typename T> constexpr int  cols_v       = detail::cols<T>::value;;
template <typename T> constexpr int  stride_v     = detail::stride<T>::value;;
template <typename T> constexpr bool transposed_v = detail::transposed<T>::value;;



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
constexpr auto for_each_constexpr(L const& lambda) {
	for_constexpr<0, rows_v<V>>([&]<auto row>() {
		for_constexpr<0, cols_v<V>>([&]<auto col>() {
			lambda.template operator()<row, col>();
		});
	});
}



template <Viewable V>
constexpr auto get(V&& v, int row, int col) -> auto& {
	if constexpr (transposed_v<V>) {
		return v.data()[row + col * stride_v<V>];
	} else {
		return v.data()[col + row * stride_v<V>];
	}
}

template <Viewable V>
constexpr auto get(V&& v, int entry) -> auto& requires(rows_v<V> == 1 or cols_v<V> == 1) {
	if constexpr ((rows_v<V> == 1 and not transposed_v<V>) or (cols_v<V> == 1 and transposed_v<V>)) {
		return v.data()[entry];
	} else {
		return v.data()[entry * stride_v<V>];
	}
}



}
