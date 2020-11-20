#pragma once

#include <type_traits>
#include <utility>

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

namespace _concept {
/*! Concept of a _concept::Matrix.
 * \shortexample _concept::Matrix

 *
 * Abstract concept of a matrix. This can be either a Matrix or a View.
 * Both can be used everywhere the _concept::Matrix concept is needed.
 */
template <typename T>
concept Matrix = is_matrix<T>::value or is_view<T>::value;

}

// value_t for finding the underlying value
namespace detail {
template <typename T>
struct value_t;

template <_concept::Matrix V>
struct value_t<V> : std::type_identity<typename V::value_t> {};
template <_concept::Matrix V>
struct value_t<V&> : std::type_identity<typename V::value_t> {};
template <_concept::Matrix V>
struct value_t<V&&> : std::type_identity<typename V::value_t> {};
template <_concept::Matrix V>
struct value_t<V const&> : std::type_identity<const typename V::value_t> {};
}

template <typename T>
using value_t = typename detail::value_t<T>::type;

// find the underlying rows, cols and stride
namespace detail {
template <typename T> struct rows;
template <_concept::Matrix V> struct rows<V> : std::integral_constant<int, V::Rows> {};
template <_concept::Matrix V> struct rows<V&> : std::integral_constant<int, V::Rows> {};
template <_concept::Matrix V> struct rows<V&&> : std::integral_constant<int, V::Rows> {};
template <_concept::Matrix V> struct rows<V const&> : std::integral_constant<int, V::Rows> {};

template <typename T> struct cols;
template <_concept::Matrix V> struct cols<V> : std::integral_constant<int, V::Cols> {};
template <_concept::Matrix V> struct cols<V&> : std::integral_constant<int, V::Cols> {};
template <_concept::Matrix V> struct cols<V&&> : std::integral_constant<int, V::Cols> {};
template <_concept::Matrix V> struct cols<V const&> : std::integral_constant<int, V::Cols> {};

template <typename T> struct stride;
template <_concept::Matrix V> struct stride<V> : std::integral_constant<int, V::Stride> {};
template <_concept::Matrix V> struct stride<V&> : std::integral_constant<int, V::Stride> {};
template <_concept::Matrix V> struct stride<V&&> : std::integral_constant<int, V::Stride> {};
template <_concept::Matrix V> struct stride<V const&> : std::integral_constant<int, V::Stride> {};

template <typename T> struct transposed;
template <_concept::Matrix V> struct transposed<V> : std::integral_constant<bool, V::Transposed> {};
template <_concept::Matrix V> struct transposed<V&> : std::integral_constant<bool, V::Transposed> {};
template <_concept::Matrix V> struct transposed<V&&> : std::integral_constant<bool, V::Transposed> {};
template <_concept::Matrix V> struct transposed<V const&> : std::integral_constant<bool, V::Transposed> {};
}

template <typename T> constexpr int  rows_v       = detail::rows<T>::value;;
template <typename T> constexpr int  cols_v       = detail::cols<T>::value;;
template <typename T> constexpr int  stride_v     = detail::stride<T>::value;;
template <typename T> constexpr bool transposed_v = detail::transposed<T>::value;;

template <_concept::Matrix T>
constexpr bool is_vector_v = (rows_v<T> == 1 or cols_v<T> == 1);

namespace _concept {
/*! Concept of a _concept::Vector.
 * \shortexample _concept::Vector
 *
 * Abstract concept of a vector. This can be either a Matrix or a View where either columns or rows is 1.
 */
template <typename T>
concept Vector = is_vector_v<T>;
}

template <_concept::Vector T>
constexpr int length_v = ((detail::rows<T>::value==1)?detail::cols<T>::value:detail::rows<T>::value);


namespace details {
// !TODO for_constexpr, is there a standard solution?
template <auto Iter, auto End, typename L> requires (std::is_same_v<void, decltype(std::declval<L>().template operator()<Iter>())>)
constexpr void for_constexpr_void(L lambda) {
	if constexpr(Iter != End) {
		lambda.template operator()<Iter>();
		for_constexpr<Iter+1, End>(lambda);
	}
}
template <auto Iter, auto End, typename L> requires (std::is_same_v<bool, decltype(std::declval<L>().template operator()<Iter>())>)
constexpr bool for_constexpr_bool(L lambda) {
	if constexpr(Iter != End) {
		auto r = lambda.template operator()<Iter>();
		if (not r) {
			return false;
		}
		return for_constexpr<Iter+1, End>(lambda);
	}
	return true;
}
}

template <auto Begin, auto End, typename L>
constexpr bool for_constexpr(L&& lambda) {
	if constexpr (Begin == End) {
		return true;
	} else {
		using R = decltype(std::declval<L>().template operator()<Begin>());
		if constexpr (std::is_same_v<R, void>) {
			details::for_constexpr_void<Begin, End>(std::forward<L>(lambda));
			return true;
		} else {
			return details::for_constexpr_bool<Begin, End>(std::forward<L>(lambda));
		}
	}
}



namespace details {
// !TODO for_each_constexpr, is there a standard solution?
template <_concept::Matrix V, typename L>
constexpr void for_each_constexpr_void(L&& lambda) {
	for_constexpr<0, rows_v<V>>([&]<auto row>() constexpr {
		for_constexpr<0, cols_v<V>>([&]<auto col>() constexpr {
			lambda.template operator()<row, col>();
		});
	});
}
template <_concept::Matrix V, typename L>
constexpr bool for_each_constexpr_bool(L&& lambda) {
	return for_constexpr<0, rows_v<V>>([&]<auto row>() constexpr {
		return for_constexpr<0, cols_v<V>>([&]<auto col>() constexpr {
			return lambda.template operator()<row, col>();
		});
	});
}

}
template <_concept::Matrix V, typename L>
constexpr bool for_each_constexpr(L&& lambda) {
	if constexpr (rows_v<V> == 0 or cols_v<V> == 0) {
		return true;
	} else {
		using R = decltype(std::declval<L>().template operator()<0, 0>());
		if constexpr (std::is_same_v<R, void>) {
			details::for_each_constexpr_void<V>(std::forward<L>(lambda));
			return true;
		} else {
			return details::for_each_constexpr_bool<V>(std::forward<L>(lambda));
		}
	}
}




template <_concept::Matrix V>
constexpr auto get(V&& v, int row, int col) -> auto& {
	if constexpr (transposed_v<V>) {
		return v.data()[row + col * stride_v<V>];
	} else {
		return v.data()[col + row * stride_v<V>];
	}
}

template <_concept::Matrix V>
constexpr auto get(V&& v, int entry) -> auto& requires(rows_v<V> == 1 or cols_v<V> == 1) {
	if constexpr ((rows_v<V> == 1 and not transposed_v<V>) or (cols_v<V> == 1 and transposed_v<V>)) {
		return v.data()[entry];
	} else {
		return v.data()[entry * stride_v<V>];
	}
}

template <int row, int col, _concept::Matrix M>
constexpr auto get(M&& m) -> auto& {
	if constexpr (transposed_v<M>) {
		return m.data()[row + col * stride_v<M>];
	} else {
		return m.data()[col + row * stride_v<M>];
	}
}

template <int entry, _concept::Matrix M>
constexpr auto get(M&& m) -> auto& requires(rows_v<M> == 1 or cols_v<M> == 1) {
	if constexpr ((rows_v<M> == 1 and not transposed_v<M>) or (cols_v<M> == 1 and transposed_v<M>)) {
		return m.data()[entry];
	} else {
		return m.data()[entry * stride_v<M>];
	}
}




}
