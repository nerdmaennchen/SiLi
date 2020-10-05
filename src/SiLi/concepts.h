#pragma once

#include <type_traits>

namespace SiLi2 {

template<int, int, typename, typename...>
class Matrix;

template<int, int, int, typename, typename...>
class MatrixView;


template <typename T>
struct is_matrix : std::false_type {};

template<int _rows, int _cols, typename T>
struct is_matrix<Matrix<_rows, _cols, T>> : std::true_type {};
template<int _rows, int _cols, typename T>
struct is_matrix<const Matrix<_rows, _cols, T>> : std::true_type {};

template <typename T>
struct is_view : std::false_type {};

template<int _rows, int _cols, int _stride, typename T>
struct is_view<MatrixView<_rows, _cols, _stride, T>> : std::true_type {};
template<int _rows, int _cols, int _stride, typename T>
struct is_view<const MatrixView<_rows, _cols, _stride, T>> : std::true_type {};

template <typename T>
concept Viewable = is_matrix<T>::value or is_view<T>::value;

}
