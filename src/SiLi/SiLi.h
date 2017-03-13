#pragma once

#include <algorithm>
#include <array>
#include <complex>
#include <iterator>
#include <type_traits>
#include <iostream>

namespace SiLi {

using DefaultType = float;

// This exception is thrown by svd() if the max iteration count is reached
struct MaxIteration {};

template<int _rowStride, bool _transposed = false, int _offset=0>
struct Properties {
	static constexpr int  rowStride  {_rowStride};
	static constexpr bool transposed {_transposed};
	static constexpr int  offset     {_offset};

	using Transposed = Properties<_rowStride, not _transposed, _offset>;
	using Diag       = Properties<_rowStride, false, 1>;
};

template<int, int, typename T = DefaultType, typename... Args>
class Matrix;

template<int trows, int tcols, typename... Args>
class MatrixView;

/** Matrix view iterator, to iterate over all elements
 * of a matrix
 */
template<typename MatrixView>
class MatrixIterator : public std::forward_iterator_tag {
private:
	int row{0};
	int col{0};
	MatrixView* mat;
public:
	MatrixIterator(MatrixView* view)
		: mat(view)
	{
		if (not view) {
			col = 0;
			row = mat->num_rows();
		}
	}

	MatrixIterator(MatrixIterator const& iter)
		: row {iter.row}
		, col {iter.col}
		, mat {iter.mat}
	{}
	auto operator=(MatrixIterator const& iter) -> MatrixIterator& {
		row = iter.row;
		col = iter.col;
		mat = iter.mat;
		return *this;
	}

	auto operator++() -> MatrixIterator& {
		col += 1;
		if (col >= mat->num_cols()) {
			col = 0;
			row += 1;
		}
		return *this;
	}
	auto operator++(int) -> MatrixIterator {
		MatrixIterator iter = *this;
		operator++();
		return iter;
	}

	bool operator==(MatrixIterator const& iter) const {
		return iter.row == row and iter.col == col;
	}
	bool operator!=(MatrixIterator const& iter) const {
		return iter.row != row or iter.col != col;
	}
	auto operator*() const -> decltype((*mat)(row, col)) {
		return (*mat)(row, col);
	}
	auto operator*() -> decltype((*mat)(row, col))&  {
		return (*mat)(row, col);
	}
};


/** Some helper class to force look up at instanciation time
 */
template <int T>
struct Number {
	static constexpr int _number {T};
};

/** SVD Results;
 */
template <int rows, int cols, typename T>
struct SVD {
	Matrix<rows, cols, T> U;
	Matrix<cols, 1, T>    S;
	Matrix<cols, cols, T> V;
};


/** Const MatrixView with all methods that are allowed to be operated on a const matrix view
 */
template<int trows, int tcols, typename prop, typename T>
class MatrixView<trows, tcols, prop, T const> {
	static_assert(trows >= 0, "Row number must be positive");
	static_assert(tcols >= 0, "Coloumn number must be positive");

public:
	using Type = T;
	using Props = prop;
protected:
	T const* const cBasePtr;

public:
	explicit MatrixView(T const* base) : cBasePtr(base) {}

	MatrixView(MatrixView const& rhs) : cBasePtr(rhs.cBasePtr) {}

	template<typename T2>
	MatrixView(MatrixView<trows, tcols, prop, T2> const& rhs) : cBasePtr(rhs.cBasePtr) {}

	constexpr int num_rows() const {return trows;}
	constexpr int num_cols() const {return tcols;}

	// read only element access not transposed
	template<bool t = prop::transposed, typename std::enable_if<not t>::type* = nullptr>
	auto operator()(int row, int col) const -> T const& {
		return *(cBasePtr + (row * prop::rowStride) + col + prop::offset*row);
	}

	// read only element access transposed
	template<bool t = prop::transposed, typename std::enable_if<t>::type* = nullptr>
	auto operator()(int row, int col) const -> T const& {
		return *(cBasePtr + (col * prop::rowStride) + row + prop::offset*col);
	}

	// read only element access not transposed
	template<int _rows = Number<trows>::_number, typename std::enable_if<_rows == 1>::type* = nullptr>
	auto operator()(int col) const -> T const& {
		return operator()(0, col);
	}

	// read only element access not transposed
	template<int _cols = Number<tcols>::_number, typename std::enable_if<_cols == 1 and trows != 1>::type* = nullptr>
	auto operator()(int row) const -> T const& {
		return operator()(row, 0);
	}

	template<int _cols = tcols, typename std::enable_if<trows == _cols and trows == 1>::type* = nullptr>
	operator T() const {
		return (*this)(0, 0);
	};



	// matrix access
	template<int subRows, int subCols>
	auto mat(int startR, int startC) const -> Matrix<subRows, subCols, T> {
		return view<subRows, subCols>(startR, startC);
	}

	// read only view access
	template<int subRows, int subCols>
	auto view(int startR, int startC) const -> MatrixView<subRows, subCols, prop, T const> {
		static_assert(subRows <= trows, "rows must be smaller or equal to the current view");
		static_assert(subCols <= tcols, "cols must be smaller or equal to the current view");

		return MatrixView<subRows, subCols, prop, T const> {&((*this)(startR, startC))};
	}

	auto view_row(int startR) const -> MatrixView<1, tcols, prop, T const> {
		return MatrixView<1, tcols, prop, T const> {&((*this)(startR, 0))};
	}

	auto view_col(int startC) const -> MatrixView<trows, 1, prop, T const> {
		return MatrixView<trows, 1, prop, T const> {&((*this)(0, startC))};
	}

	// array with view on each row
	auto rows() const -> std::array<MatrixView<1, tcols, prop, T const>, trows>;

	// array with view on each col
	auto cols() const -> std::array<MatrixView<trows, 1, prop, T const>, tcols>;


	// read only iterator
	auto begin() const -> MatrixIterator<MatrixView<trows, tcols, prop, T const> const> {
		return {this};
	}
	auto end() const -> MatrixIterator<MatrixView<trows, tcols, prop, T const> const> {
		return {nullptr};
	}

	template <int cols, typename oProp>
	auto join_rows(MatrixView<trows, cols, oProp, T const> const& _view) const -> Matrix<trows, tcols+cols, T> {
		Matrix<trows, tcols+cols, T> retValue;
		retValue.view<trows, tcols>(0, 0)    = *this;
		retValue.view<trows, cols>(0, tcols) = _view;
		return retValue;
	}

	template <int rows, typename oProp>
	auto join_cols(MatrixView<rows, tcols, oProp, T const> const& _view) const -> Matrix<trows+rows, tcols, T> {
		Matrix<trows+rows, tcols, T> retValue;
		retValue.view<trows, tcols>(0, 0)    = *this;
		retValue.view<rows,  tcols>(trows, 0) = _view;
		return retValue;
	}



	// negation
	auto operator-() const -> Matrix<trows, tcols, T> {
		Matrix<trows, tcols, T> ret;
		for (int r(0); r < trows; ++r) {
			for (int c(0); c < tcols; ++c) {
				ret(r, c) = -(*this)(r, c);
			}
		}
		return ret;
	}

	// compute inverse
	auto inv() const -> Matrix<trows, tcols, T>;

	// determinant
	auto det() const -> T;

	// read only transposed
	auto t() const -> MatrixView<tcols, trows, typename prop::Transposed, T const> {
		return MatrixView<tcols, trows, typename prop::Transposed, T const> {cBasePtr};
	}

	// read only diagonal view
	auto diag() const -> MatrixView<(tcols < trows)?tcols:trows, 1, typename prop::Diag, T const> {
		return MatrixView<(tcols < trows)?tcols:trows, 1, typename prop::Diag, T const> {cBasePtr};
	}

	/* squared frobenius norm */
	auto normSqr() const -> T {
		T ret{0.};
		for (auto const& d : (*this)) {
			ret += d*d;
		}
		return ret;
	}

	/* frobenius norm */
	auto norm() const -> T {
		return std::sqrt(normSqr());
	}

	auto normalized() const -> Matrix<trows, tcols, T> {
		return (*this) * (1./norm());
	}

	// compute svd so that *this = svd.U * make_diag(svd.S) * svd.V.t()
	// Will throw SiLi::MaxIteration if maximum of iteration is reached
	auto svd() const -> SVD<trows, tcols, T>;

	// compute abs
	auto abs() const -> Matrix<trows, tcols, T> {
		return abs(*this);
	}

	// compute sum
	auto sum() const -> T {
		return sum(*this);
	}

	auto isfinite() const -> bool {
		using std::isfinite;
		for (auto const& d : (*this)) {
			if (not isfinite(d)) {
				return false;
			}
		}
		return true;
	}
};

template<int trows, int tcols, typename prop, typename T>
class MatrixView<trows, tcols, prop, T> : public MatrixView<trows, tcols, prop, T const> {
private:
	using CView = MatrixView<trows, tcols, prop, T const>;
protected:
	T* const basePtr;

public:
	explicit MatrixView(T* base) : CView(base), basePtr(base) {}

	MatrixView(MatrixView& rhs) : CView(rhs), basePtr(rhs.basePtr) {}
	MatrixView(MatrixView&& rhs) : CView(rhs), basePtr(rhs.basePtr) {}

	// pass through
	using MatrixView<trows, tcols, prop, T const>::operator();

	// value element access
	template<bool t = prop::transposed, typename std::enable_if<not t>::type* = nullptr>
	auto operator()(int row, int col) -> T& {
		return *(basePtr + (row * prop::rowStride) + col + row * prop::offset);
	}

	// value element access
	template<bool t = prop::transposed, typename std::enable_if<t>::type* = nullptr>
	auto operator()(int row, int col) -> T& {
		return *(basePtr + (col * prop::rowStride) + row + col * prop::offset);
	}

	// access if matrix is one dimensional
	template<int _rows = Number<trows>::_number, typename std::enable_if<_rows == 1>::type* = nullptr>
	auto operator()(int col) -> T& {
		return operator()(0, col);
	}

	// access if matrix is one dimensional
	template<int _cols = Number<tcols>::_number, typename std::enable_if<_cols == 1 and trows != 1>::type* = nullptr>
	auto operator()(int row) -> T& {
		return operator()(row, 0);
	}

	// view access
	using CView::view;
	template<int subRows, int subCols>
	auto view(int startR, int startC) -> MatrixView<subRows, subCols, prop, T> {
		static_assert(subRows <= trows, "rows must be smaller or equal to the current view");
		static_assert(subCols <= tcols, "cols must be smaller or equal to the current view");

		return MatrixView<subRows, subCols, prop, T> {&((*this)(startR, startC))};
	}

	auto view_row(int startR) -> MatrixView<1, tcols, prop, T> {
		return MatrixView<1, tcols, prop, T> {&((*this)(startR, 0))};
	}

	auto view_col(int startC) -> MatrixView<trows, 1, prop, T> {
		return MatrixView<trows, 1, prop, T> {&((*this)(0, startC))};
	}

	// array with view on each row
	auto rows() -> std::array<MatrixView<1, tcols, prop, T>, trows>;

	// array with view on each col
	auto cols() -> std::array<MatrixView<trows, 1, prop, T>, tcols>;



	auto operator*=(T const& rhs) -> MatrixView& {
		for (auto& x : *this) {
			x *= rhs;
		}
		return *this;
	}
	auto operator+=(T const& rhs) -> MatrixView& {
		for (auto& x : *this) {
			x += rhs;
		}
		return *this;
	}

	template <typename oProp>
	auto operator+=(MatrixView<trows, tcols, oProp, T const> const& rhs) -> MatrixView& {
		for (int row(0); row < trows; ++row) {
			for (int col(0); col < tcols; ++col) {
				(*this)(row, col) += rhs(row, col);
			}
		}
		return (*this);
	}

	auto operator-=(T const& rhs) -> MatrixView& {
		for (auto& x : *this) {
			x -= rhs;
		}
		return *this;
	}

	template <typename oProp>
	auto operator-=(MatrixView<trows, tcols, oProp, T const> const& rhs) -> MatrixView& {
		for (int row(0); row < trows; ++row) {
			for (int col(0); col < tcols; ++col) {
				(*this)(row, col) -= rhs(row, col);
			}
		}
		return *this;
	}

	// element wise product
	template<typename oProp>
	auto operator&=(MatrixView<trows, tcols, oProp, T const> const& _view) -> MatrixView& {
		for (int row(0); row < trows; ++row) {
			for (int col(0); col < tcols; ++col) {
				(*this)(row, col) *= _view(row, col);
			}
		}
		return *this;
	}

	// element wise product
	template<typename oProp>
	auto operator&(MatrixView<trows, tcols, oProp, T const> const& _view) -> Matrix<trows, tcols, T> {
		Matrix<trows, tcols, T> retMat;
		for (int row(0); row < trows; ++row) {
			for (int col(0); col < tcols; ++col) {
				retMat(row, col) = (*this)(row, col) * _view(row, col);
			}
		}
		return retMat;
	}



	auto operator=(T const& rhs) -> MatrixView& {
		for (int r(0); r < trows; ++r) {
			for (int c(0); c < tcols; ++c) {
				(*this)(r, c) = rhs;
			}
		}
		return *this;
	}

	auto operator=(MatrixView const& rhs) const -> MatrixView& = delete;

	auto operator=(MatrixView const& rhs) -> MatrixView& {
		return operator=<prop>(rhs);
	}

	template<typename oProp>
	auto operator=(MatrixView<trows, tcols, oProp, T const> const& rhs) -> MatrixView& {
		for (int r(0); r < trows; ++r) {
			for (int c(0); c < tcols; ++c) {
				(*this)(r, c) = rhs(r, c);
			}
		}
		return *this;
	}

	auto begin() -> MatrixIterator<MatrixView<trows, tcols, prop, T>> {
		return {this};
	}
	auto end() -> MatrixIterator<MatrixView<trows, tcols, prop, T>> {
		return {nullptr};
	}

	// transposed view
	using CView::t;
	auto t() -> MatrixView<tcols, trows, typename prop::Transposed, T> {
		return MatrixView<tcols, trows, typename prop::Transposed, T> {basePtr};
	}

	// diagonal view
	using CView::diag;
	auto diag() -> MatrixView<(tcols < trows)?tcols:trows, 1, typename prop::Diag, T> {
		return MatrixView<(tcols < trows)?tcols:trows, 1, typename prop::Diag, T> {basePtr};
	}
};

template<int rows, int cols, typename T>
class Matrix<rows, cols, T> : public MatrixView<rows, cols, class Properties<cols>, T> {
	T vals[rows][cols];
public:
	using View  = MatrixView<rows, cols, Properties<cols>, T>;
	using CView = MatrixView<rows, cols, Properties<cols>, T const>;


	Matrix() : View(&(vals[0][0])) {}

	Matrix(T initVal) : View(&(vals[0][0]))  {
		((View*)(this))->operator=(initVal);
	}

	template<typename T2>
	Matrix(T2 const (&list)[rows*cols]) : View(&(vals[0][0])) {
		T2 const* iter = &(list[0]);
		for (int r(0); r < rows; ++r) {
			for (int c(0); c < cols; ++c) {
				(*this)(r, c) = T(*iter++);
			}
		}
	}

	Matrix(T const (&list)[rows][cols]) : View(&(vals[0][0])) {
		for (int r(0); r < rows; ++r) {
			for (int c(0); c < cols; ++c) {
				(*this)(r, c) = T(list[r][c]);
			}
		}
	}

	Matrix(Matrix&& other) : View(&(vals[0][0])) {
		for (int r(0); r < rows; ++r) {
			for (int c(0); c < cols; ++c) {
				(*this)(r, c) = other(r, c);
			}
		}
	}


	Matrix(Matrix const& other) : View(&(vals[0][0])) {
		for (int r(0); r < rows; ++r) {
			for (int c(0); c < cols; ++c) {
				(*this)(r, c) = other(r, c);
			}
		}
	}

	template<typename Props>
	Matrix(MatrixView<rows, cols, Props, T const> const& other) : View(&(vals[0][0])) {
		for (int r(0); r < rows; ++r) {
			for (int c(0); c < cols; ++c) {
				(*this)(r, c) = other(r, c);
			}
		}
	}

	template<typename Props>
	auto operator=(MatrixView<rows, cols, Props, T const> const& other) & -> Matrix& {
		View::operator=(other);
		return *this;
	}

	template<typename Props>
	auto operator=(MatrixView<rows, cols, Props, T const> const& other) && -> Matrix& = delete;

	auto operator=(Matrix const& other) -> Matrix& {
		for (int r(0); r < rows; ++r) {
			for (int c(0); c < cols; ++c) {
				(*this)(r, c) = other(r, c);
			}
		}

		return *this;
	}
};

// create matrix
template<int rows, int cols, typename T = DefaultType>
auto make_mat(T const (&values)[rows][cols]) -> Matrix<rows, cols, T> {
	return {values};
}

// create identity matrix
template<int rows, int cols, typename T = DefaultType>
auto make_eye() -> Matrix<rows, cols, T> {
	Matrix<rows, cols, T> retVal(0);
	retVal.diag() = 1.;
	return retVal;
}

// create matrix from vector
template<int rows, typename T, typename Prop>
auto make_diag(MatrixView<rows, 1, Prop, T const> const& _view) -> Matrix<rows, rows, T> {
	Matrix<rows, rows, T> retVal(0);
	retVal.diag() = _view;
	return retVal;
}
template<int rows, int cols, typename T, typename Prop>
auto make_diag(MatrixView<rows, 1, Prop, T const> const& _view) -> Matrix<rows, cols, T> {
	static_assert(cols >= rows, "cannot build a smaller diagonal matrix than the input vector size");
	Matrix<rows, cols, T> retVal(0);
	retVal.diag() = _view;
	return retVal;
}
template<int rows, int cols, typename T, typename Prop>
auto make_diag(MatrixView<cols, 1, Prop, T const> const& _view) -> Matrix<rows, cols, T> {
	static_assert(rows >= cols, "cannot build a smaller diagonal matrix than the input vector size");
	Matrix<rows, cols, T> retVal(0);
	retVal.diag() = _view;
	return retVal;
}

// only works for not overlaping views
template<int rows, int cols, typename P1, typename P2, typename T>
void swap(MatrixView<rows, cols, P1, T>& lhs, MatrixView<rows, cols, P2, T>& rhs) {
	for(int i = 0; i < rows; ++i) {
		for(int j = 0; j < cols; ++j) {
			std::swap(lhs(i, j), rhs(i, j));
		}
	}
}


// lu decomposition, returns L value
template <int rows, int cols, typename Props, typename T>
auto luDecomposition_L(MatrixView<rows, cols, Props, T const> const& _mat) -> Matrix<rows, cols, T> {
	Matrix<rows, cols, T> L{0};

	for (int k=0; k<rows; ++k) {
		auto kRow = L.view_row(k);
		auto kCol = L.view_col(k);
		for (int i=k;i<rows;++i) {
			L(k, i) = _mat(k, i) - kRow * L.view_col(i);
		}
		for (int i=k+1; i<rows;++i) {
			L(i, k) = (_mat(i, k) - T(L.view_row(i) * kCol)) / L(k, k);
		}
	}
	return L;
}

// compute minor matrix
template <int rows, int cols, typename Props, typename T>
auto minorMat(MatrixView<rows, cols, Props, T const> const& _view, int _row, int _col) -> Matrix<rows-1, cols-1, T> {

	Matrix<rows-1, cols-1, T> retMat;

	for(int i = 0; i < _row; ++i) {
		for(int j = 0; j < _col; ++j) {
			retMat(i, j) = _view(i, j);
		}

		for(int j = _col+1; j < cols; ++j) {
			retMat(i, j-1) = _view(i, j);
		}
	}

	for(int i = _row+1; i < rows; ++i) {
		for(int j = 0; j < cols; ++j) {
			retMat(i-1, j) = _view(i, j);
		}

		for(int j = _col+1; j < cols; ++j) {
			retMat(i-1, j-1) = _view(i, j);
		}
	}
	return retMat;
}

// compute adjugated matrix
template <int rows, int cols, typename Props, typename T>
auto adjugateMat(MatrixView<rows, cols, Props, T const> const& _view) -> Matrix<rows, cols, T> {
	Matrix<rows, cols, T> retMat;
	if (cols == 1) {
		return {1};
	}

	for(int i = 0; i < rows; ++i) {
		for(int j = 0; j < cols; ++j) {
			T sign = (i%2*2-1) * (j%2*2-1);
			auto det = minorMat(_view, i, j).det();
			retMat(i, j) = sign * det;
		}
	}
	return retMat;
}

/**
  * create all row view
  */

template<int... I>
struct index {
	template<int n>
	using append = index<I..., n>;
};

template<int N>
struct make_index {
	typedef typename make_index<N - 1>::type::template append<N - 1> type;
};

template<>
struct make_index<0> {
	typedef index<> type;
};

template<int N>
using indexer = typename make_index<N>::type;

template<int rows, int cols, typename prop, typename T, int ...nbrs>
auto createRowViews(MatrixView<rows, cols, prop, T const> const& view, index<nbrs...>) -> std::array<MatrixView<1, cols, prop, T const>, rows> {
	auto retArray = std::array<MatrixView<1, cols, prop, T const>, rows> {
		view.view_row(nbrs)...
	};
	return retArray;
}

template<int rows, int cols, typename prop, typename T, int ...nbrs>
auto createRowViews(MatrixView<rows, cols, prop, T>& view, index<nbrs...>) -> std::array<MatrixView<1, cols, prop, T>, rows> {
	auto retArray = std::array<MatrixView<1, cols, prop, T>, rows> {
		view.view_row(nbrs)...
	};
	return retArray;
}

template<int trows, int tcols, typename prop, typename T>
auto MatrixView<trows, tcols, prop, T const>::rows() const -> std::array<MatrixView<1, tcols, prop, T const>, trows> {
	return createRowViews(*this, indexer<trows>());
}

template<int trows, int tcols, typename prop, typename T>
auto MatrixView<trows, tcols, prop, T>::rows() -> std::array<MatrixView<1, tcols, prop, T>, trows> {
	return createRowViews(*this, indexer<trows>());
}

template<int rows, int cols, typename prop, typename T, int ...nbrs>
auto createColViews(MatrixView<rows, cols, prop, T const> const& view, index<nbrs...>) -> std::array<MatrixView<rows, 1, prop, T const>, cols> {
	auto retArray = std::array<MatrixView<rows, 1, prop, T const>, cols> {
		view.view_col(nbrs)...
	};
	return retArray;
}

template<int rows, int cols, typename prop, typename T, int ...nbrs>
auto createColViews(MatrixView<rows, cols, prop, T>& view, index<nbrs...>) -> std::array<MatrixView<rows, 1, prop, T>, cols> {
	auto retArray = std::array<MatrixView<rows, 1, prop, T>, cols> {
		view.view_col(nbrs)...
	};
	return retArray;
}

template<int trows, int tcols, typename prop, typename T>
auto MatrixView<trows, tcols, prop, T const>::cols() const -> std::array<MatrixView<trows, 1, prop, T const>, tcols> {
	return createColViews(*this, indexer<tcols>());
}

template<int trows, int tcols, typename prop, typename T>
auto MatrixView<trows, tcols, prop, T>::cols() -> std::array<MatrixView<trows, 1, prop, T>, tcols> {
	return createColViews(*this, indexer<tcols>());
}

/**
 * Determinant computation
 */
// compute 1x1 determinant
template<typename Props, typename T>
auto det(MatrixView<1, 1, Props, T const> const& mat) -> T {
	return mat(0, 0);
}

// compute 2x2 determinant
template<typename Props, typename T>
auto det(MatrixView<2, 2, Props, T const> const& mat) -> T {
	return mat(0, 0)*mat(1, 1) - mat(0, 1)*mat(1, 0);
}

// compute 3x3 determinant
template<typename Props, typename T>
auto det(MatrixView<3, 3, Props, T const> const& mat) -> T {
	return (mat(0, 0)*mat(1, 1)*mat(2, 2) +
			mat(0, 1)*mat(1, 2)*mat(2, 0) +
			mat(0, 2)*mat(1, 0)*mat(2, 1)) -
			(mat(0, 2)*mat(1, 1)*mat(2, 0) +
			mat(0, 1)*mat(1, 0)*mat(2, 2) +
			mat(0, 0)*mat(1, 2)*mat(2, 1));
}

// compute determinante of any other matrix
template <int rows, int cols, typename Props, typename T>
auto det(MatrixView<rows, cols, Props, T const> const& mat) -> T {
	T retValue = 1.;
	auto L = luDecomposition_L(mat);

	for (int i(0); i < rows; ++i) {
		retValue *= L(i, i);
	}
	return retValue;
}

// member implementation
template<int rows, int cols, typename Props, typename T>
auto MatrixView<rows, cols, Props, T const>::det() const -> T {
	return SiLi::det(*this);
}
/**
 * inverse computation
 */
// inverse of 1x1
template<typename Props, typename T>
auto inv(MatrixView<1, 1, Props, T const> const& mat) -> Matrix<1, 1, T> {
	return Matrix<1, 1, T>(T(1) / mat(0, 0));
}
// inverse of 2x2
template<typename Props, typename T>
auto inv(MatrixView<2, 2, Props, T const> const& mat) -> Matrix<2, 2, T> {
	const T c = T(1) / mat.det();
	return c * Matrix<2, 2, T>({
		{mat(1, 1), -mat(0, 1)},
		{-mat(1, 0), mat(0, 0)}
	});
}
//inverse of 3x3
template<typename Props, typename T>
auto inv(MatrixView<3, 3, Props, T const> const& mat) -> Matrix<3, 3, T> {
	const T c = T(1) / mat.det();
	Matrix<3, 3, T> ret;
	for (int row(0); row < 3; ++row) {
		for (int col(0); col < 3; ++col) {
			T tl = mat((row+1)%3, (col+1)%3);
			T br = mat((row+2)%3, (col+2)%3);
			T tr = mat((row+1)%3, (col+2)%3);
			T bl = mat((row+2)%3, (col+1)%3);
			ret(col, row) = (tl*br-tr*bl) * c;
		}
	}
	return ret;
}

// inverse of any square matrix
template <int rows, typename Props, typename T>
auto inv(MatrixView<rows, rows, Props, T const> const& _view) -> Matrix<rows, rows, T> {
	Matrix<rows, rows, T> retMat;

	auto detInv = T(1.) / _view.det();
	auto adj    = adjugateMat(_view);
	auto tran = adj.t();
	retMat      = tran * detInv;
	return retMat;
}

// member implementation
template<int trows, int tcols, typename props, typename T>
auto MatrixView<trows, tcols, props, T const>::inv() const -> Matrix<trows, tcols, T> {
	return SiLi::inv(*this);
}


/**
 * addition/substraction
 */
// matrix elementwise addition
template<int rows, int cols, typename P1, typename P2, typename T>
auto operator+(MatrixView<rows, cols, P1, T const > const& lhs, MatrixView<rows, cols, P2, T const> const& rhs) -> Matrix<rows, cols, T> {
	Matrix<rows, cols, T> ret;
	for (int row(0); row < rows; ++row) {
		for (int col(0); col < cols; ++col) {
			ret(row, col) = lhs(row, col) + rhs(row, col);
		}
	}
	return ret;
}
template<int rows, int cols, typename Props, typename T>
auto operator+(MatrixView<rows, cols, Props, T const> const& lhs, T const& rhs) -> Matrix<rows, cols, T> {
	Matrix<rows, cols, T> ret;
	for (int row(0); row < rows; ++row) {
		for (int col(0); col < cols; ++col) {
			ret(row, col) = lhs(row, col) + rhs;
		}
	}
	return ret;
}
template<int rows, int cols, typename Props, typename T>
auto operator+(T const& lhs, MatrixView<rows, cols, Props, T const> const& rhs) -> Matrix<rows, cols, T> {
	return rhs + lhs;
}

template<int rows, int cols, typename P1, typename P2, typename T>
auto operator-(MatrixView<rows, cols, P1, T const> const& lhs, MatrixView<rows, cols, P2, T const> const& rhs) -> Matrix<rows, cols, T> {
	Matrix<rows, cols, T> ret;
	for (int row(0); row < rows; ++row) {
		for (int col(0); col < cols; ++col) {
			ret(row, col) = lhs(row, col) - rhs(row, col);
		}
	}
	return ret;
}
template<int rows, int cols, typename Props, typename T>
auto operator-(MatrixView<rows, cols, Props, T const> const& lhs, T const& rhs) -> Matrix<rows, cols, T> {
	Matrix<rows, cols, T> ret;
	for (int row(0); row < rows; ++row) {
		for (int col(0); col < cols; ++col) {
			ret(row, col) = lhs(row, col) - rhs;
		}
	}
	return ret;
}
template<int rows, int cols, typename Props, typename T>
auto operator-(T const& lhs, MatrixView<rows, cols, Props, T const> const& rhs) -> Matrix<rows, cols, T> {
	Matrix<rows, cols, T> ret;
	for (int row(0); row < rows; ++row) {
		for (int col(0); col < cols; ++col) {
			ret(row, col) = lhs - rhs(row, col);
		}
	}
	return ret;
}

/**
 * matrix multiplication
 */
template<int lrows, int mate, int rcols, typename P1, typename P2, typename T>
auto operator*(MatrixView<lrows, mate, P1, T const> const& lhs,  MatrixView<mate, rcols, P2, T const> const& rhs) -> Matrix<lrows, rcols, T> {
	Matrix<lrows, rcols, T> ret;
	for (int oRow(0); oRow < lrows; ++oRow) {
		auto lhsRow = lhs.view_row(oRow);
		for (int oCol(0); oCol < rcols; ++oCol) {
			ret(oRow, oCol) = lhsRow * rhs.view_col(oCol);
		}
	}
	return ret;
}
template<int mate, typename P1, typename P2, typename T>
auto operator*(MatrixView<1, mate, P1, T const> const& lhs,  MatrixView<mate, 1, P2, T const> const& rhs) -> Matrix<1, 1, T> {
	T accumulator(0.);
	for (int i(0); i < mate; ++i) {
		accumulator += lhs(i) * rhs(i);
	}
	return Matrix<1, 1, T>(accumulator);
}
template<int rows, int cols, typename Props, typename T>
auto operator*(MatrixView<rows, cols, Props, T const> const& lhs, T const& rhs) -> Matrix<rows, cols, T> {
	Matrix<rows, cols, T> ret;
	for (int oRow(0); oRow < rows; ++oRow) {
		for (int oCol(0); oCol < cols; ++oCol) {
			ret(oRow, oCol) = lhs(oRow, oCol) * rhs;
		}
	}
	return ret;
}
template<int rows, int cols, typename Props, typename T>
auto operator*(T const& lhs, MatrixView<rows, cols, Props, T const> const& rhs) -> Matrix<rows, cols, T> {
	return rhs * lhs;
}

/**
 * cross product for 3x1 and 3x1 matrices
 */
// compute 1x1 determinant
template<typename P1, typename P2, typename T>
auto cross(MatrixView<3, 1, P1, T const> const& a, MatrixView<3, 1, P2, T const> const& b) -> Matrix<3, 1, T> {
	Matrix<3, 1, T> retValue;
	retValue(0) = a(1) * b(2) - a(2) * b(1);
	retValue(1) = a(2) * b(0) - a(0) * b(2);
	retValue(2) = a(0) * b(1) - a(1) * b(0);
	return retValue;
}


/*
 * Computing svd
 */
namespace detail {
template <typename T1, typename T2>
auto signCopy(T1 const& t1, T2 const& t2) -> T1 {
	using namespace std;
	return (t2 >= 0.0) ? abs(t1) : -abs(t1);
}

template <typename T>
auto pythag(T a, T b) -> T {
/* compute (a2 + b2)^1/2 without destructive underflow or overflow */
	using namespace std;
	a = abs(a);
	b = abs(b);
	if (a > b) {
		return a*sqrt(1.0 + (b/a) * (b/a));
	}
	else if (b == 0.0) {
		return 0.;
	} else {
		return b*sqrt(1.0+(a/b)*(a/b));
	}
}
}

template <int rows, int cols, typename Props, typename T>
auto svd(MatrixView<rows, cols, Props, T const> const& _view) -> SVD<rows, cols, T> {
	using namespace ::SiLi::detail;
	using namespace std;
	SVD<rows, cols, T> retValue;
	std::array<T, cols> rv1;

	auto& U = retValue.U;
	U = _view;
	auto& S = retValue.S;
	auto& V = retValue.V;
	V.diag() = 1.;

	T anorm, g, scale;
	g = scale = anorm = 0.0; /* Householder reduction to bidiagonal form */
	for (int i = 0; i < cols; ++i) {
		int l = i+1;
		rv1[i] = scale*g;
		g = scale = 0.0;
		T s = 0.;
		if (i < rows) {
			for (int k = i;k < rows; k++) scale += abs(U(k, i));
			if (scale) {
				for (int k = i; k < rows; k++) {
					U(k, i) /= scale;
					s += U(k, i)*U(k, i);
				}
				T f = U(i, i);
				g = -signCopy(sqrt(s),f);
				T h = f*g-s;
				U(i, i)=f-g;
				for (int j = l; j < cols; j++) {
					s=0;
					for (int k = i; k < rows; k++) s += U(k, i)*U(k, j);
					f=s/h;
					for (int k = i; k < rows; k++) U(k, j) += f*U(k, i);
				}
				for (int k = i; k < rows; k++) U(k, i) *= scale;
			}
		}
		S(i) = scale *g;
		g = s = scale = 0.0;
		if (i < rows && i+1 < cols) {
			for (int k = l; k < cols; k++) scale += abs(U(i, k));
			if (scale) {
				for (int k = l; k < cols; k++) {
					U(i, k) /= scale;
					s += U(i, k)*U(i, k);
				}
				T f = U(i, l);
				g = -signCopy(sqrt(s),f);
				T h = f*g-s;
				U(i, l)=f-g;
				for (int k = l; k < cols; k++) rv1[k]=U(i, k)/h;
				for (int j = l; j < rows; j++) {
					s=0.;
					for (int k = l; k < cols; k++) s += U(j, k)*U(i, k);
					for (int k = l; k < cols; k++) U(j, k) += s*rv1[k];
				}
				for (int k = l; k < cols; k++) U(i, k) *= scale;
			}
		}
		anorm = max(anorm,(abs(S(i))+abs(rv1[i])));
	}
	for (int i = cols-2; i >= 0; i--) { /* Accumulation of right-hand transformations. */
		int l = i+1;
		if (rv1[i+1]) {
			for (int j = l; j < cols; j++) /* Double division to avoid possible underflow. */
				V(j, i) = (U(i, j)/U(i, l)) / rv1[i+1];
			for (int j = l; j < cols; j++) {
				T s=0.;
				for (int k = l; k <cols; k++) s += U(i, k)*V(k, j);
				for (int k = l; k <cols; k++) V(k, j) += s*V(k, i);
			}
		}
		for (int j = l; j < cols; j++) V(i, j) = V(j, i) = 0.;
	}
	for (int i=min(rows, cols)-1; i >= 0; i--) { /* Accumulation of left-hand transformations. */
		int l = i+1;
		for (int j = l; j < cols; j++) U(i, j) = 0.;
		if (S(i)) {
			T w_inv = 1.0 / S(i);
			for (int j = l; j < cols; j++) {
				T s = 0.;
				for (int k = l; k < rows; k++) s += U(k, i) * U(k, j);
				T f = (s/U(i, i))*w_inv;
				for (int k = i; k < rows; k++) U(k, j) += f*U(k, i);
			}
			for (int j = i; j < rows; j++) U(j, i) *= w_inv;
		} else for (int j = i; j < rows; j++) U(j, i)=0.0;
		++U(i, i);
	}
	for (int k = cols-1; k >= 0; k--) { /* Diagonalization of the bidiagonal form. */
		for (int its = 1; its <=30; its++) {
			int flag = 1;
			int l;
			for (l = k; l >= 0; l--) { /* Test for splitting. */
				/* Note that rv1[0] is always zero. */
				if ((T)(abs(rv1[l])+anorm) == anorm) {
					flag=0;
					break;
				}
				if ((T)(abs(S(l-1))+anorm) == anorm) break;
			}
			if (flag) {
				T c = 0.0; /* Cancellation of rv1[l-1], if l > 1. */
				T s = 1.0;
				for (int i = l; i < k; i++) {
					T f = s*rv1[i];
					rv1[i] = c*rv1[i];
					if ((T)(std::abs(f)+anorm) == anorm) break;
					T h = pythag(f,S(i));
					S(i) = h;
					h = 1./h;
					c = S(i)*h;
					s = -f*h;
					for (int j = 0; j < rows; j++) {
						U(j, l-1) = U(j, l-1) * c + U(j, i) * s;
						U(j, i) = U(j, i) * c - U(j, l-1) * s;
					}
				}
			}
			if (l == k) { /* Convergence. */
				if (S(k) < 0.0) { /* Singular value is made nonnegative. */
					S(k) = -S(k);
					V.view_col(k) = -V.view_col(k);
				}
				break;
			}
			if (its == 30) {
				throw MaxIteration{};
			}
			T x = S(l); /* Shift from bottom 2-by-2 minor. */
			T y = S(k-1);
			T z = S(k);
			g = rv1[k-1];
			T h = rv1[k];
			T f = ((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
			g = pythag(f,1.0);
			f = ((x-z)*(x+z)+h*((y/(f+signCopy(g,f)))-h))/x;
			T c = 1.;
			T s = 1.; /* Next QR transformation: */
			for (int j = l; j < k; j++) {
				int i = j+1;
				g = rv1[i];
				y = S(i);
				h = s*g;
				g = c*g;
				z = pythag(f,h);
				rv1[j] = z;
				c = f/z;
				s = h/z;
				f = x*c+g*s;
				g = g*c-x*s;
				h = y*s;
				y *= c;
				for (int jj = 0; jj < cols; jj++) {
					x = V(jj, j);
					z = V(jj, i);
					V(jj, j) = x*c+z*s;
					V(jj, i) = z*c-x*s;
				}
				z = pythag(f,h);
				S(j) = z; /* Rotation can be arbitrary if z = 0. */
				if (z) {
					z = 1.0/z;
					c = f*z;
					s = h*z;
				}
				f = c*g+s*y;
				x = c*y-s*g;
				for (int jj = 0; jj < rows; jj++) {
					y = U(jj, j);
					z = U(jj, i);
					U(jj, j) = y*c+z*s;
					U(jj, i) = z*c-y*s;
				}
			}
			rv1[l] = 0.;
			rv1[k] = f;
			S(k) = x;
		}
	}

	// sort by singular values
	// slow sorting
	std::array<std::pair<T, int>, cols> columns;
	for (int i(0); i < cols; ++i) {
		columns[i] = std::make_pair(S(i), i);
	}
	sort(columns.begin(), columns.end(), [](std::pair<T, int> const& p1, std::pair<T, int> const& p2) {
		return p1.first > p2.first;
	});

	SiLi::SVD<rows, cols, T> sorted;
	for (int i(0); i < cols; ++i) {
		auto oldI = columns[i].second;
		sorted.U.view_col(i) = retValue.U.view_col(oldI);
		sorted.S(i)          = retValue.S(oldI);
		sorted.V.view_col(i) = retValue.V.view_col(oldI); // V is usually transposed
	}

	return sorted;
}

// member implementation
template<int trows, int tcols, typename props, typename T>
auto MatrixView<trows, tcols, props, T const>::svd() const -> SVD<trows, tcols, T> {
	return SiLi::svd(*this);
}

// compute abs global
template<int trows, int tcols, typename props, typename T>
auto abs(MatrixView<trows, tcols, props, T const> const& _view) -> Matrix<trows, tcols, T> {
	Matrix<trows, tcols, T> ret;
	for (int r(0); r < trows; ++r) {
		for (int c(0); c < tcols; ++c) {
			using std::abs;
			ret(r, c) = abs(_view(r, c));
		}
	}

	return ret;
}

template<int trows, int tcols, typename props, typename T>
auto sum(MatrixView<trows, tcols, props, T const> const& _view) -> T {
	T ret = 0.;
	for (int r(0); r < trows; ++r) {
		for (int c(0); c < tcols; ++c) {
			ret += _view(r, c);
		}
	}
	return ret;
}

template <int trows, int tcols, typename props, typename T>
auto isfinite(MatrixView<trows, tcols, props, T const> const& _view) -> bool {
	return _view.isfinite();
}

/*
 * ostream
 */
template<int rows, int cols, typename Prop, typename T>
std::ostream& operator<< (std::ostream& stream, SiLi::MatrixView<rows, cols, Prop, T const> const& view) {
	for (int i(0); i < view.num_rows(); ++i) {
		for (int j(0); j < view.num_cols(); ++j) {
			stream << view(i, j) << "\t";
		}
		stream << "\n";
	}
	return stream;
}




}
