#pragma once

#include <algorithm>
#include <array>
#include <complex>
#include <iterator>
#include <type_traits>

namespace SiLi
{

template<int _rowStride, bool _transposed = false, int _offset=0>
struct Properties {
	static constexpr int  rowStride  {_rowStride};
	static constexpr bool transposed {_transposed};
	static constexpr int  offset     {_offset};

	using Transposed = Properties<_rowStride, not _transposed, _offset>;
	using Diag       = Properties<_rowStride, false, 1>;
};

template<int, int, typename T = float, typename... Args>
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
		static_assert(subRows <= trows, "rows must be smaller or equal to the current view");
		static_assert(subCols <= tcols, "cols must be smaller or equal to the current view");

		Matrix<subRows, subCols, T> ret;
		for (int r(startR); r < trows; ++r) {
			for (int c(startC); c < tcols; ++c) {
				ret(r-startR, c-startC) = (*this)(r, c);
			}
		}
		return ret;
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
	auto begin() const -> MatrixIterator<MatrixView<trows, tcols, prop, T const>> {
		return {this};
	}
	auto end() const -> MatrixIterator<MatrixView<trows, tcols, prop, T const>> {
		return {nullptr};
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

	auto svd() const -> SVD<trows, tcols, T>;

};
template<int trows, int tcols, typename prop, typename T>
class MatrixView<trows, tcols, prop, T> : public MatrixView<trows, tcols, prop, T const> {
private:
	using Parent = MatrixView<trows, tcols, prop, T const>;
protected:
	T* const basePtr;

public:
	explicit MatrixView(T* base) : Parent(base), basePtr(base) {}

	MatrixView(MatrixView& rhs) : Parent(rhs), basePtr(rhs.basePtr) {}
	MatrixView(MatrixView&& rhs) : Parent(rhs), basePtr(rhs.basePtr) {}

	//MatrixView(MatrixView const& rhs) : basePtr(rhs.basePtr) {}

	// pass through
	using MatrixView<trows, tcols, prop, T const>::operator();
/*	auto operator()(int row, int col) const& -> T {
		return Parent::operator()(row, col);
	}*/

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
	using MatrixView<trows, tcols, prop, T const>::t;
	auto t() -> MatrixView<tcols, trows, typename prop::Transposed, T> {
		return MatrixView<tcols, trows, typename prop::Transposed, T> {basePtr};
	}

	// diagonal view
	auto diag() -> MatrixView<(tcols < trows)?tcols:trows, 1, typename prop::Diag, T> {
		return MatrixView<(tcols < trows)?tcols:trows, 1, typename prop::Diag, T> {basePtr};
	}
};

template<int rows, int cols, typename T>
class Matrix<rows, cols, T> : public MatrixView<rows, cols, class Properties<cols>, T> {
	T vals[rows][cols];
	using SuperType = MatrixView<rows, cols, Properties<cols>, T>;
public:
	using View  = SuperType;
	using CView = MatrixView<rows, cols, Properties<cols>, T const>;


	Matrix() : SuperType(&(vals[0][0])) {}

	Matrix(T initVal) : SuperType(&(vals[0][0]))  {
		((SuperType*)(this))->operator=(initVal);
	}

	template<typename T2>
	Matrix(T2 const (&list)[rows*cols]) : SuperType(&(vals[0][0])) {
		T2 const* iter = &(list[0]);
		for (int r(0); r < rows; ++r) {
			for (int c(0); c < cols; ++c) {
				(*this)(r, c) = T(*iter++);
			}
		}
	}

	Matrix(T const (&list)[rows][cols]) : SuperType(&(vals[0][0])) {
		for (int r(0); r < rows; ++r) {
			for (int c(0); c < cols; ++c) {
				(*this)(r, c) = T(list[r][c]);
			}
		}
	}

	Matrix(Matrix&& other) : SuperType(&(vals[0][0])) {
		for (int r(0); r < rows; ++r) {
			for (int c(0); c < cols; ++c) {
				(*this)(r, c) = other(r, c);
			}
		}
	}


	Matrix(Matrix const& other) : SuperType(&(vals[0][0])) {
		for (int r(0); r < rows; ++r) {
			for (int c(0); c < cols; ++c) {
				(*this)(r, c) = other(r, c);
			}
		}
	}


	template<typename Props>
	Matrix(MatrixView<rows, cols, Props, T const> const& other) : SuperType(&(vals[0][0])) {
		for (int r(0); r < rows; ++r) {
			for (int c(0); c < cols; ++c) {
				(*this)(r, c) = other(r, c);
			}
		}
	}

	template<typename Props>
	auto operator=(MatrixView<rows, cols, Props, T const> const& other) & -> Matrix& {
		SuperType::operator=(other);
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

	auto operator=(Matrix&& other) -> Matrix& {
		for (int r(0); r < rows; ++r) {
			for (int c(0); c < cols; ++c) {
				(*this)(r, c) = other(r, c);
			}
		}
		return *this;
	}
};

// create matrix
template<typename T, int rows, int cols>
auto make_mat(T const (&values)[rows][cols]) -> Matrix<rows, cols, T> {
	return {values};
}

// create identity matrix
template<typename T, int rows, int cols>
auto make_eye() -> Matrix<rows, cols, T> {
	Matrix<rows, cols, T> retVal(0);
	retVal.diag() = 1.;
	return retVal;
}

// create matrix from vector
template<typename T, typename Prop, int rows>
auto make_diag(MatrixView<rows, 1, Prop, T const> const& _view) -> Matrix<rows, rows, T> {
	Matrix<rows, rows, T> retVal(0);
	retVal.diag() = _view;
	return retVal;
}

// only works for not overlaping views
template<int rows, int cols, typename P1, typename P2, typename T>
void swap(MatrixView<rows, cols, P1, T>& lhs, MatrixView<rows, cols, P2, T>& rhs) {
	using std::swap;
	for(int i = 0; i < rows; ++i) {
		for(int j = 0; j < cols; ++j) {
			swap(lhs(i, j), rhs(i, j));
		}
	}
}


// lu decomposition, returns L value
template <int rows, int cols, typename Props, typename T>
auto luDecomposition_L(MatrixView<rows, cols, Props, T const> const& _mat) -> Matrix<rows, cols, T> {
	Matrix<rows, cols, T> L;

	for (int k=0; k<rows; ++k) {
		for (int j=k;j<rows;++j) {
			double sum =0.;
			for (int p=0; p < k; ++p) sum += L(k, p) * L(p, j);
			L(k, j) = _mat(k, j) - sum;
		}
		for (int i=k+1; i<rows;++i) {
			double sum = 0.;
			for (int p=0; p<k;++p) sum += L(i, p) * L(p, k);
			L(i, k) = (_mat(i, k) - sum) / L(k, k);
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
		return _view;
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
		for (int oCol(0); oCol < rcols; ++oCol) {
			T akkumulator = 0;
			for (int iCol(0); iCol < mate; ++iCol) {
				akkumulator += lhs(oRow, iCol) * rhs(iCol, oCol);
			}
			ret(oRow, oCol) = akkumulator;
		}
	}
	return ret;
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

/*
 * Computing svd
 */
namespace detail {
template <typename T1, typename T2>
auto signCopy(T1 const& t1, T2 const& t2) -> T1 {
	return (t2 >= 0.0) ? std::abs(t1) : -std::abs(t1);
}

template <typename T>
auto pythag(T a, T b) -> T {
	using namespace std;
/* compute (a2 + b2)^1/2 without destructive underflow or overflow */
	a = std::abs(a);
	b = std::abs(b);
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
	int const m = rows;
	int const n = cols;

	auto& a = retValue.U;
	a = _view;
	auto& w = retValue.S;
	auto& v = retValue.V;

	int j,l;
	T anorm,c,f,g,h,s,scale,x,y,z;

	std::array<T, cols> rv1;
	g=scale=anorm=0.0; /* Householder reduction to bidiagonal form */
	for (int i = 1; i <= n; ++i) {
		l=i+1;
		rv1[i-1] = scale*g;
		g = s = scale = 0.0;
		if (i <= m) {
			for (int k=i;k<=m;k++) scale += abs(a(k-1, i-1));
			if (scale) {
				for (int k=i;k<=m;k++) {
					a(k-1, i-1) /= scale;
					s += a(k-1, i-1)*a(k-1, i-1);
				}
				f=a(i-1, i-1);
				g = -signCopy(sqrt(s),f);
				h=f*g-s;
				a(i-1, i-1)=f-g;
				for (j=l;j<=n;j++) {
					s=0;
					for (int k=i;k<=m;k++) s += a(k-1, i-1)*a(k-1, j-1);
					f=s/h;
					for (int k=i;k<=m;k++) a(k-1, j-1) += f*a(k-1, i-1);
				}
				for (int k=i;k<=m;k++) a(k-1, i-1) *= scale;
			}
		}
		w(i-1)=scale *g;
		g=s=scale=0.0;
		if (i <= m && i != n) {
			for (int k=l;k<=n;k++) scale += abs(a(i-1, k-1));
			if (scale) {
				for (int k=l;k<=n;k++) {
					a(i-1, k-1) /= scale;
					s += a(i-1, k-1)*a(i-1, k-1);
				}
				f=a(i-1, l-1);
				g = -signCopy(sqrt(s),f);
				h=f*g-s;
				a(i-1, l-1)=f-g;
				for (int k=l;k<=n;k++) rv1[k-1]=a(i-1, k-1)/h;
				for (j=l;j<=m;j++) {
					s=0.;
					for (int k=l;k<=n;k++) s += a(j-1, k-1)*a(i-1, k-1);
					for (int k=l;k<=n;k++) a(j-1, k-1) += s*rv1[k-1];
				}
				for (int k=l;k<=n;k++) a(i-1, k-1) *= scale;
			}
		}
		anorm = std::max(anorm,(abs(w(i-1))+abs(rv1[i-1])));
	}
	for (int i=n;i>=1;i--) { /* Accumulation of right-hand transformations. */
		if (i < n) {
			if (g) {
				for (j=l;j<=n;j++) /* Double division to avoid possible underflow. */
					v(j-1, i-1)=(a(i-1, j-1)/a(i-1, l-1))/g;
				for (j=l;j<=n;j++) {
					s=0.;
					for (int k=l;k<=n;k++) s += a(i-1, k-1)*v(k-1, j-1);
					for (int k=l;k<=n;k++) v(k-1, j-1) += s*v(k-1, i-1);
				}
			}
			for (j=l;j<=n;j++) v(i-1, j-1)=v(j-1, i-1)=0.0;
		}
		v(i-1, i-1)=1.0;
		g=rv1[i-1];
		l=i;
	}
	for (int i=std::min(m,n);i>=1;i--) { /* Accumulation of left-hand transformations. */
		l=i+1;
		g=w(i-1);
		for (j=l;j<=n;j++) a(i-1, j-1)=0.0;
		if (g) {
			g=1.0/g;
			for (j=l;j<=n;j++) {
				s = 0.;
				for (int k=l;k<=m;k++) s += a(k-1, i-1)*a(k-1, j-1);
				f=(s/a(i-1, i-1))*g;
				for (int k=i;k<=m;k++) a(k-1, j-1) += f*a(k-1, i-1);
			}
			for (j=i;j<=m;j++) a(j-1, i-1) *= g;
		} else for (j=i;j<=m;j++) a(j-1, i-1)=0.0;
		++a(i-1, i-1);
	}
	int nm;
	for (int k=n;k>=1;k--) { /* Diagonalization of the bidiagonal form. */
		for (int its=1;its<=30;its++) {
			int flag=1;
			for (l=k;l>=1;l--) { /* Test for splitting. */
				nm=l-1; /* Note that rv1[0] is always zero. */
				if ((T)(abs(rv1[l-1])+anorm) == anorm) {
					flag=0;
					break;
				}
				if ((T)(abs(w(nm-1))+anorm) == anorm) break;
			}
			if (flag) {
				c=0.0; /* Cancellation of rv1[l-1], if l > 1. */
				s=1.0;
				for (int i=l;i<=k;i++) {
					f=s*rv1[i-1];
					rv1[i-1]=c*rv1[i-1];
					if ((T)(abs(f)+anorm) == anorm) break;
					g=w(i-1);
					h=pythag(f,g);
					w(i-1)=h;
					h=1.0/h;
					c=g*h;
					s = -f*h;
					for (j=1;j<=m;j++) {
						y=a(j-1, nm-1);
						z=a(j-1, i-1);
						a(j-1, nm-1)=y*c+z*s;
						a(j-1, i-1)=z*c-y*s;
					}
				}
			}
			z=w(k-1);
			if (l == k) { /* Convergence. */
				if (z < 0.0) { /* Singular value is made nonnegative. */
					w(k-1) = -z;
					for (j=1;j<=n;j++) v(j-1, k-1) = -v(j-1, k-1);
				}
				break;
			}
			if (its == 30) {
				throw std::runtime_error("no convergence after many svd iterations");
			}
			x=w(l-1); /* Shift from bottom 2-by-2 minor. */
			nm=k-1;
			y=w(nm-1);
			g=rv1[nm-1];
			h=rv1[k-1];
			f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
			g=pythag(f,1.0);
			f=((x-z)*(x+z)+h*((y/(f+signCopy(g,f)))-h))/x;
			c=s=1.0; /* Next QR transformation: */
			for (j=l;j<=nm;j++) {
				int i=j+1;
				g=rv1[i-1];
				y=w(i-1);
				h=s*g;
				g=c*g;
				z=pythag(f,h);
				rv1[j-1]=z;
				c=f/z;
				s=h/z;
				f=x*c+g*s;
				g = g*c-x*s;
				h=y*s;
				y *= c;
				for (int jj=1;jj<=n;jj++) {
					x=v(jj-1, j-1);
					z=v(jj-1, i-1);
					v(jj-1, j-1)=x*c+z*s;
					v(jj-1, i-1)=z*c-x*s;
				}
				z=pythag(f,h);
				w(j-1)=z; /* Rotation can be arbitrary if z = 0. */
				if (z) {
					z=1.0/z;
					c=f*z;
					s=h*z;
				}
				f=c*g+s*y;
				x=c*y-s*g;
				for (int jj=1;jj<=m;jj++) {
					y=a(jj-1, j-1);
					z=a(jj-1, i-1);
					a(jj-1, j-1)=y*c+z*s;
					a(jj-1, i-1)=z*c-y*s;
				}
			}
			rv1[l-1] = 0.;
			rv1[k-1] = f;
			w(k-1)=x;
		}
	}

	// sort by eigenvalues
	// slow sorting

	std::array<std::pair<T, int>, cols> columns;
	for (int i(0); i < cols; ++i) {
		columns[i] = std::make_pair(w(i), i);
	}
	std::sort(columns.begin(), columns.end(), [](std::pair<T, int> const& p1, std::pair<T, int> const& p2) {
		return p1.first > p2.first;
	});

	SiLi::SVD<rows, cols, T> sorted;
	for (int i(0); i < cols; ++i) {
		auto oldI = columns[i].second;
		sorted.U.view_col(i) = retValue.U.view_col(oldI);
		sorted.S(i)          = retValue.S(oldI);
		sorted.V.view_row(i) = retValue.V.view_row(oldI);
	}

	return sorted;
}

// member implementation
template<int trows, int tcols, typename props, typename T>
auto MatrixView<trows, tcols, props, T const>::svd() const -> SVD<trows, tcols, T> {
	return SiLi::svd(*this);
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
