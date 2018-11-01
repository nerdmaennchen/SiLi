#pragma once

#include <cassert>
#include <algorithm>
#include <array>
#include <complex>
#include <initializer_list>
#include <iterator>
#include <cstring>
#include <type_traits>
#include <vector>

namespace SiLi {

using DefaultType = float;

// This exception is thrown by svd() if the max iteration count is reached
struct MaxIteration {};

template<bool _dynamic, int _rowStride, bool _transposed = false, int _offset=0>
struct Properties {
	static constexpr bool dynamic    {_dynamic};
	static constexpr int  rowStride  {_rowStride};
	static constexpr bool transposed {_transposed};
	static constexpr int  offset     {_offset};

	using Transposed = Properties<_dynamic, _rowStride, not _transposed, _offset>;
	using Diag       = Properties<_dynamic, _rowStride, false, 1>;
	using Dynamic    = Properties<true,     _rowStride, _transposed, _offset>;
	using Static     = Properties<false,    _rowStride, _transposed, _offset>;
};

template<int trows, int tcols, int tstride, int toffset, bool staticMatrix = (trows >= 0 and tcols >= 0)>
class MatrixViewBase;

template<int trows, int tcols, typename... Args>
class MatrixView;

template<int, int, typename T = DefaultType, typename... Args>
class Matrix;

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
	MatrixIterator(MatrixView* view, bool end)
		: mat(view)
	{
		if (end) {
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

struct SizeMismatchError : public std::runtime_error {
	template<int trows1, int tcols1, int tstride1, int toffset1, bool staticMatrix1,
	         int trows2, int tcols2, int tstride2, int toffset2, bool staticMatrix2>
	SizeMismatchError(MatrixViewBase<trows1, tcols1, tstride1, toffset1, staticMatrix1> const& lhs,
	                  MatrixViewBase<trows2, tcols2, tstride2, toffset2, staticMatrix2> const& rhs,
	                  std::string _op)
	: std::runtime_error("size mismatch on operator " + _op + ": "
	                         "M<" + std::to_string(lhs.num_rows()) + "x" + std::to_string(lhs.num_cols()) + ">,"
	                         "M<" + std::to_string(rhs.num_rows()) + "x" + std::to_string(rhs.num_cols()) + ">")
	{}

	template<int trows1, int tcols1, int tstride1, int toffset1, bool staticMatrix1>
	SizeMismatchError(MatrixViewBase<trows1, tcols1, tstride1, toffset1, staticMatrix1> const& lhs,
	                  std::string _op)
	: std::runtime_error("size mismatch on operator " + _op + ": "
	                         "M<" + std::to_string(lhs.num_rows()) + "x" + std::to_string(lhs.num_cols()) + ">")
	{}

};
template<int trows1, int tcols1, int tstride1, int toffset1, bool staticMatrix1,
         int trows2, int tcols2, int tstride2, int toffset2, bool staticMatrix2>
auto sizeMismatchError(MatrixViewBase<trows1, tcols1, tstride1, toffset1, staticMatrix1> const& lhs,
                       MatrixViewBase<trows2, tcols2, tstride2, toffset2, staticMatrix2> const& rhs,
                       std::string _op) -> SizeMismatchError {
	return SizeMismatchError(lhs, rhs, _op);
}
template<int trows1, int tcols1, int tstride1, int toffset1, bool staticMatrix1>
auto sizeMismatchError(MatrixViewBase<trows1, tcols1, tstride1, toffset1, staticMatrix1> const& lhs,
                       std::string _op) -> SizeMismatchError {
	return SizeMismatchError(lhs, _op);
}


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

	SVD(int _rows, int _cols)
		: U(_rows, _cols, T(0.))
		, S(_cols, 1, T(0.))
		, V(_cols, _cols, T(0.))
	{}
};


namespace detail
{

template<int... I>
struct index {
	template<int n>
	using prepend = index<n, I...>;
	constexpr int size() const {
		return sizeof...(I);
	}
};

template<int start, int end, int stride, bool anchor=(start>end)>
struct make_increment {
	using type = typename make_increment<start+stride, end, stride>::type::template prepend<start>;
};

template<int start, int end, int stride>
struct make_increment<start, end, stride, true> {
	using type = index<>;
};


template<int start, int end, int stride, bool anchor=(start<end)>
struct make_decrement {
	using type = typename make_decrement<start+stride, end, stride>::type::template prepend<start>;
};

template<int start, int end, int stride>
struct make_decrement<start, end, stride, true> {
	using type = index<>;
};

template<int start, int end, int stride, typename std::enable_if<(stride>0)>::type* = nullptr>
constexpr auto indexer() -> typename make_increment<start, end, stride, (start>end or stride<0)>::type { return {}; }
template<int start, int end, int stride, typename std::enable_if<(stride<0)>::type* = nullptr>
constexpr auto indexer() -> typename make_decrement<start, end, stride, (start<end or stride>0)>::type {
	return {};
}


template<int rows, int cols, typename prop, typename T, int ...nbrs>
auto createRowViews(MatrixView<rows, cols, prop, T const> const& view, detail::index<nbrs...>) -> std::array<MatrixView<1, cols, prop, T const>, sizeof...(nbrs)> {
	auto retArray = std::array<MatrixView<1, cols, prop, T const>, sizeof...(nbrs)> {
		view.view_row(nbrs)...
	};
	return retArray;
}
template<int rows, int cols, typename prop, typename T, int ...nbrs>
auto createRowViews(MatrixView<rows, cols, prop, T>& view, detail::index<nbrs...>) -> std::array<MatrixView<1, cols, prop, T>, sizeof...(nbrs)> {
	auto retArray = std::array<MatrixView<1, cols, prop, T>, sizeof...(nbrs)> {
		view.view_row(nbrs)...
	};
	return retArray;
}
template<int rows, int cols, typename prop, typename T, int ...nbrs>
auto createColViews(MatrixView<rows, cols, prop, T const> const& view, detail::index<nbrs...>) -> std::array<MatrixView<rows, 1, prop, T const>, sizeof...(nbrs)> {
	auto retArray = std::array<MatrixView<rows, 1, prop, T const>, sizeof...(nbrs)> {
		view.view_col(nbrs)...
	};
	return retArray;
}
template<int rows, int cols, typename prop, typename T, int ...nbrs>
auto createColViews(MatrixView<rows, cols, prop, T>& view, detail::index<nbrs...>) -> std::array<MatrixView<rows, 1, prop, T>, sizeof...(nbrs)> {
	auto retArray = std::array<MatrixView<rows, 1, prop, T>, sizeof...(nbrs)> {
		view.view_col(nbrs)...
	};
	return retArray;
}
}

template<int trows, int tcols, int tstride, int toffset>
class MatrixViewBase<trows, tcols, tstride, toffset, true> {
public:
	constexpr int num_rows() const { return trows; }
	constexpr int num_cols() const { return tcols; }
	constexpr int stride() const { return tstride; }
	constexpr int offset() const { return toffset; }
};

template<int trows, int tcols, int tstride, int toffset>
class MatrixViewBase<trows, tcols, tstride, toffset, false> {
protected:
	int mRows;
	int mCols;
	int mStride;
	int mOffset;
	template<int, int, typename...> friend class MatrixView;
	template<int, int, typename, typename...> friend class Matrix;
public:
	int num_rows() const { return mRows; }
	int num_cols() const { return mCols; }
	int stride() const { return mStride; }
	int offset() const { return mOffset; }
};

constexpr auto add_size(int v1, int v2) -> int {
	return (v1 >= 0 and v2 >= 0)?(v1+v2):-1;
}


/** Const MatrixView with all methods that are allowed to be operated on a const matrix view
 */
template<int trows, int tcols, typename prop, typename T>
class MatrixView<trows, tcols, prop, T const> : public MatrixViewBase<trows, tcols, prop::rowStride, prop::offset> {
public:
	using Type = T;
	using Props = prop;

private:
	using TMatrix = Matrix<trows, tcols, T>;
	using Base    = MatrixViewBase<trows, tcols, prop::rowStride, prop::offset>;

protected:
	T const* cBasePtr;

	void changeBase(T const* base) { cBasePtr = base; }

public:
	explicit MatrixView(T const* base) : cBasePtr(base) {}

	MatrixView(MatrixView const& rhs) : cBasePtr(rhs.cBasePtr) {}
	MatrixView(MatrixView&& rhs)      : cBasePtr(rhs.cBasePtr) {}

	template<typename T2>
	explicit MatrixView(MatrixView<trows, tcols, prop, T2> const& rhs) : cBasePtr(rhs.cBasePtr) {}

	using Base::offset;
	using Base::stride;

	// read only element access not transposed
	template<bool t = prop::transposed, typename std::enable_if<not t>::type* = nullptr>
	auto operator()(int row, int col) const -> T const& {
		return *(cBasePtr + (row * stride()) + col + offset()*row);
	}

	// read only element access transposed
	template<bool t = prop::transposed, typename std::enable_if<t>::type* = nullptr>
	auto operator()(int row, int col) const -> T const& {
		return operator()<false>(col, row);
	}

	// read only element access not transposed
	template<int _rows = Number<trows>::_number, typename std::enable_if<_rows == 1 and tcols != 1>::type* = nullptr>
	auto operator()(int col) const -> T const& {
		return operator()(0, col);
	}

	// read only element access not transposed
	template<int _cols = Number<tcols>::_number, typename std::enable_if<_cols == 1>::type* = nullptr>
	auto operator()(int row) const -> T const& {
		return operator()(row, 0);
	}

	template<int _cols = tcols, typename std::enable_if<trows == _cols and trows == 1>::type* = nullptr>
	operator T() const {
		return (*this)(0, 0);
	};

	template<int _cols = tcols, typename std::enable_if<trows == _cols and trows < 0>::type* = nullptr>
	explicit operator T() const {
		if (this->num_rows() != 1 or this->num_cols() != 1) {
			throw SizeMismatchError(*this, *this, "operator T");
		}
		return (*this)(0, 0);
	};




	// matrix access
	template<int subRows, int subCols>
	auto mat(int startR, int startC) const -> Matrix<subRows, subCols, T> {
		return view<subRows, subCols>(startR, startC);
	}

	// read only view access
	template<int subRows, int subCols>
	auto view(int startR, int startC) const -> MatrixView<subRows, subCols, typename prop::Static, T const> {
		static_assert(subRows >= 0, "Row number must be positive");
		static_assert(subCols >= 0, "Coloumn number must be positive");
		static_assert(subRows <= trows, "rows must be smaller or equal to the current view");
		static_assert(subCols <= tcols, "cols must be smaller or equal to the current view");

		return MatrixView<subRows, subCols, typename prop::Static, T const> {&((*this)(startR, startC))};
	}

	// read only view access
	auto view(int startR, int startC, int rows, int cols) const -> MatrixView<-1, -1, typename prop::Dynamic, T const> {
		auto view = MatrixView<-1, -1, typename prop::Dynamic, T const> {&((*this)(startR, startC))};
		view.mRows = rows;
		view.mCols = cols;
		view.mStride = this->stride();
		view.mOffset = this->offset();
		return view;
	}

	auto operator()(int startR, int startC, int rows, int cols) const -> decltype(this->view(startR, startC, rows, cols)) {
		return view(startR, startC, rows, cols);
	}

	template <bool dynamic = prop::dynamic, typename std::enable_if<not dynamic>::type* = nullptr>
	auto view_row(int startR) const -> MatrixView<1, tcols, prop, T const> {
		return MatrixView<1, tcols, prop, T const> {&((*this)(startR, 0))};
	}

	template <bool dynamic = prop::dynamic, typename std::enable_if<not dynamic>::type* = nullptr>
	auto view_col(int startC) const -> MatrixView<trows, 1, prop, T const> {
		return MatrixView<trows, 1, prop, T const> {&((*this)(0, startC))};
	}

	template <bool dynamic = prop::dynamic, typename std::enable_if<dynamic>::type* = nullptr>
	auto view_row(int startR) const -> MatrixView<1, tcols, prop, T const> {
		auto view = MatrixView<1, tcols, prop, T const> {&((*this)(startR, 0))};
		view.mRows   = 1;
		view.mCols   = this->mCols;
		view.mStride = this->mStride;
		view.mOffset = this->offset();
		return view;
	}

	template <bool dynamic = prop::dynamic, typename std::enable_if<dynamic>::type* = nullptr>
	auto view_col(int startC) const -> MatrixView<trows, 1, prop, T const> {
		auto view = MatrixView<trows, 1, prop, T const> {&((*this)(0, startC))};
		view.mRows   = this->mRows;
		view.mCols   = 1;
		view.mStride = this->mStride;
		view.mOffset = this->offset();
		return view;
	}


	// array with view on a range of rows (inclusive end)
	template<int start=0, int end=trows-1, int stride=1>
	auto rows() const -> decltype(detail::createRowViews(*this, detail::indexer<start, end, stride>())) {
		static_assert(start >= 0 and end >= 0, "bounds must be positive");
		static_assert(start < trows and end < trows, "bounds must be smaller than the number of rows");
		static_assert(stride != 0, "stride must be different from 0");
		return createRowViews(*this, detail::indexer<start, end, stride>());
	}

	// array with view on a range of cols (inclusive end)
	template<int start=0, int end=tcols-1, int stride=1>
	auto cols() const -> decltype(detail::createColViews(*this, detail::indexer<start, end, stride>())) {
		static_assert(start >= 0 and end >= 0, "bounds must be positive");
		static_assert(start < tcols and end < tcols, "bounds must be smaller than the number of rows");
		static_assert(stride != 0, "stride must be different from 0");
		return createColViews(*this, detail::indexer<start, end, stride>());
	}

	// read only iterator
	auto begin() const -> MatrixIterator<MatrixView<trows, tcols, prop, T const> const> {
		return {this, false};
	}
	auto end() const -> MatrixIterator<MatrixView<trows, tcols, prop, T const> const> {
		return {this, true};
	}

	template <int cols, typename oProp, typename std::enable_if<not Props::dynamic and not oProp::dynamic>::type* = nullptr>
	auto join_rows(MatrixView<trows, cols, oProp, T const> const& _view) const -> Matrix<trows, add_size(tcols, cols), T> {
		Matrix<trows, add_size(tcols, cols), T> retValue(this->num_rows(), this->num_cols() + _view.num_cols());
		retValue.template view<trows, tcols>(0, 0)    = *this;
		retValue.template view<trows, cols>(0, tcols) = _view;
		return retValue;
	}

	template <int cols, typename oProp, typename std::enable_if<Props::dynamic or oProp::dynamic>::type* = nullptr>
	auto join_rows(MatrixView<trows, cols, oProp, T const> const& _view) const -> Matrix<trows, add_size(tcols, cols), T> {

		if (this->num_rows() != _view.num_rows()) {
			throw SizeMismatchError(*this, _view, "join_rows");
		}

		Matrix<trows, add_size(tcols, cols), T> ret(this->num_rows(), this->num_cols() + _view.num_cols());
		ret.view(0, 0, this->num_rows(), this->num_cols()) = *this;
		ret.view(0, this->num_cols(), _view.num_rows(), _view.num_cols()) = _view;
		return ret;
	}


	template <int rows, typename oProp, typename std::enable_if<not Props::dynamic and not oProp::dynamic>::type* = nullptr>
	auto join_cols(MatrixView<rows, tcols, oProp, T const> const& _view) const -> Matrix<add_size(trows, rows), tcols, T> {
		Matrix<add_size(trows, rows), tcols, T> retValue;
		retValue.template view<trows, tcols>(0, 0)    = *this;
		retValue.template view<rows,  tcols>(trows, 0) = _view;
		return retValue;
	}

	template <int rows, typename oProp, typename std::enable_if<Props::dynamic or oProp::dynamic>::type* = nullptr>
	auto join_cols(MatrixView<rows, tcols, oProp, T const> const& _view) const -> Matrix<add_size(trows, rows), tcols, T> {
		if (this->num_cols() != _view.num_cols()) {
			throw SizeMismatchError(*this, _view, "join_cols");
		}

		Matrix<add_size(trows, rows), tcols, T> ret(this->num_rows() + _view.num_rows(), this->num_cols());
		ret.view(0, 0, this->num_rows(), this->num_cols()) = *this;
		ret.view(this->num_rows(), 0, _view.num_rows(), _view.num_cols()) = _view;
		return ret;
	}


	// negation
	auto operator-() const -> Matrix<trows, tcols, T> {
		Matrix<trows, tcols, T> ret = *this;
		for (auto& v : ret) {
			v = -v;
		}
		return ret;
	}

	// compute inverse
	auto inv() const -> Matrix<trows, tcols, T>;

	// compute pseudo inverse
	auto pinv(T epsilon = 0.) const -> Matrix<tcols, trows, T>;

	// determinant
	auto det() const -> T;

	// read only transposed
	template <bool dynamic = prop::dynamic, typename std::enable_if<not dynamic>::type* = nullptr>
	auto t_view() const -> MatrixView<tcols, trows, typename prop::Transposed, T const> {
		return MatrixView<tcols, trows, typename prop::Transposed, T const> {cBasePtr};
	}

	template <bool dynamic = prop::dynamic, typename std::enable_if<dynamic>::type* = nullptr>
	auto t_view() const -> MatrixView<tcols, trows, typename prop::Transposed, T const> {
		auto view = MatrixView<tcols, trows, typename prop::Transposed, T const> {cBasePtr};
		view.mRows   = this->mCols;
		view.mCols   = this->mRows;
		view.mStride = this->mStride;
		view.mOffset = this->offset();
		return view;
	}


	auto t() const -> Matrix<tcols, trows, T> {
		return t_view();
	}


	// read only diagonal view
	template <bool dynamic = prop::dynamic, typename std::enable_if<not dynamic>::type* = nullptr>
	auto diag() const -> MatrixView<(tcols < trows)?tcols:trows, 1, typename prop::Diag, T const> {
		return MatrixView<(tcols < trows)?tcols:trows, 1, typename prop::Diag, T const> {cBasePtr};
	}

	template <bool dynamic = prop::dynamic, typename std::enable_if<dynamic>::type* = nullptr>
	auto diag() const -> MatrixView<-1, 1, typename prop::Diag, T const> {
		auto view = MatrixView<-1, 1, typename prop::Diag, T const> {cBasePtr};
		view.mRows   = std::min(this->num_rows(), this->num_cols());
		view.mCols   = 1;
		view.mStride = this->mStride;
		view.mOffset = 1;
		return view;
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
		using namespace std;
		return sqrt(normSqr());
	}

	auto normalized() const -> Matrix<trows, tcols, T> {
		return (*this) * (T(1.)/norm());
	}

	// compute svd so that *this = svd.U * make_diag(svd.S) * svd.V.t()
	// Will throw SiLi::MaxIteration if maximum of iteration is reached
	auto svd() const -> SVD<trows, tcols, T>;

	auto operator^(T _exponent) const -> Matrix<trows, tcols, T>;

	auto rank(T _threshold = 0.) const -> int;

	// compute abs
	auto abs() const -> Matrix<trows, tcols, T> {
		return abs(*this);
	}

	// compute sum
	auto sum() const -> T {
		return sum(*this);
	}

	// compute prod
	auto prod() const -> T {
		return prod(*this);
	}

	// compute tanh
	auto tanh() const -> T {
		Matrix<trows, tcols, T> ret = *this;
		using std::tanh;
		for (auto& e : ret) {
			e = tanh(e);
		}
		return ret;
	}

	// compute atanh
	auto atanh() const -> Matrix<trows, tcols, T> {
		Matrix<trows, tcols, T> ret = *this;
		using std::atanh;
		for (auto& e : ret) {
			e = atanh(e);
		}
		return ret;
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
	using Base  = MatrixViewBase<trows, tcols, prop::rowStride, prop::offset>;

protected:
	T* basePtr;

public:
	explicit MatrixView(T* base) : CView(base), basePtr(base) {}

	MatrixView(MatrixView& rhs) : CView(rhs), basePtr(rhs.basePtr) {}
	MatrixView(MatrixView&& rhs) : CView(rhs), basePtr(rhs.basePtr) {}

	void changeBase(T* base) { CView::changeBase(base); basePtr = base; }

	// pass through
	using CView::stride;
	using CView::offset;
	using CView::operator();

	// value element access
	auto operator()(int row, int col) -> T& {
		if constexpr (prop::transposed) {
			std::swap(row, col);
		}
		return *(basePtr + (row * stride()) + col + row * offset());
	}

	// access if matrix is one dimensional
	template<int _rows = Number<trows>::_number, typename std::enable_if<_rows == 1 and tcols != 1>::type* = nullptr>
	auto operator()(int col) -> T& {
		return operator()(0, col);
	}

	// access if matrix is one dimensional
	template<int _cols = Number<tcols>::_number, typename std::enable_if<_cols == 1>::type* = nullptr>
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

	// view access
	auto view(int startR, int startC, int rows, int cols) -> MatrixView<-1, -1, typename prop::Dynamic, T> {
		auto view = MatrixView<-1, -1, typename prop::Dynamic, T> {&((*this)(startR, startC))};
		view.mRows = rows;
		view.mCols = cols;
		view.mStride = this->stride();
		view.mOffset = this->offset();
		return view;
	}

	auto operator()(int startR, int startC, int rows, int cols) -> decltype(this->view(startR, startC, rows, cols)) {
		return view(startR, startC, rows, cols);
	}


	using CView::view_row;
	using CView::view_col;

	template <bool dynamic = prop::dynamic, typename std::enable_if<not dynamic>::type* = nullptr>
	auto view_row(int startR) -> MatrixView<1, tcols, prop, T> {
		return MatrixView<1, tcols, prop, T> {&((*this)(startR, 0))};
	}

	template <bool dynamic = prop::dynamic, typename std::enable_if<not dynamic>::type* = nullptr>
	auto view_col(int startC) -> MatrixView<trows, 1, prop, T> {
		return MatrixView<trows, 1, prop, T> {&((*this)(0, startC))};
	}

	template <bool dynamic = prop::dynamic, typename std::enable_if<dynamic>::type* = nullptr>
	auto view_row(int startR) -> MatrixView<1, tcols, prop, T> {
		auto view = MatrixView<1, tcols, prop, T> {&((*this)(startR, 0))};
		view.mRows   = 1;
		view.mCols   = this->mCols;
		view.mStride = this->mStride;
		view.mOffset = this->offset();
		return view;
	}

	template <bool dynamic = prop::dynamic, typename std::enable_if<dynamic>::type* = nullptr>
	auto view_col(int startC) -> MatrixView<trows, 1, prop, T> {
		auto view = MatrixView<trows, 1, prop, T> {&((*this)(0, startC))};
		view.mRows   = this->mRows;
		view.mCols   = 1;
		view.mStride = this->mStride;
		view.mOffset = this->offset();
		return view;
	}


	// array with view on a range of rows (inclusive end)
	template<int start=0, int end=trows-1, int stride=1>
	auto rows() -> decltype(detail::createRowViews(*this, detail::indexer<start, end, stride>())) {
		static_assert(start >= 0 and end >= 0, "bounds must be positive");
		static_assert(start < trows and end < trows, "bounds must be smaller than the number of rows");
		static_assert(stride != 0, "stride must be different from 0");
		return detail::createRowViews(*this, detail::indexer<start, end, stride>());
	}

	// array with view on a range of cols (inclusive end)
	template<int start=0, int end=tcols-1, int stride=1>
	auto cols() -> decltype(detail::createColViews(*this, detail::indexer<start, end, stride>())) {
		static_assert(start >= 0 and end >= 0, "bounds must be positive");
		static_assert(start < tcols and end < tcols, "bounds must be smaller than the number of rows");
		static_assert(stride != 0, "stride must be different from 0");
		return detail::createColViews(*this, detail::indexer<start, end, stride>());
	}



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
		if (this->num_rows() != rhs.num_rows() or
		    this->num_cols() != rhs.num_cols()) {
			throw SizeMismatchError(*this, rhs, "operator+=");
		}

		for (int row(0); row < this->num_rows(); ++row) {
			for (int col(0); col < this->num_cols(); ++col) {
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
		for (int row(0); row < this->num_rows(); ++row) {
			for (int col(0); col < this->num_cols(); ++col) {
				(*this)(row, col) -= rhs(row, col);
			}
		}
		return *this;
	}

	// element wise product
	template<typename oProp>
	auto operator&=(MatrixView<trows, tcols, oProp, T const> const& _view) -> MatrixView& {
		for (int row(0); row < this->num_rows(); ++row) {
			for (int col(0); col < this->num_cols(); ++col) {
				(*this)(row, col) *= _view(row, col);
			}
		}
		return *this;
	}

	// element wise product
	template<typename oProp>
	auto operator&(MatrixView<trows, tcols, oProp, T const> const& _view) -> Matrix<trows, tcols, T> {
		Matrix<trows, tcols, T> retMat;
		for (int row(0); row < this->num_rows(); ++row) {
			for (int col(0); col < this->num_cols(); ++col) {
				retMat(row, col) = (*this)(row, col) * _view(row, col);
			}
		}
		return retMat;
	}


	auto operator=(T const& rhs) -> MatrixView& {
		for (int r(0); r < this->num_rows(); ++r) {
			for (int c(0); c < this->num_cols(); ++c) {
				(*this)(r, c) = rhs;

			}
		}
		return *this;
	}

	//auto operator=(MatrixView const& rhs) const -> MatrixView& = delete;

	auto operator=(MatrixView const& rhs) -> MatrixView& {
		return operator=<trows, tcols, prop>(rhs);
	}

	template <int _rows, int _cols, typename oProp>
	auto operator=(MatrixView<_rows, _cols, oProp, T const> const& rhs) -> MatrixView& {
		if (this->num_rows() != rhs.num_rows()
			or this->num_cols() != rhs.num_cols()) {
			throw SizeMismatchError(*this, rhs, "operator=");
		}

		for (int r(0); r < this->num_rows(); ++r) {
			for (int c(0); c < this->num_cols(); ++c) {
				(*this)(r, c) = rhs(r, c);
			}
		}
		return *this;
	}

	auto begin() -> MatrixIterator<MatrixView<trows, tcols, prop, T>> {
		return {this, false};
	}
	auto end() -> MatrixIterator<MatrixView<trows, tcols, prop, T>> {
		return {this, true};
	}

	// transposed view
	using CView::t_view;

	template <bool dynamic = prop::dynamic, typename std::enable_if<not dynamic>::type* = nullptr>
	auto t_view() -> MatrixView<tcols, trows, typename prop::Transposed, T> {
		return MatrixView<tcols, trows, typename prop::Transposed, T> {basePtr};
	}

	template <bool dynamic = prop::dynamic, typename std::enable_if<dynamic>::type* = nullptr>
	auto t_view() -> MatrixView<tcols, trows, typename prop::Transposed, T> {
		auto view = MatrixView<tcols, trows, typename prop::Transposed, T> {basePtr};
		view.mRows   = this->mCols;
		view.mCols   = this->mRows;
		view.mStride = this->mStride;
		view.mOffset = this->offset();
;		return view;
	}


	// diagonal view
	using CView::diag;

	template <bool dynamic = prop::dynamic, typename std::enable_if<not dynamic>::type* = nullptr>
	auto diag() -> MatrixView<(tcols < trows)?tcols:trows, 1, typename prop::Diag, T> {
		return MatrixView<(tcols < trows)?tcols:trows, 1, typename prop::Diag, T> {basePtr};
	}

	template <bool dynamic = prop::dynamic, typename std::enable_if<dynamic>::type* = nullptr>
	auto diag() -> MatrixView<-1, 1, typename prop::Diag, T> {
		auto view = MatrixView<-1, 1, typename prop::Diag, T> {basePtr};
		view.mRows   = std::min(this->num_rows(), this->num_cols());
		view.mCols   = 1;
		view.mStride = this->mStride;
		view.mOffset = 1;
		return view;
	}

};

template<int rows, int cols, typename T>
class Matrix<rows, cols, T> : public MatrixView<rows, cols, class Properties<false, cols>, T> {
	std::array<std::array<T, cols>, rows> vals;
public:
	using View  = MatrixView<rows, cols, Properties<false, cols>, T>;
	using CView = MatrixView<rows, cols, Properties<false, cols>, T const>;

	template <typename Node>
	void serialize(Node& node) {
		node["values"] % vals;
	}


	Matrix() : View(&(vals[0][0])) {}

	Matrix(int /*_rows*/, int /*_cols*/) : Matrix() {
		// this function is needed for compatiblity of dynamic matrices
	}
	Matrix(int /*_rows*/, int /*_cols*/, T _initValue) : Matrix(_initValue) {
		// this function is needed for compatiblity of dynamic matrices
	}


	explicit Matrix(T initVal) : View(&(vals[0][0]))  {
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
		memcpy(&(vals[0][0]), &(list[0][0]), rows*cols*sizeof(T));
	}

	Matrix(Matrix&& other) : View(&(vals[0][0])) {
		memcpy(&(vals[0][0]), &(other.vals[0][0]), rows*cols*sizeof(T));
	}


	Matrix(Matrix const& other) : View(&(vals[0][0])) {
		memcpy(&(vals[0][0]), &(other.vals[0][0]), rows*cols*sizeof(T));
	}

	template<typename Props>
	Matrix(MatrixView<rows, cols, Props, T const> const& other) : View(&(vals[0][0])) {
		for (int r(0); r < this->num_rows(); ++r) {
			for (int c(0); c < this->num_cols(); ++c) {
				(*this)(r, c) = other(r, c);
			}
		}
	}

	using View::operator=;
	template<typename Props>
	auto operator=(MatrixView<rows, cols, Props, T const> const& other) && -> Matrix& = delete;

	auto operator=(Matrix const& other) -> Matrix& {
		memcpy(&vals[0][0], &other.vals[0][0], rows*cols*sizeof(T));
		return *this;
	}
};

template<typename T>
class Matrix<-1, -1, T> : public MatrixView<-1, -1, class Properties<true, -1>, T> {
	std::vector<T> vals;
protected:
	void resize(int _rows, int _cols, int _stride) {
		vals.resize(_rows * _cols);
		this->changeBase(vals.data());
		this->mRows = _rows;
		this->mCols = _cols;
		this->mStride = _stride;
		this->mOffset = 0;
	}
	void resize(int _rows, int _cols, int _stride, std::vector<T> values) {
		vals = std::move(values);
		this->changeBase(vals.data());
		this->mRows = _rows;
		this->mCols = _cols;
		this->mStride = _stride;
		this->mOffset = 0;
	}

public:
	using View  = MatrixView<-1, -1, Properties<true, -1>, T>;
	using CView = MatrixView<-1, -1, Properties<true, -1>, T const>;

	template <typename Node>
	void serialize(Node& node) {
		auto rows   = this->mRows;
		auto cols   = this->mCols;
		auto stride = this->mStride;
		node["rows"]   % rows;
		node["cols"]   % cols;
		node["stride"] % stride;

		if (rows != this->mRows or cols != this->mCols or stride != this->mStride) {
			resize(rows, cols, stride);
		}

		node["values"] % vals;
	}

	Matrix(int rows, int cols) : View(vals.data()) {
		resize(rows, cols, cols);
	}

	Matrix(int rows, int cols, T _initVal) : View(vals.data()) {
		resize(rows, cols, cols);
		((View*)(this))->operator=(_initVal);
	}

	Matrix() : View(vals.data()) {}

	Matrix(std::vector<T> values) : View(vals.data()) {
		int rows = values.size();
		resize(rows, 1, 1, std::move(values));
	}

	Matrix(std::array<std::vector<T>, 1> values) : View(vals.data()) {
		int cols = values.size();
		resize(1, cols, cols, std::move(values));
	}


	Matrix(Matrix&& other) : View(vals.data()) {
		resize(other.num_rows(), other.num_cols(), other.stride(), std::move(other.vals));
	}


	Matrix(Matrix const& other) : View(vals.data()) {
		resize(other.num_rows(), other.num_cols(), other.stride(), other.vals);
	}

	template<int oth_trows, int oth_tcols, typename Props>
	Matrix(MatrixView<oth_trows, oth_tcols, Props, T const> const& other) : View(vals.data()) {
		resize(other.num_rows(), other.num_cols(), other.num_cols());
		for (int r(0); r < other.num_rows(); ++r) {
			for (int c(0); c < other.num_cols(); ++c) {
				(*this)(r, c) = other(r, c);
			}
		}
	}

	auto operator=(T const& rhs) -> Matrix& {
		View::operator=(rhs);
		return *this;
	}

	template<int oth_trows, int oth_tcols, typename Props>
	auto operator=(MatrixView<oth_trows, oth_tcols, Props, T const> const& other) & -> Matrix& {
		resize(other.num_rows(), other.num_cols(), other.num_cols());
		for (int r(0); r < other.num_rows(); ++r) {
			for (int c(0); c < other.num_cols(); ++c) {
				(*this)(r, c) = other(r, c);
			}
		}
		return *this;
	}

	//template<int oth_trows, int oth_tcols, typename Props>
	//auto operator=(MatrixView<oth_trows, oth_tcols, Props, T const> const& other) && -> Matrix& = delete;

	auto operator=(Matrix const& other) -> Matrix& {
		resize(other.num_rows(), other.num_cols(), other.num_cols());
		for (int r(0); r < other.num_rows(); ++r) {
			for (int c(0); c < other.num_cols(); ++c) {
				(*this)(r, c) = other(r, c);
			}
		}
		return *this;
	}
};

template<int trows, typename T>
class Matrix<trows, -1, T> : public Matrix<-1, -1, T> {
public:
	using Matrix<-1, -1, T>::Matrix;

	using Matrix<-1, -1, T>::operator();

	// read only element access
	template<int _rows = Number<trows>::_number, typename std::enable_if<_rows == 1>::type* = nullptr>
	auto operator()(int col) const -> T const& {
		return operator()(0, col);
	}

	// element access
	template<int _rows = Number<trows>::_number, typename std::enable_if<_rows == 1>::type* = nullptr>
	auto operator()(int col) -> T& {
		return operator()(0, col);
	}

	using Matrix<-1, -1, T>::operator=;
};

template<int tcols, typename T>
class Matrix<-1, tcols, T> : public Matrix<-1, -1, T> {
public:
	using Matrix<-1, -1, T>::Matrix;

	using Matrix<-1, -1, T>::operator();

	// read only element access
	template<int _cols = Number<tcols>::_number, typename std::enable_if<_cols == 1>::type* = nullptr>
	auto operator()(int row) const -> T const& {
		return operator()(row, 0);
	}

	// element access
	template<int _cols = Number<tcols>::_number, typename std::enable_if<_cols == 1>::type* = nullptr>
	auto operator()(int row) -> T& {
		return operator()(row, 0);
	}

	using Matrix<-1, -1, T>::operator=;
};




// create matrix
template<int rows, int cols, typename T = DefaultType>
auto make_mat(T const (&values)[rows][cols]) -> Matrix<rows, cols, T> {
	return {values};
}

// create column vector
template<int rows, typename T = DefaultType>
auto make_vec(T const (&values)[rows]) -> Matrix<rows, 1, T> {
	return {values};
}

// create matrix
template<int rows, int cols, typename T = DefaultType>
auto make_mat() -> Matrix<rows, cols, T> {
	return {0.};
}

// create matrix
template<typename T = DefaultType>
auto make_mat(std::initializer_list<std::initializer_list<T>> _lists) -> Matrix<-1, -1, T> {
	int cols = 0;
	int rows = _lists.size();
	if (rows > 0) {
		cols = _lists.begin()->size();
	}
	Matrix<-1, -1, T> mat(rows, cols);
	int rowIdx = 0;
	for (auto row : _lists) {
		int colIdx = 0;
		for (auto col : row) {
			mat(rowIdx, colIdx) = col;
			++colIdx;
		}
		++rowIdx;
	}
	return mat;
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
	Matrix<rows, rows, T> retVal(_view.num_rows(), _view.num_rows());
	retVal.operator=(T(0.));
	retVal.diag() = _view;
	return retVal;
}

template<typename T, typename Prop>
auto make_diag(MatrixView<-1, -1, Prop, T const> const& _view) -> Matrix<-1, -1, T> {
	if (_view.num_cols() != 1) {
		throw SizeMismatchError(_view, "make_diag");
	}
	auto m = std::max(_view.num_rows(), _view.num_cols());
	Matrix<-1, -1, T> retVal(m, m);
	retVal.operator=(T(0.));
	retVal.diag() = _view;
	return retVal;
}

template<int rows, int cols, typename T, typename Prop>
auto make_diag(MatrixView<rows, 1, Prop, T const> const& _view) -> Matrix<rows, cols, T> {
	static_assert(cols >= rows and cols > -1 and rows > -1, "cannot build a smaller diagonal matrix than the input vector size");
	Matrix<rows, cols, T> retVal(_view.num_rows(), cols);
	retVal = T(0.);
	retVal.diag() = _view;
	return retVal;
}
template<int rows, int cols, typename T, typename Prop>
auto make_diag(MatrixView<cols, 1, Prop, T const> const& _view) -> Matrix<rows, cols, T> {
	static_assert(rows >= cols and cols > -1 and rows > -1, "cannot build a smaller diagonal matrix than the input vector size");
	Matrix<rows, cols, T> retVal(_view.num_rows(), cols);
	retVal = T(0.);
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
	Matrix<rows, cols, T> L(_mat.num_rows(), _mat.num_cols());
	L = T(0.);

	for (int k=0; k<_mat.num_rows(); ++k) {
		auto kRow = L.view_row(k);
		auto kCol = L.view_col(k);
		for (int i=k;i<_mat.num_rows();++i) {
			L(k, i) = _mat(k, i) - kRow * L.view_col(i);
		}
		for (int i=k+1; i<_mat.num_rows();++i) {
			L(i, k) = (_mat(i, k) - T(L.view_row(i) * kCol)) / L(k, k);
		}
	}
	return L;
}

// compute minor matrix
template <int rows, int cols, typename Props, typename T>
auto minorMat(MatrixView<rows, cols, Props, T const> const& _view, int _row, int _col) -> Matrix<(rows>=0)?rows-1:-1, (cols>=0)?cols-1:-1, T> {

	auto retRows = _view.num_rows()-1;
	auto retCols = _view.num_cols()-1;

	Matrix<(rows>=0)?rows-1:-1, (cols>=0)?cols-1:-1, T> retMat(retRows, retCols);

	retMat(   0,    0,         _row, _col)         = _view(     0,      0,         _row, _col);
	retMat(_row,    0, retRows-_row, _col)         = _view(_row+1,      0, retRows-_row, _col);
	retMat(   0, _col,         _row, retCols-_col) = _view(     0, _col+1,         _row, retCols-_col);
	retMat(_row, _col, retRows-_row, retCols-_col) = _view(_row+1, _col+1, retRows-_row, retCols-_col);

	return retMat;
}

// compute adjugated matrix
template <int rows, int cols, typename Props, typename T>
auto adjugateMat(MatrixView<rows, cols, Props, T const> const& _view) -> Matrix<rows, cols, T> {

	Matrix<rows, cols, T> retMat(_view.num_rows(), _view.num_cols());

	if (cols == 1) {
		retMat = T(0.);
		return retMat;
	}

	for(int i = 0; i < _view.num_rows(); ++i) {
		for(int j = 0; j < _view.num_cols(); ++j) {
			T sign = (i%2*2-1) * (j%2*2-1);
			auto det = minorMat(_view, i, j).det();
			retMat(i, j) = sign * det;
		}
	}
	return retMat;
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

	for (int i(0); i < L.num_rows(); ++i) {
		retValue *= L(i, i);
	}
	using std::isfinite;
	if (not isfinite(retValue)) {
		retValue = T{0.};
	}
	return retValue;
}

// member implementation
template<int rows, int cols, typename Props, typename T>
auto MatrixView<rows, cols, Props, T const>::det() const -> T {
	return SiLi::det(*this);
}

template<int trows, int tcols, typename Props, typename T>
auto fastInv(MatrixView<trows, tcols, Props, const T> const& _mat, T _epsilon = 0.) ->  Matrix<tcols, trows, T>{
	using std::abs;

	auto svd = _mat.svd();
	for (auto& e : svd.S) {
		if (abs(e) > _epsilon) {
			e = T(1.)/e;
		}
	}
	auto S = make_diag(svd.S);
	return svd.V * S * svd.U.t_view();
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
	if (_view.num_rows() != _view.num_cols()) {
		throw SizeMismatchError(_view, "inv");
	}
	return fastInv(_view);
/*
	Matrix<rows, rows, T> retMat(_view.num_rows(), _view.num_cols());

	auto detInv = T(1.) / _view.det();
	auto adj    = adjugateMat(_view);
	auto tran   = adj.t_view();
	retMat      = tran * detInv;
	return retMat;*/
}

// member implementation
template<int trows, int tcols, typename props, typename T>
auto MatrixView<trows, tcols, props, T const>::inv() const -> Matrix<trows, tcols, T> {
	return SiLi::inv(*this);
}

// member implementation
template<int trows, int tcols, typename props, typename T>
auto MatrixView<trows, tcols, props, T const>::pinv(T _epsilon) const -> Matrix<tcols, trows, T> {
	return fastInv(*this, _epsilon);
}

/**
 * addition/substraction
 */
// matrix elementwise addition
template<int rows, int cols, typename P1, typename P2, typename T>
auto operator+(MatrixView<rows, cols, P1, T const > const& lhs, MatrixView<rows, cols, P2, T const> const& rhs) -> Matrix<rows, cols, T> {
	Matrix<rows, cols, T> ret(lhs.num_rows(), lhs.num_cols());
	for (int row(0); row < lhs.num_rows(); ++row) {
		for (int col(0); col < lhs.num_cols(); ++col) {
			ret(row, col) = lhs(row, col) + rhs(row, col);
		}
	}
	return ret;
}
template<int rows, int cols, typename Props, typename T>
auto operator+(MatrixView<rows, cols, Props, T const> const& lhs, T const& rhs) -> Matrix<rows, cols, T> {
	Matrix<rows, cols, T> ret;(lhs.num_rows(), lhs.num_cols());
	for (int row(0); row < lhs.num_rows(); ++row) {
		for (int col(0); col < lhs.num_cols(); ++col) {
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
	Matrix<rows, cols, T> ret(lhs.num_rows(), lhs.num_cols());
	for (int row(0); row < lhs.num_rows(); ++row) {
		for (int col(0); col < lhs.num_cols(); ++col) {
			ret(row, col) = lhs(row, col) - rhs(row, col);
		}
	}
	return ret;
}
template<int rows, int cols, typename Props, typename T>
auto operator-(MatrixView<rows, cols, Props, T const> const& lhs, T const& rhs) -> Matrix<rows, cols, T> {
	Matrix<rows, cols, T> ret(lhs.num_rows(), lhs.num_cols());
	for (int row(0); row < lhs.num_rows(); ++row) {
		for (int col(0); col < lhs.num_cols(); ++col) {
			ret(row, col) = lhs(row, col) - rhs;
		}
	}
	return ret;
}
template<int rows, int cols, typename Props, typename T>
auto operator-(T const& lhs, MatrixView<rows, cols, Props, T const> const& rhs) -> Matrix<rows, cols, T> {
	Matrix<rows, cols, T> ret(rhs.num_rows(), rhs.num_cols());
	for (int row(0); row < rhs.num_rows(); ++row) {
		for (int col(0); col < rhs.num_cols(); ++col) {
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
	if (lhs.num_cols() != rhs.num_rows()) {
		throw SizeMismatchError(lhs, rhs, "multiplikation");
	}
	Matrix<lrows, rcols, T> ret(lhs.num_rows(), rhs.num_cols());
	for (int oRow(0); oRow < lhs.num_rows(); ++oRow) {
		auto lhsRow = lhs.view_row(oRow);
		for (int oCol(0); oCol < rhs.num_cols(); ++oCol) {
			ret(oRow, oCol) = T(lhsRow * rhs.view_col(oCol));
		}
	}
	return ret;
}
template<int mate, typename P1, typename P2, typename T>
auto operator*(MatrixView<1, mate, P1, T const> const& lhs,  MatrixView<mate, 1, P2, T const> const& rhs) -> Matrix<1, 1, T> {
	T accumulator(0.);
	for (int i(0); i < lhs.num_cols(); ++i) {
		accumulator += lhs(i) * rhs(i);
	}
	return Matrix<1, 1, T>(accumulator);
}
template<int rows, int cols, typename Props, typename T>
auto operator*(MatrixView<rows, cols, Props, T const> const& lhs, T const& rhs) -> Matrix<rows, cols, T> {
	Matrix<rows, cols, T> ret = lhs;
	for (auto& v : ret) {
		v *= rhs;
	}
	return ret;
}
template<int rows, int cols, typename Props, typename T>
auto operator*(T const& lhs, MatrixView<rows, cols, Props, T const> const& rhs) -> Matrix<rows, cols, T> {
	return rhs * lhs;
}
//element wise product
// matrix elementwise addition
template<int rows, int cols, typename P1, typename P2, typename T>
auto eProd(MatrixView<rows, cols, P1, T const > const& lhs, MatrixView<rows, cols, P2, T const> const& rhs) -> Matrix<rows, cols, T> {
	Matrix<rows, cols, T> ret(lhs.num_rows(), lhs.num_cols());
	for (int row(0); row < lhs.num_rows(); ++row) {
		for (int col(0); col < lhs.num_cols(); ++col) {
			ret(row, col) = lhs(row, col) * rhs(row, col);
		}
	}
	return ret;
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
		return a*sqrt(T(1.0) + (b/a) * (b/a));
	}
	else if (b == 0.0) {
		return 0.;
	} else {
		return b*sqrt(T(1.0)+(a/b)*(a/b));
	}
}
}

template <typename T, int length, bool dynamic=(length==-1)>
struct DynList : public std::array<T, length> {
	DynList(int _length) {}
};

template <typename T, int length>
struct DynList<T, length, true> : public std::vector<T> {
	DynList(int _length)
		: std::vector<T>(_length)
	{}
};


template <int rows, int cols, typename T, bool dynamic=(rows==-1) or (cols==-1)>
auto sortSVD(int _rows, int _cols, SVD<rows, cols, T> const& _svd) -> SVD<rows, cols, T> {
	// sort by singular values
	// slow sorting
	DynList<std::pair<T, int>, cols> columns(_cols);
	for (int i(0); i < _cols; ++i) {
		columns[i] = std::make_pair(_svd.S(i), i);
	}
	sort(columns.begin(), columns.end(), [](std::pair<T, int> const& p1, std::pair<T, int> const& p2) {
		return p1.first > p2.first;
	});

	SiLi::SVD<rows, cols, T> sorted(_rows, _cols);
	assert(_cols == sorted.U.num_cols());
	assert(_cols == sorted.V.num_cols());
	assert(_cols == sorted.S.num_rows());
	assert(_cols == columns.size());

	for (int i(0); i < _cols; ++i) {
		auto oldI = columns[i].second;
		sorted.U.view_col(i) = _svd.U.view_col(oldI);
		sorted.S(i)          = _svd.S(oldI);
		sorted.V.view_col(i) = _svd.V.view_col(oldI); // V is usually transposed
	}
	return sorted;
}


template <int _rows, int _cols, typename Props, typename T>
auto svd(MatrixView<_rows, _cols, Props, T const> const& _view) -> SVD<_rows, _cols, T> {
	using namespace ::SiLi::detail;
	using namespace std;
	int const rows = _view.num_rows();
	int const cols = _view.num_cols();
	SVD<_rows, _cols, T> retValue(rows, cols);
	Matrix<_cols, 1, T> rv1(rows , 1, T(0.));

	auto& U = retValue.U;
	U = _view;
	auto& S = retValue.S;
	auto& V = retValue.V;
	V.diag() = 1.;

	T anorm{0.};
	/* Householder reduction to bidiagonal form */
	for (int i = 0; i < cols and i < rows; ++i) {
		{
			auto row_view = U.view(i, i, rows-i, 1);
			auto scale = sum(abs(row_view));
			if (scale) {
				row_view *= T(1.) / scale;
				auto s = row_view.normSqr();

				auto f = U(i, i);
				auto g = -signCopy(sqrt(s), f);
				auto h = f * g - s;
				U(i, i) = f - g;
				for (int j = i+1; j < cols; j++) {
					auto row_view2 = U.view(i, j, rows-i, 1);
					auto s = T(row_view.t_view() * row_view2);
					row_view2 += row_view * (s/h);
				}
				row_view *= scale;
				S(i) = scale *g;
			}
		}
		int l = i+1;
		if (l < cols) {
			auto col_view = U.view(i, l, 1, cols-l);
			auto scale = sum(abs(col_view));
			if (scale) {
				col_view *= T(1.) / scale;
				auto s = col_view.normSqr();
				auto g = -signCopy(sqrt(s), U(i, l));
				auto h = U(i, l) * g - s;
				U(i, l) -= g;

				auto rv_view = rv1.view(l, 0, cols-l, 1);
				rv_view = col_view.t_view() * (T(1.) / h);

				for (int j = l; j < rows; j++) {
					auto col_view2 = U.view(j, l, 1, cols-l);
					auto s = T(col_view2 * col_view.t_view());
					col_view2 += rv_view.t_view() * s;
				}
				col_view *= scale;
				if (i+1 < cols) {
					rv1(i+1) = scale*g;
				}
			}
		}
		anorm = max(anorm,(abs(S(i))+abs(rv1(i))));
	}

	for (int i = cols-2; i >= 0; i--) { /* Accumulation of right-hand transformations. */
		int l = i+1;
		if (rv1(i+1)) {
			for (int j = l; j < cols; j++) /* Double division to avoid possible underflow. */
				V(j, i) = (U(i, j)/U(i, l)) / rv1(i+1);
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
			T w_inv = T(1.0) / S(i);
			for (int j = l; j < cols; j++) {
				T s = 0.;
				for (int k = l; k < rows; k++) s += U(k, i) * U(k, j);
				T f = (s/U(i, i))*w_inv;
				for (int k = i; k < rows; k++) U(k, j) += f*U(k, i);
			}
			U(i, i, rows-i, 1) *= w_inv;
		} else {
			U(i, i, rows-i, 1) = 0.;
		}
		++U(i, i);
	}
	for (int k = cols-1; k >= 0; k--) { /* Diagonalization of the bidiagonal form. */
		for (int its = 1; its <=30; its++) {
			int flag = 1;
			int l;
			for (l = k; l >= 0; l--) { /* Test for splitting. */
				/* Note that rv1[0] is always zero. */
				if ((abs(rv1(l))+anorm) == anorm) {
					flag=0;
					break;
				}
				if ((abs(S(l-1))+anorm) == anorm) break;
			}
			if (flag) {
				T c = 0.0; /* Cancellation of rv1[l-1], if l > 1. */
				T s = 1.0;
				for (int i = l; i < k; i++) {
					T f = s*rv1(i);
					rv1(i) = c*rv1(i);
					if (abs(f)+anorm == anorm) break;
					auto t1 = S(i);
					S(i) = pythag(f, t1);
					auto h = T(1.)/S(i);
					c = t1*h;
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
			T g = rv1(k-1);
			T h = rv1(k);
			T f = ((y-z)*(y+z)+(g-h)*(g+h))/(T(2.0)*h*y);
			g = pythag(f, T(1.0));
			f = ((x-z)*(x+z)+h*((y/(f+signCopy(g,f)))-h))/x;
			T c = 1.;
			T s = 1.; /* Next QR transformation: */
			for (int j = l; j < k; j++) {
				int i = j+1;
				g = rv1(i);
				y = S(i);
				h = s*g;
				g = c*g;
				z = pythag(f,h);
				rv1(j) = z;
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
					z = T(1.0)/z;
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
			rv1(l) = 0.;
			rv1(k) = f;
			S(k) = x;
		}
	}

	return sortSVD(_view.num_rows(), _view.num_cols(), retValue);
}

// member implementation
template<int trows, int tcols, typename props, typename T>
auto MatrixView<trows, tcols, props, T const>::svd() const -> SVD<trows, tcols, T> {
	return SiLi::svd(*this);
}

template<int trows, int tcols, typename props, typename T>
auto MatrixView<trows, tcols, props, T const>::operator^(T _exponent) const -> Matrix<trows, tcols, T> {
	using std::pow;

	auto tsvd = svd();
	for (auto& e : tsvd.S) {
		e = pow(e, _exponent);
	}
	auto S = make_diag(tsvd.S);
	return tsvd.U * S * tsvd.V.t_view();
}

template<int trows, int tcols, typename props, typename T>
auto MatrixView<trows, tcols, props, T const>::rank(T _threshold) const -> int {
	using std::abs;

	int r = 0;

	auto tsvd = svd();
	for (auto& e : tsvd.S) {
		if (abs(e) > _threshold) {
			r += 1;
		} else {
			break;
		}
	}
	return r;
}




// compute abs global
template<int trows, int tcols, typename props, typename T>
auto abs(MatrixView<trows, tcols, props, T const> const& _view) -> Matrix<trows, tcols, T> {
	using std::abs;
	Matrix<trows, tcols, T> ret = _view;
	for (auto& v : ret) {
		v = abs(v);
	}
	return ret;
}

template<int trows, int tcols, typename props, typename T>
auto sum(MatrixView<trows, tcols, props, T const> const& _view) -> T {
	T ret = 0.;
	for (int r(0); r < _view.num_rows(); ++r) {
		for (int c(0); c < _view.num_cols(); ++c) {
			ret += _view(r, c);
		}
	}
	return ret;
}

template<int trows, int tcols, typename props, typename T>
auto prod(MatrixView<trows, tcols, props, T const> const& _view) -> T {
	T ret = 1.;
	for (int r(0); r < _view.num_rows(); ++r) {
		for (int c(0); c < _view.num_cols(); ++c) {
			ret *= _view(r, c);
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
	if (rows != 1) {
		stream << "{\n";
	}
	for (int i(0); i < view.num_rows(); ++i) {
		stream << "{ ";
		for (int j(0); j < view.num_cols(); ++j) {
			stream << view(i, j);
			if (j+1 < view.num_cols()) {
				stream << ", ";
			}
		}
		stream << " }";
		if (i+1 < view.num_rows()) {
			stream << ",\n";
		}
	}
	if (rows != 1) {
		stream << "}";
	}
	return stream;
}

template<int rows, typename T = DefaultType>
using Vector = Matrix<rows, 1, T>;

template <int rows, int cols, typename P, typename T>
auto toVec(MatrixView<rows, cols, P, T const> const& _in) -> Matrix<rows*cols, 1, T> {
	Matrix<rows*cols, 1, T> ret;
	int i(0);
	for (auto const& e : _in) {
		ret(i) = e;
		++i;
	}
	return ret;
}
template <int rows, int cols, typename P, typename T>
auto toMat(MatrixView<rows*cols, 1, P, T const> const& _in) -> Matrix<rows, cols, T> {
	Matrix<rows, cols, T> ret;
	int i(0);
	for (auto& e : ret) {
		e = _in(i);
		++i;
	}
	return ret;
}


}
