#pragma once
#if 0

#include <array>

namespace SiLi_Mat {

template<size_t TRows, size_t TCols, typename ...Ts>
class Matrix {
public:
	static constexpr size_t Rows = TRows;
	static constexpr size_t Cols = TCols;

	static constexpr size_t rows() noexcept {
		return Rows;
	}
	static constexpr size_t cols() noexcept {
		return Cols;
	}

private:
	static_assert(sizeof...(Ts) == TRows * TCols);
	std::tuple<Ts...> values;

public:

	template <typename ...T2>
	Matrix(T2&&...) {
	}
/*	template <size_t N0, size_t ...Ns>
	Matrix(std::array<T, N0> head, std::array<T, Ns>... values)
	{
		auto a = std::array<std::array<T, TCols>, TRows> {head, values...};
	}*/

/*	constexpr Matrix(std::array<std::array<T, TCols>, TRows> values) noexcept
		: vals {values}
	{}*/
};

template <size_t LRows, size_t LCols, size_t RRows, size_t RCols, typename ...LTs, typename ...RTs>
auto operator*(Matrix<LRows, LCols, LTs...> const& l, Matrix<RRows, RCols, RTs...> const& r) {
	static_assert(LCols == RRows);
	return Matrix<LRows, RCols>{};
}

/*template <size_t Rows, size_t Cols, typename ...Ts>
Matrix(Ts...) -> Matrix<Rows, Cols, 0, Ts...>;*/

template <size_t Rows, size_t Cols, typename ...Ts>
auto make_mat(Ts&&... values) -> Matrix<Rows, Cols, Ts...> {
	return Matrix<Rows, Cols, Ts...>(std::forward<Ts>(values)...);
}



template <typename T, int Ns>
using carray = T const (&)[Ns];

/*template <typename T, int N0, int ...Ns>
Matrix(carray<T, N0>, carray<T, Ns>... list) -> Matrix<sizeof...(list)+1, N0, T>;*/


/*template<typename T, int N0, int ...Ns>
auto g(carray<T, N0>, carray<T, Ns>... list) noexcept {
	return SiLi_Mat::Matrix<sizeof...(list)+1, N0, T>{};
}*/




/*template <int TRows, int TCols, typename T>
Matrix(std::array<std::array<T, TCols>, TRows>) -> Matrix<TRows, TCols, T>;*/
/*template <int TRows, int TCols, typename T>
Matrix(std::array<std::array<T, TCols>, TRows>) -> Matrix<TRows, TCols, T>;*/


#if 0

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
	using namespace ::std;
	return (t2 >= 0.0) ? abs(t1) : -abs(t1);
}

template <typename T>
auto pythag(T a, T b) -> T {
/* compute (a2 + b2)^1/2 without destructive underflow or overflow */
	using namespace ::std;
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
	using namespace ::std;
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
#endif

}
#endif
