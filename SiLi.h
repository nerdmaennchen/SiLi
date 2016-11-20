#pragma once

#include <type_traits>

namespace SiLi
{

template<int _rowStride, bool _transposed = false, bool _readOnly = false>
struct Properties {
	static constexpr int  rowStride  {_rowStride};
	static constexpr bool transposed {_transposed};
//	static constexpr bool readOnly   {_readOnly};   // This feature doesn't work

	using Transposed = Properties<_rowStride, not _transposed>;
//	using ReadOnly   = Properties<_rowStride, _transposed, true>;
};

template<int, int, typename T = float, typename... Args>
class Matrix;

template<int rows, int cols, typename prop, typename T>
class MatrixView {
	static_assert(rows >= 0, "Row number must be positive");
	static_assert(cols >= 0, "Coloumn number must be positive");
public:
	using Propierties = prop;
protected:
	T* const basePtr;
public:
	MatrixView(T* base) : basePtr(base) {}

	MatrixView(MatrixView& rhs) : basePtr(rhs.basePtr){}
	MatrixView(MatrixView&& rhs) : basePtr(rhs.basePtr){}

	constexpr int num_rows() const {return rows;}
	constexpr int num_cols() const {return cols;}

	template<bool t = prop::transposed, typename std::enable_if<not t>::type* = nullptr>
	T& operator()(int row, int col) & {
		return *(basePtr + (row * prop::rowStride) + col);
	}

	template<bool t = prop::transposed, typename std::enable_if<not t>::type* = nullptr>
	T const& operator()(int row, int col) const& {
		return *(basePtr + (row * prop::rowStride) + col);
	}

	template<bool t = prop::transposed, typename std::enable_if<not t>::type* = nullptr>
	T operator()(int row, int col) && {
		return *(basePtr + (row * prop::rowStride) + col);
	}
	template<bool t = prop::transposed, typename std::enable_if<t>::type* = nullptr>
	T& operator()(int col, int row) & {
		return *(basePtr + (row * prop::rowStride) + col);
	}

	template<bool t = prop::transposed, typename std::enable_if<t>::type* = nullptr>
	T const& operator()(int col, int row) const& {
		return *(basePtr + (row * prop::rowStride) + col);
	}

	template<bool t = prop::transposed, typename std::enable_if<t>::type* = nullptr>
	T operator()(int col, int row) && {
		return *(basePtr + (row * prop::rowStride) + col);
	}

	template<int subRows, int subCols>
	Matrix<subRows, subCols, T> operator()(int startR, int startC) && {
		return Matrix<subRows, subCols, T>(subview<subRows, subCols>(startR, startC));
	}

	template<int subRows, int subCols>
	Matrix<subRows, subCols, T> operator()(int startR, int startC) const& {
		Matrix<subRows, subCols, T> ret;
		for (int r(0); r < rows; ++r) {
			for (int c(0); c < cols; ++c) {
				ret(r, c) = (*this)(r, c);
			}
		}
		return ret;
	}


	template<int subRows, int subCols>
	MatrixView<subRows, subCols, prop, T> subview(int startR, int startC) {
		static_assert(subRows <= rows, "rows must be smaller or equal to the current view");
		static_assert(subCols <= cols, "cols must be smaller or equal to the current view");

		return {&((*this)(startR, startC))};
	}

	template<int subRows, int subCols>
	MatrixView<subRows, subCols, prop, T> subview(int startR, int startC) const {
		static_assert(subRows <= rows, "rows must be smaller or equal to the current view");
		static_assert(subCols <= cols, "cols must be smaller or equal to the current view");

		return {&((*this)(startR, startC))};
	}

	template<int subRows, int subCols>
	Matrix<subRows, subCols, T> submat(int startR, int startC) && {
		return Matrix<subRows, subCols, T>(subview<subRows, subCols>(startR, startC));
	}

	template<int subRows, int subCols>
	Matrix<subRows, subCols, T> submat(int startR, int startC) const& {
		Matrix<subRows, subCols, T> ret;
		for (int r(0); r < rows; ++r) {
			for (int c(0); c < cols; ++c) {
				ret(r, c) = (*this)(r, c);
			}
		}
		return ret;
	}



	MatrixView& operator*=(T const& rhs) {
		for (int r(0); r < rows; ++r) {
			for (int c(0); c < cols; ++c) {
				(*this)(r, c) *= rhs;
			}
		}
		return (*this);
	}
	MatrixView& operator+=(T const& rhs) {
		for (int r(0); r < rows; ++r) {
			for (int c(0); c < cols; ++c) {
				(*this)(r, c) += rhs;
			}
		}
		return (*this);
	}

	template <typename oProp>
	MatrixView& operator+=(MatrixView<rows, cols, oProp, T> const& rhs) {
		for (int row(0); row < rows; ++row) {
			for (int col(0); col < cols; ++col) {
				(*this)(row, col) += rhs(row, col);
			}
		}
		return (*this);
	}

	MatrixView& operator-=(T const& rhs) {
		for (int r(0); r < rows; ++r) {
			for (int c(0); c < cols; ++c) {
				(*this)(r, c) -= rhs;
			}
		}
		return (*this);
	}

	template <typename oProp>
	MatrixView& operator-=(MatrixView<rows, cols, oProp, T> const& rhs) {
		for (int row(0); row < rows; ++row) {
			for (int col(0); col < cols; ++col) {
				(*this)(row, col) -= rhs(row, col);
			}
		}
		return (*this);
	}

	Matrix<rows, cols, T> operator-() const {
		Matrix<rows, cols, T> ret;
		for (int r(0); r < rows; ++r) {
			for (int c(0); c < cols; ++c) {
				ret(r, c) = -(*this)(r, c);
			}
		}
		return ret;
	}


	MatrixView& operator=(T const& rhs) {
		for (int r(0); r < rows; ++r) {
			for (int c(0); c < cols; ++c) {
				(*this)(r, c) = rhs;
			}
		}
		return (*this);
	}

	MatrixView& operator=(MatrixView const& rhs) {
		return operator=<prop>(rhs);
	}

//	template<typename oProp, bool ReadOnly = prop::readOnly, typename std::enable_if<not ReadOnly>::type* = nullptr>
	template<typename oProp>
	MatrixView& operator=(MatrixView<rows, cols, oProp, T> const& rhs) {
		for (int r(0); r < rows; ++r) {
			for (int c(0); c < cols; ++c) {
				(*this)(r, c) = rhs(r, c);
			}
		}
		return (*this);
	}

	// only works with clang, but not with gcc?
/*	template<int oRows = rows, int oCols = cols, typename std::enable_if<(oRows * oCols) == 1>::type* = nullptr>
	operator T() const {
		return (*this)(0, 0);
	}*/

	T det() const;
	Matrix<rows, cols, T> inv() const;


	MatrixView<cols, rows, typename prop::Transposed, T> t() const {
		return {basePtr};
	}
	MatrixView<cols, rows, typename prop::Transposed, T> t() {
		return {basePtr};
	}

};

template<int rows, int cols, typename T>
class Matrix<rows, cols, T> : public MatrixView<rows, cols, Properties<cols>, T> {
	T vals[rows][cols];
	using SuperType = MatrixView<rows, cols, Properties<cols>, T>;
public:


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

	template<typename T2>
	Matrix(T2 const (&list)[rows][cols]) : SuperType(&(vals[0][0])) {
		for (int r(0); r < rows; ++r) {
			for (int c(0); c < cols; ++c) {
				(*this)(r, c) = T(list[r][c]);
			}
		}
	}

	template<typename Props>
	Matrix(MatrixView<rows, cols, Props, T> const& other) : SuperType(&(vals[0][0])) {
		for (int r(0); r < rows; ++r) {
			for (int c(0); c < cols; ++c) {
				(*this)(r, c) = other(r, c);
			}
		}
	}

	template<typename Props>
	Matrix& operator=(MatrixView<rows, cols, Props, T> const& other) & {
		SuperType::operator=(other);
		return *this;
	}

	template<typename Props>
	Matrix& operator=(MatrixView<rows, cols, Props, T> const& other) && = delete;

};


template<typename Props, typename T>
T det(MatrixView<1, 1, Props, T> const& mat) {
	return mat(0, 0);
}
template<typename Props, typename T>
T det(MatrixView<2, 2, Props, T> const& mat) {
	return mat(0, 0)*mat(1, 1) - mat(0, 1)*mat(1, 0);
}
template<typename Props, typename T>
T det(MatrixView<3, 3, Props, T> const& mat) {
	return (mat(0, 0)*mat(1, 1)*mat(2, 2) +
			mat(0, 1)*mat(1, 2)*mat(2, 0) +
			mat(0, 2)*mat(1, 0)*mat(2, 1)) -
			(mat(0, 2)*mat(1, 1)*mat(2, 0) +
			mat(0, 1)*mat(1, 0)*mat(2, 2) +
			mat(0, 0)*mat(1, 2)*mat(2, 1));
}
template<int rows, int cols, typename Props, typename T>
T MatrixView<rows, cols, Props, T>::det() const {
	return SiLi::det(*this);
}

template<typename Props, typename T>
Matrix<1, 1, T> inv(MatrixView<1, 1, Props, T> const& mat) {
	return Matrix<1, 1, T>(T(1) / mat(0, 0));
}
template<typename Props, typename T>
Matrix<2, 2, T> inv(MatrixView<2, 2, Props, T> const& mat) {
	const T c = T(1) / mat.det();
	return c * Matrix<2, 2, T>({
		{mat(1, 1), -mat(0, 1)},
		{-mat(1, 0), mat(0, 0)}
	});
}
template<typename Props, typename T>
Matrix<3, 3, T> inv(MatrixView<3, 3, Props, T> const& mat) {
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
template<int rows, int cols, typename Props, typename T>
Matrix<rows, cols, T> MatrixView<rows, cols, Props, T>::inv() const {
	return SiLi::inv(*this);
}

template<int rows, int cols, typename P1, typename P2, typename T>
Matrix<rows, cols, T> operator+(MatrixView<rows, cols, P1, T> const& lhs, MatrixView<rows, cols, P2, T> const& rhs) {
	Matrix<rows, cols, T> ret;
	for (int row(0); row < rows; ++row) {
		for (int col(0); col < cols; ++col) {
			ret(row, col) = lhs(row, col) + rhs(row, col);
		}
	}
	return ret;
}
template<int rows, int cols, typename Props, typename T>
Matrix<rows, cols, T> operator+(MatrixView<rows, cols, Props, T> const& lhs, T const& rhs) {
	Matrix<rows, cols, T> ret;
	for (int row(0); row < rows; ++row) {
		for (int col(0); col < cols; ++col) {
			ret(row, col) = lhs(row, col) + rhs;
		}
	}
	return ret;
}
template<int rows, int cols, typename Props, typename T>
Matrix<rows, cols, T> operator+(T const& lhs, MatrixView<rows, cols, Props, T> const& rhs) {
	return rhs + lhs;
}

template<int rows, int cols, typename P1, typename P2, typename T>
Matrix<rows, cols, T> operator-(MatrixView<rows, cols, P1, T> const& lhs, MatrixView<rows, cols, P2, T> const& rhs) {
	Matrix<rows, cols, T> ret;
	for (int row(0); row < rows; ++row) {
		for (int col(0); col < cols; ++col) {
			ret(row, col) = lhs(row, col) - rhs(row, col);
		}
	}
	return ret;
}
template<int rows, int cols, typename Props, typename T>
Matrix<rows, cols, T> operator-(MatrixView<rows, cols, Props, T> const& lhs, T const& rhs) {
	Matrix<rows, cols, T> ret;
	for (int row(0); row < rows; ++row) {
		for (int col(0); col < cols; ++col) {
			ret(row, col) = lhs(row, col) - rhs;
		}
	}
	return ret;
}
template<int rows, int cols, typename Props, typename T>
Matrix<rows, cols, T> operator-(T const& lhs, MatrixView<rows, cols, Props, T> const& rhs) {
	return rhs - lhs;
}

template<int lrows, int mate, int rcols, typename P1, typename P2, typename T>
Matrix<lrows, rcols, T> operator*(MatrixView<lrows, mate, P1, T> const& lhs,  MatrixView<mate, rcols, P2, T> const& rhs) {
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
Matrix<rows, cols, T> operator*(MatrixView<rows, cols, Props, T> const& lhs, T const& rhs) {
	Matrix<rows, cols, T> ret;
	for (int oRow(0); oRow < rows; ++oRow) {
		for (int oCol(0); oCol < cols; ++oCol) {
			ret(oRow, oCol) = lhs(oRow, oCol) * rhs;
		}
	}
	return ret;
}
template<int rows, int cols, typename Props, typename T>
Matrix<rows, cols, T> operator*(T const& lhs, MatrixView<rows, cols, Props, T> const& rhs) {
	return rhs * lhs;
}

}
