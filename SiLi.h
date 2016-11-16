#pragma once

// to disable pretty printed std::cout of views comment out this define
#define IO_STREAM_PRESENT

namespace SiLi
{

template<int, int, typename T = float, typename... Args>
class Matrix;

template<int rows, int cols, int rowStride, typename T>
class MatrixView {
	static_assert(rows >= 0, "Row number must be positive");
	static_assert(cols >= 0, "Coloumn number must be positive");
protected:
	T* const basePtr;
public:
	MatrixView(T* base) : basePtr(base) {}
	MatrixView(MatrixView const& rhs) : basePtr(rhs.basePtr){}

	constexpr int num_rows() {return rows;}
	constexpr int num_cols() {return cols;}

	T& operator()(int row, int col) & {
		return *(basePtr + (row * rowStride) + col);
	}
	T const& operator()(int row, int col) const& {
		return *(basePtr + (row * rowStride) + col);
	}

	T operator()(int row, int col) && {
		return *(basePtr + (row * rowStride) + col);
	}


	template<int subRows, int subCols>
	MatrixView<subRows, subCols, rowStride, T> subview(int startR, int startC) & {
		return MatrixView<subRows, subCols, rowStride, T>(&((*this)(startR, startC)));
	}

	template<int subRows, int subCols>
	Matrix<subRows, subCols, T> submat(int startR, int startC) && {
		Matrix<cols, rows, T> ret;
		for (int r(0); r < subRows; ++r) {
			for (int c(0); c < subCols; ++c) {
				ret(r, c) = (*this)(r + startR, c + startC);
			}
		}
		return ret;
	}


	MatrixView<rows, cols, rowStride, T>& operator*=(T const& rhs) {
		for (int r(0); r < rows; ++r) {
			for (int c(0); c < cols; ++c) {
				(*this)(r, c) *= rhs;
			}
		}
		return (*this);
	}
	MatrixView<rows, cols, rowStride, T>& operator+=(T const& rhs) {
		for (int r(0); r < rows; ++r) {
			for (int c(0); c < cols; ++c) {
				(*this)(r, c) += rhs;
			}
		}
		return (*this);
	}

	template <int oRowStride>
	MatrixView& operator+=(MatrixView<rows, cols, oRowStride, T> const& rhs) {
		for (int row(0); row < rows; ++row) {
			for (int col(0); col < cols; ++col) {
				(*this)(row, col) += rhs(row, col);
			}
		}
		return (*this);
	}

	MatrixView<rows, cols, rowStride, T>& operator-=(T const& rhs) {
		for (int r(0); r < rows; ++r) {
			for (int c(0); c < cols; ++c) {
				(*this)(r, c) -= rhs;
			}
		}
		return (*this);
	}

	template <int oRowStride>
	MatrixView& operator-=(MatrixView<rows, cols, oRowStride, T> const& rhs) {
		for (int row(0); row < rows; ++row) {
			for (int col(0); col < cols; ++col) {
				(*this)(row, col) -= rhs(row, col);
			}
		}
		return (*this);
	}

	MatrixView<rows, cols, rowStride, T>& operator=(T const& rhs) {
		for (int r(0); r < rows; ++r) {
			for (int c(0); c < cols; ++c) {
				(*this)(r, c) = rhs;
			}
		}
		return (*this);
	}

	MatrixView<rows, cols, rowStride, T>& operator=(MatrixView<rows, cols, rowStride, T> const& rhs) {
		return operator=<rowStride>(rhs);
	}

	template<int oRowStride>
	MatrixView<rows, cols, rowStride, T>& operator=(MatrixView<rows, cols, rowStride, T> const& rhs) {
		for (int r(0); r < rows; ++r) {
			for (int c(0); c < cols; ++c) {
				(*this)(r, c) = rhs(r, c);
			}
		}
		return (*this);
	}

	T det() const;
	Matrix<rows, cols, T> inv() const;


	Matrix<cols, rows, T> t() const;
};

template<int rows, int cols, typename T>
class Matrix<rows, cols, T> : public MatrixView<rows, cols, cols, T> {
	T vals[rows][cols];
	using SuperType = MatrixView<rows, cols, cols, T>;
public:


	constexpr int num_rows() {return rows;}
	constexpr int num_cols() {return cols;}

	Matrix() : MatrixView<rows, cols, cols, T>(&(vals[0][0])) {}

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

//	template<int rowStride>
	Matrix(Matrix<rows, cols, T> const& other) : SuperType(&(vals[0][0])) {
		for (int r(0); r < rows; ++r) {
			for (int c(0); c < cols; ++c) {
				(*this)(r, c) = other(r, c);
			}
		}
	}

	template<int rowStride>
	Matrix<rows, cols, T>& operator=(MatrixView<rows, cols, rowStride, T> const& other) {
		SuperType::operator=(other);
		return *this;
	}
};

template<int rows, int cols, int rowStride, typename T>
Matrix<cols, rows, T> MatrixView<rows, cols, rowStride, T>::t() const {
	Matrix<cols, rows, T> ret;
	for (int r(0); r < rows; ++r) {
		for (int c(0); c < cols; ++c) {
			ret(c, r) = (*this)(r, c);
		}
	}
	return ret;
}

template<int rowStride, typename T>
T det(MatrixView<1, 1, rowStride, T> const& mat) {
	return mat(0, 0);
}
template<int rowStride, typename T>
T det(MatrixView<2, 2, rowStride, T> const& mat) {
	return mat(0, 0)*mat(1, 1) - mat(0, 1)*mat(1, 0);
}
template<int rowStride, typename T>
T det(MatrixView<3, 3, rowStride, T> const& mat) {
	return (mat(0, 0)*mat(1, 1)*mat(2, 2) +
			mat(0, 1)*mat(1, 2)*mat(2, 0) +
			mat(0, 2)*mat(1, 0)*mat(2, 1)) -
			(mat(0, 2)*mat(1, 1)*mat(2, 0) +
			mat(0, 1)*mat(1, 0)*mat(2, 2) +
			mat(0, 0)*mat(1, 2)*mat(2, 1));
}
template<int rows, int cols, int rowStride, typename T>
T MatrixView<rows, cols, rowStride, T>::det() const {
	return SiLi::det(*this);
}

template<int rowStride, typename T>
Matrix<1, 1, T> inv(MatrixView<1, 1, rowStride, T> const& mat) {
	return Matrix<1, 1, T>(T(1) / mat(0, 0));
}
template<int rowStride, typename T>
Matrix<2, 2, T> inv(MatrixView<2, 2, rowStride, T> const& mat) {
	const T c = T(1) / mat.det();
	return c * Matrix<2, 2, T>({
		{mat(1, 1), -mat(0, 1)},
		{-mat(1, 0), mat(0, 0)}
	});
}
template<int rowStride, typename T>
Matrix<3, 3, T> inv(MatrixView<3, 3, rowStride, T> const& mat) {
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
template<int rows, int cols, int rowStride, typename T>
Matrix<rows, cols, T> MatrixView<rows, cols, rowStride, T>::inv() const {
	return SiLi::inv(*this);
}

template<int rows, int cols, int rs1, int rs2, typename T>
Matrix<rows, cols, T> operator+(MatrixView<rows, cols, rs1, T> const& lhs, MatrixView<rows, cols, rs2, T> const& rhs) {
	Matrix<rows, cols, T> ret;
	for (int row(0); row < rows; ++row) {
		for (int col(0); col < cols; ++col) {
			ret(row, col) = lhs(row, col) + rhs(row, col);
		}
	}
	return ret;
}
template<int rows, int cols, int rs1, typename T>
Matrix<rows, cols, T> operator+(MatrixView<rows, cols, rs1, T> const& lhs, T const& rhs) {
	Matrix<rows, cols, T> ret;
	for (int row(0); row < rows; ++row) {
		for (int col(0); col < cols; ++col) {
			ret(row, col) = lhs(row, col) + rhs;
		}
	}
	return ret;
}
template<int rows, int cols, int rs1, typename T>
Matrix<rows, cols, T> operator+(T const& lhs, MatrixView<rows, cols, rs1, T> const& rhs) {
	return rhs + lhs;
}

template<int rows, int cols, int rs1, int rs2, typename T>
Matrix<rows, cols, T> operator-(MatrixView<rows, cols, rs1, T> const& lhs, MatrixView<rows, cols, rs2, T> const& rhs) {
	Matrix<rows, cols, T> ret;
	for (int row(0); row < rows; ++row) {
		for (int col(0); col < cols; ++col) {
			ret(row, col) = lhs(row, col) - rhs(row, col);
		}
	}
	return ret;
}
template<int rows, int cols, int rs1, typename T>
Matrix<rows, cols, T> operator-(MatrixView<rows, cols, rs1, T> const& lhs, T const& rhs) {
	Matrix<rows, cols, T> ret;
	for (int row(0); row < rows; ++row) {
		for (int col(0); col < cols; ++col) {
			ret(row, col) = lhs(row, col) - rhs;
		}
	}
	return ret;
}
template<int rows, int cols, int rs1, typename T>
Matrix<rows, cols, T> operator-(T const& lhs, MatrixView<rows, cols, rs1, T> const& rhs) {
	return rhs - lhs;
}

template<int lrows, int mate, int rcols, int rs1, int rs2, typename T>
Matrix<lrows, rcols, T> operator*(MatrixView<lrows, mate, rs1, T> const& lhs,  MatrixView<mate, rcols, rs2, T> const& rhs) {
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
template<int rows, int cols, int rs, typename T, typename T2>
Matrix<rows, cols, T> operator*(MatrixView<rows, cols, rs, T> const& lhs, T const& rhs) {
	Matrix<rows, cols, T> ret;
	for (int oRow(0); oRow < rows; ++oRow) {
		for (int oCol(0); oCol < cols; ++oCol) {
			ret(oRow, oCol) = lhs(oRow, oCol) * rhs;
		}
	}
	return ret;
}
template<int rows, int cols, int rs, typename T, typename T2>
Matrix<rows, cols, T> operator*(T const& lhs, MatrixView<rows, cols, rs, T> const& rhs) {
	return rhs * lhs;
}

}


#ifdef IO_STREAM_PRESENT
#include <iostream>
template<int rows, int cols, int rowStride, typename T>
std::ostream& operator<< (std::ostream& stream, SiLi::MatrixView<rows, cols, rowStride, T> const& view) {
	for (int i(0); i < view.num_rows(); ++i) {
		for (int j(0); j < view.num_cols(); ++j) {
			stream << view(i, j) << "\t";
		}
		stream << "\n";
	}
	return stream;
}
#endif

