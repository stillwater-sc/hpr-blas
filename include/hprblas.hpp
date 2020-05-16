#pragma once
// hprblas.hpp :  include file containing templated C++ interfaces to High-Performance Reproducible BLAS routines
//
// Copyright (C) 2017-2020 Stillwater Supercomputing, Inc.
//
// This file is part of the HPRBLAS project, which is released under an MIT Open Source license.
#include <boost/numeric/mtl/mtl.hpp>

namespace sw {
namespace hprblas {

// LEVEL 1 BLAS operators

// 1-norm of a vector: sum of magnitudes of the vector elements, default increment stride is 1
template<typename Vector>
typename Vector::value_type asum(size_t n, const Vector& x, size_t incx = 1) {
	typename Vector::value_type sum = 0;
	size_t ix;
	for (ix = 0; ix < n; ix += incx) {
		sum += (x[ix] < 0 ? -x[ix] : x[ix]);
	}
	return sum;
}

// sum of the vector elements, default increment stride is 1
template<typename Vector>
typename Vector::value_type sum(size_t n, const Vector& x, size_t incx = 1) {
	typename Vector::value_type sum = 0;
	size_t ix;
	for (ix = 0; ix < n; ix += incx) {
		sum += x[ix];
	}
	return sum;
}

// a time x plus y
template<typename Scalar, typename Vector>
void axpy(size_t n, Scalar a, const Vector& x, size_t incx, Vector& y, size_t incy) {
	size_t cnt, ix, iy;
	for (cnt = 0, ix = 0, iy = 0; cnt < n && ix < x.size() && iy < y.size(); ++cnt, ix += incx, iy += incy) {
		y[iy] += a * x[ix];
	}
}

// vector copy
template<typename Vector>
void copy(size_t n, const Vector& x, size_t incx, Vector& y, size_t incy) {
	size_t cnt, ix, iy;
	for (cnt = 0, ix = 0, iy = 0; cnt < n && ix < x.size() && iy < y.size(); ++cnt, ix += incx, iy += incy) {
		y[iy] = x[ix];
	}
}

// dot product: the operator vector::x[index] is limited to uint32_t, so the arguments are limited to uint32_t as well
// The library does support arbitrary posit configuration conversions, but to simplify the 
// behavior of the dot product, the element type of the vectors x and y are declared to be the same.
// TODO: investigate if the vector<> index is always a 32bit entity?
template<typename Vector>
typename Vector::value_type dot(size_t n, const Vector& x, size_t incx, const Vector& y, size_t incy) {
	typename Vector::value_type product = 0;
	size_t cnt, ix, iy;
	for (cnt = 0, ix = 0, iy = 0; cnt < n && ix < x.size() && iy < y.size(); ++cnt, ix += incx, iy += incy) {
		product += x[ix] * y[iy];
	}
	return product;
}
// specialized dot product
template<typename Vector>
typename Vector::value_type dot(const Vector& x, const Vector& y) {
	typename Vector::value_type product = 0;
	size_t cnt, ix, iy;
	for (cnt = 0, ix = 0, iy = 0; cnt < size(x); ++cnt, ++ix, ++iy) {
		product += x[ix] * y[iy];
	}
	return product;
}
///
/// fused dot product operators

// Fused dot product with quire continuation
template<typename Quire, typename Vector>
void fdp_qr(Quire& sum_of_products, size_t n, const Vector& x, size_t incx, const Vector& y, size_t incy) {
	size_t ix, iy;
	for (ix = 0, iy = 0; ix < n && iy < n; ix = ix + incx, iy = iy + incy) {
		sum_of_products += sw::unum::quire_mul(x[ix], y[iy]);
	}
}
// Resolved fused dot product, with the option to control capacity bits in the quire
template<typename Vector, size_t capacity = 10>
typename Vector::value_type fdp_stride(size_t n, const Vector& x, size_t incx, const Vector& y, size_t incy) {
	constexpr size_t nbits = Vector::value_type::nbits;
	constexpr size_t es = Vector::value_type::es;
	sw::unum::quire<nbits, es, capacity> q(0);
	size_t ix, iy;
	for (ix = 0, iy = 0; ix < n && iy < n; ix = ix + incx, iy = iy + incy) {
		q += sw::unum::quire_mul(x[ix], y[iy]);
		if (sw::unum::_trace_quire_add) std::cout << q << '\n';
	}
	typename Vector::value_type sum;
	sw::unum::convert(q.to_value(), sum);     // one and only rounding step of the fused-dot product
	return sum;
}
// Specialized resolved fused dot product that assumes unit stride and a standard vector,
// with the option to control capacity bits in the quire
template<typename Vector, size_t capacity = 10>
typename Vector::value_type fdp(const Vector& x, const Vector& y) {
	constexpr size_t nbits = Vector::value_type::nbits;
	constexpr size_t es = Vector::value_type::es;
	sw::unum::quire<nbits, es, capacity> q(0);
	size_t ix, iy, n = size(x);
	for (ix = 0, iy = 0; ix < n && iy < n; ++ix, ++iy) {
		q += sw::unum::quire_mul(x[ix], y[iy]);
	}
	typename Vector::value_type sum;
	sw::unum::convert(q.to_value(), sum);     // one and only rounding step of the fused-dot product
	return sum;
}

// rotation of points in the plane
template<typename Rotation, typename Vector>
void rot(size_t n, Vector& x, size_t incx, Vector& y, size_t incy, Rotation c, Rotation s) {
	// x_i = c*x_i + s*y_i
	// y_i = c*y_i - s*x_i
	size_t cnt, ix, iy;
	for (cnt = 0, ix = 0, iy = 0; cnt < n && ix < x.size() && iy < y.size(); ++cnt, ix += incx, iy += incy) {
		Rotation x_i = c*x[ix] + s*y[iy];
		Rotation y_i = c*y[iy] - s*x[ix];
		y[iy] = y_i;
		x[ix] = x_i;
	}
}

// compute parameters for a Givens rotation
template<typename T>
void rotg(T& a, T& b, T& c, T&s) {
	// Given Cartesian coordinates (a,b) of a point, return parameters c,s,r, and z associated with the Givens rotation.
}

// scale a vector
template<typename Scalar, typename Vector>
void scale(size_t n, Scalar a, Vector& x, size_t incx) {
	size_t cnt, ix;
	for (cnt = 0, ix = 0; cnt < n && ix < x.size(); ix += incx) {
		x[ix] *= a;
	}
}

// swap two vectors
template<typename Vector>
void swap(size_t n, Vector& x, size_t incx, Vector& y, size_t incy) {
	size_t cnt, ix, iy;
	for (cnt = 0, ix = 0, iy = 0; cnt < n && ix < x.size() && iy < y.size(); ++cnt, ix += incx, iy += incy) {
		typename Vector::value_type tmp = x[ix];
		x[ix] = y[iy];
		y[iy] = tmp;
	}
}

// find the index of the element with maximum absolute value
template<typename Vector>
size_t amax(size_t n, const Vector& x, size_t incx) {
	typename Vector::value_type running_max = -INFINITY;
	size_t ix, index;
	for (ix = 0; ix < x.size(); ix += incx) {
		if (x[ix] > running_max) {
			index = ix;
			running_max = x[ix];
		}
	}
	return index;
}

// find the index of the element with minimum absolute value
template<typename Vector>
size_t amin(size_t n, const Vector& x, size_t incx) {
	typename Vector::value_type running_min = INFINITY;
	size_t ix, index;
	for (ix = 0; ix < x.size(); ix += incx) {
		if (x[ix] < running_min) {
			index = ix;
			running_min = x[ix];
		}
	}
	return index;
}

// absolute value of a complex number
template<typename T>
T cabs(T z) {
}

// print a vector
template<typename Vector>
void strided_print(std::ostream& ostr, size_t n, Vector& x, size_t incx = 1) {
	size_t cnt, ix;
	for (cnt = 0, ix = 0; cnt < n && ix < x.size(); ++cnt, ix += incx) {
		cnt == 0 ? ostr << "[" << x[ix] : ostr << ", " << x[ix];
	}
	ostr << "]";
}


// LEVEL 2 BLAS operators

// Matrix-vector product: b = A * x
template<typename Matrix, typename Vector>
void matvec(Vector& b, const Matrix& A, const Vector& x) {
	b = A * x;
}

// Matrix-vector product: b = A * x, posit specialized
template<size_t nbits, size_t es>
void matvec(mtl::vec::dense_vector< sw::unum::posit<nbits, es> >& b, const mtl::mat::dense2D< sw::unum::posit<nbits, es> >& A, const mtl::vec::dense_vector< sw::unum::posit<nbits, es> >& x) {
	// preconditions
	assert(A.num_cols() == size(x));
	assert(size(b) == size(x));

#if HPRBLAS_TRACE_ROUNDING_EVENTS
	unsigned errors = 0;
#endif
	size_t nr = size(b);
	size_t nc = size(x);
	for (size_t i = 0; i < nr; ++i) {
		sw::unum::quire<nbits, es> q(0);
		for (size_t j = 0; j < nc; ++j) {
			q += sw::unum::quire_mul(A[i][j], x[j]);
		}
		sw::unum::convert(q.to_value(), b[i]);     // one and only rounding step of the fused-dot product
#if HPRBLAS_TRACE_ROUNDING_EVENTS
		sw::unum::quire<nbits, es> qdiff = q;
		sw::unum::quire<nbits, es> qsum = b[i];
		qdiff -= qsum;
		if (!qdiff.iszero()) {
			++errors;
			std::cout << "q    : " << q << std::endl;
			std::cout << "qsum : " << qsum << std::endl;
			std::cout << "qdiff: " << qdiff << std::endl;
			sw::unum::posit<nbits, es> roundingError;
			convert(qdiff.to_value(), roundingError);
			std::cout << "matvec b[" << i << "] = " << posit_format(b[i]) << " rounding error: " << posit_format(roundingError) << " " << roundingError << std::endl;
		}
#endif
	}
#if HPRBLAS_TRACE_ROUNDING_EVENTS
	if (errors) {
		std::cout << "HPR-BLAS: tracing found " << errors << " rounding errors in matvec operation\n";
	}
#endif
}

// A times x = b fused matrix-vector product
template<size_t nbits, size_t es>
mtl::vec::dense_vector< sw::unum::posit<nbits, es> > fmv(const mtl::mat::dense2D< sw::unum::posit<nbits, es> >& A, const mtl::vec::dense_vector< sw::unum::posit<nbits, es> >& x) {
	// preconditions
	assert(A.num_cols() == size(x));
	mtl::vec::dense_vector< sw::unum::posit<nbits, es> > b(size(x));

#if HPRBLAS_TRACE_ROUNDING_EVENTS
	unsigned errors = 0;
#endif
	size_t nr = size(b);
	size_t nc = size(x);
	for (size_t i = 0; i < nr; ++i) {
		sw::unum::quire<nbits, es> q(0);
		for (size_t j = 0; j < nc; ++j) {
			q += sw::unum::quire_mul(A[i][j], x[j]);
		}
		sw::unum::convert(q.to_value(), b[i]);     // one and only rounding step of the fused-dot product
#if HPRBLAS_TRACE_ROUNDING_EVENTS
		sw::unum::quire<nbits, es> qdiff = q;
		sw::unum::quire<nbits, es> qsum = b[i];
		qdiff -= qsum;
		if (!qdiff.iszero()) {
			++errors;
			std::cout << "q    : " << q << std::endl;
			std::cout << "qsum : " << qsum << std::endl;
			std::cout << "qdiff: " << qdiff << std::endl;
			sw::unum::posit<nbits, es> roundingError;
			convert(qdiff.to_value(), roundingError);
			std::cout << "matvec b[" << i << "] = " << posit_format(b[i]) << " rounding error: " << posit_format(roundingError) << " " << roundingError << std::endl;
		}
#endif
	}
#if HPRBLAS_TRACE_ROUNDING_EVENTS
	if (errors) {
		std::cout << "HPR-BLAS: tracing found " << errors << " rounding errors in matvec operation\n";
	}
#endif
	return b;
}

// LEVEL 3 BLAS operators

template<typename Matrix>
void matmul(Matrix& C, const Matrix& A, const Matrix& B) {
	C = A * B;
}

// C = A * B fused matrix-matrix product when posits are used
template<size_t nbits, size_t es>
void matmul(mtl::mat::dense2D< sw::unum::posit<nbits, es> >& C, const mtl::mat::dense2D< sw::unum::posit<nbits, es> >& A, const mtl::mat::dense2D< sw::unum::posit<nbits, es> >& B) {
	// precondition
	assert(A.num_cols() == B.num_rows());
	size_t nr = A.num_rows();
	size_t nc = B.num_cols();
	size_t nk = A.num_cols();
	// TODO: add asserts to make certain that C is the right size

	for (size_t i = 0; i < nr; ++i) {
		for (size_t j = 0; j < nc; ++j) {
			sw::unum::quire<nbits, es> q(0);
			for (size_t k = 0; k < nk; ++k) {
				q += sw::unum::quire_mul(A[i][k], B[k][j]);
			}
			sw::unum::convert(q.to_value(), C[i][j]);     // one and only rounding step of the fused-dot product
		}
	}
}

// C = A * B fused matrix-matrix product when posits are used
template<size_t nbits, size_t es>
mtl::mat::dense2D< sw::unum::posit<nbits, es> > fmm(const mtl::mat::dense2D< sw::unum::posit<nbits, es> >& A, const mtl::mat::dense2D< sw::unum::posit<nbits, es> >& B) {
	// precondition
	assert(A.num_cols() == B.num_rows());
	size_t nr = A.num_rows();
	size_t nc = B.num_cols();
	size_t nk = A.num_cols();
	mtl::mat::dense2D< sw::unum::posit<nbits, es> > C(nr, nc);

	for (size_t i = 0; i < nr; ++i) {
		for (size_t j = 0; j < nc; ++j) {
			sw::unum::quire<nbits, es> q(0);
			for (size_t k = 0; k < nk; ++k) {
				q += sw::unum::quire_mul(A[i][k], B[k][j]);
			}
			sw::unum::convert(q.to_value(), C[i][j]);     // one and only rounding step of the fused-dot product
		}
	}
	return C;
}

template<typename Scalar>
inline Scalar minimum(const Scalar& a, const Scalar& b) {
	return (a < b ? a : b);
}

// if you only specify one set of block parameters in the specification then by 
// the virtue of doing block matrix multiplication you generate the invariant : blockHeight == blockWidth
// Thus, no need to specify blockHeight and blockWidth: we can simplify to blockSize
// blockHeight = blockWidth = blockSize

// subBlockMM generates the partial sums of a sub-block matrix multiply
// the QuireMatrix is [blockHeight][blockWidth] submatrix
// A and B matrices are full [n][m] and [m][n] matrices
template<typename Matrix>
void subBlockMM(Matrix& C_partial, const Matrix& A, unsigned Ai, unsigned Aj, const Matrix& B, unsigned Bi, unsigned Bj) {
	assert(mtl::mat::num_rows(C_partial) == mtl::mat::num_cols(C_partial));
	using Scalar = typename Matrix::value_type;
	constexpr size_t nbits = Scalar::nbits;
	constexpr size_t es = Scalar::es;

	unsigned aRows = unsigned(mtl::mat::num_rows(A));
	unsigned aCols = unsigned(mtl::mat::num_cols(A));
	unsigned bRows = unsigned(mtl::mat::num_rows(B));
	unsigned bCols = unsigned(mtl::mat::num_cols(B));

	unsigned blockSize = unsigned(mtl::mat::num_rows(C_partial));

	unsigned aRow = Ai * blockSize;
	unsigned aCol = Aj * blockSize;
	unsigned bRow = Bi * blockSize;
	unsigned bCol = Bj * blockSize;

	// calculate the shape of submatrix A
	unsigned ar = (aRow + blockSize < aRows) ? blockSize : aRows - aRow;
	unsigned ac = (aCol + blockSize < aCols) ? blockSize : aCols - aCol;
	unsigned br = (bRow + blockSize < bRows) ? blockSize : bRows - bRow;
	unsigned bc = (bCol + blockSize < bCols) ? blockSize : bCols - bCol;
	
//	std::cout << "A=" << ar << "x" << ac << "  B=" << br << "x" << bc << std::endl;
	assert(ac == br);
	// A(ar,ac) x B(br,bc) = C(ar,bc) if ac == br
	for (unsigned i = 0; i < ar; ++i) {
		for (unsigned j = 0; j < bc; ++j) {
			for (unsigned k = 0; k < ac; ++k) {
				C_partial[i][j] += A[aRow + i][aCol + k] * B[bRow + k][bCol + j];
			}
//			std::cout << "A(" << Ai << "," << Aj << ") B(" << Bi << "," << Bj << ") C(" << i << "," << j << ")\n";
//			printMatrix(std::cout, "partial", C_partial);
		}
	}
}

// copySubBlockInto copies a subblock matrix into of the mother matrix
template<typename Matrix>
void copySubBlockInto(Matrix& A, unsigned ai, unsigned aj, const Matrix& SubBlock) {
	unsigned blockHeight = unsigned(mtl::mat::num_rows(SubBlock));
	unsigned blockWidth = unsigned(mtl::mat::num_cols(SubBlock));

	unsigned aRows = unsigned(mtl::mat::num_rows(A));
	unsigned aCols = unsigned(mtl::mat::num_cols(A));

	unsigned aRow = ai * blockHeight;
	unsigned aCol = aj * blockWidth;

	unsigned maxRow = (aRow + blockHeight < aRows) ? blockHeight : aRows - aRow;
	unsigned maxCol = (aCol + blockWidth < aCols) ? blockWidth : aCols - aCol;

	for (unsigned i = 0; i < maxRow; ++i) {
		for (unsigned j = 0; j < maxCol; ++j) {
			A[aRow + i][aCol + j] = SubBlock[i][j];
		}
	}
}

// bmm is a blocked matrix multply: C = A * B
template<typename Matrix>
Matrix bmm(const Matrix& A, const Matrix& B, unsigned blockSize) {
	// precondition
	assert(A.num_cols() == B.num_rows());
	unsigned nr = unsigned(A.num_rows());
	unsigned nc = unsigned(B.num_cols());
	unsigned nk = unsigned(A.num_cols());
	Matrix C(nr, nc);

	using Scalar = typename Matrix::value_type;
	Matrix C_partial(blockSize, blockSize);

	unsigned nrRowBlocks = nr % blockSize ? nr / blockSize + 1 : nr / blockSize;
	unsigned nrColBlocks = nc % blockSize ? nr / blockSize + 1 : nc / blockSize;
	for (unsigned bi = 0; bi < nrRowBlocks; ++bi) {			// row block index
		for (unsigned bj = 0; bj < nrColBlocks; ++bj) {		// col block index
			C_partial = Scalar(0);
			for (unsigned bk = 0; bk < nrRowBlocks; ++bk) { // block iterator
				subBlockMM(C_partial, A, bi, bk, B, bk, bj);
			}
			copySubBlockInto(C, bi, bj, C_partial);
		}
	}
	return C;
}

// subBlockMM generates the partial sums of a sub-block matrix multiply
// the QuireMatrix is a square matrix
// A and B matrices are full [n][m] and [m][n] matrices
template<typename QuireMatrix, typename Matrix>
void subBlockMM(QuireMatrix& C, const Matrix& A, unsigned Ai, unsigned Aj, const Matrix& B, unsigned Bi, unsigned Bj) {
	assert(mtl::mat::num_rows(C) == mtl::mat::num_cols(C));
	using Scalar = typename Matrix::value_type;
	constexpr size_t nbits = Scalar::nbits;
	constexpr size_t es = Scalar::es;

	unsigned aRows = unsigned(mtl::mat::num_rows(A));
	unsigned aCols = unsigned(mtl::mat::num_cols(A));
	unsigned bRows = unsigned(mtl::mat::num_rows(B));
	unsigned bCols = unsigned(mtl::mat::num_cols(B));

	unsigned blockSize = unsigned(mtl::mat::num_rows(C));

	unsigned aRow = Ai * blockSize;
	unsigned aCol = Aj * blockSize;
	unsigned bRow = Bi * blockSize;
	unsigned bCol = Bj * blockSize;

	// calculate the shape of submatrix A
	unsigned ar = (aRow + blockSize < aRows) ? blockSize : aRows - aRow;
	unsigned ac = (aCol + blockSize < aCols) ? blockSize : aCols - aCol;
	unsigned br = (bRow + blockSize < bRows) ? blockSize : bRows - bRow;
	unsigned bc = (bCol + blockSize < bCols) ? blockSize : bCols - bCol;

	//	std::cout << "A=" << ar << "x" << ac << "  B=" << br << "x" << bc << std::endl;
	assert(ac == br);
	// A(ar,ac) x B(br,bc) = C(ar,bc) if ac == br
	for (unsigned i = 0; i < ar; ++i) {
		for (unsigned j = 0; j < bc; ++j) {
			for (unsigned k = 0; k < ac; ++k) {
				C[i][j] += sw::unum::quire_mul<nbits,es>(A[aRow+i][aCol+k], B[bRow+k][bCol+j]);
			}
		}
	}
}

// subBlockRound takes a sub-block address and a QuireMatrix and rounds the partial sums

template<typename Matrix, typename QuireMatrix>
void subBlockRound(Matrix& C, unsigned ci, unsigned cj, const QuireMatrix& C_partial) {
	unsigned blockHeight = unsigned(mtl::mat::num_rows(C_partial));
	unsigned blockWidth = unsigned(mtl::mat::num_cols(C_partial));

	unsigned cRows = unsigned(mtl::mat::num_rows(C));
	unsigned cCols = unsigned(mtl::mat::num_cols(C));

	unsigned cRow = ci * blockHeight;
	unsigned cCol = cj * blockWidth;

	unsigned maxRow = (cRow + blockHeight < cRows) ? blockHeight : cRows - cRow;
	unsigned maxCol = (cCol + blockWidth < cCols) ? blockWidth : cCols - cCol;

	for (unsigned i = 0; i < maxRow; ++i) {
		for (unsigned j = 0; j < maxCol; ++j) {
			sw::unum::convert(C_partial[i][j].to_value(), C[cRow + i][cCol + j]);
		}
	}
}

// copySubBlock copies a subblock matrix out of the mother matrix
// This function is more generic than used in block matmul, as this can take non-square matrices
template<typename Matrix>
void copySubBlock(Matrix& SubBlock, const Matrix& M, unsigned bi, unsigned bj) {
	unsigned blockHeight = unsigned(mtl::mat::num_rows(SubBlock));
	unsigned blockWidth = unsigned(mtl::mat::num_cols(SubBlock));
	for (unsigned i = 0; i < blockHeight; ++i) {
		for (unsigned j = 0; j < blockWidth; ++j) {
			SubBlock[i][j] = M[bi*blockHeight + i][bj*blockWidth + j];
		}
	}
}

// bfmm is a blocked fused matrix multply: C = A * B
template<typename Matrix>
Matrix bfmm(const Matrix& A, const Matrix& B, unsigned blockSize) {
	// precondition
	assert(A.num_cols() == B.num_rows());
	unsigned nr = unsigned(A.num_rows());
	unsigned nc = unsigned(B.num_cols());
	unsigned nk = unsigned(A.num_cols());
	Matrix C(nr, nc);

	using Scalar = typename Matrix::value_type;
	constexpr size_t nbits = Scalar::nbits;
	constexpr size_t es = Scalar::es;
	using Quire = typename sw::unum::quire<nbits, es>;
	using QuireMatrix = typename mtl::mat::dense2D<Quire>;
	QuireMatrix C_partial(blockSize, blockSize);

	unsigned nrRowBlocks = nr % blockSize ? nr / blockSize + 1 : nr / blockSize;
	unsigned nrColBlocks = nc % blockSize ? nr / blockSize + 1  : nc / blockSize;
	for (unsigned bi = 0; bi < nrRowBlocks; ++bi) {			// row block index
		for (unsigned bj = 0; bj < nrColBlocks; ++bj) {		// col block index
			C_partial = Scalar(0);
			for (unsigned bk = 0; bk < nrRowBlocks; ++bk) { // block iterator
				subBlockMM(C_partial, A, bi, bk, B, bk, bj);
			}
			subBlockRound(C, bi, bj, C_partial);  // C_sub(i,j) = round(QuireMatrix)
		}
	}
	return C;
}

} // namespace hprblas
} // namespace sw
