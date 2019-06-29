#pragma once
// blas.hpp :  include file containing templated C++ interfaces to BLAS routines
//
// Copyright (C) 2017 Stillwater Supercomputing, Inc.
//
// This file is part of the universal numbers project, which is released under an MIT Open Source license.

#include <vector>

namespace sw {
	namespace hprblas {

		// LEVEL 1 BLAS operators

		// sum of magnitudes of the vector elements
		template<typename Vector>
		typename Vector::value_type asum(size_t n, const Vector& x, size_t incx) {
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
			sw::unum::quire<nbits, es, capacity> q = 0;
			size_t ix, iy;
			for (ix = 0, iy = 0; ix < n && iy < n; ix = ix + incx, iy = iy + incy) {
				q += sw::unum::quire_mul(x[ix], y[iy]);
				if (sw::unum::_trace_quire_add) std::cout << q << '\n';
			}
			typename Vector::value_type sum;
			convert(q.to_value(), sum);     // one and only rounding step of the fused-dot product
			return sum;
		}
		// Specialized resolved fused dot product that assumes unit stride and a standard vector,
		// with the option to control capacity bits in the quire
		template<typename Vector, size_t capacity = 10>
		typename Vector::value_type fdp(const Vector& x, const Vector& y) {
			constexpr size_t nbits = Vector::value_type::nbits;
			constexpr size_t es = Vector::value_type::es;
			sw::unum::quire<nbits, es, capacity> q = 0;
			size_t ix, iy, n = size(x);
			for (ix = 0, iy = 0; ix < n && iy < n; ++ix, ++iy) {
				q += sw::unum::quire_mul(x[ix], y[iy]);
			}
			typename Vector::value_type sum;
			convert(q.to_value(), sum);     // one and only rounding step of the fused-dot product
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
		void print(std::ostream& ostr, size_t n, Vector& x, size_t incx = 1) {
			size_t cnt, ix;
			for (cnt = 0, ix = 0; cnt < n && ix < x.size(); ++cnt, ix += incx) {
				cnt == 0 ? ostr << "[" << x[ix] : ostr << ", " << x[ix];
			}
			ostr << "]";
		}


		// LEVEL 2 BLAS operators

		// A times x = b matrix-vector product
		template<typename Matrix, typename Vector>
		void matvec(const Matrix& A, const Vector& x, Vector& b) {
			b = A * x;
		}

		// A times x = b fused matrix-vector product
		template<size_t nbits, size_t es>
		void matvec(const mtl::mat::dense2D< sw::unum::posit<nbits, es> >& A, const mtl::vec::dense_vector< sw::unum::posit<nbits, es> >& x, mtl::vec::dense_vector< sw::unum::posit<nbits, es> >& b) {
			// preconditions
			assert(A.num_cols() == size(x));
			assert(size(b) == A.num_rows());
			size_t nr = size(b);
			size_t nc = size(x);
			for (size_t i = 0; i < nr; ++i) {
				sw::unum::quire<nbits, es> q = 0;
				for (size_t j = 0; j < nc; ++j) {
					q += sw::unum::quire_mul(A[i][j], x[j]);
				}
				convert(q.to_value(), b[i]);     // one and only rounding step of the fused-dot product
			}
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
			assert(C.num_rows() == A.num_rows());
			assert(C.num_cols() == B.num_cols());
			size_t nr = C.num_rows();
			size_t nc = C.num_cols();
			size_t nk = A.num_cols();
			for (size_t i = 0; i < nr; ++i) {
				for (size_t j = 0; j < nc; ++j) {
					sw::unum::quire<nbits, es> q = 0;
					for (size_t k = 0; k < nk; ++k) {
						q += sw::unum::quire_mul(A[i][k], B[k][j]);
					}
					convert(q.to_value(), C[i][j]);     // one and only rounding step of the fused-dot product
				}
			}
		}

	} // namespace blas

} // namespace sw
