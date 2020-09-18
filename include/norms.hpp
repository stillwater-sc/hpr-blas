#pragma once
// norms.hpp : vector and matrix norms using the quire
//
// Copyright (C) 2017-2020 Stillwater Supercomputing, Inc.
//
// This file is part of the HPRBLAS project, which is released under an MIT Open Source license.

/*
a general vector norm | x | , sometimes written with a double bar as || x ||, 
is a nonnegative norm defined such that

1. | x | > 0 when x != 0 and | x |= 0 iff x = 0.

2. | kx |= | k || x | for any scalar k.

3. | x + y |<= |x| + |y|
*/

namespace sw {
namespace hprblas {

// L1-norm = sum of absolute value of each vector element (Manhattan Distance)
template<typename Vector>
typename mtl::traits::enable_if_vector<Vector, typename mtl::RealMagnitude<typename mtl::Collection<Vector>::value_type>::type>::type
inline l1_norm(const Vector& v) {
	using Scalar = typename Vector::value_type;

	Scalar l1 = Scalar(0);
	for (unsigned i = 0; i < size(v); ++i) {
		l1 += abs(v[i]);
	}
	return l1;
}

template<size_t nbits, size_t es>
sw::unum::posit<nbits, es> l1_norm(const mtl::dense_vector<sw::unum::posit<nbits, es> > & v) {
	using Scalar = sw::unum::posit<nbits, es>;
	sw::unum::quire<nbits, es> q(0);
	for (unsigned i = 0; i < size(v); ++i) {
		q += abs(v[i]);
	}
	Scalar l1;
	sw::unum::convert(q.to_value(), l1);// one and only rounding step of the l1-norm
	return l1;
}

template<size_t nbits, size_t es>
sw::unum::posit<nbits, es> l1_norm(const mtl::dense2D<sw::unum::posit<nbits, es> > & M) {
	using Scalar = sw::unum::posit<nbits, es>;
	sw::unum::quire<nbits, es> q(0);
	for (unsigned i = 0; i < mtl::mat::num_rows(M); ++i) {
		for (unsigned j = 0; j < mtl::mat::num_cols(M); ++j) {
			q += abs(M[i][j]);
		}
	}
	Scalar l1;
	sw::unum::convert(q.to_value(), l1);// one and only rounding step of the l1-norm
	return l1;
}

// L2-norm = Euclidean distance
template<typename Vector>
typename Vector::value_type l2_norm(const Vector& v) {
	using Scalar = typename Vector::value_type;

	Scalar l2 = Scalar(0);
	for (unsigned i = 0; i < size(v); ++i) {
		l2 += v[i] * v[i];
	}
	return sqrt(l2);
}

// L2-norm = Euclidean distance, posit specialized
template<size_t nbits, size_t es>
sw::unum::posit<nbits, es> l2_norm(const mtl::vec::dense_vector<sw::unum::posit<nbits, es> >& v) {
	using Scalar = sw::unum::posit<nbits,es>;
	sw::unum::quire<nbits, es> q(0);
	for (unsigned i = 0; i < size(v); ++i) {
		q += sw::unum::quire_mul(v[i], v[i]);
	}
	Scalar l2 = Scalar(0);
	convert(q.to_value(), l2);     // first rounding step of the l2-norm
	return sqrt(l2);               // second rounding step of the l2-norm
}

// Linfinity-norm = max of the absolute value of each vector element
template<typename Vector>
typename mtl::traits::enable_if_vector<Vector, typename mtl::RealMagnitude<typename mtl::Collection<Vector>::value_type>::type>::type
inline linf_norm(const Vector& v) {
	using Scalar = typename Vector::value_type;
	Scalar linf = Scalar(0);
	for (unsigned i = 0; i < size(v); ++i) {
		linf = abs(v[i]) > linf ? abs(v[i]) : linf;
	}
	return linf;
}

// Linfinity-norm = max of the absolute value of each matrix element
template<typename Matrix>
typename mtl::traits::enable_if_matrix<Matrix, typename mtl::RealMagnitude<typename mtl::Collection<Matrix>::value_type>::type>::type
inline linf_norm(const Matrix& M) {
	using Scalar = typename Matrix::value_type;
	Scalar linf = Scalar(0);
	for (unsigned i = 0; i < mtl::mat::num_rows(M); ++i) {
		for (unsigned j = 0; j < mtl::mat::num_cols(M); ++j) {
			linf = abs(M[i][j]) > linf ? abs(M[i][j]) : linf;
		}
	}
	return linf;
}

// Frobenius-norm = sqrt of the sum of absolute squares
template<typename Matrix>
typename Matrix::value_type frobenius_norm(const Matrix& M) {
	using Scalar = typename Matrix::value_type;
	assert(mtl::mat::num_rows(M) == mtl::mat::num_cols(M)); // assuming squareness
	int N = int(mtl::mat::num_cols(M));
	Scalar frobenius = Scalar(0);
	for (int i = 0; i < N; ++i) {
		for (int j = 0; j < N; ++j) {
			frobenius += M[i][j] * M[i][j];
		}
	}
	return sqrt(frobenius);
}

// Frobenius-norm = sqrt of the sum of absolute squares, posit specialized for dense matrices
template<size_t nbits, size_t es>
sw::unum::posit<nbits, es> frobenius_norm(const mtl::dense2D<sw::unum::posit<nbits, es> >& M) {
	assert(mtl::mat::num_rows(M) == mtl::mat::num_cols(M)); // assuming squareness
	int N = int(mtl::mat::num_cols(M));
	sw::unum::quire<nbits, es> q(0);
	for (int i = 0; i < N; ++i) {
		for (int j = 0; j < N; ++j) {
			q += sw::unum::quire_mul(M[i][j], M[i][j]);
		}
	}
	using Scalar = sw::unum::posit<nbits, es>;
	Scalar frobenius = Scalar(0);
	convert(q.to_value(), frobenius);     // first rounding step of the Frobenius-norm
	return sqrt(frobenius);               // second rounding step of the Frobenius-norm
}

// Calculate the volume of the bounding box that contains the absolute error given L-infinity norm
template<typename Scalar>
Scalar error_volume(Scalar Linfinity, unsigned dimensionality, bool measuredInULP = false) {
	Scalar ulp = std::numeric_limits<Scalar>::epsilon();
	if (measuredInULP) {
		return sw::unum::integer_power(Linfinity / ulp, dimensionality);
	}
	return sw::unum::integer_power(Linfinity, dimensionality);
}

} // namespace hprblas
} // namespace sw
