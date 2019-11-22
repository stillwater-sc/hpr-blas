#pragma once
// norms.hpp : vector and matrix norms
//
// Copyright (C) 2017-2019 Stillwater Supercomputing, Inc.
//
// This file is part of the HPR-BLAS project, which is released under an MIT Open Source license.

/*
a general vector norm | x | , sometimes written with a double bar as || x ||, 
is a nonnegative norm defined such that

1. | x | > 0 when x != 0 and | x |= 0 iff x = 0.

2. | kx |= | k || x | for any scalar k.

3. | x + y |<= |x| + |y|
*/

namespace sw {
namespace hprblas {

// L1-norm = sum of absolute value of each vector element
template<typename Vector>
typename Vector::value_type l1_norm(const Vector& v) {
	using Scalar = typename Vector::value_type;

	Scalar l1 = Scalar(0);
	for (unsigned i = 0; i < size(v); ++i) {
		l1 += abs(v[i]);
	}
	return l1;
}

// L2-norm = Manhattan distance
template<typename Vector>
typename Vector::value_type l2_norm(const Vector& v) {
	using Scalar = typename Vector::value_type;

	Scalar l2 = Scalar(0);
	for (unsigned i = 0; i < size(v); ++i) {
		l2 += v[i] * v[i];
	}
	return sqrt(l2);
}

// Linfinity-norm = max of the absolute value of each vector element
template<typename Vector>
typename Vector::value_type linf_norm(const Vector& v) {
	using Scalar = typename Vector::value_type;

	Scalar linf = Scalar(0);
	for (unsigned i = 0; i < size(v); ++i) {
		linf = abs(v[i]) > linf ? abs(v[i]) : linf;
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
