#pragma once
// gauss_jordan.hpp: Gauss-Jordan matrix inversion
//
// Copyright (C) 2017-2019 Stillwater Supercomputing, Inc.
//
// This file is part of the HPR-BLAS project, which is released under an MIT Open Source license.

#include <boost/numeric/mtl/mtl.hpp>

namespace sw {
namespace hprblas {

template<typename Matrix>
Matrix GaussJordanInversion(const Matrix& A) {
	using Scalar = typename mtl::Collection<Matrix>::value_type;

	size_t m = A.num_rows();
	size_t n = A.num_cols();
	Matrix a(n, m), inv(m, n);
	a = A; // you need a deep copy
	inv = Scalar(1);

	// Performing elementary operations 
	for (unsigned i = 0; i < m; ++i) {
		if (a[i][i] == 0) {
			unsigned c = 1;
			while (a[i + c][i] == 0 && (i + c) < n)	++c;
			if ((i + c) == n) break;

			for (unsigned j = i, k = 0; k < n; ++k) {
				std::swap(a[j][k], a[j + c][k]);
				std::cerr << "TBD" << std::endl; // need to create a permutation matrix
			}
		}
		// transform to diagonal matrix
		for (unsigned j = 0; j < m; j++) {
			if (i != j) {
				Scalar scale = a[j][i] / a[i][i];
				for (unsigned k = 0; k < n; ++k) {
					a[j][k] = a[j][k] - a[i][k] * scale;
					inv[j][k] = inv[j][k] - inv[i][k] * scale;
				}
			}
			//std::cout << i << "," << j << std::endl;
			//sw::hprblas::printMatrix(std::cout, "a", a);
			//sw::hprblas::printMatrix(std::cout, "inv", inv);
		}
	}
	// transform to identity matrix
	for (unsigned i = 0; i < m; ++i) {
		Scalar normalize = a[i][i];
		a[i][i] = Scalar(1);
		for (unsigned j = 0; j < n; ++j) {
			inv[i][j] /= normalize;
		}
	}
	//printMatrix(std::cout, "conversion", a);
	return inv;
}

} // namespace hprblas
} // namespace sw