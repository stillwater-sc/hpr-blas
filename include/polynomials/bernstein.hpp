// bernstein.hpp: evaluation of Bernstein polynomials
//
// Copyright (C) 2017-2020 Stillwater Supercomputing, Inc.
// Author: Allan Leal
//
// This file is part of the HPR-BLAS project, which is released under an MIT Open Source license.
#include <universal/functions/binomial.hpp>

// sign returns (-1)^exponent
int sign(int exponent) {
    return ( (exponent % 2) == 0 ) ? 1 : -1;
}

/*
The Bernstein matrix of order N is an NxN matrix A which can be used to
transform a vector of power basis coefficients C representing a polynomial 
P(X) to a corresponding Bernstein basis coefficient vector B:
      B = A * C
 */

// bernstein returns the Bernstein matrix
template<typename Matrix>
Matrix bernstein(int n) {
	using Scalar = typename Matrix::value_type;

	if (n < 0) {
		std::cerr << "Degree of the Bernstein matrix can't be negative: " << n << std::endl;
		return Matrix{};
	}
	Matrix A(n, n);
	for (int j = 0; j < n; ++j ) {
		for (int i = 0; i <= j; ++i ) {
			A[i][j] = sign(j - i) * 
				sw::function::binomial(n - 1 - i, j - i) * 
				sw::function::binomial(n - 1, i);
		}
	}
	return A;
}
