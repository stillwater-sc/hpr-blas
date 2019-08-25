#pragma once
// cholesky.hpp implementation of the Cholesky solver
//
// Copyright (C) 2017-2018 Stillwater Supercomputing, Inc.
//
// This file is part of the universal numbers project, which is released under an MIT Open Source license.

#include <vector>
#define POSIT_VERBOSE_OUTPUT
#define QUIRE_TRACE_ADD
#include <universal/posit/posit>


namespace sw {
	namespace hprblas {

// Cholesky requires the matrix to be symmetric positive-definite
template<typename Ty>
void Cholesky(std::vector<Ty>& S, std::vector<Ty>& D) {
	size_t d = size_t(std::sqrt(S.size()));
	assert(S.size() == d*d);
	assert(D.size() == d*d);
	for (size_t k = 0; k<d; ++k) {
		Ty sum = 0.;
		for (size_t p = 0; p<k; ++p)sum += D[k*d + p] * D[k*d + p];
		D[k*d + k] = sqrt(S[k*d + k] - sum);
		for (size_t i = k + 1; i<d; ++i) {
			Ty sum = 0.;
			for (size_t p = 0; p<k; ++p)sum += D[i*d + p] * D[k*d + p];
			D[i*d + k] = (S[i*d + k] - sum) / D[k*d + k];
		}
	}
}
// This version could be more efficient on some architectures
// Use solveCholesky for both Cholesky decompositions
template<typename Ty>
void CholeskyRow(std::vector<Ty>& S, std::vector<Ty>& D) {
	size_t d = size_t(std::sqrt(S.size()));
	assert(S.size() == d*d);
	assert(D.size() == d*d);
	for (size_t k = 0; k<d; ++k) {
		for (size_t j = 0; j<d; ++j) {
			Ty sum = 0.;
			for (size_t p = 0; p<j; ++p) sum += D[k*d + p] * D[j*d + p];
			D[k*d + j] = (S[k*d + j] - sum) / D[j*d + j];
		}
		Ty sum = 0.;
		for (size_t p = 0; p<k; ++p) sum += D[k*d + p] * D[k*d + p];
		D[k*d + k] = sqrt(S[k*d + k] - sum);
	}
}
// SolveCholesky takes an LU decomposition, LU, and a right hand side vector, b, and produces a result, x.
template<typename Ty>
void SolveCholesky(const std::vector<Ty>& LU, const std::vector<Ty>& b, std::vector<Ty>& x) {
	int d = (int)b.size();
	std::vector<Ty> y(d);
	for (int i = 0; i<d; ++i) {
		Ty sum = 0.;
		for (int k = 0; k<i; ++k) sum += LU[i*d + k] * y[k];
		y[i] = (b[i] - sum) / LU[i*d + i];
	}
	for (int i = d - 1; i >= 0; --i) {
		Ty sum = 0.;
		for (int k = i + 1; k<d; ++k) sum += LU[k*d + i] * x[k];
		x[i] = (y[i] - sum) / LU[i*d + i];
	}
}

	}  // namespace hprblas
} // namespace sw
