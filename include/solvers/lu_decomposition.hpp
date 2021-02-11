#pragma once
// lu_decomposition.hpp: LU decomposition algorithms and solvers
//
// Copyright (C) 2017-2021 Stillwater Supercomputing, Inc.
//
// This file is part of the universal numbers project, which is released under an MIT Open Source license.

#include <vector>
#define POSIT_VERBOSE_OUTPUT
#define QUIRE_TRACE_ADD
#include <universal/number/posit/posit>

namespace sw {
namespace hprblas {


// The following compact LU factorization schemes are described
// in Dahlquist, Bjorck, Anderson 1974 "Numerical Methods".
//
// S and D are d by d matrices.  However, they are stored in
// memory as 1D arrays of length d*d.  Array indices are in
// the C style such that the first element of array A is A[0]
// and the last element is A[d*d-1].
//
// These routines are written with separate source S and
// destination D matrices so the source matrix can be retained
// if desired.  However, the compact schemes were designed to
// perform in-place computations to save memory.  In
// other words, S and D can be the SAME matrix.  

// Crout implements an in-place LU decomposition, that is, S and D can be the same
// Crout uses unit diagonals for the upper triangle


/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Crout using MTL data structures

// Crout method using MTL data structures
template<typename Matrix>
void Crout(const Matrix& S, Matrix& D) {
	assert(num_rows(S) == num_rows(D));
	assert(num_cols(S) == num_cols(D));
	using value_type = typename mtl::Collection<Matrix>::value_type;
	size_t N = num_rows(S);
	for (size_t k = 0; k < N; ++k) {
		for (size_t i = k; i<N; ++i) {
			value_type sum = 0.;
			for (size_t p = 0; p<k; ++p) sum += D[i][p] * D[p][k];
			D[i][k] = S[i][k] - sum; // not dividing by diagonals
		}
		for (size_t j = k + 1; j < N; ++j) {
			value_type sum = 0.;
			for (int p = 0; p<k; ++p) sum += D[k][p] * D[p][j];
			D[k][j] = (S[k][j] - sum) / D[k][k];
		}
	}
}


// SolveCrout takes an LU decomposition, LU, and a right hand side vector, b, and produces a result, x.
template<typename Matrix, typename Vector>
void SolveCrout(const Matrix& LU, const Vector& b, Vector& x) {
	assert(num_cols(LU) == size(b));
	size_t N = size(b);
	using value_type = typename mtl::Collection<Matrix>::value_type;
	mtl::dense_vector<value_type> y(N);
	for (int i = 0; i < N; ++i) {
		value_type sum = 0.0;
		for (int k = 0; k < i; ++k) sum += LU[i][k] * y[k];
		y[i] = (b[i] - sum) / LU[i][i];

	}
	for (long i = long(N) - 1; i >= 0; --i) {
		value_type sum = 0.0;
		for (int k = i + 1; k < N; ++k) {
			//cout << "lu[] = " << LU[i][k] << " x[" << k << "] = " << x[k] << endl;
			sum += LU[i][k] * x[k];
		}
		//cout << "sum " << sum << endl;
		x[i] = (y[i] - sum); // not dividing by diagonals
	}
}


///////////////////////////////////////////////////////////////////////////////////
/// CroutFDP with MTL data structures

template<size_t nbits, size_t es, size_t capacity = 10>
void CroutFDP(mtl::dense2D< sw::universal::posit<nbits, es> >& S, mtl::dense2D< sw::universal::posit<nbits, es> >& D) {
	size_t d = num_rows(S);
	assert(size(S) == size(D));
	using namespace sw::unum;
	for (int k = 0; k < d; ++k) {
		for (int i = k; i < d; ++i) {
			quire<nbits, es, capacity> q;
			q.reset();
			//for (int p = 0; p < k; ++p) q += D[i][p] * D[p][k];   if we had expression templates for the quire
			for (int p = 0; p < k; ++p) q += quire_mul(D[i][p], D[p][k]);
			posit<nbits, es> sum;
			convert(q.to_value(), sum);     // one and only rounding step of the fused-dot product
			// TODO: can we add the difference to the quire operation?
			D[i][k] = S[i][k] - sum; // not dividing by diagonals

#if HPRBLAS_TRACE_ROUNDING_EVENTS
			quire<nbits, es, capacity> qsum(sum);
			q -= qsum;
			if (!q.iszero()) {
				sw::unum::posit<nbits, es> roundingError;
				convert(q.to_value(), roundingError);
				std::cout << "D[" << i << "," << k << "] rounding error: " << roundingError << std::endl;
			}
#endif
		}
		for (int j = k + 1; j < d; ++j) {
			quire<nbits, es, capacity> q;
			q.reset();
			//for (int p = 0; p < k; ++p) q += D[k][p] * D[p][j];   if we had expression templates for the quire
			for (int p = 0; p < k; ++p) q += quire_mul(D[k][p], D[p][j]);
			posit<nbits, es> sum;
			convert(q.to_value(), sum);   // one and only rounding step of the fused-dot product
			D[k][j] = (S[k][j] - sum) / D[k][k];

#if HPRBLAS_TRACE_ROUNDING_EVENTS
			quire<nbits, es, capacity> qsum(sum);
			q -= qsum;
			if (!q.iszero()) {
				sw::unum::posit<nbits, es> roundingError;
				convert(q.to_value(), roundingError);
				std::cout << "D[" << k << "," << j << "] rounding error: " << roundingError << std::endl;
			}
#endif

		}
	}
}

// SolveCrout takes an LU decomposition, LU, and a right hand side vector, b, and produces a result, x.
template<size_t nbits, size_t es, size_t capacity = 10>
void SolveCroutFDP(const mtl::dense2D< sw::universal::posit<nbits, es> >& LU, const mtl::dense_vector< sw::universal::posit<nbits, es> >& b, mtl::dense_vector< sw::universal::posit<nbits, es> >& x)
{
	using namespace sw::universal;

	size_t d = size(b);
	std::vector< posit<nbits, es> > y(d);
	for (long i = 0; i < d; ++i) {
		quire<nbits, es, capacity> q;
		// for (int k = 0; k < i; ++k) q += LU[i][k] * y[k];   if we had expression templates for the quire
		for (int k = 0; k < i; ++k) q += quire_mul(LU[i][k], y[k]);
		posit<nbits, es> sum;
		convert(q.to_value(), sum);   // one and only rounding step of the fused-dot product
		y[i] = (b[i] - sum) / LU[i][i];
	}
	for (long i = long(d) - 1; i >= 0; --i) {
		quire<nbits, es, capacity> q;
		// for (int k = i + 1; k < d; ++k) q += LU[i][k] * x[k];   if we had expression templates for the quire
		for (int k = i + 1; k < d; ++k) {
			//cout << "lu[] = " << LU[i][k] << " x[" << k << "] = " << x[k] << endl;
			q += quire_mul(LU[i][k], x[k]);
		}
		posit<nbits, es> sum;
		convert(q.to_value(), sum);  // one and only rounding step of the fused-dot product
									 //cout << "sum " << sum << endl;
		x[i] = (y[i] - sum); // not dividing by diagonals
	}
}


#if 0
/// DEPRECATED
// Crout generates an LU decomposition
template<typename Ty>
void Crout(std::vector<Ty>& S, std::vector<Ty>& D) {
	size_t d = size_t(std::sqrt(S.size()));
	assert(S.size() == d * d);
	assert(D.size() == d * d);
	for (size_t k = 0; k < d; ++k) {
		for (size_t i = k; i < d; ++i) {
			Ty sum = 0.;
			for (size_t p = 0; p < k; ++p) sum += D[i*d + p] * D[p*d + k];
			D[i*d + k] = S[i*d + k] - sum; // not dividing by diagonals
		}
		for (size_t j = k + 1; j < d; ++j) {
			Ty sum = 0.;
			for (int p = 0; p < k; ++p) sum += D[k*d + p] * D[p*d + j];
			D[k*d + j] = (S[k*d + j] - sum) / D[k*d + k];
		}
	}
}

/// DEPRECATED
// SolveCrout takes an LU decomposition, LU, and a right hand side vector, b, and produces a result, x.
template<typename Ty>
void SolveCrout(const std::vector<Ty>& LU, const std::vector<Ty>& b, std::vector<Ty>& x) {
	int d = (int)b.size();
	std::vector<Ty> y(d);
	for (int i = 0; i < d; ++i) {
		Ty sum = 0.0;
		for (int k = 0; k < i; ++k) sum += LU[i*d + k] * y[k];
		y[i] = (b[i] - sum) / LU[i*d + i];

	}
	for (int i = d - 1; i >= 0; --i) {
		Ty sum = 0.0;
		for (int k = i + 1; k < d; ++k) {
			//cout << "lu[] = " << LU[i*d+k] << " x[" << k << "] = " << x[k] << endl;
			sum += LU[i*d + k] * x[k];
		}
		//cout << "sum " << sum << endl;
		x[i] = (y[i] - sum); // not dividing by diagonals
	}
}

/// DEPRECATED
// CroutFDP generates an LU decomposition using fused dot products during factorization
template<size_t nbits, size_t es, size_t capacity = 10>
void CroutFDP(std::vector< sw::unum::posit<nbits, es> >& S, std::vector< sw::unum::posit<nbits, es> >& D) {
	size_t d = size_t(std::sqrt(S.size()));
	assert(S.size() == d * d);
	assert(D.size() == d * d);
	using namespace sw::unum;
	for (int k = 0; k < d; ++k) {
		for (int i = k; i < d; ++i) {
			quire<nbits, es, capacity> q = 0.0;
			//for (int p = 0; p < k; ++p) q += D[i*d + p] * D[p*d + k];   if we had expression templates for the quire
			for (int p = 0; p < k; ++p) q += quire_mul(D[i*d + p], D[p*d + k]);
			posit<nbits, es> sum;
			sum.convert(q.to_value());     // one and only rounding step of the fused-dot product
			D[i*d + k] = S[i*d + k] - sum; // not dividing by diagonals
		}
		for (int j = k + 1; j < d; ++j) {
			quire<nbits, es, capacity> q = 0.0;
			//for (int p = 0; p < k; ++p) q += D[k*d + p] * D[p*d + j];   if we had expression templates for the quire
			for (int p = 0; p < k; ++p) q += quire_mul(D[k*d + p], D[p*d + j]);
			posit<nbits, es> sum;
			sum.convert(q.to_value());   // one and only rounding step of the fused-dot product
			D[k*d + j] = (S[k*d + j] - sum) / D[k*d + k];
		}
	}
}

/// DEPRECATED
// SolveCrout takes an LU decomposition, LU, and a right hand side vector, b, and produces a result, x.
template<size_t nbits, size_t es, size_t capacity = 10>
void SolveCroutFDP(const std::vector< sw::unum::posit<nbits, es> >& LU, const std::vector< sw::unum::posit<nbits, es> >& b, std::vector< sw::unum::posit<nbits, es> >& x) {
	using namespace sw::unum;
	int d = int(b.size());
	std::vector< posit<nbits, es> > y(d);
	for (int i = 0; i < d; ++i) {
		quire<nbits, es, capacity> q = 0.0;
		// for (int k = 0; k < i; ++k) q += LU[i*d + k] * y[k];   if we had expression templates for the quire
		for (int k = 0; k < i; ++k) q += quire_mul(LU[i*d + k], y[k]);
		posit<nbits, es> sum;
		sum.convert(q.to_value());   // one and only rounding step of the fused-dot product
		y[i] = (b[i] - sum) / LU[i*d + i];
	}
	for (int i = d - 1; i >= 0; --i) {
		quire<nbits, es, capacity> q = 0.0;
		// for (int k = i + 1; k < d; ++k) q += LU[i*d + k] * x[k];   if we had expression templates for the quire
		for (int k = i + 1; k < d; ++k) {
			//cout << "lu[] = " << LU[i*d + k] << " x[" << k << "] = " << x[k] << endl;
			q += quire_mul(LU[i*d + k], x[k]);
		}
		posit<nbits, es> sum;
		sum.convert(q.to_value());   // one and only rounding step of the fused-dot product
		//cout << "sum " << sum << endl;
		x[i] = (y[i] - sum); // not dividing by diagonals
	}
}

// Doolittle uses unit diagonals for the lower triangle
template<typename Ty>
void Doolittle(std::vector<Ty>& S, std::vector<Ty>& D) {
	size_t d = size_t(std::sqrt(S.size()));
	assert(S.size() == d*d);
	assert(D.size() == d*d);
	for (size_t k = 0; k < d; ++k) {
		for (size_t j = k; j < d; ++j) {
			Ty sum = 0.;
			for (size_t p = 0; p < k; ++p) sum += D[k*d + p] * D[p*d + j];
			D[k*d + j] = (S[k*d + j] - sum); // not dividing by diagonals
		}
		for (size_t i = k + 1; i < d; ++i) {
			Ty sum = 0.;
			for (size_t p = 0; p < k; ++p) sum += D[i*d + p] * D[p*d + k];
			D[i*d + k] = (S[i*d + k] - sum) / D[k*d + k];
		}
	}
}
// SolveDoolittle takes an LU decomposition, LU, and a right hand side vector, b, and produces a result, x.
template<typename Ty>
void SolveDoolittle(const std::vector<Ty>& LU, const std::vector<Ty>& b, std::vector<Ty>& x) {
	int d = (int)b.size();
	std::vector<Ty> y(d);
	for (int i = 0; i<d; ++i) {
		Ty sum = 0.;
		for (int k = 0; k<i; ++k)sum += LU[i*d + k] * y[k];
		y[i] = (b[i] - sum); // not dividing by diagonals
	}
	for (int i = d - 1; i >= 0; --i) {
		Ty sum = 0.;
		for (int k = i + 1; k<d; ++k)sum += LU[i*d + k] * x[k];
		x[i] = (y[i] - sum) / LU[i*d + i];
	}
}
#endif 

} // namespace hprblas
} // namespace sw
