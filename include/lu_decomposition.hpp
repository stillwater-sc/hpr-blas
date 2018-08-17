// lu_decomposition.cpp example program comparing float vs posit equation solver
//
// Copyright (C) 2017 Stillwater Supercomputing, Inc.
//
// This file is part of the universal numbers project, which is released under an MIT Open Source license.

#include <vector>
#define POSIT_VERBOSE_OUTPUT
#define QUIRE_TRACE_ADD
#include <posit>


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
template<typename Ty>
void Crout(std::vector<Ty>& S, std::vector<Ty>& D) {
	size_t d = size_t(std::sqrt(S.size()));
	assert(S.size() == d*d);
	assert(D.size() == d*d);
	for (size_t k = 0; k < d; ++k) {
		for (size_t i = k; i<d; ++i) {
			Ty sum = 0.;
			for (size_t p = 0; p<k; ++p) sum += D[i*d + p] * D[p*d + k];
			D[i*d + k] = S[i*d + k] - sum; // not dividing by diagonals
		}
		for (size_t j = k + 1; j < d; ++j) {
			Ty sum = 0.;
			for (int p = 0; p<k; ++p) sum += D[k*d + p] * D[p*d + j];
			D[k*d + j] = (S[k*d + j] - sum) / D[k*d + k];
		}
	}
}


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

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Crout using MTL data structures

// Crout method using MTL data structures
template<typename Ty>
void Crout(const mtl::dense2D<Ty>& S, mtl::dense2D<Ty>& D) {
	assert(size(S) == size(D));
	size_t d = num_rows(S);
	for (size_t k = 0; k < d; ++k) {
		for (size_t i = k; i<d; ++i) {
			Ty sum = 0.;
			for (size_t p = 0; p<k; ++p) sum += D[i][p] * D[p][k];
			D[i][k] = S[i][k] - sum; // not dividing by diagonals
		}
		for (size_t j = k + 1; j < d; ++j) {
			Ty sum = 0.;
			for (int p = 0; p<k; ++p) sum += D[k][p] * D[p][j];
			D[k][j] = (S[k][j] - sum) / D[k][k];
		}
	}
}


// SolveCrout takes an LU decomposition, LU, and a right hand side vector, b, and produces a result, x.
template<typename Ty>
void SolveCrout(const mtl::dense2D<Ty>& LU, const mtl::dense_vector<Ty>& b, mtl::dense_vector<Ty>& x) {
	size_t d = size(b);
	mtl::dense_vector<Ty> y(d);
	for (int i = 0; i < d; ++i) {
		Ty sum = 0.0;
		for (int k = 0; k < i; ++k) sum += LU[i][k] * y[k];
		y[i] = (b[i] - sum) / LU[i][i];

	}
	for (long i = long(d) - 1; i >= 0; --i) {
		Ty sum = 0.0;
		for (int k = i + 1; k < d; ++k) {
			//cout << "lu[] = " << LU[i][k] << " x[" << k << "] = " << x[k] << endl;
			sum += LU[i][k] * x[k];
		}
		//cout << "sum " << sum << endl;
		x[i] = (y[i] - sum); // not dividing by diagonals
	}
}

template<size_t nbits, size_t es, size_t capacity = 10>
void CroutFDP(std::vector< sw::unum::posit<nbits, es> >& S, std::vector< sw::unum::posit<nbits, es> >& D) {
	size_t d = size_t(std::sqrt(S.size()));
	assert(S.size() == d*d);
	assert(D.size() == d*d);
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

///////////////////////////////////////////////////////////////////////////////////
/// CroutFDP with MTL data structures

template<size_t nbits, size_t es, size_t capacity = 10>
void CroutFDP(mtl::dense2D< sw::unum::posit<nbits, es> >& S, mtl::dense2D< sw::unum::posit<nbits, es> >& D) {
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
			D[i][k] = S[i][k] - sum; // not dividing by diagonals
		}
		for (int j = k + 1; j < d; ++j) {
			quire<nbits, es, capacity> q;
			q.reset();
			//for (int p = 0; p < k; ++p) q += D[k][p] * D[p][j];   if we had expression templates for the quire
			for (int p = 0; p < k; ++p) q += quire_mul(D[k][p], D[p][j]);
			posit<nbits, es> sum;
			convert(q.to_value(), sum);   // one and only rounding step of the fused-dot product
			D[k][j] = (S[k][j] - sum) / D[k][k];
		}
	}
}

// SolveCrout takes an LU decomposition, LU, and a right hand side vector, b, and produces a result, x.
template<size_t nbits, size_t es, size_t capacity = 10>
void SolveCroutFDP(const mtl::dense2D< sw::unum::posit<nbits, es> >& LU, const mtl::dense_vector< sw::unum::posit<nbits, es> >& b, mtl::dense_vector< sw::unum::posit<nbits, es> >& x) 
{
	using namespace sw::unum;

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


	}  // namespace hprblas
} // namespace sw