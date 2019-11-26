#pragma once
// cholesky.hpp implementation of the Cholesky solver
//
// Copyright (C) 2017-2020 Stillwater Supercomputing, Inc.
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


// In-place Cholesky factorization of an SPD matrix A, generating the factor UT and the diagonal
template<typename Matrix, typename Vector>
bool CholeskyFactorization(Matrix& A, Vector& diagonal) {
	using Scalar = typename Matrix::value_type;
	int N = int(mtl::mat::num_cols(A));

	for (int i = 0; i < N; ++i) {
		for (int j = i; j < N; ++j) {
			Scalar sum = A[i][j];
			for (int k = i - 1; k >= 0; k--) {
				sum -= A[i][k] * A[j][k];
			}
			if (i == j) {
				if (sum <= 0) {
					std::cerr << "Matrix is not positive definite!\n";
					return false;
				}
				diagonal[i] = sqrt(sum);
			}
			else {
				A[j][i] = sum / diagonal[i];
			}
		}
	}
	return true;
}

// Cholesky factorization of an SPD matrix A, generating the factor UT and the diagonal
template<typename Matrix, typename Vector>
bool CholeskyFactorization(const Matrix& A, Matrix& UT, Vector& diagonal) {
	using Scalar = typename Matrix::value_type;
	int N = int(mtl::mat::num_cols(A));

	UT = A;
	for (int i = 0; i < N; ++i) {
		for (int j = i; j < N; ++j) {
			Scalar sum = A[i][j];
			for (int k = i - 1; k >= 0; k--) {
				sum -= UT[i][k] * UT[j][k];
			}
			if (i == j) {
				if (sum <= 0) {
					std::cout << "Matrix is not positive definite!\n";
					return false;
				}
				diagonal[i] = sqrt(sum);
			}
			else {
				UT[j][i] = sum / diagonal[i];
			}
		}
	}
	return true;
}

// CholeskyDecomposition takes a Symmetric Positive Definite Matrix A, and returns the lower triangular Cholesky factorized matrix U transpose (UT)
template<typename Matrix>
bool CholeskyDecomposition(const Matrix& A, Matrix& UT) {
	assert(mtl::mat::num_rows(A) == mtl::mat::num_cols(A)); // assert squareness
	using Scalar = typename Matrix::value_type;
	int N = int(mtl::mat::num_cols(A));
	mtl::vec::dense_vector<Scalar> diagonal(N);

	bool success = CholeskyFactorization(A, UT, diagonal);
	if (!success) return false;

	// complete the UT by nulling out the upper triangular part
	for (int i = 0; i < N; ++i) {
		UT[i][i] = diagonal[i];
		for (int j = i + 1; j < N; ++j) {
			UT[i][j] = 0;
		}
	}

	return true;
}

// CholeskyDecomposition takes a Symmetric Positive Definite Matrix A, and returns the lower triangular Cholesky factorized matrix U transpose (UT)
template<typename Matrix>
Matrix CholeskyDecomposition(const Matrix& A) {
	using Scalar = typename Matrix::value_type;
	int N = int(mtl::mat::num_cols(A));
	Matrix UT(N, N);
	mtl::vec::dense_vector<Scalar> diagonal(N);

	UT = 0;
	bool success = CholeskyFactorization(A, UT, diagonal);
	if (!success) {
		UT = 0;
		return UT;
	}

	// complete the UT by nulling out the upper triangular part
	for (int i = 0; i < N; ++i) {
		UT[i][i] = diagonal[i];
		for (int j = i + 1; j < N; ++j) {
			UT[i][j] = 0;
		}
	}
	return UT;
}

/* -----------------------------------------------------
		 Inverse of Cholesky decomposition.

		 input    A     Symmetric positive def. matrix
		 output   Linv  inverse of lower decomposed matrix
		 uses        choldc1(Matrix,Vector)
   ----------------------------------------------------- */
template<typename Matrix>
void choldcsl(const Matrix& A, Matrix& Linv) {
	using Scalar = typename Matrix::value_type;
	int N = int(mtl::mat::num_cols(A));
	for (int i = 0; i < N; ++i) {
		for (int j = 0; j < N; ++j) {
			Linv[i][j] = A[i][j];
		}
	}

	mtl::vec::dense_vector<Scalar> p(N);
	CholeskyFactorization(Linv, p);
	for (int i = 0; i < N; ++i) {
		Linv[i][i] = 1 / p[i];
		for (int j = i + 1; j < N; ++j) {
			Scalar sum = Scalar(0);
			for (int k = i; k < j; ++k) {
				sum -= Linv[j][k] * Linv[k][i];
			}
			Linv[j][i] = sum / p[j];
		}
	}
}

// Return value of the determinant of a symmetric, positive definite matrix via a Cholesky decomposition
template<typename Matrix>
typename Matrix::value_type DeterminantSPD(const Matrix& A) {
	int N = int(mtl::mat::num_cols(A));
	Matrix C(N, N);
	CholeskyDecomposition(A, C);
	using Scalar = typename Matrix::value_type;
	Scalar d = Scalar(1);
	for (int i = 0; i < N; ++i) {
		d *= C[i][i];
	}
	return d * d;
}

/* ---------------------------------------------------
		Matrix inverse using Cholesky decomposition

		input	  A  Symmetric positive def. matrix
		output   Ainv  inverse of A
		uses        choldc1(MAT, VEC)
   --------------------------------------------------- */
template<typename Matrix>
void cholsl(const Matrix& A, Matrix& Ainv) {
	choldcsl(A, Ainv);
	//	std::cout << "first   Ainv\n" << Ainv << std::endl;

	int N = int(mtl::mat::num_cols(A));
	for (int i = 0; i < N; ++i) {
		for (int j = i + 1; j < N; ++j) {
			Ainv[i][j] = 0.0;
		}
	}
	//	std::cout << "second   Ainv\n" << Ainv << std::endl;

	for (int i = 0; i < N; ++i) {
		Ainv[i][i] *= Ainv[i][i];
		for (int k = i + 1; k < N; ++k) {
			Ainv[i][i] += Ainv[k][i] * Ainv[k][i];
		}
		for (int j = i + 1; j < N; ++j) {
			for (int k = j; k < N; ++k) {
				Ainv[i][j] += Ainv[k][i] * Ainv[k][j];
			}
		}
	}
	//	std::cout << "third  Ainv\n" << Ainv << std::endl;

	for (int i = 0; i < N; ++i) {
		for (int j = 0; j < i; ++j) {
			Ainv[i][j] = Ainv[j][i];
		}
	}
	//	std::cout << "final  Ainv\n" << Ainv << std::endl;

}

template<typename Matrix>
Matrix Inverse(const Matrix& A) {
	assert(mtl::mat::num_rows(A) == mtl::mat::num_cols(A)); // assert squareness
	using Scalar = Matrix::value_type;
	int N = int(mtl::mat::num_cols(A));
	mtl::vec::dense_vector<Scalar> diagonal(N);
	Matrix Ainv(N, N);
	Ainv = A;
	for (int i = 0; i < N; ++i) {
		for (int j = i; j < N; ++j) {
			Scalar sum = Ainv[i][j];
			for (int k = i - 1; k >= 0; k--) {
				sum -= Ainv[i][k] * Ainv[j][k];
			}
			if (i == j) {
				if (sum <= 0) {
					std::cerr << "inverse: matrix is not positive definite!\n";
					Ainv = 0;
					return Ainv;
				}
				diagonal[i] = sqrt(sum);
			}
			else {
				Ainv[j][i] = sum / diagonal[i];
			}
		}
	}
	std::cout << "step 1\n" << Ainv << std::endl;

	for (int i = 0; i < N; ++i) {
		Ainv[i][i] = 1 / diagonal[i];
		for (int j = i + 1; j < N; ++j) {
			Scalar sum = Scalar(0);
			for (int k = i; k < j; ++k) {
				sum -= Ainv[j][k] * Ainv[k][i];
			}
			Ainv[j][i] = sum / diagonal[j];
		}
	}
	std::cout << "step 2\n" << Ainv << std::endl;

	for (int i = 0; i < N; ++i) {
		for (int j = i + 1; j < N; ++j) {
			Ainv[i][j] = 0.0;
		}
	}
	std::cout << "step 3\n" << Ainv << std::endl;

	for (int i = 0; i < N; ++i) {
		Ainv[i][i] *= Ainv[i][i];
		for (int k = i + 1; k < N; ++k) {
			Ainv[i][i] += Ainv[k][i] * Ainv[k][i];
		}
		for (int j = i + 1; j < N; ++j) {
			for (int k = j; k < N; ++k) {
				Ainv[i][j] += Ainv[k][i] * Ainv[k][j];
			}
		}
	}
	std::cout << "step 4\n" << Ainv << std::endl;

	for (int i = 0; i < N; ++i) {
		for (int j = 0; j < i; ++j) {
			Ainv[i][j] = Ainv[j][i];
		}
	}
	std::cout << "step 5\n" << Ainv << std::endl;

	return Ainv;
}

}  // namespace hprblas
} // namespace sw
