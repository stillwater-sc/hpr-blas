// cholesky.cpp: Cholesky decomposition
//
// Copyright (C) 2017-2019 Stillwater Supercomputing, Inc.
//
// This file is part of the HPR-BLAS project, which is released under an MIT Open Source license.
#include <hprblas>

// CheckPositiveDefinite returns true if the Matrix A is a Positive Definite matrix
template<typename Matrix>
bool CheckPositiveDefinite(Matrix& A) {
	using Scalar = typename Matrix::value_type;
	assert(mtl::mat::num_rows(A) == mtl::mat::num_cols(A));  // got a be square
	int N = int(mtl::mat::num_rows(A));
    bool result = true;
	for (int i = 0; i < N; ++i) {
	    for (int j = i; j < N; ++j) {
              Scalar sum = A[i][j];
			  for (int k = i - 1; k >= 0; --k) {
				  sum -= A[i][k] * A[j][k];
			  }
              result = (i == j) && (sum <= 0.0) ? false : result;
	    }
    }
	return result;
}

/* 
* SAMPLE RUN:                                                    *
*                                                                *
* Inversion of a square real symetric matrix by Cholevsky method *
* (The matrix must positive definite).                           *
*                                                                *
* Size = 4                                                       *
*                                                                *
* Determinant = 432.000000                                       *
*                                                                *
* Matrix A:                                                      *
* 5.000000 -1.000000 -1.000000 -1.000000                         *
* -1.000000 5.000000 -1.000000 -1.000000                         *
* -1.000000 -1.000000 5.000000 -1.000000                         *
* -1.000000 -1.000000 -1.000000 5.000000                         *
*                                                                *
* Matrix Inv(A):                                                 *
* 0.250000 0.083333 0.083333 0.083333                            *
* 0.083333 0.250000 0.083333 0.083333                            *
* 0.083333 0.083333 0.250000 0.083333                            *
* 0.083333 0.083333 0.083333 0.250000                            *
*                                                                *
* SetupMatrix can generate the above test matrix with the call   *
* int N = 4;                                                     *
* Matrix A(N,N);                                                 *
* SetupMatrix(A, N)                                              *
*/
template<typename Matrix>
void SetupMatrix(Matrix& A, int bandwidth = 0) {
	int N = int(mtl::mat::num_rows(A));
	A = 0;
	// define lower half of symmetrical matrix
	for (int i = 0; i < N; ++i) {
		for (int j = 0; j < N; ++j) {
			if (i == j) {
				A[i][i] = 5;
			}
			else {
				int diff = (i - j) > 0 ? (i - j) : (j - i);
				if (diff <= bandwidth) {
					A[i][j] = -1;
				}
			}
		}
	}
}

// main program to demonstrate the use of function cholsl()
int main(int argc, char* argv[]) 
try {
	using namespace std;
	using namespace mtl;
	using Scalar = sw::unum::posit<32,2>;
	using Matrix = mtl::mat::dense2D<Scalar>;
	using Vector = mtl::vec::dense_vector<Scalar>;

	auto orig_precision = cout.precision();

	cout << " Inversion of a square real symmetric positive definite matrix by Cholesky method\n";
	constexpr int m = 4;
	constexpr int n = 4;
	constexpr unsigned N = 4;  m*n;
	cout << "matrix size is " << N << endl;
	Matrix A(N,N), Aorig(N, N), T(N, N);

	// intended for N = 4
	SetupMatrix(A, N);
	//mtl::mat::laplacian_setup(A, m, n);
	cout << "Original Matrix:\n" << A << endl;

	// save a copy for the verification phase, as our in-place Cholesky factorization is destructive
	Aorig = A;

	if (!CheckPositiveDefinite(A)) {
		cout << "This matrix is not positive definite !\n";
		return EXIT_FAILURE;
	}

	Scalar determinant = sw::hprblas::DeterminantSPD(A);
	cout << "Determinant = " << determinant << endl;

	/*
	Matrix UT(N, N);
	{ // syntax #1
		UT = Scalar(0);
		if (CholeskyDecomposition(A, UT)) {
			cout << "In-place Cholesky\nU^T Matrix:\n" << UT << endl;
		}
		else {
			cerr << "Couldn't factorize matrix\n";
		}
	}

	{ // syntax #2
		UT = CholeskyDecomposition(A);
		cout << "MATLAB Cholesky\nU^T Matrix:\n" << UT << endl;
	}

	

	Matrix U(N, N);
	U = trans(UT);
	T = UT * U;
	cout << "verification of the decomposition A = U^T * U = A\n" << with_format(T, 10, 3) << endl;
	{
		// non-quire norms
		cout << setprecision(30);
		Scalar oneNorm = l1_norm(T);
		cout << "l1-norm        " << oneNorm << endl;
		Scalar infNorm = linf_norm(T);
		cout << "linf-norm      " << infNorm << endl;
		Scalar frobenius = frobenius_norm(T);
		cout << "Frobenius-norm " << frobenius << endl;
		cout << setprecision(orig_precision);
	}

	// Inverse
	{
		Matrix Ainv(N, N);
		Ainv = UT;
		cholsl(A, Ainv);
		cout << "Matrix Inv(A)  :\n" << Ainv << endl;
	}
	*/
	Matrix Ainv(N, N);
	Ainv = sw::hprblas::Inverse(A);
	cout << "Matrix Inv(A)  :\n" << Ainv << endl;

	T = A * Ainv;
	cout << "Verification A * Inv(A) = I:\n" << with_format(T, 10, 3) << endl;

	{
		// non-quire norms
		cout << setprecision(30);
		Scalar oneNorm = sw::hprblas::l1_norm(T);
		cout << "l1-norm        " << oneNorm << endl;
		Scalar infNorm = sw::hprblas::linf_norm(T);
		cout << "linf-norm      " << infNorm << endl;
		Scalar frobenius = sw::hprblas::frobenius_norm(T);
		cout << "Frobenius-norm " << frobenius << endl;
		cout << setprecision(orig_precision);
	}

	
	{
		Matrix I(N, N);
		I = Scalar(1);

		Scalar oneNorm = sw::hprblas::l1_norm(I);
		cout << "l1-norm        " << oneNorm << endl;
		Scalar infNorm = sw::hprblas::linf_norm(I);
		cout << "linf-norm      " << infNorm << endl;
		Scalar frobenius = sw::hprblas::frobenius_norm(I);
		cout << "Frobenius-norm " << frobenius << endl;
		cout << setprecision(orig_precision);
	}
	
}
catch (char const* msg) {
	std::cerr << msg << std::endl;
	return EXIT_FAILURE;
}
catch (const posit_arithmetic_exception& err) {
	std::cerr << "Uncaught posit arithmetic exception: " << err.what() << std::endl;
	return EXIT_FAILURE;
}
catch (const quire_exception& err) {
	std::cerr << "Uncaught quire exception: " << err.what() << std::endl;
	return EXIT_FAILURE;
}
catch (const posit_internal_exception& err) {
	std::cerr << "Uncaught posit internal exception: " << err.what() << std::endl;
	return EXIT_FAILURE;
}
catch (std::runtime_error& err) {
	std::cerr << err.what() << std::endl;
	return EXIT_FAILURE;
}
catch (...) {
	std::cerr << "Caught unknown exception" << std::endl;
	return EXIT_FAILURE;
}
