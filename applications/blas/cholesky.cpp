/*****************************************************************
*  Inversion of a symmetric matrix by Cholesky decomposition.    *
*  The matrix must be positive definite.                         * 
* -------------------------------------------------------------- *
* REFERENCE:                                                     *
*             From a Java Library Created by Vadim Kutsyy,       *
*             "http://www.kutsyy.com".                           *
* -------------------------------------------------------------- * 
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
*                      C++ Release By Jean-Pierre Moreau, Paris. *
*                                (www.jpmoreau.fr)               *
* -------------------------------------------------------------- *
* Release 1.1 : added verification Inv(A) * A = I.               *
*****************************************************************/
#include <hprblas>
//#include <universal/posit/posit>


// In-place Cholesky factorization of an SPD matrix A, generating the factor UT and the diagonal
template<typename Matrix, typename Vector>
bool CholeskyFactorization(Matrix& A, Vector& diagonal) {
	using Scalar = typename Matrix::value_type;
	int N = int(mtl::mat::num_cols(A));

	for (int i = 0; i < N; i++) {
		for (int j = i; j < N; j++) {
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
	for (int i = 0; i < N; i++) {
		for (int j = i; j < N; j++) {
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
void CholeskyDecomposition(const Matrix& A, Matrix& UT) {
	using Scalar = typename Matrix::value_type;
	int N = int(mtl::mat::num_cols(A));
	mtl::vec::dense_vector<Scalar> diagonal(N);

	CholeskyFactorization(A, UT, diagonal);

	// complete the UT by nulling out the upper triangular part
	for (int i = 0; i < N; i++) {
		UT[i][i] = diagonal[i];
		for (int j = i + 1; j < N; j++) {
			UT[i][j] = 0;
		}
	}
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
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			Linv[i][j] = A[i][j];
		}
	}

	mtl::vec::dense_vector<Scalar> p(N);
	CholeskyFactorization(Linv, p);
	for (int i = 0; i < N; i++) {
		Linv[i][i] = 1 / p[i];
		for (int j = i + 1; j < N; j++) {
			Scalar sum = Scalar(0);
			for (int k = i; k < j; k++) {
				sum -= Linv[j][k] * Linv[k][i];
			}
			Linv[j][i] = sum / p[j];
		}
	}
}

// Return value of the determinant of a symmetric, positive definite matrix via a Cholesky decomposition
template<typename Matrix, typename Scalar>
Scalar DeterminantSPD(const Matrix& A) {	
	int N = int(mtl::mat::num_cols(A));
	Matrix C(N, N);
	CholeskyDecomposition(A, C);
	Scalar d = Scalar(1);
	for (int i = 0; i < N; i++) {
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
	for (int i = 0; i < N; i++) {
		for (int j = i + 1; j < N; j++) {
			Ainv[i][j] = 0.0;
		}
	}
//	std::cout << "second   Ainv\n" << Ainv << std::endl;

	for (int i = 0; i < N; i++) {
		Ainv[i][i] *= Ainv[i][i];
		for (int k = i + 1; k < N; k++) {
			Ainv[i][i] += Ainv[k][i] * Ainv[k][i];
		}
		for (int j = i + 1; j < N; j++) {
			for (int k = j; k < N; k++) {
				Ainv[i][j] += Ainv[k][i] * Ainv[k][j];
			}
		}
	}
//	std::cout << "third  Ainv\n" << Ainv << std::endl;

	for (int i = 0; i < N; i++) {
		for (int j = 0; j < i; j++) {
			Ainv[i][j] = Ainv[j][i];
		}
	}
//	std::cout << "final  Ainv\n" << Ainv << std::endl;

}


// CheckPositiveDefinite returns true if the Matrix A is a Positive Definite matrix
template<typename Matrix>
bool CheckPositiveDefinite(Matrix& A) {
	using Scalar = typename Matrix::value_type;
	assert(mtl::mat::num_rows(A) == mtl::mat::num_cols(A));  // got a be square
	int N = int(mtl::mat::num_rows(A));
    bool result = true;
	for (int i = 0; i < N; i++) {
	    for (int j = i; j < N; j++) {
              Scalar sum = A[i][j];
			  for (int k = i - 1; k >= 0; k--) {
				  sum -= A[i][k] * A[j][k];
			  }
              result = (i == j) && (sum <= 0.0) ? false : result;
	    }
    }
	return result;
}


// main program to demonstrate the use of function cholsl()
int main(int argc, char* argv[]) 
try {
	using namespace std;
	using namespace mtl;
	using Scalar = sw::unum::posit<32,2>;
	using Matrix = mtl::mat::dense2D<Scalar>;
	using Vector = mtl::vec::dense_vector<Scalar>;

	cout << " Inversion of a square real symmetric positive definite matrix by Cholesky method\n";
	
	constexpr unsigned N = 4;
	cout << "matrix size is " << N << endl;
	Matrix A(N,N), Ainv(N,N), A1(N,N), UT(N,N), C(N,N), Linv(N,N);
	Vector p(N);

	// define lower half of symmetrical matrix
	A[0][0]= 5;
	A[1][0]=-1; A[1][1]= 5;
	A[2][0]=-1; A[2][1]=-1; A[2][2]= 5;
	A[3][0]=-1; A[3][1]=-1; A[3][2]=-1; A[3][3]= 5;

	// define upper half by symmetry
	for (unsigned i=0; i<N; i++)
		for (unsigned j=i+1; j<N; j++)	
			A[i][j]=A[j][i];

	// save a copy
	A1 = A;
	p = 0;

	if (!CheckPositiveDefinite(A)) {
		cout << "This matrix is not positive definite !\n";
		return EXIT_FAILURE;
	}
	cout << "Original Matrix:\n" << A << endl;
	Scalar determinant = DeterminantSPD<Matrix, Scalar>(A);
	cout << "Determinant = " << determinant << endl;

	UT = Scalar(0);
	CholeskyDecomposition(A, UT);
	cout << "Cholesky U^T Matrix:\n" << UT << endl;
	Ainv = UT;
	cholsl(A, Ainv);
	cout << "Matrix Inv(A)  :\n" << Ainv << endl;

	cout << "verification of the decomposition U^T * U = A\n";
	C = A1 * Ainv;
	cout << "Verification A * Inv(A) = I:\n" << C << endl;
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
