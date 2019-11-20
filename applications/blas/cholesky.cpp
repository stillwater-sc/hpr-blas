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


/* ----------------------------------------------------
		main method for Cholesky decomposition.

		input/output  A  Symmetric positive def. matrix
		output        p  vector of resulting diag of a
		author:       <Vadum Kutsyy, kutsyy@hotmail.com>
   ----------------------------------------------------- */
template<typename Matrix, typename Vector>
void choldc1(Matrix& A, Vector& p) {
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
					std::cout << " a is not positive definite!\n";
				}
				p[i] = sqrt(sum);
			}
			else {
				A[j][i] = sum / p[i];
			}
		}
	}
}

/* -----------------------------------------------
        Cholesky decomposition.

        input    A  Symmetric positive def. matrix
        output   L  lower decomposed matrix
        uses        choldc1(Matrix&, Vector&)
   ----------------------------------------------- */
template<typename Matrix>
void choldc(const Matrix& A, Matrix& L) {
	int N = int(mtl::mat::num_cols(A));
	for (int i = 0; i < N; i++) 
		for (int j = 0; j < N; j++) 
			L[i][j] = A[i][j];

	using Scalar = typename Matrix::value_type;
	mtl::vec::dense_vector<Scalar> p(N);
	choldc1(L, p);
	for (int i = 0; i < N; i++) {
		L[i][i] = p[i];
		for (int j = i + 1; j < N; j++) {
			L[i][j] = 0;
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
	mtl::vec::dense_vector<Scalar> p(N);
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			Linv[i][j] = A[i][j];
			choldc1(Linv, p);
			for (int i = 0; i < N; i++) {
				Linv[i][i] = 1 / p[i];
				for (j = i + 1; j < N; j++) {
					Scalar sum = Scalar(0);
					for (int k = i; k < j; k++) {
						sum -= Linv[j][k] * Linv[k][i];
					}
					Linv[j][i] = sum / p[j];
				}
			}
		}
	}
}

/* -----------------------------------------------------------------------------
        Computation of Determinant of the matrix using Cholesky decomposition

        input    n  size of matrix
        input    a  Symmetric positive def. matrix
        return      det(a)
        uses        choldc(int,MAT,MAT)
   ------------------------------------------------------------------------------ */
template<typename Matrix, typename Scalar>
Scalar choldet(const Matrix& A) {
	int N = int(mtl::mat::num_cols(A));
	Matrix C(N,N); 
	choldc(A, C);
	std::cout << "choldc: \n" << C << std::endl;
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
	int N = int(mtl::mat::num_cols(A));
	for (int i = 0; i < N; i++) {
		for (int j = i + 1; j < N; j++) {
			Ainv[i][j] = 0.0;
		}
	}
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
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < i; j++) {
			Ainv[i][j] = Ainv[j][i];
		}
	}
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
              for (int k = i - 1; k >= 0; k--)
                sum -= A[i][k] * A[j][k];
              if (i == j)
                if (sum <= 0.0) result = false;
	    }
    }
	return result;
}


// main program to demonstrate the use of function cholsl()
int main(int argc, char* argv[]) 
try {
	using namespace std;
	using namespace mtl;
	using Scalar = double;
	using Matrix = mtl::mat::dense2D<Scalar>;

	cout << " Inversion of a square real symmetric positive definite matrix by Cholesky method\n";
	
	constexpr unsigned N = 4;
	cout << "matrix size is " << N << endl;
	Matrix A(N,N), A1(N,N), B(N,N), C(N,N);

	// define lower half of symmetrical matrix
	A[0][0]= 5;
	A[1][0]=-1; A[1][1]= 5;
	A[2][0]=-1; A[2][1]=-1; A[2][2]= 5;
	A[3][0]=-1; A[3][1]=-1; A[3][2]=-1; A[3][3]= 5;

	// define upper half by symmetry
	for (unsigned i=0; i<N; i++)
		for (unsigned j=i+1; j<N; j++)	
			A[i][j]=A[j][i];

	if (CheckPositiveDefinite(A)) {
		A1 = A;
		cout << "Determinant = " << choldet<Matrix,Scalar>(A) << endl;
		cout << "Original Matrix:\n" << A << endl;
		cholsl(A, B);
		cout << "Matrix Inv(A)  :\n" << B << endl;
	}
	else {
		cout << "Sorry, this matrix is not positive definite !\n";
		return EXIT_FAILURE;
	}

	cout << "verification of the decomposition U^T * U = A\n";
	C = A1 * B;
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
