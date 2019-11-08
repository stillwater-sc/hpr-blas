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
#include <stdio.h>
#include <math.h>
#include <universal/posit/posit>

#define  SIZE 25


typedef double MAT[SIZE][SIZE], VEC[SIZE];

void choldc1(int,MAT,VEC); 

/* -----------------------------------------------
        Cholesky decomposition.

        input    n  size of matrix
        input    A  Symmetric positive def. matrix
        output   a  lower deomposed matrix
        uses        choldc1(int,MAT,VEC)
   ----------------------------------------------- */
void choldc(int n,MAT A, MAT a) {
	int i,j;
	VEC p;
	for (i = 0; i < n; i++) 
		for (j = 0; j < n; j++) 
			a[i][j] = A[i][j];
	choldc1(n, a, p);
	for (i = 0; i < n; i++) {
		a[i][i] = p[i];
		for (j = i + 1; j < n; j++) {
			a[i][j] = 0;
		}
	}
}
 
/* -----------------------------------------------------
         Inverse of Cholesky decomposition.

         input    n  size of matrix
         input    A  Symmetric positive def. matrix
         output   a  inverse of lower deomposed matrix
         uses        choldc1(int,MAT,VEC)         
   ----------------------------------------------------- */
        void choldcsl(int n, MAT A, MAT a) {
	  int i,j,k; double sum;
	  VEC p;
          for (i = 0; i < n; i++) 
	    for (j = 0; j < n; j++) 
	      a[i][j] = A[i][j];
          choldc1(n, a, p);
          for (i = 0; i < n; i++) {
            a[i][i] = 1 / p[i];
            for (j = i + 1; j < n; j++) {
              sum = 0;
              for (k = i; k < j; k++) {
                sum -= a[j][k] * a[k][i];
	      }
              a[j][i] = sum / p[j];
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
        double choldet(int n, MAT a) {
	   MAT c; double d=1; int i;
           choldc(n,a,c);
           for (i = 0; i < n; i++)  d *= c[i][i];
           return d * d;
	}
 
/* ---------------------------------------------------
        Matrix inverse using Cholesky decomposition

        input    n  size of matrix
        input	  A  Symmetric positive def. matrix
        output   a  inverse of A
        uses        choldc1(MAT, VEC)
   --------------------------------------------------- */
        void cholsl(int n, MAT A, MAT a) {
	  int i,j,k;
          choldcsl(n,A,a);
          for (i = 0; i < n; i++) {
            for (j = i + 1; j < n; j++) {
              a[i][j] = 0.0;
	    }
	  }
          for (i = 0; i < n; i++) {
            a[i][i] *= a[i][i];
            for (k = i + 1; k < n; k++) {
              a[i][i] += a[k][i] * a[k][i];
	    }
            for (j = i + 1; j < n; j++) {
              for (k = j; k < n; k++) {
                a[i][j] += a[k][i] * a[k][j];
	      }
	    }
	  }
          for (i = 0; i < n; i++) {
            for (j = 0; j < i; j++) {
              a[i][j] = a[j][i];
	    }
	  }
	}

/* ----------------------------------------------------
        main method for Cholesky decomposition.

        input         n  size of matrix
        input/output  a  Symmetric positive def. matrix
        output        p  vector of resulting diag of a
        author:       <Vadum Kutsyy, kutsyy@hotmail.com>
   ----------------------------------------------------- */
        void choldc1(int n, MAT a, VEC p) {
          int i,j,k;
          double sum;

	  for (i = 0; i < n; i++) {
            for (j = i; j < n; j++) {
              sum = a[i][j];
              for (k = i - 1; k >= 0; k--) {
                sum -= a[i][k] * a[j][k];
	      }
              if (i == j) {
                if (sum <= 0) {
                  printf(" a is not positive definite!\n");
		}
                p[i] = sqrt(sum);
	      }
              else {
                a[j][i] = sum / p[i];
	      }
	    }
	  }
	}

	//print a square real matrix A of size n with caption s
	//(n items per line).
	void MatPrint(const char *s, int n, MAT A) { int i,j; printf("\n %s\n", s);
	  for (i=0; i<n; i++) {
	    for (j=0; j<n; j++) 
	      printf(" %10.6f",A[i][j]);
            printf("\n");
          }
        }

        //check if matrix A is positive definite (return 1)
        //or not positive definite (return 0) 
        int Check_Matrix(int n, MAT A) {
          int i,j,k,result; double sum;
          result=1;
	  for (i=0; i<n; i++) {
	    for (j = i; j<n; j++) {
              sum = A[i][j];
              for (k = i - 1; k>=0; k--)
                sum -= A[i][k] * A[j][k];
              if (i == j)
                if (sum <= 0.0) result=0;
	    }
          }
	  return result;
	}

/******************************************
*    MULTIPLICATION OF TWO SQUARE REAL    *                                     
*    MATRICES                             *
* --------------------------------------- *                                     
* INPUTS:    A  MATRIX N*N                *                                     
*            B  MATRIX N*N                *                                     
*            N  INTEGER                   *                                     
* --------------------------------------- *                                     
* OUTPUTS:   C  MATRIX N*N PRODUCT A*B    *                                     
*                                         *
******************************************/
void MatMult(int n, MAT A,MAT B, MAT C) {
  double SUM;
  int I,J,K;
  for (I=0; I<n; I++)                                                                  
    for (J=0; J<n; J++) {
      SUM = 0.0;                                                                
      for (K=0; K<n; K++)
       SUM += A[I][K]*B[K][J];                                               
      C[I][J]=SUM;                                                            
    }                                                                   
}

//copy MAT A in MAT A1
void MatCopy(int n, MAT A, MAT A1) {
  int i,j;
  for (i=0; i<n; i++)
    for (j=0; j<n; j++)
      A1[i][j]=A[i][j];
}

// main program to demonstrate the use of function cholsl()
int main(int argc, char* argv[]) 
try {
	MAT A, A1, B, C; int i,j, n; char answer;
	printf(" Inversion of a square real symetric matrix by Cholesky method\n");
	printf(" (The matrix must positive def.).\n");

	n = 4;
	printf("\n Size = %d\n", n);

	// define lower half of symmetrical matrix
	A[0][0]= 5;
	A[1][0]=-1; A[1][1]= 5;
	A[2][0]=-1; A[2][1]=-1; A[2][2]= 5;
	A[3][0]=-1; A[3][1]=-1; A[3][2]=-1; A[3][3]= 5;

	// define upper half by symmetry
	for (i=0; i<n; i++)
		for (j=i+1; j<n; j++)	
			A[i][j]=A[j][i];

	if (Check_Matrix(n,A)) {
		MatCopy(n,A,A1);
		printf("\n Determinant = %f\n", choldet(n,A));
		MatPrint("Matrix A:",n,A);
		cholsl(n,A,B);
		MatPrint("Matrix Inv(A):",n,B);
	}
	else {
		printf("\n Sorry, this matrix is not positive definite !\n");
		return EXIT_FAILURE;
	}

	printf("\n Do you want a verification (y/n)?: ");
	scanf("%c", &answer);
	if (answer=='y') {
		MatMult(n,A1,B,C);
		MatPrint("Verification A * Inv(A) = I:",n,C);
	}
	printf("\n");
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
