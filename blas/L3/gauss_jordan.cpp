// gauss_jordan.cpp; example program comparing float vs posit matrix inversion algorithms
//
// Copyright (C) 2017-2019 Stillwater Supercomputing, Inc.
//
// This file is part of the HPR-BLAS project, which is released under an MIT Open Source license.

#include <chrono>
// configure the posit number system behavior
#define POSIT_ROUNDING_ERROR_FREE_IO_FORMAT 0
// configure the HPR-BLAS behavior
#define HPRBLAS_TRACE_ROUNDING_EVENTS 1
#include <hprblas>
#include <mtl_extensions.hpp>
#include <matrix_utils.hpp>
#include <print_utils.hpp>

template<typename Matrix>
Matrix GaussJordanInversion(const Matrix& A) {
	using Scalar = typename mtl::Collection<Matrix>::value_type;

	size_t m = A.num_rows();
	size_t n = A.num_cols();
	Matrix a(n,m), inv(m,n);
	a = A; // you need a deep copy
	inv = Scalar(1);

	// Performing elementary operations 
	for (unsigned i = 0; i < m; ++i)	{
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
	sw::hprblas::printMatrix(std::cout, "conversion", a);
	return inv;
}

void dummy() {
	using namespace std;
    int i, j, k, n = 3;
	float a[10][10] = { 
		{ 2,3,4 },
		{ 5,6,3 },
		{ 9,8,6 }
	};
	float d;
	
 
    for (i = 1; i <= n; i++)
        for (j = 1; j <= 2 * n; j++)
            if (j == (i + n))
                a[i][j] = 1;
 
    /************** partial pivoting **************/
    for (i = n; i > 1; i--)
    {
        if (a[i - 1][1] < a[i][1])
            for (j = 1; j <= n * 2; j++)
            {
                d = a[i][j];
                a[i][j] = a[i - 1][j];
                a[i - 1][j] = d;
            }
    }
    cout << "pivoted output: " << endl;
    for (i = 1; i <= n; i++)
    {
        for (j = 1; j <= n * 2; j++)
            cout << a[i][j] << "    ";
        cout << endl;
    }
    /********** reducing to diagonal  matrix ***********/
 
    for (i = 1; i <= n; i++)
    {
        for (j = 1; j <= n * 2; j++)
            if (j != i)
            {
                d = a[j][i] / a[i][i];
                for (k = 1; k <= n * 2; k++)
                    a[j][k] -= a[i][k] * d;
            }
    }
    /************** reducing to unit matrix *************/
    for (i = 1; i <= n; i++)
    {
        d = a[i][i];
        for (j = 1; j <= n * 2; j++)
            a[i][j] = a[i][j] / d;
    }
 
    cout << "your solutions: " << endl;
    for (i = 1; i <= n; i++)
    {
        for (j = n + 1; j <= n * 2; j++)
            cout << a[i][j] << "    ";
        cout << endl;
    }
 


}

int main(int argc, char** argv)
try {
	using namespace std;
	using namespace mtl;
	using namespace sw::unum;
	using namespace sw::hprblas;

	constexpr size_t nbits = 64;
	constexpr size_t es = 3;
	constexpr size_t capacity = 10;

	using Scalar = posit<nbits, es>;
	using Matrix = mtl::mat::dense2D<Scalar>;
	size_t N = 5;
	Matrix H(N, N);
	GenerateHilbertMatrix(H);
	Matrix Hinv = GaussJordanInversion(H);
	printMatrix(cout, "Hilbert matrix order 5", H);
	printMatrix(cout, "Hilbert inverse", Hinv);
	Matrix Href(N, N);
	GenerateHilbertMatrixInverse(Href);
	printMatrix(cout, "Hilbert inverse reference", Href);

	Matrix test(N, N);
	matmul(test, H, Hinv);
	printMatrix(cout, "H * H^-1", test);

	return EXIT_SUCCESS;
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
