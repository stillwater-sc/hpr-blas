// lu_decomposition.cpp example program comparing float vs posit equation solver
//
// Copyright (C) 2017-2019 Stillwater Supercomputing, Inc.
//
// This file is part of the universal numbers project, which is released under an MIT Open Source license.

#include "common.hpp"

#include <iostream>
#include <typeinfo>

#include <boost/numeric/mtl/mtl.hpp>
#include <universal/number/posit/posit.hpp>

using namespace std;
using namespace sw::universal;

// Turn it off for now
#define USE_POSIT

int main(int argc, char** argv)
try {
	const size_t nbits = 16;
	const size_t es = 1;
	const size_t vecSize = 32;

#ifdef USE_POSIT
	using Matrix = mtl::dense2D< posit<nbits, es> >;
	using Vector = mtl::dense_vector< posit<nbits, es> >;
#else
	using Matrix = mtl::dense2D<float>;
	using Vector = mtl::dense_vector<float>;
#endif


	Matrix  A(4, 4), L(4, 4), U(4, 4), AA(4, 4);
	Vector	v(4);
	double 	c = 1.0;
	
	for (unsigned i = 0; i < 4; i++)
		for (unsigned j = 0; j < 4; j++) {
			U[i][j] = i <= j ? c * (i + j + 2) : (0);
			L[i][j] = i > j ? c * (i + j + 1) : (i == j ? (1) : (0));
		}

	std::cout << "L is:\n" << L << "U is:\n" << U;
	A = L * U;
	std::cout << "A is:\n" << A;
	AA = adjoint(A);

	for (unsigned i = 0; i < 4; i++)
		v[i] = double(i);

	Vector b(A*v), b2(adjoint(A)*v);

	Matrix LU(A);
	lu(LU);
	std::cout << "LU decomposition of A is:\n" << LU;

	Matrix B(lu_f(A));
	std::cout << "LU decomposition of A (as function result) is:\n" << B;

	Vector v1(lu_solve_straight(A, b));
	std::cout << "v1 is " << v1 << "\n";

	Vector v2(lu_solve(A, b));
	std::cout << "v2 is " << v2 << "\n";

	mtl::dense_vector<unsigned> P;
	lu(A, P);
	std::cout << "LU with pivoting is \n" << with_format(A, 5, 2) << "Permutation is " << P << "\n";
	Vector v3(lu_apply(A, P, b));
	std::cout << "v3 is " << v3 << "\n";

	Vector v4(lu_adjoint_apply(A, P, b2));
	std::cout << "v4 is " << v4 << "\n";

	Vector v5(lu_adjoint_solve(AA, b));
	std::cout << "v5 is " << v5 << "\n";

	int nrOfFailedTestCases = 0;
	return nrOfFailedTestCases > 0 ? EXIT_FAILURE : EXIT_SUCCESS;
}
catch (char const* msg) {
	cerr << msg << endl;
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
	cerr << "caught unknown exception" << endl;
	return EXIT_FAILURE;
}
