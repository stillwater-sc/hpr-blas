// linpack_unchained.cpp: system solver benchmark with exact solution
// This program completes the design goal of LINPACK, that is,
// create an idealized problem where both the input values and the correct answer 
// are expressible in the numerical vocabulary of the computing environment.
//
// Copyright (C) 2017-18 Stillwater Supercomputing, Inc.
//
// This file is part of the universal numbers project, which is released under an MIT Open Source license.

#include "common.hpp"

#include <iostream>
#include <typeinfo>

#include <boost/numeric/mtl/mtl.hpp>
#include <posit>
#include "vector_utils.hpp"
#include "matrix_utils.hpp"

using namespace std;
using namespace sw::unum;

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
	
	Matrix  A(4, 4), L(4, 4), U(4, 4), LU(4, 4);
	Vector	x(4), b(4), xx(4);
	double 	c = 1.0;
	
	{
		sw::hprblas::uniform_rand_sorted(A);
		cout << A << endl;
		sw::hprblas::printMatrix(cout, "Ordered matrix", A);
		LU = A;
		lu(LU);
		sw::hprblas::printMatrix(cout, "LU", LU);
		x = 1.0;
		sw::hprblas::printVector(cout, "x", x);
		b = A * x;
		sw::hprblas::printVector(cout, "b", b);
		xx = lu_solve(A, b);
		sw::hprblas::printVector(cout, "x", xx);
	}

	return 0;

	{
		mtl::mat::uniform_rand(A);  // uniform random with values between [0,1]
		LU = A;
		lu(LU);
		cout << A << endl;
		cout << LU << endl;
	}


	int nrOfFailedTestCases = 0;
	return nrOfFailedTestCases > 0 ? EXIT_FAILURE : EXIT_SUCCESS;
}
catch (char const* msg) {
	cerr << msg << endl;
	return EXIT_FAILURE;
}
catch (const std::runtime_error& err) {
	cerr << err.what() << endl;
	return EXIT_FAILURE;
}
catch (...) {
	cerr << "caught unknown exception" << endl;
	return EXIT_FAILURE;
}
