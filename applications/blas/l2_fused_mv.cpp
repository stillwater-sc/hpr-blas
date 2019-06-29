// l2_fused_mv.cpp example program to demonstrate BLAS L2 Reproducible Matrix-Vector product
//
// Copyright (C) 2017-2019 Stillwater Supercomputing, Inc.
//
// This file is part of the universal numbers project, which is released under an MIT Open Source license.
#include "common.hpp"
// enable posit arithmetic exceptions
#define POSIT_THROW_ARITHMETIC_EXCEPTION 1
#define MTL_WITH_INITLIST
#include <hprblas>
#include "blas_utils.hpp"
#include "vector_utils.hpp"


int main(int argc, char** argv)
try {
	using namespace std;
	using namespace mtl;
	using namespace sw::unum;
	using namespace sw::hprblas;

	int nrOfFailedTestCases = 0;

	{
		using Posit64 = posit<64, 3>;
		Posit64 array[3][3] = {
			{ 1.0, 0.0, 0.0 },
		    { 0.0, 1.0, 0.0 },
		    { 0.0, 0.0, 1.0 }
		};
		size_t N = 5;
		dense2D<Posit64> A(N,N);
		dense_vector<Posit64> x(N);
		dense_vector<Posit64> b(N);
		x = 1.0;
		Posit64 scale = 5 * 7 * 9;
		GenerateHilbertMatrix(A, scale);
		printMatrix(cout, "Hilbert Matrix", A);
		matvec(A, x, b);
		printVector(cout, "b", b);
	}

	return (nrOfFailedTestCases > 0 ? EXIT_FAILURE : EXIT_SUCCESS);
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
catch (const std::runtime_error& err) {
	std::cerr << "Uncaught runtime exception: " << err.what() << std::endl;
	return EXIT_FAILURE;
}
catch (...) {
	std::cerr << "Caught unknown exception" << std::endl;
	return EXIT_FAILURE;
}
