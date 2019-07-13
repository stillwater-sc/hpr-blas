// fused_mm.cpp example program to demonstrate BLAS L# Reproducible Matrix-Matrix product
//
// Copyright (C) 2017-2019 Stillwater Supercomputing, Inc.
//
// This file is part of the universal numbers project, which is released under an MIT Open Source license.

// enable posit arithmetic exceptions
#define POSIT_THROW_ARITHMETIC_EXCEPTION 1
#define MTL_WITH_INITLIST
#include <hprblas>
// matrix/vector helpers
#include <vector_utils.hpp>
#include <matrix_utils.hpp>
#include <print_utils.hpp>


int main(int argc, char** argv)
try {
	using namespace std;
	using namespace mtl;
	using namespace sw::unum;
	using namespace sw::hprblas;

	int nrOfFailedTestCases = 0;

	// configure the number system
	constexpr size_t nbits = 32;
	constexpr size_t es = 2;
	using Scalar = posit<nbits, es>;

	size_t N = 5;
	mtl::mat::dense2D<Scalar> H(N, N);
	GenerateHilbertMatrix(H);

	mtl::mat::dense2D<Scalar> Hinv(N, N);
	GenerateHilbertMatrixInverse(Hinv);

	mtl::mat::dense2D<Scalar> I(N, N);
	sw::hprblas::matmul(I, H, Hinv);

	Scalar lcm = HilbertScalingFactor(N);
	I = I / lcm;

	printMatrix(cout, "H * H^-1", I);

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
