// fused_mm.cpp: example program to demonstrate BLAS L# Reproducible Matrix-Matrix product
//
// Copyright (C) 2017-2020 Stillwater Supercomputing, Inc.
//
// This file is part of the HPRBLAS project, which is released under an MIT Open Source license.

// warning C4996: 'std::copy::_Unchecked_iterators::_Deprecate': Call to 'std::copy' with parameters that may be unsafe - this call relies on the caller to check that the passed values are correct. To disable this warning, use -D_SCL_SECURE_NO_WARNINGS. See documentation on how to use Visual C++ 'Checked Iterators'
#pragma warning(disable : 4996) 

// enable posit arithmetic exceptions
#define POSIT_THROW_ARITHMETIC_EXCEPTION 1
#include <hprblas>
// matrix generators
#include <generators/matrix_generators.hpp>
#include <utils/print_utils.hpp>

int main(int argc, char** argv)
try {
	using namespace std;
	using namespace mtl;
	using namespace sw::universal;
	using namespace sw::hprblas;

	int nrOfFailedTestCases = 0;

	// configure the number system
	constexpr size_t nbits = 64;
	constexpr size_t es = 2;
	using Scalar = typename sw::universal::posit<nbits, es>;
	using Matrix = typename mtl::mat::dense2D<Scalar>;

	size_t N = 10;
	Matrix H(N, N);
	GenerateHilbertMatrix(H);

	Matrix Hinv(N, N);
	GenerateHilbertMatrixInverse(Hinv);

	// using an inefficient linear form
	Matrix I1(N, N);
	I1 = sw::hprblas::fmm(H, Hinv);
	Scalar lcm = Scalar(HilbertScalingFactor(N));
	I1 = I1 / lcm;
	printMatrix(cout, "H * H^-1", I1);

	// using a blocked form
	constexpr size_t blockSize = 7;
	Matrix I2(N, N);
	I2 = sw::hprblas::bfmm(H, Hinv, blockSize);
	I2 = I2 / lcm;
	printMatrix(cout, "H * H^-1", I2);

	mtl::mat::dense2D<Scalar> eye(N, N);
	eye = Scalar(1);

	if (!isEqual(I1, eye)) ++nrOfFailedTestCases;

	if (nrOfFailedTestCases) cout << "FAIL" << endl; else cout << "PASS" << endl;

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
