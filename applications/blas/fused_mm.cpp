// fused_mm.cpp example program to demonstrate BLAS L# Reproducible Matrix-Matrix product
//
// Copyright (C) 2017-2019 Stillwater Supercomputing, Inc.
//
// This file is part of the universal numbers project, which is released under an MIT Open Source license.


// warning C4996: 'std::copy::_Unchecked_iterators::_Deprecate': Call to 'std::copy' with parameters that may be unsafe - this call relies on the caller to check that the passed values are correct. To disable this warning, use -D_SCL_SECURE_NO_WARNINGS. See documentation on how to use Visual C++ 'Checked Iterators'
#pragma warning(disable : 4996) 

// enable posit arithmetic exceptions
#define POSIT_THROW_ARITHMETIC_EXCEPTION 1
#include <hprblas>
// utilities to generate and print vectors and matrices
#include "utils/matvec.hpp"

template<typename Matrix>
bool isEqual(const Matrix& lhs, const Matrix& rhs) {
	size_t r = lhs.num_rows();
	size_t c = lhs.num_cols();
	for (size_t i = 0; i < r; ++i) {
		for (size_t j = 0; j < c; ++j) {
			if (lhs[i][j] != rhs[i][j]) return false;
		}
	}
	return true;
}

int main(int argc, char** argv)
try {
	using namespace std;
	using namespace mtl;
	using namespace sw::unum;
	using namespace sw::hprblas;

	int nrOfFailedTestCases = 0;

	// configure the number system
	constexpr size_t nbits = 64;
	constexpr size_t es = 2;
	using Scalar = typename sw::unum::posit<nbits, es>;
	using Matrix = typename mtl::mat::dense2D<Scalar>;

	size_t N = 10;
	Matrix H(N, N);
	GenerateHilbertMatrix(H);

	Matrix Hinv(N, N);
	GenerateHilbertMatrixInverse(Hinv);

	// using an inefficient linear form
	Matrix I1(N, N);
	I1 = sw::hprblas::fmm(H, Hinv);
	Scalar lcm = HilbertScalingFactor(N);
	I1 = I1 / lcm;
	printMatrix(cout, "H * H^-1", I1);

	// using a blocked form
	constexpr size_t blockHeight = 5;
	constexpr size_t blockWidth = 5;
	Matrix I2(N, N);
	I2 = sw::hprblas::bfmm(H, Hinv, blockHeight, blockWidth);
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
