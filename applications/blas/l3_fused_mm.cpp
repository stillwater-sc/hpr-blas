// l3_fused_mv.cpp example program to demonstrate BLAS L3 Reproducible Matrix-Matrix product
//
// Copyright (C) 2017-2019 Stillwater Supercomputing, Inc.
//
// This file is part of the universal numbers project, which is released under an MIT Open Source license.

#include "common.hpp"
// enable the following define to show the intermediate steps in the fused-dot product
// #define POSIT_VERBOSE_OUTPUT
#define QUIRE_TRACE_ADD
// enable posit arithmetic exceptions
#define POSIT_THROW_ARITHMETIC_EXCEPTION 1
#include <hprblas>
#include "blas_utils.hpp"
#include "vector_utils.hpp"
#include "matrix_utils.hpp"

template<typename Scalar>
void GenerateHilbertMatrixTest(size_t N) {
	using namespace std;
	using namespace mtl;
	using namespace sw::hprblas;
	cout << "Value type is " << typeid(Scalar).name() << endl;
	mtl::mat::dense2D<Scalar> A(N, N), B(N, N), C(N, N);
	GenerateHilbertMatrix(A, Scalar(3*5*7*11*13*17));
	GenerateHilbertMatrixInverse(B);
	C = 0;
	matmul(C, A, B);
	printMatrix(cout, "A matrix", A);
	printMatrix(cout, "B matrix", B);
	printMatrix(cout, "C matrix", C);
}

int main(int argc, char** argv)
try {
	using namespace std;
	using namespace sw::unum;
	using namespace sw::hprblas;

	int nrOfFailedTestCases = 0;

	constexpr size_t N = 10;
	cout << "posits\n";
	GenerateHilbertMatrixTest< posit< 32, 2> >(N);
	GenerateHilbertMatrixTest< posit< 64, 3> >(N);
//	GenerateHilbertMatrixTest< posit<128, 4> >(N);

	cout << "IEEE floating point\n";
	GenerateHilbertMatrixTest<      float>(N);
	GenerateHilbertMatrixTest<     double>(N);
//	GenerateHilbertMatrixTest<long double>(N);

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
catch (std::runtime_error& err) {
	std::cerr << err.what() << std::endl;
	return EXIT_FAILURE;
}
catch (...) {
	std::cerr << "Caught unknown exception" << std::endl;
	return EXIT_FAILURE;
}