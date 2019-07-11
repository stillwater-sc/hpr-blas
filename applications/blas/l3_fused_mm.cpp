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
// to capture all the possible bits, set this to 1 
#define POSIT_ROUNDING_ERROR_FREE_IO_FORMAT 0
#include <hprblas>
#include "print_utils.hpp"
#include "vector_utils.hpp"
#include "matrix_utils.hpp"
#include <boost/multiprecision/cpp_bin_float.hpp>

constexpr size_t bits_in_octand = 113 + 128;
typedef boost::multiprecision::number<boost::multiprecision::backends::cpp_bin_float<bits_in_octand, boost::multiprecision::backends::digit_base_2, void, boost::int16_t, -16382, 16383>, boost::multiprecision::expression_template_option::et_off> cpp_bin_float_octand;

template<typename Scalar>
void GenerateHilbertMatrixTest(size_t N, bool printA = false, bool printB = false) {
	using namespace std;
	using namespace mtl;
	using namespace sw::hprblas;
	cout << "Value type is " << typeid(Scalar).name() << endl;
	mtl::mat::dense2D<Scalar> A(N, N), B(N, N), C(N, N);
	GenerateHilbertMatrix(A);
	GenerateHilbertMatrixInverse(B);
	C = 0;
	matmul(C, A, B);
	if (printA) printMatrix(cout, "A matrix", A);
	if (printB) printMatrix(cout, "B matrix", B);
	printMatrix(cout, "C matrix", C);
}


/*
Hilbert matrices have very large condition numbers and are a compact
test case for linear algebra algorithms and number systems. 
The Hilbert matrix and its inverse are known analytically and thus
provide a perfect verification mechanism to test the computational
dynamics of the algorithm and the arithmetic used.

The condition number of a matrix is the ratio between the largest
eigenvalue and the smallest eigenvalue. Any algorithms that need
all the information in the matrix will need to be mindful of the
interplay between precision and dynamic range.

As posits are tapered floating point systems, precision at large
and small scale is limited. But the quire enables deferred rounding
removing any rounding noise. For Hilbert matrices this makes all 
the difference.

The first step to work with Hilbert matrices is to scale the input
values to numbers that the number system can represent. The Hilbert
matrix coefficient are all rational, but 'difficult' for floating
point systems. Ratios such as 1/3, 1/6, 1/7, 1/9, 1/11, 1/13, etc.
induce rounding error trying to represent them. Scaling the elements
with the least common multiple resolves this problem.
*/
int main(int argc, char** argv)
try {
	using namespace std;
	using namespace sw::unum;
	using namespace sw::hprblas;
	using sp = boost::multiprecision::cpp_bin_float_single;
	using dp = boost::multiprecision::cpp_bin_float_double;
	using qp = boost::multiprecision::cpp_bin_float_quad;
	using op = cpp_bin_float_octand;

	int nrOfFailedTestCases = 0;

	cout << "Scaling factors for a collection of Hilbert matrices\n";
	for (size_t i = 2; i < 21; ++i) {
		HilbertScalingFactor(i);
	}

	// TBD: N = 10 works, but N = 11 and up fails catastrophically... 
	// don't know where the source of the failure is
	// 1- the scaling factor still fits within a uint64, but is representable by the target number system
	// 2- does the calculation path itself not able to deal with the dynamic range
	// 3- are the binomial constants getting truncated  (the binomials fail under quad and oct precision)
	// 
	// ETLO: July 6th, 2019: The source turned out to be the calculation of the Binomial coefficients.
	// I was using a naive implementation for (n over k) = n!/(k!(n-k)!), and that was surpassing the
	// maximum value of the number system.
	constexpr size_t N = 13;

	cout << "posits\n";
//	GenerateHilbertMatrixTest< posit< 56, 3> >(N);
	GenerateHilbertMatrixTest< posit< 64, 3> >(N);
	GenerateHilbertMatrixTest< posit< 80, 3> >(N);
	GenerateHilbertMatrixTest< posit< 96, 3> >(N);
	GenerateHilbertMatrixTest< posit<128, 4> >(N);

	cout << "IEEE floating point\n";

//	GenerateHilbertMatrixTest<sp>(N, true, true);
//	GenerateHilbertMatrixTest<dp>(N);
	GenerateHilbertMatrixTest<qp>(N);
//	GenerateHilbertMatrixTest<op>(N);

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