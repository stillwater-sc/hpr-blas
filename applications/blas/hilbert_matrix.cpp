// hilbert_matrix.cpp: verification test program to study numerical behavior of very poorly conditioned matrices
//
// Copyright (C) 2017-2020 Stillwater Supercomputing, Inc.
//
// This file is part of the universal numbers project, which is released under an MIT Open Source license.
#include "common.hpp"
// Boost arbitrary precision floats
#include <boost/multiprecision/cpp_bin_float.hpp>
// Universal arbitrary precision integers
#include <universal/integer/integer>

// enable INITLISTs with MTL vectors and matrices
#ifndef MTL_WITH_INITLIST
#define MTL_WITH_INITLIST
#endif
// enable the following define to show the intermediate steps in the fused-dot product
// #define POSIT_VERBOSE_OUTPUT
#define QUIRE_TRACE_ADD
// enable posit arithmetic exceptions
#define POSIT_THROW_ARITHMETIC_EXCEPTION 1
// to capture all the possible bits, set this to 1 
#define POSIT_ROUNDING_ERROR_FREE_IO_FORMAT 0
#include <hprblas>
// matrix generators
#include <generators/matrix_generators.hpp>
#include <utils/print_utils.hpp>

// define a true 256-bit IEEE floating point type
constexpr size_t bits_in_octand = 113 + 128;
using cpp_bin_float_octand = boost::multiprecision::number<boost::multiprecision::backends::cpp_bin_float<bits_in_octand, boost::multiprecision::backends::digit_base_2, void, boost::int16_t, -16382, 16383>, boost::multiprecision::expression_template_option::et_off> ;
// define the floating point types (single, double, quad, octand)
using sp = boost::multiprecision::cpp_bin_float_single;
using dp = boost::multiprecision::cpp_bin_float_double;
using qp = boost::multiprecision::cpp_bin_float_quad;
using op = cpp_bin_float_octand;

// Generate a Hilbert matrix and inverse test cycle
template<typename Scalar>
void HilbertMatrixTest(size_t N, bool printA = false, bool printB = false) {
	using namespace std;
	using namespace mtl;
	using namespace sw::hprblas;
	cout << "Value type is " << typeid(Scalar).name() << endl;
	mtl::mat::dense2D<Scalar> A(N, N), B(N, N), C(N, N);
	GenerateHilbertMatrix(A);
	GenerateHilbertMatrixInverse(B);
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
template<size_t N = 12>     // default is set to 12 to catch any poor behavior of the (n over k) implementation
void EnumerateHilbertMatrices() {
	using namespace std;
	using namespace sw::unum;
	using namespace sw::hprblas;

	// TBD: N = 10 works, but N = 11 and up fails catastrophically... 
	// don't know where the source of the failure is
	// 1- the scaling factor still fits within a uint64, but is representable by the target number system
	// 2- does the calculation path itself not able to deal with the dynamic range
	// 3- are the binomial constants getting truncated  (the binomials fail under quad and oct precision)
	// 
	// ETLO: July 6th, 2019: The source turned out to be the calculation of the Binomial coefficients.
	// I was using a naive implementation for (n over k) = n!/(k!(n-k)!), and that was surpassing the
	// maximum value of the number system.


	cout << "posits\n";
	//	GenerateHilbertMatrixTest< posit< 56, 3> >(N);    // these are commented out as they are too small to capture the dynamic range of the Hilbert * inv(Hilbert) calculation
	//	GenerateHilbertMatrixTest< posit< 64, 3> >(N);
	HilbertMatrixTest< posit< 80, 3> >(N);
	HilbertMatrixTest< posit< 96, 3> >(N);
	HilbertMatrixTest< posit<128, 4> >(N);
	HilbertMatrixTest< posit<256, 5> >(N);

	cout << "IEEE floating point\n";
	//	GenerateHilbertMatrixTest<sp>(N, true, true);
	//	GenerateHilbertMatrixTest<dp>(N);
	HilbertMatrixTest<qp>(N);
	HilbertMatrixTest<op>(N);
}

// Generate scaling factors for a sequence of Hilbert matrix sizes constrained by [2..upperbound]
template<typename IntegerType>
void EnumerateHilbertMatrixScalingFactors(IntegerType upperbound = 30) {
/*
Scaling factors for a collection of Hilbert matrices
N is a size_t
N = 2 scaling factor = 6
N = 3 scaling factor = 60
N = 4 scaling factor = 420
N = 5 scaling factor = 2520
N = 6 scaling factor = 27720
N = 7 scaling factor = 360360
N = 8 scaling factor = 360360
N = 9 scaling factor = 12252240
N = 10 scaling factor = 232792560
N = 11 scaling factor = 232792560
N = 12 scaling factor = 5354228880
N = 13 scaling factor = 26771144400
N = 14 scaling factor = 80313433200
N = 15 scaling factor = 2329089562800
N = 16 scaling factor = 72201776446800
N = 17 scaling factor = 144403552893600
N = 18 scaling factor = 144403552893600
N = 19 scaling factor = 5342931457063200
N = 20 scaling factor = 5342931457063200
N = 21 scaling factor = 219060189739591200
N = 22 scaling factor = 9419588158802421600
N = 23 scaling factor = 8829725487644060640  <---- fail
N = 24 scaling factor = 7966566035391366368  <---- all fails from this point on
N = 25 scaling factor = 4328644540401716384
N = 26 scaling factor = 12668683009887232224
N = 27 scaling factor = 13388079202269110789
N = 28 scaling factor = 9861751895175310850
N = 29 scaling factor = 13506701862403363960
N = 30 scaling factor = 14643306287797112328

N is an integer<128>
Scaling factors for a collection of Hilbert matrices
N = 2 scaling factor = 6
N = 3 scaling factor = 60
N = 4 scaling factor = 420
N = 5 scaling factor = 2520
N = 6 scaling factor = 27720
N = 7 scaling factor = 360360
N = 8 scaling factor = 360360
N = 9 scaling factor = 12252240
N = 10 scaling factor = 232792560
N = 11 scaling factor = 232792560
N = 12 scaling factor = 5354228880
N = 13 scaling factor = 26771144400
N = 14 scaling factor = 80313433200
N = 15 scaling factor = 2329089562800
N = 16 scaling factor = 72201776446800
N = 17 scaling factor = 144403552893600
N = 18 scaling factor = 144403552893600
N = 19 scaling factor = 5342931457063200
N = 20 scaling factor = 5342931457063200
N = 21 scaling factor = 219060189739591200
N = 22 scaling factor = 9419588158802421600
N = 23 scaling factor = 9419588158802421600
N = 24 scaling factor = 442720643463713815200
N = 25 scaling factor = 3099044504245996706400
N = 26 scaling factor = 3099044504245996706400
N = 27 scaling factor = 164249358725037825439200
N = 28 scaling factor = 164249358725037825439200
N = 29 scaling factor = 164249358725037825439200
N = 30 scaling factor = 9690712164777231700912800
*/
	using namespace std;
	cout << "Scaling factors for a collection of Hilbert matrices\n";

	for (IntegerType i = 2; i <= upperbound; ++i) {
		cout << "N = " << i << " scaling factor = " << sw::hprblas::HilbertScalingFactor(i) << endl;
	}
}

// Test of Hilbert matrix inversion accuracy
// can't go higher than N = 21 as the Hilbert matrix generator is using size_t as scaling factor type
template<typename Scalar>
void HilbertInverseTest(size_t N) {
	using namespace std;
	using namespace sw::hprblas;

	size_t scale = HilbertScalingFactor(N);
	using Matrix = mtl::mat::dense2D<Scalar>;
	Matrix H(N, N), Hinv(N, N);
	GenerateHilbertMatrix(H);
	cout << H << endl;
	GenerateHilbertMatrixInverse(Hinv);
	cout << Hinv << endl;

	Matrix B(N, N), C(N, N);
	B = H * Hinv;
	cout << B << endl;
	C = B / Scalar(scale);
	// C = B / scale; this causes a compilation warning due to implicit conversion of scale to double which yields the warning loss of accuracy
	cout << C << endl;

	// Calculate the inverse via Cholesky
	Matrix Hinv2(N, N);
	Hinv2 = sw::hprblas::Inverse(H);
	B = H * Hinv2;
	cout << typeid(Scalar).name() << endl;
	cout << "Verification H * Hinv = I\n" << B << endl;
	cout << "minpos = " << std::numeric_limits<Scalar>::min() << endl;
}

int main(int argc, char** argv)
try {
	using namespace std;
	using namespace sw::unum;
	using namespace sw::hprblas;

	int nrOfFailedTestCases = 0;

	using IntegerType = sw::unum::integer<128>;
	EnumerateHilbertMatrixScalingFactors(IntegerType(30));
	
	//EnumerateHilbertMatrices();
	// N = 5 is the largest Hilbert matrix you can represent with floats
	HilbertMatrixTest< float >(5);
	// N = 9 is the largest Hilbert matrix you can represent with doubles
	HilbertMatrixTest< double >(9);

	// can't go higher than N = 21 as the Hilbert matrix generator is using size_t as scaling factor type
	// TODO: generalize to arbitrary type: for that you will need a generalized
	// integer<nbits> -> Scalar conversion layer
	size_t N = 5;
	HilbertInverseTest<double>(N);
	HilbertInverseTest<posit<32, 2>>(N);

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
