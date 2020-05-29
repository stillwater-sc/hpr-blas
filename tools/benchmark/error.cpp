// error.cpp: numerical error analysis benchmark environment
//
// Copyright (C) 2017-2020 Stillwater Supercomputing, Inc.
//
// This file is part of the HPRBLAS project, which is released under an MIT Open Source license.

// include and configure number systems
// configure the posit number system behavior
#define POSIT_ROUNDING_ERROR_FREE_IO_FORMAT 0
// configure the HPR-BLAS behavior
#define HPRBLAS_TRACE_ROUNDING_EVENTS 0
#include <hprblas>
// Boost arbitrary precision floats
#include <boost/multiprecision/cpp_bin_float.hpp>

// matrix generators
#include <generators/matrix_generators.hpp>

template<typename Scalar, typename Vector>
void GenerateNumericalAnalysisTestCase(const std::string& header, unsigned N, bool verbose = false) {
	using namespace std;
	using namespace mtl;
	using namespace sw::unum;
	using namespace sw::hprblas;

	std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n" << header << std::endl;

	// calculate the numerical error caused by the linear algebra computation
	Vector e(N), eprime(N), eabsolute(N), erelative(N), I(N);
	e = Scalar(1);
	I = Scalar(1);
	// TODO: it is not clear that for posits this would be a fused matrix-vector operation
	matvec(eprime, I, e);
	cout << "reference vector : " << e << '\n';
	cout << "error vector     : " << eprime << '\n';
	// absolute error
	eabsolute = e - eprime;
	cout << "absolute error vector : " << eabsolute << '\n';
	cout << "L1 norm   " << hex_format(l1_norm(eabsolute)) << "  " << l1_norm(eabsolute) << '\n';
	cout << "L2 norm   " << hex_format(l2_norm(eabsolute)) << "  " << l2_norm(eabsolute) << '\n';
	cout << "Linf norm " << hex_format(linf_norm(eabsolute)) << "  " << linf_norm(eabsolute) << '\n';

	// relative error
	cout << "relative error\n";
	Scalar relative_error;
	relative_error = l1_norm(eabsolute) / l1_norm(e);
	cout << "L1 norm   " << hex_format(relative_error) << "  " << relative_error << '\n';
	relative_error = l2_norm(eabsolute) / l2_norm(e);
	cout << "L2 norm   " << hex_format(relative_error) << "  " << relative_error << '\n';
	relative_error = linf_norm(eabsolute) / linf_norm(e);
	cout << "Linf norm " << hex_format(relative_error) << "  " << relative_error << '\n';

	// error volume
	cout << "error bounding box volume\n";
	cout << "Measured in Euclidean distance    : " << error_volume(linf_norm(eabsolute), N, false) << '\n';
	cout << "Measured in ULPs                  : " << error_volume(linf_norm(eabsolute), N, true) << " ulps^" << N << '\n';
	Scalar ulp = numeric_limits<Scalar>::epsilon();
	cout << "L-infinitiy norm measured in ULPs : " << linf_norm(eabsolute) / ulp << " ulps" << '\n';

	cout << endl;
}

// Benchmark Suite runner for numerical error analysis measurements
int main(int argc, char** argv)
try {
	using namespace std;
	using namespace mtl;
	using namespace sw::unum;
	using namespace sw::hprblas;

	int nrOfFailedTestCases = 0;

	// we need to enumerate along the following dimensions
	// 1- number systems: the key here is that when we have posits, we use FDP
	// 2- algorithms: different computational approaches to solve a system of linear equations
	// 3- matrices. easy, difficult, empirical

	// The measurement will always be the error to this equation:
	//    Ax = b, with a b that delivers the solution x = ones()

	// There is another measure and that is driven by the characteristic polynomial of the matrix:
	// If we have very wildy differing eigenvalues, then there are b vectors that can lift up small
	// eigenvalues compared to large eigenvalues. Those are situations in which we want to make
	// certain we don't have cancellation: a big eigenvalue multiplied by a small scaling factor
	// and a small eigenvalue multiplied by a big scaling factor.

	// the benchmark runner is structured as follows:
	// foreach test system
	//     pick a test size N
	//     foreach number system
	//         generate the test matrix A(N,N)
	//         generate the test right hand side: b(N) = A * ones()
	//         foreach algorithm
	//             generate the inverse or decomposition
	//             solve the system of equations: Ax = b
	//             measure the difference between result x and ones()
	using Scalar = sw::unum::posit<32, 2>;
	using Vector = mtl::vec::dense_vector<Scalar>;
	using Matrix = mtl::mat::dense2D<Scalar>;

	int N = 5;
	Matrix H(N, N);
	sw::hprblas::GenerateHilbertMatrix(H, false);
	Matrix Hinv = GaussJordanInversion(H);
	Matrix Href(N, N);
	GenerateHilbertMatrixInverse(Href);
	Matrix I(N, N);  // H * Hinv should yield the identity matrix
	// TODO: this is not clear that for posits this would be a fused matrix multiply
	matmul(I, H, Hinv);

	GenerateNumericalAnalysisTestCase<Scalar, Vector>("testing", 10, true);

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
