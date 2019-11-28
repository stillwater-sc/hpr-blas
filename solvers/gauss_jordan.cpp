// gauss_jordan.cpp : example program comparing float vs posit matrix inversion algorithms
//
// Copyright (C) 2017-2019 Stillwater Supercomputing, Inc.
//
// This file is part of the HPR-BLAS project, which is released under an MIT Open Source license.

//#include <chrono>
// Boost arbitrary precision floats
#include <boost/multiprecision/cpp_bin_float.hpp>

// configure the posit number system behavior
#define POSIT_ROUNDING_ERROR_FREE_IO_FORMAT 0
// configure the HPR-BLAS behavior
#define HPRBLAS_TRACE_ROUNDING_EVENTS 0
#include <hprblas>
#include <mtl_extensions.hpp>
// matrix generators
#include <generators/matrix_generators.hpp>
#include <utils/print_utils.hpp>

template<typename Scalar>
void GenerateNumericalAnalysisTestCase(const std::string& header, unsigned N, bool verbose = false) {
	using namespace std;
	using namespace mtl;
	using namespace sw::unum;
	using namespace sw::hprblas;

	std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n" << header << std::endl;
	using Vector = mtl::vec::dense_vector<Scalar>;
	using Matrix = mtl::mat::dense2D<Scalar>;

	Matrix H(N, N);
	GenerateHilbertMatrix(H, false);
	Matrix Hinv = GaussJordanInversion(H);
	Matrix Href(N, N);
	GenerateHilbertMatrixInverse(Href);
	Matrix I(N, N);  // H * Hinv should yield the identity matrix
	// TODO: this is not clear that for posits this would be a fused matrix multiply
	matmul(I, H, Hinv);

	if (verbose) {
		printMatrix(cout, "Hilbert matrix order 5", H);
		printMatrix(cout, "Hilbert inverse", Hinv);

		printMatrix(cout, "Hilbert inverse reference", Href);

		printMatrix(cout, "H * H^-1 => I", I);
	}

	// calculate the numerical error caused by the linear algebra computation
	Vector e(N), eprime(N), eabsolute(N), erelative(N);
	e = Scalar(1);
	// TODO: it is not clear that for posits this would be a fused matrix-vector operation
	matvec(eprime, I, e);
	printVector(cout, "reference vector", e);
	printVector(cout, "error vector", eprime);
	// absolute error
	eabsolute = e - eprime;
	printVector(cout, "absolute error vector", eabsolute);
	cout << "L1 norm   " << hex_format(l1_norm(eabsolute)) << "  " << l1_norm(eabsolute) << endl;
	cout << "L2 norm   " << hex_format(l2_norm(eabsolute)) << "  " << l2_norm(eabsolute) << endl;
	cout << "Linf norm " << hex_format(linf_norm(eabsolute)) << "  " << linf_norm(eabsolute) << endl;

	// relative error
	cout << "relative error\n";
	Scalar relative_error;
	relative_error = l1_norm(eabsolute) / l1_norm(e);
	cout << "L1 norm   " << hex_format(relative_error) << "  " << relative_error << endl;
	relative_error = l2_norm(eabsolute) / l2_norm(e);
	cout << "L2 norm   " << hex_format(relative_error) << "  " << relative_error << endl;
	relative_error = linf_norm(eabsolute) / linf_norm(e);
	cout << "Linf norm " << hex_format(relative_error) << "  " << relative_error << endl;

	// error volume
	cout << "error bounding box volume\n";
	cout << "Measured in Euclidean distance    : " << error_volume(linf_norm(eabsolute), N, false) << endl;
	cout << "Measured in ULPs                  : " << error_volume(linf_norm(eabsolute), N, true) << " ulps^" << N << endl;
	Scalar ulp = numeric_limits<Scalar>::epsilon();
	cout << "L-infinitiy norm measured in ULPs : " << linf_norm(eabsolute) / ulp << " ulps" << endl;
}

int main(int argc, char** argv)
try {
	using namespace std;
	using namespace mtl;
	using namespace sw::unum;
	using namespace sw::hprblas;
	using namespace boost::multiprecision;

	using sp = boost::multiprecision::cpp_bin_float_single;
	using dp = boost::multiprecision::cpp_bin_float_double;
	using qp = boost::multiprecision::cpp_bin_float_quad;


	unsigned N = 5; 
	GenerateNumericalAnalysisTestCase< posit<32, 2> >("posit<32,2>", N);
	cout << endl;
	GenerateNumericalAnalysisTestCase< float >("IEEE single precision", N);
	cout << endl;
	GenerateNumericalAnalysisTestCase< posit<64, 3> >("posit<64,3>", N);
	cout << endl;
	GenerateNumericalAnalysisTestCase< double >("IEEE double precision", N);
	cout << endl;
	GenerateNumericalAnalysisTestCase< posit<128, 4> >("posit<128,4>", N);
	cout << endl;
	GenerateNumericalAnalysisTestCase< qp >("IEEE quad precision", N);

	return EXIT_SUCCESS;
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
