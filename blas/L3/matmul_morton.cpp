// matmul_morton.cpp : example program comparing float vs posit matrix multiply algorithms
//
// Copyright (C) 2017-2019 Stillwater Supercomputing, Inc.
//
// This file is part of the HPR-BLAS project, which is released under an MIT Open Source license.

// warning C4996: 'std::copy::_Unchecked_iterators::_Deprecate': Call to 'std::copy' with parameters that may be unsafe - this call relies on the caller to check that the passed values are correct.
#pragma warning( disable : 4996)
#include "common.hpp"
#include <hprblas>
// utilities to generate and print vectors and matrices
#include <generators/matrix_generators.hpp>
#include "utils/matvec.hpp"

template<typename Scalar, size_t n = 10>
int MortonMatrixExamples()
{
	using namespace mtl; 
	using namespace mtl::mat;

	dense2D<Scalar>                            A(n, n), B(n, n);
	morton_dense<Scalar, doppled_64_row_mask>  C(n, n);

	hessian_setup(A, 3.0); 
	hessian_setup(B, 1.0);
	hessian_setup(C, 2.0);

	std::cout << A << std::endl;

	// Corresponds to A= B * B;
	mult(B, B, A);

	std::cout << A << std::endl;

	A = B * B;   // use BLAS
	A = B * C;   // use recursion + tiling from MTL4

	A += B * C;  // Increment A by the product of B and C
	A -= B * C;  // Likewise with decrement

	return 0;
}

template<typename Scalar>
void BenchmarkMorton() {
	using namespace std;
	using namespace chrono;
	using namespace mtl;
	for (int dim = 16; dim < 511; dim *= 2) {
		cout << "Matrix dimensions are: " << dim << " x " << dim << endl;

		steady_clock::time_point t1 = steady_clock::now();
		morton_dense<Scalar, doppled_64_row_mask>  A(dim, dim), B(dim, dim), C(dim, dim);
		steady_clock::time_point t2 = steady_clock::now();
		duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
		double elapsed = time_span.count();
		cout << "  Construction     " << elapsed << " seconds.\n";

		t1 = steady_clock::now();
		sw::hprblas::uniform_rand(A);
		sw::hprblas::uniform_rand(B);
		t2 = steady_clock::now();
		time_span = duration_cast<duration<double>>(t2 - t1);
		elapsed = time_span.count();
		cout << "  Random fill      " << elapsed << " seconds.\n";

		t1 = steady_clock::now();
		int N = 10;
		for (int i = 0; i < N; ++i) {
			C = A * B;
		}
		t2 = steady_clock::now();
		time_span = duration_cast<duration<double>> (t2 - t1);
		elapsed = time_span.count();
		cout << "  Matmul iteration " << elapsed << " seconds.\n";
		double flops = double(N)*dim*dim*dim / elapsed * 1.0e-9;
		cout << "  Performance:     " << flops << " GFLOPS Single Precision\n";
	}
}

template<typename Scalar>
void NumericalAttributes() {
	Scalar epsilon = std::numeric_limits<Scalar>::epsilon();
	Scalar max_exp = std::numeric_limits<Scalar>::max_exponent;
	Scalar min_exp = std::numeric_limits<Scalar>::min_exponent;
	std::cout << " epsilon " << epsilon << " min exp " << min_exp << " max exp " << max_exp << std::endl;
}

int main(int argc, char** argv)
try {
	using namespace std;
	using namespace sw::universal;
	using namespace sw::hprblas;
	using namespace mtl;

	// a 32-bit float and a <27,1> posit have the same number of significand bits around 1.0
	constexpr size_t nbits = 27;
	constexpr size_t es = 2;
	constexpr size_t capacity = 10;

	using Posit = posit<32, 2>;

	MortonMatrixExamples<float>();
//	BenchmarkMorton<Posit>();

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
