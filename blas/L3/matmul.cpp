// matmul.cpp : example program comparing float vs posit matrix multiply algorithms
//
// Copyright (C) 2017-2020 Stillwater Supercomputing, Inc.
//
// This file is part of the HPRBLAS project, which is released under an MIT Open Source license.

// warning C4996: 'std::copy::_Unchecked_iterators::_Deprecate': Call to 'std::copy' with parameters that may be unsafe - this call relies on the caller to check that the passed values are correct.
#pragma warning( disable : 4996)
#include "common.hpp"
#define POSIT_FAST_POSIT_16_1 1
#define POSIT_FAST_POSIT_32_2 1
#include <universal/number/posit/posit>
#include <hprblas>
// matrix generators
#include <generators/matrix_generators.hpp>


template<typename Ty>
void MeasureMatrixMultiplyPerformance(const std::string& metric) {
	using namespace std::chrono;
	unsigned dim = 32;
	for (unsigned i = 5; i < 8; ++i, dim *= 2) {
		std::cout << "Matrix dimensions are: " << dim << " x " << dim << std::endl;

		steady_clock::time_point t1 = steady_clock::now();
		mtl::dense2D<Ty> A(dim, dim), B(dim, dim), C(dim, dim);
		steady_clock::time_point t2 = steady_clock::now();
		duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
		double elapsed = time_span.count();
		std::cout << "  Construction     " << elapsed << " seconds.\n";

		t1 = steady_clock::now();
		sw::hprblas::uniform_rand(A, -1.0, 1.0);
		sw::hprblas::uniform_rand(B, -1.0, 1.0);
		t2 = steady_clock::now();
		time_span = duration_cast<duration<double>>(t2 - t1);
		elapsed = time_span.count();
		std::cout << "  Random fill      " << elapsed << " seconds.\n";

		t1 = steady_clock::now();
		int N = 10;
		for (int i = 0; i < N; ++i) {
			C = A * B;
		}
		t2 = steady_clock::now();
		time_span = duration_cast<duration<double>> (t2 - t1);
		elapsed = time_span.count();
		std::cout << "  Matmul iteration " << elapsed << " seconds.\n";
		double flops = double(N)*dim*dim*dim / elapsed * 1.0e-9;
		std::cout << "  Performance:     " << flops << " G" << metric << std::endl;
	}
}

int main(int argc, char** argv)
try {
	using namespace std;
	using namespace sw::universal;
	using namespace sw::hprblas;
	//using namespace mtl;  both sw::universal and mtl:: contain a size() operator, we can break that ambiguity by not pulling in mtl:: explicitely

	// a 32-bit float and a <27,1> posit have the same number of significand bits around 1.0
	constexpr size_t nbits = 16;
	constexpr size_t es = 1;
	constexpr size_t capacity = 10;

	using IEEEType=float;
	cout << dynamic_range(posit<nbits, es>()) << endl;

	float eps = std::numeric_limits<IEEEType>::epsilon();
	float epsminus = 1.0f - eps;
	float epsplus = 1.0f + eps;
	float max_fexp = std::numeric_limits<IEEEType>::max_exponent;
	float min_fexp = std::numeric_limits<IEEEType>::min_exponent;
	cout << "IEEE float: epsilon " << eps << " min exp " << min_fexp << " max exp " << max_fexp << endl;


	MeasureMatrixMultiplyPerformance<float>("spFLOPS");
	MeasureMatrixMultiplyPerformance<double>("dpFLOPS");
	MeasureMatrixMultiplyPerformance<posit<16, 1>>("p16.1POPS");
	//MeasureMatrixMultiplyPerformance<posit<16, 1>>("p27.2POPS");
	MeasureMatrixMultiplyPerformance<posit<32, 2>>("p32.2POPS");

	mtl::dense_vector<posit<32,2>> p_sixteen(16);
	cout << "size of p_sixteen is " << size(p_sixteen) << endl;

	std::vector<IEEEType> f_sixteen(16);
	cout << "size of f_sixteen is " << size(f_sixteen) << endl;

	return EXIT_SUCCESS;
}
catch (char const* msg) {
	std::cerr << msg << std::endl;
	return EXIT_FAILURE;
}
catch (const sw::universal::posit_arithmetic_exception& err) {
	std::cerr << "Uncaught posit arithmetic exception: " << err.what() << std::endl;
	return EXIT_FAILURE;
}
catch (const sw::universal::quire_exception& err) {
	std::cerr << "Uncaught quire exception: " << err.what() << std::endl;
	return EXIT_FAILURE;
}
catch (const sw::universal::posit_internal_exception& err) {
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
