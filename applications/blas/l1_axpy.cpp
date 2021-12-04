// l1_axpy.cpp: example program contrasting a BLAS L1 ?axpy routine between FLOAT and POSIT
//
// Copyright (C) 2017-2021 Stillwater Supercomputing, Inc.
//
// This file is part of the universal numbers project, which is released under an MIT Open Source license.
#include "common.hpp"
// enable posit arithmetic exceptions
#define POSIT_THROW_ARITHMETIC_EXCEPTION 1
#include <universal/number/posit/posit.hpp>

// axpy: a times x plus y
template<typename scale_T, typename vector_T>
void axpy(size_t n, scale_T a, const vector_T& x, size_t incx, vector_T& y, size_t incy) {
	size_t cnt, ix, iy;
	for (cnt = 0, ix = 0, iy = 0; cnt < n && ix < x.size() && iy < y.size(); ++cnt, ix += incx, iy += incy) {
		y[iy] += a * x[ix];
	}
}

// print a vector
template<typename vector_T>
void printStridedVector(std::ostream& ostr, size_t n, vector_T& x, size_t incx = 1) {
	size_t cnt, ix;
	for (cnt = 0, ix = 0; cnt < n && ix < x.size(); ++cnt, ix += incx) {
		cnt == 0 ? ostr << "[" << x[ix] : ostr << ", " << x[ix];
	}
	ostr << "]";
}

int main(int argc, char** argv)
try {
	using namespace sw::universal;

	constexpr size_t nbits = 16;
	constexpr size_t es = 1;
	//constexpr size_t vecSize = 32;

	int nrOfFailedTestCases = 0;

	constexpr int d = 5;
	std::vector< posit<nbits, es> > v1 = { 1.0, 2.0, 3.0, 4.0, 5.0 };
	std::vector< posit<nbits, es> > v2(d);
	posit<nbits, es> alpha(SpecificValue::minpos);

	std::cout << "AXPY is\n";
	printStridedVector(std::cout, d, v1, 1); std::cout << '\n';
	printStridedVector(std::cout, d, v2, 1); std::cout << '\n';

	axpy(d, alpha, v1, 1, v2, 1);

	printStridedVector(std::cout, d, v2, 1); std::cout << '\n';

	return (nrOfFailedTestCases > 0 ? EXIT_FAILURE : EXIT_SUCCESS);
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
catch (const std::runtime_error& err) {
	std::cerr << "Uncaught runtime exception: " << err.what() << std::endl;
	return EXIT_FAILURE;
}
catch (...) {
	std::cerr << "Caught unknown exception" << std::endl;
	return EXIT_FAILURE;
}
