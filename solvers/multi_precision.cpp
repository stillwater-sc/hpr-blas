// multi_precision.cpp : example program comparing float vs posit matrix inversion algorithms
//
// Copyright (C) 2017-2020 Stillwater Supercomputing, Inc.
//
// This file is part of the HPR-BLAS project, which is released under an MIT Open Source license.

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

using namespace sw::universal;

template<typename Scalar>
void Equation() {
	Scalar a{ 1.0 };
	Scalar b(2.0);
	Scalar c = a + b;
	std::cout << c << std::endl;
}

int main(int argc, char** argv)
try {
	using namespace std;
	using namespace mtl;
	using namespace sw::hprblas;
	using namespace boost::multiprecision;

	using sp = boost::multiprecision::cpp_bin_float_single;
	using dp = boost::multiprecision::cpp_bin_float_double;
	using qp = boost::multiprecision::cpp_bin_float_quad;



	Equation<qp>();

	return EXIT_SUCCESS;
}
catch (char const* msg) {
	std::cerr << msg << std::endl;
	return EXIT_FAILURE;
}
catch (const sw::universal::universal_arithmetic_exception& err) {
	std::cerr << "Uncaught universal arithmetic exception: " << err.what() << std::endl;
	return EXIT_FAILURE;
}
catch (const sw::universal::universal_internal_exception& err) {
	std::cerr << "Uncaught universal internal exception: " << err.what() << std::endl;
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
