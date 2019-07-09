// binomial_coefficients.cpp example program to test binomial coefficients for Hilbert matrix generation
//
// Copyright (C) 2017-2019 Stillwater Supercomputing, Inc.
//
// This file is part of the universal numbers project, which is released under an MIT Open Source license.
#include "common.hpp"
#include <boost/multiprecision/cpp_int.hpp>
#include <iostream>
// enable posit arithmetic exceptions
#define POSIT_THROW_ARITHMETIC_EXCEPTION 1
#include <hprblas>
#include <matrix_utils.hpp>

int main() 
try {
	using namespace std;
	using namespace boost::multiprecision;
	using namespace sw::hprblas;

	/*
	int128_t v = factorial<int128_t>(20);
	std::cout << v << std::endl;

	cpp_int u = factorial<cpp_int>(100);
	std::cout << u << std::endl;
	*/

	for (uint64_t n = 1; n < 10; ++n) {
		for (uint64_t k = 0; k <= n; ++k) {
			cout << "Binomial(" << n << "," << k << ") = " << BinomialCoefficient(n, k) << endl;
		}
	}
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
