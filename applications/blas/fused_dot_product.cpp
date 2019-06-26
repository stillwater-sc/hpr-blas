// fused-dot-product.cpp example program showing a fused-dot product for error free linear algebra
//
// Copyright (C) 2017 Stillwater Supercomputing, Inc.
//
// This file is part of the universal numbers project, which is released under an MIT Open Source license.

// enable the mathematical constants in cmath: old-style preprocessor magic which isn't best practice anymore
#define _USE_MATH_DEFINES
#include "common.hpp"

#include <vector>
#include <posit>

constexpr double pi = 3.14159265358979323846;  // best practice for C++

int main(int argc, char** argv)
try {
	using namespace std;
	using namespace sw::unum;

	constexpr size_t nbits = 16;
	constexpr size_t es = 1;
	constexpr size_t vecSize = 32;

	int nrOfFailedTestCases = 0;

	posit<nbits, es> p;
	vector< posit<nbits,es> > sinusoid(vecSize), cosinusoid(vecSize);

	for (int i = 0; i < vecSize; i++) {
		p = sin( (float(i) / float(vecSize)) *2.0 * pi);
		sinusoid[i] = p;
		p = cos((float(i) / float(vecSize)) *2.0 * pi);
		cosinusoid[i] = p;
	}

	// dot product
	posit<nbits, es> dot_product;
	dot_product = 0.0f;
	for (int i = 0; i < vecSize; i++) {
		dot_product += sinusoid[i] * cosinusoid[i];
	}

	cout << "Dot product is " << dot_product << endl;

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
