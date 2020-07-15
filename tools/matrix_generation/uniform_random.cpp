// uniform_random.cpp: uniform random matrix generator test
//
// Copyright (C) 2017-2020 Stillwater Supercomputing, Inc.
//
// This file is part of the HPRBLAS project, which is released under an MIT Open Source license.
#define POSIT_THROW_ARITHMETIC_EXCEPTION 1
#include <hprblas>
// matrix generators
#include "generators/matrix_generators.hpp"

int main(int argc, char** argv)
try {
	using namespace std;
	using namespace mtl;
	using namespace sw::unum;
	using namespace sw::hprblas;

	int nrOfFailedTestCases = 0;
	std::streamsize prec = std::cout.precision();

	constexpr size_t nbits = 32;
	constexpr size_t es = 2;

	using Scalar = posit<nbits, es>;
//	using Scalar = double;
	using Matrix = mtl::mat::dense2D< Scalar >;

	constexpr size_t m = 5;
	constexpr size_t n = 4;
	Matrix A(m,n);

	uniform_rand(A);
	cout << "uniform random        :\n" << A << endl;

	auto B = uniform_rand<Scalar>(m, n, -1.0, 1.0);
	cout << "uniform random(-1,1)  :\n" << B << endl;

	uniform_rand_sorted(A);
	cout << "uniform random sorted :\n" << A << endl;

	return (nrOfFailedTestCases > 0 ? EXIT_FAILURE : EXIT_SUCCESS);
}
catch (char const* msg) {
	std::cerr << "caught ad hoc exception: " << msg << std::endl;
	return EXIT_FAILURE;
}
catch (posit_arithmetic_exception& err) {
	std::cerr << "uncaught posit arithmetic exception: " << err.what() << std::endl;
	return EXIT_FAILURE;
}
catch (mtl::runtime_error& err) {
	std::cerr << "caught MTL run-time exception: " << err.what() << std::endl;
	return EXIT_FAILURE;
}
catch (...) {
	std::cerr << "caught unknown exception" << std::endl;
	return EXIT_FAILURE;
}
