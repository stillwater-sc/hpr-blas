// symmetric.cpp: symmetric matrix generation 
//
// Copyright (C) 2017-2018 Stillwater Supercomputing, Inc.
//
// This file is part of the HPR-BLAS project, which is released under an MIT Open Source license.
#include <hprblas>
#include "matrix_utils.hpp"



int main(int argc, char** argv)
try {
	using namespace std;
	using namespace mtl;
	using namespace sw::unum;
	using namespace sw::hprblas;

	int nrOfFailedTestCases = 0;

	constexpr size_t nbits = 16;
	constexpr size_t es = 1;

	using Matrix = mtl::mat::dense2D< posit<nbits, es> >;

	constexpr size_t N = 10;
	Matrix A(N, N);

	// How do I make that matrix symmetric as type?
	uniform_rand_diagonally_dominant(A);

	cout << A << endl;

	return (nrOfFailedTestCases > 0 ? EXIT_FAILURE : EXIT_SUCCESS);
}
catch (char const* msg) {
	std::cerr << "caught ad hoc exception: " << msg << std::endl;
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
