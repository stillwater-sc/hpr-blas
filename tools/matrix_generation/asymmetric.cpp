// asymmetric.cpp: asymmetric matrix generation 
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
	using namespace sw::universal;
	using namespace sw::hprblas;

	int nrOfFailedTestCases = 0;
	std::streamsize prec = std::cout.precision();

	constexpr size_t nbits = 32;
	constexpr size_t es = 2;

	using Scalar = posit<nbits, es>;
//	using Scalar = double;
	using Matrix = mtl::mat::dense2D< Scalar >;

	constexpr size_t N = 5;
	Matrix Q(N,N);

	uniform_random_orthogonal_Heiberger(Q);
	cout << Q << endl;

	Scalar pminpos;
	cout << "minpos<32,2> = " << minpos(pminpos) << endl;
	cout << "minpos<32,2> = " << setw(52) << setprecision(52) << std::fixed << pminpos << endl;

	cout << Q << endl;

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
