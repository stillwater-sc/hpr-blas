// asymmetric.cpp: asymmetric matrix generation 
//
// Copyright (C) 2017-2018 Stillwater Supercomputing, Inc.
//
// This file is part of the HPR-BLAS project, which is released under an MIT Open Source license.
#include <hprblas>

int main(int argc, char** argv)
try {
	using namespace std;
	using namespace mtl;
	using namespace sw::unum;
	using namespace sw::hprblas;

	int nrOfFailedTestCases = 0;
	std::streamsize prec = std::cout.precision();

	constexpr size_t nbits = 32;
	constexpr size_t es = 1;

	using Scalar = posit<nbits, es>;
//	using Scalar = float;
	using Matrix = mtl::mat::dense2D< Scalar >;

	constexpr size_t N = 10;
	Matrix Q(N,N);

	uniform_random_orthogonal_Heiberger(Q);
	cout << Q << endl;

	cout << "minpos<32,2> = " << minpos<32, 2>() << endl;
	cout << "minpos<32,2> = " << setw(52) << setprecision(52) << std::fixed << minpos<32, 2>() << endl;

	cout << Q[0][0] << endl;
	cout << Q[0][1] << endl;
	cout << Q[1][0] << endl;
	cout << Q[1][1] << endl;

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
