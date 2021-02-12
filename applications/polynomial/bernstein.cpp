// bernstein.cpp: evaluation of Bernstein polynomials
//
// Copyright (C) 2017-2020 Stillwater Supercomputing, Inc.
// Author:
//
// This file is part of the HPR-BLAS project, which is released under an MIT Open Source license.
#include <hprblas>
//#include <polynomials/bernstein.hpp>

int main(int argc, char** argv)
try {
	using namespace std;
	using namespace sw::universal;

	using Scalar = posit<32, 2>;
	using Matrix = mtl::mat::dense2D<Scalar>;

//	Matrix B = bernstein<Matrix>(5);
//	cout << B << endl;

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
