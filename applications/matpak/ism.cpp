// ism.cpp: test to validate that a matrix is an M matrix
//
// Copyright (C) 2017-2020 Stillwater Supercomputing, Inc.
//
// This file is part of the HPR-BLAS project, which is released under an MIT Open Source license.
#include <boost/numeric/mtl/mtl.hpp>
#include <hprblas>
#include <matpak/ism.hpp>

// Selects posits or floats
#define USE_POSIT 1

template<typename Matrix>
bool isMatrixAnMMatrix(const Matrix& A) {
	return false;
}

int main(int argc, char** argv)
try {
	using namespace std;
	using namespace mtl;
	using namespace sw::unum;
	using namespace sw::hprblas;

#if USE_POSIT
	using Ty     = sw::unum::posit<8, 0>;
	using Matrix = mtl::dense2D< Ty >;
	using Vector = mtl::dense_vector< Ty >;
#else
	using Ty     = float;
	using Matrix = mtl::dense2D<float>;
	using Vector = mtl::dense_vector<float>;
#endif

	constexpr int n = 3; // Number of states

	mtl::dense2D<Ty> A(n, n); // System dynamics matrix
    
	A = 1;  // create identity matrix

	if (sw::hprblas::ism(A)) {
		cout << "A is an M-matrix\n" << A << '\n';
	} else {
		cout << "A is not an M-Matrix\n" << A << '\n';
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
