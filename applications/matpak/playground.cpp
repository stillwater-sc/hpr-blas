// playground.cpp: test program to experiment with different matpak functions
//
// Copyright (C) 2017-2020 Stillwater Supercomputing, Inc.
//
// This file is part of the HPR-BLAS project, which is released under an MIT Open Source license.
#include <iostream>
#include <utility>
// #define MTL_WITH_INITLIST
// #define MTL_WITH_AUTO
// #define MTL_WITH_RANGEDFOR
#include <hprblas>

// check if two matrices are the same
template<typename Matrix>
std::pair<mtl::mat::dense2D<bool>, bool> eq(const Matrix& A, const Matrix& B) {
	typedef typename Matrix::value_type value_type;
	
	size_t ar = num_rows(A);
	size_t ac = num_cols(A);
	size_t br = num_rows(B);
	size_t bc = num_cols(B);

	if (ar != br || ac != bc) {
		return std::pair<mtl::mat::dense2D<bool>, bool>(mtl::mat::dense2D<bool>() , false);
	}

	mtl::mat::dense2D<bool> T(ar, ac);
	for (size_t i = 0; i < ar; ++i) {
		for (size_t j = 0; j < ac; ++j) {
			if (abs(A[i][j] - B[i][j]) > value_type(0.0001)) {
				T[i][j] = false;
			} 
			else {
				T[i][j] = true;
			}
		}
	}
	return std::pair<mtl::mat::dense2D<bool>, bool>(T, true);
}

int main(int argc, char** argv)
try {
	using namespace std;
	using namespace mtl;
	using namespace sw::unum;
	using namespace sw::hprblas;

	/// Short-cut to define parameters with unsigned and defaults otherwise
	typedef mtl::mat::parameters<col_major, index::c_index, mtl::non_fixed::dimensions, false, unsigned> column_matrix;
	using Matrix = mtl::dense2D<float, column_matrix > ;

	Matrix A(2, 2), B(2, 2);

	A.elements();
	auto result = eq(A, B);
	if (result.second) {
		cout << result.first << endl;
	}
	else {
		cout << "matrices are not same dimension\n";
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
