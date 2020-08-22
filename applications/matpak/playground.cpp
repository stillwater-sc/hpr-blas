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

/*
		this fails:
		template<typename Matrix>
std::pair<mtl::mat::dense2D<bool>, bool> eq(const Matrix& A, const Matrix& B, const typename Matrix::value::type tolerance = 0.00000000000001) {
}
with error
error C2893: Failed to specialize function template 'std::pair<mtl::mat::dense2D<bool,mtl::mat::parameters<mtl::tag::row_major,mtl::index::c_index,mtl::non_fixed::dimensions,false,size_t>>,bool> eq(const Matrix &,const Matrix &,const Matrix::value::type)'
note: With the following template arguments:
note: 'Matrix=Matrix'
error C3536: 'result': cannot be used before it is initialized

 */
// check if two matrices are the same
template<typename Matrix>
std::pair<mtl::mat::dense2D<bool>, bool> eq(const Matrix& A, const Matrix& B, const double tolerance = 0.00000000000001) {
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

	using Scalar = posit<64,2>;
	/// define a col_major matrix parameter
	typedef mtl::mat::parameters<col_major> column_order;
	using Matrix = mtl::dense2D<Scalar, column_order > ;

	constexpr size_t m = 10;
	constexpr size_t n = 10;
	Matrix A(m, n), B(m, n);
	A = 1;
	double tolerance = 0.0001;
	B = 1.0 + 10 * tolerance;
	cout << A << endl;
	cout << B << endl;

	auto column = A[iall];
	cout << typeid(column).name() << endl;
//	cout << column << endl;
	
	auto result = eq(A, B, tolerance);
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
