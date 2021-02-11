// matrix_types.cpp: examples of different matrix types
//
// Copyright (C) 2017-2019 Stillwater Supercomputing, Inc.
//
// This file is part of the HPR-BLAS project, which is released under an MIT Open Source license.

// warning C4996: 'std::copy::_Unchecked_iterators::_Deprecate': Call to 'std::copy' with parameters that may be unsafe - this call relies on the caller to check that the passed values are correct.
// \mtl4\boost/numeric/mtl/operation/update.hpp(159): warning C4244: 'argument': conversion from 'const double' to 'float', possible loss of data
#pragma warning( disable : 4996 4244)
#include "common.hpp"
#include <hprblas>

template <typename Matrix>
void fill_and_print(Matrix& A, char name)
{
	// Set values in traditional way
	A = 1.2, 3.4,
		5.6, 7.8;

	// Just print them
	std::cout << name << " is \n" << A << "\n";
}

int main(int argc, char** argv)
try {
	using namespace sw::universal;
	using namespace sw::hprblas;

#if defined(MTL_WITH_VARIADIC_TEMPLATE) && defined(MTL_WITH_TEMPLATE_ALIAS)
	using namespace mtl;

	// Compressed matrix
	matrix<float, compressed>                 A(2, 2);
	fill_and_print(A, 'A');

	// Banded matrix
	matrix<float, sparse, banded>             B(2, 2);
	fill_and_print(B, 'B');

	// Matrix in the ELLPACK format 
	matrix<double, ellpack>                   C(2, 2);
	fill_and_print(C, 'C');

	// Coordinate matrix
//	matrix<float, coordinate>                 D(2, 2);
//	fill_and_print(D, 'D');

	// Morton-order matrix with the default mask
	matrix<double, morton>                    E(2, 2);
	fill_and_print(E, 'E');

	// Matrix with a Morton mask is of course a Morton-order matrix
	matrix<double, mask<shark_z_64_row_mask>> F(2, 2);
	fill_and_print(F, 'F');
#endif

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
