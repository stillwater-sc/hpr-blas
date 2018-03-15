// axpy.cpp: example program contrasting a BLAS L1 ?axpy routine between FLOAT and POSIT
//
// Copyright (C) 2017-2018 Stillwater Supercomputing, Inc.
//
// This file is part of the HPR-BLAS project, which is released under an MIT Open Source license.

#include <vector>
#include <boost/numeric/mtl/mtl.hpp>
#include <hprblas>

/*
 An axpy operation, that is, a * X + Y, has resolution-canceling rounding error when the scales of the 
 product and the Y element are disproporitional. Reproducibility is challenged when a FMA or regular mul followed
 by an add is used because the rounding is either pre or post multiply.
 
 */

int main(int argc, char** argv)
try {
	using namespace std;
	using namespace mtl;
	using namespace sw::unum;
	using namespace sw::hprblas;

	const size_t nbits = 16;
	const size_t es = 1;
	const size_t vecSize = 32;

	int nrOfFailedTestCases = 0;

	dense_vector<posit<nbits,es>, mtl::vec::parameters<tag::row_major> > x(vecSize), y(vecSize), axpy(vecSize);
	x = posit<nbits,es>(10.0);
	y = posit<nbits,es>(-1.0);
	float a = 0.1f;
	axpy = a*x + y;

	cout << "AXPY is " << endl;
	cout << axpy << endl;

	return (nrOfFailedTestCases > 0 ? EXIT_FAILURE : EXIT_SUCCESS);
}
catch (char const* msg) {
	cerr << msg << endl;
	return EXIT_FAILURE;
}
catch (...) {
	cerr << "Caught unknown exception" << endl;
	return EXIT_FAILURE;
}
