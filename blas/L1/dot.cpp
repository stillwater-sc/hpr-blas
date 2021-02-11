// dot.cpp: example program contrasting a BLAS L1 ?dot routine between FLOAT and POSIT
//
// Copyright (C) 2017-2019 Stillwater Supercomputing, Inc.
//
// This file is part of the HPR-BLAS project, which is released under an MIT Open Source license.

#include <ratio>
#include <chrono>
#include <iostream>
#include <ctime>

////////////////////////////////////////////////////////////////////////////////////////
///  BEHAVIORAL COMPILATION SWITCHES for posit library configuration

////////////////////////////////////////////////////////////////////////////////////////
// enable/disable special posit format I/O
// POSIT_ROUNDING_ERROR_FREE_IO_FORMAT
// default is to print (long double) values
// #define POSIT_ROUNDING_ERROR_FREE_IO_FORMAT 0

////////////////////////////////////////////////////////////////////////////////////////
// enable/disable the ability to use literals in binary logic and arithmetic operators
// POSIT_ENABLE_LITERALS)
// default is to enable them
// #define POSIT_ENABLE_LITERALS 1

////////////////////////////////////////////////////////////////////////////////////////
// enable throwing specific exceptions for posit arithmetic errors
// left to application to enable
// POSIT_THROW_ARITHMETIC_EXCEPTION
// default is to use NaR as a signalling error
// #define POSIT_THROW_ARITHMETIC_EXCEPTION 0

////////////////////////////////////////////////////////////////////////////////////////
/// INCLUDE FILES posit library
#include <universal/number/posit/posit>

///////////////////////////////////////////////////////////////////////////////////////
/// useful mathematical property functions
//#include <universal/functions/functions.hpp>

///////////////////////////////////////////////////////////////////////////////////////
/// the underlying matrix/vector machinery
#ifndef MTL_WITH_INITLIST
#define MTL_WITH_INITLIST
#endif
#include <boost/numeric/mtl/mtl.hpp>

///////////////////////////////////////////////////////////////////////////////////////
/// the High-Performance Reproducible Basic Linear Algebra Subroutines
/// L1, L2, and L3 matrix/vector operations
#include <hprblas.hpp>
/// norms (l1, l2, linf, Frobenius) using HPR methods
//#include <norms.hpp>

#include <utils/matvec.hpp>

using namespace sw::universal;

int main(int argc, char** argv)
try {
	using namespace std;
	using namespace sw::hprblas;
	// configure the posit environment
	const size_t nbits = 32;
	const size_t es = 2;
	const size_t vecSize = 1024;

	int nrOfFailedTestCases = 0;
	// steady_clock example
	using namespace std::chrono;

	cout << "DOT product examples" << endl;
	vector<float> x(vecSize), y(vecSize);
	float fresult;


	randomVectorFillAroundOneEPS(vecSize, x);  //	sampleVector("x", x);
	randomVectorFillAroundOneEPS(vecSize, y);  // 	sampleVector("y", y);
	fresult = sw::hprblas::dot(vecSize, x, 1, y, 1);
	cout << "DOT product is " << setprecision(20) << fresult << endl;
#ifdef LATER
	using Posit = sw::universal::posit<nbits, es>;
	vector<Posit> px(vecSize), py(vecSize);
	Posit presult;
	randomVectorFillAroundOneEPS(vecSize, px);  //	sampleVector("px", px);
	randomVectorFillAroundOneEPS(vecSize, py);  // 	sampleVector("py", py);

	steady_clock::time_point t1 = steady_clock::now();
	presult = sw::hprblas::dot(vecSize, px, 1, py, 1);
	steady_clock::time_point t2 = steady_clock::now();
	double ops = vecSize * 2.0; // dot product is vecSize products and vecSize adds
	cout << "DOT product is " << setprecision(20) << presult << endl;
	//sampleVector("px", px);  // <-- currently shows bad conversions....

	duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
	double elapsed = time_span.count();
	std::cout << "It took " << elapsed << " seconds." << std::endl;
	std::cout << "Performance " << (uint32_t) (ops / (1000*elapsed)) << " KOPS" << std::endl;
	std::cout << std::endl;
#endif
	return (nrOfFailedTestCases > 0 ? EXIT_FAILURE : EXIT_SUCCESS);
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
