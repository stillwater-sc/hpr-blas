// mv.cpp example program to demonstrate BLAS L2 Matrix-Vector product
//
// Copyright (C) 2017-2019 Stillwater Supercomputing, Inc.
//
// This file is part of the HPR-BLAS project, which is released under an MIT Open Source license.
#include <vector>
#include <hprblas>
#include <vector_utils.hpp>
#include <print_utils.hpp>

int main(int argc, char** argv)
try {
	using namespace std;
	using namespace sw::hprblas;

	const size_t nbits = 8;
	const size_t es = 0;
	const size_t vecSize = 32;

	int nrOfFailedTestCases = 0;

	cout << "DOT product examples" << endl;
	vector<float> x(vecSize), y(vecSize);
	double result;

	randomVectorFillAroundOneEPS(vecSize, x);  //	sampleVector("x", x);
	randomVectorFillAroundOneEPS(vecSize, y);  // 	sampleVector("y", y);
	result = sw::hprblas::dot(vecSize, x, 1, y, 1);
	cout << "DOT product is " << setprecision(20) << result << endl;

	using Posit = sw::unum::posit<nbits, es>;
	vector<Posit> px(vecSize), py(vecSize);
	randomVectorFillAroundOneEPS(vecSize, px);  //	sampleVector("px", px);
	randomVectorFillAroundOneEPS(vecSize, py);  // 	sampleVector("py", py);
	result = (double)sw::hprblas::dot(vecSize, px, 1, py, 1);
	cout << "DOT product is " << setprecision(20) << result << endl;
	sampleVector("px", px);  // <-- currently shows bad conversions....

	Posit p;
	p = 1.001; cout << p << endl;
	p = 0.999; cout << p << endl;

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
