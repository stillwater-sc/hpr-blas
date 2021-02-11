// l1_asum.cpp: example program contrasting a BLAS L1 ?asum routine between FLOAT and POSIT
//
// Copyright (C) 2017-2019 Stillwater Supercomputing, Inc.
//
// This file is part of the universal numbers project, which is released under an MIT Open Source license.

// place configuration compile switches before including hprblas
#include <hprblas>

using namespace std;
using namespace sw::universal;

template<typename Real, typename Vector>
void myFunction(int N) {
	Vector v(N);
	v = Real(1.0);
	cout << 5 * v << endl;
}

int main(int argc, char** argv)
try {
	constexpr size_t nbits = 20;
	constexpr size_t es = 1;

	int nrOfFailedTestCases = 0;

	cout << "ASUM is " << 1 << endl;

	using namespace mtl;
	using Real = posit<nbits, es>;
	using Vector = mtl::dense_vector<Real>;

	myFunction<Real, Vector>(10);
	myFunction<double, Vector>(10);

	return (nrOfFailedTestCases > 0 ? EXIT_FAILURE : EXIT_SUCCESS);
}
catch (char const* msg) {
	cerr << msg << endl;
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
	cerr << "Caught unknown exception" << endl;
	return EXIT_FAILURE;
}
