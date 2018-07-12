// axpy.cpp: example program contrasting a BLAS L1 ?axpy routine between FLOAT and POSIT
//
// Copyright (C) 2017-2018 Stillwater Supercomputing, Inc.
//
// This file is part of the HPR-BLAS project, which is released under an MIT Open Source license.
#include <hprblas>

/*
 An axpy operation, that is, a * X + Y, has resolution-canceling rounding error  
 when the scales of the product and the Y element are disproporitional. 
 
 Reproducibility is challenged when an FMA or regular mul followed
 by an add is used; we have either one or two rounding events.
 
 */

template<typename Element, typename Vector>
void axpy_test(std::string tag, int vecSize, Element x_value, Element y_value)
{
	Vector X(vecSize), Y(vecSize), AXPY(vecSize);
	X = x_value;
	Y = y_value;
	double alpha = 0.1;
	AXPY = alpha*X + Y;

	std::cout << tag << std::endl;
	std::cout << AXPY << std::endl;
}


int main(int argc, char** argv)
try {
	//using namespace std;
	using namespace mtl;
	//using namespace sw::unum;
	using namespace sw::hprblas;

	int nrOfFailedTestCases = 0;

	const size_t nbits = 32;
	const size_t es = 2;
	const size_t vecSize = 32;

	{
		using Vector = mtl::dense_vector<double, mtl::vec::parameters<tag::row_major> >;
		axpy_test<double,Vector>("Double AXPY is ", vecSize, 10.0, -1.0);
	}
	{
		using Vector = mtl::dense_vector<sw::unum::posit<nbits, es>, mtl::vec::parameters<tag::row_major> >;
		axpy_test<sw::unum::posit<nbits, es>,Vector>("posit<32,2> AXPY is ", vecSize, sw::unum::posit<nbits, es>(10.0), sw::unum::posit<nbits, es>(-1.0));
	}

	return (nrOfFailedTestCases > 0 ? EXIT_FAILURE : EXIT_SUCCESS);
}
catch (char const* msg) {
	std::cerr << msg << std::endl;
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
