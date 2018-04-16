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
 by an add is used because we have either one or two rounding events.
 
 */

template<typename Ty, size_t nbits = 40>
void test() {
	Ty a(0.1), b(10.0), c(-1.0);
	Ty with_fma, without_fma;
	with_fma = std::fma(a, b, c);
	without_fma = a*b + c;
	std::cout << std::setprecision(nbits - 2);
	std::cout << std::setw(nbits) << a << " * " << std::setw(nbits) << b << " + " << std::setw(nbits) << c << " = " << std::setw(nbits+8) << without_fma << " not fused" << std::endl;
	std::cout << std::setw(nbits) << a << " * " << std::setw(nbits) << b << " + " << std::setw(nbits) << c << " = " << std::setw(nbits+8) << with_fma << "     fused" << std::endl;
}

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



 // generate specific test case that you can trace with the trace conditions in posit.h
 // for most bugs they are traceable with _trace_conversion and _trace_sub
template<size_t nbits, size_t es, typename Ty>
void GenerateFMATestCase(Ty a, Ty b, Ty c) {
	Ty with_fma, without_fma;
	sw::unum::posit<nbits, es> pa, pb, pc, pref, pfma;
	pa = a;
	pb = b;
	pc = c;
	with_fma = std::fma(a, b, c);
	without_fma = a*b + c;
	pref = with_fma;
	pfma = sw::unum::fma(pa, pb, pc);
	std::cout << "posit<" << nbits << "," << es << ">" << std::endl;
	std::cout << std::setprecision(nbits - 2);
	std::cout << std::setw(nbits) << a << " * " << std::setw(nbits) << b << " + " << std::setw(nbits) << c << " = " << std::setw(nbits) << without_fma << " not fused" << std::endl;
	std::cout << std::setw(nbits) << a << " * " << std::setw(nbits) << b << " + " << std::setw(nbits) << c << " = " << std::setw(nbits) << with_fma    << "     fused" << std::endl;
	std::cout << std::setw(nbits) << pa << " * " << std::setw(nbits) << pb << " + " << std::setw(nbits) << pc << " = " << std::setw(nbits) << pfma << std::endl;
	std::cout << pa.get() << " * " << pb.get() << " + " << pc.get() << " = " << pfma.get() << " (reference: " << pref.get() << ")  ";
	std::cout << (pref == pfma ? "PASS" : "FAIL") << std::endl << std::endl;
	std::cout << std::setprecision(5);
}

int main(int argc, char** argv)
try {
	//using namespace std;
	using namespace mtl;
	//using namespace sw::unum;
	using namespace sw::hprblas;

	int nrOfFailedTestCases = 0;
	//GenerateFMATestCase<32 ,2, double>(0.1, 10.0, -1.0);
	//GenerateFMATestCase<48, 2, double>(0.1, 10.0, -1.0);
	//GenerateFMATestCase<56, 2, double>(0.1, 10.0, -1.0);

	const size_t nbits = 32;
	const size_t es = 2;
	const size_t vecSize = 32;

	{
		using Vector = mtl::dense_vector<double, mtl::vec::parameters<tag::row_major> >;
		axpy_test<double,Vector>("Double AXPY is ", vecSize, 10.0, -1.0);
		test<double>();
	}
	{
		using Vector = mtl::dense_vector<sw::unum::posit<nbits, es>, mtl::vec::parameters<tag::row_major> >;
		axpy_test<sw::unum::posit<nbits, es>,Vector>("posit<32,2> AXPY is ", vecSize, sw::unum::posit<nbits, es>(10.0), sw::unum::posit<nbits, es>(-1.0));
		test< sw::unum::posit<32, 2> >();
	}


	return (nrOfFailedTestCases > 0 ? EXIT_FAILURE : EXIT_SUCCESS);
}
catch (char const* msg) {
	std::cerr << msg << std::endl;
	return EXIT_FAILURE;
}
catch (...) {
	std::cerr << "Caught unknown exception" << std::endl;
	return EXIT_FAILURE;
}
