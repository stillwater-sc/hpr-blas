// fam.cpp: example program contrasting fused accumulate-multiply functionality between IEEE and POSIT
//
// Copyright (C) 2017-2019 Stillwater Supercomputing, Inc.
//
// This file is part of the HPR-BLAS project, which is released under an MIT Open Source license.


#include <universal/number/posit/posit.hpp>

 // generate specific test case that you can trace with the trace conditions in posit.h
 // for most bugs they are traceable with _trace_conversion and _trace_sub
// Fused Add-Multiply = (a + b)*c
template<size_t nbits, size_t es, typename Ty>
void GenerateFAMTestCase(Ty a, Ty b, Ty c) {
#if 0
	Ty with_fam, without_fam;
	sw::universal::posit<nbits, es> pa(a), pb(b), pc(c), pref, pfam;

	++pa;         // perturb a by adding 1 machine epsilon
	a = Ty(pa);   // and update the input to reflect the new value
	pb = b;
	pc = c;
	// (a + b)*c = a*c + b*c = fma(b,c,a*c)
	with_fam = std::fmal(b, c, long double(a)*long double(c));
	without_fam = (a + b) * c;
	pref = with_fam;
	pfam = sw::universal::fam(pa, pb, pc);
	std::cout << "posit<" << nbits << "," << es << ">" << std::endl;
	std::cout << std::setprecision(nbits - 2);
	std::cout << std::setw(nbits) << a << " * " << std::setw(nbits) << b << " + " << std::setw(nbits) << c << " = " << std::setw(nbits) << without_fam << " not fused" << std::endl;
	std::cout << std::setw(nbits) << a << " * " << std::setw(nbits) << b << " + " << std::setw(nbits) << c << " = " << std::setw(nbits) << with_fam    << "     fused" << std::endl;
	std::cout << std::setw(nbits) << pa << " * " << std::setw(nbits) << pb << " + " << std::setw(nbits) << pc << " = " << std::setw(nbits) << pfam << std::endl;
	std::cout << pa.get() << " * " << pb.get() << " + " << pc.get() << " = " << pfam.get() << " (reference: " << pref.get() << ")  ";
	std::cout << (pref == pfam ? "PASS" : "FAIL") << std::endl << std::endl;
	std::cout << std::setprecision(5);
#endif

}

int main(int argc, char** argv)
try {
	using namespace std;
	//using namespace sw::universal;   it is more informative to use the namespace explicitely

	int nrOfFailedTestCases = 0;

	// NOTE: pick values that have an exact binary representation
	// otherwise you will be fighting round-off error getting into the calculation
	// which muddles the actual round-off you are trying to quantify
	GenerateFAMTestCase<32 ,2, double>(0.125, 10.0, -1.25);
	GenerateFAMTestCase<48, 2, double>(0.125, 10.0, -1.25);
	GenerateFAMTestCase<56, 2, double>(0.125, 10.0, -1.25);

	constexpr size_t nbits = 32;
	constexpr size_t es = 2;

	return (nrOfFailedTestCases > 0 ? EXIT_FAILURE : EXIT_SUCCESS);
}
catch (char const* msg) {
	std::cerr << msg << std::endl;
	return EXIT_FAILURE;
}
catch (const sw::universal::posit_arithmetic_exception& err) {
	std::cerr << "Uncaught posit arithmetic exception: " << err.what() << std::endl;
	return EXIT_FAILURE;
}
catch (const sw::universal::quire_exception& err) {
	std::cerr << "Uncaught quire exception: " << err.what() << std::endl;
	return EXIT_FAILURE;
}
catch (const sw::universal::posit_internal_exception& err) {
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
