// factorial.cpp: example program to demonstrate factorials with arbitrary precision number systems
//
// Copyright (C) 2017-2021 Stillwater Supercomputing, Inc.
//
// This file is part of the universal numbers project, which is released under an MIT Open Source license.
#include "common.hpp"
// bring in different number systems
#include <boost/multiprecision/cpp_int.hpp>
#include <universal/number/integer/integer>
#define POSIT_THROW_ARITHMETIC_EXCEPTION 1
#include <universal/number/posit/posit>
// bring in the factorial function
#include <universal/functions/factorial.hpp>

int main() 
try {
	using namespace std;
	using namespace boost::multiprecision;
	using namespace sw::function;
	using namespace sw::universal;

	int N = 30;
	int128_t v = factorial<int128_t>(N);
	std::cout << typeid(v).name() << " : \n" << v << std::endl;
	 
	cpp_int u = factorial<cpp_int>(N);
	std::cout << typeid(u).name() << " : \n" << u << std::endl;

	integer<128> w = factorial< integer<128> >(N);
	std::cout << typeid(w).name() << " : \n" << w << std::endl;

	using Posit64 = posit<64, 3>;
	Posit64 p64 = factorial<Posit64>(N);
	cout << typeid(p64).name() << " : \n" << setprecision(40) << p64 << endl;

	using Posit128 = posit<128, 4>;
	Posit128 p128 = factorial<Posit128>(N);
	cout << typeid(p128).name() << " : \n" << setprecision(40) << p128 << endl;

	using Posit256 = posit<256, 4>;
	Posit256 p256 = factorial<Posit256>(N);
	cout << typeid(p256).name() << " : \n" << setprecision(40) << p256 << endl;

	return EXIT_SUCCESS;
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
