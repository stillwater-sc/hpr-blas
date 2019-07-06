// l3_fused_mv.cpp example program to demonstrate BLAS L3 Reproducible Matrix-Matrix product
//
// Copyright (C) 2017-2019 Stillwater Supercomputing, Inc.
//
// This file is part of the universal numbers project, which is released under an MIT Open Source license.
#include "common.hpp"
#include <boost/multiprecision/cpp_int.hpp>
#define POSIT_THROW_ARITHMETIC_EXCEPTION 1
#include <posit>

template<typename Ty>
Ty factorial(unsigned fact) {
	Ty v = 1;
	for (unsigned i = 2; i <= fact; ++i) {
		v *= i;
	}
	return v;
}

int main() 
try {
	using namespace std;
	using namespace boost::multiprecision;
	using namespace sw::unum;

	int128_t v = factorial<int128_t>(20);
	std::cout << v << std::endl;

	cpp_int u = factorial<cpp_int>(20);
	std::cout << u << std::endl;

	using Posit = posit<64, 3>;
	Posit p = factorial<Posit>(20);
	cout << setprecision(27) << p << endl;

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
