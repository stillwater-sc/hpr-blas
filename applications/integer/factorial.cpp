#include "common.hpp"

#include <boost/multiprecision/cpp_int.hpp>
#include <iostream>

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
	using namespace boost::multiprecision;

	int128_t v = factorial<int128_t>(20);
	std::cout << v << std::endl;

	cpp_int u = factorial<cpp_int>(100);
	std::cout << u << std::endl;

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
