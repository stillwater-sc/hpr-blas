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

int main() {
	using namespace boost::multiprecision;

	int128_t v = factorial<int128_t>(20);
	std::cout << v << std::endl;

	cpp_int u = factorial<cpp_int>(100);
	std::cout << u << std::endl;

	return 0;
}
