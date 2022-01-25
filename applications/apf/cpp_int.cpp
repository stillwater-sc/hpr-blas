

#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <iostream>
#include <universal/native/ieee754.hpp>
#include <universal/number/cfloat/cfloat.hpp>


int main()
{
	using namespace sw::universal;
	using namespace boost::multiprecision;

	auto oldPrecision = std::cout.precision();
	constexpr std::streamsize newPrecision = 10;

	{
		std::int64_t a, b, c, d;

		a = 0x7fffffffffffffff;
		b = 1;
		c = a + b;

		std::cout << a << " + " << b << " = " << c << " " << to_binary(c) << '\n';
	}
	{
		cpp_int a, b, c;

		a = 0x7fffffffffffffff;
		b = 1;
		c = a + b;

		std::cout << a << " + " << b << " = " << c << '\n';
	}


	return EXIT_SUCCESS;
}
