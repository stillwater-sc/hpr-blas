

#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <iostream>
#include <universal/native/ieee754.hpp>
#include <universal/number/cfloat/cfloat.hpp>

const char* pistr = "3.141592653589793238462643383279";
constexpr float pi = 3.141592653589793238462643383279;
//                    3.14159265358979311599796346854
//                         B   8   4   |  = 15 digits of accuracy 

int main()
{
	using namespace sw::universal;
	using namespace boost::multiprecision;

	auto oldPrecision = std::cout.precision();
	constexpr std::streamsize newPrecision = 30;
	{
		float spnat = pi;
		cfloat<32, 8, uint32_t, true> spcf = pi;
		cpp_bin_float_single spmp = pi;
		std::cout << std::setprecision(newPrecision);
		std::cout << "single precision : " << to_binary(spnat) << '\n';
		std::cout << std::fixed << spcf << '\n';
		std::cout << "single precision : " << spmp << '\n';
		std::cout << std::left << spmp << '\n';
		std::cout << std::setw(newPrecision) << std::string(pistr) << '\n';
	}

	{
		cpp_bin_float_quad qp(pistr);
		std::cout << "quad precision   : " << qp << '\n';
		std::cout << std::setprecision(newPrecision);
		std::cout << std::setw(30) << qp << '\n';
		std::cout << std::setw(30) << std::string(pistr) << '\n';
	}
	std::cout << std::setprecision(oldPrecision);
	return 0;

	{
		// Operations at fixed precision and full numeric_limits support:
		cpp_bin_float_100 b = 2;
		std::cout << std::numeric_limits<cpp_bin_float_100>::digits << std::endl;
		std::cout << std::numeric_limits<cpp_bin_float_100>::digits10 << std::endl;
		// We can use any C++ std lib function, lets print all the digits as well:
		std::cout << std::setprecision(std::numeric_limits<cpp_bin_float_100>::max_digits10)
			<< log(b) << std::endl; // print log(2)
									// We can also use any function from Boost.Math:
		std::cout << boost::math::tgamma(b) << std::endl;
		// These even work when the argument is an expression template:
		std::cout << boost::math::tgamma(b * b) << std::endl;
		// And since we have an extended exponent range we can generate some really large 
		// numbers here (4.0238726007709377354370243e+2564):
		std::cout << boost::math::tgamma(cpp_bin_float_100(1000)) << std::endl;
	}


	return 0;
}