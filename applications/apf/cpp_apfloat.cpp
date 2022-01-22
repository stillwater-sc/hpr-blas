

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

	/* cfloat print behavior you need to emulate:
	native single precision : 0b0.11111110.00000000000000000000000 : 170141183460469231731687303715884105728.0000000000    <-- fixed format takes precision as the number of bits after the fixed point
	cfloat single precision : 0b0.11111110.00000000000000000000000 : 170141183460469231731687303715884105728.00000000000000000000000  <-- you are taking the precision of the mantissa
	boost  single precision :                                      : 170141183460469231731687303715884105728.0000000000
	native single precision : 0b0.11111110.00000000000000000000000 : 1.7014118346e+38                                      <-- scientific format takes precision as the number of bits after the decimal point
	cfloat single precision : 0b0.11111110.00000000000000000000000 : 1.70141e+38                                           <-- you are not honoring precision
	boost  single precision :                                      : 1.7014118346e+38
	*/
	auto oldPrecision = std::cout.precision();
	constexpr std::streamsize newPrecision = 10;
	{
		float v = 1.0f * pow(2.0f, 127.0f);
		float spnat = v;
		cfloat<32, 8, uint32_t, true> spcf = v;
		cpp_bin_float_single spmp = v;
		std::cout << std::setprecision(newPrecision);

		std::cout << std::fixed;
		std::cout << "native single precision : " << to_binary(spnat) << " : " << spnat << '\n';
		std::cout << "cfloat single precision : " << to_binary(spcf) << " : " << spcf << '\n';
		std::cout << "boost  single precision : " << "                                    " << " : " << spmp << '\n';
		std::cout << std::scientific;
		std::cout << "native single precision : " << to_binary(spnat) << " : " << spnat << '\n';
		std::cout << "cfloat single precision : " << to_binary(spcf) << " : " << spcf << '\n';
		std::cout << "boost  single precision : " << "                                    " << " : " << spmp << '\n';

		std::cout << spmp << '\n';
		std::cout << std::setw(newPrecision) << std::string(pistr) << '\n';
	}
	return 0;
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