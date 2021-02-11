// log.cpp: testing the use of C++ std lib with custom types
//
// Copyright (C) 2017-2020 Stillwater Supercomputing, Inc.
//
// This file is part of the universal numbers project, which is released under an MIT Open Source license.
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <hprblas>

// define a true 256-bit IEEE floating point type
constexpr size_t bits_in_octand = 113 + 128;
using cpp_bin_float_octand = boost::multiprecision::number<boost::multiprecision::backends::cpp_bin_float<bits_in_octand, boost::multiprecision::backends::digit_base_2, void, boost::int16_t, -16382, 16383>, boost::multiprecision::expression_template_option::et_off>;
// define the floating point types (single, double, quad, octand)
using sp = boost::multiprecision::cpp_bin_float_single;
using dp = boost::multiprecision::cpp_bin_float_double;
using qp = boost::multiprecision::cpp_bin_float_quad;
using op = cpp_bin_float_octand;

template<typename Real>
void dumpVector(std::vector<Real>& v) {
	for (auto value : v) {
		std::cout << value << std::endl;
	}
}

template<typename Real>
void dumpPair(std::vector<Real>& args, std::vector<Real>& nodes) {
	for (size_t i = 0; i < args.size(); ++i) {
		std::cout << args[i] << " : " << nodes[i] << std::endl;
	}
}

int main(int argc, char** argv)
try {
	using namespace std;
	using namespace sw::universal;

	using Real = qp;
	cout << setprecision(numeric_limits<qp>::digits10);

	using namespace boost::multiprecision;
	cpp_dec_float_100 b = 2;
	cout << std::numeric_limits<cpp_dec_float_100>::digits << endl;
	cout << std::numeric_limits<cpp_dec_float_100>::digits10 << endl;

	// any C++ std lib function can be used
	cpp_dec_float_100 c = log(b);
	cout << setprecision(numeric_limits<cpp_dec_float_100>::max_digits10)
		<< log(b) << endl;

	// we can also use any function from boost.math
	cout << boost::math::tgamma(b) << endl;

	// some really large numbers
	cout << boost::math::tgamma(cpp_dec_float_100(1000)) << endl;

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


/* 
 BOOST implementation of operator<<() for cpp_bin_float backend

 class number {
 ...
   //
   // String conversion functions:
   //
   std::string str(std::streamsize digits = 0, std::ios_base::fmtflags f = std::ios_base::fmtflags(0))const
   {
	  return m_backend.str(digits, f);
   }
   template<class Archive>
   void serialize(Archive & ar, const unsigned int version)
   {
   ar & m_backend;
   }
...
}

template <class Backend, expression_template_option ExpressionTemplates>
inline std::ostream& operator << (std::ostream& os, const number<Backend, ExpressionTemplates>& r)
{
	std::streamsize d = os.precision();
	std::string s = r.str(d, os.flags());
	std::streamsize ss = os.width();
	if (ss > static_cast<std::streamsize>(s.size()))
	{
		char fill = os.fill();
		if ((os.flags() & std::ios_base::left) == std::ios_base::left)
			s.append(static_cast<std::string::size_type>(ss - s.size()), fill);
		else
			s.insert(static_cast<std::string::size_type>(0), static_cast<std::string::size_type>(ss - s.size()), fill);
	}
	return os << s;
}

 template <unsigned Digits, digit_base_type DigitBase, class Allocator, class Exponent, Exponent MinE, Exponent MaxE>
std::string cpp_bin_float<Digits, DigitBase, Allocator, Exponent, MinE, MaxE>::str(std::streamsize dig, std::ios_base::fmtflags f) const
{
   if(dig == 0)
	  dig = std::numeric_limits<number<cpp_bin_float<Digits, DigitBase, Allocator, Exponent, MinE, MaxE> > >::max_digits10;

   bool scientific = (f & std::ios_base::scientific) == std::ios_base::scientific;
   bool fixed = !scientific && (f & std::ios_base::fixed);

   std::string s;

   if(exponent() <= cpp_bin_float<Digits, DigitBase, Allocator, Exponent, MinE, MaxE>::max_exponent)
   {
	  // How far to left-shift in order to demormalise the mantissa:
	  boost::intmax_t shift = (int)cpp_bin_float<Digits, DigitBase, Allocator, Exponent, MinE, MaxE>::bit_count - exponent() - 1;
	  boost::intmax_t digits_wanted = static_cast<int>(dig);
	  boost::intmax_t base10_exp = exponent() >= 0 ? static_cast<boost::intmax_t>(std::floor(0.30103 * exponent())) : static_cast<boost::intmax_t>(std::ceil(0.30103 * exponent()));
	  //
	  // For fixed formatting we want /dig/ digits after the decimal point,
	  // so if the exponent is zero, allowing for the one digit before the
	  // decimal point, we want 1 + dig digits etc.
	  //
	  if(fixed)
		 digits_wanted += 1 + base10_exp;
	  if(scientific)
		 digits_wanted += 1;
	  if(digits_wanted < -1)
	  {
		 // Fixed precision, no significant digits, and nothing to round!
		 s = "0";
		 if(sign())
			s.insert(static_cast<std::string::size_type>(0), 1, '-');
		 boost::multiprecision::detail::format_float_string(s, base10_exp, dig, f, true);
		 return s;
	  }
	  //
	  // power10 is the base10 exponent we need to multiply/divide by in order
	  // to convert our denormalised number to an integer with the right number of digits:
	  //
	  boost::intmax_t power10 = digits_wanted - base10_exp - 1;
	  //
	  // If we calculate 5^power10 rather than 10^power10 we need to move
	  // 2^power10 into /shift/
	  //
	  shift -= power10;
	  cpp_int i;
	  int roundup = 0; // 0=no rounding, 1=tie, 2=up
	  static const unsigned limb_bits = sizeof(limb_type) * CHAR_BIT;
	  //
	  // Set our working precision - this is heuristic based, we want
	  // a value as small as possible > cpp_bin_float<Digits, DigitBase, Allocator, Exponent, MinE, MaxE>::bit_count to avoid large computations
	  // and excessive memory usage, but we also want to avoid having to
	  // up the computation and start again at a higher precision.
	  // So we round cpp_bin_float<Digits, DigitBase, Allocator, Exponent, MinE, MaxE>::bit_count up to the nearest whole number of limbs, and add
	  // one limb for good measure.  This works very well for small exponents,
	  // but for larger exponents we add a few extra limbs to max_bits:
	  //
#ifdef BOOST_MP_STRESS_IO
	  boost::intmax_t max_bits = cpp_bin_float<Digits, DigitBase, Allocator, Exponent, MinE, MaxE>::bit_count + 32;
#else
	  boost::intmax_t max_bits = cpp_bin_float<Digits, DigitBase, Allocator, Exponent, MinE, MaxE>::bit_count + ((cpp_bin_float<Digits, DigitBase, Allocator, Exponent, MinE, MaxE>::bit_count % limb_bits) ? (limb_bits - cpp_bin_float<Digits, DigitBase, Allocator, Exponent, MinE, MaxE>::bit_count % limb_bits) : 0) + limb_bits;
	  if(power10)
		 max_bits += (msb(boost::multiprecision::detail::abs(power10)) / 8) * limb_bits;
#endif
	  do
	  {
		 boost::int64_t error = 0;
		 boost::intmax_t calc_exp = 0;
		 //
		 // Our integer result is: bits() * 2^-shift * 5^power10
		 //
		 i = bits();
		 if(shift < 0)
		 {
			if(power10 >= 0)
			{
			   // We go straight to the answer with all integer arithmetic,
			   // the result is always exact and never needs rounding:
			   BOOST_ASSERT(power10 <= (boost::intmax_t)INT_MAX);
			   i <<= -shift;
			   if(power10)
				  i *= pow(cpp_int(5), static_cast<unsigned>(power10));
			}
			else if(power10 < 0)
			{
			   cpp_int d;
			   calc_exp = boost::multiprecision::cpp_bf_io_detail::restricted_pow(d, cpp_int(5), -power10, max_bits, error);
			   shift += calc_exp;
			   BOOST_ASSERT(shift < 0); // Must still be true!
			   i <<= -shift;
			   cpp_int r;
			   divide_qr(i, d, i, r);
			   roundup = boost::multiprecision::cpp_bf_io_detail::get_round_mode(r, d, error, i);
			   if(roundup < 0)
			   {
#ifdef BOOST_MP_STRESS_IO
				  max_bits += 32;
#else
				  max_bits *= 2;
#endif
				  shift = (int)cpp_bin_float<Digits, DigitBase, Allocator, Exponent, MinE, MaxE>::bit_count - exponent() - 1 - power10;
				  continue;
			   }
			}
		 }
		 else
		 {
			//
			// Our integer is bits() * 2^-shift * 10^power10
			//
			if(power10 > 0)
			{
			   if(power10)
			   {
				  cpp_int t;
				  calc_exp = boost::multiprecision::cpp_bf_io_detail::restricted_pow(t, cpp_int(5), power10, max_bits, error);
				  calc_exp += boost::multiprecision::cpp_bf_io_detail::restricted_multiply(i, i, t, max_bits, error);
				  shift -= calc_exp;
			   }
			   if((shift < 0) || ((shift == 0) && error))
			   {
				  // We only get here if we were asked for a crazy number of decimal digits -
				  // more than are present in a 2^max_bits number.
#ifdef BOOST_MP_STRESS_IO
				  max_bits += 32;
#else
				  max_bits *= 2;
#endif
				  shift = (int)cpp_bin_float<Digits, DigitBase, Allocator, Exponent, MinE, MaxE>::bit_count - exponent() - 1 - power10;
				  continue;
			   }
			   if(shift)
			   {
				  roundup = boost::multiprecision::cpp_bf_io_detail::get_round_mode(i, shift - 1, error);
				  if(roundup < 0)
				  {
#ifdef BOOST_MP_STRESS_IO
					 max_bits += 32;
#else
					 max_bits *= 2;
#endif
					 shift = (int)cpp_bin_float<Digits, DigitBase, Allocator, Exponent, MinE, MaxE>::bit_count - exponent() - 1 - power10;
					 continue;
				  }
				  i >>= shift;
			   }
			}
			else
			{
			   // We're right shifting, *and* dividing by 5^-power10,
			   // so 5^-power10 can never be that large or we'd simply
			   // get zero as a result, and that case is already handled above:
			   cpp_int r;
			   BOOST_ASSERT(-power10 < INT_MAX);
			   cpp_int d = pow(cpp_int(5), static_cast<unsigned>(-power10));
			   d <<= shift;
			   divide_qr(i, d, i, r);
			   r <<= 1;
			   int c = r.compare(d);
			   roundup = c < 0 ? 0 : c == 0 ? 1 : 2;
			}
		 }
		 s = i.str(0, std::ios_base::fmtflags(0));
		 //
		 // Check if we got the right number of digits, this
		 // is really a test of whether we calculated the
		 // decimal exponent correctly:
		 //
		 boost::intmax_t digits_got = i ? static_cast<boost::intmax_t>(s.size()) : 0;
		 if(digits_got != digits_wanted)
		 {
			base10_exp += digits_got - digits_wanted;
			if(fixed)
			   digits_wanted = digits_got;  // strange but true.
			power10 = digits_wanted - base10_exp - 1;
			shift = (int)cpp_bin_float<Digits, DigitBase, Allocator, Exponent, MinE, MaxE>::bit_count - exponent() - 1 - power10;
			if(fixed)
			   break;
			roundup = 0;
		 }
		 else
			break;
	  }
	  while(true);
	  //
	  // Check whether we need to round up: note that we could equally round up
	  // the integer /i/ above, but since we need to perform the rounding *after*
	  // the conversion to a string and the digit count check, we might as well
	  // do it here:
	  //
	  if((roundup == 2) || ((roundup == 1) && ((s[s.size() - 1] - '0') & 1)))
	  {
		 boost::multiprecision::detail::round_string_up_at(s, static_cast<int>(s.size() - 1), base10_exp);
	  }

	  if(sign())
		 s.insert(static_cast<std::string::size_type>(0), 1, '-');

	  boost::multiprecision::detail::format_float_string(s, base10_exp, dig, f, false);
   }
   else
   {
	  switch(exponent())
	  {
	  case exponent_zero:
		 s = sign() ? "-0" : f & std::ios_base::showpos ? "+0" : "0";
		 boost::multiprecision::detail::format_float_string(s, 0, dig, f, true);
		 break;
	  case exponent_nan:
		 s = "nan";
		 break;
	  case exponent_infinity:
		 s = sign() ? "-inf" : f & std::ios_base::showpos ? "+inf" : "inf";
		 break;
	  }
   }
   return s;
}
 */
