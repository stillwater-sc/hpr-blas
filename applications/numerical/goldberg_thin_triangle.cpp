// goldberg_thin_triangle.cpp: example program showing the Goldberg thin triangle example
//
// Copyright (C) 2017-2019 Stillwater Supercomputing, Inc.
//
// This file is part of the HPRBLAS project, which is released under an MIT Open Source license.
#define ALIASING_ALLOWED
#include "common.hpp"


/* 
* Based on the discussion of rounding error of Golberg's thin triangle
* 
* We are following the exposition described in http://marc-b-reynolds.github.io/math/2019/02/06/Posit1.html
*/

/*

Introduction:

Goldberg’s long thin triangle
Kahan presented this problem in or prior to 1986 and the Goldberg paper (from 1991) was inspired from attended a conference by Kahan.


We are going totally compute the area A of a thin triangle using the classic form of Heron’s equation.

s = (a+b+c) / 2
A = SQRT(s(s−a)(s−b)(s−c))

The lengths are set to:
a=7
b=c=12(a+3*ulp(a))   ulp is "unit in last position"

exact  = 1000000001001111111001110×2−31
posit  = 1000000001001111111001111×2−31
IEEE   = 1001010000101001011111110×2−31

*/
namespace sw {
	namespace unum {



	
	}
}



template<typename Scalar>
Scalar GoldbergThinTriangle(const Scalar& a, const Scalar& b, const Scalar& c) {
	using namespace std;
	using namespace sw::unum;
	Scalar s, A;

	std::cout << "    a  = " << to_binary(a) << " " << to_base2_scientific(a) << " : " << std::showpos << a << std::noshowpos << std::endl;
	std::cout << "    b  = " << to_binary(b) << " " << to_base2_scientific(b) << " : " << std::showpos << std::setprecision(8) << b << std::noshowpos << std::endl;
	std::cout << "    c  = " << to_binary(c) << " " << to_base2_scientific(c) << " : " << std::showpos << c << std::noshowpos << std::endl;
	std::cout << "ulp(a) = " << to_binary(ulp<Scalar>(a)) << " " << to_base2_scientific(ulp<Scalar>(a)) << " : " << ulp<Scalar>(a) << std::endl;

	s = (a + b + c) / 2;
	std::cout << "    s  = " << to_binary(s) << " " << to_base2_scientific(s) << " : " << std::showpos << s << std::noshowpos << std::endl;

	A = sqrt(s * (s - a)*(s - b)*(s - c));
	std::cout << "    A  = " << to_binary(s) << " " << to_base2_scientific(A) << " : " << std::showpos << A << std::noshowpos << std::endl;

	return A;
}

int main(int argc, char** argv)
try {
	using namespace std;
	using namespace sw::unum;

	/*
	{
		// single precision floats
		cout << "0.125f = " << to_base2_scientific(0.125f) << endl;
		cout << "0.25f  = " << to_base2_scientific(0.25f) << endl;
		cout << "0.5f   = " << to_base2_scientific(0.5f) << endl;
		cout << "1.0f   = " << to_base2_scientific(1.0f) << endl;
		cout << "2.0f   = " << to_base2_scientific(2.0f) << endl;
		cout << "4.0f   = " << to_base2_scientific(4.0f) << endl;
		cout << "8.0f   = " << to_base2_scientific(8.0f) << endl;

		cout << "ulp(1) = " << to_base2_scientific(ulp(1.0f)) << endl;
	}


	{
		constexpr size_t nbits = 32;
		constexpr size_t es = 2;
		cout << "0.125f = " << to_base2_scientific(posit<nbits, es>(0.125f)) << endl;
		cout << "0.25f  = " << to_base2_scientific(posit<nbits, es>(0.25f)) << endl;
		cout << "0.5f   = " << to_base2_scientific(posit<nbits, es>(0.5f)) << endl;
		cout << "1.0f   = " << to_base2_scientific(posit<nbits, es>(1.0f)) << endl;
		cout << "2.0f   = " << to_base2_scientific(posit<nbits, es>(2.0f)) << endl;
		cout << "4.0f   = " << to_base2_scientific(posit<nbits, es>(4.0f)) << endl;
		cout << "8.0f   = " << to_base2_scientific(posit<nbits, es>(8.0f)) << endl;
		cout << "ulp(1) = " << to_base2_scientific(ulp(posit<nbits, es>(1.0f))) << endl;
		cout << "       = " << to_binary(posit<nbits, es>(1.0f)) << endl;

	}
	*/

	// posit<32,2>
	{
		constexpr size_t nbits = 32;
		constexpr size_t es = 2;
		constexpr size_t capacity = 10;
		using Scalar = posit<32, 2>;
		cout << "posit<32, 2>\n";
		Scalar a, b, c, A;
		a = 7;
		b = 0.5 * (a + 3 * ulp<Scalar>(a));
		c = b;
		A = GoldbergThinTriangle<Scalar>(a, b, c);
		cout << "Area = " << A << endl;
		using Quire = quire<nbits, es, capacity>;
		Scalar s = (a + b + c) / 2;
		Quire q1 = quire_mul(s, (s-a));
		Quire q2 = quire_mul((s - b), (s - c));

	}

	// IEEE single precision float
	{
		using Scalar = float;
		cout << "IEEE single precision float\n";
		Scalar a = 7;
		Scalar b = 0.5f * (a + 3 * ulp<Scalar>(a));
		Scalar c = b;
		Scalar A = GoldbergThinTriangle<Scalar>(a, b, c);
		cout << "Area = " << A << endl;
	}

	// IEEE double precision float
	{
		using Scalar = double;
		cout << "IEEE double precision float\n";
		Scalar a = 7;
		Scalar b = 0.5f * (a + 3 * ulp<Scalar>(a));
		Scalar c = b;
		Scalar A = GoldbergThinTriangle<Scalar>(a, b, c);
		cout << "Area = " << A << endl;
	}

	return EXIT_SUCCESS;
}
catch (char const* msg) {
	std::cerr << msg << std::endl;
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