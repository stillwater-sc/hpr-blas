// ism.cpp: test to validate that a matrix is an M matrix
//
// Copyright (C) 2017-2020 Stillwater Supercomputing, Inc.
//
// This file is part of the HPR-BLAS project, which is released under an MIT Open Source license.
#include <iostream>
#define MTL_WITH_INITLIST
//#define MTL_WITH_AUTO
//#define MTL_WITH_RANGEDFOR
#include <hprblas>
#include <matpak/ism.hpp>

// Selects posits or floats
#define USE_POSIT 1

int main(int argc, char** argv)
try {
	using namespace std;
	using namespace mtl;
	using namespace sw::unum;
	using namespace sw::hprblas;

	using namespace mtl;

	{
		typedef std::complex<double>  cdouble;
		dense_vector<cdouble>         u(10), v(10);
		dense_vector<double>          w(10), x(10, 4.0);

		for (unsigned i = 0; i < size(v); i++)
			v[i] = cdouble(i + 1, 10 - i), w[i] = 2 * i + 2;

		u = v + w + x;
		cout << "u is " << u << "\n";

		u -= 3 * w;
		cout << "u is " << u << "\n";

		u *= 6;
		cout << "u is " << u << "\n";

		u /= 2;
		cout << "u is " << u << "\n";

		u += dot(v, w) * w + 4.0 * v + 2 * w;
		cout << "u is " << u << "\n";

		cout << "i * w is " << cdouble(0, 1) * w << "\n";

		cout << "element sum of x = " << sum(size(x), x, 1) << '\n';
	}

	/*
	u is {10C}[(7,10),(10,9),(13,8),(16,7),(19,6),(22,5),(25,4),(28,3),(31,2),(34,1)]
	u is {10C}[(1,10),(-2,9),(-5,8),(-8,7),(-11,6),(-14,5),(-17,4),(-20,3),(-23,2),(-26,1)]
	u is {10C}[(6,60),(-12,54),(-30,48),(-48,42),(-66,36),(-84,30),(-102,24),(-120,18),(-138,12),(-156,6)]
	u is {10C}[(3,30),(-6,27),(-15,24),(-24,21),(-33,18),(-42,15),(-51,12),(-60,9),(-69,6),(-78,3)]
	u is {10C}[(1551,-810),(3090,-1697),(4629,-2584),(6168,-3471),(7707,-4358),(9246,-5245),(10785,-6132),(12324,-7019),(13863,-7906),(15402,-8793)]
	i * w is {10C}[(0,2),(0,4),(0,6),(0,8),(0,10),(0,12),(0,14),(0,16),(0,18),(0,20)]
    element sum of x = 40
	*/

	{
		using namespace mtl;
		
		vec::dense_vector<int> v{ 2, 3, 4 };
		mat::dense2D<int> A{ {1, 2}, {3, 4} };

		cout << "sum(v) = " << sum(v) << "\n";

		cout << "sum first row = " << sum(A[0][iall]) << "\n";
		cout << "sum first column = " << sum(A[iall][0]) << "\n";

	}
	/*
	sum(v) = 9
	sum first row = 3
	sum first column = 4
	*/
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
