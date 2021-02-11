// isdiagdom.cpp : Is the matrix A diagonally domainant
//
// Copyright (C) 2017-2020 Stillwater Supercomputing, Inc.
// Author: James Quinlan
//
// This file is part of the HPRBLAS project, which is released under an MIT Open Source license.

// COMMON LIBRARIES
#include <iostream>
#include <hprblas>

// DEPENDENCIES
#include <matpak/isa/isdiagdom.hpp>

// Selects posits or floats
#define USE_POSIT 1


int main ()
{
	using namespace std;
	using namespace mtl;
	using namespace sw::universal;
	using namespace sw::hprblas;
	using namespace sw::hprblas::matpak;
	cout << setprecision(5);

#if USE_POSIT
    	constexpr size_t nbits = 32;
		constexpr size_t es = 2;
		using Scalar = posit<nbits, es>;
		cout << "\nUsing POSIT<" << nbits << "," <<  es << ">\n" <<  endl;
#else	  
		using Scalar = double;
#endif

		using Matrix = mtl::mat::dense2D< Scalar >;

		Matrix A = {
				{10, 2, 3, 4},
				{5, 26, 7, 8},
				{8, 7, 30, 5},
				{4, 3, 2, 11}
			};
	cout << "The diagonal domainance of A is: " << isdiagdom(A) << "\n\n" << endl;

return 0;
}