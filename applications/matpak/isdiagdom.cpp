// isdiagdom.cpp : Is the matrix A diagonally domainant
//
// Copyright (C) 2017-2020 Stillwater Supercomputing, Inc.
// Author: James Quinlan
//
// This file is part of the HPRBLAS project, which is released under an MIT Open Source license.

// COMMON LIBRARIES
#include <iostream>
#include <hprblas>

// Selects posits or floats
#define USE_POSIT 1

int main ()
{
	using namespace std;
	// using namespace sw::hprblas::matpak;
	cout << setprecision(5);


#if USE_POSIT
    constexpr size_t nbits = 16;
	constexpr size_t es = 1;
	using Scalar = posit<nbits, es>;
	cout << "\n\nUsing POSIT<" << nbits << "," <<  es << ">\n" <<  endl;
#else	  
	using Scalar = double;
#endif

	Matrix A = {
				{1, 2, 3, 4},
				{5, 6, 7, 8},
				{8, 7, 6, 5},
				{4, 3, 2, 1}
			};
	cout <<  A << endl;

return 0;
}