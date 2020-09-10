// isdiagdom.cpp : Is the matrix A diagonally domainant
//
// Copyright (C) 2017-2020 Stillwater Supercomputing, Inc.
// Author: James Quinlan
//
// This file is part of the HPRBLAS project, which is released under an MIT Open Source license.

// COMMON LIBRARIES
#include <iostream>
#include <hprblas>


int main ()
{
	using namespace std;
	// using namespace sw::hprblas::matpak;
	cout << setprecision(5);


// Compile options
#define NUMSYS_POSIT 1
#define NUMSYS_FLOAT 0
#define NUMSYS_DOUBLE 0

// #define POSIT

#if NUMSYS_POSIT 
	{
        using namespace sw::unum;
		constexpr size_t nbits = 32;
		constexpr size_t es = 2;

		using Scalar = posit<nbits, es>;
		using Matrix = mtl::mat::dense2D< Scalar >;
		Matrix A = {
					{1, 2, 3, 4},
					{5, 6, 7, 8},
					{8, 7, 6, 5},
					{4, 3, 2, 1}
				};
		cout <<  A << endl;

		cout <<  A << endl;
	}
#endif

#if NUMSYS_DOUBLE
	{
		using Scalar = double;
		using Matrix = mtl::mat::dense2D< Scalar >;
		Matrix A = {
					{1, 2, 3, 4},
					{5, 6, 7, 8},
					{8, 7, 6, 5},
					{4, 3, 2, 1}
				};
		cout <<  A << endl;
		sw::hprblas::matpak::flipud(A);
		cout <<  A << endl;
	}
#endif

	return 0;
}