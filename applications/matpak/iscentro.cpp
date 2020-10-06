// iscentro.cpp : Determine if a matrix is centrosymmetric
//
// Copyright (C) 2017-2020 Stillwater Supercomputing, Inc.
// Author: James Quinlan
//
// This file is part of the HPRBLAS project, which is released under an MIT Open Source license.

#include <iostream>
#include <hprblas>
#include <matpak/isa/iscentro.hpp>

// Selects posits or floats
#define USE_POSIT 1

int main ()
{
	using namespace std;
	using namespace mtl;
	using namespace sw::unum;
	using namespace sw::hprblas;
	using namespace sw::hprblas::matpak;
	cout << setprecision(5);

#if USE_POSIT
    constexpr size_t nbits = 16;
	constexpr size_t es = 1;
	using Scalar = posit<nbits, es>;
	cout << "\nUsing POSIT<" << nbits << "," <<  es << ">\n" <<  endl;
#else	  
	using Scalar = double;
#endif

	using Matrix = mtl::mat::dense2D< Scalar >;
	// Matrix A = rowsto< Matrix >(5,5);   //
	Matrix A = {
				{1, 2, 3, 4},
				{5, 6, 7, 8},
				{8, 7, 6, 5},
				{4, 3, 2, 1}
			};
	std::cout <<  A << std::endl;

	int x = 0;
	if(iscentro(A)){
		x = 1;
	}
	std::cout <<  x << std::endl;

	return 0;
}
