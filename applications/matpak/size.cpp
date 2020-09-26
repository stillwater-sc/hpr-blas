// size.cpp : Calculates the size (or dimensions) of a Matrix
//            Example: size(A)
//  
// Copyright (C) 2017-2020 Stillwater Supercomputing, Inc.
// Author: James Quinlan
//
// This file is part of the HPRBLAS project, which is released under an MIT Open Source license.

// COMMON LIBRARIES
#include <iostream>
#include <hprblas>

// DEPENDENCIES
#include <matpak/rowsto.hpp>
#include <matpak/size.hpp>

// Selects posits or floats
#define USE_POSIT 0

int main ()
{
    // COMMON NAMESPACES
	using namespace std;
	using namespace mtl;
	using namespace sw::unum;
	using namespace sw::hprblas;
	using namespace sw::hprblas::matpak;
	
    cout << setprecision(5);	

#if USE_POSIT
    	constexpr size_t nbits = 32;
		constexpr size_t es = 2;
		using Scalar = posit<nbits, es>;
		using Matrix = mtl::mat::dense2D< Scalar >;
#else	  
		using Scalar = double;
		using Matrix = mtl::mat::dense2D< Scalar >;
#endif


		Matrix A = rowsto< Matrix >(5,5);   //
		cout << "Matrix A = \n" << A << endl;
		cout <<  "Size A = " << size(A) << endl;

	return 0;
}
