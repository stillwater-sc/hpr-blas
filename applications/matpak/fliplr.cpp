// fliplr.cpp : Flip matix in left/right direction.
//
// Copyright (C) 2017-2020 Stillwater Supercomputing, Inc.
// Author: James Quinlan
//
// This file is part of the HPRBLAS project, which is released under an MIT Open Source license.
#include <iostream>
#include <hprblas>
#include <matpak/rowsto.hpp>
#include <matpak/fliplr.hpp>

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

	// ---
 	using Matrix = mtl::mat::dense2D< Scalar >;
	Matrix A = rowsto< Matrix >(5,5);   //
	std::cout <<  A << std::endl;
	auto B = fliplr(A);
	std::cout <<  B << std::endl;
	 
return 0;
}
