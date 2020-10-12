// ones.cpp : Matrix of all ones
//
// Copyright (C) 2017-2020 Stillwater Supercomputing, Inc.
// Author: James Quinlan
//
// This file is part of the HPRBLAS project, which is released under an MIT Open Source license.

// COMMON LIBRARIES
#include <iostream>
#include <hprblas> // includes mtl

// DEPENDENCIES
 #include <matpak/ones.hpp>

// Selects posits or floats
#define USE_POSIT 1

int main ()
{
    // COMMON NAMESPACES
	using namespace std;
	using namespace mtl;
	using namespace sw::unum;
	using namespace sw::hprblas;
	using namespace sw::hprblas::matpak;
	
    cout << setprecision(5);

    /*
        // --- POSITS ---
        constexpr size_t nbits = 32;
		constexpr size_t es = 2;
		using Scalar = posit<nbits, es>;
		using Matrix = mtl::mat::dense2D< Scalar >;
		Matrix A = rowsto< Matrix >(5,5);   //
		std::cout <<  A << std::endl;
		fliplr(A);
		std::cout <<  A << std::endl;
	 */ 

        // --- FLOATS ---
		using Scalar = double;
		using Matrix = mtl::mat::dense2D< Scalar >;
		auto A = Ones<Matrix>(5,7);   //
		cout <<  A << endl;
	
	return 0;
}
