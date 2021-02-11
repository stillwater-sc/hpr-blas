// diam.cpp : 
//  
// Copyright (C) 2017-2020 Stillwater Supercomputing, Inc.
// Author: James Quinlan
//
// This file is part of the HPRBLAS project, which is released under an MIT Open Source license.


// Calculate the diameter of a matrix.
// Functions needed: max, norm, mtl::iall=:
#include <iostream>
#include <hprblas>
#include <matpak/ek.hpp>

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
	{
		using Scalar = posit<16,1>;
		// using Scalar = double;
		auto v = ek<Scalar>(3,8);
		std::cout <<  v << std::endl;
 	}

	return 0;
}
