// fliptranspose.cpp : Reflect Matrix across counter-identity
// 		used to create Counter-identity matrix   
//
// Copyright (C) 2017-2020 Stillwater Supercomputing, Inc.
// Author: James Quinlan
//
// This file is part of the HPRBLAS project, which is released under an MIT Open Source license.

 
#include <iostream>
#include <hprblas>
#include <matpak/rowsto.hpp>
#include <matpak/fliptranspose.hpp>

// Selects posits or floats
#define USE_POSIT 1

int main ()
{
	using namespace std;
	using namespace mtl;
	using namespace sw::unum;
	using namespace sw::hprblas;
	using namespace sw::hprblas::matpak;
	cout << setprecision(8);
	{
		using Scalar = posit<16,1>;
		// using Scalar = double;

		using Matrix = mtl::mat::dense2D< Scalar >;
		Matrix A = rowsto< Matrix >(5,5);   //
		std::cout <<  A << std::endl;
		std::cout <<  fliptranspose(A) << std::endl;
 	}

	return 0;
}
