// gt.cpp : Greater than, e.g., A>B
//
// Copyright (C) 2017-2020 Stillwater Supercomputing, Inc.
// Author: James Quinlan
//
// This file is part of the HPRBLAS project, which is released under an MIT Open Source license.


#include <iostream>
#include <hprblas>
#include <matpak/rowsto.hpp>
#include <matpak/gt.hpp>

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

	{
		using Scalar = double;
		using Matrix = mtl::mat::dense2D< Scalar >;
		Matrix A = rowsto< Matrix >(5,5);   //
		Matrix B = rowsto< Matrix >(5,5);   //

		// B = A;

		std::cout <<  gt(A, B) << std::endl;
		auto C = A > B;
		std::cout <<  C << std::endl;
	}	

	
	return 0;
}
