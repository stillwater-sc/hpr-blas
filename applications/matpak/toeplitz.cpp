// toeplitz.cpp : Generate Toeplitz matrix
//		Example: A = toeplitz(c,r); where c and r are vectors
//
// Copyright (C) 2017-2020 Stillwater Supercomputing, Inc.
// Author: James Quinlan
//
// This file is part of the HPRBLAS project, which is released under an MIT Open Source license.


#include <iostream>
#include <hprblas>
#include <matpak/toeplitz.hpp>

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
		using Scalar = double;
		using Vector = mtl::vec::dense_vector< Scalar >;
		Vector c{1, 2, 3, 4};
		Vector r{7, 5, 6, 7, 8};
		auto A = toeplitz(c,r);
		std::cout <<  A << std::endl;
 	}
	return 0;
}
