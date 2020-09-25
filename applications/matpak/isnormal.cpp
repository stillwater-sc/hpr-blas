// isnormal.cpp : Determines if a matrix is normal (i.e., A'A = AA'?)
//
// Copyright (C) 2017-2020 Stillwater Supercomputing, Inc.
// Author: James Quinlan
//
// This file is part of the HPRBLAS project, which is released under an MIT Open Source license.

#include <iostream>
#include <hprblas>
#include <matpak/rowsto.hpp>
#include <matpak/isnormal.hpp>

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
		std::cout <<  A << std::endl;
		Matrix A_transpose(mtl::mat::trans(A));
		std::cout <<  A_transpose << std::endl;
		std::cout <<  A*A_transpose << std::endl;
		std::cout <<  A_transpose*A << std::endl;

		if(isnormal(A,0.00001)){
			std::cout <<  "Matrix A is normal" << std::endl;
		}else{
			std::cout <<  "Matrix A is NOT normal" << std::endl;
		}
	}

	// -------------------- 
	{
		constexpr size_t nbits = 32;
		constexpr size_t es = 2;

		using Scalar = posit<nbits, es>;
		using Matrix = mtl::mat::dense2D< Scalar >;

	}

	return 0;
}
