// flipud.cpp : Flip matrix in up/down direction.
//
// Copyright (C) 2017-2020 Stillwater Supercomputing, Inc.
// Author: James Quinlan
//
// This file is part of the HPRBLAS project, which is released under an MIT Open Source license.
#include <iostream>

#include <hprblas>
#include <matpak/rowsto.hpp>
#include <matpak/flipud.hpp>


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
		flipud(A);
		std::cout <<  A << std::endl;
	}

	{
		constexpr size_t nbits = 32;
		constexpr size_t es = 2;

		using Scalar = posit<nbits, es>;
		using Matrix = mtl::mat::dense2D< Scalar >;
		//Matrix A = rowsto< Matrix >(5,5);   //
		//std::cout <<  A << std::endl;
		//fliplr(A);
		//std::cout <<  A << std::endl;
	}

	return 0;
}
