// hankel.cpp : Test Hankel matrix (test hankel.hpp)
//
// Copyright (C) 2017-2020 Stillwater Supercomputing, Inc.
// Author: James Quinlan
//
// This file is part of the HPRBLAS project, which is released under an MIT Open Source license.

#include <iostream>
#include <hprblas>
#include <matpak/hankel.hpp>


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
		using Vector = mtl::vec::dense_vector< Scalar >;
		Vector c{1,2, 3, 4};
		Vector r{1, 5, 6, 7, 8};
		auto A = hankel(c,r);
		std::cout <<  A << std::endl;
 	}


	return 0;
}
