// rowsto.cpp : Generates row stochastic matrix
//				
// Copyright (C) 2017-2020 Stillwater Supercomputing, Inc.
// Author: James Quinlan
//
// This file is part of the HPRBLAS project, which is released under an MIT Open Source license.

#include <iostream>
#include <hprblas>
#include <matpak/rowsto.hpp>
//#include <boost/numeric/mtl/mtl.hpp>

// Selects posits or floats
#define USE_POSIT 1


int main ()
{
	using namespace std;
	using namespace mtl;
	using namespace sw::universal;
	using namespace sw::hprblas;
	cout << setprecision(20);

#if USE_POSIT
    constexpr size_t nbits = 32;
	constexpr size_t es = 2;
	using Scalar = posit<nbits, es>;
	cout << "\nUsing POSIT<" << nbits << "," <<  es << ">\n" <<  endl;
#else	  
	using Scalar = double;
#endif

	using Matrix = mtl::mat::dense2D< Scalar >;
	Matrix A = sw::hprblas::matpak::rowsto< Matrix >(5,5);   //
	std::cout <<  A << std::endl;
	for (size_t i=0;i<5;++i){
		cout << sum(A[i][iall])<<endl;
	}
	return 0;
}
