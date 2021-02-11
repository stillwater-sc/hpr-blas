// sum.cpp : sum of elements in matrix along dimensions specified
//            S = sum(A,dim)
// Copyright (C) 2017-2020 Stillwater Supercomputing, Inc.
// Author: James Quinlan
//
// This file is part of the HPRBLAS project, which is released under an MIT Open Source license.

// COMMON LIBRARIES
#include <iostream>
#include <hprblas>

// DEPENDENCIES
#include <matpak/rowsto.hpp>
#include <matpak/sum.hpp>

// Selects posits or floats
#define USE_POSIT 1


int main ()
{
    // COMMON NAMESPACES
	using namespace std;
	using namespace mtl;
	using namespace sw::universal;
	using namespace sw::hprblas;
	using namespace sw::hprblas::matpak;
	
    cout << setprecision(5);

#if USE_POSIT
    	constexpr size_t nbits = 32;
		constexpr size_t es = 2;
		using Scalar = posit<nbits, es>;
		cout << "\n\n Using POSIT<" << nbits << "," <<  es << ">\n" <<  endl;
#else	  
		using Scalar = double;
#endif
		
		using Matrix = mtl::mat::dense2D< Scalar >;
		Matrix A = rowsto< Matrix >(5,5);   //
		
		cout <<  A << endl;
		cout << "Sum of all elements\n " << sum(A,0) << endl;
		cout <<  "Column sums = \n" << sum(A,1) << endl;
        cout <<  "Row sums = \n" << sum(A,2) << endl;
        // cout << "Tensor, Layer sums\n " << sum(A,3) << endl;
	 
	return 0;
}
