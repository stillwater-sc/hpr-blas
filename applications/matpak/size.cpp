// size.cpp : Calculates the size (or dimensions) of a Matrix
//            Example: size(A)
//  
// Copyright (C) 2017-2020 Stillwater Supercomputing, Inc.
// Author: James Quinlan
//
// This file is part of the HPRBLAS project, which is released under an MIT Open Source license.

// COMMON LIBRARIES
#include <iostream>
#include <hprblas>

// #include <boost/numeric/mtl/mtl.hpp>  // not needed for size, using for other tests.  remove in production

// DEPENDENCIES
#include <matpak/rowsto.hpp>
#include <matpak/size.hpp>
#include <boost/numeric/mtl/mtl.hpp>

#include <generators/matrix_generators.hpp>

// Selects posits or floats
#define USE_POSIT 0

int main ()
{
    // COMMON NAMESPACES
	using namespace std;
	using namespace mtl; using mtl::iall;
	using namespace sw::universal;
	using namespace sw::hprblas;
	using namespace sw::hprblas::matpak;
	
    cout << setprecision(5);	

 


#if USE_POSIT
    	constexpr size_t nbits = 32;
		constexpr size_t es = 2;
		using Scalar = posit<nbits, es>;
		using Matrix = mtl::mat::dense2D< Scalar >;
		cout << "\nUsing POSIT<" << nbits << "," <<  es << ">\n" <<  endl;
#else	  
		using Scalar = double;
		using Matrix = mtl::mat::dense2D< Scalar >;
#endif


		Matrix A = rowsto< Matrix >(7,7);   //
		cout << "Matrix A = \n" << A << endl;
		cout <<  "Size A = " << size(A) << endl;

	//	Matrix B = uniform_rand<Matrix>(6,6);
	//	cout << "Matrix B = \n" << B << endl;


	 //submatrix from matrix per irange
    using mtl::irange;
    irange row(2, 4), col(1, 7);
    dense2D<double> B1= A[row][col];
	 std::cout << "B1 is\n" << B1 << "\n";


	return 0;
}
