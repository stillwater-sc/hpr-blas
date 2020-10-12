// hadamard.cpp : Hadamard product of A.*B
//
// Copyright (C) 2017-2020 Stillwater Supercomputing, Inc.
// Author: James Quinlan
//
// This file is part of the HPRBLAS project, which is released under an MIT Open Source license.

// COMMON LIBRARIES
#include <iostream>
#include <hprblas>

// DEPENDENCIES
#include <matpak/hadamard.hpp>

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

#if USE_POSIT
    	constexpr size_t nbits = 32;
		constexpr size_t es = 2;
		using Scalar = posit<nbits, es>;
		cout << "\nUsing POSIT<" << nbits << "," <<  es << ">\n" <<  endl;
#else	  
		using Scalar = double;
#endif

		using Matrix = mtl::mat::dense2D< Scalar >;

		Matrix A = {
				{10, 2, 3, 4},
				{5, 7, 6, 2},
				{8, 7, 30, 5},
				{1, 1, 7, 41}
			};
        
        Matrix B = {
				{1, 0, 2, 4},
				{5, 2, 1, 9},
				{8, 7, 3, 5},
				{4, 3, 2, 1}
			};

    Matrix C = hadamard(A,B);
	cout << "A.*B =  \n" << A << "\n * \n" << B << "\n = \n" << C << "\n\n" << endl;

return 0;
}