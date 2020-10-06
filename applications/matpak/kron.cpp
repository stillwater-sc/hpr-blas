// kron.cpp : Kronecker Product of two matrices
//
// Copyright (C) 2017-2020 Stillwater Supercomputing, Inc.
// Author: James Quinlan
//
// This file is part of the HPRBLAS project, which is released under an MIT Open Source license.


/* From MATLAB: 
kron(X,Y) is the Kronecker tensor product of X and Y.
    The result is a large matrix formed by taking all possible
    products between the elements of X and those of Y. For
    example, if X is 2 by 3, then kron(X,Y) is
 
       [ X(1,1)*Y  X(1,2)*Y  X(1,3)*Y
         X(2,1)*Y  X(2,2)*Y  X(2,3)*Y ]
 
    If either X or Y is sparse, only nonzero elements are multiplied
    in the computation, and the result is sparse.

    A \in M(p,q) and B \in M(m,n), then kron(A,B) \in M(mxp, nxq)

*/


// COMMON LIBRARIES
#include <iostream>
#include <hprblas>
#include <matpak/kron.hpp>

// Selects posits or floats
#define USE_POSIT 0

int main ()
{
	using namespace std;
	using namespace mtl;
	using namespace sw::unum;
	using namespace sw::hprblas;
	using namespace sw::hprblas::matpak;
	cout << setprecision(5);

 #if USE_POSIT
    constexpr size_t nbits = 16;
	constexpr size_t es = 1;
	using Scalar = posit<nbits, es>;
	cout << "\nUsing POSIT<" << nbits << "," <<  es << ">\n" <<  endl;
#else	  
	using Scalar = double;
#endif

		using Matrix = mtl::mat::dense2D< Scalar >;
		Matrix A = {
				{1, 2, 3},
				{4, 5, 6}
			}; 
        Matrix B = {
				{1, 2, 3},
				{4, 5, 6},
				{7, 8, 9}
			};
	std::cout << kron(A,B) << std::endl;

	return 0;
}

