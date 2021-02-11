// checkerboard.cpp: Returns +1 / -1 checkerboard pattern matrix 
//
// Copyright (C) 2017-2020 Stillwater Supercomputing, Inc.
// Author: James Quinlan
//
// This file is part of the HPRBLAS project, which is released under an MIT Open Source license.

// COMMON LIBRARIES
#include <iostream>
#include <hprblas>
#include <matpak/checkerboard.hpp>


// Selects posits or floats
#define USE_POSIT 0

int main ()
{
	using namespace std;
	using namespace mtl;
	using namespace sw::universal;
	using namespace sw::hprblas;
	using namespace sw::hprblas::matpak;

 #if USE_POSIT
    constexpr size_t nbits = 16;
	constexpr size_t es = 1;
	using Scalar = posit<nbits, es>;
	cout << "\nUsing POSIT<" << nbits << "," <<  es << ">\n" <<  endl;
#else	  
	using Scalar = double;
#endif
	using Matrix = mtl::mat::dense2D<Scalar>;
	std::cout << checkerboard<Matrix>(8) << std::endl;

	return 0;
}

