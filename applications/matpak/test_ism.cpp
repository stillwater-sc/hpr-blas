//
// Copyright (C) 2017-2019 Stillwater Supercomputing, Inc.
//
// This file is part of the universal numbers project, which is released under an MIT Open Source license.

// enable the mathematical constants in cmath: old-style preprocessor magic which isn't best practice anymore
#include "common.hpp"



// Turn it off for now
#define USE_POSIT

int main(int argc, char** argv)
try {
	const size_t nbits = 16;
	const size_t es = 1;
	const size_t vecSize = 32;

#ifdef USE_POSIT
	using Ty     = sw::unum::posit<8, 0>;
	using Matrix = mtl::dense2D< Ty >;
	using Vector = mtl::dense_vector< Ty >;
#else
	using Ty     = float;
	using Matrix = mtl::dense2D<float>;
	using Vector = mtl::dense_vector<float>;
#endif

  int n = 3; // Number of states
  int m = 1; // Number of measurements

  double dt = 1.0/30; // Time step

  mtl::dense2D<Ty> A(n, n); // System dynamics matrix
    
  A = 1;

  if(ism(A)){
      cout << "A is an M-matrix" << A << '\n';
  }
  else{
      cout << "Failed" << '\n';
  }

    
    

  return EXIT_SUCCESS;
}
catch (char const* msg) {
	std::cerr << msg << std::endl;
	return EXIT_FAILURE;
}
catch (const posit_arithmetic_exception& err) {
	std::cerr << "Uncaught posit arithmetic exception: " << err.what() << std::endl;
	return EXIT_FAILURE;
}
catch (const quire_exception& err) {
	std::cerr << "Uncaught quire exception: " << err.what() << std::endl;
	return EXIT_FAILURE;
}
catch (const posit_internal_exception& err) {
	std::cerr << "Uncaught posit internal exception: " << err.what() << std::endl;
	return EXIT_FAILURE;
}
catch (std::runtime_error& err) {
	std::cerr << err.what() << std::endl;
	return EXIT_FAILURE;
}
catch (...) {
	std::cerr << "Caught unknown exception" << std::endl;
	return EXIT_FAILURE;
}
