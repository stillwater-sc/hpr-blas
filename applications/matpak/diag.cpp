// diag.cpp : 
//  
// Copyright (C) 2017-2020 Stillwater Supercomputing, Inc.
// Author: James Quinlan
//
// This file is part of the HPRBLAS project, which is released under an MIT Open Source license.

// Calculate the diameter of a matrix.
// Functions needed: max, norm, mtl::iall=:
#include <iostream>
#include <hprblas>
#include <matpak/rowsto.hpp>


// Selects posits or floats
#define USE_POSIT 1


/*
Returns a view of a matrix A from diagonal begin to end.

The main diagonal is numbered 0; the off-diagonal below the main one is -1.
Accordingly, the off-diagonal above the main is 1. The parameters begin and end
specify a right-open interval. For, instance bands(A, -1, 2) yields a tridiagonal matrix.

See: http://new.simunova.com/en/mtl4/
*/

int main ()
{
	using namespace std;
	using namespace mtl;
	using namespace sw::universal;
	using namespace sw::hprblas;
	using namespace sw::hprblas::matpak;
	cout << setprecision(5);
	{
		using Scalar = posit<16,1>;
		// using Scalar = double;
		using Matrix = mtl::dense2D<Scalar>;
		using Vector = mtl::dense_vector<Scalar>;
		Matrix A = sw::hprblas::matpak::rowsto< Matrix >(5,5);

		auto B = mtl::mat::bands(A, 0, 1);
		size_t m = num_rows(B);
		Vector d(m);

		for (size_t i=0;i<m;++i){
			d[i] = B[i][i];
		}

#if 0
		Matrix d(m,1);

		for (size_t i=0;i<m;++i){
			d[i,0] = B[i][i];
		}
		/*
		/Users/equinlan/Dropbox/Research/posits/forks/hpr-blas/applications/matpak/diag.cpp:37:6: warning:
      expression result unused [-Wunused-value]
                        d[i,0] = B[i][i];
                          ^
/Users/equinlan/Dropbox/Research/posits/forks/hpr-blas/applications/matpak/diag.cpp:37:11: error:
      no viable overloaded '='
                        d[i,0] = B[i][i];
                        ~~~~~~ ^ ~~~~~~~
/Users/equinlan/Dropbox/Research/posits/forks/mtl4/boost/numeric/mtl/operation/matrix_bracket.hpp:32:12: note:
      candidate function (the implicit copy assignment operator) not viable: no known conversion
      from 'typename boost::enable_if<boost::is_integral<unsigned long>, posit<16, 1> >::type'
      (aka 'sw::universal::posit<16, 1>') to 'const
      mtl::operations::bracket_proxy<mtl::mat::dense2D<sw::universal::posit<16, 1>,
      mtl::mat::parameters<mtl::tag::row_major, mtl::index::c_index, mtl::non_fixed::dimensions,
      false, unsigned long> >, mtl::mat::dense2D<sw::universal::posit<16, 1>,
      mtl::mat::parameters<mtl::tag::row_major, mtl::index::c_index, mtl::non_fixed::dimensions,
      false, unsigned long> > &, sw::universal::posit<16, 1> &>' for 1st argument
    struct bracket_proxy
	*/
#endif


		std::cout <<  d << std::endl;
 	}


	return 0;
}
