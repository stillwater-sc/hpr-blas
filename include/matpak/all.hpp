#pragma once
// all.hpp : all elements satisfy criteria?
//  
// Copyright (C) 2017-2020 Stillwater Supercomputing, Inc.
// Author: James Quinlan
//
// This file is part of the HPRBLAS project, which is released under an MIT Open Source license.


// NOTES:
/* 
1. need to have relational operator as input.  Or I want to 
 call the function as all(A>0).  Additionally, need to check per dimension
 (i.e., rows or columns or 'all').  So, call should be: all(A>0,2) 
 would check if all elements of A are greater than zero, by rows.

 2. Should be able to short-circuit.  Assume true, then if element in row
 is false, make 'all' in that row false...no need to check other row entries.

 3. lambda function - chebfun
[](){ function body}
[local variables]() { function body }
[](x) { return x*x; }

[captures]( params) { function body; }

See https://en.cppreference.com/w/cpp/language/lambda

void any(const Matrix& A, (*f)())
Matrix operator<(const Matrix& A, const Matrix& B)
A < B;

Matrix operator<(Matrix& A, Matrix& B)
Matrix operator+(Matrix& A, Matrix& B)

Matrix operator<(Matrix& A, Matrix& B)
Matrix operator+(Matrix& A, Matrix& B)
// matrix element-wise sum
template<typename Scalar>
matrix<Scalar> operator+(const matrix<Scalar>& A, const matrix<Scalar>& B) {
	matrix<Scalar> Sum(A);
	return Sum += B;
}

// matrix element-wise sum
template<typename Scalar>
matrix<Scalar> operator+(const matrix<Scalar>& A, const matrix<Scalar>& B) {
	matrix<Scalar> Sum(A);
	return Sum += B;
}


sw::unum::matrix class?
// matrix element-wise sum
	matrix& operator+=(const matrix& rhs) {
		// check if the matrices are compatible
		if (_m != rhs._m || _n != rhs._n) {
			std::cerr << "Element-wise matrix sum received incompatible matrices ("
				<< _m << ", " << _n << ") += (" << rhs._m << ", " << rhs._n << ")\n";
			return *this; // return without changing
		}
		for (size_type e = 0; e < _m * _n; ++e) {
			data[e] += rhs.data[e];
		}
		return *this;
	}

 4. Look at:  include/universal/blas/matrix.hpp
    See vector.hpp for class template.



*/

#include<matpak/relop.hpp>

namespace sw { namespace hprblas { namespace matpak {

template<typename Matrix>
Matrix all(const Matrix&A, int dim = 1) {
    auto m = num_rows(A);
    auto n = num_cols(A);

    size_t r = 0;
    size_t c = 0;
    size_t b = 0;

    switch (dim)
    {
        case 1: Matrix B(1,n) = true; c = n - 1; b = c; break;
        case 2: Matrix B(m,1) = true; r = m - 1; b = r; break;
        case 3: Matrix B(1,1) = true; break;
        default: Matrix B(1,1) = true;
    }
    
    for(size_type i = 0; i < m; ++i){
	   for(size_type j = 0; j < n; ++j){
	      	if(A[i][j]) < 0){
	            B[i][j] = false;
            }
	    }
    }
    return B;
}
}}}