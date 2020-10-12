#pragma once
// hadamard.hpp : Hadamard product of A.*B.  
//
// Copyright (C) 2017-2020 Stillwater Supercomputing, Inc.
// Author: James Quinlan
//
// This file is part of the HPRBLAS project, which is released under an MIT Open Source license.

namespace sw { namespace hprblas { namespace matpak {

template<typename Matrix>
Matrix hadamard(const Matrix&A, const Matrix&B) {
	typedef typename Matrix::value_type value_type;
	typedef typename Matrix::size_type  size_type;

	size_type m = num_rows(A);
	size_type n = num_cols(A);

    if( m!=num_rows(B) || n != num_cols(B) ){return Matrix {0};}
   
    Matrix C(m,n); 
    int k = 1;
         
    for(int i = 0; i <= num_rows(A); ++i){
        for(size_t j = 0; j <= num_cols(A); ++j){
            C[i][j] = A[i][j]*B[i][j];    
        }
    } 
    return C;
}
}}}