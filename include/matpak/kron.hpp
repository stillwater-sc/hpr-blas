#pragma once
// kron.hpp : Kronecker tensor product of A and B
// 
// The result is formed by taking all possible
// products between the elements of A and those of B.
// A in M(m,n) and B in M(p,q), then kron(A,B) in M(mxp, nxq)
//
// Copyright (C) 2017-2020 Stillwater Supercomputing, Inc.
// Author: James Quinlan
//
// This file is part of the HPRBLAS project, which is released under an MIT Open Source license.

namespace sw { namespace hprblas { namespace matpak {

template<typename Matrix>
Matrix kron(const Matrix&A, const Matrix&B) {
	typedef typename Matrix::value_type value_type;
	typedef typename Matrix::size_type  size_type;

	size_type m = num_rows(A);
	size_type n = num_cols(A);
  size_type p = num_rows(B);
	size_type q = num_cols(B);

  Matrix C(m*p,n*q);    
     
    for(size_t i = 0; i < m; ++i){
        for(size_t j = 0; j < n; ++j){
          
          Matrix X(p,q);
          X = A[i][j]*B[mtl::iall][mtl::iall];

          for(size_t r = 0; r < p; ++r){
            for(size_t c = 0; c < q; ++c){
               C[r + p*i][c  + q*j] = X[r][c];
            } // c
          } // r
        } // j
    } // i
     
    return C;
}
}}}