#pragma once
// sum.hpp :  sum of elements in matrix along dimensions specified
//            S = sum(A,dim)
//  
// Copyright (C) 2017-2020 Stillwater Supercomputing, Inc.
// Author: James Quinlan
//
// This file is part of the HPRBLAS project, which is released under an MIT Open Source license.

namespace sw { namespace hprblas { namespace matpak {

template<typename Matrix>
Matrix sum(const Matrix&A, const int& dim) {
	typedef typename Matrix::size_type  size_type;
    typedef typename Matrix::value_type value_type;

	size_type m = num_rows(A);
	size_type n = num_cols(A);
  

    if (dim == 1){ // Column sum
        Matrix S(1,n); 
        for(size_type j = 0; j < n; ++j){
            value_type temp_sum = 0;
	        for(size_type i = 0; i < m; ++i){
                temp_sum += A[i][j];
	        }
            S[0][j] = temp_sum;

        }
        return S;
    }else if(dim == 2){ // Row sums
        Matrix S(m,1); 
        for(size_type i = 0; i < m; ++i){
            value_type temp_sum = 0;
	        for(size_type j = 0; j < n; ++j){
                temp_sum += A[i][j];
	        }
            S[i][0] = temp_sum;
        }
        return S;
        }else if(dim == 0){ // Sum All elements
        Matrix S(1,1); 
        value_type temp_sum = 0;
        for(size_type i = 0; i < m; ++i){
	        for(size_type j = 0; j < n; ++j){
                temp_sum += A[i][j];
	        }
        }
        S[0][0] = temp_sum;
        return S;
    }
    return Matrix{};
}
}}}