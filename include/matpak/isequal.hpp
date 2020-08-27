#pragma once
// isequal.hpp : A == B (are two matrices equal)
//
// Copyright (C) 2017-2020 Stillwater Supercomputing, Inc.
// Author: James Quinlan
//
// This file is part of the HPRBLAS project, which is released under an MIT Open Source license.
#include <utility>

namespace sw { namespace hprblas { namespace matpak {

template<typename Matrix>
bool isequal(const Matrix&A, const Matrix&B, const typename Matrix::value_type& tolerance = 0.001) {
	typedef typename Matrix::value_type value_type;
	typedef typename Matrix::size_type  size_type;

	size_type ra = num_rows(A);
	size_type ca = num_cols(A);

	size_type rb = num_rows(B);
	size_type cb = num_cols(B);


	if (ra!=rb || ca!=cb){
        return false;
	}

    if (ra!=rb || ca!=cb){
        return false;
    }else{
	    for(size_t i = 0; i < ra; ++i){
	        for(size_t j = 0; j < ca; ++j){
	        	if(abs(A[i][j]-B[i][j]) > tolerance){
	                return false;
	            }else{
					std::cout << A[i][j] << " " << B[i][j] << " ";
					std::cout << abs(A[i][j]-B[i][j]) << std::endl;
				}
	        }
	    }
    }
    return true;
}

/*
template<typename Matrix>
std::pair<Matrix, bool> eq(const Matrix&A, const Matrix&B) {
	typedef typename Matrix::value_type value_type;
	typedef typename Matrix::size_type  size_type;

	size_type ra = num_rows(A);
	size_type ca = num_cols(A);

	size_type rb = num_rows(B);
	size_type cb = num_cols(B);

	Matrix T;


	if (ra!=rb || ca!=cb){
        return std::pair<Matrix, bool>(T,false);
	}

	T.resize(ra,ca);
	T = true;

    for(size_t i = 0; i < ra; ++i){
        for(size_t j = 0; j < ca; ++j){
        	if(abs(A[i][j]-B[i][j]) > value_type(0.000000000000001)){
                T[i][j] = false;
            }
        }
    }
    return std::pair<Matrix, bool>(T,true);
}

*/
}}} // namespace sw::hprblas::matpak
