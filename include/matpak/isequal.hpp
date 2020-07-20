#pragma once
// isequal.hpp : A == B (are two matrices equal)
//
// Copyright (C) 2017-2020 Stillwater Supercomputing, Inc.
// Author: James Quinlan
//
// This file is part of the HPRBLAS project, which is released under an MIT Open Source license.

namespace sw { namespace hprblas { namespace matpak {

template<typename Matrix>
bool isequal(const Matrix&A, const Matrix&B) {
	typedef typename Matrix::value_type value_type;
	typedef typename Matrix::size_type  size_type;

	size_type ra = num_rows(A);
	size_type ca = num_cols(A);

	size_type rb = num_rows(B);
	size_type cb = num_cols(B);

    if (ra!=rb || ca!=cb){
        return false;
    }else{
        for(size_t i = 0; i < ra; ++i){
            for(size_t j = 0; j < ca; ++j){
                if(A[i][j]!=B[i][j]){
                    // std::cout << A[i][j] <<  " " << B[i][j] <<std::endl;
                    return false;
                }
            }
        }
    }
    return true;
}
}}} // namespace sw::hprblas::matpak
