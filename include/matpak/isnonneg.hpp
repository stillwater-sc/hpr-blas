#pragma once
// isnonneg.hpp : Deteremines if A is nonnegative
//
// Copyright (C) 2017-2020 Stillwater Supercomputing, Inc.
// Author: James Quinlan
//
// This file is part of the HPRBLAS project, which is released under an MIT Open Source license.

namespace sw { namespace hprblas { namespace matpak {

template<typename Matrix>
bool isnonneg(const Matrix&A) {
	typedef typename Matrix::value_type value_type;
	typedef typename Matrix::size_type  size_type;

	size_type m = num_rows(A);
	size_type n = num_cols(A);

    for(size_t i = 0; i< m; ++i){
        for(size_t j = 0; j< n; ++j){
            if(A[i][j]<0) return false;
        }
    }
    return true;
}   
}}}