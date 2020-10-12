#pragma once
// isequal.hpp : A == B (determines if two matrices are equal)
//
// Copyright (C) 2017-2020 Stillwater Supercomputing, Inc.
// Author: James Quinlan
//
// This file is part of the HPRBLAS project, which is released under an MIT Open Source license.

namespace sw { namespace hprblas { namespace matpak {

template<typename Matrix>
bool isequal(const Matrix&A, const Matrix&B, const typename Matrix::value_type& tolerance = 0.001) {
	typedef typename Matrix::size_type  size_type;

	size_type ra = num_rows(A);
	size_type ca = num_cols(A);

	size_type rb = num_rows(B);
	size_type cb = num_cols(B);

	if (ra!=rb || ca!=cb){
		return false;
	}

	for(size_type i = 0; i < ra; ++i){
	   for(size_type j = 0; j < ca; ++j){
	      	if(abs(A[i][j]-B[i][j]) > tolerance) return false;
	   }
    }
    return true;
}
}}}
