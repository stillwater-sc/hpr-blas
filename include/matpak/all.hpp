#pragma once
// all.hpp : all elements satisfy criteria?
//  
//
// Copyright (C) 2017-2020 Stillwater Supercomputing, Inc.
// Author: James Quinlan
//
// This file is part of the HPRBLAS project, which is released under an MIT Open Source license.

#include <matpak/relop.hpp>

namespace sw { namespace hprblas { namespace matpak {

template<typename Matrix>
bool all(const Matrix&A) {
    auto m = num_rows(A);
    auto n = num_cols(A);

    for(size_type i = 0; i < m; ++i){
	   for(size_type j = 0; j < n; ++j){
	      	if(A[i][j]) > 0){
	            return false;
	   		}
	   }
    }
    return true;
}
}}}