#pragma once
// bands.hpp : Returns matrix mask M of diagonal bands
//
// Copyright (C) 2017-2020 Stillwater Supercomputing, Inc.
// Author: James Quinlan
//
// This file is part of the HPRBLAS project, which is released under an MIT Open Source license.

namespace sw { namespace hprblas { namespace matpak {

template<typename Matrix>
Matrix bands(const Matrix&A, const Matrix&v = {0}) {
    auto m = num_rows(A);
    auto n = num_cols(A);
    auto N = num_rows(v)*num_cols(v); 

    Matrix M(m,n); M = 0;

    for(size_t k = 0; k<= N; ++k){
    
        if(v[0][k] <= 0){ // Sub Diagonals (+Diagonal)
            for( size_t i = abs(v[0][k]);  i <= m ; ++i) {
                size_t j = i + v[0][k];
                M[i][j] = 1;
            }  
        }else{ // Super Diagonals
            for( size_t j = v[0][k];  j <= n ; ++j) {
                size_t i = j - v[0][k];
                M[i][j] = 1;
            }
        }
    }  
    return M;
}
}}}