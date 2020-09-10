#pragma once
// isdiagdom.hpp : Determine if matrix is diagonally dominant
//
// Copyright (C) 2017-2020 Stillwater Supercomputing, Inc.
// Author: James Quinlan
//
// This file is part of the HPRBLAS project, which is released under an MIT Open Source license.

namespace sw { namespace hprblas { namespace matpak {

template<typename Matrix>
bool isdiagdom(const Matrix&A) {
    auto m = num_rows(A);
    auto n = num_cols(A);

    if (m!=n) return false;

    int R = 0;

    for(int i=0; i<n; ++i){
        for(int j=0; j<n; ++j){
            if(j!=i){
                R = R+A[i][j];
            }
            if(A[i][i]<R) return false;
            R = 0; 
        }
    }
    return true;
}
}}}