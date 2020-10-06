#pragma once
// checkerboard.hpp: Returns +1 / -1 checkerboard pattern matrix 
//
// Definition: A is checkerboard if sign(a(i,j)) = (-1)^(i+j)
//
// Copyright (C) 2017-2020 Stillwater Supercomputing, Inc.
// Author: James Quinlan
//
// This file is part of the HPRBLAS project, which is released under an MIT Open Source license.

namespace sw { namespace hprblas { namespace matpak {

template<typename Scalar>
mtl::mat::dense2D<Scalar> checkerboard(const int&m) {
    int n;
    if ( m % 2 != 0 ){
        n = (m+1)/2;
    }else{
        n = m/2;
    }

    mtl::mat::dense2D<Scalar> C(n,n);

    for(size_t i = 0; i < n; ++i){
        for(size_t j = 0; j < n; ++j){
            C[i][j]  = pow(-1 , i+j);
        }
    }
    return C;
}
}}}