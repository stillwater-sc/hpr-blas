#pragma once
// checkerboard.hpp: Returns +1 / -1 checkerboard pattern matrix 
//
// Definition: A is checkerboard if sign(a(i,j)) = (-1)^(i+j)
//
// Copyright (C) 2017-2020 Stillwater Supercomputing, Inc.
// Author: James Quinlan
//
// This file is part of the HPRBLAS project, which is released under an MIT Open Source license.
#include<matpak/ones.hpp>

namespace sw { namespace hprblas { namespace matpak {

template<typename Matrix>
Matrix checkerboard(size_t m) {
    using value_type = typename Matrix::value_type;
    // using Vector = mtl::vec::dense_vector<value_type>;

    Matrix C(m,m);   
    for(size_t i = 0; i < m; ++i){
        for(size_t j = 0; j < m; ++j){
            C[i][j]  = pow(-1 , i+j);
        }
    }

    /*
    if ( m % 2 != 0 ){
        n = (m+1)/2;
    }else{
        n = m/2;
    }
    */


    /*
    Vector A(n);
    A = 1;
    Vector B = {1, -1};

    auto C = kron(A,B);
    C = toeplitz(C)
    */
    return C;
}
}}}