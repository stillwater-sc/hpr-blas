#pragma once
// ones.hpp: Ones matrix
//
//
// Copyright (C) 2017-2020 Stillwater Supercomputing, Inc.
// Author: James Quinlan
//
// This file is part of the HPRBLAS project, which is released under an MIT Open Source license.

namespace sw { namespace hprblas { namespace matpak {

template<typename Matrix>
Matrix Ones(size_t m, size_t n) {
    using value_type = typename Matrix::value_type;

    Matrix C(m,n);
    for (size_t i=0; i < m; ++i){
        for (size_t j = 0; j <n; ++j){
            C[i][j] = value_type(1);
        }
    }
    return C;
}
}}}