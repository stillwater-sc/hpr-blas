#pragma once
// fliptranspose.hpp : Flip and tranpose matrix A
//
// Copyright (C) 2017-2020 Stillwater Supercomputing, Inc.
// Author: James Quinlan
//
// This file is part of the HPRBLAS project, which is released under an MIT Open Source license.

#include <matpak/rot90.hpp>

namespace sw { namespace hprblas { namespace matpak {

template<typename Matrix>
Matrix fliptranspose(const Matrix&A) {
    auto m = num_rows(A);
    auto n = num_cols(A);

    if (m!=n) return Matrix {0};

    Matrix I(n,n);
    I = 1;  

    Matrix J = rot90(I); // Counter-identity
    Matrix B(m,n);
    Matrix T(mtl::mat::trans(A));
    B = J*T*J;

    return B;
}
}}}