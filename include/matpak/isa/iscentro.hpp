#pragma once
// iscentro.hpp : Determine if matrix is centrosymmetric
//
// Copyright (C) 2017-2020 Stillwater Supercomputing, Inc.
// Author: James Quinlan
//
// This file is part of the HPRBLAS project, which is released under an MIT Open Source license.

#include <matpak/rot90.hpp>
#include <matpak/isa/isequal.hpp>

namespace sw { namespace hprblas { namespace matpak {

template<typename Matrix>
bool iscentro(const Matrix&A) {
    auto m = num_rows(A);
    auto n = num_cols(A);

    if (m!=n) return false;

    Matrix I(n,n);
    I = 1; // Identity matrix

    Matrix J = rot90(I); // Counter-identity
    Matrix JAJ(m,n);
    JAJ = J*A*J;

    return isequal(JAJ,A);
}
}}}