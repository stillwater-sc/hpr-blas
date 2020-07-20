#pragma once
// flipud.hpp : Flip matrix in up/down direction.
//
// Copyright (C) 2017-2020 Stillwater Supercomputing, Inc.
// Author: James Quinlan
//
// This file is part of the HPRBLAS project, which is released under an MIT Open Source license.

namespace sw { namespace hprblas { namespace matpak {

template<typename Matrix>
Matrix flipud(const Matrix&A) {
    typedef typename Matrix::value_type value_type;
	typedef typename Matrix::size_type size_type;

    size_type r = num_rows(A);
	size_type c = num_cols(A);

    Matrix B(r,c);

    for(size_t i = 0; i < r; ++i){
        B[i][mtl::iall] = A[r - i - 1][mtl::iall];
    }
    return B;
}
}}} // namespace sw::hprblas::matpak
