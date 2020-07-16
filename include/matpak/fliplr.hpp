/*
// fliplr.hpp : Flip matix in left/right direction.
//
// Copyright (C) 2017-2020 Stillwater Supercomputing, Inc.
// Author: James Quinlan
//
// This file is part of the HPRBLAS project, which is released under an MIT Open Source license.
*/

namespace sw { namespace hprblas { namespace matpak {

template<typename Matrix>
Matrix  fliplr(const Matrix&A) {
    typedef typename Matrix::value_type value_type;

    unsigned r = num_rows(A);
    unsigned c = num_cols(A);

    Matrix B(r,c);

    for(size_t i = 0; i < c; ++i){
        B[mtl::iall][i] = A[mtl::iall][c - i - 1];
    }
    return B;
}
}}} // namespace sw::hprblas::matpak
