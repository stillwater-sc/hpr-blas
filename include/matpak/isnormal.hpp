#pragma once
// isnormal.hpp : Checks if a matrix A is normal
// 	Def: (Normal matrices)  A'A = AA'
//
// Copyright (C) 2017-2020 Stillwater Supercomputing, Inc.
// Author: James Quinlan
//
// This file is part of the HPRBLAS project, which is released under an MIT Open Source license.

#include <matpak/isequal.hpp>

namespace sw { namespace hprblas { namespace matpak {

template<typename Matrix>
bool isnormal(const Matrix&A, const typename Matrix::value_type& tolerance = 0.001) {
     Matrix A_transpose(mtl::mat::trans(A));
     return isequal(A_transpose*A, A*A_transpose, tolerance);
}

}}} // namespace sw::hprblas::matpak
