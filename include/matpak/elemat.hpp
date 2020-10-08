#pragma once
// elemat.hpp : Returns an Elementary Matrix
//              n = rank of matrix             
//              k = scalar (multiplier)
//              ri = row to be effected
//              rj = row to be scaled by k and added to ri
//
// Copyright (C) 2017-2020 Stillwater Supercomputing, Inc.
// Author: James Quinlan
//
// This file is part of the HPRBLAS project, which is released under an MIT Open Source license.

namespace sw { namespace hprblas { namespace matpak {

template<typename Matrix>
Matrix elemat(const size_t&n, const size_t&k, const size_t&ri, const size_t&rj) {
	typedef typename Matrix::value_type value_type;
	typedef typename Matrix::size_type  size_type;

    Matrix A(n,n);
    A = 1;
    A[ri][rj] = k;
    return A;
}
}}}