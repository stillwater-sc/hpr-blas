#pragma once
// size.hpp : determines size of matrix in rows and columns [m, n]
//
// Copyright (C) 2017-2020 Stillwater Supercomputing, Inc.
// Author: James Quinlan
//
// This file is part of the HPRBLAS project, which is released under an MIT Open Source license.

namespace sw { namespace hprblas { namespace matpak {

template<typename Matrix>
Matrix size(const Matrix&A) {
	typedef typename Matrix::size_type  size_type;

	size_type m = num_rows(A);
	size_type n = num_cols(A);

	// std::cout << m << n << std::endl;
    Matrix C(1,2);
	C[0][0] = m;
	C[0][1] = n;
}
}}}