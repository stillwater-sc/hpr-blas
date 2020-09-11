#pragma once
// reshape.hpp : Maps m x n --> r x c where mn = rc.
//				 reshape is row-wise.
//
// Copyright (C) 2017-2020 Stillwater Supercomputing, Inc.
// Author: James Quinlan
//
// This file is part of the HPRBLAS project, which is released under an MIT Open Source license.
#include <utility>

namespace sw { namespace hprblas { namespace matpak {

// NOTE: Matrix arg is not marked const bc the copy constructor with matrix elements
// is not const. (See Line 16)

template<typename Matrix>
Matrix reshape(Matrix&A, size_t m, size_t n) {
	typedef typename Matrix::value_type value_type;
	typedef typename Matrix::size_type  size_type;

	size_type ra = num_rows(A);
	size_type ca = num_cols(A);

	if (ra*ca!=m*n){
		return Matrix {0};
	}
	Matrix B(m, n, A.elements());
	return B;
}
}}} 
