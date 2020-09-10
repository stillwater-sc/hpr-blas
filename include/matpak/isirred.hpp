#pragma once
// isirred.hpp : Deteremines if A is an irreducible matrix
//
// Copyright (C) 2017-2020 Stillwater Supercomputing, Inc.
// Author: James Quinlan
//
// This file is part of the HPRBLAS project, which is released under an MIT Open Source license.

#include<isnonneg.hpp>
namespace sw { namespace hprblas { namespace matpak {

template<typename Matrix>
bool isirred(const Matrix&A) {
	typedef typename Matrix::value_type value_type;
	typedef typename Matrix::size_type  size_type;

	size_type m = num_rows(A);
	size_type n = num_cols(A);

    // Must be nonnegative square matrix
    if (m!=n){return false};
    if (!isnonneg(A)){return false};

    Matrix I(n,n);
    I = 1; 
    
    Matrix B(n,n);
    B = (I+A)^(n-1);

    if(isnonneg(B)) return true;
}   
}}}